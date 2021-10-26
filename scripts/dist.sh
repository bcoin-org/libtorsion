#!/bin/sh

# dist.sh - github release for packages
# Copyright (c) 2020, Christopher Jeffrey (MIT License).
# https://github.com/chjj

# Conventions:
#
#   Version tags must be in the format: v[n[.n[.n]]
#
#   CHANGELOG.md must be present (e.g.):
#
#     # Changelog
#
#     ## v0.0.2
#
#     Release notes for 0.0.2.
#
#     ## v0.0.1
#
#     Release notes for 0.0.1.
#
#   .distrc or distrc must be present (e.g.):
#
#     project_name=linux
#     repo_name=linux
#     owner_name=torvalds

set -e

global_name='dist.sh'
global_repo='dist.sh'
global_owner='chjj'
global_token=none
global_auto=no
global_dry=no
global_force=no
global_quiet=no

dist_echo() {
  if test x"$global_quiet" != x'yes'; then
    echo "$@"
  fi
}

dist_warn() {
  dist_echo "WARNING: $@" >& 2
}

dist_error() {
  dist_echo "ERROR: $@" >& 2
}

dist_archive() {
  dist_tag="$1"
  dist_ver=`echo "$dist_tag" | sed -e 's/^v//'`
  dist_name="${global_name}-${dist_ver}"

  dist_echo "Archiving $dist_tag (${dist_name}):" >& 2

  if test x"$global_force" != x'yes'; then
    if test -f "${dist_name}.tar.gz" -a -f "${dist_name}.zip"; then
      dist_echo "  SKIP ${dist_name}.tar.gz" >& 2
      dist_echo "  SKIP ${dist_name}.zip" >& 2
      return 0
    fi
  fi

  rm -f "${dist_name}.tar.gz"
  rm -f "${dist_name}.zip"

  if test x"$global_auto" = x'yes'; then
    if test x"`git ls-files --others -x '*.asc' -x sha256sums.txt`" != x; then
      dist_error 'Repository has untracked files.'
      return 1
    fi

    dist_head=`git rev-parse --abbrev-ref HEAD 2> /dev/null`

    if test x"$dist_head" = x'HEAD'; then
      dist_head=`git tag -l --points-at HEAD 2> /dev/null | sed 1q`
      if test x"$dist_head" = x; then
        dist_head=`git rev-parse HEAD 2> /dev/null`
      fi
    fi

    dist_stash=no
    dist_checkout=no

    if ! git diff --no-ext-diff --quiet --exit-code HEAD; then
      dist_stash=yes
    fi

    if test `git rev-parse "$dist_tag"` != `git rev-parse "$dist_head"`; then
      dist_checkout=yes
    fi

    if test -f ./autogen.sh; then
      dist_autogen='./autogen.sh'
    elif test -f ./configure.ac; then
      dist_autogen='autoreconf -if'
    else
      dist_error 'Cannot find autogen.sh or configure.ac.'
      return 1
    fi

    if test -f ./configure -o -f ./Makefile.in -o -f ./Makefile; then
      dist_error 'Build files already present.'
      return 1
    fi

    if test $dist_stash = yes; then
      dist_echo "  COMMAND git stash" >& 2
      git stash > /dev/null 2>& 1
    fi

    if test $dist_checkout = yes; then
      dist_echo "  COMMAND git checkout $dist_tag" >& 2
      git checkout "$dist_tag" > /dev/null 2>& 1
    fi

    set x "$dist_autogen" './configure' 'make dist-gzip' 'make dist-zip'
    shift

    dist_failed=no

    for dist_cmd in "$@"; do
      dist_echo "  COMMAND $dist_cmd" >& 2
      if ! $dist_cmd > /dev/null 2>& 1; then
        dist_failed=yes
        break
      fi
    done

    git clean -xdf -e '*.tar.gz'     \
                   -e '*.zip'        \
                   -e '*.asc'        \
                   -e sha256sums.txt > /dev/null 2>& 1

    if test $dist_checkout = yes; then
      dist_echo "  COMMAND git checkout $dist_head" >& 2
      git checkout "$dist_head" > /dev/null 2>& 1
    fi

    if test $dist_stash = yes; then
      dist_echo "  COMMAND git stash pop" >& 2
      git stash pop > /dev/null 2>& 1
    fi

    if test $dist_failed = yes; then
      dist_error 'Build failed.'
      return 1
    fi

    if ! test -f "${dist_name}.tar.gz"; then
      mv ${global_name}-*.tar.gz "${dist_name}.tar.gz"
    fi

    if ! test -f "${dist_name}.zip"; then
      mv ${global_name}-*.zip "${dist_name}.zip"
    fi
  else
    git archive -o "${dist_name}.tar.gz" --prefix "${dist_name}/" "$dist_tag"
    git archive -o "${dist_name}.zip" --prefix "${dist_name}/" "$dist_tag"
  fi

  dist_echo "  ARCHIVE ${dist_name}.tar.gz" >& 2
  dist_echo "  ARCHIVE ${dist_name}.zip" >& 2
}

dist_sign() {
  sign_name="$1"
  sign_asc="${sign_name}.asc"

  dist_echo "Signing $sign_name" >& 2

  if test x"$global_force" != x'yes'; then
    if test -f "$sign_asc"; then
      dist_echo "  SKIP $sign_asc" >& 2
      return 0
    fi
  fi

  rm -f "$sign_asc"

  if gpg --version > /dev/null 2>& 1; then
    # `ps ax` was the idiom on unix v7. Here's to
    # hoping the OS has respect for the old school.
    case "`ps ax`" in
      *gpg-agent*)
        gpg --quiet --detach-sign --armor "$sign_name" < /dev/null > /dev/null
        dist_echo "  SIGN ${sign_asc}" >& 2
      ;;
      *)
        dist_warn 'gpg-agent not running. Skipped signing.'
      ;;
    esac
  else
    dist_warn 'gpg not installed. Skipped signing.'
  fi
}

dist_sha256() {
  dist_echo "Hashing $@" >& 2

  if test x"$global_force" != x'yes'; then
    if test -f sha256sums.txt; then
      if fgrep "$1" sha256sums.txt > /dev/null 2>& 1; then
        dist_echo "  SKIP sha256sums.txt" >& 2
        return 0
      fi
    fi
  fi

  rm -f sha256sums.txt

  if echo x | sha256sum > /dev/null 2>& 1; then
    sha256sum "$@" > sha256sums.txt
    dist_echo "  HASH sha256sums.txt" >& 2
  elif echo x | sha256 > /dev/null 2>& 1; then
    sha256 "$@" > sha256sums.txt
    dist_echo "  HASH sha256sums.txt" >& 2
  else
    dist_warn 'sha256sum not installed. Skipped hashing.'
  fi
}

github_post() {
  # https://docs.github.com/en/rest
  # https://docs.github.com/en/rest/overview/resources-in-the-rest-api#oauth2-token-sent-in-a-header
  post_url="$1"
  post_type="$2"

  shift
  shift

  dist_echo "  POST $post_url" >& 2

  # Our mock Github API server.
  if test x"$global_dry" = x'yes'; then
    post_domain=`echo "$post_url" | cut -d'/' -f3`
    post_owner=`echo "$post_url" | cut -d'/' -f5`
    post_repo=`echo "$post_url" | cut -d'/' -f6`

    # https://docs.github.com/en/rest/reference/repos#create-a-release
    if test x"$post_domain" = x'api.github.com'; then
      tr -d ' \r\n' <<EOF
        {
          "upload_url": "https://uploads.github.com/repos/${post_owner}
                                                         /${post_repo}
                                                         /releases/1
                                                         /assets{?name,label}"
        }
EOF
      echo ''
      return 0
    fi

    # https://docs.github.com/en/rest/reference/repos#upload-a-release-asset
    if test x"$post_domain" = x'uploads.github.com'; then
      post_tag="$main_tag"
      post_name=`echo "$post_url" | cut -d'=' -f2 | cut -d'&' -f1`
      tr -d ' \r\n' <<EOF
        {
          "browser_download_url": "https://github.com/${post_owner}
                                                     /${post_repo}
                                                     /releases/download
                                                     /${post_tag}
                                                     /${post_name}"
        }
EOF
      echo ''
      return 0
    fi

    echo '{ "message": "Not Found" }'

    return 0
  fi

  curl -s -X POST                              \
       -H 'Accept: application/json'           \
       -H "Authorization: token $global_token" \
       -H "Content-Type: $post_type"           \
       "$@" "$post_url"
}

read_changelog() {
  log_tag="$1"
  log_ver=`echo "$log_tag" | sed -e 's/^v//'`
  log_name="${global_name}-${log_ver}"
  log_state=0

  if ! test -f "${log_name}.tar.gz"; then
    dist_error "Cannot find ${log_name}.tar.gz."
    return 1
  fi

  gunzip -c "${log_name}.tar.gz" > "${log_name}.tar"

  if ! tar tf "${log_name}.tar" | grep '^[^/]*/CHANGELOG\.md$' > /dev/null; then
    rm -f "${log_name}.tar"
    dist_error 'Cannot find CHANGELOG.md.'
    return 1
  fi

  log_file=`tar tf "${log_name}.tar" | grep '^[^/]*/CHANGELOG\.md$'`

  tar xf "${log_name}.tar" "${log_file}"

  rm -f "${log_name}.tar"

  while IFS= read -r log_line; do
    case $log_state in
      0)
        if test x"$log_line" = x"## ${log_tag}"; then
          log_state=1
        fi
      ;;
      1)
        if test x"$log_line" != x; then
          break
        fi
        log_state=2
      ;;
      2)
        if echo x"$log_line" | grep '^x## v[0-9]' > /dev/null; then
          break
        fi
        echo x"$log_line" | sed -e 's/^x//'
      ;;
    esac
  done < "$log_file"

  rm -f "$log_file"
  rmdir `dirname "$log_file"`

  if test $log_state != 2; then
    dist_error 'Could not parse changelog.'
    return 1
  fi
}

get_description() {
  desc_tag="$1"

  read_changelog "$desc_tag" | while IFS= read -r desc_line; do
    echo x"$desc_line\n" | sed -e 's/^x//' -e 's/"/\\"/g' | tr -d '\n'
  done
}

github_body() {
  body_tag="$1"
  body_commit="$2"
  body_desc="$3"

  cat <<EOF
{
  "tag_name": "${body_tag}",
  "target_commitish": "${body_commit}",
  "name": "${global_name} ${body_tag}",
  "body": "${body_desc}",
  "draft": false,
  "prerelease": false
}
EOF
}

github_release() {
  # https://docs.github.com/en/rest/reference/repos#create-a-release
  rel_tag="$1"
  rel_commit=`git rev-parse "$rel_tag"`
  rel_desc=`get_description "$rel_tag"`
  rel_url="https://api.github.com/repos/${global_owner}/${global_repo}/releases"
  rel_body=`github_body "$rel_tag" "$rel_commit" "$rel_desc"`

  dist_echo "Releasing $rel_tag (${rel_commit}):" >& 2

  github_post "$rel_url" application/json --data-raw "$rel_body" \
    | jq -r .upload_url                                          \
    | cut -d'{' -f1
}

github_upload() {
  # https://docs.github.com/en/rest/reference/repos#upload-a-release-asset
  up_name="$1"
  up_url="$2"
  up_type="$3"
  up_res=

  dist_echo "Uploading ${up_name}:" >& 2

  up_res=`github_post "${up_url}?name=${up_name}&label=${up_name}" \
                      "$up_type" --data-binary "@${up_name}"`

  up_url=`echo "$up_res" | jq -r .browser_download_url`

  dist_echo "  URL $up_url" >& 2
}

archive_repo() {
  ar_tag="$1"
  ar_ver=`echo "$ar_tag" | sed -e 's/^v//'`
  ar_name="${global_name}-${ar_ver}"

  dist_archive "$ar_tag"
  dist_sign "${ar_name}.tar.gz"
  dist_sign "${ar_name}.zip"
  dist_sha256 "${ar_name}.tar.gz" "${ar_name}.zip"
}

publish_repo() {
  pub_tag="$1"
  pub_ver=`echo "$pub_tag" | sed -e 's/^v//'`
  pub_name="${global_name}-${pub_ver}"

  dist_archive "$pub_tag"
  dist_sign "${pub_name}.tar.gz"
  dist_sign "${pub_name}.zip"
  dist_sha256 "${pub_name}.tar.gz" "${pub_name}.zip"

  if ! test -f "${pub_name}.tar.gz.asc" -a -f "${pub_name}.zip.asc"; then
    dist_error 'Files must be signed for publishing.'
    return 1
  fi

  pub_url=`github_release "$pub_tag"`

  if test x"$pub_url" = x'null'; then
    dist_error 'Github did not return an upload URL.'
    return 1
  fi

  github_upload "${pub_name}.tar.gz" "$pub_url" application/gzip
  github_upload "${pub_name}.tar.gz.asc" "$pub_url" text/plain
  github_upload "${pub_name}.zip" "$pub_url" application/zip
  github_upload "${pub_name}.zip.asc" "$pub_url" text/plain
  github_upload sha256sums.txt "$pub_url" text/plain
}

clean_repo() {
  rm -f ${global_name}-*.tar.gz
  rm -f ${global_name}-*.tar.gz.asc
  rm -f ${global_name}-*.zip
  rm -f ${global_name}-*.zip.asc
  rm -f sha256sums.txt
}

main() {
  main_action=archive
  main_deps='git'
  main_tag=''
  main_token="$GITHUB_TOKEN"
  main_auto=no
  main_dry=no
  main_force=no
  main_quiet=no

  while test $# != 0; do
    case "$1" in
      --archive|-a)
        main_action=archive
      ;;
      --publish|-p)
        main_action=publish
        main_deps="$main_deps gpg gpg-agent curl jq"

        if echo x | sha256sum > /dev/null 2>& 1; then
          : # GNU / BusyBox
        elif echo x | sha256 > /dev/null 2>& 1; then
          : # BSD
        else
          dist_error "sha256sum must be installed in order to publish."
        fi

        if ! echo x | gzip -c 2> /dev/null | gunzip -c > /dev/null 2>& 1; then
          dist_error "gzip must be installed in order to publish."
        fi
      ;;
      --log|-l)
        main_action=log
      ;;
      --clean|-c)
        main_action=clean
      ;;
      --token|-t)
        main_token="$2"
        shift
      ;;
      --autotools|-k)
        main_auto=yes
        main_deps="$main_deps autoconf"
      ;;
      --dry|-d)
        main_dry=yes
      ;;
      --force|-f)
        main_force=yes
      ;;
      --quiet|-q)
        main_quiet=yes
      ;;
      --version|-v)
        main_action=version
      ;;
      --help|-h)
        main_action=help
      ;;
      -*)
        dist_error "Invalid option '$1'."
        return 1
      ;;
      *)
        if test x"$main_tag" = x; then
          main_tag="$1"
        else
          dist_error "Invalid option '$1'."
          return 1
        fi
      ;;
    esac
    shift
  done

  if test $main_action = help; then
    dist_echo '  Usage: dist.sh [options] <tag>'
    dist_echo ''
    dist_echo '  Commands:'
    dist_echo ''
    dist_echo '    <tag>                 a git tag, commit, or branch name'
    dist_echo ''
    dist_echo '  Options:'
    dist_echo ''
    dist_echo '    -a, --archive         create archive (default)'
    dist_echo '    -p, --publish         create and publish archive'
    dist_echo '    -l, --log             read changelog from archive'
    dist_echo '    -c, --clean           clean working directory'
    dist_echo '    -t, --token <string>  set github api key'
    dist_echo '    -k, --autotools       use autotools dist'
    dist_echo '    -d, --dry             dry run for upload'
    dist_echo '    -f, --force           force dist rebuild'
    dist_echo '    -q, --quiet           silence output'
    dist_echo '    -v, --version         output version number'
    dist_echo '    -h, --help            output usage information'
    dist_echo ''
    dist_echo '  Environment Variables:'
    dist_echo ''
    dist_echo '    GITHUB_TOKEN          github api key'
    dist_echo ''
    return 0
  fi

  if test $main_action = version; then
    dist_echo 'dist.sh 0.0.0'
    return 0
  fi

  global_quiet="$main_quiet"

  if ! test -d .git; then
    dist_error 'Not a git repository.'
    return 1
  fi

  if test -f .distrc; then
    . ./.distrc || eval "`cat .distrc`"
  elif test -f distrc; then
    . ./distrc || eval "`cat distrc`"
  else
    dist_error 'Could not find distrc.'
    return 1
  fi

  if test x"$project_name" = x; then
    project_name="$repo_name"
  fi

  if test x"$repo_name" = x; then
    repo_name="$project_name"
  fi

  if test x"$project_name" = x -o x"$repo_name" = x -o x"$owner_name" = x; then
    dist_error 'Invalid distrc (must set {project,repo,owner}_name).'
    return 1
  fi

  global_name="$project_name"
  global_repo="$repo_name"
  global_owner="$owner_name"

  if test $main_action = clean; then
    clean_repo
    return 0
  fi

  for main_dep in $main_deps; do
    if ! $main_dep --version > /dev/null 2>& 1; then
      dist_error "$main_dep must be installed in order to $main_action."
      return 1
    fi
  done

  if test x"$main_tag" = x; then
    main_tag=`git tag -l 'v*' --sort v:refname | tail -n 1`
  fi

  if test x"$main_tag" = x; then
    dist_error 'Must provide a tag.'
    return 1
  fi

  if ! git rev-parse "$main_tag" > /dev/null 2>& 1; then
    dist_error "Invalid tag '$main_tag'."
    return 1
  fi

  global_auto="$main_auto"
  global_dry="$main_dry"
  global_force="$main_force"

  case "$main_action" in
    archive)
      archive_repo "$main_tag"
    ;;
    publish)
      if test x"$main_token" = x; then
        dist_echo 'GitHub Token: ' | tr -d '\n'
        if test x"$BASH_VERSION$ZSH_VERSION" != x; then
          IFS= read -rs main_token
          dist_echo ''
        else
          IFS= read -r main_token
        fi
      fi

      if test x"$main_token" = x; then
        dist_error 'Must provide a token.'
        return 1
      fi

      global_token="$main_token"

      publish_repo "$main_tag"
    ;;
    log)
      read_changelog "$main_tag"
    ;;
    *)
      dist_error "Invalid action '$main_action'."
      return 1
    ;;
  esac
}

main "$@"
