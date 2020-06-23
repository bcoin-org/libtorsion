#ifndef _TORSION_OS_H
#define _TORSION_OS_H

#if defined(TORSION_HAVE_RNG) && defined(TORSION_TEST)
#  if defined(__CloudABI__) || defined(__EMSCRIPTEN__) || defined(__wasm__)
/* Nothing. */
#  elif defined(_WIN32)
#    define TORSION_HAVE_THREADS
#    define TORSION_HAVE_WINTHREAD
#  elif defined(__linux__) || defined(__APPLE__)
#    define TORSION_HAVE_THREADS
#    define TORSION_HAVE_PTHREAD
#  endif
#endif

uint64_t
torsion_hrtime(void);

#ifdef TORSION_HAVE_THREADS
typedef struct torsion_thread_s torsion_thread_t;
typedef void torsion_thread_attr_t;
typedef void *torsion_thread_start_f(void *);

torsion_thread_t *
torsion_thread_alloc(void);

void
torsion_thread_free(torsion_thread_t *thread);

int
torsion_thread_create(torsion_thread_t *thread,
                      const torsion_thread_attr_t *attr,
                      torsion_thread_start_f *start_routine,
                      void *arg);

int
torsion_thread_join(torsion_thread_t *thread, void **retval);
#endif /* TORSION_HAVE_THREADS */

#endif /* _TORSION_OS_H */
