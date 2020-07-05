/*!
 * thread.c - thread support for libtorsion tests
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#include <stdlib.h>
#include "thread.h"

#if defined(_WIN32)
#  include <windows.h> /* CreateThread */
#  pragma comment(lib, "kernel32.lib")
#else
#  include <pthread.h>
#endif

struct torsion_thread_s {
#ifdef _WIN32
  HANDLE handle;
#else
  pthread_t handle;
#endif
};

struct torsion_thread_s *
torsion_thread_alloc(void) {
  struct torsion_thread_s *thread = malloc(sizeof(struct torsion_thread_s));

  if (thread == NULL)
    abort();

  return thread;
}

void
torsion_thread_free(struct torsion_thread_s *thread) {
  free(thread);
}

#ifdef _WIN32
typedef struct torsion_thread_args_s {
  torsion_thread_start_f *start_routine;
  void *arg;
} torsion_thread_args_t;

static DWORD WINAPI
torsion_thread_run(void *ptr) {
  torsion_thread_args_t args;

  args = *((torsion_thread_args_t *)ptr);

  free(ptr);

  (void *)args.start_routine(args.arg);

  return ERROR_SUCCESS;
}
#endif

int
torsion_thread_create(struct torsion_thread_s *thread,
                      const torsion_thread_attr_t *attr,
                      torsion_thread_start_f *start_routine,
                      void *arg) {
#ifdef _WIN32
  torsion_thread_args_t *args;

  if (attr != NULL)
    return -1;

  args = malloc(sizeof(torsion_thread_args_t));

  if (args == NULL)
    return -1;

  args->start_routine = start_routine;
  args->arg = arg;

  thread->handle = CreateThread(NULL, 0, torsion_thread_run, args, 0, NULL);

  if (thread->handle == NULL) {
    free(args);
    return -1;
  }

  return 0;
#else
  if (attr != NULL)
    return -1;

  return pthread_create(&thread->handle, NULL, start_routine, arg);
#endif
}

int
torsion_thread_join(struct torsion_thread_s *thread, void **retval) {
#ifdef _WIN32
  if (retval != NULL)
    return -1;

  WaitForSingleObject(thread->handle, INFINITE);

  if (CloseHandle(thread->handle) == FALSE)
    return -1;

  return 0;
#else
  if (retval != NULL)
    return -1;

  return pthread_join(thread->handle, NULL);
#endif
}
