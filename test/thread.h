/*!
 * thread.h - thread support for libtorsion tests
 * Copyright (c) 2020, Christopher Jeffrey (MIT License).
 * https://github.com/bcoin-org/libtorsion
 */

#ifndef TORSION_THREAD_H
#define TORSION_THREAD_H

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

#endif /* TORSION_THREAD_H */
