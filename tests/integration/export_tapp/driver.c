/* Runs the C code produced by export_tapp_generate and compares its results
 * against naive loops. */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tapp_complex.c"
#include "tapp_real.c"

enum { NOCC = 3, NVIRT = 5, MAX_ENTRIES = 16 };

/* In-memory storage backing the generated load/store callbacks */
typedef struct {
  char name[32];
  size_t bytes;
  unsigned char *data;
} entry;

typedef struct {
  entry entries[MAX_ENTRIES];
  size_t size;
  size_t elem_size;
} storage;

static entry *find(storage *st, const char *name) {
  for (size_t i = 0; i < st->size; ++i) {
    if (strcmp(st->entries[i].name, name) == 0) {
      return &st->entries[i];
    }
  }
  return NULL;
}

static int storage_load(void *user, const char *name, void *buffer,
                        size_t count) {
  storage *st = user;
  const entry *e = find(st, name);
  if (!e || e->bytes != count * st->elem_size) {
    return 1;
  }
  memcpy(buffer, e->data, e->bytes);
  return 0;
}

static int storage_store(void *user, const char *name, const void *buffer,
                         size_t count) {
  storage *st = user;
  entry *e = find(st, name);
  if (!e) {
    if (st->size == MAX_ENTRIES || strlen(name) >= sizeof e->name) {
      return 1;
    }
    e = &st->entries[st->size++];
    strcpy(e->name, name);
    e->bytes = 0;
    e->data = NULL;
  }
  unsigned char *data = realloc(e->data, count * st->elem_size);
  if (!data) {
    return 1;
  }
  e->data = data;
  e->bytes = count * st->elem_size;
  memcpy(e->data, buffer, e->bytes);
  return 0;
}

static void storage_free(storage *st) {
  for (size_t i = 0; i < st->size; ++i) {
    free(st->entries[i].data);
  }
  st->size = 0;
}

static int real_load(void *user, const char *name, double *buffer,
                     size_t count) {
  return storage_load(user, name, buffer, count);
}

static int real_store(void *user, const char *name, const double *buffer,
                      size_t count) {
  return storage_store(user, name, buffer, count);
}

static int complex_load(void *user, const char *name, double complex *buffer,
                        size_t count) {
  return storage_load(user, name, buffer, count);
}

static int complex_store(void *user, const char *name,
                         const double complex *buffer, size_t count) {
  return storage_store(user, name, buffer, count);
}

/* Deterministic, non-trivial test data */
static double value(size_t k, double shift) {
  return (double)((k * 37 + 11) % 17) / 8.0 - 1.0 + shift;
}

static double *fill(storage *st, const char *name, size_t count,
                    double shift) {
  double *data = malloc(count * sizeof *data);
  if (!data) {
    return NULL;
  }
  for (size_t k = 0; k < count; ++k) {
    data[k] = value(k, shift);
  }
  if (storage_store(st, name, data, count) != 0) {
    free(data);
    return NULL;
  }
  return data;
}

static double complex *fill_complex(storage *st, const char *name,
                                    size_t count, double shift) {
  double complex *data = malloc(count * sizeof *data);
  if (!data) {
    return NULL;
  }
  for (size_t k = 0; k < count; ++k) {
    data[k] = CMPLX(value(k, shift), value(k + 5, -shift));
  }
  if (storage_store(st, name, data, count) != 0) {
    free(data);
    return NULL;
  }
  return data;
}

static int close_to(double actual, double expected) {
  return fabs(actual - expected) <= 1e-12 * (1.0 + fabs(expected));
}

static int check(int ok, const char *what) {
  if (!ok) {
    fprintf(stderr, "FAILED: %s\n", what);
  }
  return ok ? 0 : 1;
}

/* Row-major offsets */
#define OO(i, j) ((size_t)(i) * NOCC + (j))
#define OV(i, a) ((size_t)(i) * NVIRT + (a))
#define VO(a, i) ((size_t)(a) * NOCC + (i))
#define VV(a, b) ((size_t)(a) * NVIRT + (b))
#define VVOO(a, b, i, j) ((VV(a, b) * NOCC + (i)) * NOCC + (j))
#define OOVV(i, j, a, b) ((OO(i, j) * NVIRT + (a)) * NVIRT + (b))
#define VVVV(a, b, c, d) ((VV(a, b) * NVIRT + (c)) * NVIRT + (d))

static int test_real(void) {
  const size_t n_vvoo = (size_t)NVIRT * NVIRT * NOCC * NOCC;
  int failures = 0;

  storage st = {.elem_size = sizeof(double)};
  const real_io io = {real_load, real_store, &st};
  const double s = 0.75;

  double *g = fill(&st, "g_vvvv", (size_t)NVIRT * NVIRT * NVIRT * NVIRT, 0.1);
  double *t = fill(&st, "t_vvoo", n_vvoo, 0.2);
  double *f_vv = fill(&st, "f_vv", (size_t)NVIRT * NVIRT, 0.3);
  double *v = fill(&st, "v_vvoo", n_vvoo, 0.4);
  double *w = fill(&st, "w_oovv", n_vvoo, 0.5);
  double *f_ov = fill(&st, "f_ov", (size_t)NOCC * NVIRT, 0.6);
  double *u = fill(&st, "u_ov", (size_t)NOCC * NVIRT, 0.7);
  double *p = fill(&st, "p_vo", (size_t)NVIRT * NOCC, 0.8);
  double *q = fill(&st, "q_vo", (size_t)NVIRT * NOCC, 0.9);
  if (!g || !t || !f_vv || !v || !w || !f_ov || !u || !p || !q ||
      storage_store(&st, "s", &s, 1) != 0) {
    fprintf(stderr, "Failed to set up the real test data\n");
    return 1;
  }

  const real_dims dims = {.nocc = NOCC, .nvirt = NVIRT};
  real_plans plans;
  TAPP_error error = 0;
  failures += check(real_create_plans(&plans, &dims, &error) == REAL_SUCCESS,
                    "real_create_plans");

  /* Every call must reproduce the reference, which requires the plans to be
     reusable and the results to be rebuilt from scratch */
  for (int run = 0; run < 2 && failures == 0; ++run) {
    failures += check(real_residual(&plans, &io, &error) == REAL_SUCCESS,
                      "real_residual");
    failures += check(real_singles(&plans, &io, &error) == REAL_SUCCESS,
                      "real_singles");
    failures += check(real_energy(&plans, &io, &error) == REAL_SUCCESS,
                      "real_energy");
    if (failures != 0) {
      break;
    }

    const entry *R = find(&st, "R_vvoo");
    const entry *S = find(&st, "S_vvoo");
    const entry *X = find(&st, "X_vo");
    const entry *Y = find(&st, "Y_vo");
    const entry *E = find(&st, "E");
    failures += check(R && R->bytes == n_vvoo * sizeof(double), "R stored");
    failures += check(S && S->bytes == n_vvoo * sizeof(double), "S stored");
    failures +=
        check(X && X->bytes == (size_t)NVIRT * NOCC * sizeof(double),
              "X stored");
    failures +=
        check(Y && Y->bytes == (size_t)NVIRT * NOCC * sizeof(double),
              "Y stored");
    failures += check(E && E->bytes == sizeof(double), "E stored");
    if (failures != 0) {
      break;
    }

    const double *R_data = (const double *)R->data;
    const double *S_data = (const double *)S->data;
    const double *X_data = (const double *)X->data;
    const double *Y_data = (const double *)Y->data;
    for (int a = 0; a < NVIRT; ++a) {
      for (int b = 0; b < NVIRT; ++b) {
        for (int i = 0; i < NOCC; ++i) {
          for (int j = 0; j < NOCC; ++j) {
            double ref = s * v[VVOO(a, b, i, j)];
            for (int c = 0; c < NVIRT; ++c) {
              ref += f_vv[VV(a, c)] * t[VVOO(c, b, i, j)];
              for (int d = 0; d < NVIRT; ++d) {
                ref += 0.5 * g[VVVV(a, b, c, d)] * t[VVOO(c, d, i, j)];
              }
            }
            failures += check(close_to(R_data[VVOO(a, b, i, j)], ref),
                              "R value");
            failures += check(close_to(S_data[VVOO(a, b, i, j)],
                                       R_data[VVOO(a, b, i, j)] -
                                           R_data[VVOO(b, a, i, j)]),
                              "S value");
          }
        }
      }
    }

    for (int a = 0; a < NVIRT; ++a) {
      for (int i = 0; i < NOCC; ++i) {
        double X_ref = 0.0;
        double Y_ref = 0.0;
        for (int c = 0; c < NVIRT; ++c) {
          X_ref += f_vv[VV(a, c)] * p[VO(c, i)];
          Y_ref += 2.0 * f_vv[VV(a, c)] * q[VO(c, i)];
        }
        failures += check(close_to(X_data[VO(a, i)], X_ref), "X value");
        failures += check(close_to(Y_data[VO(a, i)], Y_ref), "Y value");
      }
    }

    double E_ref = 0.0;
    for (int i = 0; i < NOCC; ++i) {
      for (int j = 0; j < NOCC; ++j) {
        for (int a = 0; a < NVIRT; ++a) {
          for (int b = 0; b < NVIRT; ++b) {
            E_ref += 0.25 * w[OOVV(i, j, a, b)] * t[VVOO(a, b, i, j)] +
                     f_ov[OV(i, a)] * t[VVOO(a, b, i, j)] * u[OV(j, b)];
          }
        }
      }
    }
    failures += check(close_to(*(const double *)E->data, E_ref), "E value");
  }

  /* A failing load must be reported and must not leak (checked by sanitizers) */
  if (failures == 0) {
    strcpy(find(&st, "u_ov")->name, "missing");
    failures += check(real_energy(&plans, &io, &error) == REAL_IO_FAILURE,
                      "real_energy reports a failed load");
  }

  real_destroy_plans(&plans);
  storage_free(&st);
  free(g);
  free(t);
  free(f_vv);
  free(v);
  free(w);
  free(f_ov);
  free(u);
  free(p);
  free(q);

  return failures;
}

static int test_complex(void) {
  const size_t n_vvoo = (size_t)NVIRT * NVIRT * NOCC * NOCC;
  int failures = 0;

  storage st = {.elem_size = sizeof(double complex)};
  const complex_io io = {complex_load, complex_store, &st};
  const double complex s = CMPLX(0.5, -1.25);

  double complex *t = fill_complex(&st, "t_vvoo", n_vvoo, 0.2);
  double complex *w = fill_complex(&st, "w_oovv", n_vvoo, 0.5);
  if (!t || !w || storage_store(&st, "s", &s, 1) != 0) {
    fprintf(stderr, "Failed to set up the complex test data\n");
    return 1;
  }

  const complex_dims dims = {.nocc = NOCC, .nvirt = NVIRT};
  complex_plans plans;
  TAPP_error error = 0;
  failures +=
      check(complex_create_plans(&plans, &dims, &error) == COMPLEX_SUCCESS,
            "complex_create_plans");
  if (failures == 0) {
    failures += check(complex_energy(&plans, &io, &error) == COMPLEX_SUCCESS,
                      "complex_energy");
  }

  const entry *E = find(&st, "E");
  failures += check(E && E->bytes == sizeof(double complex), "E stored");
  if (failures == 0) {
    double complex E_ref = 0.0;
    for (int i = 0; i < NOCC; ++i) {
      for (int j = 0; j < NOCC; ++j) {
        for (int a = 0; a < NVIRT; ++a) {
          for (int b = 0; b < NVIRT; ++b) {
            E_ref += conj(s) * w[OOVV(i, j, a, b)] * t[VVOO(a, b, i, j)];
          }
        }
      }
    }
    const double complex E_value = *(const double complex *)E->data;
    failures += check(close_to(creal(E_value), creal(E_ref)) &&
                          close_to(cimag(E_value), cimag(E_ref)),
                      "complex E value");
  }

  complex_destroy_plans(&plans);
  storage_free(&st);
  free(t);
  free(w);

  return failures;
}

int main(void) {
  const int failures = test_real() + test_complex();
  if (failures != 0) {
    fprintf(stderr, "%d check(s) failed\n", failures);
    return EXIT_FAILURE;
  }

  puts("All checks passed");
  return EXIT_SUCCESS;
}
