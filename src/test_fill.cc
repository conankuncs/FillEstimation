#include <vector>
using namespace std;

extern "C" {
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "util.h"

#define TIMEOUT 0.1

//Return time time of day as a double-precision floating point value.
double wall_time (void) {
  struct timeval t;
  gettimeofday(&t, NULL);
  return 1.0*t.tv_sec + 1.0e-6*t.tv_usec;
}

int estimate_fill_csr (size_t m,
                   size_t n,
                   size_t nnz,
                   const size_t *ptr,
                   const size_t *ind,
                   size_t B,
                   double epsilon,
                   double delta,
                   double *fill,
                   int verbose);

int estimate_fill_coo_2d (size_t m, size_t n, size_t nnz, vector<coo_2d> coo, size_t B, double epsilon, double delta, double *fill, int verbose);

int test_csr (size_t m,
          size_t n,
          size_t nnz,
          const size_t *ptr,
          const size_t *ind,
          size_t B,
          double epsilon,
          double delta,
          int trials,
          int clock,
          int results,
          int verbose) {

  double *fill = (double*)malloc(sizeof(double) * B * B * trials);
  for (size_t i = 0; i < B * B * trials; i++) {
    fill[i] = 0;
  }

  //Load problem into cache
  estimate_fill_csr(m, n, nnz, ptr, ind, B, epsilon, delta, fill, verbose);
  for (size_t i = 0; i < B * B; i++) {
    fill[i] = 0;
  }


  //Benchmark some runs
  double time = -wall_time();
  for (int t = 0; t < trials; t++){
    estimate_fill_csr(m, n, nnz, ptr, ind, B, epsilon, delta, fill + t * B * B, verbose);
  }
  time += wall_time();

  printf("{\n");
  size_t i = 0;
  if (results) {
    printf("  \"results\": [\n");
    for (int t = 0; t < trials; t++) {
      printf("    [\n");
      for (size_t b_r = 1; b_r <= B; b_r++) {
        printf("      [\n");
        for (size_t b_c = 1; b_c <= B; b_c++) {
          printf("%.*e%s", DECIMAL_DIG, fill[i], b_c <= B - 1 ? ", " : "");
          i++;
        }
        printf("      ]%s\n", b_r <= B - 1 ? "," : "");
      }
      printf("    ]%s\n", t < trials - 1 ? "," : "");
    }
    printf("  ]%s\n", clock ? "," : "");
  }
  if (clock) {
    printf("  \"time_total\": %.*e,\n", DECIMAL_DIG, time);
    printf("  \"time_mean\": %.*e%s,\n", DECIMAL_DIG, time/trials, 0 ? "," : "");
  }
  //printf("\n}\n");

  free(fill);
  return 0;
}

int test_coo_2d (size_t m,
          size_t n,
          size_t nnz,
          vector<coo_2d> &coo,
          size_t B,
          double epsilon,
          double delta,
          int trials,
          int clock,
          int results,
          int verbose) {

  double *fill = (double*)malloc(sizeof(double) * B * B * trials);
  for (size_t i = 0; i < B * B * trials; i++) {
    fill[i] = 0;
  }

  //Load problem into cache
  estimate_fill_coo_2d(m, n, nnz, coo, B, epsilon, delta, fill, verbose);
  for (size_t i = 0; i < B * B; i++) {
    fill[i] = 0;
  }


  //Benchmark some runs
  double time = -wall_time();
  for (int t = 0; t < trials; t++){
    estimate_fill_coo_2d(m, n, nnz, coo, B, epsilon, delta, fill + t * B * B, verbose);
  }
  time += wall_time();

  printf("{\n");
  size_t i = 0;
  if (results) {
    printf("  \"results\": [\n");
    for (int t = 0; t < trials; t++) {
      printf("    [\n");
      for (size_t b_r = 1; b_r <= B; b_r++) {
        printf("      [\n");
        for (size_t b_c = 1; b_c <= B; b_c++) {
          printf("%.*e%s", DECIMAL_DIG, fill[i], b_c <= B - 1 ? ", " : "");
          i++;
        }
        printf("      ]%s\n", b_r <= B - 1 ? "," : "");
      }
      printf("    ]%s\n", t < trials - 1 ? "," : "");
    }
    printf("  ]%s\n", clock ? "," : "");
  }
  if (clock) {
    printf("  \"time_total\": %.*e,\n", DECIMAL_DIG, time);
    printf("  \"time_mean\": %.*e%s,\n", DECIMAL_DIG, time/trials, 0 ? "," : "");
  }
  //printf("\n}\n");

  free(fill);
  return 0;
}

}