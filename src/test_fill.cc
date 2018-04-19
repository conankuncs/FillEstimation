#include <vector>
#include <unordered_map>
using namespace std;

//2D hashkey function
static inline int hash_key_chunks(int i,int j,int B, int n) {
  return (1 + ((int)(n/B) + (n % (int)B == 0 ? 0 : 1)) * (int)((i-1)/B) + (int)((j-1)/B));
}
static inline int local_y_chunk (int key, int B, int n) {
  int chunks_per_row = ((int)(n/B) + (n % (int)B == 0 ? 0 : 1));
  return (int)(key/chunks_per_row) + (key % chunks_per_row == 0 ? 0 : 1);
}

static inline int local_x_chunk (int key, int B, int local_y, int n) {
  int chunks_per_row = ((int)(n/B) + (n % (int)B == 0 ? 0 : 1));
  return key - chunks_per_row*(local_y-1);
}

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

int estimate_fill_coo_2d_variant (size_t m, size_t n, size_t nnz, vector<coo_2d_simplified_2_by_2> coo, size_t B, double epsilon, double delta, double *fill, int verbose);

int estimate_fill_coo_3d (int x,
                   int y,
                   int z,
                   int nnz,
                   vector<coo_3d_simplified> &coo,
                   int B,
                   double epsilon,
                   double delta,
                   double *fill,
                   int verbose);

int test_csr (size_t m,
          size_t n,
          size_t nnz,
          const size_t *ptr,
          const size_t *ind,
          int B,
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
  for (int i = 0; i < B * B; i++) {
    fill[i] = 0;
  }

  //Benchmark some runs
  double time = -wall_time();
  for (int t = 0; t < trials; t++){
    estimate_fill_coo_2d(m, n, nnz, coo, B, epsilon, delta, fill + t * B * B, verbose);
  }
  time += wall_time();

  printf("{\n");
  int i = 0;
  int exp2 = 1;
  if (results) {
    printf("  \"Original results\": [\n");
    for (int t = 0; t < trials; t++) {
      printf("    [\n");
      for (int b_r = 1; b_r <= B; b_r++) {
        if (b_r == exp2) printf("      [\n");
        for (int b_c = 1; b_c <= B; b_c++) {
          if(b_c == b_r && exp2 == b_r) printf("B=(%d,%d) %.*e%s", b_r, b_c, DECIMAL_DIG, fill[i], b_c <= B - 1 ? ", " : "");
          i++;
        }
        if (b_r == exp2) {
          printf("      ]%s\n", b_r <= B - 1 ? "," : "");
          exp2 *= 2;
        }
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

int test_coo_2d_variant (size_t m,
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

  // Hash map for list of 2x2 chunks
  unordered_map<int, coo_2d_simplified_2_by_2> mp;
  // Group every element as 2x2 chunks
  for (int i = 0; i < nnz; i++) {
    coo_2d_simplified coo2d;
    coo2d.x = coo[i].x;
    coo2d.y = coo[i].y;
    // Find the key
    int key = hash_key_chunks(coo[i].y,coo[i].x,2,n);
    unordered_map<int, coo_2d_simplified_2_by_2>::iterator it = mp.find(key);
    if(it == mp.end()) {
      // Not found. Create new chunk
      coo_2d_simplified_2_by_2 new_chunk;
      new_chunk.index=0;
      new_chunk.arr[new_chunk.index++]=coo2d;
      new_chunk.local_y = local_y_chunk(key, 2, n);
      new_chunk.local_x = local_x_chunk(key, 2, new_chunk.local_y, n);
      mp.insert(make_pair(key, new_chunk));
    } else {
      coo_2d_simplified_2_by_2 &my_chunk = it->second;
      my_chunk.arr[my_chunk.index++]=coo2d;
    }
  }

  // Make a vector of all chunks
  vector<coo_2d_simplified_2_by_2> chunk_list;
  for (auto it = mp.begin(); it != mp.end(); ++it) {
    chunk_list.push_back(it->second);
  }

  //Load problem into cache
  estimate_fill_coo_2d_variant(m/2, n/2, chunk_list.size(), chunk_list, B/2, epsilon, delta, fill, verbose);
  for (size_t i = 0; i < B * B; i++) {
    fill[i] = 0;
  }

  //Benchmark some runs
  double time = -wall_time();
  for (int t = 0; t < trials; t++){
    estimate_fill_coo_2d_variant(m/2, n/2, chunk_list.size(), chunk_list, B/2, epsilon, delta, fill + t * B/2 * B/2, verbose);
  }
  time += wall_time();

  printf("{\n");
  size_t i = 0;
  int exp2 = 1;
  if (results) {
    printf("  \"results\": [\n");
    for (int t = 0; t < trials; t++) {
      printf("    [\n");
      for (int b_r = 1; b_r <= B/2; b_r++) {
        if (b_r == exp2) printf("      [\n");
        for (int b_c = 1; b_c <= B/2; b_c++) {
          if(b_c == b_r && exp2 == b_r) printf("B=(%d,%d) %.*e%s", b_r*2, b_c*2, DECIMAL_DIG, fill[i], b_c <= B - 1 ? ", " : "");
          i++;
        }
        if (b_r == exp2) {
          printf("      ]%s\n", b_r <= B - 1 ? "," : "");
          exp2 *= 2;
        }
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


int test_coo_3d (int x,
          int y,
          int z,
          int nnz,
          vector<coo_3d_simplified> &coo,
          int B,
          double epsilon,
          double delta,
          int trials,
          int clock,
          int results,
          int verbose) {

  double *fill = (double*)malloc(sizeof(double) * B * B * B * trials);
  for (size_t i = 0; i < B * B * B * trials; i++) {
    fill[i] = 0;
  }

  //Load problem into cache
  /*estimate_fill_coo_2d(m, n, nnz, coo, B, epsilon, delta, fill, verbose);
  for (size_t i = 0; i < B * B; i++) {
    fill[i] = 0;
  }*/

  //Benchmark some runs
  double time = -wall_time();
  for (int t = 0; t < trials; t++){
    estimate_fill_coo_3d(x, y, z, nnz, coo, B, epsilon, delta, fill + t * B * B, verbose);
  }
  time += wall_time();

  printf("{\n");
  size_t i = 0;
  if (results) {
    printf("  \"results\": [\n");
    for (int t = 0; t < trials; t++) {
      printf("    [\n");
      for (int b_d = 1; b_d <= B; b_d++) {
        printf("      [\n");
        for (int b_r = 1; b_r <= B; b_r++) {
          printf("        [        ");
          for (int b_c = 1; b_c <= B; b_c++) {
            printf("%.*e%s", DECIMAL_DIG, fill[i], b_c <= B - 1 ? ", " : "");
            i++;
          }
          printf("        ]%s\n", b_r <= B - 1 ? "," : "");
        }
        printf("      ]%s\n", b_d <= B - 1 ? "," : "");
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