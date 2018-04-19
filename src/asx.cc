#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unordered_map>
#include <sys/time.h>
#include <vector>

//2D hashkey function
static inline int hash_key(int i,int j,int B, int n) {
  return (1 + ((int)(n/B) + (n % (int)B == 0 ? 0 : 1)) * (int)((i-1)/B) + (int)((j-1)/B));
}

//3D hashkey function
static inline int hash_key(int i, int j, int k, int B, int x, int y) {
  return (1 + ((int)(x/B) + (x % B == 0 ? 0 : 1))*((int)(y/B) + (y % B == 0 ? 0 : 1))*((int)((k-1)/B)) + ((int)(x/B) + (x % B == 0 ? 0 : 1))*((int)((j-1)/B)) + ((int)((i-1)/B)));
}

using namespace std;

extern "C" {

#include "util.h"

char *name () {
  return (char *)"asx";
}

//Return time time of day as a double-precision floating point value.
double wall_time_2 (void) {
  struct timeval t;
  gettimeofday(&t, NULL);
  return 1.0*t.tv_sec + 1.0e-6*t.tv_usec;
}

#define REPLACEMENT

/**
 *  Given an m by n CSR matrix A, estimates the fill ratio if the matrix were
 *  converted into b_r by b_c BCSR format. The fill ratio is b_r times b_c times
 *  the number of nonzero blocks in the BCSR format divided by the number of
 *  nonzeros. All estimates should be accurate to relative error epsilon with
 *  probability at least (1 - delta).
 *
 *  The caller supplies this routine with a maximum row and column block size B,
 *  and this routine returns the estimated fill ratios for all
 *  1 <= b_r, b_c <= B.
 *
 *  This routine assumes the CSR matrix uses full storage, and assumes that
 *  column indicies are sorted.
 *
 *  \param[in] m Logical number of matrix rows
 *  \param[in] n Logical number of matrix columns
 *  \param[in] nnz Logical number of matrix nonzeros
 *  \param[in] *ptr CSR row pointers.
 *  \param[in] *ind CSR column indices.
 *  \param[in] B Maximum desired block size
 *  \param[in] epsilon Epsilon
 *  \param[in] delta Delta
 *  \param[out] *fill Fill ratios for all specified b_r, b_c in order
 *  \param[in] verbose 0 if you should be quiet
 *
 *  Note that the fill ratios should be stored according to the following order:
 *  size_t i = 0;
 *  for (size_t b_r = 1; b_r <= B; b_r++) {
 *    for (size_t b_c = 1; b_c <= B; b_c++) {
 *      //fill[i] = fill for b_r, b_c, o_r, o_c
 *      i++;
 *    }
 *  }
 *
 *  \returns On success, returns 0. On error, returns an error code.
 */
int estimate_fill_csr (size_t m,
                   size_t n,
                   size_t nnz,
                   const size_t *ptr,
                   const size_t *ind,
                   size_t B,
                   double epsilon,
                   double delta,
                   double *fill,
                   int verbose) {
  size_t W = 2 * B;
  int Z[W][W];

  double T = 2 * log(B/delta) * B * B / (epsilon * epsilon);
  size_t s;

#ifdef REPLACEMENT
  s = T;
#else
  //s = ((T - T/nnz) + sqrt((T - T / nnz) * (T - T / nnz)  + 4 * T * (1 + T / N)))/(2 + 2 * T / nnz);
  s = ((T - T/nnz) + sqrt(T * (T + (2 * T + T / nnz) / nnz + 4)))/(2 + 2 * T / nnz);
#endif
  s = min(s, nnz);

  //Sample s items
  size_t *samples = (size_t*)malloc(s*sizeof(size_t));
  size_t *samples_i = (size_t*)malloc(s*sizeof(size_t));
  size_t *samples_j = (size_t*)malloc(s*sizeof(size_t));

#ifdef REPLACEMENT
  if (s == nnz) {
    for (size_t i = 0; i < nnz; i++) {
      samples[i] = i;
    }
  } else {
    for (size_t i = 0; i < s; i++) {
      samples[i] = random_range(0, nnz);
    }
  }
#else
  random_choose(samples, s, 0, nnz);
#endif

  //Create arrays of i and j
  sort(samples, s);
  {
    size_t i = 0;
    for (size_t t = 0; t < s; t++) {
      if (ptr[i + 1] <= samples[t]) {
        i = search_strict(ptr, i, m, samples[t]) - 1;
      }
      samples_i[t] = i;
      samples_j[t] = ind[samples[t]];
    }
  }

  for (size_t t = 0; t < s; t++) {
    size_t i = samples_i[t];
    size_t j = samples_j[t];

    //compute x for some i, j
    for (int r = 0; r < W; r++) {
      for (int c = 0; c < W; c++) {
        Z[r][c] = 0;
      }
    }

    for (size_t ii = max(i, B - 1) - (B - 1); ii <= min(i + (B - 1), m - 1); ii++) {
      int r = (B + ii) - i;
      size_t jj;
      size_t jj_min = max(j, B - 1) - (B - 1);
      size_t jj_max = min(j + (B - 1), n - 1);

      size_t scan = search(ind, ptr[ii], ptr[ii + 1], jj_min);

      while (scan < ptr[ii + 1] && (jj = ind[scan]) <= jj_max) {
        int c = (B + jj) - j;
        Z[r][c] = 1;
        scan++;
      }
    }

    for (int r = 1; r < W; r++) {
      for (int c = 1; c < W; c++) {
        Z[r][c] += Z[r][c - 1];
      }
    }

    for (int r = 1; r < W; r++) {
      for (int c = 1; c < W; c++) {
        Z[r][c] += Z[r - 1][c];
      }
    }

    int fill_index = 0;
    for (int b_r = 1; b_r <= B; b_r++) {
      int r_hi = B + b_r - 1 - (i % b_r);
      int r_lo = r_hi - b_r;
      for (int b_c = 1; b_c <= B; b_c++) {
        int c_hi = B + b_c - 1 - (j % b_c);
        int c_lo = c_hi - b_c;
        int y_0 = Z[r_hi][c_hi] - Z[r_lo][c_hi] - Z[r_hi][c_lo] + Z[r_lo][c_lo];
        fill[fill_index] += 1.0/y_0;
        fill_index++;
      }
    }
  }

  int fill_index = 0;
  for (int b_r = 1; b_r <= B; b_r++) {
    for (int b_c = 1; b_c <= B; b_c++) {
      fill[fill_index] *= b_r * b_c / (double)s;
      fill_index++;
    }
  }

  free(samples);
  free(samples_i);
  free(samples_j);
  return 0;
}

int estimate_fill_coo_2d (int m,
                   int n,
                   int nnz,
                   vector<coo_2d> &coo,
                   int B,
                   double epsilon,
                   double delta,
                   double *fill,
                   int verbose) {
  int W = 2 * B;
  int Z[W][W];

  double T = 2 * log(B/delta) * B * B / (epsilon * epsilon);
  int s;

double time = -wall_time_2();

#ifdef REPLACEMENT
  s = T;
#else
  //s = ((T - T/nnz) + sqrt((T - T / nnz) * (T - T / nnz)  + 4 * T * (1 + T / N)))/(2 + 2 * T / nnz);
  s = ((T - T/nnz) + sqrt(T * (T + (2 * T + T / nnz) / nnz + 4)))/(2 + 2 * T / nnz);
#endif
  s = min(s, nnz);

  //Sample s items
  int *samples = (int*)malloc(s*sizeof(int));

// Randomized Sampling
#ifdef REPLACEMENT
  if (s == nnz) {
    for (int i = 0; i < nnz; i++) {
      samples[i] = i;
    }
  } else {
    for (int i = 0; i < s; i++) {
      samples[i] = random_range(0, nnz);
    }
  }
#else
  random_choose(samples, s, 0, nnz);
#endif
  sort_int(samples, s);

  time += wall_time_2();
  //printf("  \"time for random sampling\": %.*e,\n", DECIMAL_DIG, time);

  unordered_map<int, vector<coo_2d_simplified>> mp;
  
  // put all elements into hash map

  time = -wall_time_2();

  for (int t = 0; t < nnz; t++) {
    int i = coo[t].y;
    int j = coo[t].x;
    coo_2d_simplified tmp;
    tmp.y = i;
    tmp.x = j;
    int key = hash_key(i,j,B,n);
    unordered_map<int, vector<coo_2d_simplified>>::iterator it = mp.find(key);
    if(it == mp.end()) {
      // if not exist, create new vector with this element.
      vector<coo_2d_simplified> vec = {tmp};
      mp.insert(make_pair(key, vec));
    } else {
      // put sample inside.
      vector<coo_2d_simplified> &vec = it->second;
      vec.push_back(tmp);
    }
  }

  time += wall_time_2();
  //printf("  \"time for insert_hashmap\": %.*e,\n", DECIMAL_DIG, time);

  int num_block_per_row = n/B + (n % B == 0 ? 0 : 1);
  
  time = -wall_time_2();
  for (int t = 0; t < s; t++) {
    int ind = samples[t];
    int i = coo[ind].y;
    int j = coo[ind].x;

    //compute x for some i, j
    for (int r = 0; r < W; r++) {
      for (int c = 0; c < W; c++) {
        Z[r][c] = 0;
      }
    }
    // nonzeroinrange

    int i_start = i-B;
    int i_end = i+B - 1;

    int j_start = j-B;
    int j_end = j+B-1;
    int start_block = hash_key(i,j,B,n);
    

    for (int r = start_block-num_block_per_row; r <= start_block + num_block_per_row; r += num_block_per_row) {
      for (int c = -1; c <=1; c++) {
        int b = r+c; // block number
        if(b <= 0) continue;

        unordered_map<int, vector<coo_2d_simplified>>::iterator it = mp.find(b);

        if(it != mp.end()) {

          vector<coo_2d_simplified> &vec = it->second;

          // iterate through all nnz element in the block if it falls in our range.
          for(int k = 0; k < vec.size(); k++) {

            
            if(j_start <= vec[k].x  && vec[k].x <= j_end && i_start <= vec[k].y && vec[k].y <= i_end) {
              
              if(vec[k].y-i_start > 0 && vec[k].x - j_start > 0) Z[vec[k].y-i_start][vec[k].x-j_start] = 1;

            }

          }
        }
      }
      
    }

    i--;j--; 

    for (int r = 1; r < W; r++) {
      for (int c = 1; c < W; c++) {
        Z[r][c] += Z[r][c - 1];
      }
    }

    for (int r = 1; r < W; r++) {
      for (int c = 1; c < W; c++) {
        Z[r][c] += Z[r - 1][c];
      }
    }

    int fill_index = 0;
    for (int b_r = 1; b_r <= B; b_r++) {
      int r_hi = B + b_r - 1 - (i % b_r);
      int r_lo = r_hi - b_r;
      for (int b_c = 1; b_c <= B; b_c++) {
        int c_hi = B + b_c - 1 - (j % b_c);
        int c_lo = c_hi - b_c;
        int y_0 = Z[r_hi][c_hi] - Z[r_lo][c_hi] - Z[r_hi][c_lo] + Z[r_lo][c_lo];
        fill[fill_index] += 1.0/y_0;
        fill_index++;
      }
    }



  }
  time += wall_time_2();
  //printf("  \"time for sampling process\": %.*e,\n", DECIMAL_DIG, time);


  time = -wall_time_2();
  int fill_index = 0;
  for (int b_r = 1; b_r <= B; b_r++) {
    for (int b_c = 1; b_c <= B; b_c++) {
      // TODO: remove b_r * b_c for block ratioi
      fill[fill_index] *= b_r * b_c / (double)s;
      fill_index++;
    }
  }

  time += wall_time_2();
  //printf("  \"time for fill calculation\": %.*e,\n", DECIMAL_DIG, time);


  free(samples);
  return 0;
}

int estimate_fill_coo_2d_variant (int m,
                   int n,
                   int nnz,
                   vector<coo_2d_simplified_2_by_2> &coo,
                   int B,
                   double epsilon,
                   double delta,
                   double *fill,
                   int verbose) {
  int W = 2 * B;
  int Z[W][W];

  double T = 2 * log(B/delta) * B * B / (epsilon * epsilon);
  int s;

double time = -wall_time_2();

#ifdef REPLACEMENT
  s = T;
#else
  //s = ((T - T/nnz) + sqrt((T - T / nnz) * (T - T / nnz)  + 4 * T * (1 + T / N)))/(2 + 2 * T / nnz);
  s = ((T - T/nnz) + sqrt(T * (T + (2 * T + T / nnz) / nnz + 4)))/(2 + 2 * T / nnz);
#endif
  s = min(s, nnz);

  //Sample s items
  int *samples = (int*)malloc(s*sizeof(int));

// Randomized Sampling
#ifdef REPLACEMENT
  if (s == nnz) {
    for (int i = 0; i < nnz; i++) {
      samples[i] = i;
    }
  } else {
    for (int i = 0; i < s; i++) {
      samples[i] = random_range(0, nnz);
    }
  }
#else
  random_choose(samples, s, 0, nnz);
#endif
  sort_int(samples, s);

  time += wall_time_2();
  //printf("  \"time for random sampling\": %.*e,\n", DECIMAL_DIG, time);

  unordered_map<int, vector<coo_2d_simplified_2_by_2>> mp;
  
  // put all elements into hash map

  time = -wall_time_2();

  for (int t = 0; t < nnz; t++) {
    int i = coo[t].local_y;
    int j = coo[t].local_x;
    int key = hash_key(i,j,B,n);
    unordered_map<int, vector<coo_2d_simplified_2_by_2>>::iterator it = mp.find(key);
    if(it == mp.end()) {
      // if not exist, create new vector with this element.
      vector<coo_2d_simplified_2_by_2> vec = {coo[t]};
      mp.insert(make_pair(key, vec));
    } else {
      // put sample inside.
      vector<coo_2d_simplified_2_by_2> &vec = it->second;
      vec.push_back(coo[t]);
    }
  }

  time += wall_time_2();
  //printf("  \"time for insert_hashmap\": %.*e,\n", DECIMAL_DIG, time);

  int num_block_per_row = n/B + (n % B == 0 ? 0 : 1);
  
  time = -wall_time_2();
  for (int t = 0; t < s; t++) {
    int ind = samples[t];
    int i = coo[ind].local_y;
    int j = coo[ind].local_x;

    //compute x for some i, j
    for (int r = 0; r < W; r++) {
      for (int c = 0; c < W; c++) {
        Z[r][c] = 0;
      }
    }
    // nonzeroinrange

    int i_start = i-B;
    int i_end = i+B - 1;

    int j_start = j-B;
    int j_end = j+B-1;
    int start_block = hash_key(i,j,B,n);
    

    for (int r = start_block-num_block_per_row; r <= start_block + num_block_per_row; r += num_block_per_row) {
      for (int c = -1; c <=1; c++) {
        int b = r+c; // block number
        if(b <= 0) continue;

        unordered_map<int, vector<coo_2d_simplified_2_by_2>>::iterator it = mp.find(b);

        if(it != mp.end()) {

          vector<coo_2d_simplified_2_by_2> &vec = it->second;

          // iterate through all nnz element in the block if it falls in our range.
          for(int k = 0; k < vec.size(); k++) {

            if(j_start <= vec[k].local_x  && vec[k].local_x <= j_end && i_start <= vec[k].local_y && vec[k].local_y <= i_end) {
              
              if(vec[k].local_y-i_start > 0 && vec[k].local_x - j_start > 0) Z[vec[k].local_y-i_start][vec[k].local_x-j_start] = 1;

            }

          }
        }
      }
      
    }

    i--;j--; 

    for (int r = 1; r < W; r++) {
      for (int c = 1; c < W; c++) {
        Z[r][c] += Z[r][c - 1];
      }
    }

    for (int r = 1; r < W; r++) {
      for (int c = 1; c < W; c++) {
        Z[r][c] += Z[r - 1][c];
      }
    }

    int fill_index = 0;
    for (int b_r = 1; b_r <= B; b_r++) {
      int r_hi = B + b_r - 1 - (i % b_r);
      int r_lo = r_hi - b_r;
      for (int b_c = 1; b_c <= B; b_c++) {
        int c_hi = B + b_c - 1 - (j % b_c);
        int c_lo = c_hi - b_c;
        int y_0 = Z[r_hi][c_hi] - Z[r_lo][c_hi] - Z[r_hi][c_lo] + Z[r_lo][c_lo];
        fill[fill_index] += 1.0/y_0;
        fill_index++;
      }
    }



  }
  time += wall_time_2();
  //printf("  \"time for sampling process\": %.*e,\n", DECIMAL_DIG, time);


  time = -wall_time_2();
  int fill_index = 0;
  for (int b_r = 1; b_r <= B; b_r++) {
    for (int b_c = 1; b_c <= B; b_c++) {
      // TODO: remove b_r * b_c for block ratioi
      fill[fill_index] *= b_r * b_c  / (double)s;
      fill_index++;
    }
  }

  time += wall_time_2();
  //printf("  \"time for fill calculation\": %.*e,\n", DECIMAL_DIG, time);


  free(samples);
  return 0;
}

int estimate_fill_coo_3d (int x,
                   int y,
                   int z,
                   int nnz,
                   vector<coo_3d_simplified> &coo,
                   int B,
                   double epsilon,
                   double delta,
                   double *fill,
                   int verbose) {
  int W = 2 * B;
  int Z[W][W][W];

  double T = 3 * log(B/delta) * B * B * B / (2 * epsilon * epsilon);
  int s;

double time = -wall_time_2();

#ifdef REPLACEMENT
  s = T;
#else
  //s = ((T - T/nnz) + sqrt((T - T / nnz) * (T - T / nnz)  + 4 * T * (1 + T / N)))/(2 + 2 * T / nnz);
  s = ((T - T/nnz) + sqrt(T * (T + (2 * T + T / nnz) / nnz + 4)))/(2 + 2 * T / nnz);
#endif
  s = min(s, nnz);

//Sample s items
int *samples = (int*)malloc(s*sizeof(int));

// Randomized Sampling
#ifdef REPLACEMENT
  if (s == nnz) {
    for (int i = 0; i < nnz; i++) {
      samples[i] = i;
    }
  } else {
    for (int i = 0; i < s; i++) {
      samples[i] = random_range(0, nnz);
    }
  }
#else
  random_choose(samples, s, 0, nnz);
#endif
  sort_int(samples, s);

  time += wall_time_2();
  printf("  \"time for random sampling\": %.*e,\n", DECIMAL_DIG, time);

  unordered_map<int, vector<coo_3d_simplified>> mp;
  
  // put all elements into hash map

  time = -wall_time_2();

  for (int t = 0; t < nnz; t++) {
    int i = coo[t].x;
    int j = coo[t].y;
    int k = coo[t].z;

    int key = hash_key(i,j,k,B,x,y);
    unordered_map<int, vector<coo_3d_simplified>>::iterator it = mp.find(key);
    if(it == mp.end()) {
      // if not exist, create new vector with this element.
      vector<coo_3d_simplified> vec = {coo[t]};
      printf("%d %d %d\n", vec[0].x, vec[0].y, vec[0].z);
      mp.insert(make_pair(key, vec));
    } else {
      // put sample inside.
      vector<coo_3d_simplified> &vec = it->second;
      vec.push_back(coo[t]);
    }
  }

  time += wall_time_2();
  printf("  \"time for insert_hashmap\": %.*e,\n", DECIMAL_DIG, time);

  int num_block_per_depth = (x/B + (x % B == 0 ? 0 : 1)) * (y/B + (y % B == 0 ? 0 : 1));
  int num_block_per_row = x/B + (x % B == 0 ? 0 : 1);
  
  time = -wall_time_2();
  for (int t = 0; t < s; t++) {
    int ind = samples[t];
    int i = coo[ind].x;
    int j = coo[ind].y;
    int k = coo[ind].z;

    for (int d = 0; d < W; d++) {
      for (int r = 0; r < W; r++) {
        for (int c = 0; c < W; c++) {
          Z[d][r][c] = 0;
        }
      }
    }
    
    // nonzeroinrange

    int i_start = i-B;
    int i_end = i+B - 1;

    int j_start = j-B;
    int j_end = j+B-1;

    int k_start = k-B;
    int k_end = k+B-1;

    int start_block = hash_key(i,j,k,B,x,y);
    
    for (int d = start_block-num_block_per_depth; d <= start_block + num_block_per_depth; d += num_block_per_depth) {
      for (int r = d-num_block_per_row; r <= d + num_block_per_row; r += num_block_per_row) {
        for (int c = r - 1; c <= r + 1; c++) {
          int b = c; // block number
          if(b <= 0) continue;
          
          unordered_map<int, vector<coo_3d_simplified>>::iterator it = mp.find(b);
          

          if(it != mp.end()) {

            vector<coo_3d_simplified> &vec = it->second;

            // iterate through all nnz element in the block if it falls in our range.
            for(int l = 0; l < vec.size(); l++) {

              
              if(j_start <= vec[l].y  && vec[l].y <= j_end && i_start <= vec[l].x && vec[l].x <= i_end && k_start <= vec[l].z && vec[l].z <= k_end) {
                
                if(vec[l].x-i_start > 0 && vec[l].y - j_start > 0 && vec[l].z - k_start > 0) {
                  Z[vec[l].z-k_start][vec[l].y-j_start][vec[l].x-i_start] = 1;
                  //printf("test2\n");
                }
              }

            }
          }
        } 
      }
    }
    
    i--;j--;k--;

    for (int d = 1; d < W; d++) {
      for (int r = 1; r < W; r++) {
        for (int c = 1; c < W; c++) {
          Z[d][r][c] += Z[d][r][c - 1];
        }
      }
    }

    for (int d = 1; d < W; d++) {
      for (int r = 1; r < W; r++) {
        for (int c = 1; c < W; c++) {
          Z[d][r][c] += Z[d][r - 1][c];
        }
      }
    }

    for (int d = 1; d < W; d++) {
      for (int r = 1; r < W; r++) {
        for (int c = 1; c < W; c++) {
          Z[d][r][c] += Z[d - 1][r][c];
        }
      }
    }

    int fill_index = 0;
    for (int b_d = 1; b_d <= B; b_d++) {
        int d_hi = B + b_d - 1 - (k % b_d);
        int d_lo = d_hi - b_d;
      for (int b_r = 1; b_r <= B; b_r++) {
        int r_hi = B + b_r - 1 - (j % b_r);
        int r_lo = r_hi - b_r;
        for (int b_c = 1; b_c <= B; b_c++) {
          int c_hi = B + b_c - 1 - (i % b_c);
          int c_lo = c_hi - b_c;
          int y_0 = Z[d_hi][r_hi][c_hi] - Z[d_hi][r_hi][c_lo] - Z[d_hi][r_lo][c_hi] - Z[d_lo][r_hi][c_hi] + Z[d_hi][r_lo][c_lo] + Z[d_lo][r_lo][c_hi] + Z[d_lo][r_hi][c_lo] - Z[d_lo][r_lo][c_lo];
          //if(y_0 ==0) continue;
          //printf("y_0 = %d / %d - %d - %d - %d + %d + %d + %d - %d\n", y_0, Z[d_hi][r_hi][c_hi], Z[d_hi][r_hi][c_lo], Z[d_hi][r_lo][c_hi],  Z[d_lo][r_hi][c_hi], Z[d_hi][r_lo][c_lo], Z[d_lo][r_lo][c_hi], Z[d_lo][r_hi][c_lo], Z[d_lo][r_lo][c_lo]);
          fill[fill_index] += 1.0/y_0;
          fill_index++;
        }
      }
    }



  }
  time += wall_time_2();
  printf("  \"time for sampling process\": %.*e,\n", DECIMAL_DIG, time);


  time = -wall_time_2();
  int fill_index = 0;
  for (int b_d = 1; b_d <= B; b_d++) {
    for (int b_r = 1; b_r <= B; b_r++) {
      for (int b_c = 1; b_c <= B; b_c++) {
        fill[fill_index] *= /*b_d * b_r * b_c */ 1 / (double)s;
        fill_index++;
      }
    }
  }

  time += wall_time_2();
  printf("  \"time for fill calculation\": %.*e,\n", DECIMAL_DIG, time);


  free(samples);
  return 0;
}

}


