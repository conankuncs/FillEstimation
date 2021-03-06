#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unordered_map>
#include <vector>

static inline int hash_key(int i,int j,int B, int n) {
  return (1 + ((int)(n/B) + (n % (int)B == 0 ? 0 : 1)) * (int)((i-1)/B) + (int)((j-1)/B));
}

using namespace std;

extern "C" {

#include "util.h"

char *name () {
  return (char *)"asx";
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

  //printf("Samples: %d (%d)\n", s, nnz);

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
        //printf("%d %d\n", r, c);
        //if(t==2) printf("+");
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

    if(t==100)
    {printf("after (%d,%d):\n", i, j);
    for (int r = 1; r < W; r++) {
      for (int c = 1; c < W; c++) {
        printf("%d ", Z[r][c]);
      }
      printf("\n");
    }}

    int fill_index = 0;
    for (int b_r = 1; b_r <= B; b_r++) {
      int r_hi = B + b_r - 1 - (i % b_r);
      int r_lo = r_hi - b_r;
      for (int b_c = 1; b_c <= B; b_c++) {
        int c_hi = B + b_c - 1 - (j % b_c);
        int c_lo = c_hi - b_c;
        int y_0 = Z[r_hi][c_hi] - Z[r_lo][c_hi] - Z[r_hi][c_lo] + Z[r_lo][c_lo];
        //printf("%d %d %d %d\n", Z[r_hi][c_hi], Z[r_lo][c_hi] , Z[r_hi][c_lo] , Z[r_lo][c_lo]);
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

#ifdef REPLACEMENT
  s = T;
#else
  //s = ((T - T/nnz) + sqrt((T - T / nnz) * (T - T / nnz)  + 4 * T * (1 + T / N)))/(2 + 2 * T / nnz);
  s = ((T - T/nnz) + sqrt(T * (T + (2 * T + T / nnz) / nnz + 4)))/(2 + 2 * T / nnz);
#endif
  s = min(s, nnz);

  //Sample s items
  int *samples = (int*)malloc(s*sizeof(int));
  int *samples_i = (int*)malloc(s*sizeof(int));
  int *samples_j = (int*)malloc(s*sizeof(int));

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
  unordered_map<int, vector<coo_2d_simplified>> mp;
  // put into hash map
  for (int t = 0; t < s; t++) {

    int ind = samples[t];
    int i = coo[ind].y;
    int j = coo[ind].x;
    coo_2d_simplified tmp;
    tmp.y = i;
    tmp.x = j;
    int key = hash_key(i,j,B,n);
    unordered_map<int, vector<coo_2d_simplified>>::iterator it = mp.find(key);
    if(it == mp.end()) {
      // if not exist, put empty vector
      vector<coo_2d_simplified> vec = {tmp};
      mp.insert(make_pair(key, vec));

      // insert
    /* it = mp.find(key);
      vector<coo_2d_simplified> &vec = it->second;
      vec.push_back(tmp);*/
    } else {
      // put sample inside.
      vector<coo_2d_simplified> &vec = it->second;
      vec.push_back(tmp);
    }
  }

  int num_block_per_row = n/B + (n % B == 0 ? 0 : 1);
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

    //int i_start = (i-B) > 1 ? i - B : 1;
    //int i_end = (i+B-1) < m ? i+B-1 : m;
    int i_start = i-B;
    int i_end = i+B - 1;
    // int j_start = (j-B) > 1 ? j-B : 1;
    // int j_end = (j+B-1) < n ? j+B-1 : n;
    int j_start = j-B;
    int j_end = j+B-1;
    int start_block = hash_key(i,j,B,n);
    

    /*printf("before:\n");
    for (int r = 1; r < W; r++) {
      for (int c = 1; c < W; c++) {
        printf("%d ", Z[r][c]);
      }
      printf("\n");
    }*/

    for (int r = start_block-num_block_per_row; r <= start_block + num_block_per_row; r += num_block_per_row) {
      
      for (int c = -1; c <=1; c++) {
        int b = r+c; // block number
        if(b < 0) continue;
        // find block in the hash map
        //printf("test3: %d %d %d-- %d %d\n", i,j, B, r,c);fflush(stdout);
        unordered_map<int, vector<coo_2d_simplified>>::iterator it = mp.find(b);
        //printf("test4 0x%x 0x%x\n", it, mp.end());fflush(stdout);*
        //printf("range: (%d,%d) / %d ~ %d , %d ~ %d\n",i,j, i_start, i_end, j_start, j_end);
        if(it != mp.end()) {

          vector<coo_2d_simplified> &vec = it->second;

          // iterate through all nnz element in the block if it falls in our range.
          for(int k = 0; k < vec.size(); k++) {

            //if(j_start <= vec[k].x  && vec[k].x <= j_end && i_start <= vec[k].y && vec[k].y <= i_end) {
              // If an element in the block is inside I-B+1 ~ I + B -1, count it
              //printf("%d %d (%d~%d, %d~%d) => %d %d\n", vec[k].y, vec[k].x, i_start, i_end,j_start, j_end,vec[k].y-i_start, vec[k].x-j_start);fflush(stdout);
              //Z[vec[k].y-i_start+1][vec[k].x-j_start+1] = 1;
              //if(t==2) printf("+");
            //}
            //printf("%d %d (%d) / %d~%d, %d~%d ", vec[k].y, vec[k].x, b, i_start, i_end, j_start, j_end);
            if(j_start < vec[k].x  && vec[k].x < j_end && i_start < vec[k].y && vec[k].y < i_end) {
              // If an element in the block is inside I-B+1 ~ I + B -1, count it
              //printf("%d %d (%d~%d, %d~%d) => %d %d\n", vec[k].y, vec[k].x, i_start, i_end,j_start, j_end,vec[k].y-i_start, vec[k].x-j_start);fflush(stdout);

              if(vec[k].y-i_start == 0 || vec[k].x-j_start == 0) {
                printf("Test\n");
              }
              Z[vec[k].y-i_start][vec[k].x-j_start] = 1;


              //printf("%d %d (%d)\n", vec[k].y-i_start+1, vec[k].x-j_start+1, b);
              //\printf(" (%d %d)", vec[k].y-i_start+1, vec[k].x-j_start+1);
              //printf(" --counted");
              //if(t==2) printf("+");
            }
          //  printf("\n");
          }
        }
      }
      printf("reached here\n");
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
if(t==100)
    {printf("after (%d,%d):\n", i, j);
    for (int r = 1; r < W; r++) {
      for (int c = 1; c < W; c++) {
        printf("%d ", Z[r][c]);
      }
      printf("\n");
    }}

    i--;j--;

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

}
