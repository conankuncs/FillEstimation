#ifndef max
    #define max(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef min
    #define min(a,b) ((a) < (b) ? (a) : (b))
#endif

//lo is inclusive, hi is exclusive

typedef char hash_t;
hash_t *hash_create ();
void hash_set (hash_t *table, size_t key, size_t value);
size_t hash_get (hash_t *table, size_t key, size_t value);
void hash_clear (hash_t *table);
void hash_destroy (hash_t *table);
size_t random_range (size_t lo, size_t hi);
double random_uniform ();
void random_choose (size_t *samples, size_t n, size_t lo, size_t hi);
void sort (size_t *stuff, size_t n);
void sort_int (int *stuff, int n);
size_t search (const size_t *stuff, size_t lo, size_t hi, size_t key);
size_t search_strict (const size_t *stuff, size_t lo, size_t hi, size_t key);

// TODO: Change name later
typedef struct coo_2d {
  int x,y;
  double val;
} coo_2d;

typedef struct coo_2d_simplified {
  int x,y;
} coo_2d_simplified;

typedef struct coo_2d_simplified_2_by_2 {
	coo_2d_simplified arr[5];
	int index;
	int local_x, local_y;
} coo_2d_simplified_2_by_2;

typedef struct coo_3d_simplified {
	int x,y,z;
} coo_3d_simplified;