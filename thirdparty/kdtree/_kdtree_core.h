typedef struct
{
    float cut_val;
    int8_t cut_dim;
    uint32_t start_idx;
    uint32_t n;
    float cut_bounds_lv;
    float cut_bounds_hv;
    struct Node_float *left_child;
    struct Node_float *right_child;
} Node_float;

typedef struct
{
    float *bbox;
    int8_t no_dims;
    uint32_t *pidx;
    struct Node_float *root; 
} Tree_float;


typedef struct
{
    double cut_val;
    int8_t cut_dim;
    uint32_t start_idx;
    uint32_t n;
    double cut_bounds_lv;
    double cut_bounds_hv;
    struct Node_double *left_child;
    struct Node_double *right_child;
} Node_double;

typedef struct
{
    double *bbox;
    int8_t no_dims;
    uint32_t *pidx;
    struct Node_double *root; 
} Tree_double;



void insert_point_float(uint32_t *closest_idx, float *closest_dist, uint32_t pidx, float cur_dist, uint32_t k);
void get_bounding_box_float(float *pa, uint32_t *pidx, int8_t no_dims, uint32_t n, float *bbox);
int partition_float(float *pa, uint32_t *pidx, int8_t no_dims, uint32_t start_idx, uint32_t n, float *bbox, int8_t *cut_dim, 
              float *cut_val, uint32_t *n_lo);
Tree_float* construct_tree_float(float *pa, int8_t no_dims, uint32_t n, uint32_t bsp);
Node_float* construct_subtree_float(float *pa, uint32_t *pidx, int8_t no_dims, uint32_t start_idx, uint32_t n, uint32_t bsp, float *bbox);
Node_float * create_node_float(uint32_t start_idx, uint32_t n, int is_leaf);
void delete_subtree_float(Node_float *root);
void delete_tree_float(Tree_float *tree);
void print_tree_float(Node_float *root, int level);
float calc_dist_float(float *point1_coord, float *point2_coord, int8_t no_dims);
float get_cube_offset_float(int8_t dim, float *point_coord, float *bbox);
float get_min_dist_float(float *point_coord, int8_t no_dims, float *bbox);
void search_leaf_float(float *restrict pa, uint32_t *restrict pidx, int8_t no_dims, uint32_t start_idx, uint32_t n, float *restrict point_coord, 
                 uint32_t k, uint32_t *restrict closest_idx, float *restrict closest_dist);
void search_splitnode_float(Node_float *root, float *pa, uint32_t *pidx, int8_t no_dims, float *point_coord, 
                      float min_dist, uint32_t k, float distance_upper_bound, float eps_fac, uint32_t *  closest_idx, float *closest_dist);
void search_tree_float(Tree_float *tree, float *pa, float *point_coords, 
                 uint32_t num_points, uint32_t k,  float distance_upper_bound, 
                 float eps, uint32_t *closest_idxs, float *closest_dists);


void insert_point_double(uint32_t *closest_idx, double *closest_dist, uint32_t pidx, double cur_dist, uint32_t k);
void get_bounding_box_double(double *pa, uint32_t *pidx, int8_t no_dims, uint32_t n, double *bbox);
int partition_double(double *pa, uint32_t *pidx, int8_t no_dims, uint32_t start_idx, uint32_t n, double *bbox, int8_t *cut_dim, 
              double *cut_val, uint32_t *n_lo);
Tree_double* construct_tree_double(double *pa, int8_t no_dims, uint32_t n, uint32_t bsp);
Node_double* construct_subtree_double(double *pa, uint32_t *pidx, int8_t no_dims, uint32_t start_idx, uint32_t n, uint32_t bsp, double *bbox);
Node_double * create_node_double(uint32_t start_idx, uint32_t n, int is_leaf);
void delete_subtree_double(Node_double *root);
void delete_tree_double(Tree_double *tree);
void print_tree_double(Node_double *root, int level);
double calc_dist_double(double *point1_coord, double *point2_coord, int8_t no_dims);
double get_cube_offset_double(int8_t dim, double *point_coord, double *bbox);
double get_min_dist_double(double *point_coord, int8_t no_dims, double *bbox);
void search_leaf_double(double *restrict pa, uint32_t *restrict pidx, int8_t no_dims, uint32_t start_idx, uint32_t n, double *restrict point_coord, 
                 uint32_t k, uint32_t *restrict closest_idx, double *restrict closest_dist);
void search_splitnode_double(Node_double *root, double *pa, uint32_t *pidx, int8_t no_dims, double *point_coord, 
                      double min_dist, uint32_t k, double distance_upper_bound, double eps_fac, uint32_t *  closest_idx, double *closest_dist);
void search_tree_double(Tree_double *tree, double *pa, double *point_coords, 
                 uint32_t num_points, uint32_t k,  double distance_upper_bound, 
                 double eps, uint32_t *closest_idxs, double *closest_dists);

