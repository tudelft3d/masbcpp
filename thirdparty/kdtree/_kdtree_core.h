
extern "C"
{
typedef struct
{
    float *bbox;
    int8_t no_dims;
    uint32_t *pidx;
    struct Node_float *root; 
} Tree_float;

typedef struct
{
    double *bbox;
    int8_t no_dims;
    uint32_t *pidx;
    struct Node_double *root; 
} Tree_double;

Tree_float* construct_tree_float(float *pa, int8_t no_dims, uint32_t n, uint32_t bsp);
void search_tree_float(Tree_float *tree, float *pa, float *point_coords, 
                 uint32_t num_points, uint32_t k,  float distance_upper_bound, 
                 float eps, uint32_t *closest_idxs, float *closest_dists);


Tree_double* construct_tree_double(double *pa, int8_t no_dims, uint32_t n, uint32_t bsp);
void search_tree_double(Tree_double *tree, double *pa, double *point_coords, 
                 uint32_t num_points, uint32_t k,  double distance_upper_bound, 
                 double eps, uint32_t *closest_idxs, double *closest_dists);

}