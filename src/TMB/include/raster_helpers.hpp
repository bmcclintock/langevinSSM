#ifndef RASTER_HELPERS_HPP
#define RASTER_HELPERS_HPP

// --- TYPE MAPPING ---
#ifdef IS_RCPP_BUILD
#define TEMPLATE_HEADER inline
#define TYPE double
#define MATRIX arma::mat
#define VECTOR Rcpp::NumericVector
#define INT_VECTOR Rcpp::IntegerVector
#define RASTER_TYPE Rcpp::NumericVector // Renamed to avoid TMB macro collision
#define AS_DOUBLE(x) (x)
#define PI_VAL M_PI
#define SET_ZERO(x) x.zeros()
// Rcpp uses standard bracket indexing
#define GET_VAL(arr, idx) arr[idx]
#define VEC_ELT(v, i) v[i]
#else
#define TEMPLATE_HEADER template<class Type>
#define TYPE Type
#define MATRIX matrix<Type>
#define VECTOR vector<Type>
#define INT_VECTOR vector<int>
#define RASTER_TYPE array<Type> // Renamed to avoid TMB macro collision
// asDouble explicitly strips AD tracking to safely evaluate integer indices
#define AS_DOUBLE(x) asDouble(x)
#define PI_VAL Type(M_PI)
#define SET_ZERO(x) x.setZero()
// TMB (Eigen) uses parentheses for array/vector indexing
#define GET_VAL(arr, idx) arr(idx)
#define VEC_ELT(v, i) v(i)
#endif


TEMPLATE_HEADER
MATRIX extract_raster_values(TYPE x, TYPE y, TYPE z,
                             RASTER_TYPE &raster_vals,
                             const VECTOR &all_z_values,
                             const INT_VECTOR &n_zvals_cov,
                             const INT_VECTOR &cov_offset,
                             const MATRIX &raster_coords,
                             const VECTOR &raster_resolution,
                             const VECTOR &raster_extent,
                             int n_covs) {

  TYPE res_x = VEC_ELT(raster_resolution, 0);
  TYPE res_y = VEC_ELT(raster_resolution, 1);

  TYPE col_raw = (x - VEC_ELT(raster_extent, 0)) / res_x - TYPE(0.5);
  TYPE row_raw = (VEC_ELT(raster_extent, 3) - y) / res_y - TYPE(0.5);

#ifdef IS_RCPP_BUILD
  int n_cols = static_cast<int>(round(AS_DOUBLE((VEC_ELT(raster_extent, 1) - VEC_ELT(raster_extent, 0)) / res_x)));
  int n_rows = static_cast<int>(round(AS_DOUBLE((VEC_ELT(raster_extent, 3) - VEC_ELT(raster_extent, 2)) / res_y)));
#else
  // TMB uses integer dimensions provided via raster_coords if needed, or derived
  int n_cols = static_cast<int>(round(AS_DOUBLE((VEC_ELT(raster_extent, 1) - VEC_ELT(raster_extent, 0)) / res_x)));
  int n_rows = static_cast<int>(round(AS_DOUBLE((VEC_ELT(raster_extent, 3) - VEC_ELT(raster_extent, 2)) / res_y)));
#endif

  int c0 = static_cast<int>(floor(AS_DOUBLE(col_raw)));
  int r0 = static_cast<int>(floor(AS_DOUBLE(row_raw)));

  MATRIX grad_values(2, n_covs);
  SET_ZERO(grad_values);

  if(c0 < 0 || c0 >= (n_cols-1) || r0 < 0 || r0 >= (n_rows-1)) return grad_values;

  TYPE dx = col_raw - TYPE(c0);
  TYPE dy = row_raw - TYPE(r0);
  int rev_row = n_rows - 1 - r0;

  for(int i = 0; i < n_covs; i++) {
    int layers = VEC_ELT(n_zvals_cov, i);
    int offset = VEC_ELT(cov_offset, i);
    int z_idx1 = 0, z_idx2 = 0;
    TYPE z_weight = 0.0;

    if (layers > 1) {
      if (z <= VEC_ELT(all_z_values, offset)) {
        z_idx1 = 0; z_idx2 = 0; z_weight = 0.0;
      } else if (z >= VEC_ELT(all_z_values, offset + layers - 1)) {
        z_idx1 = layers - 1; z_idx2 = layers - 1; z_weight = 0.0;
      } else {
        for (int k = 0; k < layers - 1; k++) {
          TYPE t1 = VEC_ELT(all_z_values, offset + k);
          TYPE t2 = VEC_ELT(all_z_values, offset + k + 1);
          if (z >= t1 && z <= t2) {
            z_idx1 = k; z_idx2 = k + 1;
            z_weight = (z - t1) / (t2 - t1);
            break;
          }
        }
      }
    }

    int idx1 = (offset + z_idx1) * (n_rows * n_cols) + rev_row * n_cols + c0;
    TYPE f1_00 = GET_VAL(raster_vals, idx1);
    TYPE f1_10 = GET_VAL(raster_vals, idx1 + 1);
    TYPE f1_01 = GET_VAL(raster_vals, idx1 - n_cols);
    TYPE f1_11 = GET_VAL(raster_vals, idx1 - n_cols + 1);

    // Optimized Analytical Gradient using normalized step lengths
    TYPE grad_x_1 = ((f1_10 - f1_00) * (TYPE(1.0) - dy) + (f1_11 - f1_01) * dy) / res_x;
    TYPE grad_y_1 = ((f1_00 - f1_01) * (TYPE(1.0) - dx) + (f1_10 - f1_11) * dx) / res_y;

    if (layers > 1 && z_idx1 != z_idx2) {
      int idx2 = (offset + z_idx2) * (n_rows * n_cols) + rev_row * n_cols + c0;
      TYPE f2_00 = GET_VAL(raster_vals, idx2);
      TYPE f2_10 = GET_VAL(raster_vals, idx2 + 1);
      TYPE f2_01 = GET_VAL(raster_vals, idx2 - n_cols);
      TYPE f2_11 = GET_VAL(raster_vals, idx2 - n_cols + 1);

      TYPE grad_x_2 = ((f2_10 - f2_00) * (TYPE(1.0) - dy) + (f2_11 - f2_01) * dy) / res_x;
      TYPE grad_y_2 = ((f2_00 - f2_01) * (TYPE(1.0) - dx) + (f2_10 - f2_11) * dx) / res_y;

      grad_values(0,i) = grad_x_1 * (TYPE(1.0) - z_weight) + grad_x_2 * z_weight;
      grad_values(1,i) = grad_y_1 * (TYPE(1.0) - z_weight) + grad_y_2 * z_weight;
    } else {
      grad_values(0,i) = grad_x_1;
      grad_values(1,i) = grad_y_1;
    }
  }
  return grad_values;
}

TEMPLATE_HEADER
MATRIX calculate_smooth_gradient(TYPE x, TYPE y, TYPE z,
                                 RASTER_TYPE &raster_vals, // Removed 'const' to satisfy TMB
                                 const VECTOR &all_z_values,
                                 const INT_VECTOR &n_zvals_cov,
                                 const INT_VECTOR &cov_offset,
                                 const MATRIX &raster_coords,
                                 const VECTOR &raster_resolution,
                                 const VECTOR &raster_extent,
                                 int n_covs,
                                 const VECTOR &weights,
                                 TYPE sigma,
                                 TYPE gamma,
                                 TYPE dt_step,
                                 TYPE zetaScale,
                                 int model) {

  TYPE zeta;
  if(model == 1) {
    zeta = zetaScale * sqrt(TYPE(2.0) * PI_VAL) * sigma;
  } else {
    zeta = zetaScale * sigma / sqrt(TYPE(2.0));
  }

  TYPE neighborhood_size = zeta * sqrt(dt_step);

  MATRIX smooth_grads = VEC_ELT(weights, 0) * extract_raster_values(x, y, z, raster_vals,
                                all_z_values, n_zvals_cov, cov_offset,
                                raster_coords, raster_resolution,
                                raster_extent, n_covs);

  int n_points = weights.size();

  for(int i = 0; i < (n_points-1); i++) {
    TYPE angle;
    if(n_points == 9) {
      angle = PI_VAL * (TYPE(2.0) * i + TYPE(1.0)) / TYPE(8.0);
    } else {
      angle = PI_VAL * (TYPE(2.0) * i + TYPE(1.0)) / TYPE(4.0);
    }

    smooth_grads += VEC_ELT(weights, i+1) * extract_raster_values(
      x + neighborhood_size * cos(angle),
      y + neighborhood_size * sin(angle),
      z,
      raster_vals,
      all_z_values,
      n_zvals_cov,
      cov_offset,
      raster_coords,
      raster_resolution,
      raster_extent,
      n_covs);
  }

  return smooth_grads;
}

TEMPLATE_HEADER
TYPE get_bilinear_val(TYPE x, TYPE y, const MATRIX &grid, const VECTOR &ext, const VECTOR &res) {

#ifdef IS_RCPP_BUILD
  int n_cols = grid.n_cols;
  int n_rows = grid.n_rows;
#else
  int n_cols = grid.cols();
  int n_rows = grid.rows();
#endif

  TYPE col_raw = (x - VEC_ELT(ext, 0)) / VEC_ELT(res, 0) - TYPE(0.5);
  TYPE row_raw = (VEC_ELT(ext, 3) - y) / VEC_ELT(res, 1) - TYPE(0.5);

  int c0 = static_cast<int>(floor(AS_DOUBLE(col_raw)));
  int r0 = static_cast<int>(floor(AS_DOUBLE(row_raw)));

  if (c0 < 0) c0 = 0; if (c0 >= n_cols - 1) c0 = n_cols - 2;
  if (r0 < 0) r0 = 0; if (r0 >= n_rows - 1) r0 = n_rows - 2;
  int c1 = c0 + 1;
  int r1 = r0 + 1;

  TYPE dx = col_raw - TYPE(c0);
  TYPE dy = row_raw - TYPE(r0);

  if(dx < TYPE(0.0)) dx = TYPE(0.0); if(dx > TYPE(1.0)) dx = TYPE(1.0);
  if(dy < TYPE(0.0)) dy = TYPE(0.0); if(dy > TYPE(1.0)) dy = TYPE(1.0);

  TYPE val = grid(r0, c0) * (TYPE(1.0) - dx) * (TYPE(1.0) - dy) +
    grid(r0, c1) * dx * (TYPE(1.0) - dy) +
    grid(r1, c0) * (TYPE(1.0) - dx) * dy +
    grid(r1, c1) * dx * dy;
  return val;
}

TEMPLATE_HEADER
void apply_barrier_penalty(TYPE x, TYPE y,
                           const MATRIX &barrier_dist,
                           const VECTOR &raster_extent,
                           const VECTOR &raster_resolution,
                           TYPE barrier_penalty,
                           TYPE &h_x, TYPE &h_y) {

  if (barrier_penalty > TYPE(0.0)) {
#ifdef IS_RCPP_BUILD
    int n_cols = barrier_dist.n_cols;
    int n_rows = barrier_dist.n_rows;
#else
    int n_cols = barrier_dist.cols();
    int n_rows = barrier_dist.rows();
#endif

    TYPE res_x = VEC_ELT(raster_resolution, 0);
    TYPE res_y = VEC_ELT(raster_resolution, 1);

    TYPE col_raw = (x - VEC_ELT(raster_extent, 0)) / res_x - TYPE(0.5);
    TYPE row_raw = (VEC_ELT(raster_extent, 3) - y) / res_y - TYPE(0.5);

    int c0 = static_cast<int>(floor(AS_DOUBLE(col_raw)));
    int r0 = static_cast<int>(floor(AS_DOUBLE(row_raw)));

    if (c0 < 0) c0 = 0; if (c0 >= n_cols - 1) c0 = n_cols - 2;
    if (r0 < 0) r0 = 0; if (r0 >= n_rows - 1) r0 = n_rows - 2;
    int c1 = c0 + 1;
    int r1 = r0 + 1;

    TYPE dx = col_raw - TYPE(c0);
    TYPE dy = row_raw - TYPE(r0);
    if(dx < TYPE(0.0)) dx = TYPE(0.0); if(dx > TYPE(1.0)) dx = TYPE(1.0);
    if(dy < TYPE(0.0)) dy = TYPE(0.0); if(dy > TYPE(1.0)) dy = TYPE(1.0);

    TYPE f00 = barrier_dist(r0, c0);
    TYPE f10 = barrier_dist(r0, c1);
    TYPE f01 = barrier_dist(r1, c0);
    TYPE f11 = barrier_dist(r1, c1);

    // Calculate distance exactly
    TYPE d_barrier = f00 * (TYPE(1.0) - dx) * (TYPE(1.0) - dy) +
      f10 * dx * (TYPE(1.0) - dy) +
      f01 * (TYPE(1.0) - dx) * dy +
      f11 * dx * dy;

    if (d_barrier <= TYPE(0.0)) {
      // Analytical gradients of the bilinear surface (no more 'eps' approximations)
      TYPE grad_sdf_x = ((f10 - f00) * (TYPE(1.0) - dy) + (f11 - f01) * dy) / res_x;
      TYPE grad_sdf_y = ((f00 - f01) * (TYPE(1.0) - dx) + (f10 - f11) * dx) / res_y;

      // Apply force = -lambda * d * grad(d)
      h_x -= barrier_penalty * d_barrier * grad_sdf_x;
      h_y -= barrier_penalty * d_barrier * grad_sdf_y;
    }
  }
}

#endif // RASTER_HELPERS_HPP
