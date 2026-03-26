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
                             RASTER_TYPE &raster_vals, // Removed 'const' to satisfy TMB
                             const VECTOR &all_z_values,
                             const INT_VECTOR &n_zvals_cov,
                             const INT_VECTOR &cov_offset,
                             const MATRIX &raster_coords,
                             const VECTOR &raster_resolution,
                             const VECTOR &raster_extent,
                             int n_covs) {

  // Use round() to prevent floating point inaccuracy truncation
  int n_cols = static_cast<int>(round(AS_DOUBLE((VEC_ELT(raster_extent, 1) - VEC_ELT(raster_extent, 0)) / VEC_ELT(raster_resolution, 0))));
  int n_rows = static_cast<int>(round(AS_DOUBLE((VEC_ELT(raster_extent, 3) - VEC_ELT(raster_extent, 2)) / VEC_ELT(raster_resolution, 1))));

  TYPE x_prop = (x - (VEC_ELT(raster_extent, 0) + VEC_ELT(raster_resolution, 0)/TYPE(2.0))) / VEC_ELT(raster_resolution, 0);
  TYPE y_prop = (y - (VEC_ELT(raster_extent, 2) + VEC_ELT(raster_resolution, 1)/TYPE(2.0))) / VEC_ELT(raster_resolution, 1);

  // CRITICAL FIX: Strip AD type using AS_DOUBLE before calling floor()
  // This calculates the cell index without crashing the CppAD tape
  int col = static_cast<int>(floor(AS_DOUBLE(x_prop)));
  int row = static_cast<int>(floor(AS_DOUBLE(y_prop)));

  MATRIX grad_values(2, n_covs);
  SET_ZERO(grad_values);

  if(col < 0 || col >= (n_cols-1) || row < 0 || row >= (n_rows-1)) return grad_values;

  TYPE x1 = VEC_ELT(raster_extent, 0) + (TYPE(col) + TYPE(0.5)) * VEC_ELT(raster_resolution, 0);
  TYPE x2 = x1 + VEC_ELT(raster_resolution, 0);
  TYPE y1 = VEC_ELT(raster_extent, 2) + (TYPE(row) + TYPE(0.5)) * VEC_ELT(raster_resolution, 1);
  TYPE y2 = y1 + VEC_ELT(raster_resolution, 1);

  int rev_row = n_rows - 1 - row;

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
            z_idx1 = k;
            z_idx2 = k + 1;
            z_weight = (z - t1) / (t2 - t1);
            break;
          }
        }
      }
    }

    int layer_idx1 = offset + z_idx1;
    int idx1 = layer_idx1 * (n_rows * n_cols) + rev_row * n_cols + col;

    TYPE f1_00 = GET_VAL(raster_vals, idx1);
    TYPE f1_10 = GET_VAL(raster_vals, idx1 + 1);
    TYPE f1_01 = GET_VAL(raster_vals, idx1 - n_cols);
    TYPE f1_11 = GET_VAL(raster_vals, idx1 - n_cols + 1);

    TYPE grad_x_1 = ((y2 - y) * (f1_10 - f1_00) + (y - y1) * (f1_11 - f1_01)) / ((y2 - y1) * (x2 - x1));
    TYPE grad_y_1 = ((x2 - x) * (f1_01 - f1_00) + (x - x1) * (f1_11 - f1_10)) / ((y2 - y1) * (x2 - x1));

    if (layers > 1 && z_idx1 != z_idx2) {
      int layer_idx2 = offset + z_idx2;
      int idx2 = layer_idx2 * (n_rows * n_cols) + rev_row * n_cols + col;

      TYPE f2_00 = GET_VAL(raster_vals, idx2);
      TYPE f2_10 = GET_VAL(raster_vals, idx2 + 1);
      TYPE f2_01 = GET_VAL(raster_vals, idx2 - n_cols);
      TYPE f2_11 = GET_VAL(raster_vals, idx2 - n_cols + 1);

      TYPE grad_x_2 = ((y2 - y) * (f2_10 - f2_00) + (y - y1) * (f2_11 - f2_01)) / ((y2 - y1) * (x2 - x1));
      TYPE grad_y_2 = ((x2 - x) * (f2_01 - f2_00) + (x - x1) * (f2_11 - f2_10)) / ((y2 - y1) * (x2 - x1));

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

#endif // RASTER_HELPERS_HPP
