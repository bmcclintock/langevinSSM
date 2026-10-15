#ifndef RASTER_HELPERS_HPP
#define RASTER_HELPERS_HPP

// ---------------------------------------------------------------------------
// Raster covariate surfaces for langevinSSM (TMB template) and simulateLangevin (Rcpp).
//
// Each covariate layer (and the barrier signed-distance grid) is represented by an
// INTERPOLATING QUINTIC B-SPLINE through the cell-centre values:
//   * it reproduces every cell-centre value exactly;
//   * it is C^4 continuous, so the habitat gradient (force) and its derivatives up to third
//     order are continuous across cell edges.  The inner Hessian needs the 2nd derivatives of
//     the surface, its derivative (used by TMB's Laplace gradient) the 3rd, and continuity of
//     that gradient the 4th.  Bilinear interpolation has a gradient that jumps at every cell
//     edge, which made the inner (latent-path) objective non-smooth.
//
// In the TMB build, surface evaluation is wrapped in a TMBad ATOMIC operator (RasterJetOp).
// The cell lookup happens inside the operator's forward pass, so it is redone at every
// evaluation of the tape.
//
// Derivatives of any order are supplied by the operator itself (reverse pass returns the
// next-order derivative operator), mirroring TMB's own tmbutils::interpol2D.
//
// Performance: one RasterJetOp per location returns ALL partials needed of a given order
// for ALL covariate layers (one cell lookup, one set of basis weights, a 6x6 window per layer,
// mixed partials computed once).  Its adjoint is a single jet operator of the next order.
// ---------------------------------------------------------------------------

#include <vector>
#include <memory>
#include <cmath>

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
#define MAT_NROW(m) ((int) (m).n_rows)
#define MAT_NCOL(m) ((int) (m).n_cols)
#else
#define TEMPLATE_HEADER template<class Type>
#define TYPE Type
#define MATRIX matrix<Type>
#define VECTOR vector<Type>
#define INT_VECTOR vector<int>
#define RASTER_TYPE array<Type> // Renamed to avoid TMB macro collision
// asDouble explicitly strips AD tracking (only used on DATA, never on parameters)
#define AS_DOUBLE(x) asDouble(x)
#define PI_VAL Type(M_PI)
#define SET_ZERO(x) x.setZero()
// TMB (Eigen) uses parentheses for array/vector indexing
#define GET_VAL(arr, idx) arr(idx)
#define VEC_ELT(v, i) v(i)
#define MAT_NROW(m) ((int) (m).rows())
#define MAT_NCOL(m) ((int) (m).cols())
#endif

namespace langevin_raster {

// ---------------------------------------------------------------------------
// B-spline bases of degree 1 (bilinear), 3 (bicubic) and 5 (quintic).
// For degree p the surface on a cell is a tensor product of degree-p polynomials built from a
// W x W window of coefficients, W = p + 1, starting LO = -(p - 1) / 2 cells from the cell's
// lower-left index: p = 1: cells i0..i0+1;  p = 3: i0-1..i0+2;  p = 5: i0-2..i0+3.
// Degree 1 with no prefilter is exactly bilinear interpolation of the cell-centre values.
// ---------------------------------------------------------------------------
static const int SPLINE_WMAX = 6;

// POLY tables: coefficient of t^d in beta_p(t - LO - k), k = 0..p (t in [0,1)).
inline const double* bspline_poly(int deg) {
  static const double P1[2][2] = { { 1.0, -1.0 }, { 0.0, 1.0 } };
  static const double P3[4][4] = {
    { 1.0/6, -3.0/6,  3.0/6, -1.0/6 },
    { 4.0/6,  0.0,   -6.0/6,  3.0/6 },
    { 1.0/6,  3.0/6,  3.0/6, -3.0/6 },
    { 0.0,    0.0,    0.0,    1.0/6 }
  };
  static const double P5[6][6] = {
    {  1.0/120,  -5.0/120,  10.0/120, -10.0/120,   5.0/120,  -1.0/120 },
    { 26.0/120, -50.0/120,  20.0/120,  20.0/120, -20.0/120,   5.0/120 },
    { 66.0/120,   0.0,     -60.0/120,   0.0,      30.0/120, -10.0/120 },
    { 26.0/120,  50.0/120,  20.0/120, -20.0/120, -20.0/120,  10.0/120 },
    {  1.0/120,   5.0/120,  10.0/120,  10.0/120,   5.0/120,  -5.0/120 },
    {  0.0,       0.0,       0.0,       0.0,       0.0,       1.0/120 }
  };
  return deg == 1 ? &P1[0][0] : (deg == 3 ? &P3[0][0] : &P5[0][0]);
}

// D[n][k][m]: coefficient of t^m in the n-th derivative of the k-th weight polynomial.
struct BsplineDerivTable {
  int deg;
  double D[SPLINE_WMAX][SPLINE_WMAX][SPLINE_WMAX];
  explicit BsplineDerivTable(int deg_) : deg(deg_) {
    const int W = deg + 1;
    const double* P = bspline_poly(deg);
    for (int n = 0; n < SPLINE_WMAX; n++) for (int k = 0; k < SPLINE_WMAX; k++) for (int m = 0; m < SPLINE_WMAX; m++) {
      D[n][k][m] = 0.0;
      int d = m + n;
      if (k >= W || d > deg) continue;
      double fac = 1.0;
      for (int q = m + 1; q <= d; q++) fac *= q;      // d! / m!
      D[n][k][m] = P[k * W + d] * fac;
    }
  }
};

// Weights w[n][k] (n = 0..K derivative order, k = 0..deg) for fractional offset t in [0,1).
// Orders above the degree are zero.
inline void bspline_weights(int deg, double t, int K, double w[SPLINE_WMAX][SPLINE_WMAX]) {
  static const BsplineDerivTable T1(1), T3(3), T5(5);
  const BsplineDerivTable& T = deg == 1 ? T1 : (deg == 3 ? T3 : T5);
  const int W = deg + 1;
  double tp[SPLINE_WMAX]; tp[0] = 1.0; for (int m = 1; m <= deg; m++) tp[m] = tp[m - 1] * t;
  for (int n = 0; n <= K; n++) {
    for (int k = 0; k < W; k++) {
      double acc = 0.0;
      if (n <= deg) { const double* c = T.D[n][k]; for (int m = 0; m <= deg - n; m++) acc += c[m] * tp[m]; }
      w[n][k] = acc;
    }
  }
}

// Whole-sample symmetric (mirror) boundary extension of an index into [0, n-1].
inline int mirror_index(int k, int n) {
  if (n == 1) return 0;
  int period = 2 * (n - 1);
  k = k % period; if (k < 0) k += period;
  return (k < n) ? k : period - k;
}

// ---- 1D interpolation prefilter (Unser 1999; Thevenaz et al. 2000), mirror boundaries ----
inline double init_causal(const double* c, int n, int stride, double z) {
  const double tol = 1e-15;
  int horizon = (int) std::ceil(std::log(tol) / std::log(std::fabs(z)));
  if (horizon < n) {
    double zn = z, sum = c[0];
    for (int k = 1; k < horizon; k++) { sum += zn * c[k * stride]; zn *= z; }
    return sum;
  }
  double zn = z, iz = 1.0 / z, z2n = std::pow(z, (double) (n - 1));
  double sum = c[0] + z2n * c[(n - 1) * stride];
  z2n *= z2n * iz;
  for (int k = 1; k <= n - 2; k++) { sum += (zn + z2n) * c[k * stride]; zn *= z; z2n *= iz; }
  return sum / (1.0 - zn * zn);
}

inline double init_anticausal(const double* c, int n, int stride, double z) {
  return (z / (z * z - 1.0)) * (z * c[(n - 2) * stride] + c[(n - 1) * stride]);
}

// Degree 1 needs no prefilter (the coefficients are the cell values).
inline void prefilter_1d(double* c, int n, int stride, int deg) {
  if (n < 2 || deg == 1) return;
  static const double poles5[2] = { -0.43057534709997379, -0.043096288203264652 };
  static const double poles3[1] = { -0.26794919243112270 };          // sqrt(3) - 2
  const double* poles = deg == 3 ? poles3 : poles5;
  const int np = deg == 3 ? 1 : 2;
  double lambda = 1.0;
  for (int p = 0; p < np; p++) lambda *= (1.0 - poles[p]) * (1.0 - 1.0 / poles[p]);
  for (int k = 0; k < n; k++) c[k * stride] *= lambda;
  for (int p = 0; p < np; p++) {
    double z = poles[p];
    c[0] = init_causal(c, n, stride, z);
    for (int k = 1; k < n; k++) c[k * stride] += z * c[(k - 1) * stride];
    c[(n - 1) * stride] = init_anticausal(c, n, stride, z);
    for (int k = n - 2; k >= 0; k--) c[k * stride] = z * (c[(k + 1) * stride] - c[k * stride]);
  }
}

// ---------------------------------------------------------------------------
// Spline table for a stack of layers on a common grid.
// Values are stored per layer in row-major order, rows from the TOP (ymax) of the raster,
// matching the existing raster_vals layout: idx = layer*(nrow*ncol) + r*ncol + c.
// ---------------------------------------------------------------------------
struct SplineTable {
  int nrow, ncol, nlayer;
  int degree;                 // 1 (bilinear), 3 (bicubic) or 5 (quintic)
  double xmin, ymax, resx, resy;
  bool clamp;                 // TRUE: clamp to the hull of cell centres (barrier);
                              // FALSE: surface and derivatives are 0 outside the hull (covariates)
  std::vector<double> coef;
  // Barrier only: neg_cell[j0 * ncol + i0] = 1 if the spline can be negative somewhere in cell
  // (i0, j0).  The B-spline basis (any degree) is nonnegative and sums to one, so the surface on a
  // cell is a convex combination of its W x W coefficient window: if every coefficient is >= 0 the
  // surface is >= 0 on the whole cell and the barrier potential min(d, 0)^2 / 2 is exactly zero,
  // with all its derivatives.  Empty = no screening.
  std::vector<unsigned char> neg_cell;
  // Barrier only: onset length l of the barrier potential (working units). l = 0: B = 1/2 min(d, 0)^2.
  // l > 0: B is the C4 potential of barrier_potential_jet(), whose force ramps in over the first l inland.
  double onset = 0.0;

  SplineTable() : nrow(0), ncol(0), nlayer(0), degree(5), xmin(0), ymax(0), resx(1), resy(1), clamp(false) {}

  int W() const { return degree + 1; }
  int LO() const { return -(degree - 1) / 2; }
  // Is the partial of order (nx, ny) structurally nonzero inside a cell?
  bool nonzero(int nx, int ny) const { return nx <= degree && ny <= degree; }

  // Fill neg_cell for layer 0 (call after prefilter).
  void build_negative_mask() {
    const int w = W(), lo = LO();
    neg_cell.assign((size_t) nrow * ncol, 0);
    std::vector<double> rmin((size_t) nrow * ncol);
    for (int r = 0; r < nrow; r++)
      for (int c = 0; c < ncol; c++) {
        double m = coef[(size_t) r * ncol + mirror_index(c + lo, ncol)];
        for (int k = 1; k < w; k++) m = std::min(m, coef[(size_t) r * ncol + mirror_index(c + lo + k, ncol)]);
        rmin[(size_t) r * ncol + c] = m;
      }
    for (int r = 0; r < nrow; r++)
      for (int c = 0; c < ncol; c++) {
        double m = rmin[(size_t) mirror_index(r + lo, nrow) * ncol + c];
        for (int k = 1; k < w; k++) m = std::min(m, rmin[(size_t) mirror_index(r + lo + k, nrow) * ncol + c]);
        neg_cell[(size_t) r * ncol + c] = m < 0.0;
      }
  }
  // False only if the surface is certainly >= 0 on the cell containing (x, y) (layer 0, clamped).
  bool maybe_negative(double x, double y) const {
    if (neg_cell.empty()) return true;
    double u = (x - xmin) / resx - 0.5, v = (ymax - y) / resy - 0.5;
    if (!(u == u) || !(v == v)) return true;
    u = std::min(std::max(u, 0.0), (double) (ncol - 1));
    v = std::min(std::max(v, 0.0), (double) (nrow - 1));
    int i0 = std::min((int) std::floor(u), ncol - 1), j0 = std::min((int) std::floor(v), nrow - 1);
    return neg_cell[(size_t) j0 * ncol + i0] != 0;
  }

  void prefilter() {
    for (int l = 0; l < nlayer; l++) {
      double* base = coef.data() + (size_t) l * nrow * ncol;
      for (int r = 0; r < nrow; r++) prefilter_1d(base + (size_t) r * ncol, ncol, 1, degree);
      for (int c = 0; c < ncol; c++) prefilter_1d(base + c, nrow, ncol, degree);
    }
  }

  static int jet_size(int k0, int k1) { int P = 0; for (int k = k0; k <= k1; k++) P += k + 1; return P; }

  // d^{nx+ny} S / dx^nx dy^ny at physical location (x, y) for one layer (reference / tests).
  double eval(int layer, double x, double y, int nx, int ny) const {
    const int k = nx + ny;
    std::vector<double> o(k + 1);
    eval_jet(&layer, 1, x, y, k, k, o.data());
    return o[ny];
  }

  // -------------------------------------------------------------------------
  // Fused "jet" evaluation (the fast path used everywhere in the model).
  // For every layer in `layers` and every derivative order k in [k0, k1], writes the k+1
  // partials d^k S / dx^(k-j) dy^j, j = 0..k, to `out` in the order
  //   out[(l * P) + off(k) + j],   P = sum_{k=k0}^{k1} (k+1),
  // sharing the cell lookup, the basis weights and the W x W window across all of them.
  // -------------------------------------------------------------------------
  void eval_jet(const int* layers, int nl, double x, double y, int k0, int k1, double* out) const {
    switch (degree) {
      case 1: eval_jet_w<2>(layers, nl, x, y, k0, k1, out); break;
      case 3: eval_jet_w<4>(layers, nl, x, y, k0, k1, out); break;
      default: eval_jet_w<6>(layers, nl, x, y, k0, k1, out); break;
    }
  }

  template<int WW>
  void eval_jet_w(const int* layers, int nl, double x, double y, int k0, int k1, double* out) const {
    const int P = jet_size(k0, k1);
    const int deg = WW - 1, lo = -(deg - 1) / 2;
    for (int i = 0; i < nl * P; i++) out[i] = 0.0;
    if (nlayer == 0 || k0 > 2 * deg) return;     // mixed partials of order > 2*deg vanish
    double u = (x - xmin) / resx - 0.5;
    double v = (ymax - y) / resy - 0.5;
    const double umax = ncol - 1, vmax = nrow - 1, tol = 1e-9;
    if (u < 0.0 && u > -tol) u = 0.0;
    if (u > umax && u < umax + tol) u = umax;
    if (v < 0.0 && v > -tol) v = 0.0;
    if (v > vmax && v < vmax + tol) v = vmax;
    bool outu = false, outv = false;
    if (clamp) {
      if (u < 0.0) { outu = true; u = 0.0; } else if (u > umax) { outu = true; u = umax; }
      if (v < 0.0) { outv = true; v = 0.0; } else if (v > vmax) { outv = true; v = vmax; }
    } else if (u < 0.0 || u > umax || v < 0.0 || v > vmax) {
      return;
    }
    const int K = k1 < deg ? k1 : deg;           // highest non-vanishing single-direction order
    int i0 = (int) std::floor(u), j0 = (int) std::floor(v);
    // Degree 1 at the last centre: use the cell to the left/top (t = 1) so the window stays in range
    if (deg == 1) { if (i0 > ncol - 2 && ncol > 1) i0 = ncol - 2; if (j0 > nrow - 2 && nrow > 1) j0 = nrow - 2; }
    double tu = u - i0, tv = v - j0;
    double wu[SPLINE_WMAX][SPLINE_WMAX], wv[SPLINE_WMAX][SPLINE_WMAX];
    bspline_weights(deg, tu, K, wu);
    bspline_weights(deg, tv, K, wv);
    double su = 1.0, sv = 1.0;
    for (int n = 0; n <= K; n++) {
      bool zu = outu && n > 0, zv = outv && n > 0;
      for (int k = 0; k < WW; k++) { wu[n][k] = zu ? 0.0 : wu[n][k] * su; wv[n][k] = zv ? 0.0 : wv[n][k] * sv; }
      su *= 1.0 / resx; sv *= -1.0 / resy;
    }
    int cu[WW], rv[WW];
    bool interior_u = (i0 + lo >= 0) && (i0 + lo + WW - 1 <= ncol - 1);
    bool interior_v = (j0 + lo >= 0) && (j0 + lo + WW - 1 <= nrow - 1);
    for (int k = 0; k < WW; k++) {
      cu[k] = interior_u ? i0 + lo + k : mirror_index(i0 + lo + k, ncol);
      rv[k] = interior_v ? j0 + lo + k : mirror_index(j0 + lo + k, nrow);
    }
    const size_t lsz = (size_t) nrow * ncol;
    for (int l = 0; l < nl; l++) {
      const double* base = coef.data() + (size_t) layers[l] * lsz;
      // row sums rs[n][a] = sum_b wu[n][b] * c(rv[a], cu[b]) for x-orders n = 0..K
      double rs[SPLINE_WMAX][WW];
      for (int a = 0; a < WW; a++) {
        const double* row = base + (size_t) rv[a] * ncol;
        double c[WW];
        if (interior_u) { const double* r = row + cu[0]; for (int b = 0; b < WW; b++) c[b] = r[b]; }
        else { for (int b = 0; b < WW; b++) c[b] = row[cu[b]]; }
        for (int n = 0; n <= K; n++) {
          const double* w = wu[n];
          double s = 0.0;
          for (int b = 0; b < WW; b++) s += w[b] * c[b];
          rs[n][a] = s;
        }
      }
      double* o = out + (size_t) l * P;
      int off = 0;
      for (int k = k0; k <= k1; k++) {
        for (int j = 0; j <= k; j++) {
          int nx = k - j, ny = j;
          if (nx <= K && ny <= K) {
            const double* w = wv[ny]; const double* r = rs[nx];
            double s = 0.0;
            for (int a = 0; a < WW; a++) s += w[a] * r[a];
            o[off + j] = s;
          }
        }
        off += k + 1;
      }
    }
  }
};

// Barrier potential B = 1/2 min(d, 0)^2 (or, if tab.onset > 0, the C4 onset potential below) of the
// spline-interpolated signed distance d (layer 0).
// The barrier force is -lambda * grad(B) = -lambda * min(d, 0) * grad(d).  Writes orders 1..k1
// of B in the jet layout with k0 = 1 (k1 <= 5).  Inside the restricted area each partial follows
// from the jet of d by the Leibniz rule:
//   d^(a,b) B = 1/2 sum_{i<=a, m<=b} C(a,i) C(b,m) d^(i,m) d^(a-i,b-m);
// in the allowed area (d >= 0) all partials are zero.
inline void barrier_potential_jet(const SplineTable& tab, double x, double y, int k1, double* out) {
  const int P = SplineTable::jet_size(1, k1);
  for (int i = 0; i < P; i++) out[i] = 0.0;
  if (!tab.maybe_negative(x, y)) return;          // d > 0 on the whole cell: B and its jet are 0
  const int lay = 0;
  double J[21];                                   // jet_size(0, 5) = 21
  tab.eval_jet(&lay, 1, x, y, 0, k1, J);
  if (!(J[0] < 0.0)) return;
  double D[6][6];
  for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) D[a][b] = 0.0;
  for (int k = 0; k <= k1; k++)
    for (int j = 0; j <= k; j++) D[k - j][j] = J[k * (k + 1) / 2 + j];
  static const double C[6][6] = {{1,0,0,0,0,0},{1,1,0,0,0,0},{1,2,1,0,0,0},
                                 {1,3,3,1,0,0},{1,4,6,4,1,0},{1,5,10,10,5,1}};
  const double l = tab.onset;
  if (l > 0.0) {
    // C4 barrier potential with onset length l. With u = -d the depth into the restricted area and
    // t = u / l, the force magnitude g(u) = dB/du ramps in with the C2 "smootherstep" slope
    // g'(u) = 6t^5 - 15t^4 + 10t^3 and continues with unit slope beyond l:
    //   0 <= u <= l : B = l^2 (t^7/7 - t^6/2 + t^5/2)
    //   u >= l      : B = d^2/2 + l d/2 + l^2/7          (g = u - l/2)
    // B and its first four derivatives are continuous at u = 0 and u = l. B is a polynomial in d on
    // each piece, so its partials follow from those of d by Horner's rule with Leibniz products.
    double c[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    if (J[0] <= -l) { c[0] = l * l / 7.0; c[1] = 0.5 * l; c[2] = 0.5; }
    else { c[5] = -0.5 / (l * l * l); c[6] = -0.5 / (l * l * l * l); c[7] = -1.0 / (7.0 * l * l * l * l * l); }
    int deg = 7; while (deg > 0 && c[deg] == 0.0) deg--;
    double E[6][6], T[6][6];
    for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) { E[a][b] = 0.0; T[a][b] = 0.0; }
    E[0][0] = c[deg];
    for (int n = deg - 1; n >= 0; n--) {          // E <- E * d + c[n]
      for (int a = 0; a <= k1; a++) for (int b = 0; a + b <= k1; b++) {
        double t = 0.0;
        for (int i = 0; i <= a; i++) for (int m = 0; m <= b; m++) t += C[a][i] * C[b][m] * E[i][m] * D[a - i][b - m];
        T[a][b] = t;
      }
      for (int a = 0; a <= k1; a++) for (int b = 0; a + b <= k1; b++) E[a][b] = T[a][b];
      E[0][0] += c[n];
    }
    int o2 = 0;
    for (int k = 1; k <= k1; k++) { for (int j = 0; j <= k; j++) out[o2 + j] = E[k - j][j]; o2 += k + 1; }
    return;
  }
  int off = 0;
  for (int k = 1; k <= k1; k++) {
    for (int j = 0; j <= k; j++) {
      const int a = k - j, b = j;
      double s = 0.0;
      for (int i = 0; i <= a; i++)
        for (int m = 0; m <= b; m++) s += C[a][i] * C[b][m] * D[i][m] * D[a - i][b - m];
      out[off + j] = 0.5 * s;
    }
    off += k + 1;
  }
}

typedef std::shared_ptr<const SplineTable> SplinePtr;

// Scalar evaluation of a single partial (reference implementation; used by tests).
inline double spline_eval(const SplinePtr& tab, int layer, double x, double y, int nx, int ny) {
  return tab->eval(layer, x, y, nx, ny);
}

// Fused jet for plain doubles: returns layers.size() * jet_size(k0, k1) values.
inline std::vector<double> spline_jet(const SplinePtr& tab, const std::vector<int>& layers,
                                      double x, double y, int k0, int k1) {
  std::vector<double> o(layers.size() * SplineTable::jet_size(k0, k1));
  tab->eval_jet(layers.data(), (int) layers.size(), x, y, k0, k1, o.data());
  return o;
}

#ifndef IS_RCPP_BUILD
// ---------------------------------------------------------------------------
// Fused jet operator (used by the model).  Inputs (x, y).  Outputs, for each layer in `layers`
// and EVERY order k in [k0, JET_KMAX], the k+1 partials d^k S / dx^(k-j) dy^j (j = 0..k).
//
// The quintic spline is a degree-5 polynomial on each cell, so all derivatives of order 6
// vanish and this set of outputs is closed under differentiation: the x- and y-derivatives of
// partial (k-j, j) are the outputs (k+1-j, j) and (k-j, j+1) of the SAME operator.  The
// reverse pass is therefore a plain linear combination of the operator's own outputs:
//   * no new raster lookup is ever created by AD, so the gradient, the inner Hessian and the
//     Laplace-gradient tapes all contain exactly ONE lookup per location;
//   * one lookup = one cell search, one set of basis weights (orders 0..JET_KMAX at once) and
//     one 6x6 coefficient window per layer, shared by all outputs.
//
// JET_KMAX is the exact minimum for this model.  The deepest chain is the Laplace gradient:
// the inner objective has a 2nd-order term, its gradient a 3rd, and TMB's implicit
// differentiation of the inner optimum a 4th, so orders 0..4 must be available.  Order 5 is
// never read (verified: setting JET_KMAX = 5 leaves the value, the inner Hessian and the
// outer gradient bit-identical, while JET_KMAX = 3 changes the gradient).  The extra order
// only costs tape (its outputs are pruned as unused) and one 6-term window row per layer, so
// stopping at 4 removes that work with no effect on any result.
// ---------------------------------------------------------------------------
static const int JET_KMAX = 4;

struct RasterJetOp : TMBad::global::DynamicOperator<2, -1> {
  SplinePtr tab;
  std::shared_ptr<const std::vector<int> > layers;
  int k0;
  bool barrier;   // false: jet of the spline layers; true: jet (k0 = 1) of the barrier potential B
  RasterJetOp(SplinePtr tab_, std::shared_ptr<const std::vector<int> > layers_, int k0_, bool barrier_ = false) :
    tab(tab_), layers(layers_), k0(k0_), barrier(barrier_) {}
  int P() const { return SplineTable::jet_size(k0, JET_KMAX); }
  // Structural sparsity of the jet: partial (nx, ny) of a degree-p tensor spline is identically
  // zero if nx > p or ny > p; for the barrier potential B = d^2 / 2 the bound is 2p, and for the
  // onset potential (a degree-7 polynomial in d on each piece) it is 7p.
  bool nonzero_partial(int nx, int ny) const {
    const int p = barrier ? (tab->onset > 0.0 ? 7 : 2) * tab->degree : tab->degree;
    return nx <= p && ny <= p;
  }
  TMBad::Index input_size() const { return 2; }
  TMBad::Index output_size() const { return (TMBad::Index) (layers->size() * P()); }
  static const bool have_input_size_output_size = true;
  static const size_t STACK_MAX = 512;
  typedef TMBad::ad_aug ad;

  void eval(double x, double y, double* out) const {
    if (barrier) barrier_potential_jet(*tab, x, y, JET_KMAX, out);
    else tab->eval_jet(layers->data(), (int) layers->size(), x, y, k0, JET_KMAX, out);
  }
  // Taped (or constant-folded) call
  std::vector<ad> operator()(ad x, ad y) const {
    size_t m = layers->size() * P();
    std::vector<ad> z(m);
    if (!CppAD::Variable(x) && !CppAD::Variable(y)) {
      std::vector<double> o(m);
      eval(asDouble(x), asDouble(y), o.data());
      for (size_t i = 0; i < m; i++) z[i] = ad(o[i]);
      return z;
    }
    std::vector<ad> xy(2); xy[0] = x; xy[1] = y;
    return TMBad::global::Complete<RasterJetOp>(*this)(xy);
  }

  void forward(TMBad::ForwardArgs<double> &args) {
    size_t m = output_size();
    double buf[STACK_MAX]; std::vector<double> big;
    double* o = buf;
    if (m > STACK_MAX) { big.resize(m); o = big.data(); }
    eval(args.x(0), args.x(1), o);
    for (size_t i = 0; i < m; i++) args.y(i) = o[i];
  }
  void forward(TMBad::ForwardArgs<ad> &args) {
    std::vector<ad> z = (*this)(args.x(0), args.x(1));
    for (size_t i = 0; i < z.size(); i++) args.y(i) = z[i];
  }
  // d out(k-j, j)/dx = out(k+1-j, j),  d out(k-j, j)/dy = out(k-j, j+1);  order-5 outputs have
  // zero derivative.  Uses only the operator's own outputs (no new lookup).
  template<class T>
  void reverse_impl(TMBad::ReverseArgs<T> &args) {
    const int Pc = P();
    T dx0 = T(0.0), dx1 = T(0.0);
    for (size_t l = 0; l < layers->size(); l++) {
      const size_t base = l * Pc;
      int off = 0;
      for (int k = k0; k < JET_KMAX; k++) {
        const int offn = off + k + 1;
        for (int j = 0; j <= k; j++) {
          // Skip partials that are identically zero for this degree (their derivatives too):
          // a degree-p spline has no x- or y-derivative above order p, and the barrier
          // potential B = d^2/2 none above 2p.  For p = 5 nothing is skipped at these orders;
          // for p = 1 and 3 this removes most of the taped products.
          if (!nonzero_partial(k + 1 - j, j) && !nonzero_partial(k - j, j + 1)) continue;
          T dy = args.dy(base + off + j);
          dx0 += args.y(base + offn + j) * dy;
          dx1 += args.y(base + offn + j + 1) * dy;
        }
        off = offn;
      }
    }
    args.dx(0) += dx0;
    args.dx(1) += dx1;
  }
  void reverse(TMBad::ReverseArgs<double> &args) { reverse_impl(args); }
  void reverse(TMBad::ReverseArgs<ad> &args) { reverse_impl(args); }
  void forward(TMBad::ForwardArgs<TMBad::Writer> &args) { Rf_error("RasterJetOp: source code writer not implemented"); }
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) { Rf_error("RasterJetOp: source code writer not implemented"); }
  const char* op_name() { return barrier ? "BJET" : "RJET"; }
};

inline std::vector<TMBad::ad_aug> spline_jet(const SplinePtr& tab, const std::vector<int>& layers,
                                             TMBad::ad_aug x, TMBad::ad_aug y, int k0, int k1) {
  // One operator computing orders [k0, 5]; return the requested orders [k0, k1] per layer,
  // in the same layout as the double version.
  RasterJetOp op(tab, std::make_shared<const std::vector<int> >(layers), k0);
  std::vector<TMBad::ad_aug> z = op(x, y);
  const int Pf = op.P(), Pr = SplineTable::jet_size(k0, k1);
  std::vector<TMBad::ad_aug> r(layers.size() * Pr);
  for (size_t l = 0; l < layers.size(); l++)
    for (int i = 0; i < Pr; i++) r[l * Pr + i] = z[l * Pf + i];
  return r;
}

// grad B for the barrier potential (B = 1/2 min(d, 0)^2, or the onset potential), as one fused
// operator: the switch between the allowed and restricted areas happens inside the operator, so
// no conditional expressions are taped and the reverse pass is the same closed linear
// combination of the operator's own outputs as for the covariate jet.
inline std::vector<TMBad::ad_aug> barrier_grad(const SplinePtr& tab, TMBad::ad_aug x, TMBad::ad_aug y) {
  RasterJetOp op(tab, std::make_shared<const std::vector<int> >(1, 0), 1, true);
  std::vector<TMBad::ad_aug> z = op(x, y);
  std::vector<TMBad::ad_aug> r(2); r[0] = z[0]; r[1] = z[1];
  return r;
}
#endif

inline std::vector<double> barrier_grad(const SplinePtr& tab, double x, double y) {
  double o[2];
  barrier_potential_jet(*tab, x, y, 1, o);
  return std::vector<double>(o, o + 2);
}

#ifndef IS_RCPP_BUILD
// Negative part min(d, 0), kept on the tape (no frozen branch).
template<class Type>
Type negative_part(Type d) { return CppAD::CondExpLe(d, Type(0.0), d, Type(0.0)); }
inline double negative_part(double d) { return d < 0.0 ? d : 0.0; }
#else
inline double negative_part(double d) { return d < 0.0 ? d : 0.0; }
#endif

} // namespace langevin_raster

// ---------------------------------------------------------------------------
// Builders (run once per template / simulator evaluation; inputs are DATA)
// ---------------------------------------------------------------------------
TEMPLATE_HEADER
langevin_raster::SplinePtr build_cov_spline(RASTER_TYPE &raster_vals,
                                            const VECTOR &raster_resolution,
                                            const VECTOR &raster_extent,
                                            int degree = 5) {
  std::shared_ptr<langevin_raster::SplineTable> tab = std::make_shared<langevin_raster::SplineTable>();
  tab->degree = degree;
  tab->resx = AS_DOUBLE(VEC_ELT(raster_resolution, 0));
  tab->resy = AS_DOUBLE(VEC_ELT(raster_resolution, 1));
  tab->xmin = AS_DOUBLE(VEC_ELT(raster_extent, 0));
  tab->ymax = AS_DOUBLE(VEC_ELT(raster_extent, 3));
  tab->ncol = static_cast<int>(std::round((AS_DOUBLE(VEC_ELT(raster_extent, 1)) - tab->xmin) / tab->resx));
  tab->nrow = static_cast<int>(std::round((tab->ymax - AS_DOUBLE(VEC_ELT(raster_extent, 2))) / tab->resy));
  tab->clamp = false;
  int n = (int) raster_vals.size();
  tab->nlayer = n / (tab->nrow * tab->ncol);
  tab->coef.resize(n);
  for (int i = 0; i < n; i++) tab->coef[i] = AS_DOUBLE(GET_VAL(raster_vals, i));
  tab->prefilter();
  return tab;
}

TEMPLATE_HEADER
langevin_raster::SplinePtr build_barrier_spline(const MATRIX &barrier_dist,
                                                const VECTOR &raster_resolution,
                                                const VECTOR &raster_extent,
                                                int degree = 5, double onset = 0.0) {
  std::shared_ptr<langevin_raster::SplineTable> tab = std::make_shared<langevin_raster::SplineTable>();
  tab->degree = degree;
  tab->resx = AS_DOUBLE(VEC_ELT(raster_resolution, 0));
  tab->resy = AS_DOUBLE(VEC_ELT(raster_resolution, 1));
  tab->xmin = AS_DOUBLE(VEC_ELT(raster_extent, 0));
  tab->ymax = AS_DOUBLE(VEC_ELT(raster_extent, 3));
  tab->nrow = MAT_NROW(barrier_dist);
  tab->ncol = MAT_NCOL(barrier_dist);
  tab->nlayer = 1;
  tab->clamp = true;
  tab->coef.resize((size_t) tab->nrow * tab->ncol);
  for (int r = 0; r < tab->nrow; r++)
    for (int c = 0; c < tab->ncol; c++)
      tab->coef[(size_t) r * tab->ncol + c] = AS_DOUBLE(barrier_dist(r, c));
  tab->prefilter();
  tab->build_negative_mask();
  tab->onset = onset;
  return tab;
}

// ---------------------------------------------------------------------------
// Gradient of every covariate surface at (x, y, z): 2 x n_covs matrix.
// Dynamic covariates are interpolated linearly in time (z is DATA, so branching on it is safe).
// ---------------------------------------------------------------------------
TEMPLATE_HEADER
MATRIX extract_raster_values(TYPE x, TYPE y, TYPE z,
                             const langevin_raster::SplinePtr &cov_spline,
                             const VECTOR &all_z_values,
                             const INT_VECTOR &n_zvals_cov,
                             const INT_VECTOR &cov_offset,
                             int n_covs) {

  MATRIX grad_values(2, n_covs);
  SET_ZERO(grad_values);
  double zd = AS_DOUBLE(z);
  std::vector<int> lay; std::vector<double> wz; std::vector<bool> two;
  lay.reserve(2 * n_covs);

  for (int i = 0; i < n_covs; i++) {
    int layers = VEC_ELT(n_zvals_cov, i);
    int offset = VEC_ELT(cov_offset, i);
    int z_idx1 = 0, z_idx2 = 0;
    double z_weight = 0.0;

    if (layers > 1) {
      double zfirst = AS_DOUBLE(VEC_ELT(all_z_values, offset));
      double zlast = AS_DOUBLE(VEC_ELT(all_z_values, offset + layers - 1));
      if (zd <= zfirst) {
        z_idx1 = 0; z_idx2 = 0;
      } else if (zd >= zlast) {
        z_idx1 = layers - 1; z_idx2 = layers - 1;
      } else {
        for (int k = 0; k < layers - 1; k++) {
          double t1 = AS_DOUBLE(VEC_ELT(all_z_values, offset + k));
          double t2 = AS_DOUBLE(VEC_ELT(all_z_values, offset + k + 1));
          if (zd >= t1 && zd <= t2) {
            z_idx1 = k; z_idx2 = k + 1;
            z_weight = (zd - t1) / (t2 - t1);
            break;
          }
        }
      }
    }

    lay.push_back(offset + z_idx1);
    if (layers > 1 && z_idx1 != z_idx2) lay.push_back(offset + z_idx2);
    wz.push_back(z_weight);
    two.push_back(layers > 1 && z_idx1 != z_idx2);
  }

  if (lay.empty()) return grad_values;
  // One fused lookup for all covariate layers: outputs (dS/dx, dS/dy) per layer.
  std::vector<TYPE> J = langevin_raster::spline_jet(cov_spline, lay, x, y, 1, 1);
  int p = 0;
  for (int i = 0; i < n_covs; i++) {
    if (two[i]) {
      double w = wz[i];
      grad_values(0, i) = J[2 * p] * TYPE(1.0 - w) + J[2 * (p + 1)] * TYPE(w);
      grad_values(1, i) = J[2 * p + 1] * TYPE(1.0 - w) + J[2 * (p + 1) + 1] * TYPE(w);
      p += 2;
    } else {
      grad_values(0, i) = J[2 * p];
      grad_values(1, i) = J[2 * p + 1];
      p += 1;
    }
  }
  return grad_values;
}

// ---------------------------------------------------------------------------
// Smoothing radius for calculate_smooth_gradient().
//
// The transition density holds the force h at its value at the start of the interval, x_0.
// In the exact process the force acts along the path x_s, s in [0, D], and enters the mean
// of the next position as  s2 * int_0^D w(s) h(x_s) ds, with kernel
//   underdamped: w(s) = (1 - exp(-gamma (D - s))) / gamma,     overdamped: w(s) = 1/2.
// To second order, and averaging over the stationary velocity (the ring is isotropic and
// cannot use the latent velocity), the w-weighted average force is
//   h_bar = h(x_0) + (1/2) m Laplacian(h)(x_0),   m = int w(s) M(s) ds / int w(s) ds,
// where M(s) is the per-coordinate mean squared displacement at lag s.  A ring of n >= 3
// equally spaced points at radius r, with centre weight c and the remaining weight shared
// equally, has per-coordinate second moment (1 - c) r^2 / 2, zero cross moment and zero odd
// moments (n = 4 and n = 8).  So the ring average reproduces h_bar to second order when
//   r* = sqrt(2 m / (1 - c)).
//
// Underdamped (dv = -gamma v dt + s2 h dt + sigma sqrt(2 gamma) dW; stationary Var(v) = s2 per
// coordinate, as in the transition and initial-velocity densities):
//   M(s) = (2 s2 / gamma^2) (gamma s - 1 + exp(-gamma s)), and with a = gamma D,
//   m    = (2 s2 / gamma^2) I(a) / W(a),
//   I(a) = a^2/2 - 2a + 3 - (3 + a) exp(-a) = sum_{k>=4} (-1)^k (k - 3) a^k / k!,
//   W(a) = a - 1 + exp(-a)                  = sum_{k>=2} (-1)^k a^k / k!.
//   Limits: a -> 0 (ballistic)  r* -> sigma D / sqrt(3 (1 - c));
//           a -> inf (diffusive) r* -> sigma sqrt(2 D / (gamma (1 - c))).
// Overdamped (dx = s2/2 h dt + sigma dW): M(s) = s2 s, uniform w, m = s2 D / 2, so
//   r* = sigma sqrt(D / (1 - c)).
// For c = 1/2 this is sqrt(2 sigma^2 D), the RMS step length used by Blackwell & Matthiopoulos
// (2024, Eqs. 29-30).
//
// Radius used (both models): r = zetaScale * r*, so zetaScale = 1 is the moment-matched radius.
// Following Blackwell & Matthiopoulos (2024), fitLangevin evaluates r* at FIXED values of sigma and
// gamma derived from the data before optimisation (error-corrected variogram estimates), passed as
// DATA, so the smoothing neighbourhood does not change with the parameters being estimated.
// (The previous rule zetaScale * sqrt(2 pi) * sigma * sqrt(D) ignored gamma and scaled as
// sqrt(D) although sigma is a velocity scale in this model; at gamma D << 1 it was 10 to 30
// times larger than r*.)
// ---------------------------------------------------------------------------

// I(a) / W(a) for the underdamped radius, accurate and AD-safe for all a > 0.  For small a the
// closed forms suffer catastrophic cancellation (I ~ a^4/24), so a truncated series is used
// below a = 0.25, where both branches are accurate to < 5e-13 relative.  The branch condition
// depends on gamma (a parameter), so in the TMB build it is taped with CondExp, not a C++ if.
TEMPLATE_HEADER
TYPE ud_smooth_IW_ratio(TYPE a) {
  // Series: I/W = a^2 P(a) / Q(a),
  //   P(a) = sum_{k=4}^{12} (-1)^k (k-3)/k! a^(k-4),   Q(a) = sum_{k=2}^{10} (-1)^k / k! a^(k-2)
  TYPE P = TYPE(0.0), Q = TYPE(0.0), apow = TYPE(1.0);
  double fact = 24.0;                                   // 4!
  for (int k = 4; k <= 12; k++) {
    double sgn = (k % 2 == 0) ? 1.0 : -1.0;
    P += TYPE(sgn * (k - 3) / fact) * apow;
    apow *= a;
    fact *= (k + 1);
  }
  apow = TYPE(1.0); fact = 2.0;                         // 2!
  for (int k = 2; k <= 10; k++) {
    double sgn = (k % 2 == 0) ? 1.0 : -1.0;
    Q += TYPE(sgn / fact) * apow;
    apow *= a;
    fact *= (k + 1);
  }
  TYPE series = a * a * P / Q;

  const double a_switch = 0.25;
  // Closed form, evaluated at max(a, a_switch) so the unused branch (and its derivative) stays finite
#ifdef IS_RCPP_BUILD
  TYPE ac = (a > a_switch) ? a : TYPE(a_switch);
#else
  TYPE ac = CppAD::CondExpGt(a, TYPE(a_switch), a, TYPE(a_switch));
#endif
  TYPE e = exp(-ac);
  TYPE I = ac * ac / TYPE(2.0) - TYPE(2.0) * ac + TYPE(3.0) - (TYPE(3.0) + ac) * e;
  TYPE W = ac - TYPE(1.0) + e;
  TYPE closed = I / W;

#ifdef IS_RCPP_BUILD
  return (a < a_switch) ? series : closed;
#else
  return CppAD::CondExpLt(a, TYPE(a_switch), series, closed);
#endif
}

// Ring radius used by calculate_smooth_gradient(); see the derivation above.
// c = centre weight (curweight).  model == 1: underdamped.
TEMPLATE_HEADER
TYPE smooth_gradient_radius(TYPE sigma, TYPE gamma, TYPE dt_step, TYPE c, TYPE zetaScale, int model) {
  if (model == 1) {
    TYPE a = gamma * dt_step;
    TYPE m = TYPE(2.0) * sigma * sigma / (gamma * gamma) * ud_smooth_IW_ratio(a);
    return zetaScale * sqrt(TYPE(2.0) * m / (TYPE(1.0) - c));
  } else {
    // overdamped: r* = sigma sqrt(D / (1 - c)); equals the Blackwell & Matthiopoulos (2024)
    // radius (RMS step length sqrt(2 sigma^2 D)) for c = 1/2
    return zetaScale * sigma * sqrt(dt_step / (TYPE(1.0) - c));
  }
}

TEMPLATE_HEADER
MATRIX calculate_smooth_gradient(TYPE x, TYPE y, TYPE z,
                                 const langevin_raster::SplinePtr &cov_spline,
                                 const VECTOR &all_z_values,
                                 const INT_VECTOR &n_zvals_cov,
                                 const INT_VECTOR &cov_offset,
                                 int n_covs,
                                 const VECTOR &weights,
                                 TYPE sigma,
                                 TYPE gamma,
                                 TYPE dt_step,
                                 TYPE zetaScale,
                                 int model) {

  // Ring radius: moment-matched to the force kernel of the transition density (see
  // smooth_gradient_radius).  weights(0) is the centre weight (curweight).  sigma and gamma here are
  // the fixed smoothing values supplied as DATA (not the parameters being estimated).
  TYPE neighborhood_size = smooth_gradient_radius(sigma, gamma, dt_step,
                                                  TYPE(VEC_ELT(weights, 0)), zetaScale, model);

  MATRIX smooth_grads = VEC_ELT(weights, 0) * extract_raster_values(x, y, z, cov_spline,
                                all_z_values, n_zvals_cov, cov_offset, n_covs);

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
      z, cov_spline, all_z_values, n_zvals_cov, cov_offset, n_covs);
  }

  return smooth_grads;
}

// Bilinear value lookup (double only); used by extract_sdf_rcpp for barrier checks in R.
TEMPLATE_HEADER
TYPE get_bilinear_val(TYPE x, TYPE y, const MATRIX &grid, const VECTOR &ext, const VECTOR &res) {

  int n_cols = MAT_NCOL(grid);
  int n_rows = MAT_NROW(grid);

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

// ---------------------------------------------------------------------------
// Barrier force: h -= lambda * grad(B), with d the spline-interpolated signed distance (positive =
// allowed) and B = 1/2 min(d, 0)^2 (gradient min(d, 0) * grad(d)) or the onset potential, evaluated by
// one fused operator (barrier_grad), so the active/inactive state is re-decided at every
// evaluation (never frozen at taping) and AD of the force uses only the operator's own outputs.
// ---------------------------------------------------------------------------
TEMPLATE_HEADER
void apply_barrier_penalty(TYPE x, TYPE y,
                           const langevin_raster::SplinePtr &barrier_spline,
                           TYPE barrier_penalty,
                           TYPE &h_x, TYPE &h_y) {
  if (AS_DOUBLE(barrier_penalty) > 0.0 && barrier_spline) {
    std::vector<TYPE> g = langevin_raster::barrier_grad(barrier_spline, x, y);
    h_x -= barrier_penalty * g[0];
    h_y -= barrier_penalty * g[1];
  }
}

#endif // RASTER_HELPERS_HPP
