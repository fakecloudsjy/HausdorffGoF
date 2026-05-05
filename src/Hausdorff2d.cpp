// Hausdorff2d.cpp
//
// H_stat_2s_2d_cpp: Hausdorff distance H(F_m^c, G_n^c) under the Chebyshev
// metric, computed via a single hsearch() call (Appendix A, paper).
//
// The algorithm is one-directional by construction:
//   projection_x = V_loc(F_m | G_n)  [locally-farthest vertices of F_m]
//   projection_y = all vertices of G_n
//
// Symmetry H(x,y) = H(y,x) follows from correct V_loc, not from two calls:
//   - Concave vertices (where F_m > G_n) capture F-above-G contributions.
//   - Convex vertices  (where G_n > F_m) capture G-above-F contributions.
// Both are in V_loc of F_m, so a single pass suffices.
//
// The key fix relative to the original code: the convex vertex condition
// is simply  G_n(ax,ay) > F_m(ax-,ay-)  (i.e. zy > z4).
// The extra guard |z4-z3| > eps that appeared in the R code is REMOVED:
// for generic continuous data z3 = z4 at every interior omnidirectional jump
// (because z3-z4 counts observations with x1 < ax AND x2 = ay exactly, which
// is zero), so the guard silently suppresses ALL interior convex vertices.
//
// Column layouts expected by hsearch() (0-indexed):
//   projection_x : [lambda0, lambda1, loc1, loc2, z]        (5 cols)
//   projection_y : [lambda0, lambda1, loc1, loc2, z, type]  (6 cols)

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <memory>
#include "Hausdorff2d.h"
#include "Hausdorffsearch.h"
#include "fastCDF.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

static MatrixXd rows_to_matrix(const std::vector<RowVectorXd> &rows, int ncols)
{
    if (rows.empty()) return MatrixXd(0, ncols);
    MatrixXd m(rows.size(), ncols);
    for (int i = 0; i < (int)rows.size(); ++i) m.row(i) = rows[i];
    return m;
}

double H_stat_2s_2d_cpp(const MatrixXd &x, const MatrixXd &y, double tol)
{
    const int mx = x.rows(), my = y.rows();

    // Sorted marginals
    VectorXd x1 = x.col(0), x2 = x.col(1);
    std::sort(x1.data(), x1.data() + mx);
    std::sort(x2.data(), x2.data() + mx);
    VectorXd y1 = y.col(0), y2 = y.col(1);
    std::sort(y1.data(), y1.data() + my);
    std::sort(y2.data(), y2.data() + my);

    const double M = std::round(std::max(x.maxCoeff(), y.maxCoeff())) + 2.5;

    // fastCDF expects (dimension x nbSim)
    ArrayXXd x_t(2, mx), y_t(2, my);
    for (int j = 0; j < mx; ++j) { x_t(0,j)=x(j,0); x_t(1,j)=x(j,1); }
    for (int j = 0; j < my; ++j) { y_t(0,j)=y(j,0); y_t(1,j)=y(j,1); }
    ArrayXd ones_x = ArrayXd::Ones(mx);
    ArrayXd ones_y = ArrayXd::Ones(my);

    // ------------------------------------------------------------------
    // Grid 1: F_m on (x1-tol, x1, M) x (x2-tol, x2, M)  size ng = 2*mx+1
    //   z_x[i0 + i1*ng] = F_m(g1[i0], g2[i1])
    //   Accessors (0-indexed xi,yi in 0..mx-1):
    //     F_m(x1[xi],   x2[yi]  ) = z_x_at(2*xi+1, 2*yi+1)   z1
    //     F_m(x1[xi],   x2[yi]-t) = z_x_at(2*xi+1, 2*yi  )   z2
    //     F_m(x1[xi]-t, x2[yi]  ) = z_x_at(2*xi,   2*yi+1)   z3
    //     F_m(x1[xi]-t, x2[yi]-t) = z_x_at(2*xi,   2*yi  )   z4
    //   Boundary:
    //     F_m(x1[xi], M)    = z_x_at(2*xi+1, 2*mx)
    //     F_m(x1[xi]-t, M)  = z_x_at(2*xi,   2*mx)
    //     F_m(M, x2[yi])    = z_x_at(2*mx,   2*yi+1)
    //     F_m(M, x2[yi]-t)  = z_x_at(2*mx,   2*yi  )
    // ------------------------------------------------------------------
    const int ng = 2*mx + 1;
    VectorXd g1(ng), g2(ng);
    for (int i=0;i<mx;++i){g1(2*i)=x1(i)-tol; g1(2*i+1)=x1(i);}
    for (int i=0;i<mx;++i){g2(2*i)=x2(i)-tol; g2(2*i+1)=x2(i);}
    g1(2*mx)=M; g2(2*mx)=M;
    auto pg1 = std::make_shared<ArrayXd>(g1.array());
    auto pg2 = std::make_shared<ArrayXd>(g2.array());
    ArrayXd z_x = fastCDF(x_t, {pg1, pg2}, ones_x);
    auto z_x_at = [&](int i0,int i1){ return z_x(i0 + i1*ng); };

    // ------------------------------------------------------------------
    // Grid 2: G_n on (x1, M) x (x2, M)  size ngc = mx+1
    //   z_yx_at(xi, yi) = G_n(x1[xi], x2[yi])
    //   z_yx_at(xi, mx) = G_n(x1[xi], M)
    //   z_yx_at(mx, yi) = G_n(M, x2[yi])
    // ------------------------------------------------------------------
    const int ngc = mx+1;
    VectorXd gc1(ngc), gc2(ngc);
    for (int i=0;i<mx;++i){gc1(i)=x1(i); gc2(i)=x2(i);}
    gc1(mx)=M; gc2(mx)=M;
    auto pgc1 = std::make_shared<ArrayXd>(gc1.array());
    auto pgc2 = std::make_shared<ArrayXd>(gc2.array());
    ArrayXd z_yx = fastCDF(y_t, {pgc1, pgc2}, ones_y);
    auto z_yx_at = [&](int i0,int i1){ return z_yx(i0 + i1*ngc); };

    // ==================================================================
    // Step 4: Build projection_x = V_loc(F_m | G_n)
    //
    // At each omnidirectional jump location (ax, ay) of F_m:
    //   Concave vertex (ax, ay, z1): in V_loc iff z1 > G_n(ax,ay)
    //   Convex  vertex (ax, ay, z4): in V_loc iff G_n(ax,ay) > z4
    //
    // NOTE: NO |z4-z3| guard.  For generic continuous data, z3=z4 at every
    // interior jump (z3-z4 = #{x1<ax, x2=ay exactly}/m = 0), so adding the
    // guard would suppress ALL interior convex vertices -- breaking symmetry.
    // ==================================================================
    std::vector<RowVectorXd> px_rows;

    // Interior: (x1[xi], x2[yi]) for xi,yi = 0..mx-1
    for (int xi=0; xi<mx; ++xi) {
        const double ax = x1(xi);
        for (int yi=0; yi<mx; ++yi) {
            const double ay = x2(yi);
            const double z1 = z_x_at(2*xi+1, 2*yi+1);
            const double z2 = z_x_at(2*xi+1, 2*yi  );
            const double z3 = z_x_at(2*xi,   2*yi+1);
            const double z4 = z_x_at(2*xi,   2*yi  );
            if (z1==z2 || z1==z3) continue;  // not omnidirectional
            const double zy = z_yx_at(xi, yi);
            if (z1 > zy) {                    // concave vertex in V_loc
                RowVectorXd r(5); r << ax+z1, ay+z1, ax, ay, z1;
                px_rows.push_back(r);
            }
            if (zy > z4 + 1e-14) {            // convex vertex in V_loc
                RowVectorXd r(5); r << ax+z4, ay+z4, ax, ay, z4;
                px_rows.push_back(r);
            }
        }
        // Boundary at (x1[xi], M)
        {
            const double z1b = z_x_at(2*xi+1, 2*mx);
            const double z4b = z_x_at(2*xi,   2*mx);
            if (z1b != z4b) {
                const double zyb = z_yx_at(xi, mx);
                if (z1b > zyb) {
                    RowVectorXd r(5); r << x1(xi)+z1b, M+z1b, x1(xi), M, z1b;
                    px_rows.push_back(r);
                }
                if (zyb > z4b + 1e-14) {
                    RowVectorXd r(5); r << x1(xi)+z4b, M+z4b, x1(xi), M, z4b;
                    px_rows.push_back(r);
                }
            }
        }
    }
    // Boundary at (M, x2[yi])
    for (int yi=0; yi<mx; ++yi) {
        const double z1b = z_x_at(2*mx, 2*yi+1);
        const double z4b = z_x_at(2*mx, 2*yi  );
        if (z1b != z4b) {
            const double zyb = z_yx_at(mx, yi);
            if (z1b > zyb) {
                RowVectorXd r(5); r << M+z1b, x2(yi)+z1b, M, x2(yi), z1b;
                px_rows.push_back(r);
            }
            if (zyb > z4b + 1e-14) {
                RowVectorXd r(5); r << M+z4b, x2(yi)+z4b, M, x2(yi), z4b;
                px_rows.push_back(r);
            }
        }
    }
    MatrixXd projection_x = rows_to_matrix(px_rows, 5);

    // ==================================================================
    // Step 5: Build projection_y = all vertices of G_n (types 1,2,3,4)
    //   Grid: (y1-tol, y1) x (y2-tol, y2)  size ngy = 2*my
    // ==================================================================
    const int ngy = 2*my;
    VectorXd gy1(ngy), gy2(ngy);
    for (int i=0;i<my;++i){gy1(2*i)=y1(i)-tol; gy1(2*i+1)=y1(i);}
    for (int i=0;i<my;++i){gy2(2*i)=y2(i)-tol; gy2(2*i+1)=y2(i);}
    auto pgy1 = std::make_shared<ArrayXd>(gy1.array());
    auto pgy2 = std::make_shared<ArrayXd>(gy2.array());
    ArrayXd z_y = fastCDF(y_t, {pgy1, pgy2}, ones_y);
    auto z_y_at = [&](int i0,int i1){ return z_y(i0 + i1*ngy); };

    std::vector<RowVectorXd> py_rows;
    for (int xi=0; xi<my; ++xi) {
        const double b1 = y1(xi);
        for (int yi=0; yi<my; ++yi) {
            const double b2 = y2(yi);
            const double z1 = z_y_at(2*xi+1, 2*yi+1);
            const double z2 = z_y_at(2*xi+1, 2*yi  );
            const double z3 = z_y_at(2*xi,   2*yi+1);
            const double z4 = z_y_at(2*xi,   2*yi  );
            if (z1==z2 || z1==z3) continue;
            { RowVectorXd r(6); r<<b1+z2,b2+z2,b1,b2,z2,2.0; py_rows.push_back(r); }
            if (z4!=z2) {
                { RowVectorXd r(6); r<<b1+z4,b2+z4,b1,b2,z4,4.0; py_rows.push_back(r); }
                if (z2!=z3) {
                    { RowVectorXd r(6); r<<b1+z3,b2+z3,b1,b2,z3,3.0; py_rows.push_back(r); }
                }
            }
            { RowVectorXd r(6); r<<b1+z1,b2+z1,b1,b2,z1,1.0; py_rows.push_back(r); }
        }
    }
    MatrixXd projection_y = rows_to_matrix(py_rows, 6);

    return hsearch(projection_x, projection_y);
}

// [[Rcpp::export]]
double H_stat_2s_2d(NumericMatrix x_r, NumericMatrix y_r, double tol = 1e-6)
{
    return H_stat_2s_2d_cpp(as<MatrixXd>(x_r), as<MatrixXd>(y_r), tol);
}



/*
// Hausdorff2d.cpp
//
// Implements H_stat_2s_2d_cpp: given two independent bivariate samples
// x (m x 2) and y (n x 2), returns the Hausdorff distance H(F^c_m, G^c_n)
// between their bivariate empirical CDF surfaces under the Chebyshev metric.
//
// The projection matrices fed to hsearch() are built entirely from direct
// bivariate ECDF counting, avoiding the fastCDF grid-ordering dependency
// present in the original R version.
//
// Column layouts (0-indexed) expected by hsearch():
//   projection_x  [proj1, proj2, loc1, loc2, z]         (5 cols)
//   projection_y  [proj1, proj2, loc1, loc2, z, type]   (6 cols)
// where proj_i = loc_i + z for both matrices.
//
// Algorithm reference: Appendix A of Dimitrova, Jia & Kaishev (2025).

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include "Hausdorff2d.h"
#include "Hausdorffsearch.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

// -----------------------------------------------------------------------------
// Private helpers - internal to this translation unit only
// -----------------------------------------------------------------------------

/// \brief Bivariate ECDF: fraction of rows of s with s[k,0]<=a AND s[k,1]<=b
static inline double becdf(const MatrixXd &s, double a, double b)
{
    int n = s.rows(), cnt = 0;
    for (int k = 0; k < n; ++k)
        cnt += (s(k, 0) <= a && s(k, 1) <= b);
    return static_cast<double>(cnt) / n;
}

/// \brief Append one row to a dynamically growing matrix
static inline void push_row(MatrixXd &mat, const RowVectorXd &row)
{
    mat.conservativeResize(mat.rows() + 1, NoChange);
    mat.bottomRows(1) = row;
}

// -----------------------------------------------------------------------------
// Core - Eigen implementation
// -----------------------------------------------------------------------------

double H_stat_2s_2d_cpp(const MatrixXd &x, const MatrixXd &y, double tol)
{
    const int mx = x.rows(), my = y.rows();

    // Sorted marginal vectors of x
    VectorXd x1 = x.col(0), x2 = x.col(1);
    std::sort(x1.data(), x1.data() + mx);
    std::sort(x2.data(), x2.data() + mx);

    // Sorted marginal vectors of y
    VectorXd y1 = y.col(0), y2 = y.col(1);
    std::sort(y1.data(), y1.data() + my);
    std::sort(y2.data(), y2.data() + my);

    // Truncation boundary (same formula as R: round(max) + 2.5)
    const double M = std::round(std::max(x.maxCoeff(), y.maxCoeff())) + 2.5;

    // -- projection_x : 5 columns ---------------------------------------------
    // Built from the omnidirectional jumps of F_m (upper curve).
    // For each jump location (ax, ay) of F_m up to two rows may be emitted:
    //   concave vertex: (ax+z1, ay+z1, ax, ay, z1)   when z1 > zy
    //   convex  vertex: (ax+z4, ay+z4, ax, ay, z4)   when zy > z4+eps
    //                                                  AND |z4-z3| > eps
    // where z1=Fm(ax,ay), z2=Fm(ax,ay-tol), z3=Fm(ax-tol,ay),
    //       z4=Fm(ax-tol,ay-tol), zy=Gn(ax,ay).
    // -------------------------------------------------------------------------
    MatrixXd projection_x(0, 5);

    for (int xi = 0; xi < mx; ++xi)
    {
        const double ax = x1[xi];

        for (int yi = 0; yi < mx; ++yi)
        {
            const double ay = x2[yi];

            const double z1 = becdf(x, ax,       ay      );
            const double z2 = becdf(x, ax,       ay - tol);
            const double z3 = becdf(x, ax - tol, ay      );
            const double z4 = becdf(x, ax - tol, ay - tol);

            // Omnidirectional jump: F_m must jump in both directions
            if (z1 == z2 || z1 == z3) continue;

            const double zy = becdf(y, ax, ay);

            if (z1 > zy) {
                RowVectorXd row(5);
                row << ax + z1, ay + z1, ax, ay, z1;
                push_row(projection_x, row);
            }
            if (zy > z4 + 1e-14 && std::abs(z4 - z3) > 1e-14) {
                RowVectorXd row(5);
                row << ax + z4, ay + z4, ax, ay, z4;
                push_row(projection_x, row);
            }
        }

        // -- Boundary at (x1[xi], M) ------------------------------------------
        // The truncation boundary M introduces additional jump vertices along
        // each boundary face (Lemma 2.3, Appendix A).  Along the face x2 = M
        // variation is only in the first coordinate, so z2 and z3 are
        // synthesised as the midpoint (no independent saddle value available).
        {
            const double z1b = becdf(x, ax,       M);
            const double z4b = becdf(x, ax - tol, M);
            if (z1b != z4b) {
                const double zyb = becdf(y, ax, M);
                if (z1b > zyb) {
                    RowVectorXd row(5);
                    row << ax + z1b, M + z1b, ax, M, z1b;
                    push_row(projection_x, row);
                }
                // |z4b - z_mid| > 0 whenever z1b != z4b (condition above)
                if (zyb > z4b + 1e-14) {
                    RowVectorXd row(5);
                    row << ax + z4b, M + z4b, ax, M, z4b;
                    push_row(projection_x, row);
                }
            }
        }
    }

    // -- Boundary at (M, x2[yi]) ----------------------------------------------
    for (int yi = 0; yi < mx; ++yi)
    {
        const double ay  = x2[yi];
        const double z1b = becdf(x, M, ay      );
        const double z4b = becdf(x, M, ay - tol);
        if (z1b != z4b) {
            const double zyb = becdf(y, M, ay);
            if (z1b > zyb) {
                RowVectorXd row(5);
                row << M + z1b, ay + z1b, M, ay, z1b;
                push_row(projection_x, row);
            }
            if (zyb > z4b + 1e-14) {
                RowVectorXd row(5);
                row << M + z4b, ay + z4b, M, ay, z4b;
                push_row(projection_x, row);
            }
        }
    }

    // -- projection_y : 6 columns ---------------------------------------------
    // Built from the omnidirectional jumps of G_n (lower curve).
    // For each jump location (b1, b2) up to four rows may be emitted:
    //   type 2 (saddle along x2,  z2): always
    //   type 4 (convex vertex,    z4): if z4 != z2
    //   type 3 (saddle along x1,  z3): if z4 != z2 AND z2 != z3
    //   type 1 (concave vertex,   z1): always
    // where z1=Gn(b1,b2), z2=Gn(b1,b2-tol), z3=Gn(b1-tol,b2),
    //       z4=Gn(b1-tol,b2-tol).
    // The type column (col 5) is read by hsearch() to select the distance
    // formula: type 1 -> |z_x - z_y|; types 2/3 -> coordinate difference.
    // -------------------------------------------------------------------------
    MatrixXd projection_y(0, 6);

    for (int xi = 0; xi < my; ++xi)
    {
        const double b1 = y1[xi];

        for (int yi = 0; yi < my; ++yi)
        {
            const double b2 = y2[yi];

            const double z1 = becdf(y, b1,       b2      );
            const double z2 = becdf(y, b1,       b2 - tol);
            const double z3 = becdf(y, b1 - tol, b2      );
            const double z4 = becdf(y, b1 - tol, b2 - tol);

            // Omnidirectional jump condition
            if (z1 == z2 || z1 == z3) continue;

            { RowVectorXd row(6); row << b1+z2, b2+z2, b1, b2, z2, 2.0;
              push_row(projection_y, row); }

            if (z4 != z2) {
                { RowVectorXd row(6); row << b1+z4, b2+z4, b1, b2, z4, 4.0;
                  push_row(projection_y, row); }
                if (z2 != z3) {
                    { RowVectorXd row(6); row << b1+z3, b2+z3, b1, b2, z3, 3.0;
                      push_row(projection_y, row); }
                }
            }

            { RowVectorXd row(6); row << b1+z1, b2+z1, b1, b2, z1, 1.0;
              push_row(projection_y, row); }
        }
    }

    return hsearch(projection_x, projection_y);
}

// -----------------------------------------------------------------------------
// Rcpp wrapper
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
double H_stat_2s_2d(NumericMatrix x_r, NumericMatrix y_r, double tol = 1e-6)
{
    return H_stat_2s_2d_cpp(as<MatrixXd>(x_r), as<MatrixXd>(y_r), tol);
}
*/
