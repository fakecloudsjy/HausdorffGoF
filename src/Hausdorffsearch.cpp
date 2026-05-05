#include <vector>
#include <memory>
#include <algorithm>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "parameter.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace std;
using namespace Eigen ;
using namespace Rcpp;
using namespace hparameter;

// Custom comparator for stable two-column descending sort
struct ColumnComparator {
  const MatrixXd& mat;
  ColumnComparator(const MatrixXd& matrix) : mat(matrix) {}
  
  bool operator()(int a, int b) const {
    if(abs(mat(a, 0)-mat(b, 0))>EPSILON) 
      return mat(a, 0) > mat(b, 0); // Descending on col 0
    return mat(a, 1) > mat(b, 1);   // Descending on col 1 as tiebreak
  }
};

MatrixXd sort_matrix(const MatrixXd& mat) {
  const int n = mat.rows();
  std::vector<int> indices(n);
  std::iota(indices.begin(), indices.end(), 0);
  std::stable_sort(indices.begin(), indices.end(), ColumnComparator(mat));
  MatrixXd sorted(mat.rows(), mat.cols());
  for(int i = 0; i < n; ++i)
    sorted.row(i) = mat.row(indices[i]);
  return sorted;
}
  

double hsearch(const MatrixXd& xmat0, const MatrixXd& ymat0) {
  // FIX H5: removed `upper` and `visit` entirely.
  // In 2D the correct dominated y-vertex for x[i+1] can have a *larger*
  // lambda[0] than ymat[upper] (i.e. it sits earlier in the sorted ymat),
  // so the upper-pointer optimisation skips it and returns a wrong vertex.
  // Always scan from j=0.
  bool type1 = false;
  double max_y = 0, h0 = 0, dist = 0, h = 0;
  int tmp = 0, tmp_l = 0;
  int xrow = xmat0.rows(), yrow = ymat0.rows();
  
  MatrixXd xmat = sort_matrix(xmat0);
  MatrixXd ymat = sort_matrix(ymat0);
  
  // frontier4: indices of all type-4 vertices on V_tilde for current x-vertex
  std::vector<int> frontier4;
  frontier4.reserve(32);

  for(int i = 0; i < xrow; i++)
  {
    tmp = -1; max_y = SMALL_init; tmp_l = 0;
    frontier4.clear();

    // ----------------------------------------------------------------
    // First pass: collect V_tilde for x-vertex A.
    // Scan dominated vertices in descending lambda[0] order, tracking
    // the running maximum lambda[1].  A vertex B is on the Pareto
    // frontier iff its lambda[1] exceeds max_y at the moment it is seen.
    //
    // Stop early on type-1 or type-2 (Proposition A.2: at most one such
    // vertex on V_tilde; Lemma A.3 applies directly -- Case 2).
    // Collect ALL type-4 frontier vertices for Case 3.
    // ----------------------------------------------------------------
    for(int j = 0; j < yrow; j++)
    {
      if( (xmat(i,0)>ymat(j,0)) && (xmat(i,1)>ymat(j,1)) )
      {
        if(ymat(j,1)>max_y)
        {
          max_y = ymat(j,1);
          tmp = j;
          if(int(ymat(j,5)+0.5)!=4)
          {
            // Type-1 or type-2 on V_tilde: Case 2 applies directly.
            type1 = true;
            tmp_l = 1;
            break;
          }
          // Type-4 on V_tilde: collect for Case 3.
          frontier4.push_back(j);
          tmp_l++;
        }
      }
    }

    if(!type1){
      // No dominated vertex at all -> Case 1.
      if(tmp_l == 0) {
        tmp = -1;
      } else {
        // ----------------------------------------------------------------
        // Second pass (Case 3): for EACH type-4 vertex on V_tilde,
        // find its two Pareto predecessors along each axis direction,
        // apply Lemma A.3 to each predecessor found, and take the
        // minimum h0 over all frontier vertices and both directions.
        // (Eq. 60 in the paper: argmin rho1 over V_tilde.)
        //
        // Predecessors need NOT be dominated by A (FIX H6): a predecessor
        // B_x has the same lambda[1] as B4 and smaller lambda[0], which
        // may still exceed lambda_A[1] or lambda_A[0].
        // ----------------------------------------------------------------
        h0 = xmat(i,4);  // Case 1 default; overwritten if any predecessor found
        bool any_pred = false;

        for(int fi = 0; fi < (int)frontier4.size(); fi++)
        {
          int k = frontier4[fi];
          double k_lam0 = ymat(k,0);
          double k_lam1 = ymat(k,1);
          bool found_x = false, found_y = false;
          int  px = -1, py = -1;
          double best_dist_x = LARGE_init, best_dist_y = LARGE_init;

          for(int j = 0; j < yrow; j++)
          {
            // Predecessor along lambda[0] axis: same lambda[1], smaller lambda[0]
            if((abs(ymat(j,1)-k_lam1)<EPSILON) && (ymat(j,0)<k_lam0-EPSILON))
            {
              double d = max(abs(xmat(i,0)-ymat(j,0)), abs(xmat(i,1)-ymat(j,1)));
              if(d < best_dist_x) { best_dist_x = d; px = j; found_x = true; }
            }
            // Predecessor along lambda[1] axis: same lambda[0], smaller lambda[1]
            if((abs(ymat(j,0)-k_lam0)<EPSILON) && (ymat(j,1)<k_lam1-EPSILON))
            {
              double d = max(abs(xmat(i,0)-ymat(j,0)), abs(xmat(i,1)-ymat(j,1)));
              if(d < best_dist_y) { best_dist_y = d; py = j; found_y = true; }
            }
          }

          // Apply Lemma A.3 to the closer predecessor (argmin rho1)
          if(found_x || found_y)
          {
            int best = (best_dist_x <= best_dist_y) ? px : py;
            double candidate = ((xmat(i,0)+ymat(best,1))<(ymat(best,0)+xmat(i,1))) ?
                               abs(xmat(i,2)-ymat(best,2)) :
                               abs(xmat(i,3)-ymat(best,3));
            if(!any_pred || candidate < h0) { h0 = candidate; }
            any_pred = true;
          }
        }

        // If no predecessor found for any frontier vertex -> Case 1
        if(!any_pred) tmp = -1;
        else          tmp = 0;  // sentinel: h0 already set, skip distance block below
      }
    }
    // ----------------------------------------------------------------
    // Compute h0 for this x-vertex
    // ----------------------------------------------------------------
    if( tmp == -1 )
    {
      // Case 1: no dominated y-vertex, or Case 3 with no predecessor found
      h0 = xmat(i,4);
    } else if(type1) {
      // Case 2: type-1 or type-2 on V_tilde -- direct Lemma A.3
      if(int(ymat(tmp,5)+0.5)==1)
      {
        h0 = abs(xmat(i,4)-ymat(tmp,4));
      } else {
        h0 = ((xmat(i,0)+ymat(tmp,1))<(ymat(tmp,0)+xmat(i,1))) ?
             abs(xmat(i,2)-ymat(tmp,2)) :
             abs(xmat(i,3)-ymat(tmp,3));
      }
    }
    // Case 3: h0 already computed in the frontier4 loop above

    h = max(h, h0);

    type1 = false;
  }
  return h;
}



/*
#include <vector>
#include <memory>
#include <algorithm>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "parameter.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace std;
using namespace Eigen ;
using namespace Rcpp;
using namespace hparameter;

// Custom comparator for stable two-column descending sort
struct ColumnComparator {
  const MatrixXd& mat;
  ColumnComparator(const MatrixXd& matrix) : mat(matrix) {}
  
  bool operator()(int a, int b) const {
    // First compare column 1 (index 1)
    if(abs(mat(a, 0)-mat(b, 0))>EPSILON) 
      return mat(a, 0) > mat(b, 0); // Descending
    // If equal, compare column 2 (index 0)
    return mat(a, 1) > mat(b, 1); // Descending
  }
};

MatrixXd sort_matrix(const MatrixXd& mat) {
  const int n = mat.rows();
  std::vector<int> indices(n);
  std::iota(indices.begin(), indices.end(), 0);
  
  // Single efficient sort using custom comparator
  std::stable_sort(indices.begin(), indices.end(), ColumnComparator(mat));
  
  // Apply permutation
  MatrixXd sorted(mat.rows(), mat.cols());
  for(int i = 0; i < n; ++i)
    sorted.row(i) = mat.row(indices[i]);
  
  return sorted;
}
  
/*
double hsearch(const MatrixXd& xmat0, const MatrixXd& ymat0) {
  int upper = 0;
  bool visit = false, type1 = false, visit_x = false, visit_y = false;
  double max_y = 0, max_y0 = 0, max_x0 = 0, min_d = 0, h0 = 0, dist = 0, h = 0;
  int tmp = 0, min_v = 0, tmp_l = 0, tmp_x = 0, tmp_y = 0;
  int xrow = xmat0.rows(), yrow = ymat0.rows();
  
  MatrixXd xmat = sort_matrix(xmat0);
  MatrixXd ymat = sort_matrix(ymat0);
  
  for(int i = 0; i < xrow; i++)
  {
    tmp = -1; max_y = SMALL_init; min_v = 0; tmp_l = 0;  min_d = LARGE_init;
    for(int j = upper; j < yrow; j++)
    {
      if((xmat(i,0)>ymat(j,0)) && (!visit))
      {
        // update the starting point of the inner loop
        visit = true;
        upper = max(j-1,0);
      }
      if( (xmat(i,0)>ymat(j,0)) && (xmat(i,1)>ymat(j,1)) )
      {
        if(ymat(j,1)>max_y)
        {
          max_y = max(max_y,ymat(j,1));
          tmp = j;
          if(int(ymat(j,5)+0.5)!=4)
          {
            // the location has been found
            type1 = true;
            tmp_l = 1;
            break;
          }
          
          dist = max(abs(xmat(i,0)-ymat(j,0)),abs(xmat(i,1)-ymat(j,1)));
          min_v = (min_d>dist? j: min_v);
          min_d = min(dist, min_d);
          
          tmp_l++; 
        }
      }
    }

    visit = false;
    if(!type1){
      // find the corresponding sucession at both directions
      tmp = min_v;max_y0 = ymat(tmp,1); max_x0 = ymat(tmp,0);
      
      for(int j = upper; j < yrow; j++)
      {
        if( (xmat(i,0)>ymat(j,0)) && (xmat(i,1)>ymat(j,1)) )
        {
          if((!visit_x)&&(ymat(j,1)==max_y0) && (ymat(j,0)<max_x0))
          {
            tmp_x = j;

              dist = max(abs(xmat(i,0)-ymat(j,0)),abs(xmat(i,1)-ymat(j,1)));
              min_v = (min_d>dist? j: min_v);
              min_d = (min_d>dist? dist: min_d);

            visit_x = true;
          }
          
          if((!visit_y)&&(ymat(j,0)==max_x0)&&(ymat(j,1)<max_y0))
          {
            tmp_y = j;
            dist = max(abs(xmat(i,0)-ymat(j,0)),abs(xmat(i,1)-ymat(j,1)));
            min_v = (min_d>dist? j: min_v);
            min_d = (min_d>dist? dist: min_d);
            visit_y = true;
          }
        }
        if(visit_x&&visit_y) break;
      }
      if(!(visit_x||visit_y)){
        tmp = -1;
      }
    }
    type1 = false;
    
    
    if( tmp == -1 )
    { 
      h0 = xmat(i,4); 
    } else {
      if(!(visit_x||visit_y))
      {
        if(int(ymat(tmp,5)+0.5)==1)
        {
          h0 = abs(xmat(i,4)-ymat(tmp,4));
        } else {
          h0 = ((xmat(i,0)+ymat(tmp,1))<(ymat(tmp,0)+xmat(i,1)))? abs(xmat(i,2)-ymat(tmp,2)): abs(xmat(i,3)-ymat(tmp,3));
        }
      } else {
        if(max(abs(xmat(i,0)-ymat(tmp_x,0)),abs(xmat(i,1)-ymat(tmp_x,1)))>max(abs(xmat(i,0)-ymat(tmp_y,0)),abs(xmat(i,1)-ymat(tmp_y,1))))
        {
          // close to y side
          h0 = min(xmat(i,0)-ymat(tmp_y,0),xmat(i,1)-ymat(tmp_y,1));
        } else {
          // close to x side
          h0 = min(xmat(i,0)-ymat(tmp_x,0),xmat(i,1)-ymat(tmp_x,1));
        }
      }
    }
    
    h  =  max(h, h0);
    
    visit_x = false;
    visit_y = false;
  }   
  return h;
}



  

double hsearch(const MatrixXd& xmat0, const MatrixXd& ymat0) {
  // FIX H5: removed `upper` and `visit` entirely.
  // In 2D the correct dominated y-vertex for x[i+1] can have a *larger*
  // lambda[0] than ymat[upper] (i.e. it sits earlier in the sorted ymat),
  // so the upper-pointer optimisation skips it and returns a wrong vertex.
  // Always scan from j=0.
  bool type1 = false, visit_x = false, visit_y = false;
  double max_y = 0, max_y0 = 0, max_x0 = 0, min_d = 0, h0 = 0, dist = 0, h = 0;
  int tmp = 0, min_v = 0, tmp_l = 0, tmp_x = 0, tmp_y = 0;
  int xrow = xmat0.rows(), yrow = ymat0.rows();
  
  MatrixXd xmat = sort_matrix(xmat0);
  MatrixXd ymat = sort_matrix(ymat0);
  
  for(int i = 0; i < xrow; i++)
  {
    tmp = -1; max_y = SMALL_init; min_v = 0; tmp_l = 0; min_d = LARGE_init;

    // ----------------------------------------------------------------
    // First pass: find dominated y-vertex with the largest lambda[1].
    // Stop early if we hit a non-type-4 (Lemma A.3 applies directly).
    // For type-4 vertices, track the closest one in min_v.
    // ----------------------------------------------------------------
    for(int j = 0; j < yrow; j++)          // FIX H5: was `j = upper`
    {
      if( (xmat(i,0)>ymat(j,0)) && (xmat(i,1)>ymat(j,1)) )
      {
        if(ymat(j,1)>max_y)
        {
          max_y = max(max_y, ymat(j,1));
          tmp = j;
          if(int(ymat(j,5)+0.5)!=4)
          {
            // Non-type-4: Lemma A.3 applies directly
            type1 = true;
            tmp_l = 1;
            break;
          }
          dist = max(abs(xmat(i,0)-ymat(j,0)), abs(xmat(i,1)-ymat(j,1)));
          min_v = (min_d>dist ? j : min_v);
          min_d = min(dist, min_d);
          tmp_l++;
        }
      }
    }

    if(!type1){
      // FIX H4: if no dominated vertex found at all, skip predecessor
      // search.  Without this guard, min_v stays 0 and the loop below
      // uses ymat(0,...) as reference - wrong.
      if(tmp_l == 0) {
        tmp = -1;
      } else {
        // ----------------------------------------------------------------
        // Second pass (Case 3): look for the two Pareto predecessors of
        // min_v along each axis direction.
        // ----------------------------------------------------------------
        tmp = min_v;
        max_y0 = ymat(tmp,1);
        max_x0 = ymat(tmp,0);
        
        // FIX H6: predecessors of a type-4 vertex need NOT be dominated
        // by A.  A predecessor B_x has lambda[1]==max_y0 which may equal
        // or exceed lambda_A[1]; similarly B_y may have lambda[0]>=lambda_A[0].
        // Restricting to dominated vertices caused the predecessor search to
        // fail silently, falling through to Case 1 (h0=z_A) instead of
        // applying Case 3 -- this was the root cause of asymmetry.
        // Scan ALL y-vertices for predecessors; no dominance filter here.
        for(int j = 0; j < yrow; j++)
        {
            // Predecessor along the lambda[0] axis
            if((!visit_x) && (ymat(j,1)==max_y0) && (ymat(j,0)<max_x0))
            {
              tmp_x = j;
              dist = max(abs(xmat(i,0)-ymat(j,0)), abs(xmat(i,1)-ymat(j,1)));
              min_v = (min_d>dist ? j : min_v);
              min_d = (min_d>dist ? dist : min_d);
              visit_x = true;
            }
            // Predecessor along the lambda[1] axis
            if((!visit_y) && (ymat(j,0)==max_x0) && (ymat(j,1)<max_y0))
            {
              tmp_y = j;
              dist = max(abs(xmat(i,0)-ymat(j,0)), abs(xmat(i,1)-ymat(j,1)));
              min_v = (min_d>dist ? j : min_v);
              min_d = (min_d>dist ? dist : min_d);
              visit_y = true;
            }
          if(visit_x && visit_y) break;
        }
        // Neither predecessor found -> treat as Case 1
        if(!(visit_x || visit_y)){
          tmp = -1;
        }
      }
    }
    type1 = false;
    
    // ----------------------------------------------------------------
    // Compute h0 for this x-vertex
    // ----------------------------------------------------------------
    if( tmp == -1 )
    {
      // Case 1: no dominated y-vertex -> h0 = z_A
      h0 = xmat(i,4);
    } else {
      if(!(visit_x || visit_y))
      {
        // Direct Lemma A.3: type-1 uses z-distance, type-2/3 uses loc-distance
        if(int(ymat(tmp,5)+0.5)==1)
        {
          h0 = abs(xmat(i,4)-ymat(tmp,4));
        } else {
          h0 = ((xmat(i,0)+ymat(tmp,1))<(ymat(tmp,0)+xmat(i,1))) ?
               abs(xmat(i,2)-ymat(tmp,2)) :
               abs(xmat(i,3)-ymat(tmp,3));
        }
      } else {
        // Case 3: select the closer predecessor (Eq. 60 argmin rho1),
        // then apply Lemma A.3 with it.
        // FIX H3: original code used min(lambda[0]-diff, lambda[1]-diff)
        // - wrong formula.  Use Lemma A.3 (loc columns) instead.
        // Guard unvisited direction so it is never selected as `best`.
        double dist_x = visit_x ?
          max(abs(xmat(i,0)-ymat(tmp_x,0)), abs(xmat(i,1)-ymat(tmp_x,1))) :
          LARGE_init;
        double dist_y = visit_y ?
          max(abs(xmat(i,0)-ymat(tmp_y,0)), abs(xmat(i,1)-ymat(tmp_y,1))) :
          LARGE_init;

        // argmin rho1(lambda_A, lambda_B) -> pick the closer predecessor
        int best = (dist_x <= dist_y) ? tmp_x : tmp_y;

        h0 = ((xmat(i,0)+ymat(best,1))<(ymat(best,0)+xmat(i,1))) ?
             abs(xmat(i,2)-ymat(best,2)) :
             abs(xmat(i,3)-ymat(best,3));
      }
    }
    
    h = max(h, h0);
    
    visit_x = false;
    visit_y = false;
  }
  return h;
}
*/
