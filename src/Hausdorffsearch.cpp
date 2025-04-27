#include <vector>
#include <memory>
#include <algorithm>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>

// [[Rcpp::depends(RcppEigen)]]

using namespace std;
using namespace Eigen ;
using namespace Rcpp;

const double small_init = -1e300, large_init = 1e300, epsilon = 1e-14;


// Custom comparator for stable two-column descending sort
struct ColumnComparator {
  const MatrixXd& mat;
  ColumnComparator(const MatrixXd& matrix) : mat(matrix) {}
  
  bool operator()(int a, int b) const {
    // First compare column 1 (index 1)
    if(abs(mat(a, 0)-mat(b, 0))>epsilon) 
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
  std::sort(indices.begin(), indices.end(), ColumnComparator(mat));
  
  // Apply permutation
  MatrixXd sorted(mat.rows(), mat.cols());
  for(int i = 0; i < n; ++i)
    sorted.row(i) = mat.row(indices[i]);
  
  return sorted;
}
  


double hsearch(const MatrixXd& xmat0, const MatrixXd& ymat0, bool boundary) {
    int upper = 0;
    bool visit = false;
    double max_y = 0, min_d = 0, h = 0, h0 = 0, dist = 0;
    int tmp = 0, min_v = 0, tmp_l = 0;
    int xrow = xmat0.rows(), yrow = ymat0.rows();
    
    MatrixXd xmat = sort_matrix(xmat0);
    MatrixXd ymat = sort_matrix(ymat0);
    for(int i = 0; i < xrow; i++)
    {
      tmp = -1; max_y = small_init; min_v = 0; tmp_l = 0; min_d = large_init;
      for(int j = upper; j < yrow; j++)
      {
        if( (xmat(i,0)>ymat(j,0)) && (xmat(i,1)>ymat(j,1)) )
        {
          if(!visit)
          {
            // update the starting point of the inner loop
            visit = true;
            upper = max(j-1,0);
          }
          if(ymat(j,1)>max_y)
          {
            max_y = max(max_y,ymat(j,1));
            if(int(ymat(j,5)+0.5)!=4)
            {
            // the location has been found
              tmp = j;
              tmp_l = 1;
              break;
            }
            if(!boundary)
            {
              dist = max(abs(xmat(i,0)-ymat(j,0)),abs(xmat(i,1)-ymat(j,1)));
              min_v = (min_d>dist? j: min_v);
              min_d = (min_d>dist? dist: min_d);
            }
            tmp_l++; 
          }
        }
      }
      if(!boundary)
      {
        tmp = (tmp_l>1? min_v: tmp);
      } else {
        tmp = ((tmp_l!=1)? -1: (int(ymat(tmp,5)+0.5)==4? -1 :tmp));
      } 
      
      visit = false;
      if( tmp == -1)
      { 
        h0 = xmat(i,4); 
      } else {
        if(int(ymat(tmp,5)+0.5)==1)
        {
          h0 = abs(xmat(i,4)-ymat(tmp,4));
        } else {
          h0 = ((xmat(i,0)+ymat(tmp,1))>(ymat(tmp,0)+xmat(i,1)))? abs(xmat(i,2)-ymat(tmp,2)): abs(xmat(i,3)-ymat(tmp,3));
        }
      }
      h = max(h,  h0);
    }   
    return h;
}





// [[Rcpp::export]]
double hsearch_Rcpp(NumericMatrix x_proj, NumericMatrix y_proj, bool boundary) {
    // Convert R matrices to Eigen matrices
    const MatrixXd xmat = Rcpp::as<MatrixXd>(x_proj);
    const MatrixXd ymat = Rcpp::as<MatrixXd>(y_proj);
    
    // Perform operation using Eigen
    double result = hsearch(xmat, ymat, boundary);
    
    // Convert back to R matrix
    return result;
}

