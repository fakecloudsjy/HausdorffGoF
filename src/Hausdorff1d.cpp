#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <map>
#include <cmath>


using namespace Rcpp;

// [[Rcpp::export]]
double Hausdorff45_Rcpp(NumericVector a, NumericVector b) {
  // Step 1: Compute frequency counts for unique elements
  std::map<double, int> a_counts, b_counts;
  for (double val : a) a_counts[val]++;
  for (double val : b) b_counts[val]++;
  
  int m_a = a.size(), m_b = b.size();
  
  // Extract sorted unique values and their proportions
  std::vector<double> a_unique, b_unique;
  std::vector<double> ay1, by1;
  
  for (auto const &pair : a_counts) {
    a_unique.push_back(pair.first);
    ay1.push_back(pair.second / static_cast<double>(m_a));
  }
  for (auto const &pair : b_counts) {
    b_unique.push_back(pair.first);
    by1.push_back(pair.second / static_cast<double>(m_b));
  }
  
  std::sort(a_unique.begin(), a_unique.end());
  std::sort(b_unique.begin(), b_unique.end());
  
  // Swap a and b if a has more unique elements
  bool swapped = false;
  if (a_unique.size() > b_unique.size()) {
    swapped = true;
    std::swap(a_unique, b_unique);
    std::swap(ay1, by1);
    std::swap(m_a, m_b);
  }
  
  // Adjust values by subtracting the minimum
  double min_val = std::min(a_unique.front(), b_unique.front());
  for (auto &val : a_unique) val -= min_val;
  for (auto &val : b_unique) val -= min_val;
  
  double max_val = std::max(a_unique.back(), b_unique.back());
  
  // Build aa and bb vectors
  int a_len = a_unique.size(), b_len = b_unique.size();
  std::vector<double> aa(2 * a_len, 0.0), bb(2 * b_len, 0.0);
  
  for (int i = 0; i < a_len; ++i) {
    double diff = (i == 0) ? a_unique[i] : (a_unique[i] - a_unique[i-1]);
    aa[2*i] = -diff;
    aa[2*i + 1] = ay1[i];
  }
  for (int i = 0; i < b_len; ++i) {
    double diff = (i == 0) ? b_unique[i] : (b_unique[i] - b_unique[i-1]);
    bb[2*i] = -diff;
    bb[2*i + 1] = by1[i];
  }
  
  // Adjust aa and bb if a's last element is max_val
  int m = 2 * a_len;
  if (!a_unique.empty() && a_unique.back() == max_val) {
    aa.resize(2 * a_len - 1);
    m = aa.size();
    if (bb.size() >= 2 * b_len && bb[2*b_len - 1] != -max_val) {
      bb.push_back(max_val);
    }
  }
  
  // Compute cumulative sums ax, ay, bx, by
  std::vector<double> ax, ay, bx, by;
  double sum_ax = 0.0, sum_ay = 0.0;
  for (double v : aa) {
    sum_ax += std::abs(v);
    sum_ay += v;
    ax.push_back(sum_ax);
    ay.push_back(sum_ay);
  }
  double sum_bx = 0.0, sum_by = 0.0;
  for (double v : bb) {
    sum_bx += std::abs(v);
    sum_by += v;
    bx.push_back(sum_bx);
    by.push_back(sum_by);
  }
  
  // Find initial j
  int j = 0;
  int begining = std::min(2, m);
  if (ay.size() > 0 && ay[0] != 0.0) begining = 1;
  double target_ax = (begining > 0) ? ax[begining - 1] : 0.0;
  
  while (j < bx.size() && bx[j] <= target_ax) j++;
  
  double max_H = 0.0;
  // Compute H_temp values
  std::vector<double> H_vec;
  for (int i = begining; i <= m; ++i) {
    if (i == 0) continue;
    int idx = i - 1;
    if (idx >= ay.size()) break;
    
    int prev_j = std::max(j - 1, 0);
    if (prev_j >= by.size()) prev_j = by.size() - 1;
    
    int parity = (j % 2 == 0) ? 1 : -1;
    double H_temp = ay[idx] - by[prev_j] + parity * (ax[idx] - bx[prev_j]);
    
    if ((H_temp < 0) == (i % 2 != 0)) {
      H_vec.push_back(H_temp);
    }
    
    if (i < m) {
      double next_ax = ax[i];
      while (j < bx.size() && bx[j] <= next_ax) j++;
    }
    max_H = std::max(abs(H_temp),max_H);
  }
  
  return max_H / 2.0;
}

// Helper function to compute ECDF values
std::vector<double> compute_ecdf(const std::vector<double>& data, 
                                const std::vector<double>& points) {
    std::vector<double> ecdf_values;
    for (double p : points) {
        int count = 0;
        for (double d : data) {
            if (d <= p) count++;
        }
        ecdf_values.push_back(static_cast<double>(count) / data.size());
    }
    return ecdf_values;
}

// [[Rcpp::export]]
double Hausdorff_Cpp(NumericVector a, NumericVector b) {
    // Convert inputs to sorted unique vectors
    std::vector<double> a_vec(a.begin(), a.end());
    std::vector<double> b_vec(b.begin(), b.end());
    
    std::sort(a_vec.begin(), a_vec.end());
    std::sort(b_vec.begin(), b_vec.end());
    
    auto last_a = std::unique(a_vec.begin(), a_vec.end());
    auto last_b = std::unique(b_vec.begin(), b_vec.end());
    a_vec.erase(last_a, a_vec.end());
    b_vec.erase(last_b, b_vec.end());

    // Create joint sorted vector
    std::vector<double> xjoint(a_vec);
    xjoint.insert(xjoint.end(), b_vec.begin(), b_vec.end());
    std::sort(xjoint.begin(), xjoint.end());
    auto last = std::unique(xjoint.begin(), xjoint.end());
    xjoint.erase(last, xjoint.end());

    // Compute ECDF values
    std::vector<double> ay1 = compute_ecdf(a_vec, xjoint);
    std::vector<double> ay2 = compute_ecdf(b_vec, xjoint);

    // Build x1/y1 and x2/y2 vectors
    std::vector<double> x1, y1, x2, y2;
    double tempy1 = 0.0, tempy2 = 0.0;
    
    for (size_t i = 0; i < xjoint.size(); ++i) {
        double current_max = std::max(ay1[i], ay2[i]);
        if (current_max > tempy1) {
            tempy1 = current_max;
            x1.push_back(xjoint[i]);
            y1.push_back(tempy1);
        }
        
        double current_min = std::min(ay1[i], ay2[i]);
        if (current_min > tempy2) {
            tempy2 = current_min;
            x2.push_back(xjoint[i]);
            y2.push_back(tempy2);
        }
    }

    // Extend x1/y1 and x2/y2
    size_t m1 = x1.size();
    size_t m2 = x2.size();
    std::vector<double> new_x1, new_y1, new_x2, new_y2;

    // Process x1/y1
    for (size_t i = 0; i < 2*m1; ++i) {
        size_t idx = (i + 1)/2 - 1;
        new_x1.push_back(x1[idx]);
        new_y1.push_back(y1[idx]);
    }
    new_x1.push_back(xjoint.back());
    new_y1.insert(new_y1.begin(), 0.0);

    // Process x2/y2
    new_x2.push_back(xjoint[0]);
    for (size_t i = 0; i < 2*m2; ++i) {
        size_t idx = (i + 1)/2 - 1;
        new_x2.push_back(x2[idx]);
        new_y2.push_back(y2[idx]);
    }
    new_y2.push_back(1.0);

    // Update sizes
    m1 = new_x1.size();
    m2 = new_x2.size();

    // Swap logic
    if (m1 > m2) {
        std::swap(new_x1, new_x2);
        std::swap(new_y1, new_y2);
        std::reverse(new_x1.begin(), new_x1.end());
        std::reverse(new_x2.begin(), new_x2.end());
        std::reverse(new_y1.begin(), new_y1.end());
        std::reverse(new_y2.begin(), new_y2.end());
        
        for (auto& x : new_x1) x = -x;
        for (auto& x : new_x2) x = -x;
        for (auto& y : new_y1) y = 1.0 - y;
        for (auto& y : new_y2) y = 1.0 - y;
        
        std::swap(m1, m2);
    }

    // Calculate parameters
    std::vector<double> par1, par2;
    for (size_t i = 1; i < m1; i += 2) {
        par1.push_back(new_x1[i] + new_y1[i]);
    }
    for (size_t i = 0; i < m2; ++i) {
        par2.push_back(new_x2[i] + new_y2[i]);
    }

    // Match parameters
    std::vector<double> combined(par1);
    combined.insert(combined.end(), par2.begin(), par2.end());
    std::sort(combined.begin(), combined.end());
    
    std::vector<int> m;
    for (double p : par1) {
        auto it = std::lower_bound(combined.begin(), combined.end(), p);
        m.push_back(it - combined.begin() + 1);
    }

    // Process matches
    std::map<int, int> counts;
    for (int val : m) counts[val]++;
    
    std::vector<int> temp2;
    for (auto& pair : counts) {
        if (pair.second > 1) temp2.push_back(pair.first);
    }
    
    bool Logic1 = !temp2.empty();
    if (Logic1) {
        std::vector<int> temp3;
        int sum = 0;
        for (auto& pair : counts) {
            sum += pair.second;
            temp3.push_back(sum);
        }
        // ... (additional processing for duplicates)
    }

    // Calculate batches and H_vec
    std::vector<int> temp;
    for (size_t i = 1; i < m.size(); ++i) {
        if (m[i] - m[i-1] > 1) temp.push_back(i);
    }
    temp.push_back(m.size());
    
    std::vector<double> H_vec;
    for (size_t i = 0; i < temp.size(); ++i) {
        int idx = temp[i];
        int ibatch = m[idx-1] - idx;
        if (ibatch % 2 == 1) {
            H_vec.push_back(new_y1[2*idx] - new_y2[ibatch+1]);
        } else {
            H_vec.push_back(new_x2[ibatch] - new_x1[2*idx]);
        }
    }

    return *std::max_element(H_vec.begin(), H_vec.end());
}
