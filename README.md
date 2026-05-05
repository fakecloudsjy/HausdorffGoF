# HausdorffGoF

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/HausdorffGoF)](https://CRAN.R-project.org/package=HausdorffGoF)
[![R-CMD-check](https://github.com/d-dimitrova/HausdorffGoF/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/d-dimitrova/HausdorffGoF/actions/workflows/R-CMD-check.yaml)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
<!-- badges: end -->

## Overview

The package `HausdorffGoF` computes the test statistic and p-values of the one-sample and two-sample Hausdorff (H) goodness-of-fit tests proposed by Dimitrova, Jia, and Kaishev (2026a, 2026b).

The H statistic measures the Hausdorff distance between the completed graphs of two cumulative distribution functions (cdfs) under the Chebyshev (l∞) metric. Geometrically, it equals the side length of the largest axis-aligned square (or hypercube in higher dimensions) that can be inscribed in the region between the two cdf curves; in the univariate case this coincides with the Lévy metric. Three testing settings are covered:

- **One-sample**: compares a random sample against a fully specified continuous null cdf *F(x)*
- **Two-sample univariate**: tests whether two independent univariate samples share the same unknown distribution
- **Two-sample bivariate**: extends the two-sample test to two-dimensional data

A key advantage of the H test is that its p-values are **not scale-invariant**. This property can be exploited to direct the test's sensitivity towards the left tail, body, or right tail of the distribution by tuning a scale parameter σ. Using this, the H test significantly outperforms classical alternatives such as Kolmogorov-Smirnov, Cramér-von Mises, and Anderson-Darling in the tails, and achieves power comparable to the Wasserstein test. This makes it particularly appealing in **finance, economics, astronomy, and extreme value theory**.

## Installation

The package is available on CRAN and can be installed with:

```r
install.packages("HausdorffGoF")
```

To install the latest development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("d-dimitrova/HausdorffGoF")
```

Building from source requires a C++ compiler. On Windows, install [Rtools](https://cran.r-project.org/bin/windows/Rtools/); on macOS, install Xcode Command Line Tools (`xcode-select --install`).

## Usage

### Specifying a null distribution (one-sample)

```r
library(HausdorffGoF)

# Bundle CDF, quantile function, and density into a NullDist object
null_exp <- distribution(
  CDF        = function(t) pexp(t, rate = 1),
  CDFinverse = function(p) qexp(p, rate = 1),
  pdf        = function(t) dexp(t, rate = 1)
)
```

### One-sample H test

```r
set.seed(1)
x <- rexp(50, rate = 1)

# Exact p-value
Hausdorff_test(x, null_exp)

# Scale-tuned: focus power on the right tail
Hausdorff_test(x, null_exp, scale_psi = c(0.99, 0.95))
```

### Two-sample univariate H test

```r
set.seed(2)
x1 <- rnorm(40)
x2 <- rnorm(40, mean = 0.5)

# Monte Carlo permutation (default)
Hausdorff_test(x1, x2)

# Exact permutation (small samples)
Hausdorff_test(rnorm(8), rnorm(8), method = "exact")

# Scale-tuned: focus on right tail
Hausdorff_test(x1, x2, scale_psi = c(0.99, 0.95))
```

### Two-sample bivariate H test

```r
set.seed(3)
xm <- matrix(rnorm(100), ncol = 2)
ym <- matrix(rnorm(100, mean = 0.5), ncol = 2)

# Monte Carlo permutation
Hausdorff_test(xm, ym, nboots = 1000)

# Order-invariant statistic
Hausdorff_test(xm, ym, nboots = 1000, invariant = TRUE)
```

### Computing the statistic only

```r
# One-sample
Hausdorff_stat(x, null_exp)

# Two-sample univariate
Hausdorff_stat(x1, x2)

# Two-sample bivariate
Hausdorff_stat(xm, ym)
```

## Functions

| Function | Description |
|---|---|
| `Hausdorff_test()` | Unified S3 generic: one-sample, two-sample univariate, and bivariate H test |
| `Hausdorff_stat()` | Unified S3 generic: computes the H statistic without a p-value |
| `distribution()` | Constructs a `NullDist` object bundling CDF, quantile function, and density |
| `H_stat_1s_1d()` | One-sample H statistic (Algorithm 1, Dimitrova et al. 2026a) |
| `H_test_1s_1d()` | One-sample H test with exact p-value |
| `H_test_c_cdf()` | Low-level rectangle-probability p-value engine |
| `H_stat_2s_1d_tr()` | Two-sample univariate H statistic via transformation method (C++) |
| `H_stat_2s_1d_p()` | Two-sample univariate H statistic via projection method (C++) |
| `H_test_2s_1d()` | Two-sample univariate H test (exact or Monte Carlo permutation) |
| `H_stat_2s_2d()` | Two-sample bivariate H statistic (projection method, C++) |
| `H_test_2s_2d()` | Two-sample bivariate H test (exact or Monte Carlo permutation) |

## Methodology

**One-sample exact p-values** are expressed as a rectangle-crossing probability for uniform order statistics (Theorem 4, Dimitrova et al. 2026a) and evaluated in O(n² log n) time using the Exact-KS-FFT method of Dimitrova, Kaishev, and Tan (2020) via the `KSgeneral` package.

**Two-sample p-values** are computed by permutation (Proposition 3.5, Dimitrova et al. 2026b): exact enumeration of all C(m+n, m) splits for small samples, or Monte Carlo permutation otherwise.

**Scale tuning** via the `scale_psi = c(ψ₁, ψ₂)` argument computes the optimal scale σ* = (ψ₁ − ψ₂) / (F⁻¹(ψ₁) − F⁻¹(ψ₂)) and rescales the data before testing, directing power towards the left tail (`c(0.05, 0.01)`), body (`c(0.70, 0.40)`), or right tail (`c(0.99, 0.95)`).

**Bivariate empirical cdf computation** uses the fast divide-and-conquer algorithm of Langrené and Warin (2021), achieving O(N log² N) complexity, implemented via C++ code from the [StOpt library](https://gitlab.com/stochastic-control/StOpt) (Copyright © 2020 EDF, LGPL-3).

## Citation

If you use `HausdorffGoF` in a publication, please cite:

**One-sample H test:**
> Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2026a). "On a One Sample Goodness-of-fit Test Based on the Hausdorff Metric". *Submitted*.

**Two-sample H test:**
> Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2026b). "On a Two-Sample Multivariate Goodness-of-Fit Test Based on the Hausdorff Metric". *Submitted to the Annals of Statistics*.

## Authors

- Dimitrina S. Dimitrova — City, University of London (<D.Dimitrova@city.ac.uk>)
- Yun Jia — maintainer (<yunjia2019@gmail.com>)
- Vladimir K. Kaishev — City, University of London (<Vladimir.Kaishev.1@city.ac.uk>)

## License

GPL (≥ 2). See [LICENSE](LICENSE) for details.

This package includes C++ source code from the [StOpt library](https://gitlab.com/stochastic-control/StOpt) (Copyright © 2020 EDF — Électricité de France), used under the GNU Lesser General Public License v3 (LGPL-3).
