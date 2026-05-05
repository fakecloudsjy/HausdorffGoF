\name{H_test_2s_2d}
\alias{H_test_2s_2d}
\title{Two-sample bivariate Hausdorff goodness-of-fit test with Monte Carlo
permutation p-values}

\description{
Tests the null hypothesis \eqn{H_{0}: F(\bm{x}) = G(\bm{x})}{} for all
\eqn{\bm{x} \in \mathbb{R}^{2}}{} against the two-sided alternative
\eqn{H_{1}: F(\bm{x}) \neq G(\bm{x})}{} for at least one \eqn{\bm{x}}{},
using the Hausdorff distance between the bivariate empirical cdfs of two
independent samples.  The p-value is obtained by Monte Carlo permutation.
The test statistic is computed by \code{\link{H_stat_2s_2d}}.
}

\details{
Given two independent random samples
\eqn{\bm{X} = \{(X_{1,1}, X_{1,2}), \ldots, (X_{m,1}, X_{m,2})\}}{} and
\eqn{\bm{Y} = \{(Y_{1,1}, Y_{1,2}), \ldots, (Y_{n,1}, Y_{n,2})\}}{} from
unknown bivariate cdfs \eqn{F}{} and \eqn{G}{} respectively, the two-sample
Hausdorff (\eqn{H}{}) test statistic is
\eqn{\hat{q} = \mathcal{H}_{m,n} = H(F^{c}_{m}, G^{c}_{n})}{},
computed by \code{\link{H_stat_2s_2d}}.  By Theorem 2.5 of Dimitrova, Jia,
and Kaishev (2026b), this equals the side length of the largest axis-aligned
hypercube that can be inscribed in the region between the two bivariate
empirical cdf surfaces \eqn{F^{c}_{m}}{} and \eqn{G^{c}_{n}}{}.

\bold{Permutation p-value.}  Because the distribution of
\eqn{\mathcal{H}_{m,n}}{} depends on the unknown \eqn{F}{} and \eqn{G}{},
p-values are computed via a permutation argument (Proposition 3.5 of
Dimitrova, Jia, and Kaishev 2026b).  Let
\eqn{Z_{m+n} = \{\bm{x}_{1},\ldots,\bm{x}_{m},\bm{y}_{1},\ldots,\bm{y}_{n}\}}{}
be the pooled sample.  The permutation Hausdorff statistic
\eqn{\mathcal{H}^{\dagger}_{m,n}}{} is defined by randomly splitting
\eqn{Z_{m+n}}{} into groups of sizes \eqn{m}{} and \eqn{n}{} and its
conditional distribution given \eqn{Z_{m+n}}{} is:
\deqn{P(\mathcal{H}^{\dagger}_{m,n} > q \mid Z_{m+n})
  = \frac{1}{B} \sum_{b=1}^{B}
    \mathbf{1}\!\left\{H(F^{c\dagger}_{m}(\pi_{b}), G^{c\dagger}_{n}(\pi_{b})) > q\right\},}{}
where \eqn{B = \code{nboots}}{} random splits are used.  If the estimated
p-value is zero it is replaced by \eqn{1/(2 \times \code{nboots})}{}.

Exact enumeration of all \eqn{\frac{(m+n)!}{m!\,n!}}{} permutations is not
supported for the bivariate case due to the computational cost of
\code{\link{H_stat_2s_2d}}.  Use the unified wrapper
\code{\link{Hausdorff_test}} if scale tuning via \code{scale_psi} is needed.

\bold{Order-invariant statistic.}  When \code{invariant = TRUE} the
order-invariant statistic
\eqn{H^{\circ}_{m,n} = \max_{1 \leq i \leq 4} H(F^{c}_{m,i}, G^{c}_{n,i})}{}
(Equation (17) of Dimitrova, Jia, and Kaishev 2026b) is used for both the
observed value and every permutation replicate.  It maximises over the
\eqn{2^{2} = 4}{} sign-flip orderings of the two coordinates, making the
test invariant to the labelling of coordinate axes.

\bold{Asymptotic equivalence (Theorems 3.2--3.4).}  Under the null
hypothesis, or under fixed or contiguous alternatives, the permutation p-value
is asymptotically equivalent to the true p-value of \eqn{\mathcal{H}_{m,n}}{}.
The permutation test controls the type I error at the nominal level for all
\eqn{k \geq 1}{}, and the power of the two tests are asymptotically equal.

\bold{Scale tuning.}  The power of \eqn{\mathcal{H}_{m,n}}{} is not
scale-invariant.  The optimal column-wise scale vector
\eqn{(\sigma^{*}_{1}, \sigma^{*}_{2})}{} can be estimated and applied
automatically by passing \code{scale_psi} to the unified wrapper
\code{\link{Hausdorff_test}}, which calls this function internally on the
already-scaled samples.
}

\usage{
H_test_2s_2d(x, y, nboots = 2000, invariant = FALSE, tol = 1e-6)
}

\arguments{
\item{x}{a two-column numeric matrix representing the first bivariate sample
of size \eqn{m}{}.}
\item{y}{a two-column numeric matrix representing the second bivariate sample
of size \eqn{n}{}.}
\item{nboots}{a positive integer giving the number of Monte Carlo permutation
replications.  Defaults to \code{2000}.}
\item{invariant}{logical.  When \code{TRUE}, the order-invariant statistic
\eqn{H^{\circ}_{m,n}}{} is used for the observed value and all permutation
replicates, maximising over the four sign-flip orderings of the two
coordinates (Equation (17) of Dimitrova, Jia, and Kaishev 2026b).
Defaults to \code{FALSE}.}
\item{tol}{a numeric value giving the tolerance for locating omnidirectional
jump vertices in the bivariate empirical cdfs, passed to
\code{\link{H_stat_2s_2d}}.  Defaults to \code{1e-6}.}
}

\value{
An object of class \code{"htest"} with the following components:
\describe{
\item{\code{statistic}}{the observed Hausdorff statistic
\eqn{\hat{q} = \mathcal{H}_{m,n} = H(F^{c}_{m}, G^{c}_{n})}{} (or
\eqn{H^{\circ}_{m,n}}{} when \code{invariant = TRUE}), named \code{"H"}.}
\item{\code{p.value}}{the Monte Carlo permutation p-value, i.e. the proportion
of permutation statistics greater than or equal to \eqn{\hat{q}}{}.  If
zero, replaced by \eqn{1/(2 \times \code{nboots})}{}.}
\item{\code{method}}{a character string describing the procedure, with the
suffix \code{", order-invariant"} appended when \code{invariant = TRUE}.}
\item{\code{alternative}}{the character string \code{"two-sided"}.}
\item{\code{data.name}}{a character string giving the names of the two data
objects.}
}
}

\references{
Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2026b). ``On a Two-Sample
Multivariate Goodness-of-Fit Test Based on the Hausdorff Metric''.
\emph{Submitted to the Annals of Statistics}.

Dimitrina S. Dimitrova, Yun Jia, Vladimir K. Kaishev (2026a). ``On a One
Sample Goodness-of-fit Test Based on the Hausdorff Metric''. \emph{Submitted}.
}

\seealso{
\code{\link{H_stat_2s_2d}},
\code{\link{Hausdorff_test}},
\code{\link{H_test_2s_1d}}.
}

\examples{
## Bivariate H test: both samples from N(0,I) (null is true)
set.seed(1)
xm <- matrix(rnorm(100), ncol = 2)
ym <- matrix(rnorm(100), ncol = 2)
H_test_2s_2d(xm, ym, nboots = 1000)

## Bivariate H test: samples from N(0,I) and N(0.5,I) (null is false)
set.seed(2)
xm2 <- matrix(rnorm(100), ncol = 2)
ym2 <- matrix(rnorm(100, mean = 0.5), ncol = 2)
H_test_2s_2d(xm2, ym2, nboots = 1000)

## Order-invariant statistic
set.seed(3)
H_test_2s_2d(xm, ym, nboots = 1000, invariant = TRUE)

## Scale-tuned test via the unified wrapper (recommended)
set.seed(4)
H_test_2s_2d(xm, ym, nboots = 1000)
Hausdorff_test(xm, ym, nboots = 1000,
               scale_psi = c(0.99, 0.95), scale_nperms = 500)
}
