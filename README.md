# GenSEM

<!-- badges: start -->
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21513948-blue.svg)](https://doi.org/10.5281/zenodo.21513948)
[![R-CMD-check](https://github.com/caloisio1/GenSEM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/caloisio1/GenSEM/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**DOI:** [10.5281/zenodo.21513948](https://doi.org/10.5281/zenodo.21513948) (concept, always the latest version) · v0.8.0: [10.5281/zenodo.21513949](https://doi.org/10.5281/zenodo.21513949)

**Multilevel multiprocess event-history models in R.** Systems of hazard equations estimated jointly, with person-level random effects that are correlated *across* equations — the models of Lillard (1993), following the estimation workflow of Bartus (2017).

## What it does

Many life-course questions involve two or more processes that cannot be assumed independent. Does hospital delivery affect subsequent fertility, or do women who deliver in hospital differ in unmeasured ways that also shape their fertility? Does leaving school raise the hazard of a first birth, or does a common unobserved disposition drive both? Estimating the two processes separately answers neither question: it attributes to the covariate what may belong to the selection.

GenSEM estimates such systems jointly. Each equation carries its own random effect; the covariance between those random effects is a free parameter, and it is what measures the endogenous selection. Duration processes are handled as piecewise-constant exponential hazards estimated as stacked Poisson regressions with log-exposure offsets (Holford, 1980), or through parametric survival families.

Two model classes are supported:

- **Nonrecursive systems** — several hazard equations with correlated random effects, where the structural effect of one process on another is recovered from the reduced form under exclusion restrictions (Bartus, 2017, eqs. 5–6).
- **Recursive systems with an endogenous exposure** — a binary or unordered categorical selection equation and a hazard equation, fitted by an explicit joint likelihood with correlated random intercepts (Bartus, 2017, eq. 7).

## Why this package

Models of this class were long the domain of aML (Lillard & Panis), now discontinued, and are today most commonly fitted with Stata's `gsem`. R has excellent tools for *shared*-frailty survival models — where one process carries one random effect — but no dedicated implementation of *multiprocess* models, where distinct equations carry distinct, correlated random effects.

GenSEM aims to make that class routine rather than artisanal:

- **A purpose-built workflow.** Exposure construction, process stacking, and the recovery of structural effects from the reduced form are functions, not hand-written data manipulation. The exposure step in particular asserts its own correctness — a mis-specified person-period exposure is the most common source of silent error in this class of model.
- **An inspectable likelihood.** The recursive fitters evaluate a joint likelihood written out explicitly and integrated by Gauss–Hermite quadrature over the random effects. Nothing is delegated to a black box: the integrand can be read, and the number of quadrature nodes is under user control.
- **Documented estimation trade-offs.** Where an approximation degrades, the package says so and ships the remedy, with the supporting simulation pinned in the test suite (see *Choosing a backend*).
- **Validation against published results.** All three worked examples of Bartus (2017) are replication targets; the scripts are shipped in `inst/validation/`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("caloisio1/GenSEM")
```

GenSEM is not yet on CRAN.

## Quick start

A recursive system: an endogenous binary exposure (`hospital`, instrumented by `distance`) and a hazard for `death`.

```r
library(GenSEM)

fit <- fit_recursive(hospital ~ distance + edu,   # selection equation
                     death ~ hospital + edu,      # hazard equation
                     data_sel = dsel, data_haz = dhaz,
                     id = "id", exposure = "dur", Q = 15)

summary(fit)
coef(fit)["haz:hospital"]  # structural effect of the endogenous exposure
fit$Sigma                  # random-effect covariance: the selection parameter
confint(fit)
```

The fitted object supports `print()`, `summary()`, `coef()`, `vcov()`, `confint()` and `logLik()`.

Identification requires an excluded instrument in the selection equation — a variable that shifts the exposure and is excluded from the hazard. The package cannot supply one, and `recover_structural()` will compute what it is asked to compute whether or not the system is identified.

A nonrecursive system follows the same shape:

```r
long <- stack_processes(spells, id = "id",
                        events = c(birth2 = "birth2", divorce = "divorce"),
                        exposure = "dur")

fit <- fit_multiproc(long, ~ hereduc3 + age + mardur, id = "id")

recover_structural(fit,    # delta-method SE, equivalent to Stata's nlcom
  "`processbirth2:hereduc3` - (`processbirth2:mardur`/`processdivorce:mardur`) * `processdivorce:hereduc3`")
```

## Features

| | |
|---|---|
| **Duration processes** | Piecewise-constant exponential hazards as stacked Poisson; `split_episodes()` for episode splitting and non-monotone duration dependence |
| **Survival families** | `weibull` (proportional hazards); `gamma`, `loglogistic`, `lognormal` (accelerated failure time); the exponential case is the piecewise-constant `poisson` route. Coefficients are not comparable across the two metrics |
| **Censoring** | Right censoring and left truncation (delayed entry) throughout |
| **Endogenous exposures** | Binary via probit selection (`fit_recursive()`); *C*-category unordered via multinomial logit with per-category loadings (`fit_recursive_mlogit()`) |
| **Multilevel structure** | Two levels by default; three levels via `fit_recursive(..., nested = "cluster_col", Qw = 7)` — spells within clusters within units, with a cluster-level frailty on the hazard |
| **Integration** | Gauss–Hermite quadrature with user-controlled nodes; nested quadrature for three-level models |

## Choosing a backend

`fit_multiproc()` offers two estimation backends:

- `backend = "glmmTMB"` — Laplace approximation. Fast, and the default.
- `backend = "GLMMadaptive"` — adaptive Gauss–Hermite quadrature; the algorithm of Rabe-Hesketh, Skrondal & Pickles (2005).

The choice matters in sparse regimes. With one spell per person, the Laplace approximation overestimates the random-effect variance; the test suite pins both the bias (seeded truth var(u) = 0.3, n = 4000) and the correction — at n = 2000 the AGHQ estimate with 11 nodes is asserted closer to the seeded truth than Laplace. Prefer AGHQ when spells per person are few or when events are strongly absorbing; note that its cost grows as *nAGQ^k* with *k* processes.

The recursive fitters are transparent rather than fast: BFGS with numerical derivatives over *Q²* nodes takes minutes, not seconds.

## Scope

GenSEM implements the multiprocess / event-history subclass of generalized structural equation modeling in the sense of Skrondal & Rabe-Hesketh (2004) — systems of simultaneous equations with latent variables (factors, random effects and frailties being the same object) over non-Gaussian responses. It is not a port of Stata's `gsem` and does not target its full family/link matrix.

Deliberately out of scope:

- **Latent measurement models.** That is lavaan's domain. The workflow bridge, when needed, is to fit the measurement model in lavaan, export factor scores and use them as GenSEM covariates — with the standard caveat that factor scores carry measurement error, so Croon-style corrections or plausible values are advisable for serious work.
- **Survey weights.** Weighted multilevel pseudo-likelihood with level-specific scaling is not implemented. Complex-survey designs require design covariates plus sensitivity analysis, or another tool.
- **Four or more levels.** A fourth level breaks the independence of top-level units and makes the integrator recursive — a different architecture rather than one more loop — while non-adaptive quadrature degrades and cost multiplies per level. At that point the honest tool is MCMC (e.g. **brms**), where additional levels are nearly free.

## Validation

- **Simulation suite:** 92 assertions, 0 failures (Windows, R 4.6.1). Seeded-truth recovery for every family, backend and nesting mode.
- **Published-results gate:** all three worked examples of Bartus (2017).
  - *Example 1* — six fixed effects replicate to a maximum absolute difference of 0.017; variance components fall inside the published confidence intervals.
  - *Example 3* — passes its canonical criteria.
  - *Example 2* — replicates the variance component within the pre-registered band. The coefficient on the endogenous exposure differs from the published point estimate (−0.788 against −0.513). Under exact-likelihood evaluation, converged for Q ≥ 61, our estimate dominates the published one by 1.33 log-likelihood units, while the published single-equation hazard model (no selection equation) replicates to within 0.0002 of its published −0.382. The published value corresponds to the optimum of the default 7-node adaptive quadrature rather than of the exact likelihood. Replication criteria therefore use published standard-error bands, not point tolerances, and the divergence is documented rather than tuned away.
- **Reproducibility note:** `children1.dta` is currently in no public archive; the Example 2 run used a verified reconstruction, documented in the script header.

To run the full suite, install with `--install-tests` and set `NOT_CRAN=true`; otherwise the heavy tests skip silently.

## Related work

For *shared*-frailty survival models — one process, one random effect — R is already well served by **survival**, **coxme**, **frailtypack**, **parfm** and **frailtyEM**; GenSEM is not an alternative to these. For multiprocess systems, a determined user can assemble something similar by stacking processes in **glmmTMB** with an unstructured process-level covariance, or by specifying a multivariate model with correlated group-level effects in **brms**. GenSEM's contribution is to make that construction a documented, validated workflow rather than a bespoke one, with the pitfalls it encountered along the way encoded as assertions.

## Citation

```r
citation("GenSEM")
```

Aloisio, C. (2026). *GenSEM: Multilevel Multiprocess Hazard Models in R*. R package version 0.8.0. https://doi.org/10.5281/zenodo.21513949

## References

Bartus, T. (2017). Multilevel multiprocess modeling with gsem. *The Stata Journal*, 17(2), 442–461. https://doi.org/10.1177/1536867X1701700211

Holford, T. R. (1980). The analysis of rates and of survivorship using log-linear models. *Biometrics*, 36(2), 299–305. https://doi.org/10.2307/2529982

Lillard, L. A. (1993). Simultaneous equations for hazards: Marriage duration and fertility timing. *Journal of Econometrics*, 56(1–2), 189–217. https://doi.org/10.1016/0304-4076(93)90106-F

Rabe-Hesketh, S., Skrondal, A., & Pickles, A. (2005). Maximum likelihood estimation of limited and discrete dependent variable models with nested random effects. *Journal of Econometrics*, 128(2), 301–323. https://doi.org/10.1016/j.jeconom.2004.08.017

Skrondal, A., & Rabe-Hesketh, S. (2004). *Generalized Latent Variable Modeling: Multilevel, Longitudinal, and Structural Equation Models*. Chapman & Hall/CRC. https://doi.org/10.1201/9780203489437

Croon, M. (2002). Using predicted latent scores in general latent structure models. In G. A. Marcoulides & I. Moustaki (Eds.), *Latent Variable and Latent Structure Models* (pp. 195–223). Lawrence Erlbaum.

Brooks, M. E., et al. (2017). glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. *The R Journal*, 9(2), 378–400. https://doi.org/10.32614/RJ-2017-066
