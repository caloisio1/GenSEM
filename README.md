# GenSEM

**GenSEM** -- the name states the model family: GENERALIZED structural
equation modeling in the Skrondal & Rabe-Hesketh (2004) sense -- systems
of simultaneous equations with latent variables (factors, random effects,
frailties are the same object) over non-Gaussian response families. Scope
honesty in the first paragraph, where it belongs: GenSEM implements the
multiprocess / event-history SUBCLASS of that family (Lillard 1993;
Bartus 2017, whose three published examples are this package's validation
gate). It is NOT a port of Stata's gsem command and does not aim at its
full family/link matrix; what it does cover, it covers with a fully
inspectable likelihood.

**Relation to lavaan, stated precisely:** complementary in SCOPE, not in
architecture. GenSEM covers the multiprocess/event-history class that
lies outside lavaan's reach; it shares no syntax, objects, or dependency
with lavaan. The one genuine workflow bridge: fit a measurement model in
lavaan, export factor scores, use them as GenSEM covariates -- with the
standard caveat that factor scores carry measurement error (attenuation
bias; use Croon-style corrections or plausible values for serious work).

Multilevel multiprocess hazard models in R: systems of piecewise-constant
exponential hazard equations with correlated person-level random effects
(Lillard 1993), estimated as stacked Poisson regressions with log-exposure
offsets, following the `gsem` workflow of Bartus (2017, *Stata Journal* 17(2)).

Four functions. The nonrecursive pipeline:

```r
long <- stack_processes(spells, id = "id",
                        events = c(birth2 = "birth2", divorce = "divorce"),
                        exposure = "dur")
fit  <- fit_multiproc(long, ~ hereduc3 + age + mardur, id = "id")
recover_structural(fit,   # nlcom equivalent, delta-method SE
  "`processbirth2:hereduc3` - (`processbirth2:mardur`/`processdivorce:mardur`) * `processdivorce:hereduc3`")
```

And the recursive mixed-family class (Bartus 2017, eq. 7: endogenous binary
exposure; probit selection + Poisson hazard with correlated random
intercepts), fit by an explicit, fully inspectable joint likelihood
integrated with 2D Gauss-Hermite quadrature:

```r
fit <- fit_recursive(hospital ~ distance + edu,   # instrument: distance
                     death ~ hospital + edu,
                     data_sel = dsel, data_haz = dhaz,
                     id = "id", exposure = "dur", Q = 15)
fit$coef["haz:hospital"]   # structural effect of the endogenous dummy
fit$Sigma                  # random-effect covariance
```

Transparent but not fast (BFGS + numerical derivatives over Q^2 nodes:
minutes, not seconds). Identification requires an excluded instrument in the
selection equation.

`fit_recursive_mlogit()` extends the same likelihood to a C-category
unordered exposure (multinomial logit selection, per-category loadings on
the selection random effect, first loading fixed to 1 as in gsem), and
`split_episodes()` provides episode splitting for piecewise-constant
duration dependence.

## Installation

```r
# install.packages("remotes")
remotes::install_github("caloisio1/GenSEM")
```

## What it does not do (by design)

- **Three-level nesting** is implemented in
  `fit_recursive(..., nested = "cluster_col", Qw = 7)`: spells nested in
  LEVEL-2 CLUSTERS nested in LEVEL-1 UNITS (`id`). A cluster-level frailty
  w ~ N(0, sw^2) enters the HAZARD in addition to the unit-level (u, v)
  pair, integrated by nested Gauss-Hermite quadrature; the inner integral
  over w is evaluated per u-node (Q x Qw kernel passes, not Q^2 x Qw),
  exploiting the fact that u depends only on the first Cholesky axis.
  The structure is domain-agnostic: children within mothers (Bartus's
  children data), episodes within youths within schools or territories,
  jobs within workers within firms. Note where each equation lives:
  selection and the (u, v) pair sit at the `id` level, the frailty at the
  `nested` level -- map your data so the endogenous exposure is an
  `id`-level choice. Scope: frailty on the hazard side only; probit fitter
  only (three-level mlogit is a follow-up). Four or more levels are out of
  scope for three written reasons: (i) structurally, a fourth level breaks
  the independence of `id` units, turning the integrator recursive -- a
  different architecture, not one more loop; (ii) numerically, non-adaptive
  quadrature degrades and costs multiply per level, and at that point the
  honest tool is MCMC (`brms`), where extra levels are nearly free; (iii)
  empirically, the Lillard/Bartus class and its published validation
  targets live at two and three levels. The boundary is where these
  reasons coincide, not a theorem: a real four-level application with
  published numbers would reopen it.
- **Survival families -- read the metric before comparing coefficients.**
  Implemented: `weibull` (PROPORTIONAL HAZARDS: positive coefficient =
  higher hazard), and `gamma` / `loglogistic` (ACCELERATED FAILURE TIME:
  positive coefficient = longer survival = lower hazard). Do not compare
  coefficients across the two metrics. Gamma here is the two-parameter
  gamma AFT with rate exp(-lp) (mean k*exp(lp)); note Stata's `streg`
  "gamma" is the three-parameter GENERALISED gamma -- verify
  parameterisations before any cross-engine numeric comparison. `lognormal`
  (AFT, sigma estimated) completes the set: all five gsem survival
  families are now covered (exponential = weibull with shape 1).
  All survival families support right censoring and left truncation
  (`entry`); for non-monotone hazards under the Poisson family, use
  `split_episodes()` + interval dummies (Bartus 2017, p. 459).
- **Latent measurement models**: permanently out of scope -- that is
  lavaan's job; duplicating it is the ambition that sinks packages.
- **Two estimation backends, no black boxes**: `backend = "glmmTMB"`
  (Laplace approximation; fast, default) or `backend = "GLMMadaptive"`
  (adaptive Gauss-Hermite quadrature, the published algorithm of
  Rabe-Hesketh, Skrondal & Pickles 2005 -- the same one gsem implements).
  In the sparse regime (one spell per person), a 24-replication Monte Carlo
  (true var(u) = 0.30, n = 2000) gives mean bias +0.076 for Laplace and
  +0.001 for AGHQ-11; the test suite pins both facts. Prefer AGHQ when
  spells per person are few; note its cost grows as nAGQ^k with k
  processes.
- **Identification is the analyst's job**: nonrecursive systems require
  excluded instruments (Bartus 2017, eq. 5-6; Maddala's condition).
  `recover_structural()` computes what you ask; it cannot make an
  unidentified system identified.
- **Survey weights**: not supported; weighted multilevel pseudo-likelihood
  with level-specific scaling is out of scope for v0.1.

## Validation status

- **Simulation suite**: 80 assertions, 0 failures on the development container
  AND on an independent Windows/R 4.6.1 machine (seeded-truth recovery for
  every family, backend, and nesting mode; Laplace sparse-regime bias pinned).
  Operational note: install with `--install-tests`, and set
  `NOT_CRAN=true` -- otherwise heavy tests silently skip and the green is
  empty.
- **Published-results gate (Bartus 2017, all three examples; JOB1 run,
  21 Jul 2026)**: Example 1 -- all six fixed effects replicate to max diff
  0.017 with spell-length exposure and the GLMMadaptive backend (Laplace is
  ruled out there: 246 absorbing events in 2121 women inflate var(V) to 61);
  variance components fall inside the published CIs (var(V) is weakly
  identified in the original itself: SE 0.379). Example 3 (mlogit) passes
  its canonical criteria. Example 2 replicates var(V) (3.53 vs 4.15) under
  Bartus's exact likelihood design, in which the probit selection equation
  enters ONCE PER SPELL (N = 2002; survivors' rows are duplicated).
- **Selection-design semantics, stated openly**: GenSEM expresses both
  designs through `data_sel` -- the package imposes neither. The docs
  recommend one row per selection unit (duplicating probit rows inflates
  the selection information and understates its SEs); the replication
  scripts use the as-published spell-level frame, declared as such.
  Replicating is not endorsing.
- **What the gate taught about integration**: the residual gaps first
  attributed to "the quadrature gap" were likelihood-design mismatches;
  non-adaptive GH was fully converged from moderate Q (Q = 21/31/41 and
  nAGQ = 11/15/21 identical to 4 decimals). The one genuine integration
  pathology found is Laplace in sparse absorbing regimes -- remedy shipped
  (`backend = "GLMMadaptive"`). An adaptive inner quadrature for the
  recursive fitters is therefore NOT on the roadmap: the evidence removed
  its justification.
- **Data caveat**: children1.dta is currently in no public archive; the
  Example 2 run used a verified reconstruction (see the script header).
  Revalidate if the original resurfaces.

## References

Bartus (2017) *Stata Journal* 17(2):442-461. Lillard (1993) *J. Econometrics*
56:189-217. Holford (1980) *Biometrics* 36:299-305. Skrondal & Rabe-Hesketh
(2004) *Generalized Latent Variable Modeling*.
