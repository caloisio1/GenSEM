# GenSEM 0.8.0
* Full S3 interface for `gensem_recursive` objects: `print()`, `summary()`,
  `coef()`, `vcov()`, `confint()`, `logLik()` (LR-comparable within data).
* Repository furniture: NEWS, CITATION, install instructions, CI workflow.

# GenSEM 0.7.2
* Example 2 validation rewritten as gates A-D after the JOB1 diagnosis:
  gsem's published point is the optimum of Stata's 7-point adaptive
  approximation, not of the exact likelihood (+1.33 nll, stable to Q = 201).
  End-to-end green on shipped scripts (suite 80/0/0; all three examples).

# GenSEM 0.7.1
* Post-gate fixes: spell-length exposure with integrity asserts (Ex. 1),
  AGHQ backend rule for sparse absorbing processes, criteria recalibrated
  to 1.96 x published SE, `shape` alias, `start =` exposed, educ/edu shim,
  children1 reconstruction caveat documented.

# GenSEM 0.7.0
* Three-level nesting (`nested =`, `Qw =`): cluster frailty in the hazard,
  inner integral evaluated per u-node. Renamed to GenSEM.

# GenSEM 0.6.0
* Lognormal AFT family; loglogistic branch made explicit.

# GenSEM 0.5.0
* Gamma and loglogistic AFT families (audited and simulation-validated).

# GenSEM 0.4.0
* Weibull PH family with estimated shape and left truncation; switchable
  hazard kernel.

# GenSEM 0.3.0
* Multinomial (mlogit) selection with per-category loadings;
  `split_episodes()` for piecewise-constant duration dependence.

# GenSEM 0.2.0
* Recursive mixed-family class (probit selection + Poisson hazard) via
  explicit 2D Gauss-Hermite likelihood.

# GenSEM 0.1.0
* Nonrecursive stacked-Poisson multiprocess class: `stack_processes()`,
  `fit_multiproc()` (glmmTMB / GLMMadaptive backends),
  `recover_structural()` (delta-method nlcom equivalent).
