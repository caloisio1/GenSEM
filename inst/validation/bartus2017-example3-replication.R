# Replication of Bartus (2017), Example 3 (children2): Poisson hazard +
# MULTINOMIAL LOGIT selection (home / hospital1 / hospital2) with correlated
# mother-level random intercepts. Preparation follows the OFFICIAL do-file
# (st0481): death <- i.hospital i.edu i.age0; hospital <- distance i.edu;
# both equations over the same spell-level rows.
# Canonical data URL: http://web.uni-corvinus.hu/bartus/stata/children2.dta
# Published gsem (p. 458): hazard 1.hospital -.5085083, 2.hospital -3.024408;
#   loading of V in 2.hospital .2732314 (1.hospital fixed at 1);
#   var(U) .2992881; var(V) 13.02677; cov(V,U) .3196824.
# WARNING: var(V) ~ 13 is an EXTREME regime for non-adaptive Gauss-Hermite.
# Use generous Q (25+); drift here measures the quadrature gap vs gsem's
# adaptive method -- report it as such.
library(haven); library(GenSEM)
tol_a <- 0.15

d <- read_dta("children2.dta")
if ("educ" %in% names(d) && !"edu" %in% names(d)) d$edu <- d$educ
d$hospital <- factor(d$hospital); d$edu <- factor(d$edu)
d$age0f <- factor(d$age0)

fit <- fit_recursive_mlogit(hospital ~ distance + edu,
                            death ~ hospital + edu + age0f,
                            data_sel = d, data_haz = d,
                            id = "id", family = "poisson",
                            exposure = "dur", Q = 25)
print(fit$coef); print(fit$Sigma)
stopifnot(fit$convergence == 0,
          abs(fit$coef["haz:hospital1"] - (-0.5085083)) < tol_a,
          abs(fit$coef["haz:hospital2"] - (-3.024408)) < 0.30)
cat("\nBartus Example 3: REPLICATED within tolerance.\n")
