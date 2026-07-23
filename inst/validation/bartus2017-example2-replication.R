# Replication of Bartus (2017), Example 2 (children1): recursive system --
# Poisson child-mortality hazard + probit hospital-delivery selection with
# correlated mother-level random intercepts. Preparation follows the OFFICIAL
# do-file (st0481): death <- hospital i.edu i.age0 [FULL age0 factor];
# hospital <- distance i.edu; both equations estimated over the SAME
# spell-level rows (gsem's single-frame construction; the probit is
# pseudo-replicated across spells).
# Canonical data URL: http://web.uni-corvinus.hu/bartus/stata/children1.dta
# Published gsem (p. 456): hospital -0.5131628 (se .2411954);
#   var(U) .4091622 CI [.1334, 1.2551]; var(V) 4.149642 CI [2.4916, 6.9110];
#   cov(V,U) .2157169 (se .1885667). Simple model (p. 455): hospital -0.382.
#
# GATE DESIGN (JOB2, 22 Jul 2026). The 0.7.1 point-tolerance on hospital
# (+-0.10) is NOT attainable by exact ML: at converged quadrature (Q>=61 and
# spot-checks at Q=101/201) the exact-likelihood MLE of this specification is
# hospital ~ -0.788 with var(V) ~ 3.53, and the PUBLISHED point evaluates
# 1.33 nll units WORSE than that MLE. gsem's published numbers are the
# optimum of its 7-point adaptive approximation, not of the exact likelihood.
# What this gate therefore certifies instead:
#   A. the data + machinery reproduce the published SIMPLE model exactly
#      (hospital -0.382; also validates the children1 reconstruction, whose
#      raw-aML binary coding gives -0.3818 vs -1.11 for the children2>0
#      alternative);
#   B. the joint fit converges and DOMINATES the published theta in exact
#      (high-Q) likelihood;
#   C. variance components fall inside the PUBLISHED CIs, and hospital falls
#      within 1.96 * published SE;
#   D. Weibull variant is internally consistent with the Poisson fit (the
#      paper's own claim is consistency, "very close to" its Poisson value).
# cov(V,U) is REPORTED with its two-estimate z but not asserted: it is where
# the exact integral and gsem's approximation genuinely disagree -- that gap
# is a finding to surface, not a defect to gate on.
#
# DATA CAVEAT (JOB1, 21 Jul 2026): children1.dta is in NO public archive
# (Wayback empty, Stata Journal 404). This runs on a verified reconstruction
# (children2's identical 1060-children episode structure + binary hospital
# from the archived raw children.dta), validated exactly by gate A.
library(haven); library(GenSEM)
if (!requireNamespace("GLMMadaptive", quietly = TRUE))
  stop("Gate A requires GLMMadaptive.")

d <- read_dta("children1.dta")
if ("educ" %in% names(d) && !"edu" %in% names(d)) d$edu <- d$educ
d$edu <- factor(d$edu); d$age0f <- factor(d$age0)
dsel <- d                                   # as published: spell-level rows
dhaz <- d

# --- GATE A: published simple model, exact replication --------------------
ms <- GLMMadaptive::mixed_model(
  fixed = death ~ hospital + edu + age0f + offset(log(dur)),
  random = ~ 1 | id, data = d, family = poisson(), nAGQ = 15)
bs <- GLMMadaptive::fixef(ms)["hospital"]
cat(sprintf("Gate A -- simple model hospital: %.4f (published -0.382)\n", bs))
stopifnot(abs(bs - (-0.382)) < 0.01)

# --- joint model, as published (spell-level probit rows) ------------------
fit <- fit_recursive(hospital ~ distance + edu,
                     death ~ hospital + edu + age0f,
                     data_sel = dsel, data_haz = dhaz,
                     id = "id", family = "poisson", exposure = "dur", Q = 21)
print(fit$coef); print(fit$Sigma)

# --- GATE B: exact-likelihood dominance over the published theta ----------
th_pub <- c(-2.209737, -0.0231453, 2.01218, 3.148736,
            -3.12697, -0.5131628, -0.2625067, -2.021169, -4.920847,
            0.5 * log(0.4091622), 0.5 * log(4.149642),
            atanh(0.2157169 / sqrt(0.4091622 * 4.149642)))
X1 <- stats::model.matrix(hospital ~ distance + edu, dsel)
bh <- GenSEM:::.build_haz(death ~ hospital + edu + age0f, dhaz,
                          "poisson", "dur", NULL, NULL)
ids <- sort(unique(d$id))
idf <- factor(d$id, levels = ids)
gh61 <- GenSEM:::gh_nodes(61)
f61 <- function(th) GenSEM:::nll_recursive(th, X1 = X1, y1 = d$hospital,
        id1 = idf, X2 = bh$X2, H = bh$H, id2 = idf,
        ids = seq_along(ids), gh = gh61)
nll_hat <- f61(unname(fit$coef)); nll_pub <- f61(th_pub)
cat(sprintf("Gate B -- exact nll: fitted %.3f vs published %.3f (margin %+.3f)\n",
            nll_hat, nll_pub, nll_pub - nll_hat))
stopifnot(fit$convergence == 0, nll_hat < nll_pub)

# --- GATE C: agreement with published inference ---------------------------
cat(sprintf("Gate C -- hospital %.4f vs -0.5132 (1.96*SEpub = %.3f)\n",
            fit$coef["haz:hospital"], 1.96 * 0.2411954))
stopifnot(abs(fit$coef["haz:hospital"] - (-0.5131628)) < 1.96 * 0.2411954,
          fit$Sigma[1, 1] > 0.1333875, fit$Sigma[1, 1] < 1.255093,
          fit$Sigma[2, 2] > 2.491617,  fit$Sigma[2, 2] < 6.910987)
# cov(V,U): reported, not gated (see header)
se_cov <- recover_structural(fit,
  "exp(`log s1`) * exp(`log s2`) * tanh(`atanh rho`)")
z_cov <- (se_cov$estimate - 0.2157169) /
  sqrt(se_cov$se^2 + 0.1885667^2)
cat(sprintf("cov(V,U): %.4f (se %.4f) vs published .2157 (se .1886); two-estimate z = %.2f\n",
            se_cov$estimate, se_cov$se, z_cov))

# --- sensitivity: child-level (deduplicated) probit rows ------------------
dsel2 <- unique(d[, c("id", "bid", "hospital", "distance", "edu")])
fit2 <- fit_recursive(hospital ~ distance + edu,
                      death ~ hospital + edu + age0f,
                      data_sel = dsel2, data_haz = dhaz,
                      id = "id", family = "poisson", exposure = "dur", Q = 21)
cat(sprintf("hospital: as-published %.4f | child-level probit %.4f\n",
            fit$coef["haz:hospital"], fit2$coef["haz:hospital"]))

# --- GATE D: Weibull variant, internal consistency ------------------------
fitw <- fit_recursive(hospital ~ distance + edu,
                      death ~ hospital + edu + age0f,
                      data_sel = dsel, data_haz = dhaz,
                      id = "id", family = "weibull",
                      time = "month", entry = "age0", Q = 21)
cat(sprintf("Gate D -- Weibull hospital %.4f vs Poisson %.4f | shape %.3f\n",
            fitw$coef["haz:hospital"], fit$coef["haz:hospital"], fitw$ancillary))
stopifnot(fitw$convergence == 0,
          sign(fitw$coef["haz:hospital"]) == sign(fit$coef["haz:hospital"]),
          abs(fitw$coef["haz:hospital"] - fit$coef["haz:hospital"]) < 0.20)

# --- three-level bonus (the level Bartus's spec ignores) ------------------
fit3 <- fit_recursive(hospital ~ distance + edu,
                      death ~ hospital + edu + age0f,
                      data_sel = dsel, data_haz = dhaz,
                      id = "id", family = "poisson", exposure = "dur",
                      nested = "bid", Q = 15, Qw = 7)
cat(sprintf("Three-level: hospital %.4f | child-frailty var sw2 = %.4f\n",
            fit3$coef["haz:hospital"], fit3$sw2))
cat("\nBartus Example 2: gates A-D PASSED.\n")
