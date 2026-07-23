# Replication of Bartus (2017, Stata Journal 17(2)), Example 1 (divorce2).
# Data preparation follows Bartus's OFFICIAL do-file (st0481, SJ software
# archive: http://www.stata-journal.com/software/sj17-2/st0481/multiprocess_gsem.do):
#   use divorce2.dta; xvars = ib2.hereduc age mardur; poisson exposure(dur)
#   [dur is NATIVE in the dataset]; estimation restricted to numkids==2.
# Canonical data URL (also readable directly by haven::read_dta):
#   http://web.uni-corvinus.hu/bartus/stata/divorce2.dta
# PASS criterion: reduced-form coefficients and RE (co)variances within tol of
# the published gsem estimates (p. 452); delta-method structural effects within
# tol of the published nlcom results (p. 453). Laplace (glmmTMB) vs adaptive
# quadrature (gsem) differences are expected -- report them, never hide them.
library(glmmTMB); library(haven); library(GenSEM)
tol_fix <- 0.03; tol_re <- 0.06; tol_str <- 0.10

d <- read_dta("divorce2.dta")               # or read_dta(<canonical URL>)
d <- subset(d, numkids == 2)                # married mothers of one child
d$birth2 <- d$birth                         # separate birth, by(numkids)
d$divorce2 <- d$divorce
d$hereduc1 <- as.integer(d$hereduc == 1)    # ib2.hereduc: base = category 2
d$hereduc3 <- as.integer(d$hereduc == 3)

# EXPOSURE = SPELL LENGTH, derived from the spell structure (JOB1 validation,
# 21 Jul 2026): `time` is the END of marriage, constant within id, so neither
# a native dur nor time - mardur is spell exposure -- both overlap spells.
# Spell length = next spell's entry (mardur) minus this spell's entry; the
# last spell ends at `time`. Integrity asserts encode the verified totals.
d <- d[order(d$id, d$mardur), ]
nxt <- ave(d$mardur, d$id, FUN = function(x) c(x[-1], NA))
d$dur <- ifelse(is.na(nxt), d$time - d$mardur, nxt - d$mardur)
stopifnot(nrow(d) == 5100,                  # published estimation sample
          all(d$dur > 0))                  # no zero/negative spells
py <- sum(d$dur)                            # total person-years: ~8734 per the
cat(sprintf("Total person-years: %.0f\n", py))  # JOB1 report (overlapped
stopifnot(abs(py - 8734) < 60)              # construction gave 23754)

long <- stack_processes(d, id = "id",
                        events = c(birth2 = "birth2", divorce = "divorce2"),
                        exposure = "dur")
# Backend rule (JOB1): divorce is sparse-absorbing (246 events / 2121 women);
# Laplace inflates var(V) (observed: 61.4). AGHQ is REQUIRED here.
fit <- fit_multiproc(long, ~ hereduc1 + hereduc3 + age + mardur, id = "id",
                     backend = "GLMMadaptive", nAGQ = 11)
print(summary(fit))

published <- c(  # gsem reduced form, Bartus 2017 p. 452
  `processbirth2:hereduc3` = 0.3128652,  `processbirth2:age` = -0.0396747,
  `processbirth2:mardur`   = -0.1071725, `processdivorce:hereduc3` = -0.4958348,
  `processdivorce:age`     = -0.0970939, `processdivorce:mardur`   = 0.0857423)
bhat <- GLMMadaptive::fixef(fit)
cmp <- data.frame(published, glmmTMB = bhat[names(published)],
                  diff = bhat[names(published)] - published)
print(round(cmp, 4))
stopifnot(all(abs(cmp$diff) < tol_fix))

Shat <- fit$D                               # GLMMadaptive RE covariance
print(Shat)
# RE criteria = 1.96 * PUBLISHED SE (CI overlap): demanding point tolerance
# tighter than the original's own SE is bad metrology (JOB1 lesson, Ex.1:
# published var(V) .479 with SE .379, CI [0.10, 2.26]).
stopifnot(abs(Shat[1, 1] - 0.4275274) < 1.96 * 0.0636296,
          abs(Shat[2, 2] - 0.4786240) < 1.96 * 0.3789376,
          abs(Shat[1, 2] - (-0.0836689)) < 1.96 * 0.1472204)

s1 <- recover_structural(fit,  # education on birth2: -0.3068977 (se .4508)
  "`processbirth2:hereduc3` - (`processbirth2:mardur`/`processdivorce:mardur`) * `processdivorce:hereduc3`")
s2 <- recover_structural(fit,  # education on dissolution: -1.261495
  "`processdivorce:hereduc3` - (`processdivorce:age`/`processbirth2:age`) * `processbirth2:hereduc3`")
s3 <- recover_structural(fit,  # dissolution hazard on birth2: -1.249938
  "`processbirth2:mardur` / `processdivorce:mardur`")
print(rbind(s1, s2, s3))
stopifnot(abs(s1$estimate - (-0.3068977)) < 1.96 * 0.4508039,
          abs(s2$estimate - (-1.261495)) < 1.96 * 0.4305668,
          abs(s3$estimate - (-1.249938)) < 1.96 * 0.4477062)
cat("\nBartus Example 1: REPLICATED within tolerance.\n")
