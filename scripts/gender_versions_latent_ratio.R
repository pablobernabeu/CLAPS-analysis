#!/usr/bin/env Rscript
# scripts/gender_versions_latent_ratio.R
# ---------------------------------------------------------------------------
# The headline numbers in reports/gender_versions_note.md: how much more precise
# the focal estimates would be if every participant rated both agent-gender
# versions of each item, measured on the latent scale of the preregistered
# cumulative-logit model, and what that precision is worth in power at N = 80.
#
# WHY A LATENT-SCALE ESTIMATE, AND WHY THIS ONE
#   scripts/assess_gender_versions.R measures the loss with linear mixed models on
#   the 1-7 ratings. Near the ceiling (82% of English and 88% of Turkish actives
#   are 7s) the rating scale compresses exactly the per-rating noise a second
#   version averages away, so rating-scale ratios understate the gain the ordinal
#   model sees. Its latent-scale projection, which sets per-rating noise to
#   pi^2/3, errs the other way by ignoring discretisation. This script estimates
#   the components directly from the pilot on the latent scale instead.
#
# THE METHOD
#   1. A penalised cumulative-logit model with one effect per verb x sentence type
#      and one per participant x sentence type (ridge SDs 3 and 5, so the fit is
#      close to maximum likelihood but always defined) gives, per verb, the latent
#      Active - Passive contrast (H1b) and the latent Passive mean (H1a).
#   2. The same model is fitted to four halves of the data per random split:
#      two participant halves that keep both versions (A, B) and two complementary
#      one-version datasets built as the recommended lists would build them
#      (H, Hc: verbs split 50/50 into woman-agent and man-agent sets, participants
#      split into two lists). Participant x verb effects cancel inside a
#      within-verb contrast, so the verb-level residual variance of
#        (A - B) / 2 estimates e + c        (both versions, half the participants)
#        (H - Hc) / 2 estimates e - c       (same participants, other version)
#      where e is the one-version noise variance of a verb contrast and c the
#      covariance of the two versions' noise within a participant. The full-data
#      residual variance minus (e + c) / 2 is the between-verb variance tau^2.
#   3. At N = 80, with Sxx the affectedness sum of squares,
#        Var(one version)  = (tau^2 + e80) / Sxx + participant-slope term
#        Var(two versions) = (tau^2 + (e80 + c80) / 2) / Sxx + participant-slope term
#      The participant-slope term is collinear with affectedness, so it cannot be
#      read off verb-level residuals. It is taken from the pilot brms model
#      (outputs/pilot_models/pilot_dgp_v2_pilot_*.rds) after rescaling by the ratio
#      of this model's threshold spacing to the brms thresholds.
#      r = SE(two versions) / SE(one version) is the square root of the ratio.
#   4. Uncertainty: a delete-one-verb jackknife on r. It leaves out participant
#      sampling, so the intervals are too narrow and are labelled as such.
#   5. A check against the generating model: the same pipeline applied to data
#      simulated from the pilot DGP with two INDEPENDENT ratings per cell.
#
# THE POWER TRANSLATION
#   A directional Savage-Dickey BF under a symmetric prior clears 10 when the
#   posterior sign probability reaches 10/11, i.e. |estimate| / SE >= qnorm(10/11)
#   under a normal posterior (R/05_hypothesis_tests.R). Assurance power at an
#   effective SE s is the mean over pilot posterior draws theta of
#   Phi(sign * theta / s - z*). s is solved so that this reproduces the brms
#   decision arm at N = 80 (outputs/design_summary_pilot/joint_power_pilot.csv),
#   and the two-version power is the same average at r * s. This is the reduction
#   scripts/gender_versions_power.R documents in full. The one-version design needs
#   N' = 80 * Q1 / Q2 participants to match two-version precision at N = 80, where
#   Q1 and Q2 are the per-participant variance terms. The verb floor tau^2 is the
#   same in both designs and so drops out, and N' can never exceed 160.
#
# OUTPUT
#   outputs/design_summary_pilot/gender_versions_latent_ratio.csv
#   outputs/design_summary_pilot/gender_versions_power_at_latent_r.csv
#
# USAGE (from design_analysis/; about 7 minutes on 3 workers)
#   Rscript scripts/gender_versions_latent_ratio.R [--splits 20] [--dgp_sims 6]
#
# Participant-level data stay in memory. Nothing but aggregate tables is written.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(parallel)
})

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--splits",     default = 20L, type = "integer",
                        help = "random list splits per language on the real data"),
  optparse::make_option("--dgp_sims",   default = 6L,  type = "integer",
                        help = "datasets simulated from the pilot DGP per language"),
  optparse::make_option("--dgp_splits", default = 4L,  type = "integer",
                        help = "list splits per simulated dataset"),
  optparse::make_option("--workers",    default = 3L,  type = "integer"),
  optparse::make_option("--outdir",     default = "outputs/design_summary_pilot")
)))

LANGS <- c("English", "Turkish", "Norwegian")
SEED  <- 20260916L

# ---------------------------------------------------------------------------
# Everything a worker needs, defined as one function so that it can be shipped
# to the cluster whole. Each language sets its own seed before drawing anything,
# which makes the result independent of how languages are spread over workers.
# ---------------------------------------------------------------------------
run_language <- function(L, splits, dgp_sims, dgp_splits, seed) {
  suppressPackageStartupMessages(library(Matrix))

  X <- read.csv("data/pilot/claps_pilot_harmonised.csv")
  d <- X[X$Language == L, c("Participant", "Verb_ID", "S_Type", "Item", "Response", "trial")]
  dgp <- readRDS(sprintf("outputs/pilot_models/pilot_dgp_v2_pilot_%s.rds", L))
  d$sem <- unname(dgp$verb_affectedness[d$Verb_ID])
  stopifnot(!anyNA(d$sem))
  d$Version <- sub(".*_", "", d$Item)
  # Norwegian watch x Active carries two rows labelled woman-agent in every pilot
  # cell; the raw sentence text shows one of them is the man-agent sentence. Any
  # duplicate label becomes the other version, so every cell has one of each.
  key <- paste(d$Participant, d$Verb_ID, d$S_Type)
  d <- d[order(key, d$trial), ]
  key <- paste(d$Participant, d$Verb_ID, d$S_Type)
  dup <- duplicated(paste(key, d$Version))
  d$Version[dup] <- ifelse(d$Version[dup] == "Woman", "Man", "Woman")
  stopifnot(all(table(paste(key, d$Version)) == 1))
  d$S_Type <- as.character(d$S_Type)

  verbs   <- names(dgp$verb_affectedness)
  s_types <- dgp$s_types

  # Penalised cumulative logit by Fisher scoring on a sparse design. Returns the
  # verb x sentence-type effects, with the fitted thresholds attached so the
  # latent scale can be compared with the brms model's.
  fit_vs <- function(dd, sd_vs = 3, sd_ps = 5, maxit = 100, tol = 1e-7) {
    yv <- as.integer(factor(dd$Response))
    K  <- max(yv)
    vs <- factor(paste(dd$Verb_ID, dd$S_Type, sep = "|"),
                 levels = as.vector(outer(verbs, s_types, paste, sep = "|")))
    ps <- factor(paste(dd$Participant, dd$S_Type, sep = "|"))
    n  <- nrow(dd)
    Xv <- sparseMatrix(i = seq_len(n), j = as.integer(vs), x = 1, dims = c(n, nlevels(vs)))
    Xp <- sparseMatrix(i = seq_len(n), j = as.integer(ps), x = 1, dims = c(n, nlevels(ps)))
    Xm <- cbind(Xv, Xp)
    p  <- ncol(Xm)
    lam <- c(rep(1 / sd_vs^2, ncol(Xv)), rep(1 / sd_ps^2, ncol(Xp)))
    cp <- cumsum(tabulate(yv, K)) / n
    theta <- qlogis(cp[-K]); beta <- rep(0, p)
    loglik <- function(theta, beta) {
      eta <- as.vector(Xm %*% beta)
      Fm <- cbind(0, plogis(outer(-eta, theta, "+")), 1)
      pr <- Fm[cbind(seq_len(n), yv + 1)] - Fm[cbind(seq_len(n), yv)]
      sum(log(pmax(pr, 1e-300))) - 0.5 * sum(lam * beta^2)
    }
    ll <- loglik(theta, beta)
    for (it in seq_len(maxit)) {
      eta <- as.vector(Xm %*% beta)
      A   <- outer(-eta, theta, "+")
      Fm  <- cbind(0, plogis(A), 1)
      fm  <- cbind(0, dlogis(A), 0)
      P   <- pmax(Fm[, -1] - Fm[, -(K + 1)], 1e-12)
      dfm <- fm[, -1] - fm[, -(K + 1)]
      py  <- P[cbind(seq_len(n), yv)]
      u_eta <- -dfm[cbind(seq_len(n), yv)] / py
      W <- rowSums(dfm^2 / P)
      s_theta <- numeric(K - 1)
      for (j in seq_len(K - 1)) {
        s_theta[j] <- sum(fm[, j + 1] * ((yv == j) - (yv == j + 1)) / py)
      }
      Itt <- matrix(0, K - 1, K - 1)
      C   <- matrix(0, n, K - 1)
      for (j in seq_len(K - 1)) {
        fj <- fm[, j + 1]
        Itt[j, j] <- sum(fj^2 * (1 / P[, j] + 1 / P[, j + 1]))
        if (j < K - 1) {
          v <- -sum(fj * fm[, j + 2] / P[, j + 1])
          Itt[j, j + 1] <- v; Itt[j + 1, j] <- v
        }
        C[, j] <- -fj * dfm[, j] / P[, j] + fj * dfm[, j + 1] / P[, j + 1]
      }
      Itb <- as.matrix(crossprod(Xm, C))
      Ibb <- as.matrix(crossprod(Xm, Xm * W)) + diag(lam)
      s_beta <- as.vector(crossprod(Xm, u_eta)) - lam * beta
      step <- solve(rbind(cbind(Itt, t(Itb)), cbind(Itb, Ibb)), c(s_theta, s_beta))
      a <- 1
      # Step halving keeps the thresholds ordered and the penalised likelihood
      # from falling, which plain Fisher scoring does not guarantee near the ceiling.
      repeat {
        th_new <- theta + a * step[1:(K - 1)]
        b_new  <- beta + a * step[-(1:(K - 1))]
        if (all(diff(th_new) > 0)) {
          ll_new <- loglik(th_new, b_new)
          if (is.finite(ll_new) && ll_new >= ll - 1e-8) break
        }
        a <- a / 2
        if (a < 1e-8) break
      }
      theta <- th_new; beta <- b_new; ll_old <- ll; ll <- ll_new
      if (max(abs(a * step)) < tol || abs(ll - ll_old) < 1e-9) break
    }
    mat <- matrix(beta[seq_len(nlevels(vs))], nrow = length(verbs),
                  dimnames = list(verbs, s_types))
    attr(mat, "theta") <- theta
    mat
  }

  # The recommended lists: verbs halved into two sets, participants into two
  # lists, the version crossed so each verb x version is seen by half the sample.
  scheme_A <- function(dd, parts) {
    vset  <- sample(verbs, length(verbs) %/% 2)
    list1 <- sample(parts, length(parts) %/% 2)
    want  <- ifelse((dd$Participant %in% list1) == (dd$Verb_ID %in% vset), "Woman", "Man")
    list(H = dd[dd$Version == want, ], Hc = dd[dd$Version != want, ])
  }

  run_replicates <- function(dd, M) {
    parts <- unique(dd$Participant)
    full  <- fit_vs(dd)
    reps  <- vector("list", M)
    for (m in seq_len(M)) {
      A  <- sample(parts, length(parts) %/% 2)
      sA <- scheme_A(dd, parts)
      reps[[m]] <- list(A  = fit_vs(dd[dd$Participant %in% A, ]),
                        B  = fit_vs(dd[!dd$Participant %in% A, ]),
                        H  = fit_vs(sA$H),
                        Hc = fit_vs(sA$Hc))
    }
    list(full = full, reps = reps, N = length(parts))
  }

  # Two conditionally independent ratings per cell from the pilot DGP: the case a
  # second version helps most, given the fitted variance components.
  simulate_two_versions <- function(N) {
    grid <- expand.grid(pi = seq_len(N), Verb = verbs, S_Type = s_types,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    sem  <- unname(dgp$verb_affectedness[grid$Verb])
    isA  <- grid$S_Type == "Active"; isPP <- grid$S_Type == "Pseudo_Passive"
    tc <- list(Intercept = rep(1, nrow(grid)), S_TypeActive = as.numeric(isA),
               S_TypePseudo_Passive = as.numeric(isPP), Semantics_scaled = sem,
               `S_TypeActive:Semantics_scaled` = isA * sem,
               `S_TypePseudo_Passive:Semantics_scaled` = isPP * sem)
    ff  <- dgp$fixef
    eta <- Reduce(`+`, lapply(names(ff), function(t) ff[[t]] * tc[[t]]))
    bP <- MASS::mvrnorm(N, rep(0, ncol(dgp$Sigma_part)), dgp$Sigma_part)
    colnames(bP) <- colnames(dgp$Sigma_part)
    eta <- eta + Reduce(`+`, lapply(colnames(bP), function(t) bP[grid$pi, t] * tc[[t]]))
    bV <- MASS::mvrnorm(length(verbs), rep(0, ncol(dgp$Sigma_verb)), dgp$Sigma_verb)
    colnames(bV) <- colnames(dgp$Sigma_verb)
    vi  <- match(grid$Verb, verbs)
    eta <- eta + Reduce(`+`, lapply(colnames(bV), function(t) bV[vi, t] * tc[[t]]))
    draw <- function() {
      Fm <- plogis(outer(-eta, dgp$thresholds, "+"))
      1L + rowSums(runif(length(eta)) > Fm)
    }
    base <- data.frame(Participant = paste0("P", grid$pi), Verb_ID = grid$Verb,
                       S_Type = grid$S_Type, sem = sem, stringsAsFactors = FALSE)
    rbind(cbind(base, Response = draw(), Version = "Woman"),
          cbind(base, Response = draw(), Version = "Man"))
  }

  set.seed(seed + match(L, c("English", "Turkish", "Norwegian")))
  real <- run_replicates(d, splits)
  sims <- lapply(seq_len(dgp_sims), function(i) run_replicates(simulate_two_versions(real$N), dgp_splits))
  list(language = L, real = real, sims = sims, sem = dgp$verb_affectedness, dgp = dgp)
}

# ---------------------------------------------------------------------------
# Variance components and r from the fitted verb contrasts
# ---------------------------------------------------------------------------
resid_var <- function(y, sem) {
  f <- lm.fit(cbind(1, sem), y)
  sum(f$residuals^2) / (length(sem) - 2)
}
contrast <- function(m, h) if (h == "h1b") m[, "Active"] - m[, "Passive"] else m[, "Passive"]

components <- function(obj, sem, h, sig_part, kscale = 1, keep = seq_along(sem), N80 = 80) {
  sem <- sem[keep]; Np <- obj$N
  g <- function(m) contrast(m, h)[keep]
  v_full <- resid_var(g(obj$full), sem)
  a  <- mean(sapply(obj$reps, function(r) resid_var(g(r$A) - g(r$B), sem) / 2))
  b  <- mean(sapply(obj$reps, function(r) resid_var(g(r$H) - g(r$Hc), sem) / 2))
  eN <- (a + b) / 2; cN <- (a - b) / 2
  tau2 <- v_full - a / 2
  Sxx  <- sum((sem - mean(sem))^2)
  e80  <- eN * Np / N80; c80 <- cN * Np / N80
  ps   <- sig_part * kscale^2 / N80
  v1   <- (tau2 + e80) / Sxx
  v2   <- (tau2 + (e80 + c80) / 2) / Sxx
  v2i  <- (tau2 + e80 / 2) / Sxx
  c(e = eN, c = cN, rho = cN / eN, tau2 = tau2, sxx = Sxx,
    r80 = sqrt((v2 + ps) / (v1 + ps)), r80_if_independent = sqrt((v2i + ps) / (v1 + ps)),
    share_verbs = (tau2 / Sxx) / (v1 + ps), share_participant_slopes = ps / (v1 + ps),
    share_rating_noise = (e80 / Sxx) / (v1 + ps), se80_one_version = sqrt(v1 + ps),
    # Per-participant variance terms, for the participant-equivalence calculation.
    q1 = sig_part * kscale^2 + eN * Np / Sxx, q2 = sig_part * kscale^2 + (eN + cN) / 2 * Np / Sxx)
}

message("[latent-r] ", length(LANGS), " languages, ", opt$splits, " splits, ",
        opt$dgp_sims, " DGP datasets x ", opt$dgp_splits, " splits, ", opt$workers, " workers")
t0 <- proc.time()[["elapsed"]]
cl <- parallel::makeCluster(min(opt$workers, length(LANGS)))
invisible(parallel::clusterCall(cl, setwd, getwd()))
fits <- parallel::parLapply(cl, LANGS, run_language, splits = opt$splits,
                            dgp_sims = opt$dgp_sims, dgp_splits = opt$dgp_splits, seed = SEED)
parallel::stopCluster(cl)
message("[latent-r] fits done in ", round(proc.time()[["elapsed"]] - t0), " s")

rows <- list()
for (R in fits) {
  L <- R$language; sem <- R$sem; th_brms <- R$dgp$thresholds
  # Ratio of this model's threshold spacing to the brms model's: the factor that
  # puts the brms participant-slope variance on this model's latent scale.
  kfun <- function(th) unname(coef(lm(th ~ th_brms))[2])
  k_real <- kfun(attr(R$real$full, "theta"))
  for (h in c("h1b", "h1a")) {
    term <- if (h == "h1b") "S_TypeActive:Semantics_scaled" else "Semantics_scaled"
    sig_part <- R$dgp$Sigma_part[term, term]
    est <- components(R$real, sem, h, sig_part, k_real)
    V <- length(sem)
    jk <- sapply(seq_len(V), function(i) components(R$real, sem, h, sig_part, k_real,
                                                     keep = setdiff(seq_len(V), i))["r80"])
    jk_se <- sqrt((V - 1) / V * sum((jk - mean(jk))^2))
    rows[[length(rows) + 1]] <- data.frame(language = L, hypothesis = h, source = "pilot",
      threshold_scale = k_real, t(est), r80_jackknife_se = jk_se, check.names = FALSE)
    if (length(R$sims)) {
      ss <- t(sapply(R$sims, function(o) components(o, sem, h, sig_part, kfun(attr(o$full, "theta")))))
      rows[[length(rows) + 1]] <- data.frame(language = L, hypothesis = h,
        source = sprintf("dgp_independent_mean_of_%d", nrow(ss)),
        threshold_scale = mean(sapply(R$sims, function(o) kfun(attr(o$full, "theta")))),
        t(colMeans(ss)), r80_jackknife_se = sd(ss[, "r80"]) / sqrt(nrow(ss)), check.names = FALSE)
    }
  }
}
latent <- do.call(rbind, rows)
latent$r80_lo <- latent$r80 - 2 * latent$r80_jackknife_se
latent$r80_hi <- pmin(1, latent$r80 + 2 * latent$r80_jackknife_se)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(latent, file.path(opt$outdir, "gender_versions_latent_ratio.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Power at the pilot estimate of r
# ---------------------------------------------------------------------------
jp <- read.csv(file.path(opt$outdir, "joint_power_pilot.csv"))
jp <- jp[jp$mode == "assurance", ]
ceil <- read.csv(file.path(opt$outdir, "pilot_params_ceilings.csv"))
zs <- qnorm(10 / 11)
power_at <- function(th, s, sign) mean(pnorm(sign * th / s - zs))
solve_s  <- function(th, p, sign) uniroot(function(s) power_at(th, s, sign) - p, c(0.01, 10), tol = 1e-10)$root

prow <- list()
for (R in fits) {
  L <- R$language
  arm80 <- jp[jp$language == L & jp$n_participants == 80, ]
  big   <- jp[jp$language == L & jp$n_participants >= 100, ]
  draws <- list(); s80 <- list(); rr <- list()
  for (h in c("h1b", "h1a")) {
    term <- if (h == "h1b") "S_TypeActive:Semantics_scaled" else "Semantics_scaled"
    sign <- if (h == "h1b") -1 else 1
    th <- as.numeric(R$dgp$fixef_draws[, term])
    p80 <- arm80[[paste0("p_", h)]]
    p_plateau <- sum(big[[paste0("p_", h)]] * big$reps) / sum(big$reps)
    lr <- latent[latent$language == L & latent$hypothesis == h & latent$source == "pilot", ]
    s  <- solve_s(th, p80, sign); spl <- solve_s(th, p_plateau, sign)
    floor_se <- ceil[ceil$language == L, if (h == "h1b") "se_floor_h1b" else "se_floor_h1a"]
    draws[[h]] <- sign * th; s80[[h]] <- s; rr[[h]] <- lr
    prow[[length(prow) + 1]] <- data.frame(language = L, hypothesis = h,
      arm_power_n80 = p80, arm_reps = arm80$reps, arm_mcse = sqrt(p80 * (1 - p80) / arm80$reps),
      s80 = s, r = lr$r80, r_lo = lr$r80_lo, r_hi = lr$r80_hi,
      power_both_versions = power_at(th, lr$r80 * s, sign),
      power_both_lo = power_at(th, lr$r80_hi * s, sign),
      power_both_hi = power_at(th, lr$r80_lo * s, sign),
      arm_power_n100plus = p_plateau, arm_reps_n100plus = sum(big$reps),
      power_both_from_n100plus = power_at(th, lr$r80 * spl, sign),
      power_limit_at_verb_floor = power_at(th, floor_se, sign),
      one_version_n_matching_two_at_80 = 80 * lr$q1 / lr$q2)
  }
  # Joint H1a and H1b, treating the two tests as independent given the generating
  # draw. The arm's own joint figure is reported beside the approximation so the
  # independence assumption can be judged; the gain is the approximation's.
  joint <- function(ra, rb) mean(pnorm(draws$h1a / (ra * s80$h1a) - zs) * pnorm(draws$h1b / (rb * s80$h1b) - zs))
  j1 <- joint(1, 1)
  prow[[length(prow) + 1]] <- data.frame(language = L, hypothesis = "joint",
    arm_power_n80 = arm80$p_joint, arm_reps = arm80$reps,
    arm_mcse = sqrt(arm80$p_joint * (1 - arm80$p_joint) / arm80$reps),
    s80 = NA, r = NA, r_lo = NA, r_hi = NA,
    power_both_versions = arm80$p_joint + joint(rr$h1a$r80, rr$h1b$r80) - j1,
    power_both_lo = arm80$p_joint + joint(rr$h1a$r80_hi, rr$h1b$r80_hi) - j1,
    power_both_hi = arm80$p_joint + joint(rr$h1a$r80_lo, rr$h1b$r80_lo) - j1,
    arm_power_n100plus = NA, arm_reps_n100plus = NA, power_both_from_n100plus = NA,
    power_limit_at_verb_floor = NA, one_version_n_matching_two_at_80 = NA)
  prow[[length(prow)]]$joint_approximation_one_version <- j1
}
power <- do.call(rbind, lapply(prow, function(x) {
  if (is.null(x$joint_approximation_one_version)) x$joint_approximation_one_version <- NA
  x
}))
power$gain_points <- 100 * (power$power_both_versions - power$arm_power_n80)
write.csv(power, file.path(opt$outdir, "gender_versions_power_at_latent_r.csv"), row.names = FALSE)

cat("\nLatent r at N = 80 (pilot rows):\n")
print(format(latent[latent$source == "pilot", c("language", "hypothesis", "rho", "r80", "r80_lo", "r80_hi",
                                                 "share_verbs", "share_participant_slopes", "share_rating_noise")],
             digits = 3), row.names = FALSE)
cat("\nPower at BF 10, N = 80:\n")
print(format(power[, c("language", "hypothesis", "arm_power_n80", "arm_mcse", "power_both_versions",
                       "gain_points", "power_both_lo", "power_both_hi", "one_version_n_matching_two_at_80")],
             digits = 3), row.names = FALSE)
cat("\nR:", R.version.string, "| elapsed", round(proc.time()[["elapsed"]] - t0), "s\n")
