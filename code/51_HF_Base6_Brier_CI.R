rm(list = ls())
gc()

library(mice)
library(survival)
library(splines)
library(dplyr)

analysis_dir <- Sys.getenv(
  "BASE6_ANALYSIS_DIR",
  unset = getwd()
)
setwd(analysis_dir)

# ============================================================
# 1. Load frozen model and MI object
# ============================================================

load("NHANES_2011_2018_MICE_FINAL.RData")
load("NHANES_Base6_FINAL_ABSOLUTE_RISK_FIXED.RData")

# ============================================================
# 2. Prepare data
# ============================================================

prepare_data <- function(d) {

  d <- d %>%
    mutate(
      age80 = pmin(ridageyr, 80),

      nwr = nwr_model,

      kappa =
        ifelse(
          female == 1,
          0.7,
          0.9
        ),

      alpha =
        ifelse(
          female == 1,
          -0.241,
          -0.302
        ),

      egfr_2021 =
        142 *
        pmin(scr_mgdl / kappa, 1)^alpha *
        pmax(scr_mgdl / kappa, 1)^(-1.200) *
        (0.9938^ridageyr) *
        ifelse(
          female == 1,
          1.012,
          1
        )
    )

  nb <- splines::ns(
    d$nwr,
    knots = nwr_knots_final,
    Boundary.knots = nwr_boundary_final,
    intercept = FALSE
  )

  d$nwr_ns1 <- nb[, 1]
  d$nwr_ns2 <- nb[, 2]
  d$nwr_ns3 <- nb[, 3]

  d
}

# ============================================================
# 3. Frozen Base6 LP
# ============================================================

get_lp <- function(d) {

  as.numeric(
    beta_shrunk_final["age80"]     * d$age80 +
      beta_shrunk_final["bmxbmi"]    * d$bmxbmi +
      beta_shrunk_final["lbxrdw"]    * d$lbxrdw +
      beta_shrunk_final["nwr_ns1"]   * d$nwr_ns1 +
      beta_shrunk_final["nwr_ns2"]   * d$nwr_ns2 +
      beta_shrunk_final["nwr_ns3"]   * d$nwr_ns3 +
      beta_shrunk_final["egfr_2021"] * d$egfr_2021 +
      beta_shrunk_final["med_cnt"]   * d$med_cnt
  )
}

# ============================================================
# 4. IPCW Brier score
# ============================================================

ipcw_brier <- function(
    time,
    status,
    pred_risk,
    t
) {

  censor_status <- 1 - status

  km_cens <- survfit(
    Surv(time, censor_status) ~ 1
  )

  G_at <- function(s) {

    out <- summary(
      km_cens,
      times = s,
      extend = TRUE
    )$surv

    as.numeric(out)
  }

  G_t <- G_at(t)

  if (
    length(G_t) == 0 ||
    !is.finite(G_t) ||
    G_t <= 0
  ) {
    return(NA_real_)
  }

  n <- length(time)

  w <- numeric(n)
  y <- numeric(n)

  for (i in seq_len(n)) {

    # Death by t
    if (
      time[i] <= t &&
      status[i] == 1
    ) {

      # G(T_i -)
      s_i <- max(
        time[i] - 1e-8,
        0
      )

      G_i <- G_at(s_i)

      if (
        length(G_i) > 0 &&
        is.finite(G_i) &&
        G_i > 0
      ) {
        w[i] <- 1 / G_i
        y[i] <- 1
      }

      # Event-free beyond t
    } else if (
      time[i] > t
    ) {

      w[i] <- 1 / G_t
      y[i] <- 0

      # Censored before or at t
    } else {

      w[i] <- 0
      y[i] <- 0
    }
  }

  mean(
    w *
      (y - pred_risk)^2,
    na.rm = TRUE
  )
}

# ============================================================
# 5. Prepare 30 MI datasets
# ============================================================

completed_list <- lapply(
  1:imp_final$m,
  function(m) {

    d <- prepare_data(
      complete(
        imp_final,
        m
      )
    )

    d$Base6_LP <- get_lp(d)

    d
  }
)

stopifnot(
  all(
    sapply(
      completed_list,
      nrow
    ) == 739
  )
)

for (m in 2:length(completed_list)) {

  stopifnot(
    identical(
      completed_list[[1]]$seqn,
      completed_list[[m]]$seqn
    )
  )
}

# ============================================================
# 6. Horizons
# ============================================================

times <- c(
  "1-year" = 12,
  "3-year" = 36,
  "5-year" = 60
)

# ============================================================
# 7. Get MI-specific baseline hazards
# ============================================================

get_H0 <- function(m, t0) {

  z <- baseline_MI[
    baseline_MI$imputation == m &
      baseline_MI$time_month == t0,
    ,
    drop = FALSE
  ]

  if (nrow(z) != 1) {
    stop(
      paste(
        "Could not uniquely identify H0 for MI",
        m,
        "time",
        t0
      )
    )
  }

  z$H0_centered
}

# ============================================================
# 8. Frozen apparent Brier
# ============================================================

apparent_Brier_MI <- matrix(
  NA_real_,
  nrow = imp_final$m,
  ncol = length(times)
)

colnames(apparent_Brier_MI) <- names(times)

for (m in seq_len(imp_final$m)) {

  d <- completed_list[[m]]

  for (j in seq_along(times)) {

    t0 <- times[j]

    H0 <- get_H0(
      m,
      t0
    )

    risk <- 1 -
      exp(
        -H0 *
          exp(
            d$Base6_LP -
              LP_center
          )
      )

    apparent_Brier_MI[m, j] <-
      ipcw_brier(
        time = d$followup_month,
        status = d$death,
        pred_risk = risk,
        t = t0
      )
  }
}

apparent_Brier_frozen <-
  apply(
    apparent_Brier_MI,
    2,
    median,
    na.rm = TRUE
  )

cat("\n============================================\n")
cat("Frozen Base6 apparent Brier\n")
cat("============================================\n")

print(apparent_Brier_frozen)

# ============================================================
# 9. Subject-level bootstrap
# ============================================================

B <- 2000

set.seed(20261003)

n <- nrow(completed_list[[1]])

boot_Brier <- matrix(
  NA_real_,
  nrow = B,
  ncol = length(times)
)

colnames(boot_Brier) <- names(times)

for (b in seq_len(B)) {

  ix <- sample.int(
    n,
    size = n,
    replace = TRUE
  )

  for (j in seq_along(times)) {

    t0 <- times[j]

    brier_MI_b <- sapply(
      seq_len(imp_final$m),
      function(m) {

        d <- completed_list[[m]][
          ix,
          ,
          drop = FALSE
        ]

        H0 <- get_H0(
          m,
          t0
        )

        risk <- 1 -
          exp(
            -H0 *
              exp(
                d$Base6_LP -
                  LP_center
              )
          )

        ipcw_brier(
          time = d$followup_month,
          status = d$death,
          pred_risk = risk,
          t = t0
        )
      }
    )

    boot_Brier[b, j] <-
      median(
        brier_MI_b,
        na.rm = TRUE
      )
  }

  if (b %% 100 == 0) {
    cat(
      "Bootstrap",
      b,
      "/",
      B,
      "\n"
    )
  }
}

# ============================================================
# 10. 95% CI
# ============================================================

Brier_summary <- data.frame(
  horizon = names(times),

  frozen_apparent =
    as.numeric(
      apparent_Brier_frozen
    ),

  CI_low =
    apply(
      boot_Brier,
      2,
      quantile,
      probs = 0.025,
      na.rm = TRUE
    ),

  CI_high =
    apply(
      boot_Brier,
      2,
      quantile,
      probs = 0.975,
      na.rm = TRUE
    ),

  bootstrap_replicates = B
)

cat("\n============================================\n")
cat("Brier bootstrap 95% CI\n")
cat("============================================\n")

print(Brier_summary)

# ============================================================
# 11. Sanity check against formal Table 2 values
# ============================================================

formal_apparent <- c(
  "1-year" = 0.06229284,
  "3-year" = 0.14022521,
  "5-year" = 0.16862330
)

comparison <- data.frame(
  horizon = names(times),
  formal_apparent =
    as.numeric(formal_apparent),
  frozen_apparent =
    as.numeric(apparent_Brier_frozen),
  difference =
    as.numeric(apparent_Brier_frozen) -
    as.numeric(formal_apparent)
)

cat("\n============================================\n")
cat("Comparison with formal Table 2 values\n")
cat("============================================\n")

print(comparison)

# ============================================================
# 12. Save
# ============================================================

write.csv(
  Brier_summary,
  "Base6_development_Brier_95CI.csv",
  row.names = FALSE
)

write.csv(
  comparison,
  "Base6_development_Brier_CI_audit.csv",
  row.names = FALSE
)

save(
  apparent_Brier_MI,
  apparent_Brier_frozen,
  boot_Brier,
  Brier_summary,
  comparison,
  file = "NHANES_Base6_Brier_95CI.RData"
)

cat("\nSTEP 51 completed\n")
