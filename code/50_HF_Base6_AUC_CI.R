rm(list = ls())
gc()

library(mice)
library(survival)
library(splines)
library(dplyr)
library(timeROC)

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
# 3. Frozen Base6 linear predictor
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
# 4. AUC function
# ============================================================

get_auc <- function(d, t0) {

  z <- tryCatch(

    timeROC(
      T = d$followup_month,
      delta = d$death,
      marker = d$Base6_LP,
      cause = 1,
      weighting = "marginal",
      times = t0,
      iid = FALSE
    ),

    error = function(e) NULL
  )

  if (is.null(z)) {
    return(NA_real_)
  }

  # timeROC automatically adds t = 0;
  # use the AUC corresponding to the requested horizon
  auc <- as.numeric(
    tail(z$AUC, 1)
  )

  if (!is.finite(auc)) {
    return(NA_real_)
  }

  auc
}

# ============================================================
# 5. Prepare all 30 imputations
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
# 7. Frozen-model apparent AUC
# ============================================================

apparent_AUC <- sapply(
  times,
  function(t0) {

    auc_mi <- sapply(
      completed_list,
      function(d) {
        get_auc(
          d,
          t0
        )
      }
    )

    mean(
      auc_mi,
      na.rm = TRUE
    )
  }
)

cat("\n============================================\n")
cat("Frozen Base6 apparent AUC\n")
cat("============================================\n")

print(apparent_AUC)

# ============================================================
# 8. Subject-level bootstrap
# ============================================================

B <- 2000

set.seed(20261003)

n <- nrow(completed_list[[1]])

boot_AUC <- matrix(
  NA_real_,
  nrow = B,
  ncol = length(times)
)

colnames(boot_AUC) <- names(times)

for (b in seq_len(B)) {

  ix <- sample.int(
    n,
    size = n,
    replace = TRUE
  )

  for (j in seq_along(times)) {

    t0 <- times[j]

    auc_mi_b <- sapply(
      completed_list,
      function(d) {

        db <- d[
          ix,
          ,
          drop = FALSE
        ]

        get_auc(
          db,
          t0
        )
      }
    )

    boot_AUC[b, j] <-
      mean(
        auc_mi_b,
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
# 9. 95% percentile CI
# ============================================================

AUC_summary <- data.frame(
  horizon = names(times),

  frozen_apparent =
    as.numeric(
      apparent_AUC
    ),

  CI_low =
    apply(
      boot_AUC,
      2,
      quantile,
      probs = 0.025,
      na.rm = TRUE
    ),

  CI_high =
    apply(
      boot_AUC,
      2,
      quantile,
      probs = 0.975,
      na.rm = TRUE
    ),

  bootstrap_replicates = B
)

cat("\n============================================\n")
cat("AUC bootstrap 95% CI\n")
cat("============================================\n")

print(AUC_summary)

# ============================================================
# 10. Save
# ============================================================

write.csv(
  AUC_summary,
  "Base6_development_AUC_95CI.csv",
  row.names = FALSE
)

save(
  apparent_AUC,
  boot_AUC,
  AUC_summary,
  file = "NHANES_Base6_AUC_95CI.RData"
)

cat("\nSTEP 50 completed\n")
