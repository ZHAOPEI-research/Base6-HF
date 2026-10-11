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
# 2. Data preparation
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
# 3. Generate frozen-model LP
# ============================================================

get_lp <- function(d) {

  as.numeric(
    beta_shrunk_final["age80"]      * d$age80 +
      beta_shrunk_final["bmxbmi"]     * d$bmxbmi +
      beta_shrunk_final["lbxrdw"]     * d$lbxrdw +
      beta_shrunk_final["nwr_ns1"]    * d$nwr_ns1 +
      beta_shrunk_final["nwr_ns2"]    * d$nwr_ns2 +
      beta_shrunk_final["nwr_ns3"]    * d$nwr_ns3 +
      beta_shrunk_final["egfr_2021"]  * d$egfr_2021 +
      beta_shrunk_final["med_cnt"]    * d$med_cnt
  )
}

# ============================================================
# 4. Harrell C
# ============================================================

get_cindex <- function(time, status, lp) {

  fit <- survival::concordance(
    Surv(time, status) ~ lp,
    reverse = TRUE
  )

  as.numeric(fit$concordance)
}

# ============================================================
# 5. Prepare all 30 imputed datasets
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

# check identical ordering across imputations
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
# 6. Apparent frozen-model C-index
# ============================================================

C_MI <- sapply(
  completed_list,
  function(d) {

    get_cindex(
      d$followup_month,
      d$death,
      d$Base6_LP
    )
  }
)

C_apparent_frozen <- mean(C_MI)

cat("\n============================================\n")
cat("Frozen Base6 apparent C\n")
cat("============================================\n")

cat(
  "Mean across 30 MI =",
  C_apparent_frozen,
  "\n"
)

cat(
  "MI range =",
  min(C_MI),
  "to",
  max(C_MI),
  "\n"
)

# ============================================================
# 7. Subject-level bootstrap
#
# Same participant indices are resampled in all 30 imputations
# ============================================================

B <- 2000

set.seed(20261003)

n <- nrow(completed_list[[1]])

boot_C <- rep(NA_real_, B)

for (b in seq_len(B)) {

  ix <- sample.int(
    n,
    size = n,
    replace = TRUE
  )

  C_b_MI <- sapply(
    completed_list,
    function(d) {

      db <- d[
        ix,
        ,
        drop = FALSE
      ]

      get_cindex(
        db$followup_month,
        db$death,
        db$Base6_LP
      )
    }
  )

  boot_C[b] <- mean(
    C_b_MI,
    na.rm = TRUE
  )

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
# 8. Percentile 95% CI
# ============================================================

C_CI <- quantile(
  boot_C,
  probs = c(
    0.025,
    0.975
  ),
  na.rm = TRUE
)

C_summary <- data.frame(
  metric = "Frozen Base6 apparent Harrell C",
  estimate = C_apparent_frozen,
  CI_low = unname(C_CI[1]),
  CI_high = unname(C_CI[2]),
  bootstrap_replicates = B
)

cat("\n============================================\n")
cat("Harrell C bootstrap CI\n")
cat("============================================\n")

print(C_summary)

# ============================================================
# 9. Save separately
# ============================================================

write.csv(
  C_summary,
  "Base6_development_HarrellC_95CI.csv",
  row.names = FALSE
)

save(
  C_MI,
  C_apparent_frozen,
  boot_C,
  C_summary,
  file = "NHANES_Base6_HarrellC_95CI.RData"
)

cat("\nSTEP 49 completed\n")
