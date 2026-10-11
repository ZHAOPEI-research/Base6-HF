# ============================================================
# Base6 1-year mortality risk prediction
#
# Final frozen model
# Predictors:
#   age, BMI, RDW, NWR, eGFR, medication count
#
# NWR = absolute neutrophil count / total white blood cell count
#
# This implementation requires all six predictors.
# No real-time imputation strategy for missing predictors
# has been validated for clinical implementation.
# ============================================================

predict_base6_1year <- function(
    age,
    bmi,
    rdw,
    nwr,
    egfr,
    med_count
) {

  # ----------------------------------------------------------
  # 1. Basic input checks
  # ----------------------------------------------------------

  inputs <- list(
    age = age,
    bmi = bmi,
    rdw = rdw,
    nwr = nwr,
    egfr = egfr,
    med_count = med_count
  )

  input_lengths <- vapply(
    inputs,
    length,
    integer(1)
  )

  if (length(unique(input_lengths)) != 1) {
    stop("All predictor inputs must have the same length.")
  }

  if (any(vapply(inputs, function(x) any(is.na(x)), logical(1)))) {
    stop(
      paste0(
        "Base6 requires all six predictors. ",
        "Predictions should not be generated when required ",
        "predictor values are missing."
      )
    )
  }

  if (any(!is.finite(unlist(inputs)))) {
    stop("All predictor values must be finite numeric values.")
  }


  # ----------------------------------------------------------
  # 2. Age handling
  #
  # NHANES public-use age is top-coded at 80 years.
  # Therefore, age is capped at 80 years for prediction.
  # ----------------------------------------------------------

  age80 <- pmin(
    age,
    80
  )


  # ----------------------------------------------------------
  # 3. NWR natural cubic spline
  #
  # Final spline specification:
  # 3 degrees of freedom
  # ----------------------------------------------------------

  nwr_knots <- c(
    0.570093457943925,
    0.652173913043478
  )

  nwr_boundary <- c(
    0.434782608695652,
    0.773607367475292
  )

  nwr_basis <- splines::ns(
    nwr,
    knots = nwr_knots,
    Boundary.knots = nwr_boundary,
    intercept = FALSE
  )

  colnames(nwr_basis) <- c(
    "nwr_ns1",
    "nwr_ns2",
    "nwr_ns3"
  )


  # ----------------------------------------------------------
  # 4. Final globally shrunken coefficients
  # ----------------------------------------------------------

  beta <- c(
    age80      =  0.035168687,
    bmi        = -0.019320732,
    rdw        =  0.187879211,
    nwr_ns1    =  0.671599209,
    nwr_ns2    =  0.163557831,
    nwr_ns3    =  0.473127774,
    egfr       = -0.008982533,
    med_count  =  0.052995926
  )


  # ----------------------------------------------------------
  # 5. Linear predictor
  # ----------------------------------------------------------

  lp <-
    beta["age80"]     * age80 +
    beta["bmi"]       * bmi +
    beta["rdw"]       * rdw +
    beta["nwr_ns1"]   * nwr_basis[, 1] +
    beta["nwr_ns2"]   * nwr_basis[, 2] +
    beta["nwr_ns3"]   * nwr_basis[, 3] +
    beta["egfr"]      * egfr +
    beta["med_count"] * med_count


  # ----------------------------------------------------------
  # 6. Frozen 1-year absolute-risk specification
  # ----------------------------------------------------------

  lp_center <- 4.415976

  H0_1year <- 0.05383994


  # ----------------------------------------------------------
  # 7. 1-year all-cause mortality risk
  # ----------------------------------------------------------

  risk_1year <-
    1 -
    exp(
      -H0_1year *
        exp(
          lp - lp_center
        )
    )


  # ----------------------------------------------------------
  # 8. Return
  # ----------------------------------------------------------

  data.frame(
    age_input = age,
    age_for_prediction = age80,
    bmi = bmi,
    rdw = rdw,
    nwr = nwr,
    egfr = egfr,
    med_count = med_count,
    linear_predictor = as.numeric(lp),
    risk_1year = as.numeric(risk_1year)
  )
}
