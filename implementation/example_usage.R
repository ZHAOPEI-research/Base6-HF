# ============================================================
# Example use of the Base6 1-year mortality model
#
# Example values below are fictional and are provided
# solely to demonstrate model implementation.
# ============================================================


source(
  "implementation/Base6_predict_1year.R"
)


example_prediction <- predict_base6_1year(
  age = 72,
  bmi = 26,
  rdw = 14.5,
  nwr = 0.65,
  egfr = 60,
  med_count = 8
)


print(
  example_prediction
)
