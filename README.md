# Base6 Heart Failure Mortality Prediction Model

This repository contains analysis code and implementation files for the study:

**Development and External Evaluation of a Parsimonious Mortality Prediction Model for Heart Failure: The Base6 Model**

Base6 is a 6-predictor Cox proportional hazards model developed to estimate the 1-year probability of all-cause mortality in adults with heart failure.

## Model predictors

Base6 uses the following six predictors:

1. Age, years
2. Body mass index (BMI), kg/m²
3. Red cell distribution width (RDW), %
4. Neutrophil-to-white-cell ratio (NWR)
5. Estimated glomerular filtration rate (eGFR), mL/min/1.73 m²
6. Medication count

NWR is defined as:

`absolute neutrophil count / total white blood cell count`

Age is capped at 80 years for prediction because age is top-coded at 80 years in the NHANES public-

## Repository structure

- `code/` — Analysis scripts used for model development, internal validation, external evaluation, recalibration, sensitivity analyses, and subgroup analyses
- `implementation/`
  - `Base6_predict_1year.R` — Frozen Base6 1-year risk implementation
  - `example_usage.R` — Example using fictional predictor values
- `metadata/`
  - Variable definitions and software/package information
- `RUN_ORDER.md` — Tested execution order for the public analysis workflow
- `LICENSE` — MIT License

## Final model implementation

The primary implementation file is:
implementation/Base6_predict_1year.R
It contains the frozen final Base6 specification, including:
- globally shrunken coefficients
- the 3-degree-of-freedom natural cubic spline for NWR
- the age cap at 80 years
- the linear-predictor centering constant
- the frozen 1-year baseline cumulative hazard
- the final 1-year absolute-risk equation
The final 1-year mortality risk is calculated as:
1 - exp[-H0(1 year) × exp(LP - LP_center)]
with:
- H0(1 year) = 0.05383994
- LP_center = 4.415976
## Example use
The file:
implementation/example_usage.R
contains a simple example using fictional predictor values.
From the repository root, run:
source("implementation/example_usage.R")

The example values are fictional and are provided solely to demonstrate model implementation.
## Input requirements
Base6 requires all six predictors for individual risk estimation.
Predictions should not be generated when required predictor values are missing or considered unreliable unless a prespecified and validated approach for handling missing inputs is available.
Direct manual calculation of the spline component is not expected. The model is intended to be implemented using software or an electronic calculator.
## Development cohort
The model was developed using NHANES 2011–2018 data linked to public-use mortality files.
The final development cohort included:
- 739 participants with heart failure
- 237 deaths
Multiple imputation, grouped penalized Cox regression, functional-form assessment, bootstrap internal validation, and global shrinkage were used during model development.
## External evaluation cohort
External evaluation was performed in an independent hospitalized heart-failure cohort from China.
The common three-model external evaluation cohort included:
- 346 participants
- 47 deaths during approximately 1 year of follow-up
Base6 was compared with MAGGIC and GWTG-HF.
## Recalibration
Post hoc recalibration was performed in the external evaluation cohort.
For Base6:
z_Base6 = log[-log(1 - p_Base6)]
where p_Base6 is the original 1-year predicted mortality risk.
Updated risks were calculated as:
1 - exp{-exp(alpha + beta × z_Base6)}
with:
- intercept-only recalibration:
  - alpha = 0.5770
  - beta = 1
- intercept-and-slope recalibration:
  - alpha = 0.6918
  - beta = 1.0519
These recalibrated models represent model updating and should not be interpreted as independently validated models.
## Data availability
NHANES data are publicly available from the National Center for Health Statistics.
Individual-level data from the external hospital cohort are not included in this repository because of institutional and participant-privacy restrictions.
Deidentified external-cohort data may be available from the corresponding author for legitimate research purposes, subject to institutional approval and applicable data-use requirements.
## Reproducibility
The public repository is intended to provide:
- the final Base6 implementation
- the main analysis workflow
- code for internal validation
- code for external evaluation
- code for recalibration
- sensitivity and subgroup analyses
Some analyses depend on source datasets that are not redistributed in this repository.
## Software
Analyses were performed using R version 4.6.0.
A complete list of relevant package versions will be provided in:
metadata/package_versions.txt
## Intended use
Base6 is intended for research and evaluation purposes.
The model should not be used for clinical decision-making without further independent multicenter validation and local assessment of calibration.
## Citation
The associated manuscript is currently under preparation/submission.
A permanent citation and Zenodo DOI will be added after the repository is archived.

## License
This repository is released under the MIT License. See `LICENSE` for details.