# Path configuration for public scripts

The public scripts avoid user-specific absolute paths.

Set these environment variables before running analyses:

```r
Sys.setenv(
  BASE6_ANALYSIS_DIR = "/path/to/your/development-analysis-working-directory",
  BASE6_EXTERNAL_DATA = "/path/to/authorized/deidentified/external_dataset.xlsx",
  BASE6_EXTERNAL_OUTPUT_DIR = "/path/to/private/external-analysis-output-directory"
)
```

`BASE6_ANALYSIS_DIR` contains NHANES-derived intermediate and model files.
`BASE6_EXTERNAL_DATA` points to the authorized deidentified hospital dataset.
`BASE6_EXTERNAL_OUTPUT_DIR` must be outside the public repository and is used for restricted external-cohort intermediate/output files.

External individual-level data are not distributed with the repository.
