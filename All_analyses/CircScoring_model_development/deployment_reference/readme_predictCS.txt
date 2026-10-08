predictCS
=========

predictCS estimates the probability that each input circRNA belongs to the P1
positive class.

Two self-contained implementations are provided:

- predictCS.py
- predictCS.R

Both implementations contain the fitted models and do not read a tuning,
reference, or training workbook at prediction time.


Models
------

Both scripts produce scores from four fitted models:

1. XGBoost CS-R using features f16-f18
2. XGBoost CS-C using features f1-f15
3. Elastic Net CS-R using features f16-f18
4. Elastic Net CS-C using features f1-f15

Each script includes:

- The two fitted XGBoost boosters
- The two fitted Elastic Net coefficient sets
- The selected model parameters
- The f1-f18 feature-name mapping


Python requirements
-------------------

Use Python 3.10 or later. Install the required packages once:

    python3 -m pip install numpy "xgboost>=3.2.1"

On Windows, this command can also be used:

    py -3 -m pip install numpy "xgboost>=3.2.1"


R requirements
--------------

Install the required R packages once:

    install.packages(c("xgboost", "base64enc"))

The glmnet package is not required for prediction because the fitted Elastic
Net coefficients are embedded in both scripts.


Input files
-----------

Both predictCS.py and predictCS.R read a tab-delimited .txt file. The first
line must contain:

- circRNA_id
- All 18 predictor columns defined by the embedded f1-f18 feature mapping

Column names must be separated by tabs. Each following line represents one
circRNA.

Requirements:

- circRNA_id values must be present and unique.
- Every predictor column must be numeric.
- Use an empty field between tabs for a missing predictor value.
- Whitespace-only fields are also treated as missing.
- Other text, such as ".", "NA", or "N/A", is rejected in predictor columns.


Missing-feature behavior
------------------------

- In both implementations, every model requires all features for its mode.
  If any required feature is missing, that model's score is written as ".".
- CS-R and CS-C are checked independently. Missing a CS-R feature does not
  prevent CS-C prediction when all CS-C features are present, and vice versa.


Run predictCS.py with runs.bat
------------------------------

The included runs.bat launches predictCS.py with:

    Input:  script/predictCS/input_predictCS.txt
    Output: script/predictCS/output_predictCS.txt

On macOS or Linux, run from the project directory:

    bash script/predictCS/runs.bat

On Windows, double-click:

    script\predictCS\runs.bat

Alternatively, run it from Command Prompt:

    script\predictCS\runs.bat


Run predictCS.py directly
-------------------------

Command:

    python3 predictCS.py INPUT_TXT OUTPUT_TXT

Example:

    python3 script/predictCS/predictCS.py \
      script/predictCS/input_predictCS.txt \
      script/predictCS/output_predictCS.txt

On Windows, replace python3 with py -3 if needed.


Run predictCS.R directly
------------------------

Command:

    Rscript predictCS.R INPUT_TXT OUTPUT_TXT

Example:

    Rscript script/predictCS/predictCS.R \
      script/predictCS/input_predictCS.txt \
      script/predictCS/output_predictCS.txt

The text files must be tab-delimited. Worksheet names do not apply to either
implementation.

Output formats
--------------

Both predictCS.py and predictCS.R write a tab-delimited text file containing
only these columns:

    circRNA_id
    XGBoost_CS-R_score
    XGBoost_CS-C_score
    ElasticNet_CS-R_score
    ElasticNet_CS-C_score

Scores range from 0 to 1 and represent the estimated probability of the P1
positive class. Input feature columns are not copied to the output.


Console messages
----------------

Both scripts print:

1. Missing-feature skip counts for all four models. Each model is reported on
   a separate line.
2. The number of circRNAs scored and the saved output path.
