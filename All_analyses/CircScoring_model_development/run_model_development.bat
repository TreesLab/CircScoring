@echo off
setlocal
cd /d "%~dp0"
Rscript CircScoring_model_development.R --input data\P1N1_training_data.tsv --output results --threads 1
if errorlevel 1 exit /b %errorlevel%
endlocal

