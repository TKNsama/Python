@echo off
setlocal enabledelayedexpansion

REM --- 设置路径 ---
set MEGACC_PATH="D:\software\ParaAT2.0\D\megacc_12.0.14_win64\megacc.exe"
set SETTINGS_MAO="D:\software\ParaAT2.0\D\megacc_12.0.14_win64\tajima_neutrality_test_nucleotide.mao"
set INPUT_DIR="D:\software\ParaAT2.0\D\filtered_out_region"
set OUTPUT_DIR="D:\software\ParaAT2.0\D\TajimaD_results_region"

if not exist %OUTPUT_DIR% (
    mkdir %OUTPUT_DIR%
)

REM --- 批量处理每个 fasta 文件 ---
for %%F in (%INPUT_DIR%\*.fasta %INPUT_DIR%\*.fa) do (
    echo Processing %%~nF ...
    "%MEGACC_PATH%" -a %SETTINGS_MAO% -d "%%F" -o "%OUTPUT_DIR%\%%~nF_TajimaD.txt"
)

echo All files processed.
pause
