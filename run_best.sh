@echo off
REM ============================================================
REM GRACE-FO POD V2.2.3 — 最佳精度运行脚本 (Windows PowerShell)
REM ============================================================
REM
REM 用法: 在 d:\prj\gnss_pod 目录下执行
REM   powershell -File run_best.ps1
REM 或在 PowerShell 中逐条复制运行
REM
REM 数据要求: 见 VERSION.md 第2节

echo ============================================================
echo V2.2.3 最佳配置 — 0.17h (10 min) 弧段
echo ============================================================

set JAVA_HOME=C:\Program Files\JetBrains\PyCharm Community Edition 2024.3.5\jbr
set OREKIT_DATA_PATH=data\orekit

py -3.12 run_sequential_pod.py ^
  --date 2024-04-29 --hours 0.17 --interval 30 --grace-id C ^
  --dynamics-mode simplified ^
  --sp3-file data\CODE\2024\COD0OPSFIN_20241200000_01D_05M_ORB.SP3 ^
  --clk-file data\CODE\2024\COD0OPSFIN_20241200000_01D_30S_CLK.CLK ^
  --dcb-file data\CODE\2024\P1P22404.DCB ^
  --antex-file data\igs14.atx ^
  --iers-c04 data\IERS\eopc04_IAU2000.txt ^
  --enable-phase-windup --enable-relativity ^
  --ar-min-epochs 6 --gravity-nmax 90

echo.
echo ============================================================
echo V2.2.3 最佳配置 — 0.5h (30 min) 弧段
echo ============================================================

py -3.12 run_sequential_pod.py ^
  --date 2024-04-29 --hours 0.5 --interval 30 --grace-id C ^
  --dynamics-mode simplified ^
  --sp3-file data\CODE\2024\COD0OPSFIN_20241200000_01D_05M_ORB.SP3 ^
  --clk-file data\CODE\2024\COD0OPSFIN_20241200000_01D_30S_CLK.CLK ^
  --dcb-file data\CODE\2024\P1P22404.DCB ^
  --antex-file data\igs14.atx ^
  --iers-c04 data\IERS\eopc04_IAU2000.txt ^
  --enable-phase-windup --enable-relativity ^
  --ar-min-epochs 6 --chi2-threshold 100 --gravity-nmax 90

echo.
echo ============================================================
echo V2.2.3 全天线评估 (6时段 × 0.5h, 生成 accuracy PNG)
echo ============================================================

py -3.12 eval_day2.py

echo.
echo 结果保存在:
echo   0.17h: results\sequential_ekf\seq_2024-04-29_C_0.17h.pkl
echo   0.5h:  results\sequential_ekf\seq_2024-04-29_C_0.5h.pkl
echo   全天线: results\accuracy_2024-04-29.png
