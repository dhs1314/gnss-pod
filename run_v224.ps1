# ============================================================
# GRACE-FO POD V2.2.4 — 运行脚本
# 用法: powershell -File run_v224.ps1
#       powershell -File run_v224.ps1 -Mode 017    (仅 0.17h)
#       powershell -File run_v224.ps1 -Mode 05     (仅 0.5h)
#       powershell -File run_v224.ps1 -Mode fullday (全天线评估)
# ============================================================
param(
    [ValidateSet("017","05","all","fullday")]
    [string]$Mode = "all"
)

$ErrorActionPreference = "Stop"

# Environment
$env:JAVA_HOME = "C:\Program Files\JetBrains\PyCharm Community Edition 2024.3.5\jbr"
$env:OREKIT_DATA_PATH = "d:\prj\gnss_pod\data\orekit"

# Common arguments
$COMMON = @(
    "--date", "2024-04-29",
    "--interval", "30",
    "--grace-id", "C",
    "--dynamics-mode", "simplified",
    "--sp3-file", "d:\prj\gnss_pod\data\CODE\2024\COD0OPSFIN_20241200000_01D_05M_ORB.SP3",
    "--clk-file", "d:\prj\gnss_pod\data\CODE\2024\COD0OPSFIN_20241200000_01D_30S_CLK.CLK",
    "--dcb-file", "d:\prj\gnss_pod\data\CODE\2024\P1P22404.DCB",
    "--antex-file", "d:\prj\gnss_pod\data\igs14.atx",
    "--iers-c04", "d:\prj\gnss_pod\data\IERS\eopc04_IAU2000.txt",
    "--enable-phase-windup", "--enable-relativity",
    "--ar-min-epochs", "6",
    "--gravity-nmax", "90"
)

function Run-017 {
    Write-Host "=" * 60 -ForegroundColor Cyan
    Write-Host "V2.2.4 短弧段 (0.17h / 10min) — 目标: 0.293m" -ForegroundColor Cyan
    Write-Host "=" * 60 -ForegroundColor Cyan
    py -3.12 d:\prj\gnss_pod\run_sequential_pod.py @COMMON --hours 0.17
}

function Run-05 {
    Write-Host "=" * 60 -ForegroundColor Cyan
    Write-Host "V2.2.4 长弧段 (0.5h / 30min) — 目标: 0.986m" -ForegroundColor Cyan
    Write-Host "=" * 60 -ForegroundColor Cyan
    py -3.12 d:\prj\gnss_pod\run_sequential_pod.py @COMMON --hours 0.5 --chi2-threshold 100
}

function Run-FullDay {
    Write-Host "=" * 60 -ForegroundColor Cyan
    Write-Host "V2.2.4 全天线评估 (6 时段 x 0.5h) — 输出 accuracy_2024-04-29.png" -ForegroundColor Cyan
    Write-Host "=" * 60 -ForegroundColor Cyan
    py -3.12 d:\prj\gnss_pod\eval_day2.py
}

Write-Host "V2.2.4 POD — GRACE-FO C, 2024-04-29" -ForegroundColor Green
Write-Host ""

switch ($Mode) {
    "017"     { Run-017 }
    "05"      { Run-05 }
    "all"     { Run-017; Write-Host ""; Run-05 }
    "fullday" { Run-FullDay }
}

Write-Host ""
Write-Host "结果保存位置:" -ForegroundColor Yellow
Write-Host "  0.17h: results\sequential_ekf\seq_2024-04-29_C_0.17h.pkl"
Write-Host "  0.5h:  results\sequential_ekf\seq_2024-04-29_C_0.5h.pkl"
Write-Host "  全天线: results\accuracy_2024-04-29.png"
