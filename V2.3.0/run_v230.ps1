# V2.3.0 POD — Run Script
# Usage: powershell -File run_v230.ps1 -Mode 017|05|all|batch|fullday

param(
    [ValidateSet("017","05","all","batch","fullday")]
    [string]$Mode = "all"
)

$BASE = @(
    "--date", "2024-04-29", "--interval", "30", "--grace-id", "C",
    "--dynamics-mode", "simplified",
    "--sp3-file", "..\data\CODE\2024\COD0OPSFIN_20241200000_01D_05M_ORB.SP3",
    "--clk-file", "..\data\CODE\2024\COD0OPSFIN_20241200000_01D_30S_CLK.CLK",
    "--dcb-file", "..\data\CODE\2024\P1P22404.DCB",
    "--antex-file", "..\data\igs14.atx",
    "--iers-c04", "..\data\IERS\eopc04_IAU2000.txt",
    "--enable-phase-windup", "--enable-relativity",
    "--ar-min-epochs", "6", "--gravity-nmax", "90"
)

Write-Host "V2.3.0 POD — GRACE-FO C, 2024-04-29" -ForegroundColor Green

switch ($Mode) {
    "017" {
        Write-Host "0.17h EKF (target: 0.293m)"
        py -3.12 run_sequential_pod.py @BASE --hours 0.17
    }
    "05" {
        Write-Host "0.5h EKF (target: 0.986m)"
        py -3.12 run_sequential_pod.py @BASE --hours 0.5 --chi2-threshold 100
    }
    "all" {
        Write-Host "0.17h EKF"
        py -3.12 run_sequential_pod.py @BASE --hours 0.17
        Write-Host ""
        Write-Host "0.5h EKF"
        py -3.12 run_sequential_pod.py @BASE --hours 0.5 --chi2-threshold 100
    }
    "batch" {
        Write-Host "0.17h + Batch solver (Framework v3)"
        py -3.12 run_sequential_pod.py @BASE --hours 0.17 --batch-lsq-v2
    }
    "fullday" {
        Write-Host "Full-day assessment"
        py -3.12 validate_v224.py
    }
}
