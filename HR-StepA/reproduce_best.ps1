param([double]$Beta=0.08,[int]$Delta=1100)
$cult=[Globalization.CultureInfo]::InvariantCulture

if (-not (Get-Command Build-HRPackageSweep -ErrorAction SilentlyContinue)) {
  throw "No encuentro Build-HRPackageSweep. Carga el módulo/ps1 donde la defines."
}

# Genera/actualiza artefactos para (Beta,Delta)
$sweepCsv = ".\out\summary_tmp_beta$($Beta.ToString($cult)).csv"
Build-HRPackageSweep -Beta $Beta -Deltas @($Delta) -CsvPath $sweepCsv | Out-Null

# Rellena lemma y compila
$betaStr=$Beta.ToString($cult)
$pw   = Get-Content ".\out\BS_PW_beta$betaStr`_D$Delta.json" | ConvertFrom-Json
$stepA= Get-Content .\out\StepA_results.json | ConvertFrom-Json

(Get-Content .\lemma_bs_insertion.tex -Raw) `
  -replace '\{BETA\}',$betaStr `
  -replace '\{DELTA\}',"$Delta" `
  -replace '\{N\}',"$Delta" `
  -replace '\{EPG\}', $pw.L1_errors.E_plus_grid.ToString('G12',$cult) `
  -replace '\{EMG\}', $pw.L1_errors.E_minus_grid.ToString('G12',$cult) `
  -replace '\{GAPp\}', $pw.valid.min_gap_majorant.ToString('G12',$cult) `
  -replace '\{GAPm\}', $pw.valid.min_gap_minorant.ToString('G12',$cult) `
  -replace '\{CBAJO\}', $stepA.constants.C_bajo.ToString('G12',$cult) `
  -replace '\{BVK\}',   $stepA.VK.BVK.ToString('G12',$cult) `
  -replace '\{BVKB\}',  $stepA.VK.beta_VK.ToString('G12',$cult) `
  -replace '\{XONE\}',  $stepA.VK.x1.ToString('G12',$cult) `
  -replace '\{CRHTOT\}',$stepA.constants_RH.C_tot_RH.ToString('G12',$cult) `
| Set-Content -Encoding UTF8 .\lemma_bs_insertion_filled.tex

pdflatex -interaction=nonstopmode .\lemma_bs_insertion_filled.tex | Out-Null
Write-Host "Repro listo para β=$betaStr, Δ=$Delta" -ForegroundColor Green
