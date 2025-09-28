param([string]$Path = 'BEST_manifest.json')

if (-not (Test-Path $Path)) { Write-Error "No existe $Path"; exit 1 }

try { $m = Get-Content -Raw $Path | ConvertFrom-Json }
catch {
  Write-Error "JSON inválido: $($_.Exception.Message)"; exit 1
}

$missing = @()
if (-not $m.assets)            { $missing += 'assets' }
if (-not $m.assets.best_zip)   { $missing += 'assets.best_zip' }
if (-not $m.assets.inputs_zip) { $missing += 'assets.inputs_zip' }
if (-not $m.assets.lemma_pdf)  { $missing += 'assets.lemma_pdf' }

if ($missing.Count) {
  Write-Error "Manifest incompleto. Faltan: $($missing -join ', ')"; exit 1
} else {
  'Manifest OK'
}
