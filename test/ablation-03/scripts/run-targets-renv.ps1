param(
  [ValidateSet('manifest', 'make', 'watch')]
  [string]$Action = 'make',
  [string]$Rscript = 'C:/R/R-4.3.1/bin/Rscript.exe'
)

$ErrorActionPreference = 'Continue'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Push-Location $project
try {
  & $Rscript --vanilla 'scripts/renv-ablation03.R' --command check
  if ($LASTEXITCODE -ne 0) { throw "ablation-03 renv check failed with exit code $LASTEXITCODE." }
  if ($Action -eq 'manifest') {
    & $Rscript --vanilla -e "source('renv/activate.R'); targets::tar_manifest(script = '_targets.R')"
  } elseif ($Action -eq 'watch') {
    & $Rscript --vanilla -e "source('renv/activate.R'); targets::tar_watch(config = '_targets.yaml', project = 'main', browse = TRUE)"
  } else {
    & $Rscript --vanilla -e "source('renv/activate.R'); targets::tar_make(script = '_targets.R')"
  }
  if ($LASTEXITCODE -ne 0) { throw "targets $Action failed with exit code $LASTEXITCODE." }
} finally { Pop-Location }
