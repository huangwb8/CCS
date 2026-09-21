param(
  [ValidateSet('manifest', 'make', 'watch')]
  [string]$Action = 'make',
  [string]$InputRds = $env:CCS_ABLATION_INPUT_RDS,
  [string]$CacheRoot = $env:CCS_ABLATION_CACHE_ROOT,
  [int]$Workers = 2,
  [string]$Rscript = 'C:/R/R-4.3.1/bin/Rscript.exe'
)

$ErrorActionPreference = 'Stop'
$project = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$repoRoot = (Resolve-Path (Join-Path $project '..\..')).Path
Push-Location $project
try {
  if ([string]::IsNullOrWhiteSpace($InputRds)) {
    throw 'InputRds is required. Formal and test runs must provide an input RDS; this is the only scientific data variation.'
  }
  if ([string]::IsNullOrWhiteSpace($CacheRoot)) {
    $CacheRoot = Join-Path $project 'tmp/targets-cache'
  }
  $inputCandidate = if (Test-Path -LiteralPath $InputRds) {
    $InputRds
  } else {
    Join-Path $repoRoot $InputRds
  }
  if (-not (Test-Path -LiteralPath $inputCandidate)) {
    throw "InputRds does not exist: $InputRds"
  }
  $env:CCS_ABLATION_INPUT_RDS = (Resolve-Path -LiteralPath $inputCandidate).Path
  $env:CCS_ABLATION_CACHE_ROOT = $CacheRoot
  $env:CCS_ABLATION_TARGET_WORKERS = [string]([Math]::Max(1, $Workers))
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
} finally {
  Pop-Location
}
