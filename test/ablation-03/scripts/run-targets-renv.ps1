param(
  [ValidateSet('manifest', 'make', 'watch')]
  [string]$Action = 'make',
  [string]$InputRds = $env:CCS_ABLATION_INPUT_RDS,
  [string]$CacheRoot = $env:CCS_ABLATION_CACHE_ROOT,
  [string]$OutputRoot = '',
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
  $outputPath = [IO.Path]::GetFullPath((Join-Path $CacheRoot '01-data/inputs.rds'))
  if ([string]::Equals($env:CCS_ABLATION_INPUT_RDS, $outputPath,
      [StringComparison]::OrdinalIgnoreCase)) {
    throw 'The input RDS must not be the generated 01-data/inputs.rds file.'
  }
  $env:CCS_ABLATION_CACHE_ROOT = [IO.Path]::GetFullPath($CacheRoot)
  $projectTmp = [IO.Path]::GetFullPath((Join-Path $project 'tmp'))
  if ([string]::IsNullOrWhiteSpace($OutputRoot)) {
    $cachePath = $env:CCS_ABLATION_CACHE_ROOT
    $inProjectTmp = [string]::Equals($cachePath, $projectTmp,
      [StringComparison]::OrdinalIgnoreCase) -or
      $cachePath.StartsWith($projectTmp + [IO.Path]::DirectorySeparatorChar,
        [StringComparison]::OrdinalIgnoreCase)
    if ($inProjectTmp) {
      throw 'OutputRoot is required for an isolated cache under the project tmp directory.'
    }
    $env:CCS_ABLATION_OUTPUT_ROOT = $project
  } else {
    $env:CCS_ABLATION_OUTPUT_ROOT = [IO.Path]::GetFullPath($OutputRoot)
  }
  $env:CCS_ABLATION_TARGET_WORKERS = [string]([Math]::Max(1, $Workers))
  $localBin = Join-Path $env:USERPROFILE '.local\bin'
  if (Test-Path -LiteralPath (Join-Path $localBin 'pandoc.exe')) {
    $env:PATH = $localBin + [IO.Path]::PathSeparator + $env:PATH
  }
  & $Rscript --vanilla 'scripts/renv-ablation03.R' --command check
  if ($LASTEXITCODE -ne 0) { throw "ablation-03 renv check failed with exit code $LASTEXITCODE." }
  # targets resolves the store before sourcing _targets.R. Synchronize the
  # config first so a new cache root cannot silently reuse the previous store.
  & $Rscript --vanilla -e "source('renv/activate.R'); targets::tar_config_set(store = normalizePath(file.path(Sys.getenv('CCS_ABLATION_CACHE_ROOT'), 'targets'), winslash = '/', mustWork = FALSE))"
  if ($LASTEXITCODE -ne 0) { throw 'Failed to configure the ablation-03 targets store.' }
  if ($Action -eq 'manifest') {
    & $Rscript --vanilla -e "source('renv/activate.R'); targets::tar_manifest(script = '_targets.R')"
  } elseif ($Action -eq 'watch') {
    & $Rscript --vanilla -e "source('renv/activate.R'); watcher <- targets::tar_watch(config = '_targets.yaml', project = 'main', browse = TRUE); watcher`$wait(); watcher`$get_result()"
  } else {
    & $Rscript --vanilla -e "source('renv/activate.R'); targets::tar_make(script = '_targets.R')"
  }
  if ($LASTEXITCODE -ne 0) { throw "targets $Action failed with exit code $LASTEXITCODE." }
} finally {
  Pop-Location
}
