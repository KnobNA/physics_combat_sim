# Fast GoogleTest suite for first-party physics units. Does not run Simulator::step.
param(
    [string]$BuildDir = "",
    [string]$Config = "Release"
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
if (-not $BuildDir) {
    $BuildDir = Join-Path $repo "build"
}

$cmake = Get-Command cmake -ErrorAction SilentlyContinue
if (-not $cmake) {
    $cmakePath = "C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
    if (-not (Test-Path $cmakePath)) {
        throw "cmake not found"
    }
} else {
    $cmakePath = $cmake.Source
}

$ctest = Get-Command ctest -ErrorAction SilentlyContinue
if (-not $ctest) {
    $ctestPath = "C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\ctest.exe"
    if (-not (Test-Path $ctestPath)) {
        throw "ctest not found"
    }
} else {
    $ctestPath = $ctest.Source
}

Write-Host "Configuring (GoogleTest FetchContent on first run)..."
& $cmakePath -S $repo -B $BuildDir
if ($LASTEXITCODE -ne 0) {
    throw "cmake configure failed with exit $LASTEXITCODE"
}

Write-Host "Building combat_sim_unit_tests ($Config)..."
& $cmakePath --build $BuildDir --config $Config --target combat_sim_unit_tests
if ($LASTEXITCODE -ne 0) {
    throw "cmake build failed with exit $LASTEXITCODE"
}

$exeDir = Join-Path $BuildDir "bin\$Config"
$env:PATH = $exeDir + ";" + $env:PATH
$env:MKL_INTERFACE_LAYER = "LP64"
$env:MKL_THREADING_LAYER = "TBB"

Write-Host "Running ctest --output-on-failure..."
& $ctestPath --test-dir $BuildDir -C $Config --output-on-failure
if ($LASTEXITCODE -ne 0) {
    throw "ctest failed with exit $LASTEXITCODE"
}
