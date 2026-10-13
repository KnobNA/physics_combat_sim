# Sequential physics-smoke sweep. MKL/Pardiso does not like two smokes at once.
param(
    [string]$BuildDir = "",
    [string]$Config = "Release",
    [int]$TimeoutSec = 180,
    [string]$Only = ""
)

$ErrorActionPreference = "Continue"
$repo = Split-Path -Parent $PSScriptRoot
if (-not $BuildDir) {
    $BuildDir = Join-Path $repo "build"
}
$exe = Join-Path $BuildDir "bin\$Config\combat_sim_physics_smoke.exe"

$cmake = Get-Command cmake -ErrorAction SilentlyContinue
if (-not $cmake) {
    $cmakePath = "C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
    if (-not (Test-Path $cmakePath)) {
        throw "cmake not found"
    }
} else {
    $cmakePath = $cmake.Source
}

Write-Host "Building combat_sim_physics_smoke ($Config)..."
& $cmakePath --build $BuildDir --config $Config --target combat_sim_physics_smoke
if ($LASTEXITCODE -ne 0) {
    throw "cmake build failed with exit $LASTEXITCODE"
}
if (-not (Test-Path $exe)) {
    throw "missing smoke binary: $exe"
}

$env:PATH = (Split-Path $exe) + ";" + $env:PATH

$stamp = Get-Date -Format "yyyyMMdd-HHmmss"
$outDir = Join-Path $repo "logs\smoke\$stamp"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null
Write-Host "Logs: $outDir"

$groups = @(
    "mesh", "motion", "ductile", "brittle", "compliant", "joints", "heat", "remesh", "editor"
)
if ($Only -and $Only.Trim().Length -gt 0) {
    $groups = @($Only.Split(',') | ForEach-Object { $_.Trim() } | Where-Object { $_ })
}

function Parse-SmokeLog {
    param([string]$text, [string]$group, [int]$exitCode, [double]$seconds)
    $pass = 0
    $fail = 0
    $failed = New-Object System.Collections.Generic.List[string]
    $dumps = New-Object System.Collections.Generic.List[string]
    $checks = New-Object System.Collections.Generic.List[hashtable]
    foreach ($line in ($text -split "`r?`n")) {
        if ($line -match '^CHECK group=(\S+) name="(.*)" result=(ok|FAIL)(?: got=(\S+) want=\[(\S+),(\S+)\])?') {
            $ok = $Matches[3] -eq "ok"
            if ($ok) { $pass++ } else { $fail++ }
            $entry = @{
                group  = $Matches[1]
                name   = $Matches[2]
                result = $Matches[3]
            }
            if ($Matches[4]) {
                $entry.got = $Matches[4]
                $entry.wantLo = $Matches[5]
                $entry.wantHi = $Matches[6]
            }
            [void]$checks.Add($entry)
            if (-not $ok) { [void]$failed.Add($line) }
        } elseif ($line -match 'FAILED' -and $line -notmatch '^CHECK ') {
            [void]$failed.Add($line)
        } elseif ($line -match '^TIMEOUT ') {
            [void]$failed.Add($line)
        } elseif ($line -match '^DUMP ') {
            [void]$dumps.Add($line)
        }
    }
    return [pscustomobject]@{
        group   = $group
        exit    = $exitCode
        seconds = [math]::Round($seconds, 2)
        pass    = $pass
        fail    = $fail
        failed  = $failed
        dumps   = $dumps
        checks  = $checks
    }
}

$results = New-Object System.Collections.Generic.List[object]
$allPath = Join-Path $outDir "ALL.log"
Set-Content -Path $allPath -Value "" -Encoding utf8

foreach ($group in $groups) {
    $logPath = Join-Path $outDir "$group.log"
    Write-Host ("`n==== {0} ====" -f $group)
    $started = Get-Date
    $limit = $TimeoutSec
    if ($group -in @("motion", "ductile", "brittle", "compliant", "editor", "heat", "joints", "remesh", "mesh")) {
        $limit = [Math]::Max($TimeoutSec, 180)
    }
    if ($group -eq "motion" -or $group -eq "ductile" -or $group -eq "editor" -or $group -eq "compliant") {
        $limit = [Math]::Max($limit, 300)
    }
    $p = Start-Process -FilePath $exe -ArgumentList "--only=$group" -WorkingDirectory (Split-Path $exe) `
        -RedirectStandardOutput $logPath -RedirectStandardError "$logPath.err" -NoNewWindow -PassThru
    $timedOut = -not $p.WaitForExit($limit * 1000)
    if ($timedOut) {
        Write-Host ("TIMEOUT {0} after {1}s - killing hung integrator" -f $group, $limit)
        try { Stop-Process -Id $p.Id -Force -ErrorAction SilentlyContinue } catch {}
        Get-Process combat_sim_physics_smoke -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue
        Start-Sleep -Milliseconds 400
        $code = 124
    } else {
        $code = $p.ExitCode
        if ($null -eq $code) { $code = 0 }
    }
    $seconds = ((Get-Date) - $started).TotalSeconds
    $stdout = ""
    $stderr = ""
    if (Test-Path $logPath) { $stdout = Get-Content -Path $logPath -Raw -ErrorAction SilentlyContinue }
    if (Test-Path "$logPath.err") { $stderr = Get-Content -Path "$logPath.err" -Raw -ErrorAction SilentlyContinue }
    $text = (($stdout + "`n" + $stderr).TrimEnd())
    if ($timedOut) {
        $text = $text + "`nTIMEOUT group=$group after ${limit}s (integrator did not return)`n"
    }
    Set-Content -Path $logPath -Value $text -Encoding utf8
    Remove-Item -Path "$logPath.err" -ErrorAction SilentlyContinue
    Add-Content -Path $allPath -Value ("`n===== GROUP {0} exit={1} wall={2:N1}s =====`n{3}" -f $group, $code, $seconds, $text) -Encoding utf8
    $parsed = Parse-SmokeLog -text $text -group $group -exitCode $code -seconds $seconds
    [void]$results.Add($parsed)
    Write-Host ("{0}: exit={1} {2:N1}s pass={3} fail={4}" -f $group, $code, $seconds, $parsed.pass, $parsed.fail)
}

$md = New-Object System.Text.StringBuilder
[void]$md.AppendLine("# Physics smoke suite")
[void]$md.AppendLine()
[void]$md.AppendLine("- stamp: $stamp")
[void]$md.AppendLine("- exe: $exe")
[void]$md.AppendLine("- meshResolution: 0.25 (coarse) except editor / motion M2-M3 (1.0, dynres on); remesh keeps dynres")
[void]$md.AppendLine()
[void]$md.AppendLine("| group | seconds | exit | pass | fail |")
[void]$md.AppendLine("|---|---:|---:|---:|---:|")
foreach ($row in $results) {
    [void]$md.AppendLine(("| {0} | {1} | {2} | {3} | {4} |" -f $row.group, $row.seconds, $row.exit, $row.pass, $row.fail))
}
[void]$md.AppendLine()
[void]$md.AppendLine("## Failures")
$anyFail = $false
foreach ($row in $results) {
    if ($row.exit -ne 0 -or $row.fail -gt 0 -or $row.failed.Count -gt 0) {
        $anyFail = $true
        [void]$md.AppendLine(("### {0} (exit {1})" -f $row.group, $row.exit))
        foreach ($line in $row.failed) {
            [void]$md.AppendLine("- $line")
        }
        [void]$md.AppendLine()
    }
}
if (-not $anyFail) {
    [void]$md.AppendLine("None.")
    [void]$md.AppendLine()
}
[void]$md.AppendLine("## dumpSim snapshots")
foreach ($row in $results) {
    if ($row.dumps.Count -eq 0) { continue }
    [void]$md.AppendLine(("### {0}" -f $row.group))
    foreach ($line in $row.dumps) {
        [void]$md.AppendLine("- $line")
    }
    [void]$md.AppendLine()
}

$summaryMd = Join-Path $outDir "SUMMARY.md"
$summaryJson = Join-Path $outDir "summary.json"
Set-Content -Path $summaryMd -Value $md.ToString() -Encoding utf8

$jsonRows = foreach ($row in $results) {
    [pscustomobject]@{
        group   = $row.group
        seconds = $row.seconds
        exit    = $row.exit
        pass    = $row.pass
        fail    = $row.fail
        failed  = @($row.failed)
        dumps   = @($row.dumps)
        checks  = @($row.checks)
    }
}
$jsonRows | ConvertTo-Json -Depth 6 | Set-Content -Path $summaryJson -Encoding utf8

Write-Host "`nWrote $summaryMd"
Write-Host "Wrote $summaryJson"
