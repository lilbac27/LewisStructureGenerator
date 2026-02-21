$ErrorActionPreference = "Stop"

$testDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$srcDir = Join-Path $testDir "..\src"
$outExe = Join-Path $testDir "lewis_engine_tests.exe"

$sources = @(
    (Join-Path $srcDir "lewis_model.c"),
    (Join-Path $srcDir "lewis_engine.c"),
    (Join-Path $testDir "lewis_engine_tests.c")
)

if (Get-Command clang -ErrorAction SilentlyContinue) {
    & clang -std=c11 -Wall -Wextra -O2 -I $srcDir @sources -o $outExe
} elseif (Get-Command gcc -ErrorAction SilentlyContinue) {
    & gcc -std=c11 -Wall -Wextra -O2 -I $srcDir @sources -o $outExe
} elseif (Get-Command zig -ErrorAction SilentlyContinue) {
    & zig cc -std=c11 -Wall -Wextra -O2 -I $srcDir @sources -o $outExe
} else {
    Write-Error "No host C compiler found (clang/gcc/zig cc)."
}

& $outExe
