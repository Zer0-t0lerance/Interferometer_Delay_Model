# build_pybind.ps1 - build the Python module ariadna (pybind11) -> ariadna.pyd (repo root).
# ASCII-only (Windows PowerShell 5.1 reads BOM-less .ps1 as ANSI).
#
# Needs: g++ (MinGW) in PATH; pybind11 (pip install pybind11); Python with dev headers.
# Python paths are taken automatically from the current `python`. libsofa.a must exist
# (build.ps1 builds it on first run).
#
# Run:
#   powershell -ExecutionPolicy Bypass -File .\python\build_pybind.ps1
# Then, from the repo root:  python -c "import ariadna; print(ariadna.__doc__)"

$ErrorActionPreference = "Stop"
Set-Location "$PSScriptRoot\.."

$PybindInc = & python -c "import pybind11; print(pybind11.get_include())"
$PyInc     = & python -c "import sysconfig; print(sysconfig.get_path('include'))"
$PyLibDir  = & python -c "import sys,os; print(os.path.join(sys.base_prefix,'libs'))"
$PyLib     = & python -c "import sysconfig; print('python'+sysconfig.get_config_var('py_version_nodot'))"

$inc = @("-I./external/eigen","-I./external","-I./external/dephem-master/include",
         "-I./external/sofa/20190722/c/src","-I$PybindInc","-I$PyInc")
$dep  = "./external/dephem-master/include/dephem/EphemerisRelease.cpp"
$sofa = "./external/sofa/20190722/c/build/libsofa.a"
$src  = (Get-ChildItem "src\*.cpp").FullName

Write-Host "Building ariadna.pyd (linking $PyLib) ..."
$gccArgs = @("-O2","-shared","-std=c++17") + $inc +
        @("python/ariadna_pybind.cpp") + $src + @($dep, $sofa,
        "-L$PyLibDir", "-l$PyLib", "-o", "ariadna.pyd")
& g++ @gccArgs
if ($LASTEXITCODE -ne 0) { Write-Host "BUILD FAILED."; exit 1 }
Write-Host "OK -> ariadna.pyd"
Write-Host "Check:  python python/example.py"
