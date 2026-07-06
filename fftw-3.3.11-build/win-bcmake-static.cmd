@echo off
echo Create MS Visual Studio projects tree
setlocal

set FFTW_MSVS="Visual Studio 15 2017"
rem set FFTW_MSVS="NMake Makefiles"
set PFMS=x64 Win32
set FFTW_VER=3.3.11

rem "Useful" mean "double"

call win-clear-static.cmd
if errorlevel 1 goto fail

for %%B in (%PFMS%) do (
    for %%A in (GEN SSE2 AVX AVX2) do (
    call :mk_cfg %%A %%B || goto fail
    )
)

exit /b 0

:fail
echo Failed.
exit /b 1

:mk_cfg
set WD=build-%1-%2-static
echo cmake %WD%
mkdir %WD% && cd %WD%
if errorlevel 1 exit /b 1
if %1 == GEN (set ARC="-DEMPTY=ON") else set ARC=-DENABLE_%1=ON
set VCMAKE=cmake -G %FFTW_MSVS% -A %2 -DFFTW_VERSION=%FFTW_VER% -DBUILD_SHARED_LIBS=OFF -DENABLE_THREADS=ON -DENABLE_RATCAT=ON %ARC% ..
echo %VCMAKE%
%VCMAKE%
if errorlevel 1 cd .. ; exit /b 1
cd ..
set BCMAKE=cmake --build %WD% --config=Release
echo %BCMAKE%
%BCMAKE%
set BCMAKE=cmake --build %WD% --config=Debug
echo %BCMAKE%
%BCMAKE%
if errorlevel 1 exit /b 1
exit /b 0

