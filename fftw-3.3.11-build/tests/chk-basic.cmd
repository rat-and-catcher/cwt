@echo off
rem Basic check one of FFTW build variant

if "%~1" == "" goto usage

if "%~2" == "" (set LOGF=NUL) else set LOGF=%~2

perl -w check.pl -r -c=30 -v %1
if errorlevel 1 goto fail

echo FFTW transforms passed basic tests
echo FFTW transforms passed basic tests >>%LOGF%

perl -w check.pl -r -c=30 -v --nthreads=2 %1
if errorlevel 1 goto fail

echo FFTW threaded transforms passed basic tests
echo FFTW threaded transforms passed basic tests >>%LOGF%

exit /b 0

:usage
echo Usage: %~nx0 bench-program.exe [log-file-name]
exit /b 0

:fail
echo FFTW transforms FAILED basic tests
echo FFTW transforms FAILED basic tests >>%LOGF%
exit /b 1
