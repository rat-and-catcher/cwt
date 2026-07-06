@echo off
rem BIG check one of FFTW build variant

if "%~1" == "" goto usage

if "%~2" == "" (set LOGF=NUL) else set LOGF=%~2

perl -w check.pl -a -v %1
if errorlevel 1 goto fail

echo FFTW transforms passed BIG tests
echo FFTW transforms passed BIG tests >>%LOGF%

perl -w check.pl -a -v --nthreads=2 %1
if errorlevel 1 goto fail
perl -w check.pl -a -v --nthreads=3 %1
if errorlevel 1 goto fail
perl -w check.pl -a -v --nthreads=10 %1
if errorlevel 1 goto fail

echo FFTW threaded transforms passed BIG tests
echo FFTW threaded transforms passed BIG tests >>%LOGF%

exit /b 0

:usage
echo Usage: %~nx0 bench-program.exe [log-file-name]
exit /b 0

:fail
echo FFTW transforms FAILED BIG tests
echo FFTW transforms FAILED BIG tests >>%LOGF%
exit /b 1
