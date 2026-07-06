@echo off
rem Paranoid check one of FFTW build variant

if "%~1" == "" goto usage

if "%~2" == "" (set LOGF=NUL) else set LOGF=%~2

perl -w check.pl -a --patient --nthreads=10 --paranoid %1
if errorlevel 1 goto fail
perl -w check.pl -a --patient --nthreads=7 --paranoid %1
if errorlevel 1 goto fail
perl -w check.pl -a --patient --nthreads=3 --paranoid %1
if errorlevel 1 goto fail
perl -w check.pl -a --patient --nthreads=2 --paranoid %1
if errorlevel 1 goto fail

echo FFTW threaded transforms passed Paranoid tests
echo FFTW threaded transforms passed Paranoid tests >>%LOGF%

perl -w check.pl -a --patient --paranoid %1
if errorlevel 1 goto fail

echo FFTW transforms passed Paranoid tests
echo FFTW transforms passed Paranoid tests >>%LOGF%

exit /b 0

:usage
echo Usage: %~nx0 bench-program.exe [log-file-name]
exit /b 0

:fail
echo FFTW transforms FAILED Paranoid tests
echo FFTW transforms FAILED Paranoid tests >>%LOGF%
exit /b 1
