rem Copy fftw3.3.4 VC build to _hlib
rem This is a part of Rat and Catcher tech. cwave tools project

set HLIB="%CD%\..\..\_hlib"
del /f /s /q "%HLIB%"
mkdir "%HLIB%"

rem x86 binaries
call :platform

rem x64 binaries
cd .\x64
call :platform
cd ..

rem API header
copy ..\api\fftw3.h "%HLIB%"\
copy ..\api\fftw3.h "%HLIB%"\fftw334.h

exit /b

rem Copy binaries for the single platform
:platform
copy Debug\*.* "%HLIB%"\
copy Release\*.* "%HLIB%"\
copy Static-Debug\*.* "%HLIB%"\
copy Static-Release\*.* "%HLIB%"\
exit /b

