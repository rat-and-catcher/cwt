@echo off
echo Completely remove MS Visual Studio projects tree
setlocal

for %%B in (x64 Win32) do (
    for %%A in (GEN SSE2 AVX AVX2) do (
    call :rm_tree %%A %%B || goto fail
    )
)

exit /b 0

:fail
echo Failed.
exit /b 1

:rm_tree
set WD=build-%1-%2
echo rmdir %WD%
:del /s /f /a:RSHAIO %WD%\*.* >NUL 2>&1
rmdir /s /q %WD% >NUL 2>&1
if errorlevel 1 exit /b 1
exit /b 0
