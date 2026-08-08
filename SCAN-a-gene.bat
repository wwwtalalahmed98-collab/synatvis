@echo off
REM ===================================================================
REM  SynAT.Vis - DRAG YOUR GENE FILE ONTO THIS ICON to scan it.
REM  (Drag a .fasta / .fa / .gb file onto SCAN-a-gene.bat)
REM  No typing needed. Requires Python 3.8+ installed.
REM ===================================================================
setlocal EnableExtensions
cd /d "%~dp0"
title SynAT.Vis - scan a gene

REM --- find Python by actually running it (robust vs App Execution Aliases) ---
set "PYEXE="
py --version >nul 2>&1 && set "PYEXE=py"
if not defined PYEXE ( python --version >nul 2>&1 && set "PYEXE=python" )
if not defined PYEXE ( python3 --version >nul 2>&1 && set "PYEXE=python3" )
if not defined PYEXE (
  echo.
  echo   [X] Python was not found. Install Python 3.8+ from
  echo       https://www.python.org/downloads/  - tick "Add to PATH" -
  echo       then try again.
  echo.
  pause
  exit /b 1
)

REM --- no file dropped? explain ---
if "%~1"=="" (
  echo.
  echo   HOW TO USE:
  echo     Drag your gene file ^(.fasta / .fa / .gb^) and drop it ONTO
  echo     this "SCAN-a-gene.bat" icon. The report opens automatically
  echo     and is also saved next to your file.
  echo.
  pause & exit /b 0
)

REM --- scan each dropped file ---
:loop
if "%~1"=="" goto done
echo.
echo === Scanning: %~nx1 ===
set "HTML=%~dpn1_SynAT.Vis_report.html"
set "TXT=%~dpn1_SynAT.Vis_report.txt"
%PYEXE% -m synatvis scan "%~1" --html --out "%HTML%"
%PYEXE% -m synatvis scan "%~1" --plain > "%TXT%" 2>&1
echo   Opening your visual report in the browser...
start "" "%HTML%"
echo   (Saved: %HTML%  and a plain-text copy: %TXT%)
shift
goto loop

:done
echo.
echo ===================================================================
echo  Done. Read the SERIOUS items first. See USER_GUIDE.pdf for more.
echo ===================================================================
echo.
pause
