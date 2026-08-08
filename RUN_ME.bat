@echo off
REM ===================================================================
REM  SynAT.Vis - double-click to run the functionality self-test.
REM  No terminal knowledge needed. Requires Python 3.8+ installed.
REM ===================================================================
setlocal
cd /d "%~dp0"
title SynAT.Vis self-test

echo Locating Python...
set "PYEXE="
py --version >nul 2>&1 && set "PYEXE=py"
if not defined PYEXE ( python --version >nul 2>&1 && set "PYEXE=python" )
if not defined PYEXE ( python3 --version >nul 2>&1 && set "PYEXE=python3" )

if not defined PYEXE (
  echo.
  echo   [X] Python was not found on this computer.
  echo       Install Python 3.8+ from https://www.python.org/downloads/
  echo       and tick "Add python.exe to PATH" during setup, then try again.
  echo.
  pause
  exit /b 1
)

echo Using: %PYEXE%
echo.
%PYEXE% selftest.py
echo.
echo ===================================================================
echo  Done. A result of "6/6 checks passed" means the tool works.
echo ===================================================================
echo.
pause
