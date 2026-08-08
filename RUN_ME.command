#!/usr/bin/env bash
# ===================================================================
#  SynAT.Vis - run the functionality self-test on macOS / Linux.
#
#  macOS: double-click this file (RUN_ME.command). If it will not open,
#         first run once in Terminal:  chmod +x RUN_ME.command
#  Linux: run:  bash RUN_ME.command   (or chmod +x then ./RUN_ME.command)
#
#  Requires Python 3.8+ installed.
# ===================================================================
cd "$(dirname "$0")" || exit 1

if command -v python3 >/dev/null 2>&1; then
  PY=python3
elif command -v python >/dev/null 2>&1; then
  PY=python
else
  echo
  echo "  [X] Python was not found on this computer."
  echo "      Install Python 3.8+ (https://www.python.org/downloads/,"
  echo "      'brew install python' on macOS, or your package manager),"
  echo "      then try again."
  echo
  read -r -p "Press Enter to close..."
  exit 1
fi

echo "Using: $PY"
echo
"$PY" selftest.py
echo
echo "==================================================================="
echo ' Done. A result of "6/6 checks passed" means the tool works.'
echo "==================================================================="
echo
read -r -p "Press Enter to close..."
