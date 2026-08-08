#!/usr/bin/env bash
# ===================================================================
#  SynAT.Vis - scan a gene file (macOS / Linux).
#
#  macOS: double-click this file. A window opens to pick your gene
#         file, then the report appears. (First time, if it won't
#         open: in Terminal run  chmod +x SCAN-a-gene.command )
#  Linux: run  bash SCAN-a-gene.command  (uses a file picker if
#         'zenity' is installed, otherwise asks you to paste a path)
#
#  No typing needed on macOS. Requires Python 3.8+ installed.
# ===================================================================
cd "$(dirname "$0")" || exit 1

# --- find Python ---
if command -v python3 >/dev/null 2>&1; then PY=python3
elif command -v python >/dev/null 2>&1; then PY=python
else
  echo "  [X] Python was not found. Install Python 3.8+ from"
  echo "      https://www.python.org/downloads/ (or 'brew install python'),"
  echo "      then try again."
  read -r -p "Press Enter to close..."
  exit 1
fi

# --- pick a file (GUI where possible) ---
FILE=""
if command -v osascript >/dev/null 2>&1; then           # macOS
  FILE=$(osascript -e 'try
    POSIX path of (choose file with prompt "Pick your gene file (.fasta / .fa / .gb)")
  end try')
elif command -v zenity >/dev/null 2>&1; then             # Linux with zenity
  FILE=$(zenity --file-selection --title="Pick your gene file" 2>/dev/null)
fi
if [ -z "$FILE" ]; then
  echo
  read -r -p "Paste the full path to your gene file (.fasta), then Enter: " FILE
fi
FILE="${FILE/#\~/$HOME}"
if [ ! -f "$FILE" ]; then
  echo "  [X] File not found: $FILE"
  read -r -p "Press Enter to close..."
  exit 1
fi

# --- scan: write a visual HTML report + a plain-text copy, then open the report ---
dir=$(dirname "$FILE"); base=$(basename "$FILE"); name="${base%.*}"
HTML="$dir/${name}_SynAT.Vis_report.html"
TXT="$dir/${name}_SynAT.Vis_report.txt"
echo
echo "=== Scanning: $base ==="
"$PY" -m synatvis scan "$FILE" --html --out "$HTML"
"$PY" -m synatvis scan "$FILE" --plain > "$TXT" 2>&1
echo "  Opening your visual report in the browser..."
( command -v open >/dev/null 2>&1 && open "$HTML" ) || \
  ( command -v xdg-open >/dev/null 2>&1 && xdg-open "$HTML" ) || \
  echo "  Open this file in your browser: $HTML"
echo "  (Saved: $HTML  and a plain-text copy: $TXT)"
echo
echo "==================================================================="
echo " Done. Read the SERIOUS items first. See USER_GUIDE.pdf for more."
echo "==================================================================="
echo
read -r -p "Press Enter to close..."
