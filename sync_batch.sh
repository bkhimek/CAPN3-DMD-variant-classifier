#!/usr/bin/env bash
# Syncs the newest CAPN3-DMD-variant-classifier batch zip from Windows into
# the WSL git repo, verifies the copy actually changed something, and runs
# the test suite -- before you ever get to `git commit`.
#
# Usage:
#   bash sync_batch.sh
#
# One-time setup (optional but recommended):
#   cp sync_batch.sh ~/projects/CAPN3-DMD-variant-classifier/sync_batch.sh
#   chmod +x ~/projects/CAPN3-DMD-variant-classifier/sync_batch.sh
# Then from then on, from inside the repo:
#   ./sync_batch.sh

set -euo pipefail

REPO_DIR="${REPO_DIR:-$HOME/projects/CAPN3-DMD-variant-classifier}"
SEARCH_ROOT="${SEARCH_ROOT:-/mnt/c/Users/krist}"
SCRATCH="/tmp/capn3_dmd_sync"

echo "== Locating newest batch zip under $SEARCH_ROOT =="
ZIP=$(find "$SEARCH_ROOT" -iname "CAPN3-DMD-variant-classifier_batch*.zip" -printf '%T@ %p\n' 2>/dev/null \
      | sort -rn | head -1 | cut -d' ' -f2-)

if [ -z "$ZIP" ]; then
  echo "ERROR: no batch zip found under $SEARCH_ROOT. Nothing to sync."
  exit 1
fi
echo "Using: $ZIP"

echo
echo "== Extracting to $SCRATCH =="
rm -rf "$SCRATCH"
mkdir -p "$SCRATCH"
unzip -q -o "$ZIP" -d "$SCRATCH"

EXTRACTED_ROOT="$SCRATCH/CAPN3-DMD-variant-classifier"
if [ ! -d "$EXTRACTED_ROOT" ]; then
  echo "ERROR: expected $EXTRACTED_ROOT after extraction, but it's not there."
  exit 1
fi

echo
echo "== Copying into $REPO_DIR =="
cp -r "$EXTRACTED_ROOT/." "$REPO_DIR/"

cd "$REPO_DIR"

echo
echo "== git status =="
git status --short

CHANGED=$(git status --porcelain | wc -l)
if [ "$CHANGED" -eq 0 ]; then
  echo
  echo "ERROR: no changes detected after copying $ZIP into the repo."
  echo "Either this zip's contents are already committed, or the copy silently"
  echo "failed. STOPPING before running tests or allowing a commit -- there is"
  echo "nothing new to commit right now, so investigate rather than proceeding."
  exit 1
fi

echo
echo "== Running tests =="
PYTHONPATH=src python3 tests/run_tests.py

echo
echo "== Diff summary =="
git diff --stat

echo
echo "Changes above are staged for your review (not yet committed)."
echo "If the diff and test results look right, commit with:"
echo "  git add -A && git commit -m \"<describe this batch>\" && git push"
