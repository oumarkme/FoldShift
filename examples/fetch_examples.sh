#!/bin/sh
# Download the structures used in docs/EXAMPLES.md from the RCSB Protein Data Bank into examples/data/.
# Usage: sh examples/fetch_examples.sh
set -eu
dir="$(dirname "$0")/data"
mkdir -p "$dir"
# id.format: the format used in the examples (some entries are only convenient as .pdb)
for f in 2BBO.cif 2BBT.cif 2BBS.cif 5UAK.cif 8EJ1.cif 1UWH.pdb 1UWJ.pdb 4OBE.pdb 5US4.pdb \
         2HHB.cif 2HBS.cif 2OCJ.cif 2J1X.cif 1L63.pdb 1L90.pdb 2B6Z.pdb 6BDD.pdb 6BDE.pdb; do
  [ -s "$dir/$f" ] || curl -fsSL -o "$dir/$f" "https://files.rcsb.org/download/$f"
  echo "$f"
done
