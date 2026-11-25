# FoldShift: methods and references

## Method

- Global sequence alignment (affine gaps, free end gaps) pairs reference and case residues.
- Superposition uses identical residues with pLDDT of at least 70 and iterative pruning of pairs beyond 2 Å (a standard outlier-pruning fit). The fit is started from the whole set and from each contiguous alignment block, so a domain
  that swings does not spoil it. Every block is also fitted alone and reported next to the global-fit RMSD.
- Both RMSDs are reported: over the fitted pairs and over all aligned pairs. The other AlphaFold models of the
  reference, compared with its best model, give the noise level of the method.
- Residues are paired with BLOSUM62 (identity pairing is an option). lDDT-C-alpha (cutoffs 0.5, 1, 2, 4 A, inclusion
  radius 15 A) does not depend on the superposition. Residues that moved more than 5 A but have pLDDT below the fit limit in
  either model are flagged "moved, low confidence". Helix, strand and loop are assigned from C-alpha geometry for both models
  and changed stretches (three or more residues) are listed; this is a simple geometric rule, not DSSP.
- Several chains are read in file order, bonds are never drawn across chains, and the chain-aware residue table is in the CSV.
- Start loss (`c.1A>G`) assumes translation from the first downstream in-frame Met; check this.
- Figures come from a built-in software renderer (ribbons from C-alpha geometry, soft ambient shading and outlines).
  Colors are Okabe-Ito and viridis.

## Develop

    cargo test --workspace --release      # uses the open data in test/data (included)
    cargo clippy --workspace --all-targets -- -D warnings && cargo fmt --check
    cargo run --release -p psplot         # --no-browser --port N for tests
    CARGO_TARGET_DIR=/tmp/target scripts/build_windows.sh   # Windows x64 and ARM64 exes into dist/
    cargo build --release -p psplot --features window       # optional native window (Windows, macOS); not tested yet

macOS (Apple Silicon) is built by `.github/workflows/release.yml` on a GitHub macOS runner (it cannot be cross-built from
Linux).

## Independence and licenses

All code is written for this project (MIT license, see `LICENSE`). The feature set was chosen by reading what PyMOL and
ChimeraX offer (superposition with pruning, per-residue distances, ribbon views). No source code, library, or script of
either program is copied, linked, called, or shipped. The embedded font is DejaVu Sans, see `THIRD_PARTY_NOTICES.md`.

### Methods and references

Each method below was implemented from its publication. Identifiers were checked against PubMed, the publisher, or the
journal record.

- Affine-gap sequence alignment (Gotoh): Gotoh, O. (1982). An improved algorithm for matching biological sequences.
  *Journal of Molecular Biology, 162*(3), 705-708. https://doi.org/10.1016/0022-2836(82)90398-9
- Least-squares superposition with unit quaternions (Horn): Horn, B. K. P. (1987). Closed-form solution of absolute
  orientation using unit quaternions. *Journal of the Optical Society of America A, 4*(4), 629-642.
  https://doi.org/10.1364/JOSAA.4.000629
- BLOSUM62 residue pairing: Henikoff, S., & Henikoff, J. G. (1992). Amino acid substitution matrices from protein blocks.
  *Proceedings of the National Academy of Sciences, 89*(22), 10915-10919. https://doi.org/10.1073/pnas.89.22.10915 (PMID 1438297)
- lDDT: Mariani, V., Biasini, M., Barbato, A., & Schwede, T. (2013). lDDT: a local superposition-free score for comparing
  protein structures and models using distance difference tests. *Bioinformatics, 29*(21), 2722-2728.
  https://doi.org/10.1093/bioinformatics/btt473 (PMID 23986568). We use the thresholds 0.5, 1, 2 and 4 Å and the 15 Å
  inclusion radius from the paper, but on C-alpha atoms only (the paper's default uses all atoms and notes a C-alpha
  option), so the value is reported as lDDT-C-alpha and is not comparable with all-atom lDDT. It is computed over aligned residue pairs only.
- TM-score: Zhang, Y., & Skolnick, J. (2004). Scoring function for automated assessment of protein structure template
  quality. *Proteins, 57*(4), 702-710. https://doi.org/10.1002/prot.20264 (PMID 15476259; erratum
  https://doi.org/10.1002/prot.21643). We use eq. 1 and eq. 5 of the paper (d0 = 1.24 (L - 15)^(1/3) - 1.8, L = reference
  length, with L set to at least 19 to keep d0 positive), but evaluate them under our own superposition, not the maximum over
  all superpositions that the paper defines. The value is therefore a lower bound of the TM-score.
- AlphaFold Server models (the input of the comparison): Abramson, J., et al. (2024). Accurate structure prediction of
  biomolecular interactions with AlphaFold 3. *Nature, 630*(8016), 493-500. https://doi.org/10.1038/s41586-024-07487-w
  (PMID 38718835)

### Our own rules, without a literature source

These are heuristics of this program and have not been validated against published benchmarks:
- Helix and strand assignment from C-alpha distances (helix: i to i+3 and i+4 distances; strand: i to i+2 above 6.2 Å;
  runs shorter than 5 and 3 residues dropped). Tested only on ideal helix and strand coordinates, not against DSSP.
- Iterative pruning of residue pairs beyond 2 Å and the pLDDT limit of 70 for the fit.
- The "moved and low confidence" flag (more than 5 Å and pLDDT below the fit limit) and the wording thresholds of the
  verdict (RMSD under 2 Å, under 4 Å).
- The color palettes follow the lab plotting rules (Okabe-Ito, viridis); no reference is cited here.
