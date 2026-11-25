# Test data

Open data only. Every file is public and may be redistributed under the license shown.

| File | Source | License |
|---|---|---|
| `AF-P69905-F1-model_v6.cif`, `AF-P69905-F1-predicted_aligned_error_v6.json` | AlphaFold Protein Structure Database, human hemoglobin subunit alpha (UniProt P69905), model v6 | CC-BY 4.0 |
| `1L63.pdb`, `1L90.pdb` | RCSB Protein Data Bank entries 1L63 and 1L90 (T4 lysozyme) | CC0 |
| `NM_000558.5.gb` | NCBI RefSeq record of the human *HBA1* mRNA | public domain (NCBI) |
| `fold_hba_reference.zip`, `fold_hba_case.zip` | AlphaFold Server results of the HBA1 protein (UniProt P69905) and the case `c.179G>A` (p.Gly60Asp), run by the maintainer, one job each | AlphaFold Server output terms (`terms_of_use.md` in each zip) |
| `hba1_ref.zip`, `hba1_del_hinge.zip`, `hba1_trunc.zip` | Built by `../make_fixtures.py` from the AlphaFold DB model above | CC-BY 4.0 (derived) |

The three `hba1_*` zips imitate the layout of an AlphaFold Server result (five models, ranking scores, predicted aligned error) with C-alpha atoms only and seeded noise of 0.2 Å. The zips are test fixtures, not predictions.

- `hba1_ref.zip`: the AlphaFold Protein Structure Database model.
- `hba1_del_hinge.zip`: residues 56-62 deleted and residues 72 to the end rotated by 25° (two rigid blocks).
- `hba1_trunc.zip`: residues 1-39 removed.
