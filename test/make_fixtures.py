#!/usr/bin/env python3
"""Build the AlphaFold Server style test zips in test/data from an open AlphaFold DB model.

Input:  data/AF-P69905-F1-model_v6.cif (human hemoglobin subunit alpha, AlphaFold DB, CC-BY 4.0).
Output: data/hba1_ref.zip       the model itself (C-alpha atoms only)
        data/hba1_del_hinge.zip residues 56-62 deleted, residues 72 to the end rotated by 25 degrees
        data/hba1_trunc.zip     residues 1-39 removed (N-terminal truncation)
Each zip holds five models (with a synthetic PAE matrix) with 0.2 A of seeded noise and ranking scores that make model 2 the best.
Run: python3 test/make_fixtures.py   (standard library only)
"""
import io, json, math, random, zipfile
from pathlib import Path

D = Path(__file__).parent / "data"
AA = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
      "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
      "TYR": "Y", "VAL": "V"}


def read_ca(path):
    out = []
    for l in path.read_text().splitlines():
        f = l.split()
        if l.startswith("ATOM") and f[3] == "CA":
            out.append((f[5], [float(f[10]), float(f[11]), float(f[12])], float(f[14])))
    return out


def rotate(p, origin, axis, deg):
    n = math.sqrt(sum(a * a for a in axis))
    k = [a / n for a in axis]
    t = math.radians(deg)
    v = [p[i] - origin[i] for i in range(3)]
    kv = [k[1] * v[2] - k[2] * v[1], k[2] * v[0] - k[0] * v[2], k[0] * v[1] - k[1] * v[0]]
    d = sum(k[i] * v[i] for i in range(3))
    r = [v[i] * math.cos(t) + kv[i] * math.sin(t) + k[i] * d * (1 - math.cos(t)) for i in range(3)]
    return [r[i] + origin[i] for i in range(3)]


def cif(res, rng, noise):
    s = ("data_model\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n"
         "_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n"
         "_atom_site.Cartn_z\n_atom_site.B_iso_or_equiv\n")
    for i, (name, xyz, b) in enumerate(res, 1):
        x, y, z = (c + rng.gauss(0, noise) for c in xyz)
        s += f"ATOM {i} CA {name} A {i} {x:.3f} {y:.3f} {z:.3f} {b:.2f}\n"
    return s + "#\n"


def write(name, res):
    seq = "".join(AA[r[0]] for r in res)
    rng = random.Random(7)
    with zipfile.ZipFile(D / f"{name}.zip", "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(f"{name}_job_request.json", json.dumps([{"name": name, "sequences": [{"proteinChain": {"sequence": seq, "count": 1}}]}]))
        for i, score in enumerate([0.80, 0.85, 0.90, 0.70, 0.60]):
            z.writestr(f"{name}_model_{i}.cif", cif(res, rng, 0.2))
            n = len(res)
            pae = [[round(1 + 0.05 * abs(a - b), 2) for b in range(n)] for a in range(n)]  # synthetic: grows with distance in sequence
            z.writestr(f"{name}_full_data_{i}.json", json.dumps({"pae": pae}))
            z.writestr(f"{name}_summary_confidences_{i}.json", json.dumps({"ranking_score": score, "ptm": 0.9}))


ca = read_ca(D / "AF-P69905-F1-model_v6.cif")
write("hba1_ref", ca)
pivot = ca[70][1]  # residue 71
moved = [(n, rotate(x, pivot, (0.3, 0.8, 0.5), 25) if k >= 71 else x, b) for k, (n, x, b) in enumerate(ca)]
write("hba1_del_hinge", [r for k, r in enumerate(moved) if not 55 <= k <= 61])
write("hba1_trunc", ca[39:])
