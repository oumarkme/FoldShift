//! Sequence alignment, superposition, per-residue C-alpha distances.
#![allow(clippy::needless_range_loop)]

use crate::ss::{assign_ss, Ss};
use crate::structure::{AfResult, Model, V3};
use serde::Serialize;

pub type M3 = [[f64; 3]; 3];

pub fn apply(r: &M3, t: &V3, p: &V3) -> V3 {
    [
        r[0][0] * p[0] + r[0][1] * p[1] + r[0][2] * p[2] + t[0],
        r[1][0] * p[0] + r[1][1] * p[1] + r[1][2] * p[2] + t[1],
        r[2][0] * p[0] + r[2][1] * p[1] + r[2][2] * p[2] + t[2],
    ]
}
pub fn dist(a: &V3, b: &V3) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// Symmetric eigen decomposition of a 4x4 matrix (cyclic Jacobi). Returns the eigenvector of the largest eigenvalue.
fn top_eigenvector(mut a: [[f64; 4]; 4]) -> [f64; 4] {
    let mut v = [[0.0; 4]; 4];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for _ in 0..60 {
        let off: f64 = (0..4)
            .flat_map(|i| (0..4).map(move |j| (i, j)))
            .filter(|(i, j)| i != j)
            .map(|(i, j)| a[i][j] * a[i][j])
            .sum();
        if off < 1e-24 {
            break;
        }
        for p in 0..3 {
            for q in p + 1..4 {
                if a[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..4 {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..4 {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for row in v.iter_mut() {
                    let (vkp, vkq) = (row[p], row[q]);
                    row[p] = c * vkp - s * vkq;
                    row[q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let best = (0..4)
        .max_by(|&i, &j| a[i][i].partial_cmp(&a[j][j]).unwrap())
        .unwrap();
    [v[0][best], v[1][best], v[2][best], v[3][best]]
}

/// Least-squares rotation R and translation t with R p + t ~ q (Horn 1987).
pub fn kabsch(p: &[V3], q: &[V3]) -> (M3, V3) {
    let n = p.len() as f64;
    let mean = |x: &[V3]| {
        let mut m = [0.0; 3];
        for a in x {
            for k in 0..3 {
                m[k] += a[k] / n;
            }
        }
        m
    };
    let (pc, qc) = (mean(p), mean(q));
    let mut s = [[0.0; 3]; 3];
    for (a, b) in p.iter().zip(q) {
        for i in 0..3 {
            for j in 0..3 {
                s[i][j] += (a[i] - pc[i]) * (b[j] - qc[j]);
            }
        }
    }
    let (sxx, sxy, sxz, syx, syy, syz, szx, szy, szz) = (
        s[0][0], s[0][1], s[0][2], s[1][0], s[1][1], s[1][2], s[2][0], s[2][1], s[2][2],
    );
    let nmat = [
        [sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
        [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
        [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
        [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
    ];
    let [q0, q1, q2, q3] = top_eigenvector(nmat);
    let r = [
        [
            q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3,
            2.0 * (q1 * q2 - q0 * q3),
            2.0 * (q1 * q3 + q0 * q2),
        ],
        [
            2.0 * (q1 * q2 + q0 * q3),
            q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3,
            2.0 * (q2 * q3 - q0 * q1),
        ],
        [
            2.0 * (q1 * q3 - q0 * q2),
            2.0 * (q2 * q3 + q0 * q1),
            q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3,
        ],
    ];
    let rp = apply(&r, &[0.0; 3], &pc);
    (r, [qc[0] - rp[0], qc[1] - rp[1], qc[2] - rp[2]])
}

pub fn rmsd(a: &[V3], b: &[V3]) -> f64 {
    if a.is_empty() {
        return f64::NAN;
    }
    (a.iter()
        .zip(b)
        .map(|(x, y)| dist(x, y).powi(2))
        .sum::<f64>()
        / a.len() as f64)
        .sqrt()
}

/// BLOSUM62 in the order ARNDCQEGHILKMFPSTWYV.
const BLOSUM62: [[i8; 20]; 20] = [
    [
        4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0,
    ],
    [
        -1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3,
    ],
    [
        -2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3,
    ],
    [
        -2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3,
    ],
    [
        0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1,
    ],
    [
        -1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2,
    ],
    [
        -1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2,
    ],
    [
        0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3,
    ],
    [
        -2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3,
    ],
    [
        -1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3,
    ],
    [
        -1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1,
    ],
    [
        -1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2,
    ],
    [
        -1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1,
    ],
    [
        -2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1,
    ],
    [
        -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2,
    ],
    [
        1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2,
    ],
    [
        0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0,
    ],
    [
        -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3,
    ],
    [
        -2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1,
    ],
    [
        0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4,
    ],
];

fn aa_index(b: u8) -> Option<usize> {
    b"ARNDCQEGHILKMFPSTWYV".iter().position(|&x| x == b)
}

fn pair_score(a: u8, b: u8, blosum: bool) -> f32 {
    if blosum {
        match (aa_index(a), aa_index(b)) {
            (Some(i), Some(j)) => BLOSUM62[i][j] as f32,
            _ => {
                if a == b {
                    1.0
                } else {
                    -1.0
                }
            }
        }
    } else if a == b {
        2.0
    } else {
        -1.0
    }
}

/// Global alignment (Gotoh, affine gaps, free end gaps). Returns 0-based (ref, case) pairs of aligned residues.
pub fn align_sequences(a: &[u8], b: &[u8], blosum: bool) -> Vec<(usize, usize)> {
    const NEG: f32 = -1e9;
    let (n, m) = (a.len(), b.len());
    let (open, ext) = if blosum {
        (-11.0f32, -1.0f32)
    } else {
        (-5.0f32, -0.5f32)
    };
    let w = m + 1;
    // pointer bytes: H: 0 diag, 1 from E, 2 from F, 3 stop; E and F: 0 open (from H), 1 extend
    let mut ph = vec![3u8; (n + 1) * w];
    let mut pe = vec![0u8; (n + 1) * w];
    let mut pf = vec![0u8; (n + 1) * w];
    let mut h_prev = vec![0.0f32; w];
    let mut e_prev = vec![NEG; w];
    let mut f_prev = vec![NEG; w];
    let mut best = (0.0f32, 0usize, 0usize);
    for i in 1..=n {
        let mut h_cur = vec![0.0f32; w];
        let mut e_cur = vec![NEG; w];
        let mut f_cur = vec![NEG; w];
        for j in 1..=m {
            // E: gap in ref (consumes case residue j): from (i, j-1)
            let (eo, ee) = (h_cur[j - 1] + open, e_cur[j - 1] + ext);
            e_cur[j] = eo.max(ee);
            pe[i * w + j] = u8::from(ee > eo);
            // F: gap in case (consumes ref residue i): from (i-1, j)
            let (fo, fe) = (h_prev[j] + open, f_prev[j] + ext);
            f_cur[j] = fo.max(fe);
            pf[i * w + j] = u8::from(fe > fo);
            let s = pair_score(a[i - 1], b[j - 1], blosum);
            let d = h_prev[j - 1] + s;
            let (mut h, mut p) = (d, 0u8);
            if e_cur[j] > h {
                h = e_cur[j];
                p = 1;
            }
            if f_cur[j] > h {
                h = f_cur[j];
                p = 2;
            }
            h_cur[j] = h;
            ph[i * w + j] = p;
            if (i == n || j == m) && h > best.0 {
                best = (h, i, j);
            }
        }
        h_prev = h_cur;
        e_prev = e_cur;
        f_prev = f_cur;
    }
    let _ = (&e_prev, &f_prev);
    let (_, mut i, mut j) = best;
    let mut pairs = Vec::new();
    let mut state = 0u8; // 0 H, 1 E, 2 F
    while i > 0 && j > 0 {
        match state {
            0 => match ph[i * w + j] {
                0 => {
                    pairs.push((i - 1, j - 1));
                    i -= 1;
                    j -= 1;
                }
                1 => state = 1,
                2 => state = 2,
                _ => break,
            },
            1 => {
                let ext_ = pe[i * w + j] == 1;
                j -= 1;
                if !ext_ {
                    state = 0;
                }
            }
            _ => {
                let ext_ = pf[i * w + j] == 1;
                i -= 1;
                if !ext_ {
                    state = 0;
                }
            }
        }
    }
    pairs.reverse();
    pairs
}

/// "1-210, 250" -> mask over `n` reference residues (1-based input).
pub fn parse_ranges(text: &str, n: usize) -> Result<Vec<bool>, String> {
    let mut mask = vec![false; n];
    for tok in text
        .split(|c: char| c == ',' || c.is_whitespace())
        .filter(|t| !t.is_empty())
    {
        let (a, b) = match tok.split_once('-') {
            Some((a, b)) => (a.parse::<usize>(), b.parse::<usize>()),
            None => (tok.parse::<usize>(), tok.parse::<usize>()),
        };
        let (Ok(a), Ok(b)) = (a, b) else {
            return Err(format!(
                "Cannot read the residue range '{tok}'. Use for example 1-210, 250."
            ));
        };
        for r in a.max(1)..=b.min(n) {
            mask[r - 1] = true;
        }
    }
    Ok(mask)
}

#[derive(Debug, Clone)]
pub struct Options {
    pub fit_region: String,
    pub min_plddt: f64,
    pub cutoff: f64,
    pub robustness: bool,
    /// Residue pairing by BLOSUM62 (default) or by exact match.
    pub blosum: bool,
}
impl Default for Options {
    fn default() -> Self {
        Options {
            fit_region: String::new(),
            min_plddt: 70.0,
            cutoff: 2.0,
            robustness: true,
            blosum: true,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct Block {
    pub ref_start: i32,
    pub ref_end: i32,
    pub case_start: i32,
    pub case_end: i32,
    pub n: usize,
    pub rmsd_local_fit: f64,
    pub rmsd_global_fit: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct Robust {
    pub comparison: String,
    pub rmsd_all_aligned: f64,
    pub fit_rmsd: f64,
    pub median_distance: f64,
    pub max_distance: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct Stats {
    pub ref_length: usize,
    pub case_length: usize,
    pub aligned_residues: usize,
    pub sequence_identity_of_aligned: f64,
    pub fit_residues: usize,
    pub fit_rmsd: f64,
    pub rmsd_all_aligned: f64,
    pub median_distance: f64,
    pub max_distance: f64,
    pub fraction_over_5a: f64,
    pub fraction_over_10a: f64,
    pub tm_score_under_this_fit: f64,
    pub ref_ranking_score: Option<f64>,
    pub case_ranking_score: Option<f64>,
    pub ref_ptm: Option<f64>,
    pub case_ptm: Option<f64>,
    pub ref_mean_plddt: f64,
    pub case_mean_plddt: f64,
    pub ref_model: usize,
    pub case_model: usize,
    pub min_plddt_for_fit: f64,
    pub fit_cutoff_a: f64,
    pub fit_region: String,
    pub lddt_ca: f64,
    pub low_confidence_shifts: usize,
    pub ss_changed_residues: usize,
    pub chains_ref: usize,
    pub chains_case: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct SsChange {
    pub ref_start: i32,
    pub ref_end: i32,
    pub from: char,
    pub to: char,
}

#[derive(Debug, Clone)]
pub struct Comparison {
    /// Residues were paired by their numbers in the files (not by sequence alignment).
    pub paired_by_number: bool,
    /// Residue numbers shown for the reference and the case (see `numbering`), one per residue.
    pub ref_num: Vec<i32>,
    pub case_num: Vec<i32>,
    pub ref_source: String,
    pub case_source: String,
    pub ref_name: String,
    pub case_name: String,
    pub ref_seq: Vec<u8>,
    pub case_seq: Vec<u8>,
    pub ref_model: Model,
    pub case_model: Model,
    pub pairs: Vec<(usize, usize)>,
    pub dist: Vec<f64>,
    pub in_fit: Vec<bool>,
    pub block: Vec<usize>,
    pub r: M3,
    pub t: V3,
    pub case_fit: Vec<V3>,
    pub stats: Stats,
    pub blocks: Vec<Block>,
    pub robustness: Vec<Robust>,
    pub missing_in_case: Vec<(usize, usize)>,
    pub extra_in_case: Vec<(usize, usize)>,
    pub lddt: Vec<f64>,
    pub low_conf: Vec<bool>,
    pub ss_ref: Vec<Ss>,
    pub ss_case: Vec<Ss>,
    pub ss_changes: Vec<SsChange>,
}

fn n_chains(m: &Model) -> usize {
    1 + m.chain.windows(2).filter(|w| w[0] != w[1]).count()
}

fn ranges(mut idx: Vec<usize>) -> Vec<(usize, usize)> {
    idx.sort_unstable();
    let mut out: Vec<(usize, usize)> = Vec::new();
    for i in idx {
        match out.last_mut() {
            Some(l) if i == l.1 + 1 => l.1 = i,
            _ => out.push((i, i)),
        }
    }
    out
}

/// Iterative outlier pruning: refit on pairs closer than `cutoff` until stable.
fn fit_core(p: &[V3], q: &[V3], cutoff: f64) -> (M3, V3, Vec<bool>) {
    let mut keep = vec![true; p.len()];
    for _ in 0..10 {
        let (pk, qk) = pick(p, q, &keep);
        let (r, t) = kabsch(&pk, &qk);
        let new: Vec<bool> = p
            .iter()
            .zip(q)
            .map(|(a, b)| dist(&apply(&r, &t, a), b) <= cutoff)
            .collect();
        if new.iter().filter(|x| **x).count() < 10 || new == keep {
            break;
        }
        keep = new;
    }
    let (pk, qk) = pick(p, q, &keep);
    let (r, t) = kabsch(&pk, &qk);
    (r, t, keep)
}
fn pick(p: &[V3], q: &[V3], keep: &[bool]) -> (Vec<V3>, Vec<V3>) {
    let it = || keep.iter().enumerate().filter(|(_, k)| **k).map(|(i, _)| i);
    (it().map(|i| p[i]).collect(), it().map(|i| q[i]).collect())
}

fn median(v: &[f64]) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = s.len();
    if n == 0 {
        f64::NAN
    } else if n % 2 == 1 {
        s[n / 2]
    } else {
        (s[n / 2 - 1] + s[n / 2]) / 2.0
    }
}
fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len().max(1) as f64
}

/// Pair residues by their numbers in the two files (same number, same chain selection). Only for plain structure
/// files whose numbers are unique and increase along the chain, and only when it agrees with the sequence alignment
/// on at least 90 % of its pairs and pairs identical residues in at least 90 % of the cases. Residues missing from
/// one file (disordered loops that differ between two crystal structures) then cannot shift the pairing.
fn pair_by_number(
    reference: &AfResult,
    case: &AfResult,
    by_sequence: &[(usize, usize)],
) -> Option<Vec<(usize, usize)>> {
    if reference.source != "structure file" || case.source != "structure file" {
        return None;
    }
    let (rm, cm) = (
        &reference.models[reference.best()],
        &case.models[case.best()],
    );
    let increasing = |v: &[i32]| v.windows(2).all(|w| w[0] < w[1]);
    if !increasing(&rm.resnums) || !increasing(&cm.resnums) {
        return None;
    }
    let pos: std::collections::HashMap<i32, usize> = cm
        .resnums
        .iter()
        .enumerate()
        .map(|(j, &n)| (n, j))
        .collect();
    let pairs: Vec<(usize, usize)> = rm
        .resnums
        .iter()
        .enumerate()
        .filter_map(|(i, n)| pos.get(n).map(|&j| (i, j)))
        .collect();
    if pairs.len() < 10 || by_sequence.is_empty() {
        return None;
    }
    let (rs, cs) = (reference.sequence.as_bytes(), case.sequence.as_bytes());
    let same = pairs.iter().filter(|&&(i, j)| rs[i] == cs[j]).count();
    let set: std::collections::HashSet<&(usize, usize)> = pairs.iter().collect();
    let agree = by_sequence.iter().filter(|p| set.contains(p)).count();
    (same as f64 >= 0.9 * pairs.len() as f64 && agree as f64 >= 0.9 * by_sequence.len() as f64)
        .then_some(pairs)
}

pub fn compare(reference: &AfResult, case: &AfResult, opt: &Options) -> Result<Comparison, String> {
    let (ri, ci) = (reference.best(), case.best());
    let by_sequence = align_sequences(
        reference.sequence.as_bytes(),
        case.sequence.as_bytes(),
        opt.blosum,
    );
    let numbered = pair_by_number(reference, case, &by_sequence);
    let paired_by_number = numbered.is_some();
    let pairs = numbered.unwrap_or(by_sequence);
    if pairs.len() < 10 {
        return Err("Fewer than 10 residues align between the two sequences. Are the two zips from the same protein?".into());
    }
    let mut c = core(
        reference,
        case,
        &reference.models[ri],
        &case.models[ci],
        pairs.clone(),
        opt,
    )?;
    c.paired_by_number = paired_by_number;
    if opt.robustness {
        let mut rows = Vec::new();
        for (label, other, same) in [("case", case, false), ("reference (self)", reference, true)] {
            for m in &other.models {
                if same && m.index == reference.models[ri].index {
                    continue;
                }
                let p: Vec<(usize, usize)> = if same {
                    (0..reference.sequence.len()).map(|i| (i, i)).collect()
                } else {
                    pairs.clone()
                };
                let other_r = if same { reference } else { case };
                let k = core(
                    reference,
                    other_r,
                    &reference.models[ri],
                    m,
                    p,
                    &Options {
                        robustness: false,
                        ..opt.clone()
                    },
                )?;
                rows.push(Robust {
                    comparison: format!(
                        "{label} model_{} vs reference model_{} (best ranked)",
                        m.index, reference.models[ri].index
                    ),
                    rmsd_all_aligned: k.stats.rmsd_all_aligned,
                    fit_rmsd: k.stats.fit_rmsd,
                    median_distance: k.stats.median_distance,
                    max_distance: k.stats.max_distance,
                });
            }
        }
        c.robustness = rows;
    }
    Ok(c)
}

/// The residue numbers shown to the user: the numbers of the file (author numbering, for example 508 for Phe508)
/// when they are unique, otherwise the position 1, 2, 3, ... in the list.
pub fn numbering(resnums: &[i32]) -> Vec<i32> {
    let mut seen = std::collections::HashSet::new();
    if resnums.iter().all(|&r| seen.insert(r)) {
        resnums.to_vec()
    } else {
        (1..=resnums.len() as i32).collect()
    }
}

fn core(
    reference: &AfResult,
    case: &AfResult,
    rm: &Model,
    cm: &Model,
    pairs: Vec<(usize, usize)>,
    opt: &Options,
) -> Result<Comparison, String> {
    let (rs, cs) = (reference.sequence.as_bytes(), case.sequence.as_bytes());
    let (rnum, cnum) = (numbering(&rm.resnums), numbering(&cm.resnums));
    let p: Vec<V3> = pairs.iter().map(|&(i, _)| rm.ca[i]).collect();
    let q: Vec<V3> = pairs.iter().map(|&(_, j)| cm.ca[j]).collect();
    let identical: Vec<bool> = pairs.iter().map(|&(i, j)| rs[i] == cs[j]).collect();
    let conf: Vec<bool> = pairs
        .iter()
        .map(|&(i, j)| rm.plddt[i] >= opt.min_plddt && cm.plddt[j] >= opt.min_plddt)
        .collect();
    let mut cand: Vec<bool> = identical.iter().zip(&conf).map(|(a, b)| *a && *b).collect();
    let region = opt.fit_region.trim();
    if !region.is_empty() {
        let mask = parse_ranges(region, rs.len())?;
        for (k, &(i, _)) in pairs.iter().enumerate() {
            cand[k] &= mask[i];
        }
        if cand.iter().filter(|x| **x).count() < 10 {
            return Err("Fewer than 10 confident, identical residues lie in the chosen region. Choose a larger region.".into());
        }
    } else if cand.iter().filter(|x| **x).count() < 10 {
        cand = identical.clone();
    }
    // alignment blocks: contiguous in both sequences
    let mut block = vec![0usize; pairs.len()];
    for k in 1..pairs.len() {
        let contiguous = pairs[k].0 == pairs[k - 1].0 + 1 && pairs[k].1 == pairs[k - 1].1 + 1;
        block[k] = block[k - 1] + usize::from(!contiguous);
    }
    let nblocks = block.last().map_or(0, |b| b + 1);
    // seed the fit from the whole set and from every block; keep the seed placing most residues within the cutoff
    let mut best: Option<(usize, M3, V3, Vec<usize>)> = None;
    let mut seeds: Vec<Vec<bool>> = vec![cand.clone()];
    for b in 0..nblocks {
        seeds.push(
            cand.iter()
                .zip(&block)
                .map(|(c, bb)| *c && *bb == b)
                .collect(),
        );
    }
    for seed in seeds {
        if seed.iter().filter(|x| **x).count() < 10 {
            continue;
        }
        let (sp, sq) = pick(&q, &p, &seed);
        let (r0, t0, _) = fit_core(&sp, &sq, opt.cutoff);
        let within: Vec<bool> = (0..q.len())
            .map(|k| cand[k] && dist(&apply(&r0, &t0, &q[k]), &p[k]) <= opt.cutoff)
            .collect();
        if within.iter().filter(|x| **x).count() < 10 {
            continue;
        }
        let (wp, wq) = pick(&q, &p, &within);
        let (r1, t1, core1) = fit_core(&wp, &wq, opt.cutoff);
        let n_ok = core1.iter().filter(|x| **x).count();
        if best.as_ref().is_none_or(|b| n_ok > b.0) {
            let idx: Vec<usize> = (0..within.len()).filter(|&k| within[k]).collect();
            let fit_idx: Vec<usize> = idx
                .iter()
                .zip(&core1)
                .filter(|(_, k)| **k)
                .map(|(i, _)| *i)
                .collect();
            best = Some((n_ok, r1, t1, fit_idx));
        }
    }
    let Some((_, r, t, fit_idx)) = best else {
        return Err("No stable superposition was found (fewer than 10 residues within the tolerance). Try a larger tolerance or choose a region.".into());
    };
    let mut in_fit = vec![false; pairs.len()];
    for k in fit_idx {
        in_fit[k] = true;
    }
    let case_fit: Vec<V3> = cm.ca.iter().map(|x| apply(&r, &t, x)).collect();
    let d: Vec<f64> = pairs
        .iter()
        .map(|&(i, j)| dist(&case_fit[j], &rm.ca[i]))
        .collect::<Vec<_>>();
    let _ = &p;
    let m = pairs.len();
    let mut lddt = vec![f64::NAN; m];
    for a in 0..m {
        let (mut sum, mut cnt) = (0.0, 0usize);
        for b in 0..m {
            if a == b {
                continue;
            }
            let dr = dist(&rm.ca[pairs[a].0], &rm.ca[pairs[b].0]);
            if dr < 15.0 {
                let diff = (dr - dist(&cm.ca[pairs[a].1], &cm.ca[pairs[b].1])).abs();
                sum += [0.5, 1.0, 2.0, 4.0].iter().filter(|t| diff < **t).count() as f64 / 4.0;
                cnt += 1;
            }
        }
        if cnt > 0 {
            lddt[a] = sum / cnt as f64;
        }
    }
    let lddt_ok: Vec<f64> = lddt.iter().cloned().filter(|x| !x.is_nan()).collect();
    let low_conf: Vec<bool> = pairs
        .iter()
        .zip(&d)
        .map(|(&(i, j), &dd)| {
            dd > 5.0 && (rm.plddt[i] < opt.min_plddt || cm.plddt[j] < opt.min_plddt)
        })
        .collect();
    let ss_ref = assign_ss(&rm.ca, &rm.breaks);
    let ss_case = assign_ss(&cm.ca, &cm.breaks);
    let mut ss_changes: Vec<SsChange> = Vec::new();
    let mut ss_changed_residues = 0;
    let mut run: Option<(usize, usize, Ss, Ss, usize)> = None; // start, end (ref idx), from, to, last pair
    for (k, &(i, j)) in pairs.iter().enumerate() {
        let changed = ss_ref[i] != ss_case[j];
        if changed {
            ss_changed_residues += 1;
        }
        match (&mut run, changed) {
            (Some(r), true) if r.2 == ss_ref[i] && r.3 == ss_case[j] && r.4 + 1 == k => {
                r.1 = i;
                r.4 = k;
            }
            _ => {
                if let Some(r) = run.take() {
                    if r.1 + 1 - r.0 >= 3 {
                        ss_changes.push(SsChange {
                            ref_start: rnum[r.0],
                            ref_end: rnum[r.1],
                            from: r.2.letter(),
                            to: r.3.letter(),
                        });
                    }
                }
                if changed {
                    run = Some((i, i, ss_ref[i], ss_case[j], k));
                }
            }
        }
    }
    if let Some(r) = run.take() {
        if r.1 + 1 - r.0 >= 3 {
            ss_changes.push(SsChange {
                ref_start: rnum[r.0],
                ref_end: rnum[r.1],
                from: r.2.letter(),
                to: r.3.letter(),
            });
        }
    }
    let mut blocks = Vec::new();
    for b in 0..nblocks {
        let sel: Vec<usize> = (0..pairs.len()).filter(|&k| block[k] == b).collect();
        if sel.len() < 20 {
            continue;
        }
        let bq: Vec<V3> = sel.iter().map(|&k| q[k]).collect();
        let bp: Vec<V3> = sel.iter().map(|&k| p[k]).collect();
        let (rb, tb) = kabsch(&bq, &bp);
        let moved: Vec<V3> = bq.iter().map(|x| apply(&rb, &tb, x)).collect();
        let glob: Vec<V3> = sel.iter().map(|&k| case_fit[pairs[k].1]).collect();
        blocks.push(Block {
            ref_start: rnum[pairs[sel[0]].0],
            ref_end: rnum[pairs[*sel.last().unwrap()].0],
            case_start: cnum[pairs[sel[0]].1],
            case_end: cnum[pairs[*sel.last().unwrap()].1],
            n: sel.len(),
            rmsd_local_fit: rmsd(&moved, &bp),
            rmsd_global_fit: rmsd(&glob, &bp),
        });
    }
    let d0 = 1.24 * ((rs.len().max(19) - 15) as f64).cbrt() - 1.8;
    let tm = d
        .iter()
        .map(|x| 1.0 / (1.0 + (x / d0).powi(2)))
        .sum::<f64>()
        / rs.len() as f64;
    let fit_p: Vec<V3> = (0..pairs.len())
        .filter(|&k| in_fit[k])
        .map(|k| p[k])
        .collect();
    let fit_q: Vec<V3> = (0..pairs.len())
        .filter(|&k| in_fit[k])
        .map(|k| apply(&r, &t, &q[k]))
        .collect();
    let all_q: Vec<V3> = pairs.iter().map(|&(_, j)| case_fit[j]).collect();
    let n = pairs.len() as f64;
    let stats = Stats {
        ref_length: rs.len(),
        case_length: cs.len(),
        aligned_residues: pairs.len(),
        sequence_identity_of_aligned: identical.iter().filter(|x| **x).count() as f64 / n,
        fit_residues: fit_p.len(),
        fit_rmsd: rmsd(&fit_q, &fit_p),
        rmsd_all_aligned: rmsd(&all_q, &p),
        median_distance: median(&d),
        max_distance: d.iter().cloned().fold(0.0, f64::max),
        fraction_over_5a: d.iter().filter(|x| **x > 5.0).count() as f64 / n,
        fraction_over_10a: d.iter().filter(|x| **x > 10.0).count() as f64 / n,
        tm_score_under_this_fit: tm,
        ref_ranking_score: rm.summary.ranking_score,
        case_ranking_score: cm.summary.ranking_score,
        ref_ptm: rm.summary.ptm,
        case_ptm: cm.summary.ptm,
        ref_mean_plddt: if rm.has_conf {
            mean(&rm.plddt)
        } else {
            f64::NAN
        },
        case_mean_plddt: if cm.has_conf {
            mean(&cm.plddt)
        } else {
            f64::NAN
        },
        ref_model: rm.index,
        case_model: cm.index,
        min_plddt_for_fit: opt.min_plddt,
        fit_cutoff_a: opt.cutoff,
        lddt_ca: lddt_ok.iter().sum::<f64>() / lddt_ok.len().max(1) as f64,
        low_confidence_shifts: low_conf.iter().filter(|x| **x).count(),
        ss_changed_residues,
        chains_ref: n_chains(rm),
        chains_case: n_chains(cm),
        fit_region: if region.is_empty() {
            "automatic (all confident, identical residues)".into()
        } else {
            region.to_string()
        },
    };
    let aligned_ref: std::collections::HashSet<usize> = pairs.iter().map(|x| x.0).collect();
    let aligned_case: std::collections::HashSet<usize> = pairs.iter().map(|x| x.1).collect();
    let missing = ranges(
        (0..rs.len())
            .filter(|i| !aligned_ref.contains(i))
            .map(|i| i + 1)
            .collect(),
    );
    let extra = ranges(
        (0..cs.len())
            .filter(|j| !aligned_case.contains(j))
            .map(|j| j + 1)
            .collect(),
    );
    Ok(Comparison {
        paired_by_number: false,
        ref_num: rnum.clone(),
        case_num: cnum.clone(),
        ref_source: reference.source.clone(),
        case_source: case.source.clone(),
        ref_name: reference.name.clone(),
        case_name: case.name.clone(),
        ref_seq: rs.to_vec(),
        case_seq: cs.to_vec(),
        ref_model: rm.clone(),
        case_model: cm.clone(),
        pairs,
        dist: d,
        in_fit,
        block,
        r,
        t,
        case_fit,
        stats,
        blocks,
        robustness: vec![],
        missing_in_case: missing,
        extra_in_case: extra,
        lddt,
        low_conf,
        ss_ref,
        ss_case,
        ss_changes,
    })
}

pub fn distance_table_csv(c: &Comparison) -> String {
    let mut s = String::from(
        "ref_resi,ref_chain,ref_chain_resi,ref_aa,case_resi,case_chain,case_chain_resi,case_aa,ref_plddt,case_plddt,ca_distance_A,lddt_ca,used_in_fit,alignment_block,ss_ref,ss_case,moved_over_5A_with_low_plddt\n",
    );
    for (k, &(i, j)) in c.pairs.iter().enumerate() {
        let l = if c.lddt[k].is_nan() {
            String::new()
        } else {
            format!("{:.2}", c.lddt[k])
        };
        s.push_str(&format!(
            "{},{},{},{},{},{},{},{},{},{},{:.2},{},{},{},{},{},{}\n",
            i + 1,
            c.ref_model.chain[i],
            c.ref_model.resnums[i],
            c.ref_seq[i] as char,
            j + 1,
            c.case_model.chain[j],
            c.case_model.resnums[j],
            c.case_seq[j] as char,
            if c.ref_model.has_conf {
                format!("{:.1}", c.ref_model.plddt[i])
            } else {
                String::new()
            },
            if c.case_model.has_conf {
                format!("{:.1}", c.case_model.plddt[j])
            } else {
                String::new()
            },
            c.dist[k],
            l,
            u8::from(c.in_fit[k]),
            c.block[k] + 1,
            c.ss_ref[i].letter(),
            c.ss_case[j].letter(),
            u8::from(c.low_conf[k])
        ));
    }
    s
}

/// Reference positions (1-based) whose amino acid is substituted in the case protein (deleted residues are not listed).
pub fn substituted_residues(reference: &str, case: &str) -> Vec<usize> {
    align_sequences(reference.as_bytes(), case.as_bytes(), false)
        .iter()
        .filter(|&&(i, j)| reference.as_bytes()[i] != case.as_bytes()[j])
        .map(|&(i, _)| i + 1)
        .collect()
}

/// One row per residue the user marked (or that a variant changed): how far it moved.
#[derive(Debug, Clone, Serialize)]
pub struct MarkRow {
    pub ref_resi: i32,
    pub ref_aa: char,
    pub case_resi: Option<i32>,
    pub case_aa: Option<char>,
    pub distance: Option<f64>,
    pub lddt_ca: Option<f64>,
    pub ref_plddt: f64,
    pub case_plddt: Option<f64>,
    pub note: String,
}

pub fn marked_rows(c: &Comparison, marks: &[usize]) -> Vec<MarkRow> {
    let by_ref: std::collections::HashMap<usize, usize> = c
        .pairs
        .iter()
        .enumerate()
        .map(|(k, &(i, _))| (i, k))
        .collect();
    marks
        .iter()
        .filter(|&&r| r >= 1 && r <= c.ref_seq.len())
        .map(|&r| {
            let i = r - 1;
            match by_ref.get(&i) {
                Some(&k) => {
                    let j = c.pairs[k].1;
                    let d = c.dist[k];
                    let note = if c.low_conf[k] {
                        "moved, but low confidence in the prediction"
                    } else if d > 5.0 {
                        "moved more than 5 \u{c5}"
                    } else if d > 2.0 {
                        "moved 2-5 \u{c5}"
                    } else {
                        "stays in place (under 2 \u{c5})"
                    };
                    MarkRow {
                        ref_resi: c.ref_num[i],
                        ref_aa: c.ref_seq[i] as char,
                        case_resi: Some(c.case_num[j]),
                        case_aa: Some(c.case_seq[j] as char),
                        distance: Some(d),
                        lddt_ca: Some(c.lddt[k]).filter(|x| !x.is_nan()),
                        ref_plddt: c.ref_model.plddt[i],
                        case_plddt: Some(c.case_model.plddt[j]),
                        note: note.to_string(),
                    }
                }
                None => MarkRow {
                    ref_resi: c.ref_num[i],
                    ref_aa: c.ref_seq[i] as char,
                    case_resi: None,
                    case_aa: None,
                    distance: None,
                    lddt_ca: None,
                    ref_plddt: c.ref_model.plddt[i],
                    case_plddt: None,
                    note: "absent in the case".into(),
                },
            }
        })
        .collect()
}

/// Plain-language summary of a comparison for people who do not read the tables.
pub fn plain_summary(c: &Comparison) -> Vec<String> {
    let s = &c.stats;
    let mut out: Vec<String> = Vec::new();
    if c.ref_source != c.case_source {
        let d = |x: &str| {
            if x == "AlphaFold Server" {
                "an AlphaFold Server result"
            } else {
                "a plain structure file"
            }
        };
        out.push(format!("The reference is {} and the case is {}; part of the difference may come from the prediction method, not from the variant.", d(&c.ref_source), d(&c.case_source)));
    }
    if !c.ref_model.has_conf || !c.case_model.has_conf {
        out.push(
            "At least one file has no confidence values, so every residue is treated as confident."
                .into(),
        );
    }
    if c.paired_by_number {
        out.push("Residues were paired by their numbers in the two files, because the sequence alignment agrees with them; residues missing from one file stay unpaired.".into());
    }
    out.extend(vec![
        format!("{} of {} reference residues have a counterpart in the case ({} residues).", s.aligned_residues, s.ref_length, s.case_length),
        format!("The two structures were lined up on {} residues that match well (RMSD {:.2} \u{c5}).", s.fit_residues, s.fit_rmsd),
        format!(
            "After this alignment the median C\u{3b1} shift of all aligned residues is {:.1} \u{c5} and the largest is {:.1} \u{c5}; {:.0}% of the residues moved more than 5 \u{c5}.",
            s.median_distance, s.max_distance, s.fraction_over_5a * 100.0
        ),
    ]);
    if let (Some(a), Some(b)) = (&c.ref_model.pae, &c.case_model.pae) {
        out.push(format!(
            "Mean predicted aligned error (PAE): {:.1} \u{c5} in the reference and {:.1} \u{c5} in the case; a higher value means the model is less sure about the relative position of its parts.",
            a.mean(),
            b.mean()
        ));
    }
    for &(a, b) in &c.missing_in_case {
        let (a, b) = (c.ref_num[a - 1], c.ref_num[b - 1]);
        out.push(format!(
            "Reference residues {a}-{b} have no counterpart in the case (deleted or cut off)."
        ));
    }
    for &(a, b) in &c.extra_in_case {
        let (a, b) = (c.case_num[a - 1], c.case_num[b - 1]);
        out.push(format!(
            "Case residues {a}-{b} have no counterpart in the reference (new sequence)."
        ));
    }
    if s.low_confidence_shifts > 0 {
        out.push(format!(
            "{} residues moved more than 5 \u{c5} but are poorly predicted (pLDDT below {:.0}) in one of the two models; treat these shifts with caution.",
            s.low_confidence_shifts, s.min_plddt_for_fit
        ));
    }
    let word = |l: char| match l {
        'H' => "helix",
        'E' => "strand",
        _ => "loop",
    };
    for ch in c.ss_changes.iter().take(4) {
        out.push(format!(
            "Reference residues {}-{} are a {} in the reference but a {} in the case.",
            ch.ref_start,
            ch.ref_end,
            word(ch.from),
            word(ch.to)
        ));
    }
    if c.ss_changes.len() > 4 {
        out.push(format!(
            "{} more secondary-structure changes are listed in the residue table.",
            c.ss_changes.len() - 4
        ));
    }
    let selfs: Vec<f64> = c
        .robustness
        .iter()
        .filter(|r| r.comparison.starts_with("reference"))
        .map(|r| r.rmsd_all_aligned)
        .collect();
    if !selfs.is_empty() {
        let lo = selfs.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = selfs.iter().cloned().fold(0.0, f64::max);
        out.push(format!(
            "Noise level: the other AlphaFold models of the reference itself differ from its best model by {lo:.1}-{hi:.1} \u{c5} RMSD. Differences of this size or smaller are not evidence of a real effect."
        ));
    }
    if c.blocks.len() > 1 {
        out.push("The case has more than one rigid part. Each part may fit the reference well on its own while their relative angle differs; see the table of regions and choose which one to line up on.".into());
    }
    out
}
