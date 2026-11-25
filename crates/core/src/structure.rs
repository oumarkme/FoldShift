//! Reading AlphaFold Server result zips: C-alpha coordinates, pLDDT, confidences, sequence.

use serde::Serialize;
use std::io::{Read, Seek};

pub type V3 = [f64; 3];

#[derive(Debug, Clone, Default, Serialize)]
pub struct Summary {
    pub ranking_score: Option<f64>,
    pub ptm: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct Model {
    pub index: usize,
    /// Predicted aligned error of this model, when the file or an extra PAE file provides it.
    pub pae: Option<std::sync::Arc<crate::pae::Pae>>,
    /// False when the file has no confidence values (experimental structure, or B-factors that are not pLDDT).
    pub has_conf: bool,
    /// The model file is PDB text (otherwise mmCIF).
    pub pdb: bool,
    /// Residue number within its chain.
    pub resnums: Vec<i32>,
    pub chain: Vec<String>,
    pub breaks: Vec<bool>,
    pub ca: Vec<V3>,
    pub plddt: Vec<f64>,
    pub summary: Summary,
    pub cif: String,
}

#[derive(Debug, Clone)]
pub struct AfResult {
    /// "AlphaFold Server" for its result zips, "structure file" for any other .cif/.pdb input.
    pub source: String,
    pub name: String,
    pub sequence: String,
    pub models: Vec<Model>,
}

impl AfResult {
    /// Keep only the residues of the listed chains (author chain letters, for example "A" or "A, B"; empty = all).
    /// Residues of different chains are never joined: a break is set wherever the kept residues are not neighbors.
    pub fn keep_chains(&mut self, spec: &str) -> Result<(), String> {
        let want: Vec<&str> = spec
            .split([',', ';', ' '])
            .map(str::trim)
            .filter(|c| !c.is_empty())
            .collect();
        if want.is_empty() {
            return Ok(());
        }
        let mut seq = Vec::new();
        let first_chain = self.models[0].chain.clone();
        for (k, m) in self.models.iter_mut().enumerate() {
            let idx: Vec<usize> = (0..m.chain.len())
                .filter(|&i| want.contains(&m.chain[i].as_str()))
                .collect();
            if idx.is_empty() {
                let mut have: Vec<&str> = Vec::new();
                for c in &first_chain {
                    if !have.contains(&c.as_str()) {
                        have.push(c);
                    }
                }
                return Err(format!(
                    "Chain {} was not found. This file has the chains: {}.",
                    want.join(", "),
                    have.join(", ")
                ));
            }
            if k == 0 {
                seq = self.sequence.bytes().collect::<Vec<u8>>();
                seq = idx.iter().filter_map(|&i| seq.get(i).copied()).collect();
            }
            let breaks: Vec<bool> = (0..idx.len())
                .map(|t| t + 1 < idx.len() && (idx[t + 1] != idx[t] + 1 || m.breaks[idx[t]]))
                .collect();
            m.resnums = idx.iter().map(|&i| m.resnums[i]).collect();
            m.chain = idx.iter().map(|&i| m.chain[i].clone()).collect();
            m.ca = idx.iter().map(|&i| m.ca[i]).collect();
            m.plddt = idx.iter().map(|&i| m.plddt[i]).collect();
            m.breaks = breaks;
            if let Some(p) = &m.pae {
                if p.n == first_chain.len() {
                    m.pae = Some(std::sync::Arc::new(p.subset(&idx)));
                } else {
                    m.pae = None;
                }
            }
        }
        self.sequence = String::from_utf8(seq).unwrap_or_default();
        Ok(())
    }

    /// Attach a PAE file (json, npy, npz) to the best model; its size must equal the residue count.
    pub fn attach_pae(&mut self, bytes: &[u8]) -> Result<(), String> {
        let p = crate::pae::parse_pae(bytes)?;
        let b = self.best();
        let m = &mut self.models[b];
        if p.n != m.ca.len() {
            return Err(format!(
                "The PAE matrix is {n} x {n} but the structure has {r} residues.",
                n = p.n,
                r = m.ca.len()
            ));
        }
        m.pae = Some(std::sync::Arc::new(p));
        Ok(())
    }

    /// Treat every residue as fully confident (files without confidence values).
    pub fn drop_confidence(&mut self) {
        for m in &mut self.models {
            m.has_conf = false;
            m.plddt.iter_mut().for_each(|x| *x = 100.0);
        }
    }

    /// Index (into `models`) of the best-ranked model: highest ranking_score, ties and missing scores -> lowest index.
    pub fn best(&self) -> usize {
        let mut best = 0;
        for (i, m) in self.models.iter().enumerate() {
            let s = m.summary.ranking_score.unwrap_or(f64::NEG_INFINITY);
            let b = self.models[best]
                .summary
                .ranking_score
                .unwrap_or(f64::NEG_INFINITY);
            if s > b {
                best = i;
            }
        }
        best
    }
}

fn aa1(r: &str) -> char {
    match r {
        "ALA" => 'A',
        "ARG" => 'R',
        "ASN" => 'N',
        "ASP" => 'D',
        "CYS" => 'C',
        "GLN" => 'Q',
        "GLU" => 'E',
        "GLY" => 'G',
        "HIS" => 'H',
        "ILE" => 'I',
        "LEU" => 'L',
        "LYS" => 'K',
        "MET" => 'M',
        "PHE" => 'F',
        "PRO" => 'P',
        "SER" => 'S',
        "THR" => 'T',
        "TRP" => 'W',
        "TYR" => 'Y',
        "VAL" => 'V',
        "MSE" => 'M',
        _ => 'X',
    }
}

struct AtomCols {
    atom: usize,
    comp: usize,
    asym: usize,
    seq: usize,
    auth_asym: Option<usize>,
    auth_seq: Option<usize>,
    alt: Option<usize>,
    model: Option<usize>,
    x: usize,
    b: usize,
    n: usize,
    first_row: usize,
}

fn atom_cols(lines: &[&str]) -> Option<AtomCols> {
    let start = lines
        .iter()
        .position(|l| l.trim_end() == "_atom_site.group_PDB")?;
    let mut names = Vec::new();
    let mut i = start;
    while i < lines.len() && lines[i].starts_with("_atom_site.") {
        names.push(lines[i]["_atom_site.".len()..].trim());
        i += 1;
    }
    let f = |n: &str| names.iter().position(|x| *x == n);
    // some tools write only the author columns
    let g = |a: &str, b: &str| f(a).or_else(|| f(b));
    Some(AtomCols {
        atom: g("label_atom_id", "auth_atom_id")?,
        comp: g("label_comp_id", "auth_comp_id")?,
        asym: g("label_asym_id", "auth_asym_id")?,
        seq: g("label_seq_id", "auth_seq_id")?,
        auth_asym: f("auth_asym_id"),
        auth_seq: f("auth_seq_id"),
        alt: f("label_alt_id"),
        model: f("pdbx_PDB_model_num"),
        x: f("Cartn_x")?,
        b: f("B_iso_or_equiv")?,
        n: names.len(),
        first_row: i,
    })
}

/// (residue names, resnums, CA coords, pLDDT) of the first chain.
/// C-alpha data of all protein chains of one model, in file order.
pub struct CifData {
    pub aa: Vec<char>,
    pub resnums: Vec<i32>,
    pub chain: Vec<String>,
    pub ca: Vec<V3>,
    pub plddt: Vec<f64>,
    /// `breaks[i]`: no bond between residue i and i + 1 (chain change, numbering gap, or too far apart).
    pub breaks: Vec<bool>,
}

pub fn parse_cif(text: &str) -> Result<CifData, String> {
    let lines: Vec<&str> = text.lines().collect();
    let c = atom_cols(&lines).ok_or("The model file has no atom list.")?;
    let mut d = CifData {
        aa: vec![],
        resnums: vec![],
        chain: vec![],
        ca: vec![],
        plddt: vec![],
        breaks: vec![],
    };
    let mut lseq: Vec<i32> = Vec::new();
    for l in &lines[c.first_row..] {
        if l.starts_with('#') || l.starts_with("loop_") || l.starts_with('_') {
            break;
        }
        let f: Vec<&str> = l.split_whitespace().collect();
        if f.len() < c.n || f[c.atom] != "CA" {
            continue;
        }
        // first model only; the first alternative location only
        if c.model.is_some_and(|m| f[m] != "1") {
            continue;
        }
        if c.alt
            .is_some_and(|a| !matches!(f[a], "." | "?" | "A" | "1"))
        {
            continue;
        }
        let p = |k: usize| {
            f[k].parse::<f64>()
                .map_err(|_| "A coordinate could not be read.".to_string())
        };
        let aa = aa1(f[c.comp]);
        if aa == 'X' && f[c.comp] != "UNK" {
            continue; // ions, ligands and nucleotides are not drawn
        }
        let label_seq = f[c.seq]
            .parse::<i32>()
            .map_err(|_| "A residue number could not be read.".to_string())?;
        // the author's chain letter and residue number are what a person knows (F508, chain A)
        let auth_num = c.auth_seq.and_then(|k| f[k].parse::<i32>().ok());
        d.aa.push(aa);
        d.resnums.push(auth_num.unwrap_or(label_seq));
        lseq.push(label_seq);
        d.chain.push(f[c.auth_asym.unwrap_or(c.asym)].to_string());
        d.ca.push([p(c.x)?, p(c.x + 1)?, p(c.x + 2)?]);
        d.plddt.push(p(c.b)?);
    }
    if d.ca.is_empty() {
        return Err("No C-alpha atoms were found in the model.".into());
    }
    let n = d.ca.len();
    d.breaks = (0..n)
        .map(|i| {
            i + 1 < n
                && (d.chain[i] != d.chain[i + 1]
                    || lseq[i + 1] != lseq[i] + 1
                    || crate::compare::dist(&d.ca[i], &d.ca[i + 1]) > 4.6)
        })
        .collect();
    Ok(d)
}

/// C-alpha data of the first model of a PDB file (fixed columns). B-factor = confidence column.
pub fn parse_pdb(text: &str) -> Result<CifData, String> {
    let mut d = CifData {
        aa: vec![],
        resnums: vec![],
        chain: vec![],
        ca: vec![],
        plddt: vec![],
        breaks: vec![],
    };
    let mut last: Option<(String, i32, char)> = None;
    for l in text.lines() {
        if l.starts_with("ENDMDL") {
            break;
        }
        if !(l.starts_with("ATOM  ") || l.starts_with("HETATM"))
            || l.len() < 54
            || !l.is_char_boundary(54)
        {
            continue;
        }
        let g = |a: usize, b: usize| l.get(a..b.min(l.len())).unwrap_or("").trim();
        if g(12, 16) != "CA" {
            continue;
        }
        let alt = l.as_bytes()[16] as char;
        if alt != ' ' && alt != 'A' && alt != '1' {
            continue;
        }
        let aa = aa1(g(17, 20));
        if aa == 'X' && g(17, 20) != "UNK" {
            continue;
        }
        let chain = g(21, 22).to_string();
        let seq = g(22, 26)
            .parse::<i32>()
            .map_err(|_| "A residue number could not be read.".to_string())?;
        let icode = l.as_bytes().get(26).map_or(' ', |b| *b as char);
        let key = (chain.clone(), seq, icode);
        if last.as_ref() == Some(&key) {
            continue; // a second C-alpha of the same residue
        }
        // an insertion code repeats the residue number: number such residues on
        let mut num = seq;
        if let (Some(p), Some(&prev)) = (&last, d.resnums.last()) {
            if p.0 == chain && num <= prev && (p.1 == seq || icode != ' ') {
                num = prev + 1;
            }
        }
        last = Some(key);
        let p = |a: usize, b: usize| {
            g(a, b)
                .parse::<f64>()
                .map_err(|_| "A coordinate could not be read.".to_string())
        };
        d.aa.push(aa);
        d.resnums.push(num);
        d.chain.push(chain);
        d.ca.push([p(30, 38)?, p(38, 46)?, p(46, 54)?]);
        d.plddt.push(g(60, 66).parse::<f64>().unwrap_or(0.0));
    }
    if d.ca.is_empty() {
        return Err("No C-alpha atoms were found in the model.".into());
    }
    let n = d.ca.len();
    d.breaks = (0..n)
        .map(|i| {
            i + 1 < n
                && (d.chain[i] != d.chain[i + 1]
                    || d.resnums[i + 1] != d.resnums[i] + 1
                    || crate::compare::dist(&d.ca[i], &d.ca[i + 1]) > 4.6)
        })
        .collect();
    Ok(d)
}

/// True when the file says the coordinates come from an experiment (X-ray, EM, NMR): PDB `EXPDTA` line or
/// mmCIF `_exptl.method`. Predicted models carry no such line or say "THEORETICAL MODEL". The B-factor of an
/// experimental structure measures motion and is high where the model is uncertain, so it must never be read as pLDDT.
pub fn is_experimental(text: &str) -> bool {
    for l in text.lines().take(4000) {
        let up = l.trim_start().to_ascii_uppercase();
        if up.starts_with("EXPDTA") || up.starts_with("_EXPTL.METHOD") {
            return !up.contains("THEORETICAL") && !up.contains("PREDICTED");
        }
    }
    false
}

/// Confidence of a model file: the B-factor column read as pLDDT, except for experimental structures (none).
fn confidence(text: &str, plddt: &mut [f64]) -> bool {
    if is_experimental(text) {
        plddt.iter_mut().for_each(|x| *x = 100.0);
        return false;
    }
    normalize_conf(plddt)
}

/// Bring the B-factor column to pLDDT units (0-100). Returns false when it cannot be a confidence
/// (all equal, or above 100 as in many experimental structures); the values are then set to 100.
fn normalize_conf(v: &mut [f64]) -> bool {
    let (lo, hi) = v
        .iter()
        .fold((f64::MAX, f64::MIN), |a, &x| (a.0.min(x), a.1.max(x)));
    if hi <= 1.0 + 1e-9 && hi > 0.0 && lo != hi {
        v.iter_mut().for_each(|x| *x *= 100.0); // 0-1 scale
        return true;
    }
    if hi > 100.0 + 1e-6 || lo == hi || lo < 0.0 {
        v.iter_mut().for_each(|x| *x = 100.0);
        return false;
    }
    true
}

/// PDB text with every ATOM and HETATM position replaced by R x + t.
pub fn transform_pdb(text: &str, r: &[[f64; 3]; 3], t: &V3) -> String {
    let mut out = String::with_capacity(text.len());
    for l in text.lines() {
        if (l.starts_with("ATOM  ") || l.starts_with("HETATM"))
            && l.len() >= 54
            && l.is_char_boundary(54)
        {
            let q = |a: usize, b: usize| l[a..b].trim().parse::<f64>().unwrap_or(0.0);
            let p = [q(30, 38), q(38, 46), q(46, 54)];
            let n: Vec<f64> = (0..3)
                .map(|k| r[k][0] * p[0] + r[k][1] * p[1] + r[k][2] * p[2] + t[k])
                .collect();
            out.push_str(&format!(
                "{}{:8.3}{:8.3}{:8.3}{}\n",
                &l[..30],
                n[0],
                n[1],
                n[2],
                &l[54..]
            ));
        } else {
            out.push_str(l);
            out.push('\n');
        }
    }
    out
}

/// Rewrite atom coordinates with x' = R x + t (all atoms), returning a new mmCIF text.
pub fn transform_cif(text: &str, r: &[[f64; 3]; 3], t: &V3) -> String {
    let lines: Vec<&str> = text.lines().collect();
    let Some(c) = atom_cols(&lines) else {
        return text.to_string();
    };
    let mut out = String::with_capacity(text.len());
    for (i, l) in lines.iter().enumerate() {
        if i >= c.first_row && (l.starts_with("ATOM") || l.starts_with("HETATM")) {
            let mut f: Vec<String> = l.split_whitespace().map(String::from).collect();
            if f.len() >= c.n {
                let p = [
                    f[c.x].parse::<f64>().unwrap_or(0.0),
                    f[c.x + 1].parse().unwrap_or(0.0),
                    f[c.x + 2].parse().unwrap_or(0.0),
                ];
                for k in 0..3 {
                    f[c.x + k] = format!(
                        "{:.3}",
                        r[k][0] * p[0] + r[k][1] * p[1] + r[k][2] * p[2] + t[k]
                    );
                }
                out.push_str(&f.join(" "));
                out.push('\n');
                continue;
            }
        }
        out.push_str(l);
        out.push('\n');
    }
    out
}

pub fn load_zip<R: Read + Seek>(src: R) -> Result<AfResult, String> {
    let mut z =
        zip::ZipArchive::new(src).map_err(|e| format!("This is not a readable zip file ({e})."))?;
    let names: Vec<String> = (0..z.len())
        .filter_map(|i| z.by_index(i).ok().map(|f| f.name().to_string()))
        .collect();
    let model_idx = |n: &str| -> Option<usize> {
        let stem = n.strip_suffix(".cif")?;
        stem.rsplit_once("_model_")?.1.parse().ok()
    };
    let mut cifs: Vec<(usize, String)> = names
        .iter()
        .filter(|n| !n.contains("template"))
        .filter_map(|n| model_idx(n).map(|i| (i, n.clone())))
        .collect();
    cifs.sort();
    let mut source = "AlphaFold Server";
    if cifs.is_empty() {
        // any other structure file(s): .cif, .mmcif, .pdb, .ent; "rank_N" in the name orders the models
        source = "structure file";
        let is_model = |n: &str| {
            let l = n.to_ascii_lowercase();
            !l.contains("__macosx")
                && !l.rsplit('/').next().unwrap_or("").starts_with('.')
                && [".cif", ".mmcif", ".pdb", ".ent"]
                    .iter()
                    .any(|e| l.ends_with(e))
        };
        let rank = |n: &str| -> usize {
            let l = n.to_ascii_lowercase();
            l.find("rank_")
                .map(|k| {
                    l[k + 5..]
                        .chars()
                        .take_while(char::is_ascii_digit)
                        .collect::<String>()
                })
                .and_then(|d| d.parse::<usize>().ok())
                .unwrap_or(usize::MAX)
        };
        let mut v: Vec<(usize, String)> = names
            .iter()
            .filter(|n| is_model(n))
            .map(|n| (rank(n), n.clone()))
            .collect();
        v.sort();
        cifs = v
            .into_iter()
            .enumerate()
            .map(|(k, (_, n))| (k + 1, n))
            .collect();
    }
    if cifs.is_empty() {
        return Err(
            "No model file found. Use an AlphaFold Server result zip, or a .cif or .pdb file."
                .into(),
        );
    }
    let mut read = |name: &str| -> Option<String> {
        let mut f = z.by_name(name).ok()?;
        let mut s = String::new();
        f.read_to_string(&mut s).ok()?;
        Some(s)
    };
    let file0 = cifs[0].1.rsplit('/').next().unwrap_or("").to_string();
    let base = match file0.rsplit_once("_model_") {
        Some(x) if source == "AlphaFold Server" => x.0.to_string(),
        _ => file0
            .rsplit_once('.')
            .map_or(file0.clone(), |x| x.0.to_string()),
    };
    let mut job_seq = String::new();
    if let Some(jr) = names.iter().find(|n| n.ends_with("_job_request.json")) {
        if let Some(t) = read(jr) {
            if let Ok(v) = serde_json::from_str::<serde_json::Value>(&t) {
                if let Some(s) = v[0]["sequences"].as_array().and_then(|a| {
                    a.iter()
                        .find_map(|s| s["proteinChain"]["sequence"].as_str())
                }) {
                    job_seq = s.to_string();
                }
            }
        }
    }
    let mut models = Vec::new();
    let mut cif_seq = String::new();
    let mut multi = false;
    for (i, n) in &cifs {
        let text = read(n).ok_or("A model file could not be read.")?;
        let pdb = {
            let l = n.to_ascii_lowercase();
            l.ends_with(".pdb") || l.ends_with(".ent")
        };
        let mut d = if pdb {
            parse_pdb(&text)?
        } else {
            parse_cif(&text)?
        };
        let has_conf = confidence(&text, &mut d.plddt);
        if cif_seq.is_empty() {
            cif_seq = d.aa.iter().collect();
            multi = d.chain.iter().any(|c| *c != d.chain[0]);
        }
        let sm = n.replace(
            &format!("_model_{i}.cif"),
            &format!("_summary_confidences_{i}.json"),
        );
        let summary = read(&sm)
            .and_then(|t| serde_json::from_str::<serde_json::Value>(&t).ok())
            .map(|v| Summary {
                ranking_score: v["ranking_score"].as_f64(),
                ptm: v["ptm"].as_f64(),
            })
            .unwrap_or_default();
        models.push(Model {
            index: *i,
            pae: None,
            has_conf,
            pdb,
            resnums: d.resnums,
            chain: d.chain,
            breaks: d.breaks,
            ca: d.ca,
            plddt: d.plddt,
            summary,
            cif: text,
        });
    }
    // one chain: the sequence submitted to AlphaFold Server; several chains: the residues of the model, chain after chain
    let sequence = if job_seq.is_empty() || multi || job_seq.len() != cif_seq.len() {
        cif_seq
    } else {
        job_seq
    };
    let mut res = AfResult {
        source: source.to_string(),
        name: base,
        sequence,
        models,
    };
    if source == "AlphaFold Server" {
        let b = res.best();
        let idx = res.models[b].index;
        let fname = cifs.iter().find(|(i, _)| *i == idx).map(|x| {
            x.1.replace(
                &format!("_model_{idx}.cif"),
                &format!("_full_data_{idx}.json"),
            )
        });
        if let Some(t) = fname.and_then(|f| read(&f)) {
            if let Ok(p) = crate::pae::parse_pae(t.as_bytes()) {
                if p.n == res.models[b].ca.len() {
                    res.models[b].pae = Some(std::sync::Arc::new(p));
                }
            }
        }
    }
    Ok(res)
}

/// Load an uploaded file: an AlphaFold Server zip, or a single .cif or .pdb model (detected from the content),
/// with an optional PAE file (.json, .npy, .npz) for it. `name` is the original file name.
pub fn load_any(bytes: &[u8], name: &str, pae: Option<&[u8]>) -> Result<AfResult, String> {
    let mut res = if bytes.starts_with(b"PK") {
        load_zip(std::io::Cursor::new(bytes))?
    } else {
        let text = String::from_utf8(bytes.to_vec())
            .map_err(|_| "This file is not a text structure file (.cif or .pdb).".to_string())?;
        let pdb = !text.trim_start().starts_with("data_") && text.contains("\nATOM  ")
            || text.starts_with("ATOM  ")
            || text.starts_with("HEADER");
        let mut d = if pdb {
            parse_pdb(&text)?
        } else {
            parse_cif(&text)?
        };
        let has_conf = confidence(&text, &mut d.plddt);
        let sequence: String = d.aa.iter().collect();
        let stem = name.rsplit(['/', '\\']).next().unwrap_or(name);
        let base = stem.rsplit_once('.').map_or(stem, |x| x.0).to_string();
        AfResult {
            source: "structure file".into(),
            name: base,
            sequence,
            models: vec![Model {
                index: 1,
                pae: None,
                has_conf,
                pdb,
                resnums: d.resnums,
                chain: d.chain,
                breaks: d.breaks,
                ca: d.ca,
                plddt: d.plddt,
                summary: Summary::default(),
                cif: text,
            }],
        }
    };
    if let Some(p) = pae {
        res.attach_pae(p)?;
    }
    Ok(res)
}
