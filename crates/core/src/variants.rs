//! HGVS c. and p. variants applied to a transcript, giving the mutant protein.

use std::fmt;

#[derive(Debug, Clone)]
pub struct Transcript {
    pub accession: String,
    pub gene: String,
    pub definition: String,
    /// Full transcript (DNA, upper case). Empty when only a protein was given.
    pub mrna: Vec<u8>,
    pub cds_start: usize,
    pub cds_end: usize,
    pub protein: String,
}

#[derive(Debug)]
pub struct VariantError(pub String);
impl fmt::Display for VariantError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for VariantError {}
type Res<T> = Result<T, VariantError>;
fn err<T>(m: impl Into<String>) -> Res<T> {
    Err(VariantError(m.into()))
}

#[derive(Debug, Clone)]
pub struct VariantResult {
    pub variants: Vec<String>,
    pub protein: String,
    pub notes: Vec<String>,
    pub start_loss: bool,
    pub changed: bool,
}

const THREE: [(&str, char); 22] = [
    ("Ala", 'A'),
    ("Arg", 'R'),
    ("Asn", 'N'),
    ("Asp", 'D'),
    ("Cys", 'C'),
    ("Gln", 'Q'),
    ("Glu", 'E'),
    ("Gly", 'G'),
    ("His", 'H'),
    ("Ile", 'I'),
    ("Leu", 'L'),
    ("Lys", 'K'),
    ("Met", 'M'),
    ("Phe", 'F'),
    ("Pro", 'P'),
    ("Ser", 'S'),
    ("Thr", 'T'),
    ("Trp", 'W'),
    ("Tyr", 'Y'),
    ("Val", 'V'),
    ("Ter", '*'),
    ("Sec", 'U'),
];

fn three(c: char) -> &'static str {
    THREE
        .iter()
        .find(|(_, o)| *o == c)
        .map(|(t, _)| *t)
        .unwrap_or("?")
}

/// Translate from `start` until the first stop (excluded). Returns (protein, stop_found).
pub fn translate(seq: &[u8], start: usize) -> (String, bool) {
    const T: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    let idx = |b: u8| match b {
        b'T' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    };
    let mut out = String::new();
    let mut i = start;
    while i + 3 <= seq.len() {
        let aa = match (idx(seq[i]), idx(seq[i + 1]), idx(seq[i + 2])) {
            (Some(a), Some(b), Some(c)) => T[a * 16 + b * 4 + c] as char,
            _ => 'X',
        };
        if aa == '*' {
            return (out, true);
        }
        out.push(aa);
        i += 3;
    }
    (out, false)
}

struct Edit {
    start: usize,
    end: usize,
    repl: Vec<u8>,
    label: String,
}

fn num(s: &str) -> usize {
    s.parse().unwrap_or(0)
}

/// Split "c.123A>G" (optionally with "NM_x.y(GENE):" prefix) into ('c', "123A>G").
fn clean(v: &str) -> Res<(char, String)> {
    let mut v = v.trim();
    if let Some(i) = v.find(':') {
        v = &v[i + 1..];
    }
    let mut ch = v.chars();
    let kind = ch.next().unwrap_or(' ');
    if (kind != 'c' && kind != 'p') || ch.next() != Some('.') {
        return err(format!(
            "'{v}': start with c. or p. (for example c.1A>G or p.Val56_Lys62del)."
        ));
    }
    let body = v[2..]
        .trim()
        .trim_start_matches('(')
        .trim_end_matches(')')
        .to_string();
    Ok((kind, body))
}

/// Parsed c. body: (start, optional end, kind, trailing text).
fn split_c(body: &str) -> Option<(usize, Option<usize>, &'static str, String)> {
    let b = body.as_bytes();
    let mut i = 0;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    if i == 0 {
        return None;
    }
    let a = num(&body[..i]);
    let mut e = None;
    let mut rest = &body[i..];
    if let Some(r) = rest.strip_prefix('_') {
        let j = r.bytes().take_while(|c| c.is_ascii_digit()).count();
        if j == 0 {
            return None;
        }
        e = Some(num(&r[..j]));
        rest = &r[j..];
    }
    let dna = |s: &str| !s.is_empty() && s.bytes().all(|c| matches!(c, b'A' | b'C' | b'G' | b'T'));
    for (kw, kind) in [
        ("delins", "delins"),
        ("del", "del"),
        ("dup", "dup"),
        ("ins", "ins"),
    ] {
        if let Some(t) = rest.strip_prefix(kw) {
            if (kind == "delins" || kind == "ins") && !dna(t) {
                return None;
            }
            if (kind == "del" || kind == "dup") && !t.is_empty() && !dna(t) {
                return None;
            }
            return Some((a, e, kind, t.to_string()));
        }
    }
    let rb = rest.as_bytes();
    if e.is_none() && rb.len() == 3 && rb[1] == b'>' && dna(&rest[0..1]) && dna(&rest[2..3]) {
        return Some((a, None, "sub", rest.to_string()));
    }
    None
}

fn parse_c(body: &str, tx: &Transcript) -> Res<Edit> {
    let bytes = body.as_bytes();
    let has_offset = body.contains('*')
        || body.starts_with('-')
        || bytes
            .windows(2)
            .any(|w| w[0].is_ascii_digit() && (w[1] == b'+' || w[1] == b'-'));
    if has_offset {
        return err(format!(
            "c.{body}: only positions inside the coding sequence are supported (no UTR or intron offsets)."
        ));
    }
    let Some((a, e, kind, tail)) = split_c(body) else {
        return err(format!(
            "c.{body}: not understood. Supported: A>G, del, dup, ins, delins."
        ));
    };
    let cds_len = tx.cds_end - tx.cds_start;
    let idx = |p: usize| -> Res<usize> {
        if p < 1 || p > cds_len {
            return err(format!(
                "c.{body}: position {p} is outside the CDS (1-{cds_len})."
            ));
        }
        Ok(tx.cds_start + p - 1)
    };
    let label = format!("c.{body}");
    let i0 = idx(a)?;
    let i1 = match e {
        Some(e) => idx(e)?,
        None => i0,
    };
    if i1 < i0 {
        return err(format!("c.{body}: the end is before the start."));
    }
    let refseq = &tx.mrna[i0..=i1];
    match kind {
        "sub" => {
            let (r, alt) = (tail.as_bytes()[0], tail.as_bytes()[2]);
            if tx.mrna[i0] != r {
                return err(format!(
                    "c.{body}: reference base at c.{a} is {}, not {}.",
                    tx.mrna[i0] as char, r as char
                ));
            }
            Ok(Edit {
                start: i0,
                end: i0 + 1,
                repl: vec![alt],
                label,
            })
        }
        "del" | "dup" => {
            if !tail.is_empty() && refseq != tail.as_bytes() {
                return err(format!(
                    "c.{body}: reference at that range is {}, not {tail}.",
                    String::from_utf8_lossy(refseq)
                ));
            }
            if kind == "del" {
                Ok(Edit {
                    start: i0,
                    end: i1 + 1,
                    repl: vec![],
                    label,
                })
            } else {
                Ok(Edit {
                    start: i1 + 1,
                    end: i1 + 1,
                    repl: refseq.to_vec(),
                    label,
                })
            }
        }
        "ins" => {
            if e != Some(a + 1) {
                return err(format!(
                    "c.{body}: the insertion flanks must be adjacent positions."
                ));
            }
            Ok(Edit {
                start: i1,
                end: i1,
                repl: tail.into_bytes(),
                label,
            })
        }
        _ => Ok(Edit {
            start: i0,
            end: i1 + 1,
            repl: tail.into_bytes(),
            label,
        }),
    }
}

fn check_overlap(edits: &[Edit]) -> Res<()> {
    // edits sorted by start descending
    for w in edits.windows(2) {
        if w[1].end > w[0].start {
            return err(format!(
                "{} and {} overlap; list one variant per position.",
                w[1].label, w[0].label
            ));
        }
    }
    Ok(())
}

fn apply_c(tx: &Transcript, mut edits: Vec<Edit>) -> Res<VariantResult> {
    edits.sort_by_key(|e| std::cmp::Reverse(e.start));
    check_overlap(&edits)?;
    let mut mrna = tx.mrna.clone();
    for e in &edits {
        mrna.splice(e.start..e.end, e.repl.iter().copied());
    }
    let mut notes = Vec::new();
    let mut start = tx.cds_start;
    let mut start_loss = false;
    if &mrna[start..start + 3] != b"ATG" {
        start_loss = true;
        let k = (start + 3..mrna.len().saturating_sub(2))
            .step_by(3)
            .find(|&j| &mrna[j..j + 3] == b"ATG");
        let Some(k) = k else {
            return err("The start codon is lost and there is no downstream in-frame ATG, so no protein can be predicted.");
        };
        notes.push(format!(
            "Start codon lost: translation is assumed to begin at the first downstream in-frame Met (residue {} of the reference). Check this assumption; the real start site may differ.",
            (k - start) / 3 + 1
        ));
        start = k;
    }
    let (protein, stopped) = translate(&mrna, start);
    if !stopped {
        notes.push("No stop codon was found before the end of the transcript.".into());
    }
    edits.reverse();
    Ok(VariantResult {
        variants: edits.iter().map(|e| e.label.clone()).collect(),
        changed: protein != tx.protein,
        protein,
        notes,
        start_loss,
    })
}

fn aa(tok: &str) -> Res<char> {
    if let Some((_, c)) = THREE.iter().find(|(t, _)| *t == tok) {
        return Ok(*c);
    }
    let mut ch = tok.chars();
    match (ch.next(), ch.next()) {
        (Some(c), None) if c.is_ascii_uppercase() || c == '*' => Ok(c),
        _ => err(format!("Unknown amino acid '{tok}'.")),
    }
}

/// Read one amino acid token at the start of `s` (3-letter, or 1-letter). Returns (char, rest).
fn take_aa(s: &str) -> Option<(char, &str)> {
    if s.len() >= 3 && s.is_char_boundary(3) {
        if let Some((_, c)) = THREE.iter().find(|(t, _)| *t == &s[..3]) {
            return Some((*c, &s[3..]));
        }
    }
    let c = s.chars().next()?;
    if c.is_ascii_uppercase() || c == '*' {
        Some((c, &s[1..]))
    } else {
        None
    }
}

fn aa_seq(s: &str) -> Res<String> {
    let mut out = String::new();
    let mut rest = s;
    while !rest.is_empty() {
        match take_aa(rest) {
            Some((c, r)) => {
                out.push(c);
                rest = r;
            }
            None => return err(format!("Cannot read the amino acid string '{s}'.")),
        }
    }
    Ok(out)
}

fn take_pos(s: &str) -> Option<(char, usize, &str)> {
    let (c, r) = take_aa(s)?;
    let n = r.bytes().take_while(|b| b.is_ascii_digit()).count();
    if n == 0 {
        return None;
    }
    Some((c, num(&r[..n]), &r[n..]))
}

fn apply_p(tx: &Transcript, bodies: &[String]) -> Res<VariantResult> {
    let mut prot: Vec<char> = tx.protein.chars().collect();
    let mut notes = Vec::new();
    let mut edits: Vec<(usize, usize, Vec<char>, String)> = Vec::new();
    let mut start_loss = false;
    for body in bodies {
        let label = format!("p.{body}");
        if body.contains("fs") || body.contains("ext") {
            return err(format!("{label}: a frameshift or extension needs the nucleotide change; enter the c. variant instead."));
        }
        if let Some((_, 1, "?")) = take_pos(body) {
            start_loss = true;
            let Some(k) = prot.iter().skip(1).position(|&c| c == 'M').map(|k| k + 1) else {
                return err("There is no downstream Met after the start loss.");
            };
            notes.push(format!(
                "Start loss: translation is assumed from the first downstream Met (residue {}).",
                k + 1
            ));
            edits.push((0, k, vec![], label));
            continue;
        }
        let Some((a_aa, a, rest)) = take_pos(body) else {
            return err(format!(
                "{label}: not understood. Supported: substitution, del, dup, ins, delins."
            ));
        };
        let (b_aa, b, rest) = if let Some(r) = rest.strip_prefix('_') {
            match take_pos(r) {
                Some((c, n, r2)) => (c, n, r2),
                None => return err(format!("{label}: not understood.")),
            }
        } else {
            (a_aa, a, rest)
        };
        let check = |c: char, pos: usize| -> Res<()> {
            if pos < 1 || pos > prot.len() || prot[pos - 1] != c {
                let have = if pos >= 1 && pos <= prot.len() {
                    three(prot[pos - 1]).to_string()
                } else {
                    "none".into()
                };
                return err(format!(
                    "{label}: residue {pos} of {} is {have}, not {}. Check the isoform numbering.",
                    tx.accession,
                    three(c)
                ));
            }
            Ok(())
        };
        check(a_aa, a)?;
        check(b_aa, b)?;
        let (s, e) = (a - 1, b);
        if rest == "del" {
            edits.push((s, e, vec![], label));
        } else if rest == "dup" {
            edits.push((e, e, prot[s..e].to_vec(), label));
        } else if let Some(t) = rest.strip_prefix("delins") {
            edits.push((s, e, aa_seq(t)?.chars().collect(), label));
        } else if let Some(t) = rest.strip_prefix("ins") {
            if b != a + 1 {
                return err(format!("{label}: the insertion flanks must be adjacent."));
            }
            edits.push((a, a, aa_seq(t)?.chars().collect(), label));
        } else if a == b && !rest.is_empty() {
            let new = aa(rest)?;
            edits.push((s, e, vec![new], label));
        } else {
            return err(format!(
                "{label}: not understood. Supported: substitution, del, dup, ins, delins."
            ));
        }
    }
    edits.sort_by_key(|x| std::cmp::Reverse(x.0));
    for w in edits.windows(2) {
        if w[1].1 > w[0].0 {
            return err(format!("{} and {} overlap.", w[1].3, w[0].3));
        }
    }
    for (s, e, r, _) in &edits {
        prot.splice(*s..*e, r.iter().copied());
    }
    if let Some(i) = prot.iter().position(|&c| c == '*') {
        prot.truncate(i);
        notes.push("A stop codon was introduced; the protein is cut there.".into());
    }
    let protein: String = prot.into_iter().collect();
    let variants = edits.iter().rev().map(|e| e.3.clone()).collect();
    Ok(VariantResult {
        variants,
        changed: protein != tx.protein,
        protein,
        notes,
        start_loss,
    })
}

/// Apply all variants to one allele. c. and p. cannot be mixed.
pub fn apply_variants(tx: &Transcript, variants: &[String]) -> Res<VariantResult> {
    let vs: Vec<&str> = variants
        .iter()
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .collect();
    if vs.is_empty() {
        return Ok(VariantResult {
            variants: vec![],
            protein: tx.protein.clone(),
            notes: vec![],
            start_loss: false,
            changed: false,
        });
    }
    let parsed: Vec<(char, String)> = vs.iter().map(|v| clean(v)).collect::<Res<_>>()?;
    let first = parsed[0].0;
    if parsed.iter().any(|(k, _)| *k != first) {
        return err("Do not mix c. and p. variants on one allele.");
    }
    if first == 'c' {
        if tx.mrna.is_empty() {
            return err("A protein sequence was given, so c. variants cannot be applied. Use p. notation, or paste the DNA sequence.");
        }
        let edits = parsed
            .iter()
            .map(|(_, b)| parse_c(b, tx))
            .collect::<Res<Vec<_>>>()?;
        apply_c(tx, edits)
    } else {
        apply_p(tx, &parsed.into_iter().map(|(_, b)| b).collect::<Vec<_>>())
    }
}
