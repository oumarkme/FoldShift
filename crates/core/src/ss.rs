//! Secondary structure (helix, strand, loop) from C-alpha geometry: a light DSSP-like rule, adequate for
//! drawing and for reporting gained or lost helices and strands between two predictions.

use crate::structure::V3;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Ss {
    Helix,
    Strand,
    Loop,
}

impl Ss {
    pub fn letter(self) -> char {
        match self {
            Ss::Helix => 'H',
            Ss::Strand => 'E',
            Ss::Loop => 'L',
        }
    }
}

fn d(a: &V3, b: &V3) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

/// `breaks[i]` is true when there is no bond between residue i and i + 1.
pub fn assign_ss(ca: &[V3], breaks: &[bool]) -> Vec<Ss> {
    let n = ca.len();
    let ok = |i: usize, k: usize| i + k < n && !(i..i + k).any(|j| breaks[j]);
    let mut helix = vec![false; n];
    for i in 0..n.saturating_sub(4) {
        if ok(i, 4) {
            let (d3, d4, d14) = (
                d(&ca[i], &ca[i + 3]),
                d(&ca[i], &ca[i + 4]),
                d(&ca[i + 1], &ca[i + 4]),
            );
            if (4.5..=5.7).contains(&d3) && (5.4..=7.2).contains(&d4) && (4.5..=5.7).contains(&d14)
            {
                for h in helix.iter_mut().take(i + 5).skip(i) {
                    *h = true;
                }
            }
        }
    }
    let mut strand = vec![false; n];
    for i in 0..n.saturating_sub(2) {
        if ok(i, 2) && d(&ca[i], &ca[i + 2]) > 6.2 && !helix[i] {
            let mut run = 0;
            while i + run + 2 < n && ok(i + run, 2) && d(&ca[i + run], &ca[i + run + 2]) > 6.2 {
                run += 1;
            }
            if run >= 3 {
                for j in i..i + run + 2 {
                    if !helix[j] {
                        strand[j] = true;
                    }
                }
            }
        }
    }
    drop_short(&mut helix, 5);
    drop_short(&mut strand, 3);
    (0..n)
        .map(|i| {
            if helix[i] {
                Ss::Helix
            } else if strand[i] {
                Ss::Strand
            } else {
                Ss::Loop
            }
        })
        .collect()
}

fn drop_short(mask: &mut [bool], min: usize) {
    let n = mask.len();
    let mut i = 0;
    while i < n {
        if mask[i] {
            let s = i;
            while i < n && mask[i] {
                i += 1;
            }
            if i - s < min {
                for m in mask.iter_mut().take(i).skip(s) {
                    *m = false;
                }
            }
        } else {
            i += 1;
        }
    }
}
