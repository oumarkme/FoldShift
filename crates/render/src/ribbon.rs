//! Secondary structure from C-alpha geometry, and a ribbon (cartoon) mesh swept along a spline.

use crate::raster::*;

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Ss {
    Helix,
    Strand,
    Loop,
}

fn d(a: V3, b: V3) -> f32 {
    dot(sub(a, b), sub(a, b)).sqrt()
}

/// Helix/strand/loop from C-alpha distances (a light DSSP-like rule; adequate for drawing).
pub fn assign_ss(ca: &[V3], breaks: &[bool]) -> Vec<Ss> {
    let n = ca.len();
    let mut ss = vec![Ss::Loop; n];
    let ok = |i: usize, k: usize| i + k < n && !(i..i + k).any(|j| breaks[j]);
    let mut helix = vec![false; n];
    for i in 0..n.saturating_sub(4) {
        if ok(i, 4) {
            let (d3, d4) = (d(ca[i], ca[i + 3]), d(ca[i], ca[i + 4]));
            if (4.5..=5.7).contains(&d3)
                && (5.4..=7.2).contains(&d4)
                && (4.5..=5.7).contains(&d(ca[i + 1], ca[i + 4]))
            {
                for j in i..=i + 4 {
                    helix[j] = true;
                }
            }
        }
    }
    let mut strand = vec![false; n];
    for i in 0..n.saturating_sub(2) {
        if ok(i, 2) && d(ca[i], ca[i + 2]) > 6.2 && !helix[i] {
            let mut run = 0;
            while i + run + 2 < n && ok(i + run, 2) && d(ca[i + run], ca[i + run + 2]) > 6.2 {
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
    // drop helix runs shorter than 5 residues and strand runs shorter than 3
    let drop_short = |mask: &mut Vec<bool>, min: usize| {
        let mut i = 0;
        while i < n {
            if mask[i] {
                let s = i;
                while i < n && mask[i] {
                    i += 1;
                }
                if i - s < min {
                    for j in s..i {
                        mask[j] = false;
                    }
                }
            } else {
                i += 1;
            }
        }
    };
    drop_short(&mut helix, 5);
    drop_short(&mut strand, 3);
    for i in 0..n {
        ss[i] = if helix[i] {
            Ss::Helix
        } else if strand[i] {
            Ss::Strand
        } else {
            Ss::Loop
        };
    }
    ss
}

fn catmull(p0: V3, p1: V3, p2: V3, p3: V3, t: f32) -> V3 {
    let t2 = t * t;
    let t3 = t2 * t;
    let f = |a: f32, b: f32, c: f32, e: f32| {
        0.5 * ((2.0 * b)
            + (-a + c) * t
            + (2.0 * a - 5.0 * b + 4.0 * c - e) * t2
            + (-a + 3.0 * b - 3.0 * c + e) * t3)
    };
    [
        f(p0[0], p1[0], p2[0], p3[0]),
        f(p0[1], p1[1], p2[1], p3[1]),
        f(p0[2], p1[2], p2[2], p3[2]),
    ]
}

pub struct Section {
    pub pos: V3,
    pub tan: V3,
    pub nrm: V3,
    pub width: f32,
    pub thick: f32,
    pub res: usize,
    pub frac: f32,
}

/// Ribbon cross-sections along each unbroken segment of the chain.
pub fn sections(ca: &[V3], breaks: &[bool], ss: &[Ss], sub_n: usize) -> Vec<Vec<Section>> {
    let n = ca.len();
    let mut out = Vec::new();
    let mut s = 0;
    while s < n {
        let mut e = s;
        while e + 1 < n && !breaks[e] {
            e += 1;
        }
        if e > s {
            out.push(segment(&ca[s..=e], &ss[s..=e], s, sub_n));
        }
        s = e + 1;
    }
    out
}

fn dims(s: Ss) -> (f32, f32) {
    match s {
        Ss::Helix => (1.7, 0.32),
        Ss::Strand => (1.7, 0.32),
        Ss::Loop => (0.34, 0.34),
    }
}

fn segment(ca: &[V3], ss: &[Ss], off: usize, sub_n: usize) -> Vec<Section> {
    let n = ca.len();
    let at = |i: isize| ca[i.clamp(0, n as isize - 1) as usize];
    // per-residue reference normal (towards the helix axis for helices; plane normal for strands)
    let mut nrm: Vec<V3> = (0..n as isize)
        .map(|i| {
            let (a, b, c) = (at(i - 1), at(i), at(i + 1));
            let curv = sub(add(a, c), mul(b, 2.0));
            if dot(curv, curv) > 1e-6 {
                norm(curv)
            } else {
                [0.0, 0.0, 1.0]
            }
        })
        .collect();
    for i in 0..n {
        if ss[i] == Ss::Strand {
            let (a, b, c) = (at(i as isize - 1), at(i as isize), at(i as isize + 1));
            let pn = cross(sub(b, a), sub(c, b));
            if dot(pn, pn) > 1e-6 {
                nrm[i] = norm(pn);
            }
        }
        if i > 0 && dot(nrm[i], nrm[i - 1]) < 0.0 && ss[i] == Ss::Strand {
            nrm[i] = mul(nrm[i], -1.0);
        }
    }
    // smooth the normals over neighbours
    let sm: Vec<V3> = (0..n as isize)
        .map(|i| {
            let g = |k: isize| nrm[k.clamp(0, n as isize - 1) as usize];
            norm(add(add(g(i - 1), mul(g(i), 2.0)), g(i + 1)))
        })
        .collect();
    let mut out = Vec::new();
    for i in 0..n - 1 {
        for k in 0..sub_n {
            let t = k as f32 / sub_n as f32;
            let pos = catmull(
                at(i as isize - 1),
                at(i as isize),
                at(i as isize + 1),
                at(i as isize + 2),
                t,
            );
            let pos2 = catmull(
                at(i as isize - 1),
                at(i as isize),
                at(i as isize + 1),
                at(i as isize + 2),
                t + 0.02,
            );
            let tan = norm(sub(pos2, pos));
            let nn = norm(add(mul(sm[i], 1.0 - t), mul(sm[i + 1], t)));
            // make the normal perpendicular to the tangent
            let nperp = norm(sub(nn, mul(tan, dot(nn, tan))));
            let (w0, h0) = dims(ss[i]);
            let (w1, h1) = dims(ss[i + 1]);
            // smooth width change across the transition (over the whole residue step)
            let s = t * t * (3.0 - 2.0 * t);
            out.push(Section {
                pos,
                tan,
                nrm: nperp,
                width: w0 + (w1 - w0) * s,
                thick: h0 + (h1 - h0) * s,
                res: off + i,
                frac: t,
            });
        }
    }
    let (w, h) = dims(ss[n - 1]);
    let last = ca[n - 1];
    let tan = norm(sub(last, ca[n - 2]));
    out.push(Section {
        pos: last,
        tan,
        nrm: norm(sub(sm[n - 1], mul(tan, dot(sm[n - 1], tan)))),
        width: w,
        thick: h,
        res: off + n - 1,
        frac: 0.0,
    });
    out
}

/// Sweep an elliptical cross-section along the sections. Calls `emit(verts, residue)` for each triangle.
pub fn sweep(secs: &[Section], sides: usize, mut emit: impl FnMut([V3; 3], [V3; 3], usize)) {
    if secs.len() < 2 {
        return;
    }
    let ring = |s: &Section| -> Vec<(V3, V3)> {
        let bin = norm(cross(s.tan, s.nrm));
        (0..sides)
            .map(|k| {
                let a = k as f32 / sides as f32 * std::f32::consts::TAU;
                let (c, sn) = (a.cos(), a.sin());
                // ellipse with half axes width (along bin) and thick (along nrm)
                let p = add(s.pos, add(mul(bin, s.width * c), mul(s.nrm, s.thick * sn)));
                let n = norm(add(
                    mul(bin, c / s.width.max(1e-3)),
                    mul(s.nrm, sn / s.thick.max(1e-3)),
                ));
                (p, n)
            })
            .collect()
    };
    let mut prev = ring(&secs[0]);
    for w in secs.windows(2) {
        let cur = ring(&w[1]);
        let res = if w[0].frac > 0.5 {
            w[0].res + 1
        } else {
            w[0].res
        };
        for k in 0..sides {
            let k2 = (k + 1) % sides;
            let (a, b, c, e) = (prev[k], prev[k2], cur[k], cur[k2]);
            emit([a.0, b.0, c.0], [a.1, b.1, c.1], res);
            emit([b.0, e.0, c.0], [b.1, e.1, c.1], res);
        }
        prev = cur;
    }
}
