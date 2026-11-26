//! Builds the two structure views (overlay and distance map) from a Comparison.

use crate::raster::*;
use crate::ribbon::{assign_ss, sections, sweep};
use psplot_core::compare::Comparison;

pub const REF: [u8; 3] = [0xBD, 0xBD, 0xBD];
pub const REF_DARK: [u8; 3] = [0x4D, 0x4D, 0x4D];
pub const CASE: [u8; 3] = [0xE6, 0x9F, 0x00];
pub const EXTRA: [u8; 3] = [0xD5, 0x5E, 0x00];
pub const JUNCTION: [u8; 3] = [0x00, 0x72, 0xB2];
pub const MARK: [u8; 3] = [0x00, 0x9E, 0x73];
pub const UNALIGNED: [u8; 3] = [0x75, 0x75, 0x75];
pub const FAINT: f32 = 0.28;
/// Empty margin around the structure in Angstrom.
const PAD: f32 = 3.0;

const VIRIDIS: [[u8; 3]; 11] = [
    [0x44, 0x01, 0x54],
    [0x48, 0x25, 0x76],
    [0x41, 0x44, 0x87],
    [0x35, 0x60, 0x8D],
    [0x2A, 0x78, 0x8E],
    [0x21, 0x90, 0x8C],
    [0x22, 0xA8, 0x84],
    [0x44, 0xBF, 0x70],
    [0x7A, 0xD1, 0x51],
    [0xBB, 0xDF, 0x27],
    [0xFD, 0xE7, 0x25],
];
pub fn viridis(t: f64) -> [u8; 3] {
    let x = t.clamp(0.0, 1.0) * 10.0;
    let i = (x.floor() as usize).min(9);
    let f = x - i as f64;
    let m =
        |k: usize| (VIRIDIS[i][k] as f64 * (1.0 - f) + VIRIDIS[i + 1][k] as f64 * f).round() as u8;
    [m(0), m(1), m(2)]
}

#[derive(Clone, Copy, PartialEq)]
pub enum Kind {
    Overlay,
    Distance,
    /// Case colored by its pLDDT (confidence), reference faded.
    Plddt,
}

/// Everything the user can tune about a picture; `Style::default()` is the automatic choice.
#[derive(Clone, Debug)]
pub struct Style {
    /// Distance (Angstrom) at which the color scale saturates.
    pub cap: f64,
    /// 0 = automatic, 1 and 2 = turned 90 degrees.
    pub view: u8,
    /// Reference residues (1-based) to mark with green spheres.
    pub marks: Vec<usize>,
    /// Reference residues (1-based) to zoom on; the rest is faded. Empty = whole protein.
    pub focus: Vec<usize>,
    /// Also make the confidence (pLDDT) figure.
    pub plddt_panel: bool,
    /// Rotation about the vertical axis in degrees, on top of `view` (used by the interactive viewer).
    pub yaw: f32,
    pub pitch: f32,
    /// Rotation in the plane of the picture in degrees, counterclockwise, applied after yaw and pitch.
    pub roll: f32,
    /// Magnification about the middle of the structure; 1 shows all of it.
    pub zoom: f32,
    /// Framing of the interactive viewer: a square whose size and middle do not change with the angle.
    pub steady: bool,
}

impl Default for Style {
    fn default() -> Self {
        Style {
            cap: 10.0,
            view: 0,
            marks: vec![],
            focus: vec![],
            plddt_panel: true,
            yaw: 0.0,
            pitch: 0.0,
            roll: 0.0,
            zoom: 1.0,
            steady: false,
        }
    }
}

/// One chain to draw: coordinates, per-residue color and opacity.
struct Chain {
    ca: Vec<V3>,
    col: Vec<[u8; 3]>,
    alpha: Vec<f32>,
    breaks: Vec<bool>,
}

fn v3(a: [f64; 3]) -> V3 {
    [a[0] as f32, a[1] as f32, a[2] as f32]
}

/// Principal axes (columns e1 >= e2 >= e3) of a point set, as a right-handed rotation.
fn principal_axes(pts: &[V3]) -> [V3; 3] {
    let n = pts.len() as f32;
    let c = pts.iter().fold([0.0; 3], |a, p| add(a, mul(*p, 1.0 / n)));
    let mut m = [[0.0f64; 3]; 3];
    for p in pts {
        let d = sub(*p, c);
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] += d[i] as f64 * d[j] as f64;
            }
        }
    }
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..50 {
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if m[p][q].abs() < 1e-12 {
                continue;
            }
            let th = (m[q][q] - m[p][p]) / (2.0 * m[p][q]);
            let t = th.signum() / (th.abs() + (th * th + 1.0).sqrt());
            let t = if th == 0.0 { 1.0 } else { t };
            let cs = 1.0 / (t * t + 1.0).sqrt();
            let sn = t * cs;
            for k in 0..3 {
                let (a, b) = (m[k][p], m[k][q]);
                m[k][p] = cs * a - sn * b;
                m[k][q] = sn * a + cs * b;
            }
            for k in 0..3 {
                let (a, b) = (m[p][k], m[q][k]);
                m[p][k] = cs * a - sn * b;
                m[q][k] = sn * a + cs * b;
            }
            for row in v.iter_mut() {
                let (a, b) = (row[p], row[q]);
                row[p] = cs * a - sn * b;
                row[q] = sn * a + cs * b;
            }
        }
    }
    let mut idx = [0usize, 1, 2];
    idx.sort_by(|&a, &b| m[b][b].partial_cmp(&m[a][a]).unwrap());
    let e = |k: usize| {
        [
            v[0][idx[k]] as f32,
            v[1][idx[k]] as f32,
            v[2][idx[k]] as f32,
        ]
    };
    let (e1, e2) = (e(0), e(1));
    [e1, e2, cross(e1, e2)]
}

pub struct Camera {
    pub axes: [V3; 3],
    pub center: V3,
    /// pixels per Angstrom at the given output width
    pub scale: f32,
    pub w: usize,
    pub h: usize,
    pub minxy: [f32; 2],
}

/// View: 0 = widest face (automatic), 1 = turned 90 degrees about the vertical axis, 2 = turned about the horizontal axis.
pub fn camera(c: &Comparison, width: usize, st: &Style) -> Camera {
    let (mut axes, center, pts) = base_view(c, st);
    if st.yaw != 0.0 || st.pitch != 0.0 {
        let (sy, cy) = st.yaw.to_radians().sin_cos();
        let (sp, cp) = st.pitch.to_radians().sin_cos();
        let (x, y, z) = (axes[0], axes[1], axes[2]);
        let x2 = add(mul(x, cy), mul(z, -sy));
        let z2 = add(mul(x, sy), mul(z, cy));
        let y2 = add(mul(y, cp), mul(z2, sp));
        let z3 = add(mul(z2, cp), mul(y, -sp));
        axes = [x2, y2, z3];
    }
    if st.roll != 0.0 {
        let (s, c) = st.roll.to_radians().sin_cos();
        let (x, y) = (axes[0], axes[1]);
        axes[0] = add(mul(x, c), mul(y, -s));
        axes[1] = add(mul(x, s), mul(y, c));
    }
    if st.steady {
        steady_camera(axes, &pts, width, st.zoom)
    } else {
        fit_camera(axes, center, &pts, width, st.zoom)
    }
}

/// Middle and half side (Angstrom, at zoom 1) of the square viewer frame: a ball that holds every point
/// (Badoiu-Clarkson iteration), plus the padding. It holds the structure at every angle.
fn view_window(pts: &[V3]) -> (V3, f32) {
    let dist2 = |a: V3, b: V3| {
        let d = sub(a, b);
        dot(d, d)
    };
    let mut c = pts
        .iter()
        .fold([0.0; 3], |a, p| add(a, mul(*p, 1.0 / pts.len() as f32)));
    for i in 1..=100 {
        let far = pts
            .iter()
            .fold(c, |f, &p| if dist2(p, c) > dist2(f, c) { p } else { f });
        c = add(c, mul(sub(far, c), 1.0 / (i as f32 + 1.0)));
    }
    let r2 = pts.iter().fold(0.0f32, |m, &p| m.max(dist2(p, c)));
    (c, r2.sqrt() + PAD)
}

/// The square viewer frame: its middle stays at the picture's middle at every angle.
fn steady_camera(axes: [V3; 3], pts: &[V3], width: usize, zoom: f32) -> Camera {
    let (center, half) = view_window(pts);
    let scale = width as f32 / (2.0 * half) * zoom.max(0.1);
    let min = -(width as f32) / (2.0 * scale);
    Camera {
        axes,
        center,
        scale,
        w: width,
        h: width,
        minxy: [min, min],
    }
}

/// Axes (before yaw, pitch and roll), center and the points that set the frame.
fn base_view(c: &Comparison, st: &Style) -> ([V3; 3], V3, Vec<V3>) {
    let in_focus = focus_mask(c, &st.focus);
    let mut pts: Vec<V3> = Vec::new();
    let mut all: Vec<V3> = Vec::new();
    for (i, p) in c.ref_model.ca.iter().enumerate() {
        all.push(v3(*p));
        if in_focus.as_ref().is_none_or(|m| m.0[i]) {
            pts.push(v3(*p));
        }
    }
    for (j, p) in c.case_fit.iter().enumerate() {
        all.push(v3(*p));
        if in_focus.as_ref().is_none_or(|m| m.1[j]) {
            pts.push(v3(*p));
        }
    }
    if pts.is_empty() {
        pts = all.clone();
    }
    // the axes come from the whole protein so a zoom keeps the same orientation
    let mut axes = principal_axes(&all);
    match st.view {
        1 => axes = [axes[2], axes[1], mul(axes[0], -1.0)],
        2 => axes = [axes[0], axes[2], mul(axes[1], -1.0)],
        _ => {}
    }
    let center = pts
        .iter()
        .fold([0.0; 3], |a, p| add(a, mul(*p, 1.0 / pts.len() as f32)));
    (axes, center, pts)
}

/// Frame of a figure: cropped to the structure. With a zoom, it shows what the square viewer frame shows,
/// with the empty margins cut off.
fn fit_camera(axes: [V3; 3], center: V3, pts: &[V3], width: usize, zoom: f32) -> Camera {
    let project = |p: V3| {
        let d = sub(p, center);
        [dot(axes[0], d), -dot(axes[1], d)]
    };
    let (mut lo, mut hi) = ([f32::MAX; 2], [f32::MIN; 2]);
    for &p in pts {
        let q = project(p);
        for k in 0..2 {
            lo[k] = lo[k].min(q[k]);
            hi[k] = hi[k].max(q[k]);
        }
    }
    let (mut bw, mut bh) = (hi[0] - lo[0] + 2.0 * PAD, hi[1] - lo[1] + 2.0 * PAD);
    let mut minxy = [lo[0] - PAD, lo[1] - PAD];
    if zoom != 1.0 {
        let (middle, half) = view_window(pts);
        let (middle, half) = (project(middle), half / zoom.max(0.1));
        let mut size = [bw, bh];
        for k in 0..2 {
            let start = minxy[k].max(middle[k] - half);
            size[k] = (hi[k] + PAD).min(middle[k] + half) - start;
            minxy[k] = start;
        }
        [bw, bh] = size;
    }
    let scale = width as f32 / bw;
    let h = ((bh * scale).round() as usize).clamp(width / 3, width * 3 / 2);
    // fit height too if the clamp shrank it
    let scale = scale.min(h as f32 / bh);
    // a zoomed frame that the height limits made larger than its content keeps the content in the middle
    if zoom != 1.0 {
        minxy[0] -= (width as f32 / scale - bw) / 2.0;
        minxy[1] -= (h as f32 / scale - bh) / 2.0;
    }
    Camera {
        axes,
        center,
        scale,
        w: width,
        h,
        minxy,
    }
}

impl Camera {
    fn project(&self, p: V3, ss: f32) -> V3 {
        let d = sub(p, self.center);
        let s = self.scale * ss;
        let x = (dot(self.axes[0], d) - self.minxy[0]) * s;
        let y = (-dot(self.axes[1], d) - self.minxy[1]) * s;
        [x, y, dot(self.axes[2], d) * s]
    }
    fn dir(&self, n: V3) -> V3 {
        [
            dot(self.axes[0], n),
            dot(self.axes[1], n),
            dot(self.axes[2], n),
        ]
    }
}

/// Which reference and case residues lie in the focus region (None: no focus).
fn focus_mask(c: &Comparison, focus: &[usize]) -> Option<(Vec<bool>, Vec<bool>)> {
    if focus.is_empty() {
        return None;
    }
    let mut r = vec![false; c.ref_seq.len()];
    let mut k = vec![false; c.case_seq.len()];
    for &x in focus {
        if x >= 1 && x <= r.len() {
            r[x - 1] = true;
        }
    }
    for &(i, j) in &c.pairs {
        if r[i] {
            k[j] = true;
        }
    }
    Some((r, k))
}

fn chains(c: &Comparison, kind: Kind, st: &Style) -> (Chain, Chain) {
    let cap = st.cap;
    let ref_ca: Vec<V3> = c.ref_model.ca.iter().map(|p| v3(*p)).collect();
    let case_ca: Vec<V3> = c.case_fit.iter().map(|p| v3(*p)).collect();
    let mut ref_col = vec![REF; ref_ca.len()];
    let mut ref_alpha = vec![1.0f32; ref_ca.len()];
    let mut case_col = vec![CASE; case_ca.len()];
    let mut case_alpha = vec![1.0f32; case_ca.len()];
    let mut aligned_ref = vec![false; ref_ca.len()];
    let mut aligned_case = vec![false; case_ca.len()];
    for (k, &(i, j)) in c.pairs.iter().enumerate() {
        aligned_ref[i] = true;
        aligned_case[j] = true;
        if kind == Kind::Distance {
            case_col[j] = if c.low_conf[k] {
                UNALIGNED
            } else {
                viridis(c.dist[k] / cap)
            };
        }
        if kind == Kind::Plddt {
            case_col[j] = viridis(c.case_model.plddt[j] / 100.0);
        }
        if !c.in_fit[k] && kind == Kind::Overlay {
            case_alpha[j] = FAINT;
        }
    }
    for i in 0..ref_ca.len() {
        if !aligned_ref[i] && kind == Kind::Overlay {
            ref_col[i] = REF_DARK;
        }
        if kind != Kind::Overlay {
            ref_alpha[i] = 0.35;
        }
    }
    for j in 0..case_ca.len() {
        if !aligned_case[j] {
            case_col[j] = EXTRA;
            case_alpha[j] = FAINT;
        }
    }
    if kind == Kind::Plddt {
        for j in 0..case_ca.len() {
            case_col[j] = viridis(c.case_model.plddt[j] / 100.0);
            case_alpha[j] = 1.0;
        }
    }
    if let Some((fr, fc)) = focus_mask(c, &st.focus) {
        for i in 0..ref_alpha.len() {
            if !fr[i] {
                ref_alpha[i] = ref_alpha[i].min(0.18);
            }
        }
        for j in 0..case_alpha.len() {
            if !fc[j] {
                case_alpha[j] = case_alpha[j].min(0.18);
            }
        }
    }
    (
        Chain {
            breaks: c.ref_model.breaks.clone(),
            ca: ref_ca,
            col: ref_col,
            alpha: ref_alpha,
        },
        Chain {
            breaks: c.case_model.breaks.clone(),
            ca: case_ca,
            col: case_col,
            alpha: case_alpha,
        },
    )
}

/// Draw a chain's ribbon into the frame. `translucent` selects which pass.
fn draw_chain(f: &mut Frame, cam: &Camera, ch: &Chain, ss: f32, translucent: bool) {
    let mut list: Vec<Tri> = Vec::new();
    let sst = assign_ss(&ch.ca, &ch.breaks);
    let secs = sections(&ch.ca, &ch.breaks, &sst, 8);
    for seg in secs {
        sweep(&seg, 10, |p, n, res| {
            let a = ch.alpha[res.min(ch.alpha.len() - 1)];
            if (a < 1.0) != translucent {
                return;
            }
            let col = ch.col[res.min(ch.col.len() - 1)];
            let v = |k: usize| Vert {
                p: cam.project(p[k], ss),
                n: cam.dir(n[k]),
            };
            list.push(Tri {
                v: [v(0), v(1), v(2)],
                col: [col; 3],
                alpha: a,
            });
        });
    }
    f.tris(&list);
}

pub struct Marks {
    pub flank_ref: Vec<usize>,
    pub flank_case: Vec<usize>,
    pub mark_ref: Vec<usize>,
    pub mark_case: Vec<usize>,
}

pub fn marks(c: &Comparison, user_marks: &[usize]) -> Marks {
    let mut m = Marks {
        flank_ref: vec![],
        flank_case: vec![],
        mark_ref: vec![],
        mark_case: vec![],
    };
    let map: std::collections::HashMap<usize, usize> = c.pairs.iter().cloned().collect();
    for &(a, b) in &c.missing_in_case {
        if a > 1 && b < c.ref_seq.len() {
            for r in [a - 2, b] {
                m.flank_ref.push(r);
                if let Some(&j) = map.get(&r) {
                    m.flank_case.push(j);
                }
            }
        }
    }
    for &r in user_marks {
        if r >= 1 && r <= c.ref_seq.len() {
            m.mark_ref.push(r - 1);
            if let Some(&j) = map.get(&(r - 1)) {
                m.mark_case.push(j);
            }
        }
    }
    m
}

pub fn render(c: &Comparison, cam: &Camera, kind: Kind, st: &Style) -> Image {
    let user_marks = &st.marks[..];
    // 2x supersampling up to 300 dpi; at higher resolution the pixels are already fine enough
    let ss: usize = if cam.w >= 1500 { 1 } else { 2 };
    let ssf = ss as f32;
    let kk = cam.w as f32 / 1039.0; // outline width follows the resolution
    let mut f = Frame::new(cam.w * ss, cam.h * ss);
    let (rc, cc) = chains(c, kind, st);
    // opaque pass, then outlines and ambient occlusion, then translucent parts
    draw_chain(&mut f, cam, &rc, ssf, false);
    draw_chain(&mut f, cam, &cc, ssf, false);
    f.post(cam.scale * ssf, 1.6 * ssf * kk);
    draw_chain(&mut f, cam, &rc, ssf, true);
    draw_chain(&mut f, cam, &cc, ssf, true);
    if kind == Kind::Overlay {
        let m = marks(c, user_marks);
        let r = 1.15 * cam.scale * ssf;
        let sph = |f: &mut Frame, ca: &Vec<V3>, ids: &[usize], col: [u8; 3]| {
            for &i in ids {
                if i < ca.len() {
                    let p = cam.project(ca[i], ssf);
                    f.sphere([p[0], p[1], p[2] + 1.5 * cam.scale * ssf], r, col, 1.0);
                }
            }
        };
        sph(&mut f, &rc.ca, &m.flank_ref, JUNCTION);
        sph(&mut f, &cc.ca, &m.flank_case, CASE);
        sph(&mut f, &rc.ca, &m.mark_ref, MARK);
        sph(&mut f, &cc.ca, &m.mark_case, MARK);
    }
    if ss == 1 {
        f.img
    } else {
        f.img.downsample(ss)
    }
}

/// What the page needs to turn the structure at full frame rate while the mouse is down:
/// the same axes, frame and colors as the server picture of the viewer (drawn as a thin trace).
pub struct ViewData {
    /// Axes before yaw, pitch and roll.
    pub axes: [V3; 3],
    /// Middle of the square frame, and its half side in Angstrom at zoom 1.
    pub center: V3,
    pub half: f32,
    /// The reference, then the case.
    pub chains: Vec<ViewChain>,
}

pub struct ViewChain {
    pub ca: Vec<V3>,
    pub col: Vec<[u8; 3]>,
    pub alpha: Vec<f32>,
    pub breaks: Vec<bool>,
}

/// The viewer data of one figure kind; `st` sets the starting axes and the focus, not the angle.
pub fn view_data(c: &Comparison, kind: Kind, st: &Style) -> ViewData {
    let (axes, _, pts) = base_view(c, st);
    let (center, half) = view_window(&pts);
    let (rc, cc) = chains(c, kind, st);
    let view_chain = |ch: Chain| ViewChain {
        ca: ch.ca,
        col: ch.col,
        alpha: ch.alpha,
        breaks: ch.breaks,
    };
    ViewData {
        axes,
        center,
        half,
        chains: vec![view_chain(rc), view_chain(cc)],
    }
}
