//! 2D drawing (text, lines, axes) and the figure layouts. Layout is written in units of 300 dpi pixels
//! (text of 7 pt = 29 units); every drawing call multiplies by the resolution factor of the figure
//! being drawn (`dpi / 300`), so any dpi gives the same layout. Okabe-Ito colors, no title.

use crate::raster::Image;
use crate::scene::*;
use ab_glyph::{point, Font, FontRef, PxScale, ScaleFont};
use psplot_core::compare::Comparison;

/// Resolution of the downloadable figures and of the on-screen preview.
pub const PRINT_DPI: f32 = 600.0;
pub const PREVIEW_DPI: f32 = 300.0;
const DPI: f32 = 300.0;
thread_local! { static KF: std::cell::Cell<f32> = const { std::cell::Cell::new(1.0) }; }
fn kf() -> f32 {
    KF.with(|k| k.get())
}
fn set_dpi(dpi: f32) {
    KF.with(|k| k.set(dpi / DPI));
}
pub fn mm(v: f32) -> usize {
    (v / 25.4 * DPI).round() as usize
}
const PT7: f32 = 7.0 / 72.0 * DPI;
const PT8: f32 = 8.0 / 72.0 * DPI;
const BLACK: [u8; 3] = [0, 0, 0];
const BLUE: [u8; 3] = [0x00, 0x72, 0xB2];

fn font() -> FontRef<'static> {
    FontRef::try_from_slice(include_bytes!("../assets/DejaVuSans.ttf")).expect("embedded font")
}

pub fn text_width(s: &str, px: f32) -> f32 {
    let f = font();
    let sf = f.as_scaled(PxScale::from(px));
    s.chars().map(|c| sf.h_advance(f.glyph_id(c))).sum()
}

/// Draw text with its baseline at `y` and left edge at `x`.
fn text_raw(img: &mut Image, x: f32, y: f32, px: f32, s: &str, col: [u8; 3]) {
    let f = font();
    let sf = f.as_scaled(PxScale::from(px));
    let mut cx = x;
    for ch in s.chars() {
        let g = f
            .glyph_id(ch)
            .with_scale_and_position(PxScale::from(px), point(cx, y));
        if let Some(o) = f.outline_glyph(g) {
            let b = o.px_bounds();
            o.draw(|gx, gy, cov| {
                img.blend(
                    b.min.x as i64 + gx as i64,
                    b.min.y as i64 + gy as i64,
                    col,
                    cov,
                )
            });
        }
        cx += sf.h_advance(f.glyph_id(ch));
    }
}

pub fn text_center(img: &mut Image, cx: f32, y: f32, px: f32, s: &str, col: [u8; 3]) {
    text(img, cx - text_width(s, px) / 2.0, y, px, s, col);
}

/// Text rotated 90 degrees counter-clockwise, centred at (cx, cy).
fn text_rot_raw(img: &mut Image, cx: f32, cy: f32, px: f32, s: &str, col: [u8; 3]) {
    let w = text_width(s, px).ceil() as usize + 4;
    let h = (px * 1.4) as usize;
    let mut tmp = Image::new(w, h, [255, 255, 255, 255]);
    text_raw(&mut tmp, 2.0, px, px, s, [0, 0, 0]);
    for y in 0..h {
        for x in 0..w {
            let a = 1.0 - tmp.px(x, y)[0] as f32 / 255.0;
            if a > 0.0 {
                // (x, y) in text space -> rotated: x' = y, y' = (w - x)
                let dx = cx - h as f32 / 2.0 + y as f32;
                let dy = cy - w as f32 / 2.0 + (w - 1 - x) as f32;
                img.blend(dx as i64, dy as i64, col, a);
            }
        }
    }
}

fn rect_raw(img: &mut Image, x0: f32, y0: f32, x1: f32, y1: f32, col: [u8; 3], alpha: f32) {
    for y in y0.floor() as i64..y1.ceil() as i64 {
        for x in x0.floor() as i64..x1.ceil() as i64 {
            let cx = (x1.min(x as f32 + 1.0) - x0.max(x as f32)).max(0.0);
            let cy = (y1.min(y as f32 + 1.0) - y0.max(y as f32)).max(0.0);
            img.blend(x, y, col, cx * cy * alpha);
        }
    }
}

/// Anti-aliased thick line segment.
fn line_raw(img: &mut Image, a: (f32, f32), b: (f32, f32), w: f32, col: [u8; 3]) {
    let (minx, maxx) = (a.0.min(b.0) - w - 1.0, a.0.max(b.0) + w + 1.0);
    let (miny, maxy) = (a.1.min(b.1) - w - 1.0, a.1.max(b.1) + w + 1.0);
    let (dx, dy) = (b.0 - a.0, b.1 - a.1);
    let l2 = (dx * dx + dy * dy).max(1e-9);
    for y in miny.floor() as i64..=maxy.ceil() as i64 {
        for x in minx.floor() as i64..=maxx.ceil() as i64 {
            let (px, py) = (x as f32 + 0.5, y as f32 + 0.5);
            let t = (((px - a.0) * dx + (py - a.1) * dy) / l2).clamp(0.0, 1.0);
            let (qx, qy) = (a.0 + t * dx, a.1 + t * dy);
            let d = ((px - qx).powi(2) + (py - qy).powi(2)).sqrt();
            img.blend(x, y, col, (w / 2.0 + 0.5 - d).clamp(0.0, 1.0));
        }
    }
}

fn disc_raw(img: &mut Image, c: (f32, f32), r: f32, col: [u8; 3]) {
    for y in (c.1 - r - 1.0).floor() as i64..=(c.1 + r + 1.0).ceil() as i64 {
        for x in (c.0 - r - 1.0).floor() as i64..=(c.0 + r + 1.0).ceil() as i64 {
            let d = ((x as f32 + 0.5 - c.0).powi(2) + (y as f32 + 0.5 - c.1).powi(2)).sqrt();
            img.blend(x, y, col, (r + 0.5 - d).clamp(0.0, 1.0));
        }
    }
}

/// Text with its baseline at (x, y) and size `px`, all in layout units.
pub fn text(img: &mut Image, x: f32, y: f32, px: f32, s: &str, col: [u8; 3]) {
    let k = kf();
    text_raw(img, x * k, y * k, px * k, s, col);
}

/// Text rotated 90 degrees counter-clockwise, centred at (cx, cy).
pub fn text_rot(img: &mut Image, cx: f32, cy: f32, px: f32, s: &str, col: [u8; 3]) {
    let k = kf();
    text_rot_raw(img, cx * k, cy * k, px * k, s, col);
}

pub fn rect(img: &mut Image, x0: f32, y0: f32, x1: f32, y1: f32, col: [u8; 3]) {
    let k = kf();
    rect_raw(img, x0 * k, y0 * k, x1 * k, y1 * k, col, 1.0);
}

fn rect_a(img: &mut Image, x0: f32, y0: f32, x1: f32, y1: f32, col: [u8; 3], alpha: f32) {
    let k = kf();
    rect_raw(img, x0 * k, y0 * k, x1 * k, y1 * k, col, alpha);
}

/// Anti-aliased thick line segment.
pub fn line(img: &mut Image, a: (f32, f32), b: (f32, f32), w: f32, col: [u8; 3]) {
    let k = kf();
    line_raw(img, (a.0 * k, a.1 * k), (b.0 * k, b.1 * k), w * k, col);
}

pub fn disc(img: &mut Image, c: (f32, f32), r: f32, col: [u8; 3]) {
    let k = kf();
    disc_raw(img, (c.0 * k, c.1 * k), r * k, col);
}

pub fn blit(dst: &mut Image, src: &Image, x0: usize, y0: usize) {
    for y in 0..src.h {
        for x in 0..src.w {
            if x0 + x < dst.w && y0 + y < dst.h {
                let p = src.px(x, y);
                let i = ((y0 + y) * dst.w + x0 + x) * 4;
                dst.rgba[i..i + 4].copy_from_slice(&p);
            }
        }
    }
}

fn nice_step(range: f64, target: f64) -> f64 {
    let raw = range / target;
    let mag = 10f64.powf(raw.log10().floor());
    for m in [1.0, 2.0, 2.5, 5.0, 10.0] {
        if raw <= m * mag {
            return m * mag;
        }
    }
    10.0 * mag
}

fn fmt_tick(v: f64, step: f64) -> String {
    if step >= 1.0 {
        format!("{}", v.round() as i64)
    } else {
        format!("{:.*}", if step >= 0.1 { 1 } else { 2 }, v)
    }
}

/// The per-residue distance profile into `img` inside the box (x0,y0,x1,y1) including axis labels.
pub fn profile(img: &mut Image, bx: (usize, usize, usize, usize), c: &Comparison) {
    let (left, bottom, top, right) = (
        mm(15.0) as f32,
        mm(11.0) as f32,
        mm(2.0) as f32,
        mm(2.0) as f32,
    );
    let (x0, x1) = (bx.0 as f32 + left, bx.2 as f32 - right);
    let (y0, y1) = (bx.1 as f32 + top, bx.3 as f32 - bottom);
    let n = c.ref_seq.len() as f64;
    let ymax_raw = c.dist.iter().cloned().fold(1.0, f64::max) * 1.08;
    let ystep = nice_step(ymax_raw, 5.0);
    let ymax = (ymax_raw / ystep).ceil() * ystep;
    // x axis in residue numbers of the file when they increase along the chain (gaps then show as gaps)
    let nums = &c.ref_num;
    let numeric = nums.windows(2).all(|w| w[0] < w[1]);
    let coord = |idx1: f64| -> f64 {
        if numeric {
            nums[idx1 as usize - 1] as f64
        } else {
            idx1
        }
    };
    let lo = coord(1.0) - 0.5;
    let span = coord(n) - coord(1.0) + 1.0;
    let sx = |v: f64| x0 + ((v - lo) / span) as f32 * (x1 - x0);
    let sy = |v: f64| y1 - (v / ymax) as f32 * (y1 - y0);
    let mut last_end = [f32::MIN; 2];
    for &(a, b) in &c.missing_in_case {
        let (xa, xb) = (sx(coord(a as f64) - 0.5), sx(coord(b as f64) + 0.5));
        rect(img, xa, y0, xb, y1, [0xBD, 0xBD, 0xBD]);
        let label = if a == b {
            format!("absent {}", c.ref_num[a - 1])
        } else {
            format!("absent {}-{}", c.ref_num[a - 1], c.ref_num[b - 1])
        };
        let w = text_width(&label, PT7);
        let cx = ((xa + xb) / 2.0).clamp(x0 + w / 2.0, x1 - w / 2.0);
        // two rows so that neighbouring labels do not overlap
        let row = if cx - w / 2.0 > last_end[0] + 8.0 {
            0
        } else {
            1
        };
        last_end[row] = cx + w / 2.0;
        text_center(
            img,
            cx,
            y0 + PT7 * (1.15 + 1.5 * row as f32),
            PT7,
            &label,
            BLACK,
        );
    }
    // curve: do not connect across residues missing in the case
    let pts: Vec<(usize, f32, f32)> = c
        .pairs
        .iter()
        .enumerate()
        .map(|(k, &(i, _))| (i, sx(coord(i as f64 + 1.0)), sy(c.dist[k])))
        .collect();
    for w in pts.windows(2) {
        if w[1].0 == w[0].0 + 1 {
            line(img, (w[0].1, w[0].2), (w[1].1, w[1].2), 3.3, BLUE);
        }
    }
    for (k, p) in pts.iter().enumerate() {
        if c.in_fit[k] {
            disc(img, (p.1, p.2), 3.0, BLACK);
        }
    }
    // axes and ticks
    line(img, (x0, y0), (x0, y1), 2.1, BLACK);
    line(img, (x0, y1), (x1, y1), 2.1, BLACK);
    let mut v = 0.0;
    while v <= ymax + 1e-9 {
        let y = sy(v);
        line(img, (x0 - 12.0, y), (x0, y), 2.1, BLACK);
        let s = fmt_tick(v, ystep);
        text(
            img,
            x0 - 18.0 - text_width(&s, PT7),
            y + PT7 * 0.35,
            PT7,
            &s,
            BLACK,
        );
        v += ystep;
    }
    // ticks at round residue numbers
    let mut ticks: Vec<(f64, String)> = Vec::new();
    let (tlo, thi) = (coord(1.0), coord(n));
    let step = nice_step((thi - tlo + 1.0).max(1.0), 6.0);
    let mut t = (tlo / step).ceil() * step;
    while t <= thi {
        ticks.push((t, format!("{}", t as i64)));
        t += step;
    }
    for (pos, lab) in ticks {
        let x = sx(pos);
        line(img, (x, y1), (x, y1 + 12.0), 2.1, BLACK);
        text_center(img, x, y1 + 18.0 + PT7 * 0.8, PT7, &lab, BLACK);
    }
    text_center(
        img,
        (x0 + x1) / 2.0,
        y1 + 18.0 + PT7 * 1.6 + PT8 * 1.25,
        PT8,
        "Reference residue number",
        BLACK,
    );
    text_rot(
        img,
        bx.0 as f32 + PT8 * 0.8,
        (y0 + y1) / 2.0,
        PT8,
        "C\u{3b1} distance (\u{c5})",
        BLACK,
    );
}

/// 88 mm wide distance profile (single-column figure).
pub fn profile_figure(c: &Comparison, dpi: f32) -> Image {
    set_dpi(dpi);
    let k = kf();
    let (w, h) = (mm(FIG_W_MM), mm(62.0));
    let mut img = Image::new(
        (w as f32 * k).round() as usize,
        (h as f32 * k).round() as usize,
        [255, 255, 255, 255],
    );
    profile(&mut img, (0, 0, w, h - mm(4.0)), c);
    if c.in_fit.iter().any(|&b| b) {
        let y = h as f32 - mm(1.5) as f32;
        disc(&mut img, (mm(18.0) as f32, y - PT7 * 0.35), 5.0, BLACK);
        text(
            &mut img,
            mm(18.0) as f32 + 16.0,
            y,
            PT7,
            "Residues used to line up the two structures",
            BLACK,
        );
    }
    img
}

const FIG_W_MM: f32 = 88.0;

/// A legend entry: a square swatch (or a dot for marker spheres) and a label.
struct Item {
    col: [u8; 3],
    alpha: f32,
    dot: bool,
    label: String,
}

fn item(col: [u8; 3], alpha: f32, dot: bool, label: &str) -> Item {
    Item {
        col,
        alpha,
        dot,
        label: label.to_string(),
    }
}

/// Legend entries in rows that wrap inside [x0, x1]; returns the y below the last row.
fn legend(img: &mut Image, items: &[Item], x0: f32, x1: f32, y: f32) -> f32 {
    let s = PT7;
    let row_h = PT7 * 1.9;
    let (mut x, mut row) = (x0, 0.0f32);
    for it in items {
        let need = s + 10.0 + text_width(&it.label, PT7);
        if x > x0 && x + need > x1 {
            x = x0;
            row += row_h;
        }
        let base = y + row + s;
        if it.dot {
            disc(img, (x + s / 2.0, base - s * 0.4), s * 0.5, it.col);
        } else {
            rect_a(
                img,
                x,
                base - s * 0.85,
                x + s,
                base + s * 0.15,
                it.col,
                it.alpha,
            );
        }
        text(img, x + s + 10.0, base, PT7, &it.label, BLACK);
        x += need + 36.0;
    }
    y + row + row_h
}

/// A viridis colorbar from x0 to x1 with ticks, tick labels and a label; returns the y below it.
fn colorbar(
    img: &mut Image,
    (x0, x1, y): (f32, f32, f32),
    ticks: &[(f64, String)],
    label: &str,
) -> f32 {
    let bh = PT7 * 0.9;
    for xi in x0 as usize..x1 as usize {
        let col = viridis((xi as f64 - x0 as f64) / (x1 - x0) as f64);
        rect(img, xi as f32, y, xi as f32 + 1.0, y + bh, col);
    }
    for (f, lab) in ticks {
        let x = x0 + *f as f32 * (x1 - x0);
        line(img, (x, y + bh), (x, y + bh + 8.0), 2.1, BLACK);
        text_center(img, x, y + bh + 10.0 + PT7 * 0.85, PT7, lab, BLACK);
    }
    let yl = y + bh + 10.0 + PT7 * 2.2;
    text_center(img, (x0 + x1) / 2.0, yl, PT7, label, BLACK);
    yl + PT7 * 1.0
}

fn crop_height(img: &mut Image, h: usize) {
    img.rgba.truncate(img.w * h * 4);
    img.h = h;
}

/// One stand-alone figure, 88 mm wide: the structure view with its own legend or colorbar below.
/// All figures of one run use the same camera, so they show the same orientation.
pub fn figure(c: &Comparison, kind: Kind, st: &Style, dpi: f32) -> Image {
    set_dpi(dpi);
    let k = kf();
    let w = mm(FIG_W_MM);
    let pw = (w as f32 * k).round() as usize;
    let cam = camera(c, pw, st);
    let view = render(c, &cam, kind, st);
    let margin = mm(3.0) as f32;
    let mut key = Image::new(pw, (mm(40.0) as f32 * k) as usize, [255, 255, 255, 255]);
    let (x0, x1) = (margin, w as f32 - margin);
    let mut y = PT7 * 0.6;
    let mut items: Vec<Item> = Vec::new();
    let has_junction = c
        .missing_in_case
        .iter()
        .any(|&(s, e)| s > 1 && e < c.ref_seq.len());
    match kind {
        Kind::Overlay => {
            items.push(item(REF, 1.0, false, "Reference"));
            items.push(item(CASE, 1.0, false, "Case (fitted)"));
            items.push(item(CASE, FAINT, false, "Case (not fitted)"));
            if !c.missing_in_case.is_empty() {
                items.push(item(REF_DARK, 1.0, false, "Absent in case"));
            }
            if has_junction {
                items.push(item(JUNCTION, 1.0, true, "Deletion junction"));
            }
            if !c.extra_in_case.is_empty() {
                items.push(item(EXTRA, FAINT, false, "No counterpart in reference"));
            }
            if !st.marks.is_empty() {
                items.push(item(MARK, 1.0, true, "Marked residue"));
            }
        }
        Kind::Distance => {
            let cap = st.cap;
            let step = if cap <= 12.0 {
                2.0
            } else {
                nice_step(cap, 5.0)
            };
            let mut ticks = Vec::new();
            let mut t = 0.0;
            while t <= cap + 1e-9 {
                ticks.push((t / cap, fmt_tick(t, step)));
                t += step;
            }
            y = colorbar(
                &mut key,
                (x0 + mm(6.0) as f32, x1 - mm(6.0) as f32, y),
                &ticks,
                &format!(
                    "C\u{3b1} distance to reference (\u{c5}; {} and above saturated)",
                    fmt_tick(cap, 1.0)
                ),
            ) + PT7 * 0.8;
            items.push(item(REF, 0.35, false, "Reference"));
            if c.stats.low_confidence_shifts > 0 {
                items.push(item(UNALIGNED, 1.0, false, "Moved, low confidence"));
            }
            if !c.extra_in_case.is_empty() {
                items.push(item(EXTRA, FAINT, false, "No counterpart in reference"));
            }
        }
        Kind::Plddt => {
            let ticks: Vec<(f64, String)> = [0.0, 0.5, 1.0]
                .iter()
                .map(|f| (*f, format!("{:.0}", f * 100.0)))
                .collect();
            y = colorbar(
                &mut key,
                (x0 + mm(6.0) as f32, x1 - mm(6.0) as f32, y),
                &ticks,
                "Case pLDDT (AlphaFold confidence, 0-100)",
            ) + PT7 * 0.8;
            items.push(item(REF, 0.35, false, "Reference"));
        }
    }
    let y_end = legend(&mut key, &items, x0, x1, y);
    crop_height(&mut key, ((y_end + PT7 * 0.4) * k) as usize);
    let mut img = Image::new(pw, view.h + key.h, [255, 255, 255, 255]);
    blit(&mut img, &view, 0, 0);
    blit(&mut img, &key, 0, view.h);
    img
}

pub fn encode_jpeg(img: &Image, dpi: f32) -> Vec<u8> {
    let mut out = Vec::new();
    let mut e = jpeg_encoder::Encoder::new(&mut out, 95);
    e.set_density(jpeg_encoder::Density::Inch {
        x: dpi as u16,
        y: dpi as u16,
    });
    let rgb: Vec<u8> = img
        .rgba
        .chunks_exact(4)
        .flat_map(|p| [p[0], p[1], p[2]])
        .collect();
    e.encode(
        &rgb,
        img.w as u16,
        img.h as u16,
        jpeg_encoder::ColorType::Rgb,
    )
    .expect("jpeg");
    out
}
