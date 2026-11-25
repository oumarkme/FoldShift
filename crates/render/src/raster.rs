//! A small z-buffered triangle rasterizer (orthographic), supersampled, with soft shading and outlines.

pub type V3 = [f32; 3];

/// Run `f(first_row, rows)` on bands of whole rows of `data` (`row_len` values per row) in parallel threads.
/// Every row is handled by exactly one call and the bands do not overlap, so the result equals a serial loop.
pub fn par_rows<T: Send>(data: &mut [T], row_len: usize, f: impl Fn(usize, &mut [T]) + Sync) {
    let rows = data.len() / row_len.max(1);
    let threads = std::thread::available_parallelism()
        .map_or(1, |n| n.get())
        .min(8)
        .min(rows.max(1));
    if threads <= 1 || rows < 64 {
        f(0, data);
        return;
    }
    // many small bands handed out one by one: the model fills the middle rows, so equal bands would be unequal work
    let per = rows.div_ceil(threads * 8).max(4);
    let queue = std::sync::Mutex::new(data.chunks_mut(per * row_len).enumerate());
    std::thread::scope(|sc| {
        for _ in 0..threads {
            sc.spawn(|| loop {
                let next = queue.lock().unwrap().next();
                match next {
                    Some((i, band)) => f(i * per, band),
                    None => break,
                }
            });
        }
    });
}

#[derive(Clone)]
pub struct Image {
    pub w: usize,
    pub h: usize,
    pub rgba: Vec<u8>,
}

impl Image {
    pub fn new(w: usize, h: usize, bg: [u8; 4]) -> Image {
        let mut rgba = vec![0u8; w * h * 4];
        for p in rgba.chunks_exact_mut(4) {
            p.copy_from_slice(&bg);
        }
        Image { w, h, rgba }
    }
    pub fn px(&self, x: usize, y: usize) -> [u8; 4] {
        let i = (y * self.w + x) * 4;
        [
            self.rgba[i],
            self.rgba[i + 1],
            self.rgba[i + 2],
            self.rgba[i + 3],
        ]
    }
    pub fn blend(&mut self, x: i64, y: i64, c: [u8; 3], a: f32) {
        if x < 0 || y < 0 || x >= self.w as i64 || y >= self.h as i64 || a <= 0.0 {
            return;
        }
        let i = (y as usize * self.w + x as usize) * 4;
        for k in 0..3 {
            self.rgba[i + k] =
                (c[k] as f32 * a + self.rgba[i + k] as f32 * (1.0 - a)).round() as u8;
        }
        self.rgba[i + 3] = 255;
    }
    /// Box downsample by an integer factor.
    pub fn downsample(&self, f: usize) -> Image {
        let (w, h) = (self.w / f, self.h / f);
        let mut out = Image::new(w, h, [255, 255, 255, 255]);
        par_rows(&mut out.rgba, w * 4, |y0, band| {
            for (r, row) in band.chunks_exact_mut(w * 4).enumerate() {
                let y = y0 + r;
                for x in 0..w {
                    let mut s = [0u32; 3];
                    for dy in 0..f {
                        for dx in 0..f {
                            let p = self.px(x * f + dx, y * f + dy);
                            for k in 0..3 {
                                s[k] += p[k] as u32;
                            }
                        }
                    }
                    let n = (f * f) as u32;
                    row[x * 4] = (s[0] / n) as u8;
                    row[x * 4 + 1] = (s[1] / n) as u8;
                    row[x * 4 + 2] = (s[2] / n) as u8;
                }
            }
        });
        out
    }
}

pub fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub fn mul(a: V3, s: f32) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub fn dot(a: V3, b: V3) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub fn norm(a: V3) -> V3 {
    let l = dot(a, a).sqrt();
    if l < 1e-12 {
        [0.0, 0.0, 1.0]
    } else {
        mul(a, 1.0 / l)
    }
}

#[derive(Clone, Copy)]
pub struct Vert {
    pub p: V3, // x,y in pixels, z = depth (larger = nearer to the viewer)
    pub n: V3, // view-space normal
}

pub struct Frame {
    pub w: usize,
    pub h: usize,
    pub img: Image,
    pub depth: Vec<f32>,
    /// Light direction (towards the light, view space), ambient and diffuse strength.
    pub light: V3,
    pub ambient: f32,
    pub diffuse: f32,
}

const NEAR: f32 = -1e30;

/// One triangle to draw: vertices, per-vertex colors and opacity.
#[derive(Clone, Copy)]
pub struct Tri {
    pub v: [Vert; 3],
    pub col: [[u8; 3]; 3],
    pub alpha: f32,
}

/// Light settings of a frame.
#[derive(Clone, Copy)]
struct Look {
    light: V3,
    ambient: f32,
    diffuse: f32,
}

impl Look {
    fn shade(&self, n: V3, c: [u8; 3]) -> [u8; 3] {
        let n = norm(n);
        // hemispheric ambient (sky from above) plus a soft key light, normals facing the viewer
        let hemi = 0.5 + 0.5 * n[1];
        let key = dot(n, self.light).max(0.0);
        let rim = (1.0 - n[2].abs()).powi(2) * 0.0;
        let l = self.ambient * (0.72 + 0.28 * hemi) + self.diffuse * key + rim;
        [
            (c[0] as f32 * l).min(255.0) as u8,
            (c[1] as f32 * l).min(255.0) as u8,
            (c[2] as f32 * l).min(255.0) as u8,
        ]
    }
}

/// Draw one triangle into the rows that start at row `y0` (`depth` and `rgba` hold only those rows).
fn raster_tri(
    look: &Look,
    w: usize,
    h: usize,
    y0: usize,
    depth: &mut [f32],
    rgba: &mut [u8],
    t: &Tri,
) {
    let (v, col, alpha) = (&t.v, &t.col, t.alpha);
    let rows = depth.len() / w;
    let (a, b, c) = (v[0].p, v[1].p, v[2].p);
    let area = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    if area.abs() < 1e-9 {
        return;
    }
    let minx = a[0].min(b[0]).min(c[0]).floor().max(0.0) as usize;
    let maxx = (a[0].max(b[0]).max(c[0]).ceil().min(w as f32 - 1.0)).max(0.0) as usize;
    let miny = (a[1].min(b[1]).min(c[1]).floor().max(0.0) as usize).max(y0);
    let maxy = (a[1].max(b[1]).max(c[1]).ceil().min(h as f32 - 1.0)).max(0.0) as usize;
    let maxy = maxy.min(y0 + rows - 1);
    if miny > maxy {
        return;
    }
    let inv = 1.0 / area;
    for y in miny..=maxy {
        for x in minx..=maxx {
            let (px, py) = (x as f32 + 0.5, y as f32 + 0.5);
            let w0 = ((b[0] - px) * (c[1] - py) - (b[1] - py) * (c[0] - px)) * inv;
            let w1 = ((c[0] - px) * (a[1] - py) - (c[1] - py) * (a[0] - px)) * inv;
            let w2 = 1.0 - w0 - w1;
            if w0 < 0.0 || w1 < 0.0 || w2 < 0.0 {
                continue;
            }
            let z = w0 * a[2] + w1 * b[2] + w2 * c[2];
            let idx = (y - y0) * w + x;
            if z <= depth[idx] {
                continue;
            }
            let n = [
                w0 * v[0].n[0] + w1 * v[1].n[0] + w2 * v[2].n[0],
                w0 * v[0].n[1] + w1 * v[1].n[1] + w2 * v[2].n[1],
                w0 * v[0].n[2] + w1 * v[1].n[2] + w2 * v[2].n[2],
            ];
            let n = if n[2] < 0.0 { [-n[0], -n[1], -n[2]] } else { n }; // two-sided lighting
            let mix = |k: usize| {
                (w0 * col[0][k] as f32 + w1 * col[1][k] as f32 + w2 * col[2][k] as f32) as u8
            };
            let rgb = look.shade(n, [mix(0), mix(1), mix(2)]);
            let i = idx * 4;
            if alpha >= 1.0 {
                depth[idx] = z;
                rgba[i..i + 3].copy_from_slice(&rgb);
            } else if alpha > 0.0 {
                for k in 0..3 {
                    rgba[i + k] =
                        (rgb[k] as f32 * alpha + rgba[i + k] as f32 * (1.0 - alpha)).round() as u8;
                }
                rgba[i + 3] = 255;
            }
        }
    }
}

impl Frame {
    pub fn new(w: usize, h: usize) -> Frame {
        Frame {
            w,
            h,
            img: Image::new(w, h, [255, 255, 255, 255]),
            depth: vec![NEAR; w * h],
            light: norm([-0.45, 0.55, 0.7]),
            ambient: 0.55,
            diffuse: 0.55,
        }
    }

    fn shade(&self, n: V3, c: [u8; 3]) -> [u8; 3] {
        self.look().shade(n, c)
    }

    fn look(&self) -> Look {
        Look {
            light: self.light,
            ambient: self.ambient,
            diffuse: self.diffuse,
        }
    }

    /// Draw a triangle. `alpha` < 1 blends over what is there (translucent pass: depth test only, no depth write).
    pub fn tri(&mut self, v: [Vert; 3], col: [[u8; 3]; 3], alpha: f32) {
        let (w, h, look) = (self.w, self.h, self.look());
        raster_tri(
            &look,
            w,
            h,
            0,
            &mut self.depth,
            &mut self.img.rgba,
            &Tri { v, col, alpha },
        );
    }

    /// Draw many triangles in the given order. The picture is cut into bands of rows that are drawn on
    /// separate threads; every band sees all triangles in the same order, so each pixel gets exactly the
    /// result of drawing them one after the other.
    pub fn tris(&mut self, list: &[Tri]) {
        let (w, h, look) = (self.w, self.h, self.look());
        let threads = std::thread::available_parallelism()
            .map_or(1, |n| n.get())
            .min(8);
        if threads <= 1 || h < 64 || list.len() < 2000 {
            for t in list {
                raster_tri(&look, w, h, 0, &mut self.depth, &mut self.img.rgba, t);
            }
            return;
        }
        let per = h.div_ceil(threads * 8).max(4);
        let queue = std::sync::Mutex::new(
            self.depth
                .chunks_mut(per * w)
                .zip(self.img.rgba.chunks_mut(per * w * 4))
                .enumerate(),
        );
        std::thread::scope(|sc| {
            for _ in 0..threads {
                sc.spawn(|| loop {
                    let next = queue.lock().unwrap().next();
                    let Some((b, (depth, rgba))) = next else {
                        break;
                    };
                    for t in list {
                        raster_tri(&look, w, h, b * per, depth, rgba, t);
                    }
                });
            }
        });
    }

    /// Shaded sphere (orthographic impostor). `r` in pixels.
    pub fn sphere(&mut self, c: V3, r: f32, col: [u8; 3], alpha: f32) {
        let (x0, x1) = (
            (c[0] - r).floor().max(0.0) as usize,
            (c[0] + r).ceil().min(self.w as f32 - 1.0) as usize,
        );
        let (y0, y1) = (
            (c[1] - r).floor().max(0.0) as usize,
            (c[1] + r).ceil().min(self.h as f32 - 1.0) as usize,
        );
        for y in y0..=y1 {
            for x in x0..=x1 {
                let (dx, dy) = (x as f32 + 0.5 - c[0], y as f32 + 0.5 - c[1]);
                let d2 = dx * dx + dy * dy;
                if d2 > r * r {
                    continue;
                }
                let nz = (r * r - d2).sqrt();
                let z = c[2] + nz;
                let idx = y * self.w + x;
                if z <= self.depth[idx] {
                    continue;
                }
                let rgb = self.shade([dx / r, -dy / r, nz / r], col);
                if alpha >= 1.0 {
                    self.depth[idx] = z;
                    let i = idx * 4;
                    self.img.rgba[i..i + 3].copy_from_slice(&rgb);
                } else {
                    self.img.blend(x as i64, y as i64, rgb, alpha);
                }
            }
        }
    }

    /// Soft ambient occlusion (screen space) and dark outlines at depth discontinuities.
    /// `scale` = pixels per Angstrom. Each pixel depends only on the depth buffer as it was before this
    /// call, so the rows are processed in parallel with the same arithmetic as a serial loop.
    pub fn post(&mut self, scale: f32, outline_px: f32) {
        let (w, h) = (self.w, self.h);
        let src = self.depth.clone();
        let rad = (3.5 * scale).max(2.0);
        let taps: Vec<(f32, f32)> = (0..12)
            .map(|i| {
                let a = i as f32 * 2.399963; // golden angle
                let r = rad * ((i as f32 + 1.0) / 12.0).sqrt();
                (a.cos() * r, a.sin() * r)
            })
            .collect();
        let jump = 1.2 * scale; // depth jump (in pixels = scaled Angstrom) that counts as an edge
        let ow = outline_px.ceil() as i64;
        let out_col = [38u8, 38, 38];
        let (near_dz, far_dz) = (0.4 * scale, 6.0 * scale);
        let ntaps = taps.len() as f32;
        // Sample positions per tap, worked out once with the same float arithmetic as `x + dx`:
        // column index per x and row start per y, or NONE when the sample falls outside the picture.
        const NONE: usize = usize::MAX;
        let cols: Vec<Vec<usize>> = taps
            .iter()
            .map(|&(dx, _)| {
                (0..w)
                    .map(|x| {
                        let sx = x as f32 + dx;
                        if sx < 0.0 || sx >= w as f32 {
                            NONE
                        } else {
                            sx as usize
                        }
                    })
                    .collect()
            })
            .collect();
        let rows_of: Vec<Vec<usize>> = taps
            .iter()
            .map(|&(_, dy)| {
                (0..h)
                    .map(|y| {
                        let sy = y as f32 + dy;
                        if sy < 0.0 || sy >= h as f32 {
                            NONE
                        } else {
                            sy as usize * w
                        }
                    })
                    .collect()
            })
            .collect();
        par_rows(&mut self.img.rgba, w * 4, |y0, band| {
            for (r, row) in band.chunks_exact_mut(w * 4).enumerate() {
                let y = y0 + r;
                for x in 0..w {
                    let z = src[y * w + x];
                    let p = &mut row[x * 4..x * 4 + 3];
                    if z == NEAR {
                        // background pixel next to the model gets the outline (outer silhouette)
                        'o: for dy in -ow..=ow {
                            for dx in -ow..=ow {
                                let (nx, ny) = (x as i64 + dx, y as i64 + dy);
                                if nx >= 0
                                    && ny >= 0
                                    && nx < w as i64
                                    && ny < h as i64
                                    && src[ny as usize * w + nx as usize] != NEAR
                                    && dx * dx + dy * dy <= ow * ow
                                {
                                    p.copy_from_slice(&out_col);
                                    break 'o;
                                }
                            }
                        }
                        continue;
                    }
                    let mut occ = 0.0;
                    for t in 0..taps.len() {
                        let (cx, ry) = (cols[t][x], rows_of[t][y]);
                        if cx == NONE || ry == NONE {
                            continue;
                        }
                        let zs = src[ry + cx];
                        if zs != NEAR {
                            let dz = zs - z; // neighbor nearer to the viewer than this point
                            if dz > near_dz && dz < far_dz {
                                occ += 1.0 - (dz / far_dz);
                            }
                        }
                    }
                    let ao = (1.0 - 0.55 * occ / ntaps).max(0.45);
                    // outline: this pixel is much farther than a neighbor
                    let mut edge = false;
                    for (dx, dy) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)] {
                        let (nx, ny) = (x as i64 + dx * ow, y as i64 + dy * ow);
                        if nx < 0 || ny < 0 || nx >= w as i64 || ny >= h as i64 {
                            continue;
                        }
                        let zn = src[ny as usize * w + nx as usize];
                        if zn - z > jump && zn != NEAR {
                            edge = true;
                        }
                    }
                    if edge {
                        p.copy_from_slice(&out_col);
                    } else {
                        for v in p.iter_mut() {
                            *v = (*v as f32 * ao) as u8;
                        }
                    }
                }
            }
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The banded multi-thread drawing must give exactly the picture of drawing the triangles one by one.
    #[test]
    fn banded_triangles_equal_serial() {
        let mut seed = 12345u64;
        let mut rnd = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((seed >> 33) as f32) / (1u64 << 31) as f32
        };
        let (w, h) = (300usize, 260usize);
        let list: Vec<Tri> = (0..6000)
            .map(|k| {
                let c = [rnd() * w as f32, rnd() * h as f32];
                let vert = |r: &mut dyn FnMut() -> f32| Vert {
                    p: [
                        c[0] + r() * 40.0 - 20.0,
                        c[1] + r() * 40.0 - 20.0,
                        r() * 10.0,
                    ],
                    n: norm([r() - 0.5, r() - 0.5, r() + 0.1]),
                };
                Tri {
                    v: [vert(&mut rnd), vert(&mut rnd), vert(&mut rnd)],
                    col: [[(k % 250) as u8, 90, 200]; 3],
                    alpha: if k % 5 == 0 { 0.4 } else { 1.0 },
                }
            })
            .collect();
        let mut a = Frame::new(w, h);
        let mut b = Frame::new(w, h);
        for t in &list {
            a.tri(t.v, t.col, t.alpha);
        }
        b.tris(&list);
        assert!(a.img.rgba == b.img.rgba && a.depth == b.depth);
        // the same picture after the screen-space pass, which also runs on several threads
        a.post(3.0, 2.0);
        let mut c = Frame::new(w, h);
        c.tris(&list);
        c.post(3.0, 2.0);
        assert!(a.img.rgba == c.img.rgba);
    }
}
