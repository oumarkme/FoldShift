# Third-party notices

## DejaVu Sans (`crates/render/assets/DejaVuSans.ttf`)

Embedded for figure text. DejaVu fonts are based on Bitstream Vera, copyright Bitstream, Inc.; the DejaVu changes are in the
public domain. The font license allows use and redistribution, including embedding in software, with the copyright notice
kept. Full text: https://dejavu-fonts.github.io/License.html

## Rust crates

The program statically links open-source Rust crates (for example `tiny_http`, `ureq`, `zip`, `serde`, `ab_glyph`,
`jpeg-encoder`). Checked with `cargo metadata` (170 packages including optional ones): all are MIT, Apache-2.0, BSD, ISC, Zlib or similar permissive licenses; two certificate-data crates use CDLA-Permissive-2.0, and one crate offers LGPL only as an extra option beside MIT. `Cargo.lock` fixes the versions.
