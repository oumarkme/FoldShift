//! Software renderer and figure composer.
#![allow(clippy::needless_range_loop, clippy::chunks_exact_to_as_chunks)]
pub mod figure;
pub mod raster;
pub mod ribbon;
pub mod scene;
pub use figure::{encode_jpeg, figure, profile_figure, PREVIEW_DPI, PRINT_DPI};
pub use raster::Image;
pub use scene::{view_data, Kind, Style};
