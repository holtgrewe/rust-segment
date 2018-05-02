use super::haarseglib::*;

use std::ops::Range;

/// The segment with corresponding value.
pub struct HaarSegment {
    /// Range in input vector of segment.
    pub segment: Range<usize>,
    /// Value of the signal segment.
    pub value: f64,
}

/// Segmentation result.
pub struct HaarSegResult {
    /// Resulting segmented signal.
    pub segments: Vec<HaarSegment>,
    /// Segmented values, same length as `intensities` from input.
    pub seg_values: Vec<f64>,
}

/// Rust interface to the original `haarSeq()` implementation from R.
///
/// This is a reimplementation of the R part of `haarSeg()` using the original
/// C implementation of the Haar Wavelet bundled with the R package.
///
/// # Args
///
/// - `intensities` - corresponds to `I`
/// - `weights` - corresponds to `W`
/// - `raw_intensities` - corresponds to `rI`
/// - `chrom_pos` - corresponds to `chromPos`
/// - `breaks_fdr_q` - corresponds to `breaksFdrQ`
/// - `haar_start_level` - corresponds to `haarStartLevel`
/// - `haar_end_level` - corresponds to `haarEndLevel`
///
/// # Result
///
/// The segmentation result.
pub fn seg_haar(
    intensities: &[f64],
    weights: Option<&[f64]>,
    raw_intensities: Option<&[f64]>,
    chrom_pos: Vec<Range<usize>>,
    breaks_fdr_q: f64,
    haar_start_level: usize,
    haar_end_level: usize,
) -> HaarSegResult {
    HaarSegResult {
        segments: Vec::new(),
        seg_values: Vec::new(),
    }
}
