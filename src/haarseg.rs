use super::haarseglib;

use stats::Stats;

use std::ops::Range;
use std::ptr;

/// The segment with corresponding value.
#[derive(Debug)]
pub struct HaarSegment {
    /// Range in input vector of segment.
    pub range: Range<usize>,
    /// Value of the signal segment.
    pub value: f64,
}

/// Segmentation result.
#[derive(Debug)]
pub struct HaarSegResult {
    /// Resulting segmented signal.
    pub segments: Vec<HaarSegment>,
    /// Segmented values, same length as `intensities` from input.
    pub seg_values: Vec<f64>,
}

impl HaarSegResult {
    /// Add a segment with the given `range` and `value`.
    fn add_segment(&mut self, range: Range<usize>, value: f64) {
        let old_len = self.seg_values.len();
        self.seg_values.resize(old_len + range.len(), value);
        self.segments.push(HaarSegment { range, value });
    }
}

/// Value to scale MAD with (assuming normality).
const NORMAL_MAD_SCALE: f64 = 0.6745;

// Cumulative density function of the normal distribution
fn pnorm(mean: f64, sd: f64) -> f64 {
    0.5 * (1.0 + unsafe { haarseglib::erfl(mean / sd / 2.0_f64.sqrt()) })
}

/// FDR thresholding.
fn fdr_thresh(x: &[f64], q: f64, sigma: f64) -> f64 {
    let n = x.len();
    if n < 2 {
        return 0.0;
    }

    let m = (1..(n + 1))
        .map(|x| (x as f64) / (n as f64))
        .collect::<Vec<f64>>();
    // Sort x *decreasingly*
    let mut x_sorted = Vec::from(x.clone());
    x_sorted.sort_by(|a, b| b.partial_cmp(a).unwrap());

    // pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)

    let p = x_sorted
        .iter()
        .map(|&q| 2.0 * (1.0 - pnorm(q, sigma)))
        .collect::<Vec<f64>>();
    let k = p.iter()
        .enumerate()
        .filter(|(ref i, &x)| x < m[*i] * q)
        .map(|(ref i, &_x)| *i)
        .collect::<Vec<usize>>();
    if k.last().is_none() {
        const EPSILON: f64 = 1e-16;
        x_sorted.first().unwrap() + EPSILON
    } else {
        x_sorted[*k.last().unwrap()]
    }
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
/// - `breaks_fdr_q` - corresponds to `breaksFdrQ` (R default: `0.001`)
/// - `haar_start_level` - corresponds to `haarStartLevel` (R default: `1`)
/// - `haar_end_level` - corresponds to `haarEndLevel` (R default: `5`)
///
/// # Result
///
/// The segmentation result.
pub fn seg_haar(
    intensities: &[f64],
    _weights: Option<&[f64]>,
    _raw_intensities: Option<&[f64]>,
    chrom_pos: &[Range<usize>],
    breaks_fdr_q: f64,
    haar_start_level: u32,
    haar_end_level: u32,
) -> HaarSegResult {
    let num_probes = intensities.len();
    let mut result = HaarSegResult {
        segments: Vec::new(),
        seg_values: Vec::new(),
    };

    trace!("Starting segmentation using seg_haar()");

    // Perform initial convolution for estimation of peak sigma.
    trace!("Initial convolution for estimating peak sigma");
    let peak_sigma_est = {
        let mut conv_result = vec![0.0_f64; num_probes];
        unsafe {
            haarseglib::HaarConv(
                intensities.as_ptr(),
                ptr::null(),
                num_probes as i32,
                1_i32,
                conv_result.as_mut_ptr(),
            );
        }
        // No non-stationary variance for now, estimate sigma.
        conv_result
            .iter()
            .map(|x| x.abs())
            .collect::<Vec<f64>>()
            .as_slice()
            .median() / NORMAL_MAD_SCALE
    };

    // Perform chromosome-wise segmentation.
    trace!("Performing chromosome-wise segmentation");
    for (i, chrom_range) in chrom_pos.iter().enumerate() {
        trace!("{}-th chromosome from entries {:?}", i, chrom_range);
        let chrom_intensities = &intensities[chrom_range.clone()];
        let chrom_num_probes = chrom_range.len();

        let mut uni_peak_loc = vec![-1_i32; 1];
        for level in haar_start_level..(haar_end_level + 1) {
            trace!("Performing segmentation level {}", level);
            let step_half_size = 2_i32.pow(level);

            // Perform Haar convolution and call peak candidates.
            let mut conv_result = vec![0.0_f64; chrom_num_probes];
            let mut peak_loc_for_c = vec![-1_i32; chrom_num_probes];
            unsafe {
                // `rConvAndPeak(...)`
                haarseglib::HaarConv(
                    chrom_intensities.as_ptr(),
                    ptr::null(),
                    chrom_num_probes as i32,
                    step_half_size,
                    conv_result.as_mut_ptr(),
                );
                haarseglib::FindLocalPeaks(
                    conv_result.as_ptr(),
                    chrom_num_probes as i32,
                    peak_loc_for_c.as_mut_ptr(),
                );
            }

            // Get easier handle on the peak locations.
            let num_peaks = peak_loc_for_c.iter().position(|x| *x == -1).unwrap();
            let peak_loc = Vec::from(&peak_loc_for_c[0..num_peaks]);

            // Perform FDR-based thresholding.
            let conv_result_at_peak_loc = peak_loc
                .iter()
                .map(|&i| conv_result[i as usize])
                .collect::<Vec<f64>>();
            let thresh = fdr_thresh(&conv_result_at_peak_loc, breaks_fdr_q, peak_sigma_est);

            // Apply thresholds and perform unification.
            let unify_win = 2_i32.pow(level - 1);
            let tmp_peak_loc = uni_peak_loc;

            uni_peak_loc = vec![-1_i32; chrom_num_probes];
            unsafe {
                // `rThresAndUnify(...)`
                haarseglib::HardThreshold(
                    conv_result.as_ptr(),
                    thresh,
                    peak_loc_for_c.as_mut_ptr(),
                );
                haarseglib::UnifyLevels(
                    tmp_peak_loc.as_ptr(),
                    peak_loc_for_c.as_ptr(),
                    unify_win,
                    num_peaks as i32,
                    uni_peak_loc.as_mut_ptr(),
                );
            }
        }

        // Collect contig-wise breakpoints and extend result.
        let num_breakpoints = uni_peak_loc.iter().position(|x| *x == -1).unwrap();
        let breakpoints = {
            let mut breakpoints = vec![0_usize];
            breakpoints.extend(uni_peak_loc[0..num_breakpoints].iter().map(|&x| x as usize));
            breakpoints.push(chrom_intensities.len());
            breakpoints
        };
        for i in 1..(breakpoints.len()) {
            let range = Range {
                start: breakpoints[i - 1],
                end: breakpoints[i],
            };
            let value = chrom_intensities[range.clone()].mean();
            result.add_segment(range, value);
        }
    }

    trace!("Done with segmentation using seg_haar()");

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn haar_small_example() {
        let intensities = vec![0.0, 0.1, 0.2, 0.8, 0.9, 1.1, 0.7, 0.0, 0.0];
        let chrom_pos = vec![0..(intensities.len())];

        let result = seg_haar(&intensities, None, None, &chrom_pos, 0.001, 1, 5);

        assert_eq!(result.segments.len(), 3);
        assert_eq!(result.segments[0].range, 0..3);
        assert_eq!(result.segments[1].range, 3..7);
        assert_eq!(result.segments[2].range, 7..9);

        assert_approx_eq!(result.segments[0].value, 0.1);
        assert_approx_eq!(result.segments[1].value, 0.875);
        assert_approx_eq!(result.segments[2].value, 0.0);
    }
}
