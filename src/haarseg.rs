use super::haarseglib;
pub use super::shared::{Segment, Segmentation};

use stats::Stats;

use std::ops::Range;
use std::ptr;

/// Value to scale MAD with (assuming normality).
const NORMAL_MAD_SCALE: f64 = 0.6745;

// Cumulative density function of the normal distribution
fn pnorm(x: f64, mean: f64, sd: f64) -> f64 {
    0.5 * unsafe { haarseglib::erfl((x - mean) / sd / 2.0_f64.sqrt()) }
}

/// Compute p value from distribution in `mu = 0` and given `sigma` for the given `x`.
fn pvalue(x: f64, sigma: f64) -> f64 {
    2.0 * (1.0 - pnorm(x, 0.0, sigma))
    // use statrs::distribution::{Univariate, Poisson};
    // let dist = Poisson::new(sigma * sigma).unwrap();
    // return 1.0 - dist.cdf(x);
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

    let p = x_sorted
        .iter()
        .map(|&q| pvalue(q, sigma))
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
/// - `values_log2` - corresponds to `I`
/// - `weights` - corresponds to `W`
/// - `raw_values` - corresponds to `rI`
/// - `chrom_pos` - corresponds to `chromPos`
/// - `breaks_fdr_q` - corresponds to `breaksFdrQ` (R default: `0.001`)
/// - `haar_start_level` - corresponds to `haarStartLevel` (R default: `1`)
/// - `haar_end_level` - corresponds to `haarEndLevel` (R default: `5`)
///
/// # Result
///
/// The segmentation result.
pub fn seg_haar(
    values_log2: &[f64],
    _weights: Option<&[f64]>,
    _raw_values: Option<&[f64]>,
    chrom_pos: &[Range<usize>],
    breaks_fdr_q: f64,
    haar_start_level: u32,
    haar_end_level: u32,
) -> Segmentation {
    let num_probes = values_log2.len();
    let mut result = Segmentation::new();

    trace!("Starting segmentation using seg_haar()");

    // Perform initial convolution for estimation of peak sigma.
    trace!("Initial convolution for estimating peak sigma");
    let peak_sigma_est = {
        let mut conv_result = vec![0.0_f64; num_probes];
        unsafe {
            haarseglib::HaarConv(
                values_log2.as_ptr(),
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
            .filter(|&x| x > 1e-5)
            .collect::<Vec<f64>>()
            .as_slice()
            .median() / NORMAL_MAD_SCALE * 5.0
    };

    // Perform chromosome-wise segmentation.
    trace!("Performing chromosome-wise segmentation");
    for (i, chrom_range) in chrom_pos.iter().enumerate() {
        trace!("{}-th chromosome from entries {:?}", i, chrom_range);
        let chrom_values = &values_log2[chrom_range.clone()];
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
                    chrom_values.as_ptr(),
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
            breakpoints.push(chrom_values.len());
            breakpoints
        };
        let offset = result.values_log2.len();
        for i in 1..(breakpoints.len()) {
            let range_here = Range {
                start: breakpoints[i - 1],
                end: breakpoints[i],
            };

            let range = Range {
                start: offset + range_here.start,
                end: offset + range_here.end,
            };

            let mean_log2 = chrom_values[range_here.clone()].mean();
            let std_dev_log2 = chrom_values[range_here.clone()].std_dev();
            let values = chrom_values[range_here.clone()]
                .iter()
                .map(|x| 2_f64.powf(*x))
                .collect::<Vec<f64>>();
            let mean = values.mean();
            let std_dev = values.std_dev();

            result.extend(
                &Segment {
                    range,
                    mean,
                    std_dev,
                    mean_log2,
                    std_dev_log2,
                },
                1e-5,
            );
        }
    }

    trace!("Done with segmentation using seg_haar()");

    result
}

/// Refine segmentation by moving breaks left or right.
///
/// # Returns
///
/// A pair (`n`, `res`) with
///
/// - `n` -- The number of breaks that were adjusted.
/// - `res` -- The updated segmentation result.
pub fn adjust_breaks(segmentation: &Segmentation, values: &[f64]) -> (usize, Segmentation) {
    let mut peak_locs = segmentation.segments[1..segmentation.segments.len()]
        .iter()
        .map(|x| x.range.start as i32)
        .collect::<Vec<i32>>();
    peak_locs.push(-1);
    let peak_locs = peak_locs;
    let mut new_peak_locs = vec![-1_i32; peak_locs.len()];
    let num_adjusted = unsafe {
        haarseglib::AdjustBreaks(
            values.as_ptr(),
            values.len() as i32,
            peak_locs.as_ptr(),
            new_peak_locs.as_mut_ptr(),
        )
    };

    let mut segments: Vec<Segment> = Vec::new();
    let mut seg_values: Vec<f64> = Vec::new();
    let mut prev: usize = 0;
    for (i, end) in new_peak_locs.iter().take_while(|i| **i >= 0).enumerate() {
        if prev != *end as usize {
            let segment = segmentation.segments[i].clone();
            segments.push(Segment {
                range: Range {
                    start: prev,
                    end: *end as usize,
                },
                ..segment
            });
            seg_values.resize(values.len() + *end as usize - prev, segment.mean);
            prev = *end as usize;
        }
    }

    if prev != seg_values.len() {
        let segment = segmentation.segments.last().unwrap().clone();
        segments.push(Segment {
            range: Range {
                start: prev,
                end: seg_values.len(),
            },
            ..segment
        });
        let new_len = values.len() + seg_values.len() - prev;
        seg_values.resize(new_len, segment.mean);
    }

    let mut segmentation = Segmentation::with_values(segments, seg_values);
    segmentation.update_stats(values);

    (num_adjusted as usize, segmentation)
}

/// Implement simple approach for aberrant interval detection as proposed in HaarSeg paper.
///
/// Segments with a log2 value above `m * sigma` will be accepted, the others rejected.
pub fn reject_nonaberrant(segmentation: &Segmentation, values: &[f64], m: f64) -> Segmentation {
    // Estimate sigma.
    let est_sigma = segmentation
        .values
        .iter()
        .map(|x| x.log2())
        .zip(values.iter())
        .map(|(y, x)| (y - *x).abs())
        .collect::<Vec<f64>>()
        .as_slice()
        .median() / 0.6745;

    // Rebuild new `Segmentation`.
    let mut segments: Vec<Segment> = Vec::new();
    let mut seg_values: Vec<f64> = Vec::new();
    let mut curr: Option<Segment> = None;
    for seg in &segmentation.segments {
        let range = seg.range.clone();
        let (mean, mean_log2) = if seg.mean_log2.abs() > m * est_sigma {
            (seg.mean, seg.mean_log2)
        } else {
            (1.0, 0.0)
        };

        curr = match curr {
            Some(curr) => {
                if mean_log2 != curr.mean_log2 {
                    // Start new segment.
                    segments.push(curr);
                    Some(Segment {
                        range,
                        mean,
                        mean_log2,
                        std_dev: 0.0,
                        std_dev_log2: 0.0,
                    })
                } else {
                    // Extend old segment.
                    Some(Segment {
                        range: Range {
                            end: seg.range.end,
                            ..curr.range
                        },
                        ..curr
                    })
                }
            }
            None => Some(Segment {
                range,
                mean,
                mean_log2,
                std_dev: 0.0,
                std_dev_log2: 0.0,
            }),
        };

        seg_values.extend_from_slice(&vec![mean; seg.range.len()]);
    }
    // Add last segment, if any, values are already expanded.
    if let Some(curr) = curr {
        segments.push(curr);
    }

    assert_eq!(values.len(), values.len());
    assert!(segments.len() <= segmentation.segments.len());

    let mut segmentation = Segmentation::with_values(segments, seg_values);
    segmentation.update_stats(values);

    segmentation
}

/// Implement simple approach for aberrant interval detection as proposed in CNVnator paper.
///
/// Segments with a value above `m * sigma` will be accepted, the others rejected.
pub fn reject_nonaberrant_pvalue(
    segmentation: &Segmentation,
    values: &[f64],
    p_value_threshold: f64,
) -> Segmentation {
    // Rebuild new `Segmentation`.
    let mut segments: Vec<Segment> = Vec::new();
    let mut seg_values: Vec<f64> = Vec::new();
    let mut curr: Option<Segment> = None;
    for seg in &segmentation.segments {
        let range = seg.range.clone();
        let (mean, mean_log2) = if range.len() >= 2
            && seg.p_value_significant_student(p_value_threshold) < p_value_threshold
        {
            (seg.mean, seg.mean_log2)
        } else {
            (1.0, 0.0)
        };

        curr = match curr {
            Some(curr) => {
                if mean != curr.mean {
                    // Start new segment.
                    segments.push(curr);
                    Some(Segment {
                        range,
                        mean,
                        mean_log2,
                        std_dev: 0.0,
                        std_dev_log2: 0.0,
                    })
                } else {
                    // Extend old segment.
                    Some(Segment {
                        range: Range {
                            end: seg.range.end,
                            ..curr.range
                        },
                        ..curr
                    })
                }
            }
            None => Some(Segment {
                range,
                mean,
                mean_log2,
                std_dev: 0.0,
                std_dev_log2: 0.0,
            }),
        };

        seg_values.extend_from_slice(&vec![mean; seg.range.len()]);
    }
    // Add last segment, if any, seg_values are already expanded.
    if let Some(curr) = curr {
        segments.push(curr);
    }

    assert_eq!(seg_values.len(), values.len());
    assert!(segments.len() <= segmentation.segments.len());

    let mut segmentation = Segmentation::with_values(segments, seg_values);
    segmentation.update_stats(values);

    segmentation
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn haar_small_example() {
        let values = vec![0.0, 0.1, 0.2, 0.8, 0.9, 1.1, 0.7, 0.0, 0.0];
        let chrom_pos = vec![0..(values.len())];

        let result = seg_haar(&values, None, None, &chrom_pos, 0.001, 1, 5);

        assert_eq!(result.segments.len(), 3);
        assert_eq!(result.segments[0].range, 0..3);
        assert_eq!(result.segments[1].range, 3..7);
        assert_eq!(result.segments[2].range, 7..9);

        assert_approx_eq!(result.segments[0].value, 0.1);
        assert_approx_eq!(result.segments[1].value, 0.875);
        assert_approx_eq!(result.segments[2].value, 0.0);
    }
}
