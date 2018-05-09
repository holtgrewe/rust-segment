// Shared code for segmentation.

use std::ops::Range;

use statrs::distribution::{StudentsT, Univariate};

use super::stats::*;

/// A segment with an optional score.
#[derive(Debug, Clone)]
pub struct Segment {
    /// Range of the segment in the signal array.
    pub range: Range<usize>,

    /// Mean of the segment in normal number space.
    pub mean: f64,
    /// Standard deviation in normal number space.
    pub std_dev: f64,

    /// Mean value of the segment in log2 space.
    pub mean_log2: f64,
    /// Standard deviation in log2 space.
    pub std_dev_log2: f64,
}

/// Helper functions.
impl Segment {
    /// Segment length as `f64`.
    pub fn length(&self) -> f64 {
        self.range.len() as f64
    }

    /// Update the means and standard values from the raw values.
    pub fn update_stats(&mut self, all_values: &[f64]) {
        let values = &all_values[self.range.clone()];
        self.mean = values.mean();
        self.std_dev = values.std_dev();

        let values = values.iter().map(|x| x.log2()).collect::<Vec<f64>>();
        self.mean_log2 = values.mean();
        self.std_dev_log2 = values.std_dev();
    }
}

/// Implementation of CNVnator-style tests on `Segment`.
impl Segment {
    /// Return P-value for the segment, given the signal, using Student's t-test.
    ///
    /// # Arguments
    ///
    /// - `p_thresh` - Threshold on the P-value to use for multiple testing correction.
    pub fn p_value_significant_student(&self, p_thresh: f64) -> f64 {
        // Compute test statistics and obtain distribution.
        let t = (1.0 - self.mean as f64) / self.std_dev * self.length().sqrt();
        let dist =
            StudentsT::new(0.0, 1.0, self.length() - 1.0).expect("Could not find t distribution");
        // Compute p value.
        let p = 1.0 - dist.cdf(t);
        // P-value with Bonferroni correction.
        p * p_thresh / self.length()
    }

    /// Return P-value from Welch's test with Bonferroni correction.
    ///
    /// # Arguments
    ///
    /// - `other` - Other segment to compare to.
    pub fn p_value_difference_welch(&self, other: &Self) -> f64 {
        let x_1 = self.mean;
        let x_2 = other.mean;
        let s_1 = self.std_dev;
        let s_2 = other.std_dev;
        let n_1 = self.length();
        let n_2 = other.length();

        // Helper for computing `x^2`.
        fn sq(x: f64) -> f64 {
            x * x
        }
        /// Helper for computing `x^4`.
        fn q(x: f64) -> f64 {
            x * x * x * x
        }

        // Compute test statistic and degrees of freedom for t Distribution, then create
        // distribution.
        let t = (x_1 - x_2) / (sq(s_1) / n_1 + sq(s_2) / n_2).sqrt();
        let nu = sq(sq(s_1) / n_1 + sq(s_2) / n_2)
            / ((q(s_1) / sq(n_1) / (n_1 - 1.0)) + q(s_2) / sq(n_2) / (n_2 - 1.0));
        let dist = StudentsT::new(0.0, 1.0, nu).unwrap();

        // P-value.
        1.0 - dist.cdf(t)
    }
}

/// A segmentation consists of the segments and the underlying signals.
#[derive(Debug, Clone)]
pub struct Segmentation {
    /// The segments.
    pub segments: Vec<Segment>,
    /// The signals/intensities (from the segment) in number space.
    pub values: Vec<f64>,
    /// The signals/intensities (from the segment) in log2 space.
    pub values_log2: Vec<f64>,
}

// Constructors
impl Segmentation {
    // Construct empty.
    pub fn new() -> Self {
        Segmentation {
            segments: Vec::new(),
            values: Vec::new(),
            values_log2: Vec::new(),
        }
    }

    // Construct with number space values.
    pub fn with_values(segments: Vec<Segment>, values: Vec<f64>) -> Self {
        let values_log2 = values.iter().map(|x| x.log2()).collect::<Vec<f64>>();
        Segmentation {
            segments,
            values,
            values_log2,
        }
    }

    // Construct with number space values.
    pub fn with_values_log2(segments: Vec<Segment>, values_log2: Vec<f64>) -> Self {
        let values = values_log2
            .iter()
            .map(|x| 2_f64.powf(*x))
            .collect::<Vec<f64>>();
        Segmentation {
            segments,
            values,
            values_log2,
        }
    }
}

// Extending a `Segmentation`.
impl Segmentation {
    /// Extend segmentation with a segment. If the value of both segments is within `eps`
    /// then the last segment will be extended. Note that the `p_value` and standard deviation
    // members will not be updated by this!
    pub fn extend(&mut self, segment: &Segment, eps: f64) {
        let last = self.segments.last().map(|s| s.clone());
        match last {
            Some(last) => {
                let last_index = self.segments.len() - 1;
                if (last.mean - segment.mean).abs() < eps {
                    let new_len = self.values.len() + segment.range.len();
                    self.values.resize(new_len, last.mean);
                    self.values_log2.resize(new_len, last.mean_log2);
                    self.segments[last_index] = Segment {
                        range: Range {
                            start: last.range.start,
                            end: segment.range.end,
                        },
                        ..last.clone()
                    }
                } else {
                    let new_len = self.values.len() + segment.range.len();
                    self.values.resize(new_len, segment.mean);
                    self.values_log2.resize(new_len, segment.mean_log2);
                    self.segments.push(segment.clone())
                }
            }
            None => {
                let new_len = self.values.len() + segment.range.len();
                self.values.resize(new_len, segment.mean);
                self.values_log2.resize(new_len, segment.mean_log2);
                self.segments.push(segment.clone());
            }
        }
    }

    /// Update statistics of all segments.
    pub fn update_stats(&mut self, all_values: &[f64]) {
        for ref mut seg in self.segments.iter_mut() {
            seg.update_stats(all_values);
        }
    }

    // Update statistics of all segments with log2 values.
    pub fn update_stats_log2(&mut self, all_values_log2: &[f64]) {
        let all_values = all_values_log2
            .iter()
            .map(|x| 2_f64.powf(*x))
            .collect::<Vec<f64>>();
        for ref mut seg in self.segments.iter_mut() {
            seg.update_stats(&all_values);
        }
    }
}
