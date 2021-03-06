#[macro_use]
extern crate log;
extern crate env_logger;

#[cfg(test)]
#[macro_use]
extern crate test_logger;

#[cfg(test)]
#[macro_use]
extern crate assert_approx_eq;

extern crate statrs;

pub mod haarseg;
pub mod haarseglib;
pub mod shared;
mod stats;

pub use haarseg::{adjust_breaks, reject_nonaberrant, reject_nonaberrant_pvalue, seg_haar};

pub use haarseg::{Segment, Segmentation};
