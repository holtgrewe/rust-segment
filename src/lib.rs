#[macro_use]
extern crate log;
extern crate env_logger;

#[cfg(test)]
#[macro_use]
extern crate test_logger;

#[cfg(test)]
#[macro_use]
extern crate assert_approx_eq;

pub mod haarseg;
pub mod haarseglib;
mod stats;

pub use haarseg::{seg_haar, HaarSegResult};
