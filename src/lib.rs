pub mod haarseg;
pub mod haarseglib;

pub use haarseg::{seg_haar, HaarSegResult};

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
