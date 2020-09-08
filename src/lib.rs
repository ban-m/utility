//!A tiny implementation of utility package such as file input output
#![crate_name = "utility"]
#![crate_type = "lib"]
#![warn(missing_docs)]

extern crate bio;
extern crate dtw;
extern crate fast5wrapper;
extern crate histogram_minimizer;
extern crate knn_predictor;
extern crate rand;
extern crate rayon;
extern crate squiggler;
extern crate num;
/// Utilities I'm using.
pub mod utilities;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
