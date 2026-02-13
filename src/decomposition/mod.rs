pub mod lu;
pub mod plu;
pub mod qr;

pub use lu::lu_decomposition;
pub use plu::lu_pivot_decomposition;

use crate::Matrix2DError;

#[derive(Debug)]
pub enum DecompositionError {
    NonSquareMatrix,
    SingularMatrix,
    InvalidVector(Matrix2DError),
    InvalidBounds,
}

impl From<Matrix2DError> for DecompositionError {
    fn from(err: Matrix2DError) -> Self {
        DecompositionError::InvalidVector(err)
    }
}
