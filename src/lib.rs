pub mod core;
pub mod decomposition;

pub use crate::core::Matrix2DError;
pub use crate::core::matrix2D::Matrix2D;
pub use crate::core::matrix2D::Rounding;

pub mod substitution {
    pub use crate::core::substitution::back_substitution;
    pub use crate::core::substitution::forward_substitution;
}
