#[allow(non_snake_case)]
pub mod matrix2D;
pub mod substitution;

#[derive(Debug)]
pub enum Matrix2DError {
    InconsistentRowLengths,
    NonSquareMatrix,
    SingularMatrix,
    InvalidReshape {
        size: usize,
        new_height: usize,
    },
    InvalidShape {
        input_size: usize,
        output_size: usize,
    },
    InvalidDotShape {
        lhs: usize,
        rhs: usize,
    },
    ConversionFailed {
        from: &'static str,
        to: &'static str,
    },
}
