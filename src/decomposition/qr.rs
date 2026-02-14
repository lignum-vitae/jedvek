use crate::decomposition::DecompositionError;
use crate::{Matrix2D, Matrix2DError};

fn norm(matrix: &Matrix2D<f64>) -> f64 {
    // Assumes the matrix is 1 dimensional, i.e. a vector
    assert!(matrix.width == 1 || matrix.height == 1);
    let mut sum: f64 = 0.0;
    for m in 0..matrix.height {
        for n in 0..matrix.width {
            sum += matrix[m][n].powi(2)
        }
    }

    sum.sqrt()
}

fn submatrix(
    matrix: &Matrix2D<f64>,
    hstart: usize,
    hend: usize,
    wstart: usize,
    wend: usize,
) -> Result<Matrix2D<f64>, Matrix2DError> {
    // Uses 0-indexing
    if hstart > hend || wstart > wend || hend >= matrix.height || wend >= matrix.width {
        return Err(Matrix2DError::OutOfBounds);
    }

    // Refactor matrix data
    let mut inner = Vec::<f64>::new();
    for h in hstart..(hend + 1) {
        for w in wstart..(wend + 1) {
            inner.push(matrix[h][w]);
        }
    }

    Matrix2D::from_flat(inner, 0.0, hend - hstart + 1, wend - wstart + 1)
}

fn sign(entry: f64) -> f64 {
    if entry < 0.0 {
        return -1.0;
    }
    1.0
}

pub fn qr_factorization<M>(matrix: M) -> Result<(Matrix2D<f64>, Matrix2D<f64>), DecompositionError>
where
    M: TryInto<Matrix2D<f64>, Error = Matrix2DError>,
{
    // Returns the non-compressed form of the QR factorization, less space efficient
    // but returns a direct result at the end
    let mut matrix: Matrix2D<f64> = matrix.try_into()?;
    if matrix.height < matrix.width {
        return Err(DecompositionError::InvalidBounds);
    }

    let mut q: Matrix2D<f64> = Matrix2D::identity(matrix.height);

    for j in 0..matrix.width {
        let hvec = submatrix(&matrix, j, matrix.width - 1, j, j).unwrap();
        let normx = norm(&hvec);
        let s = -sign(matrix[j][j]);
        let u1 = matrix[j][j] - s * normx;
        let mut w = hvec / u1;

        // I can't currently think of a better way to mutate the first element
        // of an Matrix2D which doesn't have a statically defined size at compile time
        let w_iter = w.rows_mut();
        if let Some(row) = w_iter.into_iter().next()
            && let Some(elem) = row.iter_mut().next()
        {
            *elem = 1_f64;
        }

        let tau = -s * u1 / normx;

        let r_end = submatrix(&matrix, j, matrix.height - 1, 0, matrix.width - 1).unwrap();
        let q_end = submatrix(&q, 0, matrix.height - 1, j, matrix.width - 1).unwrap();

        let mut tau_w = w.clone();
        let tau_w_iter = tau_w.rows_mut();
        for row in tau_w_iter.into_iter() {
            for elem in row {
                *elem *= tau;
            }
        }
        let r_sub_by = (tau_w.clone()) * (w.transpose() * &r_end);
        let q_sub_by = (q_end * w) * tau_w.transpose();

        let r_sub_by_shape = r_sub_by.shape();
        let q_sub_by_shape = q_sub_by.shape();

        for row in j..matrix.height {
            for col in 0..matrix.width {
                if row - j < r_sub_by_shape.0 && col < r_sub_by_shape.1 {
                    matrix[row][col] -= r_sub_by[row - j][col];
                }
                if row - j < q_sub_by_shape.1 && col < q_sub_by_shape.0 {
                    q[col][row] -= q_sub_by[col][row - j];
                }
            }
        }
    }

    Ok((q, matrix))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn check_orthogonal(matrix: &Matrix2D<f64>) -> bool {
        let mat_tol_f64 = 10.0_f64.powf(-12.0);
        for col in 0..matrix.width {
            for other_col in (col + 1)..matrix.width {
                let mut product = 0_f64;
                for row in 0..matrix.height {
                    product += matrix[row][col] * matrix[row][other_col];
                }

                if product.abs() > mat_tol_f64 {
                    return false;
                }
            }
        }
        true
    }

    fn check_upper_triangular(matrix: &Matrix2D<f64>) -> bool {
        let mat_tol_f64 = 10.0_f64.powf(-12.0);
        for row in 0..matrix.height {
            for col in 0..matrix.width {
                if row > col && matrix[row][col].abs() > mat_tol_f64 {
                    return false;
                } 
            }
        }
        true
    }

    fn equal_within_tol(lhs: &Matrix2D<f64>, rhs: &Matrix2D<f64>) -> bool {
        let mat_tol_f64 = 10.0_f64.powf(-12.0);
        if lhs.height != rhs.height && lhs.width != rhs.width {
            return false;
        }
        if lhs.height == 0 || rhs.height == 0 {
            return true;
        }

        for r in 0..lhs.height {
            for c in 0..lhs.width {
                if (lhs[r][c] - rhs[r][c]).abs() > mat_tol_f64 {
                    return false;
                }
            }
        }
        true
    }

    #[test]
    fn basic_norm() {
        let vec: Matrix2D<f64> = Matrix2D::from(&[[-1.0, 3.0, 5.0]]);
        let result = norm(&vec);
        assert_eq!(result, 5.916079783099616);
    }

    #[test]
    fn basic_submatrix() {
        let vec = Matrix2D::from(&[[-1.0, 2.0, 3.0], [6.0, 5.0, 4.0]]);
        assert_eq!(vec.shape(), (2, 3));
        let result = submatrix(&vec, 0, 1, 1, 2).unwrap();
        let expected: Matrix2D<f64> = Matrix2D::from(&[[2.0, 3.0], [5.0, 4.0]]);
        assert_eq!(result, expected);
    }

    #[test]
    fn basic_qr() {
        let matr: Matrix2D<f64> = Matrix2D::from(&[[-1.0, 3.0], [1.0, 5.0]]);
        let (q, r) = qr_factorization(&matr).unwrap();
        let expected = &q * &r;
        assert!(check_orthogonal(&q));
        assert!(check_upper_triangular(&r));
        assert!(equal_within_tol(&matr, &expected));
    }

    #[test]
    fn qr_test_3x3() {
        let matr: Matrix2D<f64> =
            Matrix2D::from(&[[2.0, -1.0, -2.0], [-4.0, 6.0, -3.0], [-4.0, -2.0, 8.0]]);
        let (q, r) = qr_factorization(&matr).unwrap();
        let expected = &q * &r;
        assert!(check_orthogonal(&q));
        assert!(check_upper_triangular(&r));
        assert!(equal_within_tol(&matr, &expected));
    }

    #[test]
    fn qr_test_4x4() {
        let matr: Matrix2D<f64> = Matrix2D::from(&[
            [1.0, -1.0, 0.0, 0.0],
            [1.0, 1.0, 1.0, 0.0],
            [0.0, 1.0, 1.0, 1.0],
            [0.0, 0.0, 1.0, 1.0],
        ]);
        let (q, r) = qr_factorization(&matr).unwrap();
        let expected = &q * &r;
        assert!(check_orthogonal(&q));
        assert!(check_upper_triangular(&r));
        assert!(equal_within_tol(&matr, &expected));
    }
}
