use crate::decomposition::DecompositionError;
use crate::{Matrix2D, Matrix2DError};

fn norm(matrix: &Matrix2D<f64>) -> f64 {
    // Assumes the matrix is 1 dimensional, i.e. a vector
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
    if hstart > hend
        || wstart > wend
        || hend >= matrix.height
        || wend >= matrix.width
    {
        return Err(Matrix2DError::InconsistentRowLengths);
    }

    // refactor matrix data
    let mut inner = Vec::<f64>::new();
    for h in hstart..(hend + 1){
        for w in wstart..(wend + 1) {
            inner.push(matrix[h][w]);
        }
    }

    let result = Matrix2D::from_flat(inner, 0.0, hend - hstart + 1, wend - wstart + 1);
    if let Ok(result) = result {
        return Ok(result);
    } else {
        println!("Error when creating result");
        return Err(Matrix2DError::InconsistentRowLengths);
    }
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
    // Returns the compressed form of the QR factorization
    let mut matrix: Matrix2D<f64> = matrix.try_into()?;
    if matrix.height < matrix.width {
        return Err(DecompositionError::NonSquareMatrix);
    }

    let mut q: Matrix2D<f64> = Matrix2D::identity(matrix.height);

    for j in 0..matrix.width {
        let hvec = submatrix(&matrix, j, matrix.width - 1, j, j).unwrap();
        println!("hvec: {}", hvec);
        let normx = norm(&hvec);
        println!("normx: {}", &normx);
        let s = -sign(matrix[j][j]);
        println!("s: {}", &s);
        let u1 = matrix[j][j] - s * normx;
        println!("u1 {}", &u1);
        let mut w = hvec / u1;

        // I can't currently think of a better way to mutate the first element
        // of an Matrix2D which doesn't have a statically defined size at compile time
        let w_iter = w.rows_mut();
        for row in w_iter.into_iter() {
            for elem in row {
                *elem = 1_f64;
                break;
            }
            break;
        }

        println!("w: {}", &w);

        //matrix[j+1:end, j] = w(2:end)

        //matrix[j][j] = s * normx;

        let tau = -s * u1 / normx;
        println!("tau: {}", tau);


        // Now modify the original matrix
        let r_end = submatrix(&matrix, j, matrix.height - 1, 0, matrix.width - 1).unwrap();
        let q_end = submatrix(&q, 0, matrix.height - 1, j, matrix.width - 1).unwrap();
        println!("r_end: {}", r_end);
        println!("q_end: {}", q_end);

        let mut tau_w = w.clone();
        let tau_w_iter = tau_w.rows_mut();
        for row in tau_w_iter.into_iter() {
            for elem in row {
                *elem = *elem * tau;
            }
        }
        
        println!("tau_w: {}", tau_w);
        println!("Getting r_end");
        let r_sub_by = (tau_w.clone()) * (w.transpose() * &r_end);
        println!("tau_w shape: {:?}", tau_w.shape());
        println!("w transpose: {:?}", w.transpose());
        println!("w_transpose shape: {:?}", w.transpose().shape());
        println!("r_end shape: {:?}", r_end.shape());
        println!("w' * r_end: {}", (w.transpose() * &r_end));
        println!("r sub by: {}", &r_sub_by);

        println!("q_end shape: {:?}", q_end.shape());
        println!("w shape: {:?}", w.shape());
        let q_sub_by = (q_end * w) * tau_w.transpose();
        println!("q sub by: {}", &q_sub_by);
        
        let r_sub_by_shape = r_sub_by.shape();
        let q_sub_by_shape = q_sub_by.shape();

        println!("r_sub_by_shape: {:?}", r_sub_by_shape);
        println!("q_sub_by_shape: {:?}", q_sub_by_shape);
        for row in j..matrix.height {
            for col in 0..matrix.width {
                if row - j < r_sub_by_shape.0 && col < r_sub_by_shape.1 {
                    matrix[row][col] = matrix[row][col] - r_sub_by[row - j][col];
                }
                if row - j < q_sub_by_shape.1 && col < q_sub_by_shape.0 {
                    q[col][row] = q[col][row] - q_sub_by[col][row - j];
                }
            }
        }

        println!("q: {}", &q);
        println!("r: {}", &matrix);
    }

    Ok((q, matrix))
}

#[cfg(test)]
mod tests {
    use super::*;

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
        println!("q: {} \n r: {}", q, r);
        let expected = q * r;
        assert_eq!(matr, expected);
    }

    #[test]
    fn qr_test_3x3() {
        let matr: Matrix2D<f64> = Matrix2D::from(&[[2.0, -1.0, -2.0], [-4.0, 6.0, -3.0], [-4.0, -2.0, 8.0]]);
        let (q, r) = qr_factorization(&matr).unwrap();
        println!("q: {} \n r: {}", q, r);
        let expected = q * r;
        assert_eq!(matr, expected);
    }
}