use jedvek::{Matrix2D, Rounding};

fn main() {
    let empty_mat: Matrix2D<f64> = Matrix2D::new();
    println!(
        "The Matrix2D::new() function creates a {:?} matrix that is empty: {}\n",
        empty_mat.shape(),
        empty_mat.is_empty()
    );

    let mat = Matrix2D::try_from(vec![vec![2, 1, 0], vec![1, 2, 1], vec![0, 1, 2]]).unwrap();
    println!("{mat}");
    println!("{:?}", mat.shape());
    println!(
        "Matrix size: {}\nMax value: {}\nMin value: {}\n",
        mat.size(),
        mat.max_unchecked(),
        mat.min_unchecked()
    );

    let mat_inverse = mat.inverse().unwrap();
    println!("{}", mat_inverse.round_to_decimal(2)); // Rounds every value in matrix to specified number of decimal places

    let mut singular_mat = Matrix2D::from(&[[1, 2, 3], [4, 5, 6], [8, 10, 12]]);
    let singular_mat_inverse = singular_mat.inverse();
    println!("{singular_mat_inverse:?}");
    println!("Original Matrix:\n{singular_mat}");
    singular_mat.reshape(1).unwrap(); // Reshape sets height to input and auto adjusts width
    println!("1x9 Matrix:\n{singular_mat}");
    singular_mat.reshape(9).unwrap();
    println!("9x1 Matrix:\n{singular_mat}");
}
