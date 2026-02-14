use jedvek::{Matrix2D, Rounding};

fn main() {
    let empty_mat: Matrix2D<f64> = Matrix2D::new();
    println!(
        "The Matrix2D::new() function creates a {:?} matrix that is empty: {}\n",
        empty_mat.shape(),
        empty_mat.is_empty()
    );

    let mat = Matrix2D::try_from(vec![vec![2, 1, 6], vec![1, 2, 1], vec![0, 1, 2]]).unwrap();
    println!("Original Matrix:\n{mat}\n");
    let mut mat_t = mat.transpose();
    println!("Transposed Matrix:\n{mat_t}\n");
    mat_t.swap_rows(1, 2); // Rows are 0-indexed
    println!("First and Second row swap:\n{mat_t}\n");
    println!("{:?}", mat.shape());
    println!(
        "Matrix size: {}\nMax value: {}\nMin value: {}\n",
        mat.size(),
        mat.max_unchecked(),
        mat.min_unchecked()
    );

    let mat_inverse = mat.inverse().unwrap();
    println!("Matrix Inverse:\n{}\n", mat_inverse.round_to_decimal(2)); // Rounds every value in matrix to specified number of decimal places

    let mut singular_mat = Matrix2D::from(&[[1, 2, 3], [4, 5, 6], [8, 10, 12]]);
    println!("Uninverted singular matrix:\n{singular_mat}");
    let singular_mat_inverse = singular_mat.inverse();
    println!("{singular_mat_inverse:?}\n");
    println!("Original Matrix:\n{singular_mat}");
    singular_mat.reshape(1).unwrap(); // Reshape sets height to input and auto adjusts width
    println!("1x9 Matrix:\n{singular_mat}");
    singular_mat.reshape(9).unwrap();
    println!("9x1 Matrix:\n{singular_mat}\n");

    let flat_array = [1, 3, 5, 6, 7, 8];
    // Parameters for from_flat are AsRef<[T]>, default_val, height, width
    let matrix = Matrix2D::from_flat(flat_array, 0, 3, 3).unwrap_or(Matrix2D::full(0, 3, 3));
    println!("Flat array:\n{flat_array:?}\n2D Matrix from flat:\n{matrix}\n");

    println!("Identity matrix:\n{}\n", Matrix2D::<f64>::identity(3));

    let filled_matrix = Matrix2D::<i32>::full(4, 3, 3);
    println!("Filled matrix:\n{filled_matrix}\n");
    println!(
        "Filled matrix map subtract 2:\n{}\n",
        filled_matrix.map(|x| x - 2)
    );
}
