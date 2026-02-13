use jedvek::Matrix2D;

fn main() {
    let mat_1 = Matrix2D::from(&[[1, 2], [3, 5]]);
    let mat_2 = Matrix2D::from(&[[2, 5], [7, 9]]);
    println!("mat_1:\n{mat_1}\nmat_2:\n{mat_2}\n");
    println!("mat_1 * mat_2:\n{}\n", &mat_1 * &mat_2);

    let scalar = 2;
    let mat_mul_scalar = &mat_1 * scalar;
    println!("mat_1 * {scalar}:\n{mat_mul_scalar}");
    let mat_to_scalar = Matrix2D::from_flat([2], 0, 1, 1)
        .unwrap()
        .as_scalar_unchecked();
    println!("{}\n", &mat_1 * mat_to_scalar);

    println!("mat_1 * {scalar} / {scalar}:\n{}", &mat_mul_scalar / scalar);
}
