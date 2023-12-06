use num::complex::{Complex, Complex64};

type Cvecf = Vec<Complex<f64>>;

pub fn corr_init_vec_ada_order_2__2e_2(init_vec: &mut Cvecf) {
	init_vec[0] = Complex::new(0.0, 0.0); // 0 ada
	init_vec[1] = Complex::new(0.0, 0.0); // 1 spa
	init_vec[2] = Complex::new(0.0, 0.0); // 2 smsp
	init_vec[3] = Complex::new(1.0, 0.0); // 3 sz
	init_vec[4] = Complex::new(1.0, 0.0); // 4 szsz
	init_vec[5] = Complex::new(0.0, 0.0); // 5 szad
	init_vec[6] = Complex::new(0.0, 0.0); // 6 szsp
	init_vec[7] = Complex::new(0.0, 0.0); // 7 sp
	init_vec[8] = Complex::new(0.0, 0.0); // 8 a
	init_vec[9] = Complex::new(0.0, 0.0); // 9 szsm
	init_vec[10] = Complex::new(0.0, 0.0); // 10 spad
	init_vec[11] = Complex::new(0.0, 0.0); // 11 spsp
	init_vec[12] = Complex::new(0.0, 0.0); // 12 adad
}