use nalgebra::{ClosedAdd, ComplexField, SVector, Scalar, Vector3};
use num::complex::{Complex, Complex64};
use num_traits::Zero;
use std::io::Write;
use std::{fs::File, ops::Mul};
use std::{io, result};
mod corr_init_vec_ada_order_2__2e_2;
mod corr_sys_ada_order_2;

type Cvecf = Vec<Complex<f64>>;

const A_COEFF: [f64; 6] = [0.0, 1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 1.0 / 2.0];
const B1_COEFF: [f64; 1] = [1.0 / 4.0];
const B2_COEFF: [f64; 2] = [3.0 / 32.0, 9.0 / 32.0];
const B3_COEFF: [f64; 3] = [1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0];
const B4_COEFF: [f64; 4] = [439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0];
const B5_COEFF: [f64; 5] = [
    -8.0 / 27.0,
    2.0,
    -3544.0 / 2565.0,
    1859.0 / 4104.0,
    -11.0 / 40.0,
];
const CH_COEFF: [f64; 6] = [
    16.0 / 135.0,
    0.0,
    6656.0 / 12825.0,
    28561.0 / 56430.0,
    -9.0 / 50.0,
    2.0 / 55.0,
];
const CT_COEFF: [f64; 6] = [
    -1.0 / 360.0,
    0.0,
    128.0 / 4275.0,
    2197.0 / 75240.0,
    -1.0 / 50.0,
    -2.0 / 55.0,
];
const MIN_INTERVALS: usize = 10;
const MAX_INTERVALS: usize = 1000000000;
const TOL: f64 = 0.0001;

fn rk4(
    t_start: f64,
    t_end: f64,
    num_intervals: usize,
    func: fn(f64, &Cvecf, &mut Cvecf),
    init_vec: &Cvecf,
) -> Vec<Cvecf> {
    let vec_len = init_vec.len();
    let h = (t_end - t_start) / num_intervals as f64;
    let mut result = vec![vec![Complex::new(0.0, 0.0); num_intervals + 1]; vec_len];
    let mut curr_vec = init_vec.clone();
    let mut curr_t = t_start;
    let mut k1 = vec![Complex::new(0.0, 0.0); vec_len];
    let mut k2 = vec![Complex::new(0.0, 0.0); vec_len];
    let mut k3 = vec![Complex::new(0.0, 0.0); vec_len];
    let mut k4 = vec![Complex::new(0.0, 0.0); vec_len];
    let mut temp = vec![Complex::new(0.0, 0.0); vec_len];
    for i in 0..(num_intervals + 1) {
        //println!("{} {}", i, curr_t);

        // k1
        func(curr_t, &curr_vec, &mut k1);

        // yn + h*k1/2
        for j in 0..init_vec.len() {
            temp[j] = curr_vec[j] + h * k1[j] / 2.0;
        }

        // k2
        func(curr_t + h / 2.0, &temp, &mut k2);

        // yn + h*k2/2
        for j in 0..init_vec.len() {
            temp[j] = curr_vec[j] + h * k2[j] / 2.0;
        }

        // k3
        func(curr_t + h / 2.0, &temp, &mut k3);

        // yn + h*k3
        for j in 0..init_vec.len() {
            temp[j] = curr_vec[j] + h * k3[j];
        }

        // k4
        func(curr_t + h, &temp, &mut k4);

        // update
        for j in 0..vec_len {
            result[j][i] = curr_vec[j].clone();
            //println!("{}", h / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]));
            curr_vec[j] += h / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
        }

        //panic!("haha");
        curr_t += h;
    }
    return result;
}

fn get_norm(vec: &Cvecf) -> f64 {
    let mut res_add = 0.0;
    for elem in vec {
        res_add += (elem * elem.conj()).re;
        //println!("{} {}", res_add, elem)
    }
    return res_add.sqrt();
}

fn rkf45(
    t_start: f64,
    t_end: f64,
    func: fn(f64, &Cvecf, &mut Cvecf),
    init_vec: &Cvecf,
) -> (Vec<f64>, Vec<Cvecf>) {
    let vec_len = init_vec.len();
    let h_max = (t_end - t_start) / MIN_INTERVALS as f64;
    let h_min = (t_end - t_start) / MAX_INTERVALS as f64;
    let mut h = h_min;
    let mut result = vec![Cvecf::new(); vec_len];
    let mut result_time = Vec::<f64>::new();
    let mut curr_vec = init_vec.clone();
    let mut curr_t = t_start;
    let mut k1 = vec![Complex::new(0.0, 0.0); vec_len];
    let mut k2 = vec![Complex::new(0.0, 0.0); vec_len];
    let mut k3 = vec![Complex::new(0.0, 0.0); vec_len];
    let mut k4 = vec![Complex::new(0.0, 0.0); vec_len];
    let mut k5 = vec![Complex::new(0.0, 0.0); vec_len];
    let mut k6 = vec![Complex::new(0.0, 0.0); vec_len];
    let mut temp = vec![Complex::new(0.0, 0.0); vec_len];
    let mut te_vec = vec![Complex::new(0.0, 0.0); vec_len];
    while curr_t <= t_end {
        //println!("{} {}", curr_t, h);
        func(curr_t + A_COEFF[0] * h, &curr_vec, &mut k1);
        for j in 0..vec_len {
            k1[j] *= h;
            temp[j] = curr_vec[j] + B1_COEFF[0] * k1[j];
        }
        func(curr_t + A_COEFF[1] * h, &temp, &mut k2);
        for j in 0..vec_len {
            k2[j] *= h;
            temp[j] = curr_vec[j] + B2_COEFF[0] * k1[j] + B2_COEFF[1] * k2[j];
        }
        func(curr_t + A_COEFF[2] * h, &temp, &mut k3);
        for j in 0..vec_len {
            k3[j] *= h;
            temp[j] = curr_vec[j] + B3_COEFF[0] * k1[j] + B3_COEFF[1] * k2[j] + B3_COEFF[2] * k3[j];
        }
        func(curr_t + A_COEFF[3] * h, &temp, &mut k4);
        for j in 0..vec_len {
            k4[j] *= h;
            temp[j] = curr_vec[j]
                + B4_COEFF[0] * k1[j]
                + B4_COEFF[1] * k2[j]
                + B4_COEFF[2] * k3[j]
                + B4_COEFF[3] * k4[j];
        }
        func(curr_t + A_COEFF[4] * h, &temp, &mut k5);
        for j in 0..vec_len {
            k5[j] *= h;
            temp[j] = curr_vec[j]
                + B4_COEFF[0] * k1[j]
                + B5_COEFF[1] * k2[j]
                + B5_COEFF[2] * k3[j]
                + B5_COEFF[3] * k4[j]
                + B5_COEFF[4] * k5[j];
        }
        func(curr_t + A_COEFF[5] * h, &temp, &mut k6);
        for j in 0..vec_len {
            k6[j] *= h;
            temp[j] = curr_vec[j]
                + CH_COEFF[0] * k1[j]
                + CH_COEFF[1] * k2[j]
                + CH_COEFF[2] * k3[j]
                + CH_COEFF[3] * k4[j]
                + CH_COEFF[4] * k5[j]
                + CH_COEFF[5] * k6[j];
        }
        for j in 0..vec_len {
            te_vec[j] = CT_COEFF[0] * k1[j]
                + CT_COEFF[1] * k2[j]
                + CT_COEFF[2] * k3[j]
                + CT_COEFF[3] * k4[j]
                + CT_COEFF[4] * k5[j]
                + CT_COEFF[5] * k6[j];
        }
        println!("{}", temp[0]);
        let te = get_norm(&te_vec);
        if te <= TOL {
            for j in 0..vec_len {
                result[j].push(curr_vec[j].clone());
                curr_vec[j] = temp[j].clone();
            }
            result_time.push(curr_t.clone());
            curr_t += h
        }
        h = (0.9 * h * (TOL / te).powf(1.0 / 5.0)).max(h_min).min(h_max);
        println!("h: {} te: {}", h, te);
        if h == 0.0 {
            panic!("haha");
        }
    }
    return (result_time, result);
}

fn save_data(path: &str, result_time: &Vec<f64>, result: &Vec<Cvecf>) {
    let mut output = File::create(path).unwrap();
    let mut time_output = String::new();
    let mut count = 0;
    for t in result_time {
        time_output += &t.to_string();
        time_output += " ";
        count += 1;
    }
    writeln!(output, "{}", time_output).unwrap();
    for ind_cvecf in result {
        //println!("{}", ind_cvecf.len());
        let mut line = String::new();
        let mut count_2 = 0;
        for ind_num in ind_cvecf {
            count_2 += 1;
            let mut add_sign = "+";
            if ind_num.im < 0.0 {
                add_sign = "";
            }
            line += &(ind_num.re.to_string() + add_sign + &ind_num.im.to_string() + "j");
            line += " ";
        }
        writeln!(output, "{}", line).unwrap();
        println!("{}, {}", count, count_2);
    }
}

fn diff_equ_test(time: f64, yn: &Cvecf, fty: &mut Cvecf) {
    fty[0] = time.sin().powi(2) * yn[0];
}

fn diff_equ_test_2(time: f64, yn: &Cvecf, fty: &mut Cvecf) {
    fty[0] = 1.0 + yn[0] * yn[0];
}

fn main() {
    //let a: Complex<i32> = Complex::new(10, 20);
    //println!("Hello, world! {}", &(a.re.to_string() + "+" + &a.im.to_string() + "j"));
    println!("Computing...");
    //let init_vec = vec![Complex64::new(1.0, 0.0)];
    //let t_start = 0.0;
    //let t_end: f64 = 5.0;
    //let num_intervals = 24;
    //let result = rk4(t_start, t_end, num_intervals, diff_equ_test, &init_vec);

    let mut init_vec = vec![Complex64::new(0.0, 0.0); 13];
    corr_init_vec_ada_order_2__2e_2::corr_init_vec_ada_order_2__2e_2(&mut init_vec);
    let t_start = 0.0;
    let t_end: f64 = 10.0;
    let results = rkf45(t_start, t_end, corr_sys_ada_order_2::corr_system_ada_order_2, &init_vec);
    //let init_vec = vec![Complex64::new(1.0, 0.0)];
    //let t_start = 0.0;
    //let t_end: f64 = 10.0;
    //let results = rkf45(t_start, t_end, diff_equ_test, &init_vec);
    println!("Done");
    println!("Saving data...");
    save_data("corr_ada_2.dat", &results.0, &results.1);
    println!("Done");
}
