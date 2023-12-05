use num::complex::{Complex, Complex64};
use std::fs::File;
use std::io::Write;

type Cvecf = Vec<Complex<f64>>;

fn rk4(t_start: f64, t_end: f64, num_intervals: usize, func: fn(f64, &Cvecf, &mut Cvecf), init_vec: &Cvecf) -> Vec<Cvecf> {
    let h = (t_end - t_start) / num_intervals as f64;
    let mut result = vec![vec![Complex::new(0.0, 0.0); num_intervals + 1]; init_vec.len()];
    let mut curr_vec = init_vec.clone();
    let mut curr_t = t_start;
    let mut k1 = vec![Complex::new(0.0, 0.0); init_vec.len()];
    let mut k2 = vec![Complex::new(0.0, 0.0); init_vec.len()];
    let mut k3 = vec![Complex::new(0.0, 0.0); init_vec.len()];
    let mut k4 = vec![Complex::new(0.0, 0.0); init_vec.len()];
    let mut temp = vec![Complex::new(0.0, 0.0); init_vec.len()];
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
        for j in 0..init_vec.len() {
            result[j][i] = curr_vec[j];
            //println!("{}", h / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]));
            curr_vec[j] += h / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
        }
        curr_t += h;
    }
    return result;
}

fn save_data(path: &str, t_start: f64, t_end: f64, num_intervals: usize, result: &Vec<Cvecf>) {
    let mut output =  File::create(path).unwrap();
    let mut time_output = String::new();
    let h = (t_end - t_start) / num_intervals as f64;
    for i in 0..(num_intervals + 1) {
        time_output += &(t_start + i as f64 * h).to_string();
        time_output += " ";
    }
    writeln!(output, "{}", time_output).unwrap();
    for ind_cvecf in result {
        //println!("{}", ind_cvecf.len());
        let mut line = String::new();
        for i in 0..(num_intervals + 1) {
            line += &(ind_cvecf[i].re.to_string() + "+" + &ind_cvecf[i].im.to_string() + "j");
            line += " ";
        }
        writeln!(output, "{}", line).unwrap();
    }
}

fn diff_equ_test(time: f64, yn: &Cvecf, fty: &mut Cvecf) {
    fty[0] = time.sin().powi(2) * yn[0];
}

fn main() {
    let a: Complex<i32> = Complex::new(10, 20);
    //println!("Hello, world! {}", &(a.re.to_string() + "+" + &a.im.to_string() + "j"));
    print!("Computing...");
    let init_vec = vec![Complex64::new(1.0, 0.0)];
    let t_start = 0.0;
    let t_end: f64 = 5.0;
    let num_intervals = 24;
    let result = rk4(t_start, t_end, num_intervals, diff_equ_test, &init_vec);
    println!("Done");
    print!("Saving data...");
    save_data("diff_equ_test.dat", t_start, t_end, num_intervals, &result);
    println!("Done");
}
