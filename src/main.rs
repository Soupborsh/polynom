use std::io::stdin;

use generic_array::arr;
use num::complex::c32;
use polynom::{Polynomial, cr};

fn main() {
    println!("Type 4 coefficients, space separated:");
    let mut s = String::new();
    stdin().read_line(&mut s).unwrap();
    let a = s.split_whitespace();
    let mut b = a.map(|x| x.parse::<f32>().unwrap());
    let mut p = arr![cr(0.0);4];
    for i in 0..4 {
        p[i].re = b.next().unwrap();
    }
    let p = Polynomial::new(p);
    let d = p.der();
    println!("x1: {}", p.newton_method(&d, c32(0.0, 1.0_f32), 32))
}
