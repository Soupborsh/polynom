#![no_std]
use core::ops::Sub;

use generic_array::{
    ArrayLength, GenericArray,
    typenum::{B1, Sub1, U3},
};
use libm::sqrtf;
use num::Complex;

#[derive(Debug, PartialEq)]
pub struct Polynomial<N: ArrayLength>(GenericArray<Complex<f32>, N>);

impl<N: ArrayLength> Polynomial<N> {
    pub fn new(p: GenericArray<Complex<f32>, N>) -> Self {
        Self(p)
    }
    pub fn eval(&self, x: Complex<f32>) -> Complex<f32> {
        let mut s = self.0[0];
        let mut xp = x;
        for k in self.0.iter().skip(1) {
            s += xp * k;
            xp *= x;
        }
        s
    }
}
// #[cfg(test)]
// extern crate std;

impl<N: ArrayLength> Polynomial<N>
where
    N: Sub<B1>,
    Sub1<N>: ArrayLength,
{
    pub fn der(&self) -> Polynomial<Sub1<N>> {
        let mut d = GenericArray::<Complex<f32>, Sub1<N>>::default();
        for (i, k) in d.iter_mut().enumerate() {
            *k = self.0[i + 1] * (i + 1) as f32;
        }
        Polynomial(d)
    }
    // Gorner
    // pub fn div_by_root(&self, x: Complex<f32>) -> Polynomial<Sub1<N>> {
    //     let mut d = GenericArray::<f32, Sub1<N>>::default();
    //     *d.last_mut().unwrap() = *self.0.last().unwrap();
    //     for (i, k) in d.iter_mut().rev().skip(1).enumerate() {
    //         *k = self.0[i] + d[i - 1] * x;
    //     }
    //     // for (i, k) in d.iter_mut().enumerate() {
    //     //     *k = self.0[i + 1] * (i + 1) as f32;
    //     // }
    //     Polynomial(d)
    // }
    pub fn newton_method(
        &self,
        d: &Polynomial<Sub1<N>>,
        mut x: Complex<f32>,
        n: usize,
    ) -> Complex<f32> {
        for _ in 0..n {
            // #[cfg(test)]
            // std::dbg!(x);
            x -= self.eval(x) / d.eval(x);
        }
        x
    }
}

impl Polynomial<U3> {
    pub fn roots(&self) -> [Complex<f32>; 2] {
        let d: Complex<f32> = self.0[1] * self.0[1] - 4.0 * self.0[0] * self.0[2];
        let sqrt_d = d.sqrt();
        [
            (-self.0[1] + sqrt_d) / (2.0 * self.0[0]),
            (-self.0[1] - sqrt_d) / (2.0 * self.0[0]),
        ]
        // if d < 0.0 {
        //     [
        //         (-self.0[1] + c32(0.0, sqrt_d)) / (2.0 * self.0[0]),
        //         (-self.0[1] - c32(0.0, sqrt_d)) / (2.0 * self.0[0]),
        //     ]
        // } else {
        //     [
        //         c32((-self.0[1] + sqrt_d) / (2.0 * self.0[0]), 0.0),
        //         c32((-self.0[1] - sqrt_d) / (2.0 * self.0[0]), 0.0),
        //     ]
        // }
    }
}

pub fn c32(re: f32, im: f32) -> Complex<f32> {
    Complex { re, im }
}
pub fn cr(re: f32) -> Complex<f32> {
    Complex { re, im: 0.0 }
}
pub fn ci(im: f32) -> Complex<f32> {
    Complex { re: 0.0, im }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use generic_array::arr;

    use crate::{Polynomial, c32, cr};
    extern crate std;
    use std::dbg;

    #[test]
    fn test_der() {
        let p = Polynomial::new(arr![cr(1.0), cr(1.0), cr(1.0)]);
        assert_eq!(p.der(), Polynomial(arr![cr(1.0), cr(2.0)]));
    }
    // #[test]
    // fn test_eval_re() {
    //     let p = Polynomial(arr![1.0, 1.0, 1.0]);
    //     assert_eq!(p.eval_re(0.0), 1.0);
    //     assert_eq!(p.eval_re(1.0), 3.0);
    //     assert_eq!(p.eval_re(2.0), 7.0);
    //     let p = Polynomial(arr![1.0, 1.0, 1.0, 1.0]);
    //     assert_eq!(p.eval_re(-1.0), 0.0);
    // }
    #[test]
    fn test_roots_quadratic() {
        let p = Polynomial::new(arr![cr(1.0), cr(3.0), cr(2.0)]);
        let r = p.roots();
        dbg!(r);
    }
    #[test]
    fn test_newton_method() {
        let p = Polynomial::new(arr![cr(1.0), cr(1.0), cr(1.0), cr(1.0)]);
        let d = p.der();
        assert_relative_eq!(
            p.newton_method(&d, c32(1.0, -2.0f32), 32),
            c32(0.0, -1.0f32)
        );
    }
}
