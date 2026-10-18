#![no_std]
use core::ops::Sub;

use generic_array::{
    ArrayLength, GenericArray,
    typenum::{B1, Sub1, U3},
};
use libm::sqrtf;
use num::{Complex, complex::c32};

#[derive(Debug, PartialEq)]
pub struct Polynomial<N: ArrayLength>(GenericArray<f32, N>);

impl<N: ArrayLength> Polynomial<N>
where
    N: Sub<B1>,
    Sub1<N>: ArrayLength,
{
    pub fn new(p: GenericArray<f32, N>) -> Self {
        Self(p)
    }
    pub fn der(&self) -> Polynomial<Sub1<N>> {
        let mut d = GenericArray::<f32, Sub1<N>>::default();
        for (i, k) in d.iter_mut().enumerate() {
            *k = self.0[i] * (i + 1) as f32;
        }
        Polynomial(d)
    }
    pub fn eval_re(&self, mut x: f32) -> f32 {
        let mut s = self.0[0];
        for k in self.0.iter().skip(1) {
            s += x * k;
            x *= x;
        }
        s
    }
    pub fn eval(&self, mut x: Complex<f32>) -> Complex<f32> {
        let mut s = self.0[0].into();
        for k in self.0.iter().skip(1) {
            s += x * k;
            x *= x;
        }
        s
    }
    pub fn newton_method(
        &self,
        d: &Polynomial<Sub1<N>>,
        mut x: Complex<f32>,
        n: usize,
    ) -> Complex<f32> {
        for _ in 0..n {
            x -= self.eval(x) / d.eval(x);
        }
        x
    }
}

impl Polynomial<U3> {
    pub fn roots(&self) -> [Complex<f32>; 2] {
        let d = self.0[1] * self.0[1] - 4.0 * self.0[0] * self.0[2];
        let sqrt_d = sqrtf(d.abs());
        if d < 0.0 {
            [
                (-self.0[1] + c32(0.0, sqrt_d)) / (2.0 * self.0[0]),
                (-self.0[1] - c32(0.0, sqrt_d)) / (2.0 * self.0[0]),
            ]
        } else {
            [
                c32((-self.0[1] + sqrt_d) / (2.0 * self.0[0]), 0.0),
                c32((-self.0[1] - sqrt_d) / (2.0 * self.0[0]), 0.0),
            ]
        }
    }
}

#[cfg(test)]
mod tests {
    use generic_array::arr;

    use crate::Polynomial;
    extern crate std;
    use std::dbg;

    #[test]
    fn test_der() {
        let p = Polynomial::new(arr![1.0, 1.0, 1.0]);
        assert_eq!(p.der(), Polynomial(arr![1.0, 2.0]));
    }
    #[test]
    fn test_eval_re() {
        let p = Polynomial(arr![1.0, 1.0, 1.0]);
        assert_eq!(p.eval_re(0.0), 1.0);
        assert_eq!(p.eval_re(1.0), 3.0);
        assert_eq!(p.eval_re(2.0), 7.0);
    }
    #[test]
    fn test_roots_quadratic() {
        let p = Polynomial::new(arr![1.0, 3.0, 2.0]);
        // assert_eq!(p.der(), Polynomial(arr![1.0, 2.0]));
        let r = p.roots();
        dbg!(r);
    }
}
