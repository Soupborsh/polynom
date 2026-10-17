use std::ops::Sub;

use generic_array::{
    ArrayLength, GenericArray,
    typenum::{B1, Sub1},
};

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
    pub fn eval(&self, mut x: f32) -> f32 {
        let mut s = self.0[0];
        for k in self.0.iter().skip(1) {
            s += x * k;
            x *= x;
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use generic_array::arr;

    use crate::Polynomial;

    #[test]
    pub(crate) fn test_der() {
        let p = Polynomial::new(arr![1.0, 1.0, 1.0]);
        assert_eq!(p.der(), Polynomial(arr![1.0, 2.0]));
    }
    #[test]
    pub(crate) fn test_eval() {
        let p = Polynomial(arr![1.0, 1.0, 1.0]);
        assert_eq!(p.eval(0.0), 1.0);
        assert_eq!(p.eval(1.0), 3.0);
        assert_eq!(p.eval(2.0), 7.0);
    }
}
