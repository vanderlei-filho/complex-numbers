//! Simple complex numbers library written for educational purposes.
//! Instead of using generic types, this library is made to be used
//! only with f64 primitives.
use std::ops::{Add, Sub, Mul};
use std::cmp::PartialEq;
use std::f64::consts::PI;


/// A complex number in Cartesian form.
#[derive(Copy, Clone, Debug)]
pub struct ComplexNumber {
    /// Real part of the complex number
    pub re: f64,
    /// Imaginary part of the complex number  
    pub im: f64
}

impl ComplexNumber {
    /// Returns the argument of a complex number
    /// in radians. The atan2() function is used
    /// to handle results based on quadrants
    fn argument(&self) -> f64 {
        self.im.atan2(self.re)
    }
    /// Returns the magnitude of a complex number
    fn magnitude(&self) -> f64 {
        //(self.re.powf(2f64) + self.im.powf(2f64)).sqrt()
        self.re.hypot(self.im)
    }
    /// Returns a ComplexNumber struct representing
    /// the conjugate complex of the input argument.
    fn conjugate(&self) -> ComplexNumber {
        ComplexNumber {re: self.re, im: -self.im,}
    }
    /// Returns a ComplexNumber struct representing 
    /// the imaginary number "i"
    fn i() -> ComplexNumber {
        ComplexNumber {re: 0f64, im: 1f64}
    }
    // Converts self into a PolarComplex struct
    /*fn to_polar(&self) -> PolarComplex {
        PolarComplex { arg: self.argument(), mag: self.magnitude() }
    }*/
}

impl Add for ComplexNumber {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }
}

/// Add trait for ComplexNumber to support
/// addition with the f64 primitive (or "real"
/// numbers)
impl Add<f64> for ComplexNumber {
    type Output = Self;

    fn add(self, other: f64) -> Self {
        Self {
            re: self.re + other,
            im: self.im
        }
    }
}

impl Sub for ComplexNumber {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }
}

/// Add trait for ComplexNumber to support
/// subtraction with the f64 primitive
impl Sub<f64> for ComplexNumber {
    type Output = Self;

    fn sub(self, other: f64) -> Self {
        Self {
            re: self.re - other,
            im: self.im
        }
    }
}

/// A complex number is equal to another complex
/// number if and only if their real and imaginary
/// parts are equal
impl PartialEq for ComplexNumber {
    fn eq(&self, other: &Self) -> bool {
        self.re == other.re && 
        self.im == other.im
    }
}

impl Mul for ComplexNumber {
    type Output = Self;
    /// complex numbers multiplication work based on
    /// distributive rules for common numbers:
    /// z1 = a + i*b; z2 = c + i*d;
    /// z1*z2 = (a+i*b)*(c+i*d)
    ///       = a*c + i*a*d + i*b*c + (i^2)*b*d
    ///       = a*c + i*(a*d+b*c) - b*d
    ///       = (a*c-b*d) + i*(a*d+b*c)
    fn mul(self, other: Self) -> Self {
        Self {
            re: self.re*other.re - self.im*other.im,
            im: self.re*other.im + self.im*other.re
        }
    }
}

impl Mul<f64> for ComplexNumber {
    type Output = Self;
    /// Multiplying a complex number by a real
    /// number is an operation that is also called
    /// "scaling" the complex number. It works the
    /// same way as the multiplication between two
    /// complex numbers, but resulting in a simpler
    /// distribution rule:
    /// z1 = a + i*b; r1 = c;
    /// z1*r1 = c*(a+i*b)
    ///       = (c*a) + i*(c*b)
    fn mul(self, other: f64) -> Self {
        Self {
            re: self.re*other,
            im: self.im*other
        }
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cartesian() {
        // Assert new function
        let z: ComplexNumber = ComplexNumber::i();
        // Assert argument function
        dbg!(z.argument());
        // Assert magnitude function
        dbg!(z.magnitude());
        // Assert conjugate function
        assert!(z.conjugate().re == 0f64);
        assert!(z.conjugate().im == -1f64);
        // Assert Add and Sub operator
        let z1 = ComplexNumber { re: 0f64, im: 1f64 };
        let z2 = ComplexNumber { re: 1f64, im: -1f64 };
        let z3 = z1 + z2;
        assert_eq!(z3, ComplexNumber { re: 1f64, im: 0f64 }); 
        assert!(z1 - z2 == ComplexNumber { re: -1f64, im: 2f64 })
    }

    #[test]
    fn test_sum_traits() {
        let z = ComplexNumber::i();
        dbg!(z.re);
        dbg!(z.im);
        assert!(z.re == 0.0f64);
        assert!(z.im == 1.0f64);
        let sum = z + 3f64;
        dbg!(sum);
        assert!(sum == ComplexNumber { re:3.0f64, im:1.0f64 })
    }

    #[test]
    fn test_mul_traits() {
        let z1 = ComplexNumber::i();
        assert_eq!(z1*4f64, ComplexNumber { re:0.0f64, im:4.0f64 });
        let z2 = ComplexNumber { re: 2.0f64.sqrt(), im:-2.0f64 };
        assert_eq!(z1*z2, ComplexNumber { re: 2.0f64, im: 2.0f64.sqrt() });
        assert_eq!(z1*z2, z2*z1)
    }
}
