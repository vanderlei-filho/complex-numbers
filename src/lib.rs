//! Simple complex numbers library written for educational purposes.
//! Instead of using generic types, this library is made to be used
//! only with f64 primitives in a seamless manner.
use std::ops::{Add, Sub, Mul, Div};
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
    /// the conjugate complex of the input argument
    fn conjugate(&self) -> ComplexNumber {
        ComplexNumber {re: self.re, im: -self.im,}
    }

    /// Returns a ComplexNumber struct representing 
    /// the imaginary unit number "i"
    fn i() -> ComplexNumber {
        ComplexNumber {re: 0f64, im: 1f64}
    }

    /// Checks if the ComplexNumber is zero
    fn is_zero(&self) -> bool {
        if self.re == 0f64 && self.im == 0f64 {
            true
        } else {
            false
        }
    }
    // Converts self into a PolarComplex struct
    /*fn to_polar(&self) -> PolarComplex {
        PolarComplex { arg: self.argument(), mag: self.magnitude() }
    }*/
}

impl Add<Self> for ComplexNumber {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}

/// Add traits for ComplexNumber to support
/// addition with the f64 primitive
impl Add<f64> for ComplexNumber {
    type Output = Self;

    fn add(self, rhs: f64) -> Self {
        Self {
            re: self.re + rhs,
            im: self.im
        }
    }
}

impl Add<ComplexNumber> for f64 {
    type Output = ComplexNumber;

    fn add(self, rhs: ComplexNumber) -> ComplexNumber {
        ComplexNumber {
            re: self + rhs.re,
            im: rhs.im
        }
    }
}

impl Sub<Self> for ComplexNumber {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self {
            re: self.re - rhs.re,
            im: self.im - rhs.im,
        }
    }
}

/// Sub trait for ComplexNumber to support
/// subtraction with the f64 primitive
impl Sub<f64> for ComplexNumber {
    type Output = Self;

    fn sub(self, rhs: f64) -> Self {
        Self {
            re: self.re - rhs,
            im: self.im
        }
    }
}

impl Sub<ComplexNumber> for f64 {
    type Output = ComplexNumber;

    fn sub(self, rhs: ComplexNumber) -> ComplexNumber {
        ComplexNumber {
            re: self - rhs.re,
            im: (-rhs.im)
        }
    }
}

/// A complex number is equal to anrhs complex
/// number if and only if their real and imaginary
/// parts are equal
impl PartialEq<Self> for ComplexNumber {
    fn eq(&self, rhs: &Self) -> bool {
        self.re == rhs.re && 
        self.im == rhs.im
    }
}

impl PartialEq<ComplexNumber> for f64 {
    fn eq(&self, rhs: &ComplexNumber) -> bool {
        self == &rhs.re && rhs.im == 0f64
    }
}

impl PartialEq<f64> for ComplexNumber {
    fn eq(&self, rhs: &f64) -> bool {
        &self.re == rhs && self.im == 0f64
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
    fn mul(self, rhs: Self) -> Self {
        Self {
            re: self.re*rhs.re - self.im*rhs.im,
            im: self.re*rhs.im + self.im*rhs.re
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
    fn mul(self, rhs: f64) -> Self {
        Self {
            re: self.re*rhs,
            im: self.im*rhs
        }
    }
}

impl Mul<ComplexNumber> for f64 {
    type Output = ComplexNumber;
    fn mul(self, rhs: ComplexNumber) -> ComplexNumber {
        ComplexNumber {
            re: self*rhs.re,
            im: self*rhs.im
        }
    }
}

impl Div<Self> for ComplexNumber {
    type Output = Self;
    /// The division operation between complex numbers
    /// is implemented based on the following derived
    /// formula: w = u + vi; z = x + yi;
    /// => w/y = (1/(x^2+y^2)) * ((ux + vy) + (vx - uy)i)
    /// See: https://en.wikipedia.org/wiki/Complex_number#Reciprocal_and_division
    fn div(self, rhs: Self) -> Self {
        let norm_sqr = 1f64/(rhs.re.powf(2f64)+rhs.im.powf(2f64));
        Self {
            re: norm_sqr*(self.re*rhs.re + self.im*rhs.im),
            im: norm_sqr*(self.im*rhs.re - self.re*rhs.im)
        }
    }
}

impl Div<f64> for ComplexNumber {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Self {
            re: self.re / rhs,
            im: self.im / rhs
        }
    }
}

impl Div<ComplexNumber> for f64 {
    type Output = ComplexNumber;
    fn div(self, rhs: ComplexNumber) -> ComplexNumber {
        let norm_sqr = 1f64/(rhs.re.powf(2f64)+rhs.im.powf(2f64));
        ComplexNumber {
            re: norm_sqr*(self*rhs.re),
            im: norm_sqr*(-self*rhs.im)
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
        assert_eq!(z1 - z2, ComplexNumber { re: -1f64, im: 2f64 })
    }

    #[test]
    fn test_add_traits() {
        let z = ComplexNumber::i();
        assert!(z.re == 0.0f64);
        assert!(z.im == 1.0f64);
        let sum = z + 3f64;
        assert_eq!(sum, ComplexNumber { re:3.0f64, im:1.0f64 });
        assert_eq!(z + 3f64, 3f64 + z);
        let z2 = ComplexNumber::i() + 3.0f64;
        assert_eq!(z2.re, 3.0f64);
    }

    #[test]
    fn test_sub_traits() {
        let z = ComplexNumber::i();
        assert_eq!(z - 3f64, ComplexNumber { re:-3.0f64, im:1.0f64 });
        assert_eq!(3f64 - z, ComplexNumber { re: 3f64, im:-1.0f64 });
        let z1 = ComplexNumber { re: 4f64, im: 2f64 };
        let z2 = ComplexNumber { re: 2f64, im: -2f64 };
        assert_eq!(z1 - z2, ComplexNumber { re: 2f64, im: 4f64 });
        assert_eq!(z2 - z1, ComplexNumber { re: -2f64, im: -4f64 })
    }

    #[test]
    fn test_mul_traits() {
        let z1 = ComplexNumber::i();
        assert_eq!(z1*4f64, ComplexNumber { re:0.0f64, im:4.0f64 });
        assert_eq!(4f64*z1, ComplexNumber { re:0.0f64, im:4.0f64 });
        let z2 = ComplexNumber { re: 2.0f64.sqrt(), im:-2.0f64 };
        assert_eq!(z1*z2, ComplexNumber { re: 2.0f64, im: 2.0f64.sqrt() });
        assert_eq!(z1*z2, z2*z1);
        assert_eq!(z1*z2*3f64, 3f64*z2*z1);
        // composed arithmetic
        assert_eq!(z1*(3f64 + 2f64), ComplexNumber { re: 0f64, im: 5f64 });
        assert_eq!((z1+z2)*z2, ComplexNumber {re: 2.0f64.sqrt(), im: -1f64} * z2)
    }

    #[test]
    fn test_div_traits() {
        let z1 = ComplexNumber::i();
        assert_eq!(z1/4f64, ComplexNumber { re:0.0f64, im:0.25f64 });
        dbg!(z1/z1);
        assert_eq!(z1/z1, 1f64);
        let z2 = ComplexNumber { re: 2.0f64, im: 2.0f64.sqrt() };
        let z3 = ComplexNumber { re: 1.0f64, im: 3.0f64.sqrt() };
        dbg!(z3/z2);
        dbg!(4f64/z1);
        assert!(true)
    }

}
