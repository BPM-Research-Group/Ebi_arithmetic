pub trait EbiNumber: Zero + One + Infinite + Normal + Round + Clone {}

pub trait Zero: Sized {
    fn zero() -> Self;

    fn set_zero(&mut self) {
        *self = Zero::zero();
    }

    fn is_zero(&self) -> bool;
}

pub trait One: Sized {
    fn one() -> Self;

    fn set_one(&mut self) {
        *self = One::one();
    }

    fn is_one(&self) -> bool;
}

pub trait Signed: Sized {
    fn abs(&self) -> Self;

    /// Returns true if the number is positive and false if the number is zero or negative.
    fn is_positive(&self) -> bool;

    /// Returns true if the number is negative and false if the number is zero or positive.
    fn is_negative(&self) -> bool;

    /// For exact arithmetic: Returns true if the number is positive or zero.
    /// For approximate arithmetic: returns true if the number is larger than -epsilon
    fn is_not_negative(&self) -> bool {
        !self.is_negative()
    }

    /// For exact arithmetic: Returns true if the number is negative or zero.
    /// For approximate arithmetic: returns true if the number is smaller than epsilon
    fn is_not_positive(&self) -> bool {
        !self.is_positive()
    }
}

pub trait Infinite: Signed {
    /// Returns false if this value is positive infinity or negative infinity, and true otherwise.
    fn is_finite(&self) -> bool {
        !self.is_infinite()
    }

    /// Returns true if this value is positive infinity or negative infinity, and false otherwise.
    fn is_infinite(&self) -> bool;

    fn is_positive_infinite(&self) -> bool {
        self.is_infinite() && self.is_positive()
    }

    fn is_negative_infinite(&self) -> bool {
        self.is_infinite() && self.is_negative()
    }
}

pub trait Normal: Infinite + Signed + Zero + Sized {
    /// Returns true if the number is neither zero, infinite, subnormal, or NaN.
    fn is_normal(&self) -> bool {
        !self.is_zero() && !self.is_infinite() && !self.is_nan()
    }

    fn is_nan(&self) -> bool;
}

pub trait Round: Sized {
    /// Returns the largest integer less than or equal to `self`.
    fn floor(self) -> Self;

    /// Returns the smallest integer greater than or equal to `self`.
    fn ceil(self) -> Self;
}

pub trait Fractional: Sized {
    /// Takes the reciprocal (inverse) of a number, `1/x`.
    fn recip(&self) -> Self;
}
