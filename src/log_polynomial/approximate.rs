use crate::{
    fraction::approximate::Approximate,
    log_polynomial::{
        log_polynomial_enum::LogPolynomialEnum, log_polynomial_exact::LogPolynomialExact,
        log_polynomial_f64::LogPolynomialF64,
    },
};
use anyhow::{Result, anyhow};

impl Approximate for LogPolynomialExact {
    fn approximate(self) -> Result<f64> {
        let mut sum = 0.0;
        for (argument, coefficient) in self.argument2coefficient {
            sum += (argument.approx_ln() / 0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964186875420014810205706857336855202) * coefficient.approximate()?;
        }
        Ok(sum)
    }
}

impl Approximate for LogPolynomialF64 {
    fn approximate(self) -> Result<f64> {
        Ok(self.0)
    }
}

impl Approximate for LogPolynomialEnum {
    fn approximate(self) -> Result<f64> {
        match self {
            LogPolynomialEnum::Approx(log_polynomial_f64) => log_polynomial_f64.approximate(),
            LogPolynomialEnum::Exact(log_polynomial_exact) => log_polynomial_exact.approximate(),
            LogPolynomialEnum::CannotCombineExactAndApprox => {
                Err(anyhow!("Cannot combine approximate and exact arithmetic."))
            }
        }
    }
}
