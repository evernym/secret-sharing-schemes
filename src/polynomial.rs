use amcl_wrapper::field_elem::{FieldElement, FieldElementVector};
use std::collections::HashSet;

pub struct Polynomial(FieldElementVector);

impl Polynomial {
    /// Return a randomly chosen polynomial (each coefficient is randomly chosen) of degree `degree`.
    pub fn random(degree: usize) -> Self {
        Self(FieldElementVector::random(degree + 1)) // +1 for constant term
    }

    pub fn degree(&self) -> usize {
        self.0.len() - 1
    }

    /// Return coefficients starting from lowest degree term
    pub fn coefficients(&self) -> &FieldElementVector {
        &self.0
    }

    // Evaluate polynomial at given `x`
    pub fn eval(&self, x: &FieldElement) -> FieldElement {
        if x.is_zero() {
            self.coefficients()[0].clone()
        } else {
            let exp = FieldElementVector::new_vandermonde_vector(x, self.degree() + 1);
            self.0.inner_product(&exp).unwrap()
        }
    }

    /// Return the Lagrange basis polynomial at x = 0 given the x coordinates
    pub fn lagrange_basis_at_0(x_coords: HashSet<usize>, i: usize) -> FieldElement {
        let mut numerator = FieldElement::one();
        let mut denominator = FieldElement::one();
        let i_as_field_elem = FieldElement::from(i as u64);
        let neg_i = -i_as_field_elem; // -i
        for x in x_coords {
            if x == i {
                continue;
            }
            // numerator = numerator * x
            let x_as_field_elem = FieldElement::from(x as u64);
            numerator = &numerator * &x_as_field_elem;
            let x_minus_i = &x_as_field_elem + &neg_i;
            // denominator = denominator * (x - i)
            denominator = &denominator * &x_minus_i;
        }
        denominator.inverse_mut();
        // (x_coords[0]) * (x_coords[1]) * ... / ((x_coords[0] - i) * (x_coords[1] - i) * ...)
        numerator * denominator
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_poly() {
        for _ in 0..10 {
            let degree = 10;
            let poly = Polynomial::random(degree);
            assert_eq!(poly.degree(), degree);
            let coeffs = poly.coefficients();

            // Evaluation at 0 results in coefficient of constant term
            assert_eq!(poly.eval(&FieldElement::zero()), coeffs[0]);

            // Evaluation at 1 results in sum of all coefficients
            assert_eq!(poly.eval(&FieldElement::one()), coeffs.sum());
        }
    }
}