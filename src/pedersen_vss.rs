// Pedersen Verifiable secret sharing

use amcl_wrapper::group_elem::{GroupElement, GroupElementVector};
use amcl_wrapper::group_elem_g1::{G1, G1Vector};
use amcl_wrapper::field_elem::{FieldElement, FieldElementVector};
use std::collections::HashMap;
use crate::shamir_secret_sharing::get_shared_secret_with_polynomial;

// Pedersen Verifiable secret sharing. Based on the paper "Non-interactive and information-theoretic
// secure verifiable secret sharing", section 4. https://www.cs.cornell.edu/courses/cs754/2001fa/129.PDF.
/* The basic idea is the following
    Dealer wants to share a secret s in k-of-n manner with n participants
    Dealer commits to secret s with randomness t so C_0 = C(s, t) = g^s.h^t
    Create polynomial F(x) = s + F_1.x + F_2.x^2 + ... F_{k-1}.x^{k-1} such that F(0) = s.
    Create polynomial G(x) = t + G_1.x + G_2.x^2 + ... G_{k-1}.x^{k-1} such that G(0) = t.
    Commits to coefficients as C_1 = C(F_1, G_1), C_2 = C(F_2, G_2),... till C_k = C(F_k, G_k), broadcast to all n participants
    Dealer sends (F(i), G(i)) to participant i
    Each participant verifies C(F(i), G(i)) = C_0 * C_1^i * C_2^{i^2} * ... C_{k-1}^{k-1}
*/
pub struct PedersenVSS {}

impl PedersenVSS {
    /// Generators used for commitment.
    // TODO: Wrap `gens` in a struct
    pub fn gens(label: &[u8]) -> (G1, G1) {
        // For NUMS.
        let g = G1::from_msg_hash(&[label, " : g".as_bytes()].concat());
        let h = G1::from_msg_hash(&[label, " : h".as_bytes()].concat());
        (g, h)
    }

    /// Executed by dealer. Output secret, blinding, commitment to coefficients of both polynomials
    /// and shares for each participant. Each participant has access to all commitments to coefficients
    /// but only to its own share.
    pub fn deal(
        threshold: usize,
        total: usize,
        g: &G1,
        h: &G1,
    ) -> (
        FieldElement,                 // secret
        FieldElement,                 // blinding
        HashMap<usize, G1>,           // commitment to coefficients
        HashMap<usize, FieldElement>, // shares for secret
        HashMap<usize, FieldElement>, // shares for blinding
    ) {
        let (s, s_shares, s_poly) = get_shared_secret_with_polynomial(threshold, total);
        let (t, t_shares, t_poly) = get_shared_secret_with_polynomial(threshold, total);
        // map of i -> g^s_poly.coefficients[i] * h^t_poly.coefficients[i]
        let commitment_coeffs = (0..threshold)
            .map(|i| {
                (
                    i,
                    // g^s_poly.coefficients[i] * h^t_poly.coefficients[i]
                    g.binary_scalar_mul(&h, &s_poly.coefficients()[i], &t_poly.coefficients()[i]),
                )
            })
            .collect::<HashMap<usize, G1>>();
        (s, t, commitment_coeffs, s_shares, t_shares)
    }

    /// Executed by each participant to verify its share received from the dealer.
    pub fn verify_share(
        threshold: usize,
        id: usize,
        share: (&FieldElement, &FieldElement),
        commitment_coeffs: &HashMap<usize, G1>,
        g: &G1,
        h: &G1,
    ) -> bool {
        assert!(commitment_coeffs.len() >= threshold);
        // Check commitment_coeffs[0] * commitment_coeffs[1]^id * commitment_coeffs[2]^{id^2} * ... commitment_coeffs[threshold-1]^{id^threshold-1} == g^share.0 * h^share.1
        // => commitment_coeffs[0] * commitment_coeffs[1]^id * commitment_coeffs[2]^{id^2} * ... commitment_coeffs[threshold-1]^{id^threshold-1} * {g^share.0 * h^share.1}^-1 == 1

        // exp will be [1, id, id^2, ... id^threshold-1]
        let mut exp =
            FieldElementVector::new_vandermonde_vector(&FieldElement::from(id as u64), threshold);

        // add share.0 and share.1 to exp
        exp.push(share.0.clone());
        exp.push(share.1.clone());

        let mut bases = G1Vector::with_capacity(threshold + 2);
        for i in 0..threshold {
            bases.push(commitment_coeffs[&i].clone())
        }

        // g^share.0 and h^share.1 will need to be inverted. To do one multi-scalar multiplication,invert g and h
        bases.push(g.negation());
        bases.push(h.negation());

        bases.multi_scalar_mul_var_time(&exp).unwrap().is_identity()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shamir_secret_sharing::reconstruct_secret;

    #[test]
    fn test_Pedersen_VSS() {
        let threshold = 5;
        let total = 10;
        let (g, h) = PedersenVSS::gens("test".as_bytes());
        let (secret, _, comm_coeffs, s_shares, t_shares) =
            PedersenVSS::deal(threshold, total, &g, &h);
        assert_eq!(s_shares.len(), total);
        assert_eq!(t_shares.len(), total);
        assert_eq!(comm_coeffs.len(), threshold);
        for i in 1..=total {
            assert!(PedersenVSS::verify_share(
                threshold,
                i,
                (&s_shares[&i], &t_shares[&i]),
                &comm_coeffs,
                &g,
                &h
            ));
        }
        let recon_secret = reconstruct_secret(
            threshold,
            s_shares
                .into_iter()
                .take(threshold)
                .collect::<HashMap<usize, FieldElement>>(),
        );
        assert_eq!(secret, recon_secret);
    }
}