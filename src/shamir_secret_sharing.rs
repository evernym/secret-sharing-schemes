// Shamir secret sharing

use amcl_wrapper::field_elem::FieldElement;
use std::collections::{HashMap, HashSet};
use crate::polynomial::Polynomial;

/// Generate a random polynomial with the secret at the polynomial evaluation at 0.
pub fn get_shared_secret_with_polynomial(
    threshold: usize,
    total: usize,
) -> (FieldElement, HashMap<usize, FieldElement>, Polynomial) {
    let random_poly = Polynomial::random(threshold - 1);
    let secret = random_poly.eval(&FieldElement::zero());
    let shares = (1..=total)
        .map(|x| (x, random_poly.eval(&FieldElement::from(x as u64))))
        .collect::<HashMap<usize, FieldElement>>();
    (secret, shares, random_poly)
}

/// Generate a secret with its shares according to Shamir secret sharing.
/// Returns the secret and a map of share_id -> share
pub fn get_shared_secret(
    threshold: usize,
    total: usize,
) -> (FieldElement, HashMap<usize, FieldElement>) {
    let (secret, shares, _) = get_shared_secret_with_polynomial(threshold, total);
    (secret, shares)
}

pub fn reconstruct_secret(threshold: usize, shares: HashMap<usize, FieldElement>) -> FieldElement {
    assert!(shares.len() >= threshold);
    let mut secret = FieldElement::zero();
    let share_ids = shares
        .iter()
        .take(threshold)
        .map(|(i, _)| *i)
        .collect::<HashSet<usize>>();
    for id in share_ids.clone() {
        let share = shares.get(&id).unwrap();
        let l = Polynomial::lagrange_basis_at_0(share_ids.clone(), id);
        secret += &(&l * share)
    }
    secret
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_secret_sharing() {
        for _ in 0..10 {
            let threshold = 5;
            let total = 10;
            let (secret, shares) = get_shared_secret(threshold, total);
            assert_eq!(shares.len(), total);
            let recon_secret = reconstruct_secret(
                threshold,
                shares
                    .into_iter()
                    .take(threshold)
                    .collect::<HashMap<usize, FieldElement>>(),
            );
            assert_eq!(secret, recon_secret);
        }
    }

    #[test]
    fn test_secret_sharing_1() {
        {
            let threshold = 5;
            let total = 10;
            let (secret, shares) = get_shared_secret(threshold, total);
            let mut some_shares = HashMap::<usize, FieldElement>::new();
            for i in vec![1, 3, 4, 7, 9] {
                some_shares.insert(i, shares.get(&i).unwrap().clone());
            }
            let recon_secret = reconstruct_secret(threshold, some_shares);
            assert_eq!(secret, recon_secret);
        }

        {
            let threshold = 3;
            let total = 5;
            let (secret, shares) = get_shared_secret(threshold, total);
            let mut some_shares = HashMap::<usize, FieldElement>::new();
            for i in vec![1, 2, 4] {
                some_shares.insert(i, shares.get(&i).unwrap().clone());
            }
            let recon_secret = reconstruct_secret(threshold, some_shares);
            assert_eq!(secret, recon_secret);
        }

        {
            let threshold = 2;
            let total = 5;
            let (secret, shares) = get_shared_secret(threshold, total);
            let mut some_shares = HashMap::<usize, FieldElement>::new();
            for i in vec![1, 4] {
                some_shares.insert(i, shares.get(&i).unwrap().clone());
            }
            let recon_secret = reconstruct_secret(threshold, some_shares);
            assert_eq!(secret, recon_secret);
        }

        {
            let threshold = 3;
            let total = 5;
            let (secret, shares) = get_shared_secret(threshold, total);
            let mut some_shares = HashMap::<usize, FieldElement>::new();
            for i in vec![1, 2, 4, 5] {
                some_shares.insert(i, shares.get(&i).unwrap().clone());
            }
            let recon_secret = reconstruct_secret(threshold, some_shares);
            assert_eq!(secret, recon_secret);
        }
    }
}