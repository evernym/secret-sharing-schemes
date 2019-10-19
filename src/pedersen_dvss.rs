// Pedersen Decentralized Verifiable secret sharing

use amcl_wrapper::group_elem::GroupElement;
use amcl_wrapper::group_elem_g1::G1;
use amcl_wrapper::field_elem::{FieldElement, FieldElementVector};
use std::collections::HashMap;
use crate::pedersen_vss::PedersenVSS;


// Pedersen Decentralized Verifiable secret sharing. Based on the paper "Non-interactive and information-theoretic
// secure verifiable secret sharing", section 5. https://www.cs.cornell.edu/courses/cs754/2001fa/129.PDF
// Does not involve a trusted third party but assumes that all participants (and not just threshold) participate till the end.
// Even if one participant aborts, the protocol needs to be restarted. A workaround is for each participant to ignore the
// faulty participant's share essentially making it such that the faulty participant was never there.
/*
    n participants want to generate a shared secret s k-of-n manner
    Each of the n participants chooses a secret and runs a VSS for that secret in k-of-n manner. Say participant i chooses a secret s_i_0
    The shared secret s the becomes sum of secrets chosen by all n participants so s = s_1_0 + s_2_0 + s_3_0 + ... s_n_0
    After each of the n participants has successfully runs a VSS, they generate their corresponding share of s by adding
    their shares of each s_i_0 for i in 1 to n.
*/
// TODO: Model the code as state machine
pub struct PedersenDVSSParticipant {
    pub id: usize,
    pub secret: FieldElement,
    pub comm_coeffs: HashMap<usize, G1>,
    pub s_shares: HashMap<usize, FieldElement>,
    pub t_shares: HashMap<usize, FieldElement>,
    all_comm_coeffs: HashMap<usize, HashMap<usize, G1>>,
    all_shares: HashMap<usize, (FieldElement, FieldElement)>,
    // XXX: Should be in a different struct if the protocol is modelled as a state machine
    pub final_comm_coeffs: HashMap<usize, G1>,
    pub secret_share: FieldElement,
}

impl PedersenDVSSParticipant {
    /// Generates a new secret and verifiable shares of that secret for every participant
    pub fn new(id: usize, threshold: usize, total: usize, g: &G1, h: &G1) -> Self {
        let (secret, _, comm_coeffs, s_shares, t_shares) =
            PedersenVSS::deal(threshold, total, &g, &h);
        // TODO: As mentioned in the paper, there should be a signature from the participant for non-repudiation
        Self {
            id,
            secret,
            comm_coeffs,
            s_shares,
            t_shares,
            all_comm_coeffs: HashMap::new(),
            all_shares: HashMap::new(),
            final_comm_coeffs: HashMap::new(),
            secret_share: FieldElement::new(),
        }
    }

    /// Called by a participant when it receives a share from another participant with id `sender_id`
    pub fn received_share(
        &mut self,
        sender_id: usize,
        comm_coeffs: HashMap<usize, G1>,
        share: (FieldElement, FieldElement),
        threshold: usize,
        total: usize,
        g: &G1,
        h: &G1,
    ) {
        assert!(sender_id <= total);
        assert!(!self.all_comm_coeffs.contains_key(&sender_id));
        assert!(!self.all_shares.contains_key(&sender_id));
        // Verify received share
        assert!(PedersenVSS::verify_share(
            threshold,
            self.id,
            (&share.0, &share.1),
            &comm_coeffs,
            &g,
            &h
        ));
        self.all_comm_coeffs.insert(sender_id, comm_coeffs);
        self.all_shares.insert(sender_id, share);
    }

    /// Called by a participant when it has received shares from all participants. Computes the final
    /// share of the distributed secret
    pub fn compute_final_comm_coeffs_and_shares(
        &mut self,
        threshold: usize,
        total: usize,
        g: &G1,
        h: &G1,
    ) {
        assert_eq!(self.all_comm_coeffs.len(), total - 1);
        assert_eq!(self.all_shares.len(), total - 1);

        // Compute own share and commitment to coefficients of the distributed secret.
        for i in 0..threshold {
            // cm is the sum of coefficients of each signer's polynomial's ith degree term
            let mut cm = G1::identity();
            for j in 1..=total {
                if j != self.id {
                    cm += self.all_comm_coeffs[&j].get(&i).unwrap();
                } else {
                    cm += self.comm_coeffs.get(&i).unwrap();
                }
            }
            self.final_comm_coeffs.insert(i, cm);
        }

        let mut final_s_share = FieldElement::zero();
        let mut final_t_share = FieldElement::zero();
        for i in 1..=total {
            let (s, t) = if i != self.id {
                let tpl = &self.all_shares[&i];
                (&tpl.0, &tpl.1)
            } else {
                (&self.s_shares[&i], &self.t_shares[&i])
            };
            final_s_share += s;
            final_t_share += t;
        }

        // Verify computed share of the distributed secret
        assert!(PedersenVSS::verify_share(
            threshold,
            self.id,
            (&final_s_share, &final_t_share),
            &self.final_comm_coeffs,
            &g,
            &h
        ));

        self.secret_share = final_s_share;
    }
}

/// Create participants that take part in a decentralized secret sharing and perform the secret sharing.
#[cfg(test)]
pub fn share_secret_for_testing(
    threshold: usize,
    total: usize,
    g: &G1,
    h: &G1,
) -> Vec<PedersenDVSSParticipant> {
    let mut participants = vec![];

    // Each participant generates a new secret and verifiable shares of that secret for everyone.
    for i in 1..=total {
        let p = PedersenDVSSParticipant::new(i, threshold, total, g, h);
        participants.push(p);
    }

    // Every participant gives shares of its secret to others
    for i in 0..total {
        for j in 0..total {
            if i == j {
                continue;
            }
            let (id, comm_coeffs, (s, t)) = (
                participants[j].id.clone(),
                participants[j].comm_coeffs.clone(),
                (
                    participants[j].s_shares[&(i + 1)].clone(),
                    participants[j].t_shares[&(i + 1)].clone(),
                ),
            );

            let recv_p = &mut participants[i];
            recv_p.received_share(id, comm_coeffs, (s, t), threshold, total, g, h);
        }
    }

    // Every participant computes its share to the distributed secret.
    for i in 0..total {
        participants[i].compute_final_comm_coeffs_and_shares(threshold, total, g, h);
    }
    participants
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shamir_secret_sharing::reconstruct_secret;

    #[test]
    fn test_Pedersen_DVSS() {
        let threshold = 5;
        let total = 10;
        let (g, h) = PedersenVSS::gens("test".as_bytes());
        let participants = share_secret_for_testing(threshold, total, &g, &h);

        let mut expected_shared_secret = FieldElement::zero();
        for p in &participants {
            expected_shared_secret += &p.secret;
        }
        let mut shares = HashMap::new();
        for i in 0..threshold {
            shares.insert(participants[i].id, participants[i].secret_share.clone());
        }

        // Verify that the secret can be recomputed.
        let recon_secret = reconstruct_secret(threshold, shares);

        assert_eq!(expected_shared_secret, recon_secret);
    }
}