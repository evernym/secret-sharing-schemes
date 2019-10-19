# Secret sharing
<strong>The implemented schemes are prototypes and should not be used in production. 
They are used to assist in demonstrating threshold signatures (in other repos)</strong>.

1. [Shamir secret sharing (Requires trusted third party)](src/shamir_secret_sharing.rs)
1. [Pedersen verifiable secret sharing (Requires trusted third party)](src/pedersen_vss.rs)
1. [Pedersen decentralized verifiable secret sharing (Does not require a trusted party)](src/pedersen_dvss.rs)


# Pending:
1. More tests
2. Model Pedersen DVSS as state machine