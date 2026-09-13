# Endpoint-conditioned forward-query certificates

These rational witnesses predate the stronger full-controller certificate in
`certificates/forward_queries/robust_VUVUV.json`. They restrict controllers to
perfect ideal Bell/Pauli recovery; do not import the newer score-selected scope
into these older results.

Run `python history/endpoint_forward_queries/verify_certificates.py` from any
working directory. All four certificates and their supplied reason-specific
mutations are checked using the standard library. Discovery arrays and solver
scripts remain in the original Forward_Query_Classification handoff; they are
not dependencies of this historical verifier.
