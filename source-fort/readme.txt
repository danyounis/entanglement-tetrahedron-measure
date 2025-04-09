tetra_main: Double-Nelder-Mead optimization, Max[x,y,z] and Min[11θ,7ϕ].

tetra_mini: Parses a binary file containing an array defining the (x,y,z) values, named
    Xeval-xxxxxxxxxx.dat where the x's in {0..9} serve as a unique identifier (UID),
    and evaluates the Nelder-Mead interior minimization over (11θ,7ϕ).
    Compile and move to: ../source-mixed; it is intended to be driven by Python.

tetra_adhoc1: Here, we've made the ansatz z->Infinity, in which case Δ_ij is no longer
    dependent on (y,z). Thus, we evaluate Min[11θ,7ϕ] of σ_12(x) for a discrete set of x points.
