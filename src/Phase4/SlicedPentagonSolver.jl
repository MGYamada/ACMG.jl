"""
    SlicedPentagonSolver

Placeholder module. Pentagon slice generation will be rewired once the
Kitaev χ³ ansatz is validated by `KitaevComplex.verify_homotopy`.

Currently exports only the pentagon system helper (same as before).
"""
module SlicedPentagonSolver

using Oscar
using ..KitaevComplex
using ..PentagonEquations
using ..PentagonSolver

export get_pentagon_system_only

"""
    get_pentagon_system_only(Nijk, r) -> (R, pent_eqs, n)

Trivial wrapper around `PentagonEquations.get_pentagon_system`.
"""
function get_pentagon_system_only(Nijk::Array{Int,3}, r::Int)
    R, eqs, n = PentagonEquations.get_pentagon_system(Nijk, r)
    return R, eqs, n
end

end # module SlicedPentagonSolver
