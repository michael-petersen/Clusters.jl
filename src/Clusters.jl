module Clusters

# add any `using` statements here
using Optim           # for zero-finding I
using Roots           # for zero finder II
using Interpolations  # for basic interpolation

using HDF5            # To be able to dump/read HDF5 files


# define any global constants here
const astronomicalG = 1

# add libraries here,
# and expose whatever we'd like to see externally.

# basic utility libraries for functions
include("utils/coordinates.jl")         # tools for mapping (rperi,rapo)<->(a,e)
include("utils/conserved.jl")           # tools for mapping (rperi,rapo)<->(E,L)
export E_from_rpra_pot,L_from_rpra_pot
include("utils/bisect.jl")              # tools for stable bisection
export extremise_function
include("utils/mappings.jl")            # tools for mapping (alpha,beta)<->(u,v)

# functions for numerical computation of frequencies
include("frequencies/epicycle.jl")      # frequencies (Ω₁,Ω₂) for circular orbits
export Omega1_circular,Omega2_circular
include("frequencies/frequencies.jl")   # frequencies for generic orbits (a,e)<->(Ω₁,Ω₂)
export compute_frequencies_henon_ae

# functions to return basis elements
include("basis/basis.jl")
export potential_function,density_function,accumulate_ln
include("basis/rawCB73.jl")
export tabUlnpCB73!

# functions to analyse N-body runs
include("nbody/nbody.jl")
export return_particles,return_density_centre,find_rbary,perturbation_density,discrete_fourier_k

# the isochrone model
include("models/isochrone.jl")
export isochrone_psi,isochrone_dpsi_dr,isochrone_ddpsi_ddr

# the plummer model
include("models/plummer.jl")
export plummer_psi

# some utilities
include("utils/potentials.jl")
export dpotential_numerical,extremise_function

end # module
