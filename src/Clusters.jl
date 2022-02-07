module Clusters

# add any `using` statements here
using Optim           # for zero-finding I
using Roots           # for zero finder II
using Interpolations  # for basic interpolation

# define any global constants here
const astronomicalG = 1

# add libraries here,
# and expose whatever we'd like to see externally.

# basic libraries for functions
include("utils/coordinates.jl")
include("utils/conserved.jl")
export E_from_rpra_pot,L_from_rpra_pot

# functions for numerical computation of frequencies
# first for circular orbits
include("frequencies/epicycle.jl")
# no exports needed? can do beta_c for fun visualisations
export make_betac

# then for more generic orbits
include("frequencies/frequencies.jl")
export compute_frequencies_henon_ae

# functions to return basis elements
include("basis/basis.jl")
export potential_function,density_function,accumulate_ln

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
