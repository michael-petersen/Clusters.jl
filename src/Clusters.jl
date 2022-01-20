module Clusters

# add any `using` statements here
using Optim

# define any global constants here
const astronomicalG = 1

# add libraries here,
# and expose whatever we'd like to see externally.

# basic libraries for functions
include("utils/coordinates.jl")

# functions for numerical computation of frequencies
include("frequencies/frequencies.jl")
export compute_frequencies,compute_frequencies_EK

# functions to return basis elements
include("basis/basis.jl")
export potential_function,density_function,accumulate_ln

# functions to analyse N-body runs
include("nbody/nbody.jl")
export return_particles,return_density_centre,find_rbary,perturbation_density,discrete_fourier_k

# the isochrone model
include("models/isochrone.jl")
export isochrone_psi,isochrone_Omega_1_2_ae

# the plummer model
include("models/plummer.jl")
export plummer_psi

# some utilities
include("utils/potentials.jl")
export dpotential_numerical,extremise_function

end # module
