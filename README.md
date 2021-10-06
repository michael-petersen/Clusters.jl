
# Clusters.jl

`Clusters.jl` is a package written in Julia to analyse spherical
clusters, from analytical descriptions to _n_-body simulations.

-----------------------------

## Analysing _n_-body runs
#### reader abilities

`return_particles(filename::String)` reads in a delimited file to x,y,z,vx,vy,vz

`return_density_centre(filename::String)` reads in a centre file

#### centring support

`find_rbary` computes the barycentre and returns the radius for a given x,y,z series

#### power spectra computation

`perturbation_density` computes eq. 3 of Heggie (2020)

`discrete_fourier_k` computes eq. 4 of Heggie (2020)

## Computing numerical frequencies


## Basis element support


## Author

Mike Petersen -  @michael-petersen - petersen@iap.fr



