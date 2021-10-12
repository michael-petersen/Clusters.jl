
# Clusters.jl

`Clusters.jl` is a package written in Julia to analyse spherical
clusters, from analytical descriptions to _n_-body simulations.

-----------------------------

## Analysing _n_-body runs
#### reader abilities

`return_particles(filename)` reads in a delimited file to x,y,z,vx,vy,vz

`return_density_centre(filename)` reads in a centre file



#### centring support

`find_rbary(x,y,z)` computes the barycentre and returns the radius for a given x,y,z series
`find_rbary(x,y,z,xc,yc,zc)` offsets series by xc,yc,zc and returns the radius for a given x,y,z series

`crosscorrelation(x,y)` computes the cross correlation of two density centre series à la Heggie.

#### power spectra computation

`perturbation_density(x,y,z,rmax,m=1)` computes eq. 3 of Heggie (2020)

`discrete_fourier_k(x,k)` computes eq. 4 of Heggie (2020)

#### escaper analysis

`plummer_potential(r,a=1)` return the simplest possible Plummer potential

`escaper_energy(x,y,z,vx,vy,vz)` return the naive total energy (½v² + φ)

`spherical_coordinates(x,y,z,vx,vy,vz)` return the spherical positions and velocities

## Computing numerical frequencies




## Basis element support


## Author

Mike Petersen -  @michael-petersen - petersen@iap.fr



