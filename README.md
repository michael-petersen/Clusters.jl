
# Clusters.jl

`Clusters.jl` is a package written in Julia to analyse spherical
clusters, from analytical descriptions to _n_-body simulations.

-----------------------------

## Quick activate

In the main directory where you the package lives, enter the Julia environment (`julia`), then the package manager (`]`), then activate (`activate .`), then return to the Julia interpreter (`[backspace]`), then you are good to go with the latest version of the package! Import by typing `using Clusters` into the Julia interpreter.

-----------------------------

## Computing numerical frequencies

`compute_frequencies_henon_ae(φ,dφ/dr,d²φ/dr²,a,e)` return frequencies Ω₁ and Ω₂ for a given a,e (semi-major axis and eccentricity) + potential model using Henon anomaly. Will fall back to epicycle approximation when too close to circular.


## Basis element support

`potential_function(r,l,n,bc)` return the Clutton-Brock potential basis element at (l,n)  
`density_function(r,l,n,bc)` return the Clutton-Brock density basis element at (l,n)  

## Specific models

#### Isochrone
`isochrone_psi(r,M=1,bc=1)` return the isochrone potential  
`isochrone_dpsi_dr(r,M=1,bc=1)` return the isochrone potential first derivative  
`isochrone_ddpsi_ddr(r,M=1,bc=1)` return the isochrone potential second derivative  

#### Plummer
`plummer_psi(r,M,bc=1)` return the plummer potential  


## Analysing _n_-body runs
#### reader abilities

`return_particles(filename)` reads in a delimited file to x,y,z,vx,vy,vz  
`return_density_centre(filename)` reads in a centre file

#### centring support

`find_rbary(x,y,z)` computes the barycentre and returns the radius for a given x,y,z series  
`find_rbary(x,y,z,xc,yc,zc)` offsets series by xc,yc,zc and returns the radius for a given x,y,z series


#### power spectra computation

`crosscorrelation(x,y)` computes the cross correlation of two density centre series à la Heggie.  
`perturbation_density(x,y,z,rmax,m=1)` computes eq. 3 of Heggie (2020)  
`discrete_fourier_k(x,k)` computes eq. 4 of Heggie (2020)

#### escaper analysis

`plummer_potential(r,a=1)` return the simplest possible Plummer potential  
`escaper_energy(x,y,z,vx,vy,vz)` return the naive total energy (½v² + φ)  
`spherical_coordinates(x,y,z,vx,vy,vz)` return the spherical positions and velocities


-----------------------------

## Author

Mike Petersen -  @michael-petersen - petersen@iap.fr
