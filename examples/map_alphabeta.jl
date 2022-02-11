
using Plots                  # Plotting capabilities
using FastGaussQuadrature    # To have access to the nodes and weights of the G-L (Gauss-Legendre) quadrature

using Clusters

const K_u = 50
const tabuGLquad = zeros(Float64,K_u) # Table of the nodes   for the G-L quadrature
const tabwGLquad = zeros(Float64,K_u) # Table of the weights for the G-L quadrature

function tabuwGLquad!()
    tabuGLquad[:], tabwGLquad[:] = gausslegendre(K_u) # Computing the nodes (u) and weights (w) of the G-L quadrature
end

tabuwGLquad!()

plot(tabuGLquad,tabwGLquad)

savefig("figures/gausstest.png")
# savefig(plot_ref, fn) # save the fig referenced by plot_ref as fn
