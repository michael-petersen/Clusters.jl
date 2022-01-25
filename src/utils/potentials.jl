#=
some potential utilities

all could be improved!

=#

function dpotential_numerical(potential::Function,r::Float64,eps::Float64=0.0000001)
    # rough numerical derivative, but cheap
    # making eps any smaller incurs numerical noise
    return (potential(r+eps)-potential(r))/eps
    # tests suggest that this will have error of approximately 1.e-8
end
