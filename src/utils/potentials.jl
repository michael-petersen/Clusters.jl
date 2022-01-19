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

function extremise_function(func::Function,neps::Int64=32)
    #=
    Find the single extremum of a function between 0 and 1

    Accuracy will be set by deps, 1/2^neps:
      can never do better than evaluating the function on a grid this fine.

    Requires exactly 2*neps evaluations of the function

    =#
    if neps > 53
        # can't be for double precision
        neps = 53
        print("extremise_function- reached maximum precision.\n")
    end
    
    deps = 1/2^(neps)

    # check the endpoint derivatives
    left_end_derivative = func(deps)-func(0.)
    right_end_derivative = func(1.) - func(1-deps)

    if left_end_derivative*right_end_derivative > 0
        print("Monotonic")
        return -1
    end

    left_derivative  = left_end_derivative
    right_derivative = right_end_derivative
    leftmin = 0.
    rightmax = 1.

    for iter=1:neps
        # assign a new midpoint
        midpoint = leftmin + 1/(2^iter)
        midfunc = func(midpoint)

        # compute finite difference derivatives on either side of the midpoint
        mid_left_derivative  = midfunc - func(midpoint-deps)
        mid_right_derivative = func(midpoint+deps) - midfunc

        # consider a special case of equal derivatives to block extra evaluations
        #print(mid_left_derivative,'=',mid_right_derivative)
        if mid_left_derivative == -mid_right_derivative
            # we are directly on the extremum!
            #print("Directly centred at ",midpoint)
            return midpoint
        end

        # check if on left side of the interval
        if mid_left_derivative*left_derivative < 0
            # set the right side (maximum) to be the midpoint
            rightmax = midpoint
            # set the right derivative to be the old midpoint left derivative
            right_derivative = mid_left_derivative
        end
        #else
        if mid_right_derivative*right_derivative < 0 # if not on left, must be on the right side
            # set the left side (minimum) to be the midpoint
            leftmin = midpoint

            left_derivative = mid_right_derivative
        end

        # watch the convergence fly by...
        #print(iter,' ',midpoint,' ',leftmin,' ',rightmax,'\n')

        # can add some fancy check here for oscillatory convergence, if we want

        # otherwise, return the value at the maximum number of iterations
        # (with an advancement to the next midpoint)
        if iter == neps
            return midpoint + 1/(2^(neps+1))
        end
    end
end
