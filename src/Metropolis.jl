######################
# Metropolis Algorithm
######################
# Julien Pascal
#
# based on:
#  - Nichola Metropolis, Arianna W. Rosenbluth, Marshall N. Rosenbluth, Augusta H. Teller, and Edward Telle (1953) "Equation of State Calculations by Fast Computing Machines"
#  J. Chem. Phys. 21, 1087 (1953); doi: 10.1063/1.1699114
# - W.K. Hastings (1970), "Monte Carlo Sampling Methods Using Markov Chains and Their Applications", Biometrika, Vol. 57, No. 1. (Apr., 1970), pp. 97-109
#
# This script: define the function used in the algorithm

# function proportional to the beta distribution:
function prop_Beta(x, alpha, beta)
  if (x <= 1) & (x>=0)
    f_x = (x^(alpha-1))*(1-x)^(beta-1)
    return f_x
  else
    return 0
  end
end

# function the generate a "proposal" draw from the state space (based on the previous draw)
function generate_proposal(from, to, sigma, type_generation)
  if type_generation == "normal"
    #gaussian:
    if size(sigma) == () #we are in the scalar case
      return rand(Normal(from, sigma))
    elseif size(sigma)[1] > 1
      return rand(MvNormal(from, sigma))
    else
      error("Check the dimension of sigma")
    end
  #As in "An MCMC approach to classical estimation" Victor Chernozhukov and Han Hong (2003)
  elseif type_generation == "cauchy"
    if size(sigma) == () #we are in the scalar case
      return rand(Cauchy(from))
    elseif size(sigma)[1] > 1
      error("Multivariate Cauchy not implemented in Julia so far")
    else
      error("Check the dimension of sigma")
    end
  else
    error("Select a correct type for the proposal matrix")
  end

end
