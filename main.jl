###################################################################################
# LTE Estimation of the Least Absolute Deviations for the Censored Regression model
###################################################################################
# Julien Pascal
#
# based on:
# - J. L. Powell (1984), "Least Absolute Deviations estimation for the censored regression model", Journal of Econometrics 25, 303-325
# - V. Chernozhukov and H. Hong (2003) "An MCMC approach to classical estimation" Journal of Econometrics 115, 293 - 346
#
#
# This script: Estimate the LAD for the censored regression model with a LTE approach
#
# To do: Finish the adaptive MCMC algorithm
#
using PyPlot
using Distributions
using LaTeXStrings
using Distributions
using GLM
using DataFrames

srand(123)

#######################
# path to main.jl file
#######################
# You may want to adapt this path to you setting:
path_to_main = "/home/julien/MCMC_Approach_Classical_Estimation"

include(string(path_to_main,"/src/censored_regression.jl"))
include(string(path_to_main,"/src/Metropolis.jl"))

######################################
# Decide what type of algorithm to use
######################################
algo_type = "metropolis"
#algo_type = "adaptive_metropolis"

################
# Generate data
###############
# true parameter
# produces approximately 40% of censoring rate:
beta = transpose([1 3 3 3]) #beta as a column vector
number_regressors = Int(size(beta)[1])
M = 1000 #number of people in the sample

#Generate a 3*M matrix of observations:
X = transpose(rand(MvNormal(zeros(3), eye(3)), M)) #transpose to have the usual matrix representation

#Add a first column of one for the constant term
X = hcat(ones(M),X)

#Generate error terms:
u = rand(Normal(0,1), M).*((X[:, 2]).^2)
#u = rand(Normal(0,1), M).*abs(X[:, 2])

#Generate Y star:
Y_star = X*beta[1:number_regressors] + u

#Generate truncation:
Y = max(0, Y_star)

#Censoring rate:
println("Censoring rate = $(sum(Y .== 0)/M)")

fig = figure("distribution of Y star",figsize=(10,10))
subplot(121)
title("distribution of Y star")
h1 = plt[:hist](Y_star, 30, normed="True", color="blue") # Histogram

subplot(122)
title("distribution of Y")
h2 = plt[:hist](Y, 30,  normed="True", color = "blue") # Histogram
savefig(string(path_to_main,"/output/censoring.png"))

#################
# MCMC parameter:
#################

# # 1-D Search
# sigma = 0.5 #"jump" parameter for the proposal matrix.

# Multidimensional
sigma = 0.15*eye(number_regressors) #"jump" parameter for the proposal matrix.

#A = rand(number_regressors,number_regressors)
#sigma = 0.001*(transpose(A)*A + eye(number_regressors))


N = 10000000 #number of times to iterate the markov chain
burnin = Int((1/10)*N)
##########################
# initial draw for MCMC is
##########################
# Use the OLS estimates
data = convert(DataFrame, hcat(Y,X))
# Update names:
rename!(data, {:x1=>:y})
rename!(data, {:x2=>:x0})
rename!(data, {:x3=>:x1})
rename!(data, {:x4=>:x2})
rename!(data, {:x5=>:x3})
OLS = glm(y ~ x1 + x2 + x3, data, Normal(), IdentityLink())

#OLS = glm(x1 ~ x2 + x3 + x4,  convert(DataFrame, hcat(Y_star,X)), Normal(), IdentityLink())

#intialize theta:
theta = zeros(number_regressors, N)
theta0 = coef(OLS)

# 1-D Search
# true value, except for beta1 which takes value from the OLS coeff:
#theta[1:4,1] = [beta[1]; theta0[2]; beta[3]; beta[4]] #initial value of the markov chain

# Multidimensional Search
# Initialize the first draw using OLS coefficients:
theta[1:4,1] = coef(OLS)

# Define the parameter space:
theta_upper_limit = theta[1:4,1] + 10
theta_lower_limit = theta[1:4,1] - 10

# To keep track of acceptance and rejection rates:
Accept = 0
Reject = 0

if algo_type == "metropolis"

  #############################
  # Simple Metropolis Algorithm
  #############################   
  for m in 2:N
    #sample proposed from a symmetric jumping distribution:

    # Multidimensional Search
    theta_trial = generate_proposal(theta[1:4, m-1], 100, sigma, "normal")

    # 1-D Search
    #theta0_trial = generate_proposal(theta[2, m-1], 100, sigma, "normal")
    

    #draw from a uniform on [0,1]:
    u = rand(Uniform(0,1))

    if u < (exp(Ln(theta_trial, Y, X, M))*prior(theta_trial, theta_upper_limit, theta_lower_limit, "normal")/(exp(Ln(theta[1:4, m-1], Y, X, M))*prior(theta[1:4, m-1], theta_upper_limit, theta_lower_limit, "normal")))
      #accept the new draw
      Accept = Accept + 1
      theta[1:4, m] = theta_trial
    else
      #keep the previous draw
      Reject = Reject + 1
      theta[1:4, m] = theta[1:4, m-1]
    end

  end
##########################################
elseif algo_type == "adaptive_metropolis"
  #########################################
  # Adaptive MCMC with vanishing adaptation
  #########################################
  #########################################
  # WORK IN PROGRESS : the algorith is not 
  # operational so far
  #########################################
  # based on "Adaptive Markov chain Monte Carlo sampling
  # and estimation in Mata", J Baker 

  # Initialization of extra parameters (follow the paper):
  # the initial draw (X0) was already set above
  mu_0 = theta0 + 1 #need to initialize in a way. 
  Sigma_0 = 1.0*eye(number_regressors) 
  #
  #A = rand(number_regressors,number_regressors)
  #Sigma_0 = 0.1*(transpose(A)*A + eye(number_regressors))

  d = number_regressors
  lambda_0 = (2.38^2)/d #recommended initial scaling parameter (according to the paper quoted above)

  alpha_star = 0.234 #optimal acceptance rate
  delta = 2/3

  # Initialize vectors in which I am storing the "tuning" parameters as the iteration moves forward:
  mu_t = zeros(number_regressors, N)*NaN
  lambda_t = zeros(N)*NaN
  gamma_t = zeros(N)*NaN
  Sigma_t = zeros(number_regressors,number_regressors,N)*NaN #mutlidimensional array 
  Var_covariance_proposal_t = zeros(number_regressors,number_regressors,N)*NaN #mutlidimensional array 

  # Store the initial values:
  mu_t[:, 1] = mu_0
  lambda_t[1] = lambda_0
  gamma_t[1] = 1/((1+1)^delta)
  Sigma_t[:,:,1] = Sigma_0
  Var_covariance_proposal_t[:,:,1] = lambda_0*Sigma_0

  tol = 0.01

  for m in 2:N

    # Draw a candidate: 
    # Multidimensional Search:
    theta_trial = generate_proposal(theta[1:4, m-1], 100, Var_covariance_proposal_t[:,:,(m-1)], "normal")

    # Calculate the "pseudo-likelihood" ratio:
    pseudo_likelihood_ratio = (exp(Ln(theta_trial, Y, X, M))*prior(theta_trial, theta_upper_limit, theta_lower_limit, "normal")/(exp(Ln(theta[1:4, m-1], Y, X, M))*prior(theta[1:4, m-1], theta_upper_limit, theta_lower_limit, "normal")))
    
    # will be useful for later steps:
    # have to take care of NaN that can occur (division by zero):
    if isnan(pseudo_likelihood_ratio) == true
      alpha_Y_X = 0
    elseif pseudo_likelihood_ratio > 1 
      alpha_Y_X  = 1
    else 
      alpha_Y_X = pseudo_likelihood_ratio
    end
    

    #draw from a uniform on [0,1]:
    u = rand(Uniform(0,1))

    if u < alpha_Y_X 
      #accept the new draw
      Accept = Accept + 1
      theta[1:4, m] = theta_trial
    else
      #keep the previous draw
      Reject = Reject + 1
      theta[1:4, m] = theta[1:4, m-1]
    end

    # [Step 6 in the paper] Compute weighting parameter gamma_t:
    gamma_t[m-1] = 1/((1+(m-1))^delta)

    # when gamma_t is small enough, stop the adujstement phase:
    if gamma_t[m-1] > tol
      # [Step 7 in the paper] update  lambda_t_+1:
      lambda_t[m] = exp(gamma_t[m-1]*(alpha_Y_X - alpha_star))*lambda_t[m-1]

      # [Step 8 in the paper] update mu_t_+1:
      #precalculate an element:
      A = theta[1:4, m] - mu_t[:,(m-1)]
      mu_t[:,m] = mu_t[:,(m-1)] + gamma_t[m-1]*(A)

      # [Step 9 in the paper] update Sigma_t_+1:
      #precalculate an element:
      B = A*transpose(A) - Sigma_t[:,:,(m-1)]
      if isnan(gamma_t[m-1]*(B)) == false #when gamma_t goes to zero, there is a moment when it appears as "NaN"
        Sigma_t[:,:,m] = Sigma_t[:,:,(m-1)] + gamma_t[m-1]*(B)
      else
        Sigma_t[:,:,m] = Sigma_t[:,:,(m-1)]
      end
      # [Extra step not explicit in the paper] update the variance-covariance for the proposal distribution:
      Var_covariance_proposal_t[:,:,m]  = lambda_t[m]*Sigma_t[:,:,m]
    else # when the adjustement phase is over, keep the previous variance-covariance:
      Var_covariance_proposal_t[:,:,m]  = Var_covariance_proposal_t[:,:,(m-1)]
    end
    
    
  end
#################
else 
  error("Specify a correct algorithm type")
end


############
# Plot Draws
############
fig = figure("Draws of theta",figsize=(10,10))
subplot(221)
title("theta0")
#Use the function "squeeze" to go from Array{Float64,2} to Array{Float64,1}: necessary for plotting
plot(1:N, squeeze(theta[1,:], 1), color=:grey, alpha=0.5)
ylabel("theta")
xlabel("period")


subplot(222)
title("theta1")
plot(1:N, squeeze(theta[2,:],1), color=:grey, alpha=0.5)
ylabel("theta")
xlabel("period")

subplot(223)
title("theta2")
plot(1:N, squeeze(theta[3,:],1), color=:grey, alpha=0.5)
ylabel("theta")
xlabel("period")

subplot(224)
title("theta3")
plot(1:N, squeeze(theta[4,:],1), color=:grey, alpha=0.5)
ylabel("theta")
xlabel("period")
savefig(string(path_to_main,"/output/draws_censored_regression.png"))

#################
# Plot Histograms
#################
fig = figure("Histogram for theta",figsize=(10,10))
subplot(221)
title("theta0")
#Use the function "squeeze" to go from Array{Float64,2} to Array{Float64,1}: necessary for plotting
h0 = plt[:hist](squeeze(theta[1, burnin:N],1), 50, normed="True", color="blue") # Histogram

subplot(222)
title("theta1")
h1 = plt[:hist](squeeze(theta[2, burnin:N],1), 50, normed="True", color="blue") # Histogram


subplot(223)
title("theta2")
h2 = plt[:hist](squeeze(theta[3, burnin:N],1), 50, normed="True", color="blue") # Histogram

subplot(224)
h3 = plt[:hist](squeeze(theta[4, burnin:N],1), 50, normed="True", color="blue") # Histogram
title("theta3")
savefig(string(path_to_main,"/output/histograms_censored_regression.png"))

################################
# Use the sample to find moments
################################
result_beta0 = quantile(squeeze(theta[1, burnin:N],1),[0.10, 0.5, 0.90])
result_beta1 = quantile(squeeze(theta[2, burnin:N],1),[0.10, 0.5, 0.90])
result_beta2 = quantile(squeeze(theta[3, burnin:N],1),[0.10, 0.5, 0.90])
result_beta3 = quantile(squeeze(theta[4, burnin:N],1),[0.10, 0.5, 0.90])

results = DataFrame(variable = ["P10"; "Median"; "P90"], beta0 = result_beta0, beta1 = result_beta1, beta2 = result_beta2, beta3 = result_beta3)
writetable(string(path_to_main,"/output/output_table.csv"), results)


println("Rejection rate = $(Reject/N)")
println("Acception rate = $(Accept/N)")



