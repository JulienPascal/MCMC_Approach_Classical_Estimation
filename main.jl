#------------------------------------------------------------------------------
# LTE Estimation of the Least Absolute Deviations for the Censored Regression model
#------------------------------------------------------------------------------
# Julien Pascal
# https://julienpascal.github.io/
#
# based on:
# - J. L. Powell (1984), "Least Absolute Deviations estimation for the censored regression model", Journal of Econometrics 25, 303-325
# - V. Chernozhukov and H. Hong (2003) "An MCMC approach to classical estimation" Journal of Econometrics 115, 293 - 346
#
# This script: Estimate the LAD for the censored regression model with a LTE approach
#------------------------------------------------------------------------------
using Distributed
using RobustPmap
using Test
using Random
using Statistics
using Plots

using CSV
using GLM
using DataFrames

addprocs(3)
@everywhere using AffineInvariantMCMC
@everywhere using Distributions
@everywhere using LinearAlgebra
#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------
path_to_main = pwd()

#------------------------------------------------------------------------------
# Generate data
#------------------------------------------------------------------------------
# true parameters
# produces approximately 40% of censoring rate:
beta = transpose([1 2 3 4]) #beta as a column vector
number_regressors = size(beta,1)
M = 1000 #number of people in the sample

#Generate a 3*M matrix of observations:
d = MvNormal(zeros(3), I(3))
X = transpose(rand(d, M)) #transpose to have the usual matrix representation

#Add a first column of one for the constant term
X = hcat(ones(M),X)
#Generate error terms:
u = rand(Normal(0,1), M).*((X[:, 2]).^2)
#Generate Y star:
Y_star = X*beta .+ u
#Generate truncation:
Y = max.(0.0, Y_star)

#Censoring rate:
println("Censoring rate = $(100*sum(Y .== 0)/M) %")

h1 = histogram(Y_star) # Histogram
title!(h1, "Non-Truncated Distribution", legend=:none)
h2 = histogram(Y) # Histogram
title!(h2, "Truncated Distribution", legend=:none)
h = plot(h1, h2)
savefig(h, joinpath(path_to_main,"output/censoring.png"))


#-------------------------------------------------------------------------------
# MCMC Estimation
#-------------------------------------------------------------------------------
@everywhere begin
	numdims = 4
	numwalkers = 100
	thinning = 10
	numsamples_perwalker = 1000
	burnin = Int((1/10)*numsamples_perwalker)
	lb = 0 .* ones(numdims) #lower bound
	ub = 5 .* ones(numdims) #upper bound
	# Uniform prior
	d_prior = Product(Uniform.(lb, ub))
	# Normal
	# d_prior = MvNormal(zeros(numdims), 0.1 .* I(numdims))
end

# Log-likelihood
@everywhere function Ln(theta, Y, X)
	Ln_theta = - sum(abs.(Y .- max.(0.0, X*theta)))
	return Ln_theta
end

# Log quasi-posterior: Log(likelihood) + log(prior)
@everywhere function quasi_posterior(theta, Y, X, d_prior)
	 return Ln(theta, Y, X) + log(pdf(d_prior, theta))
end

# initial draw for MCMC is OLS
data = DataFrame(y = Y[:,1], x1 = X[:,2], x2 = X[:,3], x3 = X[:,4])
OLS = glm(@formula(y ~ x1 + x2 + x3), data, Normal(), IdentityLink())
theta0 = coef(OLS)

# Option 1: slightly perturb the initial draws for the walkers
x0 = ones(numdims, numwalkers).*theta0 .+ rand(numdims, numwalkers) .* 1.0
# Option 2: introduce a bias to see how fast the chains converge
#x0 = ones(numdims, numwalkers).*theta0 .+ ones(numdims, numwalkers) .* 1.0

# See https://github.com/madsjulia/AffineInvariantMCMC.jl
#=
chain, llhoodvals = AffineInvariantMCMC.sample(x -> Ln(x, Y, X), numwalkers, x0, burnin, 1)
chain, llhoodvals = AffineInvariantMCMC.sample(x -> Ln(x, Y, X), numwalkers, chain[:, :, end], numsamples_perwalker, thinning)
flatchain, flatllhoodvals = AffineInvariantMCMC.flattenmcmcarray(chain, llhoodvals)
=#

chain, llhoodvals = AffineInvariantMCMC.sample(x -> quasi_posterior(x, Y, X, d_prior), numwalkers, x0, burnin, 1)
chain, llhoodvals = AffineInvariantMCMC.sample(x -> quasi_posterior(x, Y, X, d_prior), numwalkers, chain[:, :, end], numsamples_perwalker, thinning)
flatchain, flatllhoodvals = AffineInvariantMCMC.flattenmcmcarray(chain, llhoodvals)

#-------------------------------------------------------------------------------
# Plot Draws
#-------------------------------------------------------------------------------
#Use the function "squeeze" to go from Array{Float64,2} to Array{Float64,1}: necessary for plotting
p1 = plot(flatchain[1,:], ylabel="theta0", xlabel="T", legend=:none)
p2 = plot(flatchain[2,:], ylabel="theta1", xlabel="T", legend=:none)
p3 = plot(flatchain[3,:], ylabel="theta2", xlabel="T", legend=:none)
p4 = plot(flatchain[4,:], ylabel="theta3", xlabel="T", legend=:none)
p5 = plot(p1, p2, p3, p4)
savefig(p5, joinpath(path_to_main,"output/chains.png"))
display(p5)

#-------------------------------------------------------------------------------
# Plot Histograms
#-------------------------------------------------------------------------------
hh1 = histogram(flatchain[1,burnin:end], title="theta0", legend=:none)
hh2 = histogram(flatchain[2,burnin:end], title="theta1", legend=:none)
hh3 = histogram(flatchain[3,burnin:end], title="theta2", legend=:none)
hh4 = histogram(flatchain[4,burnin:end], title="theta3", legend=:none)
hh5 = plot(hh1, hh2, hh3, hh4)
savefig(hh5, joinpath(path_to_main,"output/histograms.png"))
display(hh5)

#-------------------------------------------------------------------------------
# Use the sample to find quasi-posterior summary statistics
#-------------------------------------------------------------------------------
# append true value as well
result_beta0 = append!(quantile(flatchain[1,burnin:end],[0.10, 0.5, 0.90]), beta[1])
result_beta1 = append!(quantile(flatchain[2,burnin:end],[0.10, 0.5, 0.90]), beta[2])
result_beta2 = append!(quantile(flatchain[3,burnin:end],[0.10, 0.5, 0.90]), beta[3])
result_beta3 = append!(quantile(flatchain[4,burnin:end],[0.10, 0.5, 0.90]), beta[4])
results = DataFrame(variable = ["P10"; "Median"; "P90"; "True value"], beta0 = result_beta0,
						beta1 = result_beta1, beta2 = result_beta2, beta3 = result_beta3)

CSV.write(joinpath(path_to_main,"output/output_table.csv"), results)

#-------------------------------------------------------------------------------
# Example with GMM
#-------------------------------------------------------------------------------
# Let's assume now that we observe Y_star
@everywhere begin
	numdims = 4
	numwalkers = 100
	thinning = 10
	numsamples_perwalker = 2000
	burnin = Int((1/10)*numsamples_perwalker)
	lb = 0 .* ones(numdims) #lower bound
	ub = 5 .* ones(numdims) #upper bound
	# Uniform prior
	d_prior = Product(Uniform.(lb, ub))
end


# Log-likelihood
@everywhere function Ln_GMM(theta, Y, X)
	m = (1/sqrt(size(X,1))).*(transpose(X)*(Y .- X*theta))
	Ln_theta = - 0.5*transpose(m)*I(size(theta,1))*m
	return Ln_theta[1]
end

chain, llhoodvals = AffineInvariantMCMC.sample(x -> Ln_GMM(x, Y_star, X), numwalkers, x0, burnin, 1)
chain, llhoodvals = AffineInvariantMCMC.sample(x -> Ln_GMM(x, Y_star, X), numwalkers, chain[:, :, end], numsamples_perwalker, thinning)
flatchain, flatllhoodvals = AffineInvariantMCMC.flattenmcmcarray(chain, llhoodvals)

#-------------------------------------------------------------------------------
# Plot Draws
#-------------------------------------------------------------------------------
#Use the function "squeeze" to go from Array{Float64,2} to Array{Float64,1}: necessary for plotting
p1 = plot(flatchain[1,:], ylabel="theta0", xlabel="T", legend=:none)
p2 = plot(flatchain[2,:], ylabel="theta1", xlabel="T", legend=:none)
p3 = plot(flatchain[3,:], ylabel="theta2", xlabel="T", legend=:none)
p4 = plot(flatchain[4,:], ylabel="theta3", xlabel="T", legend=:none)
p5 = plot(p1, p2, p3, p4)
savefig(p5, joinpath(path_to_main,"output/chains_GMM.png"))
display(p5)


hh1 = histogram(flatchain[1,burnin:end], title="theta0", legend=:none)
hh2 = histogram(flatchain[2,burnin:end], title="theta1", legend=:none)
hh3 = histogram(flatchain[3,burnin:end], title="theta2", legend=:none)
hh4 = histogram(flatchain[4,burnin:end], title="theta3", legend=:none)
hh5 = plot(hh1, hh2, hh3, hh4)
savefig(hh5, joinpath(path_to_main,"output/histograms_GMM.png"))
display(hh5)

# Compare to traditional OLS
data = DataFrame(y = Y_star[:,1], x1 = X[:,2], x2 = X[:,3], x3 = X[:,4])
GMM = glm(@formula(y ~ x1 + x2 + x3), data, Normal(), IdentityLink())
stds_GMM = stderror(GMM)
coefs_GMM = coef(GMM)

result_beta0 = append!(quantile(flatchain[1,burnin:end],[0.05, 0.10, 0.5, 0.90, 0.95]), std(flatchain[1,burnin:end]), beta[1])
result_beta1 = append!(quantile(flatchain[2,burnin:end],[0.05, 0.10, 0.5, 0.90, 0.95]), std(flatchain[2,burnin:end]), beta[2])
result_beta2 = append!(quantile(flatchain[3,burnin:end],[0.05, 0.10, 0.5, 0.90, 0.95]), std(flatchain[3,burnin:end]), beta[3])
result_beta3 = append!(quantile(flatchain[4,burnin:end],[0.05, 0.10, 0.5, 0.90, 0.95]), std(flatchain[4,burnin:end]), beta[4])
results = DataFrame(variable = ["P5"; "P10"; "Median"; "P90"; "P95"; "std"; "True value"], beta0 = result_beta0,
						beta1 = result_beta1, beta2 = result_beta2, beta3 = result_beta3)

CSV.write(joinpath(path_to_main,"output/output_table_GMM.csv"), results)
