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
# This script: 

#####################
# Objective function
####################
# multidimension:
function Ln(theta, Y, X, M::Int64)

    Ln_theta = - sum(abs(Y- max(0, X*theta[:])))

  return Ln_theta

end


####################
# Prior distribution
####################
# define a uniform prior
function prior(theta, theta_upper_limit, theta_lower_limit, type_prior)

  if ((theta .<= theta_upper_limit) & (theta .>= theta_lower_limit)) == [true; true; true; true]

    if type_prior == "uniform"
      # 1-D case
      # Here it's for finding beta1
      # Adjust indexes if necessary ! 
      if size(theta) == () 
       pdf_value = pdf(Uniform(theta_lower_limit[2], theta_upper_limit[2]), theta[2])
      
      # Mutlidimensional case
      elseif size(theta)[1] > 1
        pdf_value = 1 #initialization

        for i=1:4
          pdf_value = pdf_value*pdf(Uniform(theta_lower_limit[i], theta_upper_limit[i]), theta[i])
        end
      else
        error("Check the dimension of theta")
      end

      return pdf_value

    elseif type_prior == "normal"
      #Select a large variance:
        sigma = 10.0

      if size(theta) == () #we are in the scalar case
        return pdf(Normal(0, sigma))
      elseif size(theta)[1] > 1
        return pdf(MvNormal(zeros(size(theta)[1]), sigma*eye(size(theta)[1])), theta)
      else
        error("Prior: check the dimension in theta")
      end
    end
  else  #if outside of the parameter set
      return 0
  end
end

