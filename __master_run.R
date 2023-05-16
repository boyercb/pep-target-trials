# Target trials for post-exposure vaccination 
# Author: Christopher Boyer

N_SIMS <- 1000

# simulations -------------------------------------------------------------

# demonstrations of the advantages of adopting a target trial approach for the
# estimation of efficacy of postexposure vaccines.

source("1_code/0_packages.R")

source("1_code/1_functions.R")

source("1_code/2_sims.R")

source("1_code/3_plots.R")


# application (NYC Mpox) --------------------------------------------------

# NOTE: this code was shipped to NYC DOH so that they could run the analysis
# onsite. For debugging, I created a dummy dataset based on marginal
# distributions of important variables.

#source("1_code/mpox.R")