
set.seed(8761275)

# scenario 1: VE = 0% -----------------------------------------------------

dgp0 <- gendata(
  sims = N_SIMS,
  Y_dist = function(n, k, Z, X) { 
    alpha0 <- logit(bhaz(k, 2.1, 0.59))
    VE <- 0
    rbinom(n, 1, plogis(alpha0 + log(1 - VE) * Z * I(X < k & X < 21)))
  }
) 

with_progress({
  p <- progressor(steps = N_SIMS)
  p("scenario 1: VE = 0%", class = "sticky")
  
  sim0 <- dgp0[, .(
    VE = unlist(runsim(.SD, p)),
    truth = 0,
    type = c("leave",
             "move",
             "cox",
             "target_trial",
             "leave_rr",
             "move_rr",
             "target_trial_rr")
  ), by = sim]
})

sim0_results <-
  sim0[, .(est = mean(VE),
           bias = mean(VE - mean(truth)),
           sd = sd(VE),
           rmse = sqrt(mean((VE - mean(truth))^2))
           ),
       by = list(type)]


# scenario 2: VE = 40% ----------------------------------------------------

dgp40 <- gendata(sims = N_SIMS) 

with_progress({
  p <- progressor(steps = N_SIMS)
  p("scenario 2: VE = 40%", class = "sticky")
  
  sim40 <- dgp40[, .(
    VE = unlist(runsim(.SD, p)),
    truth = c(rep(0.4, 4), rep(calc_true_rr(0.4), 3)),
    type = c("leave",
             "move",
             "cox",
             "target_trial",
             "leave_rr",
             "move_rr",
             "target_trial_rr")
  ), by = sim]
})

sim40_results <- 
  sim40[, .(est = mean(VE),
            bias = mean(VE - mean(truth)),
            sd = sd(VE),
            rmse = sqrt(mean((VE - mean(truth))^2))
            ),
        by = list(type)]


# scenario 3: VE(t) -------------------------------------------------------

VE_func <- function(x) {
  0.8 / (1 + exp(0.75 * (x - 4)))
}

dgpV <- gendata(
  sims = N_SIMS,
  Y_dist = function(n, k, Z, X) {   # outcome distribution
    alpha0 <- logit(bhaz(k, 2.1, 0.59))
    rbinom(
      n = n, 
      size = 1, 
      prob = plogis(alpha0 + log(1 - VE_func(X)) * Z * I(X < k & X < 21))
    )
  }
)  

with_progress({
  p <- progressor(steps = N_SIMS)
  p("scenario 3: VE(t)", class = "sticky")
  
  simV <- dgpV[, .(
    VE = unlist(runsim(.SD, p, hetero = TRUE)),
    X = rep(0:7, 4),
    truth = rep(VE_func(0:7), 4),
    type = rep(c("leave", "move",  "cox", "target_trial"), each = 8)
  ), by = sim]
})

simV_results <-
  simV[, .(est = mean(VE),
           bias = mean(VE - mean(truth)),
           sd = sd(VE),
           rmse = sqrt(mean((VE - mean(truth))^2))
           ),
       by = list(type, X)]


# scenario 4: overlap -----------------------------------------------------

meanlogs <- c(log(3), log(9), log(15))
  
dgp_overlap <- 
  map(meanlogs, function(mu) {
    dt <- gendata(
      sims = N_SIMS,
      Y_dist = function(n, k, Z, X) { 
        alpha0 <- logit(bhaz(k, mu, 0.59))
        VE <- 0
        rbinom(n, 1, plogis(alpha0 + log(1 - VE) * Z * I(X < k & X < 21)))
      }
    ) 
    dt$meanlog <- mu
    
    dt
  })
  
dgp_overlap <- rbindlist(dgp_overlap)

with_progress({ 
  p <- progressor(steps = N_SIMS * length(meanlogs))
  p("scenario 4: overlap", class = "sticky")
  
  sim_overlap <- dgp_overlap[, .(
    VE = unlist(runsim(.SD, p)),
    truth = 0,
    type = c("leave",
             "move",
             "cox",
             "target_trial",
             "leave_rr",
             "move_rr",
             "target_trial_rr")
  ), by = list(meanlog, sim)]
})

sim_overlap_results <-
  sim_overlap[, .(est = mean(VE),
                  bias = abs(mean(VE) - mean(truth)),
                  sd = sd(VE),
                  rmse = sqrt(mean((VE - mean(truth))^2))
                  ),
              by = list(meanlog, type)]

