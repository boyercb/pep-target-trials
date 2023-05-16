
# helper functions --------------------------------------------------------

# logit and expit functions
expit <- function(x) exp(x) / (1 + exp(x))
logit <- function(x) log(x / (1 - x))

# truncated lognormal distribution functions
rlnormt <- function(n, range, meanlog, sdlog) {
  a <- plnorm(range[1], meanlog, sdlog)
  b <- plnorm(range[2], meanlog, sdlog)
  
  u <- runif(n, a, b)
  
  floor(qlnorm(u, meanlog, sdlog))
}

dlnormt <- function(q, range, meanlog, sdlog) {
  a <- plnorm(range[1], meanlog, sdlog)
  b <- plnorm(range[2], meanlog, sdlog)
  
  dlnorm(q, meanlog, sdlog) / (b - a)
}

plnormt <- function(q, range, meanlog, sdlog) {
  a <- plnorm(range[1], meanlog, sdlog)
  b <- plnorm(range[2], meanlog, sdlog)
  
  (plnorm(q, meanlog, sdlog) - a) / (b - a)
}


# data generation function ------------------------------------------------

# baseline hazard 
bhaz<- function(k, meanlog, sdlog) {
  f <- plnorm(k, meanlog, sdlog) - plnorm(k - 1, meanlog, sdlog) 
  S <- 1 - plnorm(k - 1, meanlog, sdlog) 
  0.25 * f/S
}

# data generation function
gendata <-
  function(N = 1000,                         # sample size
           sims = 5000,                      # number of simulations
           mu = 2.1,                         # shape parameter
           sigma = 0.59,                     # scale parameter
           X_dist = function(n) rpois(n, 5), # distribution of vaccine administration
           Y_dist = function(n, k, Z, X) {   # outcome distribution
             alpha0 <- logit(bhaz(k, 2.1, 0.59))
             VE <- 0.4
             rbinom(n, 1, plogis(alpha0 + log(1 - VE) * Z * I(X < k & X < 21)))
           }
           
  ) {
    
    n <- N * sims
    
    Z <- rbinom(n, 1, 0.5)
    W <- rbinom(n, 1, 0.25)
    X <- X_dist(n)
    
    Y <- matrix(nrow = n, ncol = 21)
    for (k in 1:21) {
      Y[, k] <- Y_dist(n, k, Z, X)
      
      if (k > 1) {
        Y[, k] <- replace(Y[, k], Y[, k - 1] == 1, 1L)
      }
    }
    
    T <- 22 - rowSums(Y)
    T <- replace(T, T > 21, 21)
    Y <- Y[, 21]

    dt <- data.table(
      sim = rep(1:sims, each = N),
      id = rep(1:N, sims),
      Z = Z,
      W = W,
      X = X * Z + (1 - Z) * 21,
      T = T * Y + (1 - Y) * 21,
      Y = Y
    ) 
    
    dt$A_leave <- 
      as.numeric(dt$X < dt$T | (dt$X >= dt$T & dt$W == 1 & dt$Z == 1))
    
    dt$A_move <- 
      as.numeric(dt$X < dt$T)

    setkey(dt, sim)
    
    return(dt)
  }


# simulation function -----------------------------------------------------

runsim <- function(d, p, hetero = FALSE) {
  p()
  
  d_long <- d |>
    uncount(pmin(X + 1, T), .id = "trial") |>
    group_by(id) |>
    mutate(
      trial = trial - 1,
      A = as.numeric(I(X == trial & T != trial & X != 21)),
      Y = Y,
      censor = pmin(T - trial, (1 - A) * (X - trial) + A * 21),
      fup = T - trial,
      pid = paste(trial, id, sep = "_")
    ) |>
    ungroup() |> 
    uncount(censor, .id = "day") |>
    group_by(pid) |>
    mutate(
      day = day - 1,
      Y = if_else(row_number() == n() & (day + trial + 1) == T,
                  Y,
                  0L),
      censor = if_else(row_number() == n() &
                         Y == 0 & A == 0 &
                         (day + trial + 1) == X & X != 21,
                       1L,
                       0L)
    ) |>
    ungroup()
  
  if (max(d_long$censor > 0)) {
    w <- ipcw(
      numerator = NULL,
      denominator = censor ~ bs(day + trial, 3), 
      data = d_long,
      treatment = "A",
      id = "pid",
      time = "day"
    )
  } else {
    w$weights <- 1:nrow(d_long)
  }
  
  if (hetero) {
    suppressWarnings({
      fits <- list(
        leave =
          glm(Y ~ A_leave + A_leave:bs(X, 4),
              data = d, 
              family = poisson(link = "log"), 
              offset = log(T)),
        move =
          glm(Y ~ A_move + A_move:bs(X, 4),
              data = d, 
              family = poisson(link = "log"), 
              offset = log(T)),
        cox = 
          coxph(Surv(start, stop, Y) ~ A + A:bs(X, 4), 
                data = d |>
                  uncount(A_move * I(X > 0) + 1, .id = "split") |>
                  mutate(start = 0, stop = as.numeric(T)) |>
                  mutate(
                    start = if_else(split == 2 & A_move == 1 & X > 0, X, start),
                    stop = if_else(split == 1 & A_move == 1 & X > 0, X, stop),
                    A = replace(A_move, split == 1 & A_move == 1 & X > 0, 0),
                    Y = replace(Y, split == 1 & A_move == 1 & X > 0, 0)
                  )
          ),
        target_trial =
          glm(Y ~ A + A:bs(trial, 4) + bs(day + trial, 3),
              data = d_long,
              family = binomial(link = "logit"),
              weights = w$weights)
        # leave_rr =
        #   glm(Y ~ A_leave + A_leave:bs(X, 4),
        #       data = d,
        #       family = poisson(link = "log")),
        # move_rr =
        #   glm(Y ~ A_move + A_move:bs(X, 4),
        #       data = d,
        #       family = poisson(link = "log")),
        # target_trial_rr =
        #   glm(Y ~ A + A:bs(trial, 4) + bs(day + trial, 3),
        #       data = d_long,
        #       family = binomial(link = "logit"),
        #       weights = w$weights) |>
        #   calc_log_rr()
      )
    })
    
    unlist(map(fits, 
        function(f, idx) {
          lp1 <- predict(
            f,
            newdata = tibble(
              X = 0:7,
              trial = 0:7,
              day = 0,
              Z = 1,
              A = 1,
              A_move = 1,
              A_leave = 1
            ),
            type = "term"
          ) |> rowSums() 
          
          lp0 <- predict(
            f,
            newdata = tibble(
              X = 0:7,
              trial = 0:7,
              day = 0,
              Z = 0,
              A = 0,
              A_move = 0,
              A_leave = 0
            ),
            type = "term"
          ) |> rowSums() 
          
          1 - exp(lp1 - lp0)
        }))
    
  } else {
    suppressWarnings({
      fits <- list(
        leave =
          glm(Y ~ A_leave,
              data = d, 
              family = poisson(link = "log"), 
              offset = log(T)),
        move =
          glm(Y ~ A_move,
              data = d, 
              family = poisson(link = "log"), 
              offset = log(T)),
        cox = 
          coxph(Surv(start, stop, Y) ~ A, 
                data = d |>
                  uncount(A_move * I(X > 0) + 1, .id = "split") |>
                  mutate(start = 0, stop = as.numeric(T)) |>
                  mutate(
                    start = if_else(split == 2 & A_move == 1 & X > 0, X, start),
                    stop = if_else(split == 1 & A_move == 1 & X > 0, X, stop),
                    A = replace(A_move, split == 1 & A_move == 1 & X > 0, 0),
                    Y = replace(Y, split == 1 & A_move == 1 & X > 0, 0)
                  )
          ),
        target_trial =
          glm(Y ~ A + bs(day + trial, 3),
              data = d_long,
              family = binomial(link = "logit"),
              weights = w$weights),
        leave_rr =
          glm(Y ~ A_leave,
              data = d,
              family = poisson(link = "log")),
        move_rr =
          glm(Y ~ A_move,
              data = d,
              family = poisson(link = "log")),
        target_trial_rr =
          calc_log_rr(
            glm(Y ~ A + bs(day + trial, 3),
              data = d_long,
              family = binomial(link = "logit"),
              weights = w$weights),
            df = d_long
          )
      )
    })
    
    pmap_dbl(
      list(
        f = fits,
        idx = c(2, 2, 1, 2, 2, 2, 1), 
        c = c(rep(TRUE, 6), FALSE)
      ),
      function(f, idx, c) {
        if (c) {
          1 - exp(coef(f)[idx])
        } else {
          1 - exp(f)
        }
      }
    )
  }
}


# standardized survival curves --------------------------------------------

calc_log_rr <- function(fit, df) {
  newdf <-
    df |>
    filter(day == 0) |>
    mutate(remaining = 21 - trial) |>
    uncount(remaining, .id = "time") |>
    mutate(day = time - 1)
  
  newdf$surv1 <- 1 - predict(
    fit,
    newdata = mutate(newdf, A = 1),
    type = "response"
  ) 
  
  newdf$surv0 <- 1 - predict(
    fit,
    newdata = mutate(newdf, A = 0),
    type = "response"
  ) 

  newdf <- 
    newdf |>
    group_by(pid, trial) |>
    mutate(
      p1 = 1 - cumprod(surv1),
      p0 = 1 - cumprod(surv0)
    ) |>
    ungroup() |>
    filter(day == 20)
  
  # print(head(select(newdf, trial, day, p1, p0)))
  log(mean(newdf$p1 / newdf$p0))
}

# calc true RR
calc_true_rr <- function(HR, meanlog = 2.1, sdlog = 0.59) {
  surv0 <- vector(length = length(1:21))
  surv1 <- vector(length = length(1:21))
  for (k in 1:21) {
    alpha0 <- logit(bhaz(k, meanlog, sdlog))
    surv0[k] <- 1 - plogis(alpha0)
    surv1[k] <- 1 - plogis(alpha0 + log(1 - HR))
  }
  
  df <- tibble(
    time = 1:21,
    risk1 = (1 - cumprod(surv1)),
    risk0 = (1 - cumprod(surv0))
  )
  
  1 - df$risk1[21] / df$risk0[21]
}


# ipcw --------------------------------------------------------------------

# function for creating inverse probability of censoring weights
ipcw <-
  function(numerator = NULL,
           denominator,
           data,
           treatment = "treatment",
           regimes = c(0, 1),
           id = "id",
           time = "time") {
    # identify rows for always and never take regimes
    rows <- lapply(regimes, function(x)
      data[[treatment]] == x)
    names(rows) <- as.character(regimes)
    
    # extract censoring indicator
    censor <- all.vars(denominator)[1]
    
    # fit denominator model(s)
    denominator_model <- lapply(rows,
                                function(x) {
                                  if (sum(data[[censor]][x]) > 0) {
                                    glm(
                                      formula = denominator,
                                      family = binomial(link = "logit"),
                                      data = data[x,]
                                    )
                                  } else {
                                    NULL
                                  }
                                })
    
    # get denominator predictions
    p_den <- lapply(denominator_model,
                    function(x) {
                      if (!is.null(x)) {
                        1 - predict(x, newdata = data, type = "response")
                      } else {
                        1
                      }
                    })
    names(p_den) <- as.character(regimes)
    
    # fit numerator model
    if (!is.null(numerator)) {
      stopifnot(all.vars(numerator)[1] == all.vars(denominator)[1])
      
      numerator_model <- lapply(rows,
                                function(x) {
                                  if (sum(data[[censor]][x]) > 0) {
                                    glm(
                                      formula = numerator,
                                      family = binomial(link = "logit"),
                                      data = data[x,]
                                    )
                                  } else {
                                    glm(
                                      formula = reformulate("1", censor),
                                      family = binomial(link = "logit"),
                                      data = data[x,]
                                    )
                                  }
                                })
      
      # get predictions
      p_num <- lapply(numerator_model,
                      function(x) {
                        if (!is.null(x)) {
                          1 - predict(x, newdata = data, type = "response")
                        } else {
                          1
                        }
                      })
      names(p_num) <- as.character(regimes)
      
      # add to dataset
      data <- dplyr::mutate(
        data,
        p_num = rowSums(dplyr::bind_cols(rows) * dplyr::bind_cols(p_num)),
        p_den = rowSums(dplyr::bind_cols(rows) * dplyr::bind_cols(p_den))
      )
      
      # calculate cumulative product
      data <- dplyr::group_by(data, .data[[id]])
      
      data <- dplyr::mutate(
        data,
        p_num = cumprod(p_num),
        # p_num = dplyr::if_else(
        #   .data[[censor]] == 1,
        #   (1 - p_num) * cumprod(dplyr::lag(p_num, 1, default = 1)),
        #   p_num * cumprod(dplyr::lag(p_num, 1, default = 1))
        # ),
        p_den = cumprod(p_den),
        # p_den = dplyr::if_else(
        #   .data[[censor]] == 1,
        #   (1 - p_den) * cumprod(dplyr::lag(p_den, 1, default = 1)),
        #   p_den * cumprod(dplyr::lag(p_den, 1, default = 1))
        # ),
        ipcw = p_num / p_den
      )
      
      data <- dplyr::ungroup(data)
      
    } else {
      numerator_model <- NULL
      
      # add to dataset
      data <- dplyr::mutate(
        data,
        p_num = 1,
        p_den = rowSums(dplyr::bind_cols(rows) * dplyr::bind_cols(p_den))
      )
      
      # calculate cumulative product
      data <- dplyr::group_by(data, .data[[id]])
      
      data <- dplyr::mutate(data,
                            p_den = cumprod(p_den),
                            # p_den = dplyr::if_else(
                            #   .data[[censor]] == 1,
                            #   (1 - p_den) * cumprod(dplyr::lag(p_den, 1, default = 1)),
                            #   p_den * cumprod(dplyr::lag(p_den, 1, default = 1))
                            # ),
                            ipcw = p_num / p_den)
      
      data <- dplyr::ungroup(data)
      
    }
    
    ret <- list(
      "numerator_model" = numerator_model,
      "denominator_model" = denominator_model,
      "p_num" = data$p_num,
      "p_den" = data$p_den,
      "weights" = data$ipcw
    )
  }

