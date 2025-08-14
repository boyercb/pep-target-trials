
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

colSds <- function(x, na.rm=TRUE) {
  if (na.rm) {
    n <- colSums(!is.na(x)) # thanks @flodel
  } else {
    n <- nrow(x)
  }
  colVar <- colMeans(x*x, na.rm=na.rm) - (colMeans(x, na.rm=na.rm))^2
  return(sqrt(colVar * n/(n-1)))
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

bootfun <- function(d, i, hetero, compare_regimes, msg = FALSE) {
  db <- d[i, ]
  db$id <- 1:nrow(db)
  
  res <- runsim(db, NULL, hetero, compare_regimes, se = FALSE)
  
  if (msg) {
    print(str(res))
  }
  
  return(res)
}

runsim <- function(d, p, hetero = FALSE, compare_regimes = FALSE, se = FALSE, R = 1000) {
  if (!is.null(p)) {
    p()
  }

  if (compare_regimes) {
    # 7-day fixed enrollment period
    d_fixed <- 
      d |>
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
      ungroup() |>
      filter(trial <= 6)
    
    
    # day zero trial with following regimes
    #   0: never vaccinate
    #   1: vaccine on day 0
    #   2: vaccine on day 1
    #   3: vaccine on day 2
    #   4: vaccine on day 3
    #   5: vaccine on day 4
    #   6: vaccine on day 5
    #   7: vaccine on day 6
    d_dayzero <- d |>
      uncount(8, .id = "regime") |>
      group_by(id) |>
      mutate(
        Y = Y,
        censor = case_when(
          regime == 1 ~ pmin(T, I(X == 0) * 1 + I(X > 0) * (X + 1)),
          regime == 2 ~ pmin(T, I(X == 0) * T + I(X != 0) * 1),
          regime == 3 ~ pmin(T, I(X == 1) * T + I(X < 1) * (X + 1) + I(X > 1) * 1),
          regime == 4 ~ pmin(T, I(X == 2) * T + I(X < 2) * (X + 1) + I(X > 2) * 2),
          regime == 5 ~ pmin(T, I(X == 3) * T + I(X < 3) * (X + 1) + I(X > 3) * 3),
          regime == 6 ~ pmin(T, I(X == 4) * T + I(X < 4) * (X + 1) + I(X > 4) * 4),
          regime == 7 ~ pmin(T, I(X == 5) * T + I(X < 5) * (X + 1) + I(X > 5) * 5),
          regime == 8 ~ pmin(T, I(X == 6) * T + I(X < 6) * (X + 1) + I(X > 6) * 6)
        ),
        fup = T,
        pid = paste(id, regime, sep = "_")
      ) |>
      ungroup() |> 
      uncount(censor, .id = "day") |>
      group_by(pid) |>
      mutate(
        day = day - 1,
        Y = if_else(row_number() == n() & (day + 1) == T,
                    Y,
                    0L),
        censor = if_else(row_number() == n() &
                           Y == 0 &
                           (regime >= 2 & day != 20) | 
                           #(regime == 1 & day == I(X == 0) * (0) + I(X > 0) * (X - 1) & X != 21),
                           (regime == 1 & day == X & X != 21),
                         1L,
                         0L),
        A = as.numeric(regime >= 2),
        Y = if_else(censor == 1 & ((regime > 2 & X == 0) | regime == 2 | (regime == 1 & day == 0)), NA, Y)
      ) |>
      ungroup() 
    
    # 7-day grace period trial with following regimes
    #   0: never vaccinate
    #   1: vaccine within 7 days
    d_grace <- d |>
      uncount(2, .id = "regime") |>
      group_by(id) |>
      mutate(
        Y = Y,
        censor = case_when(
          regime == 1 ~ pmin(T, X + 1),
          regime == 2 ~ pmin(T, I(X <= 6) * T + I(X > 6) * 7)
        ),
        fup = T,
        pid = paste(id, regime, sep = "_")
      ) |>
      ungroup() |> 
      uncount(censor, .id = "day") |>
      group_by(pid) |>
      mutate(
        day = day - 1,
        Y = if_else(row_number() == n() & (day + 1) == T,
                    Y,
                    0L),
        censor = if_else(row_number() == n() &
                           Y == 0 &
                           (regime == 2 & day != 20) | 
                           (regime == 1 & day == X & X != 21),
                         1L,
                         0L),
        censor = if_else(regime == 2 & row_number() != n(), NA, censor)
      ) |>
      ungroup() 
    
    if (max(d_fixed$censor) > 0) {
      w_fixed <- ipcw(
        numerator = NULL,
        denominator = censor ~ bs(day + trial, 3), 
        data = d_fixed,
        treatment = "A",
        id = "pid",
        time = "day"
      )
      
    } else {
      w_fixed$weights <- rep(1, 1:nrow(d_fixed))
    }
    
    if (max(d_dayzero$censor) > 0) {
      w_dayzero <- ipcw(
        numerator = NULL,
        denominator = list(
          "1" = censor ~ bs(day, 3), 
          "2" = censor ~ 1, 
          "3" = censor ~ -1 + I(day == 0),
          "4" = censor ~ -1 + I(day == 0) + I(day == 1), 
          "5" = censor ~ -1 + I(day == 0) + I(day == 1) + I(day == 2), 
          "6" = censor ~ -1 + I(day == 0) + I(day == 1) + I(day == 2) + I(day == 3), 
          "7" = censor ~ -1 + I(day == 0) + I(day == 1) + I(day == 2) + I(day == 3) + I(day == 4), 
          "8" = censor ~ -1 + I(day == 0) + I(day == 1) + I(day == 2) + I(day == 3) + I(day == 4) + I(day == 5)
        ),
        data = d_dayzero,
        treatment = "regime",
        regimes = 1:8,
        id = "pid",
        time = "day",
        stop = as.list(
          c(21,0:6)
        ),
        dayzero = TRUE
      )
      cutoff <- quantile(w_dayzero$weights, 0.99)
      w_dayzero$weights[w_dayzero$weights > cutoff] <- cutoff
    } else {
      w_dayzero$weights <- rep(1, 1:nrow(d_dayzero))
    }

    if (max(d_grace$censor, na.rm = TRUE) > 0) {
      w_grace <- ipcw(
        numerator = NULL,
        denominator = list(
          "1" = censor ~ bs(day, 3), 
          "2" = censor ~ -1 + I(day == 6)
        ),
        data = d_grace,
        treatment = "regime",
        regimes = 1:2,
        id = "pid",
        time = "day", 
        grace = TRUE
      )
      # print(summary(w_grace$denominator_model[[2]]))
    } else {
      w_grace$weights <- rep(1, 1:nrow(d_grace))
    }
  
    suppressWarnings({
      fits <- list(
        target_trial_fixed =
          glm(Y ~ A + A:bs(trial, 4) + bs(day + trial, 3),
              data = d_fixed,
              family = binomial(link = "logit"),
              weights = w_fixed$weights),
        target_trial_dayzero =
          glm(Y ~ A + A:bs(regime - 1, 4) + bs(day, 4),
              data = d_dayzero,
              family = binomial(link = "logit"),
              weights = w_dayzero$weights
              ),
        target_trial_grace =
          glm(Y ~ I(regime == 2) * bs(day, 4),
              data = d_grace,
              family = binomial(link = "logit"),
              weights = w_grace$weights
              )
      )
    })

    ret <- list(
      calc_log_rr(
        fit = fits[[1]], 
        df = d_fixed, 
        pooled = TRUE,
        natural = TRUE
      ),
      calc_log_rr(
        fit = fits[[1]], 
        df = d_fixed, 
        pooled = FALSE
      ),
      sapply(
        2:8,
        function(x)
          calc_log_rr_regime(
            fit = fits[[2]],
            df = d,
            regimes = c(1, x)
          )
      ),
      # calc_hr_regime(fits[[2]]),
      calc_log_rr_regime(
        fit = fits[[3]], 
        df = d,
        regimes = c(1, 2)
      )
    )
    
    ret <- unlist(lapply(ret, function(x) x))
  
  } else {
    # 21-day fixed enrollment period design
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
        )
      })
      
      ret <- unlist(map(fits, 
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
                   
                   lp1 - lp0
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
      
      ret <- pmap_dbl(
        list(
          f = fits,
          idx = c(2, 2, 1, 2, 2, 2, 1), 
          c = c(rep(TRUE, 6), FALSE)
        ),
        function(f, idx, c) {
          if (c) {
            coef(f)[idx]
          } else {
            f
          }
        }
      )
    }
  }
  
  if (se) {
    
    bt <-
      boot(
        data = d,
        statistic = bootfun,
        R = R,
        parallel = "multicore",
        ncpus = 1, #detectCores() - 1,
        hetero = hetero,
        compare_regimes = compare_regimes,
        msg = FALSE
      )
    
    #print(runsim(d[sample(.N, .N, replace = TRUE)][, id := 1:.N], NULL, hetero, compare_regimes, se = FALSE))
    # bt <- mclapply(
    #   1:R,
    #   function(iter, d, hetero, compare_regimes) {
    #     db <- d[sample(.N, .N, replace = TRUE)][, id := 1:.N]
    #     res <- runsim(db, NULL, hetero, compare_regimes, se = FALSE)
    #     return(res)
    #   },
    #   d = d,
    #   hetero = hetero,
    #   compare_regimes = compare_regimes,
    #   mc.cores = detectCores() - 1
    # )

    ret_se <- colSds(bt$t)
    list("est" = ret, "se" = ret_se)
  } else {
    list("est" = ret)
  }
}


# standardized survival curves --------------------------------------------

calc_log_rr <- function(fit, df, pooled = TRUE, natural = FALSE) {
  newdf <-
    df |>
    filter(day == 0) |>
    mutate(remaining = 21 - trial) |>
    uncount(remaining, .id = "time") |>
    mutate(day = time - 1)
  
  if (natural) {
    newdf <- filter(newdf, A == 1)
  } 
    
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
    filter(day == 20 - trial)

  if (pooled) {
    log(mean(newdf$p1 / newdf$p0))
  } else {
    newdf |>
      group_by(trial) |>
      summarise(
        log_rr = log(mean(p1/p0))
      ) |>
      ungroup() |>
      pull(log_rr)
  }
}

calc_log_rr_regime <- function(fit, df, regimes, print = FALSE) {
  newdf <-
    df |>
    mutate(remaining = 21) |>
    uncount(remaining, .id = "time") |>
    mutate(day = time - 1)
  
  newdf$surv1 <- 1 - predict(
    fit,
    newdata = mutate(newdf, regime = regimes[2], A = 1),
    type = "response"
  ) 
  
  newdf$surv0 <- 1 - predict(
    fit,
    newdata = mutate(newdf, regime = regimes[1], A = 0),
    type = "response"
  ) 
  
  newdf <- 
    newdf |>
    group_by(id) |>
    mutate(
      p1 = 1 - cumprod(surv1),
      p0 = 1 - cumprod(surv0)
    ) |>
    filter(day == 20) |>
    ungroup()
  
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

calc_true_rr_X <- function(VE_func, X, conditional = FALSE, meanlog = 2.1, sdlog = 0.59) {
  if (conditional) {
    K <- (X + 1):21
  } else {
    K <- 1:21
  }
  
  surv0 <- vector(length = length(K))
  surv1 <- vector(length = length(K))
  
  for (i in seq_along(K)) {
    k <- K[i]
    alpha0 <- logit(bhaz(k, meanlog, sdlog))
    surv0[i] <- 1 - plogis(alpha0)
    surv1[i] <- 1 - plogis(alpha0 + log(1 - VE_func(X)) * I(X < k & X < 21))
  }
  
  df <- tibble(
    time = K,
    risk1 = (1 - cumprod(surv1)),
    risk0 = (1 - cumprod(surv0))
  )
  
  1 - rev(df$risk1)[1] / rev(df$risk0)[1]
}


calc_hr_regime <- function(fit) {
  lp1 <- predict(
    fit,
    newdata = tibble(
      regime = 2:8,
      day = 20,
      A = 1
    ),
    type = "term"
  ) |> rowSums()
  
  lp0 <- predict(
    fit,
    newdata = tibble(
      regime = 2:8,
      day = 20,
      A = 0
    ),
    type = "term"
  ) |> rowSums()

  lp1 - lp0
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
           time = "time",
           stop = NULL, 
           grace = FALSE,
           dayzero = FALSE) {
    # identify rows for always and never take regimes
    rows <- lapply(regimes, function(x) data[[treatment]] == x)
    names(rows) <- as.character(regimes)
    
    if (is.null(stop)) {
      stop <- rep(Inf, length(rows))
    }
    
    # extract censoring indicator
    if (is.list(denominator)) {
      censor <- all.vars(denominator[[1]])[1]
    } else {
      censor <- all.vars(denominator)[1]
    }
    
    # fit denominator model(s)
    denominator_model <- lapply(1:length(rows),
                                function(x) {
                                  if (sum(data[[censor]][rows[[x]]], na.rm = TRUE) > 0) {
                                    glm(
                                      formula = if (is.list(denominator)) {
                                        denominator[[x]]
                                      } else {
                                        denominator
                                      },
                                      family = binomial(link = "logit"),
                                      data = data[rows[[x]],]
                                    )
                                  } else {
                                    NULL
                                  }
                                })
    
    # get denominator predictions
    p_den <- lapply(1:length(rows),
                    function(x) {
                      if (grace & x == 2) {
                        ifelse(data[[time]] == 7 & data[["X"]] == 6, 
                               1 / (0.5 * (1 - 0.1124712) * (1 - ppois(5, 5)) * dpois(6, 5)),
                               1#1 - predict(denominator_model[[x]], newdata = data, type = "response")
                        )
                      } else {
                        if (!is.null(denominator_model[[x]])) {
                          # ifelse(data[[time]] <= stop[[x]], 
                          #        1,
                          #        1 - predict(denominator_model[[x]], newdata = data, type = "response")
                          # )
                          if (dayzero) {
                            ifelse(data[[time]] <= stop[[x]] & data[[time]] > 0 & x != 2, 
                                   1 - predict(denominator_model[[x]], newdata = mutate(data, day = ifelse(day == 0, 0, day - 1)), type = "response"),
                                   ifelse(data[[time]] <= stop[[x]] & x == 2,
                                          1 - predict(denominator_model[[x]], newdata = data, type = "response"),
                                          1
                                   )
                              )
                          } else {
                            1 - predict(denominator_model[[x]], newdata = data, type = "response")
                          }
                        } else {
                          1
                        }
                      }
                    })
    names(p_den) <- as.character(regimes)
    
    # fit numerator model
    if (!is.null(numerator)) {
      stopifnot(all.vars(numerator)[1] == all.vars(denominator)[1])
      
      numerator_model <- lapply(rows,
                                function(x) {
                                  if (sum(data[[censor]][x], na.rm = TRUE) > 0) {
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
      p_num <- lapply(1:length(rows),
                      function(x) {
                        if (grace & x == 2) {
                          ifelse(data[[time]] != 7, 
                                 1,
                                 1 - predict(numerator_model[[x]], newdata = data, type = "response")
                          )
                        } else {
                          if (!is.null(numerator_model[[x]])) {
                            # ifelse(data[[time]] <= stop[[x]], 
                            #        1,
                            #        1 - predict(numerator_model[[x]], newdata = data, type = "response")
                            # )
                            1 - predict(denominator_model[[x]], newdata = data, type = "response")
                            
                          } else {
                            1
                          }
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

