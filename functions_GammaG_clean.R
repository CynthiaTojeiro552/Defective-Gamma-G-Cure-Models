################################################################################
# functions_GammaG_clean.R
#
# Core functions for the paper:
# Defective Gamma-G family for cure fraction models
#
# This file keeps only the functions needed for the defective Gompertz,
# defective Dagum, Gamma-Gompertz and Gamma-Dagum models used in the manuscript.
################################################################################

################################################################################
# 1. Survival and density functions of the baseline defective distributions
################################################################################

# Defective Dagum survival function
s_dd <- function(y, alpha, beta, p0, theta = NULL, log.p = FALSE) {
  if (any(alpha < 0)) stop("alpha must be positive")
  if (any(beta < 0)) stop("beta must be positive")
  if (any(y <= 0)) stop("y must be greater than 0")

  theta <- 1 - p0

  surv <- (beta + theta * y^(-alpha) - theta * beta) /
    (beta + theta * y^(-alpha))

  if (log.p) {
    surv <- log(surv)
  }

  return(surv)
}

# Defective Dagum density function
d_dd <- function(y, alpha, beta, p0, theta = NULL) {
  if (any(alpha < 0)) stop("alpha must be positive")
  if (any(beta < 0)) stop("beta must be positive")
  if (any(y <= 0)) stop("y must be greater than 0")

  theta <- 1 - p0

  density <- (alpha * beta * theta^2 * y^(-alpha - 1)) /
    (beta + theta * y^(-alpha))^2

  return(density)
}

# Defective Gompertz survival function
s_dg <- function(y, alpha, beta, p0 = NULL, theta = NULL, log.p = FALSE) {
  if (any(beta < 0)) stop("beta must be positive")
  if (any(y <= 0)) stop("y must be greater than 0")

  surv <- exp(-(beta / alpha) * (exp(alpha * y) - 1))

  if (log.p) {
    surv <- log(surv)
  }

  return(surv)
}

# Defective Gompertz density function
d_dg <- function(y, alpha, beta, p0 = NULL, theta = NULL) {
  if (any(beta < 0)) stop("beta must be positive")
  if (any(y <= 0)) stop("y must be greater than 0")

  density <- beta * exp(alpha * y - (beta / alpha) * (exp(alpha * y) - 1))

  return(density)
}


################################################################################
# 2. Gamma-G survival and density functions
################################################################################

# Gamma-G survival function
surv_gammag <- function(y, mu, surv_func) {
  func_result <- -log(surv_func(y))
  out <- 1 - pgamma(func_result, shape = mu, lower.tail = TRUE)
  return(out)
}

# Gamma-G density function
pdf_gammag <- function(y, mu, distribution, alpha, beta, p0 = NULL, theta = NULL) {

  pdf_result <- switch(
    distribution,
    "dagum" = d_dd(y, alpha, beta, p0, theta),
    "gompertz" = d_dg(y, alpha, beta, p0, theta),
    stop("Invalid distribution. Use 'dagum' or 'gompertz'.")
  )

  surv_result <- switch(
    distribution,
    "dagum" = s_dd(y, alpha, beta, p0, theta, log.p = FALSE),
    "gompertz" = s_dg(y, alpha, beta, p0, theta, log.p = FALSE),
    stop("Invalid distribution. Use 'dagum' or 'gompertz'.")
  )

  out <- (1 / gamma(mu)) * pdf_result * (-log(surv_result))^(mu - 1)

  return(out)
}


################################################################################
# 3. Log-likelihood functions
################################################################################

# Log-likelihood for Gamma-G models
log_like <- function(
    par,
    dados,
    delta,
    distribution,
    pred_alpha,
    pred_beta,
    pred_p0 = NULL,
    log = TRUE
) {

  surv_func <- switch(
    distribution,
    "dagum" = s_dd,
    "gompertz" = s_dg,
    stop("Invalid distribution. Use 'dagum' or 'gompertz'.")
  )

  delta <- get(delta, dados)
  t <- as.numeric(stats::model.extract(stats::model.frame(pred_alpha, dados), "response"))

  # Model matrices
  x_alpha <- stats::model.matrix(pred_alpha, dados)
  x_beta  <- stats::model.matrix(pred_beta, dados)

  J_alpha <- dim(x_alpha)[2]
  J_beta  <- dim(x_beta)[2]

  # Parameter vector:
  # par = (log(mu), alpha coefficients, beta coefficients, p0 coefficients if Dagum)
  mu <- exp(par[1])

  aa <- par[2:(J_alpha + 1)]
  be <- par[(J_alpha + 2):(J_alpha + J_beta + 1)]

  # Linear predictors
  alpha <- tcrossprod(x_alpha, t(aa))
  if (identical(distribution, "dagum")) {
    alpha <- exp(alpha)
  }

  beta <- exp(tcrossprod(x_beta, t(be)))

  if (!is.null(pred_p0)) {
    x_p0 <- stats::model.matrix(pred_p0, dados)
    J_p0 <- dim(x_p0)[2]
    pi <- par[(J_alpha + J_beta + 2):(J_alpha + J_beta + J_p0 + 1)]
    p0 <- exp(tcrossprod(x_p0, t(pi))) /
      (1 + exp(tcrossprod(x_p0, t(pi))))
  } else {
    p0 <- NULL
  }

  surv_func_partial <- purrr::partial(
    .f = surv_func,
    alpha = alpha,
    beta = beta,
    p0 = p0,
    theta = NULL
  )

  part_1 <- delta * log(
    pdf_gammag(
      y = t,
      mu = mu,
      distribution = distribution,
      alpha = alpha,
      beta = beta,
      p0 = p0,
      theta = NULL
    )
  )

  part_2 <- (1 - delta) * log(
    surv_gammag(
      y = t,
      mu = mu,
      surv_func = surv_func_partial
    )
  )

  out <- sum(part_1 + part_2)

  if (log) {
    lik <- out
  } else {
    lik <- exp(out)
  }

  return(lik)
}


# Log-likelihood for baseline defective models
log_like_baseline <- function(
    par,
    dados,
    delta,
    distribution,
    pred_alpha,
    pred_beta,
    pred_p0 = NULL,
    log = TRUE
) {

  surv_func <- switch(
    distribution,
    "dagum" = s_dd,
    "gompertz" = s_dg,
    stop("Invalid distribution. Use 'dagum' or 'gompertz'.")
  )

  pdf_func <- switch(
    distribution,
    "dagum" = d_dd,
    "gompertz" = d_dg,
    stop("Invalid distribution. Use 'dagum' or 'gompertz'.")
  )

  delta <- get(delta, dados)
  t <- as.numeric(stats::model.extract(stats::model.frame(pred_alpha, dados), "response"))

  # Model matrices
  x_alpha <- stats::model.matrix(pred_alpha, dados)
  x_beta  <- stats::model.matrix(pred_beta, dados)

  J_alpha <- dim(x_alpha)[2]
  J_beta  <- dim(x_beta)[2]

  # Parameter vector:
  # par = (alpha coefficients, beta coefficients, p0 coefficients if Dagum)
  aa <- par[1:J_alpha]
  be <- par[(J_alpha + 1):(J_alpha + J_beta)]

  # Linear predictors
  alpha <- tcrossprod(x_alpha, t(aa))
  if (identical(distribution, "dagum")) {
    alpha <- exp(alpha)
  }

  beta <- exp(tcrossprod(x_beta, t(be)))

  if (!is.null(pred_p0)) {
    x_p0 <- stats::model.matrix(pred_p0, dados)
    J_p0 <- dim(x_p0)[2]
    pi <- par[(J_alpha + J_beta + 1):(J_alpha + J_beta + J_p0)]
    p0 <- exp(tcrossprod(x_p0, t(pi))) /
      (1 + exp(tcrossprod(x_p0, t(pi))))
  } else {
    p0 <- NULL
  }

  part_1 <- delta * log(
    pdf_func(y = t, alpha = alpha, beta = beta, p0 = p0, theta = NULL)
  )

  part_2 <- (1 - delta) * log(
    surv_func(y = t, alpha = alpha, beta = beta, p0 = p0, theta = NULL)
  )

  out <- sum(part_1 + part_2)

  if (log) {
    lik <- out
  } else {
    lik <- exp(out)
  }

  return(lik)
}


################################################################################
# 4. Cure fraction functions
################################################################################

# Baseline cure fraction
calc_p0 <- function(alpha, beta, p0 = NULL, distribution) {
  switch(
    distribution,
    "dagum" = p0,
    "gompertz" = exp(beta / alpha),
    stop("Invalid distribution. Use 'dagum' or 'gompertz'.")
  )
}

# Gamma-G cure fraction
p0_gammag <- function(mu, alpha, beta, p0 = NULL, distribution) {
  p0_base <- calc_p0(
    alpha = alpha,
    beta = beta,
    p0 = p0,
    distribution = distribution
  )

  p0_gamma_g <- 1 - pgamma(-log(p0_base), shape = mu, lower.tail = TRUE)

  return(p0_gamma_g)
}


################################################################################
# 5. Random generation from the Gamma-G distribution
################################################################################

# Vectorized generation from Gamma-G distribution using inverse transform
generate_t_gammag_vector <- function(u1, mu, distribution, alpha, beta, p0 = NULL) {

  fnl_gammag <- function(t, u1, mu, distribution, alpha, beta, p0 = NULL) {
    surv_func <- switch(
      distribution,
      "dagum" = s_dd(t, alpha, beta, p0, log.p = FALSE),
      "gompertz" = s_dg(t, alpha, beta, p0, log.p = FALSE),
      stop("Invalid distribution. Use 'dagum' or 'gompertz'.")
    )

    pgamma(-log(surv_func), shape = mu, lower.tail = TRUE) - u1
  }

  solve_t <- function(single_u1) {
    tryCatch({
      lower <- 0.01
      upper <- 10000

      f_lower <- fnl_gammag(lower, single_u1, mu, distribution, alpha, beta, p0)
      f_upper <- fnl_gammag(upper, single_u1, mu, distribution, alpha, beta, p0)

      if (f_lower * f_upper > 0) {
        return(NA)
      }

      uniroot(
        fnl_gammag,
        c(lower, upper),
        tol = 0.0001,
        u1 = single_u1,
        mu = mu,
        distribution = distribution,
        alpha = alpha,
        beta = beta,
        p0 = p0
      )$root
    }, error = function(e) {
      return(NA)
    })
  }

  t_out <- sapply(u1, solve_t)

  return(t_out)
}


# Generation from Gamma-G distribution for vectors of alpha, beta and p0
generate_t_gammag <- function(mu, distribution, alpha, beta, p0 = NULL) {

  fnl_gammag <- function(t, single_u1, mu, distribution, single_alpha, single_beta, single_p0) {
    surv_func <- switch(
      distribution,
      "dagum" = s_dd(t, single_alpha, single_beta, single_p0, log.p = FALSE),
      "gompertz" = s_dg(t, single_alpha, single_beta, single_p0, log.p = FALSE),
      stop("Invalid distribution. Use 'dagum' or 'gompertz'.")
    )

    pgamma(-log(surv_func), shape = mu, lower.tail = TRUE) - single_u1
  }

  solve_t <- function(single_alpha, single_beta, single_p0) {
    max_attempts <- 1000
    attempt <- 1
    valid_result <- FALSE

    p0_base <- calc_p0(
      alpha = single_alpha,
      beta = single_beta,
      p0 = single_p0,
      distribution = distribution
    )

    p0_gamma_g <- 1 - pgamma(-log(p0_base), shape = mu, lower.tail = TRUE)

    while (!valid_result && attempt <= max_attempts) {
      single_u1 <- runif(1, min = 0, max = 1 - p0_gamma_g)

      tryCatch({
        lower <- 0.01
        upper <- 10000

        f_lower <- fnl_gammag(
          lower, single_u1, mu, distribution,
          single_alpha, single_beta, single_p0
        )

        f_upper <- fnl_gammag(
          upper, single_u1, mu, distribution,
          single_alpha, single_beta, single_p0
        )

        if (f_lower * f_upper > 0) {
          attempt <<- attempt + 1
          next
        }

        result <- uniroot(
          fnl_gammag,
          c(lower, upper),
          tol = 0.0001,
          single_u1 = single_u1,
          mu = mu,
          distribution = distribution,
          single_alpha = single_alpha,
          single_beta = single_beta,
          single_p0 = single_p0
        )$root

        valid_result <<- TRUE
        return(result)
      }, error = function(e) {
        attempt <<- attempt + 1
      })
    }

    if (!valid_result) {
      warning("Max attempts reached, returning NA")
      return(NA)
    }
  }

  # For Gompertz, p0 can be NULL. Use a vector of NULL-like values.
  if (is.null(p0)) {
    p0 <- rep(NA, length(alpha))
  }

  t_out <- mapply(solve_t, alpha, beta, p0)

  return(t_out)
}


################################################################################
# 6. Monte Carlo data generation
################################################################################

data_generate <- function(
    n,
    pc = 0.2,
    distribution,
    mu,
    a0,
    a1,
    b0,
    b1,
    v0 = NULL,
    v1 = NULL
) {

  compute_lambda <- function(a, b, v = NULL) {
    alpha <- if (distribution == "dagum") exp(a) else a
    beta <- exp(b)

    p0_base <- if (distribution == "dagum") {
      exp(v) / (1 + exp(v))
    } else {
      NULL
    }

    p0 <- calc_p0(
      alpha = alpha,
      beta = beta,
      p0 = p0_base,
      distribution = distribution
    )

    p0_gamma_g <- 1 - pgamma(-log(p0), shape = mu, lower.tail = TRUE)

    u1 <- runif(10000, min = 0, max = 1 - p0_gamma_g)

    simulated_times <- generate_t_gammag_vector(
      u1 = u1,
      mu = mu,
      distribution = distribution,
      alpha = alpha,
      beta = beta,
      p0 = p0_base
    )

    mean(simulated_times, na.rm = TRUE) / pc
  }

  lambda0 <- compute_lambda(a0, b0, v0)
  lambda1 <- compute_lambda(a0 + a1, b0 + b1, if (!is.null(v0)) v0 + v1 else NULL)

  cond <- 1

  while (cond > 0.8) {
    x <- rbinom(n, 1, 0.5)
    X <- stats::model.matrix(~ 1 + x)

    alpha <- X %*% c(a0, a1)
    if (distribution == "dagum") {
      alpha <- exp(alpha)
    }

    beta <- exp(X %*% c(b0, b1))

    p0_base <- if (distribution == "dagum") {
      exp(X %*% c(v0, v1)) / (1 + exp(X %*% c(v0, v1)))
    } else {
      NULL
    }

    p0 <- calc_p0(
      alpha = alpha,
      beta = beta,
      p0 = p0_base,
      distribution = distribution
    )

    p0_gamma_g <- 1 - pgamma(-log(p0), shape = mu, lower.tail = TRUE)

    # Generate susceptible survival times
    ts <- generate_t_gammag(
      mu = mu,
      distribution = distribution,
      alpha = alpha,
      beta = beta,
      p0 = p0_base
    )

    # Generate censoring times
    lambda <- ifelse(x == 0, lambda0, lambda1)
    cens <- runif(n, 0, lambda)

    # Long-term survivors
    u_star <- runif(n)
    w <- ifelse(u_star < p0_gamma_g, Inf, ts)

    t <- pmin(w, cens)
    d <- as.integer(w == t)

    # Check if the generated Kaplan-Meier plateaus are close enough
    mKM <- survival::survfit(survival::Surv(t, d) ~ x, se.fit = FALSE)

    p00 <- min(mKM[1]$surv)
    p01 <- min(mKM[2]$surv)

    if (a1 < 0) {
      cond <- max(abs(c(min(p0_gamma_g), max(p0_gamma_g)) - c(p00, p01)))
    } else {
      cond <- max(abs(c(max(p0_gamma_g), min(p0_gamma_g)) - c(p00, p01)))
    }
  }

  data.frame(t = t, d = d, x = x)
}


################################################################################
# 7. Inference utilities
################################################################################

# Numerical variance-covariance matrix, standard errors and confidence intervals
variances <- function(
    lik,
    est,
    dados,
    delta,
    distribution,
    pred_alpha,
    pred_beta,
    pred_p0 = NULL,
    log = TRUE
) {

  hes_numDeriv <- numDeriv::hessian(
    func = lik,
    x = est,
    dados = dados,
    delta = delta,
    distribution = distribution,
    pred_alpha = pred_alpha,
    pred_beta = pred_beta,
    pred_p0 = pred_p0,
    log = log
  )

  aux_numDeriv <- MASS::ginv(-hes_numDeriv)
  var_numDeriv <- diag(aux_numDeriv)

  cont_neg <- sapply(var_numDeriv, function(x) ifelse(x < 0, 1, 0))

  if (any(var_numDeriv < 0)) {
    var_numDeriv <- abs(var_numDeriv)
    diag(aux_numDeriv) <- var_numDeriv
  }

  li <- est - 1.96 * sqrt(var_numDeriv)
  ls <- est + 1.96 * sqrt(var_numDeriv)

  out <- list()
  out$covar <- aux_numDeriv
  out$ic <- data.frame(var = var_numDeriv, li = li, ls = ls, cont_neg = cont_neg)

  return(out)
}

# Univariate delta method
delta_uni <- function(f, est_mv, vari) {
  var_delta <- vari * (numDeriv::grad(f, est_mv))^2

  li <- f(est_mv) - 1.96 * sqrt(var_delta)
  ls <- f(est_mv) + 1.96 * sqrt(var_delta)

  out <- data.frame(var = var_delta, li = li, ls = ls)

  return(out)
}

# Multivariate delta method
delta_multi <- function(f, est_mv, vari) {
  grad_f <- numDeriv::grad(f, est_mv)
  var_delta <- t(grad_f) %*% vari %*% grad_f

  li <- f(est_mv) - 1.96 * sqrt(var_delta)
  ls <- f(est_mv) + 1.96 * sqrt(var_delta)

  out <- data.frame(var = var_delta, li = li, ls = ls)

  return(out)
}

# Safe wrapper for optim
try_optim <- function(...) {
  tryCatch(
    expr = stats::optim(...),
    error = function(e) NULL
  )
}
