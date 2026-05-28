################################################################################
# 01_simulation_gammagompertz_clean.R
#
# Monte Carlo simulation for the defective Gamma-Gompertz model used in the
# manuscript. The script generates the numerical summaries for Table 1 and the
# figures related to bias/RMSE, coverage, and information-criterion selection.
################################################################################

rm(list = ls())

################################################################################
# 1. Packages and functions
################################################################################

required_packages <- c(
  "numDeriv", "MASS", "survival", "ggplot2", "dplyr", "purrr",
  "tibble", "tidyr", "readr"
)

invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Package not installed: ", pkg))
  }
  library(pkg, character.only = TRUE)
}))

source("R/functions_GammaG_clean.R")

if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)

################################################################################
# 2. Simulation settings
################################################################################

# Full settings used for the manuscript.
# To make a quick test, set quick_run <- TRUE.


# IMPORTANT:
# Set quick_run = TRUE only to check that the script runs correctly.
# This quick version uses a reduced number of Monte Carlo replications and
# does not reproduce exactly the numerical results reported in the manuscript.
# To reproduce the full Monte Carlo study used in the paper, set quick_run = FALSE.
# The full simulation may take several hours depending on the computer.


quick_run <- TRUE#FALSE

n_mc <- if (quick_run) 20 else 1000
ns <- if (quick_run) c(100, 500) else c(100, 500, 1000, 2500, 5000)

mu_values <- c(1.0, 1.3)
main_mu <- 1.3

pc <- 0.20
sig <- 0.05
distribution <- "gompertz"

# Parameter values motivated by the melanoma application
a0_true <- -0.13
a1_true <- -0.23
b0_true <- -2.70
b1_true <-  2.20

################################################################################
# 3. One Monte Carlo replication
################################################################################

true_values_gompertz <- function(mu_true) {
  alpha0_true <- a0_true
  alpha1_true <- a0_true + a1_true

  beta0_true <- exp(b0_true)
  beta1_true <- exp(b0_true + b1_true)

  p00_true <- p0_gammag(
    mu = mu_true,
    alpha = alpha0_true,
    beta = beta0_true,
    distribution = distribution
  )

  p01_true <- p0_gammag(
    mu = mu_true,
    alpha = alpha1_true,
    beta = beta1_true,
    distribution = distribution
  )

  c(
    mu = mu_true,
    zeta0 = a0_true,
    zeta1 = a1_true,
    eta0 = b0_true,
    eta1 = b1_true,
    p00 = as.numeric(p00_true),
    p01 = as.numeric(p01_true)
  )
}

one_step_gompertz <- function(m, n, mu_true) {

  param_true <- true_values_gompertz(mu_true)

  for (attempt in seq_len(500)) {

    set.seed(10000 * m + attempt + round(1000 * mu_true))

    dados <- data_generate(
      n = n,
      pc = pc,
      distribution = distribution,
      mu = mu_true,
      a0 = a0_true,
      a1 = a1_true,
      b0 = b0_true,
      b1 = b1_true,
      v0 = NULL,
      v1 = NULL
    )

    pred_alpha <- t ~ x
    pred_beta <- t ~ x
    pred_p0 <- NULL

    x_alpha <- stats::model.matrix(pred_alpha, dados)
    x_beta <- stats::model.matrix(pred_beta, dados)

    J_alpha <- ncol(x_alpha)
    J_beta <- ncol(x_beta)

    # Gamma-Gompertz fit
    fit_gammag <- suppressWarnings(try_optim(
      par = c(log(mu_true), a0_true, a1_true, b0_true, b1_true),
      fn = function(par) log_like(
        par = par,
        dados = dados,
        delta = "d",
        distribution = distribution,
        pred_alpha = pred_alpha,
        pred_beta = pred_beta,
        pred_p0 = pred_p0,
        log = TRUE
      ),
      control = list(fnscale = -1, maxit = 3000),
      method = "BFGS"
    ))

    # Baseline defective Gompertz fit
    fit_baseline <- suppressWarnings(try_optim(
      par = c(a0_true, a1_true, b0_true, b1_true),
      fn = function(par) log_like_baseline(
        par = par,
        dados = dados,
        delta = "d",
        distribution = distribution,
        pred_alpha = pred_alpha,
        pred_beta = pred_beta,
        pred_p0 = pred_p0,
        log = TRUE
      ),
      control = list(fnscale = -1, maxit = 3000),
      method = "BFGS"
    ))

    if (is.null(fit_gammag) || is.null(fit_baseline)) next
    if (fit_gammag$convergence != 0L || fit_baseline$convergence != 0L) next

    var_aux <- try(
      variances(
        lik = log_like,
        est = fit_gammag$par,
        dados = dados,
        delta = "d",
        distribution = distribution,
        pred_alpha = pred_alpha,
        pred_beta = pred_beta,
        pred_p0 = pred_p0,
        log = TRUE
      ),
      silent = TRUE
    )

    if (inherits(var_aux, "try-error")) next
    if (any(!is.finite(var_aux$ic$var))) next
    if (sum(var_aux$ic$cont_neg) > 0) next

    est <- fit_gammag$par

    estimates <- c(
      mu = exp(est[1]),
      zeta0 = est[2],
      zeta1 = est[3],
      eta0 = est[4],
      eta1 = est[5]
    )

    # Delta-method CI for mu
    f_mu <- function(x) exp(x)
    out_mu <- delta_uni(
      f = f_mu,
      est_mv = est[1],
      vari = var_aux$ic[1, "var"]
    )

    ci <- var_aux$ic[, c("var", "li", "ls")]
    ci[1, ] <- out_mu

    # Cure fractions for x = 0 and x = 1
    p0_fun <- function(pars, x_value) {
      x_vec <- c(1, x_value)
      mu <- exp(pars[1])
      alpha <- sum(x_vec * pars[2:3])
      beta <- exp(sum(x_vec * pars[4:5]))

      p0_gammag(
        mu = mu,
        alpha = alpha,
        beta = beta,
        distribution = distribution
      )
    }

    p0_estimates <- matrix(NA_real_, nrow = 2, ncol = 3)

    for (j in 1:2) {
      x_value <- j - 1
      f_j <- function(pars) p0_fun(pars, x_value = x_value)

      aux <- delta_multi(
        f = f_j,
        est_mv = est,
        vari = var_aux$covar[1:5, 1:5]
      )

      p0_estimates[j, ] <- as.numeric(aux[1, c("var", "li", "ls")])
    }

    p00_hat <- p0_fun(est, 0)
    p01_hat <- p0_fun(est, 1)

    estimates <- c(estimates, p00 = as.numeric(p00_hat), p01 = as.numeric(p01_hat))

    ci_full <- rbind(
      ci,
      data.frame(
        var = p0_estimates[, 1],
        li = p0_estimates[, 2],
        ls = p0_estimates[, 3],
        row.names = c("p00", "p01")
      )
    )

    rownames(ci_full) <- names(estimates)

    coverage <- as.integer(
      param_true[names(estimates)] >= ci_full[names(estimates), "li"] &
        param_true[names(estimates)] <= ci_full[names(estimates), "ls"]
    )

    # Model comparison
    ll_gammag <- log_like(
      par = fit_gammag$par,
      dados = dados,
      delta = "d",
      distribution = distribution,
      pred_alpha = pred_alpha,
      pred_beta = pred_beta,
      pred_p0 = pred_p0,
      log = TRUE
    )

    ll_baseline <- log_like_baseline(
      par = fit_baseline$par,
      dados = dados,
      delta = "d",
      distribution = distribution,
      pred_alpha = pred_alpha,
      pred_beta = pred_beta,
      pred_p0 = pred_p0,
      log = TRUE
    )

    k_gammag <- length(fit_gammag$par)
    k_baseline <- length(fit_baseline$par)

    aic_gammag <- -2 * ll_gammag + 2 * k_gammag
    aic_baseline <- -2 * ll_baseline + 2 * k_baseline

    aicc_gammag <- aic_gammag + (2 * k_gammag * (k_gammag + 1)) / (n - k_gammag - 1)
    aicc_baseline <- aic_baseline + (2 * k_baseline * (k_baseline + 1)) / (n - k_baseline - 1)

    bic_gammag <- -2 * ll_gammag + k_gammag * log(n)
    bic_baseline <- -2 * ll_baseline + k_baseline * log(n)

    hqic_gammag <- -2 * ll_gammag + 2 * k_gammag * log(log(n))
    hqic_baseline <- -2 * ll_baseline + 2 * k_baseline * log(log(n))

    caic_gammag <- -2 * ll_gammag + k_gammag * (log(n) + 1)
    caic_baseline <- -2 * ll_baseline + k_baseline * (log(n) + 1)

    return(tibble::tibble(
      id_mc = m,
      attempt = attempt,
      n = n,
      mu_setting = mu_true,
      parameter = names(estimates),
      true = as.numeric(param_true[names(estimates)]),
      estimate = as.numeric(estimates),
      variance = as.numeric(ci_full[names(estimates), "var"]),
      li = as.numeric(ci_full[names(estimates), "li"]),
      ls = as.numeric(ci_full[names(estimates), "ls"]),
      coverage = coverage,
      ll_gammag = ll_gammag,
      ll_baseline = ll_baseline,
      aic_gammag = aic_gammag,
      aic_baseline = aic_baseline,
      aicc_gammag = aicc_gammag,
      aicc_baseline = aicc_baseline,
      bic_gammag = bic_gammag,
      bic_baseline = bic_baseline,
      hqic_gammag = hqic_gammag,
      hqic_baseline = hqic_baseline,
      caic_gammag = caic_gammag,
      caic_baseline = caic_baseline
    ))
  }

  warning(paste("No valid fit obtained for replication", m, "n =", n, "mu =", mu_true))

  tibble::tibble()
}

################################################################################
# 4. Run Monte Carlo simulation
################################################################################

simulation_grid <- expand.grid(
  mu_setting = mu_values,
  n = ns
)

raw_results <- purrr::map_dfr(seq_len(nrow(simulation_grid)), function(i) {
  mu_i <- simulation_grid$mu_setting[i]
  n_i <- simulation_grid$n[i]

  message("Running Gamma-Gompertz simulation: mu = ", mu_i, ", n = ", n_i)

  purrr::map_dfr(seq_len(n_mc), function(m) {
    one_step_gompertz(m = m, n = n_i, mu_true = mu_i)
  })
})

readr::write_csv(raw_results, "outputs/simulation_gammagompertz_raw.csv")
saveRDS(raw_results, "outputs/simulation_gammagompertz_raw.rds")

################################################################################
# 5. Numerical summaries for manuscript tables
################################################################################

summary_results <- raw_results |>
  dplyr::group_by(mu_setting, n, parameter, true) |>
  dplyr::summarise(
    mean = mean(estimate, na.rm = TRUE),
    bias = mean(estimate - true, na.rm = TRUE),
    rmse = sqrt(mean((estimate - true)^2, na.rm = TRUE)),
    coverage = mean(coverage, na.rm = TRUE),
    n_valid = dplyr::n_distinct(id_mc),
    .groups = "drop"
  )

table_gammagompertz <- summary_results |>
  dplyr::filter(mu_setting == main_mu) |>
  dplyr::arrange(factor(parameter, levels = c("mu", "zeta0", "zeta1", "eta0", "eta1", "p00", "p01")), n)

readr::write_csv(summary_results, "outputs/simulation_gammagompertz_summary_all_mu.csv")
readr::write_csv(table_gammagompertz, "outputs/simulation_gammagompertz_table1.csv")

################################################################################
# 6. Information-criterion selection summaries
################################################################################

criteria_results <- raw_results |>
  dplyr::distinct(
    id_mc, n, mu_setting,
    aic_gammag, aic_baseline,
    aicc_gammag, aicc_baseline,
    bic_gammag, bic_baseline,
    hqic_gammag, hqic_baseline,
    caic_gammag, caic_baseline
  ) |>
  dplyr::mutate(
    AIC = as.integer(aic_gammag < aic_baseline),
    AICc = as.integer(aicc_gammag < aicc_baseline),
    BIC = as.integer(bic_gammag < bic_baseline),
    HQIC = as.integer(hqic_gammag < hqic_baseline),
    CAIC = as.integer(caic_gammag < caic_baseline)
  ) |>
  tidyr::pivot_longer(
    cols = c(AIC, AICc, BIC, HQIC, CAIC),
    names_to = "criterion",
    values_to = "selected_gamma_g"
  ) |>
  dplyr::group_by(mu_setting, n, criterion) |>
  dplyr::summarise(
    selection_proportion = mean(selected_gamma_g, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(criteria_results, "outputs/simulation_gammagompertz_selection.csv")

################################################################################
# 7. Figures
################################################################################

plot_data <- summary_results |>
  dplyr::filter(parameter != "mu") |>
  dplyr::mutate(
    parameter = factor(parameter, levels = c("zeta0", "zeta1", "eta0", "eta1", "p00", "p01")),
    mu_label = paste0("mu = ", mu_setting)
  )

p_bias_rmse <- ggplot(plot_data, aes(x = n, y = bias, linetype = mu_label, shape = mu_label)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = bias - rmse, ymax = bias + rmse), width = 0) +
  geom_point(size = 2) +
  facet_wrap(~ parameter, scales = "free_y") +
  scale_x_continuous(breaks = ns) +
  labs(x = "n", y = "Bias", linetype = expression(mu), shape = expression(mu)) +
  theme_minimal()

ggsave("outputs/figure_bias_rmse_gammagompertz.pdf", p_bias_rmse, width = 9, height = 6)

p_coverage <- ggplot(plot_data, aes(x = n, y = coverage, linetype = mu_label, shape = mu_label)) +
  geom_hline(yintercept = 0.95) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ parameter) +
  scale_x_continuous(breaks = ns) +
  coord_cartesian(ylim = c(0.6, 1.0)) +
  labs(x = "n", y = "Coverage", linetype = expression(mu), shape = expression(mu)) +
  theme_minimal()

ggsave("outputs/figure_coverage_gammagompertz.pdf", p_coverage, width = 9, height = 6)

selection_plot_data <- criteria_results |>
  dplyr::filter(mu_setting == main_mu)

p_selection <- ggplot(selection_plot_data, aes(x = n, y = selection_proportion, linetype = criterion, shape = criterion)) +
  geom_hline(yintercept = 0.5) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = ns) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    x = "n",
    y = "Observed selection proportions based on information criteria",
    linetype = "Criterion",
    shape = "Criterion"
  ) +
  theme_minimal()

ggsave("outputs/figure_selection_gammagompertz.pdf", p_selection, width = 8, height = 5)

################################################################################
# 8. Print main summaries
################################################################################

cat("\nGamma-Gompertz simulation completed.\n")
cat("Raw results saved to: outputs/simulation_gammagompertz_raw.csv\n")
cat("Table 1 summary saved to: outputs/simulation_gammagompertz_table1.csv\n")
cat("Figures saved in: outputs/\n")
