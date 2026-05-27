################################################################################
# Application: Melanoma dataset
# Defective Gamma--Gompertz model with one binary covariate (clinical stage)
#
# This script reproduces the melanoma application reported in the manuscript.
# It fits:
#   1) Defective Gamma--Gompertz model
#   2) Defective Gompertz baseline model
#
# Covariate:
#   x = EC_cat
#       0: Stage I--II / non-metastatic
#       1: Stage III--IV / metastatic
################################################################################

rm(list = ls())

################################################################################
# Packages
################################################################################

required_packages <- c(
  "numDeriv",
  "MASS",
  "survival",
  "survminer",
  "ggplot2",
  "broom",
  "dplyr",
  "purrr",
  "scales",
  "gridExtra"
)

missing_packages <- required_packages[!(required_packages %in% rownames(installed.packages()))]

if (length(missing_packages) > 0) {
  stop(
    "Please install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

################################################################################
# Paths
################################################################################

# This script assumes that it is run from the root folder of the repository.
# Expected structure:
#   R/functions_GammaG_clean.R
#   data/melanoma.txt
#   outputs/

functions_file <- "R/functions_GammaG_clean.R"
data_file <- "data/melanoma.txt"
output_dir <- "outputs"

if (!file.exists(functions_file)) {
  stop("Functions file not found: ", functions_file)
}

if (!file.exists(data_file)) {
  stop(
    "Data file not found: ", data_file,
    "\nPlease place the melanoma dataset in the data/ folder with the name melanoma.txt."
  )
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

source(functions_file)

################################################################################
# Data
################################################################################

dados <- read.table(data_file, header = TRUE)

# Required columns:
#   tempo_anos: observed survival/censoring time
#   status: event indicator, 1 = event, 0 = censored
#   EC_cat: clinical stage group, 0 = I--II, 1 = III--IV

required_columns <- c("tempo_anos", "status", "EC_cat")
missing_columns <- required_columns[!(required_columns %in% names(dados))]

if (length(missing_columns) > 0) {
  stop(
    "The following required columns are missing from the dataset: ",
    paste(missing_columns, collapse = ", ")
  )
}

dados <- dados |>
  dplyr::mutate(
    t = tempo_anos,
    x = EC_cat
  )

################################################################################
# Kaplan--Meier curve by clinical stage
################################################################################

fit_km <- survival::survfit(survival::Surv(tempo_anos, status) ~ EC_cat, data = dados)

km_plot_with_risk_table <- survminer::ggsurvplot(
  fit_km,
  data = dados,
  risk.table = TRUE,
  risk.table.col = "strata",
  surv.median.line = "none",
  conf.int = FALSE,
  palette = c("#0072B2", "#D55E00"),
  xlab = "Time (years)",
  ylab = "Survival probability (%)",
  legend.title = "Stage group",
  legend.labs = c("I--II", "III--IV"),
  risk.table.height = 0.25,
  risk.table.y.text.col = TRUE,
  ggtheme = ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
)

grDevices::pdf(file.path(output_dir, "km_melanoma_risk_table.pdf"), width = 8, height = 9)
gridExtra::grid.arrange(
  km_plot_with_risk_table$plot,
  km_plot_with_risk_table$table,
  ncol = 1,
  heights = c(2, 1)
)
grDevices::dev.off()

################################################################################
# Defective Gamma--Gompertz model
################################################################################

distribution <- "gompertz"

pred_alpha <- t ~ x
pred_beta <- t ~ x
pred_p0 <- NULL

x_alpha <- stats::model.matrix(pred_alpha, dados)
x_beta <- stats::model.matrix(pred_beta, dados)

J_alpha <- ncol(x_alpha)
J_beta <- ncol(x_beta)

# Initial values based on the final estimates used in the manuscript.
# Order: log(mu), alpha intercept, alpha x, beta intercept, beta x.
par0_gammagompertz <- c(
  log(1.295),
  -0.131,
  -0.226,
  -2.680,
  2.197
)

if (length(par0_gammagompertz) != 1 + J_alpha + J_beta) {
  stop("The length of par0_gammagompertz is not compatible with the model matrices.")
}

emv_gammagompertz <- suppressWarnings(try_optim(
  par0_gammagompertz,
  fn = \(par) log_like(
    par = par,
    dados = dados,
    delta = "status",
    distribution = distribution,
    pred_alpha = pred_alpha,
    pred_beta = pred_beta,
    pred_p0 = pred_p0,
    log = TRUE
  ),
  control = list(fnscale = -1, maxit = 5000),
  method = "BFGS"
))

if (is.null(emv_gammagompertz) || emv_gammagompertz$convergence != 0) {
  stop("The defective Gamma--Gompertz model did not converge.")
}

var_gammagompertz <- try(
  variances(
    lik = log_like,
    est = emv_gammagompertz$par,
    dados = dados,
    delta = "status",
    distribution = distribution,
    pred_alpha = pred_alpha,
    pred_beta = pred_beta,
    pred_p0 = pred_p0,
    log = TRUE
  ),
  silent = TRUE
)

if (inherits(var_gammagompertz, "try-error")) {
  stop("Could not compute the covariance matrix for the defective Gamma--Gompertz model.")
}

################################################################################
# Parameter estimates: Gamma--Gompertz
################################################################################

idx_mu <- 1
idx_alpha <- 2:(1 + J_alpha)
idx_beta <- (max(idx_alpha) + 1):(max(idx_alpha) + J_beta)

hat_mu <- exp(emv_gammagompertz$par[idx_mu])
hat_alpha <- emv_gammagompertz$par[idx_alpha]
hat_beta <- emv_gammagompertz$par[idx_beta]

estimates_gammagompertz <- c(hat_mu, hat_alpha, hat_beta)

names_alpha <- c("zeta0", "zeta1")
names_beta <- c("eta0", "eta1")

names(estimates_gammagompertz) <- c("mu", names_alpha, names_beta)

efe_mu <- function(mu) exp(mu)

out_mu <- cbind(
  emv = efe_mu(emv_gammagompertz$par[idx_mu]),
  delta_uni(
    f = efe_mu,
    est_mv = emv_gammagompertz$par[idx_mu],
    vari = var_gammagompertz$ic[idx_mu, idx_mu]
  )
)

estimates_table_gammagompertz <- cbind(
  emv = estimates_gammagompertz,
  var_gammagompertz$ic[, 1:3]
)

estimates_table_gammagompertz[1, ] <- out_mu
rownames(estimates_table_gammagompertz) <- names(estimates_gammagompertz)

################################################################################
# Cure fractions: Gamma--Gompertz
################################################################################

# Profiles:
#   x = 0: Stage I--II
#   x = 1: Stage III--IV

p0_pars_gammagompertz <- emv_gammagompertz$par
comb_covars <- as.matrix(expand.grid(intercept = 1, x = 0:1))

efe_p0_gammagompertz <- function(pars) {
  mu <- exp(pars[idx_mu])
  alpha <- as.numeric(t(comb_c) %*% pars[idx_alpha])
  beta <- exp(as.numeric(t(comb_c) %*% pars[idx_beta]))

  if (alpha >= 0) {
    return(0)
  }

  p0_baseline <- exp(beta / alpha)
  p0_gammag <- 1 - pgamma(-log(p0_baseline), shape = mu, lower.tail = TRUE)

  return(p0_gammag)
}

p0_estimates <- matrix(NA, nrow = nrow(comb_covars), ncol = 3)
p0_est <- numeric(nrow(comb_covars))

for (i in seq_len(nrow(comb_covars))) {
  comb_c <- comb_covars[i, ]

  aux <- delta_multi(
    f = efe_p0_gammagompertz,
    est_mv = p0_pars_gammagompertz,
    vari = var_gammagompertz$covar[
      seq_along(p0_pars_gammagompertz),
      seq_along(p0_pars_gammagompertz)
    ]
  )

  p0_estimates[i, 1:3] <- as.matrix(aux[1:3])
  p0_est[i] <- efe_p0_gammagompertz(p0_pars_gammagompertz)
}

p0_table_gammagompertz <- data.frame(
  emv = p0_est,
  var = p0_estimates[, 1],
  li = p0_estimates[, 2],
  ls = p0_estimates[, 3]
)

rownames(p0_table_gammagompertz) <- c("p00", "p01")

estimates_table_gammagompertz <- rbind(
  estimates_table_gammagompertz,
  p0_table_gammagompertz
)

################################################################################
# Log-likelihood, AIC, and BIC: Gamma--Gompertz
################################################################################

loglik_gammagompertz <- log_like(
  par = emv_gammagompertz$par,
  dados = dados,
  delta = "status",
  distribution = distribution,
  pred_alpha = pred_alpha,
  pred_beta = pred_beta,
  pred_p0 = pred_p0,
  log = TRUE
)

k_gammagompertz <- length(emv_gammagompertz$par)
n_obs <- nrow(dados)

aic_gammagompertz <- -2 * loglik_gammagompertz + 2 * k_gammagompertz
bic_gammagompertz <- -2 * loglik_gammagompertz + k_gammagompertz * log(n_obs)

################################################################################
# Estimated correlation matrix: Gamma--Gompertz
################################################################################

V_raw <- var_gammagompertz$covar[seq_len(k_gammagompertz), seq_len(k_gammagompertz)]

J <- diag(k_gammagompertz)
J[1, 1] <- exp(emv_gammagompertz$par[1])

V_nat <- J %*% V_raw %*% t(J)

corr_gammagompertz <- cov2cor(V_nat)
rownames(corr_gammagompertz) <- colnames(corr_gammagompertz) <- c(
  "mu", "zeta0", "zeta1", "eta0", "eta1"
)

kappa_cov_gammagompertz <- kappa(V_nat)
eig_cov_gammagompertz <- eigen(V_nat, symmetric = TRUE)$values

################################################################################
# Defective Gompertz baseline model
################################################################################

par0_gompertz <- c(
  -0.065,
  -0.215,
  -3.499,
  2.555
)

if (length(par0_gompertz) != J_alpha + J_beta) {
  stop("The length of par0_gompertz is not compatible with the model matrices.")
}

emv_gompertz <- suppressWarnings(try_optim(
  par0_gompertz,
  fn = \(par) log_like_baseline(
    par = par,
    dados = dados,
    delta = "status",
    distribution = distribution,
    pred_alpha = pred_alpha,
    pred_beta = pred_beta,
    pred_p0 = pred_p0,
    log = TRUE
  ),
  control = list(fnscale = -1, maxit = 5000),
  method = "BFGS"
))

if (is.null(emv_gompertz) || emv_gompertz$convergence != 0) {
  stop("The defective Gompertz baseline model did not converge.")
}

var_gompertz <- try(
  variances(
    lik = log_like_baseline,
    est = emv_gompertz$par,
    dados = dados,
    delta = "status",
    distribution = distribution,
    pred_alpha = pred_alpha,
    pred_beta = pred_beta,
    pred_p0 = pred_p0,
    log = TRUE
  ),
  silent = TRUE
)

if (inherits(var_gompertz, "try-error")) {
  stop("Could not compute the covariance matrix for the defective Gompertz baseline model.")
}

################################################################################
# Parameter estimates: defective Gompertz baseline
################################################################################

estimates_gompertz <- emv_gompertz$par
names(estimates_gompertz) <- c("zeta0", "zeta1", "eta0", "eta1")

estimates_table_gompertz <- cbind(
  emv = estimates_gompertz,
  var_gompertz$ic[, 1:3]
)

rownames(estimates_table_gompertz) <- names(estimates_gompertz)

################################################################################
# Cure fractions: defective Gompertz baseline
################################################################################

p0_pars_gompertz <- emv_gompertz$par
comb_covars <- as.matrix(expand.grid(intercept = 1, x = 0:1))

efe_p0_gompertz <- function(pars) {
  alpha <- as.numeric(t(comb_c) %*% pars[1:J_alpha])
  beta <- exp(as.numeric(t(comb_c) %*% pars[(J_alpha + 1):(J_alpha + J_beta)]))

  if (alpha >= 0) {
    return(0)
  }

  p0 <- exp(beta / alpha)

  return(p0)
}

p0_estimates_gompertz <- matrix(NA, nrow = nrow(comb_covars), ncol = 3)
p0_est_gompertz <- numeric(nrow(comb_covars))

for (i in seq_len(nrow(comb_covars))) {
  comb_c <- comb_covars[i, ]

  aux <- delta_multi(
    f = efe_p0_gompertz,
    est_mv = p0_pars_gompertz,
    vari = var_gompertz$covar[
      seq_along(p0_pars_gompertz),
      seq_along(p0_pars_gompertz)
    ]
  )

  p0_estimates_gompertz[i, 1:3] <- as.matrix(aux[1:3])
  p0_est_gompertz[i] <- efe_p0_gompertz(p0_pars_gompertz)
}

p0_table_gompertz <- data.frame(
  emv = p0_est_gompertz,
  var = p0_estimates_gompertz[, 1],
  li = p0_estimates_gompertz[, 2],
  ls = p0_estimates_gompertz[, 3]
)

rownames(p0_table_gompertz) <- c("p00", "p01")

estimates_table_gompertz <- rbind(
  estimates_table_gompertz,
  p0_table_gompertz
)

################################################################################
# Log-likelihood, AIC, and BIC: defective Gompertz
################################################################################

loglik_gompertz <- log_like_baseline(
  par = emv_gompertz$par,
  dados = dados,
  delta = "status",
  distribution = distribution,
  pred_alpha = pred_alpha,
  pred_beta = pred_beta,
  pred_p0 = pred_p0,
  log = TRUE
)

k_gompertz <- length(emv_gompertz$par)

aic_gompertz <- -2 * loglik_gompertz + 2 * k_gompertz
bic_gompertz <- -2 * loglik_gompertz + k_gompertz * log(n_obs)

################################################################################
# Likelihood ratio test: defective Gompertz vs defective Gamma--Gompertz
################################################################################

lrt_mu1 <- -2 * (loglik_gompertz - loglik_gammagompertz)
p_value_mu1 <- pchisq(lrt_mu1, df = 1, lower.tail = FALSE)

################################################################################
# Save numerical outputs
################################################################################

model_comparison <- data.frame(
  Model = c("Defective Gamma-Gompertz", "Defective Gompertz"),
  logLik = c(loglik_gammagompertz, loglik_gompertz),
  AIC = c(aic_gammagompertz, aic_gompertz),
  BIC = c(bic_gammagompertz, bic_gompertz),
  n_parameters = c(k_gammagompertz, k_gompertz)
)

write.csv(
  estimates_table_gammagompertz,
  file.path(output_dir, "melanoma_estimates_gamma_gompertz.csv"),
  row.names = TRUE
)

write.csv(
  estimates_table_gompertz,
  file.path(output_dir, "melanoma_estimates_gompertz.csv"),
  row.names = TRUE
)

write.csv(
  model_comparison,
  file.path(output_dir, "melanoma_model_comparison.csv"),
  row.names = FALSE
)

write.csv(
  round(corr_gammagompertz, 6),
  file.path(output_dir, "melanoma_correlation_matrix_gamma_gompertz.csv"),
  row.names = TRUE
)

diagnostics_gammagompertz <- data.frame(
  kappa_covariance = kappa_cov_gammagompertz,
  min_eigenvalue = min(eig_cov_gammagompertz),
  max_eigenvalue = max(eig_cov_gammagompertz)
)

write.csv(
  diagnostics_gammagompertz,
  file.path(output_dir, "melanoma_diagnostics_gamma_gompertz.csv"),
  row.names = FALSE
)

################################################################################
# Plot: Kaplan--Meier and fitted Gamma--Gompertz curves
################################################################################

km_data <- broom::tidy(fit_km) |>
  dplyr::mutate(
    group = dplyr::case_when(
      strata == "EC_cat=0" ~ "I--II",
      strata == "EC_cat=1" ~ "III--IV",
      TRUE ~ as.character(strata)
    )
  )

time_grid <- seq(0.01, max(dados$tempo_anos), length.out = 300)

mu_hat <- estimates_table_gammagompertz["mu", "emv"]
a0_hat <- estimates_table_gammagompertz["zeta0", "emv"]
a1_hat <- estimates_table_gammagompertz["zeta1", "emv"]
b0_hat <- estimates_table_gammagompertz["eta0", "emv"]
b1_hat <- estimates_table_gammagompertz["eta1", "emv"]

curve_gamma_gompertz <- function(t, group_value) {
  alpha <- a0_hat + a1_hat * group_value
  beta <- exp(b0_hat + b1_hat * group_value)
  surv_gammag(
    y = t,
    mu = mu_hat,
    surv_func = purrr::partial(s_dg, alpha = alpha, beta = beta)
  )
}

fitted_curves <- data.frame(
  time = rep(time_grid, 2),
  group = factor(rep(c(0, 1), each = length(time_grid)), labels = c("I--II", "III--IV")),
  survival = c(
    curve_gamma_gompertz(time_grid, 0),
    curve_gamma_gompertz(time_grid, 1)
  )
)

plot_km_fitted <- ggplot2::ggplot() +
  ggplot2::geom_step(
    data = km_data,
    ggplot2::aes(x = time, y = estimate, color = group),
    linewidth = 1.2
  ) +
  ggplot2::geom_line(
    data = fitted_curves,
    ggplot2::aes(x = time, y = survival, color = group),
    linetype = "dashed",
    linewidth = 1.2
  ) +
  ggplot2::labs(
    x = "Time (years)",
    y = "Survival probability (%)",
    color = "Stage group"
  ) +
  ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
  ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00")) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(panel.grid = ggplot2::element_blank())

grDevices::pdf(file.path(output_dir, "km_melanoma_gamma_gompertz_fitted.pdf"), width = 8, height = 6)
print(plot_km_fitted)
grDevices::dev.off()

################################################################################
# Plot: estimated hazard functions under Gamma--Gompertz
################################################################################

hazard_gamma_gompertz <- function(t, group_value) {
  alpha <- a0_hat + a1_hat * group_value
  beta <- exp(b0_hat + b1_hat * group_value)

  f_t <- pdf_gammag(
    y = t,
    mu = mu_hat,
    distribution = "gompertz",
    alpha = alpha,
    beta = beta
  )

  s_t <- surv_gammag(
    y = t,
    mu = mu_hat,
    surv_func = purrr::partial(s_dg, alpha = alpha, beta = beta)
  )

  f_t / s_t
}

hazard_data <- data.frame(
  time = rep(time_grid, 2),
  group = factor(rep(c(0, 1), each = length(time_grid)), labels = c("I--II", "III--IV")),
  hazard = c(
    hazard_gamma_gompertz(time_grid, 0),
    hazard_gamma_gompertz(time_grid, 1)
  )
)

plot_hazard <- ggplot2::ggplot(
  hazard_data,
  ggplot2::aes(x = time, y = hazard, color = group)
) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::labs(
    x = "Time (years)",
    y = "Hazard function h(t)",
    color = "Stage group"
  ) +
  ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00")) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(panel.grid = ggplot2::element_blank())

grDevices::pdf(file.path(output_dir, "hazard_melanoma_gamma_gompertz.pdf"), width = 8, height = 6)
print(plot_hazard)
grDevices::dev.off()

################################################################################
# Print summary
################################################################################

cat("\nDefective Gamma--Gompertz estimates:\n")
print(round(estimates_table_gammagompertz, 4))

cat("\nDefective Gompertz estimates:\n")
print(round(estimates_table_gompertz, 4))

cat("\nModel comparison:\n")
print(model_comparison)

cat("\nLikelihood ratio test for Gamma--Gompertz vs Gompertz:\n")
cat("LRT =", round(lrt_mu1, 4), " p-value =", signif(p_value_mu1, 4), "\n")

cat("\nCorrelation matrix of Gamma--Gompertz MLEs:\n")
print(round(corr_gammagompertz, 3))

cat("\nCovariance diagnostics:\n")
print(diagnostics_gammagompertz)

cat("\nOutputs saved in folder:", output_dir, "\n")
