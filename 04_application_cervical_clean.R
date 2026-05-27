################################################################################
# 04_application_cervical_clean.R
#
# Cervical cancer application for the paper:
# Defective Gamma-G family for cure fraction models
#
# Models fitted:
#   1. Defective Gamma-Dagum model
#   2. Defective Dagum baseline model
#
# Covariate:
#   Clinical stage grouped as I-II vs III-IV
#
# Expected data file:
#   data/cervical.txt
#
# Required columns:
#   tempo   = observed survival/censoring time
#   status  = event indicator (1 = event, 0 = censored)
#   ECGRUP  = clinical stage, with values including I, II, III, IV
################################################################################

#rm(list = ls())



################################################################################
# 1. Packages and functions
################################################################################

required_packages <- c(
  "numDeriv", "MASS", "survival", "survminer", "ggplot2",
  "dplyr", "purrr", "tibble", "broom", "scales"
)

invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Package not installed: ", pkg))
  }
  library(pkg, character.only = TRUE)
}))

# Load clean Gamma-G functions
source("R/functions_GammaG_clean.R")

if (!dir.exists("outputs")) {
  dir.create("outputs", recursive = TRUE)
}

################################################################################
# 2. Data
################################################################################

data_path <- "data/cervical.txt"

if (!file.exists(data_path)) {
  stop(
    paste0(
      "Data file not found: ", data_path, "\n",
      "Please place the cervical cancer dataset in the data/ folder with the name cervical.txt.\n",
      "The file must contain at least the columns: tempo, status, ECGRUP."
    )
  )
}

dados <- read.table(data_path, header = TRUE, sep = "\t")

required_cols <- c("tempo", "status", "ECGRUP")
missing_cols <- setdiff(required_cols, names(dados))

if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# Load clean Gamma-G functions
source("R/functions_GammaG_clean.R")


# Group stages I-II and III-IV
dados <- dados |>
  dplyr::mutate(
    ECGRUPAG = ifelse(ECGRUP %in% c("I", "II"), "I_II", "III_IV"),
    ECGRUPAG = factor(ECGRUPAG, levels = c("I_II", "III_IV")),
    x = ECGRUPAG,
    t = tempo
  ) |>
  dplyr::filter(!is.na(t), !is.na(status), !is.na(x), t > 0)

################################################################################
# 3. Kaplan-Meier curve by clinical stage
################################################################################

fit_km <- survival::survfit(survival::Surv(tempo, status) ~ ECGRUPAG, data = dados)

km_plot <- survminer::ggsurvplot(
  fit_km,
  data = dados,
  risk.table = TRUE,
  risk.table.col = "strata",
  surv.median.line = "none",
  conf.int = FALSE,
  palette = c("#1F77B4", "#D62728"),
  xlab = "Time (days)",
  ylab = "Survival probability (%)",
  legend.title = "Stage group",
  legend.labs = c("I-II", "III-IV"),
  risk.table.height = 0.25,
  risk.table.y.text.col = TRUE,
  ggtheme = ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
)

pdf("outputs/km_cervical_stage_risk_table.pdf", width = 8, height = 9)
print(km_plot)
dev.off()

################################################################################
# 4. Model specification
################################################################################

distribution <- "dagum"

pred_alpha <- t ~ x
pred_beta  <- t ~ x
pred_p0    <- t ~ x

x_alpha <- stats::model.matrix(pred_alpha, dados)
x_beta  <- stats::model.matrix(pred_beta, dados)
x_p0    <- stats::model.matrix(pred_p0, dados)

J_alpha <- ncol(x_alpha)
J_beta  <- ncol(x_beta)
J_p0    <- ncol(x_p0)

################################################################################
# 5. Defective Gamma-Dagum model
################################################################################

emv_dagum <- suppressWarnings(try_optim(
  par = c(2, rep(-0.5, J_alpha), rep(-0.5, J_beta), rep(-0.5, J_p0)),
  fn = function(par) {
    log_like(
      par = par,
      dados = dados,
      delta = "status",
      distribution = distribution,
      pred_alpha = pred_alpha,
      pred_beta = pred_beta,
      pred_p0 = pred_p0,
      log = TRUE
    )
  },
  control = list(fnscale = -1, maxit = 5000),
  method = "BFGS"
))

if (is.null(emv_dagum) || emv_dagum$convergence != 0) {
  stop("The defective Gamma-Dagum model did not converge.")
}

var_gamma_dagum <- try(
  variances(
    lik = log_like,
    est = emv_dagum$par,
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

if (inherits(var_gamma_dagum, "try-error")) {
  stop("Variance-covariance matrix could not be computed for the Gamma-Dagum model.")
}

# Parameter estimates
hat_mu <- exp(emv_dagum$par[1L])
hat_a0 <- emv_dagum$par[2L]
hat_a1 <- emv_dagum$par[3L]
hat_b0 <- emv_dagum$par[4L]
hat_b1 <- emv_dagum$par[5L]
hat_v0 <- emv_dagum$par[6L]
hat_v1 <- emv_dagum$par[7L]

estimates_gamma_dagum <- c(
  hat_mu, hat_a0, hat_a1, hat_b0, hat_b1, hat_v0, hat_v1
)

efe_mu <- function(mu) exp(mu)

out_mu <- cbind(
  emv = efe_mu(emv_dagum$par[1L]),
  delta_uni(
    f = efe_mu,
    est_mv = emv_dagum$par[1L],
    vari = var_gamma_dagum$ic[1, 1]
  )
)

estimativas_gamma_dagum <- cbind(
  emv = estimates_gamma_dagum,
  var_gamma_dagum$ic[, 1:3]
)

estimativas_gamma_dagum[1, ] <- out_mu

rownames(estimativas_gamma_dagum) <- c(
  "mu", "a0", "a1", "b0", "b1", "v0", "v1"
)

# Cure fractions under the Gamma-Dagum model
p0_pars_gamma_dagum <- c(emv_dagum$par[1L], hat_v0, hat_v1)
comb_covars <- as.matrix(expand.grid(inter = 1, x = 0:1))

efe_p0_gamma_dagum <- function(pars) {
  pars_p0 <- pars[2:3]
  p0_baseline <- exp(t(comb_c) %*% pars_p0) /
    (1 + exp(t(comb_c) %*% pars_p0))
  p0_gamma <- 1 - pgamma(-log(p0_baseline), shape = exp(pars[1]), lower.tail = TRUE)
  return(as.numeric(p0_gamma))
}

p0_estimates_gamma <- matrix(NA, nrow = nrow(comb_covars), ncol = 3)
p0_est_gamma <- numeric(nrow(comb_covars))

for (i in seq_len(nrow(comb_covars))) {
  comb_c <- comb_covars[i, ]
  aux <- delta_multi(
    f = efe_p0_gamma_dagum,
    est_mv = p0_pars_gamma_dagum,
    vari = var_gamma_dagum$covar[c(1, 6, 7), c(1, 6, 7)]
  )
  p0_estimates_gamma[i, 1:3] <- as.matrix(aux[1:3])
  p0_est_gamma[i] <- efe_p0_gamma_dagum(p0_pars_gamma_dagum)
}

p0_gamma_dagum <- data.frame(
  emv = p0_est_gamma,
  var = p0_estimates_gamma[, 1],
  li = p0_estimates_gamma[, 2],
  ls = p0_estimates_gamma[, 3]
)

rownames(p0_gamma_dagum) <- c("p00", "p01")

estimativas_gamma_dagum <- rbind(estimativas_gamma_dagum, p0_gamma_dagum)

loglik_gamma_dagum <- log_like(
  par = emv_dagum$par,
  dados = dados,
  delta = "status",
  distribution = distribution,
  pred_alpha = pred_alpha,
  pred_beta = pred_beta,
  pred_p0 = pred_p0,
  log = TRUE
)

k_gamma_dagum <- length(emv_dagum$par)
aic_gamma_dagum <- -2 * loglik_gamma_dagum + 2 * k_gamma_dagum
bic_gamma_dagum <- -2 * loglik_gamma_dagum + k_gamma_dagum * log(nrow(dados))

################################################################################
# 6. Defective Dagum baseline model
################################################################################

dagum_baseline <- suppressWarnings(try_optim(
  par = c(rep(-0.5, J_alpha), rep(-0.5, J_beta), rep(-0.5, J_p0)),
  fn = function(par) {
    log_like_baseline(
      par = par,
      dados = dados,
      delta = "status",
      distribution = distribution,
      pred_alpha = pred_alpha,
      pred_beta = pred_beta,
      pred_p0 = pred_p0,
      log = TRUE
    )
  },
  control = list(fnscale = -1, maxit = 5000),
  method = "BFGS"
))

if (is.null(dagum_baseline) || dagum_baseline$convergence != 0) {
  stop("The defective Dagum baseline model did not converge.")
}

var_dagum <- try(
  variances(
    lik = log_like_baseline,
    est = dagum_baseline$par,
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

if (inherits(var_dagum, "try-error")) {
  stop("Variance-covariance matrix could not be computed for the Dagum baseline model.")
}

estimativas_dagum_baseline <- cbind(
  emv = dagum_baseline$par,
  var_dagum$ic[, 1:3]
)

rownames(estimativas_dagum_baseline) <- c("a0", "a1", "b0", "b1", "v0", "v1")

# Cure fractions under the defective Dagum baseline model
p0_pars_dagum <- dagum_baseline$par[5:6]

efe_p0_dagum <- function(pars) {
  p0 <- exp(t(comb_c) %*% pars) / (1 + exp(t(comb_c) %*% pars))
  return(as.numeric(p0))
}

p0_estimates_dagum <- matrix(NA, nrow = nrow(comb_covars), ncol = 3)
p0_est_dagum <- numeric(nrow(comb_covars))

for (i in seq_len(nrow(comb_covars))) {
  comb_c <- comb_covars[i, ]
  aux <- delta_multi(
    f = efe_p0_dagum,
    est_mv = p0_pars_dagum,
    vari = var_dagum$covar[c(5, 6), c(5, 6)]
  )
  p0_estimates_dagum[i, 1:3] <- as.matrix(aux[1:3])
  p0_est_dagum[i] <- efe_p0_dagum(p0_pars_dagum)
}

p0_dagum <- data.frame(
  emv = p0_est_dagum,
  var = p0_estimates_dagum[, 1],
  li = p0_estimates_dagum[, 2],
  ls = p0_estimates_dagum[, 3]
)

rownames(p0_dagum) <- c("p00", "p01")

estimativas_dagum_baseline <- rbind(estimativas_dagum_baseline, p0_dagum)

loglik_dagum <- log_like_baseline(
  par = dagum_baseline$par,
  dados = dados,
  delta = "status",
  distribution = distribution,
  pred_alpha = pred_alpha,
  pred_beta = pred_beta,
  pred_p0 = pred_p0,
  log = TRUE
)

k_dagum <- length(dagum_baseline$par)
aic_dagum <- -2 * loglik_dagum + 2 * k_dagum
bic_dagum <- -2 * loglik_dagum + k_dagum * log(nrow(dados))

################################################################################
# 7. Model comparison and LRT
################################################################################

fit_statistics <- data.frame(
  Model = c("Defective Gamma-Dagum", "Defective Dagum"),
  logLik = c(loglik_gamma_dagum, loglik_dagum),
  AIC = c(aic_gamma_dagum, aic_dagum),
  BIC = c(bic_gamma_dagum, bic_dagum),
  n_parameters = c(k_gamma_dagum, k_dagum)
)

lrt <- -2 * (loglik_dagum - loglik_gamma_dagum)
p_value <- pchisq(lrt, df = 1, lower.tail = FALSE)

################################################################################
# 8. Correlation matrix and covariance diagnostics for Gamma-Dagum
################################################################################

V_raw_GD <- var_gamma_dagum$covar[1:7, 1:7]

J_GD <- diag(7)
J_GD[1, 1] <- exp(emv_dagum$par[1])

V_nat_GD <- J_GD %*% V_raw_GD %*% t(J_GD)

corr_gammadagum <- cov2cor(V_nat_GD)

rownames(corr_gammadagum) <- colnames(corr_gammadagum) <- c(
  "mu", "zeta0", "zeta1", "eta0", "eta1", "nu0", "nu1"
)

covariance_diagnostics <- data.frame(
  kappa_covariance = kappa(V_nat_GD),
  min_eigenvalue = min(eigen(V_nat_GD, symmetric = TRUE)$values),
  max_eigenvalue = max(eigen(V_nat_GD, symmetric = TRUE)$values)
)

################################################################################
# 9. Fitted survival curves: Kaplan-Meier and Gamma-Dagum
################################################################################

surv_gamma_dagum_curve <- function(t, group_value) {
  x_val <- ifelse(group_value == "I_II", 0, 1)

  alpha <- exp(hat_a0 + hat_a1 * x_val)
  beta  <- exp(hat_b0 + hat_b1 * x_val)
  eta_p0 <- hat_v0 + hat_v1 * x_val
  p0 <- exp(eta_p0) / (1 + exp(eta_p0))
  mu <- hat_mu

  s_dagum <- s_dd(t, alpha = alpha, beta = beta, p0 = p0)
  s_gamma_dagum <- 1 - pgamma(-log(s_dagum), shape = mu, lower.tail = TRUE)

  return(as.vector(s_gamma_dagum))
}

time_grid <- seq(0.01, max(dados$tempo), length.out = 300)

df_gamma <- dplyr::bind_rows(
  data.frame(
    time = time_grid,
    survival = surv_gamma_dagum_curve(time_grid, "I_II"),
    group = "I_II",
    model = "Gamma-Dagum"
  ),
  data.frame(
    time = time_grid,
    survival = surv_gamma_dagum_curve(time_grid, "III_IV"),
    group = "III_IV",
    model = "Gamma-Dagum"
  )
)

df_km <- survminer::surv_summary(fit_km, data = dados) |>
  dplyr::rename(time = time, survival = surv, group = strata) |>
  dplyr::mutate(
    group = gsub("ECGRUPAG=", "", group),
    model = "Kaplan-Meier"
  )

survival_plot <- ggplot2::ggplot() +
  ggplot2::geom_step(
    data = df_km,
    ggplot2::aes(x = time, y = survival, color = group),
    linewidth = 1.4
  ) +
  ggplot2::geom_line(
    data = df_gamma,
    ggplot2::aes(x = time, y = survival, color = group),
    linetype = "dashed",
    linewidth = 1.4
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  ggplot2::scale_color_manual(
    name = "Stage group",
    values = c("I_II" = "#1F77B4", "III_IV" = "#D62728"),
    labels = c("I-II", "III-IV")
  ) +
  ggplot2::labs(
    title = "",
    x = "Time (days)",
    y = "Survival probability (%)"
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(panel.grid = ggplot2::element_blank())

pdf("outputs/km_cervical_gamma_dagum.pdf", width = 8, height = 6)
print(survival_plot)
dev.off()

################################################################################
# 10. Estimated hazard functions for Gamma-Dagum
################################################################################

hazard_gamma_dagum_curve <- function(t, group_value) {
  x_val <- ifelse(group_value == "I_II", 0, 1)

  alpha <- exp(hat_a0 + hat_a1 * x_val)
  beta  <- exp(hat_b0 + hat_b1 * x_val)
  eta_p0 <- hat_v0 + hat_v1 * x_val
  p0 <- exp(eta_p0) / (1 + exp(eta_p0))
  mu <- hat_mu

  s_dagum <- s_dd(t, alpha = alpha, beta = beta, p0 = p0)
  f_gamma <- pdf_gammag(
    y = t,
    mu = mu,
    distribution = "dagum",
    alpha = alpha,
    beta = beta,
    p0 = p0
  )
  s_gamma <- 1 - pgamma(-log(s_dagum), shape = mu, lower.tail = TRUE)

  return(as.vector(f_gamma / s_gamma))
}

df_hazard <- dplyr::bind_rows(
  data.frame(
    time = time_grid,
    hazard = hazard_gamma_dagum_curve(time_grid, "I_II"),
    group = "I_II"
  ),
  data.frame(
    time = time_grid,
    hazard = hazard_gamma_dagum_curve(time_grid, "III_IV"),
    group = "III_IV"
  )
)

hazard_plot <- ggplot2::ggplot(df_hazard, ggplot2::aes(x = time, y = hazard, color = group)) +
  ggplot2::geom_line(linewidth = 1.6) +
  ggplot2::scale_color_manual(
    name = "Stage group",
    values = c("I_II" = "#1F77B4", "III_IV" = "#D62728"),
    labels = c("I-II", "III-IV")
  ) +
  ggplot2::labs(
    title = "",
    x = "Time (days)",
    y = "Hazard function h(t)"
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(panel.grid = ggplot2::element_blank())

pdf("outputs/hazard_cervical_gamma_dagum.pdf", width = 8, height = 6)
print(hazard_plot)
dev.off()

################################################################################
# 11. Save outputs
################################################################################

write.csv(
  estimativas_gamma_dagum,
  "outputs/estimates_gamma_dagum_cervical.csv",
  row.names = TRUE
)

write.csv(
  estimativas_dagum_baseline,
  "outputs/estimates_dagum_baseline_cervical.csv",
  row.names = TRUE
)

write.csv(
  fit_statistics,
  "outputs/model_comparison_cervical.csv",
  row.names = FALSE
)

write.csv(
  round(corr_gammadagum, 6),
  "outputs/correlation_matrix_gamma_dagum_cervical.csv",
  row.names = TRUE
)

write.csv(
  covariance_diagnostics,
  "outputs/covariance_diagnostics_gamma_dagum_cervical.csv",
  row.names = FALSE
)

################################################################################
# 12. Print summary
################################################################################

cat("\nDefective Gamma--Dagum estimates:\n")
print(round(estimativas_gamma_dagum, 4))

cat("\nDefective Dagum estimates:\n")
print(round(estimativas_dagum_baseline, 4))

cat("\nModel comparison:\n")
print(fit_statistics)

cat("\nLikelihood ratio test for Gamma--Dagum vs Dagum:\n")
cat("LRT =", round(lrt, 4), " p-value =", signif(p_value, 4), "\n")

cat("\nCorrelation matrix of Gamma--Dagum MLEs:\n")
print(round(corr_gammadagum, 3))

cat("\nCovariance diagnostics:\n")
print(covariance_diagnostics)

cat("\nOutputs saved in folder: outputs\n")
