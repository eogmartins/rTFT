# check-prior.R

# Functions to check whether the priors of a fitted failure-rate model pull its
# estimates. A prior that is narrow beside what the data say about a parameter draws
# the posterior towards its centre, and nothing in the convergence diagnostics shows
# it. The check scales the priors by a power, one block of them at a time, and
# measures how far each estimate moves (power-scaling sensitivity, after Kallioinen,
# Paananen, Burkner and Vehtari 2024): an estimate that moves when one block is
# scaled is pulled by that block.
#
# brms keeps the log prior of each draw as one sum, lprior. The check needs it by
# block: one block per parameter of the model (its intercept and its slopes), one per
# co-stressor (its main effect and its interaction with temperature together), and
# one for every prior on a group-level term. The priors of the population-level
# coefficients are read from the fit's own table of priors and evaluated at each
# draw; the group-level block is what is left of lprior, so that its priors
# (exponential, half Student-t, LKJ) need no evaluation of their own. Where the model
# has no group-level terms, what is left must be zero, which tests that the
# coefficients' priors were read as brms used them. The check also needs the
# log-likelihood of each draw summed over the rows of the data, which is computed in
# batches of draws so that no table of draws by rows is held at once.

# ===========================================================================
# 1. The log prior and the log-likelihood of each draw
# ===========================================================================

# .coef_prior() / .prior_lpdf() ------------------------------------------------

# The prior in force on one population-level coefficient, `coef` of the parameter
# `nlpar`, read from `pt`, the fit's table of priors: the prior set on that
# coefficient when there is one, else the prior set on all coefficients of the
# parameter, else the flat prior "". A bound on the coefficients is refused, as brms
# then adds to lprior a correction for the truncation that the check does not
# evaluate. Returns a character string.
.coef_prior <- function(pt, nlpar, coef) {
  rows <- pt$class == "b" & pt$nlpar == nlpar & pt$group == ""
  own  <- pt[rows & pt$coef == coef, , drop = FALSE]
  wide <- pt[rows & pt$coef == "", , drop = FALSE]
  bounds <- c(own$lb, own$ub, wide$lb, wide$ub)
  if (any(!is.na(bounds) & nzchar(bounds)))
    cli::cli_abort(c(
      "The coefficients of {.field {nlpar}} have a bound, which {.fn check_prior} cannot take into account.",
      "i" = "A bounded prior adds a correction for its truncation to the log prior of the model."
    ), class = "rTFT_prior_unsupported")
  if (nrow(own) && nzchar(own$prior[1])) return(own$prior[1])
  if (nrow(wide) && nzchar(wide$prior[1])) return(wide$prior[1])
  ""
}

# The log density of the prior string `prior` at the draws `x` of a coefficient, as
# brms adds it to lprior: a normal(mu, sigma), a student_t(nu, mu, sigma), or zero
# for the flat prior "", which brms leaves out of lprior. Any other prior is refused.
# `coef` names the coefficient in the message. Returns a numeric vector the length of
# `x`.
.prior_lpdf <- function(x, prior, coef) {
  if (identical(prior, "")) return(numeric(length(x)))
  p <- .prior_args(prior)
  a <- p$args
  if (!is.null(p) && p$family == "normal" && length(a) == 2L && a[[2]] > 0)
    return(stats::dnorm(x, a[[1]], a[[2]], log = TRUE))
  if (!is.null(p) && p$family == "student_t" && length(a) == 3L && all(a[-2] > 0))
    return(stats::dt((x - a[[2]]) / a[[3]], df = a[[1]], log = TRUE) - log(a[[3]]))
  cli::cli_abort(c(
    "{.fn check_prior} cannot evaluate the prior {.val {prior}} on {.field {coef}}.",
    "i" = "It evaluates normal and Student-t priors (and the flat prior) on the population-level coefficients."
  ), class = "rTFT_prior_unsupported")
}


# .prior_blocks() --------------------------------------------------------------

# The prefixes of the draws of terms other than linear and group-level ones (splines,
# monotonic effects, Gaussian processes), whose priors the check does not read.
.other_term_prefixes <- c("bs", "bsp", "bcs", "simo", "sds", "sdgp", "lscale")

# The log prior of each draw by block, from `draws`, the draws of `fit` as a data
# frame with its `lprior` column. Each population-level coefficient (a column
# `b_<parameter>_<coefficient>`) is assigned to the block of its parameter, or for a
# co-stressor's parameters (`c<name>` and `dT<name>`) to the block `c<name>`, and its
# prior (.coef_prior()) is evaluated at its draws (.prior_lpdf()). The block `sd` is
# the rest of lprior, present when the model has group-level terms (columns `sd_` or
# `cor_`); without them the rest must be zero to within 1e-6, or the priors were not
# read as brms used them and the check stops. Returns a list of `blocks` (a named
# list of numeric vectors, one value per draw) and `coefs` (a data frame with one
# row per coefficient: `variable`, `nlpar`, `coef`, `block` and `prior`).
.prior_blocks <- function(fit, draws) {
  vars <- names(draws)
  if (!"lprior" %in% vars)
    cli::cli_abort(c(
      "The draws of {.arg fit} hold no log prior ({.field lprior}).",
      "i" = "brms keeps it with the draws of every model it fits; draws taken apart or rebuilt may lose it. Refit the model with {.fn fit_hazard}."
    ), class = "rTFT_no_lprior")
  other <- vars[sub("_.*", "", vars) %in% .other_term_prefixes]
  if (length(other))
    cli::cli_abort(c(
      "{.fn check_prior} reads the priors of linear and group-level terms only.",
      "x" = "The model has {cli::qty(length(other))}{?a term/terms} of another kind: {.field {other}}."
    ), class = "rTFT_prior_unsupported")

  pt <- as.data.frame(fit$prior)
  nlpars <- unique(pt$nlpar[pt$class == "b" & nzchar(pt$nlpar)])
  costressor_block <- unlist(lapply(.fit_settings(fit)$costressors, \(cs) {
    dp <- c(paste0("c", cs$name), if (cs$interaction) paste0("dT", cs$name))
    stats::setNames(rep(paste0("c", cs$name), length(dp)), dp)
  }))
  coefs <- do.call(rbind, lapply(nlpars, \(nl) {
    v <- vars[startsWith(vars, paste0("b_", nl, "_"))]
    if (!length(v)) return(NULL)
    cf <- substring(v, nchar(nl) + 4L)
    data.frame(variable = v, nlpar = nl, coef = cf,
               block = if (nl %in% names(costressor_block)) costressor_block[[nl]] else nl,
               prior = vapply(cf, \(x) .coef_prior(pt, nl, x), character(1)),
               row.names = NULL)
  }))
  stray <- setdiff(vars[startsWith(vars, "b_")], coefs$variable)
  if (length(stray))
    cli::cli_abort(c(
      "{.fn check_prior} found no prior for {cli::qty(length(stray))}{?a coefficient/coefficients}: {.field {stray}}.",
      "i" = "Each coefficient should belong to a parameter of the fit's table of priors."
    ), class = "rTFT_prior_mismatch")

  lp <- lapply(seq_len(nrow(coefs)), \(j)
    .prior_lpdf(draws[[coefs$variable[j]]], coefs$prior[j], coefs$variable[j]))
  blocks <- lapply(split(seq_len(nrow(coefs)), factor(coefs$block, unique(coefs$block))),
                   \(j) Reduce(`+`, lp[j]))
  rest <- draws$lprior - Reduce(`+`, blocks)
  if (any(startsWith(vars, "sd_") | startsWith(vars, "cor_"))) {
    blocks$sd <- rest
  } else if (max(abs(rest)) > 1e-6) {
    cli::cli_abort(c(
      "The priors read from the fit do not add up to its log prior.",
      "x" = "They differ by up to {signif(max(abs(rest)), 3)} in a model with no group-level terms, where they should agree.",
      "i" = "A prior on a coefficient may have been changed after fitting, or set in a way the check does not read."
    ), class = "rTFT_prior_mismatch")
  }
  list(blocks = blocks, coefs = coefs)
}


# .summed_log_lik() ------------------------------------------------------------

# The log-likelihood of each draw of `fit`, summed over the rows of its data,
# computed `batch` draws at a time so that a table of draws by rows is never held
# whole (a fit of 4000 draws to 14,000 rows would hold 56 million values). Returns a
# numeric vector, one value per draw, in the order of the draws.
.summed_log_lik <- function(fit, batch = 250L) {
  n  <- brms::ndraws(fit)
  ll <- numeric(n)
  for (ids in split(seq_len(n), ceiling(seq_len(n) / batch)))
    ll[ids] <- rowSums(brms::log_lik(fit, draw_ids = ids))
  ll
}
