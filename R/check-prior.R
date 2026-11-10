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
# `cor_`); without them the rest must be zero to within what the rounding of the
# stored draws allows, or the priors were not read as brms used them and the check
# stops. Returns a list of `blocks` (a named
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
  # The draws of a cmdstanr fit are stored to eight significant figures, so the
  # rebuilt log prior differs from lprior by rounding, in proportion to its size.
  tol  <- 1e-5 * max(1, max(abs(draws$lprior)))
  if (any(startsWith(vars, "sd_") | startsWith(vars, "cor_"))) {
    blocks$sd <- rest
  } else if (max(abs(rest)) > tol) {
    cli::cli_abort(c(
      "The priors read from the fit do not add up to its log prior.",
      "x" = "They differ by up to {signif(max(abs(rest)), 3)}, more than the {signif(tol, 3)} that rounding of the stored draws allows, in a model with no group-level terms, where they should agree.",
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


# ===========================================================================
# 2. The check
# ===========================================================================

# .has_priorsense() ------------------------------------------------------------

# Checks if priorsense is installed: a function of its own so that a test can stand
# in for a library without it.
.has_priorsense <- function() requireNamespace("priorsense", quietly = TRUE)


# .check_quantities() ----------------------------------------------------------

# The quantities whose sensitivity is measured, as one vector of draws each, in the
# order of the draws of `fit` (`draws` as a data frame): every population-level
# coefficient in `coefs` (from .prior_blocks()), and with `newdata` the reported
# value of every parameter at each of its rows (.param_draws(), with random effects
# at zero), with sCTmax at fraction `p` within `t_ref` as well when `t_ref` is
# given. `levels` sets the co-stressors' levels for those. Returns a list of
# `values` (a list of numeric vectors) and `index` (a tibble with one row per
# quantity: `quantity`, its name, and `row`, its row of `newdata` or NA for a
# coefficient).
.check_quantities <- function(fit, draws, coefs, newdata, p, t_ref, levels) {
  values <- lapply(coefs$variable, \(v) draws[[v]])
  index  <- tibble::tibble(quantity = coefs$variable, row = NA_integer_)
  if (is.null(newdata)) return(list(values = values, index = index))

  st   <- .fit_settings(fit)
  cs   <- st$costressors
  pars <- c(M_catalogue[[st$m_form]]$dpars, h0_catalogue[[st$baseline]]$dpars,
            .costressor_dpars(cs))
  dr <- .fit_draws(fit, pars, newdata)
  k  <- if (M_catalogue[[st$m_form]]$k_sw != 0) st$k_sw
  .rate_setup(st$m_form, purrr::map(dr, \(m) m[, 1]), st$Tr, k, cs, levels)
  cs_names <- vapply(cs, `[[`, character(1), "name")
  levels   <- if (length(cs)) levels[cs_names]
  z_at     <- stats::setNames(purrr::map(cs, \(x) .costressor_z(x, levels[[x$name]])),
                              cs_names)
  H0 <- H0_of(st$baseline)
  for (r in seq_len(ncol(dr[[1]]))) {
    d  <- purrr::map(dr, \(m) m[, r])
    at <- .param_draws(st$m_form, st$baseline, d, cs, z_at, st$tau_ref)
    if (!is.null(t_ref)) {
      bl <- .baseline_draws(st$baseline, d, length(d[[1]]), st$tau_ref)
      target <- -log1p(-p) / H0(t_ref, bl$shape, bl$scale)
      at$sCTmax <- .sctmax_draws(st$m_form, d, target, st$Tr, k, cs, levels)
    }
    values <- c(values, unname(at))
    index  <- dplyr::bind_rows(index, tibble::tibble(quantity = names(at), row = r))
  }
  list(values = values, index = index)
}


# .sensitivity_table() / .prior_table() ----------------------------------------

# The power-scaling sensitivity of each quantity in `values` to each block of
# `blocks` (from .prior_blocks()) and to all the priors together (the block `all`),
# with the sensitivity to the likelihood, from priorsense at its defaults (alpha 0.99
# and 1.01, the cumulative Jensen-Shannon distance). `lprior` and `log_lik` are the
# log prior and the summed log-likelihood of each draw, and `draws` carries the
# chain, iteration and draw of each. A quantity with one or more infinite draws (an
# sCTmax that some draws never reach) gets NA, with a message, since a distance
# between distributions with infinite values is undefined. A block whose priors are
# all flat is the same at every draw, so there is nothing to scale: it is left out,
# with a message, and so is `all` when every prior is flat. Each scaling is an
# importance-sampling estimate, and the Pareto k of its weights says whether they
# can be trusted (Kallioinen et al. 2024, section 2.3): the larger k of a block's
# two scalings, and of the likelihood's two, is returned, with a warning when one is
# above priorsense's threshold, min(1 - 1 / log10(S), 0.7) for S draws. Returns a
# tibble with one row per block and quantity: `block`, `q` (the quantity's position
# in `values`), `prior`, `likelihood`, `diagnosis`, `k_prior` and `k_likelihood`.
.sensitivity_table <- function(values, blocks, lprior, log_lik, draws) {
  finite <- vapply(values, \(x) all(is.finite(x)), logical(1))
  if (!all(finite))
    cli::cli_inform(c(
      "i" = "Some draws of {sum(!finite)} quantit{?y/ies} are infinite (an sCTmax never reached), so {?its/their} sensitivity is given as NA."
    ), class = "rTFT_unreached")
  flat <- vapply(blocks, \(x) all(x == x[1]), logical(1))
  if (any(flat))
    cli::cli_inform(c(
      "i" = "The priors of the block{?s} {.field {names(blocks)[flat]}} are all flat, so there is nothing to scale and {?it is/they are} left out."
    ), class = "rTFT_flat_block")
  scaled <- c(if (!all(lprior == lprior[1])) "all", names(blocks)[!flat])
  empty  <- tibble::tibble(block = character(), q = integer(), prior = numeric(),
                           likelihood = numeric(), diagnosis = character(),
                           k_prior = numeric(), k_likelihood = numeric())
  if (!length(scaled)) return(empty)
  q_names <- paste0("q", seq_along(values))[finite]
  tab <- data.frame(stats::setNames(values[finite], q_names), lprior = lprior,
                    stats::setNames(blocks, paste0("lprior_", names(blocks))),
                    log_lik = log_lik, .chain = draws$.chain,
                    .iteration = draws$.iteration, .draw = draws$.draw,
                    check.names = FALSE)
  x   <- posterior::as_draws_df(tab)
  psd <- priorsense::create_priorsense_data(
    posterior::subset_draws(x, variable = q_names),
    log_prior = posterior::subset_draws(x, variable = c("lprior", paste0("lprior_", names(blocks)))),
    log_lik   = posterior::subset_draws(x, variable = "log_lik"))
  # the larger Pareto k of the two scalings of one component (the block `selection`
  # of the prior, or the likelihood), with priorsense's threshold for it
  k_of <- function(component, selection = NULL) {
    d <- lapply(c(0.99, 1.01), \(a) attr(priorsense::powerscale(
      psd, component = component, alpha = a, selection = selection),
      "powerscaling")$diagnostics)
    c(k = max(d[[1]]$khat, d[[2]]$khat), threshold = d[[1]]$khat_threshold)
  }
  k_lik <- k_of("likelihood")
  out <- dplyr::bind_rows(lapply(scaled, \(b) {
    sel <- if (b != "all") b
    s <- as.data.frame(priorsense::powerscale_sensitivity(psd, prior_selection = sel))
    out <- tibble::tibble(block = b, q = seq_along(values), prior = NA_real_,
                          likelihood = NA_real_, diagnosis = NA_character_,
                          k_prior = k_of("prior", sel)[["k"]],
                          k_likelihood = k_lik[["k"]])
    at <- match(s$variable, paste0("q", seq_along(values)))
    out$prior[at]      <- s$prior
    out$likelihood[at] <- s$likelihood
    out$diagnosis[at]  <- as.character(s$diagnosis)
    out
  }))
  threshold <- k_lik[["threshold"]]
  bad <- c(unique(out$block[out$k_prior > threshold]),
           if (k_lik[["k"]] > threshold) "likelihood")
  if (length(bad))
    cli::cli_warn(c(
      "!" = "The power-scaling of {.field {bad}} is unreliable: {?its/their} Pareto k is above {signif(threshold, 2)}.",
      "i" = "The sensitivities of {cli::qty(length(bad))}{?that/those} scaling{?s} rest on importance weights that a few draws dominate. Read them as a sign at most, and test a flagged block by refitting with its priors widened."
    ), class = "rTFT_unreliable_scaling")
  out
}

# The prior of each population-level coefficient in `coefs` (from .prior_blocks())
# beside its posterior in `draws`: the centre and scale of a normal or Student-t
# prior (NA for the flat prior), the posterior median and standard deviation, the
# distance of the median from the prior's centre in units of the prior's scale, and
# the ratio of the posterior standard deviation to that scale. Returns a tibble with
# one row per coefficient: `parameter`, `coefficient`, `prior`, `prior_centre`,
# `prior_scale`, `median`, `sd`, `distance` and `sd_ratio`.
.prior_table <- function(coefs, draws) {
  centre_scale <- vapply(coefs$prior, \(pr) {
    p <- .prior_args(pr)
    if (is.null(p)) return(c(NA_real_, NA_real_))
    switch(p$family, normal = p$args[1:2], student_t = p$args[2:3], c(NA_real_, NA_real_))
  }, numeric(2), USE.NAMES = FALSE)
  med <- vapply(coefs$variable, \(v) stats::median(draws[[v]]), numeric(1), USE.NAMES = FALSE)
  sds <- vapply(coefs$variable, \(v) stats::sd(draws[[v]]), numeric(1), USE.NAMES = FALSE)
  tibble::tibble(parameter = coefs$nlpar, coefficient = coefs$coef, prior = coefs$prior,
                 prior_centre = centre_scale[1, ], prior_scale = centre_scale[2, ],
                 median = med, sd = sds,
                 distance = (med - centre_scale[1, ]) / centre_scale[2, ],
                 sd_ratio = sds / centre_scale[2, ])
}


# check_prior() ----------------------------------------------------------------

#' Check whether the priors of a fitted model pull its estimates
#'
#' Measures how much each estimate of a model fitted with [fit_hazard()] depends on
#' its priors, block by block, by power-scaling (Kallioinen et al. 2024): each block
#' of priors is raised to a power slightly above and below 1, the posterior is
#' reweighted to match, and the change in each estimate's distribution is measured.
#' An estimate that moves when one block is scaled is pulled by that block. A narrow
#' prior can pull an estimate with no sign in the convergence diagnostics. That is
#' most likely in a design whose groups fail far from the reference temperature or
#' far from each other. log a is the log failure rate at the reference temperature,
#' so a group that fails some degrees away has its log a extrapolated from its
#' failures over that distance, a change of about b times the distance. A prior
#' narrow beside that change draws log a back towards the prior's centre, and the
#' group's other parameters move to fit its failures.
#'
#' The model is not refitted. The check needs the log prior of each draw, which
#' brms keeps, and the log-likelihood of each draw, which it computes once over all
#' draws (about half a minute for 4000 draws of a model fitted to 14,000 rows).
#'
#' @section Blocks of priors:
#' The priors are scaled in blocks: one per parameter of the model (`loga`, `logb`,
#' `Tc`, the baseline's parameters), holding its intercept and its slopes; one per
#' co-stressor, named `c<name>` after its main effect, holding the main effect and
#' the interaction with temperature; `sd`, holding every prior on a group-level term
#' (standard deviations and correlations); and `all`, every prior together. The
#' priors of the population-level coefficients are read from the fit and evaluated
#' at each draw, and `sd` is the rest of the log prior brms kept. Without group-level
#' terms that rest must be zero, which checks that the priors were read as brms used
#' them. Normal and Student-t priors (and flat ones) on the coefficients can be
#' read; any other is refused, as are bounds on the coefficients and terms other than
#' linear and group-level ones (splines, monotonic effects). A block whose priors are
#' all flat has nothing to scale and is left out, with a message.
#'
#' @section Reading the sensitivity:
#' priorsense reads the two sensitivities together, against a threshold of 0.05.
#' When both are at or above it, the prior and the data disagree, and it reports a
#' "potential prior-likelihood conflict". When only the sensitivity to the prior is
#' above it, the prior rather than the data sets the estimate, and it reports a
#' "potential strong prior / weak likelihood"; every other row reads "-". In a model
#' with group-level terms, scaling the likelihood is like changing the number of
#' individuals within each group, not the number of groups, so a quantity that the
#' groups inform, such as their standard deviation, can show little likelihood
#' sensitivity without the likelihood being weak (Kallioinen et al. 2024, section
#' 2.5). The threshold marks a sensitivity that can be detected, not one that
#' matters. In simulations of a design with three species whose failures lay some
#' degrees apart, a sensitivity of the Weibull shape of 0.18 to 0.28 to the block
#' `loga` came with a bias in the shape of about two and a half posterior standard
#' deviations, while sensitivities up to 0.09 came with biases under a third of one;
#' no fit fell between 0.09 and 0.18. Nor is the size of a sensitivity the size of a
#' bias: in the same simulations sCTmax reached sensitivities of 0.06 to 0.12 while
#' its bias was at most 0.08 degrees. The other blocks never exceeded 0.03.
#'
#' Each scaling is an importance-sampling estimate, and the Pareto k of its weights
#' says whether they can be trusted: `k_prior` is the larger k of the two scalings of
#' the block and `k_likelihood` that of the likelihood's two. A k above priorsense's
#' threshold (min(1 - 1 / log10(S), 0.7) for S draws; 0.7 from 4,000 draws) warns,
#' and the sensitivities of that scaling are a sign at most. A block whose prior lies
#' far from the posterior is where the weights degrade, and that is the case the
#' check is for, so the warning is to be read rather than silenced.
#'
#' The test of a flagged block is a refit with that block's priors widened (through
#' `intercept_priors` and `slope_prior` of [build_hazard_prior()]), compared with
#' the first fit. The table `priors` shows which coefficients lie far out in their
#' prior or have a posterior almost as wide as it.
#'
#' @param fit A `brmsfit`, a model fitted with [fit_hazard()].
#' @param newdata A data frame of covariate values, one row per setting at which to
#'   check the reported parameters (and sCTmax, with `t_ref`), e.g. one row per
#'   species. It must hold every covariate in the model's parameter formulas, as in
#'   [derive_tdt_ttl()], and the quantities are computed with the random effects at
#'   zero, as there. `NULL` (the default) checks the coefficients only.
#' @param t_ref A single positive number, the reference exposure time for sCTmax, in
#'   the time unit of the data the model was fitted to. `NULL` (the default) leaves
#'   sCTmax out; it needs `newdata`.
#' @param p A single number strictly between 0 and 1, the fraction failing within
#'   `t_ref` that defines sCTmax.
#' @param levels A named numeric vector giving the level of each co-stressor, on its
#'   measurement scale, at which the parameters and sCTmax are checked; required
#'   with `newdata` when the model has co-stressors, as in [derive_tdt_ttl()], and
#'   ignored, with a warning, without `newdata`.
#' @return A list of two tibbles:
#'   * `sensitivity`, one row per block and quantity: `block`; `quantity`, a
#'     coefficient as brms names its draws (`b_loga_Intercept`) or, for each row of
#'     `newdata`, a reported parameter (`a`, `b`, `shape`, ...) or `sCTmax`; the
#'     columns of `newdata`, NA for a coefficient; `prior` and `likelihood`, the
#'     sensitivity to scaling the block and to scaling the likelihood; `diagnosis`,
#'     priorsense's reading of the two, in its own words; and `k_prior` and
#'     `k_likelihood`, the Pareto k of the scaling of the block and of the
#'     likelihood (see *Reading the sensitivity*).
#'   * `priors`, one row per population-level coefficient: `parameter`,
#'     `coefficient`, `prior`; `prior_centre` and `prior_scale` (the mean and
#'     standard deviation of a normal, the location and scale of a Student-t, NA for
#'     a flat prior); the posterior `median` and `sd`; `distance`, the median's
#'     distance from the prior's centre in units of its scale; and `sd_ratio`, the
#'     posterior standard deviation over the prior's scale.
#' @references
#' Kallioinen, N., Paananen, T., Burkner, P.-C. and Vehtari, A. (2024). Detecting
#' and diagnosing prior and likelihood sensitivity with power-scaling. *Statistics
#' and Computing* 34, 57.
#' @examples
#' \dontrun{
#' # a model fitted to a small ramping trial: 20 individuals warmed from 20 to 35
#' # degrees over 60 minutes and checked every minute, in two trials (see fit_hazard())
#' path   <- data.frame(path_id = "A", time = c(0, 60), temp = c(20, 35))
#' checks <- data.frame(path_id = "A", time = seq(0, 60, by = 1))
#' ind    <- data.frame(id = 1:20, path_id = "A", event = 1, mass = runif(20, 8, 12),
#'                      trial = rep(c("a", "b"), each = 10),
#'                      t_upper = sample(seq(10, 60, 10), 20, replace = TRUE))
#' long   <- build_hazard_data(ind, path, "interval", covariates = c("mass", "trial"),
#'                             checks = checks)
#' fit    <- fit_hazard(long, "exp", dpar_formulas = list(loga = ~ 1 + mass))
#'
#' chk <- check_prior(fit, newdata = data.frame(mass = c(8, 12)), t_ref = 60)
#' chk$sensitivity
#' chk$priors
#' }
#' @export
check_prior <- function(fit, newdata = NULL, t_ref = NULL, p = 0.5, levels = NULL) {
  .fit_settings(fit)
  if (!.has_priorsense())
    cli::cli_abort(c(
      "{.fn check_prior} needs the package {.pkg priorsense}.",
      "i" = "Install it with {.code install.packages(\"priorsense\")}."
    ), class = "rTFT_needs_priorsense")
  if (!is.numeric(p) || length(p) != 1L || is.na(p) || p <= 0 || p >= 1)
    cli::cli_abort(c(
      "{.arg p} must be a single fraction of individuals failing, strictly between 0 and 1.",
      "x" = "You supplied {.val {p}}."
    ), class = "rTFT_bad_p")
  if (!is.null(t_ref)) {
    if (is.null(newdata))
      cli::cli_abort(c(
        "Give {.arg newdata} with {.arg t_ref}: sCTmax is checked at each of its rows.",
        "i" = "For a model without covariates, {.code newdata = data.frame(row.names = 1)} is one row."
      ), class = "rTFT_needs_newdata")
    if (!is.numeric(t_ref) || length(t_ref) != 1L || !is.finite(t_ref) || t_ref <= 0)
      cli::cli_abort(c(
        "{.arg t_ref} must be a single positive exposure time, in the time unit of the data the model was fitted to.",
        "x" = "You supplied {.val {t_ref}}."
      ), class = "rTFT_bad_t_ref")
  }
  if (is.data.frame(newdata)) {
    taken <- intersect(names(newdata), c("block", "quantity", "prior", "likelihood",
                                         "diagnosis", "k_prior", "k_likelihood"))
    if (length(taken))
      cli::cli_abort(c(
        "{.arg newdata} has {cli::qty(length(taken))}{?a column/columns} named {.field {taken}}, {?a name/names} the table of sensitivities uses.",
        "i" = "Rename {cli::qty(length(taken))}{?it/them} in {.arg newdata}."
      ), class = "rTFT_newdata_names")
    if (anyNA(newdata))
      cli::cli_abort(c(
        "{.arg newdata} has missing values.",
        "i" = "Each row sets the values at which the parameters are checked, so every value is needed."
      ), class = "rTFT_missing_values")
  } else if (is.null(newdata) && !is.null(levels)) {
    cli::cli_warn(c(
      "!" = "Ignoring {.arg levels}: without {.arg newdata} only the coefficients are checked, which no level of a co-stressor enters.",
      "i" = "Give {.arg newdata} to check the parameters and sCTmax at those levels."
    ), class = "rTFT_unused_levels")
  }

  draws <- as.data.frame(brms::as_draws_df(fit))
  pb    <- .prior_blocks(fit, draws)
  qs    <- .check_quantities(fit, draws, pb$coefs, newdata, p, t_ref, levels)
  sens  <- .sensitivity_table(qs$values, pb$blocks, draws$lprior, .summed_log_lik(fit),
                              draws)

  at <- qs$index[sens$q, ]
  sensitivity <- tibble::tibble(block = sens$block, quantity = at$quantity)
  if (!is.null(newdata))
    sensitivity <- dplyr::bind_cols(sensitivity,
                                    tibble::as_tibble(newdata[at$row, , drop = FALSE]))
  sensitivity <- dplyr::bind_cols(
    sensitivity, sens[c("prior", "likelihood", "diagnosis", "k_prior", "k_likelihood")])
  list(sensitivity = sensitivity, priors = .prior_table(pb$coefs, draws))
}
