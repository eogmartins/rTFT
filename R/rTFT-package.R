# rTFT-package.R - package-level documentation

#' @keywords internal
#' @importFrom utils globalVariables
"_PACKAGE"

# Column names referenced inside dplyr non-standard evaluation, declared so R CMD
# check does not read them as undefined globals. Only the names the reviewed
# functions actually use are listed; more are added as further functions return.
# build_hazard_data(): `id` (select/relocate), `tau` and `exposure` (mutate).
globalVariables(c("exposure", "id", "tau"))
