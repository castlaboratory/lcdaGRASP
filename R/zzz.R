# R/zzz.R

# Column names referenced inside ggplot2::aes() in R/plot_leaders.R. Declared
# here so R CMD check does not flag them as undefined globals (ggplot2 is a
# Suggested dependency, so we avoid pulling in rlang just for `.data`).
utils::globalVariables(c("x", "y", "xend", "yend", "intra", "community",
                         "label", "is_leader", "scope"))

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    cli::format_inline(
      "{.pkg lcdaGRASP} {utils::packageVersion('lcdaGRASP')} - ",
      "GRASP / Reactive-GRASP for joint community-leader detection."
    ), "\n",
    cli::format_inline(
      "Reference: Ospina et al. (2026), preprint. ",
      "See {.code vignette('lcdaGRASP-intro')}."
    )
  )
}
