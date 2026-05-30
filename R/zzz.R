# R/zzz.R

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
