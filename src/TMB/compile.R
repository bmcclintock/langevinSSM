tmb_name <- "langevinSSM_TMBExports"
args <- commandArgs(trailingOnly = TRUE)
tmb_flags <- if (length(args) >= 1) args[1] else ""
# Value of $(SHLIB_OPENMP_CXXFLAGS) as expanded by make: empty when R has no OpenMP flag for C++ (e.g. CRAN's
# macOS build of R).  OpenMP can then still be enabled by the user's Makevars (e.g. 'CPPFLAGS += -Xclang -fopenmp'
# and 'LDFLAGS += -lomp' on macOS), which TMB::compile() passes on to the compiler as well.
omp_flags <- if (length(args) >= 2) trimws(args[2]) else ""
use_openmp <- nzchar(omp_flags)

# The user's Makevars, looked up as R CMD INSTALL does (TMB::compile() itself only looks at R_MAKEVARS_USER and
# ~/.R/Makevars), so that the TMB library is compiled with the same flags as the rest of the package.
if (is.na(Sys.getenv("R_MAKEVARS_USER", NA))) {
  cand <- if (.Platform$OS.type == "windows") {
    c(if (Sys.getenv("R_ARCH") == "/x64") c("~/.R/Makevars.ucrt", "~/.R/Makevars.win64"), "~/.R/Makevars.win", "~/.R/Makevars")
  } else {
    c(paste0("~/.R/Makevars-", Sys.getenv("R_PLATFORM")), "~/.R/Makevars")
  }
  cand <- path.expand(cand)
  cand <- cand[file.exists(cand)]
  if (length(cand)) Sys.setenv(R_MAKEVARS_USER = cand[1])
}

if(file.exists(paste0(tmb_name, ".cpp"))) {
  options(tmb.ad.framework = "TMBad")
  # TMB::compile(openmp = TRUE) writes PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS) to its Makevars, but the
  # PKG_CXXFLAGS given here is written after it and takes precedence, so the OpenMP flag is added here too.
  # PKG_LIBS gets $(SHLIB_OPENMP_CXXFLAGS) from TMB::compile(openmp = TRUE).
  TMB::compile(file = paste0(tmb_name, ".cpp"),
               PKG_CXXFLAGS = paste(tmb_flags, if (use_openmp) "$(SHLIB_OPENMP_CXXFLAGS)"),
               openmp = use_openmp,
               safebounds = FALSE, safeunload = FALSE)
  file.copy(from = paste0(tmb_name, .Platform$dynlib.ext), to = "..", overwrite = TRUE)
}
