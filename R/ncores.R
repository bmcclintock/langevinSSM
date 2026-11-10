# --- OpenMP / multi-threading helpers ---
#
# langevinSSM can use several threads (OpenMP) in fitLangevin(), predLangevin(), simLangevin() and cAIC() via their
# 'ncores' argument.  The default is getOption("langevinSSM.ncores", 1), i.e. a single thread unless the user
# asks for more.  The package has two compiled parts, which are checked separately because they can, in principle,
# be built with different flags:
#   * the TMB library 'langevinSSM_TMBExports' (likelihood: fitLangevin, predLangevin, cAIC, simLangevin.fitLangevin)
#   * the Rcpp library 'langevinSSM' (track simulator: simLangevin)
# Whether a part has OpenMP is read from the compiled code itself (was _OPENMP defined when it was compiled?), so the
# answer is correct however OpenMP was enabled (R's $(SHLIB_OPENMP_CXXFLAGS), or e.g. 'CPPFLAGS += -Xclang -fopenmp'
# in ~/.R/Makevars on macOS).  A part without OpenMP runs serially and a request for more threads gives a message.

.tmb_dll <- "langevinSSM_TMBExports"

# Was the TMB library (likelihood) compiled with OpenMP?
.tmb_has_openmp <- function() {
  fw <- tryCatch(.Call("getFramework", PACKAGE = .tmb_dll), error = function(e) NULL)
  isTRUE(attr(fw, "openmp"))
}

# Was the Rcpp library (simulator) compiled with OpenMP?
.rcpp_has_openmp <- function() lssm_omp_num_procs() > 0L

# Validate 'ncores' and return the number of threads that will actually be used (an integer >= 1).
# 'has_openmp': does the part of the package that will use the threads have OpenMP support?
.check_ncores <- function(ncores, has_openmp = .tmb_has_openmp()) {
  if (is.null(ncores)) ncores <- 1L
  if (length(ncores) != 1 || !is.numeric(ncores) || is.na(ncores) || ncores < 1 || ncores != round(ncores))
    stop("'ncores' must be a single positive integer.")
  ncores <- as.integer(ncores)
  if (ncores == 1L) return(1L)
  if (!has_openmp) {
    if (!isTRUE(getOption("langevinSSM.ncores.warned"))) {
      message("langevinSSM was compiled without OpenMP support; 'ncores' is ignored and a single thread is used. ",
              "See ?has_omp.")
      options(langevinSSM.ncores.warned = TRUE)
    }
    return(1L)
  }
  nprocs <- lssm_omp_num_procs()
  if (nprocs > 0L && ncores > nprocs) {
    warning("'ncores' (", ncores, ") exceeds the number of available processors (", nprocs, "); using ", nprocs, ".")
    ncores <- nprocs
  }
  ncores
}

# Set the number of threads of the TMB library and return the previous setting (to be restored with
# .restore_tmb_threads(), typically via on.exit()).  Must be called before TMB::MakeADFun: the objective is taped in
# 'ncores' parts (one per thread) when the object is created, and evaluated with that many threads.  (An object taped
# in several parts remains valid when the number of threads is changed afterwards; it is then evaluated with the new
# number of threads.)  TMB's automatic parallelisation ('autopar') is not used, so the split is the one defined by
# the parallel_accumulator in the template.  Returns NULL if the TMB library has no OpenMP support.
.set_tmb_threads <- function(ncores) {
  if (!.tmb_has_openmp()) return(NULL)
  old <- TMB::openmp(DLL = .tmb_dll)
  TMB::openmp(n = ncores, autopar = FALSE, DLL = .tmb_dll)
  old
}

.restore_tmb_threads <- function(old) {
  if (!is.null(old)) TMB::openmp(n = as.integer(old), autopar = isTRUE(attr(old, "autopar")), DLL = .tmb_dll)
  invisible(NULL)
}

#' Is OpenMP multi-threading available in langevinSSM?
#'
#' @description
#' Checks whether the installed copy of \pkg{langevinSSM} was compiled with OpenMP support, i.e. whether the
#' \code{ncores} argument of \code{\link{fitLangevin}}, \code{\link{predLangevin}}, \code{\link{simLangevin}} and
#' \code{\link{cAIC}} will actually use several threads.
#'
#' @param verbose Logical; print a short report? Default: \code{FALSE}.
#'
#' @details
#' The package has two compiled parts: the TMB library that evaluates the likelihood (used by \code{fitLangevin},
#' \code{predLangevin}, \code{cAIC} and \code{simLangevin} for a fitted model) and the Rcpp library that contains
#' the track simulator (\code{simLangevin}). For each part, the answer is read from the compiled code itself (was
#' OpenMP enabled when it was compiled?), so it is correct however OpenMP was enabled: through R's
#' \code{SHLIB_OPENMP_CXXFLAGS} (Linux, Windows) or, on macOS, through the compiler flags in \file{~/.R/Makevars},
#' e.g. (see \url{https://mac.r-project.org/openmp/})
#' \preformatted{CPPFLAGS += -Xclang -fopenmp
#' LDFLAGS += -lomp}
#' The package must be installed from source after setting these flags (CRAN's macOS binary is built without
#' OpenMP).
#'
#' The number of threads is set per call with \code{ncores}, or for the whole session with
#' \code{options(langevinSSM.ncores = n)}. It is never larger than the number of processors (\code{nprocs}).
#' \code{\link{residuals.fitLangevin}} parallelises over tracks with separate R processes instead (its own
#' \code{ncores} argument), which does not need OpenMP.
#'
#' @return \code{TRUE} if both compiled parts have OpenMP support, \code{FALSE} otherwise, with attributes
#' \item{likelihood}{Logical; OpenMP support of the TMB library (likelihood).}
#' \item{simulator}{Logical; OpenMP support of the track simulator.}
#' \item{nprocs}{Number of processors available to OpenMP (\code{NA} without OpenMP).}
#' \item{ncores}{Default number of threads, \code{getOption("langevinSSM.ncores", 1)}.}
#'
#' @examples
#' has_omp()
#' has_omp(verbose = TRUE)
#'
#' @export
has_omp <- function(verbose = FALSE) {
  lik <- .tmb_has_openmp()
  sim <- .rcpp_has_openmp()
  np <- lssm_omp_num_procs()
  out <- lik && sim
  attr(out, "likelihood") <- lik
  attr(out, "simulator") <- sim
  attr(out, "nprocs") <- if (np > 0L) np else NA_integer_
  attr(out, "ncores") <- getOption("langevinSSM.ncores", 1L)
  if (verbose) {
    yn <- function(x) if (x) "yes" else "no"
    message("langevinSSM OpenMP support:\n",
            "  likelihood (fitLangevin, predLangevin, cAIC): ", yn(lik), "\n",
            "  simulator (simLangevin):                      ", yn(sim), "\n",
            "  processors available:                         ", if (np > 0L) np else "NA", "\n",
            "  default ncores (option langevinSSM.ncores):   ", attr(out, "ncores"),
            if (!out) "\n  To enable OpenMP, see ?has_omp." else "")
  }
  out
}
