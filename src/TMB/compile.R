tmb_name <- "langevinSSM_TMBExports"
tmb_flags <- commandArgs(trailingOnly = TRUE)

if(file.exists(paste0(tmb_name, ".cpp"))) {
  if(length(tmb_flags) == 0) tmb_flags <- ""

  tmb_flags <- paste(tmb_flags, "-g0")

  # ON WINDOWS ONLY: Force the compilation flag to -O2
  # This stops GCC 14's hyper-aggressive array boundary warnings
  # while keeping SIMD vectorization fully functional.
  if (.Platform$OS.type == "windows") {
    tmb_flags <- paste(tmb_flags, "-O2")
  }

  options(tmb.ad.framework = "TMBad")

  TMB::compile(file = paste0(tmb_name, ".cpp"),
               PKG_CXXFLAGS = tmb_flags,
               safebounds = FALSE, safeunload = FALSE)

  file.copy(from = paste0(tmb_name, .Platform$dynlib.ext),
            to = "..", overwrite = TRUE)
}
