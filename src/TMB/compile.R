tmb_name <- "langevinSSM_TMBExports"
tmb_flags <- paste(commandArgs(trailingOnly = TRUE), collapse = " ")

if(file.exists(paste0(tmb_name, ".cpp"))) {

  options(tmb.ad.framework = "TMBad")

  TMB::compile(file = paste0(tmb_name, ".cpp"),
               PKG_CXXFLAGS = tmb_flags,
               safebounds = FALSE, safeunload = FALSE)

  file.copy(from = paste0(tmb_name, .Platform$dynlib.ext),
            to = "..", overwrite = TRUE)
}
