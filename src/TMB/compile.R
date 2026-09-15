tmb_name <- "langevinSSM_TMBExports"
tmb_flags <- commandArgs(trailingOnly = TRUE)

if(file.exists(paste0(tmb_name, ".cpp"))) {
  if(length(tmb_flags) == 0) tmb_flags <- ""

  # Strip literal quotes passed by Windows cmd.exe
  tmb_flags <- gsub("['\"]", "", tmb_flags)
  tmb_flags <- paste(tmb_flags, collapse = " ")

  # Append optimization and vectorization suppression flags
  tmb_flags <- paste(tmb_flags, "-g0 -DEIGEN_DONT_VECTORIZE")

  options(tmb.ad.framework = "TMBad")

  TMB::compile(file = paste0(tmb_name, ".cpp"),
               PKG_CXXFLAGS = tmb_flags,
               safebounds = FALSE, safeunload = FALSE)

  file.copy(from = paste0(tmb_name, .Platform$dynlib.ext),
            to = "..", overwrite = TRUE)
}
