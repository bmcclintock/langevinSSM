tmb_name <- "langevinSSM_TMBExports"
tmb_flags <- commandArgs(trailingOnly = TRUE)

if(file.exists(paste0(tmb_name, ".cpp"))) {
  if(length(tmb_flags) == 0) tmb_flags <- ""

  # Strip literal quotes passed by Windows cmd.exe
  tmb_flags <- gsub("['\"]", "", tmb_flags)
  tmb_flags <- paste(tmb_flags, collapse = " ")

  tmb_flags <- paste(tmb_flags, "-g0")

  options(tmb.ad.framework = "TMBad")

  # TMB::compile uses R CMD SHLIB, which reads Makeconf (-Wall) *after* PKG_CXXFLAGS.
  # To prevent -Wall from re-enabling the array-bounds warning, we must append our
  # suppression directly to CXX17FLAGS using a local Makevars.win file.
  if (.Platform$OS.type == "windows") {
    writeLines("CXX17FLAGS += -Wno-array-bounds\n", "Makevars.win")
  }

  TMB::compile(file = paste0(tmb_name, ".cpp"),
               PKG_CXXFLAGS = tmb_flags,
               safebounds = FALSE, safeunload = FALSE)

  # Cleanup
  if (.Platform$OS.type == "windows" && file.exists("Makevars.win")) {
    file.remove("Makevars.win")
  }

  file.copy(from = paste0(tmb_name, .Platform$dynlib.ext),
            to = "..", overwrite = TRUE)
}
