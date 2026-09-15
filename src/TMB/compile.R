tmb_name <- "langevinSSM_TMBExports"
tmb_flags <- commandArgs(trailingOnly = TRUE)

if(file.exists(paste0(tmb_name, ".cpp"))) {
  if(length(tmb_flags) == 0) tmb_flags <- ""

  # Strip literal quotes passed by Windows cmd.exe
  tmb_flags <- gsub("['\"]", "", tmb_flags)
  tmb_flags <- paste(tmb_flags, collapse = " ")

  tmb_flags <- paste(tmb_flags, "-g0")

  # Force append suppression flag to CXX17FLAGS on Windows only.
  # This ensures it evaluates AFTER R's default -Wall, successfully silencing the GCC bug.
  if (.Platform$OS.type == "windows") {
    cxx17 <- Sys.getenv("CXX17FLAGS")
    Sys.setenv(CXX17FLAGS = paste(cxx17, "-Wno-array-bounds"))
  }

  options(tmb.ad.framework = "TMBad")

  TMB::compile(file = paste0(tmb_name, ".cpp"),
               PKG_CXXFLAGS = tmb_flags,
               safebounds = FALSE, safeunload = FALSE)

  file.copy(from = paste0(tmb_name, .Platform$dynlib.ext),
            to = "..", overwrite = TRUE)
}
