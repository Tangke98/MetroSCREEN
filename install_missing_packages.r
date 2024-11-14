if (!requireNamespace("kpcalg", quietly = TRUE)) {
  install.packages("kpcalg", repos = "http://cran.us.r-project.org")
}
if (!requireNamespace("nichenetr", quietly = TRUE)) {
  devtools::install_github("saeyslab/nichenetr")  # nichenetr 是 GitHub 包
}
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
devtools::install(".", upgrade="never")