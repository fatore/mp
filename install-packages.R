## Create the personal library if it doesn't exist. Ignore a warning if the directory already exists.
dir.create(Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)

install.packages(c("roxygen2", "RcppArmadillo", "proxy"), Sys.getenv("R_LIBS_USER"), repos = "http://cran.case.edu" )

## Install a package that you have copied to the remote system.
## install.packages("file_name.tar.gz", Sys.getenv("R_LIBS_USER")

