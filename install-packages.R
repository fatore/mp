packages = c(
                "roxygen2",
                "Rcpp",
                "RcppArmadillo",
                "proxy"
            )


# Create the personal library if it doesn't exist. Ignore a warning if the directory already exists.
dir.create(Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)

# Install listed packages
sapply(packages, install.packages, Sys.getenv("R_LIBS_USER"))

