# List of packages to check and install if necessary
packages <- c("nlme", 
              # "rms",  # Uncomment if you decide to use this library
              "phytools", "geiger", "caper", "tidyverse", 
              "cowplot", "ggrepel", "ggsci", "patchwork", 
              "poolr", "parallel", "MASS", "filelock", "R.utils")

# Function to check and install missing packages
check_and_install <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Apply the function to each package
sapply(packages, check_and_install)

# Message indicating the process is complete
cat("All specified packages are installed and loaded.\n")
