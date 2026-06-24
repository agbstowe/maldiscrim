
<!-- README.md is generated from README.Rmd. Please edit that file -->

# maldiscrim

<!-- badges: start -->
<!-- badges: end -->

`maldiscrim` is an R package designed for functional data
discrimination, specifically using fPLS-DA methods and deep learning
approach to classify and analyze spectroscopic data.

## Installation

You can install the development version of maldiscrim with:

``` r
# Option 1: Using devtools (Standard)
install.packages("devtools")
devtools::install_github("agbstowe/maldiscrim")

# Option 2: Using pak (Fast/Modern)
# install.packages("pak")
pak::pak("agbstowe/maldiscrim")
```

## Create Python environnement for deep learning approach (FNN method)

To use the deep learning features (FNN method) of this package, a
dedicated Python virtual environment named `maldiscrim-env` must be
configured. This background environment isolates mandatory Python
dependencies—such as TensorFlow and Keras—preventing version conflicts
with other Python projects on your computer.

**This setup only needs to be performed once**. You can automatically
build, update, audit, or manually enforce the use of this environment
using the following instructions:

``` r
library(maldiscrim)
```

``` r
## Initialize and build the mandatory Python environment (Required once)
maldiscrim_install_python()

## Force-reinstall or update the environment
# If your environment becomes unstable, corrupted, or if you apply systemic modifications that disrupt the Python layer, you can pass 'force = TRUE'. This will cleanly wipe the existing 'maldiscrim-env' and reconstruct it from scratch.
maldiscrim_install_python(force = TRUE)

## List available virtual environments on your machine
# Use this reticulate helper function to audit your local system. It allows you to 
# visually confirm if 'maldiscrim-env' is correctly registered among your virtual environments.
reticulate::virtualenv_list()

## Explicitly select and lock the required environment before execution
# In case of environment ambiguity, multi-version conflicts, or to guarantee that your code executes strictly within the intended sandbox, call this function to force the current R session to bind directly to 'maldiscrim-env' before starting any training.
reticulate::use_virtualenv("maldiscrim-env", required = TRUE)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.
