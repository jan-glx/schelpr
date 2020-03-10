# scHelpR
random collection of helper functions that have yet to make it into appropriate packages


## Installation

```R
# install.packages("remotes") # lightway package that provides the implementation of "install_*" functions for devtools
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true") # to avoid warnings in Seurat cause the installation to fail
remotes::install_github("jan-glx/schelpr")
```
