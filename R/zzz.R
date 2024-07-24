.onLoad <- function(libname, pkgname) {
    ## Attempt to set the correct python
    tryCatch({
          reticulate::py_set_seed(1234)
          reticulate::use_python("/usr/bin/python3")
        },
        error=function(cond) {
          message("Reticulate failed to detect a usable python environment. Defaulting to Docker path.")
          message("If this error persists, please set the correct python path via `reticulate::use_python` before running scattch-mapping.")
        },
        finally={
          ## Run once more to produce an error if Docker path also fails.
          reticulate::py_set_seed(1234)
          if (py_config()$python != "/usr/bin/python3") {
            message("Reticulate failed to detect a usable python environment. Please set the correct python path via `reticulate::use_python(\"/usr/bin/python3\")`.")
          }
        }
    )
}