
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Attaching peaksat version ",
                        packageDescription("peaksat")$Version, ".")
  #When adding new options here, also add them to the "names" setMethod below
  PS_OPTIONS <<- new("PS_OPTIONS")
  PS_STATS <<- valid_stats
}

setClass("PS_OPTIONS", representation = list(
  is_valid = "logical"
  ))

setMethod("names", "PS_OPTIONS",
          function(x)
          {
            c(
              "PS_MACS2_PATH",
              "PS_OUTDIR",
              "PS_JOB_IDS"
            )
          })


setMethod("$", "PS_OPTIONS",
          function(x, name)
          {
            getOption(name)
          })

setReplaceMethod("$", "PS_OPTIONS",
                 function(x, name, value)
                 {
                   value = list(value)
                   names(value) = name
                   do.call("options", value)
                   x
                 })



setMethod("show", "PS_OPTIONS",
          function(object)
          {
            message("Use the $ accessor (i.e. PS_OPTIONS$PS_MACS2_PATH) to get/set PS relevant options.")
            message("Use names(PS_OPTIONS) to view all options.")
          })


