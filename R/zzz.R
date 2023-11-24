



.onLoad <- function(libname, pkgname) {

  invisible(suppressPackageStartupMessages(library("ggplot2")))
  invisible(suppressPackageStartupMessages(library("IOBR")))

  invisible(suppressPackageStartupMessages(
    sapply(c("crayon", "randomForest", "randomForestSRC", "survminer", "ggplot2"),
           requireNamespace, quietly = TRUE)
    ))
}



##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0("==========================================================================\n",
                " ", pkgname, " v", pkgVersion, "  ",

                "  For help: https://github.com/LiaoWJLab/uaiscore/issues", "\n\n")

  citation <-paste0(" If you use ", pkgname, " in published research, please cite:\n",
                    " DQ Zeng, XT Huang, …, WJ Liao*.\n",
                    " UAIscore enables precise stratification of adjuvant immunotherapy in urothelial carcinoma \n",
                    " Under Review. \n",
                    # " XXXX, 2020", "\n",
                    # " DOI: 10.3389/bioRxiv.2021.687975\n",
                    # " PMID:  ","\n",
                    "==========================================================================")

  packageStartupMessage(paste0(msg, citation))
}



