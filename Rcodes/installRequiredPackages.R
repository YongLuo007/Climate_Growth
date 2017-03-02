rm(list = ls())
#checkPackageVersion on your default
install.packages("httr")
library(httr)
workPath <- "~/GitHub/Climate_Growth"
packageTable <- data.frame(rbind(c("SpaDES", "1.3.0"),
                                 c("nlme", "3.1-128"),
                                 c("raster", "2.5-2"),
                                 c("ade4", "1.7-4"),
                                 c("data.table", "1.10.0"),
                                 c("ggplot2", "2.1.0"),
                                 c("dplyr", "0.4.3"),
                                 c("MuMIn", "1.15.1"),
                                 c("gridExtra", "2.2.0"),
                                 c("colorRamps", "2.2"),
                                 c("maptools", "0.8-41"),
                                 c("rgeos", "0.3-20")),
                           stringsAsFactors = FALSE)
names(packageTable) <- c("PackageName", "requiedVersion")

libPath <- file.path(workPath, "RequiredRPackages")
if(!dir.exists(libPath)){
  dir.create(libPath)
}
# Download all the packages and install them


for (i in 1:nrow(packageTable)){
  packageurl <- paste("https://cran.r-project.org/src/contrib/Archive/",
                      packageTable$PackageName[i], "/", packageTable$PackageName[i],
                      "_", packageTable$requiedVersion[i], ".tar.gz", sep = "")
  if(http_error(packageurl)){
    packageurl <- paste("https://cran.r-project.org/src/contrib/",
                        packageTable$PackageName[i],
                        "_", packageTable$requiedVersion[i], ".tar.gz", sep = "")
  }
  if(!dir.exists(file.path(libPath, packageTable$PackageName[i]))){
    install.packages(packageurl, repos=NULL, type="source", lib = libPath, 
                     INSTALL_opts = c("--no-lock"))
  }
}
