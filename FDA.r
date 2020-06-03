library(fda)

scriptsDir <- system.file('scripts', package='fda')
Rscripts <- dir(scriptsDir, full.names=TRUE, pattern='R$')
fdarm <- grep('fdarm', Rscripts, value=TRUE)
chapters <- length(fdarm)
# NOTE:  If R fails in any of these scripts,
# this for loop will not end normally,
# and the abnormal termination will be displayed:
for(ch in 1:chapters){
  cat('Running', fdarm[ch], '\n')
  invisible(source(fdarm[ch]))
}
find.package('fda')
help(fda)
fdarm
