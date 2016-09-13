library(rmarkdown) 
argv  = commandArgs (TRUE)
# rmarkdown::run(argv[1])
rmarkdown::render(argv[1])

