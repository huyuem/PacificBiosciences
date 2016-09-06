library(rmarkdown) 
argv  = commandArgs (TRUE)
rmarkdown::render(argv[1])
