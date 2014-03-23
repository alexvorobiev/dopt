suppressWarnings(suppressMessages(require('RBGL')))

args <- commandArgs(trailingOnly = TRUE)

file <- substr(args, 7, nchar(args[1]))
#file <- 'data/gc_50_3'

data <- read.table(file)

nodes <- as.character(0:(data[1, 1] - 1))

g <- new("graphNEL", nodes = nodes)

g <- addEdge(as.character(data[2:nrow(data), 1]),
            as.character(data[2:nrow(data), 2]), g)


res <- sequential.vertex.coloring(g)

cat(paste0(res$`no. of colors needed`, ' 0', '\n'))
cat(res$`colors of nodes`)
cat('\n')
