args <- commandArgs(trailingOnly = TRUE)

data <- read.table(substr(args, 7, nchar(args[1])))

n <- data[1, 1]
K <- data[1, 2]
v <- data[-1, 1]
w <- data[-1, 2]

## suppressWarnings(suppressMessages (require ('adagio')))
## suppressWarnings(suppressMessages (require ('Rsymphony')))
suppressWarnings(suppressMessages (require ('Rglpk')))

## knapsack function
##res <- knapsack(w = w, p = v, cap = K)
##cat(paste0(res$profit, ' 0\n'))
##cat(paste0(ifelse(1:n %in% res$indices, '1', '0')))

mat <- matrix(w, nrow = 1)
dir <- c('<=')
rhs <- c(K)
max <- TRUE
types <- rep('B', n)

#res <-Rsymphony_solve_LP(v, mat, dir, rhs, max = max, types = types)
#cat(paste0(res$objval, ' ', abs(1 - res$status), '\n'))

res <-Rglpk_solve_LP(v, mat, dir, rhs, max = max, types = types)

cat(paste0(res$optimum, ' ', abs(1 - res$status), '\n'))
cat(paste0(res$solution, collapse = ' '))
cat('\n')


