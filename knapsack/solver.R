args <- commandArgs(trailingOnly = TRUE)

data <- read.table(substr(args, 7, nchar(args[1])))

data <- read.table('data/ks_100_0')

n <- data[1, 1]
K <- data[1, 2]

v <- data[-1, 1]
w <- data[-1, 2]

suppressWarnings(suppressMessages (require ('Rsymphony')))
suppressWarnings(suppressMessages (require ('Rglpk')))


knapsack <- function (w, p, cap) 
{
    n <- length(w)
    x <- logical(n)
    F <- matrix(0, nrow = cap + 1, ncol = n)
    G <- matrix(0, nrow = cap + 1, ncol = 1)
    for (k in 1:n) {
        F[, k] <- G
        H <- c(numeric(w[k]), G[1:(cap + 1 - w[k]), 1] + p[k])
        G <- pmax(G, H)
    }
    fmax <- G[cap + 1, 1]
    f <- fmax
    j <- cap + 1
    for (k in n:1) {
        if (F[j, k] < f) {
            x[k] <- TRUE
            j <- j - w[k]
            f <- F[j, k]
        }
    }
    inds <- which(x)
    wght <- sum(w[inds])
    prof <- sum(p[inds])
    return(list(capacity = wght, profit = prof, indices = inds))
}

res <- knapsack(w = w, p = v, cap = K)

mat <- matrix(w, nrow = 1)

dir <- c('<=')

rhs <- c(K)

max <- TRUE

types <- rep('B', n)

res <-Rsymphony_solve_LP(v, mat, dir, rhs, max = max, types = types)

res <-Rglpk_solve_LP(v, mat, dir, rhs, max = max, types = types,
                     control = list("verbose" = TRUE))
    
cat(paste0(res$profit, ' 0\n'))
cat(paste0(ifelse(1:n %in% res$indices, '1', '0')))
cat('\n')


