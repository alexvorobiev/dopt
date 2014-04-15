args <- commandArgs(trailingOnly = TRUE)

file <- substr(args, 7, nchar(args[1]))

                                        #file <- 'data/fl_25_2'
#file <- 'data/fl_test'

params <- read.table(file, nrows = 1)
names(params) <- c('N', 'M')

## N Facilities
f <- read.table(file, skip = 1, nrows = params$N)
names(f) <- c('s', 'cap', 'x', 'y')

## M Customers
c <- read.table(file, skip = 1 + params$N, nrows = params$M)
names(c) <- c('d', 'x', 'y')

## Coordinates as lists of pairs x, y
## fs <- split(f[, c('x', 'y')], rownames(f))
## cs <- split(c[, c('x', 'y')], rownames(c))

## fs <- split(f, row(f))[, c('x', 'y')]
## cs <- split(c, row(c))[, c('x', 'y')]

## fs <- t(apply(f, 1, function(x) x[c(3, 4)]))
## cs <- t(apply(c, 1, function(x) x[c(3, 4)]))

## There is got to be a better way...
width.f <- max(nchar(rownames(f)))
width.c <- max(nchar(rownames(c)))

fs <- split(f[, c('x', 'y')], formatC(rownames(f), width = width.f, format = 'd', flag = '0'))
cs <- split(c[, c('x', 'y')], formatC(rownames(c), width = width.c, format = 'd', flag = '0'))


## N by M matrix of distances
t <- outer(fs, cs, FUN = Vectorize(function(a, b) dist(rbind(a, b))))

## Objective function: all facilities, then distance matrix columnwise
## (so that variables for each customer go together)
obj <- c(f$s, as.vector(t))

## Matrix of constraints
## Assign customer to facility only if it is open: y_{w, c} <= x_{w}

lhs1 <- cbind(
    do.call(rbind, lapply(1:params$M,
                          function(x) diag(-1, nrow = params$N))),
    diag(1, nrow = params$N * params$M))

rhs1 <- rep(0, params$N * params$M)
rel1 <- rep('<=', params$N * params$M)

## One facility for a customer: Sum_{w} y_{w, c} = 1

lhs2 <- cbind(
    matrix(0, nrow = params$M, ncol = params$N),
    do.call(rbind, lapply(1:params$M,
                          function(w) {
                              c(rep(0, params$N * (w - 1)),
                                rep(1, params$N),
                                rep(0, params$N * (params$M - w)))
                          })))

rhs2 <- rep(1, params$M)
rel2 <- rep('==', params$M)

## Capacity constraint 
lhs3 <- cbind(
    matrix(0, nrow = params$N, ncol = params$N),
    do.call(cbind, lapply(1:params$M,
                          function(w) diag(params$N) * c$d[[w]]
                          )))
rhs3 <- f$cap
rel3 <- rep('<=', params$N)

suppressWarnings(suppressMessages(require('Rglpk')))

mat <- rbind(lhs1, lhs2, lhs3)
dir <- c(rel1, rel2, rel3)
rhs <- c(rhs1, rhs2, rhs3)

max <- FALSE

types <- rep('B', params$N * (1 + params$M))

res <- Rglpk_solve_LP(obj, mat, dir, rhs, max = max, types = types)

cat(paste0(round(res$optimum), ' ', abs(1 - res$status), '\n'))

cw <- sapply(
    split(res$solution[-(1:params$N)],
          sapply(1:params$M, function(c) rep(c, params$N))),
    function(c) which.max(c))

cat(paste0(cw - 1, collapse = ' '))
cat('\n')


