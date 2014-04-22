width.f <- max(nchar(rownames(f)))
width.c <- max(nchar(rownames(c)))


fs <- split(f[, c('x', 'y')], formatC(rownames(f), width = width.f, format = 'd', flag = '0'))
cs <- split(c[, c('x', 'y')], formatC(rownames(c), width = width.c, format = 'd', flag = '0'))


## N by M matrix of distances
t <- outer(fs, cs, FUN = Vectorize(function(a, b) dist(rbind(a, b))))

## Objective function: all facilities, then distance matrix columnwise
## (so that variables for each customer go together)
obj <- c(f$s, as.vector(t))

suppressWarnings(suppressMessages(require ('Matrix')))

## Matrix of constraints
## Assign customer to facility only if it is open: y_{w, c} <= x_{w}

## lhs1 <- cbind(
##     do.call(rbind, lapply(1:params$M,
##                           function(x) diag(-1, nrow = params$N))),
##     diag(1, nrow = params$N * params$M))

m1 <- Matrix(0, nrow = params$N * params$M, ncol = params$N)

m1[c(sapply(1:params$N,
            function(w) {
                (w - 1) * (params$N * params$M + 1) +
                    seq(1, params$N * params$M, params$N)
            }))] <- -1


lhs1 <- cBind(
    m1, 
    .sparseDiagonal(x = 1, n = params$N * params$M))

rhs1 <- rep(0, params$N * params$M)
rel1 <- rep('<=', params$N * params$M)

## One facility for a customer: Sum_{w} y_{w, c} = 1

## lhs2 <- cbind(
##     matrix(0, nrow = params$M, ncol = params$N),
##     do.call(rbind, lapply(1:params$M,
##                           function(w) {
##                               c(rep(0, params$N * (w - 1)),
##                                 rep(1, params$N),
##                                 rep(0, params$N * (params$M - w)))
##                           })))


bandM <- Matrix(0, nrow = params$M, ncol = params$N * params$M)

bandM[c(sapply(1:params$M,
               function(c) {
                   (c - 1) * (params$N * params$M + 1) +
                       seq(1, params$N * params$M, params$M)
               }))] <- 1

lhs2 <- cBind(
    Matrix(0, nrow = params$M, ncol = params$N),
    bandM)


rhs2 <- rep(1, params$M)

rel2 <- rep('==', params$M)

                                        #rel2 <- rep('=', params$M)

## Capacity constraint 

m3 <- Matrix(0, nrow = params$N, ncol = params$N * params$M)

for (i in 1:params$M) {
    sq <- params$N * params$N
    m3[((i - 1) * sq) + seq(1, sq, params$N + 1)] <- c$d[[i]]
}

## lhs3 <- cBind(
##     Matrix(0, nrow = params$N, ncol = params$N),
##     do.call(cBind, lapply(1:params$M,
##                           function(w) Diagonal(params$N) * c$d[[w]]
##                           )))

lhs3 <- cBind(
    Matrix(0, nrow = params$N, ncol = params$N),
    m3)

rhs3 <- f$cap
rel3 <- rep('<=', params$N)

                                        #suppressWarnings(suppressMessages(require('Rglpk')))
###suppressWarnings(suppressMessages(require('gurobi')))

problem <- list(sense = "min")

problem$A <- rBind(lhs1, lhs2, lhs3)

problem$c <- obj

dir <- c(rel1, rel2, rel3)

rhs <- c(rhs1, rhs2, rhs3)

## lower and upper constraints
blc <- rep(-Inf, length(dir))
idx <- which(dir %in% c(">=", "=="))
blc[idx] <- rhs[idx]

buc <- rep(Inf, length(dir))
idx <- which(dir %in% c("<=", "=="))
buc[idx] <- rhs[idx]

problem$bc <- rbind(blc, buc)

## variable constraints ([0, 1])
problem$bx <- rbind(rep(0, params$N * (1 + params$M)),
                    rep(1, params$N * (1 + params$M)))

problem$intsub <- 1:(params$N * (1 + params$M))
    
print(params$N)
print(params$M)
print(dim(m3))
print(dim(lhs3))
print(length(rhs3))

res <- mosek(problem)

## model  <- list()

## model$A      <- mat
## model$obj    <- obj

## model$sense  <- dir

## model$modelsense <- 'min'

## model$rhs    <- rhs
## model$vtype  <- "B"

##                                     #gparams <- list(Presolve=2, TimeLimit=100.0)
## gparams <- list(Presolve=2, TuneOutput = 0,
##                 OutputFlag = 0, TimeLimit = 200.,
##                 ConcurrentMIP =  4)

## res <- gurobi(model, gparams)


## res <- gurobi(obj, mat, dir, rhs, max = max, types = types)

                                        #cat(paste0(round(res$optimum), ' ', abs(1 - res$status), '\n'))
cw <- sapply(
    split(res$x[-(1:params$N)],
          sapply(1:params$M, function(c) rep(c, params$N))),
    function(c) which.max(c))

list(res$objval, cw)
