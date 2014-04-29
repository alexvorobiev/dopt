args <- commandArgs(trailingOnly = TRUE)

file <- substr(args, 7, nchar(args[1]))


##file <- 'data/fl_test'

## #4

## 6
file <- 'data/fl_500_7'

## 7
file <- 'data/fl_1000_2'

params <- read.table(file, nrows = 1)
names(params) <- c('N', 'M')

## N Facilities
f <- read.table(file, skip = 1, nrows = params$N)
names(f) <- c('s', 'cap', 'x', 'y')

## M Customers
c <- read.table(file, skip = 1 + params$N, nrows = params$M)
names(c) <- c('d', 'x', 'y')

## break the plane into rectangular areas based on x- and y- ranges of facilities

x.N <- 4
y.N <- 4

x.q <- quantile(f$x, probs = seq(0, 1, 1./x.N))
y.q <- quantile(f$y, probs = seq(0, 1, 1./y.N))

x.breaks <- c(-Inf, x.q[c(-1, -length(x.q))], Inf)
y.breaks <- c(-Inf, y.q[c(-1, -length(y.q))], Inf)

f$x.rect <- cut(f$x, breaks = x.breaks, right = TRUE)
f$y.rect <- cut(f$y, breaks = y.breaks, right = TRUE)

## apply to customers
c$x.rect <- cut(c$x, breaks = x.breaks, right = TRUE)
c$y.rect <- cut(c$y, breaks = y.breaks, right = TRUE)

require(foreach)
require(doParallel)

require(Matrix)

require('Rglpk')

## 32bit R only
require('gurobi')


require(Rmosek)

cl <- makeCluster()

registerDoParallel(2)

## See if the spatial distribution of facilities and customers is as
## uniform as they claim...
foreach(i = levels(c$x.rect), .combine = rbind) %:% 
    foreach(j = levels(c$y.rect), .combine = c) %dopar% {
        nrow(subset(f, x.rect == i & y.rect == j))
    }

foreach(i = levels(c$x.rect), .combine = rbind) %:% 
    foreach(j = levels(c$y.rect), .combine = c) %dopar% {
        nrow(subset(c, x.rect == i & y.rect == j))
    }

fOptim <- function(f, c, params, solver, opt) {
    ## Coordinates as lists of pairs x, y

    ## There is got to be a better way...
    width.f <- max(nchar(rownames(f)))
    width.c <- max(nchar(rownames(c)))


    fs <- split(f[, c('x', 'y')],
                formatC(rownames(f), width = width.f,
                        format = 'd', flag = '0'))
    cs <- split(c[, c('x', 'y')],
                formatC(rownames(c), width = width.c,
                        format = 'd', flag = '0'))


    ## N by M matrix of distances
    t <- outer(fs, cs, FUN = Vectorize(function(a, b) dist(rbind(a, b))))

    ## Objective function: all facilities, then distance matrix columnwise
    ## (so that variables for each customer go together)
    obj <- c(f$s, as.vector(t))

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

    print(dim(lhs1))
    ## One facility for a customer: Sum_{w} y_{w, c} = 1

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

    lhs3 <- cBind(
        Matrix(0, nrow = params$N, ncol = params$N),
        m3)

    rhs3 <- f$cap
    rel3 <- rep('<=', params$N)

    mat <- rBind(lhs1, lhs2, lhs3)

    dir <- c(rel1, rel2, rel3)

    rhs <- c(rhs1, rhs2, rhs3)

    max <- FALSE

    types <- rep('B', params$N * (1 + params$M))

    solver(obj, mat, dir, rhs, max, types, opt, params)
}
solver.Rglpk <- function(obj, mat, dir, rhs, max, type, opt, params) {
    res <- Rglpk_solve_LP(obj, mat, dir, rhs, max = max, types = types)

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
}
solver.gurobi <- function(obj, mat, dir, rhs, max, type, opt, params) {
    model  <- list()
    
    model$A      <- mat
    model$obj    <- obj

    model$sense  <- ifelse(dir == '==', '=', dir)

    model$modelsense <- 'min'

    model$rhs    <- rhs
    model$vtype  <- "B"
    
                                        #gparams <- list(Presolve=2, TimeLimit=100.0)
    gparams <- opt
    
    res <- gurobi(model, gparams)


    ## res <- gurobi(obj, mat, dir, rhs, max = max, types = types)

                                        #cat(paste0(round(res$optimum), ' ', abs(1 - res$status), '\n'))
    cw <- sapply(
        split(res$x[-(1:params$N)],
              sapply(1:params$M, function(c) rep(c, params$N))),
        function(c) which.max(c))
    
    list(res$objval, cw)
}
solver.mosek <- function(obj, mat, dir, rhs, max, type, opt) {

    problem <- list(sense = "min")

    problem$A <- mat

    problem$c <- obj

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

    problem$dparam <- opt$dparam
    
    res <- mosek(problem, opts = list(soldetail = 1))

    ## cw <- sapply(
    ##     split(res$x[-(1:params$N)],
    ##           sapply(1:params$M, function(c) rep(c, params$N))),
    ##     function(c) which.max(c))

                                        #list(res$pobjval, res$sol)
    return(res)
}

capture.output.helper <- function (file, code) {
    capture.output(file = file, res <- code(), print ("Done"))
    return (res)
}

optimRes <- setRefClass('optimRes',
                        fields = c('c.idx', 'f.idx', 'res'))

## mosek
res <- foreach(i = levels(c$x.rect), .combine = list,
                .packages = c('Rglpk', 'Matrix')) %:% 
    foreach(j = levels(c$y.rect), .combine = list,
            .packages = c('Rglpk', 'Matrix')) %dopar% {
        c.idx <- which(c$x.rect == i & c$y.rect == j)
        f.idx <- which(f$x.rect == i & f$y.rect == j)

        p <- list(N = length(f.idx), M = length(c.idx))
        
        list(c.idx=c.idx, f.idx=f.idx,
             res = fOptim(f = f[f.idx, ], c = c[c.idx, ], p))
    }


opt.mosek <- list(dparam = list(OPTIMIZER_MAX_TIME = 3600))

res <- fOptim(f = f, c = c, params = params,
              solver = solver.mosek, opt = opt.mosek)


## gurobi
opt.gurobi <- list(Presolve=2,
                   ##TuneOutput = 0,
                   ##OutputFlag = 0,
                   TimeLimit = 600.,
                   ConcurrentMIP =  4)
    
res <- fOptim(f = f, c = c, params = params,
              solver = solver.gurobi, opt = opt.gurobi)

res <- foreach(i = levels(c$x.rect),
               .combine = list, .multicombine = F) %:% 

    foreach(j = levels(c$y.rect), .combine = list, .multicombine = F,
            .packages = c('gurobi', 'Matrix')) %dopar%

        capture.output.helper(
            file = paste0('output_', i, '_', j),
            code = function() {
                c.idx <- which(c$x.rect == i & c$y.rect == j)
                f.idx <- which(f$x.rect == i & f$y.rect == j)
                
                p <- list(N = length(f.idx), M = length(c.idx))
                optimRes <- setRefClass('optimRes',
                        fields = c('c.idx', 'f.idx', 'res'))
                optimRes(c.idx = c.idx, f.idx = f.idx,
                     res = fOptim(f = f[f.idx, ], c = c[c.idx, ],
                         params = p,
                         solver = solver.gurobi,
                         opt = opt.gurobi))
            })

res.flat <- unlist(res)

res.full <- rep(0, params$M)
res.obj <- 0

for(i in 1:length(res.flat)) {
    res.full[res.flat[[i]]$c.idx] <-
        res.flat[[i]]$f.idx[res.flat[[i]]$res[[2]]]

    res.obj <- res.obj + res.flat[[i]]$res[[1]]
}

conn <- file('sol7.txt', 'w')
writeLines(paste0(round(res.obj), ' ', 0), conn)
writeLines(paste0(res.full - 1, collapse = ' '), conn)
close(conn)


stopCluster(cl)

cat(paste0(round(res[[1]]), ' ', 0, '\n'))

## cw <- sapply(
##     split(res$solution[-(1:params$N)],
##           sapply(1:params$M, function(c) rep(c, params$N))),
##     function(c) which.max(c))

cat(paste0(res[[2]] - 1, collapse = ' '))
cat('\n')

conn <- file('sol7.txt', 'w')
writeLines(paste0(round(res[[1]]), ' ', 1), conn)
writeLines(paste0(res[[2]] - 1, collapse = ' '), conn)
close(conn)

