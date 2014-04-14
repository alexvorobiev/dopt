# Return full path to problem data file.
dataFileFullPath <- function(fname) {
  paste("data", fname, sep="/")
}

# Returns solution for a given problem.
computeSolutionFor <- function(fname) {
  # This is for the standard, python solver
  sol <- system(paste('python solver.py', dataFileFullPath(fname)), intern=T)
  #sol <- system(paste('./facility', dataFileFullPath(fname)), intern=T)
  cost <- as.numeric(substr(sol[1], 1, nchar(sol[1])-1))
  solution <- as.numeric(unlist(strsplit(sol[2], " ")))
  list(cost=cost, solution=solution)
}

# Create a color mapping function between the given vector (which may contain duplicate values)
# and a color for each unique value in the vector.
makeColorMapper <- function(vals) {
  uniq <- sort(unique(vals))
  colors <- rainbow(length(uniq))
  function(val) {
    colors[which(uniq == val)]
  }
}

# Loads one data file and returns a list with facilities (fx) and customers (cx).
loadDataFile <- function(fname) {
  rawData <- read.csv(dataFileFullPath(fname), header=F, sep=' ')
  facCount <- rawData$V1[1]; custCount <-rawData$V2[1]

  facilities <- rawData[2:(1+facCount), 1:4]
  colnames(facilities) <- c('setupCost', 'capacity', 'x', 'y')

  customers <- rawData[(2+facCount):(1+facCount+custCount), 1:3]
  colnames(customers) <- c('demand', 'x', 'y')

  list(fx = facilities, cx = customers)
}

# Plots one instance of the Facility Location problem defined in fname, including
# the solution (obtained via python solver.py fname).
plotInstance <- function(fname, plotLegend=F, plotSol=T) {
  z <- loadDataFile(fname)
  facilities <- z$fx
  customers <- z$cx

  plot(customers$y ~ customers$x, col="red", main=paste("Facility Location",fname), xlab="x", ylab="y",
       xlim=c(min(c(min(customers$x), min(facilities$x))),max(c(max(customers$x), max(facilities$x)))),
       ylim=c(min(c(min(customers$y), min(facilities$y))),max(c(max(customers$y), max(facilities$y)))))
  points(facilities$y ~ facilities$x, col="blue")

  if (plotLegend) {
    colors <- c("red", "blue", "blue")
    legend("topleft", legend=c("customer", "facility", "facility open"), col=colors,
           text.col=colors, y.intersp=0.85, pch=c(1,1,16))
  }

  if (plotSol) {
    plotSolution(fname, customers, facilities)
  }
}

# Plots solution for one instance of the problem, defined in fname.
plotSolution <- function(fname, customers=NULL, facilities=NULL) {
  # Lazy load customers & facilities
  if (is.null(customers) || is.null(facilities)) {
    z <- loadDataFile(fname)
    facilities <- z$fx
    customers <- z$cx
  }

  z <- computeSolutionFor(fname)
  # converting from 0 based arrays of the output to 1 based arrays of R
  z$solution <- z$solution + rep(1,length(z$solution))
  facSel <- facilities[sort(unique(z$solution)),]
  colors <- makeColorMapper(z$solution)
  mtext(formatC(z$cost, format="d", big.mark=","), line=0.2, outer=F, col="#CC2211")
  points(facSel$y ~ facSel$x, col="blue", pch=16)

  ci <- 0
  for (fi in z$solution) {
    ci <- ci + 1
    f <- facilities[fi,]
    c <- customers[ci, ]
    segments(f$x, f$y, c$x, c$y, col=colors(fi))
  }
}

# Plots multiple data instances (all found in _metadata by default).
# If only specific instances need to be plotted, then limit can be passed a vector
# of their IDs (e.g. c(1,3,5) would only plot problems #1, #3 and #5).
plotMulti <- function(limit=NULL, plotSol=T, fname = './_metadata') {
  z <- read.table(fname, header=F, sep=' ', skip=3)
  names <- lapply(as.character(z$V2), function (a) substr(a,8,nchar(a)-1))
  if (!is.null(limit)) {
    names <- names[limit]
  }

  if (length(names) > 1) {
    par(mfrow=c(2,length(names)/2))
  } else {
    par(mfrow=c(1,1))
  }

  lapply(names, function (a) plotInstance(a, a == names[[1]], plotSol))
}

#png(width=1600, height=900, res=96, antialias="subpixel", filename="out.png")
#plotMulti()
#dev.off()

# Manually plotting selected problems:
#plotInstance('fl_25_2')
