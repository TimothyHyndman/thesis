

library(deconvolve)
library(ggplot2)
library(tictoc)
library(R.matlab)


set.seed(5000)
n <- 500
NSR <- 0.2

# X
dchi <- 3
X <- rchisq(n, df = dchi)
varX <- 2 * dchi
X <- X / sqrt(varX)
varX <- 1
cvar <- sqrt(2 * dchi)

truedens <- function(xx) {
	cvar*dchisq(xx * cvar, df = dchi)
}

# U
sigU <- sqrt(NSR * varX)
U <- rnorm(n, mean = 0, sd = sigU)

# W
W <- X + U
# Same as matlab
W <- readMat('../thesis/Figures/Deconvolution/Code/Wexample.mat')
W <- as.vector(W$W)

fWEF <- readMat('../thesis/Figures/Deconvolution/Code/naiveexample.mat')
fWEF <- as.vector(fWEF$fWEF)
h <- 0.221444788809846

tic()
decon_object <- deconvolve(W, show_diagnostics = TRUE, m = 20)
toc()

xx <- decon_object$x
yy <- truedens(xx)
df1 <- data.frame(x = xx, y = yy)
df2 <- data.frame(x = decon_object$support, y = decon_object$probweights)
df3 <- data.frame(x = xx, y = fWEF)

truedenscol <- rgb(102, 177, 255, maxColorValue = 255)
pointcol <- rgb(150, 150, 150, maxColorValue = 255)

plot(decon_object) +
    geom_line(data = df3, aes(x, y), colour = truedenscol, size = 1, linetype = "dashed") +
	geom_line(data = df1, aes(x, y), colour = truedenscol, size = 1) +
	geom_point(data = df2, aes(x, y), colour = pointcol)


save_image <- FALSE

if (save_image){
	ggsave("../thesis/Figures/Deconvolution/R_example.png")
}
