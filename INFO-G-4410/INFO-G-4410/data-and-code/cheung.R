



#load and inspect data
load("ng-cheung-2007.Rda")
ls()
dim(e)
hist(e)
origin
table(origin)

colors <- rep("black", length(origin))
colors[origin=="CHB"] <- "red"
colors[origin=="CHLA"] <- "blue"
colors[origin=="JPT"] <- "green"


##HC
d <- as.dist(1-cor(e))

sclus <- hclust(d, method = "complete")
plot(sclus, main="complete", xlab="", sub="", cex=0.8, col=colors)

sclus <- hclust(d, method = "average")
plot(sclus, main="average", xlab="", sub="", cex=0.8, col=colors)

sclus <- hclust(d, method = "ward")
plot(sclus, main="ward", xlab="", sub="", cex=0.8, col=colors)


##PCA
pca <- prcomp(t(e))

plot(pca$x[, 1:2], col=colors)
plot(pca$x[, 3:4], col=colors)
plot(pca$x[, 5:6], col=colors)

plot(pca)
x <- 100*pca$sdev^2/sum(pca$sdev^2)
plot(x)
plot(cumsum(x))


##MDS
library(MASS)
d <- as.dist(1-cor(e))
mds <- isoMDS(d, k=2, maxit=1000, tol=1e-20)
plot(mds$points, col=colors)
mds$stress


##Look at correlation matrix
c <- 1-cor(e)
plot(density(c), xlab="distance")
lines(density(c[origin=="CEU", origin=="CEU"]), col="red")
lines(density(c[origin=="CEU", origin %in% c("CHB", "CHLA", "JPT")]), col="green")
lines(density(c[origin=="CHLA", origin=="CHLA"]), col="blue")

##Same with histograms
par(mfrow=c(3,1))
par(mar=c(1,1,1,1)) #Line added
hist(c, xlab="distance")
hist(c[origin=="CEU", origin=="CEU"], col="red")
hist(c[origin=="CEU", origin %in% c("CHB", "CHLA", "JPT")], xlim=c(0, 0.15), col="green")
par(mfrow=c(1,1))


#SAM Europeans vs. Asians
library(samr)

cl <- rep(1, length(origin))
cl[origin=="CHB"] <- 2
cl[origin=="CHLA"] <- 2
cl[origin=="JPT"] <- 2

data <- list(x=e,
             y=cl,
             geneid=rownames(e),
             genenames=rownames(e),
             logged2=TRUE)

sam.d <- samr(data, resp.type="Two class unpaired", nperms=100) ##better results with nperm=1000, but slow
delta.table <- samr.compute.delta.table(sam.d)
invisible(edit(delta.table))

delta <- 0
samr.plot(sam.d, delta)

table <- samr.compute.siggenes.table(sam.d, delta, data, delta.table, all.genes=TRUE)
x <- edit(table$genes.up)
colnames(table$genes.up)
hist(as.numeric(table$genes.up[table$genes.up[,8] <=5, 7]), breaks=50)









x <- cbind(rnorm(100), rnorm(100))
png("pca-2d-before.png")
plot(x %*% cbind(c(1, 1.5), c(0.5, 2.5)), xlab="x", ylab="y", pch=19)
dev.off()

x <- prcomp(x)$x
png("pca-2d-after.png")
plot(x, pch=19)
dev.off()
