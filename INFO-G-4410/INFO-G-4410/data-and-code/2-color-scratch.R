##investigate slide 1, scratch

library(marray) ##load library for dual channel array preprocessing
data(swirl) ##load the SWIRL data set


norm <- maNorm(swirl)
maImage(norm[, 1], x="maM") ##here it is again


norm <- maNorm(swirl, norm="loess")
norm <- maNorm(norm, norm="twoD")
maPlot(norm[,1])
maImage(norm[, 1], x="maM") ##here it is again
maImage(norm[, 1], x="maMloc") ##here it is again

maImage(norm[, 3], x="maMloc") 

plot(log10(maGb(swirl[,1])), log10(maRb(swirl[,1]))) ##funny spots here
abline(h=2.4, col="red")
x <- log10(maRb(swirl[,1]))>2.4
maImage(swirl[, 1], x = "maM", subset = x) ##ok, easy to filter them out...


plot(log10(maGf(swirl[,1])),log10(maRf(swirl[,1])))
points(log10(maGf(swirl[x,1])),log10(maRf(swirl[x,1])), col="red")

maPlot(norm[,1])
points(norm[x,1], col="red")


##average dye swap replicate
M <- cbind((maM(norm)[,2]-maM(norm)[,1])/2, (maM(norm)[,4]-maM(norm)[,3])/2)
plot(M) 
points(M[x,], col="red")
abline(0,1, col="red")

y <- x & (M[,1] < -1)
maImage(swirl[, 1], x = "maM", subset = y) 

##so 'useful' data is more like this
plot(M[!y,]) 
abline(0,1, col="red")

getClassDef("marrayNorm")
getClassDef("marrayInfo")
bmp2 <- grep("^BMP2$", norm@maGnames@maLabels)
bmp2
points(M[bmp2,], col="green", pch=19)
points(M[y,], col="orange", pch="+")

ctl <- norm@maGnames@maInfo[, 1]=="control"
points(M[ctl,], col="green", pch="+")

down.genes <- M[,1] < -1 & M[,2] < -1
unique(norm@maGnames@maLabels[down.genes])

up.genes <- M[,1] > 1 & M[,2] > 1
unique(norm@maGnames@maLabels[up.genes])

M[bmp2,]
