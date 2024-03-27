#####################################################################################################
## The experiment:
## Swirl is a point 
## mutant in the BMP2 gene that aﬀects the dorsal/ventral body axis. Ventral fates such as blood 
## are reduced, whereas dorsal structures such as somites and notochord are expanded. A goal of 
## the Swirl experiment is to identify genes with altered expression in the swirl mutant compared to 
## wild–type zebraﬁsh. Two sets of dye–swap experiments were performed, for a total of four replicate 
## hybridizations. For each of these hybridizations, target cDNA from the swirl mutant was labeled 
## using one of the Cy3 or Cy5 dyes and the target cDNA wild–type mutant was labeled using the 
## other dye. Target cDNA was hybridized to microarrays containing 8,448 cDNA probes, including 
## 768 controls spots (e.g. negative, positive, and normalization controls spots). Microarrays were 
## printed using 4 
## × 4 print–tips and are thus partitioned into a 4 × 4 grid matrix. Each grid consists 
## of a 22 
## × 24 spot matrix that was printed with a single print–tip. Here, spot row and plate coor- 
## dinates should coincide, as each row of spots corresponds to probe sequences from the same 384 
## well–plate. 


library(marray) ##load library for dual channel array preprocessing
data(swirl) ##load the SWIRL data set

summary(swirl) ##overview of the swirl object

##look at array #1
maImage(swirl[, 1], x="maGf")
maImage(swirl[, 1], x="maGb") ##see dark line

maImage(swirl[, 1], x="maRf") ##line show again!
maImage(swirl[, 1], x="maRb") ##here it is again

maImage(swirl[, 1], x="maM") ##here it is again


##look at array #3
maImage(swirl[, 3], x="maGf")
maImage(swirl[, 3], x="maGb") ##see the left to right trend 

maImage(swirl[, 3], x="maRf") ##block 3,3 seems more faint
maImage(swirl[, 3], x="maRb")

maImage(swirl[, 3], x="maM") ##block 3,3 is yellowish, more blue in the left
maImage(swirl[, 3], x = "maM", subset = maTop(maM(swirl[, 3]), h = 0.1, l = 0.1))
boxplot(swirl[,3], yvar = "maM") ##box plots confirm this

hist(maGf(swirl)) #need to use log
plot(log2(maGf(swirl[,1])),log2(maRf(swirl[,1])))


##look at intensity biases
maPlot(swirl[,1])
maPlot(swirl[,2])
maPlot(swirl[,3])
maPlot(swirl[,4])


##normalization
norm <- maNorm(swirl)

maPlot(norm[,1])
maPlot(norm[,2])
maPlot(norm[,3])
maPlot(norm[,4])


##check slide #3
maImage(norm[, 3], x="maM") ##look at block 3-3 on slide 3
boxplot(norm[,3], yvar = "maM") ##box plots confirm this


##average dye swap replicate
summary(norm)
M <- cbind((maM(norm)[,1]-maM(norm)[,2])/2, (maM(norm)[,3]-maM(norm)[,4])/2)
head(M)
plot(M) 
abline(0,1, col="red")
##Not huge signal. What's about the right hand M1-M2 tail?



