##http://homepages.ulb.ac.be/~vdetours/thyroid-1.Rda

setwd("../selec-de-genes/")
##load thyroid cancer data (MAS 5.0 normalization)
load("thyroid-2.Rda")
ls()

##quick look at them
dim(M)
head(M)
dim(ga)
head(ga)
cbind(ga[1:6,], M[1:6, 1:6])
colnames(M)
hist(M) ##that's log2 data

##'T' is tumor, 'N' is normal
##let's compute fold change T/N
fc <- M[, seq(2, 24, 2)] - M[, seq(1, 23, 2)]
hist(fc)
dim(fc)

##now we select differentially expressed genes, ie 2-fold in 2/3 of 12 patients
de <- apply(fc, 1, function(x) 8 <= sum(x >= 1)) |
      apply(fc, 1, function(x) 8 <= sum(x <= -1))
sum(de)
med <- apply(fc, 1, median)
hist(med[de])

##let's have look at the top genes
cbind(med[de], ga[de, c("SYMBOL", "GENENAME")])

##same ordered
cbind(med[de], ga[de, c("SYMBOL", "GENENAME")])[order(med[de]),]


##run random, negative, control
z <- c()
for (i in (1:100)) {
  M <- M[, sample(1:ncol(M))]
  fc <- M[, seq(2, 24, 2)] - M[, seq(1, 23, 2)]
  de <- apply(fc, 1, function(x) 8 <= sum(x >= 1)) | apply(fc, 1, function(x) 8 <= sum(x <= -1))
  z <- c(z, sum(de))
}

summary(z)

hist(z)
abline(v=28, col="red")
sum(z >= 28)/100






































































































































load("/Users/detours/data/genriskt/compare-gliwice-vs-iribhm/genriskt-gcrma-gene-wise-average.Rda")
M <- M[,1:24]
save(M, ga, file="/Users/detours/BIOL-F-523/selec-de-genes/thyroid-2.Rda", compress=T)

c <- colnames(M)
M <- M[,sample(1:24)]
colnames(M) <- c
save(M, ga, file="/Users/detours/BIOL-F-523/selec-de-genes/thyroid-1.Rda", compress=T)
