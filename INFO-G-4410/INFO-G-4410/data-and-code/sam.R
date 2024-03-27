library(samr)

##load scrambled data
load("thyroid-1.Rda")
ls()


##create SAM input
cl <- as.vector(sapply(1:12, function(x) c(x,-x)))
data <- list(x=M,
             y=cl,
             geneid=rownames(M),
             genenames=ga[, "SYMBOL"],
             logged2=TRUE)

##run SAM
sam.d <- samr(data, resp.type="Two class paired", nperms=100)
delta <- 0
samr.plot(sam.d, delta)



##load real data
load("thyroid-2.Rda")

##create SAM input
data <- list(x=M, y=cl, geneid=rownames(M), genenames=ga[, "SYMBOL"], logged2=TRUE)

##run SAM
sam.d <- samr(data, resp.type="Two class paired", nperms=100) ##better results with nperm=1000, but slow

delta.table <- samr.compute.delta.table(sam.d)
edit(delta.table)

delta <- 1.15
samr.plot(sam.d, delta)

table <- samr.compute.siggenes.table(sam.d, delta, data, delta.table, all.genes=TRUE)
x <- edit(table$genes.up)
dim(x)
