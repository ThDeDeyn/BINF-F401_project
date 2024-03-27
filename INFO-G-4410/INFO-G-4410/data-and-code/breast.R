##http://homepages.ulb.ac.be/~vdetours/breast.Rda

load("breast.Rda")
ls()
names(tab)
sapply(tab, dim)
attach(tab)
colnames(sa)
attach(sa)
head(sa)
hist(sa[, "time.disease.reccurence"], breaks=20)
table(sa[, "event.disease.reccurence"])

median(sa[sa[, "event.disease.reccurence"]==0, "time.disease.reccurence"])
median(sa[sa[, "event.disease.reccurence"]==1, "time.disease.reccurence"])

sa[1:10, c("time.disease.reccurence", "event.disease.reccurence")]

library(survival)

s <- Surv(time.disease.related.death, event.disease.related.death)
event.disease.related.death[1:50]
s[1:50]
plot(survfit(s ~ 1)) ##plot KM curve

hist(time.disease.related.death[event.disease.related.death==1], breaks=20)
hist(time.disease.related.death[event.disease.related.death==0])

table(ER)
plot(survfit(s ~ ER))
summary(coxph(s ~ ER))
summary(coxph(s ~ ER))$logtest["pvalue"]

table(histologic.grade)
plot(survfit(s ~ histologic.grade))
summary(coxph(s ~ histologic.grade))

x <- histologic.grade
x[x==3] <- 2
summary(coxph(s ~ x))

table(lymph.node.positive)
plot(survfit(s ~ lymph.node.positive))
summary(coxph(s ~ lymph.node.positive))


table(treatment)
plot(survfit(s ~ treatment))
summary(coxph(s ~ treatment))
sum(treatment=="tamoxifen" & ER==1)

hist(tumor.size)
x <- tumor.size <= 1
plot(survfit(s ~ x))
summary(coxph(s ~ tumor.size))

summary(coxph(s ~ ER + tumor.size + lymph.node.positive))

summary(coxph(s ~ ER + x + lymph.node.positive))

##signature
##http://homepages.ulb.ac.be/~vdetours/sigs.Rda
load("sigs.Rda")
ls()
names(sigs.ltr)
ggi <- sigs.ltr[["sotiriou.2006"]]
ggi
length(ggi)

head(M)
x <- M[intersect(rownames(M), ggi),]
dim(x)
c <- 1 - cor(x)
clust <- hclust(as.dist(c), method="ward")
pred <- cutree(clust, 2)
pred[1:50]

table(pred)
table(pred, histologic.grade)
plot(survfit(s ~ pred))
summary(coxph(s ~ pred))

summary(coxph(s ~ pred + ER + tumor.size + lymph.node.positive))

#### A control
x <- M[sample(rownames(M), 100),]
dim(x)
c <- 1 - cor(x)
clust <- hclust(as.dist(c), method="ward")
pred <- cutree(clust, 2)
pred[1:50]

table(pred)
plot(survfit(s ~ pred))
summary(coxph(s ~ pred))

