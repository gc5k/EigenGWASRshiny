# lambda1 mixture distribution / genetic drift only
M=10000
f=0.01 # genetic drift

P=runif(M, 0.2, 0.8) # ancestral
Sz=c(300, 300) #subgroup sample size
Fq=matrix(0, 2, M) 
Fq[1,]=rbeta(M, P*(1-f)/f, (1-P)*(1-f)/f)
Fq[2,]=rbeta(M, P*(1-f)/f, (1-P)*(1-f)/f)

G=matrix(0, sum(Sz), M)
for(i in 1:Sz[1]) {
  G[i,] = rbinom(M, 2, Fq[1,])
}

for(i in (Sz[1]+1):(Sz[1]+Sz[2])) {
  G[i,] = rbinom(M, 2, Fq[2,])
}

fq1 = colMeans(G[1:Sz[1],])/2 # observed allele freq, in subgroup 1
fq2 = colMeans(G[(Sz[1]+1):(Sz[1]+Sz[2]),])/2 # observed allele freq, in subgroup 2
H1 = 2*fq1*(1-fq1)
H2 = 2*fq2*(1-fq2)
Hs = (H1+H2)/2
HT = 2*((fq1+fq2)/2)*(1-(fq1+fq2)/2)
#Fst = sum((HT-Hs)/HT)
Fst = mean((fq1-fq2)^2/(HT))*(Sz[1]+Sz[2])

Gs=apply(G, 2, scale)
GG=Gs %*% t(Gs) / M
EigenG=eigen(GG)
layout(matrix(1:2, 1, 2))
barplot(EigenG$values[1:20])
plot(main=paste(M, "markers"), EigenG$vectors[,1], EigenG$vectors[,2], bty='n', xlab="PC 1", ylab="PC 2", pch=16, col=c(rep("red", Sz[1]), rep("blue", Sz[2])))

# get lambda GC
eigenspace1 = EigenG$vectors[,1]
para = vector(length = M)
for (i in 1:M){
  model = lm(eigenspace1~Gs[,i])
  para[i]=summary(model)$coefficients[2,4] #p value
}
lambdaGC = qchisq(median(para), 1, lower.tail = F)/qchisq(0.5, 1, lower.tail = F)
lambda1 = EigenG$values[1]
layout(matrix(1:1, 1))
a = c(Fst, lambda1, lambdaGC)
barplot(a)



# lambda1 mixture distribution / genetic drift and selection
m=c(8000) #backgroud marker
mq=c(2000) #loci under selection
f1=0.01 # genetic drift
fq=0.05 #selection
f=c(rep(f1, m), rep(fq, mq))
M=length(f)

P=runif(M, 0.2, 0.8) #ancestral
Sz=c(300, 300) #subgroup sample size
Fq=matrix(0, 2, M)
Fq[1,]=rbeta(M, P*(1-f)/f, (1-P)*(1-f)/f) #freq for subgroup 1
Fq[2,]=rbeta(M, P*(1-f)/f, (1-P)*(1-f)/f) #freq for subgroup 2

G=matrix(0, sum(Sz), M)
for(i in 1:Sz[1]) {
  G[i,] = rbinom(M, 2, Fq[1,])
}

for(i in (Sz[1]+1):(Sz[1]+Sz[2])) {
  G[i,] = rbinom(M, 2, Fq[2,])
}

fq1 = colMeans(G[1:Sz[1],])/2
fq2 = colMeans(G[(Sz[1]+1):(Sz[1]+Sz[2]),])/2
H1 = 2*fq1*(1-fq1)
H2 = 2*fq2*(1-fq2)
Hs = (H1+H2)/2
HT = 2*((fq1+fq2)/2)*(1-(fq1+fq2)/2)
#Fst = sum((HT-Hs)/HT)
Fst = mean((fq1-fq2)^2/(HT))*(Sz[1]+Sz[2])

Gs=apply(G, 2, scale)
GG=Gs %*% t(Gs)/M
EigenG=eigen(GG)
layout(matrix(1:3, 1, 3))
plot(colMeans(G[1:Sz[1],1:m])/2, colMeans(G[(Sz[1]+1):(Sz[1]+Sz[2]), 1:m]/2))
points(colMeans(G[1:Sz[1],(m+1):M])/2, colMeans(G[(Sz[1]+1):(Sz[1]+Sz[2]), (m+1):M]/2), pch=16, cex=2, col="red")
abline(a=0, b=1, col="green")
barplot(EigenG$values[1:20])
plot(main=paste(M, "markers"), EigenG$vectors[,1], EigenG$vectors[,2], bty='n', xlab="PC 1", ylab="PC 2", pch=16, col=c(rep("red", Sz[1]), rep("blue", Sz[2])))

eigenspace1 = EigenG$vectors[,1]
para = vector(length = M)
for (i in 1:M){
  model = lm(eigenspace1~Gs[,i])
  para[i]=summary(model)$coefficients[2,4] #p value
}
lambdaGC = qchisq(median(para), 1, lower.tail = F)/qchisq(0.5, 1, lower.tail = F)
lambda1 = EigenG$values[1]

layout(matrix(1:1, 1))
b = c(Fst, lambda1, lambdaGC)

data = cbind(a,b)
colnames(data) = c('Genetic drift only', 'Genetic drift and selection')
color = cm.colors(3)
color = grey.colors(3)
color = terrain.colors(3)
color = c("#0072B2","#000000","#999999")
barplot(as.matrix(data),
        beside=T,
        legend=c(expression(F[st]),expression(lambda[1]),expression(lambda[GC])),
        col = color,
        width = 6,
        ylim = c(0,12))




