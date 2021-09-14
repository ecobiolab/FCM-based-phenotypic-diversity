rm(list=ls())

{
library(vegan)
library(ggplot2)  
library(caret)
library(sampling)
library(plyr)
library(ggpubr)
library(dplyr)
library(lmtest)
library(nparcomp)
}

# Add datasets=====

data = read.csv("TEST.data.csv", header=T)

# Normalize the datasets to the [0,1] range=====

data.pro1 = preProcess(data[ ,2:5], method = "range")

data.nor = predict(data.pro1, data[,2:5])
data.nor$key = data$key
data.nor$stat = data$stat

summary(data.nor)

# Randomly resample to the lowest sample size=====

min(table(data.nor$key)) # selected the number of sample size > ex) 3383

data.nor.re = ddply(data.nor,.(key),function(x) x[sample(nrow(x),3383),])

table(data.nor.re$key)

# data binning=====

DAT = data.nor.re
DAT

DAT$fg  = NA
xbin = range(DAT$FSC.A)
ybin = range(DAT$SSC.A)

Fbin = seq(round_any(xbin[1], 1, trunc), round_any(xbin[2], 1, ceiling), by = 0.1)
Sbin = seq(round_any(ybin[1], 1, trunc), round_any(ybin[2], 1, ceiling), by = 0.1)

for (i in 1:(length(Fbin)-1)){
  for (j in 1:(length(Sbin)-1)){
    DAT[DAT$FSC.A >= Fbin[i] & DAT$FSC.A < Fbin[i+1] & 
          DAT$SSC.A >= Sbin[j] & DAT$SSC.A < Sbin[j+1], "fg"] = paste0("F", i, "S", j)
  }
}

ID.fg = unique(DAT$fg)
ID.name = unique(DAT$key)
figname = unlist(strsplit(unique(DAT$key), ".csv"))
col.jk = rainbow(length(ID.fg))

freq = as.data.frame(unique(DAT$fg)); colnames(freq) = "Var1"

fn = 1; i = 1
for (fn in 1:length(ID.name)){
  slct = DAT[DAT$key == ID.name[fn], ]
  fre = as.data.frame(table(slct$fg))
  freq = merge(freq, fre, by = "Var1", all.x = T)
  colnames(freq)[fn+1] = ID.name[fn]
}  

data.freq = t(freq)


## calculate the alpha-diversity====

data.freq[is.na(data.freq)] = 0
colnames(data.freq) = data.freq[1,]
data.freq = data.freq[-1,-67]
data.freq <- as.data.frame(data.freq)
for (i in 1:length(data.freq)){
  data.freq[,i] <- as.numeric(data.freq[,i])
  }
str(data.freq)


diversity(data.freq, index="shannon") # calculation of shannon-Wiener index
exp(diversity(data.freq, index="shannon")) # calculation of effective (true) shannon-Wiener index
diversity(data.freq, index="simpson") # calculation of index
abundance = apply(data.freq>0,1,sum) 
diversity(data.freq, index="simpson")/log(abundance) # calculation of pilou evenness 


## visualization of diversity index
data.sha =  data.frame(diversity(data.freq, index="shannon"))
data.sha.e =  data.frame(exp(diversity(data.freq, index="shannon")))
data.sim =  data.frame(diversity(data.freq, index="simpson"))
data.eve =  data.frame(diversity(data.freq, index="simpson")/log(abundance))

data.div=do.call("cbind", list(data.sha, data.sha.e, data.sim, data.eve))
colnames(data.div) = c("shannon","shannon.e","simpson","evenness")

#data.div$label= data.freq$label
data.div$label[1:10] = c("Clean")
data.div$label[11:15] = c("Contaminated")
head(data.div)


ggplot(data=data.div, aes(x=label, y=shannon, fill=label))+
  geom_boxplot()+
  scale_fill_manual(values=c("#0388A6", "#F2C12E"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))

## calculate the beta-diversity =====
## PCA ====

pca = prcomp(data.freq, center=T, retx=TRUE)
summary(pca)
names(pca)
a=as.data.frame(pca$x)
a=a[,c(1,2)]
a$label = data.div$label

data.freq$label = data.div$label

ggplot(data=a, aes(x=PC1, y=PC2, col=label))+
  geom_point()+
  geom_point(aes(fill = label), shape = 21, color="black", size=6)+
  geom_hline(yintercept=0, lty=2, col="grey")+
  geom_vline(xintercept=0, lty=2, col="grey")+
  scale_fill_manual(values=c("#0388A6", "#F2C12E"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12))+
  xlab("PC1 (57.9 %)")+ylab("PC2 (14.9 %)")

dist <- vegdist(data.freq[3:40], method = "euclidean")
adonis(dist ~ label, data = data.freq)
