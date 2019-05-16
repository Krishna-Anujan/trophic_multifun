#bef_multifun_biorXiv.R is the code used in the paper titled "Trophic complexity alters the diversity-multifunctionality relationship in experimental grassland mesocosms"

bef<-read.csv("bef_troph_bioRxiv.csv")

library(wesanderson)
library("png")
library("corrplot")

head(bef)

bef$rslnc<-bef$ABOVETOT02/bef$TOTBIOM01
bef$abvbiom<-bef$TOTBIOM01-bef$ROOT01
bef$biomstd<-(bef$abvbiom-min(bef$abvbiom))/(max(bef$abvbiom)-min(bef$abvbiom))
bef$flostd<-(bef$FLOWTIME -min(bef$FLOWTIME))/(max(bef$FLOWTIME)-min(bef$FLOWTIME))
bef$recstd<-(bef$rslnc -min(bef$rslnc))/(max(bef$rslnc)-min(bef$rslnc))
bef$rtstd<-(bef$ROOT01 -min(bef$ROOT01))/(max(bef$ROOT01)-min(bef$ROOT01))
bef$INTERACT<-paste(bef$INS, bef$LIT, sep="")
bef$INTERACT[bef$INTERACT=="00"]<-"NONE"
bef$INTERACT[bef$INTERACT=="01"]<-"LIT"
bef$INTERACT[bef$INTERACT=="10"]<-"INS"
bef$INTERACT[bef$INTERACT=="11"]<-"BOTH"

bef$INDEX<-rep(0, nrow(bef))
bef$INDEX[bef$INTERACT=="NONE"]<-1
bef$INDEX[bef$INTERACT=="LIT"]<-2
bef$INDEX[bef$INTERACT=="INS"]<-3
bef$INDEX[bef$INTERACT=="BOTH"]<-4

col<-wes_palette("Royal1", 4, type="continuous")

#Multifunctionality using maximum as the average of the top 5-------------------------
biomax<-mean(bef$TOTBIOM01[order(bef$abvbiom, decreasing=T)][1:5])
recmax<-mean(bef$rslnc[order(bef$rslnc, decreasing=T)][1:5])
rtmax<-mean(bef$ROOT01[order(bef$ROOT01, decreasing=T)][1:5])
flomax<-mean(bef$FLOWTIME[order(bef$FLOWTIME, decreasing=T)][1:5])

bef$biomstd<-(bef$abvbiom-min(bef$abvbiom))/(biomax-min(bef$abvbiom))
bef$flostd<-(bef$FLOWTIME -min(bef$FLOWTIME))/(flomax-min(bef$FLOWTIME))
bef$recstd<-(bef$rslnc -min(bef$rslnc))/(recmax-min(bef$rslnc))
bef$rtstd<-(bef$ROOT01 -min(bef$ROOT01))/(rtmax-min(bef$ROOT01))


mult<-bef[,c(2,5, 71, 72, 73, 74, 76)]
mult$l90<-rep(0, nrow(mult))
mult$l75<-rep(0, nrow(mult))
mult$l50<-rep(0, nrow(mult))
mult$l25<-rep(0, nrow(mult))
levs<-c(90, 75, 50, 25)

for(i in 1:nrow(mult)){
  for (j in 1:length(levs)){
    lev<-levs[j]
    count<-0
    for (k in 3:6){
      if (mult[i, k]>=0.01*lev) {count<-count+1}
    }
    mult[i, j+7]<-count
  }
}

mult$INDEX<-bef$INDEX
mult$fdis<-bef$fdis
mult$mpd<-bef$mpd

mult$sp.std<-rescale(mult$SPP)
mult$fdis.std<-rescale(mult$fdis)
mult$mpd.std<-rescale(mult$mpd)

fit.90<-lm(l90~SPP+INDEX+SPP*INDEX, mult)
fit.90<-lm(l90~sp.std+INDEX+sp.std*INDEX, mult)
fit.90<-lm(l90~sp.std+INDEX, mult)
summary(fit.90)

fit.75<-lm(l75~SPP+INDEX+SPP*INDEX, mult)
fit.75<-lm(l75~sp.std+INDEX+sp.std*INDEX, mult)
fit.75<-lm(l75~sp.std+INDEX, mult)
summary(fit.75)

fit.50<-lm(l50~SPP+INDEX+SPP*INDEX, mult)
fit.50<-lm(l50~sp.std+INDEX+sp.std*INDEX, mult)
fit.50<-lm(l50~sp.std+INDEX, mult)
summary(fit.50)

fit.25<-lm(l25~SPP+INDEX+SPP*INDEX, mult)
fit.25<-lm(l25~sp.std+INDEX+sp.std*INDEX, mult)
fit.25<-lm(l25~sp.std+INDEX, mult)
summary(fit.25)

library(lme4)
lme.90<-lmer(l90~sp.std+(1|INDEX), mult)
lme.75<-lmer(l75~sp.std+(1|INDEX), mult)
lme.50<-lmer(l50~sp.std+(1|INDEX), mult)
summary(lme.50)

fit.90<-lm(l90~fdis.std+INDEX+fdis.std*INDEX, mult)
summary(fit.90)

fit.75<-lm(l75~fdis.std+INDEX+fdis.std*INDEX, mult)
summary(fit.75)

fit.50<-lm(l50~fdis.std+INDEX+fdis.std*INDEX, mult)
summary(fit.50)

fit.25<-lm(l25~fdis.std+INDEX+fdis.std*INDEX, mult)
summary(fit.25)

#Figure 3----------------------------------------------------------

png("thresh_4_final_nee.png", width=8, height=8, units="in", res=300)
par(xpd=F)

par(mfrow=c(2,2), oma=c(2,0,0,0), mai=c(1, 1, 0.1, 0.1), xpd=F)
#par(mfrow=c(1,1), xpd=F)
plot(jitter(l90)~SPP, mult, pch=16, #main="90% threshold", 
     xlab="Number of species", ylab="Number of functions \nabove threshold", cex.lab=1.3,
     ylim=c(-0.5,4.5), 
     col="grey")
abline(lm(l90~SPP, mult), lwd=3, lty=2)
legend("topright", "90%", cex=1.5, bty="n")
for (i in 1:4){
  sub<-subset(mult, INDEX==i)
  abline(lm(l90~SPP, sub), col=col[i], lwd=3, lty=3)
}


plot(jitter(l75)~SPP, mult, pch=16, #main="75% threshold", 
     xlab="Number of species", ylab="Number of functions \nabove threshold", cex.lab=1.3,
     ylim=c(-0.5, 4.5), 
     col="grey")
abline(lm(l75~SPP, mult), lwd=3, lty=2)
legend("topright", "75%", cex=1.5, bty="n")
for (i in 1:4){
  sub<-subset(mult, INDEX==i)
  abline(lm(l75~SPP, sub), col=col[i], lwd=3, lty=3)
}


plot(jitter(l50)~SPP, mult, pch=16, #main="50% threshold", 
     xlab="Number of species", ylab="Number of functions \nabove threshold", cex.lab=1.3,
     ylim=c(-0.5, 4.5), 
     col="grey")
abline(lm(l50~SPP, mult), lwd=3, lty=2)
legend("topright", "50%", cex=1.5, bty="n")
for (i in 1:4){
  sub<-subset(mult, INDEX==i)
  abline(lm(l50~SPP, sub), col=col[i], lwd=3, lty=3)
}


plot(jitter(l25)~SPP, mult, pch=16, #main="25% threshold", 
     xlab="Number of species", ylab="Number of functions \nabove threshold", cex.lab=1.3, ylim=c(-0.5, 4.5), col="grey")
abline(lm(l25~SPP, mult), lwd=3, lty=2)
legend("topright", "25%", cex=1.5, bty="n")
for (i in 1:4){
  sub<-subset(mult, INDEX==i)
  abline(lm(l25~SPP, sub), col=col[i], lwd=3, lty=3)
}
par(xpd=NA)
legend(x=-22, y=-2.4, bty="n", legend=c("NONE", "INS", "LIT", "BOTH"), col=col, lty=3, lwd=3, horiz = T, cex=1.3)
dev.off()


#Threshold as x-axis and slope on the y-axis (Figure 4)------------------------------------

thresh<-bef[,c(2,5,69,72,73,74,77)]

x<-seq(0, 100)
plot(1, type="n", xlim=c(0, 100), ylim=c(-0.1, 0.1), xlab="Threshold", ylab="Diversity-multifunctionality effect", cex.lab=1.5)
slope.store<-numeric(0)
ind<-numeric(0)
slope.ind<-numeric(0)
thresh.ind<-numeric(0)
slope.sigma<-numeric(0)
ind.sigma<-numeric(0)

for(i in 1:length(x)){
  t<-x[i]
  thresh.1<-thresh
  thresh.1$val<-rep(t, nrow(thresh.1))
  thresh.1$mult<-rep(0, nrow(thresh.1))
  for(j in 1:nrow(thresh.1)){
    count<-0
    for (k in 3:6){
      if (thresh.1[j, k]>=0.01*t) {count<-count+1}
    }
    thresh.1$mult[j]<-count
  }
  main.lm<-lm(mult~SPP, thresh.1)
  slope<-main.lm$coefficients[2]
  sig<-confint(main.lm)[2,]
  slope.store<-c(slope.store, slope)
  slope.sigma<-rbind(slope.sigma, sig)
  for (l in 1:4){
    sub<-subset(thresh.1, INDEX==l)
    ind.lm<-lm(mult~SPP, sub)
    slp<-ind.lm$coefficients[2]
    sig.ind<-confint(ind.lm)[2,]
    ind<-c(ind, l)
    thresh.ind<-c(thresh.ind, t)
    slope.ind<-c(slope.ind, slp)
    ind.sigma<-rbind(ind.sigma, sig.ind)

  }
}

col<-wes_palette("Royal1", 4, type="continuous")

wes_palette("BottleRocket")

png("threshold_bef_abv_nee.png")
par(mfrow=c(1,1), xpd=F, mai=c(1,1,0.5,0.5))
plot(x, slope.store, ylim=c(-0.075, 0.075), xlab="Threshold", ylab="Diversity-multifunctionality effect",
     cex=0.75, pch=16, col="grey", cex.lab=2)

lines(lowess(x, slope.store, f=0.15), lwd=3)

polygon(c(x, rev(x)), c(slope.sigma[,2], rev(slope.sigma[,1])),border=NA, col=paste(rgb(t(col2rgb("grey70"))/255), "AA", sep=""))
ind.all<-cbind(thresh.ind, ind, slope.ind, ind.sigma)

for(i in 1:4){
  sub<-ind.all[ind.all[,2]==i,]
  points(x, sub[,3], cex=0.75, col=col[i], pch=16)
  lines(lowess(x, sub[,3], f=0.15), lwd=2, col=col[i])
}
abline(h=0)
legend("topright", c("PLANT ONLY", "INS", "LIT", "BOTH", "POOLED"), col=c(col,"black"), pch=16, cex=0.75)
dev.off()  

#Plots to see which functions perform differently (Figure 2) ---------------------

head(mult)
x<-seq(0, 100)

mult_id<-as.data.frame(x)
mult_id$biom_num<-numeric(length(x))
mult_id$flo_num<-numeric(length(x))
mult_id$rec_num<-numeric(length(x))
mult_id$rt_num<-numeric(length(x))

mult_id$biom_none<-numeric(length(x))
mult_id$flo_none<-numeric(length(x))
mult_id$rec_none<-numeric(length(x))
mult_id$rt_none<-numeric(length(x))

mult_id$biom_ins<-numeric(length(x))
mult_id$flo_ins<-numeric(length(x))
mult_id$rec_ins<-numeric(length(x))
mult_id$rt_ins<-numeric(length(x))

mult_id$biom_lit<-numeric(length(x))
mult_id$flo_lit<-numeric(length(x))
mult_id$rec_lit<-numeric(length(x))
mult_id$rt_lit<-numeric(length(x))

mult_id$biom_both<-numeric(length(x))
mult_id$flo_both<-numeric(length(x))
mult_id$rec_both<-numeric(length(x))
mult_id$rt_both<-numeric(length(x))

for(i in 1:nrow(mult_id)){
  mult_id$biom_none[i]<-sum(mult[mult$INTERACT=="NONE",]$biomstd>mult_id$x[i]*0.01)/nrow(mult[mult$INTERACT=="NONE",])
  mult_id$flo_none[i]<-sum(mult[mult$INTERACT=="NONE",]$flostd>mult_id$x[i]*0.01)/nrow(mult[mult$INTERACT=="NONE",])
  mult_id$rec_none[i]<-sum(mult[mult$INTERACT=="NONE",]$recstd>mult_id$x[i]*0.01)/nrow(mult[mult$INTERACT=="NONE",])
  mult_id$rt_none[i]<-sum(mult[mult$INTERACT=="NONE",]$rtstd>mult_id$x[i]*0.01)/nrow(mult[mult$INTERACT=="NONE",])
  
  mult_id$biom_ins[i]<-sum(mult[mult$INTERACT=="INS",]$biomstd>mult_id$x[i]*0.01)/nrow(mult[mult$INTERACT=="INS",])
  mult_id$flo_ins[i]<-sum(mult[mult$INTERACT=="INS",]$flostd>mult_id$x[i]*0.01)/nrow(mult[mult$INTERACT=="INS",])
  mult_id$rec_ins[i]<-sum(mult[mult$INTERACT=="INS",]$recstd>mult_id$x[i]*0.01)/nrow(mult[mult$INTERACT=="INS",])
  mult_id$rt_ins[i]<-sum(mult[mult$INTERACT=="INS",]$rtstd>mult_id$x[i]*0.01)/nrow(mult[mult$INTERACT=="INS",])
  
  mult_id$biom_lit[i]<-sum(mult[mult$INTERACT=="LIT",]$biomstd>mult_id$x[i]*0.01)/nrow(mult[mult$INTERACT=="LIT",])
  mult_id$flo_lit[i]<-sum(mult[mult$INTERACT=="LIT",]$flostd>mult_id$x[i]*0.01)/nrow(mult[mult$INTERACT=="LIT",])
  mult_id$rec_lit[i]<-sum(mult[mult$INTERACT=="LIT",]$recstd>mult_id$x[i]*0.01)/nrow(mult[mult$INTERACT=="LIT",])
  mult_id$rt_lit[i]<-sum(mult[mult$INTERACT=="LIT",]$rtstd>mult_id$x[i]*0.01)/nrow(mult[mult$INTERACT=="LIT",])
  
  mult_id$biom_both[i]<-sum(mult[mult$INTERACT=="BOTH",]$biomstd>mult_id$x[i]*0.01)/nrow(mult[mult$INTERACT=="BOTH",])
  mult_id$flo_both[i]<-sum(mult[mult$INTERACT=="BOTH",]$flostd>mult_id$x[i]*0.01)/nrow(mult[mult$INTERACT=="BOTH",])
  mult_id$rec_both[i]<-sum(mult[mult$INTERACT=="BOTH",]$recstd>mult_id$x[i]*0.01)/nrow(mult[mult$INTERACT=="BOTH",])
  mult_id$rt_both[i]<-sum(mult[mult$INTERACT=="BOTH",]$rtstd>mult_id$x[i]*0.01)/nrow(mult[mult$INTERACT=="BOTH",])
  
}

png("/media/krishna/New Volume/Columbia/StatMod/cedar_creek/multifun_id_prop.png", width=8, height=8, units="in", res=300)
par(mfrow=c(2,2), oma=c(2,0,0,0), mai=c(1, 1, 0.1, 0.1), xpd=F)
plot(x, mult_id$biom_none, type="n", xlab="Threshold", cex.lab=1.6, ylab=NA)
title(ylab="Proportion above threshold", line=2.2, cex.lab=1.3)
lines(x, mult_id$biom_none, lwd=4, col=col[1])
lines(x, mult_id$biom_ins, lwd=4, col=col[2])
lines(x, mult_id$biom_lit, lwd=4, col=col[3])
lines(x, mult_id$biom_both, lwd=4, col=col[4])
legend("topright", "Aboveground biomass", bty="n", cex=1.3)

plot(x, mult_id$biom_none, type="n", xlab="Threshold", cex.lab=1.6, ylab=NA)
title(ylab="Proportion above threshold", line=2.2, cex.lab=1.3)
lines(x, mult_id$flo_none, lwd=4, col=col[1])
lines(x, mult_id$flo_ins, lwd=4, col=col[2])
lines(x, mult_id$flo_lit, lwd=4, col=col[3])
lines(x, mult_id$flo_both, lwd=4, col=col[4])
legend("topright", "Flow time", bty="n", cex=1.3)

plot(x, mult_id$biom_none, type="n", xlab="Threshold", cex.lab=1.6, ylab=NA)
title(ylab="Proportion above threshold", line=2.2, cex.lab=1.3)
lines(x, mult_id$rec_none, lwd=4, col=col[1])
lines(x, mult_id$rec_ins, lwd=4, col=col[2])
lines(x, mult_id$rec_lit, lwd=4, col=col[3])
lines(x, mult_id$rec_both, lwd=4, col=col[4])
legend("topright", "Biomass recovery", bty="n", cex=1.3)

plot(x, mult_id$biom_none, type="n", xlab="Threshold", cex.lab=1.6, ylab=NA)
title(ylab="Proportion above threshold", line=2.2, cex.lab=1.3)
lines(x, mult_id$rt_none, lwd=4, col=col[1])
lines(x, mult_id$rt_ins, lwd=4, col=col[2])
lines(x, mult_id$rt_lit, lwd=4, col=col[3])
lines(x, mult_id$rt_both, lwd=4, col=col[4])
legend("topright", "Belowground biomass", bty="n", cex=1.3)
par(xpd=NA)
legend(x=-150, y=-0.4, bty="n", legend=c("NONE", "INS", "LIT", "BOTH"), col=col, lty=3, lwd=3, horiz = T, cex=1.3)
dev.off()

#Binomial curves---------------------------------------------

png("/media/krishna/New Volume/Columbia/StatMod/cedar_creek/multifun_id_binom.png", width=8, height=8, units="in", res=300)
par(mfrow=c(2,2), oma=c(2,0,0,0), mai=c(1, 1, 0.1, 0.1), xpd=F)
plot(x, mult_id$biom_none, type="n", xlab="Threshold", cex.lab=1.6, ylab=NA)
title(ylab="Proportion above threshold", line=2.2, cex.lab=1.3)
lines(x, mult_id$biom_none, lwd=4, col=col[1])
lines(x, mult_id$biom_ins, lwd=4, col=col[2])
lines(x, mult_id$biom_lit, lwd=4, col=col[3])
lines(x, mult_id$biom_both, lwd=4, col=col[4])
legend("topright", "Aboveground biomass", bty="n", cex=1.3)

plot(x, mult_id$biom_none, type="n", xlab="Threshold", cex.lab=1.6, ylab=NA)
title(ylab="Proportion above threshold", line=2.2, cex.lab=1.3)
lines(x, mult_id$rt_none, lwd=4, col=col[1])
lines(x, mult_id$rt_ins, lwd=4, col=col[2])
lines(x, mult_id$rt_lit, lwd=4, col=col[3])
lines(x, mult_id$rt_both, lwd=4, col=col[4])
legend("topright", "Belowground biomass", bty="n", cex=1.3)

plot(x, mult_id$biom_none, type="n", xlab="Threshold", cex.lab=1.6, ylab=NA)
title(ylab="Proportion above threshold", line=2.2, cex.lab=1.3)
lines(x, mult_id$flo_none, lwd=4, col=col[1])
lines(x, mult_id$flo_ins, lwd=4, col=col[2])
lines(x, mult_id$flo_lit, lwd=4, col=col[3])
lines(x, mult_id$flo_both, lwd=4, col=col[4])
legend("topright", "Water retention", bty="n", cex=1.3)

plot(x, mult_id$biom_none, type="n", xlab="Threshold", cex.lab=1.6, ylab=NA)
title(ylab="Proportion above threshold", line=2.2, cex.lab=1.3)
lines(x, mult_id$rec_none, lwd=4, col=col[1])
lines(x, mult_id$rec_ins, lwd=4, col=col[2])
lines(x, mult_id$rec_lit, lwd=4, col=col[3])
lines(x, mult_id$rec_both, lwd=4, col=col[4])
legend("topright", "Biomass recovery", bty="n", cex=1.3)


par(xpd=NA)
legend(x=-150, y=-0.4, bty="n", legend=c("NONE", "INS", "LIT", "BOTH"), col=col, lty=3, lwd=3, horiz = T, cex=1.3)
dev.off()

#Function for transparent colours------------------------
t_col <- function(color, percent = 50, name = NULL) {
  
  #	  color = color name
  #	percent = % transparency
  #	   name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  
  ## Save the color
  invisible(t.col)
  
}

#Spider plots (Figure S1)-----------------------------------------
#first make mult, the dataframe
library(plyr)
mult_sum<-ddply(mult, .(SPP,INTERACT), summarize, meanbiom=mean(biomstd), meanrt=mean(rtstd), meanflo=mean(flostd), meanrec=mean(recstd))

colnames(mult_sum)<-c("SPP", "INTERACT", "Aboveground biomass", "Belowground biomass", "Water retention", "Biomass recovery")

mult_sum_all<-ddply(mult, .(INTERACT), summarize, meanbiom=mean(biomstd), meanrt=mean(rtstd), meanflo=mean(flostd), meanrec=mean(recstd))

colnames(mult_sum_all)<-c("SPP", "Aboveground biomass", "Belowground biomass", "Water retention", "Biomass recovery")


# Library
library(fmsb)

# Create data: note in High school for several students
set.seed(99)
data=as.data.frame(matrix( sample( 0:20 , 15 , replace=F) , ncol=5))
colnames(data)=c("math" , "english" , "biology" , "music" , "R-coding" )
rownames(data)=paste("mister" , letters[1:3] , sep="-")

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
data=rbind(rep(20,5) , rep(0,5) , data)

mult_sum_1<-mult_sum[mult_sum$SPP==1,3:6]
mult_sum_1<-rbind(rep(1, 4), rep(0, 4), mult_sum_1)

mult_sum_2<-mult_sum[mult_sum$SPP==2,3:6]
mult_sum_2<-rbind(rep(1, 4), rep(0, 4), mult_sum_2)

mult_sum_4<-mult_sum[mult_sum$SPP==4,3:6]
mult_sum_4<-rbind(rep(1, 4), rep(0, 4), mult_sum_4)

mult_sum_8<-mult_sum[mult_sum$SPP==8,3:6]
mult_sum_8<-rbind(rep(1, 4), rep(0, 4), mult_sum_8)

mult_sum_16<-mult_sum[mult_sum$SPP==16,3:6]
mult_sum_16<-rbind(rep(1, 4), rep(0, 4), mult_sum_16)

mult_sum_all<-rbind(rep(1, 5), rep(0, 5), mult_sum_all)

colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )

colors<-wes_palette("Royal1", 4, type="continuous")
colors_border<-c(colors[4], colors[2], colors[3], colors[1])

colors_in<-c(t_col(colors_border[1], 40), t_col(colors_border[2], 40), t_col(colors_border[3], 40),
             t_col(colors_border[4], 40))


par(mar=c(1, 1, 1, 1))
radarchart( mult_sum_all[,2:5]  , axistype=1 ,
            pcol=colors_border , 
            plwd=4 , plty=1,
            cglcol="grey", cglty=1, axislabcol="grey45", caxislabels=seq(0,1,5), cglwd=0.8,
            vlabels=NA
)

xpos<-c(0, -1.05, 0, 1.05)
ypos<-c(1.15, 0, -1.15, 0)
labs<-c("Aboveground\nbiomass", "Belowground\nbiomass", 
        "Water\nretention", "Biomass\nrecovery")
text(xpos, ypos, labs, cex=1.3)

?text
colnames(mult_sum_all)
?radarchart

png("spider_plot_nofill_2.png", width=8, height=12, units="in", res=300)
par(mfrow=c(3,2), mar=c(1,1,1,1))
xpos<-c(0, -1, 0, 1.05)
ypos<-c(1.15, 0, -1.15, 0)
labs<-c("Aboveground\nbiomass", "Belowground\nbiomass", 
        "Water\nretention", "Biomass\nrecovery")

radarchart( mult_sum_all[,2:5]  , axistype=1 ,
            pcol=colors_border , 
            plwd=4 , plty=1,
            cglcol="grey", cglty=1, axislabcol="grey45", caxislabels=seq(0,1,5), cglwd=0.8,
            vlabels=NA
)
legend("topleft", "a) All", cex=1.3, bty="n")
text(xpos, ypos, labs, cex=1.3)

radarchart( mult_sum_1  , axistype=1 ,
            pcol=colors_border , 
            plwd=4 , plty=1,
            cglcol="grey", cglty=1, axislabcol="grey45", caxislabels=seq(0,1,5), cglwd=0.8,
            vlabels=NA
)
legend("topleft", "b) #sp = 1", cex=1.3, bty="n")
text(xpos, ypos, labs, cex=1.3)

radarchart( mult_sum_2  , axistype=1 ,
            pcol=colors_border , 
            plwd=4 , plty=1,
            cglcol="grey", cglty=1, axislabcol="grey45", caxislabels=seq(0,1,5), cglwd=0.8,
            vlabels=NA
)
legend("topleft", "c) #sp = 2", cex=1.3, bty="n")
text(xpos, ypos, labs, cex=1.3)

radarchart( mult_sum_4  , axistype=1 ,
            pcol=colors_border , 
            plwd=4 , plty=1,
            cglcol="grey", cglty=1, axislabcol="grey45", caxislabels=seq(0,1,5), cglwd=0.8,
            vlabels=NA
)
legend("topleft", "d) #sp = 4", cex=1.3, bty="n")
text(xpos, ypos, labs, cex=1.3)

radarchart( mult_sum_8  , axistype=1 ,
            pcol=colors_border , 
            plwd=4 , plty=1,
            cglcol="grey", cglty=1, axislabcol="grey45", caxislabels=seq(0,1,5), cglwd=0.8,
            vlabels=NA
)
legend("topleft", "e) #sp = 8", cex=1.3, bty="n")
text(xpos, ypos, labs, cex=1.3)

radarchart( mult_sum_16  , axistype=1 ,
            pcol=colors_border , 
            plwd=4 , plty=1,
            cglcol="grey", cglty=1, axislabcol="grey45", caxislabels=seq(0,1,5), cglwd=0.8,
            vlabels=NA
)
legend("topleft", "f) #sp = 16", cex=1.3, bty="n")
text(xpos, ypos, labs, cex=1.3)


par(xpd=NA)
legend(x=-1.7, y=5.2, bty="n", legend=c("NONE", "INS", "LIT", "BOTH"), col=colors, lty=1, lwd=3, horiz = F, cex=1.3)
dev.off()
