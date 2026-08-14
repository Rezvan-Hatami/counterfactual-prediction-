# 20_counterfactual_analysis.R
# Author: Rezvan Hatami
# Date: 2026-08-01
# This script performs the counterfactual analysis for discharge-related variables and plots the resulting response profiles.

rm(list = ls())

# ----------------------------------------------------------------------
# 1. Package setup
# ----------------------------------------------------------------------
required_pkgs <- c(
  "fs",
  "here",
  "vegan",
  "car",
  "lattice",
  "ecodist",
  "BiodiversityR",
  "latticeExtra"
)

missing_required_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_required_pkgs) > 0) {
  stop(
    "These packages must be installed before running the script: ",
    paste(missing_required_pkgs, collapse = ", ")
  )
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

options(stringsAsFactors = FALSE)
options(scipen = 1)

# ----------------------------------------------------------------------
# 2. Paths
# ----------------------------------------------------------------------
project_dir <- here::here()
output_dir <- fs::path(project_dir, "output")
figures_dir <- fs::path(project_dir, "figures")
script_figure_dir <- fs::path(figures_dir, "20")

fs::dir_create(figures_dir)
fs::dir_create(script_figure_dir)

wangbug_file <- fs::path(output_dir, "01_wangbug_raw.rds")
wangenv_file <- fs::path(output_dir, "01_wangenv_prepped.rds")
script3_inputs_file <- fs::path(output_dir, "03_environmental_script_inputs.rds")

if (!fs::file_exists(wangbug_file)) {
  stop("Missing input file: ", wangbug_file)
}
if (!fs::file_exists(wangenv_file)) {
  stop("Missing input file: ", wangenv_file)
}
if (!fs::file_exists(script3_inputs_file)) {
  stop("Missing input file: ", script3_inputs_file)
}

# ----------------------------------------------------------------------
# 3. Read inputs
# ----------------------------------------------------------------------
wangbug <- readRDS(wangbug_file)
wangenv <- readRDS(wangenv_file)
script3_inputs <- readRDS(script3_inputs_file)

alldata <- script3_inputs$alldata

# ----------------------------------------------------------------------
# 4. Counterfactual analysis
# ----------------------------------------------------------------------# ----------------------------------------------------------------------
######””””””””””#######

###””””” PCO (Principal coordinate analysis)””””” #######
wangbug.BC<-vegdist(sqrt(wangbug)) #Bray-Curtis dissimilarity matrix of square-root transformed abundances
n<-dim(wangbug)[1]
p<-n-1
#cmdscale is multidimensional scaling, also known as principal coordinates analysis (Gower, 1966)
wangbug.mds<-cmdscale(wangbug.BC, k = p, eig = TRUE, add = TRUE, x.ret = FALSE)# PCoA using the Bray-Curtis dissimilarity measure on square-root transformed abundances, and a correction for negative eigenvalues 
# add=TRUE, Logical indicating if an additive constant should be computed, and added to the non-diagonal dissimilarities such that all eigenvalues are non-negative in the underlying Principal Co-ordinates Analysis (see cmdscale for details)
pco.varpercent<-round(wangbug.mds$eig/sum(wangbug.mds$eig)*100,digits=1) # Percentage of variation explained by each successive PCO axis
round(cumsum(wangbug.mds$eig)/sum(wangbug.mds$eig)*100,digits=2)
# #view(wangbug.mds$points)
pwangbug<-wangbug.mds$points
dim(pwangbug)
colnames(pwangbug)<-c(paste("pco",sep="",1:p))
#view(pwangbug)

bugenv<-cbind(wangenv,pwangbug)# combine environmental data and PCO matrix
##view(bugenv)
dim(bugenv)
##view(pwangbug)

wangenv$eff<-ifelse(wangenv$dist<4,0,1)
wangenv$time<-as.factor(wangenv$time)
alldata<-cbind(wangenv,pwangbug)#

#””#*****#””#*****#””#*****counterfactual analysis for conductivity discharge load#””#

tiff(file = fs::path(script_figure_dir, "20_figure_2_conductivity.jpg"), width = 8, height = 10, units = "in", pointsize = 12, bg = "transparent", res = 800, compression = "lzw")
par(mfrow=c(4,1), mar=c(5,4.5,2.5,3),mai = c(0.6, 1.2, 0.3, 0.5),cex=1.5,cex.axis=0.9,las=1,cex.main=1,cex.lab=0.8)

plot(bugenv[bugenv$day==1,"dist"],bugenv[bugenv$day==1,"cond"],type = "b",ylim=c(0,800),pch=19,xlab="",ylab="Conductivity (mg/L)",main="",lwd=3,lty=1)
points(bugenv[bugenv$day==126,"dist"],bugenv[bugenv$day==126,"cond"],type="b",col="gold2",pch=15,lwd=3,cex=0.9,lty=2)
points(bugenv[bugenv$day==260,"dist"],bugenv[bugenv$day==260,"cond"],type="b",col="blue",pch=8,cex=0.9,lwd=3,lty=3)
points(bugenv[bugenv$day==336,"dist"],bugenv[bugenv$day==336,"cond"],type="b",col="green3",pch=17,cex=0.9,lwd=3,lty=4)
points(bugenv[bugenv$day==518,"dist"],bugenv[bugenv$day==518,"cond"],type="b",col="red",pch=18,cex=0.9,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.5,title="months")
abline(v=4,lty=2)

#by dividing all the discharge values in half
## First, alkalinity value is needed for TP model (calculating again considering missing data)
alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dalk<-alldata8$dalk*0.5
alkpred<-exp(predict(tempo,newdata)+residuals(tempo))
envdata<-(cbind(alkpred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(tp~alk+ph2+rain2+rain3+eff:I(dtp/dayflow):time+eff:I(dtp/dayflow):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dtp<-alldata8$dtp*0.5
newdata$alk<-alkpred
tppred<-predict(tempo,newdata)+residuals(tempo)
envdata1<-(cbind(tppred,envdata))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(cond)~tp+alk+time+dist:time:eff,data=alldata8)
newdata<-alldata8
newdata<-newdata
newdata$tp<-tppred
newdata$alk<-alkpred
condcounter<-exp(predict(tempo,newdata)+residuals(tempo))
envdata<-(cbind(condcounter,alldata8))
plot(envdata[envdata$day==126,"dist"],envdata[envdata$day==126,"condcounter"],type = "l",ylim=c(0,800),pch=15,col="gold2",xlab="",ylab="Conductivity (mg/L)",main="",lwd=3,lty=2)
#points(envdata[envdata$day==126,"dist"],envdata[envdata$day==126,"condcounter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(envdata[envdata$day==260,"dist"],envdata[envdata$day==260,"condcounter"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(envdata[envdata$day==336,"dist"],envdata[envdata$day==336,"condcounter"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(envdata[envdata$day==518,"dist"],envdata[envdata$day==518,"condcounter"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
abline(v=4,lty=2)

#by dividing all the discharge values in 10
alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dalk<-alldata8$dalk*0.1
alkpred<-exp(predict(tempo,newdata)+residuals(tempo))
envdata<-(cbind(alkpred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(tp~alk+ph2+rain2+rain3+eff:I(dtp/dayflow):time+eff:I(dtp/dayflow):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dtp<-alldata8$dtp*0.1
newdata$alk<-alkpred
tppred<-predict(tempo,newdata)+residuals(tempo)
envdata1<-(cbind(tppred,envdata))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(cond)~tp+alk+time+dist:time:eff,data=alldata8)
newdata<-alldata8
newdata<-newdata
newdata$tp<-tppred
newdata$alk<-alkpred
condcounter<-exp(predict(tempo,newdata)+residuals(tempo))
envdata<-(cbind(condcounter,alldata8))
plot(envdata[envdata$day==126,"dist"],envdata[envdata$day==126,"condcounter"],type = "l",ylim=c(0,800),pch=15,col="gold2",xlab="",ylab="Conductivity (mg/L)",main="",lwd=3,lty=2)
#points(envdata[envdata$day==126,"dist"],envdata[envdata$day==126,"condcounter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(envdata[envdata$day==260,"dist"],envdata[envdata$day==260,"condcounter"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(envdata[envdata$day==336,"dist"],envdata[envdata$day==336,"condcounter"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(envdata[envdata$day==518,"dist"],envdata[envdata$day==518,"condcounter"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
abline(v=4,lty=2)

#by dividing all the discharge values in 100
alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dalk<-alldata8$dalk*0.01
alkpred<-exp(predict(tempo,newdata)+residuals(tempo))
envdata<-(cbind(alkpred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(tp~alk+ph2+rain2+rain3+eff:I(dtp/dayflow):time+eff:I(dtp/dayflow):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dtp<-alldata8$dtp*0.01
newdata$alk<-alkpred
tppred<-predict(tempo,newdata)+residuals(tempo)
envdata1<-(cbind(tppred,envdata))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(cond)~tp+alk+time+dist:time:eff,data=alldata8)
newdata<-alldata8
newdata<-newdata
newdata$tp<-tppred
newdata$alk<-alkpred
condcounter<-exp(predict(tempo,newdata)+residuals(tempo))
envdata<-(cbind(condcounter,alldata8))
plot(envdata[envdata$day==126,"dist"],envdata[envdata$day==126,"condcounter"],type = "l",ylim=c(0,800),pch=15,col="gold2",xlab="Distance(km)",ylab="Conductivity (mg/L)",main="",lwd=3,lty=2)
#points(envdata[envdata$day==126,"dist"],envdata[envdata$day==126,"condcounter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(envdata[envdata$day==260,"dist"],envdata[envdata$day==260,"condcounter"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(envdata[envdata$day==336,"dist"],envdata[envdata$day==336,"condcounter"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(envdata[envdata$day==518,"dist"],envdata[envdata$day==518,"condcounter"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
abline(v=4,lty=2)
dev.off()

#***** Figure 3 *****#””#counterfactual analysis for TOC discharge load

tiff(file = fs::path(script_figure_dir, "20_figure_3_toc.jpg"), width = 8, height = 10, units = "in", pointsize = 12, bg = "transparent", res = 800, compression = "lzw")
par(mfrow=c(4,1), mar=c(5,4.5,2.5,3),mai = c(0.6, 1.2, 0.3, 0.5),cex=1.5,cex.axis=0.9,las=1,cex.main=1,cex.lab=0.8)

#toc plotted against spatial position
plot(bugenv[bugenv$day==1,"dist"],bugenv[bugenv$day==1,"toc"],type = "b",ylim=c(0,10),pch=19,cex=0.9,xlab="",ylab="TOC (mg/L)",main="",lwd=3,lty=1)
points(bugenv[bugenv$day==126,"dist"],bugenv[bugenv$day==126,"toc"],type="b",col="gold2",pch=15,lwd=3,cex=1,lty=2)
points(bugenv[bugenv$day==260,"dist"],bugenv[bugenv$day==260,"toc"],type="b",col="blue",pch=8,cex=0.8,lwd=3,lty=3)
points(bugenv[bugenv$day==336,"dist"],bugenv[bugenv$day==336,"toc"],type="b",col="green3",pch=17,cex=1,lwd=3,lty=4)
points(bugenv[bugenv$day==518,"dist"],bugenv[bugenv$day==518,"toc"],type="b",col="red",pch=18,cex=1.2,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.5,title="months")
abline(v=4,lty=2)

# what if toc had been half
tempo<-lm(toc~dayflow+rain2+rain3+eff:I((dtoc/dayflow)):time+eff:I((dtoc/dayflow)):I(dist==4.08):time,data=alldata)# Checked RH5
anova(tempo)
newdata<-alldata
newdata$dtoc<-alldata$dtoc*0.5
toccounter<-predict(tempo,newdata)+residuals(tempo)
wangenv1<-cbind(toccounter,alldata)
plot(wangenv1[wangenv1$day==1,"dist"],wangenv1[wangenv1$day==1,"toccounter"],type = "l",ylim=c(0,10),pch=19,col="black",xlab="",ylab="TOC (mg/L)",main="",lwd=3,lty=1)
points(wangenv1[wangenv1$day==126,"dist"],wangenv1[wangenv1$day==126,"toccounter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(wangenv1[wangenv1$day==260,"dist"],wangenv1[wangenv1$day==260,"toccounter"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(wangenv1[wangenv1$day==336,"dist"],wangenv1[wangenv1$day==336,"toccounter"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(wangenv1[wangenv1$day==518,"dist"],wangenv1[wangenv1$day==518,"toccounter"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.5,title="months")
abline(v=4,lty=2)

## if toc had been 1/10
tempo<-lm((toc)~dayflow+rain2+rain3+eff:I((dtoc/dayflow)):time+eff:I((dtoc/dayflow)):I(dist==4.08):time,data=alldata)# Checked RH5
anova(tempo)
newdata<-alldata
newdata$dtoc<-alldata$dtoc*0.1
toccounter<-predict(tempo,newdata=newdata)+residuals(tempo)
wangenv1<-(cbind(toccounter,alldata))
plot(wangenv1[wangenv1$day==1,"dist"],wangenv1[wangenv1$day==1,"toccounter"],type = "l",ylim=c(0,10),pch=19,col="black",xlab="",ylab="TOC (mg/L)",main="",lwd=3,lty=1)
points(wangenv1[wangenv1$day==126,"dist"],wangenv1[wangenv1$day==126,"toccounter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(wangenv1[wangenv1$day==260,"dist"],wangenv1[wangenv1$day==260,"toccounter"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(wangenv1[wangenv1$day==336,"dist"],wangenv1[wangenv1$day==336,"toccounter"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(wangenv1[wangenv1$day==518,"dist"],wangenv1[wangenv1$day==518,"toccounter"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.5,title="months")
abline(v=4,lty=2)

## if toc had been 1/100
tempo<-lm((toc)~dayflow+rain2+rain3+eff:I((dtoc/dayflow)):time+eff:I((dtoc/dayflow)):I(dist==4.08):time,data=alldata)# Checked RH5
newdata<-alldata
newdata$dtoc<-alldata$dtoc*0.01
toccounter<-predict(tempo,newdata)+residuals(tempo)
wangenv1<-(cbind(toccounter,alldata))
plot(wangenv1[wangenv1$day==1,"dist"],wangenv1[wangenv1$day==1,"toccounter"],type = "l",ylim=c(0,10),pch=19,col="black", xlab="",ylab="TOC (mg/L)",main="",lwd=3,lty=1)
points(wangenv1[wangenv1$day==126,"dist"],wangenv1[wangenv1$day==126,"toccounter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(wangenv1[wangenv1$day==260,"dist"],wangenv1[wangenv1$day==260,"toccounter"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(wangenv1[wangenv1$day==336,"dist"],wangenv1[wangenv1$day==336,"toccounter"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(wangenv1[wangenv1$day==518,"dist"],wangenv1[wangenv1$day==518,"toccounter"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.5,title="months")
abline(v=4,lty=2)
dev.off()


#***** Figure 4 *****###””##counterfactual analysis for chlorophyll a discharge load 

tiff(file = fs::path(script_figure_dir, "20_figure_4_chlorophyll_a.jpg"), width = 8, height = 10, units = "in", pointsize = 12, bg = "transparent", res = 800, compression = "lzw")
par(mfrow=c(4,1), mar=c(5,4.5,2.5,3),mai = c(0.6, 1.2, 0.3, 0.5),cex=1.5,cex.axis=0.9,las=1,cex.main=1,cex.lab=0.8)

#Chlorophyll A plotted against spatial position
plot(bugenv[bugenv$day==1,"dist"],bugenv[bugenv$day==1,"chla"],type = "b",ylim=c(0,40),cex=0.9, pch=19,xlab="",ylab="Chlorophyll A (mg/L)",main="",lwd=3,lty=1)
points(bugenv[bugenv$day==126,"dist"],bugenv[bugenv$day==126,"chla"],type="b",col="gold2",pch=15,lwd=3,cex=1,lty=2)
points(bugenv[bugenv$day==260,"dist"],bugenv[bugenv$day==260,"chla"],type="b",col="blue",pch=8,cex=0.8,lwd=3,lty=3)
points(bugenv[bugenv$day==336,"dist"],bugenv[bugenv$day==336,"chla"],type="b",col="green3",pch=17,cex=1,lwd=3,lty=4)
points(bugenv[bugenv$day==518,"dist"],bugenv[bugenv$day==518,"chla"],type="b",col="red",pch=18,cex=1.2,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.5,title="months")
abline(v=4,lty=2)

### dchla, dtp, dno3, dalk divided by half 
## First, TOC and alkalinity value are needed for TP model
alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dalk<-alldata8$dalk*0.5
alkpred<-exp(predict(tempo,newdata)+residuals(tempo))
envdata<-(cbind(alkpred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(no3)~temp+rain2+rain3+eff:I(log(dno3)/dayflow):time+eff:I(log(dno3)/dayflow):I(dist==4.08):time,data=alldata8)##RH5
newdata<-alldata8
newdata$dno3<-alldata8$dno3*0.5
no3pred<-exp(predict(tempo,newdata)+residuals(tempo))
envdata<-(cbind(no3pred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(tp~alk+ph2+rain2+rain3+eff:I(dtp/dayflow):time+eff:I(dtp/dayflow):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dtp<-alldata8$dtp*0.5
newdata$alk<-alkpred
tppred<-predict(tempo,newdata)+residuals(tempo)
envdata1<-(cbind(tppred,envdata))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(chla)~solar+temp+tp+no3+turb+solar:temp+eff:log(I(dchla/dayflow)),data=alldata8)# Checked RH
newdata<-alldata8
newdata$dchla<-alldata8$dchla*0.5
newdata$no3<-no3pred
newdata$tp<-tppred
chlcounter<-exp(predict(tempo,newdata)+residuals(tempo))
envdata2<-(cbind(chlcounter,envdata1))
##view(envdata2)
plot(envdata2[envdata2$day==126,"dist"],envdata2[envdata2$day==126,"chlcounter"],type = "l",ylim=c(0,40),pch=15,col="gold2",xlab="",ylab="Chlorophyll A (mg/L)",main="",lwd=3,lty=2)
#points(envdata2[envdata2$day==126,"dist"],envdata2[envdata2$day==126,"chlcounter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(envdata2[envdata2$day==260,"dist"],envdata2[envdata2$day==260,"chlcounter"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(envdata2[envdata2$day==336,"dist"],envdata2[envdata2$day==336,"chlcounter"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(envdata2[envdata2$day==518,"dist"],envdata2[envdata2$day==518,"chlcounter"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
abline(v=4,lty=2)

### dchla, dtp, dno3, dalk divided by 10
## First, TOC and alkalinity value are needed for TP model
alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dalk<-alldata8$dalk*0.1
alkpred<-exp(predict(tempo,newdata)+residuals(tempo))
envdata<-(cbind(alkpred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(no3)~temp+rain2+rain3+eff:I(log(dno3)/dayflow):time+eff:I(log(dno3)/dayflow):I(dist==4.08):time,data=alldata8)##RH5
newdata<-alldata8
newdata$dno3<-alldata8$dno3*0.1
no3pred<-exp(predict(tempo,newdata)+residuals(tempo))
envdata<-(cbind(no3pred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(tp~alk+ph2+rain2+rain3+eff:I(dtp/dayflow):time+eff:I(dtp/dayflow):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dtp<-alldata8$dtp*0.1
newdata$alk<-alkpred
tppred<-predict(tempo,newdata)+residuals(tempo)
envdata1<-(cbind(tppred,envdata))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(chla)~solar+temp+tp+no3+turb+solar:temp+eff:log(I(dchla/dayflow)),data=alldata8)# Checked RH
newdata<-alldata8
newdata$dchla<-alldata8$dchla*0.1
newdata$no3<-no3pred
newdata$tp<-tppred
chlcounter<-exp(predict(tempo,newdata)+residuals(tempo))
envdata2<-(cbind(chlcounter,envdata1))
##view(envdata2)
plot(envdata2[envdata2$day==126,"dist"],envdata2[envdata2$day==126,"chlcounter"],type = "l",ylim=c(0,40),pch=15,col="gold2",xlab="",ylab="Chlorophyll A (mg/L)",lwd=3,lty=2)
#points(envdata2[envdata2$day==126,"dist"],envdata2[envdata2$day==126,"chlcounter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(envdata2[envdata2$day==260,"dist"],envdata2[envdata2$day==260,"chlcounter"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(envdata2[envdata2$day==336,"dist"],envdata2[envdata2$day==336,"chlcounter"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(envdata2[envdata2$day==518,"dist"],envdata2[envdata2$day==518,"chlcounter"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
abline(v=4,lty=2)

### dchla, dtp, dno3, dalk divided by 100
## First, TOC and alkalinity value are needed for TP model
alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dalk<-alldata8$dalk*0.01
alkpred<-exp(predict(tempo,newdata)+residuals(tempo))
envdata<-(cbind(alkpred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(no3)~temp+rain2+rain3+eff:I(log(dno3)/dayflow):time+eff:I(log(dno3)/dayflow):I(dist==4.08):time,data=alldata8)##RH5
newdata<-alldata8
newdata$dno3<-alldata8$dno3*0.01
no3pred<-exp(predict(tempo,newdata)+residuals(tempo))
envdata<-(cbind(no3pred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(tp~alk+ph2+rain2+rain3+eff:I(dtp/dayflow):time+eff:I(dtp/dayflow):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dtp<-alldata8$dtp*0.01
newdata$alk<-alkpred
tppred<-predict(tempo,newdata)+residuals(tempo)
envdata1<-(cbind(tppred,envdata))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(chla)~solar+temp+tp+no3+turb+solar:temp+eff:log(I(dchla/dayflow)),data=alldata8)# Checked RH
newdata<-alldata8
newdata$dchla<-alldata8$dchla*0.01
newdata$no3<-no3pred
newdata$tp<-tppred
chlcounter<-exp(predict(tempo,newdata)+residuals(tempo))
envdata2<-(cbind(chlcounter,envdata1))
##view(envdata2)
plot(envdata2[envdata2$day==126,"dist"],envdata2[envdata2$day==126,"chlcounter"],type = "l",cex=0.9,ylim=c(0,40),pch=15,col="gold2",xlab="",ylab="Chlorophyll A (mg/L)",main="",lwd=3,lty=2)
#points(envdata2[envdata2$day==126,"dist"],envdata2[envdata2$day==126,"chlcounter"],type="l",col="gold2",pch=15,lwd=3,cex=0.9,lty=2)
points(envdata2[envdata2$day==260,"dist"],envdata2[envdata2$day==260,"chlcounter"],type="l",col="blue",pch=8,cex=0.9,lwd=3,lty=3)
points(envdata2[envdata2$day==336,"dist"],envdata2[envdata2$day==336,"chlcounter"],type="l",col="green3",pch=17,cex=0.9,lwd=3,lty=4)
points(envdata2[envdata2$day==518,"dist"],envdata2[envdata2$day==518,"chlcounter"],type="l",col="red",pch=18,cex=0.9,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
abline(v=4,lty=2)
dev.off()


#””#*****#””#*****#””#*****manipulating chlorophyll a discharge load

tiff(file = fs::path(script_figure_dir, "20_figure_5a_chlorophyll_a_isolated_and_joint_interventions.jpg"), width = 8, height = 10, units = "in", pointsize = 12, bg = "transparent", res = 800, compression = "lzw")
par(mfrow=c(4,1), mar=c(5,4.5,2.5,3),mai = c(0.6, 1.2, 0.3, 0.5),cex=1.5,cex.axis=0.9,las=1,cex.main=1,cex.lab=0.8)

#prediction with environmental variables
tempo<-lm(log(chla)~solar+temp+tp+no3+turb+solar:temp+eff:log(I(dchla/dayflow)),data=alldata[complete.cases(alldata[,c("temp","no3","tp")]),])# 
summary(tempo)
anova(tempo)
tempo1<-cbind(cbind(exp(fitted(tempo)),exp(residuals(tempo))),alldata[complete.cases(alldata[,c("temp","no3","tp")]),])
names(tempo1[,1:2])<-c("fitted","residual")
tempo1[,"2"]
plot(tempo1[tempo1$day==126,"dist"],tempo1[tempo1$day==126,"1"],type = "l",ylim=c(0,40),pch=15,col="gold2",xlab="",ylab="Chlorophyll A (mg/L)",main="",lwd=3,lty=2)
#points(tempo1[tempo1$day==126,"dist"],tempo1[tempo1$day==126,"1"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(tempo1[tempo1$day==260,"dist"],tempo1[tempo1$day==260,"1"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(tempo1[tempo1$day==336,"dist"],tempo1[tempo1$day==336,"1"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(tempo1[tempo1$day==518,"dist"],tempo1[tempo1$day==518,"1"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
abline(v=4,lty=2)

#########chla in discharge divided by half
tempo<-lm(log(chla)~solar+temp+tp+no3+turb+solar:temp+eff:log(I(dchla/dayflow)),data=alldata[complete.cases(alldata[,c("temp","no3","tp")]),])# 
##view(alldata)
newdata<-alldata
newdata$dchla<-alldata$dchla*0.5
chlpred<-exp(predict(tempo,newdata=newdata))
envdata<-(cbind(chlpred,alldata))
plot(envdata[envdata$day==126,"dist"],envdata[envdata$day==126,"chlpred"],type = "l",ylim=c(0,40),pch=15,col="gold2",xlab="",ylab="Chlorophyll A (mg/L)",main="",lwd=3,lty=2)
#points(envdata[envdata$day==126,"dist"],envdata[envdata$day==126,"chlpred"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(envdata[envdata$day==260,"dist"],envdata[envdata$day==260,"chlpred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(envdata[envdata$day==336,"dist"],envdata[envdata$day==336,"chlpred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(envdata[envdata$day==518,"dist"],envdata[envdata$day==518,"chlpred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
abline(v=4,lty=2)

#########chla in discharge divided by 10
tempo<-lm(log(chla)~solar+temp+tp+no3+turb+solar:temp+eff:log(I(dchla/dayflow)),data=alldata[complete.cases(alldata[,c("no3","tp")]),])# 
##view(alldata)
newdata<-alldata
newdata$dchla<-alldata$dchla*0.1
chlpred<-exp(predict(tempo,newdata=newdata))
envdata<-(cbind(chlpred,alldata))
plot(envdata[envdata$day==126,"dist"],envdata[envdata$day==126,"chlpred"],type = "l",ylim=c(0,40),pch=15,col="gold2",xlab="",ylab="Chlorophyll A (mg/L)",main="",lwd=3,lty=2)
#points(envdata[envdata$day==126,"dist"],envdata[envdata$day==126,"chlpred"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(envdata[envdata$day==260,"dist"],envdata[envdata$day==260,"chlpred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(envdata[envdata$day==336,"dist"],envdata[envdata$day==336,"chlpred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(envdata[envdata$day==518,"dist"],envdata[envdata$day==518,"chlpred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
abline(v=4,lty=2)

#########chla in discharge divided by 100
tempo<-lm(log(chla)~solar+temp+tp+no3+turb+solar:temp+eff:log(I(dchla/dayflow)),data=alldata[complete.cases(alldata[,c("no3","tp")]),])# 
##view(alldata)
newdata<-alldata
newdata$dchla<-alldata$dchla*0.01
chlpred<-exp(predict(tempo,newdata=newdata))
envdata<-(cbind(chlpred,alldata))
plot(envdata[envdata$day==126,"dist"],envdata[envdata$day==126,"chlpred"],type = "l",ylim=c(0,40),pch=15,col="gold2",xlab="Distance(km)",ylab="Chlorophyll A (mg/L)",main="",lwd=3,lty=2)
#points(envdata[envdata$day==126,"dist"],envdata[envdata$day==126,"chlpred"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(envdata[envdata$day==260,"dist"],envdata[envdata$day==260,"chlpred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(envdata[envdata$day==336,"dist"],envdata[envdata$day==336,"chlpred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(envdata[envdata$day==518,"dist"],envdata[envdata$day==518,"chlpred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
abline(v=4,lty=2)

dev.off()

##*****Figure 5b *****##

tiff(file = fs::path(script_figure_dir, "20_figure_5b_chlorophyll_a_isolated_and_joint_interventions.jpg"), width = 8, height = 10, units = "in", pointsize = 12, bg = "transparent", res = 800, compression = "lzw")
par(mfrow=c(4,1), mar=c(5,4.5,2.5,3),mai = c(0.6, 1.2, 0.3, 0.5),cex=1.5,cex.axis=0.9,las=1,cex.main=1,cex.lab=0.8)

#prediction with environmental variables
tempo<-lm(log(chla)~solar+temp+tp+no3+turb+solar:temp+eff:log(I(dchla/dayflow)),data=alldata[complete.cases(alldata[,c("temp","no3","tp")]),])# Checked RH
summary(tempo)
anova(tempo)
tempo1<-cbind(cbind(exp(fitted(tempo)),exp(residuals(tempo))),alldata[complete.cases(alldata[,c("temp","no3","tp")]),])
names(tempo1[,1:2])<-c("fitted","residual")
tempo1[,"2"]
plot(tempo1[tempo1$day==126,"dist"],tempo1[tempo1$day==126,"1"],type = "l",ylim=c(0,40),pch=15,col="gold2",xlab="",ylab="Chlorophyll A (mg/L)",main="",lwd=3,lty=2)
#points(tempo1[tempo1$day==126,"dist"],tempo1[tempo1$day==126,"1"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(tempo1[tempo1$day==260,"dist"],tempo1[tempo1$day==260,"1"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(tempo1[tempo1$day==336,"dist"],tempo1[tempo1$day==336,"1"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(tempo1[tempo1$day==518,"dist"],tempo1[tempo1$day==518,"1"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
abline(v=4,lty=2)

### dchla, dtp, dno3, dalk divided by half 
## First, TOC and alkalinity value are needed for TP model
alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata8[complete.cases(alldata[,c("tp")]),]) 
newdata<-alldata8
newdata$dalk<-alldata8$dalk*0.5
alkpred<-exp(predict(tempo,newdata=newdata))
envdata<-(cbind(alkpred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(no3)~temp+rain2+rain3+eff:I(log(dno3)/dayflow):time+eff:I(log(dno3)/dayflow):I(dist==4.08):time,data=alldata)##RH5
newdata<-alldata8
newdata$dno3<-alldata8$dno3*0.5
no3pred<-exp(predict(tempo,newdata=newdata))
envdata<-(cbind(no3pred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(tp~alk+ph2+rain2+rain3+eff:I(dtp/dayflow):time+eff:I(dtp/dayflow):I(dist==4.08):time,data=alldata8[complete.cases(alldata[,c("tp")]),])
newdata<-alldata8
newdata$dtp<-alldata8$dtp*0.5
newdata$alk<-alkpred
tppred<-(predict(tempo,newdata=newdata))
envdata1<-(cbind(tppred,envdata))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(chla)~solar+temp+tp+no3+turb+solar:temp+eff:log(I(dchla/dayflow)),data=alldata8)# Checked RH
newdata<-alldata8
newdata$dchla<-alldata8$dchla*0.5
newdata$no3<-no3pred
newdata$tp<-tppred
chlpred<-exp(predict(tempo,newdata=newdata))
envdata2<-(cbind(chlpred,envdata1))
##view(envdata2)
plot(envdata2[envdata2$day==126,"dist"],envdata2[envdata2$day==126,"chlpred"],type = "l",ylim=c(0,40),pch=15,col="gold2",xlab="",ylab="Chlorophyll A (mg/L)",main="",lwd=3,lty=2)
#points(envdata2[envdata2$day==126,"dist"],envdata2[envdata2$day==126,"chlpred"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(envdata2[envdata2$day==260,"dist"],envdata2[envdata2$day==260,"chlpred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(envdata2[envdata2$day==336,"dist"],envdata2[envdata2$day==336,"chlpred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(envdata2[envdata2$day==518,"dist"],envdata2[envdata2$day==518,"chlpred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
abline(v=4,lty=2)

### dchla, dtp, dno3, dalk divided by 10
alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata8[complete.cases(alldata[,c("tp")]),]) 
newdata<-alldata8
newdata$dalk<-alldata8$dalk*0.1
alkpred<-exp(predict(tempo,newdata=newdata))
envdata<-(cbind(alkpred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(no3)~temp+rain2+rain3+eff:I(log(dno3)/dayflow):time+eff:I(log(dno3)/dayflow):I(dist==4.08):time,data=alldata)##RH5
newdata<-alldata8
newdata$dno3<-alldata8$dno3*0.1
no3pred<-exp(predict(tempo,newdata=newdata))
envdata<-(cbind(no3pred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(tp~alk+ph2+rain2+rain3+eff:I(dtp/dayflow):time+eff:I(dtp/dayflow):I(dist==4.08):time,data=alldata8[complete.cases(alldata[,c("tp")]),])
newdata<-alldata8
newdata$dtp<-alldata8$dtp*0.1
newdata$alk<-alkpred
tppred<-(predict(tempo,newdata=newdata))
envdata1<-(cbind(tppred,envdata))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(chla)~solar+temp+tp+no3+turb+solar:temp+eff:log(I(dchla/dayflow)),data=alldata8)# Checked RH
newdata<-alldata8
newdata$dchla<-alldata8$dchla*0.1
newdata$no3<-no3pred
newdata$tp<-tppred
chlpred<-exp(predict(tempo,newdata=newdata))
envdata2<-(cbind(chlpred,envdata1))
##view(envdata2)
plot(envdata2[envdata2$day==126,"dist"],envdata2[envdata2$day==126,"chlpred"],type = "l",ylim=c(0,40),pch=15,col="gold2",xlab="",ylab="Chlorophyll A (mg/L)",main="",lwd=3,lty=2)
#points(envdata2[envdata2$day==126,"dist"],envdata2[envdata2$day==126,"chlpred"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(envdata2[envdata2$day==260,"dist"],envdata2[envdata2$day==260,"chlpred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(envdata2[envdata2$day==336,"dist"],envdata2[envdata2$day==336,"chlpred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(envdata2[envdata2$day==518,"dist"],envdata2[envdata2$day==518,"chlpred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
abline(v=4,lty=2)

### dchla, dtp, dno3, dalk divided by 100
alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata8[complete.cases(alldata[,c("tp")]),]) 
newdata<-alldata8
newdata$dalk<-alldata8$dalk*0.01
alkpred<-exp(predict(tempo,newdata=newdata))
envdata<-(cbind(alkpred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(no3)~temp+rain2+rain3+eff:I(log(dno3)/dayflow):time+eff:I(log(dno3)/dayflow):I(dist==4.08):time,data=alldata)##RH5
newdata<-alldata8
newdata$dno3<-alldata8$dno3*0.01
no3pred<-exp(predict(tempo,newdata=newdata))
envdata<-(cbind(no3pred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(tp~alk+ph2+rain2+rain3+eff:I(dtp/dayflow):time+eff:I(dtp/dayflow):I(dist==4.08):time,data=alldata8[complete.cases(alldata[,c("tp")]),])
newdata<-alldata8
newdata$dtp<-alldata8$dtp*0.01
newdata$alk<-alkpred
tppred<-(predict(tempo,newdata=newdata))
envdata1<-(cbind(tppred,envdata))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(chla)~solar+temp+tp+no3+turb+solar:temp+eff:log(I(dchla/dayflow)),data=alldata8)#
newdata<-alldata8
newdata$dchla<-alldata8$dchla*0.01
newdata$no3<-no3pred
newdata$tp<-tppred
chlpred<-exp(predict(tempo,newdata=newdata))
envdata2<-(cbind(chlpred,envdata1))
##view(envdata2)
plot(envdata2[envdata2$day==126,"dist"],envdata2[envdata2$day==126,"chlpred"],type = "l",ylim=c(0,40),pch=15,col="gold2",xlab="Distance(km)",ylab="Chlorophyll A (mg/L)",main="",lwd=3,lty=2)
#points(envdata2[envdata2$day==126,"dist"],envdata2[envdata2$day==126,"chlpred"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(envdata2[envdata2$day==260,"dist"],envdata2[envdata2$day==260,"chlpred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(envdata2[envdata2$day==336,"dist"],envdata2[envdata2$day==336,"chlpred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(envdata2[envdata2$day==518,"dist"],envdata2[envdata2$day==518,"chlpred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
abline(v=4,lty=2)
dev.off()


##***** Figure 6 *****###””###counterfactual analysis for TP discharge load 

tiff(file = fs::path(script_figure_dir, "20_figure_6_total_phosphorus.jpg"), width = 8, height = 10, units = "in", pointsize = 12, bg = "transparent", res = 800, compression = "lzw")
par(mfrow=c(4,1), mar=c(5,4.5,2.5,3),mai = c(0.6, 1.2, 0.3, 0.5),cex=1.5,cex.axis=0.9,las=1,cex.main=1,cex.lab=0.8)

##total phosphorus plotted against spatial position
plot(bugenv[bugenv$day==1,"dist"],bugenv[bugenv$day==1,"tp"],type = "b",ylim=c(0,0.8),pch=19,xlab="",cex=0.9,ylab="TP (mg/L)",main="",lwd=3,lty=1)
points(bugenv[bugenv$day==126,"dist"],bugenv[bugenv$day==126,"tp"],type="b",col="gold2",pch=15,lwd=3,cex=1,lty=2)
points(bugenv[bugenv$day==260,"dist"],bugenv[bugenv$day==260,"tp"],type="b",col="blue",pch=8,cex=0.8,lwd=3,lty=3)
points(bugenv[bugenv$day==336,"dist"],bugenv[bugenv$day==336,"tp"],type="b",col="green3",pch=17,cex=1,lwd=3,lty=4)
points(bugenv[bugenv$day==518,"dist"],bugenv[bugenv$day==518,"tp"],type="b",col="red",pch=18,cex=1.2,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.6,title="months")
abline(v=4,lty=2)

##TP divided by half

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata8[complete.cases(alldata[,c("tp")]),]) 
newdata<-alldata8
newdata$dalk<-alldata8$dalk*0.5
alkpred<-exp(predict(tempo,newdata)+residuals(tempo))
envdata<-(cbind(alkpred,alldata8))

tempo<-lm(log(tp)~ph2+alk+time+eff:I(dtp/dayflow)+eff:(I(dtp/dayflow)):time,data=alldata8)##checked model2+revision
newdata<-alldata8
newdata<-alldata8
newdata$dtp<-alldata8$dtp*0.5
newdata$alk<-alkpred
tpcounter<-exp(predict(tempo,newdata)+residuals(tempo))
envdata<-(cbind(tpcounter,alldata8))
plot(envdata[envdata$day==126,"dist"],envdata[envdata$day==126,"tpcounter"],type = "l",ylim=c(0,0.8),pch=15,col="gold2",xlab="",ylab="TP (mg/L)",main="",lwd=3,lty=2)
#points(envdata[envdata$day==126,"dist"],envdata[envdata$day==126,"tpcounter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(envdata[envdata$day==260,"dist"],envdata[envdata$day==260,"tpcounter"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(envdata[envdata$day==336,"dist"],envdata[envdata$day==336,"tpcounter"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(envdata[envdata$day==518,"dist"],envdata[envdata$day==518,"tpcounter"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.6,title="months")
abline(v=4,lty=2)

##TP divided by 10
alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dalk<-alldata8$dalk*0.1
alkpred<-exp(predict(tempo,newdata)+residuals(tempo))
envdata<-(cbind(alkpred,alldata8))

tempo<-lm(log(tp)~ph2+alk+time+eff:I(dtp/dayflow)+eff:(I(dtp/dayflow)):time,data=alldata8)
newdata<-alldata8
newdata$dtp<-alldata8$dtp*0.1
newdata$alk<-alkpred
tpcounter<-exp(predict(tempo,newdata)+residuals(tempo))
envdata<-(cbind(tpcounter,alldata8))
plot(envdata[envdata$day==126,"dist"],envdata[envdata$day==126,"tpcounter"],type = "l",ylim=c(0,0.8),pch=15,col="gold2",xlab="",ylab="TP (mg/L)",main="",lwd=3,lty=2)
#points(envdata[envdata$day==126,"dist"],envdata[envdata$day==126,"tpcounter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(envdata[envdata$day==260,"dist"],envdata[envdata$day==260,"tpcounter"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(envdata[envdata$day==336,"dist"],envdata[envdata$day==336,"tpcounter"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(envdata[envdata$day==518,"dist"],envdata[envdata$day==518,"tpcounter"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.6,title="months")
abline(v=4,lty=2)

##TP divided by 100
alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dalk<-alldata8$dalk*0.01
alkpred<-exp(predict(tempo,newdata)+residuals(tempo))
envdata<-(cbind(alkpred,alldata8))

tempo<-lm(log(tp)~ph2+alk+time+eff:I(dtp/dayflow)+eff:(I(dtp/dayflow)):time,data=alldata8)
newdata$dtp<-alldata8$dtp*0.01
newdata$alk<-alkpred
tpcounter<-exp(predict(tempo,newdata)+residuals(tempo))
envdata<-(cbind(tpcounter,alldata8))
plot(envdata[envdata$day==126,"dist"],envdata[envdata$day==126,"tpcounter"],type = "l",ylim=c(0,0.8),pch=15,col="gold2",xlab="Distance (km)",ylab="TP (mg/L)",main="",lwd=3,lty=2)
#points(envdata[envdata$day==126,"dist"],envdata[envdata$day==126,"tpcounter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(envdata[envdata$day==260,"dist"],envdata[envdata$day==260,"tpcounter"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(envdata[envdata$day==336,"dist"],envdata[envdata$day==336,"tpcounter"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(envdata[envdata$day==518,"dist"],envdata[envdata$day==518,"tpcounter"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.6,title="months")
abline(v=4,lty=2)
dev.off()

##***** Figure 8 *****###””#counterfactual analysis for nitrate discharge load  #””# 
tiff(file = fs::path(script_figure_dir, "20_figure_8_nitrate.jpg"), width = 8, height = 10, units = "in", pointsize = 12, bg = "transparent", res = 800, compression = "lzw")
par(mfrow=c(4,1), mar=c(5,4.5,2.5,3),mai = c(0.6, 1.2, 0.3, 0.5),cex=1.5,cex.axis=0.9,las=1,cex.main=1,cex.lab=0.8)

plot(bugenv[bugenv$day==1,"dist"],bugenv[bugenv$day==1,"no3"],type = "b",ylim=c(0,0.8),pch=19,cex=0.9,xlab="",ylab="Nitrate (mg/L)",main="",lwd=3,lty=1)
points(bugenv[bugenv$day==126,"dist"],bugenv[bugenv$day==126,"no3"],type="b",col="gold2",pch=15,lwd=3,cex=1,lty=2)
points(bugenv[bugenv$day==260,"dist"],bugenv[bugenv$day==260,"no3"],type="b",col="blue",pch=8,cex=0.8,lwd=3,lty=3)
points(bugenv[bugenv$day==336,"dist"],bugenv[bugenv$day==336,"no3"],type="b",col="green3",pch=17,cex=1,lwd=3,lty=4)
points(bugenv[bugenv$day==518,"dist"],bugenv[bugenv$day==518,"no3"],type="b",col="red",pch=18,cex=1.2,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.5,title="months")
abline(v=4,lty=2)

## if no3 had been half
tempo<-lm(log(no3)~temp+rain2+rain3+eff:I(log(dno3)/dayflow):time+eff:I(log(dno3)/dayflow):I(dist==4.08):time,data=alldata)##RH5
newdata<-alldata
newdata$dno3<-alldata$dno3*0.5
no3counter<-exp(predict(tempo,newdata)+(residuals(tempo)))
wangenv1<-(cbind(no3counter,alldata))
plot(wangenv1[wangenv1$day==1,"dist"],wangenv1[wangenv1$day==1,"no3counter"],type = "l",ylim=c(0,0.8),pch=19,col="black",xlab="",ylab="Nitrate (mg/L)",main="",lwd=3,lty=1)
points(wangenv1[wangenv1$day==126,"dist"],wangenv1[wangenv1$day==126,"no3counter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(wangenv1[wangenv1$day==260,"dist"],wangenv1[wangenv1$day==260,"no3counter"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(wangenv1[wangenv1$day==336,"dist"],wangenv1[wangenv1$day==336,"no3counter"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(wangenv1[wangenv1$day==518,"dist"],wangenv1[wangenv1$day==518,"no3counter"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.5,title="months")
abline(v=4,lty=2)

## if no3 had been 1/10
tempo<-lm(log(no3)~temp+rain2+rain3+eff:I(log(dno3)/dayflow):time+eff:I(log(dno3)/dayflow):I(dist==4.08):time,data=alldata)##RH5
newdata<-alldata
newdata$dno3<-alldata$dno3*0.1
no3counter<-exp(predict(tempo,newdata)+residuals(tempo))
wangenv1<-(cbind(no3counter,alldata))
plot(wangenv1[wangenv1$day==1,"dist"],wangenv1[wangenv1$day==1,"no3counter"],type = "l",ylim=c(0,0.8),pch=19,col="black",xlab="",ylab="Nitrate (mg/L)",main="",lwd=3,lty=1)
points(wangenv1[wangenv1$day==126,"dist"],wangenv1[wangenv1$day==126,"no3counter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(wangenv1[wangenv1$day==260,"dist"],wangenv1[wangenv1$day==260,"no3counter"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(wangenv1[wangenv1$day==336,"dist"],wangenv1[wangenv1$day==336,"no3counter"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(wangenv1[wangenv1$day==518,"dist"],wangenv1[wangenv1$day==518,"no3counter"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.5,title="months")
abline(v=4,lty=2)

## if no3 had been 1/100
tempo<-lm(log(no3)~temp+rain2+rain3+eff:I(log(dno3)/dayflow):time+eff:I(log(dno3)/dayflow):I(dist==4.08):time,data=alldata)##RH5
newdata<-alldata
newdata$dno3<-alldata$dno3*0.01
no3counter<-exp(predict(tempo,newdata)+residuals(tempo))
wangenv1<-(cbind(no3counter,alldata))
plot(wangenv1[wangenv1$day==1,"dist"],wangenv1[wangenv1$day==1,"no3counter"],type = "l",ylim=c(0,0.8),pch=19,col="black",xlab="",ylab="Nitrate (mg/L)",main="",lwd=3,lty=1)
points(wangenv1[wangenv1$day==126,"dist"],wangenv1[wangenv1$day==126,"no3counter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(wangenv1[wangenv1$day==260,"dist"],wangenv1[wangenv1$day==260,"no3counter"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(wangenv1[wangenv1$day==336,"dist"],wangenv1[wangenv1$day==336,"no3counter"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(wangenv1[wangenv1$day==518,"dist"],wangenv1[wangenv1$day==518,"no3counter"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.5,title="months")
abline(v=4,lty=2)

dev.off()


##***** Figure 9 *****###””# prediction the effect of an intervention (dalk)
tiff(file = fs::path(script_figure_dir, "20_figure_9_alkalinity.jpg"), width = 8, height = 10, units = "in", pointsize = 12, bg = "transparent", res = 800, compression = "lzw")
par(mfrow=c(4,1), mar=c(5,4.5,2.5,3),mai = c(0.6, 1.2, 0.3, 0.5),cex=1.5,cex.axis=0.9,las=1,cex.main=1,cex.lab=0.8)

#alkalinity plotted against spatial position
plot(alldata[alldata$day==1,"dist"],alldata[alldata$day==1,"alk"],type = "b",ylim=c(0,160),pch=19,cex=0.9,xlab="",ylab="Alkalinity (mg/L)",main="",lwd=3,lty=1)
points(alldata[alldata$day==126,"dist"],alldata[alldata$day==126,"alk"],type="b",col="gold2",pch=15,lwd=3,cex=1,lty=2)
points(alldata[alldata$day==260,"dist"],alldata[alldata$day==260,"alk"],type="b",col="blue",pch=8,cex=0.8,lwd=3,lty=3)
points(alldata[alldata$day==336,"dist"],alldata[alldata$day==336,"alk"],type="b",col="green3",pch=17,cex=1,lwd=3,lty=4)
points(alldata[alldata$day==518,"dist"],alldata[alldata$day==518,"alk"],type="b",col="red",pch=18,cex=1.2,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.5,title="months")
abline(v=4,lty=2)

## if alkalinity had been half
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata) # Checked RH5
newdata<-alldata
newdata$dalk<-alldata$dalk*0.5
alkcounter<-exp(predict(tempo,newdata)+residuals(tempo))
wangenv1<-(cbind(alkcounter,newdata))
plot(wangenv1[wangenv1$day==1,"dist"],wangenv1[wangenv1$day==1,"alkcounter"],type = "l",ylim=c(0,160),pch=19,col="black",xlab="",ylab="Alkalinity(mg/L)",main="",lwd=3,lty=1)
points(wangenv1[wangenv1$day==126,"dist"],wangenv1[wangenv1$day==126,"alkcounter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(wangenv1[wangenv1$day==260,"dist"],wangenv1[wangenv1$day==260,"alkcounter"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(wangenv1[wangenv1$day==336,"dist"],wangenv1[wangenv1$day==336,"alkcounter"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(wangenv1[wangenv1$day==518,"dist"],wangenv1[wangenv1$day==518,"alkcounter"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
abline(v=4,lty=2)

## if alkalinity had been 1/10
# first, Toc is needed for alkalinity plot
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata) # Checked RH5
newdata<-alldata
newdata$dalk<-alldata$dalk*0.1
##view(newdata)
alkcounter<-exp(predict(tempo,newdata)+residuals(tempo))
wangenv1<-(cbind(alkcounter,alldata))
plot(wangenv1[wangenv1$day==1,"dist"],wangenv1[wangenv1$day==1,"alkcounter"],type = "l",ylim=c(0,160),pch=19,col="black",xlab="",ylab="Alkalinity(mg/L)",main="",lwd=3,lty=1)
points(wangenv1[wangenv1$day==126,"dist"],wangenv1[wangenv1$day==126,"alkcounter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(wangenv1[wangenv1$day==260,"dist"],wangenv1[wangenv1$day==260,"alkcounter"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(wangenv1[wangenv1$day==336,"dist"],wangenv1[wangenv1$day==336,"alkcounter"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(wangenv1[wangenv1$day==518,"dist"],wangenv1[wangenv1$day==518,"alkcounter"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
abline(v=4,lty=2)

## if alkalinity had been 1/100
# first, Toc is needed for alkalinity plot
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata) # Checked RH5
newdata<-alldata
newdata$dalk<-alldata$dalk*0.01
##view(newdata)
alkcounter<-exp(predict(tempo,newdata)+residuals(tempo))
wangenv1<-(cbind(alkcounter,alldata))
plot(wangenv1[wangenv1$day==1,"dist"],wangenv1[wangenv1$day==1,"alkcounter"],type = "l",ylim=c(0,160),pch=19,col="black",xlab="",ylab="Alkalinity(mg/L)",main="",lwd=3,lty=1)
points(wangenv1[wangenv1$day==126,"dist"],wangenv1[wangenv1$day==126,"alkcounter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(wangenv1[wangenv1$day==260,"dist"],wangenv1[wangenv1$day==260,"alkcounter"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
points(wangenv1[wangenv1$day==336,"dist"],wangenv1[wangenv1$day==336,"alkcounter"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
points(wangenv1[wangenv1$day==518,"dist"],wangenv1[wangenv1$day==518,"alkcounter"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.6,title="months")
abline(v=4,lty=2)
dev.off()


#””#””#””#””#*****Figure 10 *****###””#””#””#””#””############Counterfactual prediction for bugs ##

tiff(file = fs::path(script_figure_dir, "20_figure_10_macroinvertebrate_pco1.jpg"), width = 8, height = 14, units = "in", pointsize = 12, bg = "transparent", res = 800, compression = "lzw")
par(mfrow=c(4,1), mar=c(4,4.5,2.5,3),cex=1.2,cex.axis=0.9,las=1,cex.main=1,cex.lab=0.8)

## pco plot against distance
plot(bugenv[bugenv$day==1,"dist"],bugenv[bugenv$day==1,"pco1"],type = "b",ylim=c(-.6,1),pch=19,col="black",cex=0.9,xlab="",ylab="PCO1 observed values",main="",lwd=3,lty=2)
points(bugenv[bugenv$day==126,"dist"],bugenv[bugenv$day==126,"pco1"],type="b",col="gold2",pch=15,lwd=3,cex=1,lty=1)
points(bugenv[bugenv$day==260,"dist"],bugenv[bugenv$day==260,"pco1"],type="b",col="blue",pch=8,cex=1,lwd=3,lty=3)
points(bugenv[bugenv$day==336,"dist"],bugenv[bugenv$day==336,"pco1"],type="b",col="green3",pch=17,cex=1,lwd=3,lty=4)
points(bugenv[bugenv$day==518,"dist"],bugenv[bugenv$day==518,"pco1"],type="b",col="red",pch=18,cex=1,lwd=3,lty=5)
#legend("topright",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(2,1,3,4,5), pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.6,title="months")
abline(v=4,lty=2)


##””””## all the discharge variables divided by half
##””””## create new data by dividing discharge variables by half


alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(toc~dayflow+rain2+rain3+eff:I((dtoc/dayflow)):time+eff:I((dtoc/dayflow)):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dtoc<-alldata8$dtoc*0.5
tocpred<-predict(tempo,newdata)

tocpred<-predict(tempo,newdata)
wangintervention<-(cbind(tocpred,alldata8))
#plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"tocpred"],type = "l",ylim=c(0,10),pch=15,col="gold2",xlab="",ylab="TOC (mg/L)",main="If Toc in discharge had been devided by half",lwd=3,lty=2)
#points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"tocpred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
#points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"tocpred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
#points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"tocpred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
#abline(v=4,lty=2)

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(no3)~temp+rain2+rain3+eff:I(log(dno3)/dayflow):time+eff:I(log(dno3)/dayflow):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dno3<-alldata8$dno3*0.5
newdata$no3<-no3pred
no3pred<-exp(predict(tempo,newdata))
wangintervention<-(cbind(no3pred,alldata8))
#plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"no3pred"],type = "l",ylim=c(0,0.7),pch=15,col="gold2",xlab="",ylab="NO3 (mg/L)",main="If NO3 in discharge had been devided by half",lwd=3,lty=2)
#points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"no3pred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
#points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"no3pred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
#points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"no3pred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
#abline(v=4,lty=2)

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dalk<-alldata8$dalk*0.5
newdata$alk<-alkpred
alkpred<-exp(predict(tempo,newdata))
wangintervention<-(cbind(alkpred,alldata8))
#plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"alkpred"],type = "l",ylim=c(0,150),pch=15,col="gold2",xlab="",ylab="Alkalinity (mg/L)",main="If alkalinity in discharge had been devided by half",lwd=3,lty=2)
#points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"alkpred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
#points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"alkpred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
#points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"alkpred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
#abline(v=4,lty=2)

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(tp~alk+ph2+rain2+rain3+eff:I(dtp/dayflow):time+eff:I(dtp/dayflow):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dtp<-alldata8$dtp*0.5
newdata$alk<-alkpred
tppred<-predict(tempo,newdata)
wangintervention<-(cbind(tppred,alldata8))
#plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"tppred"],type = "l",ylim=c(0,0.8),pch=15,col="gold2",xlab="TP",ylab="TP (mg/L)",main="If alkalinity in discharge had been devided by half",lwd=3,lty=2)
#points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"tppred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
#points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"tppred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
#points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"tppred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
#abline(v=4,lty=2)

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(chla)~solar+temp+tp+no3+turb+solar:temp+eff:log(I(dchla/dayflow)),data=alldata8)# 
newdata<-alldata8
newdata$dchla<-alldata8$dchla*0.5
newdata$no3<-no3pred
newdata$tp<-tppred
chlpred<-exp(predict(tempo,newdata))
wangintervention<-(cbind(chlpred,envdata1))
#plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"chlpred"],type = "l",ylim=c(0,40),pch=15,col="gold2",xlab="",ylab="Chlorophyll A (mg/L)",main="Prediction the effect of an intervention \n (when chlorophyll a, TP, and nitrate in discharge were divided by half)",lwd=3,lty=2)
#points(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"chlpred"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
#points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"chlpred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
#points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"chlpred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
#points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"chlpred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
#abline(v=4,lty=2)

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dalk<-alldata8$dalk*0.5
alkpred<-exp(predict(tempo,newdata))
wangintervention<-(cbind(alkpred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(tp~alk+ph2+rain2+rain3+eff:I(dtp/dayflow):time+eff:I(dtp/dayflow):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dtp<-alldata8$dtp*0.5
newdata$alk<-alkpred
tppred<-predict(tempo,newdata)
wangintervention<-(cbind(tppred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(cond)~tp+alk+time+dist:time:eff,data=alldata[complete.cases(alldata[,c("cond","tp")]),])##
newdata<-alldata8
newdata$tp<-tppred
newdata$alk<-alkpred
condpred<-exp(predict(tempo,newdata))
wangintervention<-(cbind(condpred,alldata8))
#plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"condpred"],type = "l",ylim=c(0,700),pch=15,col="gold2",xlab="",ylab="Conductivity (mg/L)",main=" Prediction the effect of an intervention \n (when TP and alkalinity were divided by half)",lwd=3,lty=2)
#points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"condpred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
#points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"condpred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
#points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"condpred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
#abline(v=4,lty=2)


### PCO1 prediction plot with newdata
alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(pco1~temp+cond+temp:cond+log(zn)+toc+temp:toc+log(chla),data=alldata8)##final model
newdata<-alldata8
newdata$toc<-tocpred
newdata$tp<-tppred
newdata$chla<-chlpred
newdata$cond<-condpred
newdata$alk<-alkpred
pco1counter<-predict(tempo,newdata)+residuals(tempo)
wangintervention<-(cbind(pco1counter,alldata8))
plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"pco1counter"],type = "l",ylim=c(-.6,1),pch=15,col="gold2",cex=0.9,xlab="",ylab="PCO1 predicted values",main="",lwd=3,lty=1)
#points(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"pco1counter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"pco1counter"],type="l",col="blue",pch=8,cex=1,lwd=3,lty=3)
points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"pco1counter"],type="l",col="green3",pch=17,cex=1,lwd=3,lty=4)
points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"pco1counter"],type="l",col="red",pch=18,cex=1,lwd=3,lty=5)
#legend("topright",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.6,title="months")
abline(v=4,lty=2)

##””””## all the discharge variables divided by 0.1
##””””## create new data by dividing discharge variables by 0.1

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(toc~dayflow+rain2+rain3+eff:I((dtoc/dayflow)):time+eff:I((dtoc/dayflow)):I(dist==4.08):time,data=alldata8)# Checked RH5
newdata<-alldata8
newdata$dtoc<-alldata8$dtoc*0.1
newdata$toc<-tocpred
tocpred<-predict(tempo,newdata)
wangintervention<-(cbind(tocpred,alldata8))
#plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"tocpred"],type = "l",ylim=c(0,10),pch=15,col="gold2",xlab="",ylab="TOC (mg/L)",main="If Toc in discharge had been devided by half",lwd=3,lty=2)
#points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"tocpred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
#points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"tocpred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
#points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"tocpred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
#abline(v=4,lty=2)

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(no3)~temp+rain2+rain3+eff:I(log(dno3)/dayflow):time+eff:I(log(dno3)/dayflow):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dno3<-alldata8$dno3*0.1
newdata$no3<-no3pred
no3pred<-exp(predict(tempo,newdata))
wangintervention<-(cbind(no3pred,alldata8))
#plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"no3pred"],type = "l",ylim=c(0,0.7),pch=15,col="gold2",xlab="",ylab="NO3 (mg/L)",main="If NO3 in discharge had been devided by half",lwd=3,lty=2)
#points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"no3pred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
#points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"no3pred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
#points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"no3pred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
#abline(v=4,lty=2)

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dalk<-alldata8$dalk*0.1
newdata$alk<-alkpred
alkpred<-exp(predict(tempo,newdata))
wangintervention<-(cbind(alkpred,alldata8))
#plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"alkpred"],type = "l",ylim=c(0,150),pch=15,col="gold2",xlab="",ylab="Alkalinity (mg/L)",main="If alkalinity in discharge had been devided by half",lwd=3,lty=2)
#points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"alkpred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
#points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"alkpred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
#points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"alkpred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
#abline(v=4,lty=2)

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(tp~alk+ph2+rain2+rain3+eff:I(dtp/dayflow):time+eff:I(dtp/dayflow):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dtp<-alldata8$dtp*0.1
newdata$alk<-alkpred
tppred<-predict(tempo,newdata)
wangintervention<-(cbind(tppred,alldata8))
#plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"tppred"],type = "l",ylim=c(0,0.8),pch=15,col="gold2",xlab="TP",ylab="TP (mg/L)",main="If alkalinity in discharge had been devided by half",lwd=3,lty=2)
#points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"tppred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
#points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"tppred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
#points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"tppred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
#abline(v=4,lty=2)

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(chla)~solar+temp+tp+no3+turb+solar:temp+eff:log(I(dchla/dayflow)),data=alldata8)# Checked RH
newdata<-alldata8
newdata$dchla<-alldata8$dchla*0.1
newdata$no3<-no3pred
newdata$tp<-tppred
chlpred<-exp(predict(tempo,newdata))
wangintervention<-(cbind(chlpred,envdata1))
#plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"chlpred"],type = "l",ylim=c(0,40),pch=15,col="gold2",xlab="",ylab="Chlorophyll A (mg/L)",main="Prediction the effect of an intervention \n (when chlorophyll a, TP, and nitrate in discharge were divided by half)",lwd=3,lty=2)
#points(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"chlpred"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
#points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"chlpred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
#points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"chlpred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
#points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"chlpred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
#abline(v=4,lty=2)

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dalk<-alldata8$dalk*0.1
alkpred<-exp(predict(tempo,newdata))
wangintervention<-(cbind(alkpred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(tp~alk+ph2+rain2+rain3+eff:I(dtp/dayflow):time+eff:I(dtp/dayflow):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dtp<-alldata8$dtp*0.1
newdata$alk<-alkpred
tppred<-predict(tempo,newdata)
wangintervention<-(cbind(tppred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(cond)~tp+alk+time+dist:time:eff,data=alldata[complete.cases(alldata[,c("cond","tp")]),])##final model Checked RH3
newdata<-alldata8
newdata$tp<-tppred
newdata$alk<-alkpred
condpred<-exp(predict(tempo,newdata))
wangintervention<-(cbind(condpred,alldata8))
#plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"condpred"],type = "l",ylim=c(0,700),pch=15,col="gold2",xlab="",ylab="Conductivity (mg/L)",main=" Prediction the effect of an intervention \n (when TP and alkalinity were divided by half)",lwd=3,lty=2)
#points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"condpred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
#points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"condpred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
#points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"condpred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
#abline(v=4,lty=2)


### PCO1 prediction plot with newdata
alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(pco1~temp+cond+temp:cond+log(zn)+toc+temp:toc+log(chla),data=alldata8)##final model
newdata<-alldata8
newdata$toc<-tocpred
newdata$tp<-tppred
newdata$chla<-chlpred
newdata$cond<-condpred
newdata$alk<-alkpred
pco1counter<-predict(tempo,newdata)+residuals(tempo)
wangintervention<-(cbind(pco1counter,alldata8))
plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"pco1counter"],type = "l",ylim=c(-.6,1),pch=15,col="gold2",cex=0.9,xlab="",ylab="PCO1 predicted values",main="",lwd=3,lty=1)
#points(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"pco1counter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"pco1counter"],type="l",col="blue",pch=8,cex=1,lwd=3,lty=3)
points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"pco1counter"],type="l",col="green3",pch=17,cex=1,lwd=3,lty=4)
points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"pco1counter"],type="l",col="red",pch=18,cex=1,lwd=3,lty=5)
#legend("topright",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.6,title="months")
abline(v=4,lty=2)


##””””## all the discharge variables divided by 0.01
##””””## create new data by dividing discharge variables by 0.01

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(toc~dayflow+rain2+rain3+eff:I((dtoc/dayflow)):time+eff:I((dtoc/dayflow)):I(dist==4.08):time,data=alldata8)# Checked RH5
newdata<-alldata8
newdata$dtoc<-alldata8$dtoc*0.01
newdata$toc<-tocpred
tocpred<-predict(tempo,newdata)
wangintervention<-(cbind(tocpred,alldata8))
#plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"tocpred"],type = "l",ylim=c(0,10),pch=15,col="gold2",xlab="",ylab="TOC (mg/L)",main="If Toc in discharge had been devided by half",lwd=3,lty=2)
#points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"tocpred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
#points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"tocpred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
#points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"tocpred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
#abline(v=4,lty=2)

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(no3)~temp+rain2+rain3+eff:I(log(dno3)/dayflow):time+eff:I(log(dno3)/dayflow):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dno3<-alldata8$dno3*0.01
newdata$no3<-no3pred
no3pred<-exp(predict(tempo,newdata))
wangintervention<-(cbind(no3pred,alldata8))
#plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"no3pred"],type = "l",ylim=c(0,0.7),pch=15,col="gold2",xlab="",ylab="NO3 (mg/L)",main="If NO3 in discharge had been devided by half",lwd=3,lty=2)
#points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"no3pred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
#points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"no3pred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
#points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"no3pred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
#abline(v=4,lty=2)

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dalk<-alldata8$dalk*0.01
newdata$alk<-alkpred
alkpred<-exp(predict(tempo,newdata))
wangintervention<-(cbind(alkpred,alldata8))
#plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"alkpred"],type = "l",ylim=c(0,150),pch=15,col="gold2",xlab="",ylab="Alkalinity (mg/L)",main="If alkalinity in discharge had been devided by half",lwd=3,lty=2)
#points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"alkpred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
#points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"alkpred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
#points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"alkpred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
#abline(v=4,lty=2)

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(tp~alk+ph2+rain2+rain3+eff:I(dtp/dayflow):time+eff:I(dtp/dayflow):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dtp<-alldata8$dtp*0.01
newdata$alk<-alkpred
tppred<-predict(tempo,newdata)
wangintervention<-(cbind(tppred,alldata8))
#plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"tppred"],type = "l",ylim=c(0,0.8),pch=15,col="gold2",xlab="TP",ylab="TP (mg/L)",main="If alkalinity in discharge had been devided by half",lwd=3,lty=2)
#points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"tppred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
#points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"tppred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
#points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"tppred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),pch=c(19,15,8,17,18),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
#abline(v=4,lty=2)

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(chla)~solar+temp+tp+no3+turb+solar:temp+eff:log(I(dchla/dayflow)),data=alldata8)# 
newdata<-alldata8
newdata$dchla<-alldata8$dchla*0.01
newdata$no3<-no3pred
newdata$tp<-tppred
chlpred<-exp(predict(tempo,newdata))
wangintervention<-(cbind(chlpred,envdata1))
#plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"chlpred"],type = "l",ylim=c(0,40),pch=15,col="gold2",xlab="",ylab="Chlorophyll A (mg/L)",main="Prediction the effect of an intervention \n (when chlorophyll a, TP, and nitrate in discharge were divided by 0.01)",lwd=3,lty=2)
#points(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"chlpred"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
#points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"chlpred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
#points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"chlpred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
#points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"chlpred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
#abline(v=4,lty=2)

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(alk)~rain3+dayflow+eff:I(log(dalk/dayflow)):time+eff:I(log(dalk/dayflow)):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dalk<-alldata8$dalk*0.01
alkpred<-exp(predict(tempo,newdata))
wangintervention<-(cbind(alkpred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(tp~alk+ph2+rain2+rain3+eff:I(dtp/dayflow):time+eff:I(dtp/dayflow):I(dist==4.08):time,data=alldata8)
newdata<-alldata8
newdata$dtp<-alldata8$dtp*0.01
newdata$alk<-alkpred
tppred<-predict(tempo,newdata)
wangintervention<-(cbind(tppred,alldata8))

alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(log(cond)~tp+alk+time+dist:time:eff,data=alldata[complete.cases(alldata[,c("cond","tp")]),])##
newdata<-alldata8
newdata$tp<-tppred
newdata$alk<-alkpred
condpred<-exp(predict(tempo,newdata))
wangintervention<-(cbind(condpred,alldata8))
#plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"condpred"],type = "l",ylim=c(0,700),pch=15,col="gold2",xlab="",ylab="Conductivity (mg/L)",main=" Prediction the effect of an intervention \n (when TP and alkalinity were divided by half)",lwd=3,lty=2)
#points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"condpred"],type="l",col="blue",pch=8,cex=1.2,lwd=3,lty=3)
#points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"condpred"],type="l",col="green3",pch=17,cex=1.2,lwd=3,lty=4)
#points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"condpred"],type="l",col="red",pch=18,cex=1.4,lwd=3,lty=5)
#legend("topleft",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=2,horiz=FALSE,cex=0.6,title="months")
#abline(v=4,lty=2)


### PCO1 prediction plot with newdata
alldata8<-alldata[complete.cases(alldata[,c("cond","tp")]),]
tempo<-lm(pco1~temp+cond+temp:cond+log(zn)+toc+temp:toc+log(chla),data=alldata8)##
newdata<-alldata8
newdata$toc<-tocpred
newdata$tp<-tppred
newdata$chla<-chlpred
newdata$cond<-condpred
newdata$alk<-alkpred
pco1counter<-predict(tempo,newdata)+residuals(tempo)
wangintervention<-(cbind(pco1counter,alldata8))
plot(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"pco1counter"],type = "l",ylim=c(-.6,1),pch=19,col="gold2",cex=0.9,xlab="",ylab="PCO1 predicted values",lwd=3,lty=1)
#points(wangintervention[wangintervention$day==126,"dist"],wangintervention[wangintervention$day==126,"pco1counter"],type="l",col="gold2",pch=15,lwd=3,cex=1.5,lty=2)
points(wangintervention[wangintervention$day==260,"dist"],wangintervention[wangintervention$day==260,"pco1counter"],type="l",col="blue",pch=8,cex=1,lwd=3,lty=3)
points(wangintervention[wangintervention$day==336,"dist"],wangintervention[wangintervention$day==336,"pco1counter"],type="l",col="green3",pch=17,cex=1,lwd=3,lty=4)
points(wangintervention[wangintervention$day==518,"dist"],wangintervention[wangintervention$day==518,"pco1counter"],type="l",col="red",pch=18,cex=1,lwd=3,lty=5)
#legend("topright",inset=c(0,0),legend=c("Dec 13", "April 14", "Aug 2014", "Nov 2014","May 2015"),lty=c(1,5),lwd=2,col=c("black","gold2","blue","green3","red"),ncol=3,horiz=FALSE,cex=0.6,title="months")
abline(v=4,lty=2)
dev.off()



figure_dir <- script_figure_dir

set.seed(123)

#------------------------------------------------------------
# Figure 2. Counterfactual predictions of conductivity
#------------------------------------------------------------

alldata8 <- alldata[complete.cases(alldata[, c("cond", "tp")]), ]
alldata8 <- alldata8[order(alldata8$day, alldata8$dist), ]
alldata8$time <- factor(alldata8$time, levels = levels(alldata$time))

m_alk0 <- lm(
  log(alk) ~ rain3 + dayflow +
    eff:I(log(dalk / dayflow)):time +
    eff:I(log(dalk / dayflow)):I(dist == 4.08):time,
  data = alldata8
)

m_tp0 <- lm(
  tp ~ alk + ph2 + rain2 + rain3 +
    eff:I(dtp / dayflow):time +
    eff:I(dtp / dayflow):I(dist == 4.08):time,
  data = alldata8
)

m_cond0 <- lm(
  log(cond) ~ tp + alk + time + dist:time:eff,
  data = alldata8
)

fit_residual_boot_chain <- function(base_data, m_alk0, m_tp0, m_cond0) {
  boot_data <- base_data
  boot_data$time <- factor(boot_data$time, levels = levels(base_data$time))
  
  fitted_logalk <- predict(m_alk0, newdata = boot_data)
  res_alk <- sample(residuals(m_alk0), size = nrow(boot_data), replace = TRUE)
  boot_data$alk <- exp(fitted_logalk + res_alk)
  
  fitted_tp <- predict(m_tp0, newdata = boot_data)
  res_tp <- sample(residuals(m_tp0), size = nrow(boot_data), replace = TRUE)
  boot_data$tp <- fitted_tp + res_tp
  
  fitted_logcond <- predict(m_cond0, newdata = boot_data)
  res_cond <- sample(residuals(m_cond0), size = nrow(boot_data), replace = TRUE)
  boot_data$cond <- exp(fitted_logcond + res_cond)
  
  boot_data$time <- factor(boot_data$time, levels = levels(base_data$time))
  
  m_alk_b <- lm(
    log(alk) ~ rain3 + dayflow +
      eff:I(log(dalk / dayflow)):time +
      eff:I(log(dalk / dayflow)):I(dist == 4.08):time,
    data = boot_data
  )
  
  m_tp_b <- lm(
    tp ~ alk + ph2 + rain2 + rain3 +
      eff:I(dtp / dayflow):time +
      eff:I(dtp / dayflow):I(dist == 4.08):time,
    data = boot_data
  )
  
  m_cond_b <- lm(
    log(cond) ~ tp + alk + time + dist:time:eff,
    data = boot_data
  )
  
  list(m_alk = m_alk_b, m_tp = m_tp_b, m_cond = m_cond_b)
}

predict_cond_chain <- function(reduction_factor, base_data, m_alk_b, m_tp_b, m_cond_b) {
  newdata <- base_data
  newdata$time <- factor(newdata$time, levels = levels(base_data$time))
  newdata$dalk <- base_data$dalk * reduction_factor
  newdata$dtp  <- base_data$dtp  * reduction_factor
  
  alkpred <- exp(predict(m_alk_b, newdata = newdata))
  newdata$alk <- alkpred
  
  tppred <- predict(m_tp_b, newdata = newdata)
  newdata$tp <- tppred
  
  exp(predict(m_cond_b, newdata = newdata))
}

simulate_cond_chain_residual <- function(reduction_factor, base_data, m_alk0, m_tp0, m_cond0, n_boot = 500L) {
  cond_draws <- matrix(NA_real_, nrow = nrow(base_data), ncol = n_boot)
  
  for (b in seq_len(n_boot)) {
    res <- try({
      boot_models <- fit_residual_boot_chain(base_data, m_alk0, m_tp0, m_cond0)
      cond_draws[, b] <- predict_cond_chain(
        reduction_factor = reduction_factor,
        base_data = base_data,
        m_alk_b = boot_models$m_alk,
        m_tp_b = boot_models$m_tp,
        m_cond_b = boot_models$m_cond
      )
    }, silent = TRUE)
    
    if (inherits(res, "try-error")) cond_draws[, b] <- NA_real_
  }
  
  out <- base_data[, c("day", "dist")]
  out$lower  <- apply(cond_draws, 1, quantile, probs = 0.025, na.rm = TRUE)
  out$median <- apply(cond_draws, 1, median,   na.rm = TRUE)
  out$upper  <- apply(cond_draws, 1, quantile, probs = 0.975, na.rm = TRUE)
  out
}

get_band_fill <- function(line_col) {
  switch(
    line_col,
    "black"  = "grey70",
    "gold2"  = "khaki1",
    "blue"   = "lightskyblue1",
    "green3" = "palegreen2",
    "red"    = "mistyrose2",
    line_col
  )
}

get_band_edge <- function(line_col) {
  switch(
    line_col,
    "black"  = "grey35",
    "gold2"  = "goldenrod3",
    "blue"   = "blue",
    "green3" = "green3",
    "red"    = "red3",
    line_col
  )
}

draw_band_lines <- function(df_day, line_col, line_lty, line_pch = NA,
                            line_lwd = 2.8, point_cex = 1.0) {
  df_day <- df_day[order(df_day$dist), ]
  if (nrow(df_day) == 0) return(invisible(NULL))
  
  fill_col <- get_band_fill(line_col)
  edge_col <- get_band_edge(line_col)
  
  polygon(
    x = c(df_day$dist, rev(df_day$dist)),
    y = c(df_day$lower, rev(df_day$upper)),
    col = fill_col,
    border = NA
  )
  
  lines(df_day$dist, df_day$lower, col = edge_col, lwd = 1.0, lty = 1)
  lines(df_day$dist, df_day$upper, col = edge_col, lwd = 1.0, lty = 1)
  lines(df_day$dist, df_day$median, col = line_col, lwd = line_lwd, lty = line_lty)
  
  if (!is.na(line_pch)) {
    points(df_day$dist, df_day$median, col = line_col, pch = line_pch, cex = point_cex)
  }
}

panel_style <- function() {
  par(
    mar = c(3.3, 4.8, 1.5, 0.8),
    mgp = c(2.35, 0.75, 0),
    cex = 1.15,
    cex.axis = 1.0,
    cex.lab = 1.15,
    font.lab = 2,
    las = 1,
    bty = "o",
    xaxs = "r",
    yaxs = "r"
  )
}

draw_panel_tag <- function(tag) {
  usr <- par("usr")
  text(
    x = usr[1] + 0.035 * (usr[2] - usr[1]),
    y = usr[4] - 0.085 * (usr[4] - usr[3]),
    labels = tag,
    font = 2,
    cex = 1.25,
    xpd = NA,
    adj = c(0, 1)
  )
}

n_boot <- 500L

cond_half      <- simulate_cond_chain_residual(0.50, alldata8, m_alk0, m_tp0, m_cond0, n_boot = n_boot)
cond_tenth     <- simulate_cond_chain_residual(0.10, alldata8, m_alk0, m_tp0, m_cond0, n_boot = n_boot)
cond_hundredth <- simulate_cond_chain_residual(0.01, alldata8, m_alk0, m_tp0, m_cond0, n_boot = n_boot)

obs_days  <- c(1, 126, 260, 336, 518)
cols_show <- c("black", "gold2", "blue", "green3", "red")
ltys_show <- c(1, 2, 3, 4, 5)
pchs_show <- c(19, 15, 8, 17, 18)
cexs_show <- c(1.0, 1.0, 1.05, 1.0, 1.0)

outfile_tif <- file.path(figure_dir, "Figure_2_counterfactual_conductivity.tif")

plot_figure_2 <- function(device_file) {
  tiff(
    file = device_file,
    width = 10,
    height = 13.5,
    units = "in",
    pointsize = 12,
    bg = "white",
    res = 800,
    compression = "lzw"
  )
  
  layout(matrix(c(1, 2, 3, 4, 5), ncol = 1), heights = c(1.15, 1.15, 1.15, 1.15, 0.55))
  
  panel_style()
  plot(
    bugenv[bugenv$day == 1, "dist"],
    bugenv[bugenv$day == 1, "cond"],
    type = "b",
    ylim = c(0, 600),
    xlim = range(alldata8$dist),
    pch = 19,
    cex = 1.0,
    xlab = "",
    ylab = "Conductivity (mg/L)",
    lwd = 2.8,
    lty = 1
  )
  points(bugenv[bugenv$day == 126, "dist"], bugenv[bugenv$day == 126, "cond"], type = "b", col = "gold2",  pch = 15, lwd = 2.8, cex = 1.0,  lty = 2)
  points(bugenv[bugenv$day == 260, "dist"], bugenv[bugenv$day == 260, "cond"], type = "b", col = "blue",   pch = 8,  lwd = 2.8, cex = 1.05, lty = 3)
  points(bugenv[bugenv$day == 336, "dist"], bugenv[bugenv$day == 336, "cond"], type = "b", col = "green3", pch = 17, lwd = 2.8, cex = 1.0,  lty = 4)
  points(bugenv[bugenv$day == 518, "dist"], bugenv[bugenv$day == 518, "cond"], type = "b", col = "red",    pch = 18, lwd = 2.8, cex = 1.0,  lty = 5)
  abline(v = 4, lty = 2, lwd = 1.2)
  box(lwd = 1.2)
  draw_panel_tag("(a)")
  
  panel_style()
  plot(
    cond_half[cond_half$day == 126, "dist"],
    cond_half[cond_half$day == 126, "median"],
    type = "n",
    xlim = range(alldata8$dist),
    ylim = c(0, 600),
    xlab = "",
    ylab = "Conductivity (mg/L)"
  )
  for (i in 2:5) {
    draw_band_lines(cond_half[cond_half$day == obs_days[i], ], cols_show[i], ltys_show[i], pchs_show[i], point_cex = cexs_show[i])
  }
  abline(v = 4, lty = 2, lwd = 1.2)
  box(lwd = 1.2)
  draw_panel_tag("(b)")
  
  panel_style()
  plot(
    cond_tenth[cond_tenth$day == 126, "dist"],
    cond_tenth[cond_tenth$day == 126, "median"],
    type = "n",
    xlim = range(alldata8$dist),
    ylim = c(0, 600),
    xlab = "",
    ylab = "Conductivity (mg/L)"
  )
  for (i in 2:5) {
    draw_band_lines(cond_tenth[cond_tenth$day == obs_days[i], ], cols_show[i], ltys_show[i], pchs_show[i], point_cex = cexs_show[i])
  }
  abline(v = 4, lty = 2, lwd = 1.2)
  box(lwd = 1.2)
  draw_panel_tag("(c)")
  
  panel_style()
  plot(
    cond_hundredth[cond_hundredth$day == 126, "dist"],
    cond_hundredth[cond_hundredth$day == 126, "median"],
    type = "n",
    xlim = range(alldata8$dist),
    ylim = c(0, 600),
    xlab = "Distance(km)",
    ylab = "Conductivity (mg/L)"
  )
  for (i in 2:5) {
    draw_band_lines(cond_hundredth[cond_hundredth$day == obs_days[i], ], cols_show[i], ltys_show[i], pchs_show[i], point_cex = cexs_show[i])
  }
  abline(v = 4, lty = 2, lwd = 1.2)
  box(lwd = 1.2)
  draw_panel_tag("(d)")
  
  par(mar = c(0.4, 1.5, 0.2, 1.5))
  plot.new()
  legend(
    "center",
    legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
    lty = c(1, 2, 3, 4, 5),
    pch = c(19, 15, 8, 17, 18),
    lwd = 2.5,
    pt.cex = c(1.0, 1.0, 1.05, 1.0, 1.0),
    col = c("black", "gold2", "blue", "green3", "red"),
    ncol = 5,
    cex = 1.0,
    bty = "o",
    box.lwd = 1.1,
    x.intersp = 0.8,
    y.intersp = 1.0
  )
  
  dev.off()
}

plot_figure_2(outfile_tif)

#------------------------------------------------------------
# Figure 3. Counterfactual predictions of TOC
#------------------------------------------------------------

toc_vars <- c("toc", "dayflow", "rain2", "rain3", "dtoc", "time", "dist", "eff")
toc_base <- alldata[complete.cases(alldata[, toc_vars]), ]
toc_base <- toc_base[order(toc_base$day, toc_base$dist), ]
toc_base$time <- factor(toc_base$time, levels = levels(alldata$time))

m_toc0 <- lm(
  toc ~ dayflow + rain2 + rain3 +
    eff:I((dtoc / dayflow)):time +
    eff:I((dtoc / dayflow)):I(dist == 4.08):time,
  data = toc_base
)

fit_residual_boot_toc <- function(base_data, m_toc0) {
  boot_data <- base_data
  boot_data$time <- factor(boot_data$time, levels = levels(base_data$time))
  
  fitted_toc <- predict(m_toc0, newdata = boot_data)
  res_toc <- sample(residuals(m_toc0), size = nrow(boot_data), replace = TRUE)
  boot_data$toc <- fitted_toc + res_toc
  boot_data$time <- factor(boot_data$time, levels = levels(base_data$time))
  
  lm(
    toc ~ dayflow + rain2 + rain3 +
      eff:I((dtoc / dayflow)):time +
      eff:I((dtoc / dayflow)):I(dist == 4.08):time,
    data = boot_data
  )
}

predict_toc_residual <- function(reduction_factor, base_data, m_toc_b) {
  newdata <- base_data
  newdata$time <- factor(newdata$time, levels = levels(base_data$time))
  newdata$dtoc <- base_data$dtoc * reduction_factor
  predict(m_toc_b, newdata = newdata)
}

simulate_toc_residual <- function(reduction_factor, base_data, m_toc0, n_boot = 500L) {
  toc_draws <- matrix(NA_real_, nrow = nrow(base_data), ncol = n_boot)
  
  for (b in seq_len(n_boot)) {
    res <- try({
      m_toc_b <- fit_residual_boot_toc(base_data, m_toc0)
      toc_draws[, b] <- predict_toc_residual(
        reduction_factor = reduction_factor,
        base_data = base_data,
        m_toc_b = m_toc_b
      )
    }, silent = TRUE)
    
    if (inherits(res, "try-error")) toc_draws[, b] <- NA_real_
  }
  
  out <- base_data[, c("day", "dist")]
  out$lower  <- apply(toc_draws, 1, quantile, probs = 0.025, na.rm = TRUE)
  out$median <- apply(toc_draws, 1, median,   na.rm = TRUE)
  out$upper  <- apply(toc_draws, 1, quantile, probs = 0.975, na.rm = TRUE)
  out
}

toc_half      <- simulate_toc_residual(0.50, toc_base, m_toc0, n_boot = n_boot)
toc_tenth     <- simulate_toc_residual(0.10, toc_base, m_toc0, n_boot = n_boot)
toc_hundredth <- simulate_toc_residual(0.01, toc_base, m_toc0, n_boot = n_boot)

outfile_tif <- file.path(figure_dir, "Figure_3_counterfactual_TOC.tif")

plot_figure_3 <- function(device_file) {
  tiff(
    file = device_file,
    width = 10,
    height = 13.5,
    units = "in",
    pointsize = 12,
    bg = "white",
    res = 800,
    compression = "lzw"
  )
  
  layout(matrix(c(1, 2, 3, 4, 5), ncol = 1), heights = c(1.15, 1.15, 1.15, 1.15, 0.55))
  
  panel_style()
  plot(
    bugenv[bugenv$day == 1, "dist"],
    bugenv[bugenv$day == 1, "toc"],
    type = "b",
    ylim = c(0, 10.5),
    xlim = range(toc_base$dist),
    pch = 19,
    cex = 1.0,
    xlab = "",
    ylab = "TOC (mg/L)",
    lwd = 2.8,
    lty = 1
  )
  points(bugenv[bugenv$day == 126, "dist"], bugenv[bugenv$day == 126, "toc"], type = "b", col = "gold2",  pch = 15, lwd = 2.8, cex = 1.0,  lty = 2)
  points(bugenv[bugenv$day == 260, "dist"], bugenv[bugenv$day == 260, "toc"], type = "b", col = "blue",   pch = 8,  lwd = 2.8, cex = 1.05, lty = 3)
  points(bugenv[bugenv$day == 336, "dist"], bugenv[bugenv$day == 336, "toc"], type = "b", col = "green3", pch = 17, lwd = 2.8, cex = 1.0,  lty = 4)
  points(bugenv[bugenv$day == 518, "dist"], bugenv[bugenv$day == 518, "toc"], type = "b", col = "red",    pch = 18, lwd = 2.8, cex = 1.0,  lty = 5)
  abline(v = 4, lty = 2, lwd = 1.2)
  box(lwd = 1.2)
  draw_panel_tag("(a)")
  
  panel_style()
  plot(
    toc_half[toc_half$day == 1, "dist"],
    toc_half[toc_half$day == 1, "median"],
    type = "n",
    xlim = range(toc_base$dist),
    ylim = c(0, 10.5),
    xlab = "",
    ylab = "TOC (mg/L)"
  )
  for (i in seq_along(obs_days)) {
    draw_band_lines(toc_half[toc_half$day == obs_days[i], ], cols_show[i], ltys_show[i], pchs_show[i], point_cex = cexs_show[i])
  }
  abline(v = 4, lty = 2, lwd = 1.2)
  box(lwd = 1.2)
  draw_panel_tag("(b)")
  
  panel_style()
  plot(
    toc_tenth[toc_tenth$day == 1, "dist"],
    toc_tenth[toc_tenth$day == 1, "median"],
    type = "n",
    xlim = range(toc_base$dist),
    ylim = c(0, 10.5),
    xlab = "",
    ylab = "TOC (mg/L)"
  )
  for (i in seq_along(obs_days)) {
    draw_band_lines(toc_tenth[toc_tenth$day == obs_days[i], ], cols_show[i], ltys_show[i], pchs_show[i], point_cex = cexs_show[i])
  }
  abline(v = 4, lty = 2, lwd = 1.2)
  box(lwd = 1.2)
  draw_panel_tag("(c)")
  
  panel_style()
  plot(
    toc_hundredth[toc_hundredth$day == 1, "dist"],
    toc_hundredth[toc_hundredth$day == 1, "median"],
    type = "n",
    xlim = range(toc_base$dist),
    ylim = c(0, 10.5),
    xlab = "Distance(km)",
    ylab = "TOC (mg/L)"
  )
  for (i in seq_along(obs_days)) {
    draw_band_lines(toc_hundredth[toc_hundredth$day == obs_days[i], ], cols_show[i], ltys_show[i], pchs_show[i], point_cex = cexs_show[i])
  }
  abline(v = 4, lty = 2, lwd = 1.2)
  box(lwd = 1.2)
  draw_panel_tag("(d)")
  
  par(mar = c(0.4, 1.5, 0.2, 1.5))
  plot.new()
  legend(
    "center",
    legend = c("Dec 13", "April 14", "Aug 2014", "Nov 2014", "May 2015"),
    lty = c(1, 2, 3, 4, 5),
    pch = c(19, 15, 8, 17, 18),
    lwd = 2.5,
    pt.cex = c(1.0, 1.0, 1.05, 1.0, 1.0),
    col = c("black", "gold2", "blue", "green3", "red"),
    ncol = 5,
    cex = 1.0,
    bty = "o",
    box.lwd = 1.1,
    x.intersp = 0.8,
    y.intersp = 1.0
  )
  
  dev.off()
}

plot_figure_3(outfile_tif)

#------------------------------------------------------------
# Figure 4. Counterfactual predictions of chlorophyll a
#------------------------------------------------------------

chla_days_obs  <- c(1, 126, 260, 336, 518)
chla_days_show <- c(126, 260, 336, 518)
chla_cols_obs  <- c("black", "gold2", "blue", "green3", "red")
chla_cols_show <- c("gold2", "blue", "green3", "red")
chla_ltys_obs  <- c(1, 2, 3, 4, 5)
chla_ltys_show <- c(2, 3, 4, 5)
chla_pchs_obs  <- c(19, 15, 8, 17, 18)
chla_pchs_show <- c(15, 8, 17, 18)
chla_cexs_obs  <- c(0.95, 1.0, 1.0, 1.0, 1.0)
chla_cexs_show <- c(1.0, 1.0, 1.0, 1.0)

iso_vars <- c("chla", "solar", "temp", "tp", "no3", "turb", "dchla", "dayflow", "eff", "dist", "time")
iso_base <- alldata[complete.cases(alldata[, iso_vars]), ]
iso_base <- iso_base[order(iso_base$day, iso_base$dist), ]
iso_base$time <- factor(iso_base$time, levels = levels(alldata$time))

m_chla_iso <- lm(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)),
  data = iso_base
)

fit_residual_boot_iso <- function(base_data, m_chla_iso) {
  boot_data <- base_data
  boot_data$time <- factor(boot_data$time, levels = levels(base_data$time))
  
  fitted_logchla <- predict(m_chla_iso, newdata = boot_data)
  res_chla <- sample(residuals(m_chla_iso), size = nrow(boot_data), replace = TRUE)
  boot_data$chla <- exp(fitted_logchla + res_chla)
  
  boot_data$time <- factor(boot_data$time, levels = levels(base_data$time))
  
  lm(
    log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)),
    data = boot_data
  )
}

predict_iso_chla <- function(reduction_factor, base_data, m_chla_b) {
  newdata <- base_data
  newdata$time <- factor(newdata$time, levels = levels(base_data$time))
  newdata$dchla <- base_data$dchla * reduction_factor
  exp(predict(m_chla_b, newdata = newdata))
}

simulate_iso_chla_residual <- function(reduction_factor, base_data, m_chla_iso, n_boot = 500L) {
  chla_draws <- matrix(NA_real_, nrow = nrow(base_data), ncol = n_boot)
  
  for (b in seq_len(n_boot)) {
    res <- try({
      m_chla_b <- fit_residual_boot_iso(base_data, m_chla_iso)
      chla_draws[, b] <- predict_iso_chla(reduction_factor, base_data, m_chla_b)
    }, silent = TRUE)
    
    if (inherits(res, "try-error")) chla_draws[, b] <- NA_real_
  }
  
  out <- base_data[, c("day", "dist")]
  out$lower  <- apply(chla_draws, 1, quantile, probs = 0.025, na.rm = TRUE)
  out$median <- apply(chla_draws, 1, median,   na.rm = TRUE)
  out$upper  <- apply(chla_draws, 1, quantile, probs = 0.975, na.rm = TRUE)
  out
}

joint_vars <- c(
  "chla", "tp", "alk", "no3", "rain3", "dayflow", "dalk", "temp", "rain2",
  "dno3", "ph2", "dtp", "solar", "turb", "dchla", "time", "dist", "eff"
)

joint_base <- alldata[complete.cases(alldata[, joint_vars]), ]
joint_base <- joint_base[order(joint_base$day, joint_base$dist), ]
joint_base$time <- factor(joint_base$time, levels = levels(alldata$time))

m_alk_joint <- lm(
  log(alk) ~ rain3 + dayflow +
    eff:I(log(dalk / dayflow)):time +
    eff:I(log(dalk / dayflow)):I(dist == 4.08):time,
  data = joint_base
)

m_no3_joint <- lm(
  log(no3) ~ temp + rain2 + rain3 +
    eff:I(log(dno3) / dayflow):time +
    eff:I(log(dno3) / dayflow):I(dist == 4.08):time,
  data = joint_base
)

m_tp_joint <- lm(
  tp ~ alk + ph2 + rain2 + rain3 +
    eff:I(dtp / dayflow):time +
    eff:I(dtp / dayflow):I(dist == 4.08):time,
  data = joint_base
)

m_chla_joint <- lm(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)),
  data = joint_base
)

fit_residual_boot_joint <- function(base_data, m_alk_joint, m_no3_joint, m_tp_joint, m_chla_joint) {
  boot_data <- base_data
  boot_data$time <- factor(boot_data$time, levels = levels(base_data$time))
  
  fitted_logalk <- predict(m_alk_joint, newdata = boot_data)
  res_alk <- sample(residuals(m_alk_joint), size = nrow(boot_data), replace = TRUE)
  boot_data$alk <- exp(fitted_logalk + res_alk)
  
  fitted_logno3 <- predict(m_no3_joint, newdata = boot_data)
  res_no3 <- sample(residuals(m_no3_joint), size = nrow(boot_data), replace = TRUE)
  boot_data$no3 <- exp(fitted_logno3 + res_no3)
  
  fitted_tp <- predict(m_tp_joint, newdata = boot_data)
  res_tp <- sample(residuals(m_tp_joint), size = nrow(boot_data), replace = TRUE)
  boot_data$tp <- fitted_tp + res_tp
  
  fitted_logchla <- predict(m_chla_joint, newdata = boot_data)
  res_chla <- sample(residuals(m_chla_joint), size = nrow(boot_data), replace = TRUE)
  boot_data$chla <- exp(fitted_logchla + res_chla)
  
  boot_data$time <- factor(boot_data$time, levels = levels(base_data$time))
  
  m_alk_b <- lm(
    log(alk) ~ rain3 + dayflow +
      eff:I(log(dalk / dayflow)):time +
      eff:I(log(dalk / dayflow)):I(dist == 4.08):time,
    data = boot_data
  )
  
  m_no3_b <- lm(
    log(no3) ~ temp + rain2 + rain3 +
      eff:I(log(dno3) / dayflow):time +
      eff:I(log(dno3) / dayflow):I(dist == 4.08):time,
    data = boot_data
  )
  
  m_tp_b <- lm(
    tp ~ alk + ph2 + rain2 + rain3 +
      eff:I(dtp / dayflow):time +
      eff:I(dtp / dayflow):I(dist == 4.08):time,
    data = boot_data
  )
  
  m_chla_b <- lm(
    log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)),
    data = boot_data
  )
  
  list(m_alk = m_alk_b, m_no3 = m_no3_b, m_tp = m_tp_b, m_chla = m_chla_b)
}

predict_joint_chla <- function(reduction_factor, base_data, m_alk_b, m_no3_b, m_tp_b, m_chla_b) {
  newdata <- base_data
  newdata$time <- factor(newdata$time, levels = levels(base_data$time))
  newdata$dalk  <- base_data$dalk  * reduction_factor
  newdata$dno3  <- base_data$dno3  * reduction_factor
  newdata$dtp   <- base_data$dtp   * reduction_factor
  newdata$dchla <- base_data$dchla * reduction_factor
  
  alkpred <- exp(predict(m_alk_b, newdata = newdata))
  newdata$alk <- alkpred
  
  no3pred <- exp(predict(m_no3_b, newdata = newdata))
  newdata$no3 <- no3pred
  
  tppred <- predict(m_tp_b, newdata = newdata)
  newdata$tp <- tppred
  
  exp(predict(m_chla_b, newdata = newdata))
}

simulate_joint_chla_residual <- function(reduction_factor, base_data, m_alk_joint, m_no3_joint, m_tp_joint, m_chla_joint, n_boot = 500L) {
  chla_draws <- matrix(NA_real_, nrow = nrow(base_data), ncol = n_boot)
  
  for (b in seq_len(n_boot)) {
    res <- try({
      boot_models <- fit_residual_boot_joint(base_data, m_alk_joint, m_no3_joint, m_tp_joint, m_chla_joint)
      chla_draws[, b] <- predict_joint_chla(
        reduction_factor = reduction_factor,
        base_data = base_data,
        m_alk_b = boot_models$m_alk,
        m_no3_b = boot_models$m_no3,
        m_tp_b = boot_models$m_tp,
        m_chla_b = boot_models$m_chla
      )
    }, silent = TRUE)
    
    if (inherits(res, "try-error")) chla_draws[, b] <- NA_real_
  }
  
  out <- base_data[, c("day", "dist")]
  out$lower  <- apply(chla_draws, 1, quantile, probs = 0.025, na.rm = TRUE)
  out$median <- apply(chla_draws, 1, median,   na.rm = TRUE)
  out$upper  <- apply(chla_draws, 1, quantile, probs = 0.975, na.rm = TRUE)
  out
}

iso_unchanged  <- simulate_iso_chla_residual(1.00, iso_base,   m_chla_iso,   n_boot = n_boot)
iso_half       <- simulate_iso_chla_residual(0.50, iso_base,   m_chla_iso,   n_boot = n_boot)
iso_tenth      <- simulate_iso_chla_residual(0.10, iso_base,   m_chla_iso,   n_boot = n_boot)
iso_hundredth  <- simulate_iso_chla_residual(0.01, iso_base,   m_chla_iso,   n_boot = n_boot)

joint_half      <- simulate_joint_chla_residual(0.50, joint_base, m_alk_joint, m_no3_joint, m_tp_joint, m_chla_joint, n_boot = n_boot)
joint_tenth     <- simulate_joint_chla_residual(0.10, joint_base, m_alk_joint, m_no3_joint, m_tp_joint, m_chla_joint, n_boot = n_boot)
joint_hundredth <- simulate_joint_chla_residual(0.01, joint_base, m_alk_joint, m_no3_joint, m_tp_joint, m_chla_joint, n_boot = n_boot)

outfile_tif <- file.path(figure_dir, "Figure_4_counterfactual_chlorophyll_a.tif")

plot_observed_panel_chla <- function(tag, xlim_range, ylim_range, xlab_text = "") {
  panel_style()
  
  plot(
    bugenv[bugenv$day == chla_days_obs[1], "dist"],
    bugenv[bugenv$day == chla_days_obs[1], "chla"],
    type = "b",
    xlim = xlim_range,
    ylim = ylim_range,
    pch = chla_pchs_obs[1],
    cex = chla_cexs_obs[1],
    col = chla_cols_obs[1],
    lwd = 2.2,
    lty = chla_ltys_obs[1],
    xlab = xlab_text,
    ylab = "Chlorophyll a (mg/L)"
  )
  
  for (i in 2:length(chla_days_obs)) {
    dd <- chla_days_obs[i]
    points(
      bugenv[bugenv$day == dd, "dist"],
      bugenv[bugenv$day == dd, "chla"],
      type = "b",
      pch = chla_pchs_obs[i],
      cex = chla_cexs_obs[i],
      col = chla_cols_obs[i],
      lwd = 2.2,
      lty = chla_ltys_obs[i]
    )
  }
  
  abline(v = 4, lty = 2, lwd = 1.0, col = "grey40")
  box(lwd = 1.1)
  draw_panel_tag(tag)
}

plot_scenario_panel_chla <- function(panel_df, tag, xlim_range, ylim_range, xlab_text = "") {
  panel_style()
  
  plot(
    panel_df[panel_df$day == chla_days_show[1], "dist"],
    panel_df[panel_df$day == chla_days_show[1], "median"],
    type = "n",
    xlim = xlim_range,
    ylim = ylim_range,
    xlab = xlab_text,
    ylab = "Chlorophyll a (mg/L)"
  )
  
  for (i in seq_along(chla_days_show)) {
    draw_band_lines(
      panel_df[panel_df$day == chla_days_show[i], ],
      line_col = chla_cols_show[i],
      line_lty = chla_ltys_show[i],
      line_pch = chla_pchs_show[i],
      point_cex = chla_cexs_show[i]
    )
  }
  
  abline(v = 4, lty = 2, lwd = 1.0, col = "grey40")
  box(lwd = 1.1)
  draw_panel_tag(tag)
}

plot_figure_4 <- function(device_file) {
  tiff(
    file = device_file,
    width = 12.5,
    height = 16,
    units = "in",
    pointsize = 12,
    bg = "white",
    res = 800,
    compression = "lzw"
  )
  
  ylim_range <- c(0, 40)
  xlim_range <- range(alldata$dist, na.rm = TRUE)
  
  layout(
    matrix(c(1, 2,
             3, 4,
             5, 6,
             7, 8,
             9, 9), nrow = 5, byrow = TRUE),
    heights = c(1, 1, 1, 1, 0.42),
    widths = c(1, 1)
  )
  
  plot_observed_panel_chla("(a)", xlim_range, ylim_range, "")
  plot_scenario_panel_chla(iso_unchanged,  "(b)", xlim_range, ylim_range, "")
  plot_scenario_panel_chla(iso_half,       "(c)", xlim_range, ylim_range, "")
  plot_scenario_panel_chla(iso_tenth,      "(d)", xlim_range, ylim_range, "")
  plot_scenario_panel_chla(iso_hundredth,  "(e)", xlim_range, ylim_range, "")
  plot_scenario_panel_chla(joint_half,     "(f)", xlim_range, ylim_range, "")
  plot_scenario_panel_chla(joint_tenth,    "(g)", xlim_range, ylim_range, "Distance (km)")
  plot_scenario_panel_chla(joint_hundredth,"(h)", xlim_range, ylim_range, "Distance (km)")
  
  par(mar = c(0.3, 1.0, 0.2, 1.0))
  plot.new()
  legend(
    "center",
    legend = c("Dec 2013", "April 2014", "Aug 2014", "Nov 2014", "May 2015"),
    lty = c(1, 2, 3, 4, 5),
    pch = c(19, 15, 8, 17, 18),
    lwd = 2.4,
    pt.cex = c(0.95, 1.0, 1.0, 1.0, 1.0),
    col = c("black", "gold2", "blue", "green3", "red"),
    ncol = 5,
    cex = 1.0,
    bty = "o",
    box.lwd = 1.1,
    x.intersp = 0.8,
    y.intersp = 1.0
  )
  
  dev.off()
}

plot_figure_4(outfile_tif)

#------------------------------------------------------------
# Figure 5. Counterfactual predictions of total phosphorus
#------------------------------------------------------------

tp_days_obs  <- c(1, 126, 260, 336, 518)
tp_days_show <- c(126, 260, 336, 518)

tp_cols_obs  <- c("black", "gold2", "blue", "green3", "red")
tp_cols_show <- c("gold2", "blue", "green3", "red")

tp_ltys_obs  <- c(1, 2, 3, 4, 5)
tp_ltys_show <- c(2, 3, 4, 5)

tp_pchs_obs  <- c(19, 15, 8, 17, 18)
tp_pchs_show <- c(15, 8, 17, 18)

tp_cexs_obs  <- c(0.95, 1.0, 1.0, 1.0, 1.0)
tp_cexs_show <- c(1.00, 1.0, 1.0, 1.0)

tp_vars <- c("tp", "alk", "ph2", "rain2", "rain3", "dtp", "dayflow", "dalk", "time", "dist", "eff")
tp_base <- alldata[complete.cases(alldata[, tp_vars]), ]
tp_base <- tp_base[order(tp_base$day, tp_base$dist), ]
tp_base$time <- factor(tp_base$time, levels = levels(alldata$time))

m_tp_iso0 <- lm(
  tp ~ alk + ph2 + rain2 + rain3 +
    eff:I(dtp / dayflow):time +
    eff:I(dtp / dayflow):I(dist == 4.08):time,
  data = tp_base
)

m_alk_tp0 <- lm(
  log(alk) ~ rain3 + dayflow +
    eff:I(log(dalk / dayflow)):time +
    eff:I(log(dalk / dayflow)):I(dist == 4.08):time,
  data = tp_base
)

m_tp_joint0 <- lm(
  tp ~ alk + ph2 + rain2 + rain3 +
    eff:I(dtp / dayflow):time +
    eff:I(dtp / dayflow):I(dist == 4.08):time,
  data = tp_base
)

fit_tp_iso_boot <- function(base_data, m_tp0) {
  boot_data <- base_data
  boot_data$time <- factor(boot_data$time, levels = levels(base_data$time))
  
  fitted_tp <- predict(m_tp0, newdata = boot_data)
  res_tp <- sample(residuals(m_tp0), size = nrow(boot_data), replace = TRUE)
  boot_data$tp <- fitted_tp + res_tp
  
  boot_data$time <- factor(boot_data$time, levels = levels(base_data$time))
  
  lm(
    tp ~ alk + ph2 + rain2 + rain3 +
      eff:I(dtp / dayflow):time +
      eff:I(dtp / dayflow):I(dist == 4.08):time,
    data = boot_data
  )
}

predict_tp_iso <- function(reduction_factor, base_data, m_tp_b) {
  newdata <- base_data
  newdata$time <- factor(newdata$time, levels = levels(base_data$time))
  newdata$dtp <- base_data$dtp * reduction_factor
  predict(m_tp_b, newdata = newdata)
}

simulate_tp_iso <- function(reduction_factor, base_data, m_tp0, n_boot = 500L) {
  tp_draws <- matrix(NA_real_, nrow = nrow(base_data), ncol = n_boot)
  
  for (b in seq_len(n_boot)) {
    res <- try({
      m_tp_b <- fit_tp_iso_boot(base_data, m_tp0)
      tp_draws[, b] <- predict_tp_iso(reduction_factor, base_data, m_tp_b)
    }, silent = TRUE)
    
    if (inherits(res, "try-error")) tp_draws[, b] <- NA_real_
  }
  
  out <- base_data[, c("day", "dist")]
  out$lower  <- apply(tp_draws, 1, quantile, probs = 0.025, na.rm = TRUE)
  out$median <- apply(tp_draws, 1, median,   na.rm = TRUE)
  out$upper  <- apply(tp_draws, 1, quantile, probs = 0.975, na.rm = TRUE)
  out
}

fit_tp_joint_boot <- function(base_data, m_alk0, m_tp0) {
  boot_data <- base_data
  boot_data$time <- factor(boot_data$time, levels = levels(base_data$time))
  
  fitted_logalk <- predict(m_alk0, newdata = boot_data)
  res_alk <- sample(residuals(m_alk0), size = nrow(boot_data), replace = TRUE)
  boot_data$alk <- exp(fitted_logalk + res_alk)
  
  fitted_tp <- predict(m_tp0, newdata = boot_data)
  res_tp <- sample(residuals(m_tp0), size = nrow(boot_data), replace = TRUE)
  boot_data$tp <- fitted_tp + res_tp
  
  boot_data$time <- factor(boot_data$time, levels = levels(base_data$time))
  
  m_alk_b <- lm(
    log(alk) ~ rain3 + dayflow +
      eff:I(log(dalk / dayflow)):time +
      eff:I(log(dalk / dayflow)):I(dist == 4.08):time,
    data = boot_data
  )
  
  m_tp_b <- lm(
    tp ~ alk + ph2 + rain2 + rain3 +
      eff:I(dtp / dayflow):time +
      eff:I(dtp / dayflow):I(dist == 4.08):time,
    data = boot_data
  )
  
  list(m_alk = m_alk_b, m_tp = m_tp_b)
}

predict_tp_joint <- function(reduction_factor, base_data, m_alk_b, m_tp_b) {
  newdata <- base_data
  newdata$time <- factor(newdata$time, levels = levels(base_data$time))
  newdata$dalk <- base_data$dalk * reduction_factor
  newdata$dtp  <- base_data$dtp  * reduction_factor
  
  alkpred <- exp(predict(m_alk_b, newdata = newdata))
  newdata$alk <- alkpred
  
  predict(m_tp_b, newdata = newdata)
}

simulate_tp_joint <- function(reduction_factor, base_data, m_alk0, m_tp0, n_boot = 500L) {
  tp_draws <- matrix(NA_real_, nrow = nrow(base_data), ncol = n_boot)
  
  for (b in seq_len(n_boot)) {
    res <- try({
      boot_models <- fit_tp_joint_boot(base_data, m_alk0, m_tp0)
      tp_draws[, b] <- predict_tp_joint(
        reduction_factor = reduction_factor,
        base_data = base_data,
        m_alk_b = boot_models$m_alk,
        m_tp_b  = boot_models$m_tp
      )
    }, silent = TRUE)
    
    if (inherits(res, "try-error")) tp_draws[, b] <- NA_real_
  }
  
  out <- base_data[, c("day", "dist")]
  out$lower  <- apply(tp_draws, 1, quantile, probs = 0.025, na.rm = TRUE)
  out$median <- apply(tp_draws, 1, median,   na.rm = TRUE)
  out$upper  <- apply(tp_draws, 1, quantile, probs = 0.975, na.rm = TRUE)
  out
}

n_boot <- 500L

tp_unchanged       <- simulate_tp_iso(1.00, tp_base, m_tp_iso0,   n_boot = n_boot)
tp_half            <- simulate_tp_iso(0.50, tp_base, m_tp_iso0,   n_boot = n_boot)
tp_tenth           <- simulate_tp_iso(0.10, tp_base, m_tp_iso0,   n_boot = n_boot)
tp_hundredth       <- simulate_tp_iso(0.01, tp_base, m_tp_iso0,   n_boot = n_boot)

tp_joint_half      <- simulate_tp_joint(0.50, tp_base, m_alk_tp0, m_tp_joint0, n_boot = n_boot)
tp_joint_tenth     <- simulate_tp_joint(0.10, tp_base, m_alk_tp0, m_tp_joint0, n_boot = n_boot)
tp_joint_hundredth <- simulate_tp_joint(0.01, tp_base, m_alk_tp0, m_tp_joint0, n_boot = n_boot)

plot_observed_panel_tp <- function(tag, xlim_range, ylim_range, xlab_text = "") {
  panel_style()
  
  plot(
    bugenv[bugenv$day == tp_days_obs[1], "dist"],
    bugenv[bugenv$day == tp_days_obs[1], "tp"],
    type = "b",
    xlim = xlim_range,
    ylim = ylim_range,
    pch = tp_pchs_obs[1],
    cex = tp_cexs_obs[1],
    col = tp_cols_obs[1],
    lwd = 2.2,
    lty = tp_ltys_obs[1],
    xlab = xlab_text,
    ylab = "TP (mg/L)"
  )
  
  for (i in 2:length(tp_days_obs)) {
    dd <- tp_days_obs[i]
    points(
      bugenv[bugenv$day == dd, "dist"],
      bugenv[bugenv$day == dd, "tp"],
      type = "b",
      pch = tp_pchs_obs[i],
      cex = tp_cexs_obs[i],
      col = tp_cols_obs[i],
      lwd = 2.2,
      lty = tp_ltys_obs[i]
    )
  }
  
  abline(v = 4, lty = 2, lwd = 1.0, col = "grey40")
  box(lwd = 1.1)
  draw_panel_tag(tag)
}

plot_scenario_panel_tp <- function(panel_df, tag, xlim_range, ylim_range, xlab_text = "") {
  panel_style()
  
  plot(
    panel_df[panel_df$day == tp_days_show[1], "dist"],
    panel_df[panel_df$day == tp_days_show[1], "median"],
    type = "n",
    xlim = xlim_range,
    ylim = ylim_range,
    xlab = xlab_text,
    ylab = "TP (mg/L)"
  )
  
  for (i in seq_along(tp_days_show)) {
    draw_band_lines(
      panel_df[panel_df$day == tp_days_show[i], ],
      line_col = tp_cols_show[i],
      line_lty = tp_ltys_show[i],
      line_pch = tp_pchs_show[i],
      point_cex = tp_cexs_show[i]
    )
  }
  
  abline(v = 4, lty = 2, lwd = 1.0, col = "grey40")
  box(lwd = 1.1)
  draw_panel_tag(tag)
}

outfile_tif <- file.path(figure_dir, "Figure_5_counterfactual_total_phosphorus.tif")

plot_figure_5 <- function(device_file) {
  tiff(
    file = device_file,
    width = 12.5,
    height = 16,
    units = "in",
    pointsize = 12,
    bg = "white",
    res = 800,
    compression = "lzw"
  )
  
  ylim_range <- c(0, 0.8)
  xlim_range <- range(tp_base$dist)
  
  layout(
    matrix(c(1, 2,
             3, 4,
             5, 6,
             7, 8,
             9, 9), nrow = 5, byrow = TRUE),
    heights = c(1, 1, 1, 1, 0.42),
    widths = c(1, 1)
  )
  
  plot_observed_panel_tp("(a)", xlim_range, ylim_range, "")
  plot_scenario_panel_tp(tp_unchanged,       "(b)", xlim_range, ylim_range, "")
  plot_scenario_panel_tp(tp_half,            "(c)", xlim_range, ylim_range, "")
  plot_scenario_panel_tp(tp_tenth,           "(d)", xlim_range, ylim_range, "")
  plot_scenario_panel_tp(tp_hundredth,       "(e)", xlim_range, ylim_range, "")
  plot_scenario_panel_tp(tp_joint_half,      "(f)", xlim_range, ylim_range, "")
  plot_scenario_panel_tp(tp_joint_tenth,     "(g)", xlim_range, ylim_range, "Distance (km)")
  plot_scenario_panel_tp(tp_joint_hundredth, "(h)", xlim_range, ylim_range, "Distance (km)")
  
  par(mar = c(0.3, 1.0, 0.2, 1.0))
  plot.new()
  legend(
    "center",
    legend = c("Dec 2013", "April 2014", "Aug 2014", "Nov 2014", "May 2015"),
    lty = c(1, 2, 3, 4, 5),
    pch = c(19, 15, 8, 17, 18),
    lwd = 2.4,
    pt.cex = c(0.95, 1.0, 1.0, 1.0, 1.0),
    col = c("black", "gold2", "blue", "green3", "red"),
    ncol = 5,
    cex = 1.0,
    bty = "o",
    box.lwd = 1.1,
    x.intersp = 0.8,
    y.intersp = 1.0
  )
  
  dev.off()
}

plot_figure_5(outfile_tif)

#------------------------------------------------------------
# Figure 6. Counterfactual predictions of nitrate
#------------------------------------------------------------

no3_days_obs <- c(1, 126, 260, 336, 518)

no3_cols_obs <- c("black", "gold2", "blue", "green3", "red")
no3_ltys_obs <- c(1, 2, 3, 4, 5)
no3_pchs_obs <- c(19, 15, 8, 17, 18)
no3_cexs_obs <- c(0.95, 1.0, 1.0, 1.0, 1.0)

no3_vars <- c(
  "day",
  "no3",
  "temp",
  "rain2",
  "rain3",
  "dno3",
  "dayflow",
  "time",
  "dist",
  "eff"
)

no3_base <- alldata[complete.cases(alldata[, no3_vars]), ]
no3_base <- no3_base[order(no3_base$day, no3_base$dist), ]
no3_base$time <- factor(no3_base$time, levels = levels(alldata$time))

no3_observed <- no3_base[, c("day", "dist", "no3")]

m_no30 <- lm(
  log(no3) ~ temp + rain2 + rain3 +
    eff:I(log(dno3) / dayflow):time +
    eff:I(log(dno3) / dayflow):I(dist == 4.08):time,
  data = no3_base
)

fit_no3_boot <- function(base_data, fitted_model) {
  boot_data <- base_data
  boot_data$time <- factor(boot_data$time, levels = levels(base_data$time))
  
  fitted_logno3 <- predict(fitted_model, newdata = boot_data)
  resampled_residuals <- sample(
    residuals(fitted_model),
    size = nrow(boot_data),
    replace = TRUE
  )
  
  boot_data$no3 <- exp(fitted_logno3 + resampled_residuals)
  boot_data$time <- factor(boot_data$time, levels = levels(base_data$time))
  
  lm(
    log(no3) ~ temp + rain2 + rain3 +
      eff:I(log(dno3) / dayflow):time +
      eff:I(log(dno3) / dayflow):I(dist == 4.08):time,
    data = boot_data
  )
}

predict_no3 <- function(reduction_factor, base_data, bootstrap_model) {
  newdata <- base_data
  newdata$time <- factor(newdata$time, levels = levels(base_data$time))
  newdata$dno3 <- base_data$dno3 * reduction_factor
  
  exp(predict(bootstrap_model, newdata = newdata))
}

simulate_no3 <- function(reduction_factor, base_data, fitted_model, n_boot = 500L) {
  no3_draws <- matrix(NA_real_, nrow = nrow(base_data), ncol = n_boot)
  
  for (b in seq_len(n_boot)) {
    result <- try({
      bootstrap_model <- fit_no3_boot(
        base_data = base_data,
        fitted_model = fitted_model
      )
      
      no3_draws[, b] <- predict_no3(
        reduction_factor = reduction_factor,
        base_data = base_data,
        bootstrap_model = bootstrap_model
      )
    }, silent = TRUE)
    
    if (inherits(result, "try-error")) {
      no3_draws[, b] <- NA_real_
    }
  }
  
  output <- base_data[, c("day", "dist")]
  output$lower  <- apply(no3_draws, 1, quantile, probs = 0.025, na.rm = TRUE)
  output$median <- apply(no3_draws, 1, median,   na.rm = TRUE)
  output$upper  <- apply(no3_draws, 1, quantile, probs = 0.975, na.rm = TRUE)
  output
}

n_boot <- 500L

no3_half <- simulate_no3(
  reduction_factor = 0.50,
  base_data = no3_base,
  fitted_model = m_no30,
  n_boot = n_boot
)

no3_tenth <- simulate_no3(
  reduction_factor = 0.10,
  base_data = no3_base,
  fitted_model = m_no30,
  n_boot = n_boot
)

no3_hundredth <- simulate_no3(
  reduction_factor = 0.01,
  base_data = no3_base,
  fitted_model = m_no30,
  n_boot = n_boot
)

plot_observed_panel_no3 <- function(panel_df, tag, xlim_range, ylim_range, xlab_text = "") {
  panel_style()
  
  plot(
    NA_real_,
    NA_real_,
    type = "n",
    xlim = xlim_range,
    ylim = ylim_range,
    xlab = xlab_text,
    ylab = "Nitrate (mg/L)"
  )
  
  for (i in seq_along(no3_days_obs)) {
    df_day <- panel_df[panel_df$day == no3_days_obs[i], ]
    df_day <- df_day[order(df_day$dist), ]
    
    lines(
      x = df_day$dist,
      y = df_day$no3,
      col = no3_cols_obs[i],
      lwd = 2.4,
      lty = no3_ltys_obs[i]
    )
    
    points(
      x = df_day$dist,
      y = df_day$no3,
      col = no3_cols_obs[i],
      pch = no3_pchs_obs[i],
      cex = no3_cexs_obs[i]
    )
  }
  
  abline(v = 4, lty = 2, lwd = 1.0, col = "grey40")
  box(lwd = 1.1)
  draw_panel_tag(tag)
}

plot_scenario_panel_no3 <- function(panel_df, tag, xlim_range, ylim_range, xlab_text = "") {
  panel_style()
  
  plot(
    NA_real_,
    NA_real_,
    type = "n",
    xlim = xlim_range,
    ylim = ylim_range,
    xlab = xlab_text,
    ylab = "Nitrate (mg/L)"
  )
  
  for (i in seq_along(no3_days_obs)) {
    draw_band_lines(
      df_day = panel_df[panel_df$day == no3_days_obs[i], ],
      line_col = no3_cols_obs[i],
      line_lty = no3_ltys_obs[i],
      line_pch = no3_pchs_obs[i],
      point_cex = no3_cexs_obs[i]
    )
  }
  
  abline(v = 4, lty = 2, lwd = 1.0, col = "grey40")
  box(lwd = 1.1)
  draw_panel_tag(tag)
}

outfile_tif <- file.path(
  figure_dir,
  "Figure_6_counterfactual_nitrate.tif"
)

plot_figure_6 <- function(device_file) {
  tiff(
    file = device_file,
    width = 10,
    height = 13.5,
    units = "in",
    pointsize = 12,
    bg = "white",
    res = 800,
    compression = "lzw"
  )
  
  xlim_range <- range(no3_base$dist, na.rm = TRUE)
  ylim_range <- c(0, 0.7)
  
  layout(
    matrix(c(1, 2, 3, 4, 5), ncol = 1),
    heights = c(1, 1, 1, 1, 0.42)
  )
  
  plot_observed_panel_no3(no3_observed, "(a)", xlim_range, ylim_range, "")
  plot_scenario_panel_no3(no3_half, "(b)", xlim_range, ylim_range, "")
  plot_scenario_panel_no3(no3_tenth, "(c)", xlim_range, ylim_range, "")
  plot_scenario_panel_no3(no3_hundredth, "(d)", xlim_range, ylim_range, "Distance(km)")
  
  par(mar = c(0.3, 1.0, 0.2, 1.0))
  plot.new()
  legend(
    "center",
    legend = c("Dec 2013", "April 2014", "Aug 2014", "Nov 2014", "May 2015"),
    lty = c(1, 2, 3, 4, 5),
    pch = c(19, 15, 8, 17, 18),
    lwd = 2.4,
    pt.cex = c(0.95, 1.0, 1.0, 1.0, 1.0),
    col = c("black", "gold2", "blue", "green3", "red"),
    ncol = 5,
    cex = 1.0,
    bty = "o",
    box.lwd = 1.1,
    x.intersp = 0.8,
    y.intersp = 1.0
  )
  
  dev.off()
}

plot_figure_6(outfile_tif)

#------------------------------------------------------------
# Figure 7. Counterfactual predictions of alkalinity
#------------------------------------------------------------

alk_days_obs <- c(1, 126, 260, 336, 518)

alk_cols_obs <- c("black", "gold2", "blue", "green3", "red")
alk_ltys_obs <- c(1, 2, 3, 4, 5)
alk_pchs_obs <- c(19, 15, 8, 17, 18)
alk_cexs_obs <- c(0.95, 1.0, 1.0, 1.0, 1.0)

alk_vars <- c(
  "day",
  "alk",
  "rain3",
  "dayflow",
  "dalk",
  "time",
  "dist",
  "eff"
)

alk_base <- alldata[complete.cases(alldata[, alk_vars]), ]
alk_base <- alk_base[order(alk_base$day, alk_base$dist), ]
alk_base$time <- factor(alk_base$time, levels = levels(alldata$time))

alk_observed <- alk_base[, c("day", "dist", "alk")]

m_alk0 <- lm(
  log(alk) ~ rain3 + dayflow +
    eff:I(log(dalk / dayflow)):time +
    eff:I(log(dalk / dayflow)):I(dist == 4.08):time,
  data = alk_base
)

fit_alk_boot <- function(base_data, fitted_model) {
  boot_data <- base_data
  boot_data$time <- factor(boot_data$time, levels = levels(base_data$time))
  
  fitted_logalk <- predict(fitted_model, newdata = boot_data)
  resampled_residuals <- sample(
    residuals(fitted_model),
    size = nrow(boot_data),
    replace = TRUE
  )
  
  boot_data$alk <- exp(fitted_logalk + resampled_residuals)
  boot_data$time <- factor(boot_data$time, levels = levels(base_data$time))
  
  lm(
    log(alk) ~ rain3 + dayflow +
      eff:I(log(dalk / dayflow)):time +
      eff:I(log(dalk / dayflow)):I(dist == 4.08):time,
    data = boot_data
  )
}

predict_alk <- function(reduction_factor, base_data, bootstrap_model) {
  newdata <- base_data
  newdata$time <- factor(newdata$time, levels = levels(base_data$time))
  newdata$dalk <- base_data$dalk * reduction_factor
  
  exp(predict(bootstrap_model, newdata = newdata))
}

simulate_alk <- function(reduction_factor, base_data, fitted_model, n_boot = 500L) {
  alk_draws <- matrix(NA_real_, nrow = nrow(base_data), ncol = n_boot)
  
  for (b in seq_len(n_boot)) {
    result <- try({
      bootstrap_model <- fit_alk_boot(
        base_data = base_data,
        fitted_model = fitted_model
      )
      
      alk_draws[, b] <- predict_alk(
        reduction_factor = reduction_factor,
        base_data = base_data,
        bootstrap_model = bootstrap_model
      )
    }, silent = TRUE)
    
    if (inherits(result, "try-error")) {
      alk_draws[, b] <- NA_real_
    }
  }
  
  output <- base_data[, c("day", "dist")]
  output$lower  <- apply(alk_draws, 1, quantile, probs = 0.025, na.rm = TRUE)
  output$median <- apply(alk_draws, 1, median,   na.rm = TRUE)
  output$upper  <- apply(alk_draws, 1, quantile, probs = 0.975, na.rm = TRUE)
  output
}

n_boot <- 500L

alk_half <- simulate_alk(
  reduction_factor = 0.50,
  base_data = alk_base,
  fitted_model = m_alk0,
  n_boot = n_boot
)

alk_tenth <- simulate_alk(
  reduction_factor = 0.10,
  base_data = alk_base,
  fitted_model = m_alk0,
  n_boot = n_boot
)

alk_hundredth <- simulate_alk(
  reduction_factor = 0.01,
  base_data = alk_base,
  fitted_model = m_alk0,
  n_boot = n_boot
)

plot_observed_panel_alk <- function(panel_df, tag, xlim_range, ylim_range, xlab_text = "") {
  panel_style()
  
  plot(
    NA_real_,
    NA_real_,
    type = "n",
    xlim = xlim_range,
    ylim = ylim_range,
    xlab = xlab_text,
    ylab = "Alkalinity (mg/L)"
  )
  
  for (i in seq_along(alk_days_obs)) {
    df_day <- panel_df[panel_df$day == alk_days_obs[i], ]
    df_day <- df_day[order(df_day$dist), ]
    
    lines(
      x = df_day$dist,
      y = df_day$alk,
      col = alk_cols_obs[i],
      lwd = 2.4,
      lty = alk_ltys_obs[i]
    )
    
    points(
      x = df_day$dist,
      y = df_day$alk,
      col = alk_cols_obs[i],
      pch = alk_pchs_obs[i],
      cex = alk_cexs_obs[i]
    )
  }
  
  abline(v = 4, lty = 2, lwd = 1.0, col = "grey40")
  box(lwd = 1.1)
  draw_panel_tag(tag)
}

plot_scenario_panel_alk <- function(panel_df, tag, xlim_range, ylim_range, xlab_text = "") {
  panel_style()
  
  plot(
    NA_real_,
    NA_real_,
    type = "n",
    xlim = xlim_range,
    ylim = ylim_range,
    xlab = xlab_text,
    ylab = "Alkalinity (mg/L)"
  )
  
  for (i in seq_along(alk_days_obs)) {
    draw_band_lines(
      df_day = panel_df[panel_df$day == alk_days_obs[i], ],
      line_col = alk_cols_obs[i],
      line_lty = alk_ltys_obs[i],
      line_pch = alk_pchs_obs[i],
      point_cex = alk_cexs_obs[i]
    )
  }
  
  abline(v = 4, lty = 2, lwd = 1.0, col = "grey40")
  box(lwd = 1.1)
  draw_panel_tag(tag)
}

outfile_tif <- file.path(
  figure_dir,
  "Figure_7_counterfactual_alkalinity.tif"
)

plot_figure_7 <- function(device_file) {
  tiff(
    file = device_file,
    width = 10,
    height = 13.5,
    units = "in",
    pointsize = 12,
    bg = "white",
    res = 800,
    compression = "lzw"
  )
  
  xlim_range <- range(alk_base$dist, na.rm = TRUE)
  ylim_range <- c(0, 160)
  
  layout(
    matrix(c(1, 2, 3, 4, 5), ncol = 1),
    heights = c(1, 1, 1, 1, 0.42)
  )
  
  plot_observed_panel_alk(alk_observed, "(a)", xlim_range, ylim_range, "")
  plot_scenario_panel_alk(alk_half, "(b)", xlim_range, ylim_range, "")
  plot_scenario_panel_alk(alk_tenth, "(c)", xlim_range, ylim_range, "")
  plot_scenario_panel_alk(alk_hundredth, "(d)", xlim_range, ylim_range, "Distance(km)")
  
  par(mar = c(0.3, 1.0, 0.2, 1.0))
  plot.new()
  legend(
    "center",
    legend = c("Dec 2013", "April 2014", "Aug 2014", "Nov 2014", "May 2015"),
    lty = c(1, 2, 3, 4, 5),
    pch = c(19, 15, 8, 17, 18),
    lwd = 2.4,
    pt.cex = c(0.95, 1.0, 1.0, 1.0, 1.0),
    col = c("black", "gold2", "blue", "green3", "red"),
    ncol = 5,
    cex = 1.0,
    bty = "o",
    box.lwd = 1.1,
    x.intersp = 0.8,
    y.intersp = 1.0
  )
  
  dev.off()
}

plot_figure_7(outfile_tif)

#------------------------------------------------------------
# Figure 8. Counterfactual predictions of pco
#------------------------------------------------------------

pco_days_obs <- c(1, 126, 260, 336, 518)
pco_days_show <- c(126, 260, 336, 518)

pco_cols_obs <- c("black", "gold2", "blue", "green3", "red")
pco_cols_show <- c("gold2", "blue", "green3", "red")

pco_ltys_obs <- c(2, 1, 3, 4, 5)
pco_ltys_show <- c(1, 3, 4, 5)

pco_pchs_obs <- c(19, 15, 8, 17, 18)
pco_pchs_show <- c(15, 8, 17, 18)

pco_cexs_obs <- c(0.95, 1.0, 1.0, 1.0, 1.0)
pco_cexs_show <- c(1.0, 1.0, 1.0, 1.0)

pco_vars <- c(
  "pco1", "toc", "dtoc", "no3", "dno3", "alk", "dalk", "tp", "dtp",
  "chla", "dchla", "cond", "temp", "rain2", "rain3", "dayflow", "ph2",
  "solar", "turb", "zn", "time", "dist", "eff"
)

pco_base <- alldata[complete.cases(alldata[, pco_vars]), ]
pco_base <- pco_base[order(pco_base$day, pco_base$dist), ]
pco_base$time <- factor(pco_base$time, levels = levels(alldata$time))

m_toc_pco0 <- lm(
  toc ~ dayflow + rain2 + rain3 +
    eff:I(dtoc / dayflow):time +
    eff:I(dtoc / dayflow):I(dist == 4.08):time,
  data = pco_base
)

m_no3_pco0 <- lm(
  log(no3) ~ temp + rain2 + rain3 +
    eff:I(log(dno3) / dayflow):time +
    eff:I(log(dno3) / dayflow):I(dist == 4.08):time,
  data = pco_base
)

m_alk_pco0 <- lm(
  log(alk) ~ rain3 + dayflow +
    eff:I(log(dalk / dayflow)):time +
    eff:I(log(dalk / dayflow)):I(dist == 4.08):time,
  data = pco_base
)

m_tp_pco0 <- lm(
  tp ~ alk + ph2 + rain2 + rain3 +
    eff:I(dtp / dayflow):time +
    eff:I(dtp / dayflow):I(dist == 4.08):time,
  data = pco_base
)

m_chla_pco0 <- lm(
  log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)),
  data = pco_base
)

m_cond_pco0 <- lm(
  log(cond) ~ tp + alk + time + dist:time:eff,
  data = pco_base
)

m_pco10 <- lm(
  pco1 ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla),
  data = pco_base
)

fit_residual_boot_pco_chain <- function(base_data, m_toc0, m_no30, m_alk0, m_tp0, m_chla0, m_cond0, m_pco10) {
  boot_data <- base_data
  boot_data$time <- factor(boot_data$time, levels = levels(base_data$time))
  
  fitted_toc <- predict(m_toc0, newdata = boot_data)
  res_toc <- sample(residuals(m_toc0), size = nrow(boot_data), replace = TRUE)
  boot_data$toc <- fitted_toc + res_toc
  
  fitted_logno3 <- predict(m_no30, newdata = boot_data)
  res_no3 <- sample(residuals(m_no30), size = nrow(boot_data), replace = TRUE)
  boot_data$no3 <- exp(fitted_logno3 + res_no3)
  
  fitted_logalk <- predict(m_alk0, newdata = boot_data)
  res_alk <- sample(residuals(m_alk0), size = nrow(boot_data), replace = TRUE)
  boot_data$alk <- exp(fitted_logalk + res_alk)
  
  fitted_tp <- predict(m_tp0, newdata = boot_data)
  res_tp <- sample(residuals(m_tp0), size = nrow(boot_data), replace = TRUE)
  boot_data$tp <- fitted_tp + res_tp
  
  fitted_logchla <- predict(m_chla0, newdata = boot_data)
  res_chla <- sample(residuals(m_chla0), size = nrow(boot_data), replace = TRUE)
  boot_data$chla <- exp(fitted_logchla + res_chla)
  
  fitted_logcond <- predict(m_cond0, newdata = boot_data)
  res_cond <- sample(residuals(m_cond0), size = nrow(boot_data), replace = TRUE)
  boot_data$cond <- exp(fitted_logcond + res_cond)
  
  fitted_pco1 <- predict(m_pco10, newdata = boot_data)
  res_pco1 <- sample(residuals(m_pco10), size = nrow(boot_data), replace = TRUE)
  boot_data$pco1 <- fitted_pco1 + res_pco1
  
  boot_data$time <- factor(boot_data$time, levels = levels(base_data$time))
  
  m_toc_b <- lm(
    toc ~ dayflow + rain2 + rain3 +
      eff:I(dtoc / dayflow):time +
      eff:I(dtoc / dayflow):I(dist == 4.08):time,
    data = boot_data
  )
  
  m_no3_b <- lm(
    log(no3) ~ temp + rain2 + rain3 +
      eff:I(log(dno3) / dayflow):time +
      eff:I(log(dno3) / dayflow):I(dist == 4.08):time,
    data = boot_data
  )
  
  m_alk_b <- lm(
    log(alk) ~ rain3 + dayflow +
      eff:I(log(dalk / dayflow)):time +
      eff:I(log(dalk / dayflow)):I(dist == 4.08):time,
    data = boot_data
  )
  
  m_tp_b <- lm(
    tp ~ alk + ph2 + rain2 + rain3 +
      eff:I(dtp / dayflow):time +
      eff:I(dtp / dayflow):I(dist == 4.08):time,
    data = boot_data
  )
  
  m_chla_b <- lm(
    log(chla) ~ solar + temp + tp + no3 + turb + solar:temp + eff:log(I(dchla / dayflow)),
    data = boot_data
  )
  
  m_cond_b <- lm(
    log(cond) ~ tp + alk + time + dist:time:eff,
    data = boot_data
  )
  
  m_pco1_b <- lm(
    pco1 ~ temp + cond + temp:cond + log(zn) + toc + temp:toc + log(chla),
    data = boot_data
  )
  
  list(
    m_toc = m_toc_b,
    m_no3 = m_no3_b,
    m_alk = m_alk_b,
    m_tp = m_tp_b,
    m_chla = m_chla_b,
    m_cond = m_cond_b,
    m_pco1 = m_pco1_b
  )
}

predict_pco_chain <- function(reduction_factor, base_data, m_toc_b, m_no3_b, m_alk_b, m_tp_b, m_chla_b, m_cond_b, m_pco1_b) {
  newdata <- base_data
  newdata$time <- factor(newdata$time, levels = levels(base_data$time))
  
  newdata$dtoc  <- base_data$dtoc  * reduction_factor
  newdata$dno3  <- base_data$dno3  * reduction_factor
  newdata$dalk  <- base_data$dalk  * reduction_factor
  newdata$dtp   <- base_data$dtp   * reduction_factor
  newdata$dchla <- base_data$dchla * reduction_factor
  
  newdata$toc <- predict(m_toc_b, newdata = newdata)
  newdata$no3 <- exp(predict(m_no3_b, newdata = newdata))
  newdata$alk <- exp(predict(m_alk_b, newdata = newdata))
  newdata$tp  <- predict(m_tp_b, newdata = newdata)
  newdata$chla <- exp(predict(m_chla_b, newdata = newdata))
  newdata$cond <- exp(predict(m_cond_b, newdata = newdata))
  
  predict(m_pco1_b, newdata = newdata)
}

simulate_pco1_residual <- function(reduction_factor, base_data, m_toc0, m_no30, m_alk0, m_tp0, m_chla0, m_cond0, m_pco10, n_boot = 500L) {
  pco_draws <- matrix(NA_real_, nrow = nrow(base_data), ncol = n_boot)
  
  for (b in seq_len(n_boot)) {
    res <- try({
      boot_models <- fit_residual_boot_pco_chain(
        base_data = base_data,
        m_toc0 = m_toc0,
        m_no30 = m_no30,
        m_alk0 = m_alk0,
        m_tp0 = m_tp0,
        m_chla0 = m_chla0,
        m_cond0 = m_cond0,
        m_pco10 = m_pco10
      )
      
      pco_draws[, b] <- predict_pco_chain(
        reduction_factor = reduction_factor,
        base_data = base_data,
        m_toc_b = boot_models$m_toc,
        m_no3_b = boot_models$m_no3,
        m_alk_b = boot_models$m_alk,
        m_tp_b = boot_models$m_tp,
        m_chla_b = boot_models$m_chla,
        m_cond_b = boot_models$m_cond,
        m_pco1_b = boot_models$m_pco1
      )
    }, silent = TRUE)
    
    if (inherits(res, "try-error")) {
      pco_draws[, b] <- NA_real_
    }
  }
  
  out <- base_data[, c("day", "dist")]
  out$lower  <- apply(pco_draws, 1, quantile, probs = 0.025, na.rm = TRUE)
  out$median <- apply(pco_draws, 1, median,   na.rm = TRUE)
  out$upper  <- apply(pco_draws, 1, quantile, probs = 0.975, na.rm = TRUE)
  out
}

n_boot <- 500L

pco_half <- simulate_pco1_residual(
  reduction_factor = 0.50,
  base_data = pco_base,
  m_toc0 = m_toc_pco0,
  m_no30 = m_no3_pco0,
  m_alk0 = m_alk_pco0,
  m_tp0 = m_tp_pco0,
  m_chla0 = m_chla_pco0,
  m_cond0 = m_cond_pco0,
  m_pco10 = m_pco10,
  n_boot = n_boot
)

pco_tenth <- simulate_pco1_residual(
  reduction_factor = 0.10,
  base_data = pco_base,
  m_toc0 = m_toc_pco0,
  m_no30 = m_no3_pco0,
  m_alk0 = m_alk_pco0,
  m_tp0 = m_tp_pco0,
  m_chla0 = m_chla_pco0,
  m_cond0 = m_cond_pco0,
  m_pco10 = m_pco10,
  n_boot = n_boot
)

pco_hundredth <- simulate_pco1_residual(
  reduction_factor = 0.01,
  base_data = pco_base,
  m_toc0 = m_toc_pco0,
  m_no30 = m_no3_pco0,
  m_alk0 = m_alk_pco0,
  m_tp0 = m_tp_pco0,
  m_chla0 = m_chla_pco0,
  m_cond0 = m_cond_pco0,
  m_pco10 = m_pco10,
  n_boot = n_boot
)

plot_observed_panel_pco <- function(tag, xlim_range, ylim_range, xlab_text = "") {
  panel_style()
  
  plot(
    bugenv[bugenv$day == pco_days_obs[1], "dist"],
    bugenv[bugenv$day == pco_days_obs[1], "pco1"],
    type = "b",
    xlim = xlim_range,
    ylim = ylim_range,
    pch = pco_pchs_obs[1],
    cex = pco_cexs_obs[1],
    col = pco_cols_obs[1],
    lwd = 2.8,
    lty = pco_ltys_obs[1],
    xlab = xlab_text,
    ylab = "PCO1"
  )
  
  for (i in 2:length(pco_days_obs)) {
    dd <- pco_days_obs[i]
    points(
      bugenv[bugenv$day == dd, "dist"],
      bugenv[bugenv$day == dd, "pco1"],
      type = "b",
      col = pco_cols_obs[i],
      pch = pco_pchs_obs[i],
      cex = pco_cexs_obs[i],
      lwd = 2.8,
      lty = pco_ltys_obs[i]
    )
  }
  
  abline(v = 4, lty = 2, lwd = 1.2)
  box(lwd = 1.2)
  draw_panel_tag("(a)")
}

plot_scenario_panel_pco <- function(panel_df, tag, xlim_range, ylim_range, xlab_text = "") {
  panel_style()
  
  plot(
    panel_df[panel_df$day == pco_days_show[1], "dist"],
    panel_df[panel_df$day == pco_days_show[1], "median"],
    type = "n",
    xlim = xlim_range,
    ylim = ylim_range,
    xlab = xlab_text,
    ylab = "PCO1"
  )
  
  for (i in seq_along(pco_days_show)) {
    draw_band_lines(
      df_day = panel_df[panel_df$day == pco_days_show[i], ],
      line_col = pco_cols_show[i],
      line_lty = pco_ltys_show[i],
      line_pch = pco_pchs_show[i],
      point_cex = pco_cexs_show[i]
    )
  }
  
  abline(v = 4, lty = 2, lwd = 1.2)
  box(lwd = 1.2)
  draw_panel_tag(tag)
}

outfile_tif <- file.path(
  figure_dir,
  "Figure_8_counterfactual_pco.tif"
)

plot_figure_8 <- function(device_file) {
  tiff(
    file = device_file,
    width = 10,
    height = 13.5,
    units = "in",
    pointsize = 12,
    bg = "white",
    res = 800,
    compression = "lzw"
  )
  
  xlim_range <- range(pco_base$dist)
  ylim_range <- c(-0.6, 1.0)
  
  layout(
    matrix(c(1, 2, 3, 4, 5), ncol = 1),
    heights = c(1.15, 1.15, 1.15, 1.15, 0.55)
  )
  
  plot_observed_panel_pco("(a)", xlim_range, ylim_range, "")
  plot_scenario_panel_pco(pco_half, "(b)", xlim_range, ylim_range, "")
  plot_scenario_panel_pco(pco_tenth, "(c)", xlim_range, ylim_range, "")
  plot_scenario_panel_pco(pco_hundredth, "(d)", xlim_range, ylim_range, "Distance(km)")
  
  par(mar = c(0.4, 1.5, 0.2, 1.5))
  plot.new()
  legend(
    "center",
    legend = c("Dec 2013", "April 2014", "Aug 2014", "Nov 2014", "May 2015"),
    lty = c(2, 1, 3, 4, 5),
    pch = c(19, 15, 8, 17, 18),
    lwd = 2.5,
    pt.cex = c(0.95, 1.0, 1.0, 1.0, 1.0),
    col = c("black", "gold2", "blue", "green3", "red"),
    ncol = 5,
    cex = 1.0,
    bty = "o",
    box.lwd = 1.1,
    x.intersp = 0.8,
    y.intersp = 1.0
  )
  
  dev.off()
}

plot_figure_8(outfile_tif)
