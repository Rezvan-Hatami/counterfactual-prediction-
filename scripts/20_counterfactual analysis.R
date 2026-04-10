# 20_counterfactual_analysis.R
# Author: Rezvan Hatami
# Date: 2026-04-01
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


