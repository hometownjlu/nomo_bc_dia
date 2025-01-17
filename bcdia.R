library(rms)
library(foreign)
library(survival)
library(glmnet)
library(Hmisc)

setwd("D:/R test")

mydata<-read.csv("wbc-data2.csv") 
mydata<-read.csv("wbc-data2-random.csv") 
mydata<-read.csv("lungp.csv") 
lungp.csv
# -P  cla~clt+bn+bc
summary(mydata$BN)
sd(c(mydata$BN))
summary(mydata)
# mydata$BN=impute(mydata$BN,"random")
 mydata$BN=impute(mydata$BN,mean)
#mydata$BN=impute(mydata$BN,median)
#View(mydata)
names(mydata)

summary(mydata)

dd <- datadist(mydata)
options(datadist='dd')
#mydata=mydata[,2-10]
#names(mydata)
#fit<-lrm(cla~id+uncs+uncsh+maad+siecs+bn,data=mydata,x=T,y=T)
#fmla1<-as.formula(cla~clt+uncsh+maad+bn+bc)
#fmla1<-as.formula(cla~clt+bc)
fmla1<-as.formula(Dia~CT+BN+BC+MA)
fmla2<-as.formula(Dia~BC)
fmla3<-as.formula(Dia~BN)
fmla4<-as.formula(Dia~CT)
fmla5<-as.formula(Dia~MA)
fmla1<-as.formula(price~.)
fit<-lrm(fmla1,data=mydata,x=T,y=T)
#pit<-psm(cla~+uncs+uncsh+maad+siecs+bn,data=mydata,dist = 'lognormal')

fit
summary(fit)
nomomodelA <- nomogram(fit,
                       lp=F, 
                       fun=function(x)1/(1+exp(-x)),
                       fun.at=seq(0.1,1,by=0.2),
                       funlabel="Diagnostic possibility")
pdf("bc-nomo.pdf")
plot(nomomodelA,xfrac = 0.15,cex.axis = 1,cex.var = 1.1)
nomomodelA
dev.off()

#dca curve
library(rmda)
modul=decision_curve(fmla1,data=mydata,
                     family=binomial(link='logit'),
                     thresholds=seq(0,1,by=0.01),
                     confidence.intervals=0.95)
modul2=decision_curve(fmla2,data=mydata,
                     family=binomial(link='logit'),
                     thresholds=seq(0,1,by=0.01),
                     confidence.intervals=0.95)

modul3=decision_curve(fmla3,data=mydata,
                      family=binomial(link='logit'),
                      thresholds=seq(0,1,by=0.01),
                      confidence.intervals=0.95)
modul4=decision_curve(fmla4,data=mydata,
                      family=binomial(link='logit'),
                      thresholds=seq(0,1,by=0.01),
                      confidence.intervals=0.95)
modul5=decision_curve(fmla5,data=mydata,
                      family=binomial(link='logit'),
                      thresholds=seq(0,1,by=0.01),
                      confidence.intervals=0.95)
List=list(modul,modul2,modul3,modul4,modul5)
pdf("bc-dca.pdf")
plot_decision_curve(List,
                    #curve.names = "Malignant Prediction Nomogram",
                    curve.names=c('Nomogram','BC','BN','CT','MA'),
                    xlab="Threshold Probability",ylab="Net Benifit",
                    cost.benefit.axis=FALSE,
                    #col=c('red','blue','green',#223344,#334455),
                    confidence.intervals=FALSE,
                    standardize=FALSE)
dev.off()
modul
modul2
modul3
modul4
modul5

# com 1-100;bc 7-100;bn


#内部 随机抽样验证
# modul  the start is the first value of FPR less than 1,
#end is the value of last colom(snB) is less than 0
set.seed(500)
myc=validate(fit,method = "b",B=500,pr=T,dxy=T)
c_index=(myc[1,5]+1)/2
c_index

#1.4-二分类logistic回归-涉及nomogram列线图
#https://www.bilibili.com/video/BV1AV411d7rv
c=rcorrcens(Dia~predict(fit,newdata=mydata),data=mydata)
c[1,1]
c[1,1]-1.96*c[1,4]/2
c[1,1]+1.96*c[1,4]/2
summary(c)

# coxpe<-predict(fit)
# c_index=1-rcorr.cens(coxpe,mydata$Dia)
# c_index
library(ROCR)
pre_rate=predict(fit)
roc1=prediction(pre_rate,mydata$Dia)
roc2=performance(roc1,"tpr","fpr")
auc=performance(roc1,"auc")

auc

pdf("roc-bc.pdf")
plot(roc2,col="blue",xlab="False Positive Rate",ylab="True Positive Rate",)
abline(0,1,lty=2,lwd=3)
dev.off()
summary(auc)





#缂佹ê鍩楅弲顕�鈧艾鍨痪鍨禈


library(forestplot)

#fun.at=c(.01,.05,seq(.1,.9,by=.1),.95,.99)

#c(0.0001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.9999)
#fun.at=seq(0.1,1,by=0.1)
#缂佹ê鍩楅弲顕�鈧艾鍨痪鍨禈

plot(nomomodelA,xfrac=.45)
library(regplot)
# 输出模型参数相关信息
nomomodelA

#<- calibrate(fit,  method = "KM", B = 1000)
 
cal <- calibrate(fit,  method = "boot", B = 1000)
pdf("bc-cali.pdf")
plot(cal, xlab = "Nomogram Predicted Malignant", ylab = "Actual Malignant",main = "Calibration Curve")
dev.off()
#cal <- calibrate(fit,cmethod="KM",method = 'boot',x=T,y=T)
#cal <- calibrate(fit,method ='boot'),x=TRUE,y=TRU
#plot(xlim = c(0,1.0),ylim = C(0,1.0),xlab = "prob dia",ylab = "ob dia")




#cal <- calibrate(fit,method = 'boot',x=TRUE,y=TRUE)
#cal <- calibrate(fit,method ='boot')
#plot(xlim = c(0,1.0),ylim = C(0,1.0),xlab = "prob dia",ylab = "ob dia")
#plot(cal)

#妫灄鍥?
#http://www.iikx.com/news/statistics/1985.html1:2]
rs_forest<-read.csv("slt2.csv")
#View(dat)
names(rs_forest)
library(forestplot)
pdf("bc-forest.pdf")
tiff('bc-forest.tiff',height = 600,width = 700,res= 60)
forestplot(as.matrix(rs_forest[,1:2]), 
           rs_forest$v3, rs_forest$v4, rs_forest$v5, 
           zero = 1, clip = c(0.4,15), graph.pos = 2, 
           xticks = c(0.1,1,10,100), 
           boxsize=0.05, xlog=TRUE)
dev.off()

#https://www.meiwen.com.cn/subject/jjbmyftx.html
#https://www.plob.org/article/22371.html


forestplot(labeltext = as.matrix(rs_forest[,0:2]),
           
           #设置用于文本展示的列，此处我们用数据的前四列作为文本，在图中展示
           
          rs_forest$v2,rs_forest$v3,rs_forest$v4,)
           #mean = rs_forest$V2, #设置均值
           
           #lower = rs_forest$V3, #设置均值的lowlimits限
           
           #upper = rs_forest$V4, #设置均值的uplimits限
           
           #is.summary=c(T,T,T,F,F,T,F,F,T,F,F,T,F,F,F,T,F,F,T,F,F,T,F,F,T,F,F),
           
           #该参数接受一个逻辑向量，用于定义数据中每一行是否是汇总值，若是，则在对应位置设置为TRUE，若否，则设置为FALSE；设置为TRUE的行则以粗体出现
           
           #zero = 1, #设置参照值，此处我们展示的是HR值，故参照值是1，而不是0
           
          # boxsize = 0.4, #设置点估计的方形大小
           
           #lineheight = unit(8,'mm'),#设置图形中的行距
           
           #colgap = unit(2,'mm'),#设置图形中的列间距
           
           #lwd.zero = 2,#设置参考线的粗细
           
           #lwd.ci = 2,#设置区间估计线的粗细
           
           #col=fpColors(box='#458B00',summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
           
           #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
           
          # xlab="The estimates",#设置x轴标签
           
           #lwd.xaxis=2,#设置X轴线的粗细
           
          # lty.ci = "solid",
           
           #graph.pos = 1)#设置森林图的位置，此处设置为4，则出现在第四列




dev.off()

library(regplot)

#glm model 好用，再查什么意思，广义线性回归模型，线性预测

modelC <- glm(fmla1, data = mydata, family = binomial(link="logit"))
summary(modelC)
modelC

pdf("glmmodel.pdf")

regplot(modelC)
dev.off()

cbind(coef= coef(modelC),confint(modelC))
exp(cbind(OR= coef(modelC),confint(modelC)))

mydata$predmodelC<- predict(newdata=mydata,modelC,"response")

regplot(modelC,observation=mydata[1,]) 
regplot(modelC)

Lusurv<-Surv(time=Lung$time,event = Lung$status)
Lufit <- survfit(Lusurv~Lung$sex)
survdiff(Lusurv~Lung$sex)




dev = rawdata[rawdata$category==0,]
vad = rawdata[rawdata$category==1,]

#閹垫挸瀵橀弫鐗堝祦

ddist <- datadist(dev)
options(datadist='ddist')

#閺嬪嫬缂撴稉澶夐嚋閸ョ偛缍婂Ο鈥崇�?

modelA2 <- lrm(MN ~lnuPCX,data=dev)
modelB2 <- lrm(MN ~ageper10 + eGFRper10 + DM,data=dev)
modelC2 <- lrm(MN ~ageper10 + eGFRper10 + DM + lnuPCX,data=dev)

#鐠佸墽鐤嗛崚妤冨殠閸ユ儳寮弫?
#缁楊兛绔寸悰瀹畂delA鐏忚鲸妲搁崚姘logistic閸ョ偛缍婇惃鍕侀崹瀣倳缁夎埇鈧康p闁瀚═rue閹存湉alse閿涘本妲搁崥锔芥▔缁�铏瑰殠閹囶暕濞村娼楅弽鍥风礄linear predictor閿涘绱漟un閺勵垵顩﹂懛顏勭箒鐠佸彞绔存稉顏勫毐閺佸府绱濈�电p鏉╂稖顢戞潪顒佸床閿涘苯鑻熷铏圭彌娑撯偓娑擃亝鏌婇崸鎰垼鏉炴番鈧倹顒濇径鍕皑閻⑩暔ogit閸欐ɑ宕查惃鍕冀閸戣姤鏆熼敍灞界殺lp鏉烆剚宕叉稉鐑樺灉娴狀剛鍟涢幃澶屾畱妞嬪酣娅撳鍌滃芳-閵嗕咖unction(x) 1/(1+exp(-x))鏉╂瑤绔存稉璇х礉閸楀厖濞囬悽鈺nction()閺嬪嫬缂撴稉鈧稉顏囧殰鐎规矮绠熼崙鑺ユ殶閿涘本瀚崣铚傝厬閻ㄥ墡娴犲窅p閻ㄥ嫯瀵栭崶缈犺厬閸欐牕鈧》绱濇禒锝呭弳1/(1+exp(-x))娑擃叀绻嶇粻妞尖偓?
#fun.at閸掓瑦妲哥紒娆愭煀閻ㄥ嫬娼楅弽鍥叡鐠佸墽鐤嗛懠鍐ㄦ纯閵嗕咖unlabel閸掓瑦妲哥紒娆庣瑐闂堛垼娴嗛幑銏犮偨閻ㄥ嫭鏌婇崸鎰垼鏉炵鎹ｆ稉顏勬倳鐎涙绱滵iagnostic possibility閵嗗倸鍙剧�圭偞婀佹禍鍡氱箹閺夆�虫綏閺嶅洩閰遍敍灞肩瑐闂堫晵p闁綁鍣锋稊鐔峰讲娴犮儴顔曟稉绡庨敍灞肩瑝閺勫墽銇氭禍鍡愨偓?

nomomodelA <- nomogram(modelA2,lp=F, 
                       fun=function(x)1/(1+exp(-x)),
                       fun.at=seq(0.1,1,by=0.1),
                       funlabel="Diagnostic possibility")

nomomodelB <- nomogram(modelB2,lp=F, 
                       fun=function(x)1/(1+exp(-x)),
                       fun.at=seq(0.1,1,by=0.1),
                       funlabel="Diagnostic possibility")

nomomodelC <- nomogram(modelC2,lp=F, 
                       fun=function(x)1/(1+exp(-x)),
                       fun.at=seq(0.1,1,by=0.1),
                       funlabel="Diagnostic possibility")


#缂佹ê鍩楅弲顕�鈧艾鍨痪鍨禈

plot(nomomodelA)
plot(nomomodelB)
plot(nomomodelC)

#缂佹ê鍩楁禍銈勭鞍瀵繐鍨痪鍨禈鐎瑰顥婄粙瀣碍閸栧崨nstall.packages("regplot")

library(regplot)

#娴溿倓绨板蹇撳灙缁惧灝娴樿箛鍛淬�忛悽鈺m閸戣姤鏆?

modelC <- glm(MN ~ageper10 + eGFRper10 + DM + lnuPCX, data = dev, family = binomial(link="logit"))
summary(modelC)
regplot(modelC) 



cbind(coef= coef(modelC),confint(modelC))
exp(cbind(OR= coef(modelC),confint(modelC)))

mydata$predmodelC<- predict(newdata=dev,modelC,"response")
regplot(modelC,observation=mydata[10,]) 




#LASSO鍒嗘瀽
v1<-as.matrix(mydata[,c(3:11)])
v2<-mydata[,2]

#myfit<-glmnet(v1,v2,)

myfit = glmnet(v1,v2,family = "binomial")

pdf("lambda.pdf")
plot(myfit,xvar="lambda",label=TRUE)
dev.off()

myfit2=cv.glmnet(v1,v2,family="binomial")
pdf("min.pdf")
plot(myfit2)
abline(v=log(c(myfit2$lambda.min,myfit2$lambda.1se)),lty="dashed")
dev.off()
myfit2$lambda.min

coe=coef(myfit,s=myfit2$lambda.min)
act_index=which(coe!=0)
act_coe= coe[act_index]
row.names(coe)[act_index]
coe
