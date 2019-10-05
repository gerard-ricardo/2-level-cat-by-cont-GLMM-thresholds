#EC10's for 1 continuous treatment and 1 categorical treatment


#1) Import data
data1 <-  read.table("https://pastebin.com/raw/ev7vYmZ0", header= TRUE,dec=",")  #spe.rec.sur
head(data1)
options(scipen = 999)  # turn off scientific notation

#2)Organising and wrangling
str(data1)
data1$raw.x <- as.numeric(as.character(data1$raw.x))
data1$factor <- as.factor(as.character(data1$factor))
data1$tank <- as.factor(as.character(data1$tank))
data1$suc <- as.numeric(as.character(data1$suc))
data1$tot <- as.numeric(as.character(data1$tot))
data1$obs <- factor(formatC(1:nrow(data1), flag="0", width = 3))# unique tank ID for later on
data1$prop <- data1$suc/data1$tot  #creates a proportion for binomial data
nrow(data1)
str(data1)

#3)Pre-model Data exploration
hist(data1$tot)  #unbalanaced in response lowest grouping factor. GLMs can allow for some unblance in this factor
table(data1$raw.x)  #unblanced in predictor (x)

#3a)Visulise data - plot data split at every factor
library(ggplot2)
p0 = ggplot()+geom_point(data1, mapping = aes(x = raw.x, y = prop),size = (data1$tot*0.3))+facet_wrap(~factor)+scale_x_log10(name ="dep sed")
p0
#is the control healthy/ or representative of field data health?
#how does the control relates to the first treatment?
#are data evenly spaced across the logx scale?
#is the effect linear or non-linear?
#is there any non-monotonic (bell shape) effects?
#are there any outliers
##keep all this information in mind for post-model data exploration later


#4) Fit the model
library(lme4)
library(splines)
md3 <- glmer(cbind(suc,(tot - suc)) ~ scale(raw.x) * factor + (1|obs) ,family = binomial (link = logit),data = data1) 
summary(md3)
library(RVAideMemoire) #GLMM overdispersion test
overdisp.glmer(md3) #Overdispersion for GLMM
#is there an interactive effect?
#is it overdispersed?


#5)Predict the model across a range of x
md3 <- glmer(cbind(suc,(tot - suc)) ~ raw.x * factor + (1|obs) ,family = binomial (link = logit),data = data1) #dont scale predictors to get prediction
range = range(data1$raw.x)
my.predict <- expand.grid(raw.x = seq(range[1], range[2], length = 100),
                          factor    = c("1", '2'))
mm.m <- model.matrix(~raw.x * factor, data = my.predict)
betas = fixef(md3)
eta <- mm.m %*% betas
my.predict$prediction  <- as.vector(exp(eta) / (1 + exp(eta)))
se    <- sqrt(diag(mm.m %*% vcov(md3) %*% t(mm.m)))
my.predict$upper  <- as.vector(exp(eta + 1.96 *se) /(1 + exp(eta  + 1.96 *se)))
my.predict$lower  <- as.vector(exp(eta - 1.96 *se) /(1 + exp(eta  - 1.96 *se)))
my.predict
p1 = ggplot()+geom_point(data1, mapping = aes(x = raw.x, y = prop),size = (data1$tot*0.3))+facet_wrap(~factor)+scale_x_log10(name ="dep sed")
p1 = p1 + geom_line(data = my.predict, aes(x = raw.x, y = prediction), color = 'red', size=1) #mean
p1 = p1 + geom_line(data = my.predict, aes(x = raw.x, y = upper), color = 'grey', size=1) #upper 95% CI
p1 = p1 + geom_line(data = my.predict, aes(x = raw.x, y = lower), color = 'grey', size=1) #lower 95% CI
p1

#6) Getting ECx's
ec = 10 #put in your ecx here
library(dplyr)
group.fac <-my.predict %>% group_by(factor)%>%summarise(estimate = max(prediction))%>%as.data.frame()
top1 = group.fac$estimate[1]  #the modelled control of factor 1
top2 = group.fac$estimate[2]  #the modelled control of factor 2
inhib1 = top1 -((ec/100) * top1)  #an x% decrease from the control for factor 1
inhib2 = top2 -((ec/100) * top2)  #an x% decrease from the control for factor 2
library(VGAM)
eta1 <- logit(inhib1) 
eta2 <- logit(inhib2)
md3.1 <- glmer(cbind(suc,(tot - suc)) ~ raw.x * factor + (1|obs) ,family = binomial (link = logit),data = data1) 
data1$factor <- relevel(data1$factor, ref = "2") #set reference levels for GLMs   
md3.2 <- glmer(cbind(suc,(tot - suc)) ~ raw.x * factor + (1|obs) ,family = binomial (link = logit),data = data1) 
data1$factor <- relevel(data1$factor, ref = "1") #set reference levels for GLMs   
betas1 = fixef(md3.1)[1:2]  #intercept and slope for ref 1
betas2 = fixef(md3.2)[1:2]  #intercept and slope for ref 2
ecx1 <- (eta1 - betas1[1])/betas1[2] 
ecx2 <- (eta2 - betas2[1])/betas2[2]
pd1 <- -cbind(1, ecx1)/betas1[2] 
pd2 <- -cbind(1, ecx2)/betas2[2] 
ff1 = as.matrix(vcov(md3.1)[1:2,1:2])
ff2 = as.matrix(vcov(md3.2)[1:2,1:2])
ec.se1 <- sqrt(((pd1 %*% ff1 )* pd1) %*% c(1, 1))
ec.se2 <- sqrt(((pd2 %*% ff2 )* pd2) %*% c(1, 1))
upper1 = (ecx1+ec.se1*1.96)
lower1 = (ecx1-ec.se1*1.96)
upper2 = (ecx2+ec.se2*1.96)
lower2 = (ecx2-ec.se2*1.96)
ec.df1 = data.frame(ecx1, lower1, upper1)
ec.df2 = data.frame(ecx2, lower2, upper2)
ec.df1  #this is your factor 1 ecx values
ec.df2  #this is your factor 1 ecx values

#Visualising the ECX
ecx.all = bind_cols(data.frame(factor    = c("1", '2')), data.frame(ecx = c(ecx1, ecx2)), data.frame(inhib = c(inhib1, inhib2)))
upper.all = bind_cols(data.frame(factor    = c("1", '2')), data.frame(upper = c(upper1, upper2)), data.frame(inhib = c(inhib1, inhib2)))
lower.all = bind_cols(data.frame(factor    = c("1", '2')), data.frame(lower = c(lower1, lower2)), data.frame(inhib = c(inhib1, inhib2)))
p1 = p1 + geom_point(data = upper.all, aes(x = upper, y = inhib), color = 'red', size=2) 
p1 = p1 + geom_point(data = ecx.all, aes(x = ecx, y = inhib), color = 'red', size=2) 
p1 = p1 + geom_point(data = lower.all, aes(x = lower, y = inhib), color = 'red', size=2)
p1


