library(nycflights13)
library(tidyverse)
library(modelr)
library(glmnetUtils)

fitfun=function(dF){
  cv.glmnet(arr_delay~distance+air_time+dep_time,data=dF)
}
gnetr2=function(model,datavals){
  yvar=all.vars(formula(model)[[2]])
  # print(paste('y variable:',yvar))
  # print('observations')
  # print(str(as.data.frame(datavals)[[yvar]]))
  # print('predictions')
  # print(str(predict(object=model,newdata=datavals)))
  stats::cor(stats::predict(object=model, newdata=datavals), as.data.frame(datavals)[[yvar]], use='complete.obs')^2
}


flights %>%
  group_by(carrier) %>%
  do({
    crossv_mc(.,4) %>%
      mutate(mdl=map(train,fitfun),
             r2=map2_dbl(mdl,test,gnetr2))
  })

f2=flights %>% slice(1:150000)
cl=parallel::makePSOCKcluster(4)
cvalpha=cvAlpha.glmnet(arr_delay~distance+air_time+dep_time,data=flights,outerParallel=cl)
parallel::stopCluster(cl)

cv0=cv.glmnet(arr_delay~dis0tance+air_time+dep_time,data=f2,alpha=0)
predict(cv0,newdata=head(f2))

maxlambda=100
alpha=cvalpha$alpha
ialpha=1
cvdf=data_frame(alpha) %>%
  bind_cols(data.frame(matrix(NA,ncol=maxlambda,nrow=length(alpha)))) %>%
  bind_cols(as.data.frame(matrix(NA,ncol=maxlambda,nrow=length(alpha))))#
for(ialpha in seq_along(alpha)){
  cvlambdas=cvalpha$modlist[[ialpha]]$lambda
  numlambda=length(cvlambdas)
  cvdf[ialpha,2:(numlambda+1)] <- cvlambdas
  cvdf[ialpha,maxlambda+1+(1:numlambda)] <- cvalpha$modlist[[ialpha]]$cvm
}
cvresults <-
  cvdf %>%
  mutate(alpha=alpha) %>%
  gather(lambdaid,lambdaval,num_range('X',1:100)) %>%
  gather(mseid,mseval,num_range('V',1:100)) %>%
  mutate(lambdaid=extract_numeric(lambdaid),
         mseid=extract_numeric(mseid))

ggplot(cvresults)+
  geom_boxplot(aes(x=as.factor(alpha),y=mseval))

ggplot(cvresults)+
  geom_point(aes(x=lambdaid,y=alpha,color=mseval))+
  viridis::scale_color_viridis()

med_mse <-
  cvresults %>%
  group_by(alpha) %>%
  summarise(medmse=median(mseval,na.rm=T))

med_se <- med_mse %>%
  summarise(
    se=sd(medmse,na.rm=T)/sqrt(n())
  ) %>% as.numeric

alphaind=which.min(abs(med_mse$medmse-(min(med_mse$medmse,na.rm=T)+med_se)))
bestalpha=alpha[alphaind]
predict(cvalpha,newdata=head(flights),which=alphaind)

ggplot(cvresults)+
  geom_tile(aes(x=lambdaid,y=alpha,fill=mseval))+
  viridis::scale_fill_viridis()

range(cvresults$mseval,na.rm=T)
median(cvresults$mseval,na.rm=T)
ggplot(cvresults)+
  geom_point(aes(x=lambdaid,y=as.numeric(alpha),fill=mseval))+
  scale_fill_continuous(limits=c(1820,1890))

predict(cvalpha,newdata=head(f2),alpha=0)

coef(cvalpha,which=alphaind)

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)

  dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

im=list(x=alpha,y=1:maxlambda)
image(im,z=Z,col=viridis::plasma(20))
color.bar(viridis::plasma(length(alpha)*maxlambda),1800,1900)

