options(stringsAsFactors=FALSE,warn=1,survey.lonely.psu='adjust')
invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
.libPaths(c('C:/Users/Public/CodexRLib42Copy',.libPaths()))
suppressPackageStartupMessages({library(gamlss);library(gamlss.dist);library(splines);library(survey);library(ggplot2)})
wq<-function(x,w,p){o<-order(x);cw<-cumsum(w[o])/sum(w);vapply(p,function(pp)x[o[which(cw>=pp)[1]]],numeric(1))}
wsd<-function(x,w)sqrt(weighted.mean((x-weighted.mean(x,w))^2,w))
ns_spec<-function(x,k){a<-ns(x,df=k);list(knots=attr(a,'knots'),Boundary.knots=attr(a,'Boundary.knots'),intercept=FALSE)}
ns_at<-function(x,s)do.call(ns,c(list(x=x),s))
make_spec<-function(d,mu){list(age=ns_spec(d$age,4),height=ns_spec(d$height,if(mu=='M2')4 else 3),sigma_height=ns_spec(d$height,3),mu=mu)}
make_frame<-function(d,sp){
 A<-ns_at(d$age,sp$age);H<-ns_at(d$height,sp$height);SH<-ns_at(d$height,sp$sigma_height)
 colnames(A)<-paste0('A',seq_len(ncol(A)));colnames(H)<-paste0('H',seq_len(ncol(H)));colnames(SH)<-paste0('SH',seq_len(ncol(SH)))
 out<-data.frame(y=d$hgs,A,H,SH,interaction=(d$age-50)*(d$height-165)/100)
 out
}
formulas<-function(sp,scale){
 mn<-c(paste0('A',1:4),paste0('H',seq_len(length(sp$height$knots)+1)),if(sp$mu=='M3')'interaction')
 sn<-switch(scale,S0=character(),S1=paste0('SH',1:3),S2=c(paste0('A',1:4),paste0('SH',1:3)))
 list(mu=reformulate(mn,response='y'),sigma=if(length(sn))reformulate(sn) else ~1)
}
score_model<-function(b,d){
 if(b$type=='traditional'){
  h<-d$hgs
  y<-switch(b$method,Absolute=h,BMI=h/d$bmi,Height2=h/d$height_m^2,Allometry=h/d$height_m^b$b,Allometry_ageadj=h/d$height_m^b$b)
  A<-ns_at(d$age,b$age_spec);X<-cbind(1,A);mu<-as.numeric(X%*%b$coefficients);z<-(y-mu)/b$sigma
  return(list(z=z,mu=mu,sigma=rep(b$sigma,nrow(d)),clipped=rep(FALSE,nrow(d)),logdens=dnorm(y,mu,b$sigma,log=TRUE)))
 }
 f<-make_frame(d,b$spec);mu<-as.numeric(model.matrix(delete.response(terms(b$formulas$mu)),f)%*%b$coefs$mu)
 sigma<-exp(as.numeric(model.matrix(b$formulas$sigma,f)%*%b$coefs$sigma))
 stopifnot(all(is.finite(mu)),all(is.finite(sigma)),all(sigma>0))
 if(b$family!='NO')stopifnot(all(mu>0))
 args<-list(q=d$hgs,mu=mu,sigma=sigma)
 if(!is.null(b$coefs$nu))args$nu<-rep(unname(b$coefs$nu),nrow(d))
 if(!is.null(b$coefs$tau))args$tau<-rep(exp(unname(b$coefs$tau)),nrow(d))
 p<-do.call(get(paste0('p',b$family),asNamespace('gamlss.dist')),args)
 stopifnot(all(is.finite(p)),all(p>=0&p<=1))
 clipped<-p<=1e-12 | p>=1-1e-12;z<-qnorm(pmin(pmax(p,1e-12),1-1e-12))
 args$x<-args$q;args$q<-NULL;args$log<-TRUE
 ld<-do.call(get(paste0('d',b$family),asNamespace('gamlss.dist')),args)
 stopifnot(all(is.finite(ld)))
 list(z=z,mu=mu,sigma=sigma,clipped=clipped,logdens=ld)
}
fit_model<-function(d,family,mu,scale){
 sp<-make_spec(d,mu);fr<-make_frame(d,sp);ff<-formulas(sp,scale);w<-d$weight/mean(d$weight)
 fr$w<-w;init<-lm(ff$mu,fr,weights=w);m0<-as.numeric(fitted(init));s0<-sqrt(weighted.mean(residuals(init)^2,w))
 warning_messages<-character();last_error<-''
 for(attempt in 1:2){
  ans<-tryCatch(withCallingHandlers({
   a<-list(formula=ff$mu,sigma.formula=ff$sigma,nu.formula=~1,tau.formula=~1,family=get(family,asNamespace('gamlss.dist'))(),data=fr,weights=w,mu.start=m0,sigma.start=rep(if(family=='NO')s0 else s0/weighted.mean(d$hgs,w),nrow(d)),control=gamlss.control(n.cyc=if(attempt==1)150 else 400,c.crit=1e-5,trace=FALSE,mu.step=if(attempt==1)1 else .5,sigma.step=if(attempt==1)1 else .5,nu.step=if(attempt==1)1 else .5,tau.step=if(attempt==1)1 else .5))
   if(family!='NO')a$nu.start<-rep(1,nrow(d))
   if(family %in% c('BCT','BCPE'))a$tau.start<-rep(if(family=='BCT')10 else 2,nrow(d))
   g<-do.call(gamlss,a)
   if(!isTRUE(g$converged))stop('GAMLSS nonconvergence')
   b<-list(type='distribution',family=family,mu=mu,scale=scale,spec=sp,formulas=ff,coefs=list(mu=g$mu.coefficients,sigma=g$sigma.coefficients,nu=g$nu.coefficients,tau=g$tau.coefficients),n=nrow(d),df=g$df.fit,GD=g$G.deviance,BIC=g$G.deviance+log(nrow(d))*g$df.fit,GAIC=g$G.deviance+2*g$df.fit,iterations=g$iter,attempt=attempt)
   pred<-score_model(b,d);stopifnot(max(abs(pred$mu-g$mu.fv))<1e-6,max(abs(pred$sigma-g$sigma.fv))<1e-6)
   b
  },warning=function(e){warning_messages<<-c(warning_messages,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e){last_error<<-conditionMessage(e);NULL})
  if(!is.null(ans)){ans$warnings<-paste(unique(warning_messages),collapse=' | ');return(ans)}
 }
 stop(last_error)
}
fit_traditional<-function(d,method){
 b<-NA_real_;age_spec<-ns_spec(d$age,4);A<-ns_at(d$age,age_spec);fr<-data.frame(logh=log(d$hgs),loght=log(d$height_m),A);w<-d$weight/mean(d$weight)
 if(method=='Allometry')b<-unname(coef(lm(logh~loght,fr,weights=w))[2])
 if(method=='Allometry_ageadj')b<-unname(coef(lm(logh~.,fr,weights=w))['loght'])
 y<-switch(method,Absolute=d$hgs,BMI=d$hgs/d$bmi,Height2=d$hgs/d$height_m^2,Allometry=d$hgs/d$height_m^b,Allometry_ageadj=d$hgs/d$height_m^b)
 X<-cbind(1,A);g<-lm.wfit(X,y,w);sig<-sqrt(weighted.mean(g$residuals^2,w))
 list(type='traditional',method=method,b=b,age_spec=age_spec,coefficients=g$coefficients,sigma=sig,n=nrow(d),df=length(g$coefficients)+1)
}
metrics<-function(d,z,edges){
 w<-d$weight;center<-z-weighted.mean(z,w);sd<-wsd(z,w);h<-d$height/10;hc<-h-weighted.mean(h,w)
 slope<-sum(w*hc*center)/sum(w*hc^2);q<-cut(d$height,c(-Inf,edges,Inf),labels=paste0('Q',1:5),include.lowest=TRUE)
 tab<-do.call(rbind,lapply(paste0('Q',1:5),function(lev){i<-q==lev;zz<-z[i];ww<-w[i];data.frame(quintile=lev,n=sum(i),mean_z=weighted.mean(zz,ww),sd_z=wsd(zz,ww),iqr_z=diff(wq(zz,ww,c(.25,.75))),p5=weighted.mean(zz< -1.645,ww),p10=weighted.mean(zz< -1.282,ww),p90=weighted.mean(zz>1.282,ww),p95=weighted.mean(zz>1.645,ww))}))
 loss<-mean(((tab$p10-.1)^2/.09+(tab$p5-.05)^2/.0475)/2)
 overview<-data.frame(n=nrow(d),mean_z=weighted.mean(z,w),sd_z=sd,median_z=wq(z,w,.5),slope_per10cm=slope,weighted_correlation=slope*wsd(h,w)/sd,p5=weighted.mean(z< -1.645,w),p10=weighted.mean(z< -1.282,w),skewness=weighted.mean(center^3,w)/sd^3,excess_kurtosis=weighted.mean(center^4,w)/sd^4-3,tail_loss=loss,sd_min=min(tab$sd_z),sd_max=max(tab$sd_z),sd_range=diff(range(tab$sd_z)),p5_drift_pp=100*(tab$p5[1]-tab$p5[5]),p10_drift_pp=100*(tab$p10[1]-tab$p10[5]),p5_range_pp=100*diff(range(tab$p5)),p10_range_pp=100*diff(range(tab$p10)))
 overview$guardrail<-abs(overview$mean_z)<=.10 & overview$sd_z>=.90 & overview$sd_z<=1.10 & abs(slope)<=.10 & overview$sd_min>=.80 & overview$sd_max<=1.20
 list(summary=overview,quintiles=tab)
}
