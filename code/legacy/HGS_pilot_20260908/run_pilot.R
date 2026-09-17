options(stringsAsFactors=FALSE, survey.lonely.psu='adjust')
invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
.libPaths(c('C:/Users/Public/CodexRLib42Copy',.libPaths()))
library(haven); library(survey); library(splines); library(ggplot2)
out <- normalizePath('新一篇/HGS_pilot_20260908',winslash='/')
raw <- normalizePath('../数据',winslash='/')
wc <- function(x,n) write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='')
wq <- function(x,w,p) {o<-order(x); approx(cumsum(w[o])/sum(w),x[o],xout=p,method='constant',rule=2)$y}
readcycle <- function(cy) {
 d<-read_xpt(file.path(raw,'人口学',paste0('DEMO_',cy,'.xpt')))
 b<-read_xpt(file.path(raw,'身体测量',paste0('BMX_',cy,'.xpt')))
 g<-read_xpt(file.path(raw,'握力',paste0('MGX_',cy,'.xpt')))
 x<-merge(merge(d,b,by='SEQN',all.x=TRUE),g,by='SEQN',all.x=TRUE)
 stopifnot(!anyDuplicated(x$SEQN))
 trials<-grep('^MGXH[12]T[123]$',names(x),value=TRUE)
 stopifnot(length(trials)==6)
 valid<-sapply(trials,function(t) ifelse(!is.na(x[[paste0(t,'E')]]) & x[[paste0(t,'E')]]==1 & x[[t]]>=0,x[[t]],NA_real_))
 mx<-function(a) {z<-apply(a,1,function(r) if(all(is.na(r))) NA_real_ else max(r,na.rm=TRUE));z[!is.finite(z)]<-NA;z}
 x$hgs<-mx(valid);x$n_valid<-rowSums(!is.na(valid))
 x$max_hand1<-mx(valid[,grepl('H1',trials),drop=FALSE]);x$max_hand2<-mx(valid[,grepl('H2',trials),drop=FALSE])
 x$max_right<-ifelse(x$MGATHAND==1,x$max_hand1,ifelse(x$MGATHAND==2,x$max_hand2,NA))
 x$max_left<-ifelse(x$MGATHAND==2,x$max_hand1,ifelse(x$MGATHAND==1,x$max_hand2,NA))
 x$bilateral_sum<-x$max_right+x$max_left;x$bilateral_mean<-x$bilateral_sum/2
 x$dominant_max<-ifelse(x$MGD130==1,x$max_right,ifelse(x$MGD130==2,x$max_left,NA))
 x$age<-as.numeric(x$RIDAGEYR);x$sex<-ifelse(x$RIAGENDR==1,'Male','Female')
 x$height<-as.numeric(x$BMXHT);x$height10<-x$height/10;x$bmi<-as.numeric(x$BMXBMI);x$weight<-as.numeric(x$WTMEC2YR)
 x$cycle<-cy
 x$eligible<-with(x,age>=20 & age<=79 & is.finite(hgs) & hgs>0 & is.finite(height) & height>0 & is.finite(bmi) & bmi>0 & weight>0 & !is.na(SDMVPSU) & !is.na(SDMVSTRA) & MGDEXSTS %in% c(1,2))
 x$eligible[is.na(x$eligible)]<-FALSE
 keep<-unique(c('SEQN','cycle','age','sex','height','BMXWT','bmi','BMXWAIST','weight','SDMVPSU','SDMVSTRA','eligible','hgs','n_valid','max_right','max_left','bilateral_sum','bilateral_mean','dominant_max',grep('^MG',names(x),value=TRUE)))
 x[,keep]
}
x<-rbind(readcycle('G'),readcycle('H')); saveRDS(x,file.path(out,'harmonized_raw_qc.rds'))
wc(do.call(rbind,lapply(split(x,x$cycle),function(a) data.frame(cycle=a$cycle[1],all=nrow(a),age20_79=sum(a$age>=20&a$age<=79),valid_hgs_age20_79=sum(a$age>=20&a$age<=79&is.finite(a$hgs)&a$hgs>0),analysis=sum(a$eligible)))),'sample_flow')
methods<-c('Absolute','BMI','Height2','Allometry','Conditional')
res<-list();qt<-list();curves<-list();models<-list();exps<-list();scores<-list()
for(s in c('Female','Male')) {
 dev<-subset(x,cycle=='G' & eligible & sex==s);val<-subset(x,cycle=='H' & eligible & sex==s)
 af<-lm(log(hgs)~log(height/100),data=dev,weights=weight/mean(weight));b<-coef(af)[2]
 exps[[s]]<-data.frame(sex=s,b=b,n_development=nrow(dev))
 transform_y<-function(a,m) switch(m,Absolute=a$hgs,BMI=a$hgs/a$bmi,Height2=a$hgs/(a$height/100)^2,Allometry=a$hgs/(a$height/100)^b,Conditional=a$hgs)
 breaks<-c(-Inf,wq(val$height,val$weight,c(.2,.4,.6,.8)),Inf)
 val$height10<-val$height/10
 val$q<-cut(val$height,breaks,labels=paste0('Q',1:5),include.lowest=TRUE)
 for(m in methods) {
  dev$y<-transform_y(dev,m);val$y<-transform_y(val,m)
  f<-if(m=='Conditional') y~ns(age,df=4)+ns(height,df=3) else y~ns(age,df=4)
  fit<-lm(f,data=dev,weights=weight/mean(weight));sigma<-sqrt(weighted.mean(residuals(fit)^2,dev$weight))
  val$z<-(val$y-predict(fit,newdata=val))/sigma
  val$p10<-as.numeric(val$z< -1.282);val$p5<-as.numeric(val$z< -1.645)
  full<-subset(x,cycle=='H'); ii<-match(paste(val$SEQN),paste(full$SEQN));full$z<-full$p10<-full$p5<-NA_real_;full$q<-NA_character_
  full$z[ii]<-val$z;full$p10[ii]<-val$p10;full$p5[ii]<-val$p5;full$q[ii]<-as.character(val$q);full$height10<-full$height/10
  design<-svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~weight,data=subset(full,weight>0),nest=TRUE)
  design<-subset(design,!is.na(z))
  slope<-svyglm(z~height10,design);sc<-summary(slope)$coef;ci<-confint(slope)[2,]
  v<-as.matrix(svyvar(~z+height10,design));cor<-v[1,2]/sqrt(v[1,1]*v[2,2])
  rr<-data.frame(sex=s,method=m,n_dev=nrow(dev),n_val=nrow(val),mean_z=weighted.mean(val$z,val$weight),sd_z=sqrt(v[1,1]),median_z=wq(val$z,val$weight,.5),slope10=sc[2,1],slope_lo=ci[1],slope_hi=ci[2],slope_p=sc[2,4],weighted_cor=cor)
  for(t in c('p10','p5')) {
   ff<-svyglm(as.formula(paste(t,'~height10')),design,family=quasibinomial());cc<-confint(ff)[2,]
   rr[[paste0(t,'_overall')]]<-weighted.mean(val[[t]],val$weight)
   rr[[paste0(t,'_OR10')]]<-exp(coef(ff)[2]);rr[[paste0(t,'_ORlo')]]<-exp(cc[1]);rr[[paste0(t,'_ORhi')]]<-exp(cc[2]);rr[[paste0(t,'_trend_p')]]<-summary(ff)$coef[2,4]
  }
  qq<-do.call(rbind,lapply(levels(val$q),function(q) {a<-val[val$q==q,];data.frame(sex=s,method=m,quintile=q,n=nrow(a),height_min=min(a$height),height_max=max(a$height),mean_z=weighted.mean(a$z,a$weight),sd_z=NA_real_,iqr=diff(wq(a$z,a$weight,c(.25,.75))),p10=weighted.mean(a$p10,a$weight),p5=weighted.mean(a$p5,a$weight))}))
  # Weighted descriptive population SD within each height quintile.
  qq$sd_z<-sapply(levels(val$q),function(q) {a<-val[val$q==q,];sqrt(weighted.mean((a$z-weighted.mean(a$z,a$weight))^2,a$weight))})
  rr$p10_drift<-qq$p10[1]-qq$p10[5];rr$p5_drift<-qq$p5[1]-qq$p5[5]
  sf<-svyglm(z~ns(height10,df=3),design);grid<-data.frame(height10=seq(wq(val$height10,val$weight,.01),wq(val$height10,val$weight,.99),length.out=100))
  grid$z<-as.numeric(predict(sf,newdata=grid));grid$sex<-s;grid$method<-m;grid$height<-grid$height10*10
  key<-paste(s,m);res[[key]]<-rr;qt[[key]]<-qq;curves[[key]]<-grid;models[[key]]<-list(fit=fit,sigma=sigma,b=b);scores[[key]]<-data.frame(SEQN=val$SEQN,sex=s,method=m,z=val$z)
 }
}
summary_table<-do.call(rbind,res);quintiles<-do.call(rbind,qt);curve<-do.call(rbind,curves)
summary_table$slope_reduction<-unlist(lapply(split(summary_table,summary_table$sex),function(a) 1-abs(a$slope10)/abs(a$slope10[a$method=='Absolute'])),use.names=FALSE)
wc(summary_table,'summary');wc(quintiles,'height_quintiles');wc(do.call(rbind,exps),'allometric_exponents');wc(do.call(rbind,scores),'frozen_validation_scores');wc(curve,'spline_curves')
saveRDS(models,file.path(out,'frozen_US_2011_2012_models.rds'))
p<-ggplot(curve,aes(height,z,color=method))+geom_hline(yintercept=0,linetype=2,color='grey60')+geom_line(linewidth=.8)+facet_wrap(~sex)+theme_bw(base_size=12)+labs(x='Height (cm)',y='Frozen validation Z (survey-weighted spline)',color='Method')
ggsave(file.path(out,'Figure1_Z_vs_height.png'),p,width=11,height=5,dpi=180)
long<-rbind(transform(quintiles,tail='P10',prevalence=p10),transform(quintiles,tail='P5',prevalence=p5))
p<-ggplot(long,aes(quintile,prevalence*100,color=method,group=method))+geom_line()+geom_point()+facet_grid(tail~sex)+theme_bw(base_size=12)+labs(x='Sex-specific weighted height quintile (validation)',y='Low-strength prevalence (%)',color='Method')
ggsave(file.path(out,'Figure2_low_prevalence.png'),p,width=11,height=7,dpi=180)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
print(summary_table,digits=3)

capture.output(warnings(),file=file.path(out,"warnings.txt"))

