invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
root<-getwd();phase<-file.path(dirname(root),'KNHANES_HGS_phase2_20260910');source(file.path(phase,'core.R'))
wc<-function(x,n)write.csv(x,file.path(root,paste0(n,'.csv')),row.names=FALSE,na='',fileEncoding='UTF-8')
old<-file.path(dirname(root),'KNHANES_HGS_pilot_20260909');dev<-readRDS(file.path(old,'development_data.rds'));val<-readRDS(file.path(old,'validation_data.rds'));kr<-rbind(dev,val)
kr$ageband<-cut(kr$age,c(20,30,40,50,60,70,80),right=FALSE,labels=paste0(seq(20,70,10),'-',seq(29,79,10)))
annual<-list();groups<-list();qc<-list();distributions<-list()
for(s in c('Female','Male')){
 ref<-subset(dev,eligible&sex==s);edges<-wq(ref$height,ref$weight,c(.2,.4,.6,.8));full<-kr
 full$heightband<-cut(full$height,c(-Inf,edges,Inf),labels=paste0('Q',1:5),include.lowest=TRUE)
 design<-svydesign(ids=~psu_pool,strata=~strata_pool,weights=~weight,data=subset(full,design_ok),nest=TRUE)
 for(y in 2014:2019){
  des<-subset(design,eligible&sex==s&year==y);dd<-des$variables
  r<-data.frame(sex=s,year=y,n=nrow(dd))
  for(var in c('hgs','height','bmi','age')){a<-svymean(reformulate(var),des);r[[paste0(var,'_mean')]]<-coef(a)[1];r[[paste0(var,'_se')]]<-SE(a)[1];q<-wq(dd[[var]],dd$weight,c(.1,.25,.5,.75,.9));distributions[[paste(s,y,var)]]<-data.frame(sex=s,year=y,variable=var,mean=coef(a)[1],sd=wsd(dd[[var]],dd$weight),p10=q[1],p25=q[2],median=q[3],p75=q[4],p90=q[5])}
  for(var in c('ageband','heightband')){
   lev<-levels(full[[var]]);refgroup<-if(var=='ageband')cut(ref$age,c(20,30,40,50,60,70,80),right=FALSE,labels=lev) else cut(ref$height,c(-Inf,edges,Inf),labels=lev,include.lowest=TRUE)
   pw<-as.numeric(tapply(ref$weight,refgroup,sum));pw<-pw/sum(pw)
   fit<-svyglm(reformulate(var,intercept=FALSE,response='hgs'),des);stopifnot(length(coef(fit))==length(pw),all(table(dd[[var]])>0))
   est<-sum(pw*coef(fit));se<-sqrt(as.numeric(t(pw)%*%vcov(fit)%*%pw));r[[paste0(var,'_standardized_hgs')]]<-est;r[[paste0(var,'_standardized_se')]]<-se
   groups[[paste(s,y,var)]]<-data.frame(sex=s,year=y,standardization=var,group=lev,n=as.numeric(table(dd[[var]])),reference_proportion=pw,year_proportion=as.numeric(tapply(dd$weight,dd[[var]],sum))/sum(dd$weight),hgs=coef(fit))
  }
  annual[[paste(s,y)]]<-r
  a<-subset(full,year==y&sex==s&age>=20&age<=79&design_ok);w<-a$weight
  q<-data.frame(sex=s,year=y,n_age20_79_design=nrow(a),n_eligible=sum(a$eligible),weighted_eligible_pct=100*weighted.mean(a$eligible,w),weighted_missing_hgs_pct=100*weighted.mean(!is.finite(a$hgs),w),six_valid_pct=100*weighted.mean(a$n_valid==6,w),any_valid_pct=100*weighted.mean(a$n_valid>0,w))
  for(nm in grep('^gs_mea_[rl]_[123]$',names(a),value=TRUE)){q[[paste0(nm,'_missing_pct')]]<-100*weighted.mean(is.na(a[[nm]]),w);q[[paste0(nm,'_mean')]]<-weighted.mean(a[[nm]],w,na.rm=TRUE)}
  qc[[paste(s,y)]]<-q
 }
}
ann<-do.call(rbind,annual);wc(ann,'korea_annual_hgs');wc(do.call(rbind,groups),'korea_standardization_groups');wc(do.call(rbind,distributions),'korea_annual_distributions');wc(do.call(rbind,qc),'korea_annual_qc')
long<-rbind(data.frame(ann[,c('sex','year')],method='Raw',hgs=ann$hgs_mean,se=ann$hgs_se),data.frame(ann[,c('sex','year')],method='Age standardized',hgs=ann$ageband_standardized_hgs,se=ann$ageband_standardized_se),data.frame(ann[,c('sex','year')],method='Height standardized',hgs=ann$heightband_standardized_hgs,se=ann$heightband_standardized_se))
p<-ggplot(long,aes(year,hgs,color=method,group=method))+geom_line()+geom_point()+geom_errorbar(aes(ymin=hgs-1.96*se,ymax=hgs+1.96*se),width=.1)+facet_wrap(~sex,scales='free_y')+scale_x_continuous(breaks=2014:2019)+theme_bw(base_size=12)+theme(legend.position='bottom')+labs(title='Korea: annual grip strength',subtitle='Direct standardization to fixed sex-specific 2014-2016 distributions',x='Year',y='Maximum handgrip strength (kg)',caption='Age: 10-year bands; height: development quintiles. Separate adjustments; reference shares treated as fixed. 95% normal intervals.')
ggsave(file.path(root,'korea_annual_hgs.png'),p,width=12,height=6,dpi=180)
# Fully frozen cross-population stress test, no refit or recalibration.
frozen<-file.path(phase,'locked_development_models.rds');hash0<-unname(tools::md5sum(frozen));bundles<-readRDS(frozen);us<-readRDS(file.path(dirname(root),'HGS_pilot_20260908','harmonized_raw_qc.rds'))
us$uid<-paste0('US_',us$SEQN);us$height_m<-us$height/100;us$height10<-us$height/10;us$weight<-us$weight/2
us$strata_pool<-interaction(us$cycle,us$SDMVSTRA,drop=TRUE);us$psu_pool<-interaction(us$cycle,us$SDMVSTRA,us$SDMVPSU,drop=TRUE);us$design_ok<-is.finite(us$weight)&us$weight>0&!is.na(us$psu_pool)&!is.na(us$strata_pool)
stopifnot(!anyDuplicated(us$uid),all(subset(us,eligible)$BMXWT>0),sum(us$eligible)==9334)
results<-list();qt<-list();scores<-list();support<-list();cuts<-list()
for(s in c('Female','Male')){
 a<-subset(us,eligible&sex==s);ref<-subset(dev,eligible&sex==s);edge<-wq(a$height,a$weight,c(.2,.4,.6,.8));cuts[[s]]<-data.frame(sex=s,probability=c(.2,.4,.6,.8),height_cm=edge)
 support[[s]]<-data.frame(sex=s,n=nrow(a),korea_height_min=min(ref$height),korea_height_max=max(ref$height),us_height_min=min(a$height),us_height_max=max(a$height),outside_height_pct=100*weighted.mean(a$height<min(ref$height)|a$height>max(ref$height),a$weight),outside_age_pct=100*weighted.mean(a$age<min(ref$age)|a$age>max(ref$age),a$weight))
 bs<-list(Height2=readRDS(file.path(phase,'candidates',paste0(s,'_Height2.rds')))$model,Allometry_ageadj=bundles[[s]]$age_adjusted_allometry,Gaussian=bundles[[s]]$simplified,BCT=bundles[[s]]$primary)
 for(m in names(bs)){
  pr<-score_model(bs[[m]],a);mm<-metrics(a,pr$z,edge);full<-us;full$z<-NA_real_;full$z[match(a$uid,full$uid)]<-pr$z
  de<-svydesign(ids=~psu_pool,strata=~strata_pool,weights=~weight,data=subset(full,design_ok),nest=TRUE);de<-subset(de,!is.na(z));f<-svyglm(z~height10,de);ci<-confint(f)[2,]
  rr<-cbind(sex=s,method=m,mm$summary);rr$slope_ci_low<-ci[1];rr$slope_ci_high<-ci[2];rr$clipped<-sum(pr$clipped);results[[paste(s,m)]]<-rr;qt[[paste(s,m)]]<-cbind(sex=s,method=m,mm$quintiles)
  scores[[paste(s,m)]]<-data.frame(uid=a$uid,cycle=a$cycle,sex=s,method=m,z=pr$z,clipped=pr$clipped)
 }
}
ss<-do.call(rbind,results);qq<-do.call(rbind,qt);wc(ss,'nhanes_frozen_summary');wc(qq,'nhanes_height_quintiles');wc(do.call(rbind,support),'nhanes_training_support');wc(do.call(rbind,cuts),'nhanes_height_cutpoints');saveRDS(do.call(rbind,scores),file.path(root,'nhanes_frozen_scores.rds'))
stopifnot(unname(tools::md5sum(frozen))==hash0,nrow(ss)==8,nrow(qq)==40,all(ss$clipped==0));wc(data.frame(model_md5=hash0,scored_at=as.character(Sys.time()),refit=FALSE,recalibrated=FALSE,n_US=sum(us$eligible)),'nhanes_freeze_manifest')
ll<-rbind(transform(qq,tail='P10',pct=100*p10,target=10),transform(qq,tail='P5',pct=100*p5,target=5))
p<-ggplot(ll,aes(quintile,pct,color=method,group=method))+geom_line()+geom_point()+geom_hline(data=unique(ll[,c('tail','target')]),aes(yintercept=target),linetype=2)+facet_grid(tail~sex)+theme_bw(base_size=12)+theme(legend.position='bottom')+labs(title='Korea 2014-2016 models frozen in NHANES 2011-2014',x='US sex-specific weighted height quintile',y='Low-strength prevalence (%)')
ggsave(file.path(root,'nhanes_tail_stress_test.png'),p,width=12,height=8,dpi=180)
print(ann,row.names=FALSE);print(ss[,c('sex','method','mean_z','sd_z','slope_per10cm','p10','p5','tail_loss')],row.names=FALSE)
capture.output(sessionInfo(),file=file.path(root,'sessionInfo.txt'));cat('ANNUAL DIAGNOSTICS AND FROZEN US STRESS TEST COMPLETE\n')