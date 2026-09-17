invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
root<-getwd();source(file.path(dirname(root),'KNHANES_HGS_phase2_20260910','core.R'))
wc<-function(x,n)write.csv(x,paste0(n,'.csv'),row.names=FALSE,na='',fileEncoding='UTF-8')
old<-file.path(dirname(root),'KNHANES_HGS_pilot_20260909');full<-rbind(readRDS(file.path(old,'development_data.rds')),readRDS(file.path(old,'validation_data.rds')))
# Original files use annual weight / 3; six-year average uses annual weight / 6.
full$weight<-full$weight/2;stopifnot(!anyDuplicated(full$uid),all(full$year %in% 2014:2019),sum(full$eligible)==32311)
models<-list();scores<-list();overall<-list();quints<-list();years<-list();diagnostics<-list();coefficients<-list();worms<-list();cutpoints<-list()
stats<-function(d,z)data.frame(N=nrow(d),mean_z=weighted.mean(z,d$weight),sd_z=wsd(z,d$weight),P10=100*weighted.mean(z< -1.282,d$weight),P5=100*weighted.mean(z< -1.645,d$weight))
for(s in c('Female','Male')){
 d<-subset(full,eligible&sex==s);cat('FITTING',s,nrow(d),'\n');flush.console()
 b<-fit_model(d,'BCT','M1','S2');models[[s]]<-b;p<-score_model(b,d);z<-p$z;stopifnot(all(is.finite(z)),!any(p$clipped),b$df==18)
 edges<-wq(d$height,d$weight,c(.2,.4,.6,.8));q<-cut(d$height,c(-Inf,edges,Inf),labels=paste0('Q',1:5),include.lowest=TRUE);cutpoints[[s]]<-data.frame(sex=s,prob=c(.2,.4,.6,.8),height=edges)
 scores[[s]]<-data.frame(uid=d$uid,sex=s,year=d$year,height=d$height,weight=d$weight,quintile=q,z=z,mu=p$mu,sigma=p$sigma,clipped=p$clipped)
 f<-full;f$z<-NA_real_;f$z[match(d$uid,f$uid)]<-z;f$P10<-as.numeric(f$z< -1.282);f$P5<-as.numeric(f$z< -1.645)
 de<-svydesign(ids=~psu_pool,strata=~strata_pool,weights=~weight,data=subset(f,design_ok),nest=TRUE);de<-subset(de,!is.na(z))
 fit<-svyglm(z~height10,de);ci<-confint(fit)[2,];r<-cbind(sex=s,stats(d,z),slope_10cm=unname(coef(fit)[2]),slope_CI_low=ci[1],slope_CI_high=ci[2])
 for(t in c('P10','P5')){fit<-svyglm(reformulate('height10',response=t),de,family=quasibinomial());ci<-confint(fit)[2,];r[[paste0('OR_',t)]]<-exp(coef(fit)[2]);r[[paste0('OR_',t,'_CI_low')]]<-exp(ci[1]);r[[paste0('OR_',t,'_CI_high')]]<-exp(ci[2])}
 qq<-do.call(rbind,lapply(levels(q),function(lev){i<-q==lev;cbind(sex=s,quintile=lev,stats(d[i,],z[i]),P90=100*weighted.mean(z[i]>1.282,d$weight[i]),P95=100*weighted.mean(z[i]>1.645,d$weight[i]))}))
 r$P10_drift_pp<-qq$P10[1]-qq$P10[5];r$P5_drift_pp<-qq$P5[1]-qq$P5[5]
 for(nm in c('mean_z','sd_z','P10','P5','P90','P95'))r[[paste0(nm,'_quintile_range')]]<-diff(range(qq[[nm]]))
 overall[[s]]<-r;quints[[s]]<-qq
 years[[s]]<-do.call(rbind,lapply(2014:2019,function(y){i<-d$year==y;cbind(sex=s,year=y,stats(d[i,],z[i]))}))
 nu<-unname(b$coefs$nu);tau<-exp(unname(b$coefs$tau));m<-metrics(d,z,edges)$summary
 diagnostics[[s]]<-data.frame(sex=s,converged=TRUE,iterations=b$iterations,attempt=b$attempt,df=b$df,all_coefficients_finite=all(is.finite(unlist(b$coefs))),nu=nu,tau=tau,mu_min=min(p$mu),mu_max=max(p$mu),sigma_min=min(p$sigma),sigma_max=max(p$sigma),clipped=sum(p$clipped),residual_skewness=m$skewness,residual_excess_kurtosis=m$excess_kurtosis,warnings=b$warnings)
 coefficients[[s]]<-do.call(rbind,lapply(names(b$coefs),function(nm)data.frame(sex=s,parameter=nm,term=names(b$coefs[[nm]]),coefficient=unname(b$coefs[[nm]]))))
 pp<-seq(.01,.99,.01)
 for(lev in c('Overall',levels(q))){i<-if(lev=='Overall')rep(TRUE,nrow(d))else q==lev;worms[[paste(s,lev)]]<-data.frame(sex=s,group=lev,theory=qnorm(pp),observed=wq(z[i],d$weight[i],pp),difference=wq(z[i],d$weight[i],pp)-qnorm(pp))}
 cat('DONE',s,'iterations',b$iterations,'\n');flush.console()
}
ss<-do.call(rbind,overall);qq<-do.call(rbind,quints);yy<-do.call(rbind,years);sc<-do.call(rbind,scores);ww<-do.call(rbind,worms)
wc(ss,'pooled_BCT_overall');wc(qq,'pooled_BCT_height_quintiles');wc(yy,'pooled_BCT_by_year');wc(do.call(rbind,diagnostics),'model_diagnostics');wc(do.call(rbind,coefficients),'model_coefficients');wc(do.call(rbind,cutpoints),'height_cutpoints');wc(ww,'worm_plot_data')
saveRDS(models,'pooled_BCT_models.rds');saveRDS(sc,'pooled_BCT_scores.rds')
wc(data.frame(file='pooled_BCT_models.rds',md5=unname(tools::md5sum('pooled_BCT_models.rds')),fit_years='2014-2019',evaluation='same fitting sample',created=as.character(Sys.time())),'model_manifest')
th<-theme_bw(base_size=12)+theme(legend.position='bottom')
p<-ggplot(yy,aes(year,mean_z,color=sex))+geom_hline(yintercept=0,linetype=2)+geom_line()+geom_point(size=2)+scale_x_continuous(breaks=2014:2019)+th+labs(title='Pooled BCT: yearly mean Z',subtitle='One sex-specific 2014-2019 model; in-sample year diagnostics',x='Year',y='Weighted mean Z')
ggsave('yearly_Z.png',p,width=10,height=5.5,dpi=180)
yl<-rbind(transform(yy,tail='P10',pct=P10,target=10),transform(yy,tail='P5',pct=P5,target=5))
p<-ggplot(yl,aes(year,pct,color=sex))+geom_hline(data=unique(yl[,c('tail','target')]),aes(yintercept=target),linetype=2)+geom_line()+geom_point()+facet_wrap(~tail)+scale_x_continuous(breaks=2014:2019)+th+labs(title='Pooled BCT: yearly lower-tail proportions',x='Year',y='Proportion (%)')
ggsave('yearly_P10_P5.png',p,width=11,height=5.5,dpi=180)
p<-ggplot(subset(ww,group=='Overall'),aes(theory,observed))+geom_abline(intercept=0,slope=1,linetype=2)+geom_line()+facet_wrap(~sex)+th+labs(title='Pooled BCT: overall quantile residual QQ',x='Normal quantile',y='Weighted residual quantile')
ggsave('overall_RQR.png',p,width=10,height=5.5,dpi=180)
p<-ggplot(subset(ww,group=='Overall'),aes(theory,difference))+geom_hline(yintercept=0,linetype=2)+geom_line()+facet_wrap(~sex)+th+labs(title='Pooled BCT: overall worm plot',x='Normal quantile',y='Observed minus normal quantile')
ggsave('overall_worm.png',p,width=10,height=5.5,dpi=180)
for(s in c('Female','Male')){p<-ggplot(subset(ww,sex==s&group!='Overall'),aes(theory,difference))+geom_hline(yintercept=0,linetype=2)+geom_line()+facet_wrap(~group,ncol=3)+th+labs(title=paste(s,'pooled BCT: height-quintile worm plots'),x='Normal quantile',y='Observed minus normal quantile',caption='Weighted descriptive detrended QQ; no independence-based confidence bands.')
ggsave(paste0(s,'_height_worm.png'),p,width=12,height=7,dpi=180)}
stopifnot(sum(ss$N)==32311,nrow(qq)==10,nrow(yy)==12,all(is.finite(sc$z)),!anyDuplicated(sc$uid))
for(s in c('Female','Male')){stopifnot(sum(subset(qq,sex==s)$N)==subset(ss,sex==s)$N,sum(subset(yy,sex==s)$N)==subset(ss,sex==s)$N)}
wc(data.frame(check=c('Only BCT M1 S2 fitted with 18 parameters per sex','Same 32311 eligible records; unique IDs','All model predictions finite with no CDF clipping','Height and year tables sum to sex-specific sample sizes','Frozen prediction agrees with GAMLSS fitted location and scale'),pass=TRUE),'integrity_checks')
capture.output(sessionInfo(),file='sessionInfo.txt');print(ss,row.names=FALSE);print(yy,row.names=FALSE);cat('POOLED BCT COMPLETE\n')