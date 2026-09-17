invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
root<-getwd();source(file.path(root,'core.R'))
wc<-function(x,n)write.csv(x,file.path(root,paste0(n,'.csv')),row.names=FALSE,na='',fileEncoding='UTF-8')
old<-file.path(dirname(root),'KNHANES_HGS_pilot_20260909')
frozenpath<-file.path(root,'locked_development_models.rds');stopifnot(file.exists(frozenpath))
hash0<-unname(tools::md5sum(frozenpath));manifest<-read.csv(file.path(root,'model_lock_manifest.csv'));stopifnot(hash0==manifest$md5)
bundles<-readRDS(frozenpath);v0<-readRDS(file.path(old,'validation_data.rds'));stopifnot(all(v0$year %in% 2017:2019))
edges<-read.csv(file.path(old,'height_quintile_boundaries.csv'))
scorefile<-file.path(root,'temporal_scores.rds')
if(file.exists(scorefile)){sc<-readRDS(scorefile)}else{
 prior<-readRDS(file.path(dirname(root),'KNHANES_HGS_sigma_comparison_20260910','validation_scores_seven_methods.rds'))
 news<-list()
 for(s in c('Female','Male'))for(role in c('primary','simplified','age_adjusted_allometry')){
  v<-subset(v0,eligible & sex==s);p<-score_model(bundles[[s]][[role]],v)
  q<-cut(v$height,c(-Inf,edges$height_cm[edges$sex==s],Inf),labels=paste0('Q',1:5),include.lowest=TRUE)
  news[[paste(s,role)]]<-data.frame(uid=v$uid,year=v$year,sex=s,method=role,height_cm=v$height,weight=v$weight,quintile=q,z=p$z,p5=as.numeric(p$z< -1.645),p10=as.numeric(p$z< -1.282),p90=as.numeric(p$z>1.282),p95=as.numeric(p$z>1.645))
 }
 sc<-rbind(prior,do.call(rbind,news));saveRDS(sc,scorefile)
 wc(data.frame(scored_at=as.character(Sys.time()),model_md5=hash0,n=nrow(sc)),'temporal_scoring_manifest')
}
methods<-unique(sc$method);rows<-list();quints<-list();curves<-list()
for(s in c('Female','Male')) for(m in methods){
 a<-sc[sc$sex==s & sc$method==m,];full<-v0;idx<-match(a$uid,full$uid)
 stopifnot(!anyNA(idx),!anyDuplicated(a$uid))
 for(n in c('z','p5','p10','p90','p95')){full[[n]]<-NA_real_;full[[n]][idx]<-a[[n]]};full$q<-NA_character_;full$q[idx]<-as.character(a$quintile)
 des<-svydesign(ids=~psu_pool,strata=~strata_pool,weights=~weight,data=subset(full,design_ok),nest=TRUE);des<-subset(des,!is.na(z))
 f<-svyglm(z~height10,des);co<-summary(f)$coef;ci<-confint(f)[2,]
 corr<-weighted.mean((a$z-weighted.mean(a$z,a$weight))*(a$height_cm-weighted.mean(a$height_cm,a$weight)),a$weight)/(wsd(a$z,a$weight)*wsd(a$height_cm,a$weight))
 r<-data.frame(sex=s,method=m,n_development=bundles[[s]]$primary$n,n_validation=nrow(a),mean_z=weighted.mean(a$z,a$weight),sd_z=wsd(a$z,a$weight),median_z=wq(a$z,a$weight,.5),slope_per10cm=co[2,1],slope_ci_low=ci[1],slope_ci_high=ci[2],slope_p=co[2,4],weighted_correlation=corr)
 for(t in c('p5','p10','p90','p95'))r[[paste0(t,'_overall_pct')]]<-100*weighted.mean(a[[t]],a$weight)
 for(t in c('p5','p10')){
  f2<-svyglm(as.formula(paste(t,'~height10')),des,family=quasibinomial());c2<-confint(f2)[2,]
  r[[paste0(t,'_OR_per10cm')]]<-exp(coef(f2)[2]);r[[paste0(t,'_OR_ci_low')]]<-exp(c2[1]);r[[paste0(t,'_OR_ci_high')]]<-exp(c2[2]);r[[paste0(t,'_trend_p')]]<-summary(f2)$coef[2,4]
 }
 qq<-do.call(rbind,lapply(paste0('Q',1:5),function(lev){
  aa<-a[a$quintile==lev,];r0<-data.frame(sex=s,method=m,quintile=lev,n=nrow(aa),weighted_n=sum(aa$weight),mean_z=weighted.mean(aa$z,aa$weight),sd_z=wsd(aa$z,aa$weight),iqr_z=diff(wq(aa$z,aa$weight,c(.25,.75))))
  for(t in c('p5','p10','p90','p95')){r0[[paste0(t,'_pct')]]<-100*weighted.mean(aa[[t]],aa$weight);r0[[paste0(t,'_events')]]<-sum(aa[[t]])}
  r0
 }))
 for(t in c('p5','p10','p90','p95'))r[[paste0(t,'_Q1minusQ5_pp')]]<-qq[[paste0(t,'_pct')]][1]-qq[[paste0(t,'_pct')]][5]
 r$sd_Q1<-qq$sd_z[1];r$sd_Q5<-qq$sd_z[5];r$sd_Q5_over_Q1<-qq$sd_z[5]/qq$sd_z[1]
 r$sd_quintile_range<-diff(range(qq$sd_z));r$p10_quintile_range_pp<-diff(range(qq$p10_pct));r$p5_quintile_range_pp<-diff(range(qq$p5_pct))
 sf<-svyglm(z~ns(height10,df=3),des);g<-data.frame(height10=seq(wq(a$height_cm/10,a$weight,.01),wq(a$height_cm/10,a$weight,.99),length.out=120));g$z<-as.numeric(predict(sf,g));g$height_cm<-g$height10*10;g$sex<-s;g$method<-m
 key<-paste(s,m);rows[[key]]<-r;quints[[key]]<-qq;curves[[key]]<-g
}
ss<-do.call(rbind,rows);qq<-do.call(rbind,quints);gg<-do.call(rbind,curves)
wc(ss,'summary_validation');wc(qq,'height_quintiles');wc(gg,'spline_curve_data')

stopifnot(identical(hash0,unname(tools::md5sum(frozenpath))))
shown<-c('primary','simplified','age_adjusted_allometry','Height2','Conditional_heightSD')
qq<-subset(qq,method %in% shown);gg<-subset(gg,method %in% shown)
th<-theme_bw(base_size=11)+theme(legend.position='bottom',legend.title=element_blank())
p<-ggplot(gg,aes(height_cm,z,color=method))+geom_hline(yintercept=0,linetype=2)+geom_line()+facet_wrap(~sex,scales='free_x')+th+labs(title='Frozen temporal validation: mean dependence',x='Height (cm)',y='Survey spline of Z')
ggsave(file.path(root,'temporal_mean.png'),p,width=12,height=6,dpi=180)
long<-rbind(transform(qq,tail='P10',pct=p10_pct,target=10),transform(qq,tail='P5',pct=p5_pct,target=5))
p<-ggplot(long,aes(quintile,pct,color=method,group=method))+geom_line()+geom_point()+geom_hline(data=unique(long[,c('tail','target')]),aes(yintercept=target),linetype=2)+facet_grid(tail~sex)+th+labs(title='Frozen temporal validation: lower-tail coverage',y='Prevalence (%)',x='Height quintile')
ggsave(file.path(root,'temporal_tails.png'),p,width=12,height=8,dpi=180)
p<-ggplot(qq,aes(quintile,sd_z,color=method,group=method))+geom_line()+geom_point()+geom_hline(yintercept=1,linetype=2)+facet_wrap(~sex)+th+labs(title='Frozen temporal validation: scale',y='Weighted SD(Z)',x='Height quintile')
ggsave(file.path(root,'temporal_scale.png'),p,width=12,height=6,dpi=180)
print(subset(ss,method %in% shown)[,c('sex','method','mean_z','sd_z','slope_per10cm','p10_overall_pct','p10_Q1minusQ5_pp','p10_quintile_range_pp')],row.names=FALSE)