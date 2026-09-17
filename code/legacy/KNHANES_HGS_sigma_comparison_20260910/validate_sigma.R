options(stringsAsFactors=FALSE,survey.lonely.psu='adjust',warn=1)
invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
.libPaths(c('C:/Users/Public/CodexRLib42Copy',.libPaths()))
suppressPackageStartupMessages({library(survey);library(splines);library(ggplot2)})
script_path<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
root<-dirname(normalizePath(script_path,winslash='/'));old<-file.path(dirname(root),'KNHANES_HGS_pilot_20260909')
wc<-function(x,n) write.csv(x,file.path(root,paste0(n,'.csv')),row.names=FALSE,na='',fileEncoding='UTF-8')
wq<-function(x,w,p){o<-order(x);cw<-cumsum(w[o])/sum(w);vapply(p,function(pp)x[o[which(cw>=pp)[1]]],numeric(1))}
wsd<-function(x,w)sqrt(weighted.mean((x-weighted.mean(x,w))^2,w))
frozenpath<-file.path(root,'frozen_Gaussian_location_scale_models.rds');hash0<-unname(tools::md5sum(frozenpath))
bundles<-readRDS(frozenpath);v0<-readRDS(file.path(old,'validation_data.rds'))
prior<-file.path(dirname(root),'KNHANES_HGS_phase15_20260910')
old_scores<-readRDS(file.path(prior,'validation_scores_six_methods.rds'));old_summary<-read.csv(file.path(prior,'summary_validation.csv'))
edges<-read.csv(file.path(old,'height_quintile_boundaries.csv'));new_scores<-list();score_params<-list();rows<-list();quints<-list();curves<-list()
for(s in c('Female','Male')){
 v<-subset(v0,eligible & sex==s);b<-bundles[[s]]
 stopifnot(!any(v$uid %in% b$development_uids),all(v$year %in% 2017:2019))
 X<-model.matrix(b$terms,data=v);mu<-as.numeric(X%*%b$mu_coefficients);sig<-exp(as.numeric(X%*%b$log_sigma_coefficients));z<-(v$hgs-mu)/sig
 stopifnot(all(is.finite(z)),all(is.finite(sig)),all(sig>0))
 q<-cut(v$height,c(-Inf,edges$height_cm[edges$sex==s],Inf),labels=paste0('Q',1:5),include.lowest=TRUE)
 new_scores[[s]]<-data.frame(uid=v$uid,year=v$year,sex=s,method='Conditional_heightSD',height_cm=v$height,weight=v$weight,quintile=q,z=z,p5=as.numeric(z< -1.645),p10=as.numeric(z< -1.282),p90=as.numeric(z>1.282),p95=as.numeric(z>1.645))
 score_params[[s]]<-data.frame(uid=v$uid,sex=s,mu=mu,sigma=sig,z=z)
}
sc<-rbind(old_scores,do.call(rbind,new_scores));saveRDS(sc,file.path(root,'validation_scores_seven_methods.rds'));saveRDS(do.call(rbind,score_params),file.path(root,'validation_mu_sigma.rds'))
methods<-c('Absolute','BMI','Height2','Allometry','Conditional','Conditional_heightSD','Conditional_LS')
for(s in c('Female','Male')) for(m in methods){
 a<-sc[sc$sex==s & sc$method==m,];full<-v0;idx<-match(a$uid,full$uid)
 stopifnot(!anyNA(idx),!anyDuplicated(a$uid))
 for(n in c('z','p5','p10','p90','p95')){full[[n]]<-NA_real_;full[[n]][idx]<-a[[n]]};full$q<-NA_character_;full$q[idx]<-as.character(a$quintile)
 des<-svydesign(ids=~psu_pool,strata=~strata_pool,weights=~weight,data=subset(full,design_ok),nest=TRUE);des<-subset(des,!is.na(z))
 f<-svyglm(z~height10,des);co<-summary(f)$coef;ci<-confint(f)[2,]
 corr<-weighted.mean((a$z-weighted.mean(a$z,a$weight))*(a$height_cm-weighted.mean(a$height_cm,a$weight)),a$weight)/(wsd(a$z,a$weight)*wsd(a$height_cm,a$weight))
 r<-data.frame(sex=s,method=m,n_development=bundles[[s]]$n_development,n_validation=nrow(a),mean_z=weighted.mean(a$z,a$weight),sd_z=wsd(a$z,a$weight),median_z=wq(a$z,a$weight,.5),slope_per10cm=co[2,1],slope_ci_low=ci[1],slope_ci_high=ci[2],slope_p=co[2,4],weighted_correlation=corr)
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
delta<-do.call(rbind,lapply(c('Female','Male'),function(s){
 c0<-ss[ss$sex==s & ss$method=='Conditional',];c1<-ss[ss$sex==s & ss$method=='Conditional_heightSD',]
 do.call(rbind,lapply(c('mean_z','sd_z','slope_per10cm','p5_overall_pct','p10_overall_pct','p5_Q1minusQ5_pp','p10_Q1minusQ5_pp','sd_Q5_over_Q1','sd_quintile_range','p5_quintile_range_pp','p10_quintile_range_pp'),function(metric)data.frame(sex=s,metric=metric,constant_sigma=c0[[metric]],varying_sigma=c1[[metric]],change=c1[[metric]]-c0[[metric]])))
}))
names(delta)[names(delta)=='varying_sigma']<-'height_sigma'
wc(delta,'height_sigma_vs_baseline')
three<-ss[ss$method %in% c('Conditional','Conditional_heightSD','Conditional_LS'),]
wc(three,'three_sigma_models')
# All original outputs must remain numerically unchanged for their common fields.
common<-intersect(names(ss),names(old_summary));common_num<-common[vapply(ss[common],is.numeric,logical(1))]
errs<-sapply(seq_len(nrow(old_summary)),function(i){j<-which(ss$sex==old_summary$sex[i]&ss$method==old_summary$method[i]);max(abs(unlist(ss[j,common_num])-unlist(old_summary[i,common_num])))})
stopifnot(max(errs)<1e-9,identical(hash0,unname(tools::md5sum(frozenpath))))
wc(data.frame(check=c('Original six methods unchanged','Frozen height-sigma model unchanged during validation'),pass=c(max(errs)<1e-9,TRUE),max_error=c(max(errs),0)),'validation_checks')
# Focus figures on the prespecified competitors and old/new conditional models.
shown<-c('Conditional','Conditional_heightSD','Conditional_LS');pal<-c(Conditional='#0072B2',Conditional_heightSD='#009E73',Conditional_LS='#D55E00')
labels<-c(Conditional='Baseline: constant SD',Conditional_heightSD='Model 1: SD(height)',Conditional_LS='Model 2: SD(age + height)')
qq<-subset(qq,method %in% shown);gg<-subset(gg,method %in% shown)
qq$method<-factor(qq$method,levels=shown);gg$method<-factor(gg$method,levels=shown)
th<-theme_bw(base_size=12)+theme(legend.position='bottom',legend.title=element_blank(),panel.grid.minor=element_blank(),strip.background=element_rect(fill='#F0F2F5'),plot.title=element_text(face='bold'),plot.caption=element_text(hjust=0,size=9))
cs<-scale_color_manual(values=pal,labels=labels)
p<-ggplot(gg,aes(height_cm,z,color=method))+geom_hline(yintercept=0,linetype=2,color='grey55')+geom_line(linewidth=.85)+facet_wrap(~sex,scales='free_x')+cs+th+labs(title='Phase 1.5: mean dependence',subtitle='Frozen KNHANES 2017–2019 validation; development 2014–2016',x='Height (cm)',y='Z (survey-weighted spline)',caption='No validation-set location or scale recalibration. Curves are validation diagnostics.')
ggsave(file.path(root,'Z_vs_height.png'),p,width=11.5,height=5.5,dpi=200)
long<-rbind(transform(qq,tail='P10',pct=p10_pct,target=10),transform(qq,tail='P5',pct=p5_pct,target=5))
p<-ggplot(long,aes(quintile,pct,color=method,group=method))+geom_hline(data=unique(long[,c('tail','target')]),aes(yintercept=target),linetype=2,color='grey55')+geom_line(linewidth=.8)+geom_point(size=2)+facet_grid(tail~sex)+cs+th+labs(title='Phase 1.5: lower-tail classification',x='Frozen sex-specific height quintile',y='Low-strength prevalence (%)',caption='P10: Z < −1.282; P5: Z < −1.645. Dashed lines: nominal coverage. Q1 shortest; Q5 tallest.')
ggsave(file.path(root,'P10_P5_by_height_quintile.png'),p,width=11.5,height=7.4,dpi=200)
p<-ggplot(qq,aes(quintile,sd_z,color=method,group=method))+geom_hline(yintercept=1,linetype=2,color='grey55')+geom_line(linewidth=.8)+geom_point(size=2)+facet_wrap(~sex)+cs+th+labs(title='Phase 1.5: scale dependence',x='Frozen sex-specific height quintile',y='Weighted SD(Z)',caption='Descriptive within-quintile SD. Q1 shortest; Q5 tallest. All models use the same participants.')
ggsave(file.path(root,'SDZ_by_height_quintile.png'),p,width=11.5,height=5.5,dpi=200)
capture.output(sessionInfo(),file=file.path(root,'validation_sessionInfo.txt'))
print(three[,c('sex','method','mean_z','sd_z','slope_per10cm','p10_overall_pct','p5_overall_pct','p10_Q1minusQ5_pp','p5_Q1minusQ5_pp','sd_Q1','sd_Q5')],row.names=FALSE,digits=4)
