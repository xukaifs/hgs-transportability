invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
root<-getwd();base<-dirname(root);source(file.path(base,'KNHANES_HGS_phase2_20260910','core.R'))
mp<-file.path(base,'KNHANES_pooled_BCT_20260910','pooled_BCT_models.rds');hash0<-unname(tools::md5sum(mp));b<-readRDS(mp)
kr<-readRDS(file.path(base,'KNHANES_pooled_BCT_20260910','pooled_BCT_scores.rds'))
old<-file.path(base,'KNHANES_HGS_pilot_20260909');kd<-rbind(readRDS(file.path(old,'development_data.rds')),readRDS(file.path(old,'validation_data.rds')));kd<-subset(kd,eligible)
us<-readRDS(file.path(base,'HGS_pilot_20260908','harmonized_raw_qc.rds'))
us$height_m<-us$height/100;us$height10<-us$height/10;us$weight<-us$weight/2
us$strata_pool<-interaction(us$cycle,us$SDMVSTRA,drop=TRUE);us$psu_pool<-interaction(us$cycle,us$SDMVSTRA,us$SDMVPSU,drop=TRUE)
us$design_ok<-is.finite(us$weight)&us$weight>0&!is.na(us$strata_pool)&!is.na(us$psu_pool)
# Preserve existing basic QC and NHANES valid-effort trial screening.
stopifnot(sum(us$eligible)==9334,all(subset(us,eligible)$BMXWT>0),all(subset(us,eligible)$age>=20&subset(us,eligible)$age<=79))
stats<-function(d,z)data.frame(N=nrow(d),mean_z=weighted.mean(z,d$weight),sd_z=wsd(z,d$weight),p10=100*weighted.mean(z< -1.282,d$weight),p5=100*weighted.mean(z< -1.645,d$weight))
rows<-list();quints<-list()
for(s in c('Female','Male')){
 d<-subset(us,eligible&sex==s);ref<-subset(kd,sex==s);stopifnot(nrow(ref)==b[[s]]$n)
 p<-score_model(b[[s]],d);z<-p$z;stopifnot(all(is.finite(z)),!any(p$clipped))
 edges<-wq(d$height,d$weight,c(.2,.4,.6,.8));q<-cut(d$height,c(-Inf,edges,Inf),labels=paste0('Q',1:5),include.lowest=TRUE)
 qq<-do.call(rbind,lapply(levels(q),function(k){i<-q==k;cbind(sex=s,quintile=k,stats(d[i,],z[i]))}));quints[[s]]<-qq
 f<-us;f$z<-NA_real_;f$z[match(d$SEQN,f$SEQN)]<-z;f$p10<-as.numeric(f$z< -1.282);f$p5<-as.numeric(f$z< -1.645)
 de<-svydesign(ids=~psu_pool,strata=~strata_pool,weights=~weight,data=subset(f,design_ok),nest=TRUE);de<-subset(de,!is.na(z))
 fit<-svyglm(z~height10,de);ci<-confint(fit)[2,];r<-cbind(sex=s,stats(d,z),slope_10cm=unname(coef(fit)[2]),slope_ci_low=ci[1],slope_ci_high=ci[2])
 for(t in c('p10','p5')){fit<-svyglm(reformulate('height10',response=t),de,family=quasibinomial());ci<-confint(fit)[2,];r[[paste0('OR_',t)]]<-exp(coef(fit)[2]);r[[paste0('OR_',t,'_ci_low')]]<-exp(ci[1]);r[[paste0('OR_',t,'_ci_high')]]<-exp(ci[2]);r[[paste0(t,'_Q1_pct')]]<-qq[[t]][1];r[[paste0(t,'_Q5_pct')]]<-qq[[t]][5];r[[paste0(t,'_drift_pp')]]<-qq[[t]][1]-qq[[t]][5];r[[paste0(t,'_range_pp')]]<-diff(range(qq[[t]]))}
 pct<-wq(ref$height,ref$weight,c(.01,.99));r$korea_age_min<-min(ref$age);r$korea_age_max<-max(ref$age);r$korea_height_min<-min(ref$height);r$korea_height_max<-max(ref$height);r$korea_height_p01<-pct[1];r$korea_height_p99<-pct[2]
 flags<-list(age_outside=d$age<min(ref$age)|d$age>max(ref$age),height_below_min=d$height<min(ref$height),height_above_max=d$height>max(ref$height),height_below_p01=d$height<pct[1],height_above_p99=d$height>pct[2])
 for(nm in names(flags)){r[[paste0(nm,'_N')]]<-sum(flags[[nm]]);r[[paste0(nm,'_weighted_pct')]]<-100*weighted.mean(flags[[nm]],d$weight)}
 for(i in 1:4)r[[paste0('US_height_cut_',i)]]<-edges[i]
 r$model_md5<-hash0;r$CDF_clipped_N<-sum(p$clipped);rows[[s]]<-r
 stopifnot(sum(qq$N)==nrow(d),!anyDuplicated(d$SEQN))
}
ss<-do.call(rbind,rows);qq<-do.call(rbind,quints);stopifnot(sum(ss$N)==9334,nrow(qq)==10,unname(tools::md5sum(mp))==hash0)
write.csv(ss,'US_transport_overall.csv',row.names=FALSE,na='',fileEncoding='UTF-8');write.csv(qq,'US_transport_height_quintiles.csv',row.names=FALSE,na='',fileEncoding='UTF-8')
long<-rbind(data.frame(qq[,c('sex','quintile')],metric='Mean Z',value=qq$mean_z,target=0),data.frame(qq[,c('sex','quintile')],metric='P10 (%)',value=qq$p10,target=10),data.frame(qq[,c('sex','quintile')],metric='P5 (%)',value=qq$p5,target=5));long$metric<-factor(long$metric,levels=c('Mean Z','P10 (%)','P5 (%)'))
p<-ggplot(long,aes(quintile,value,color=sex,group=sex))+geom_hline(data=unique(long[,c('metric','target')]),aes(yintercept=target),linetype=2,color='grey50')+geom_line()+geom_point()+facet_wrap(~metric,ncol=1,scales='free_y')+theme_bw(base_size=12)+theme(legend.position='bottom')+labs(title='Frozen pooled Korean BCT in NHANES 2011-2014',subtitle='Korea 2014-2019 parameters unchanged; no recalibration',x='US sex-specific weighted height quintile',y=NULL,caption='P10: Z < -1.282; P5: Z < -1.645. Dashed lines show reference targets.')
ggsave('US_transport_by_height.png',p,width=9,height=10,dpi=160)
print(ss[,c('sex','N','mean_z','sd_z','slope_10cm','slope_ci_low','slope_ci_high','p10_Q1_pct','p10_Q5_pct','p5_Q1_pct','p5_Q5_pct')],row.names=FALSE)
cat('PASS: fixed model hash; same US cohort; finite CDF; support counted; quintile totals.\n')