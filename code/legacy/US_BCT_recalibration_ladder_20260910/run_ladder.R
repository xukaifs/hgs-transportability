invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
root<-getwd();base<-dirname(root);source(file.path(base,'KNHANES_HGS_phase2_20260910','core.R'))
wc<-function(x,n)write.csv(x,paste0(n,'.csv'),row.names=FALSE,na='',fileEncoding='UTF-8')
mp<-file.path(base,'KNHANES_pooled_BCT_20260910','pooled_BCT_models.rds');hash0<-unname(tools::md5sum(mp));models<-readRDS(mp)
us<-readRDS(file.path(base,'HGS_pilot_20260908','harmonized_raw_qc.rds'));us$height_m<-us$height/100;us$height10<-us$height/10
us$strata_pool<-factor(us$SDMVSTRA);us$psu_pool<-interaction(us$SDMVSTRA,us$SDMVPSU,drop=TRUE);us$design_ok<-is.finite(us$weight)&us$weight>0&!is.na(us$strata_pool)&!is.na(us$psu_pool)
# Fit calibration only in G (2011-2012). Do not inspect H outcomes until parameters saved.
params<-list()
for(s in c('Female','Male')){
 a<-subset(us,eligible&cycle=='G'&sex==s);z<-score_model(models[[s]],a)$z;w<-a$weight
 center<-weighted.mean(z,w);scale<-wsd(z,w);h0<-weighted.mean(a$height10,w);hc<-a$height10-h0
 slope<-sum(w*hc*(z-center))/sum(w*hc^2);res<-z-center-slope*hc;rs<-wsd(res,w)
 stopifnot(scale>0,rs>0,abs(weighted.mean(res,w))<1e-10)
 params[[s]]<-data.frame(sex=s,N_adaptation=nrow(a),center=center,scale=scale,height_center_10cm=h0,slope_z_per10cm=slope,residual_scale=rs,korean_model_md5=hash0)
}
pa<-do.call(rbind,params);wc(pa,'adaptation_parameters');saveRDS(pa,'frozen_adaptation_parameters.rds');ph<-unname(tools::md5sum('frozen_adaptation_parameters.rds'))
rows<-list();quints<-list();scores<-list();cuts<-list()
for(s in c('Female','Male')){
 d<-subset(us,eligible&cycle=='H'&sex==s);a<-subset(pa,sex==s);z0<-score_model(models[[s]],d)$z;edges<-wq(d$height,d$weight,c(.2,.4,.6,.8));q<-cut(d$height,c(-Inf,edges,Inf),labels=paste0('Q',1:5),include.lowest=TRUE);cuts[[s]]<-data.frame(sex=s,prob=c(.2,.4,.6,.8),height=edges)
 zz<-list(Frozen=z0,Location=z0-a$center,Location_scale=(z0-a$center)/a$scale,Location_scale_height=(z0-a$center-a$slope_z_per10cm*(d$height10-a$height_center_10cm))/a$residual_scale)
 for(method in names(zz)){
 z<-zz[[method]];stopifnot(all(is.finite(z)));f<-subset(us,cycle=='H');f$z<-NA_real_;f$z[match(d$SEQN,f$SEQN)]<-z
 de<-svydesign(ids=~psu_pool,strata=~strata_pool,weights=~weight,data=subset(f,design_ok),nest=TRUE);de<-subset(de,!is.na(z));fit<-svyglm(z~height10,de);ci<-confint(fit)[2,]
 m<-metrics(d,z,edges);r<-data.frame(sex=s,method=method,N=nrow(d),mean_z=m$summary$mean_z,sd_z=m$summary$sd_z,slope_10cm=unname(coef(fit)[2]),slope_ci_low=ci[1],slope_ci_high=ci[2],p10=100*m$summary$p10,p5=100*m$summary$p5,p10_drift_pp=m$summary$p10_drift_pp,p5_drift_pp=m$summary$p5_drift_pp,p10_range_pp=m$summary$p10_range_pp,p5_range_pp=m$summary$p5_range_pp)
 qt<-m$quintiles;qt<-data.frame(sex=s,method=method,quintile=qt$quintile,N=qt$n,mean_z=qt$mean_z,sd_z=qt$sd_z,p10=100*qt$p10,p5=100*qt$p5)
 rows[[paste(s,method)]]<-r;quints[[paste(s,method)]]<-qt;scores[[paste(s,method)]]<-data.frame(SEQN=d$SEQN,sex=s,method=method,z=z)
 }
}
ss<-do.call(rbind,rows);qq<-do.call(rbind,quints);wc(ss,'US_ladder_heldout_overall');wc(qq,'US_ladder_heldout_height_quintiles');wc(do.call(rbind,cuts),'heldout_height_cutpoints');saveRDS(do.call(rbind,scores),'heldout_scores.rds')
stopifnot(sum(pa$N_adaptation)==4436,sum(ss$N[ss$method=='Frozen'])==4898,nrow(ss)==8,nrow(qq)==40,unname(tools::md5sum(mp))==hash0,unname(tools::md5sum('frozen_adaptation_parameters.rds'))==ph,!any(subset(us,cycle=='G')$SEQN %in% subset(us,cycle=='H')$SEQN))
# Affine transformations must match their implied held-out slopes and SDs.
for(s in c('Female','Male')){r<-subset(ss,sex==s);a<-subset(pa,sex==s);f<-subset(r,method=='Frozen');l<-subset(r,method=='Location');ls<-subset(r,method=='Location_scale');h<-subset(r,method=='Location_scale_height');stopifnot(abs(f$slope_10cm-l$slope_10cm)<1e-10,abs(f$sd_z-l$sd_z)<1e-10,abs(ls$slope_10cm-f$slope_10cm/a$scale)<1e-10,abs(h$slope_10cm-(f$slope_10cm-a$slope_z_per10cm)/a$residual_scale)<1e-10)}
wc(data.frame(check=c('2011-2012 adaptation only; 4436 participants','2013-2014 evaluation only; 4898 participants','Disjoint participant IDs','Korean and calibration parameter hashes unchanged','All four methods use same evaluation participants and cutpoints','Affine slope and scale identities verified'),pass=TRUE),'integrity_checks')
l<-rbind(data.frame(qq[,c('sex','method','quintile')],metric='Mean Z',value=qq$mean_z,target=0),data.frame(qq[,c('sex','method','quintile')],metric='P10 (%)',value=qq$p10,target=10),data.frame(qq[,c('sex','method','quintile')],metric='P5 (%)',value=qq$p5,target=5));l$metric<-factor(l$metric,levels=c('Mean Z','P10 (%)','P5 (%)'))
p<-ggplot(l,aes(quintile,value,color=method,group=method))+geom_hline(data=unique(l[,c('metric','target')]),aes(yintercept=target),linetype=2)+geom_line()+geom_point()+facet_grid(metric~sex,scales='free_y')+theme_bw(base_size=11)+theme(legend.position='bottom')+labs(title='US recalibration ladder: held-out 2013-2014',subtitle='Adaptation in 2011-2012 only; pooled Korean BCT unchanged',x='Held-out US weighted height quintile',y=NULL)
ggsave('US_ladder_heldout.png',p,width=12,height=10,dpi=160)
print(pa,row.names=FALSE);print(ss,row.names=FALSE);cat('LADDER COMPLETE; ALL CHECKS PASSED\n')