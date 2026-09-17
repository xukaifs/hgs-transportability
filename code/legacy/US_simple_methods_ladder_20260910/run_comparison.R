invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
root<-getwd();base<-dirname(root);source(file.path(base,'KNHANES_HGS_phase2_20260910','core.R'));wc<-function(x,n)write.csv(x,paste0(n,'.csv'),row.names=FALSE,na='',fileEncoding='UTF-8')
prior<-file.path(base,'US_BCT_recalibration_ladder_20260910');old<-file.path(base,'KNHANES_HGS_pilot_20260909')
kr<-rbind(readRDS(file.path(old,'development_data.rds')),readRDS(file.path(old,'validation_data.rds')));kr$weight<-kr$weight/2;kr<-subset(kr,eligible)
models<-list();params<-list();us<-readRDS(file.path(base,'HGS_pilot_20260908','harmonized_raw_qc.rds'));us$height_m<-us$height/100;us$height10<-us$height/10
us$strata_pool<-factor(us$SDMVSTRA);us$psu_pool<-interaction(us$SDMVSTRA,us$SDMVPSU,drop=TRUE);us$design_ok<-is.finite(us$weight)&us$weight>0&!is.na(us$strata_pool)&!is.na(us$psu_pool)
for(s in c('Female','Male'))for(m in c('Height2','Allometry_ageadj')){
 k<-paste(s,m);b<-fit_traditional(subset(kr,sex==s),m);models[[k]]<-b
 a<-subset(us,eligible&cycle=='G'&sex==s);z<-score_model(b,a)$z;w<-a$weight;center<-weighted.mean(z,w);scale<-wsd(z,w);h0<-weighted.mean(a$height10,w);hc<-a$height10-h0;slope<-sum(w*hc*(z-center))/sum(w*hc^2);rs<-wsd(z-center-slope*hc,w)
 params[[k]]<-data.frame(sex=s,normalization=m,N_Korea=b$n,allometry_b=b$b,N_adaptation=nrow(a),center=center,scale=scale,height_center_10cm=h0,slope_z_per10cm=slope,residual_scale=rs)
}
pa<-do.call(rbind,params);wc(pa,'simple_adaptation_parameters');saveRDS(models,'frozen_pooled_Korean_simple_models.rds');saveRDS(pa,'frozen_simple_adaptation.rds')
cutpoints<-read.csv(file.path(prior,'heldout_height_cutpoints.csv'));rows<-list();quints<-list();scores<-list()
for(s in c('Female','Male'))for(m in c('Height2','Allometry_ageadj')){
 d<-subset(us,eligible&cycle=='H'&sex==s);a<-subset(pa,sex==s&normalization==m);z0<-score_model(models[[paste(s,m)]],d)$z;edges<-subset(cutpoints,sex==s)$height
 zz<-list(Frozen=z0,Location=z0-a$center,Location_scale=(z0-a$center)/a$scale,Location_scale_height=(z0-a$center-a$slope_z_per10cm*(d$height10-a$height_center_10cm))/a$residual_scale)
 for(method in names(zz)){
 z<-zz[[method]];mm<-metrics(d,z,edges);f<-subset(us,cycle=='H');f$z<-NA_real_;f$z[match(d$SEQN,f$SEQN)]<-z
 de<-svydesign(ids=~psu_pool,strata=~strata_pool,weights=~weight,data=subset(f,design_ok),nest=TRUE);de<-subset(de,!is.na(z));fit<-svyglm(z~height10,de);ci<-confint(fit)[2,]
 r<-data.frame(sex=s,normalization=m,method=method,N=nrow(d),mean_z=mm$summary$mean_z,sd_z=mm$summary$sd_z,slope_10cm=unname(coef(fit)[2]),slope_ci_low=ci[1],slope_ci_high=ci[2],p10=100*mm$summary$p10,p5=100*mm$summary$p5,p10_drift_pp=mm$summary$p10_drift_pp,p5_drift_pp=mm$summary$p5_drift_pp,p10_range_pp=mm$summary$p10_range_pp,p5_range_pp=mm$summary$p5_range_pp)
 qt<-mm$quintiles;rows[[paste(s,m,method)]]<-r;quints[[paste(s,m,method)]]<-data.frame(sex=s,normalization=m,method=method,quintile=qt$quintile,N=qt$n,mean_z=qt$mean_z,sd_z=qt$sd_z,p10=100*qt$p10,p5=100*qt$p5);scores[[paste(s,m,method)]]<-data.frame(SEQN=d$SEQN,sex=s,normalization=m,method=method,z=z)
 }
}
bct<-read.csv(file.path(prior,'US_ladder_heldout_overall.csv'));bct$normalization<-'BCT';ss<-do.call(rbind,rows);ss<-rbind(ss,bct[,names(ss)])
bq<-read.csv(file.path(prior,'US_ladder_heldout_height_quintiles.csv'));bq$normalization<-'BCT';qq<-do.call(rbind,quints);qq<-rbind(qq,bq[,names(qq)])
wc(ss,'heldout_normalization_ladder_comparison');wc(qq,'heldout_normalization_height_quintiles');saveRDS(do.call(rbind,scores),'simple_heldout_scores.rds')
stopifnot(nrow(ss)==24,nrow(qq)==120,all(ss$N[ss$sex=='Female']==2527),all(ss$N[ss$sex=='Male']==2371),!any(subset(us,cycle=='G')$SEQN %in% subset(us,cycle=='H')$SEQN))
for(s in c('Female','Male'))for(m in c('Height2','Allometry_ageadj')){r<-subset(ss,sex==s&normalization==m);a<-subset(pa,sex==s&normalization==m);f<-subset(r,method=='Frozen');l<-subset(r,method=='Location');ls<-subset(r,method=='Location_scale');h<-subset(r,method=='Location_scale_height');stopifnot(abs(f$slope_10cm-l$slope_10cm)<1e-10,abs(f$sd_z-l$sd_z)<1e-10,abs(ls$slope_10cm-f$slope_10cm/a$scale)<1e-10,abs(h$slope_10cm-(f$slope_10cm-a$slope_z_per10cm)/a$residual_scale)<1e-10)}
# BCT equality is checked numerically below.
# Numeric comparison confirms BCT values were copied without changes.
for(i in seq_len(nrow(bct))){r<-subset(ss,normalization=='BCT'&sex==bct$sex[i]&method==bct$method[i]);nn<-names(bct)[vapply(bct,is.numeric,logical(1))];stopifnot(max(abs(as.numeric(r[1,nn])-as.numeric(bct[i,nn])))<1e-12)}
wc(data.frame(check=c('Only two Korean simple references added','Allometry exponent fit on Korea only','US adaptation and evaluation IDs disjoint','Same 4898 held-out participants and original BCT cutpoints','Original BCT results numerically unchanged','Affine recalibration identities verified'),pass=TRUE),'integrity_checks')
print(subset(ss,method=='Location_scale')[,c('sex','normalization','mean_z','sd_z','slope_10cm','p10','p5','p10_range_pp','p5_range_pp')],row.names=FALSE);cat('COMPARISON COMPLETE\n')