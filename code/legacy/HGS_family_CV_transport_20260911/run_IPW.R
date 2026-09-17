invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
root<-getwd();base<-dirname(root);source(file.path(base,'KNHANES_HGS_phase2_20260910','core.R'));wc<-function(x,n)write.csv(x,paste0(n,'.csv'),row.names=FALSE,na='',fileEncoding='UTF-8')
old<-file.path(base,'KNHANES_HGS_pilot_20260909');all<-rbind(readRDS(file.path(old,'development_data.rds')),readRDS(file.path(old,'validation_data.rds')));all$weight<-all$weight/2
frame<-subset(all,age>=20&age<=79&sex %in% c('Female','Male')&design_ok&is.finite(height)&height>0&is.finite(bmi)&bmi>0&is.finite(body_weight)&body_weight>0)
frame$measured<-as.numeric(is.finite(frame$hgs)&frame$hgs>0&frame$n_valid>0);stopifnot(sum(frame$measured)==32311,identical(frame$measured==1,frame$eligible))
frame$sex<-factor(frame$sex);frame$year_f<-factor(frame$year)
de<-svydesign(ids=~psu_pool,strata=~strata_pool,weights=~weight,data=frame,nest=TRUE)
fit<-svyglm(measured~ns(age,4)+sex+ns(height,3)+ns(bmi,3)+year_f,de,family=quasibinomial());frame$prob<-as.numeric(predict(fit,type='response'));stopifnot(isTRUE(fit$converged),all(is.finite(frame$prob)),all(frame$prob>0&frame$prob<1))
scpath<-file.path(base,'KNHANES_pooled_BCT_20260910','pooled_BCT_scores.rds');mh<-unname(tools::md5sum(scpath));scores<-readRDS(scpath);frame$z<-scores$z[match(frame$uid,scores$uid)];stopifnot(all(is.finite(frame$z[frame$measured==1])))
curves<-list();diag<-list();weightrows<-list();flow<-list()
for(s in c('Female','Male'))for(y in 2014:2019){
 a<-subset(frame,sex==s&year==y);o<-subset(a,measured==1);resp<-weighted.mean(a$measured,a$weight);factor<-resp/o$prob;cap<-wq(factor,o$weight,.99);capped<-pmin(factor,cap)
 wd<-data.frame(uid=o$uid,sex=s,year=y,prob=o$prob,stabilized_factor=factor,capped_factor=capped,original_weight=o$weight);weightrows[[paste(s,y)]]<-wd
 flags<-list(Observed=rep(1,nrow(o)),IPW=factor,IPW_cap99=capped)
 for(m in names(flags)){w<-o$weight*flags[[m]];curves[[paste(s,y,m)]]<-data.frame(sex=s,year=y,method=m,N_observed=nrow(o),mean_HGS=weighted.mean(o$hgs,w),mean_z=weighted.mean(o$z,w),mean_age=weighted.mean(o$age,w),mean_height=weighted.mean(o$height,w),Kish_effective_N=sum(w)^2/sum(w^2))}
 diag[[paste(s,y)]]<-data.frame(sex=s,year=y,N_base=nrow(a),N_observed=nrow(o),weighted_response_pct=100*resp,prediction_min=min(a$prob),prediction_p01=wq(a$prob,a$weight,.01),prediction_median=wq(a$prob,a$weight,.5),prediction_max=max(a$prob),IPW_factor_max=max(factor),IPW_factor_p99=cap,cap_N=sum(factor>cap),cap_weighted_pct=100*weighted.mean(factor>cap,o$weight))
 target<-subset(all,age>=20&age<=79&sex==s&year==y);flow[[paste(s,y)]]<-data.frame(sex=s,year=y,N_age20_79=nrow(target),N_base=nrow(a),excluded_design_or_other_covariates=nrow(target)-nrow(a),N_observed=nrow(o),N_missing_or_invalid_HGS=sum(a$measured==0))
}
cu<-do.call(rbind,curves);di<-do.call(rbind,diag);wc(cu,'IPW_yearly_curves');wc(di,'IPW_overlap_weight_diagnostics');wc(do.call(rbind,flow),'IPW_sample_flow');saveRDS(do.call(rbind,weightrows),'IPW_observed_weights.rds');saveRDS(fit,'HGS_response_model.rds')
cf<-summary(fit)$coef;wc(data.frame(term=rownames(cf),coefficient=cf[,1],SE=cf[,2]),'HGS_response_coefficients')
prior<-read.csv(file.path(base,'HGS_transportability_20260910','korea_annual_hgs.csv'));for(i in seq_len(nrow(prior))){a<-subset(cu,sex==prior$sex[i]&year==prior$year[i]&method=='Observed');stopifnot(abs(a$mean_HGS-prior$hgs_mean[i])<1e-10)}
stopifnot(unname(tools::md5sum(scpath))==mh,nrow(cu)==36)
changes<-do.call(rbind,lapply(split(cu,list(cu$sex,cu$method),drop=TRUE),function(a){v<-function(y,nm)a[a$year==y,nm];data.frame(sex=a$sex[1],method=a$method[1],HGS_2014_to2018=v(2018,'mean_HGS')-v(2014,'mean_HGS'),HGS_2018_to2019=v(2019,'mean_HGS')-v(2018,'mean_HGS'),Z_2014_to2018=v(2018,'mean_z')-v(2014,'mean_z'),Z_2018_to2019=v(2019,'mean_z')-v(2018,'mean_z'))}));wc(changes,'IPW_year_change_summary')
l<-rbind(data.frame(cu[,c('sex','year','method')],metric='HGS (kg)',value=cu$mean_HGS),data.frame(cu[,c('sex','year','method')],metric='Pooled BCT mean Z',value=cu$mean_z))
p<-ggplot(l,aes(year,value,color=method,group=method))+geom_line()+geom_point()+facet_wrap(~sex+metric,scales='free_y',ncol=2)+scale_x_continuous(breaks=2014:2019)+theme_bw(base_size=11)+theme(legend.position='bottom')+labs(title='Korean HGS participation sensitivity',subtitle='Same observed participants and frozen pooled BCT; reweighting only',x='Year',y=NULL,caption='Response model: age, sex, height, BMI and year. Stabilized IPW and weighted 99th-percentile cap within sex/year.\nDescriptive estimates; no uncertainty propagation for response estimation. Does not address unmeasured selection or instrument changes.')
ggsave('IPW_yearly_sensitivity.png',p,width=12,height=8,dpi=180,bg='white')
wc(data.frame(check=c('Observed IPW cohort identical to original 32311 complete HGS cases','Response model converged; fitted probabilities strictly between zero and one','Original annual HGS means reproduced within 1e-10','Pooled BCT score file unchanged','Untrimmed and capped weights both retained'),pass=TRUE),'IPW_integrity_checks')
writeLines(c('IPW sensitivity is conditional on measured age/sex/height/BMI/year and assumes response exchangeability given these variables. It cannot establish absence of selection bias, causality, secular decline or device stability.',
'Base population: 20-79 with valid survey design, sex, height, bodyweight and BMI. Participants missing these covariates are excluded and counted, not imputed.',
'Outcome is valid observed HGS under existing QC. Stabilized numerator is sex-year weighted response proportion; denominator is pooled weighted logistic prediction. Cap is weighted 99th percentile of stabilized factor among observed in each sex/year; not a cap on survey weights.',
'No HGS reference fit/recalibration, no race analysis, no reciprocal transport. BCT remains locked primary.',
'Plots are descriptive without confidence intervals; response estimation and IPW uncertainty are not propagated.'),'IPW_notes.txt')
print(changes,row.names=FALSE);print(di,row.names=FALSE);cat('IPW SENSITIVITY COMPLETE\n')