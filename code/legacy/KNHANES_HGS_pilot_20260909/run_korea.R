options(stringsAsFactors=FALSE, survey.lonely.psu='adjust', warn=1)
invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
.libPaths(c('C:/Users/Public/CodexRLib42Copy',.libPaths()))
suppressPackageStartupMessages({library(haven);library(survey);library(splines);library(ggplot2)})
root<-normalizePath('新一篇/KNHANES_HGS_pilot_20260909',winslash='/')
wc<-function(x,n) write.csv(x,file.path(root,paste0(n,'.csv')),row.names=FALSE,na='',fileEncoding='UTF-8')
wq<-function(x,w,p){stopifnot(all(is.finite(x)),all(is.finite(w)),all(w>0));o<-order(x);cw<-cumsum(w[o])/sum(w);vapply(p,function(pp) x[o[which(cw>=pp)[1]]],numeric(1))}
wsd<-function(x,w) sqrt(weighted.mean((x-weighted.mean(x,w))^2,w))
rowmax<-function(a) apply(a,1,function(v) if(all(is.na(v))) NA_real_ else max(v,na.rm=TRUE))
trials<-paste0('gs_mea_',rep(c('r','l'),each=3),'_',rep(1:3,2))
needed<-c('id','year','sex','age','he_ht','he_wt','he_bmi','he_wc',trials,'gs_use','wt_itvex','psu','kstrata')
units<-c(id='identifier',year='calendar year',sex='code',age='years',he_ht='cm',he_wt='kg',he_bmi='kg/m^2',he_wc='cm',gs_use='code',wt_itvex='population weight',psu='identifier',kstrata='identifier')
units<-c(units,setNames(rep('kg',6),trials))
dict<-list();flow<-list();qc<-list();raws<-list();missing_fields<-list();allmeta<-list()
for(yr in 2014:2019){
 path<-file.path(root,'raw',paste0('hn',yr-2000,'_all.sas7bdat'))
 x<-read_sas(path);orig<-names(x);names(x)<-tolower(orig)
 stopifnot(all(needed %in% names(x)),all(x$year==yr),!anyDuplicated(x$id))
 labs<-vapply(x,function(v){a<-attr(v,'label');if(is.null(a)) '' else a},character(1))
 allmeta[[as.character(yr)]]<-data.frame(year=yr,variable=orig,label=labs)
 for(n in needed){
  v<-x[[n]];obs<-v[!is.na(v)];missing<-if(n=='gs_use') 'NA; observed 8/9 retained as raw codes, not used in analysis' else 'SAS missing read as NA; no labelled numeric missing codes supplied'
  dict[[paste(yr,n)]]<-data.frame(year=yr,variable=orig[match(n,names(x))],harmonized_name=n,label=labs[n],unit=units[n],n=nrow(x),n_missing=sum(is.na(v)),tagged_missing=if(is.numeric(v)) sum(is_tagged_na(v)) else 0,observed_min=as.character(min(obs)),observed_max=as.character(max(obs)),observed_codes=if(n %in% c('sex','gs_use')) paste(sort(unique(obs)),collapse=';') else '',missing_code_note=missing,source=basename(path))
 }
 # Audit every name/label for test-level QC, not just a selected column list.
 matched<-names(x)[grepl('^gs|grip|effort',names(x),ignore.case=TRUE)|grepl('악력|손목|오른손|왼손|grip|effort',labs,ignore.case=TRUE)]
 missing_fields[[as.character(yr)]]<-data.frame(year=yr,matching_variables=paste(matched,collapse=';'),effort_flag_available=FALSE,completion_flag_available=FALSE,cannot_test_reason_available=FALSE)
 d<-as.data.frame(x[,needed]);raws[[as.character(yr)]]<-d
 a<-as.matrix(d[,trials]);bad<-!is.na(a)&(!is.finite(a)|a<0|a>100);a[bad]<-NA_real_
 d$hgs<-rowmax(a);d$n_valid<-rowSums(!is.na(a));d$max_right<-rowmax(a[,1:3,drop=FALSE]);d$max_left<-rowmax(a[,4:6,drop=FALSE]);d$bilateral_sum<-d$max_right+d$max_left;d$bilateral_mean<-d$bilateral_sum/2
 d$height<-as.numeric(d$he_ht);d$height10<-d$height/10;d$height_m<-d$height/100;d$bmi<-as.numeric(d$he_bmi);d$body_weight<-as.numeric(d$he_wt)
 d$sex_code<-d$sex;d$sex<-ifelse(d$sex_code==1,'Male',ifelse(d$sex_code==2,'Female',NA))
 d$weight<-as.numeric(d$wt_itvex)/3;d$psu_pool<-paste(d$year,d$psu,sep='_');d$strata_pool<-paste(d$year,d$kstrata,sep='_')
 d$design_ok<-is.finite(d$weight)&d$weight>0&!is.na(d$psu)&nzchar(d$psu)&!is.na(d$kstrata)
 ageok<-!is.na(d$age)&d$age>=20&d$age<=79
 hgsok<-ageok&is.finite(d$hgs)&d$hgs>0
 anthropok<-hgsok&is.finite(d$height)&d$height>0&is.finite(d$body_weight)&d$body_weight>0&is.finite(d$bmi)&d$bmi>0
 d$eligible<-anthropok&!is.na(d$sex)&d$design_ok
 d$split<-if(yr<=2016) 'development' else 'validation';d$uid<-paste(d$year,d$id,sep='_')
 d$age20_69<-d$age>=20&d$age<=69;d$age50_79<-d$age>=50&d$age<=79
 flow[[as.character(yr)]]<-data.frame(year=yr,raw_n=nrow(d),age20_79_n=sum(ageok),hgs_valid_n=sum(hgsok),height_weight_bmi_complete_n=sum(anthropok),final_n=sum(d$eligible),final_female_n=sum(d$eligible&d$sex=='Female',na.rm=TRUE),final_male_n=sum(d$eligible&d$sex=='Male',na.rm=TRUE))
 qc[[as.character(yr)]]<-data.frame(year=yr,invalid_trial_count=sum(bad),nonmissing_trials=sum(!is.na(as.matrix(d[,trials]))),positive_below5_trials=sum(a>0&a<5,na.rm=TRUE),final_one_hand_n=sum(d$eligible&xor(is.na(d$max_right),is.na(d$max_left))),final_six_valid_n=sum(d$eligible&d$n_valid==6),final_hgs_min=min(d$hgs[d$eligible]),final_hgs_max=max(d$hgs[d$eligible]),waist_missing_final=sum(d$eligible&is.na(d$he_wc)))
 raws[[as.character(yr)]]<-d
 cat('Prepared',yr,'eligible',sum(d$eligible),'\n')
}
dat<-do.call(rbind,raws);rownames(dat)<-NULL;stopifnot(!anyDuplicated(dat$uid))
dat$psu_pool<-factor(dat$psu_pool);dat$strata_pool<-factor(dat$strata_pool)
wc(do.call(rbind,dict),'variable_dictionary');wc(do.call(rbind,flow),'sample_flow');wc(do.call(rbind,qc),'qc_audit');wc(do.call(rbind,missing_fields),'qc_field_availability')
saveRDS(allmeta,file.path(root,'all_variable_metadata.rds'))
saveRDS(dat,file.path(root,'harmonized_data.rds'))
saveRDS(subset(dat,split=='development'),file.path(root,'development_data.rds'))
saveRDS(subset(dat,split=='validation'),file.path(root,'validation_data.rds'))

# Fit phase receives development records only; do not load validation here.
rm(dat,raws,x,d,a)
dev_all<-readRDS(file.path(root,'development_data.rds'))
dev_design<-svydesign(ids=~psu_pool,strata=~strata_pool,weights=~weight,data=subset(dev_all,design_ok),nest=TRUE)
methods<-c('Absolute','BMI','Height2','Allometry','Conditional')
y_for<-function(a,m,b) switch(m,Absolute=a$hgs,BMI=a$hgs/a$bmi,Height2=a$hgs/a$height_m^2,Allometry=a$hgs/a$height_m^b,Conditional=a$hgs)
models<-list();exponents<-list()
for(s in c('Female','Male')){
 dd<-subset(dev_design,eligible & sex==s);d<-dd$variables
 af<-svyglm(log(hgs)~log(height_m),design=dd);b<-unname(coef(af)[2]);ac<-confint(af)[2,]
 exponents[[s]]<-data.frame(sex=s,n=nrow(d),b=b,se=sqrt(vcov(af)[2,2]),ci_low=ac[1],ci_high=ac[2],design_df=af$df.residual)
 for(m in methods){
  d$y<-y_for(d,m,b)
  ff<-if(m=='Conditional') y~ns(age,df=4)+ns(height,df=3) else y~ns(age,df=4)
  fit<-lm(ff,data=d,weights=weight/mean(weight))
  sig<-sqrt(weighted.mean(residuals(fit)^2,d$weight))
  models[[paste(s,m,sep='_')]]<-list(model=fit,sigma=sig,b=b,sex=s,method=m,development_years=2014:2016,development_n=nrow(d),development_age_range=range(d$age),development_height_range=range(d$height),development_uids=d$uid)
 }
}
wc(do.call(rbind,exponents),'allometry_exponents')
saveRDS(models,file.path(root,'frozen_Korea_2014_2016_models.rds'));rm(models,d,dd,af,dev_design,dev_all)
model_md5_before<-unname(tools::md5sum(file.path(root,'frozen_Korea_2014_2016_models.rds')))

# Validation phase: reload frozen file and validation records only.
models<-readRDS(file.path(root,'frozen_Korea_2014_2016_models.rds'))
val_all<-readRDS(file.path(root,'validation_data.rds'))
validation_scores<-list();summary_rows<-list();quintile_rows<-list();curve_rows<-list();boundary_rows<-list();design_rows<-list()
for(s in c('Female','Male')){
 v<-subset(val_all,eligible & sex==s)
 edges<-wq(v$height,v$weight,c(.2,.4,.6,.8));stopifnot(all(diff(edges)>0))
 v$q<-cut(v$height,c(-Inf,edges,Inf),labels=paste0('Q',1:5),include.lowest=TRUE)
 boundary_rows[[s]]<-data.frame(sex=s,p=c(.2,.4,.6,.8),height_cm=edges)
 for(m in methods){
  bundle<-models[[paste(s,m,sep='_')]]
  stopifnot(!any(v$uid %in% bundle$development_uids),all(v$year %in% 2017:2019))
  v$z<-(y_for(v,m,bundle$b)-as.numeric(predict(bundle$model,newdata=v)))/bundle$sigma
  stopifnot(all(is.finite(v$z)))
  v$p5<-as.numeric(v$z< -1.645);v$p10<-as.numeric(v$z< -1.282);v$p90<-as.numeric(v$z>1.282);v$p95<-as.numeric(v$z>1.645)
  full<-val_all;ii<-match(v$uid,full$uid)
  for(n in c('z','p5','p10','p90','p95')){full[[n]]<-NA_real_;full[[n]][ii]<-v[[n]]}
  full$q<-NA_character_;full$q[ii]<-as.character(v$q)
  des_full<-svydesign(ids=~psu_pool,strata=~strata_pool,weights=~weight,data=subset(full,design_ok),nest=TRUE)
  des<-subset(des_full,!is.na(z))
  slope<-svyglm(z~height10,design=des);ct<-summary(slope)$coef;ci<-confint(slope)[2,]
  covzh<-weighted.mean((v$z-weighted.mean(v$z,v$weight))*(v$height-weighted.mean(v$height,v$weight)),v$weight)
  rr<-data.frame(sex=s,method=m,n_development=bundle$development_n,n_validation=nrow(v),mean_z=weighted.mean(v$z,v$weight),sd_z=wsd(v$z,v$weight),median_z=wq(v$z,v$weight,.5),slope_per10cm=ct[2,1],slope_ci_low=ci[1],slope_ci_high=ci[2],slope_p=ct[2,4],weighted_correlation=covzh/(wsd(v$z,v$weight)*wsd(v$height,v$weight)),sigma_development=bundle$sigma)
  for(t in c('p5','p10','p90','p95')) rr[[paste0(t,'_overall_pct')]]<-100*weighted.mean(v[[t]],v$weight)
  for(t in c('p5','p10')){
   lf<-svyglm(as.formula(paste(t,'~height10')),design=des,family=quasibinomial());lc<-confint(lf)[2,]
   rr[[paste0(t,'_OR_per10cm')]]<-exp(coef(lf)[2]);rr[[paste0(t,'_OR_ci_low')]]<-exp(lc[1]);rr[[paste0(t,'_OR_ci_high')]]<-exp(lc[2]);rr[[paste0(t,'_trend_p')]]<-summary(lf)$coef[2,4]
  }
  qq<-do.call(rbind,lapply(levels(v$q),function(qqq){
   a<-v[as.character(v$q)==qqq,];quant<-wq(a$z,a$weight,c(.25,.75))
   row<-data.frame(sex=s,method=m,quintile=qqq,n=nrow(a),weighted_n=sum(a$weight),height_min_cm=min(a$height),height_max_cm=max(a$height),mean_z=weighted.mean(a$z,a$weight),sd_z=wsd(a$z,a$weight),iqr_z=diff(quant))
   for(t in c('p5','p10','p90','p95')){row[[paste0(t,'_pct')]]<-100*weighted.mean(a[[t]],a$weight);row[[paste0(t,'_events')]]<-sum(a[[t]])}
   row
  }))
  for(t in c('p5','p10','p90','p95')) rr[[paste0(t,'_Q1minusQ5_pp')]]<-qq[[paste0(t,'_pct')]][1]-qq[[paste0(t,'_pct')]][5]
  rr$sd_Q1<-qq$sd_z[1];rr$sd_Q5<-qq$sd_z[5];rr$sd_Q5_over_Q1<-qq$sd_z[5]/qq$sd_z[1]
  sf<-svyglm(z~ns(height10,df=3),design=des)
  grid<-data.frame(height10=seq(wq(v$height10,v$weight,.01),wq(v$height10,v$weight,.99),length.out=120))
  gp<-predict(sf,newdata=grid,se.fit=TRUE);gs<-sqrt(as.numeric(attr(gp,'var')));crit<-qt(.975,df=sf$df.residual)
  grid$z<-as.numeric(gp);grid$ci_low<-grid$z-crit*gs;grid$ci_high<-grid$z+crit*gs;grid$height_cm<-10*grid$height10;grid$sex<-s;grid$method<-m
  key<-paste(s,m,sep='_');summary_rows[[key]]<-rr;quintile_rows[[key]]<-qq;curve_rows[[key]]<-grid
  validation_scores[[key]]<-data.frame(uid=v$uid,year=v$year,sex=s,method=m,height_cm=v$height,weight=v$weight,quintile=v$q,z=v$z,p5=v$p5,p10=v$p10,p90=v$p90,p95=v$p95)
  cat('Validated',key,'slope',rr$slope_per10cm,'P10 drift pp',rr$p10_Q1minusQ5_pp,'\n')
 }
}
ss<-do.call(rbind,summary_rows);qq<-do.call(rbind,quintile_rows);gg<-do.call(rbind,curve_rows);sc<-do.call(rbind,validation_scores)
ss$residual_slope_reduction<-vapply(seq_len(nrow(ss)),function(i) 1-abs(ss$slope_per10cm[i])/abs(ss$slope_per10cm[ss$sex==ss$sex[i]&ss$method=='Absolute']),numeric(1))
wc(ss,'summary_validation');wc(qq,'height_quintiles');wc(gg,'spline_curve_data');wc(do.call(rbind,boundary_rows),'height_quintile_boundaries')
saveRDS(sc,file.path(root,'frozen_validation_scores.rds'))
stopifnot(identical(model_md5_before,unname(tools::md5sum(file.path(root,'frozen_Korea_2014_2016_models.rds')))))
colors<-c(Absolute='#D55E00',BMI='#009E73',Height2='#CC79A7',Allometry='#B39A00',Conditional='#0072B2')
gg$method<-factor(gg$method,levels=methods);qq$method<-factor(qq$method,levels=methods)
th<-theme_bw(base_size=12)+theme(legend.position='bottom',legend.title=element_blank(),panel.grid.minor=element_blank(),strip.background=element_rect(fill='#F0F2F5'),plot.title=element_text(face='bold'),plot.caption=element_text(hjust=0,size=9))
p1<-ggplot(gg,aes(height_cm,z,color=method))+geom_hline(yintercept=0,color='grey55',linetype=2)+geom_line(linewidth=.85)+facet_wrap(~sex,scales='free_x')+scale_color_manual(values=colors)+th+labs(title='Frozen temporal validation: mean dependence',subtitle='KNHANES 2014–2016 development / 2017–2019 validation; ages 20–79',x='Height (cm)',y='Z (survey-weighted spline)',caption='Curves are validation diagnostics, not refitted scoring models. Display: weighted height P1–P99.')
ggsave(file.path(root,'Z_vs_height.png'),p1,width=11,height=5.5,dpi=200)
long<-rbind(transform(qq,tail='P10',prevalence=p10_pct,target=10),transform(qq,tail='P5',prevalence=p5_pct,target=5))
p2<-ggplot(long,aes(quintile,prevalence,color=method,group=method))+geom_hline(data=unique(long[,c('tail','target')]),aes(yintercept=target),linetype=2,color='grey55')+geom_line(linewidth=.8)+geom_point(size=2)+facet_grid(tail~sex)+scale_color_manual(values=colors)+th+labs(title='Frozen temporal validation: lower-tail classification',x='Sex-specific weighted height quintile',y='Low-strength prevalence (%)',caption='P10: Z < −1.282; P5: Z < −1.645. Q1 shortest; Q5 tallest. Dashed lines: nominal coverage.')
ggsave(file.path(root,'P10_P5_by_height_quintile.png'),p2,width=11,height=7.4,dpi=200)
p3<-ggplot(qq,aes(quintile,sd_z,color=method,group=method))+geom_hline(yintercept=1,linetype=2,color='grey55')+geom_line(linewidth=.8)+geom_point(size=2)+facet_wrap(~sex)+scale_color_manual(values=colors)+th+labs(title='Frozen temporal validation: scale dependence',x='Sex-specific weighted height quintile',y='Weighted SD(Z)',caption='Descriptive population SD within each quintile, not a standard error. Q1 shortest; Q5 tallest.')
ggsave(file.path(root,'SDZ_by_height_quintile.png'),p3,width=11,height=5.5,dpi=200)
wc(data.frame(file=basename(list.files(file.path(root,'raw'),pattern='sas7bdat$',full.names=TRUE)),md5=unname(tools::md5sum(list.files(file.path(root,'raw'),pattern='sas7bdat$',full.names=TRUE)))),'raw_file_manifest')
capture.output(sessionInfo(),file=file.path(root,'sessionInfo.txt'))
print(ss[,c('sex','method','n_development','n_validation','slope_per10cm','p10_Q1minusQ5_pp','p5_Q1minusQ5_pp','sd_Q1','sd_Q5')],row.names=FALSE,digits=4)
