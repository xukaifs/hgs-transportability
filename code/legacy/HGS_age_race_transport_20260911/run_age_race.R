invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
root<-getwd();base<-dirname(root);source(file.path(base,'KNHANES_HGS_phase2_20260910','core.R'));suppressPackageStartupMessages(library(haven))
wc<-function(x,n)write.csv(x,paste0(n,'.csv'),row.names=FALSE,na='',fileEncoding='UTF-8')
phase<-file.path(base,'KNHANES_HGS_phase2_20260910');ladder<-file.path(base,'US_BCT_recalibration_ladder_20260910');simple<-file.path(base,'US_simple_methods_ladder_20260910')
inputs<-c(file.path(phase,'temporal_scores.rds'),file.path(ladder,'heldout_scores.rds'),file.path(simple,'simple_heldout_scores.rds'))
hashes<-unname(tools::md5sum(inputs));ks<-readRDS(inputs[1]);bs<-readRDS(inputs[2]);ss<-readRDS(inputs[3]);bs$normalization<-'BCT';uscores<-rbind(ss,bs[,names(ss)])
kr<-readRDS(file.path(base,'KNHANES_HGS_pilot_20260909','validation_data.rds'))
us<-readRDS(file.path(base,'HGS_pilot_20260908','harmonized_raw_qc.rds'));us<-subset(us,cycle=='H');us$height10<-us$height/10;us$strata_pool<-factor(us$SDMVSTRA);us$psu_pool<-interaction(us$SDMVSTRA,us$SDMVPSU,drop=TRUE);us$design_ok<-is.finite(us$weight)&us$weight>0&!is.na(us$strata_pool)&!is.na(us$psu_pool)
# Official NHANES RIDRETH3: 1/2 Hispanic, 3 NH White, 4 NH Black, 6 NH Asian, 7 Other including multiracial.
demo<-read_xpt('D:/DXA身体成分分析/数据/人口学/DEMO_H.xpt');stopifnot(!anyDuplicated(demo$SEQN));ix<-match(us$SEQN,demo$SEQN);stopifnot(!anyNA(ix),all(us$age==demo$RIDAGEYR[ix]),all(us$height>0|!us$eligible,na.rm=TRUE))
us$race_code<-as.numeric(demo$RIDRETH3[ix]);race_map<-c('1'='Hispanic','2'='Hispanic','3'='Non-Hispanic White','4'='Non-Hispanic Black','6'='Non-Hispanic Asian','7'='Other / multiracial');us$race<-unname(race_map[as.character(us$race_code)]);stopifnot(!anyNA(subset(us,eligible)$race))
wc(data.frame(code=names(race_map),group=unname(race_map)),'race_code_mapping')
getdes<-function(full,score,idcol){stopifnot(!anyDuplicated(score[[idcol]]));ii<-match(score[[idcol]],full[[idcol]]);stopifnot(!anyNA(ii));full$z<-NA_real_;full$z[ii]<-score$z;full$age10<-full$age/10;full$age_c<-(full$age-50)/10;full$height_c<-(full$height-165)/10;full$age_group<-cut(full$age,c(20,40,60,80),right=FALSE,labels=c('20-39','40-59','60-79'));full$p10<-as.numeric(full$z< -1.282);full$p5<-as.numeric(full$z< -1.645)
 de<-svydesign(ids=~psu_pool,strata=~strata_pool,weights=~weight,data=subset(full,design_ok),nest=TRUE);subset(de,!is.na(z))}
stat<-function(de){d<-de$variables;meanfit<-svymean(~z,de);ci<-confint(meanfit);data.frame(N=nrow(d),mean_z=weighted.mean(d$z,d$weight),mean_ci_low=ci[1],mean_ci_high=ci[2],sd_z=wsd(d$z,d$weight),p10=100*weighted.mean(d$p10,d$weight),p5=100*weighted.mean(d$p5,d$weight),p10_events=sum(d$p10),p5_events=sum(d$p5))}
slopes<-list();agegroups<-list();interactions<-list();checks<-list();designs<-list()
for(country in c('Korea_temporal','US_heldout'))for(s in c('Female','Male'))for(m in c('BCT','Height2','Allometry_ageadj'))for(cal in if(country=='Korea_temporal')'Frozen'else c('Frozen','Location_scale')){
 if(country=='Korea_temporal'){mid<-switch(m,BCT='primary',Height2='Height2',Allometry_ageadj='age_adjusted_allometry');sc<-subset(ks,sex==s&method==mid);de<-getdes(kr,sc,'uid')}else{sc<-subset(uscores,sex==s&normalization==m&method==cal);de<-getdes(us,sc,'SEQN')}
 key<-paste(country,s,m,cal);designs[[key]]<-de;meta<-data.frame(scenario=country,sex=s,normalization=m,calibration=cal)
 fit<-svyglm(z~age10,de);ci<-confint(fit)[2,];slopes[[key]]<-cbind(meta,stat(de),age_slope_per10y=unname(coef(fit)[2]),age_slope_ci_low=ci[1],age_slope_ci_high=ci[2])
 agegroups[[key]]<-do.call(rbind,lapply(c('20-39','40-59','60-79'),function(gr)cbind(meta,age_group=gr,stat(subset(de,age_group==gr)))))
 ff<-svyglm(z~age_c*height_c,de);cf<-coef(ff);cc<-confint(ff);termnames<-c('age_c','height_c','age_c:height_c')
 interactions[[key]]<-do.call(rbind,lapply(termnames,function(nm)cbind(meta,term=nm,estimate=unname(cf[nm]),ci_low=cc[nm,1],ci_high=cc[nm,2],p_value=summary(ff)$coef[nm,4])))
 # Independent join-and-score integrity: new summaries must match existing source tables.
 if(country=='Korea_temporal'){tb<-read.csv(file.path(phase,'summary_validation.csv'));ref<-subset(tb,sex==s&method==mid)}else{tb<-read.csv(file.path(simple,'heldout_normalization_ladder_comparison.csv'));ref<-subset(tb,sex==s&normalization==m&method==cal)}
 stopifnot(nrow(ref)==1,abs(weighted.mean(de$variables$z,de$variables$weight)-ref$mean_z)<1e-10,abs(wsd(de$variables$z,de$variables$weight)-ref$sd_z)<1e-10)
 checks[[key]]<-cbind(meta,N=nrow(de$variables),score_join_unchanged=TRUE)
}
a<-do.call(rbind,slopes);ag<-do.call(rbind,agegroups);it<-do.call(rbind,interactions);wc(a,'age_transport_overall');wc(ag,'age_transport_groups');wc(it,'age_height_interaction_diagnostics')
cat('AGE DIAGNOSTICS COMPLETE\n');print(subset(a,normalization=='BCT')[,c('scenario','sex','calibration','age_slope_per10y','age_slope_ci_low','age_slope_ci_high')],row.names=FALSE);flush.console()
races<-list()
for(s in c('Female','Male'))for(m in c('BCT','Height2','Allometry_ageadj'))for(cal in c('Frozen','Location_scale')){
 de<-designs[[paste('US_heldout',s,m,cal)]]
 for(gr in c('Non-Hispanic Asian','Non-Hispanic White','Non-Hispanic Black','Hispanic','Other / multiracial')){
 dd<-subset(de,race==gr);fit<-svyglm(z~height10,dd);ci<-confint(fit)[2,]
 races[[paste(s,m,cal,gr)]]<-cbind(data.frame(sex=s,normalization=m,calibration=cal,race=gr),stat(dd),height_slope_per10cm=unname(coef(fit)[2]),height_slope_ci_low=ci[1],height_slope_ci_high=ci[2],PSU_N=length(unique(dd$variables$psu_pool)),mean_age=weighted.mean(dd$variables$age,dd$variables$weight),mean_height=weighted.mean(dd$variables$height,dd$variables$weight))
 }
}
rr<-do.call(rbind,races);wc(rr,'US_race_ethnicity_transport');wc(do.call(rbind,checks),'score_integrity_checks')
for(s in c('Female','Male'))for(m in c('BCT','Height2','Allometry_ageadj'))for(cal in c('Frozen','Location_scale'))stopifnot(sum(subset(rr,sex==s&normalization==m&calibration==cal)$N)==if(s=='Female')2527 else 2371)
stopifnot(nrow(a)==18,nrow(ag)==54,nrow(it)==54,nrow(rr)==60,identical(unname(tools::md5sum(inputs)),hashes));wc(data.frame(file=inputs,md5=hashes),'input_provenance')
cat('RACE DIAGNOSTICS COMPLETE; SOURCE SCORES UNCHANGED\n');print(subset(rr,normalization=='BCT'&calibration=='Frozen')[,c('sex','race','N','mean_z','mean_ci_low','mean_ci_high','height_slope_per10cm','p10','p5','p10_events','p5_events')],row.names=FALSE)
writeLines(c('Scope: sections III and IV only; no HGS model or recalibration parameters refitted.',
'Korea temporal: 2014-2016 trained references, evaluated 2017-2019. US held-out: pooled 2014-2019 Korean references, evaluated 2013-2014; Location_scale adapted using US 2011-2012 only.',
'Age groups: 20-39,40-59,60-79; age slope marginal per 10 years. Interaction model: Z ~ (age-50)/10 * (height-165)/10. Age main effect at height 165cm, height main effect at age 50; interaction per 10 years x 10cm. These are linear diagnostics, not new normalization models.',
'All mean/age/height/interactions use survey design and domain analysis; CI conditional on existing scores and calibration parameters. No multiple-testing adjustment; exploratory diagnostics, not confirmatory proof.',
'US RIDRETH3 joined by SEQN from DEMO_H: 1+2 Hispanic,3 NH White,4 NH Black,6 NH Asian,7 Other/multiracial. All five groups retained.',
'Official coding source: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/DEMO_H.htm',
'Race/ethnicity categories are self-reported broad population categories, not genetic ancestry, and NH Asian is not equivalent to Korean. Subgroup differences do not identify a causal country or ancestry effect. Age/height composition is not standardized across race groups here.',
'p10 and p5 are percentages. Raw event counts supplied; sparse tails should not be interpreted as stable estimates. P-values are exploratory. No subgroup-specific recalibration.',
'Age and race grouped sample counts reproduce original samples; original mean and SD verified unchanged; source score file hashes unchanged.'),'analysis_notes.txt')
# Minimal descriptive figures emphasize BCT while tables retain all three methods.
b<-subset(ag,normalization=='BCT');b$setting<-paste(b$scenario,b$calibration,sep=' / ')
p<-ggplot(b,aes(age_group,mean_z,color=setting,group=setting))+geom_hline(yintercept=0,linetype=2)+geom_line()+geom_point()+geom_errorbar(aes(ymin=mean_ci_low,ymax=mean_ci_high),width=.1)+facet_wrap(~sex)+theme_bw(base_size=11)+theme(legend.position='bottom')+labs(title='Age transportability: existing BCT scores',x='Age group',y='Weighted mean Z (95% CI)')
ggsave('BCT_age_transport.png',p,width=12,height=6,dpi=180,bg='white')
b<-subset(rr,normalization=='BCT');p<-ggplot(b,aes(mean_z,race,color=calibration))+geom_vline(xintercept=0,linetype=2)+geom_point(position=position_dodge(width=.4))+geom_errorbar(aes(xmin=mean_ci_low,xmax=mean_ci_high),position=position_dodge(width=.4),width=.2,orientation='y')+facet_wrap(~sex)+theme_bw(base_size=11)+theme(legend.position='bottom')+labs(title='US 2013-2014: BCT race/ethnicity sensitivity',x='Weighted mean Z (95% CI)',y=NULL)
ggsave('BCT_race_ethnicity.png',p,width=12,height=6,dpi=180,bg='white')
cat('ALL REQUESTED AGE AND RACE OUTPUTS COMPLETE\n')