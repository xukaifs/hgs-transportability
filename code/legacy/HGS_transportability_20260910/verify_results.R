invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
root<-getwd();phase<-file.path(dirname(root),'KNHANES_HGS_phase2_20260910');source(file.path(phase,'core.R'))
a<-read.csv('korea_annual_hgs.csv');q<-read.csv('korea_annual_qc.csv');u<-read.csv('nhanes_frozen_summary.csv');g<-read.csv('korea_standardization_groups.csv');z<-readRDS('nhanes_frozen_scores.rds');support<-read.csv('nhanes_training_support.csv')
stopifnot(nrow(a)==12,sum(a$n)==32311,nrow(u)==8,nrow(z)==4*9334,all(is.finite(z$z)),!any(z$clipped),!anyDuplicated(z[,c('uid','method')]))
for(s in c('Female','Male'))for(y in 2014:2019)for(k in c('ageband','heightband')){
 gg<-subset(g,sex==s&year==y&standardization==k);rr<-subset(a,sex==s&year==y)
 stopifnot(all(gg$n>0),abs(sum(gg$reference_proportion)-1)<1e-10,abs(sum(gg$reference_proportion*gg$hgs)-rr[[paste0(k,'_standardized_hgs')]])<1e-9,abs(sum(gg$year_proportion*gg$hgs)-rr$hgs_mean)<1e-9)
}
locked<-read.csv(file.path(phase,'model_lock_manifest.csv'));now<-unname(tools::md5sum(file.path(phase,'locked_development_models.rds')));stopifnot(now==locked$md5,read.csv('nhanes_freeze_manifest.csv')$model_md5==now)
checks<-data.frame(check=c('Korea 12 sex-year groups and 32311 eligible records','All direct-standardization cells nonempty','Standardized estimates independently match weighted group sums','Annual raw means match group-weighted means','US 9334 unique participants scored by all four methods','US scores finite and no CDF clipping','Frozen Korean model MD5 matches original Phase 2 lock'),pass=TRUE)
write.csv(checks,'integrity_checks.csv',row.names=FALSE)

cat('NUMERICAL CHECKS COMPLETE\n')
