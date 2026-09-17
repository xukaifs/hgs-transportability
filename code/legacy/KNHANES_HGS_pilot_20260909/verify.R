options(stringsAsFactors=FALSE)
invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
.libPaths(c('C:/Users/Public/CodexRLib42Copy',.libPaths()))
suppressPackageStartupMessages({library(survey);library(splines)})
root<-'新一篇/KNHANES_HGS_pilot_20260909'
d<-readRDS(file.path(root,'harmonized_data.rds'));mods<-readRDS(file.path(root,'frozen_Korea_2014_2016_models.rds'));scores<-readRDS(file.path(root,'frozen_validation_scores.rds'))
ss<-read.csv(file.path(root,'summary_validation.csv'));qq<-read.csv(file.path(root,'height_quintiles.csv'));flow<-read.csv(file.path(root,'sample_flow.csv'))
stopifnot(nrow(d)==47309,nrow(ss)==10,nrow(qq)==50,sum(flow$final_n)==sum(d$eligible),!anyDuplicated(d$uid))
checks<-list();add<-function(n,pass,detail=''){stopifnot(pass);checks[[length(checks)+1]]<<-data.frame(check=n,pass=pass,detail=detail)}
add('Raw N and six-year flow reconcile',TRUE,paste(nrow(d),'records'))
trials<-paste0('gs_mea_',rep(c('r','l'),each=3),'_',rep(1:3,2))
mx<-apply(as.matrix(d[,trials]),1,function(a) {a<-a[is.finite(a)&a>=0&a<=100];if(length(a))max(a) else NA_real_})
add('HGSmax independently recomputed',isTRUE(all.equal(mx,d$hgs,check.attributes=FALSE)))
add('Development-validation separation',all(d$year[d$split=='development'] %in% 2014:2016)&&all(d$year[d$split=='validation'] %in% 2017:2019))
for(s in c('Female','Male')) for(m in c('Absolute','BMI','Height2','Allometry','Conditional')){
 key<-paste(s,m,sep='_');bundle<-mods[[key]];v<-subset(d,eligible & split=='validation' & sex==s)
 scored<-scores[scores$sex==s & scores$method==m,];scored<-scored[match(v$uid,scored$uid),]
 y<-switch(m,Absolute=v$hgs,BMI=v$hgs/v$bmi,Height2=v$hgs/v$height_m^2,Allometry=v$hgs/v$height_m^bundle$b,Conditional=v$hgs)
 design_matrix<-model.matrix(delete.response(terms(bundle$model)),data=v)
 manual<-(y-as.numeric(design_matrix%*%coef(bundle$model)))/bundle$sigma
 add(paste(key,'frozen matrix scoring'),max(abs(manual-scored$z))<1e-10)
 ix<-c(1,round(nrow(v)/2),nrow(v));single<-sapply(ix,function(i) (y[i]-as.numeric(predict(bundle$model,newdata=v[i,,drop=FALSE])))/bundle$sigma)
 add(paste(key,'single versus batch scores'),max(abs(single-scored$z[ix]))<1e-10)
 r<-ss[ss$sex==s & ss$method==m,];h<-v$height10;z<-scored$z;w<-v$weight
 slope<-sum(w*(h-weighted.mean(h,w))*(z-weighted.mean(z,w)))/sum(w*(h-weighted.mean(h,w))^2)
 add(paste(key,'slope independently recomputed'),abs(slope-r$slope_per10cm)<1e-10)
 add(paste(key,'same sample and no leakage'),nrow(v)==r$n_validation&&!any(v$uid %in% bundle$development_uids))
 for(q in paste0('Q',1:5)){
  a<-scored[scored$quintile==q,];target<-qq[qq$sex==s & qq$method==m & qq$quintile==q,]
  sd1<-sqrt(sum(a$weight*(a$z-sum(a$weight*a$z)/sum(a$weight))^2)/sum(a$weight))
  rates<-sapply(list(a$z< -1.645,a$z< -1.282,a$z>1.282,a$z>1.645),function(k)100*sum(a$weight[k])/sum(a$weight))
  add(paste(key,q,'SD and four tails'),nrow(a)==target$n&&abs(sd1-target$sd_z)<1e-10&&max(abs(rates-unlist(target[,c('p5_pct','p10_pct','p90_pct','p95_pct')])))<1e-9)
 }
}
add('Finite primary summary and intervals',all(vapply(ss[sapply(ss,is.numeric)],function(v) all(is.finite(v)),logical(1))))
add('Ordered confidence limits',all(ss$slope_ci_low<=ss$slope_per10cm & ss$slope_per10cm<=ss$slope_ci_high)&all(ss$p5_OR_ci_low<=ss$p5_OR_per10cm&ss$p5_OR_per10cm<=ss$p5_OR_ci_high)&all(ss$p10_OR_ci_low<=ss$p10_OR_per10cm&ss$p10_OR_per10cm<=ss$p10_OR_ci_high))
design_check<-do.call(rbind,lapply(c('development','validation'),function(sp){a<-subset(d,split==sp & design_ok);counts<-tapply(as.character(a$psu_pool),droplevels(a$strata_pool),function(p)length(unique(p)));data.frame(split=sp,n_full_design=nrow(a),psu_n=length(unique(a$psu_pool)),strata_n=length(counts),min_psu_per_stratum=min(counts),singleton_strata=sum(counts==1))}))
write.csv(design_check,file.path(root,'design_audit.csv'),row.names=FALSE)
write.csv(do.call(rbind,checks),file.path(root,'verification.csv'),row.names=FALSE)
cat(length(checks),'checks passed\n');print(design_check)
