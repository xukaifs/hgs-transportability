options(stringsAsFactors=FALSE)
invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
library(splines)
script_path<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
root<-dirname(normalizePath(script_path,winslash='/'));old<-file.path(dirname(root),'KNHANES_HGS_pilot_20260909')
dev<-readRDS(file.path(old,'development_data.rds'));val<-readRDS(file.path(old,'validation_data.rds'))
baseline<-readRDS(file.path(old,'frozen_Korea_2014_2016_models.rds'));models<-readRDS(file.path(root,'frozen_Gaussian_location_scale_models.rds'))
sc<-readRDS(file.path(root,'validation_scores_seven_methods.rds'));params<-readRDS(file.path(root,'validation_mu_sigma.rds'))
ss<-read.csv(file.path(root,'summary_validation.csv'));qq<-read.csv(file.path(root,'height_quintiles.csv'))
checks<-list();add<-function(name,pass,error=0){stopifnot(pass);checks[[length(checks)+1]]<<-data.frame(check=name,pass=pass,error=error)}
for(s in c('Female','Male')){
 b<-models[[s]];v<-subset(val,eligible & sex==s);d<-subset(dev,eligible & sex==s);oldfit<-baseline[[paste0(s,'_Conditional')]]$model
 # Independent basis assembly from the original stored ns objects.
 Xv<-cbind(1,predict(oldfit$model[[2]],newx=v$age),predict(oldfit$model[[3]],newx=v$height))
 Xd<-cbind(1,predict(oldfit$model[[2]],newx=d$age),predict(oldfit$model[[3]],newx=d$height))
 mu<-as.numeric(Xv%*%b$mu_coefficients);sigma<-exp(as.numeric(Xv%*%b$log_sigma_coefficients));z<-(v$hgs-mu)/sigma
 a<-sc[sc$sex==s & sc$method=='Conditional_heightSD',];a<-a[match(v$uid,a$uid),]
 error<-max(abs(z-a$z));add(paste(s,'independent ns scoring'),error<1e-10,error)
 ix<-c(1,round(nrow(v)/2),nrow(v));zs<-sapply(ix,function(i){XX<-model.matrix(b$terms,v[i,,drop=FALSE]);(v$hgs[i]-as.numeric(XX%*%b$mu_coefficients))/exp(as.numeric(XX%*%b$log_sigma_coefficients))})
 error<-max(abs(zs-z[ix]));add(paste(s,'single batch consistency'),error<1e-10,error)
 add(paste(s,'development validation separation'),!any(v$uid %in% b$development_uids)&&identical(sort(d$uid),sort(b$development_uids)))
 w<-d$weight/mean(d$weight);p<-ncol(Xd);S<-Xd[,b$sigma_column_indices,drop=FALSE];q<-ncol(S);theta<-c(b$mu_coefficients,b$log_sigma_coefficients[b$sigma_column_indices])
 add(paste(s,'no age in sigma'),all(b$log_sigma_coefficients[-b$sigma_column_indices]==0))
 fn<-function(t){mu<-as.numeric(Xd%*%t[seq_len(p)]);ls<-as.numeric(S%*%t[p+seq_len(q)]);sum(w*(ls+.5*(d$hgs-mu)^2*exp(-2*ls)))}
 H<-optimHess(theta,fn);ev<-eigen(H,symmetric=TRUE,only.values=TRUE)$values
 add(paste(s,'positive objective curvature at solution'),min(ev)>0,min(ev))
 r<-ss[ss$sex==s & ss$method=='Conditional_heightSD',];slope<-sum(a$weight*(a$height_cm/10-weighted.mean(a$height_cm/10,a$weight))*(a$z-weighted.mean(a$z,a$weight)))/sum(a$weight*(a$height_cm/10-weighted.mean(a$height_cm/10,a$weight))^2)
 add(paste(s,'independent slope'),abs(slope-r$slope_per10cm)<1e-10,abs(slope-r$slope_per10cm))
 for(q in paste0('Q',1:5)){
  aa<-a[a$quintile==q,];rr<-qq[qq$sex==s & qq$method=='Conditional_heightSD' & qq$quintile==q,]
  sd<-sqrt(sum(aa$weight*(aa$z-sum(aa$weight*aa$z)/sum(aa$weight))^2)/sum(aa$weight))
  rates<-sapply(list(aa$z< -1.645,aa$z< -1.282,aa$z>1.282,aa$z>1.645),function(flag)100*sum(aa$weight[flag])/sum(aa$weight))
  error<-max(abs(rates-unlist(rr[,c('p5_pct','p10_pct','p90_pct','p95_pct')])),abs(sd-rr$sd_z))
  add(paste(s,q,'independent SD tails'),error<1e-9&&nrow(aa)==rr$n,error)
 }
}
add('All predictions and primary results finite',all(is.finite(params$mu))&&all(is.finite(params$sigma))&&all(params$sigma>0)&&all(vapply(ss[sapply(ss,is.numeric)],function(a)all(is.finite(a)),logical(1))))
write.csv(do.call(rbind,checks),file.path(root,'independent_verification.csv'),row.names=FALSE)
cat(length(checks),'independent checks passed\n')
