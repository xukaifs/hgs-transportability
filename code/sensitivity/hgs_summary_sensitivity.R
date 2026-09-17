invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
args<-commandArgs(TRUE)
base<-if(length(args))normalizePath(args[1],winslash='/')else dirname(getwd())
out<-if(length(args)>1)args[2]else getwd()
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source(file.path(base,'KNHANES_HGS_phase2_20260910','core.R'))
wc<-function(x,n)write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='',fileEncoding='UTF-8')
writeLines(c('BCT M1/S2 only, unchanged df; refit source for bilateral mean outcome.', 'Each hand requires at least one valid trial under original country QC. Existing eligibility also retained.', 'US 2011-2014 and Korea 2014-2019 source models frozen before target evaluation.', 'US adaptation uses 2011-2012 only; location and population-weighted SD. Evaluate 2013-2014 only.', 'No CV, family selection, or target height/age recalibration. Slopes are marginal, not mutually adjusted.'),file.path(out,'analysis_lock.txt'))
kp<-file.path(base,'KNHANES_HGS_pilot_20260909')
k<-rbind(readRDS(file.path(kp,'development_data.rds')),readRDS(file.path(kp,'validation_data.rds')));k$weight<-k$weight/2
u<-readRDS(file.path(base,'HGS_pilot_20260908','harmonized_raw_qc.rds'));u$weight<-u$weight/2
mx<-function(x)apply(x,1,function(a)if(all(is.na(a)))NA_real_ else max(a,na.rm=TRUE))
kt<-paste0('gs_mea_',rep(c('r','l'),each=3),'_',rep(1:3,2));a<-as.matrix(k[,kt]);a[!is.na(a)&(!is.finite(a)|a<0|a>100)]<-NA_real_
k$hand1<-mx(a[,1:3]);k$hand2<-mx(a[,4:6]);stopifnot(isTRUE(all.equal(mx(a),k$hgs,check.attributes=FALSE)))
ut<-grep('^MGXH[12]T[123]$',names(u),value=TRUE)
a<-sapply(ut,function(t)ifelse(!is.na(u[[paste0(t,'E')]])&u[[paste0(t,'E')]]==1&u[[t]]>=0,u[[t]],NA_real_))
u$hand1<-mx(a[,grepl('H1',ut),drop=FALSE]);u$hand2<-mx(a[,grepl('H2',ut),drop=FALSE]);stopifnot(isTRUE(all.equal(mx(a),u$hgs,check.attributes=FALSE)))
flow<-list()
for(country in c('Korea','US')){
 d<-if(country=='Korea')k else u
 d$bilateral_ok<-d$eligible&is.finite(d$hand1)&is.finite(d$hand2)&(d$hand1+d$hand2)>0
 for(s in c('Female','Male'))for(period in if(country=='US')c('G','H')else '2014-2019'){
  ii<-d$sex==s;if(country=='US')ii<-ii&d$cycle==period
  ii[is.na(ii)]<-FALSE;flow[[paste(country,s,period)]]<-data.frame(country=country,sex=s,period=period,N_original=sum(d$eligible&ii),N_bilateral=sum(d$bilateral_ok&ii),N_excluded=sum(d$eligible&ii&!d$bilateral_ok))
 }
 d$hgs_original<-d$hgs;d$hgs<-(d$hand1+d$hand2)/2;d$height10<-d$height/10;d$height_m<-d$height/100
 d<-subset(d,bilateral_ok);stopifnot(all(d$hgs<=d$hgs_original+1e-10))
 if(country=='Korea')k<-d else u<-d
}
wc(do.call(rbind,flow),'sensitivity_sample_flow')
oldpaths<-c(file.path(base,'KNHANES_pooled_BCT_20260910','pooled_BCT_models.rds'),file.path(base,'BCT_reciprocal_transport_20260911','final_US_BCT_models.rds'));oldhash<-unname(tools::md5sum(oldpaths))
mods<-list();diagnostics<-list()
for(country in c('Korea','US'))for(s in c('Female','Male')){
 d<-subset(if(country=='Korea')k else u,sex==s);key<-paste(country,s,sep='_');b<-fit_model(d,'BCT','M1','S2');mods[[key]]<-b
 diagnostics[[key]]<-data.frame(country=country,sex=s,N=nrow(d),converged=TRUE,iterations=b$iterations,nu=unname(b$coefs$nu),tau=exp(unname(b$coefs$tau)),warnings=b$warnings)
 cat('FIT COMPLETE',key,'\n');flush.console()
}
mp<-file.path(out,'frozen_bilateral_BCT_models.rds');saveRDS(mods,mp);mh<-unname(tools::md5sum(mp));mods<-readRDS(mp);wc(do.call(rbind,diagnostics),'sensitivity_fit_diagnostics')
adapt<-list()
for(s in c('Female','Male')){d<-subset(u,sex==s&cycle=='G');p<-score_model(mods[[paste('Korea',s,sep='_')]],d);stopifnot(!any(p$clipped));adapt[[s]]<-data.frame(sex=s,N=nrow(d),center=weighted.mean(p$z,d$weight),scale=wsd(p$z,d$weight))}
ad<-do.call(rbind,adapt);wc(ad,'sensitivity_adaptation_parameters');saveRDS(ad,file.path(out,'frozen_bilateral_adaptation.rds'))
rows<-list()
for(scene in c('Korea_to_US_frozen','US_to_Korea_frozen','US_heldout_LS'))for(s in c('Female','Male')){
 d<-subset(if(scene=='US_to_Korea_frozen')k else u,sex==s);if(scene=='US_heldout_LS')d<-subset(d,cycle=='H')
 b<-mods[[paste(if(scene=='US_to_Korea_frozen')'US'else'Korea',s,sep='_')]];p<-score_model(b,d);z<-p$z
 if(scene=='US_heldout_LS'){aa<-subset(ad,sex==s);z<-(z-aa$center)/aa$scale}
 sl<-function(x)sum(d$weight*(x-weighted.mean(x,d$weight))*(z-weighted.mean(z,d$weight)))/sum(d$weight*(x-weighted.mean(x,d$weight))^2)
 rows[[paste(scene,s)]]<-data.frame(scenario=scene,sex=s,definition='bilateral_mean',N=nrow(d),mean_z=weighted.mean(z,d$weight),sd_z=wsd(z,d$weight),age_slope_per10y=sl(d$age/10),height_slope_per10cm=sl(d$height/10),p10=100*weighted.mean(z< -1.282,d$weight),p5=100*weighted.mean(z< -1.645,d$weight),n_cdf_clipped=sum(p$clipped))
}
r<-do.call(rbind,rows);wc(r,'sensitivity_hgs_definition')
# Existing HGSmax results only; no extra reference fits.
ref<-read.csv(file.path(base,'BCT_reciprocal_transport_20260911','reciprocal_overall.csv'));ref$scenario<-ifelse(ref$direction=='US to Korea','US_to_Korea_frozen','Korea_to_US_frozen')
h<-read.csv(file.path(base,'HGS_family_CV_transport_20260911','US_heldout_raw_diagnostics.csv'));h<-subset(h,family=='BCT');h$scenario<-'US_heldout_LS';h$height_slope_per10cm<-h$slope_per10cm;h$p10<-100*h$p10;h$p5<-100*h$p5;h$n_cdf_clipped<-0
cols<-setdiff(names(r),'definition');ref<-rbind(ref[,cols],h[,cols]);ref$definition<-'HGSmax';ref<-ref[,names(r)]
wc(rbind(ref,r),'sensitivity_vs_primary')
stopifnot(identical(oldhash,unname(tools::md5sum(oldpaths))),mh==unname(tools::md5sum(mp)),all(r$n_cdf_clipped==0),nrow(r)==6)
wc(data.frame(file=c(basename(oldpaths),basename(mp)),md5=c(oldhash,mh),unchanged=TRUE),'sensitivity_model_integrity')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
print(r,row.names=FALSE);cat('SENSITIVITY COMPLETE\n')
