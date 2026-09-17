invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
root<-getwd();base<-dirname(root);source(file.path(base,'KNHANES_HGS_phase2_20260910','core.R'))
wc<-function(x,n)write.csv(x,paste0(n,'.csv'),row.names=FALSE,na='',fileEncoding='UTF-8')
# Reuse only the existing metric function; never execute the prior analysis.
e<-parse(file.path(base,'HGS_family_CV_transport_20260911','run_US_transport.R'))
for(x in e)if(is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name('profile')))eval(x)
writeLines(c('BCT only; fixed M1/S2; sex-specific US 2011-2014 full-source fitting.', 'Freeze before scoring Korea 2014-2019. No recalibration, CV, or selection.', 'Both directions use full pooled target samples and target-sex weighted height quintiles.', 'Age groups 20-39/40-59/60-79. Percentages and ranges in percentage points.', 'Descriptive survey-weighted point estimates; no model-parameter uncertainty or paired tests.'),'analysis_lock.txt')
us<-readRDS(file.path(base,'HGS_pilot_20260908','harmonized_raw_qc.rds'));us$height_m<-us$height/100;us$height10<-us$height/10
old<-file.path(base,'KNHANES_HGS_pilot_20260909');kr<-rbind(readRDS(file.path(old,'development_data.rds')),readRDS(file.path(old,'validation_data.rds')));kr$weight<-kr$weight/2
us<-subset(us,eligible);kr<-subset(kr,eligible);stopifnot(nrow(us)==9334,nrow(kr)==32311)
kpath<-file.path(base,'KNHANES_pooled_BCT_20260910','pooled_BCT_models.rds');kh<-unname(tools::md5sum(kpath));km<-readRDS(kpath)
um<-list();diag<-list()
for(s in c('Female','Male')){d<-subset(us,sex==s);b<-fit_model(d,'BCT','M1','S2');um[[s]]<-b;diag[[s]]<-data.frame(sex=s,N=nrow(d),converged=TRUE,iterations=b$iterations,nu=unname(b$coefs$nu),tau=exp(unname(b$coefs$tau)),BIC=b$BIC,warnings=b$warnings);saveRDS(um,'final_US_BCT_models.rds');cat(s,'US FIT COMPLETE\n');flush.console()}
wc(do.call(rbind,diag),'US_fit_diagnostics');uh<-unname(tools::md5sum('final_US_BCT_models.rds'));um<-readRDS('final_US_BCT_models.rds')
out<-list();pr<-list();ag<-list();hg<-list();sc<-list();support<-list()
for(direction in c('US to Korea','Korea to US'))for(s in c('Female','Male')){
 d<-subset(if(direction=='US to Korea')kr else us,sex==s);tr<-subset(if(direction=='US to Korea')us else kr,sex==s);b<-if(direction=='US to Korea')um[[s]]else km[[s]]
 p<-score_model(b,d);stopifnot(all(is.finite(p$z)));m<-profile(d,p$z,wq(d$height,d$weight,c(.2,.4,.6,.8)));v<-m$summary;meta<-data.frame(direction=direction,sex=s);key<-paste(direction,s)
 out[[key]]<-cbind(meta,data.frame(N=nrow(d),mean_z=v$mean_z,sd_z=v$sd_z,age_slope_per10y=v$age_slope_per10y,height_slope_per10cm=v$slope_per10cm,p10=100*v$p10,p5=100*v$p5,p10_age_range_pp=m$profile$p10_age_range_pp,p10_height_range_pp=v$p10_range_pp,p5_height_range_pp=v$p5_range_pp,n_cdf_clipped=sum(p$clipped)))
 pr[[key]]<-cbind(meta,m$profile);ag[[key]]<-cbind(meta,m$age);h<-m$height;h[,c('p5','p10','p90','p95')]<-100*h[,c('p5','p10','p90','p95')];hg[[key]]<-cbind(meta,h)
 sc[[key]]<-data.frame(direction=direction,sex=s,id=if(direction=='US to Korea')as.character(d$uid)else as.character(d$SEQN),z=p$z,clipped=p$clipped)
 q<-wq(tr$height,tr$weight,c(.01,.99));support[[key]]<-cbind(meta,data.frame(N=nrow(d),source_height_min=min(tr$height),source_height_max=max(tr$height),source_height_p01=q[1],source_height_p99=q[2],age_outside=sum(d$age<min(tr$age)|d$age>max(tr$age)),height_below_min=sum(d$height<min(tr$height)),height_above_max=sum(d$height>max(tr$height)),height_below_p01=sum(d$height<q[1]),height_above_p99=sum(d$height>q[2])))
}
a<-do.call(rbind,out);pp<-do.call(rbind,pr);wc(a,'reciprocal_overall');wc(pp,'reciprocal_nine_metric_profile');wc(do.call(rbind,ag),'reciprocal_age_groups');wc(do.call(rbind,hg),'reciprocal_height_quintiles');wc(do.call(rbind,support),'reciprocal_support');saveRDS(do.call(rbind,sc),'frozen_transport_scores.rds')
stopifnot(kh==unname(tools::md5sum(kpath)),uh==unname(tools::md5sum('final_US_BCT_models.rds')))
wc(data.frame(source=c('Korea','US'),md5=c(kh,uh),unchanged_after_scoring=TRUE),'freeze_integrity')
met<-names(pp)[-(1:2)];l<-do.call(rbind,lapply(met,function(nm)data.frame(pp[,1:2],metric=nm,value=pp[[nm]])));l$metric<-factor(l$metric,levels=met);mx<-tapply(l$value,l$metric,max);l$intensity<-l$value/pmax(as.numeric(mx[as.character(l$metric)]),1e-12);l$label<-ifelse(l$metric %in% met[1:4],sprintf('%.3f',l$value),sprintf('%.2f',l$value))
p<-ggplot(l,aes(metric,direction,fill=intensity))+geom_tile(color='white')+geom_text(aes(label=label,color=intensity>.6),size=4,show.legend=FALSE)+scale_color_manual(values=c('FALSE'='#172A3A','TRUE'='white'))+scale_fill_gradient(low='#F2F6FA',high='#125384',name='Within-metric\nrelative magnitude')+scale_x_discrete(labels=c('|mean Z|','|SD-1|','|age slope|\nper 10y','|height slope|\nper 10cm','|P10-10%|\npp','|P5-5%|\npp','P10 age\nrange pp','P10 height\nrange pp','P5 height\nrange pp'),position='top')+facet_wrap(~sex,ncol=1)+theme_minimal(base_size=12)+theme(panel.grid=element_blank(),axis.title=element_blank(),plot.caption=element_text(hjust=0))+labs(title='Frozen BCT: bidirectional population transport',subtitle='US 2011-2014 and Korea 2014-2019; same fixed structure; no target recalibration',caption='Full pooled source fit -> full pooled target evaluation in both directions. Target-specific weighted height quintiles.\nSmaller absolute deviations/ranges indicate closer calibration; colors scaled per metric, with no composite score.\nDescriptive point estimates; source countries and survey periods differ. Neither direction is an untouched validation dataset.')
ggsave('reciprocal_transport_profile.png',p,width=16,height=6,dpi=180,bg='white')
print(a,row.names=FALSE);cat('RECIPROCAL COMPLETE\n')