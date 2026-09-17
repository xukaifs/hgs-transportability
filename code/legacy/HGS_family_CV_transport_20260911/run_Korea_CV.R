invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
root<-getwd();base<-dirname(root);source(file.path(base,'KNHANES_HGS_phase2_20260910','core.R'))
wc<-function(x,n)write.csv(x,paste0(n,'.csv'),row.names=FALSE,na='',fileEncoding='UTF-8');dir.create('cache',showWarnings=FALSE)
old<-file.path(base,'KNHANES_HGS_pilot_20260909');full<-rbind(readRDS(file.path(old,'development_data.rds')),readRDS(file.path(old,'validation_data.rds')));full$weight<-full$weight/2;d0<-subset(full,eligible);stopifnot(nrow(d0)==32311,!anyDuplicated(d0$uid))
oldb<-readRDS(file.path(base,'KNHANES_pooled_BCT_20260910','pooled_BCT_models.rds'));oldscores<-readRDS(file.path(base,'KNHANES_pooled_BCT_20260910','pooled_BCT_scores.rds'));stopifnot(setequal(oldscores$uid,d0$uid))
ps<-unique(d0[,c('psu_pool','strata_pool')]);set.seed(20260911);folds<-list()
for(rr in 1:2){v<-integer(nrow(ps));for(st in unique(ps$strata_pool)){i<-which(ps$strata_pool==st);off<-sample(0:4,1);v[i[sample.int(length(i))]]<-((seq_along(i)-1+off)%%5)+1};folds[[rr]]<-data.frame(psu_pool=ps$psu_pool,repeat_id=rr,fold=v)}
wc(do.call(rbind,folds),'PSU_fold_assignment');saveRDS(folds,'PSU_fold_assignment.rds')
profile<-function(d,z,edges){m<-metrics(d,z,edges);w<-d$weight;age<-d$age/10;aslope<-sum(w*(age-weighted.mean(age,w))*(z-weighted.mean(z,w)))/sum(w*(age-weighted.mean(age,w))^2);gr<-cut(d$age,c(20,40,60,80),right=FALSE,labels=c('20-39','40-59','60-79'));ag<-do.call(rbind,lapply(levels(gr),function(g){i<-gr==g;data.frame(age_group=g,N=sum(i),mean_z=weighted.mean(z[i],w[i]),sd_z=wsd(z[i],w[i]),p10=100*weighted.mean(z[i]< -1.282,w[i]),p5=100*weighted.mean(z[i]< -1.645,w[i]))}));s<-m$summary
 out<-data.frame(abs_mean_z=abs(s$mean_z),abs_sd_minus1=abs(s$sd_z-1),abs_age_slope_per10y=abs(aslope),abs_height_slope_per10cm=abs(s$slope_per10cm),abs_p10_minus10_pp=abs(100*s$p10-10),abs_p5_minus5_pp=abs(100*s$p5-5),p10_age_range_pp=diff(range(ag$p10)),p10_height_range_pp=s$p10_range_pp,p5_height_range_pp=s$p5_range_pp)
 list(profile=out,summary=cbind(s,age_slope_per10y=aslope),age=ag,height=m$quintiles)
}
results<-list();fitlogs<-list();models<-list();diag<-list();worms<-list()
for(s in c('Female','Male'))for(fam in c('NO','BCCG','BCT')){
 key<-paste(s,fam,sep='_');path<-file.path('cache',paste0(key,'.rds'))
 if(file.exists(path)){r<-readRDS(path);cat('RESUME',key,'\n')}else{
 d<-subset(d0,sex==s);edges<-wq(d$height,d$weight,c(.2,.4,.6,.8));alloo<-list();logs<-list()
 for(rr in 1:2){fld<-folds[[rr]]$fold[match(d$psu_pool,folds[[rr]]$psu_pool)];z<-rep(NA_real_,nrow(d));clip<-rep(FALSE,nrow(d))
 for(k in 1:5){tr<-d[fld!=k,];te<-d[fld==k,];stopifnot(!any(tr$psu_pool %in% te$psu_pool));t<-proc.time()[3]
 ans<-tryCatch({b<-fit_model(tr,fam,'M1','S2');p<-score_model(b,te);list(b=b,p=p)},error=function(e)e);ok<-!inherits(ans,'error')
 logs[[length(logs)+1]]<-data.frame(sex=s,family=fam,repeat_id=rr,fold=k,N_train=nrow(tr),N_test=nrow(te),converged=ok,seconds=proc.time()[3]-t,message=if(ok)ans$b$warnings else conditionMessage(ans),iterations=if(ok)ans$b$iterations else NA)
 if(ok){z[fld==k]<-ans$p$z;clip[fld==k]<-ans$p$clipped};cat(key,'repeat',rr,'fold',k,if(ok)'OK'else conditionMessage(ans),'\n');flush.console()
 }
 alloo[[rr]]<-data.frame(uid=d$uid,repeat_id=rr,fold=fld,z=z,clipped=clip)
 }
 logs<-do.call(rbind,logs);wc(logs,paste0('fit_log_',key));stopifnot(all(logs$converged),all(vapply(alloo,function(a)all(is.finite(a$z))&&!any(a$clipped),logical(1))))
 b<-if(fam=='BCT')oldb[[s]]else fit_model(d,fam,'M1','S2');p<-score_model(b,d);stopifnot(!any(p$clipped),b$n==nrow(d))
 if(fam=='BCT'){prev<-subset(oldscores,sex==s);stopifnot(max(abs(p$z-prev$z[match(d$uid,prev$uid)]))<1e-12)}
 r<-list(sex=s,family=fam,model=b,oof=do.call(rbind,alloo),edges=edges,fitlog=logs);saveRDS(r,path)
 }
 d<-subset(d0,sex==s);b<-r$model;p<-score_model(b,d);models[[key]]<-b;results[[key]]<-r;fitlogs[[key]]<-r$fitlog
 prof<-profile(d,p$z,r$edges);diag[[key]]<-data.frame(sex=s,family=fam,N=nrow(d),df=b$df,BIC=b$BIC,GAIC=b$GAIC,converged=TRUE,iterations=b$iterations,skewness=prof$summary$skewness,excess_kurtosis=prof$summary$excess_kurtosis,full_fit_reused=fam=='BCT')
 for(rr in 0:2){z<-if(rr==0)p$z else subset(r$oof,repeat_id==rr)$z;pp<-seq(.01,.99,.01);worms[[paste(key,rr)]]<-data.frame(sex=s,family=fam,sample=if(rr==0)'Full source fit'else paste('OOF repeat',rr),theory=qnorm(pp),observed=wq(z,d$weight,pp),difference=wq(z,d$weight,pp)-qnorm(pp))}
 cat('FAMILY COMPLETE',key,'\n');flush.console()
}
repmetrics<-list();raw<-list();age<-list();height<-list()
for(key in names(results)){r<-results[[key]];d<-subset(d0,sex==r$sex)
 for(rr in 1:2){o<-subset(r$oof,repeat_id==rr);stopifnot(identical(o$uid,d$uid));p<-profile(d,o$z,r$edges);meta<-data.frame(sex=r$sex,family=r$family,repeat_id=rr);repmetrics[[paste(key,rr)]]<-cbind(meta,p$profile);raw[[paste(key,rr)]]<-cbind(meta,p$summary);age[[paste(key,rr)]]<-cbind(meta,p$age);height[[paste(key,rr)]]<-cbind(meta,p$height)}
}
rm<-do.call(rbind,repmetrics);mn<-aggregate(rm[,4:12],rm[,c('sex','family')],mean);wc(mn,'Korea_OOF_nine_metrics');wc(rm,'Korea_OOF_repeat_metrics');wc(do.call(rbind,raw),'Korea_OOF_raw_diagnostics');wc(do.call(rbind,age),'Korea_OOF_age_groups');wc(do.call(rbind,height),'Korea_OOF_height_groups');wc(do.call(rbind,diag),'source_fit_diagnostics');wc(do.call(rbind,fitlogs),'CV_convergence');wc(do.call(rbind,worms),'overall_residual_quantiles')
saveRDS(models,'frozen_pooled_family_models.rds');wc(data.frame(file='frozen_pooled_family_models.rds',md5=unname(tools::md5sum('frozen_pooled_family_models.rds')),frozen_at=as.character(Sys.time()),primary='BCT unchanged'),'freeze_manifest')
w<-do.call(rbind,worms);p<-ggplot(w,aes(theory,difference,color=sample))+geom_hline(yintercept=0,linetype=2)+geom_line()+facet_grid(sex~family)+theme_bw(base_size=11)+theme(legend.position='bottom')+labs(title='Korean source: weighted overall worm diagnostics',x='Standard normal quantile',y='Observed minus normal quantile');ggsave('Korea_overall_worm.png',p,width=13,height=7,dpi=180,bg='white')
p<-ggplot(w,aes(theory,observed,color=sample))+geom_abline(intercept=0,slope=1,linetype=2)+geom_line()+facet_grid(sex~family)+theme_bw(base_size=11)+theme(legend.position='bottom')+labs(title='Korean source: weighted quantile residual QQ',x='Standard normal quantile',y='Observed residual quantile');ggsave('Korea_overall_QQ.png',p,width=13,height=7,dpi=180,bg='white')
stopifnot(nrow(mn)==6,nrow(rm)==12,nrow(do.call(rbind,fitlogs))==60);cat('KOREA CV COMPLETE; BCT PRIMARY UNCHANGED; NO US DATA LOADED\n');print(mn,row.names=FALSE)