invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
script_path<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]);root<-dirname(normalizePath(script_path,winslash='/'))
source(file.path(root,'core.R'))
old<-file.path(dirname(root),'KNHANES_HGS_pilot_20260909')
wc<-function(x,n)write.csv(x,file.path(root,paste0(n,'.csv')),row.names=FALSE,na='',fileEncoding='UTF-8')
dir.create(file.path(root,'candidates'),showWarnings=FALSE)
all<-readRDS(file.path(old,'development_data.rds'));stopifnot(all(all$year %in% 2014:2016))
dat<-subset(all,eligible);set.seed(20260910)
ps<-unique(dat[,c('psu_pool','strata_pool')]);folds<-list()
for(rep in 1:2){
 v<-integer(nrow(ps))
 for(st in unique(ps$strata_pool)){i<-which(ps$strata_pool==st);offset<-sample(0:4,1);v[i[sample.int(length(i))]]<-((seq_along(i)-1+offset)%%5)+1}
 folds[[rep]]<-data.frame(psu_pool=ps$psu_pool,repeat_id=rep,fold=v)
}
foldtab<-do.call(rbind,folds);wc(foldtab,'psu_cv_assignments');saveRDS(foldtab,file.path(root,'psu_cv_assignments.rds'))
totals<-list();details<-list();quints<-list();fitlog<-list();saved<-list();errors<-list();gate<-list()
run_candidate<-function(sex,family='NO',mu='M1',scale='S1',traditional=NULL){
 id<-if(is.null(traditional))paste(family,mu,scale,sep='_') else traditional;key<-paste(sex,id,sep='_')
 cache<-file.path(root,'candidates',paste0(key,'.rds'))
 if(file.exists(cache)){r<-readRDS(cache);if(isTRUE(r)){cat('RESUME',key,'\n');return(r)}}
 d<-dat[dat$sex==sex,];edges<-wq(d$height,d$weight,c(.2,.4,.6,.8));oo<-list();logs<-list();valid<-TRUE
 for(rr in 1:2){
  fld<-folds[[rr]]$fold[match(d$psu_pool,folds[[rr]]$psu_pool)];z<-rep(NA_real_,nrow(d));ld<-z;clip<-rep(FALSE,nrow(d))
  for(k in 1:5){
   train<-d[fld!=k,];test<-d[fld==k,];stopifnot(!any(train$psu_pool %in% test$psu_pool));t0<-proc.time()[3]
   ans<-tryCatch({b<-if(is.null(traditional))fit_model(train,family,mu,scale) else fit_traditional(train,traditional);p<-score_model(b,test);list(b=b,p=p)},error=function(e)e)
   ok<-!inherits(ans,'error');msg<-if(ok)if(is.null(ans$b$warnings))'' else ans$b$warnings else conditionMessage(ans)
   logs[[length(logs)+1]]<-data.frame(sex=sex,candidate=id,repeat_id=rr,fold=k,n_training=nrow(train),n_testing=nrow(test),ok=ok,seconds=proc.time()[3]-t0,message=msg)
   if(ok){z[fld==k]<-ans$p$z;ld[fld==k]<-ans$p$logdens;clip[fld==k]<-ans$p$clipped}else valid<-FALSE
   cat(key,'repeat',rr,'fold',k,if(ok)'OK' else msg,'\n');flush.console()
  }
  oo[[rr]]<-data.frame(uid=d$uid,repeat_id=rr,fold=fld,z=z,logdensity=ld,clipped=clip)
 }
 full<-tryCatch(if(is.null(traditional))fit_model(d,family,mu,scale) else fit_traditional(d,traditional),error=function(e)e)
 fullok<-!inherits(full,'error');valid<-valid&&fullok
 sm<-list();qt<-list()
 if(valid){for(rr in 1:2){m<-metrics(d,oo[[rr]]$z,edges);sm[[rr]]<-cbind(sex=sex,candidate=id,repeat_id=rr,m$summary,clipped=sum(oo[[rr]]$clipped),mean_nll=-weighted.mean(oo[[rr]]$logdensity,d$weight));qt[[rr]]<-cbind(sex=sex,candidate=id,repeat_id=rr,m$quintiles)}
  fp<-score_model(full,d);fm<-metrics(d,fp$z,edges)$summary
  ds<-data.frame(sex=sex,candidate=id,family=if(is.null(traditional))family else 'traditional',mu=mu,scale=scale,valid=TRUE,df=full$df,BIC=if(is.null(full$BIC))NA_real_ else full$BIC,GAIC=if(is.null(full$GAIC))NA_real_ else full$GAIC,tail_loss=mean(vapply(sm,function(a)a$tail_loss,numeric(1))),guardrail=all(vapply(sm,function(a)a$guardrail&&a$clipped==0,logical(1))),development_skewness=fm$skewness,development_excess_kurtosis=fm$excess_kurtosis)
 }else ds<-data.frame(sex=sex,candidate=id,family=if(is.null(traditional))family else 'traditional',mu=mu,scale=scale,valid=FALSE,df=NA_real_,BIC=NA_real_,GAIC=NA_real_,tail_loss=NA_real_,guardrail=FALSE,development_skewness=NA_real_,development_excess_kurtosis=NA_real_)
 r<-list(summary=ds,repeats=if(length(sm))do.call(rbind,sm)else NULL,quintiles=if(length(qt))do.call(rbind,qt)else NULL,fitlog=do.call(rbind,logs),oof=do.call(rbind,oo),model=if(fullok)full else NULL,full_error=if(fullok)'' else conditionMessage(full),edges=edges)
 saveRDS(r,cache);cat('CANDIDATE DONE',key,'valid',valid,'tail loss',ds$tail_loss,'\n');flush.console();r
}
choose<-function(rs){
 t<-do.call(rbind,lapply(rs,`[[`,'summary'));t<-t[t$valid,];stopifnot(nrow(t)>0)
 accepted<-t[t$guardrail,];if(!nrow(accepted))accepted<-t
 best<-min(accepted$tail_loss);near<-accepted[accepted$tail_loss<=best*1.05+1e-12,];near<-near[order(near$df,near$BIC,near$candidate),];near$candidate[1]
}
final<-list();selection<-list();allresults<-list()
for(s in c('Male','Female')){
 rs<-list()
 for(fam in c('NO','BCCG','BCT'))for(sc in c('S0','S1','S2')){r<-run_candidate(s,fam,'M1',sc);rs[[r$summary$candidate]]<-r}
 trigger<-TRUE;gaterows<-list()
 for(fam in c('BCCG','BCT')){
  candidates<-rs[vapply(rs,function(a)a$summary$family==fam&&a$summary$valid,logical(1))]
  if(!length(candidates)){trigger<-FALSE;gaterows[[fam]]<-data.frame(sex=s,family=fam,best=NA_character_,kurtosis_r1=NA_real_,kurtosis_r2=NA_real_,thick=FALSE)}else{
   best<-choose(candidates);ku<-candidates[[best]]$repeats$excess_kurtosis;thick<-all(ku>.5);trigger<-trigger&&thick;gaterows[[fam]]<-data.frame(sex=s,family=fam,best=best,kurtosis_r1=ku[1],kurtosis_r2=ku[2],thick=thick)
  }
 }
 gate[[s]]<-do.call(rbind,gaterows);gate[[s]]$BCPE_triggered<-trigger
 if(trigger)for(sc in c('S0','S1','S2')){r<-run_candidate(s,'BCPE','M1',sc);rs[[r$summary$candidate]]<-r}
 provisional<-choose(rs)
 if(s=='Female'){
  selected<-rs[[provisional]]$summary
  for(mu in c('M2','M3')){r<-run_candidate(s,selected$family,mu,selected$scale);rs[[r$summary$candidate]]<-r}
 }
 winner<-choose(rs);simple<-if(winner=='NO_M1_S1')'NO_M1_S0' else 'NO_M1_S1'
 stopifnot(rs[[simple]]$summary$valid)
 final[[s]]<-list(primary=rs[[winner]]$model,simplified=rs[[simple]]$model)
 selection[[s]]<-data.frame(sex=s,role=c('primary','simplified'),candidate=c(winner,simple),internal_guardrail=c(rs[[winner]]$summary$guardrail,rs[[simple]]$summary$guardrail),tail_loss=c(rs[[winner]]$summary$tail_loss,rs[[simple]]$summary$tail_loss),df=c(rs[[winner]]$summary$df,rs[[simple]]$summary$df))
 for(method in c('Absolute','BMI','Height2','Allometry','Allometry_ageadj')){r<-run_candidate(s,traditional=method);rs[[method]]<-r;if(method=='Allometry_ageadj')final[[s]]$age_adjusted_allometry<-r$model}
 allresults[[s]]<-rs
 saveRDS(allresults,file.path(root,'development_results_checkpoint.rds'))
 cat('SEX LOCKED',s,winner,'simplified',simple,'\n')
}
flat<-unlist(allresults,recursive=FALSE)
wc(do.call(rbind,lapply(flat,`[[`,'summary')),'development_candidate_summary')
wc(do.call(rbind,lapply(flat,`[[`,'repeats')),'cv_repeat_metrics')
wc(do.call(rbind,lapply(flat,`[[`,'quintiles')),'cv_height_quintiles')
wc(do.call(rbind,lapply(flat,`[[`,'fitlog')),'cv_fit_log')
wc(do.call(rbind,gate),'BCPE_gate')
wc(do.call(rbind,selection),'locked_selection')
# Development-only survey exponent inference.
design<-svydesign(ids=~psu_pool,strata=~strata_pool,weights=~weight,data=subset(all,design_ok),nest=TRUE)
exps<-list()
for(s in c('Female','Male')){de<-subset(design,eligible & sex==s);fit<-svyglm(log(hgs)~ns(age,4)+log(height_m),de);nm<-'log(height_m)';ci<-confint(fit)[nm,];exps[[s]]<-data.frame(sex=s,b=coef(fit)[nm],SE=sqrt(vcov(fit)[nm,nm]),ci_low=ci[1],ci_high=ci[2]);stopifnot(abs(coef(fit)[nm]-final[[s]]$age_adjusted_allometry$b)<1e-8)}
wc(do.call(rbind,exps),'age_adjusted_allometry_exponents')
saveRDS(final,file.path(root,'locked_development_models.rds'))
wc(data.frame(file='locked_development_models.rds',md5=unname(tools::md5sum(file.path(root,'locked_development_models.rds'))),locked_at=as.character(Sys.time()),training_years='2014-2016',temporal_loaded=FALSE),'model_lock_manifest')
saveRDS(allresults,file.path(root,'development_results.rds'))
capture.output(sessionInfo(),file=file.path(root,'selection_sessionInfo.txt'))
cat('DEVELOPMENT SELECTION COMPLETE. TEMPORAL DATA NOT LOADED.\n')
