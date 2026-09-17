invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
root<-getwd();source('core.R');d0<-readRDS(file.path(dirname(root),'KNHANES_HGS_pilot_20260909','development_data.rds'));d0<-subset(d0,eligible)
rs<-readRDS('development_results.rds');dir.create('diagnostics',showWarnings=FALSE);checks<-list()
for(s in names(rs))for(id in names(rs[[s]])){
 r<-rs[[s]][[id]];if(!r$summary$valid)next
 d<-subset(d0,sex==s);pp<-seq(.01,.99,.01);worm<-list();hist<-list()
 for(rr in 0:2){
  if(rr==0){z<-score_model(r$model,d)$z}else{a<-subset(r$oof,repeat_id==rr);stopifnot(!anyDuplicated(a$uid),setequal(a$uid,d$uid));z<-a$z[match(d$uid,a$uid)]}
  q<-cut(d$height,c(-Inf,r$edges,Inf),labels=paste0('Q',1:5),include.lowest=TRUE)
  for(lev in c('Overall',paste0('Q',1:5))){i<-if(lev=='Overall')rep(TRUE,nrow(d))else q==lev;worm[[length(worm)+1]]<-data.frame(group=lev,theory=qnorm(pp),difference=wq(z[i],d$weight[i],pp)-qnorm(pp),sample=if(rr==0)'Development fit' else paste('OOF repeat',rr))}
 }
 a<-do.call(rbind,worm);write.csv(a,file.path('diagnostics',paste0(s,'_',id,'_worm.csv')),row.names=FALSE)
 p<-ggplot(a,aes(theory,difference,color=sample))+geom_hline(yintercept=0,linetype=2)+geom_line()+facet_wrap(~group,ncol=3)+theme_bw(base_size=11)+theme(legend.position='bottom')+labs(title=paste(s,id,'weighted detrended QQ'),subtitle='Continuous CDF residuals; descriptive curves without independent-observation bands',x='Standard normal quantile',y='Observed minus normal quantile')
 ggsave(file.path('diagnostics',paste0(s,'_',id,'_worm.png')),p,width=12,height=7,dpi=150)
 checks[[paste(s,id)]]<-data.frame(sex=s,candidate=id,oof_rows=nrow(r$oof),expected=2*nrow(d),finite=all(is.finite(r$oof$z)),full_model=TRUE)
}
ch<-do.call(rbind,checks);stopifnot(all(ch$oof_rows==ch$expected),all(ch$finite));write.csv(ch,'development_integrity_checks.csv',row.names=FALSE)
fold<-readRDS('psu_cv_assignments.rds');stopifnot(!anyDuplicated(fold[,c('psu_pool','repeat_id')]),all(fold$fold %in% 1:5))
cat('DIAGNOSTICS AND INTEGRITY CHECKS COMPLETE\n')