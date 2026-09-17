invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
.libPaths(c('C:/Users/Public/CodexRLib42Copy',.libPaths()));suppressPackageStartupMessages(library(ggplot2))
root<-getwd();base<-dirname(root);phase<-file.path(base,'KNHANES_HGS_phase2_20260910');usdir<-file.path(base,'US_simple_methods_ladder_20260910')
inputs<-c(file.path(phase,'cv_repeat_metrics.csv'),file.path(phase,'summary_validation.csv'),file.path(usdir,'heldout_normalization_ladder_comparison.csv'))
cv<-read.csv(inputs[1]);tv<-read.csv(inputs[2]);us<-read.csv(inputs[3]);wc<-function(x,n)write.csv(x,paste0(n,'.csv'),row.names=FALSE,na='',fileEncoding='UTF-8')
metric_names<-c('abs_mean_z','abs_sd_minus1','abs_slope_per10cm','abs_p10_minus10_pp','abs_p5_minus5_pp','p10_height_range_pp','p5_height_range_pp')
calc<-function(mean,sd,slope,p10,p5,r10,r5)data.frame(abs_mean_z=abs(mean),abs_sd_minus1=abs(sd-1),abs_slope_per10cm=abs(slope),abs_p10_minus10_pp=abs(p10-10),abs_p5_minus5_pp=abs(p5-5),p10_height_range_pp=r10,p5_height_range_pp=r5)
rows<-list();repeatrows<-list()
for(s in c('Female','Male'))for(m in c('Height2','Allometry_ageadj','BCT')){
 id<-if(m=='BCT')'BCT_M1_S2'else m
 a<-subset(cv,sex==s&candidate==id);stopifnot(nrow(a)==2,setequal(a$repeat_id,1:2))
 aa<-calc(a$mean_z,a$sd_z,a$slope_per10cm,100*a$p10,100*a$p5,a$p10_range_pp,a$p5_range_pp)
 repeatrows[[paste(s,m)]]<-cbind(sex=s,method=m,repeat_id=a$repeat_id,aa)
 rows[[paste(s,m,'A')]]<-cbind(sex=s,scenario='A',method=m,N=a$n[1],source_reference='Korea 2014-2016, training-fold fits',evaluation='Korea 2014-2016, mean of 2 PSU-grouped OOF repeats',adaptation='None',as.data.frame(as.list(colMeans(aa))))
 tid<-switch(m,BCT='primary',Allometry_ageadj='age_adjusted_allometry',Height2='Height2');a<-subset(tv,sex==s&method==tid);stopifnot(nrow(a)==1)
 rows[[paste(s,m,'B')]]<-cbind(sex=s,scenario='B',method=m,N=a$n_validation,source_reference='Korea 2014-2016, frozen full development fit',evaluation='Korea 2017-2019 temporal evaluation',adaptation='None',calc(a$mean_z,a$sd_z,a$slope_per10cm,a$p10_overall_pct,a$p5_overall_pct,a$p10_quintile_range_pp,a$p5_quintile_range_pp))
 a<-subset(us,sex==s&normalization==m&method=='Location_scale');stopifnot(nrow(a)==1)
 rows[[paste(s,m,'C')]]<-cbind(sex=s,scenario='C',method=m,N=a$N,source_reference='Korea 2014-2019, frozen pooled fit',evaluation='NHANES 2013-2014 held-out',adaptation='NHANES 2011-2012 location + scale only',calc(a$mean_z,a$sd_z,a$slope_10cm,a$p10,a$p5,a$p10_range_pp,a$p5_range_pp))
}
pf<-do.call(rbind,rows);pf<-pf[order(pf$sex,pf$scenario,match(pf$method,c('Height2','Allometry_ageadj','BCT'))),];rownames(pf)<-NULL
wc(pf,'performance_profile');wc(do.call(rbind,repeatrows),'OOF_repeat_profiles');wc(data.frame(file=inputs,md5=unname(tools::md5sum(inputs))),'input_provenance')
long<-do.call(rbind,lapply(metric_names,function(nm)data.frame(pf[,c('sex','scenario','method')],metric=nm,value=pf[[nm]])))
long$metric<-factor(long$metric,levels=metric_names);mx<-tapply(long$value,long$metric,max);long$within_metric_fraction<-long$value/as.numeric(mx[as.character(long$metric)]);long$display<-ifelse(long$metric %in% metric_names[1:3],sprintf('%.3f',long$value),sprintf('%.2f',long$value))
long$display[long$metric %in% metric_names[1:3] & long$value>0 & long$value<.0005]<-"<0.001";wc(long,"performance_profile_heatmap_data")
method_label<-c(Height2='Height squared',Allometry_ageadj='Age-adjusted allometry',BCT='BCT')
long$row<-paste(long$scenario,method_label[long$method],sep=' | ');rowlevels<-unlist(lapply(c('A','B','C'),function(sc)paste(sc,method_label,sep=' | ')));long$row<-factor(long$row,levels=rev(rowlevels))
labels<-c('|mean Z|','|SD - 1|','|height slope|\nZ / 10 cm','|P10 - 10%|\npp','|P5 - 5%|\npp','P10 range\npp','P5 range\npp')
p<-ggplot(long,aes(metric,row,fill=within_metric_fraction))+geom_tile(color='white',linewidth=.7)+geom_text(aes(label=display,color=within_metric_fraction>.58),size=3.7,show.legend=FALSE)+scale_color_manual(values=c('FALSE'='#172A3A','TRUE'='white'))+scale_fill_gradient(low='#F2F6FA',high='#125384',limits=c(0,1),name='Within-metric\nrelative magnitude',breaks=c(0,.5,1),labels=c('0','0.5','1'))+scale_x_discrete(labels=labels,position='top')+facet_wrap(~sex,ncol=1)+theme_minimal(base_size=12)+theme(panel.grid=element_blank(),axis.title=element_blank(),axis.text.y=element_text(size=11),axis.text.x=element_text(size=10),strip.text=element_text(face='bold',size=13),legend.position='right',plot.title=element_text(face='bold',size=17),plot.caption=element_text(hjust=0,size=10),panel.spacing=grid::unit(1,'lines'))+labs(title='Handgrip normalization: multidimensional performance profile',subtitle='A: Korean OOF   |   B: Korean temporal evaluation   |   C: US held-out after location + scale adaptation',caption=paste('A/B reference: Korea 2014-2016. C reference: Korea 2014-2019; adaptation: US 2011-2012; evaluation: US 2013-2014.',
'Cells show absolute deviations or height-quintile ranges in original units; smaller is closer to the stated target.',
'Color is scaled separately for each metric across both sexes and all scenarios; no composite score or clinical thresholds.',
'A averages the two repeat-specific deviations/ranges. Source years, calibration and height cutpoints differ across scenarios.',sep='\n'))
ggsave('performance_profile_heatmap.png',p,width=14,height=10,dpi=240,bg='white');ggsave('performance_profile_heatmap.pdf',p,width=14,height=10,bg='white')
stopifnot(nrow(pf)==18,nrow(long)==126,all(is.finite(as.matrix(pf[,metric_names]))),all(as.matrix(pf[,metric_names])>=0),!anyDuplicated(pf[,c('sex','scenario','method')]),identical(unname(tools::md5sum(inputs)),read.csv('input_provenance.csv')$md5))
writeLines(c('No models fitted. All values copied or algebraically transformed from existing results.',
'A: repeat-specific absolute deviations calculated first, then averaged over two OOF repeats; participants not doubled.',
'B: 2014-2016 frozen references evaluated on 2017-2019. C: pooled 2014-2019 references adapted only with 2011-2012 US location and scale, evaluated in 2013-2014.',
'C is uniformly location + scale, not each method best-performing calibration level.',
'A/B/C do not estimate an isolated country effect: training years, recalibration, sample composition and height quintiles differ.',
'OOF performance is not nested model-selection validation. Temporal and US samples have been viewed previously; no untouched-validation claim.',
'Color=value / metric-specific maximum across the 18 cells. Never compare color intensity between different metrics as an absolute effect size.',
'No composite score, clinical cutoffs, significance tests, or overall winner. No claims of equivalence from similar point estimates.',
'7 metrics, 3 methods, 3 scenarios, 2 sexes: 126 numerical cells; input hashes unchanged.'),'figure_notes.txt')
cat('PROFILE COMPLETE: 18 rows, 126 cells; no refitting; input hashes unchanged.\n');print(pf[,c('sex','scenario','method',metric_names)],row.names=FALSE)
