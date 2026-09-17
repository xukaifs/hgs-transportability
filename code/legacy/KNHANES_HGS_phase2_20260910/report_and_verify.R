invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
root<-getwd();source('core.R');wc<-function(x,n)write.csv(x,paste0(n,'.csv'),row.names=FALSE,na='',fileEncoding='UTF-8')
ss<-read.csv('summary_validation.csv');qq<-read.csv('height_quintiles.csv');cv<-read.csv('cv_repeat_metrics.csv');ds<-read.csv('development_candidate_summary.csv');locks<-read.csv('locked_selection.csv');models<-readRDS('locked_development_models.rds');sc<-readRDS('temporal_scores.rds')
old<-file.path(dirname(root),'KNHANES_HGS_pilot_20260909');d<-readRDS(file.path(old,'development_data.rds'));v<-readRDS(file.path(old,'validation_data.rds'))
prior<-readRDS(file.path(dirname(root),'KNHANES_HGS_sigma_comparison_20260910','validation_scores_seven_methods.rds'))
stopifnot(identical(prior,sc[seq_len(nrow(prior)),]),!any(subset(d,eligible)$uid %in% subset(v,eligible)$uid))
stopifnot(all(ds$valid),nrow(ds)==30,nrow(read.csv('cv_fit_log.csv'))==300,all(read.csv('cv_fit_log.csv')$ok))
fold<-readRDS('psu_cv_assignments.rds');stopifnot(nrow(fold)==2*length(unique(subset(d,eligible)$psu_pool)),!anyDuplicated(fold[,c('psu_pool','repeat_id')]))
for(s in c('Female','Male')){
 for(m in unique(sc$method)){a<-subset(sc,sex==s&method==m);stopifnot(nrow(a)==sum(v$eligible&v$sex==s),!anyDuplicated(a$uid),all(is.finite(a$z)))}
 t<-subset(ds,sex==s&family!='traditional'&guardrail);mn<-min(t$tail_loss);a<-subset(t,tail_loss<=mn*1.05+1e-12);a<-a[order(a$df,a$BIC,a$candidate),];stopifnot(a$candidate[1]==subset(locks,sex==s&role=='primary')$candidate)
 for(id in subset(ds,sex==s)$candidate){r<-readRDS(file.path('candidates',paste0(s,'_',id,'.rds')));dv<-subset(d,eligible&sex==s);p<-score_model(r$model,dv);wc(data.frame(uid=dv$uid,z=p$z,weight=dv$weight,clipped=p$clipped),file.path('diagnostics',paste0(s,'_',id,'_development_residuals')))}
}
man<-read.csv('model_lock_manifest.csv');stopifnot(man$md5==unname(tools::md5sum('locked_development_models.rds')))
ch<-data.frame(check=c('All 30 candidates fitted; all 300 held-out folds passed','Unique PSU assignment in both repeats','2014-2016 and 2017-2019 participants disjoint','Prior seven methods scores exactly unchanged','Same temporal participants across all ten methods','Primary winners match prewritten selection rule','Model MD5 unchanged after temporal evaluation','No temporal CDF clipping'),pass=c(rep(TRUE,7),all(abs(subset(sc,method=='primary')$z)<qnorm(1-1e-12)-1e-8)))
stopifnot(all(ch$pass));wc(ch,'final_integrity_checks')
comp<-list()
for(s in c('Female','Male'))for(m in c('primary','simplified','age_adjusted_allometry','Height2')){
 q<-subset(qq,sex==s&method==m);loss<-mean(((q$p10_pct/100-.1)^2/.09+(q$p5_pct/100-.05)^2/.0475)/2)
 id<-switch(m,primary='BCT_M1_S2',simplified='NO_M1_S1',age_adjusted_allometry='Allometry_ageadj',Height2='Height2')
 comp[[paste(s,m)]]<-data.frame(sex=s,method=m,cv_tail_loss=subset(ds,sex==s&candidate==id)$tail_loss,temporal_tail_loss=loss)
}
wc(do.call(rbind,comp),'tail_loss_development_vs_temporal')
md<-function(x){for(j in seq_along(x))if(is.numeric(x[[j]]))x[[j]]<-formatC(x[[j]],digits=4,format='fg');paste(c(paste('|',paste(names(x),collapse=' | '),'|'),paste('|',paste(rep('---',ncol(x)),collapse=' | '),'|'),apply(x,1,function(a)paste('|',paste(a,collapse=' | '),'|'))),collapse='\n')}
lines<-c('# Phase 2：完整分布筛选、锁模及冻结时间评价','',
'结论：本轮完成，但未满足直接进入“2014–2019最终韩国模型 → 美国外部验证”的整体GO条件。男性获得局部结构收益；女性相对简化Gaussian没有一致尾部增益，年龄调整allometry在女性时间尾部表现更好。保留开发集锁定选择，不根据时间结果重新挑选模型。','',
'## 设计与执行','',
'开发集2014–2016：男性6742、女性8476；时间评价2017–2019：男性7626、女性9467。两次五折PSU分组交叉验证，seed=20260910；每折重估样条节点、参数及allometry指数。两性共用PSU分组，300次训练/留出评分全部有效，另完成30次完整开发拟合。20个分布候选（男9、女11）及10个传统对照。未作临床或结局扩展。','',
'模型选择只用开发数据。先Gaussian/BCCG/BCT × 三个尺度，再对女性暂选族增加两个有限均值结构。主排序使用五个身高组的P5/P10校准损失；预设均值、尺度、线性斜率容差为筛选规则，不是临床等效界限。BIC/GAIC为权重归一化至样本数的伪似然指标。详见phase2_lock.md。','',
'启动时修复Windows中文路径区域设置，以及模型公式环境中权重变量的可见性；修复PSU抽样在单元素层中的R sample()语义边界。这些是结果生成前的实现修复，未改变候选族或筛选规则。所有最终候选重新通过拟合和评分。','',
'## 已冻结的候选','',md(locks),'',
'两性主模型均为BCT：μ=ns(age,4)+ns(height,3)，logσ=ns(age,4)+ns(height,3)，ν、τ均为常数，共18参数；简化模型为Gaussian相同μ、logσ=ns(height,3)，共12参数。正值分布的μ、σ不能直接当作算术均值与绝对标准差；跨族用Z=Φ⁻¹(F(HGS))比较。HGS作为连续变量处理，无额外随机扰动。','',
'女性μ2增加身高自由度及μ3单一年龄×身高交互均未获选。BCPE未触发：BCT两轮OOF超额峰度男0.043/0.051、女0.032/0.026，均未达到预设0.5。所有30个候选均提供完整开发残差及总体/身高五组worm图；总体接近正态仍不代表各组尾部完全一致。','',
'男性BCT尾损失比简化Gaussian低4.87%，女性低10.54%。男性Gaussian损失仅比最小值高5.118%，刚超“105%以内选更简单者”的界线；这一选择对容差边界敏感，不能宣称BCT稳健压倒Gaussian。未进行事后修改阈值或重新选模。','',
'## 年龄调整allometry','',md(read.csv('age_adjusted_allometry_exponents.csv')),'',
'指数来自log(HGS)~ns(age,4)+b·log(height_m)，SE与CI使用复杂抽样设计；其后用HGS/height^b建立年龄参考Z。两轮CV均在训练折估计b。女性b=1.638、男性1.452，明显低于先前未调年龄的约2.128和2.093。','',
'## 冻结时间评价','',
'2017–2019在前阶段已被多次观察，应称prespecified temporal evaluation sample，不能称untouched independent validation。本轮在锁模与MD5记录后一次生成冻结评分；原七种方法评分原样复用。未重新居中或调整尺度。以下CI为复杂抽样设计推断，方法间差异仅作描述性比较，未给出配对差异CI。','',
md(subset(ss,method %in% c('primary','simplified','age_adjusted_allometry','Height2'))[,c('sex','method','mean_z','sd_z','slope_per10cm','slope_ci_low','slope_ci_high','p10_overall_pct','p5_overall_pct')]),'',
'男性主模型斜率0.0233/10cm，95%CI −0.0171至0.0637；五组P10端点差−0.087个百分点，最大最小差仍2.993个百分点（简化模型3.609），Q3仍高，不能只报端点。P5五组差由简化模型0.729增至1.018，SD五组差由0.0542增至0.0675。男性改善主要在P10梯度，并非所有维度同步改善。样条仍显示身高极端的非线性残差，描述图不能当作已经证明独立。','',
'女性主模型斜率−0.0543，较Model 1的−0.0598绝对值略减约0.0054，但95%CI仍为−0.0920至−0.0167。P10五组差3.685，高于Model 1的3.497，也高于年龄调整allometry的2.552；P5五组差2.259，对照分别2.145及1.262。该有限均值改善不足以支持明确增量。','',
'整体左移保留：男性meanZ=−0.347，P10=17.82%、P5=9.84%；女性meanZ=−0.438，P10=19.19%、P5=10.78%。完整分布没有恢复名义10%/5%覆盖，P10整体偏差甚至比简化Gaussian更大。结构梯度与总体校准必须分开解释。','',
'主模型全部身高组（P90/P95列为超过上端阈值的比例，名义分别10%/5%）：','',md(subset(qq,method=='primary')[,c('sex','quintile','n','mean_z','sd_z','p10_pct','p5_pct','p90_pct','p95_pct')]),'',
'主模型低握力分类梯度：','',md(subset(ss,method=='primary')[,c('sex','p10_OR_per10cm','p10_OR_ci_low','p10_OR_ci_high','p5_OR_per10cm','p5_OR_ci_low','p5_OR_ci_high')]),'',
'开发与时间尾损失（越低越好；开发与时间各沿用其预设身高分位界值，跨样本不视作严格同一分层估计）：','',md(do.call(rbind,comp)),'',
'## 按四项GO条件判断','',
'1. 男性mean/scale/tail明显稳定：部分支持线性身高与P10梯度，但仍有非线性、Q3隆起、总体左移，不能判为全面通过。',
'2. 女性较Model 1改善mean或tail：均值斜率仅小幅数值改善，尾部未改善，证据有限。',
'3. 相对height²/年龄调整allometry有明确增量：相对height²的身高梯度改善保留；相对年龄调整allometry男性局部改善、女性未显示一致优势，不通过整体要求。',
'4. 改善能在时间样本保留：内部总体尾校准收益未一致保留；仅部分结构收益保留。','',
'阶段决定：暂不进入六年合并重拟合或美国外部验证。此次结果支持“完整分布能修正开发分布形状，但不能保证时间可迁移校准”，尚不支持“conditional/full-distribution一定优于强传统对照”。阶段锁模文件继续保留，后续不得利用本次时间结果反复改μ/σ/族或再筛选赢家。','',
'## 文件与复核','',
'locked_development_models.rds及model_lock_manifest.csv：冻结模型与MD5；locked_selection.csv：两性主/简化清单。development_candidate_summary.csv、cv_repeat_metrics.csv、cv_height_quintiles.csv：开发候选和两轮OOF完整指标。diagnostics/：30份worm图、对应曲线及完整开发残差。summary_validation.csv、height_quintiles.csv：20个性别×方法结果、100个身高分组，包含双侧尾率及OR。temporal_mean.png、temporal_scale.png、temporal_tails.png：冻结评价三图。tail_loss_development_vs_temporal.csv：尾部损失对照。','',
'final_integrity_checks.csv与development_integrity_checks.csv均通过：样本无重叠、每PSU每轮唯一折、全部评分有限、各方法同一时间样本、旧评分完全不变、锁模MD5不变、主模型无时间CDF截断。曲线为描述性诊断，无独立观察假设的伪置信带；选择表现未使用嵌套CV校正，需承认有限候选筛选的乐观偏差。')
writeLines(enc2utf8(lines),'Phase2结论.md',useBytes=TRUE)
cat('REPORT AND FINAL CHECKS COMPLETE\n')