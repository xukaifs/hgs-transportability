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
md<-function(x){for(j in seq_along(x))if(is.numeric(x[[j]]))x[[j]]<-formatC(x[[j]],digits=4,format='fg');paste(c(paste('|',paste(names(x),collapse=' | '),'|'),paste('|',paste(rep('---',ncol(x)),collapse=' | '),'|'),apply(x,1,function(a)paste('|',paste(a,collapse=' | '),'|'))),collapse='\n')}
v<-u[,c('sex','method','mean_z','sd_z','slope_per10cm','slope_ci_low','slope_ci_high','p10','p5','tail_loss')];v$p10<-100*v$p10;v$p5<-100*v$p5;names(v)[names(v)=='p10']<-'P10_pct';names(v)[names(v)=='p5']<-'P5_pct'
lines<-c('# 临时结果更新：韩国逐年漂移与冻结美国压力测试','',
'日期：2026-09-10。按当前授权只完成描述性时间诊断与四种冻结策略的美国压力测试。统计结果及基本核验已完成；韩国逐年官方测量协议、设备校准与执行变更记录尚未核实完整，因此本报告为临时解释，不作漂移原因定论。','',
'目前最有依据的论文主线是：更复杂的分布拟合并不保证更好的时间或跨人群可迁移性。本次美国结果进一步支持“没有稳定、全维度的复杂模型优势”，但不能改写为“复杂模型必然更差”或“简单方法已普遍获胜”。没有根据新结果重新选模。','',
'## 一、韩国逐年结果：2018年低谷，2019年回升','',
'同一完整病例规则、男女分别、20–79岁。下表握力单位kg。raw为复杂抽样加权均值；年龄标准化采用20–29至70–79六组，身高标准化采用开发样本五分位，各自固定为2014–2016同性别参考构成。两种标准化分开完成，不是联合调整。','',
md(a[,c('sex','year','n','hgs_mean','ageband_standardized_hgs','heightband_standardized_hgs')]),'',
'女性年龄标准化握力：26.22 → 25.74 → 24.81 → 24.27 → 23.07 → 24.31 kg。男性：43.38 → 42.96 → 42.29 → 41.42 → 40.10 → 41.21 kg。两性2014–2018均下降，2019均回升；下降在2017年前已存在，不能将整个三年窗口差异归于2017单一断点。','',
'年龄/身高标准化后形态仍在，说明本次粗分组所捕获的年龄或身高构成变化不能解释全部差异。分组内残余构成、其他特征、选择进入测量的人群变化和测量执行等仍可能参与；不能据此认定为真实人群握力世代下降。年龄标准化参考比例视作固定，图中95%正态区间不包含参考比例估计误差。','',
'完整年龄、身高、BMI的均值、SD、P10/P25/中位数/P75/P90已导出korea_annual_distributions.csv。年度均值和SE见korea_annual_hgs.csv；所有标准化分组的样本数、当年及参考比例见korea_standardization_groups.csv，所有分组非空。','',
'## 二、QC发现：2017年前后可测量比例改变','',
md(q[,c('sex','year','n_age20_79_design','n_eligible','weighted_missing_hgs_pct','six_valid_pct')]),'',
'分母为20–79岁且抽样设计有效者。女性握力缺失率2016年的8.48%降至2017年的3.04%；男性4.67%降至1.83%。六次有效率分别由88.51%升至94.64%、92.45%升至95.94%。这提示测量完成或进入测量样本的构成存在年份变化，但不能说明新纳入者一定更弱，也不能定量解释左移。','',
'六次左右手原始读数的年度均值和缺失率均已保留。它们可用于定位试次顺序或测量完成的变化；不同试次可测者未必完全相同，不能将这些均值差直接当作同一个人的疲劳/学习效应。既往数据核查显示韩国各年均有六次握力字段，但未找到可与NHANES effort/完成原因完全对应的公开标记。字段稳定不是协议或设备校准稳定的证据。','',
'美国官方2011–2012与2013–2014文档均描述双手各三次交替测试、同手间至少60秒、定期监督和设备校准；本次使用有效试次最大值，而非官方派生的左右手最大值之和。[CDC 2011–2012文档](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/MGX_G.htm)、[CDC 2013–2014文档](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/MGX_H.htm)。','',
'本地韩国数据目录未发现对应六年的测量手册；官方站检索尚未取得足以逐年核实2014–2019设备型号、校准/更换、握距、试次休息、开始手、鼓励语及不能测试规则的完整原始记录。此次公开资料核查不能判定“2017换设备”，也不能判定“六年协议完全相同”。临时报告保留这一缺口，后续优先核实2017测量缺失变化与2018低谷/2019回升的执行背景。','',
'## 三、冻结模型在NHANES的表现','',
'直接使用韩国2014–2016的Height²、年龄调整allometry、简化Gaussian、BCT。未重新估计指数、年龄参考、样条、μ、σ、ν、τ；没有位置或尺度重校准。美国2011–2014合并9334人（女性4718、男性4616），每种方法同一批人；WTMEC2YR/2，按周期组合PSU及分层建立设计后做domain分析。保留先前NHANES有效试次QC规则。美国数据在项目早期已被用于另一个美国内部试跑，因此不称从未接触的样本；本次没有使用美国数据选择韩国模型。','',
'身高分组采用美国合并样本同性别加权五分位，固定后所有方法共用。这是本地身高梯度比较，不是韩国同一绝对身高分组；跨国家的tail loss不能解释为同一人群的配对损失差。斜率CI条件于冻结模型，未包含开发拟合或方法间差异不确定性。','',
md(v),'',
'BCT到美国的总体meanZ为女性+0.739、男性+0.469，与韩国时间评价的左移方向相反。实际P10女性3.14%、男性6.10%，P5女性1.77%、男性3.30%，仍未达到名义覆盖。不能将两个国家的方向差异直接解释为族群生物学差异，人口构成、日历时间和测量执行均可能参与。','',
'男性BCT身高斜率+0.250（95%CI +0.200至+0.301）Z/10cm，P10五组最大差6.32个百分点；Gaussian为+0.231、5.78个百分点，年龄调整allometry为+0.156、5.24个百分点，Height²为−0.053、1.51个百分点。这次男性BCT在韩国时间样本的相对平坦结构并未跨国保留。','',
'同时，美国BCT的tail loss仍低于Gaussian（男0.01530对0.01824，女0.03757对0.04047），且BCT的分组SD差小于Gaussian。这些局部收益必须保留，不能为了新的论文叙事把BCT全部写成失败。Height²本次tail loss点估计最低（男0.01095、女0.03630），但没有进行方法间差异推断或据此改选主模型，不能宣称已证明Height²普遍最佳。','',
'训练支持度：','',md(support),'',
'女性0.37%、男性1.45%的美国加权人群身高超出韩国同性别训练最小/最大值；年龄均在开发范围内。外推者原样保留，无CDF数值截断。该比例仅衡量边界外推，不代表边界内协变量分布完全重叠，也不构成稳定性证明。','',
'## 四、当前可以写到什么程度','',
'可以写：开发内部更好的分布校准，没有转化为稳定的时间和跨人群整体优势；表现排序依赖评价人群、校准维度和简单对照的定义。韩国存在年龄/身高粗标准化后仍保留的年度握力变化，同时可测量比例也随年份变化。','',
'暂不能写：握力真实世代下降已被证实；2017发生了设备更换；BCT永远差于简单方法；allometry在所有人群都更稳健；复杂模型只需重新居中就能解决运输问题。','',
'项目建议：继续围绕transportability开展解释与写作，保持目前模型冻结。眼下未解决的问题是韩国年度测量与样本选择的来源核实，而不是再增加模型。研究标题可暂用“Transportability of Handgrip-Strength Normalization Methods Across Time and Populations”，避免承诺条件分布框架优越。','',
'## 五、交付与复核','',
'先行版临时结论_先行版.md保留Phase 2原始判断；本报告加入逐年与美国结果。korea_annual_hgs.png和nhanes_tail_stress_test.png已检查可读性。nhanes_frozen_summary.csv、nhanes_height_quintiles.csv提供四种方法全部结果；nhanes_frozen_scores.rds保存逐人评分。模型MD5仍与Phase 2原始锁定一致；7项integrity_checks.csv均通过。代码run_diagnostics.R与report_verify.R可复现当前计算。报告中的协议来源缺口仍未完成核实。')
writeLines(enc2utf8(lines),'临时结果_韩国逐年与美国压力测试.md',useBytes=TRUE)
cat('REPORT AND SEVEN CHECKS COMPLETE\n')