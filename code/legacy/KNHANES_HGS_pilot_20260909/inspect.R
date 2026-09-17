invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
.libPaths(c('C:/Users/Public/CodexRLib42Copy',.libPaths()))
library(haven)
root<-'新一篇/KNHANES_HGS_pilot_20260909'
for(y in 14:19){
 x<-read_sas(file.path(root,'raw',paste0('hn',y,'_all.sas7bdat')))
 labs<-vapply(x,function(v) {l<-attr(v,'label');if(is.null(l)) '' else l},character(1))
 cat('\nALL GRIP/HAND LABEL MATCHES\n'); print(labs[grepl('악력|손목|오른손|왼손|grip|effort',labs,ignore.case=TRUE)])
 nm<-grep('^ID$|^age$|^sex$|^year$|^HE_(ht|wt|BMI|wc)$|^GS|^wt_|^psu$|^kstrata$',names(x),value=TRUE,ignore.case=TRUE)
 cat('\nYEAR',y,'N',nrow(x),'\n')
 for(n in nm){v<-x[[n]];cat(n,'LABEL:',attr(v,'label'),' RANGE:',paste(range(v,na.rm=TRUE),collapse=' / '),' NA:',sum(is.na(v)),'\n');if(grepl('^GS',n))print(head(sort(table(v),decreasing=TRUE),12))}
}
z<-read_sav(file.path(root,'raw','HN14_all.sav'),user_na=TRUE)
for(n in c('sex','GS_use',grep('^GS_mea',names(z),value=TRUE))){cat('\nSPSS ATTR',n,'\n');print(attributes(z[[n]]))}
