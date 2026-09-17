options(stringsAsFactors=FALSE,warn=1)
invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
library(splines)
script_path<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
root<-dirname(normalizePath(script_path,winslash='/'));old<-file.path(dirname(root),'KNHANES_HGS_pilot_20260909')
d0<-readRDS(file.path(old,'development_data.rds'));base<-readRDS(file.path(old,'frozen_Korea_2014_2016_models.rds'))
bundles<-list();diagnostics<-list()
for(s in c('Female','Male')){
 d<-subset(d0,eligible & sex==s);b0<-base[[paste0(s,'_Conditional')]]
 termsX<-delete.response(terms(b0$model));X<-model.matrix(termsX,data=d)
 sigma_idx<-c(1,which(grepl('ns(height',colnames(X),fixed=TRUE)));S<-X[,sigma_idx,drop=FALSE]
 stopifnot(ncol(X)==8,ncol(S)==4)
 y<-d$hgs;w<-d$weight/mean(d$weight);p<-ncol(X);q<-ncol(S)
 nll<-function(theta){mu<-as.numeric(X%*%theta[seq_len(p)]);ls<-as.numeric(S%*%theta[p+seq_len(q)]);sum(w*(ls+.5*(y-mu)^2*exp(-2*ls)))}
 grad<-function(theta){mu<-as.numeric(X%*%theta[seq_len(p)]);ls<-as.numeric(S%*%theta[p+seq_len(q)]);r<-y-mu;iv<-exp(-2*ls);c(-crossprod(X,w*r*iv),crossprod(S,w*(1-r^2*iv)))}
 start<-c(coef(b0$model),log(b0$sigma),rep(0,q-1))
 fit<-optim(start,nll,grad,method='BFGS',control=list(maxit=2000,reltol=1e-12,parscale=c(rep(10,p),rep(1,q))))
 normgrad<-max(abs(grad(fit$par)))/sum(w);stopifnot(fit$convergence==0,normgrad<1e-5,fit$value<=nll(start)+1e-6)
 mu<-as.numeric(X%*%fit$par[seq_len(p)]);sigma<-exp(as.numeric(S%*%fit$par[p+seq_len(q)]));stopifnot(all(is.finite(sigma)),all(sigma>0))
 probe<-start+seq_along(start)*1e-5;eps<-1e-5
 ng<-vapply(seq_along(probe),function(j){a<-b<-probe;a[j]<-a[j]+eps;b[j]<-b[j]-eps;(nll(a)-nll(b))/(2*eps)},numeric(1))
 ge<-max(abs(ng-grad(probe)))/max(1,max(abs(ng)));stopifnot(ge<1e-5)
 H<-optimHess(fit$par,nll,grad);mineig<-min(eigen(H,symmetric=TRUE,only.values=TRUE)$values);stopifnot(mineig>0)
 # Zero-pad age terms explicitly so scoring uses the same frozen 8-column basis.
 gamma<-setNames(rep(0,p),colnames(X));gamma[sigma_idx]<-fit$par[p+seq_len(q)]
 bundles[[s]]<-list(sex=s,method='Conditional_heightSD',family='Gaussian',mu_coefficients=setNames(fit$par[seq_len(p)],colnames(X)),log_sigma_coefficients=gamma,sigma_column_indices=sigma_idx,terms=termsX,development_uids=d$uid,development_years=2014:2016,n_development=nrow(d),convergence=fit$convergence,mean_gradient=normgrad)
 diagnostics[[s]]<-data.frame(sex=s,n=nrow(d),parameters=length(start),convergence=fit$convergence,max_gradient_per_weight=normgrad,finite_difference_gradient_relative_error=ge,min_hessian_eigenvalue=mineig,nll_constant_sigma=nll(start),nll_height_sigma=fit$value,development_sigma_min=min(sigma),development_sigma_max=max(sigma))
 cat(s,'converged; mean gradient',normgrad,'\n')
}
saveRDS(bundles,file.path(root,'frozen_Gaussian_location_scale_models.rds'))
write.csv(do.call(rbind,diagnostics),file.path(root,'fit_diagnostics.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(root,'fit_sessionInfo.txt'))
