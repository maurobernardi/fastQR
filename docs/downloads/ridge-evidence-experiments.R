# Run after sourcing ridge-evidence.R; all stochastic settings are explicit.
evidence_results <- local({
 set.seed(912);X<-matrix(rnorm(100*40),100,40);beta<-rnorm(40,sd=sqrt(4/40));y<-drop(X%*%beta+rnorm(100))
 fit<-ridge_evidence_fit(X,y);p<-fit$prepared
 # Independent dense covariance and augmented-QR checks at several penalties.
 checks<-list();j<-0L
 for(lambda in c(.001,.1,2)){
  alpha<-nrow(X)*lambda;C<-diag(nrow(X))+tcrossprod(X)/alpha;u<-chol(C)
  q<-sum(forwardsolve(t(u),y)^2);ld<-2*sum(log(diag(u)))
  t<-ridge_evidence_terms(p,log(lambda))
  A<-rbind(X,sqrt(alpha)*diag(ncol(X)));f<-fastQR::qr_lapack(A,'thin')
  b<-backsolve(f$R,crossprod(f$Q,c(y,rep(0,ncol(X)))))
  j<-j+1L;checks[[j]]<-data.frame(lambda=lambda,quadratic_error=abs(q-t$q)/max(1,q),
   determinant_error=abs(ld-t$logdet)/max(1,abs(ld)),
   coefficient_error=max(abs(b-ridge_evidence_coefficients(p,lambda))),
   derivative_error=abs(ridge_evidence_score_derivative(p,log(lambda))-(ridge_evidence_terms(p,log(lambda)+1e-5)$score-ridge_evidence_terms(p,log(lambda)-1e-5)$score)/2e-5),
   penalized_identity_error=abs(sum((y-X%*%b)^2)+alpha*sum(b^2)-t$q)/max(1,t$q))
 }
 checks<-do.call(rbind,checks);stopifnot(max(as.matrix(checks[,c("quadratic_error","determinant_error","coefficient_error","penalized_identity_error")]))<1e-10, max(checks$derivative_error)<1e-6)
 # Compare the efficient CV calculation with fastQR's foldwise procedure.
 set.seed(913);xc<-matrix(rnorm(24*32),24,32);yc<-rnorm(24);lg<-c(.05,.3,2)
 native<-fastQR::qrridge_cv(yc,xc,lg,k=3L,seed=19L,type='R')
 spectral<-ridge_cv_spectral(xc,yc,lg,native$foldid)
 cv_error<-max(abs(native$mse-spectral$mse));stopifnot(cv_error<1e-10)
 time_method<-function(fun){invisible(fun());calls<-1L
  repeat{dt<-system.time(for(i in seq_len(calls))z<-fun())[['elapsed']]
   if(dt>=.03||calls>=64L)break;calls<-calls*2L}
  median(replicate(3L,system.time(for(i in seq_len(calls))z<-fun())[['elapsed']]/calls))
 }
 simulations<-list();idx<-0L;grid<-exp(seq(log(1e-4),log(1e2),length.out=31L))
 for(shape in list(c(100L,40L),c(60L,100L)))for(signal in c('Dense Gaussian','Sparse coefficients'))for(rep in 1:5){
  set.seed(10000+100*shape[1]+10*rep+match(signal,c('Dense Gaussian','Sparse coefficients')))
  n<-shape[1];pp<-shape[2];xx<-matrix(rnorm(n*pp),n,pp)
  bb<-if(signal=='Dense Gaussian')rnorm(pp,sd=sqrt(4/pp)) else c(rep(2,3),rep(0,pp-3))
  yy<-drop(xx%*%bb+rnorm(n));xt<-matrix(rnorm(400*pp),400,pp);mu<-drop(xt%*%bb);yt<-mu+rnorm(400)
  ids<-sample(rep(1:5,length.out=n));eb<-ridge_evidence_fit(xx,yy);cv<-ridge_cv_spectral(xx,yy,grid,ids)
  timings<-c(evidence=time_method(function()ridge_evidence_fit(xx,yy)),CV=time_method(function()ridge_cv_spectral(xx,yy,grid,ids)))
  for(method in c('evidence','CV')){
   f<-if(method=='evidence')eb else cv;pred<-drop(xt%*%f$coeff)
   idx<-idx+1L;simulations[[idx]]<-data.frame(n=n,p=pp,signal=signal,repetition=rep,method=method,
    lambda=f$lambda,test_mse=mean((yt-pred)^2),mean_prediction_mse=mean((mu-pred)^2),
    selection_ms=1000*timings[method],boundary=if(method=='evidence')f$status else if(f$boundary)'CV grid boundary' else 'interior')
  }
 }
 list(demo=list(X=X,y=y,truth=beta,fit=fit),checks=checks,cv_error=cv_error,simulations=do.call(rbind,simulations))
})
evidence_summary<-aggregate(cbind(test_mse,mean_prediction_mse,selection_ms)~n+p+signal+method,evidence_results$simulations,median)
evidence_root<-if(dir.exists('../results'))'..' else '.'
write.csv(evidence_results$checks,file.path(evidence_root,'results/ridge-evidence-validation.csv'),row.names=FALSE)
write.csv(data.frame(spectral_native_CV_error=evidence_results$cv_error),file.path(evidence_root,'results/ridge-evidence-cv-validation.csv'),row.names=FALSE)
write.csv(evidence_results$simulations,file.path(evidence_root,'results/ridge-evidence-simulations.csv'),row.names=FALSE)
write.csv(evidence_summary,file.path(evidence_root,'results/ridge-evidence-summary.csv'),row.names=FALSE)
