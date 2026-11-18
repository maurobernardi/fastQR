# Guide-level empirical Bayes helpers; no changes to fastQR's public API.
ridge_evidence_prepare <- function(X,y) {
 if(!is.matrix(X)||!is.numeric(X)||!all(is.finite(X))||!nrow(X)||!ncol(X))stop('X must be a nonempty finite numeric matrix.')
 if(!is.numeric(y)||length(y)!=nrow(X)||any(!is.finite(y))||sum(y^2)<=0)stop('y must be finite, nonzero, and match X.')
 s<-svd(X,nu=min(dim(X)),nv=min(dim(X)))
 z<-drop(crossprod(s$u,y));perp<-y-drop(s$u%*%z)
 list(n=nrow(X),p=ncol(X),d=s$d,z=z,v=s$v,perp2=sum(perp^2),y2=sum(y^2))
}
ridge_evidence_terms <- function(prepared,log_lambda) {
 a<-prepared$n*exp(log_lambda)
 q<-prepared$perp2+sum(prepared$z^2*(a/(prepared$d^2+a)))
 logdet<-sum(log1p(prepared$d^2/a))
 list(score=.5*(prepared$n*log(q/prepared$n)+logdet),
      q=q,logdet=logdet,sigma2=q/prepared$n)
}
ridge_evidence_score_derivative <- function(prepared, log_lambda) {
 a <- prepared$n * exp(log_lambda)
 q <- ridge_evidence_terms(prepared, log_lambda)$q
 dq <- sum(prepared$z^2 * a * prepared$d^2 / (prepared$d^2 + a)^2)
 df <- sum(prepared$d^2 / (prepared$d^2 + a))
 .5 * (prepared$n * dq / q - df)
}
ridge_evidence_coefficients <- function(prepared,lambda) {
 if(is.infinite(lambda))return(rep(0,prepared$p))
 drop(prepared$v%*%(prepared$d*prepared$z/(prepared$d^2+prepared$n*lambda)))
}
ridge_evidence_fit <- function(X,y,log_bounds=c(-12,12),grid_length=121L) {
 if(length(log_bounds)!=2L||any(!is.finite(log_bounds))||log_bounds[1]>=log_bounds[2]||max(abs(log_bounds))>600)stop('Supply finite, increasing log-lambda bounds of moderate magnitude.')
 if(length(grid_length)!=1L||!is.finite(grid_length)||grid_length<3L||grid_length!=as.integer(grid_length))stop("grid_length must be an integer of at least three.")
 p<-ridge_evidence_prepare(X,y)
 objective<-function(t)ridge_evidence_terms(p,t)$score
 grid<-seq(log_bounds[1],log_bounds[2],length.out=grid_length)
 values<-vapply(grid,objective,numeric(1));candidate<-c(1L,grid_length)
 local_min<-which(values[2:(grid_length-1)]<=values[1:(grid_length-2)] & values[2:(grid_length-1)]<=values[3:grid_length])+1L
 candidates<-data.frame(log_lambda=grid[candidate],score=values[candidate])
 for(i in local_min){o<-optimize(objective,c(grid[i-1],grid[i+1]),tol=1e-8)
  candidates<-rbind(candidates,data.frame(log_lambda=o$minimum,score=o$objective))}
 best<-candidates[which.min(candidates$score),];null_score<-.5*p$n*log(p$y2/p$n)
 if(null_score<=best$score){lambda<-Inf;status<-'null-model limit';sigma2<-p$y2/p$n;score<-null_score}
 else {lambda<-exp(best$log_lambda);terms<-ridge_evidence_terms(p,best$log_lambda);sigma2<-terms$sigma2;score<-best$score
  status<-if(abs(best$log_lambda-log_bounds[1])<1e-6)'lower search boundary' else if(abs(best$log_lambda-log_bounds[2])<1e-6)'upper search boundary' else 'interior'}
 list(lambda=lambda,alpha=p$n*lambda,sigma2=sigma2,coeff=ridge_evidence_coefficients(p,lambda),
  score=score,status=status,prepared=p,curve=data.frame(lambda=exp(grid),score=values))
}
ridge_cv_spectral <- function(X,y,lambda_grid,foldid) {
 if(any(!is.finite(lambda_grid))||any(lambda_grid<=0))stop('Use a finite positive CV grid.')
 mse<-numeric(length(lambda_grid))
 for(f in sort(unique(foldid))){tr<-which(foldid!=f);te<-which(foldid==f)
  p<-ridge_evidence_prepare(X[tr,,drop=FALSE],y[tr])
  pred<-vapply(lambda_grid,function(lambda)drop(X[te,,drop=FALSE]%*%ridge_evidence_coefficients(p,lambda)),numeric(length(te)))
  mse<-mse+colMeans((pred-y[te])^2)
 }
 mse<-mse/length(unique(foldid));selected<-lambda_grid[which.min(mse)]
 p<-ridge_evidence_prepare(X,y)
 list(lambda=selected,coeff=ridge_evidence_coefficients(p,selected),mse=mse,lambda_grid=lambda_grid,
      boundary=which.min(mse)%in%c(1L,length(mse)))
}
