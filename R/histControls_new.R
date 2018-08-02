hct<-function(data,estimate,standardError,N,iter=2000,rseed=NA,
              silent=TRUE,constantStderr=TRUE) {

if(constantStderr){
 stanmodelcode <- "
  data {int<lower=0> N;real y[N];real sig[N];}
  parameters {real mu;real u[N];real<lower=0> vsig;}
  model {
  mu ~ normal(0, 100);
  vsig~uniform(0,100);
  u ~ normal(mu, vsig);
  y~normal(u,sig);#observed y, #observed sig(STDERR Y)
  }
  generated quantities{real y_pred;
  y_pred = normal_rng(mu,vsig);
  }"
 dat=list(N=dim(data)[1],y=data[,estimate],sig=data[,standardError])
 }
else{
stanmodelcode <- "
  data {int<lower=0> N;real y[N];real sig[N];real ssize[N];}
 transformed data {
 real<lower=0> sig2[N];
 for(i in 1:N){
 sig2[i]=sig[i]*sig[i];}
 }
 parameters {real mu;real u[N];real<lower=0> esigc[N];real mus;real<lower=0> sigmas;real<lower=0> vsig;}
 model {
 for(i in 1:N){
 esigc[i]~lognormal(mus,sigmas);
 u[i] ~ normal(mu, vsig);
 sig2[i]~gamma((ssize[i]-1)/2,ssize[i]*(ssize[i]-1)/(2*esigc[i]));
 y[i]~normal(u[i],sqrt(esigc[i]/ssize[i]));//observed y, #observed sig(STDERR Y)
 }
 }
 generated quantities{real y_pred;
  y_pred = normal_rng(mu,vsig);
}"
dat <- list(N = dim(data)[1], y = data[,estimate],sig=data[,standardError],ssize=data[,N])
}

if(silent) sink(file='messages.txt')
                fit <- rstan::stan(model_code =stanmodelcode,
                model_name = "Prior treatment effect",
                data =dat,
                iter = iter,verbose = FALSE,seed=rseed,
                control=list(adapt_delta=.998))
  sampleStudy=extract(fit,c("y_pred"),permuted=FALSE, inc_warmup=FALSE)
  #Samples from posterior of mu and sigma.
  if(silent) {file.remove('messages.txt');sink()}
  VAR1=sum(data[,standardError]^2*data[,N]*(data[,N]-1))/(sum(data[,N])-1)
  power=function(t,meanValue,se,df=NULL){
    if(is.null(df)) mean(pnorm(t,meanValue+sampleStudy,se,lower.tail=FALSE))
    else mean(pt((t-sampleStudy)/se,df,lower.tail=FALSE))
    }
  criteria=function(p,se,df=NULL){
  lowValue=qnorm(1-p,min(sampleStudy),se)
  highValue=qnorm(1-p,max(sampleStudy),se)
  return(uniroot(function(t) power(t,0,se,df)-p,c(lowValue,highValue))$root)
  }
  outp=list(criteria=criteria,power=power,effective.SD=sqrt(VAR1),fit=fit)
  class(outp)<-append(class(outp),'hct')
  return(outp)
  }
summary.hct<-function(object,...) list(effective.SD=object[[3]],
                                       prior.distribution=summary(object[[4]]))
print.hct<-function(x,...) print(summary(x,...))

