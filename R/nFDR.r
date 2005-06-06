#########################################
###
### nFDR: Nonparameric Estimate of FDR
###   Based on Bernstein Polynomials
###
#########################################

# Function Bf returns Bernstein Polynomial
# Approximation of f(x) at x. Vector f contains
# sampled values of f(x) at 0, 1/k, ..., (k-1)/k, 1
Bf<-function(f, k, x){
    if(length(f)!=k+1) stop("length of f must be k+1.")
    m<-length(x)
    bf<-NULL
    for(i in 1:m)  bf[i]<-sum(f*dbinom(0:k, k, x[i]))
    return(bf)
}


#######################
##       h_k(r)       #
#######################

hk<-function(R,k){
    hr<-NULL
    if(max(R)>k) stop("R exceed k")
    LR<-length(R)
    for(i in 1:LR){
        hr[i]<-0
        for(j in 0:(k-1)){
            hr[i]<-hr[i]+(mean(dbinom(j,k-1,1-(1:R[i])/k)))^2
        }
    }
    hr
}

F.tilde<-function(tt,x,n)
{
    Ftilde<-NULL
    EDF<-ecdf(x)
    temp<-0:n
    for(i in 1:length(tt))
        Ftilde[i]<-sum(EDF(temp/n)*dbinom(temp, n, tt[i]))
    return(Ftilde)
}

f.tilde<-function(tt,x,n, Smooth = TRUE)
{
    ftilde<-NULL
    temp<-0:(n-1)
    if(Smooth) f<-n*(F.tilde((temp+1)/n,x,n)-F.tilde(temp/n,x,n))
    else{
        EDF<-ecdf(x)
        f<-n*(EDF((temp+1)/n)-EDF(temp/n))
    }
    for(i in 1:length(tt))
        ftilde[i]<-sum(f*dbinom(temp, n-1,  tt[i]))
    return(ftilde)
}

#########################################
##  Function nFDR returns estmates 
##  of pi_0, confidence intervals for
##  FDR in increasing order of p-values
nFDR<-function(x, r0 = floor(length(x)/100), k0 = floor(length(x)/10), K  = 2*floor(length(x)/10), 
    alpha = 0.05, Trial.r = 5, Trial.k = 10, Method = "approx", Smooth = TRUE)
{
    n<-length(x)
    res<-c(r0, k0)
    rk<-res
    if(Trial.r >= 1 && Trial.k >=1){
        cat("It may take a while to search the optimal r and k.\n")
        cat("Please be patient.\n")
        rk<-.C("est_rk", as.double(x), as.integer(n), as.integer(r0), as.integer(k0), as.integer(K), res = as.integer(c(r0, k0)),
            as.integer(Trial.r), as.integer(Trial.k), as.integer(Method == "approx"), as.integer(Smooth == TRUE))$res
        r<-rk[1]; k<-rk[2]
        #cat("rk=",rk,"\n")
    }
    else{ r<-r0; k<-k0  }
    #cat("r*=", r, " k*=", k, "\n")
    tt<-(0:k)/k
    EDF<-ecdf(x)
    edf<-EDF(tt)
    Ftilde<-Bf(edf, k, tt)
    if(Smooth) f<-k*(Ftilde[2:(k+1)]-Ftilde[1:k])
    else f<-k*(edf[2:(k+1)]-edf[1:k])
    ftilde<-Bf(f, k-1, tt)
    PI0<-mean(ftilde[(k-r+1):(k+1)])
    d<-qnorm(1-alpha/2,0,1)*sqrt(k*hk(r,k)*PI0/n)
    cint.pi0<-c(max(0,PI0-d), min(1,PI0+d))
    #CI<-ci4fdr(x, RK[2], RK[1], alpha, Smooth = TRUE)
    p<-sort(x)
    tt<-(0:n)/n
    edf<-EDF(tt)
    if(Smooth) Fhat<-Bf(edf, n, p)
    else  Fhat<-EDF(p)
    Cint.fdr<-rbind(p*max(0,PI0-d)/Fhat, p*min(1,PI0+d)/Fhat) 
    #d<-qnorm(1-alpha/4,0,1)*sqrt(k*hk(r,k)*PI0/n)
    #dF<-qnorm(1-alpha/4,0,1)*sqrt(Fhat*(1-Fhat)/n)
    #F.L<-Fhat-dF
    #F.U<-Fhat+dF
    #cint.fdr<-rbind(p*max(0,PI0-d)/F.U, p*min(1,PI0+d)/F.L) 
    FDR<-p*PI0/Fhat
    FNR<-apply(matrix(1-(1-p)*PI0/(1-Fhat), ncol = 1), 1, function(x) max(0,x))
    qvalue<-sort(cummin(FDR[n:1]))
    return(list(p = p, r = r, k = k, PI0 = PI0, cint.pi0 = cint.pi0, FDR = FDR, FNR = FNR, 
        Cint.fdr = Cint.fdr, qvalue = qvalue))
 }
