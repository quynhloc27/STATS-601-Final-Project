library(pracma)
library(MASS)
library(usmap)
library(ggplot2)
library(data.table)
library(dplyr)

covid.states <- read.csv(url("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv"))
covid.counties <- read.csv(url("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv"))
fips.states <- unique(covid.states$fips)
counties <- unique(covid.counties$fips)

dates <- unique(covid.states$date)
no_dates <- nlevels(dates)
data.states <- matrix(0, nrow=no_dates,ncol=length(fips.states))
colnames(data.states) <- fips.states
row.names(data.states) <- dates

for (row in dates){
  data.states[row,c(as.character(covid.states$fips[covid.states$date == row]))] <- covid.states$cases[covid.states$date == row]
}
data.states = data.states[,-c(51,52,53,55,56)]
# Cumulative Sums

fips.states <- colnames(data.states)
dates <- row.names(data.states)
covid = t(data.states)
start_dates = matrix(0,nrow = 1,ncol = dim(covid)[1])
for (ii in c(1:dim(covid)[1])) {
  start_date_ii = 1
  for (tt in c(1:dim(covid)[2])) {
    if(covid[ii,tt] == 0){
      start_date_ii = start_date_ii + 1
    }
  }
  start_dates[,ii] = start_date_ii
}

q = dim(covid)[1]
T = dim(covid)[2]
vars = matrix(0, nrow = q, ncol = 1)
for (i in c(1:q)) {
  if(start_dates[1,i] == T){
    vars[i,] = 0
  } else {
    vars[i,] = var(covid[i,(start_dates[1,i]:T)])
  }
}

covid = covid/sqrt(max(vars))
y = matrix(0, nrow=q,ncol=T)
for (i in c(1:q)) {
  if(start_dates[1,i] == T){
    y[i,start_dates[1,i]] = covid[i,start_dates[1,i]]
  } else {
    y[i,(start_dates[1,i]:T)] = covid[i,(start_dates[1,i]:T)]
  }
}

p = dim(y)[1]
N = dim(y)[2]
# In this model:
# y_i = \Theta\xi(x_i)\eta_i + \epsilon_i
# \eta_i = \psi(x_i) + \xi_i
# \epsilon_i \sim N(0,\Sigma_0),    \xi_i \sim N(0,I).
x = c(1:N)/N

c = 0.1
d = 1
r = 1e-5
K = matrix(0, nrow=N, ncol=N)
for(ii in c(1:N)){
  for(jj in c(1:N)){
    dist_ii_jj = abs(x[ii]-x[jj])
    K[ii,jj] = d*exp(-c*(dist_ii_jj^2))
  }
}
K = K + diag(r,nrow=N)
invK = solve(K)

prior_params = c(K.c_prior = 1, sig.a_sig = 1, sig.b_sig = 0.1, hypers.a_phi = 1.5, hypers.b_phi = 1.5, hypers.a1 = 2, hypers.a2 = 2); 

settings = c(
  L = 8, # truncation level for dictionary of latent GPs
  k = 12, # latent factor dimension
  Niter = 1000) # number of Gibbs iterations to run

# invK, K, inds_y are important function parameters
inds_y = matrix(1, nrow = p, ncol = N)
inds_y[y==0 ] = 0
inds_y = inds_y > 0

BNP_covreg <- function(y, prior_params, settings, invK, K, inds_y){
  p = dim(y)[1] # number of regions
  N = dim(y)[2] # number of observations
  sample_K_flag = settings["sample_K_flag"]
  k = settings["k"]
  L = settings["L"]
  Niter = settings["Niter"]
  
  # Sample hyperparams from prior:
  delta = c(rep(0, times=L))
  delta[1] = rgamma(1, shape=prior_params["hypers.a1"], scale = 1)
  delta[2:L] = rgamma(length(c(2:L)),shape=prior_params["hypers.a2"], scale=1)
  tau = exp(cumsum(log(delta)))
  phi = matrix(rgamma(n=p*L,shape=prior_params["hypers.a_phi"],scale=1),nrow = p, ncol = L)
  
  # Sample theta, eta, and Sigma initially as prior draws:
  theta = matrix(0, nrow=p, ncol=L)
  for(pp in c(1:p)){
    theta[pp,] = t(chol(diag(1/(phi[pp,]*tau))))%*%as.matrix(rnorm(L))
  }
  
  xi = matrix(rnorm(k*N),nrow=k,ncol=N)
  psi = matrix(0, nrow=k, ncol=N)
  eta = psi + xi
  
  # invSig_vec represents the diagonal elements of \Sigma_0^{-1}:
  invSig_vec = rgamma(n=p,shape = prior_params["sig.a_sig"],scale=1)/prior_params["sig.b_sig"]
  
  
  # Sample zeta_i using initialization scheme based on data and other
  # sampled params:
  zeta = array(rep(0,times = L*k*N), c(L,k,N))
  zeta = sample_zeta(y,theta,eta,invSig_vec,zeta,invK,inds_y)
  num_iters = 1
  nstart = 1
  
  washtenaw_list = list()
  snohomish_list = list()
  lincoln_list = list()
  for (nn in c(nstart:Niter)) {
    # Sample latent factor model additive Gaussian noise covariance
    # (diagonal matrix) given the data, weightings matrix, latent GP
    # functions, and latent factors:
    invSig_vec = sample_sig(y,theta,eta,zeta,prior_params,inds_y)
    #print("invSig_vec")
    #print(invSig_vec)
    
    # Sample weightings matrix hyperparameters given weightings matrix:
    phi_tau = sample_hypers(theta,phi,tau,prior_params)
    phi <- phi_tau$phi
    tau <- phi_tau$tau
    #print("phi")
    #print(colSums(phi))
    #print("tau: ")
    #print(tau)
    
    # Sample weightings matrix given the data, latent factors, latent GP
    # functions, noise parameters, and hyperparameters:
    theta = sample_theta(y,eta,invSig_vec,zeta,phi,tau,inds_y)
    #print("theta")
    #print(colSums(theta))
    
    # sample the latent mean GPs \psi(x) marginalizing \xi
    psi = sample_psi_margxi(y,theta,invSig_vec,zeta,psi,invK,inds_y,nn)
    #psi = matrix(0, nrow=k, ncol=N)
    #print("psi")
    #print(colSums(psi))
    
    # Sample latent factors given the data, weightings matrix, latent GP
    # functions, and noise parameters.
    xi = sample_xi(y,theta,invSig_vec,zeta,psi,inds_y)
    #print("xi")
    #print(colSums(xi))
    
    eta = psi + xi
    
    # Sample latent GP functions zeta_i given the data, weightings matrix,
    # latent factors, noise params, and GP cov matrix (hyperparameter):
    zeta = sample_zeta(y,theta,eta,invSig_vec,zeta,invK,inds_y)
    #print("zeta")
    #print(rowSums(colSums(zeta)))
    
    print(nn)
    
    
    
    if(mod(nn,1000)==0){
      cov_est = matrix(0, nrow = p, ncol = N)
      michigan_list = list()
      washington_list = list()
      nebraska_list = list()
      for (tt in c(1:N)) {
        SigmaX = theta %*% zeta[, ,tt] %*% t(zeta[, ,tt]) %*% t(theta) + diag((1/invSig_vec))
        cov_est[,tt] = diag(SigmaX)
        mat = cov2cor(SigmaX)
        colnames(mat) <- fips.states
        # Michigan
        michigan = as.matrix(mat[,"26"])
        row.names(michigan) <- fips.states
        michigan <- tibble::rownames_to_column(data.frame(michigan),"fips")
        colnames(michigan) <- c("fips","val")
        p1 <- plot_usmap(regions="states",data=michigan,values="val",labels=TRUE) + labs(title="Michigan (1000 Gibbs Iterations)", subtitle=dates[tt]) + scale_fill_continuous(name="time-varying spatial correlation", low ="white", high = "blue", label = scales::comma) + theme(legend.position = "right")
        michigan_list[[tt]] = p1
        
        # Washington
        washington = as.matrix(mat[,"53"])
        row.names(washington) <- fips.states
        washington<-tibble::rownames_to_column(data.frame(washington),"fips")
        colnames(washington) <- c("fips","val")
        p2 <- plot_usmap(regions="states",data=washington,values="val",labels=TRUE) + labs(title="Washington (1000 Gibbs Iterations)", subtitle=dates[tt]) + scale_fill_continuous(name="time-varying spatial correlation",low ="white", high = "red", label = scales::comma) + theme(legend.position = "right")
        washington_list[[tt]] = p2
        
        # Nebraska
        nebraska = as.matrix(mat[,"31"])
        row.names(nebraska) <- fips.states
        nebraska<-tibble::rownames_to_column(data.frame(nebraska),"fips")
        colnames(nebraska) <- c("fips","val")
        p3 <- plot_usmap(regions="states",data=nebraska,values="val",labels=TRUE) + labs(title="Nebraska (1000 Gibbs Iterations)", subtitle=dates[tt]) + scale_fill_continuous(name="time-varying spatial correlation",low ="white", high = "green", label = scales::comma) + theme(legend.position = "right")
        nebraska_list[[tt]] = p3
      }
      
      pdf("michigan.pdf")
      for (i in c(1:N)) {
        print(michigan_list[[i]])
      }
      dev.off()
      pdf("washington.pdf")
      for (i in c(1:N)) {
        print(washington_list[[i]])
      }
      dev.off()
      pdf("nebraska.pdf")
      for (i in c(1:N)) {
        print(nebraska_list[[i]])
      }
      dev.off()
      return()
    }
  }
}

sample_zeta <- function(y, theta, eta, invSig_vec,zeta,invK,inds_y){
  
  p = dim(theta)[1]
  L = dim(theta)[2]
  k = dim(eta)[1]
  N = dim(eta)[2]
  
  mu_tot = matrix(0, nrow=p,ncol=N)
  error_nn = matrix(0, nrow=L,ncol=N)
  for (nn in c(1:N)) {
    mu_tot[,nn] = theta%*%zeta[, ,nn]%*%eta[,nn]
    error_nn[,nn] = t(theta^2) %*% (as.matrix(invSig_vec) * as.matrix(1-inds_y[,nn]))
  }
  for (ll in c(1:L)) {
    theta_ll = theta[,ll]
    for (kk in randperm(k,k)) {
      eta_kk = eta[kk,]
      zeta_ll_kk = t(drop(zeta[ll,kk,]))
      mu_tot = mu_tot - matrix(theta_ll,nrow=p,ncol=N,byrow = FALSE) * matrix(eta_kk,nrow=p,ncol=N,byrow = TRUE) * matrix(zeta_ll_kk,nrow=p,ncol=N,byrow = TRUE)
      
      # Using standard Gaussian identities, form posterior of
      # zeta_{ll,kk} using information form of Gaussian prior and likelihood:
      A_lk_invSig_A_lk = (eta[kk,]^2)*(t(theta[,ll]^2)%*%as.matrix(invSig_vec) - error_nn[ll,])
      theta_tmp = t(theta[,ll])*t(invSig_vec)
      ytilde = (y - mu_tot)*inds_y # normalize data by subtracting mean of tilde(eps)
      theta_lk = as.matrix(eta[kk,])*t(theta_tmp%*%ytilde)
      # Transform information parameters
      cholSig_lk_trans = mldivide(chol(invK + diag(A_lk_invSig_A_lk,nrow = N)),diag(nrow=N))
      m_lk = cholSig_lk_trans%*%(t(cholSig_lk_trans)%*%theta_lk)
      
      # Sample zeta_{ll,kk} from posterior Gaussian
      zeta[ll,kk,] = m_lk + cholSig_lk_trans%*%as.matrix(rnorm(N,1))
      zeta_ll_kk = t(drop(zeta[ll,kk,]))
      mu_tot = mu_tot + matrix(theta_ll,nrow=p,ncol=N,byrow=FALSE) * matrix(eta_kk,nrow=p,ncol=N,byrow=TRUE) * matrix(zeta_ll_kk,nrow=p,ncol=N,byrow=TRUE)
    }
  }
  return(zeta)  
}

sample_sig <- function(y,theta, eta, zeta,prior_params,inds_y){
  p = dim(y)[1]
  N = dim(y)[2]
  a_sig = prior_params["sig.a_sig"]
  b_sig = prior_params["sig.b_sig"]
  
  inds_vec = seq(1,N)
  
  invSig_vec = seq(1,p)
  for (pp in seq(1,p)) {
    sq_err = 0
    for (nn in inds_vec[inds_y[pp,]]) {
      sq_err = sq_err + (y[pp,nn] - theta[pp,]%*%zeta[ , ,nn]%*%eta[,nn])^2
    }
    invSig_vec[pp] = rgamma(n=1, shape = (a_sig + 0.5*rowSums(inds_y)[pp]),scale=1)/(b_sig + 0.5*sq_err)
  }
  return(invSig_vec)
}

sample_hypers <- function(theta, phi, tau, prior_params){
  p = dim(theta)[1]
  L = dim(theta)[2]
  
  a1 = prior_params["hypers.a1"]
  a2 = prior_params["hypers.a2"]
  a_phi = prior_params["hypers.a_phi"]
  b_phi = prior_params["hypers.b_phi"]
  
  a = c(a1,rep(1,times=L-1))
  delta = c(exp(log(tau[1])),exp(diff(log(tau))))
  for (numIter in c(1:50)) {
    phi = matrix(rgamma(n=p*L,shape=(a_phi + 0.5),scale=1),nrow=p,ncol=L) / (b_phi + 0.5*matrix(tau,nrow=p,ncol=L,byrow = TRUE)*(theta^2))
    sum_phi_theta = colSums(phi*(theta^2))
    for (hh in c(1:L)) {
      tau_hh = exp(cumsum(log(delta)))*c(rep(0,times=hh-1),rep(1,times=L-hh+1)/delta[hh])
      delta[hh] = rgamma(n=1,shape=a[hh] + 0.5*p*(L-hh+1),scale = 1)/(1+0.5*sum(tau_hh*sum_phi_theta))
    }
    tau = exp(cumsum(log(delta)))
  }
  phi_tau <- list("phi"=phi, "tau"=tau)
  return(phi_tau)
}

sample_theta <- function(y,eta,invSig_vec,zeta,phi,tau,inds_y){
  p = dim(y)[1]
  N = dim(y)[2]
  L = dim(zeta)[1]
  theta = matrix(0, nrow = p, ncol = L)
  
  eta_tilde = matrix(0, nrow = L, ncol = N)
  for (nn in c(1:N)) {
    eta_tilde[,nn] = zeta[, ,nn]%*%eta[,nn]
  }
  eta_tilde = t(eta_tilde)
  
  for (pp in c(1:p)) {
    inds_y_p = t(inds_y[pp,])
    eta_tilde_p = eta_tilde*matrix(inds_y_p, nrow=N, ncol=L,byrow = FALSE)
    chol_Sig_theta_p_trans = mldivide(chol(diag(phi[pp,]*tau,nrow=L)+ invSig_vec[pp]*(t(eta_tilde_p)%*%eta_tilde_p)),diag(nrow=L))
    m_theta_p = invSig_vec[pp]*(t(chol_Sig_theta_p_trans)%*%chol_Sig_theta_p_trans)%*%(t(eta_tilde_p)%*%y[pp,])
    theta[pp,] = m_theta_p + chol_Sig_theta_p_trans%*%as.matrix(rnorm(L,1))
  }
  return(theta)
} 

sample_xi <- function(y,theta,invSig_vec,zeta,psi,inds_y){
  # Sample latent factors eta_i using standard Gaussian identities based on
  # the fact that:
  
  # y_i = (theta*zeta)*eta_i + eps_i,   eps_i \sim N(0,Sigma_0),   
  # eta_i \sim N(0,I) 
  
  # and using the information form of the Gaussian likelihood and prior.
  p = dim(y)[1]
  N = dim(y)[2]
  L = dim(zeta)[1]
  k = dim(zeta)[2]
  
  xi = matrix(0, nrow=k, ncol=N)
  for (nn in c(1:N)) {
    theta_zeta_n = theta%*%zeta[, ,nn]
    y_tilde_n = y[,nn]-theta_zeta_n%*%psi[,nn]
    invSigMat = invSig_vec*t(inds_y[,nn])
    invSigMat = matrix(invSigMat, nrow=k,ncol=p,byrow=TRUE)
    zeta_theta_invSig = t(theta_zeta_n)*invSigMat
    cholSig_xi_n_trans = mldivide(chol(diag(nrow=k) + zeta_theta_invSig%*%theta_zeta_n),diag(nrow=k))
    m_xi_n = cholSig_xi_n_trans%*%(t(cholSig_xi_n_trans)%*%(zeta_theta_invSig%*%y_tilde_n))
    xi[,nn] = m_xi_n + cholSig_xi_n_trans%*%as.matrix(rnorm(k,1))
  }
  return(xi)
}

sample_psi_margxi <- function(y,theta,invSig_vec,zeta,psi,invK,inds_y,iter){
  p = dim(theta)[1]
  L = dim(theta)[2]
  k = dim(psi)[1]
  N = dim(psi)[2]
  Sigma_0 = diag(1/invSig_vec)
  
  # We derive the sequential sampling of the zeta_{ll,kk} by reformulating
  # the full regression problem as one that is solely in terms of
  # zeta_{ll,kk} having conditioned on the other latent GP functions:
  
  # y_i = eta(i,m)*theta(:,ll)*zeta_{ll,kk}(x_i) + tilde(eps)_i,
  
  # Initialize the structure for holding the conditional mean for additive
  # Gaussian noise term tilde(eps)_i and add values based on previous zeta:
  mu_tot = matrix(0,nrow = p, ncol = N)
  Omega = array(rep(0,times = p*k*N), c(p,k,N))
  OmegaInvOmegaOmegaSigma0 = array(rep(0,times = k*p*N), c(k,p,N))
  for (nn in 1:N) {
    Omega[inds_y[,nn], ,nn] = theta[inds_y[,nn],]%*%zeta[, ,nn]
    if(dim(as.matrix(Omega[inds_y[,nn], ,nn]))[2]==1){
      temp = mldivide((t(Omega[inds_y[,nn], ,nn]) %*% as.matrix(Omega[inds_y[,nn], ,nn]) + Sigma_0[inds_y[,nn],inds_y[,nn]]) , diag(sum(inds_y[,nn])))
      OmegaInvOmegaOmegaSigma0[,inds_y[,nn],nn] = as.matrix(Omega[inds_y[,nn], ,nn]) %*% temp
    } else {
      temp = mldivide((as.matrix(Omega[inds_y[,nn], ,nn]) %*% t(Omega[inds_y[,nn], ,nn]) + Sigma_0[inds_y[,nn],inds_y[,nn]]), diag(sum(inds_y[,nn])))
      OmegaInvOmegaOmegaSigma0[,inds_y[,nn],nn] = t(Omega[inds_y[,nn], ,nn])%*%temp
    }
    mu_tot[,nn] = Omega[, ,nn]%*%psi[,nn]
  }
  if (sum(sum(mu_tot))==0){ # if this is a call to initialize psi
    numToIters = 50
  } else {
    numToIters = 5
  }
  for (numIter in c(1:numToIters)) {
    for (kk in randperm(k,k)) { # create random ordering for kk in sampling zeta_{ll,kk} 
      Omega_kk = drop(Omega[,kk,])
      psi_kk = psi[kk,]
      mu_tot = mu_tot - Omega_kk*matrix(psi_kk,nrow=p,ncol=N,byrow=TRUE)
      
      theta_k = as.matrix(diag(t(drop(OmegaInvOmegaOmegaSigma0[kk, ,]))%*%(y-mu_tot)))
      Ak_invSig_Ak = diag(t(drop(OmegaInvOmegaOmegaSigma0[kk, ,]))%*%Omega_kk)
      cholSig_k_trans = mldivide(chol(invK + diag(Ak_invSig_Ak)),diag(nrow=N))
      
      # Transform information parameters
      m_k = cholSig_k_trans%*%(t(cholSig_k_trans)%*%theta_k)
      
      # Sample zeta_{ll,kk} from posterior Gaussian
      psi[kk,] = m_k+cholSig_k_trans%*%matrix(rnorm(N),nrow=N,ncol=1)
      
      psi_kk = psi[kk,]
      mu_tot = mu_tot + Omega_kk*matrix(psi_kk,nrow=p,ncol=N,byrow=TRUE)
    }
  }
  return(psi)
}

BNP_covreg(y, prior_params, settings, invK, K, inds_y)
