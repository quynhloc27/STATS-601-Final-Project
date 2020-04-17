library(pracma)

covid.states <- read.csv(url("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv"))
covid.counties <- read.csv(url("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv"))
states <- unique(covid.states$state)
counties <- unique(covid.counties$county)

dates <- unique(covid.states$date)
no_dates <- nlevels(dates)
data.states <- matrix(0, nrow=no_dates,ncol=nlevels(states))
colnames(data.states) <- states
row.names(data.states) <- unique(covid.states$date)

for (row in dates){
  data.states[row,c(as.character(covid.states$state[covid.states$date == row]))] <- covid.states$cases[covid.states$date == row]
}
data.states = data.states[,-56]
# Cumulative Sums

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
  vars[i,]=var(covid[i,start_dates[1,i]:T])
}

covid = covid/sqrt(max(vars))
y = matrix(0, nrow=q,ncol=T)
for (i in c(1:q)) {
  y[i,start_dates[1,i]:T] = covid[i,start_dates[1,i]:T]
}
y = 1.75*y

p = dim(y)[1]
N = dim(y)[2]
# In this model:
# y_i = \Theta\xi(x_i)\eta_i + \epsilon_i
# \eta_i = \psi(x_i) + \xi_i
# \epsilon_i \sim N(0,\Sigma_0),    \xi_i \sim N(0,I).
x = c(1:N)/N

c = 100  # Gaussian process length scale parameter.  Smaller c will lead to slower changes in correlations.
d = 1
r = 1e-5
K = matrix(0, nrow=N, ncol=N)
for(ii in c(1:N)){
  for(jj in c(1:N)){
    dist_ii_jj = abs(x[ii]-x[jj])
    K[ii,jj] = d*exp(-c*(dist_ii_jj^2))
  }
}
K = K + diag(r*rep(1,times=N),nrow = N, ncol = N)
invK = solve(K)
logdetK = 2*sum(log(diag(chol(K))))

prior_params = c(K.c_prior = 1, K.logdetK = logdetK, sig.a_sig = 1, sig.b_sig = 0.1, hypers.a_phi = 1.5, hypers.b_phi = 1.5, hypers.a1 = 2, hypers.a2 = 2); 

settings = c(
  L = 5, # truncation level for dictionary of latent GPs
  k = 4, # latent factor dimension
  Niter = 2) # number of Gibbs iterations to run

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
    theta[pp,] = t(chol(diag(1/(phi[pp,]*tau)), pivot = TRUE))%*%matrix(rnorm(L),nrow=L,ncol=1)
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
  for (nn in c(nstart:Niter)) {
    # Sample latent factor model additive Gaussian noise covariance
    # (diagonal matrix) given the data, weightings matrix, latent GP
    # functions, and latent factors:
    invSig_vec = sample_sig(y,theta,eta,zeta,prior_params,inds_y)
    
    # Sample weightings matrix hyperparameters given weightings matrix:
    phi_tau = sample_hypers(theta,phi,tau,prior_params)
    phi <- phi_tau$phi
    tau <- phi_tau$tau
    
    # Sample weightings matrix given the data, latent factors, latent GP
    # functions, noise parameters, and hyperparameters:
    theta = sample_theta(y,eta,invSig_vec,zeta,phi,tau,inds_y)
    
    # sample the latent mean GPs \psi(x) marginalizing \xi
    #psi = sample_psi_margxi(y,theta,invSig_vec,zeta,psi,invK,inds_y,nn)
    psi = matrix(0, nrow=k, ncol=N)
    
    # Sample latent factors given the data, weightings matrix, latent GP
    # functions, and noise parameters.
    xi = sample_xi(y,theta,invSig_vec,zeta,psi,inds_y)
    
    eta = psi + xi
    
    # Sample latent GP functions zeta_i given the data, weightings matrix,
    # latent factors, noise params, and GP cov matrix (hyperparameter):
    zeta = sample_zeta(y,theta,eta,invSig_vec,zeta,invK,inds_y)
    
    cov_est = matrix(0, nrow = p, ncol = N)
    print(nn)
    for (tt in c(1:N)) {
      cov_est[,tt] = diag(theta %*% zeta[, ,tt] %*% t(zeta[, ,tt]) %*% t(theta) + diag(1/invSig_vec))
    }
    print(cov_est)
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
      cholSig_lk_trans = chol2inv(invK + diag(A_lk_invSig_A_lk,nrow = N))
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
    eta_tilde[,nn] = zeta[, ,nn]%*%as.matrix(eta[,1])
  }
  eta_tilde = t(eta_tilde)
  
  for (pp in c(1:p)) {
    inds_y_p = t(inds_y[pp,])
    eta_tilde_p = eta_tilde*matrix(inds_y_p, nrow=N, ncol=L)
    chol_Sig_theta_p_trans = chol2inv(diag(phi[pp,]*tau) + invSig_vec[pp]*(t(eta_tilde_p)%*%eta_tilde_p))
    
    m_theta_p = invSig_vec[pp]*(t(chol_Sig_theta_p_trans)%*%chol_Sig_theta_p_trans)%*%(t(eta_tilde_p)%*%as.matrix(y[pp,]))
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
    cholSig_xi_n_trans = chol2inv(diag(nrow=k)+zeta_theta_invSig%*%theta_zeta_n)
    m_xi_n = cholSig_xi_n_trans%*%(t(cholSig_xi_n_trans)%*%(zeta_theta_invSig%*%y_tilde_n))
    xi[,nn] = m_xi_n + cholSig_xi_n_trans%*%as.matrix(rnorm(k,1))
  }
  return(xi)
}
