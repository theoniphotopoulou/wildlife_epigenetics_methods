# functions used for probage

invlogit <- function(x){exp(x)/(exp(x) + 1)}
logit <- function(x){log(x/(1-x))}

beta_mean = function(t, eta0, omega, p){ eta0 + exp(-omega * t) * (p - eta0)}

beta_var = function(t, eta0, omega, p, N, cc){ 
  eta1 = eta0
  varterm0 = eta0 * eta1
  varterm1 =  (1-p)*eta0^2 + p*eta1^2
  varterm0 / N + 
    exp(-omega * t) * (varterm1 - varterm0) / N + 
    exp(-2*omega*t) * ( (cc/N^2) - varterm1 / N)
}

# only valid beta if var <= mu(1-mu)
beta_a = function(mean, var){ mean^2 * (1 - mean) / var - mean }
beta_b = function(mean, var){ (mean^2 * (1 - mean) / var - mean) * (1/mean - 1) }

beta_ab = function(mean, var){ 
  k = (mean*(1-mean)/var)-1
  a = mean*k
  b = (1-mean)*k
  return(list(a, b))
}

llh = function(pars = c(eta0, omega, p, N, cc), x, ages, weights = NULL){
  if(is.null(weights)){weights = rep(1, length(ages))}
  eta0 = invlogit(pars[1])
  omega = exp(pars[2])
  p = invlogit(pars[3])
  N = exp(pars[4])
  cc = exp(pars[5])
  mu_j = beta_mean(ages, eta0, omega, p)
  sigma2_j = beta_var(ages, eta0, omega, p, N, cc)
  valid_beta = (sigma2_j <= mu_j * (1 - mu_j))
  pen = ifelse(valid_beta, 0, 1e6)
  a = beta_a(mu_j, sigma2_j)
  b = beta_b(mu_j, sigma2_j)
  val = mapply(function(x, a, b) dbeta(x, a, b, log = TRUE), x = x, a = a, b = b)
  val = val[!is.nan(val)] 
  wts = weights[!is.nan(val)] 
  val = -sum(wts * val) + sum(pen)
  # val = sum(weights * mapply(function(x, a, b) dbeta(x, a, b, log = TRUE), x = x, a = a, b = b))
  # val = -val + sum(weights * pen)
}

llh2 = function(pars = c(aa, bb), x, ages, eta0, omega, p0, N, cc){
  aa = pars[1]
  bb = pars[2]
  new_eta0 = invlogit(logit(eta0) + bb)
  new_p = invlogit(logit(p0) + bb)
  new_omega = exp(log(omega) + aa)
  # mu_j = sapply(ages, beta_mean, eta0 = new_eta0, omega = new_omega, p = new_p)
  mu_j = beta_mean(ages, new_eta0, new_omega, new_p)
  # sigma2_j = sapply(ages, beta_var, eta0 = new_eta0, omega = new_omega, p = new_p, N = N, cc = cc)
  sigma2_j = beta_var(ages, new_eta0, new_omega, new_p, N, cc)
  valid_beta = (sigma2_j <= mu_j * (1 - mu_j))
  pen = ifelse(valid_beta, 0, 1e6)
  a = beta_a(mu_j, sigma2_j)
  b = beta_b(mu_j, sigma2_j)
  val = mapply(function(x, a, b) dbeta(x, a, b, log = TRUE), x = x, a = a, b = b)
  #val[is.nan(val)] = -1e6
  val = val[!is.nan(val)] 
  val = -sum(val) + sum(pen)
  val
}

probage = function(betas, age, n_cpg, weights = NULL){
  
  cors = cor(betas, age)
  inds = order(abs(cors), decreasing = TRUE)[1:n_cpg]
  betas = betas[, inds]

  pars_per_site <- list()
  for(i in 1:ncol(betas)){
    tt <- optim(c(0, 0, 0, 1, 1), fn = llh, x = betas[,i], ages = age, weights = weights)
    pars_per_site[[i]] <- data.frame(eta0 = invlogit(tt$par[1]), omega = exp(tt$par[2]), p0 = invlogit(tt$par[3]), N = exp(tt$par[4]), cc = exp(tt$par[5]), ll = tt$value)
  }
  pars_per_site <- do.call(rbind, pars_per_site)
  
  # drop any sites with huge omega or with eta0 or p0 close to 1 (saturating sites)
  include_sites1 = (pars_per_site$omega < 1000)
  include_sites2 = (pars_per_site$eta0 < 0.95)
  include_sites3 = (pars_per_site$p0 < 0.95)
  include_sites = include_sites1 & include_sites2 & include_sites3 
  pars_per_site = pars_per_site[include_sites, ]
  betas = betas[, include_sites]
  
  pars_per_animal <- list()
  for(i in 1:nrow(betas)){
    tt2 <- try(optim(c(0, 0), fn = llh2, x = betas[i, ], 
                     ages = age[i], eta0 = pars_per_site$eta0, omega = pars_per_site$omega, 
                     p0 = pars_per_site$p, N = pars_per_site$N, cc = pars_per_site$cc), silent = TRUE)
    if(class(tt2) != "try-error"){
      pars_per_animal[[i]] <- data.frame(
        acc = tt2$par[1],
        bias = tt2$par[2],
        ll = tt2$value,
        old_eta0_s2 = pars_per_site$eta0[2],
        old_p_s2 = pars_per_site$p0[2],
        old_omega_s2 = pars_per_site$omega[2],
        new_eta0_s2 = invlogit(logit(pars_per_site$eta0[2]) + tt2$par[2]), 
        new_p_s2 = invlogit(logit(pars_per_site$p0[2]) + tt2$par[2]),
        new_omega_s2 = exp(log(pars_per_site$omega[2]) + tt2$par[1]))
    } else {
      pars_per_animal[[i]] <- data.frame(acc = NA, bias = NA, old_eta0_s2 = NA, old_p_s2 = NA, old_omega_s2 = NA, 
                                         new_eta0_s2 = NA, new_p_s2 = NA, new_omega_s2 = NA)
    }
  }
  pars_per_animal <- do.call(rbind, pars_per_animal)
  
  pars_per_animal
  
}

# ww <- runif(n = nrow(x))
# xx <- probage(betas = x, age = y, n_cpg = 500, weights = ww)

probage_fit_sites = function(betas, age, n_cpg, weights = NULL){
  
  cors = cor(betas, age)
  inds = order(abs(cors), decreasing = TRUE)[1:n_cpg]
  betas = betas[, inds]
  
  pars_per_site <- list()
  for(i in 1:ncol(betas)){
    tt <- optim(c(0, 0, 0, 1, 1), fn = llh, x = betas[,i], ages = age, weights = weights)
    pars_per_site[[i]] <- data.frame(eta0 = invlogit(tt$par[1]), omega = exp(tt$par[2]), p0 = invlogit(tt$par[3]), N = exp(tt$par[4]), cc = exp(tt$par[5]), ll = tt$value)
  }
  pars_per_site <- do.call(rbind, pars_per_site)
  
  # drop any sites with huge omega or with eta0 or p0 close to 1 (saturating sites)
  include_sites1 = (pars_per_site$omega < 1000)
  include_sites2 = (pars_per_site$eta0 < 0.95)
  include_sites3 = (pars_per_site$p0 < 0.95)
  include_sites = include_sites1 & include_sites2 & include_sites3 
  pars_per_site = pars_per_site[include_sites, ]
  betas = betas[, include_sites]
  
  return(list(age_cor_cpgs = inds, noerror_cpgs = include_sites, pars_per_site = pars_per_site))

}

probage_fit_animals = function(betas, age, site_results){
  
  inds = site_results$age_cor_cpgs
  include_sites = site_results$noerror_cpgs
  pars_per_site = site_results$pars_per_site
  
  betas = betas[, inds]
  betas = betas[, include_sites]

  pars_per_animal <- list()
  for(i in 1:nrow(betas)){
    tt2 <- try(optim(c(0, 0), fn = llh2, x = betas[i, ], 
                     ages = age[i], eta0 = pars_per_site$eta0, omega = pars_per_site$omega, 
                     p0 = pars_per_site$p, N = pars_per_site$N, cc = pars_per_site$cc), silent = TRUE)
    if(class(tt2) != "try-error"){
      pars_per_animal[[i]] <- data.frame(
        acc = tt2$par[1],
        bias = tt2$par[2],
        ll = tt2$value,
        old_eta0_s2 = pars_per_site$eta0[2],
        old_p_s2 = pars_per_site$p0[2],
        old_omega_s2 = pars_per_site$omega[2],
        new_eta0_s2 = invlogit(logit(pars_per_site$eta0[2]) + tt2$par[2]), 
        new_p_s2 = invlogit(logit(pars_per_site$p0[2]) + tt2$par[2]),
        new_omega_s2 = exp(log(pars_per_site$omega[2]) + tt2$par[1]))
    } else {
      pars_per_animal[[i]] <- data.frame(acc = NA, bias = NA, old_eta0_s2 = NA, old_p_s2 = NA, old_omega_s2 = NA, 
                                         new_eta0_s2 = NA, new_p_s2 = NA, new_omega_s2 = NA)
    }
  }
  pars_per_animal <- do.call(rbind, pars_per_animal)
  
  pars_per_animal
  
}

# xx <- probage_fit_sites(betas = x, age = y, n_cpg = 500)
# xx2 <- probage_fit_animals(betas = x, age = y, site_results = xx)

