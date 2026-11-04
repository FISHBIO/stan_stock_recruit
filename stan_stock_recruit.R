# File-Name:       stan_stock_recruit
# Version:         v.1
# Date:            2026-11-02
# Author:          Rob Eckelbecker
# Email(s):        Roberteckelbecker@fishbio.com
# Description:     Statistical methods for manuscript. Create stock recruit relationship between Oakdale and Caswell with covariates.
#                  Simulation modeling to estimate abundance under different water management scenarios.
#                  Ver 1. Used for Drought and deluge manuscript
# Functions:       NONE

## Oakdale juvenile production and Oakdale-Caswell survival ##
#### Packages and file loading ####
# install.packages(c('readxl', 'jagsUI', 'coda'))
library(readxl)
library(jagsUI)
library(coda)

#JAGS has an error on my computer and named different
#work around to rename the version to allow jags to work
#likely this section is not necessary for others
original_jags_version <- rjags::jags.version
unlockBinding("jags.version", asNamespace("rjags"))

assign(
  "jags.version",
  function() package_version("4.3.2"),
  envir = asNamespace("rjags")
)
lockBinding("jags.version", asNamespace("rjags"))


#### Files and model settings #####
# Import and prepare data
file = 'SR-data_master_v9.xlsx'
output.dir = 'OakCas_model_comparison'
dir.create(output.dir, showWarnings = FALSE)

stan = as.data.frame(read_excel(file, sheet = 'Data'))
covs = as.data.frame(read_excel(file, sheet = 'Covariates'))
idx = match(stan$`Brood Year`, covs$Brood_Year)
if(anyNA(idx)) stop('Missing brood years in Covariates sheet.')
covs = covs[idx, , drop = FALSE]
stopifnot(all(stan$`Brood Year` == covs$Brood_Year))

Nyrs = nrow(stan)
years = stan$`Brood Year`
Fspawn = as.numeric(stan$Females)

## Approximate female SD from reported total-spawner SD.
sex.ratio = Fspawn / as.numeric(stan$`Total Spawners`)
Fspawn_sd = as.numeric(stan$`Spawner SD`) * sex.ratio
if(any(!is.finite(Fspawn_sd)) || any(Fspawn_sd < 0)) stop('Invalid spawner SD')
Fspawn_sd[Fspawn_sd == 0] = 0.001
Fspawn_tau = 1 / Fspawn_sd^2

OakJuv = as.numeric(stan$`Oakdale Juveniles`)
CasJuv = as.numeric(stan$`Caswell Juveniles`)
OakJuv_sd = as.numeric(stan$`Oakdale SD`)
CasJuv_sd = as.numeric(stan$`Cas Juv SD`)

## Manuscript describes log-space SD as CV (SD / mean).
## Observed uncertainty on log-abundance scale
OakJuv_cv = OakJuv_sd / OakJuv
CasJuv_cv = CasJuv_sd / CasJuv
lgobsR_Oak = log(OakJuv)
lgobsR_Cas = log(CasJuv)
if(any(!is.finite(c(Fspawn, Fspawn_tau, lgobsR_Oak, lgobsR_Cas,
                    OakJuv_cv, CasJuv_cv))) ||
   any(Fspawn <= 0) || any(OakJuv_cv <= 0) || any(CasJuv_cv <= 0)) {
  stop('Missing, zero, or invalid observations or uncertainties.')
}

muJuv = mean(OakJuv)
sdJuv = sd(OakJuv)

## Manuscript AM2 (temperature date) = workbook AM4.
## Manuscript AM3 (reservoir storage) = workbook AM5.
## Workbook AM2 and AM3 describe different biological covariates.
covs$AM2_paper = covs$AM4
covs$AM3_paper = covs$AM5
Oak.vars = c('AM1', 'AM2_paper', 'AM3_paper',
             'IN1', 'IN2', 'IN3', 'ER1', 'ER2', 'ER3')
Cas.vars = c('EM1', 'EM2', 'EM3', 'EM4', 'EM5',
             'LM1', 'LM2', 'LM3', 'LM4', 'LM5')

## Check covariate completeness and variation
all.vars = unique(c(Oak.vars, Cas.vars))
for(v in all.vars) {
  z = as.numeric(covs[[v]])
  if(any(!is.finite(z)) || sd(z) == 0) stop(paste('Invalid covariate:', v))
}

## Standardize using the same years for all candidate models
cov.std = covs
for(v in all.vars) cov.std[[v]] = as.numeric(scale(as.numeric(covs[[v]])))

## Candidate model table
## null + Oakdale-only + Caswell-only + all Oakdale x  Caswell interactions
## NA indicates no covariate at that location.
models = data.frame(
  Oak = c(
    NA_character_,
    Oak.vars,
    rep(NA_character_, length(Cas.vars))),
  Cas = c(
    NA_character_,
    rep(NA_character_, length(Oak.vars)),
    Cas.vars),stringsAsFactors = FALSE)

## Add every Oakdale x Caswell covariate combination
pairs = expand.grid(
  Oak = Oak.vars,
  Cas = Cas.vars,
  stringsAsFactors = FALSE)

models = rbind(models, pairs)

## Should always be 110 candidate models
stopifnot(nrow(models) == 110)
models$Model = ifelse(is.na(models$Oak) & is.na(models$Cas), 'Null',
                      ifelse(is.na(models$Oak), paste0('Cas_', models$Cas),
                             ifelse(is.na(models$Cas), paste0('Oak_', models$Oak),
                                    paste0(models$Oak, '_', models$Cas))))
models = models[, c('Model', 'Oak', 'Cas')]
models$Model = gsub('AM2_paper', 'AM2', models$Model, fixed = TRUE)
models$Model = gsub('AM3_paper', 'AM3', models$Model, fixed = TRUE)
print(models)

## MCMC settings 
use.density = FALSE
n.chains = 3
n.iter = 200000
n.burnin = 100000
n.thin = 10

#### Fit candidate models ####
## Write only parameters used by each candidate model. This avoids unconstrained 
#covariate coefficients in the null models.
results = vector('list', nrow(models))
fits = vector('list', nrow(models))
names(fits) = models$Model

for(k in 1:nrow(models)) {
  Oak.cov = models$Oak[k]
  Cas.cov = models$Cas[k]
  has.oak = !is.na(Oak.cov)
  has.cas = !is.na(Cas.cov)
  cat('\nModel', k, 'of', nrow(models), ':', models$Model[k], '\n')
  
  X_Oak = if(has.oak) cov.std[[Oak.cov]] else rep(0, Nyrs)
  X_Cas = if(has.cas) cov.std[[Cas.cov]] else rep(0, Nyrs)
  
  oak.prior = if(has.oak) 'gOak ~ dnorm(0, 1.0E-3)' else ''
  cas.prior = if(has.cas) 'gCas ~ dnorm(0, 1.0E-3)' else ''
  oak.term = if(has.oak) ' + gOak * X_Oak[i]' else ''
  cas.term = if(has.cas) ' + gCas * X_Cas[i]' else ''
  den.prior = if(use.density) 'bden ~ dnorm(0, 1.0E-3)' else ''
  den.term = if(use.density) ' + bden * stdJuv[i]' else ''
  den.line = if(use.density) 'stdJuv[i] <- (pred_R_Oak[i] - muJuv) / sdJuv' else ''
  
  model.txt = paste0('model {
',
  '  alpha ~ dnorm(0, 1.0E-3)
  beta ~ dnorm(0, 1.0E-3)',
  '  ', oak.prior, '
  ', cas.prior, '
  ', den.prior, '',
'  b0 ~ dnorm(0, 1.0E-3)',
'  tau_pro_Oak ~ dgamma(0.01, 0.01)',
'  tau_pro_Cas ~ dgamma(0.01, 0.01)',
'  sd_pro_Oak <- pow(tau_pro_Oak, -0.5)',
'  sd_pro_Cas <- pow(tau_pro_Cas, -0.5)',
'  for(i in 1:Nyrs) {
',
'    SP[i] ~ dnorm(Fspawn[i], Fspawn_tau[i]) T(0,)',
'    lgpred_R_Oak[i] <- log(SP[i]) + alpha + beta * SP[i]', oak.term, '',
'    pred_R_Oak[i] <- exp(lgpred_R_Oak[i])',
'    Tautot_Oak[i] <- 1 / (pow(OakJuv_cv[i], 2) + pow(sd_pro_Oak, 2))',
'    lgobsR_Oak[i] ~ dnorm(lgpred_R_Oak[i], Tautot_Oak[i])',
'    ', den.line, '',
'    logit(surv[i]) <- b0', cas.term, den.term, '',
'    lgpred_R_Cas[i] <- lgpred_R_Oak[i] + log(surv[i])',
'    pred_R_Cas[i] <- exp(lgpred_R_Cas[i])',
'    Tautot_Cas[i] <- 1 / (pow(CasJuv_cv[i], 2) + pow(sd_pro_Cas, 2))',
'    lgobsR_Cas[i] ~ dnorm(lgpred_R_Cas[i], Tautot_Cas[i])',
'  }
}
')
  model.file = file.path(output.dir, paste0(models$Model[k], '.jags'))
  writeLines(model.txt, model.file)
  
  jags.data = list(Nyrs = Nyrs, Fspawn = Fspawn, Fspawn_tau = Fspawn_tau,
                   X_Oak = X_Oak, X_Cas = X_Cas,
                   OakJuv_cv = OakJuv_cv, CasJuv_cv = CasJuv_cv,
                   lgobsR_Oak = lgobsR_Oak, lgobsR_Cas = lgobsR_Cas)
  if(use.density) {
    jags.data$muJuv = muJuv
    jags.data$sdJuv = sdJuv
  }
  params = c('alpha','beta','b0','sd_pro_Oak','sd_pro_Cas',
             'surv','pred_R_Oak','pred_R_Cas',
             'lgpred_R_Oak','lgpred_R_Cas','deviance')
  if(has.oak) params = c(params, 'gOak')
  if(has.cas) params = c(params, 'gCas')
  if(use.density) params = c(params, 'bden')
  
  inits = function() {
    ans = list(alpha = -2, beta = -0.0001, b0 = -2,
               tau_pro_Oak = 1, tau_pro_Cas = 1,
               SP = pmax(Fspawn, 0.001))
    if(has.oak) ans$gOak = 0
    if(has.cas) ans$gCas = 0
    if(use.density) ans$bden = 0
    ans
  }
  
  #  set.seed(123 + k)
  mod = jags(data = jags.data, inits = inits,
             parameters.to.save = params, model.file = model.file,
             n.chains = n.chains, n.iter = n.iter,
             n.burnin = n.burnin, n.thin = n.thin,
             DIC = TRUE, parallel = FALSE)
  fits[[k]] = mod
  p = as.matrix(mod$samples)
  s = mod$summary
  get.mean = function(v) if(v %in% rownames(s)) s[v, 'mean'] else NA_real_
  get.sd = function(v) if(v %in% rownames(s)) s[v, 'sd'] else NA_real_
  get.prob = function(v) if(v %in% colnames(p)) mean(p[, v] > 0) else NA_real_
  get.rhat = function() {
    x = s[, 'Rhat']
    if(all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  }
  results[[k]] = data.frame(Model = models$Model[k], Oak = Oak.cov,
                            Cas = Cas.cov, alpha_mean = get.mean('alpha'), beta_mean = get.mean('beta'),
                            b0_mean = get.mean('b0'), gOak_mean = get.mean('gOak'),
                            gOak_SD = get.sd('gOak'), p_gOak_gt_0 = get.prob('gOak'),
                            gCas_mean = get.mean('gCas'), gCas_SD = get.sd('gCas'),
                            p_gCas_gt_0 = get.prob('gCas'),
                            bden_mean = get.mean('bden'),
                            sd_pro_Oak = get.mean('sd_pro_Oak'),
                            sd_pro_Cas = get.mean('sd_pro_Cas'),
                            DIC = mod$DIC, max_Rhat = get.rhat())
  cat('DIC:', round(mod$DIC, 2), '\n')
  saveRDS(mod, file.path(output.dir, paste0(models$Model[k], '.rds')))
}

## 5. Compare models -------------------------------------------
comparison.table = do.call(rbind, results)
comparison.table = comparison.table[order(comparison.table$DIC), ]
comparison.table$DeltaDIC = comparison.table$DIC - min(comparison.table$DIC)
comparison.table = comparison.table[, c('Model','Oak','Cas','DIC','DeltaDIC',
                                        'alpha_mean','beta_mean','b0_mean','gOak_mean','gOak_SD','p_gOak_gt_0',
                                        'gCas_mean','gCas_SD','p_gCas_gt_0',
                                        'sd_pro_Oak','sd_pro_Cas','max_Rhat')]
## R-squared: squared Pearson correlation on log observed/predicted abundance.
comparison.table$Oak_R2 = NA_real_
comparison.table$Cas_R2 = NA_real_
for(k in seq_len(nrow(comparison.table))) {
  m = fits[[comparison.table$Model[k]]]
  oak.pred = m$summary[paste0('lgpred_R_Oak[', seq_len(Nyrs), ']'), 'mean']
  cas.pred = m$summary[paste0('lgpred_R_Cas[', seq_len(Nyrs), ']'), 'mean']
  comparison.table$Oak_R2[k] = cor(lgobsR_Oak, oak.pred)^2
  comparison.table$Cas_R2[k] = cor(lgobsR_Cas, cas.pred)^2
}
print(comparison.table, row.names = FALSE, digits = 3)
write.csv(print(comparison.table, row.names = FALSE, digits = 3), "comparison.table.csv")

## Save all fits so you can inspect or simulate without rerunning
saveRDS(list(models = models, comparison = comparison.table, fits = fits,
             years = years, standardized_covariates = cov.std),
        file.path(output.dir, 'all_model_fits.rds'))



##### Projected Simulations under different management scenarios ####
## Read the projection data
projdat <- as.data.frame(
  read_excel(
    "SR_covariate_projections_v2.xlsx",
    sheet = "Covariates"
  )
)
#Restrict projecions to the same brood years ued to fit model
projdat <- projdat[projdat$`Brood Year` %in% stan$`Brood Year`,,drop = FALSE]
projdat <- projdat[order(projdat$`Brood Year`),,drop = FALSE]
## Verify year alignment
stopifnot(
  nrow(projdat) == nrow(stan),
  all(projdat$`Brood Year` == stan$`Brood Year`)
)
projdat$`Brood Year`

## Posterior from the fitted model
sim.mod <- fits[["AM3_EM2"]]

p <- as.matrix(sim.mod$samples)

dim(p)
colnames(p)

### Projection Inputs
muX_oak <- mean(covs$AM5)
sdX_oak <- sd(covs$AM5)

muX_cas <- mean(covs$EM2)
sdX_cas <- sd(covs$EM2)

Nyrs <- nrow(projdat)

mean_Fspawn <- mean(Fspawn)
Fspawn_proj <- rep(mean_Fspawn,Nyrs)

## Define Scenarios
scenarios <- list(
  Historical = list(Oak = projdat$Obs_AM5, Cas = projdat$Obs_EM2),
  SRP = list(Oak = projdat$SRP_AM5, Cas = projdat$SRP_EM2),
  P40 = list(Oak = projdat$P40_AM5, Cas = projdat$P40_EM2)
)


### Simulation
Nmcmc <- nrow(p)
Nyrs <- nrow(projdat)
Nscens <- length(scenarios)

## Storage
OakStats <- array(NA_real_, dim = c(Nscens, Nyrs, 3),
                  dimnames = list(names(scenarios), projdat$`Brood Year`, c("p10", "mean", "p90")))
CasStats <- OakStats


## Generate process deviations
pro_dev_oak <- matrix(NA_real_, nrow = Nmcmc, ncol = Nyrs)
pro_dev_cas <- matrix(NA_real_, nrow = Nmcmc, ncol = Nyrs)

for(j in seq_len(Nmcmc)) {
  #  pro_dev_oak[j, ] <- rnorm(Nyrs, mean = 0, sd = p[j, "sd_pro_Oak"])
  #  pro_dev_cas[j, ] <- rnorm(Nyrs, mean = 0, sd = p[j, "sd_pro_Cas"])
  sigma_oak <- p[j, "sd_pro_Oak"]
  sigma_cas <- p[j, "sd_pro_Cas"]
  
  pro_dev_oak[j, ] <- rnorm(
    Nyrs,
    mean = -0.5 * sigma_oak^2,
    sd = sigma_oak
  )
  
  pro_dev_cas[j, ] <- rnorm(
    Nyrs,
    mean = -0.5 * sigma_cas^2,
    sd = sigma_cas
  )
}


for(iscen in seq_along(scenarios)) {
  scen <- scenarios[[iscen]]
  
  ## Standardize using FITTING-data mean and SD
  X_oak <- (scen$Oak - muX_oak) / sdX_oak
  X_cas <- (scen$Cas - muX_cas) / sdX_cas
  
  for(iyr in seq_len(Nyrs)) {
    
    ## Oakdale juvenile abundance
    log_JuvOak <- log(Fspawn_proj[iyr]) + p[, "alpha"] + p[, "beta"] * Fspawn_proj[iyr] +
      p[, "gOak"] * X_oak[iyr] + pro_dev_oak[, iyr]
    JuvOak <- exp(log_JuvOak)
    
    ## Oakdale -> Caswell survival
    logit_surv <- p[, "b0"] + p[, "gCas"] * X_cas[iyr]
    surv <- plogis(logit_surv)
    if(!exists("SurvStats")) {
      SurvStats <- matrix(
        NA_real_,
        nrow = Nscens,
        ncol = Nyrs,
        dimnames = list(
          names(scenarios),
          projdat$`Brood Year`)
      )
    }
    SurvStats[iscen, iyr] <- mean(surv)
    
    ## Caswell abundance
    log_JuvCas <- log(JuvOak * surv) + pro_dev_cas[, iyr]
    JuvCas <- exp(log_JuvCas)

    ## Save posterior summaries
    OakStats[iscen, iyr, "p10"] <- quantile(JuvOak, 0.10)
    OakStats[iscen, iyr, "mean"] <- mean(JuvOak)
    OakStats[iscen, iyr, "p90"] <- quantile(JuvOak, 0.90)
    
    CasStats[iscen, iyr, "p10"] <- quantile(JuvCas, 0.10)
    CasStats[iscen, iyr, "mean"] <- mean(JuvCas)
    CasStats[iscen, iyr, "p90"] <- quantile(JuvCas, 0.90)
  }
}

#Run output
simulation_results <- do.call(
  rbind,
  lapply(
    names(scenarios),
    function(s) {
      data.frame(
        BroodYear = projdat$`Brood Year`,
        Scenario = s,
        Oak_p10 = OakStats[s, , "p10"],
        Oak_mean = OakStats[s, , "mean"],
        Oak_p90 = OakStats[s, , "p90"],
        Cas_p10 = CasStats[s, , "p10"],
        Cas_mean = CasStats[s, , "mean"],
        Cas_p90 = CasStats[s, , "p90"])
    }
  )
)

rownames(simulation_results) <- NULL

simulation_results

scenario_table <- do.call(
  rbind,lapply(
    names(scenarios),
    function(s) {
      x <- simulation_results[
        simulation_results$Scenario == s,]
      
      data.frame(
        Scenario = s,
        Oak_mean = mean(x$Oak_mean),
        Oak_SD = sd(x$Oak_mean),
        Cas_mean = mean(x$Cas_mean),
        Cas_SD = sd(x$Cas_mean))
    }
  )
)

scenario_table
scenario_table[, 2:5] <- scenario_table[, 2:5] / 1000
scenario_table

write.csv(simulation_results, "scenario_table.csv")
