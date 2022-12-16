#############################################
#### this is the script to use the NLHPC ####
#############################################

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# read the files to use the NLHPC
sample_dir = readRDS("sample_dir.rds")
control = readRDS("control.rds")
grid = readRDS("grid.rds")
grid_2 = readRDS("grid_2.rds")
beta_margin_params = readRDS("beta_margin.rds")
weights = readRDS("weights.rds")

require(tidyr)
require(dplyr)
require(maxLik)

# it handles the marginal betas decomposing them onto nutritional
# statuses and considering an extra one in order to include those
# non-observed values using a mixture model computing it via
# expectation-maximization algorithm

# estimate the starting distributions via method of
# moments (actually the parameters in control object)
init_estim = function(x = control, sex, tran, ns) {
  params = control %>%
    filter(NS == ns,
           Sex == sex,
           Transition == ifelse(ns == "Underweight" & tran == "Decrease",
                                "Remain", ifelse(ns == "Obese" & tran == "Increase",
                                                 "Remain", tran))) %>%
    dplyr::select(c("Mean", "SD")) %>%
    as.matrix() %>% as.numeric()
  # params[1]: mean; params[2]: sd
  aux = params[1]*(1-params[1])/(params[2]^2) - 1
  return(list(alpha = params[1]*aux, beta = (1-params[1])*aux))
}

for (i in 1:nrow(grid)) {
  grid[i, "alpha"] = init_estim(sex = grid$sex[i],
                                tran = grid$tran[i],
                                ns = grid$ns[i])$alpha
  grid[i, "beta"] = init_estim(sex = grid$sex[i],
                               tran = grid$tran[i],
                               ns = grid$ns[i])$beta
}

# mixture density
beta_mixture_dens = function(x, mixtures, alphas, betas, log = T) {
  len = length(mixtures)
  dens = 0
  for (i in 1:len) {
    dens = dens + mixtures[i]*dbeta(x, alphas[i], betas[i])
  }
  if (log) { return(log(dens)) }
  else {return(dens)}
}

# likelihood
likelihood = function(theta, x) {
  if (any(theta <= 0)) {
    NA
  } else {
    len = length(theta)
    alphas = theta[seq(1, len, by = 2)]
    betas = theta[seq(2, len, by = 2)]
    beta_mixture_dens(x, mixture_params, alphas, betas, log = T)
  }
}

# take the marginal beta density function
beta_goal = function(x, id) {
  shapes = beta_margin_params %>%
    filter(Sex == grid_2$sex[id],
           Transition == grid_2$tran[id]) %>%
    dplyr::select(c("shape1", "shape2")) %>%
    as.numeric()
  
  dbeta(x, shape1 = shapes[1], shape2 = shapes[2])
}

# compute the minimum in order to get the intersection between
# marginal beta and the mixture model
area = function(x, id_1, mixtures, alphas, betas, log = FALSE) {
  pmin(beta_goal(x, id_1),
       beta_mixture_dens(x, mixtures, alphas, betas, log))
}

# --------- expectation maximization algorithm ---------
# get the initial parameters obtained from the method of moments
# there are ordered by nutritional status and where first four
# are the shape1 and second four are the shape2
# while the last two are the shape1 and shape2 for non-observed
shape_ns = grid %>%
  filter(sex == grid_2$sex[as.numeric(args[1])],
         tran == grid_2$tran[as.numeric(args[1])]) %>%
  dplyr::select(c("ns", "alpha", "beta")) %>%
  gather(key = "param", value = "value",
         -ns) %>%
  dplyr::select(c("value")) %>% pull()
# if we add non-observed (uniform)
# c(., rep(1, 2))

# states the initial mixture weights
# if we add non-observed
# rep(1/5, 5)
# mixture_params = rep(1/4, 4)
mixture_params = weights %>%
  filter(Sex == grid_2$sex[as.numeric(args[1])],
         Transition == grid_2$tran[as.numeric(args[1])]) %>%
  dplyr::select(c("weight")) %>%
  pull()

if(grid_2$w_inverse[as.numeric(args[1])] == "Yes") {
  print("------------ Transforming weights ------------")
  mixture_params[1] = 1 - mixture_params[1]
  mixture_params[2] = 1 - mixture_params[2]
  mixture_params[3] = 1 - mixture_params[3]
  mixture_params[4] = 1 - mixture_params[4]
  total = mixture_params[1] + mixture_params[2] + mixture_params[3] + mixture_params[4]
  mixture_params[1] = mixture_params[1]/total
  mixture_params[2] = mixture_params[2]/total
  mixture_params[3] = mixture_params[3]/total
  mixture_params[4] = mixture_params[4]/total
}

print("------- starting weights -------")
print(paste0("Underweight:", mixture_params[1]))
print(paste0("Normal:", mixture_params[2]))
print(paste0("Overweight:", mixture_params[3]))
print(paste0("Obese:", mixture_params[4]))

# subset the dirichlet random numbers
data = sample_dir %>%
  filter(Sex == grid_2$sex[as.numeric(args[1])],
         Transition == grid_2$tran[as.numeric(args[1])]) %>%
  dplyr::select(c("prob")) %>% pull() %>%
  as.numeric()

# ---- maximization step ----
tol = 100 # starting tolerance criteria
itmax = 1

# counter in case we use the uniform transformation
count_retry = 0
count_retry_max = 0 # maximum times of retrying until transform
count_trans = 0
count_trans_max = 0 # maximum times of transformation allowed

# parameters to optimize in each density
component_parameters = c(shape_ns[1], shape_ns[5], # underweight
                         shape_ns[2], shape_ns[6], # normal
                         shape_ns[3], shape_ns[7], # overweight
                         shape_ns[4], shape_ns[8]) # obese
                       # shape_ns[9], shape_ns[10]) # non-observed

max_lik = 0 # starting loglik

set.seed(12345)

res = data.frame(id = as.numeric(args[1]),
                 Iteration = 0,
                 Loglik = 0,
                 area = 0,
                 Tol = 0,
                 Method = 0,
                 Time = 0,
                 w_inverse = 0,
                 w_un = mixture_params[1],
                 alpha_un = component_parameters[1],
                 beta_un = component_parameters[2],
                 w_no = mixture_params[2],
                 alpha_no = component_parameters[3],
                 beta_no = component_parameters[4],
                 w_ov = mixture_params[3],
                 alpha_ov = component_parameters[5],
                 beta_ov = component_parameters[6],
                 w_ob = mixture_params[4],
                 alpha_ob = component_parameters[7],
                 beta_ob = component_parameters[8])

t0 = proc.time()[3] 
t_it = proc.time()[3]-t0
while (tol >= 1e-10 & itmax <= 600 & t_it <= 200) {
  if(grid_2$algo[as.numeric(args[1])] == "NM") {
    # Nelder-Mead
    print("------------ Nelder-Mead algorithm ------------")
    method = "NM"
    new_component_params = maxLik::maxLik(likelihood,
                                          start = component_parameters,
                                          method = "NM",
                                          x = data,
                                          control = list(iterlim = 500,   # 500
                                                         tol = 1e-10,
                                                         nm_alpha = 1,    # 1
                                                         nm_beta = 0.5,   # 0.5
                                                         nm_gamma = 2))   # 2
  } else {
    if(grid_2$algo[as.numeric(args[1])] == "SANN") {
      # Simulated Annealing
      print("------------ Simulated Annealing ------------")
      method = "SANN"
      new_component_params = maxLik::maxLik(likelihood,
                                            start = component_parameters,
                                            method = "SANN",
                                            x = data,
                                            control = list(iterlim = 10000, # 10000
                                                           tol = 1e-10,
                                                           # sann_cand = 2.5,        # Gaussian
                                                           sann_temp = 10,           # 10
                                                           sann_tmax = 10,           # 10
                                                           sann_randomSeed = 12345)) # 123
    } else {
      # use Broyden, Fletcher, Goldfarb and Shanno (quasi-Newton)
      print("------------ BFGS ------------")
      method = "BFGS"
      new_component_params = maxLik::maxLik(likelihood,
                                            start = component_parameters,
                                            method = "BFGS",
                                            x = data,
                                            control = list(iterlim = 200, # 200
                                                           tol = 1e-10)) # 123
    }
  }
  
  print(paste0("Max Loglik = ", new_component_params$maximum))
  
  # if the optimization does not get a value stop the process
  if(is.na(new_component_params$maximum)) {
    print("---------- All algorithms failed, then stop ----------")
    break }
  
  print(paste0("Current loglik value = ", max_lik))
  
  tol = (new_component_params$maximum - max_lik)/max_lik
  
  if(tol > 0) {
    # tol = tol
    # r_alpha = rep(0, 8)
    t_it = proc.time()[3]-t0
    print(paste0(grid_2$sex[as.numeric(args[1])], ", ",
                 grid_2$tran[as.numeric(args[1])], ". tol = ",
                 tol, ", time = ", t_it, ", iteration = ", itmax))
    
    # reset the counter
    count_retry = 0
  } else {
    
    print("")
    print("Algorithm did not get a better solution")
    t_it = proc.time()[3]-t0
    print(paste0(grid_2$sex[as.numeric(args[1])], ", ",
                 grid_2$tran[as.numeric(args[1])], ". tol = ",
                 tol, ", time = ", t_it, ", iteration = ", itmax))
    print("")
    print("------------------ Retrying ------------------")
    print("")
    
    # increase the counter
    count_retry = count_retry + 1
    
    # if the retrying is repeated too much times break
    if(count_retry == 50) {
      print("------ Maximum iterations withouth improvement achieved ------")
      break
    }
    
    tol = 100
    # component_parameters = runif(8, 0.9, 1.1)*component_parameters
    new_params = new_component_params$estimate
    sum_loglik = new_component_params$maximum + max_lik
    prop_new = new_component_params$maximum/sum_loglik
    prop_old = max_lik/sum_loglik
    component_parameters = runif(1, 0.85, 1.15)*(prop_new*new_params + prop_old*component_parameters)
    next
    }
  
  # update the maximum likelihood value
  max_lik = new_component_params$maximum

  # updating the parameters
  component_parameters = new_component_params$estimate
  
  dbeta_un = mixture_params[1]*dbeta(data, component_parameters[1],
                                     component_parameters[5])
  dbeta_no = mixture_params[2]*dbeta(data, component_parameters[2],
                                     component_parameters[6])
  dbeta_ov = mixture_params[3]*dbeta(data, component_parameters[3],
                                     component_parameters[7])
  dbeta_ob = mixture_params[4]*dbeta(data, component_parameters[4],
                                     component_parameters[8])
  # dbeta_na = mixture_params[5]*dbeta(data, component_parameters[9],
  #                                    component_parameters[10])
  # updating mixture weights
  responsibility_denom = dbeta_un + dbeta_no + dbeta_ov + dbeta_ob #+ dbeta_na
  un_responsibility = dbeta_un / responsibility_denom
  no_responsibility = dbeta_no / responsibility_denom
  ov_responsibility = dbeta_ov / responsibility_denom
  ob_responsibility = dbeta_ob / responsibility_denom
  # na_responsibility = dbeta_na / responsibility_denom
    
  mixture_param1 = sum(un_responsibility)/length(data)
  mixture_param2 = sum(no_responsibility)/length(data)
  mixture_param3 = sum(ov_responsibility)/length(data)
  mixture_param4 = sum(ob_responsibility)/length(data)
  # mixture_param5 = sum(na_responsibility)/length(data)
  mixture_params = c(mixture_param1, mixture_param2, mixture_param3,
                     mixture_param4)#, mixture_param5
  
  # if some weight is too low then break
  # if(any(mixture_params < 0.05)) break
  
  # print("Here 1")
  
  # compute the common area
  integral = integrate(area, id_1 = as.numeric(args[1]),
                       mixtures = mixture_params,
                       alphas = component_parameters[seq(1, length(component_parameters), by = 2)],
                       betas = component_parameters[seq(2, length(component_parameters), by = 2)],
                       lower = 0, upper = 1, subdivisions = 10000)$value
  
  print(paste0("Common area: ", integral))
  
  res_it = data.frame(id = as.numeric(args[1]),
                      Iteration = itmax,
                      Loglik = max_lik,
                      area = integral,
                      Tol = ifelse(is.infinite(tol), 0, tol),
                      Method = method,
                      Time = t_it,
                      w_inverse = grid_2$w_inverse[as.numeric(args[1])],
                      w_un = mixture_params[1],
                      alpha_un = component_parameters[1],
                      beta_un = component_parameters[2],
                      w_no = mixture_params[2],
                      alpha_no = component_parameters[3],
                      beta_no = component_parameters[4],
                      w_ov = mixture_params[3],
                      alpha_ov = component_parameters[5],
                      beta_ov = component_parameters[6],
                      w_ob = mixture_params[4],
                      alpha_ob = component_parameters[7],
                      beta_ob = component_parameters[8])
  
  res = rbind(res, res_it)
  # print("Here 2")
  itmax = itmax + 1
  }

# take the optimized values from the EM algorithm
em_params = data.frame(Sex = grid_2$sex[as.numeric(args[1])],
                       Transition = grid_2$tran[as.numeric(args[1])],
                       NS = rep(c("Underweight", "Normal",
                                  "Overweight", "Obese"), 2),#, "NA"), 2),
                       param = rep(c("alpha", "beta"), each = 4),#5),
                       value = component_parameters) %>%
  spread(param, value) %>%
  mutate(weight = mixture_params,
         id = as.numeric(args[1]))

# this part saves the files in NLHPC
saveRDS(em_params, paste0("Experiments/Par_",
                          grid_2$sex[as.numeric(args[1])], "_",
                          grid_2$tran[as.numeric(args[1])], "_",
                          as.numeric(args[1]), ".rds"))

saveRDS(res, paste0("Experiments/Res_",
                    grid_2$sex[as.numeric(args[1])], "_",
                    grid_2$tran[as.numeric(args[1])], "_",
                    as.numeric(args[1]), ".rds"))
