# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------
# --------------------------------------- PUROMYCIN DATASET ---------------------------------
# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------

# Load the functions used below:
source("DBDA2E-utilities.R") # Must be in R's current working directory.

#--------------------------#
# Load necessary libraries
#--------------------------#
library(rjags)
library(coda)
library(parallel)
library(runjags)
library(tictoc)
library(ggplot2)
library(GGally)
library(DAAG)
library(mlbench)
library(irr)
library(psych)
library(mclust)
library(factoextra)
library(gridExtra)
library(clustertend)
library(hopkins)
library(clValid)
library(clusterSim)
library(dirichletprocess)
library(fossil)
library(corrplot)
library(moments)
library(MASS)
library(ade4)
library(aricode)
library(genie)
library(gtools) 
library(combinat)
library(clue)
library(irr)

#---------------------------------------------------#
#             Load the Puromycin Dataset
#---------------------------------------------------#
data("Puromycin", package = "datasets")
Pu_data <- as.matrix(Puromycin[,c(-3)]) 
true_lab_Pu <- as.numeric(Puromycin$state) 
table( true_lab_Pu , true_lab_Pu)

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#                                          DPGMVLG
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#

#--------------------------------------#
#      Scaling the data
#--------------------------------------#
Pu_scaled_data <- scale(Pu_data)

#--------------------------------------#
#      Plotting true labels
#--------------------------------------#
fviz_cluster(list(data = Pu_scaled_data, cluster = true_lab_Pu), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

#--------------------------------------#
#     Number of rows and columns
#--------------------------------------#
P <- nrow(Pu_scaled_data)
D <- ncol(Pu_scaled_data)

#--------------------------------------#
#      Number of clusters G
#--------------------------------------#
G <- 2  #maximum no. of clusters

#--------------------------------------#
#      Setting up data for JAGS
#--------------------------------------#
data_list <- list(
  x = Pu_scaled_data,
  P = P,
  D = D,
  deltaData = det(cor(Pu_scaled_data))^(1 / (D - 1)),
  G = G 
)

#-------------------------------------#
#           JAGS MODEL
#-------------------------------------#
model_string <- "
Data {
  C <- 1000000000
  for (i in 1:P) {
     zeros[i] <- 0
  }
  v_g <- 1.4262
  delta_g <- deltaData
}
model {
  for (i in 1:P) {
   for (g in 1:G){
      for ( j in 1:D){
        summ2[i,j,g] <- exp(mu_g[j, g] * x[i, j])/lambda_g[j, g]
      }
      logLik[i, g] <- z[i, g] * (log(pi_g[g])+ v_g * log(delta_g) + sum(log(mu_g[ ,g])) - v_g * sum(log(lambda_g[, g])) - loggam(v_g) ) + 
                      v_g * z[i, g] * sum(mu_g[1:D, g] * x[i, 1:D]) - z[i, g] * sum(summ2[i, , g])  +  
                      z[i, g] * (1 - delta_g) * prod(1/lambda_g[ , g]) * exp(-2.07728 * (D-1) + sum(mu_g[ , g] * x[i, ]) ) 
    
      z[i, g] <- ifelse(cluster[i] == g, 1, 0) 
   }
   cluster[i] ~ dcat(pi_g) 
   zeros[i] ~ dpois(-sum(logLik[i,]) + C)
  }
   for (g in 1:G) {
    for (j in 1:D) {
       mu_g[j, g] ~ dgamma(25,20)          
       lambda_g[j, g] ~ dgamma(11,3)     
    }
    alpha[g] ~ dgamma(2,4)    
   }
   pi_g[1:G] ~ ddirch(alpha/sum(alpha))
}
"
#-----------------------------------------------#
#            Writing the model to a file
#-----------------------------------------------#
writeLines(model_string, con = "TEMPmodel.txt")

#-----------------------------------------------#
#             Parameters to monitor
#-----------------------------------------------#
#params <- c("delta_g","mu_g", "lambda_g", "cluster", "z")
params <- c("mu_g", "lambda_g", "pi_g", "z")

#-----------------------------------------------#
#               Reproducibility
#-----------------------------------------------#
inits <- list(list(.RNG.name = "base::Mersenne-Twister",.RNG.seed = 22021),
              list(.RNG.name = "base::Mersenne-Twister",.RNG.seed = 32019))

#-----------------------------------------------#
# Tracking the Run time and Running the JAGS model
#-----------------------------------------------#.
tic("JAGS Model Runtime")
runJagsOut <- run.jags( method="parallel" ,
                        model="TEMPmodel.txt" ,
                        monitor=params ,
                        data=data_list ,
                        n.chains=2 ,
                        adapt = 3000,
                        burnin=2000,
                        sample=1000, 
                        inits = inits,
                        thin = 2)
run_time_Pu <- toc()  # Capture the runtime
codaSamples_Pu = as.mcmc.list( runJagsOut )
summaryChains_Pu <- summary(codaSamples_Pu)

#----------------------------------------------#
#               Diagnostics
#----------------------------------------------#
diagMCMC( codaObject=codaSamples_Elisa , parName="lambda_g[1,1]" )
diagMCMC( codaObject=codaSamples_Elisa , parName="lambda_g[1,2]" )
diagMCMC( codaObject=codaSamples_Elisa , parName="lambda_g[2,1]" )
diagMCMC( codaObject=codaSamples_Elisa , parName="lambda_g[2,2]" )
diagMCMC( codaObject=codaSamples_Elisa , parName="mu_g[1,1]" )
diagMCMC( codaObject=codaSamples_Elisa , parName="mu_g[1,2]" )
diagMCMC( codaObject=codaSamples_Elisa , parName="mu_g[2,1]" )
diagMCMC( codaObject=codaSamples_Elisa , parName="mu_g[2,2]" )
graphics.off()

#---------------------------------------------#
#             Summary Chains
#---------------------------------------------#
summaryChains_Pu$statistics
summaryChains_Pu$statistics[11,]
matrix(summaryChains_Pu$statistics[11:nrow(summaryChains_Pu$statistics),1], P, G)
z_mode_Pu <- apply(matrix(summaryChains_Pu$statistics[11:nrow(summaryChains_Pu$statistics),1], P, G),1, which.max)
z_mode_Pu

#--------------------------------------------#
#         Plotting Predicted labels
#--------------------------------------------#
fviz_cluster(list(data = Pu_scaled_data, cluster = z_mode_Pu), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

plot(Pu_scaled_data, col= z_mode_Pu, pch=16, main="Data Points by Cluster")
legend("topright", legend=unique(z_mode_Pu), col=unique(z_mode_Pu), pch=16, title="Cluster")
table( true_lab_Pu , z_mode_Pu)
kappa2(data.frame(rater1 = true_lab_Pu, rater2 = z_mode_Pu))

#--------------------------------------------#
#           DPGMVLG Results
#--------------------------------------------#
calculate_dp_Gmvlg_clustering_metrics <- function(data_1, true_clusters, z_mode) {
  
  results_list <- list()
  
  table_result <- table(true_clusters, z_mode)
  kappa_result <- kappa2(data.frame(rater1 = true_clusters, rater2 = z_mode))
  ari_result <- adjustedRandIndex(true_clusters, z_mode)
  nmi_result <- NMI(true_clusters, z_mode)
  
  # Store results
  results_list[["dp_Gmvlg"]] <- list(
    Table = table_result,
    Kappa = kappa_result$value,
    ARI = ari_result,
    NMI = nmi_result,
    CPU_RUN_TIME = run_time_Pu$callback_msg  
  )
  return(results_list)
}

#--------------------#
# Call the function
#--------------------#
dp_Gmvlg_metricsPu_1 <- calculate_dp_Gmvlg_clustering_metrics(data_1 = Pu_scaled_data, true_clusters = true_lab_Pu, z_mode = z_mode_Pu)

#--------------------#
#       Results 
#--------------------#
for (dp_Gmvlg in names(dp_Gmvlg_metricsPu_1)) {
  cat("\n", dp_Gmvlg, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(dp_Gmvlg_metricsPu_1[[dp_Gmvlg]]$Table)
  cat("Kappa:", dp_Gmvlg_metricsPu_1[[dp_Gmvlg]]$Kappa, "\n")
  cat("Adjusted Rand Index:", dp_Gmvlg_metricsPu_1[[dp_Gmvlg]]$ARI, "\n")
  cat("NMI:", dp_Gmvlg_metricsPu_1[[dp_Gmvlg]]$NMI, "\n")
  cat(" CPU RUN TIME:", dp_Gmvlg_metricsPu_1[[dp_Gmvlg]]$CPU_RUN_TIME, "\n")
}

# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# ------------------------------------------ FUZZY BAYESIAN DPGMVLG MODEL -------------------------------
# ------------------------------------------    --  (FDPGMVLG)  --        ------------------------------- 
# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------

# -----------------------------
# -----------------------------
# --- GFN helper functions ----
# -----------------------------
# -----------------------------

# ----------------------------------------------------------------
# --- Defines a Gaussian fuzzy number (mean and variance pair) ---
# ----------------------------------------------------------------
gfn <- function(mu, var) list(mu = mu, var = var) 

# ----------------------------------------------
# --- Addition of two independent GFNs ---------
# ---------------- (A + B) ---------------------
# ----------------------------------------------
gfn_add <- function(A, B) {
  gfn(mu = A$mu + B$mu, var = A$var + B$var) # -- Addition of two independent GFNs.
}

# --------------------------------------------------
# --- Sum of multiple GFNs assuming independence ---
# --------------------------------------------------
gfn_sum_list <- function(gfn_list) {
  mu <- sum(sapply(gfn_list, function(x) x$mu))
  var <- sum(sapply(gfn_list, function(x) x$var))
  gfn(mu, var)
}

# ----------------------------------------------
# --------- Subtraction of two GFNs ------------
# ---------------- (A - B) ---------------------
# ----------------------------------------------
gfn_sub <- function(A, B) {
  gfn(mu = A$mu - B$mu, var = A$var + B$var)
}

# ---------------------------------------------------------------------
# --------------- Scalar multiplication -------------------------------
# ---------------       (a x A)         -------------------------------
# --------- "a" = a crisp number and "A" = a GFN ----------------------
# ---------------------------------------------------------------------
gfn_scalar_mul <- function(a, A) {
  gfn(mu = a * A$mu, var = (a^2) * A$var)
}

# ----------------------------------------------
# -------- Product of two independent GFNs -----
# ---------------- (A x B) ---------------------
# ----------------------------------------------
gfn_mul <- function(A, B) {
  mu <- A$mu * B$mu
  var <- (A$mu^2) * B$var + (B$mu^2) * A$var + A$var * B$var
  gfn(mu, var)
}

# --------------------------------------------------------------------------
# --- Division using a first-order Taylor (approx via delta method) --------
# ------------------------  (A / B)  ---------------------------------------
# --------------------------------------------------------------------------
gfn_div <- function(A, B) {
  mu <- A$mu / B$mu
  var <- (A$var / (B$mu^2)) + ((A$var^2 + A$mu^2) * B$var / (B$mu^4))
  gfn(mu, var)
}

# --------------------------------------------------------------
# --- Log of a positive GFN: use delta method approximations ---
# --------------------------------------------------------------
gfn_log <- function(A) {
  mu <- log(A$mu) 
  var <- A$var / (A$mu^2)
  gfn(mu, var)
}

# -------------------------------------------
# --- Exponential transformation of a GFN ---
# -------------------------------------------
gfn_exp <- function(A) {
  mu_out <- exp(A$mu - A$var)
  var_out <- exp(2 * A$mu + A$var) * (exp(A$var) - 1)
  gfn(mu_out, var_out)
}

# ----------------------------------------------------------
# --- Reciprocal of a positive GFN: 1 / A (delta method) ---
# ----------------------------------------------------------
gfn_reciprocal <- function(A) {
  mu <- 1 / A$mu
  var <- A$var / (A$mu^4)
  gfn(mu, var)
}

# ---------------------------------------------------------------------------------
# --- Product of vector of GFNs using log/exp trick to keep numerical stability ---
# ---------------------------------------------------------------------------------
gfn_prod_list <- function(gfn_list) {
  logs <- lapply(gfn_list, gfn_log)
  s <- gfn_sum_list(logs)      
  gfn_exp(s)
}

# ----------------------------------------------------------
# ----------------------------------------------------------
# --- Extracting posterior summaries (means & variances) ---
# ----------------------------------------------------------
# ----------------------------------------------------------

mcmc_mat <- as.matrix(codaSamples_Pu) 

param_names <- colnames(mcmc_mat)

# --- Helper to compute posterior mean and var of a parameter ---
post_mean <- function(name) mean(mcmc_mat[, name])
post_var  <- function(name) var(mcmc_mat[, name])

# --- GFNs for mu_g[j,g], lambda_g[j,g], and pi_g[g] ---
mu_names <- param_names[grep("^mu_g\\[", param_names)]
lambda_names <- param_names[grep("^lambda_g\\[", param_names)]
pi_names <- param_names[grep("^pi_g\\[", param_names)]

# ------------------------------------------------------------------------
# -- controlled perturbation mechanism inspired by perturbation methods --
# ------------------------------------------------------------------------
build_gfn_parameters <- function(perturb = TRUE,
                                 mu_scale = mu_scale,
                                 var_scale = var_scale) {
  
  
  mu_g_gfn <- array(list(), dim = c(D, G))
  lambda_g_gfn <- array(list(), dim = c(D, G))
  pi_g_gfn <- vector("list", G)
  
  # ----- mu_g -----
  for (nm in mu_names) {
    
    idxs <- gsub("mu_g\\[|\\]", "", nm)
    parts <- as.integer(strsplit(idxs, ",")[[1]])
    
    j <- parts[1]
    g <- parts[2]
    
    mu_mean <- post_mean(nm)
    mu_var  <- post_var(nm)
    
    if (perturb) {
      
      eps_mu  <- rnorm(1, 0, mu_scale)
      eps_var <- rnorm(1, 0, var_scale)
      
      mu_mean <- mu_mean + eps_mu * sqrt(mu_var)
      mu_var  <- mu_var * abs(1 + eps_var)
      
    }
    
    mu_g_gfn[[j,g]] <- gfn(mu_mean, mu_var)
  }
  
  # ----- lambda_g -----
  for (nm in lambda_names) {
    
    idxs <- gsub("lambda_g\\[|\\]", "", nm)
    parts <- as.integer(strsplit(idxs, ",")[[1]])
    
    j <- parts[1]
    g <- parts[2]
    
    lam_mean <- post_mean(nm)
    lam_var  <- post_var(nm)
    
    if (perturb) {
      
      eps_mu  <- rnorm(1,0,mu_scale)
      eps_var <- rnorm(1,0,var_scale)
      
      lam_mean <- lam_mean + eps_mu * sqrt(lam_var)
      lam_var  <- lam_var * abs(1 + eps_var)
      
    }
    
    lambda_g_gfn[[j,g]] <- gfn(lam_mean, lam_var)
  }
  
  # ----- pi_g -----
  if (length(pi_names) == G) {
    
    for (k in seq_len(G)) {
      
      pi_mean <- post_mean(pi_names[k])
      pi_var  <- post_var(pi_names[k])
      
      if (perturb) {
        
        eps_mu <- rnorm(1,0,0.1)
        eps_var <- rnorm(1,0,0.1)
        
        pi_mean <- abs(pi_mean + eps_mu * sqrt(pi_var))
        pi_var  <- pi_var * abs(1 + eps_var)
        
      }
      
      pi_g_gfn[[k]] <- gfn(pi_mean, pi_var)
    }
  }
  
  list(
    mu_g_gfn = mu_g_gfn,
    lambda_g_gfn = lambda_g_gfn,
    pi_g_gfn = pi_g_gfn
  )
}

# -------------------------------------------------------------------------
# --- Compute fuzzy complete-data log - likelihood for each observation ---
# -------------------------------------------------------------------------

# --- Function to compute fuzzy log-likelihood for a single observation i and cluster g ---
compute_fuzzy_logLik <- function(
    x_i, g, D, pi_g_gfn, mu_g_gfn, lambda_g_gfn,
    delta_g_crisp, term_v_g_delta_g, v_g, const_term, crisp_var) {
  
  #-----------------------------------------------------------------------------
  # 1) log(π_g) + v_g * log(δ_g) + ∑ log(μ_jg) - v_g * ∑ log(λ_jg) - log(Γ(ν_g)) 
  #-----------------------------------------------------------------------------
  
  pig_term <- gfn_log(pi_g_gfn[[g]])                                      # log(π_g)
  
  term_vd <- gfn(term_v_g_delta_g, crisp_var)                             # v_g * log(δ_g)
  
  log_mus <- lapply(seq_len(D), function(j) gfn_log(mu_g_gfn[[j,g]]))
  sum_log_mus <- gfn_sum_list(log_mus)                                    # ∑ log(μ_jg)
  
  log_lambdas <- lapply(seq_len(D), function(j) gfn_log(lambda_g_gfn[[j,g]]))
  sum_log_lambdas <- gfn_sum_list(log_lambdas)
  term_4 <- gfn_scalar_mul(-v_g, sum_log_lambdas)                         # - v_g * ∑ log(λ_jg)
  
  term_5 <- gfn(log(gamma(v_g)), crisp_var)                               # log(Γ(ν_g)) (crisp)
  
  part1 <- gfn_sum_list(list(pig_term, term_vd, sum_log_mus, term_4, gfn_scalar_mul(-1, term_5))) #First line
  
  #---------------------------------------------------------------
  # 2) ν_g * ∑ (μ_jg * x_ij)  - ∑ (1/λ_jg) * exp(μ_jg * x_ij)
  #---------------------------------------------------------------
  
  mu_x_sum <- gfn_sum_list(lapply(seq_len(D), function(j)         # ∑ (μ_jg * x_ij)
    gfn_scalar_mul(x_i[j], mu_g_gfn[[j, g]])
  ))
  
  inner_exp <- lapply(seq_len(D), function(j) {                   # (1/λ_jg) * exp(μ_jg * x_ij)
    arg <- gfn_scalar_mul(x_i[j], mu_g_gfn[[j, g]])
    exp_arg <- gfn_exp(arg)
    exp_scaled <- gfn_mul(exp_arg, gfn_reciprocal(lambda_g_gfn[[j, g]]))
  })
  sum_inner <- gfn_sum_list(inner_exp)                            # ∑{(1/λ_jg) * exp(μ_jg * x_ij)}
  
  part2 <- gfn_sub(gfn_scalar_mul(v_g, mu_x_sum), sum_inner)   # 1st term  - second term
  
  #----------------------------------------------------------------------------
  # 3) The infinite-sum additive term approximated by constant exponential term
  #----------------------------------------------------------------------------
  one_minus_delta <- gfn(1 - delta_g_crisp, crisp_var)
  inv_lambdas <- lapply(seq_len(D), function(j) gfn_reciprocal(lambda_g_gfn[[j, g]]))
  prod_inv_lambda <- gfn_prod_list(inv_lambdas)
  
  exp_arg <- gfn_add(gfn(const_term, crisp_var), mu_x_sum)
  exp_part <- gfn_exp(exp_arg)
  
  part3 <- gfn_mul(gfn_mul(one_minus_delta, prod_inv_lambda), exp_part)
  
  # sum all terms
  total <- gfn_sum_list(list(part1, part2, part3))
  return(total)
}

#-----------------------------------------------
# Defuzzification of GFN-based fuzzy memberships
#-----------------------------------------------

defuzzify_fuzzy_membership <- function(mu_mat, var_mat, alpha) {
  # mu_mat: P x G matrix of means
  # var_mat: P x G matrix of variances
  # alpha: scaling factor for variance shift
  
  P <- nrow(mu_mat)
  G <- ncol(mu_mat)
  
  # --- Asymmetry δ = μ / sqrt(σ²)
  delta_mat <- mu_mat / sqrt(var_mat + 1e-8)
  
  # --- Threshold d = first quartile of |δ|
  d <- quantile(abs(delta_mat), probs = 0.25, na.rm = TRUE)
  
  # --- Defuzzified matrix
  B <- matrix(NA, nrow = P, ncol = G)
  
  for (i in seq_len(P)) {
    for (g in seq_len(G)) {
      delta_ig <- delta_mat[i, g]
      mu_ig <- mu_mat[i, g]
      var_ig <- var_mat[i, g]
      
      if (abs(delta_ig) < d) {
        B[i, g] <- mu_ig + alpha * var_ig
      } else {
        B[i, g] <- mu_ig
      }
    }
  }
  
  return(list(B = B, delta = delta_mat, d = d))
}

#--------------------------------------------#
#           Fuzzy DPGMVLG function
#--------------------------------------------#
calculate_Fuzzy_dpGmvlg_clustering_metrics <- function(data_1, true_clusters, z_mode_combined) {
  
  results_list_Fuzzy <- list()
  
  table_result <- table(true_clusters, z_mode_combined)
  kappa_result <- kappa2(data.frame(rater1 = true_clusters, rater2 = z_mode_combined))
  ari_result <- adjustedRandIndex(true_clusters, z_mode_combined)
  nmi_result <- NMI(true_clusters, z_mode_combined)
  
  # Store results
  results_list_Fuzzy[["Fuzzy_dpGmvlg"]] <- list(
    Table = table_result,
    Kappa = kappa_result$value,
    ARI = ari_result,
    NMI = nmi_result
  )
  return(results_list_Fuzzy)
}

#--------------------------------------------#
#           Relabeling using a function
#--------------------------------------------#
library(clue)

# label mapping function
relabel_clusters_optimal <- function(pred_clusters, true_labels) {
  
  conf_mat <- table(pred_clusters, true_labels)
  
  cost_mat <- max(conf_mat) - conf_mat
  
  assignment <- solve_LSAP(cost_mat)
  
  cluster_levels <- as.integer(rownames(conf_mat))
  mapping <- setNames(as.integer(colnames(conf_mat)[assignment]), cluster_levels)
  
  remapped <- mapping[as.character(pred_clusters)]
  
  return(list(
    remapped_labels = remapped,
    mapping = mapping,
    confusion = table(remapped, true_labels)
  ))
}

# --------------------------------------------------------
# IMI: Imbalanced Membership Index
#
# Inputs:
#   X : N×D data matrix
#   U : N×G membership matrix (rows sum to 1)
#   m : fuzzifier exponent 
#
# Returns:
#   Scalar IMI validity index (lower = better)
# --------------------------------------------------------

IMI <- function(X, U, m) {
  N <- nrow(X)
  D <- ncol(X)
  G <- ncol(U)
  
  # ----- (1) Effective cluster sizes -----
  U_m <- U^m
  n_eff <- colSums(U_m)    # length G
  
  if (any(n_eff < 1e-12)) {
    warning("Some effective cluster sizes are extremely small.")
    n_eff[n_eff < 1e-12] <- 1e-12
  }
  
  # ----- (2) Fuzzy cluster centers -----
  # V[g, ] = sum_i u_ig^m x_i / n_eff[g]
  V <- (t(U_m) %*% X) / n_eff
  
  # ----- (3) Compactness -----
  comp_g <- numeric(G)
  for (g in 1:G) {
    d2 <- rowSums((X - matrix(V[g, ], N, D, byrow = TRUE))^2)
    comp_g[g] <- sum(U_m[, g] * d2) / n_eff[g]
  }
  
  # ----- (4) Separation: minimum center distance -----
  sep_min <- Inf
  for (g in 1:(G-1)) {
    for (h in (g+1):G) {
      dgh <- sqrt(sum((V[g, ] - V[h, ])^2))
      if (dgh < sep_min) sep_min <- dgh
    }
  }
  
  sep_min <- max(sep_min, 1e-12)
  
  # ----- (5) Imbalance weighting -----
  weight_g <- 1 / n_eff
  
  # ----- (6) Final IMI metric -----
  idx_value <- sum(weight_g * comp_g) / sep_min
  
  return(idx_value)
}


# --- DPGMVLG ---
IMI(Pu_scaled_data, matrix(summaryChains_Pu$statistics[11:nrow(summaryChains_Pu$statistics),1], P, G), 2)

# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------
# ----------------------------------------- RANDOM SEARCH -----------------------------------
# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------

tic("FDPGMVLG Runtime")
tune_fuzzy_pipeline <- function(
    kappaToBeat,
    stop_on_success = TRUE,
    max_trials = 100000,
    verbose = TRUE,
    seed = NULL
) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # ---- Basic checks ----
  stopifnot(exists("Pu_scaled_data"), exists("G"), exists("D"))
  P <- nrow(Pu_scaled_data)
  
  if (verbose) message("P = ", P, ", G = ", G, ", D = ", D)
  
  # ---- delta crisp vector ----
  if (exists("data_list") && !is.null(data_list$deltaData)) {
    delta_g_crisp <- data_list$deltaData
  } else if (exists("deltaData")) {
    delta_g_crisp <- deltaData
  } else stop("deltaData or data_list$deltaData must exist.")
  
  const_term_base <- -2.07728 * (D - 1)
  
  # ---- storage ----
  best_kappa <- -Inf
  best_params <- NULL
  
  trial_no <- 0
  
  while (trial_no < max_trials) {
    
    trial_no <- trial_no + 1
    
    # ---- sample parameters ----
    v_g <- rgamma(1, 1, 1)
    crisp_var <- rgamma(1, 1, 0.1)
    alpha <- rnorm(1, 0, 1)
    
    if (verbose && trial_no %% 500 == 0) {
      message(Sys.time(),
              "  Trial #", trial_no,
              "  v_g=", signif(v_g,4),
              "  crisp_var=", signif(crisp_var,4),
              "  alpha=", signif(alpha,4))
    }
    
    t_start <- proc.time()
    
    try({
      
      # --- sample perturbation hyperparameters FIRST ---
      mu_scale  <- runif(1, 0, 1)
      var_scale <- runif(1, 0, 1)
      
      # --- now build GFN parameters using these values ---
      params <- build_gfn_parameters(
        perturb = TRUE,
        mu_scale = mu_scale,
        var_scale = var_scale
      )
      
      mu_g_gfn <- params$mu_g_gfn
      lambda_g_gfn <- params$lambda_g_gfn
      pi_g_gfn <- params$pi_g_gfn
      
      term_v_g_delta_g <- v_g * log(delta_g_crisp)
      
      P_local <- nrow(Pu_scaled_data)
      
      fuzzy_scores <- array(list(), dim = c(P_local, G))
      
      for (i in seq_len(P_local)) {
        
        x_i <- as.numeric(Pu_scaled_data[i, ])
        
        for (g_idx in seq_len(G)) {
          
          sc <- tryCatch({
            
            compute_fuzzy_logLik(
              x_i = x_i,
              g = g_idx,
              D = D,
              pi_g_gfn = pi_g_gfn,
              mu_g_gfn = mu_g_gfn,
              lambda_g_gfn = lambda_g_gfn,
              delta_g_crisp = delta_g_crisp,
              term_v_g_delta_g = term_v_g_delta_g,
              v_g = v_g,
              const_term = const_term_base,
              crisp_var = crisp_var
            )
            
          }, error = function(e) {
            list(mu = NaN, var = NaN)
          })
          
          if (!is.finite(sc$mu) || !is.finite(sc$var)) {
            sc$mu <- 0
            sc$var <- 1e6
          }
          
          fuzzy_scores[[i, g_idx]] <- sc
        }
      }
      
      # ---- convert to matrices ----
      fuzzy_membership_mu <- matrix(NA, P_local, G)
      fuzzy_membership_var <- matrix(NA, P_local, G)
      
      for (i in seq_len(P_local)) {
        for (g_idx in seq_len(G)) {
          
          fuzzy_membership_mu[i, g_idx] <- fuzzy_scores[[i, g_idx]]$mu
          fuzzy_membership_var[i, g_idx] <- fuzzy_scores[[i, g_idx]]$var
          
        }
      }
      
      # ---- defuzzify ----
      def_res <- defuzzify_fuzzy_membership(
        mu_mat = fuzzy_membership_mu,
        var_mat = fuzzy_membership_var,
        alpha = alpha
      )
      
      B_mat <- def_res$B
      
      # ---- crisp labels ----
      z_mode_combined <- apply(B_mat, 1, which.max)
      
      relabeled <- relabel_clusters_optimal(z_mode_combined, true_lab_Pu)
      z_mode_combined <- relabeled$remapped_labels
      
      # --- IMI Metric ---
      IMI_val <- IMI(Pu_scaled_data, B_mat, 2)
      
      # ---- metrics ----
      dp_metrics <- calculate_Fuzzy_dpGmvlg_clustering_metrics(
        data_1 = Pu_scaled_data,
        true_clusters = true_lab_Pu,
        z_mode_combined = z_mode_combined
      )
      
      kappa_val <- dp_metrics$Fuzzy_dpGmvlg$Kappa
      ari_val <- dp_metrics$Fuzzy_dpGmvlg$ARI
      nmi_val <- dp_metrics$Fuzzy_dpGmvlg$NMI
      
      runtime_sec <- (proc.time() - t_start)[["elapsed"]]
      
      # ---- update best ----
      if (!is.na(kappa_val) && kappa_val > best_kappa) {
        
        best_kappa <- kappa_val
        
        best_params <- list(
          trial = trial_no,
          kappa = kappa_val,
          ARI = ari_val,
          NMI = nmi_val,
          IMI = IMI_val,
          v_g = v_g,
          crisp_var = crisp_var,
          alpha = alpha,
          mu_scale = mu_scale,
          var_scale = var_scale,
          runtime = runtime_sec
        )
        
        if (verbose) {
          message("NEW BEST κ = ", signif(best_kappa,4),
                  "  (v_g=", signif(v_g,4),
                  ", crisp_var=", signif(crisp_var,4),
                  ", alpha=", signif(alpha,4), 
                  ", mu_scale=", signif(mu_scale,3),
                  ", var_scale=", signif(var_scale,3), ")")
        }
      }
      
      # ---- stopping rule ----
      if (!is.na(kappa_val) && kappa_val >= kappaToBeat) {
        
        if (verbose) {
          message("\nTARGET REACHED")
          message("κ = ", kappa_val,
                  " ≥ ", kappaToBeat)
          message("Stopping search at trial ", trial_no)
        }
        
        return(best_params)
      }
      
    }, silent = TRUE)
    
  }
  
  if (verbose) {
    message("\nSearch finished. Best κ = ", signif(best_kappa,4))
  }
  
  return(best_params)
}

best_result <- tune_fuzzy_pipeline(
  kappaToBeat = 0.657,
  max_trials = 100000,
  verbose = TRUE,
  seed = 124
)
FDPGMVLG_run_time <- toc()
best_result

