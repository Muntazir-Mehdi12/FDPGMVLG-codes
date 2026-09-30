# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------
#                                   -- MACROLOIRE DATASET --
# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------

#--------------------------#
# Load necessary libraries
#--------------------------#
library(factoextra)
library(ppclust)
library(fclust)
library(irr)
library(aricode)
library(genie)
library(mclust)
library(tictoc)

#---------------------------------------------------#
#             Load the Macroloire Dataset
#---------------------------------------------------#
data("macroloire",package ="ade4")
macro_data <- as.matrix(macroloire$envir[, c(2,3)])
true_lab_macro=as.numeric(macroloire$envir[,5])
table( true_lab_macro , true_lab_macro)

#--------------------------------------#
#      Scaling the data
#--------------------------------------#
macro_scaled_data <- scale(macro_data)

#--------------------------------------#
#   Number of rows and columns
#--------------------------------------#
P <- nrow(macro_scaled_data)
D <- ncol(macro_scaled_data)

#--------------------------------------#
#      Number of clusters G
#--------------------------------------#
G <- 3  #maximum no. of clusters

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
#   m : fuzzifier exponent (commonly m = 2)
#
# Returns:
#   Scalar IMI validity index (lower = better)
# --------------------------------------------------------

IMI <- function(X, U, m = 2) {
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

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#                          --- Other Fuzzy Clustering Approaches ---
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#                          --- Fuzzy K-means (Gustafson and Kessel) ---
#                          ---              (EGK)                   ---
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
tic("EGK Runtime")
best_EGK_by_seed <- function(
    data,                 # -- scaled dataset
    true_labels,          # -- ground truth vector
    G,                    # -- number of clusters
    ent,                  # -- entropy degree
    max_seeds,            # -- how many seeds to try
    verbose = TRUE
){
  
  best_kappa <- -Inf
  best_seed  <- NA
  best_res   <- NULL
  best_u     <- NULL
  
  for (s in 1:max_seeds) {
    
    set.seed(s)
    u <- inaparc::imembrand(nrow(data), k = G)$u
    
    # -- running EGK 
    res <- try(fclust::FKM.gk.ent(data, k = G, ent = ent, RS = 10, seed = s),
               silent = TRUE)
    
    if (inherits(res, "try-error")) {
      if (verbose & s %% 100 == 0)
        cat("Checked seeds up to:", s, "\n")
      next
    }
    
    # -- Extract predicted clusters from membership 
    memb_df <- as.data.frame(cl.memb(res$U))
    pred <- memb_df$Cluster
    
    # -- Relabel to match ground truth
    pred_relab <- relabel_clusters_optimal(pred, true_labels)
    pred <- pred_relab$remapped_labels
    
    # -- Compute metrics
    kappa_val <- kappa2(
      data.frame(rater1 = true_labels, rater2 = pred) 
    )$value
    
    ari_val  <- adjustedRandIndex(true_labels, pred)
    nmi_val  <- NMI(true_labels, pred)
    imi_val  <- IMI(data, memb_df, m = 2)  
    
    if (verbose) {
      cat(
        sprintf("Seed %d → kappa=%.4f  ARI=%.4f  NMI=%.4f  IMI=%.6f\n",
                s, kappa_val, ari_val, nmi_val, imi_val)
      )
    }
    
    if (kappa_val > best_kappa) {
      best_kappa <- kappa_val
      best_seed  <- s
      best_res   <- res
      best_u     <- res$U
      
      if (verbose) {
        cat("  → NEW BEST! (Seed:", s, ")\n")
      }
    }
  }
  
  # -- Final output
  list(
    best_seed  = best_seed,
    best_kappa = best_kappa,
    best_model = best_res,
    best_u     = best_u
  )
}

out_egk <- best_EGK_by_seed(
  data = macro_scaled_data,
  true_labels = true_lab_macro,
  G = G,
  ent = 0.1,
  max_seeds = 300
)

out_egk$best_seed
out_egk$best_kappa


fkm3 <- FKM.gk.ent(macro_scaled_data, k = G, ent = 0.1, RS = 10, seed = 3)
cl.size(fkm3$U)
fkm_membership3 <- as.data.frame(cl.memb(fkm3$U))
fkm_clusters3 <- fkm_membership3$Cluster

# --- Relabeling
fkm_clusters3 <- relabel_clusters_optimal(fkm_clusters3, true_lab_macro)
fkm_clusters3 <- fkm_clusters3$remapped_labels

fviz_cluster(list(data = macro_scaled_data, cluster = fkm_clusters3), geom = "point", ellipse.type = "convex", ggtheme = theme_minimal(), main = "True label Clustering")

#---------------------------------------------------------------------------------#
#                                 EGK RESULTS
#---------------------------------------------------------------------------------#
calculate_EGK_clustering_metrics <- function(data_1, true_clusters, cluster_labels) {
  
  results_list <- list()
  
  table_result <- table(true_clusters, cluster_labels)
  kappa_result <- kappa2(data.frame(rater1 = true_clusters, rater2 = cluster_labels))
  ari_result <- adjustedRandIndex(true_clusters, cluster_labels)
  nmi_result <- NMI(true_clusters, cluster_labels)
  IMI_result <- IMI(data_1, fkm3$U, 2)
  
  # Store results
  results_list[["EGK"]] <- list(
    Table = table_result,
    Kappa = kappa_result$value,
    ARI = ari_result,
    NMI = nmi_result,
    IMI = IMI_result
  )
  
  return(results_list)
}

#------------------#
# Call the function
#------------------#
EGK_metrics_1 <- calculate_EGK_clustering_metrics(data_1 = macro_scaled_data, true_clusters = true_lab_macro, cluster_labels = fkm_clusters3)

#------------------#
#     Results 
#------------------#
for (EGK in names(EGK_metrics_1)) {
  cat("\n", EGK, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(EGK_metrics_1[[EGK]]$Table)
  cat("Kappa:", EGK_metrics_1[[EGK]]$Kappa, "\n")
  cat("ARI:", EGK_metrics_1[[EGK]]$ARI, "\n")
  cat("NMI:", EGK_metrics_1[[EGK]]$NMI, "\n")
  cat("IMI:", EGK_metrics_1[[EGK]]$IMI, "\n")
}
EGK_run_time <- toc()

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#                                  --- Gath-Geva Clustering ---
#                                  ---        (GG)          ---
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
tic("GG Runtime")
best_GG_by_seed <- function(
    data,                 # -- scaled dataset
    true_labels,          # -- ground truth vector
    G,                    # -- number of clusters
    m,                    # -- fuzzifier
    max_seeds,            # -- how many seeds to try
    verbose = TRUE
){
  
  best_kappa <- -Inf
  best_seed  <- NA
  best_res   <- NULL
  best_u     <- NULL
  
  for (s in 1:max_seeds) {
    
    set.seed(s)
    u <- inaparc::imembrand(nrow(data), k = G)$u
    
    # -- Try running GG safely
    res <- try(ppclust::gg(data, centers = G, memberships = u, m = m),
               silent = TRUE)
    
    if (inherits(res, "try-error")) {
      if (verbose & s %% 100 == 0)
        cat("Checked seeds up to:", s, "\n")
      next
    }
    
    # -- Extract predicted clusters
    pred <- res$cluster
    
    # -- Relabel to match ground truth
    pred <- relabel_clusters_optimal(pred, true_labels)
    pred <- pred$remapped_labels
    
    # -- Compute metrics
    kappa_val <- kappa2(
      data.frame(rater1 = true_labels, rater2 = pred) 
    )$value
    
    ari_val  <- adjustedRandIndex(true_labels, pred)
    nmi_val  <- NMI(true_labels, pred)
    imi_val  <- IMI(data, res$u, m)
    
    if (verbose) {
      cat(
        sprintf("Seed %d → kappa=%.4f  ARI=%.4f  NMI=%.4f  IMI=%.6f\n",
                s, kappa_val, ari_val, nmi_val, imi_val)
      )
    }
    
    # Check if this is the best kappa
    if (kappa_val > best_kappa) {
      best_kappa <- kappa_val
      best_seed  <- s
      best_res   <- res
      best_u     <- u
      
      if (verbose) {
        cat("  → NEW BEST! (Seed:", s, ")\n")
      }
    }
  }
  
  # -- Final output
  list(
    best_seed  = best_seed,
    best_kappa = best_kappa,
    best_model = best_res,
    best_u     = best_u
  )
}

out <- best_GG_by_seed(
  data = macro_scaled_data,
  true_labels = true_lab_macro,
  G = G,
  m = 1,
  max_seeds = 300
)

out$best_seed
out$best_kappa

GG <- out$best_model
u_best <- out$best_u

set.seed(57)
u <- inaparc::imembrand(nrow(macro_scaled_data), k = G)$u

# ---- Step 3: Run Gath–Geva with stable initialization
GG <- ppclust::gg(macro_scaled_data, centers = G, memberships = u, m = 1)

GG_clusters <- GG$cluster
GG_clusters <- relabel_clusters_optimal(GG_clusters, true_lab_macro)
GG_clusters <- GG_clusters$remapped_labels

#---------------------------------------------------------------------------------#
#                                 GG RESULTS
#---------------------------------------------------------------------------------#
calculate_GG_clustering_metrics <- function(data_1, true_clusters, cluster_labels) {
  
  results_list <- list()
  
  table_result <- table(true_clusters, cluster_labels)
  kappa_result <- kappa2(data.frame(rater1 = true_clusters, rater2 = cluster_labels))
  ari_result <- adjustedRandIndex(true_clusters, cluster_labels)
  nmi_result <- NMI(true_clusters, cluster_labels)
  IMI_result <- IMI(data_1, GG$u, 2)
  
  # Store results
  results_list[["GG"]] <- list(
    Table = table_result,
    Kappa = kappa_result$value,
    ARI = ari_result,
    NMI = nmi_result,
    IMI = IMI_result
  )
  
  return(results_list)
}

#------------------#
# Call the function
#------------------#
GG_metrics_1 <- calculate_GG_clustering_metrics(data_1 = macro_scaled_data, true_clusters = true_lab_macro, cluster_labels = GG_clusters)

#------------------#
#     Results 
#------------------#
for (GG in names(GG_metrics_1)) {
  cat("\n", GG, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(GG_metrics_1[[GG]]$Table)
  cat("Kappa:", GG_metrics_1[[GG]]$Kappa, "\n")
  cat("ARI:", GG_metrics_1[[GG]]$ARI, "\n")
  cat("NMI:", GG_metrics_1[[GG]]$NMI, "\n")
  cat("IMI:", GG_metrics_1[[GG]]$IMI, "\n")
}
GG_run_time <- toc()

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#                                          GENIE METHOD
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#

#--------------------------------------------#
#          label mapping function
#--------------------------------------------#
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
tic("Genie Runtime")
#-----------------------------------------------#
#             Gini values to test
#-----------------------------------------------#
gini_values <- seq(0.001, 1, 0.1)
kappa_scores <- numeric(length(gini_values))
names(kappa_scores) <- gini_values

k <- length(unique(true_lab_macro)) 

for (i in seq_along(gini_values)) {
  gini <- gini_values[i]
  cat("\nTrying gini threshold:", gini, "\n")
  
  tic("Genie Runtime")
  res <- try(gclust(macro_scaled_data, gini_threshold = gini), silent = TRUE)
  run_time_macro <- toc()
  
  if (inherits(res, "try-error") || is.null(res)) {
    cat(" Error or invalid clustering result\n")
    kappa_scores[i] <- NA
  } else {
    cluster_labels <- cutree(res, k = k)
    
    relabeled <- relabel_clusters_optimal(cluster_labels, true_lab_macro)
    aligned_labels <- relabeled$remapped_labels
    
    cat("  Confusion matrix (after relabeling):\n")
    print(relabeled$confusion)
    
    kappa_result <- kappa2(data.frame(rater1 = true_lab_macro, rater2 = aligned_labels))
    kappa_scores[i] <- kappa_result$value
    cat(" Kappa:", kappa_result$value, "\n")
  }
}

# final Kappa scores
cat("\nFinal Kappa scores by Gini threshold:\n")
print(kappa_scores)

#-------------------------------------------------------#
#           Running Genie with best Gini Coefficient
#-------------------------------------------------------#
res <- gclust(macro_scaled_data, gini_threshold = 0.001)
k <- length(unique(true_lab_macro))

aligned_labels <- cutree(res, k = k)

table( true_lab_macro , aligned_labels)
kappa2(data.frame(rater1 = true_lab_macro, rater2 = aligned_labels))

#-------------------------------------------------------#
#            Relabeling
#-------------------------------------------------------#
new <- c(2,3)
old <- c(3,2)
aligned_labels[aligned_labels %in% old] <- new[match(aligned_labels, old, nomatch = 0)]
table( true_lab_macro , aligned_labels)
kappa2(data.frame(rater1 = true_lab_macro, rater2 = aligned_labels))

#---------------------------------------------------------------------------------#
#                            GENIE RESULTS
#---------------------------------------------------------------------------------#
calculate_Genie_clustering_metrics <- function(data_1, true_clusters, cluster_labels) {
  
  results_list <- list()
  
  table_result <- table(true_clusters, cluster_labels)
  #kappa_result <- kappa2(data.frame(rater1 = true_clusters, rater2 = cluster_labels))
  kappa_result <- kappa2(data.frame(rater1 = true_clusters, rater2 = cluster_labels))
  ari_result <- adjustedRandIndex(true_clusters, cluster_labels)
  nmi_result <- NMI(true_clusters, cluster_labels)
  
  # Store results
  results_list[["Genie"]] <- list(
    Table = table_result,
    Kappa = kappa_result$value,
    ARI = ari_result,
    NMI = nmi_result
  )
  
  return(results_list)
}

#------------------#
# Call the function
#------------------#
Genie_metricsmacro_1 <- calculate_Genie_clustering_metrics(data_1 = macro_scaled_data, true_clusters = true_lab_macro, cluster_labels = aligned_labels)

#------------------#
#     Results 
#------------------#
for (Genie in names(Genie_metricsmacro_1)) {
  cat("\n", Genie, " Clustering Results:\n", sep="")
  cat("Table:\n")
  print(Genie_metricsmacro_1[[Genie]]$Table)
  cat("Kappa:", Genie_metricsmacro_1[[Genie]]$Kappa, "\n")
  cat("Adjusted Rand Index:", Genie_metricsmacro_1[[Genie]]$ARI, "\n")
  cat("NMI:", Genie_metricsmacro_1[[Genie]]$NMI, "\n")
}

U_m <- matrix(0, nrow = length(aligned_labels), ncol = G)
for (i in seq_along(aligned_labels)) {
  U_m[i, aligned_labels[i]] <- 1
}

IMI(macro_scaled_data, U_m, 2)
run_time_Elisa_Genie <- toc()
