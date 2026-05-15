if (file.exists("code/00_packages.R")) {
  source("code/00_packages.R")
} else if (file.exists("00_packages.R")) {
  source("00_packages.R")
} else if (file.exists("replication_package/code/00_packages.R")) {
  source("replication_package/code/00_packages.R")
} else {
  stop("Could not locate replication_package bootstrap.")
}

# Convert AdjMats to Edgelists ------------------------------------------------------
mat_to_df = function(x, dims) {
  
  ## dims = villages to process
  y = eval(as.name(x))
  
  for (i in dims) {
    
    colnames(y[[i]]) = 1:nrow(y[[i]])
    rownames(y[[i]]) = 1:nrow(y[[i]])
    
    ind = which(upper.tri(y[[i]]), arr.ind = TRUE)
    nm = dimnames(y[[i]])
    
    y[[i]] = data.frame(id1 = nm[[1]][ind[, 1]],
                        id2 = nm[[2]][ind[, 2]],
                        link = y[[i]][ind],
                        village = i) %>% 
      arrange(as.numeric(id1), as.numeric(id2))
    
    names(y[[i]]) = c("id1", "id2", x, "village")
  }
  
  y = y[!sapply(y, is.null)]
  return(bind_rows(y))
}

# bootstrap function ----------------------------------------------------------------
Boot_eigen = function(dat, i) {
  X2 = as.data.frame(dat[i, ])
  X2 = as.matrix(X2)
  X2_cov = cov(X2)
  X2_eig = eigen(X2_cov)
  return(X2_eig)
}

# Regression model output -----------------------------------------------------------
mod_out = function(reg_input) {
  
  modelsummary(reg_input, statistic = c("std.error","p.value"),
               gof_omit = '[^Num.Obs.]',
               coef_omit = '[num_hh_random]')
  
}

# Giant Component Construction ------------------------------------------------------
giant_comp_const = function(net_type) {
  
  y = eval(as.name(net_type))
  gc = graph_from_adjacency_matrix(y[[i]], mode = "max")
  verts_gc = as.vector(
    V(gc)[igraph::components(gc)$membership] == which.max(igraph::components(gc)$csize)
  )
  gc = y[[i]][verts_gc, verts_gc]
  d = rowSums(gc)
  
  T_gc = matrix(data = NA, nrow = nrow(gc), ncol = nrow(gc))
  
  for (j in 1:nrow(T_gc)) {
    for (k in 1:nrow(T_gc)){
      
      if (j != k) {
        T_gc[j, k] = gc[j,k]/(d[j] + 1)
      } else {
        T_gc[j, k] = 1/(d[j] + 1)
      }
    }
  }
  
  return(Mod(sort(eigen(T_gc)$values, decreasing = TRUE)[2]))
}

# Diffusion Centrality (Only seeds (Random)) ----------------------------------------
DC_out_seed = function(net_type, k, seed_data, lambda_data) {

  y = eval(as.name(net_type))
  
  ## set q
  lambda_1 = lambda_data[lambda_data$village == k, grep(net_type, names(lambda_data))]
  q = 1/lambda_1
  
  ## set T
  g = graph_from_adjacency_matrix(y[[k]], mode = "max")

  diam = diameter(g, directed = FALSE)
  T_diam = diam

  DC = list()
  
  for (i in 0:T_diam) {
    DC[[i + 1]] = (q * y[[k]]) %^% i
  }
  
  DC = Reduce('+', DC) %*% matrix(rep(1, nrow(y[[k]])), ncol = 1)
  seeds = seed_data %>% filter(villageid == k) %>% select(vertex) %>% pull()
  DC = DC[seeds]
  
  return(sum(DC, na.rm = TRUE))
}

loading_plot = function(dat, labsize) { 
  out = factoextra::fviz_pca_var(dat, col.var = "red",
                                 axes = c(i, j),
                                 labelsize = labsize,
                                 arrowsize = 0.4,
                                 repel = TRUE) +
    theme(text = element_text(size = 7.5),
          axis.title = element_text(size = 7.5),
          axis.text = element_text(size = 7.5))
  
  return(out)
  
}

# BEV Function (part of laddle plot formation) --------------------------------------
bev_out = function(n_boots, X_mat, d) {
  
  boot_out = list()
  for (i in seq_len(n_boots)) {
    
    X2 = X_mat[sample(seq_len(nrow(X_mat)), nrow(X_mat), replace = TRUE), ]
    X2_cov = cov(X2)
    X2_eig = eigen(X2_cov)
    
    B_mat = list()
    for (k in 1:d) {
      B_mat[[k]] = as.matrix(X2_eig$vectors[, 1:k])
    }
    
    boot_out[[i]] = list(eigen(X2_cov), B_mat)
  }
  
  ## with the original data
  B_main = list()
  for (k in 1:d) {
    B_main[[k]] = as.matrix(eigen(cov(X_mat))$vectors[, 1:k])
  }
  
  bev = list()
  for (k in 1:(d - 1)) {
    
    res_out = list()
    for (i in 1:n_boots) {
      res_out[[i]] = (1 - abs(det(t(B_main[[k]]) %*% boot_out[[i]][[2]][[k]])))
    }
    
    tot_out = Reduce('+', res_out)
    bev[[k]] = tot_out/n_boots
  }
  
  df_bev = data.frame(bev = c(0, unlist(bev)),
                      k = 0:(d - 1))
  
  return(df_bev)
}


# Making all matrices symmetric in rct data -----------------------------------------
adj_mat_adjustment = function(m, dims, data) {
  
  y = list()
  for (i in dims) {
    y[[i]] = as.matrix((forceSymmetric(as.matrix(data[[1]][[i]][[1]][[m]][[1]]),
                                       uplo = "U") +
                          forceSymmetric(as.matrix(data[[1]][[i]][[1]][[m]][[1]]),
                                         uplo = "L") > 0)*1)
    diag(y[[i]]) = 0
  }
  return(y)
}

# Calculate lambda_1's --------------------------------------------------------------
get_lambda_1 = function(x, dims) {
  
  lambda_list = list()
  for (i in dims) {
    lambda_list[[i]] = max(abs(eigen(x[[i]])$values))
  }
  return(lambda_list)
}

# Regression Function ---------------------------------------------------------------
model_out = function(dep_var, exog_vars, controls, dat) {
  return(lm_robust(reformulate(response = dep_var, termlabels = c(exog_vars, controls) %>% unlist),
                   data = dat))
}

# To Calculate Network Stats --------------------------------------------------------
network_stats_avg = function(adj_mat, dims) {
  
  Y = eval(as.name(adj_mat))
  netout = list()
  
  for (i in dims) {
    
    A = Y[[i]]  
    n = dim(A)[1]
    M = A * (1 - ((A %^% 2) > 0 ))
    
    ## link stats
    ulinks_g = sum(rowSums(M))/2
    triangles_g = sum(diag(A %^% 3))/6
    degree_g = mean(rowSums(A))
    degree_sd = sd(rowSums(A))
    
    ## clustering
    clustering_g = transitivity(graph_from_adjacency_matrix(A, mode = "max"), type="global")
    
    ## density
    density_g = edge_density(graph_from_adjacency_matrix(A, mode = "max"))
    
    netout[[i]] =   c("degree" = degree_g,
                      "degree S.D." = degree_sd,
                      "density" = density_g,
                      "triangles" = triangles_g,
                      "clustering" = clustering_g
                      )
  }
  
  netout = bind_rows(netout)
  return(apply(netout, 2, mean, na.rm=T))

}

network_stats = function(adj_mat, dims) {
  
  Y = eval(as.name(adj_mat))
  netout = list()
  
  for (i in dims) {
    
    A = Y[[i]]  
    n = dim(A)[1]
    M = A * (1 - ((A %^% 2) > 0 ))
    
    ## link stats
    ulinks_g = sum(rowSums(M))/2
    triangles_g = sum(diag(A %^% 3))/6
    degree_g = mean(rowSums(A))
    degree_sd = sd(rowSums(A))
    
    ## clustering
    clustering_g = transitivity(graph_from_adjacency_matrix(A, mode = "max"), type="global")
    
    ## density
    density_g = edge_density(graph_from_adjacency_matrix(A, mode = "max"))
    
    netout[[i]] =   c("degree" = degree_g,
                      "degree S.D." = degree_sd,
                      "density" = density_g,
                      "triangles" = triangles_g,
                      "clustering" = clustering_g
                      )
  }
  
  netout = bind_rows(netout)
  return(netout)
}

# For Puffer N Transformation -------------------------------------------------------
puffer_N_transform = function(data, index_dep, index_exog) {
  
  inp_data = data[, c(index_dep, index_exog)]
  inp_data = na.omit(inp_data)
  
  X = inp_data[, -1]
  y = inp_data[, 1]
  
  names_X = paste0(names(X))
  names_y = paste0(names(y))
  
  X = as.matrix(X)
  y = as.matrix(y)
  
  V = solve(t(X) %*% X, tol = 1e-20)
  V = diag(diag(V))
  N = sqrt(V)
  
  X_N = X %*% N
  X_N.svd = svd(X_N)
  D = diag(X_N.svd$d)
  
  Fn = X_N.svd$u %*% solve(D) %*% t(X_N.svd$u)
  
  X_new = Fn %*% X_N
  y_new = Fn %*% y
  
  dat = cbind(y_new, X_new) %>% as.data.frame() %>% set_names(c(names_y, names_X)) %>% as.matrix()
  return(dat)
  
}

# Post Puffer_N Lasso and OLS -------------------------------------------------------
post_lasso_ols_glmnet = function(exog_selection, lambda_quantile, data, var_names) {
  
  puffer_data = puffer_N_transform(data = data, index_dep = 1, index_exog = grep(exog_selection, var_names))
  
  lasso_glmnet = glmnet(x = puffer_data[, -1], y = puffer_data[, 1])
  
  lambda = quantile(lasso_glmnet$lambda, lambda_quantile)
  
  coefs = coef(lasso_glmnet ,s= lambda)[, 1]
  support_puffer = coefs[coefs != 0] %>% names()
  
  support_puffer = support_puffer[-1] ## remove intercept
  
  post_ols = lm_robust(reformulate(response = "CallsReceived",
                                   termlabels = support_puffer),
                       data = data)
  return(post_ols)
}

lasso_plot_vals = function(s) {
  
  coefs = coef(lasso_glmnet ,s= s)[, 1]
  support_puffer = coefs[coefs != 0] %>% names()
  
  support_puffer = support_puffer[-1] ## remove intercept
  grep("sum|num_hh|dummy", support_puffer, value = T)
  
}

network_stats_diff = function(layer1, layer2, data, netlist_in){

  diff_list = list()
  diff_list[["degree"]] = t.test(stats_out[[which(netlist_in == layer1)]][["degree"]],
    stats_out[[which(netlist_in == layer2)]][["degree"]])$p.value

  diff_list[["degree S.D."]] = t.test(stats_out[[which(netlist_in == layer1)]][["degree S.D."]],
    stats_out[[which(netlist_in == layer2)]][["degree S.D."]])$p.value

  diff_list[["density"]] = t.test(stats_out[[which(netlist_in == layer1)]][["density"]],
    stats_out[[which(netlist_in == layer2)]][["density"]])$p.value

  diff_list[["triangles"]] = t.test(stats_out[[which(netlist_in == layer1)]][["triangles"]],
    stats_out[[which(netlist_in == layer2)]][["triangles"]])$p.value

  diff_list[["clustering"]] = t.test(stats_out[[which(netlist_in == layer1)]][["clustering"]],
    stats_out[[which(netlist_in == layer2)]][["clustering"]])$p.value

  return(diff_list)
}

## F-test table

F_calc = function(r1, r2, data) {

  k1 = length(r1)
  k2 = length(r2)
  lm_1 = lm(reformulate(response = "CallsReceived", termlabels = r1), data = data)
  lm_2 = lm(reformulate(response = "CallsReceived", termlabels = r2), data = data)

  if (k1==k2) {
    a = anova(lm_1)
    F = a$`F value`[1]
    p_val = a$`Pr(>F)`[1]
    return(c("F" = F, "p_val" = p_val))
  }

  a = anova(lm_1, lm_2)
  F = a$F[2]
  p_val = a$`Pr(>F)`[2]
  return(c("F" = F, "p_val" = p_val))

}
