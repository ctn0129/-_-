# ============================================================================
# Bayesian Constrained Ordinal Probit Variable Selection
# Real Data 版本（PKD 資料集）：train_pkd.csv / test_pkd.csv
# ============================================================================

library(MASS)
library(mvtnorm)
library(parallel)
library(truncnorm)
library(coda)
library(mcmcse)

# ----------------------------------------------------------------------------
# 約束條件處理
# ----------------------------------------------------------------------------

compute_gamma <- function(eta, phi, delta, xi) {
  p <- nrow(phi)
  g <- ncol(phi)
  gamma <- numeric(p)
  for (j in 1:p) gamma[j] <- prod(eta^phi[j, ])
  max_iter <- 100
  for (iter in 1:max_iter) {
    gamma_old <- gamma
    for (j in 1:p) {
      if (gamma[j] > 0) {
        hierarchical_part <- 1
        for (i in 1:p) {
          if (i != j && delta[i, j] == 1) hierarchical_part <- hierarchical_part * gamma[i]
        }
        anti_hierarchical_part <- 1
        for (i in 1:p) {
          if (i != j && xi[i, j] == 1) anti_hierarchical_part <- anti_hierarchical_part * (1 - gamma[i])
        }
        gamma[j] <- prod(eta^phi[j, ]) * hierarchical_part * anti_hierarchical_part
      }
    }
    if (all(gamma == gamma_old)) break
  }
  return(gamma)
}

# ----------------------------------------------------------------------------
# 計算邊際概似
# ----------------------------------------------------------------------------

calc_marginal_likelihood_ordinal <- function(X_sub, Y, tau_sq, sigma_sq = 1) {
  n <- length(Y)
  if (is.null(X_sub) || ncol(X_sub) == 0) {
    return(-0.5 * n * log(2 * pi * sigma_sq) - 0.5 * sum(Y^2) / sigma_sq)
  }
  p_sub <- ncol(X_sub)
  if (!is.matrix(X_sub)) X_sub <- as.matrix(X_sub)
  XtX <- crossprod(X_sub)
  Xty <- crossprod(X_sub, Y)
  B <- diag(p_sub) + (tau_sq / sigma_sq) * XtX
  log_det_B <- determinant(B, logarithm = TRUE)$modulus
  A <- XtX / sigma_sq + diag(1/tau_sq, p_sub)
  tryCatch({
    L <- chol(A)
    Sigma_post <- chol2inv(L)
  }, error = function(e) Sigma_post <<- solve(A))
  mu_post <- Sigma_post %*% (Xty / sigma_sq)
  quad_term <- as.numeric(0.5 * crossprod(mu_post, A %*% mu_post))
  log_mlik <- -0.5 * n * log(2 * pi * sigma_sq) - 0.5 * sum(Y^2) / sigma_sq -
    0.5 * log_det_B + quad_term
  return(log_mlik)
}

# ----------------------------------------------------------------------------
# 單條 MCMC 核心抽樣函數
# ----------------------------------------------------------------------------

run_single_chain <- function(C, X_scaled, phi, delta, xi,
                             K = 3,
                             b_sq = 100,
                             a_tau = 1,
                             w = 0.5,
                             n_iter,
                             init_beta = NULL,
                             init_tau = NULL,
                             init_eta = NULL,
                             init_Y = NULL,
                             chain_seed = NULL) {
  
  if (!is.null(chain_seed)) set.seed(chain_seed)
  
  n <- length(C)
  p <- ncol(X_scaled)
  g <- ncol(phi)
  
  beta <- if (!is.null(init_beta)) init_beta else rep(0, p)
  eta  <- if (!is.null(init_eta))  init_eta  else rep(1, g)
  tau  <- if (!is.null(init_tau))  init_tau  else seq(-1, 1, length.out = K - 1)
  Y    <- if (!is.null(init_Y))    init_Y    else rnorm(n, X_scaled %*% beta, 1)
  
  gamma <- compute_gamma(eta, phi, delta, xi)
  
  eta_samples   <- matrix(0, nrow = n_iter, ncol = g)
  gamma_samples <- matrix(0, nrow = n_iter, ncol = p)
  beta_samples  <- matrix(0, nrow = n_iter, ncol = p)
  tau_samples   <- matrix(0, nrow = n_iter, ncol = K - 1)
  
  for (iter in 1:n_iter) {
    
    gamma <- compute_gamma(eta, phi, delta, xi)
    
    active <- which(gamma == 1)
    beta   <- rep(0, p)
    if (length(active) > 0) {
      X_act        <- X_scaled[, active, drop = FALSE]
      A            <- crossprod(X_act) + diag(1/b_sq, length(active))
      L            <- tryCatch(chol(A), error = function(e) NULL)
      Sigma        <- if (!is.null(L)) chol2inv(L) else solve(A)
      mu           <- Sigma %*% crossprod(X_act, Y)
      beta[active] <- mvrnorm(1, as.vector(mu), Sigma)
    }
    
    order_k <- sample(1:g)
    for (k in order_k) {
      eta_0 <- eta; eta_0[k] <- 0
      eta_1 <- eta; eta_1[k] <- 1
      
      gamma_0 <- compute_gamma(eta_0, phi, delta, xi)
      gamma_1 <- compute_gamma(eta_1, phi, delta, xi)
      
      affected_vars <- which(gamma_0 != gamma_1)
      
      if (length(affected_vars) == 0) {
        eta[k] <- rbinom(1, 1, w)
      } else {
        vars_H0     <- which(gamma_0 == 1)
        log_mlik_H0 <- if (length(vars_H0) > 0)
          calc_marginal_likelihood_ordinal(X_scaled[, vars_H0, drop=FALSE], Y, b_sq)
        else -0.5 * sum(Y^2) - 0.5 * n * log(2 * pi)
        
        vars_H1     <- which(gamma_1 == 1)
        log_mlik_H1 <- if (length(vars_H1) > 0)
          calc_marginal_likelihood_ordinal(X_scaled[, vars_H1, drop=FALSE], Y, b_sq)
        else -0.5 * sum(Y^2) - 0.5 * n * log(2 * pi)
        
        log_BF  <- log_mlik_H1 - log_mlik_H0
        log_p1  <- log(w) + log_BF
        log_p0  <- log(1 - w)
        log_max <- max(log_p1, log_p0)
        prob_1  <- exp(log_p1 - log_max) / (exp(log_p1 - log_max) + exp(log_p0 - log_max))
        prob_1  <- max(1e-10, min(1 - 1e-10, prob_1))
        eta[k]  <- rbinom(1, 1, prob_1)
      }
    }
    
    tau_star <- tau
    for (k in 1:(K-1)) {
      lower_bound <- if (k == 1) -Inf else tau_star[k-1]
      upper_bound <- if (k == K-1) Inf else tau[k+1]
      tau_star[k] <- rtruncnorm(1, a = lower_bound, b = upper_bound,
                                mean = tau[k], sd = a_tau)
    }
    
    log_ratio <- 0
    for (i in 1:n) {
      mu_i <- X_scaled[i, ] %*% beta
      if (C[i] == 1) {
        prob_star <- pnorm(tau_star[1] - mu_i)
        prob_old  <- pnorm(tau[1]      - mu_i)
      } else if (C[i] == K) {
        prob_star <- 1 - pnorm(tau_star[K-1] - mu_i)
        prob_old  <- 1 - pnorm(tau[K-1]      - mu_i)
      } else {
        prob_star <- pnorm(tau_star[C[i]] - mu_i) - pnorm(tau_star[C[i]-1] - mu_i)
        prob_old  <- pnorm(tau[C[i]]      - mu_i) - pnorm(tau[C[i]-1]      - mu_i)
      }
      log_ratio <- log_ratio + log(max(prob_star, 1e-10)) - log(max(prob_old, 1e-10))
    }
    log_ratio <- log_ratio - 0.5 / (a_tau^2) * (sum(tau_star^2) - sum(tau^2))
    if (log(runif(1)) < log_ratio) tau <- tau_star
    
    mu <- X_scaled %*% beta
    for (i in 1:n) {
      if (C[i] == 1) {
        Y[i] <- rtruncnorm(1, a = -Inf,        b = tau[1],     mean = mu[i], sd = 1)
      } else if (C[i] == K) {
        Y[i] <- rtruncnorm(1, a = tau[K-1],    b = Inf,        mean = mu[i], sd = 1)
      } else {
        Y[i] <- rtruncnorm(1, a = tau[C[i]-1], b = tau[C[i]], mean = mu[i], sd = 1)
      }
    }
    
    eta_samples[iter, ]   <- eta
    gamma_samples[iter, ] <- gamma
    beta_samples[iter, ]  <- beta
    tau_samples[iter, ]   <- tau
  }
  
  return(list(
    eta_samples   = eta_samples,
    gamma_samples = gamma_samples,
    beta_samples  = beta_samples,
    tau_samples   = tau_samples,
    last_beta = beta,
    last_tau  = tau,
    last_eta  = eta,
    last_Y    = Y
  ))
}

# ----------------------------------------------------------------------------
# 模型選擇與預測
# ----------------------------------------------------------------------------

select_median_model <- function(gamma_samples) {
  posterior_prob <- colMeans(gamma_samples)
  median_model   <- as.numeric(posterior_prob > 0.5)
  return(list(median_model = median_model, posterior_prob = posterior_prob))
}

predict_ordinal_probit <- function(X_new, beta, tau) {
  n_new <- nrow(X_new)
  K     <- length(tau) + 1
  pred  <- numeric(n_new)
  mu    <- X_new %*% beta
  for (i in 1:n_new) {
    probs    <- numeric(K)
    probs[1] <- pnorm(tau[1] - mu[i])
    if (K > 2) for (k in 2:(K-1)) probs[k] <- pnorm(tau[k] - mu[i]) - pnorm(tau[k-1] - mu[i])
    probs[K] <- 1 - pnorm(tau[K-1] - mu[i])
    pred[i]  <- which.max(probs)
  }
  return(pred)
}

# ----------------------------------------------------------------------------
# 讀取 train_pkd.csv 及 test_pkd.csv
# ----------------------------------------------------------------------------

load_pkd_data <- function(train_path,
                          test_path,
                          target_col     = "Y",
                          cor_threshold  = 0.93,
                          non_image_vars = c("Sex", "Age")) {
  
  train_df <- read.csv(train_path, check.names = FALSE)
  test_df  <- read.csv(test_path,  check.names = FALSE)
  
  C_train <- as.integer(train_df[[target_col]])
  C_test  <- as.integer(test_df[[target_col]])
  K <- length(unique(c(C_train, C_test)))
  
  feat_cols <- setdiff(names(train_df), target_col)
  
  int_prefix_pattern <- paste0("^(", paste(non_image_vars, collapse = "|"), ")_")
  is_int_col  <- grepl(int_prefix_pattern, feat_cols)
  main_cols   <- feat_cols[!is_int_col]
  int_cols_raw <- feat_cols[is_int_col]
  
  image_vars <- setdiff(main_cols, non_image_vars)
  
  missing_niv <- setdiff(non_image_vars, main_cols)
  if (length(missing_niv) > 0)
    stop(sprintf("找不到非影像變數：%s", paste(missing_niv, collapse = ", ")))
  
  rename_int_col <- function(col_name) {
    for (niv in non_image_vars) {
      prefix <- paste0(niv, "_")
      if (startsWith(col_name, prefix)) {
        suffix <- substring(col_name, nchar(prefix) + 1)
        return(paste0(niv, ":", suffix))
      }
    }
    return(col_name)
  }
  int_names <- vapply(int_cols_raw, rename_int_col, character(1), USE.NAMES = FALSE)
  
  X_train_main <- as.matrix(train_df[, main_cols,    drop = FALSE])
  X_test_main  <- as.matrix(test_df[,  main_cols,    drop = FALSE])
  X_train_int  <- as.matrix(train_df[, int_cols_raw, drop = FALSE])
  X_test_int   <- as.matrix(test_df[,  int_cols_raw, drop = FALSE])
  colnames(X_train_int) <- int_names
  colnames(X_test_int)  <- int_names
  
  X_train <- cbind(X_train_main, X_train_int)
  X_test  <- cbind(X_test_main,  X_test_int)
  var_names <- colnames(X_train)
  p_main    <- length(main_cols)
  p_total   <- ncol(X_train)
  
  cor_mat  <- cor(X_train_main, use = "pairwise.complete.obs")
  assigned <- rep(FALSE, p_main)
  groups   <- list()
  
  for (i in seq_len(p_main)) {
    if (assigned[i]) next
    current_group <- i
    for (j in seq_len(p_main)) {
      if (j == i || assigned[j]) next
      if (all(!is.na(cor_mat[j, current_group])) &&
          all(abs(cor_mat[j, current_group]) > cor_threshold)) {
        current_group <- c(current_group, j)
      }
    }
    groups[[length(groups) + 1]] <- current_group
    assigned[current_group] <- TRUE
  }
  
  main_group_assign <- character(p_main)
  for (grp in groups) {
    gname <- main_cols[grp[1]]
    for (idx in grp) main_group_assign[idx] <- gname
  }
  names(main_group_assign) <- main_cols
  
  group_sizes         <- table(main_group_assign)
  multi_member_groups <- names(group_sizes[group_sizes > 1])
  
  return(list(
    X_train             = X_train,
    X_test              = X_test,
    C_train             = C_train,
    C_test              = C_test,
    var_names           = var_names,
    main_names          = main_cols,
    int_names           = int_names,
    main_group_assign   = main_group_assign,
    multi_member_groups = multi_member_groups,
    non_image_vars      = non_image_vars,
    image_vars          = image_vars,
    p_main              = p_main,
    p_total             = p_total,
    K                   = K,
    n_train             = nrow(X_train),
    n_test              = nrow(X_test),
    cor_threshold       = cor_threshold
  ))
}

# ----------------------------------------------------------------------------
# 建立 phi、delta、xi 約束矩陣
# ----------------------------------------------------------------------------

build_constraint_matrices_pkd <- function(var_names, main_names, int_names,
                                          main_group_assign,
                                          multi_member_groups,
                                          group_interactions = FALSE) {
  
  p      <- length(var_names)
  p_main <- length(main_names)
  
  if (!group_interactions) {
    int_group_assign <- int_names
  } else {
    int_group_assign <- character(length(int_names))
    names(int_group_assign) <- int_names
    
    for (int_name in int_names) {
      parts <- strsplit(int_name, ":", fixed = TRUE)[[1]]
      niv   <- parts[1]
      imv   <- parts[2]
      img_group_leader <- main_group_assign[imv]
      
      if (img_group_leader %in% multi_member_groups) {
        int_group_assign[int_name] <- paste0(niv, "__grp__", img_group_leader)
      } else {
        int_group_assign[int_name] <- int_name
      }
    }
  }
  
  all_group_assign <- c(main_group_assign, int_group_assign)
  
  group_names <- unique(all_group_assign)
  g           <- length(group_names)
  
  phi <- matrix(0L, nrow = p, ncol = g,
                dimnames = list(var_names, group_names))
  for (j in seq_len(p)) {
    k         <- which(group_names == all_group_assign[j])
    phi[j, k] <- 1L
  }
  
  delta <- matrix(0L, nrow = p, ncol = p,
                  dimnames = list(var_names, var_names))
  
  for (int_name in int_names) {
    parts   <- strsplit(int_name, ":", fixed = TRUE)[[1]]
    parent1 <- parts[1]
    parent2 <- parts[2]
    idx_int <- which(var_names == int_name)
    idx_p1  <- which(var_names == parent1)
    idx_p2  <- which(var_names == parent2)
    if (length(idx_int) == 1 && length(idx_p1) == 1)
      delta[idx_p1, idx_int] <- 1L
    if (length(idx_int) == 1 && length(idx_p2) == 1)
      delta[idx_p2, idx_int] <- 1L
  }
  
  xi <- matrix(0L, nrow = p, ncol = p,
               dimnames = list(var_names, var_names))
  
  return(list(
    phi              = phi,
    delta            = delta,
    xi               = xi,
    group_names      = group_names,
    all_group_assign = all_group_assign,
    g                = g
  ))
}

# ----------------------------------------------------------------------------
# Cross-Validation 選取最佳 b_sq
# ----------------------------------------------------------------------------

cv_select_tau1sq <- function(X_train,
                             C_train,
                             phi,
                             delta,
                             xi,
                             K             = 3,
                             b_sq_grid  = c(1, 2, 3, 5, 10),
                             n_folds       = 5,
                             cv_burnin     = 2000,
                             cv_post       = 1000,
                             w             = 0.5,
                             a_tau         = 1,
                             seed          = 42,
                             verbose       = TRUE) {
  
  set.seed(seed)
  n <- nrow(X_train)
  
  fold_ids <- sample(rep(1:n_folds, length.out = n))
  
  cv_acc <- numeric(length(b_sq_grid))
  
  for (ti in seq_along(b_sq_grid)) {
    b_sq_val <- b_sq_grid[ti]
    fold_acc    <- numeric(n_folds)
    
    for (fold in 1:n_folds) {
      val_idx   <- which(fold_ids == fold)
      train_idx <- which(fold_ids != fold)
      
      X_sub_train <- X_train[train_idx, , drop = FALSE]
      X_sub_val   <- X_train[val_idx,   , drop = FALSE]
      C_sub_train <- C_train[train_idx]
      C_sub_val   <- C_train[val_idx]
      
      res_burnin <- run_single_chain(
        C = C_sub_train, X_scaled = X_sub_train,
        phi = phi, delta = delta, xi = xi,
        K = K, b_sq = b_sq_val, a_tau = a_tau, w = w,
        n_iter = cv_burnin,
        chain_seed = seed * 42 + fold
      )
      
      res_post <- run_single_chain(
        C = C_sub_train, X_scaled = X_sub_train,
        phi = phi, delta = delta, xi = xi,
        K = K, b_sq = b_sq_val, a_tau = a_tau, w = w,
        n_iter    = cv_post,
        init_beta = res_burnin$last_beta,
        init_tau  = res_burnin$last_tau,
        init_eta  = res_burnin$last_eta,
        init_Y    = res_burnin$last_Y,
        chain_seed = seed * 42 + fold + 50000
      )
      
      post_prob <- colMeans(res_post$gamma_samples)
      selected  <- as.numeric(post_prob > 0.5)
      beta_est  <- colMeans(res_post$beta_samples)
      tau_est   <- colMeans(res_post$tau_samples)
      beta_est[selected == 0] <- 0
      
      pred_val       <- predict_ordinal_probit(X_sub_val, beta_est, tau_est)
      fold_acc[fold] <- mean(pred_val == C_sub_val)
    }
    
    cv_acc[ti] <- mean(fold_acc)
    
    if (verbose) {
      cat(sprintf("  tau1_sq = %s,平均 ACC = %.4f\n",
                  as.character(b_sq_val), cv_acc[ti]))
    }
  }
  
  best_idx    <- which.max(cv_acc)
  best_tau1sq <- b_sq_grid[best_idx]
  
  return(list(
    best_tau1sq  = best_tau1sq,
    cv_acc       = setNames(cv_acc, as.character(b_sq_grid)),
    b_sq_grid = b_sq_grid
  ))
}

# ----------------------------------------------------------------------------
# 找收斂期
# ----------------------------------------------------------------------------

find_burnin_realdata <- function(
    X_train,
    C_train,
    phi,
    delta,
    xi,
    K              = 3,
    b_sq        = 100,
    w              = 0.5,
    n_chains       = 5,
    max_burnin     = 20000,
    check_every    = 1000,
    grd_threshold  = 1.1,
    mcse_threshold = 0.06,
    verbose        = TRUE) {
  
  n       <- nrow(X_train)
  p_total <- ncol(X_train)
  g       <- ncol(phi)
  
  X_scaled <- X_train
  
  chain_states <- vector("list", n_chains)
  for (ch in 1:n_chains) {
    tau_init <- switch(ch,
                       seq(-2,   2,   length.out = K - 1),
                       seq(-3,   3,   length.out = K - 1),
                       seq(-1,   1,   length.out = K - 1),
                       seq(-2.5, 1.5, length.out = K - 1),
                       seq(-1.5, 2.5, length.out = K - 1))
    beta_init <- rnorm(p_total, 0, 0.1) * ch
    chain_states[[ch]] <- list(
      beta = beta_init,
      tau  = tau_init,
      eta  = rep(1, g),
      Y    = rnorm(n, X_scaled %*% beta_init, 1)
    )
  }
  
  chain_tau_history   <- lapply(1:n_chains, function(x) matrix(nrow = 0, ncol = K - 1))
  chain_gamma_history <- lapply(1:n_chains, function(x) matrix(nrow = 0, ncol = p_total))
  
  converged   <- FALSE
  burnin_used <- max_burnin
  n_blocks    <- max_burnin %/% check_every
  final_rhat  <- NA
  final_mcse  <- NA
  
  for (block in 1:n_blocks) {
    for (ch in 1:n_chains) {
      st <- chain_states[[ch]]
      res_ch <- run_single_chain(
        C = C_train, X_scaled = X_scaled,
        phi = phi, delta = delta, xi = xi,
        K = K, b_sq = b_sq, a_tau = 1, w = w,
        n_iter    = check_every,
        init_beta = st$beta, init_tau = st$tau,
        init_eta  = st$eta,  init_Y   = st$Y,
        chain_seed = ch * 1000 + block
      )
      chain_states[[ch]] <- list(
        beta = res_ch$last_beta, tau = res_ch$last_tau,
        eta  = res_ch$last_eta,  Y   = res_ch$last_Y
      )
      chain_tau_history[[ch]]   <- rbind(chain_tau_history[[ch]],   res_ch$tau_samples)
      chain_gamma_history[[ch]] <- rbind(chain_gamma_history[[ch]], res_ch$gamma_samples)
    }
    
    steps_so_far <- block * check_every
    
    max_tau_rhat <- NA
    if (block > 2) {
      tryCatch({
        mcmc_list_tau <- lapply(1:n_chains, function(ch) coda::mcmc(chain_tau_history[[ch]]))
        mc_tau        <- coda::mcmc.list(mcmc_list_tau)
        max_tau_rhat  <- coda::gelman.diag(mc_tau, multivariate = TRUE)$mpsrf
      }, error = function(e) { max_tau_rhat <<- NA })
    }
    
    mcse_per_chain <- rep(NA, n_chains)
    for (ch in 1:n_chains) {
      gamma_mat <- chain_gamma_history[[ch]]
      n_ch      <- nrow(gamma_mat)
      if (n_ch < 10) next
      tryCatch({
        mcse_result        <- mcmcse::mcse.mat(gamma_mat, size = "cuberoot")
        mcse_per_chain[ch] <- max(mcse_result[, 2], na.rm = TRUE)
      }, error = function(e) { mcse_per_chain[ch] <<- NA })
    }
    max_mcse_stat <- max(mcse_per_chain, na.rm = TRUE)
    
    final_rhat <- max_tau_rhat
    final_mcse <- max_mcse_stat
    
    tau_ok  <- !is.na(max_tau_rhat) && max_tau_rhat < grd_threshold
    mcse_ok <- !is.na(max_mcse_stat) && max_mcse_stat < mcse_threshold
    
    if (tau_ok && mcse_ok) {
      burnin_used <- steps_so_far
      converged   <- TRUE
      break
    }
  }
  
  return(list(
    burnin     = burnin_used,
    converged  = converged,
    final_rhat = final_rhat,
    final_mcse = final_mcse
  ))
}

# ----------------------------------------------------------------------------
# main function
# ----------------------------------------------------------------------------

run_pkd_analysis <- function(
    train_path,
    test_path,
    target_col     = "Y",
    cor_threshold  = 0.93,
    non_image_vars = c("Sex", "Age"),
    group_interactions = FALSE,
    b_sq        = NULL,
    b_sq_grid   = c(1, 2, 3, 5, 10),
    cv_folds       = 5,
    cv_burnin      = 2000,
    cv_post        = 1000,
    w              = 0.5,
    a_tau          = 1,
    n_chains       = 5,
    max_burnin     = 20000,
    check_every    = 1000,
    grd_threshold  = 1.1,
    mcse_threshold = 0.06,
    n_post_burnin  = 2000,
    verbose        = TRUE,
    preset_burnin  = NULL,
    data_seed      = 42) {
  
  dat <- load_pkd_data(
    train_path     = train_path,
    test_path      = test_path,
    target_col     = target_col,
    cor_threshold  = cor_threshold,
    non_image_vars = non_image_vars
  )
  X_train             <- dat$X_train
  X_test              <- dat$X_test
  C_train             <- dat$C_train
  C_test              <- dat$C_test
  var_names           <- dat$var_names
  main_names          <- dat$main_names
  int_names           <- dat$int_names
  main_group_assign   <- dat$main_group_assign
  multi_member_groups <- dat$multi_member_groups
  K                   <- dat$K
  
  cons <- build_constraint_matrices_pkd(var_names, main_names, int_names,
                                        main_group_assign, multi_member_groups, group_interactions = group_interactions)
  phi         <- cons$phi
  delta       <- cons$delta
  xi          <- cons$xi
  group_names <- cons$group_names
  g           <- cons$g
  
  if (is.null(b_sq)) {
    cv_result <- cv_select_tau1sq(
      X_train      = X_train,
      C_train      = C_train,
      phi          = phi,
      delta        = delta,
      xi           = xi,
      K            = K,
      b_sq_grid = b_sq_grid,
      n_folds      = cv_folds,
      cv_burnin    = cv_burnin,
      cv_post      = cv_post,
      w            = w,
      a_tau        = a_tau,
      seed         = data_seed,
      verbose      = verbose
    )
    b_sq_used <- cv_result$best_tau1sq
    cat(sprintf("\nCV 選定 tau1_sq = %g\n", b_sq_used))
  } else {
    b_sq_used <- b_sq
    cv_result    <- NULL
  }
  
  if (!is.null(preset_burnin)) {
    fixed_burnin     <- preset_burnin
    burnin_converged <- NA
  } else {
    burnin_info <- find_burnin_realdata(
      X_train        = X_train,
      C_train        = C_train,
      phi            = phi,
      delta          = delta,
      xi             = xi,
      K              = K,
      b_sq        = b_sq_used,
      w              = w,
      n_chains       = n_chains,
      max_burnin     = max_burnin,
      check_every    = check_every,
      grd_threshold  = grd_threshold,
      mcse_threshold = mcse_threshold,
      verbose        = verbose
    )
    fixed_burnin     <- burnin_info$burnin
    burnin_converged <- burnin_info$converged
  }
  cat(sprintf("確定 burn-in = %d 步\n", fixed_burnin))
  
  burnin_result <- run_single_chain(
    C = C_train, X_scaled = X_train,
    phi = phi, delta = delta, xi = xi,
    K = K, b_sq = b_sq_used, a_tau = a_tau, w = w,
    n_iter     = fixed_burnin,
    chain_seed = data_seed
  )
  
  post_result <- run_single_chain(
    C = C_train, X_scaled = X_train,
    phi = phi, delta = delta, xi = xi,
    K = K, b_sq = b_sq_used, a_tau = a_tau, w = w,
    n_iter     = n_post_burnin,
    init_beta  = burnin_result$last_beta,
    init_tau   = burnin_result$last_tau,
    init_eta   = burnin_result$last_eta,
    init_Y     = burnin_result$last_Y,
    chain_seed = data_seed + 99999
  )
  
  model_sel      <- select_median_model(post_result$gamma_samples)
  selected_model <- model_sel$median_model
  posterior_prob <- model_sel$posterior_prob
  M              <- sum(selected_model)
  
  selected_vars <- var_names[selected_model == 1]
  sel_main <- selected_vars[selected_vars %in% main_names]
  sel_int  <- selected_vars[!(selected_vars %in% main_names)]
  
  cat(sprintf("選中變數數 M = %d\n", M))
  if (length(sel_main) > 0) cat("  主效應:",   paste(sel_main, collapse = ", "), "\n")
  if (length(sel_int)  > 0) cat("  交互作用:", paste(sel_int,  collapse = ", "), "\n")
  
  prob_df <- data.frame(
    Variable  = var_names,
    Type      = ifelse(var_names %in% main_names, "主效應", "交互作用"),
    Group     = cons$all_group_assign,
    Post_Prob = round(posterior_prob, 4),
    Selected  = selected_model == 1,
    stringsAsFactors = FALSE
  )
  prob_df <- prob_df[order(-prob_df$Post_Prob), ]
  
  prob_df_show <- prob_df[prob_df$Post_Prob > 0.05, ]
  cat(sprintf("各變數後驗選入機率(僅顯示 > 0.05,共 %d 個):\n", nrow(prob_df_show)))
  print(prob_df_show, row.names = FALSE)
  cat("\n")
  
  beta_est <- colMeans(post_result$beta_samples)
  tau_est  <- colMeans(post_result$tau_samples)
  beta_est[selected_model == 0] <- 0
  
  pred_train <- predict_ordinal_probit(X_train, beta_est, tau_est)
  pred_test  <- predict_ordinal_probit(X_test,  beta_est, tau_est)
  ACC_train  <- mean(pred_train == C_train)
  ACC_test   <- mean(pred_test  == C_test)
  conf_train <- table(Predicted = pred_train, Actual = C_train)
  conf_test  <- table(Predicted = pred_test,  Actual = C_test)
  
  cat(sprintf("  訓練集 ACC = %.4f  (n = %d)\n", ACC_train, dat$n_train))
  cat(sprintf("  測試集 ACC = %.4f  (n = %d)\n\n", ACC_test,  dat$n_test))
  cat("Confusion Matrix【訓練集】:\n"); print(conf_train); cat("\n")
  cat("Confusion Matrix【測試集】:\n");  print(conf_test);  cat("\n")
  
  return(list(
    selected_vars       = selected_vars,
    selected_main       = sel_main,
    selected_int        = sel_int,
    selected_model      = selected_model,
    posterior_prob      = posterior_prob,
    prob_df             = prob_df,
    beta_est            = beta_est,
    tau_est             = tau_est,
    ACC_train           = ACC_train,
    ACC_test            = ACC_test,
    conf_train          = conf_train,
    conf_test           = conf_test,
    M                   = M,
    b_sq_used        = b_sq_used,
    cv_result           = cv_result,
    burnin_used         = fixed_burnin,
    burnin_converged    = burnin_converged,
    post_gamma_samples  = post_result$gamma_samples,
    post_beta_samples   = post_result$beta_samples,
    post_tau_samples    = post_result$tau_samples,
    phi                 = phi,
    delta               = delta,
    xi                  = xi,
    group_names         = group_names,
    var_names           = var_names,
    main_names          = main_names,
    int_names           = int_names,
    main_group_assign   = main_group_assign,
    multi_member_groups = multi_member_groups
  ))
}

# ----------------------------------------------------------------------------
# 執行
# ----------------------------------------------------------------------------

results <- run_pkd_analysis(
  train_path     = "train_pkd.csv",
  test_path      = "test_pkd.csv",
  target_col     = "Y",
  cor_threshold  = 0.92,
  non_image_vars = c("Sex", "Age"),
  group_interactions = FALSE,                #FALSE代表將交互作用項視為singleton，TRUE代表將交互作用項分成群組
  b_sq        = NULL,                        #這裡可以填入指定的b_sq，即可跳過交叉驗證
  b_sq_grid   = c(0.01, 0.04, 0.09, 0.25, 1, 4, 9, 25, 100),
  cv_folds       = 5,
  cv_burnin      = 5000,
  cv_post        = 2000,
  w              = 0.5,
  a_tau          = 1,
  n_chains       = 5,
  max_burnin     = 20000,
  check_every    = 1000,
  grd_threshold  = 1.1,
  mcse_threshold = 0.06,
  n_post_burnin  = 9000,
  verbose        = TRUE,
  preset_burnin  = NULL
)