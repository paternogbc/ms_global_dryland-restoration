# Function to check common column names
cc <- function(x, y) {
  names(x)[names(x) %in% names(y)]
  
}

# Check overdispersion and zero-inflantion
zinf <- function(x, xlab = "variable") {
  dd <- data.frame(x = x)
  dd %>%
    summarise(
      n_obs    = n(),
      n_na     = sum(is.na(x), na.rm = T),
      n_zero   = sum(x == 0, na.rm = T),
      p_zero   = round((n_zero / n_obs),digits = 3),
      p_na     = round((n_na / n_obs),digits = 3),
      p_va     = round(((n_obs - n_na - n_zero) / n_obs), digits = 3),
      mean     = mean(x, na.rm = T),
      median   = median(x, na.rm = T),
      sd       = sd(x, na.rm = T),
      var      = var(x, na.rm = T),
      disper   = var / mean,
      zip      = 1 + log(p_zero) / mean
    ) -> out
  
  print(out)
  
  g1 <-
    ggplot(dd %>% drop_na(x), aes(x = x)) +
    geom_histogram(alpha = .7, fill = "red", color = "white") +
    labs(x = xlab,
         y = "Frequency",
         subtitle = "raw format") +
    theme_cowplot(font_size = 16);
  
  g2 <-
    ggplot(dd %>% drop_na(x), aes(x = sqrt(x))) +
    geom_histogram(alpha = .7, fill = "red", color = "white") +
    labs(x = xlab,
         y = "Frequency",
         subtitle = "sqrt(x)");
  
  g3 <-
    ggplot(dd %>% drop_na(x), aes(x = log(x + 1))) +
    geom_histogram(alpha = .7, fill = "red", color = "white") +
    labs(x = xlab,
         y = "Frequency",
         subtitle = "log(x + 1)");
  
  
  g4 <-
    ggplot(dd %>% drop_na(x), aes(x = log(x))) +
    geom_histogram(alpha = .7, fill = "red", color = "white") +
    labs(x = xlab,
         y = "Frequency",
         subtitle = "log(x) | droping zeros");
  
  g5 <-
    ggplot(dd %>% drop_na(x), aes(x = x ^ (2/3))) +
    geom_histogram(alpha = .7, fill = "red", color = "white") +
    labs(x = xlab,
         y = "Frequency",
         subtitle = "x ^ 2/3");
  g6 <-
    ggplot(dd %>% drop_na(x), aes(x = x ^ 2)) +
    geom_histogram(alpha = .7, fill = "red", color = "white") +
    labs(x = xlab,
         y = "Frequency",
         subtitle = "x^2");
  
  # title
  tit_text <- paste(
    xlab, "has", 
    out$n_obs, "observations |",
    paste(out$n_na), "NA`s |",
    paste(out$p_zero * 100, "%"), "of zeros | ",
    paste(out$p_va * 100, "%"), "of non-zero-na")
  tit_text
  plots <- plot_grid(g1, g2, g3, g4, g5, g6)
  tit   <- ggdraw() + 
    draw_label(tit_text, fontface = 'bold')
  
  plots <- plot_grid(tit, plots, ncol = 1, rel_heights = c(0.1, 1))
  
  return(list(stats = out,
              plots = plots))
}

# Beta-rescale------------------------------------------------------------------
transform_beta <- function(x) {
  (x * (length(x) - 1) + 0.5) / (length(x))
}

# Function to match seedrates by genus and lifeform-----------------------------
genuscheck <- function(miss_seed, full_data) {
  for(ms in miss_seed) {
    
    if (substr(full_data$speciesid[ms],
               start = 1, stop = 2) == "G_") {
      
      if(full_data$projectid[ms] %in% proj$studyseedrates$projectid) {
        
        genus <- full_data$genus[ms]
        rates <- proj$studyseedrates %>%
          filter(projectid == full_data$projectid[ms]) %>%
          left_join(proj$species, by = "speciesid")
        
        if(genus %in% rates$genus) {
          
          full_data$seededrate[ms] <- 
            sum(rates$seededrate[rates$genus == genus])
          
        }
        
      } else {
        
        if(full_data$treatmentid[ms] %in% proj$trtseedrates$treatmentid) {
          
          genus <- full_data$genus[ms]
          rates <- proj$trtseedrates %>%
            filter(treatmentid == full_data$treatmentid[ms]) %>%
            left_join(proj$species, by = "speciesid")
          
          if(genus %in% rates$genus) {
            
            full_data$seededrate[ms] <- 
              sum(rates$seededrate[rates$genus == genus])
            
          }      
        }
      }
      
    } else next
  }
  
  return(full_data)
}


# Function to match seedrates by lifeform
lifecheck <- function(miss_seed, full_data) {
  for(ms in miss_seed) {
    
    if (substr(full_data$speciesid[ms], 
               start = 1, stop = 2) == "L_") {
      
      if(full_data$projectid[ms] %in% proj$studyseedrates$projectid) {
        
        lifeform <- full_data$lifeform[ms]
        rates <- proj$studyseedrates %>%
          filter(projectid == full_data$projectid[ms]) %>%
          left_join(proj$species, by = "speciesid")
        
        if(lifeform %in% rates$lifeform) {
          
          full_data$seededrate[ms] <- 
            sum(rates$seededrate[rates$lifeform == lifeform])
          
        }
        
      } else {
        
        if(full_data$treatmentid[ms] %in% proj$trtseedrates$treatmentid) {
          
          lifeform <- full_data$lifeform[ms]
          rates <- proj$trtseedrates %>%
            filter(treatmentid == full_data$treatmentid[ms]) %>%
            left_join(proj$species, by = "speciesid")
          
          if(lifeform %in% rates$lifeform) {
            
            full_data$seededrate[ms] <- 
              sum(rates$seededrate[rates$lifeform == lifeform])
            
          }      
        }
      }
      
    } else next
  }
  
  return(full_data)
}

# Function to organize glmmTMB summary table------------------------------------
summary_glmmTMB <- function(x) {
  if (!inherits(x, "glmmTMB")) {
    stop("x must be an object of class `glmmTMB`")
  }
  msum <- summary(x)
  econd   <- msum$coefficients$cond
  ezi   <- msum$coefficients$zi
  
  cond <- cbind(data.frame(model = "Conditional",
                           term = rownames(econd)),
                econd)
  if (!is.null(ezi)) {
    zinf <- cbind(data.frame(model = "Zero-inflation",
                             term = rownames(ezi)),
                  ezi)
    estim <- rbind(cond, zinf)
    
    
  } else {
    estim <- cond
  }
  
  estim <-
    estim %>%
    rename(
      estimate = Estimate,
      std.error = `Std. Error`,
      statistic = `z value`,
      p.value = `Pr(>|z|)`
    ) %>%
    as.tibble()
  return(estim)
}

diagnose_vcov <- function(model, tol=1e-5, digits=2, analyze_hessian=FALSE) {
  vv <- vcov(model, full=TRUE)
  nn <- rownames(vv)
  if (!all(is.finite(vv))) {
    if (missing(analyze_hessian)) warning("analyzing Hessian, not vcov")
    if (!analyze_hessian) stop("can't analyze vcov")
    analyze_hessian <- TRUE
  }
  if (analyze_hessian) {
    par.fixed <- model$obj$env$last.par.best
    r <- model$obj$env$random
    if (!is.null(r)) par.fixed <- par.fixed[-r]
    vv <- optimHess(par.fixed, fn=model$obj$fn, gr=model$obj$gr)
    ## note vv is now HESSIAN, not vcov
  }
  ee <- eigen(vv)
  if (all(ee$values>tol)) {message("var-cov matrix OK"); return(invisible(NULL))}
  ## find negative or small-positive eigenvalues (flat/wrong curvature)
  bad_evals <- which(ee$values<tol)
  ## order worst to best
  bad_evals <- bad_evals[order(-ee$values[bad_evals])]
  ret <- lapply(bad_evals,
                function(i) {
                  ## extract loadings
                  v <- setNames(ee$vectors[,i], nn)
                  ## order in decreasing magnitude & round
                  list(val=ee$values[i],vec=round(v[order(-abs(v))],digits))
                })
  return(ret)
}

## print table; add space (From Ben Bolker)
pxt <- function(x,title) {
  cat(sprintf("{\n\n\\textbf{%s}\n\\ \\\\\\vspace{2pt}\\ \\\\\n",title))
  print(xtable(x), floating=FALSE); cat("\n\n")
  cat("\\ \\\\\\vspace{5pt}\\ \\\\\n")
}

# Un-scale predictions 
unzt <- function(formula, data, xu) {
  mod <- lm(formula, data = data)
  coefs <- coef(mod)
  xl <- log(xu)
  xz <- (xl - coefs[1]) / coefs[2]
  data.frame(xu, xl, xz)
}
