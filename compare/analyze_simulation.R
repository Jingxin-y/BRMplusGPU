
### input param, event, hypothesis, n, and R used in 'run_simulation'

library(reshape2)
library(ggplot2)
library(ggsci)
library(dplyr)
library(tidyr)
library(matrixStats)

est.result <- function(df,para.true){
  est <-  df$estimate
  bias <- mean(est)-para.true
  
  se <- df$se
  se.est<- mean(se)/sqrt(length(se))
  
  sd.est<- mean(se)
  sd.mc <- sd(est)
  acc <- sd.est/sd.mc
  
  low <- df$low
  up <- df$up
  cov <- mean((low<para.true)*(up>para.true))
  
  p <- mean(df$p<=0.05)
  
  estm.ml <- c(bias,se.est,acc,cov,p)
  return(estm.ml)
}

plot_metric <- function(df, metric, title, value, ymax = 1, file = NULL) {
  # 
  df <- df %>%
    mutate(y_plot = pmin(!!sym(metric), ymax))   
  
  p <- ggplot(df, aes(x = reorder(method, !!sym(metric)), y = y_plot)) +
    geom_col(fill = "steelblue") +
    geom_hline(yintercept = value, linetype = "dashed", color = "grey60") +
    
    
    geom_text(
      data = subset(df, !!sym(metric) > ymax),
      aes(label = sprintf("%.2f", !!sym(metric)), y = ymax),
      vjust = -0.3, size = 3.5, color = "red"
    ) +
    
    coord_cartesian(ylim = c(0, ymax * 1.05)) +   
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) +
    labs(x = "Method", y = metric, title = title)
  
  if (!is.null(file))
    ggsave(file, p, width = 8, height = 6, dpi = 300)
  return(p)
}

## read data
data <- read.csv(paste0("simulation_results_",param,"_",event,"_",hypothesis,"_n_", n, "_R_", R,".csv"))


if(param == "RR"){
  data.brm <- data[(seq_len(nrow(data)) %% 14) == 1,]
  data.brm.ad <- data[(seq_len(nrow(data)) %% 14) == 2,]
  data.CMH <- data[(seq_len(nrow(data)) %% 14) == 3,]
  data.lb <- data[(seq_len(nrow(data)) %% 14) == 4,]
  data.lp <- data[(seq_len(nrow(data)) %% 14) == 5,]
  data.rlp <- data[(seq_len(nrow(data)) %% 14) == 6,]
  data.firth <- data[(seq_len(nrow(data)) %% 14) == 7,]
  data.exact <- data[(seq_len(nrow(data)) %% 14) == 8,]
  data.exact.ad <- data[(seq_len(nrow(data)) %% 14) == 9,]
  data.GC <- data[(seq_len(nrow(data)) %% 14) == 10,]
  data.GC.BR <- data[(seq_len(nrow(data)) %% 14) == 11,]
  data.GC.FC <- data[(seq_len(nrow(data)) %% 14) == 12,]
  data.GC.FC.BR1 <- data[(seq_len(nrow(data)) %% 14) == 13,]
  data.GC.FC.BR2 <- data[(seq_len(nrow(data)) %% 14) == 0,]
  num_cols <- sapply(data.CMH, is.numeric)
  df.CMH <- data.CMH[ apply(data.CMH[ , num_cols], 1, function(x) all(is.finite(x))), ]
  num_cols <- sapply(data.GC, is.numeric)
  df.GC <- data.GC[ apply(data.GC[ , num_cols], 1, function(x) all(is.finite(x))), ]
  num_cols <- sapply(data.GC.BR, is.numeric)
  df.GC.BR <- data.GC.BR[ apply(data.GC.BR[ , num_cols], 1, function(x) all(is.finite(x))), ]
}else{
  data.brm <- data[(seq_len(nrow(data)) %% 14) == 1,]
  data.brm.ad <- data[(seq_len(nrow(data)) %% 14) == 2,]
  data.bayes <- data[(seq_len(nrow(data)) %% 14) == 3,]
  data.glm <- data[(seq_len(nrow(data)) %% 14) == 4,]
  data.lpm <- data[(seq_len(nrow(data)) %% 14) == 5,]
  data.MN <- data[(seq_len(nrow(data)) %% 14) == 6,]
  data.firth <- data[(seq_len(nrow(data)) %% 14) == 7,]
  data.exact <- data[(seq_len(nrow(data)) %% 14) == 8,]
  data.exact.ad <- data[(seq_len(nrow(data)) %% 14) == 9,]
  data.GC <- data[(seq_len(nrow(data)) %% 14) == 10,]
  data.GC.BR <- data[(seq_len(nrow(data)) %% 14) == 11,]
  data.GC.FC <- data[(seq_len(nrow(data)) %% 14) == 12,]
  data.GC.FC.BR1 <- data[(seq_len(nrow(data)) %% 14) == 13,]
  data.GC.FC.BR2 <- data[(seq_len(nrow(data)) %% 14) == 0,]
  
  num_cols <- sapply(data.GC, is.numeric)
  df.GC <- data.GC[ apply(data.GC[ , num_cols], 1, function(x) all(is.finite(x))), ]
  num_cols <- sapply(data.GC.BR, is.numeric)
  df.GC.BR <- data.GC.BR[ apply(data.GC.BR[ , num_cols], 1, function(x) all(is.finite(x))), ]
  num_cols <- sapply(data.lpm, is.numeric)
  df.lpm <- data.lpm[ apply(data.lpm[ , num_cols], 1, function(x) all(is.finite(x))), ]
}

## true value
if(param == "RR"){
  if (event == "common"){
    if (hypothesis == "null"){
      alpha.true <- 0
      beta.true  <- c(1.5, 0.6)
      gamma.true <- c(0.2, -0.5)
    }else{
      alpha.true <- 0.3
      beta.true  <- c(1.65, 0.5)
      gamma.true <- c(0.2, -0.5)
    }
  }else{
    if (hypothesis == "null"){
      alpha.true <- 0
      beta.true  <- c(-4.7, 0.5)
      gamma.true <- c(0.2, -0.5)
    }else{
      alpha.true <- 0.7
      beta.true  <- c(-5.5, 0.5)
      gamma.true <- c(0.2, -0.5)
    }
  }
}else{
  if (event == "common"){
    if (hypothesis == "null"){
      alpha.true = 0
      beta.true   = c(0.9,0.5)
      gamma.true  = c(0.2,-0.5)
    }else{
      alpha.true = 0.1
      beta.true   = c(0.9,0.2)
      gamma.true  = c(0.2,-0.5)
    }
  }else{
    if (hypothesis == "null"){
      alpha.true = 0
      beta.true   = c(-4.5,0.5)
      gamma.true  = c(0.2,-0.5)
    }else{
      alpha.true = 0.05
      beta.true   = c(-5.5,0.2)
      gamma.true  = c(0.2,-0.5)# rare
    }
  }
}

## results 
if(param == "RR"){
  result.brm <- est.result(data.brm,alpha.true)
  result.brm_b <- est.result(data.brm.ad,alpha.true)
  result.CMH <- est.result(df.CMH,alpha.true)
  result.LB <- est.result(data.lb,alpha.true)
  result.LP <- est.result(data.lp,alpha.true)
  result.RLP <- est.result(data.rlp,alpha.true)
  result.brm.FC <- est.result(data.firth,alpha.true)
  result.brm.BC <- est.result(data.exact,alpha.true)
  result.brm_b.BC <- est.result(data.exact.ad,alpha.true)
  result.GC <- est.result(data.GC,alpha.true)
  result.GC.BR <- est.result(data.GC.BR,alpha.true)
  result.GC.FC <- est.result(data.GC.FC,alpha.true)
  result.GC.FC.BR1 <- est.result(data.GC.FC.BR1,alpha.true)
  result.GC.FC.BR2 <- est.result(data.GC.FC.BR2,alpha.true)
  
  result <- cbind(result.brm, result.CMH, result.LB, result.LP,
                  result.RLP, result.brm.FC, result.brm.BC,result.brm_b,
                  result.brm_b.BC, result.GC,result.GC.BR,result.GC.FC
                  ,result.GC.FC.BR1,result.GC.FC.BR2)
  
  rownames(result) <- c("bias", "se", "acc", "coverage", "p")
  colnames(result) <- c("brm","CMH","LB","LP","RLP","brm-FC","brm-BC",
                        "brm_b","brm_b-BC","GC","GC-BR","GC-FC","GC-FC-BR1","GC-FC-BR2")
}else{
  result.brm <- est.result(data.brm,alpha.true)
  result.brm_b <- est.result(data.brm.ad,alpha.true)
  result.bayesian <- est.result(data.bayes,alpha.true)
  result.GLM <- est.result(data.glm,alpha.true)
  result.LPM <- est.result(df.lpm,alpha.true)
  result.MN<- est.result(data.MN,alpha.true)
  result.brm.FC <- est.result(data.firth,alpha.true)
  result.brm.BC <- est.result(data.exact,alpha.true)
  result.brm_b.BC <- est.result(data.exact.ad,alpha.true)
  result.GC <- est.result(df.GC,alpha.true)
  result.GC.BR <- est.result(df.GC.BR,alpha.true)
  result.GC.FC <- est.result(data.GC.FC,alpha.true)
  result.GC.FC.BR1 <- est.result(data.GC.FC.BR1,alpha.true)
  result.GC.FC.BR2 <- est.result(data.GC.FC.BR2,alpha.true)
  
  result <- cbind(result.brm,result.bayesian, result.GLM, result.LPM,
                  result.MN, result.brm.FC,result.brm.BC,result.brm_b,result.brm_b.BC,
                  result.GC,result.GC.BR,result.GC.FC,result.GC.FC.BR1,result.GC.FC.BR2)
  
  rownames(result) <- c("bias", "se", "acc", "coverage", "p")
}

write.csv(result, paste0("result_", param, "_", n, "_", event, "_",hypothesis,".csv"))
## plots

if (param == "RR"){
  data_names <- c("data.brm", "data.CMH", "data.lb", "data.lp",
                  "data.rlp", "data.firth", "data.exact","data.brm.ad",
                  "data.exact.ad", "data.GC","data.GC.BR","data.GC.FC"
                  ,"data.GC.FC.BR1","data.GC.FC.BR2")
  
  data_list <- mget(data_names)
  
  est <- do.call(cbind, lapply(data_list, `[[`, "estimate"))
  colnames(est) <- c("brm","CMH","LB","LP","RLP","brm-FC","brm-BC",
                     "brm_b","brm_b-BC","GC","GC-BR","GC-FC","GC-FC-BR1","GC-FC-BR2")
  
  se <- do.call(cbind, lapply(data_list, `[[`, "se"))
  colnames(se) <- c("brm","CMH","LB","LP","RLP","brm-FC","brm-BC",
                    "brm_b","brm_b-BC","GC","GC-BR","GC-FC","GC-FC-BR1","GC-FC-BR2")
  
}else{
  data_names <- c("data.brm", "data.brm.ad", "data.bayes", "data.glm",
                  "data.lpm", "data.MN", "data.firth", "data.exact",
                  "data.exact.ad", "data.GC", "data.GC.BR" , "data.GC.FC",
                  "data.GC.FC.BR1","data.GC.FC.BR2")
  data_list <- mget(data_names)
  
  est <- do.call(cbind, lapply(data_list, `[[`, "estimate"))
  colnames(est) <- c("brm","brm_b","bayesian","GLM","LPM","MN","brm-FC","brm-BC",
                     "brm_b-BC","GC","GC-BR","GC-FC","GC-FC-BR1","GC-FC-BR2")
  
  se <- do.call(cbind, lapply(data_list, `[[`, "se"))
  colnames(se) <- c("brm","brm_b","bayesian","GLM","LPM","MN","brm-FC","brm-BC",
                     "brm_b-BC","GC","GC-BR","GC-FC","GC-FC-BR1","GC-FC-BR2")
}


est_df <- as.data.frame(est) 
est_long <- melt(est, variable.name = "method", value.name = "estimate")
colnames(est_long) <- c("number","method","estimate")

removed_counts <- est_long %>%
  group_by(method) %>%
  summarise(removed = sum(is.na(estimate) | abs(estimate) >= 5))

est_long_small <- est_long %>% filter(!is.na(estimate), abs(estimate) < 5)


p1 = ggplot(est_long_small, aes(x = method, y = estimate, fill = method)) +
  #geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(color = "grey30",width = 0.6, outlier.size = 0.5, alpha = 0.9) +
  ggsci::scale_fill_d3(palette = "category20") +
  theme_minimal(base_size = 14) +
  geom_hline(yintercept = alpha.true, linetype = "dashed", color = "steelblue") + 
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1.5, fill = "white", color = "black") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = paste0("Monte Carlo when n = ", n),
       x = "Method", y = "Estimate")

stats <- est_long_small %>%
  group_by(method) %>%
  summarise(ypos = max(estimate, na.rm = TRUE))


removed_counts <- left_join(removed_counts, stats, by = "method")


p2 = p1 + geom_text(data = removed_counts,
                    aes(x = method, y = ypos + 0.2, 
                        label = paste0("(", removed,")")),
                    inherit.aes = FALSE,
                    size = 3.5, color = "steelblue")

if(sum(removed_counts$removed)==0){
  ggsave(paste0("est_",param,"_",event,"_",hypothesis,"_n_", n,".png"), p1, width = 12, height = 6, dpi = 300)
}else{
  ggsave(paste0("est_",param,"_",event,"_",hypothesis,"_n_", n,".png"), p2, width = 12, height = 6, dpi = 300)
}

### plot of SE

se_df <- as.data.frame(se) 
se_long <- melt(se, variable.name = "method", value.name = "SE")
colnames(se_long) <- c("number","method","SE")

removed_counts <- se_long %>%
  group_by(method) %>%
  summarise(removed = sum(is.na(SE) | SE >= 5))


se_long_small <- se_long %>% filter(!is.na(SE), SE < 5)

p3 <- ggplot(se_long_small, aes(x = method, y = SE, fill = method)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = paste0("Boxplot of RR for n = ", n),
       x = "Method", y = "SE")

stats <- se_long_small %>%
  group_by(method) %>%
  summarise(ypos = max(SE, na.rm = TRUE))

removed_counts <- left_join(removed_counts, stats, by = "method")

p4 = p3 + geom_text(data = removed_counts,
                   aes(x = method, y = ypos + 0.2, 
                       label = paste0("(", removed,")")),
                   inherit.aes = FALSE,
                   size = 3.5, color = "steelblue")

if(sum(removed_counts$removed)==0){
  ggsave(paste0("se_",param,"_",event,"_",hypothesis,"_n_", n,".png"), p3, width = 12, height = 6, dpi = 300)
}else{
  ggsave(paste0("se_",param,"_",event,"_",hypothesis,"_n_", n,".png"), p4, width = 12, height = 6, dpi = 300)
}

### bar for accuracy, coverage, and p-value
df <- as.data.frame(t(result))
df$method <- sub("^result\\.", "", rownames(df))
rownames(df) <- NULL

num_cols <- c("bias","se","acc","coverage","p")
df <- df %>% mutate(across(all_of(num_cols), as.numeric))

p_acc <- plot_metric(df, "acc",      paste0("Barplot of accuarcy at n = ",n), 1, 1.5,
                     paste0("accuracy ", param," ", event," ",hypothesis, " n = ",n, ".png"))
p_cov <- plot_metric(df, "coverage", paste0("Barplot of coverage at n = ",n), 0.95,1.1,
                     paste0("coverage ", param," ", event," ",hypothesis, " n = ",n, ".png"))


