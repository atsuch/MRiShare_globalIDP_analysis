library(stats)
library(ggplot2)
library(ggpubr)
library(heplots)

# Set working directory. Change this according to where your folder is.
setwd("PATH_TO_WD")

# Data for all 1831* subjects (Note that one subject withdrew from the study since
# the writing of associated publication, and the data had to be dropped)
master_dat <- read.table("data/MRiShare_global_IDPs_BSAF2021.csv", header=TRUE, sep=",", dec=".")
head(master_dat)
summary(master_dat)


#### Analysis plan ###########################################################
# Since only 5% of our data are >26 yrs old, we will
# perform and report analyses for those under <26 so that we can focus 
# interpretation of the results for those between 18 to 26.
#
# In the following analyses, we will test the effect of age separately for
# males and females, followed by the combined analyses to also test the effect
# of sex.
#
# In each case, both linear and quadratic age effects are fitted, with or w/o
# eTIV or WM mask vol for DTI/NODDI metrics.
#
# The best model is evaluated using BIC (but AIC is also 
# computed for comparison). We will focus on linear age effect models, but 
# BIC/AIC for quadratic models are also computed.
##############################################################################
#### Data ####

## Create a list of df for those under 26
u26_dat <- list("male" = master_dat[which(master_dat$Age <= 26 & master_dat$Sex == 'M'), ],
                "female" = master_dat[which(master_dat$Age <= 26 & master_dat$Sex == 'F'), ],
                "both" = master_dat[which(master_dat$Age <= 26), ])

## Add centered age, squared age and eTIV
u26_dat <- lapply(u26_dat, function(df) {
  df$Age_c = as.numeric(scale(df$Age, scale=F))
  df$Sq_Age_c = as.numeric(df$Age_c^2)
  df$eTIV_c = as.numeric(scale(df$eTIV, scale=F))
  return(df)
})

## For DTI/NODDI, some subjects are missing data. Create dtinoddi_u26_dat.
# We also need to recompute Age_c, Sq_Age_c in each df
dtinoddi_dat <- lapply(u26_dat, function(df) {
  new_df <- df[!is.na(df$cerebralWM_mean_FA), !(names(df) %in% c("Age_c", "eTIV_c", "Sq_Age_c"))]
  new_df$Age_c = as.numeric(scale(new_df$Age, scale=F))
  new_df$Sq_Age_c = as.numeric(new_df$Age_c^2)
  new_df$WMmask_c = as.numeric(scale(new_df$FS6_cerebralWM_mask_vol, scale=F))
  return(new_df)
})


## Create df for prediction that contain hypothetical data for male/female,
## with age range of 18-26 and constant eTIV or WMmask volume (mean of the group)
u26_pred_dat <- lapply(u26_dat, function(df) {
  pred_df <- data.frame(Age=seq(18.5, 25.5, 0.5),
                        eTIV_c=mean(df$eTIV_c))
  pred_df$Age_c = pred_df$Age - mean(df$Age)
  pred_df$Sq_Age_c = pred_df$Age_c^2
  if (nrow(unique(df["Sex"])) == 2) {
    pred_df_m <- data.frame(pred_df)
    pred_df_m$Sex = factor("M", levels=c("F", "M"))
    pred_df_f <- data.frame(pred_df)
    pred_df_f$Sex = factor("F", levels=c("F", "M"))
    pred_df <- rbind(pred_df_m, pred_df_f)
  } else {
    pred_df$Sex = factor(unique(df["Sex"])[[1]], levels=c("F", "M"))
  }
  pred_df
})

dtinoddi_pred_dat <- lapply(dtinoddi_dat, function(df) {
  pred_df <- data.frame(Age=seq(18.5, 25.5, 0.5),
                        WMmask_c=mean(df$WMmask_c))
  pred_df$Age_c = pred_df$Age - mean(df$Age)
  pred_df$Sq_Age_c = pred_df$Age_c^2
  if (nrow(unique(df["Sex"])) == 2) {
    pred_df_m <- data.frame(pred_df)
    pred_df_m$Sex = factor("M", levels=c("F", "M"))
    pred_df_f <- data.frame(pred_df)
    pred_df_f$Sex = factor("F", levels=c("F", "M"))
    pred_df <- rbind(pred_df_m, pred_df_f)
  } else {
    pred_df$Sex = factor(unique(df["Sex"])[[1]], levels=c("F", "M"))
  }
  pred_df
})


## Save outputs in one folder
dir.create("Results_test")
setwd("Result_test/")

################################################################################
### Functions that will be used in this script
## Return outlier upper and lower values
getOutlierVals <- function(data, iqr_thr=1.5) {
  lowerq = unname(quantile(data)[2])
  upperq = unname(quantile(data)[4])
  iqr = upperq - lowerq #Or use IQR(data)
  threshold.upper = (iqr * iqr_thr) + upperq
  threshold.lower = lowerq - (iqr * iqr_thr)
  return(list("upper"=threshold.upper, "lower"=threshold.lower))
}

## Get numbe of outliers for a given data 
numOutliers <- function(data, iqr_thr=1.5) {
  thr <- getOutlierVals(data, iqr_thr)
  result <- length(which(data > thr$upper | data < thr$lower))
}

## Write descriptive stats for a given variables
writeDStxt <- function(df_list, column, fname) {
  sink(fname)
  cat("Descriptive stats for", column, "\n")
  cat("\nMale stats:\n")
  cat("N:", nrow(df_list$male), "\n")
  cat("Mean:", mean(df_list$male[, column]), "\n")
  cat("SD:", sd(df_list$male[, column]), "\n")
  cat("Range:", range(df_list$male[, column]), "\n")
  cat("\nFemale stats:\n")
  cat("N:", nrow(df_list$female), "\n")
  cat("Mean:", mean(df_list$female[, column]), "\n")
  cat("SD:", sd(df_list$female[, column]), "\n")
  cat("Range:", range(df_list$female[, column]), "\n")
  cat("\nCombined stats:\n")
  cat("N:", nrow(df_list$both), "\n")
  cat("Mean:", mean(df_list$both[, column]), "\n")
  cat("SD:", sd(df_list$both[, column]), "\n")
  cat("Range:", range(df_list$both[, column]), "\n")
  sink()
}

## Get all linear models for male, female, combined data, linear/
# quadratic age model, with or without eTIV
getLMmodels <- function(df_list, metric_name) {
  # Sex specific models
  s_lin_noVol_model_eq <- substitute(n ~ Age_c,
                                      list(n=as.name(metric_name)))
  lin_noVol_m <- lm(s_lin_noVol_model_eq, data=df_list$male)
  lin_noVol_f <- lm(s_lin_noVol_model_eq, data=df_list$female)
  
  s_quad_noVol_model_eq <- substitute(n ~ Age_c + Sq_Age_c,
                                      list(n=as.name(metric_name)))
  quad_noVol_m <- lm(s_quad_noVol_model_eq, data=df_list$male)
  quad_noVol_f <- lm(s_quad_noVol_model_eq, data=df_list$female)
  
  # Combined group models
  b_lin_noVol_model_eq <- substitute(n ~ Sex + Age_c + Sex*Age_c,
                                      list(n=as.name(metric_name)))
  lin_noVol_b <- lm(b_lin_noVol_model_eq, data=df_list$both, contrasts = list(Sex="contr.sum"))
  lin_noVol_b_cdef <- lm(b_lin_noVol_model_eq, data=df_list$both)
  
  b_quad_noVol_model_eq <- substitute(n ~ Sex + Age_c + Sq_Age_c +
                                         Sex*Age_c + Sex*Sq_Age_c,
                                      list(n=as.name(metric_name)))
  quad_noVol_b <- lm(b_quad_noVol_model_eq, data=df_list$both, contrasts = list(Sex="contr.sum"))
  quad_noVol_b_cdef <- lm(b_quad_noVol_model_eq, data=df_list$both)
  
  if (grepl("cerebralWM", metric_name, fixed=TRUE) | metric_name != "eTIV") {
    if (grepl("cerebralWM", metric_name, fixed=TRUE)) {
      s_lin_Vol_model_eq <- substitute(n ~ Age_c + WMmask_c,
                                      list(n=as.name(metric_name)))
      s_quad_Vol_model_eq <- substitute(n ~ Age_c + Sq_Age_c + WMmask_c,
                                       list(n=as.name(metric_name)))
      b_lin_Vol_model_eq <- substitute(n ~ Sex + Age_c + WMmask_c +
                                        Sex*Age_c + Sex*WMmask_c,
                                      list(n=as.name(metric_name)))
      b_quad_Vol_model_eq <- substitute(n ~ Sex + Age_c + Sq_Age_c + WMmask_c +
                                         Sex*Age_c + Sex*Sq_Age_c + Sex*WMmask_c,
                                       list(n=as.name(metric_name)))
    } else if (metric_name != "eTIV") {
      s_lin_Vol_model_eq <- substitute(n ~ Age_c + eTIV_c,
                                     list(n=as.name(metric_name)))
      s_quad_Vol_model_eq <- substitute(n ~ Age_c + Sq_Age_c + eTIV_c,
                                      list(n=as.name(metric_name)))
      b_lin_Vol_model_eq <- substitute(n ~ Sex + Age_c + eTIV_c +
                                       Sex*Age_c + Sex*eTIV_c,
                                     list(n=as.name(metric_name)))
      b_quad_Vol_model_eq <- substitute(n ~ Sex + Age_c + Sq_Age_c + eTIV_c +
                                        Sex*Age_c + Sex*Sq_Age_c + Sex*eTIV_c,
                                      list(n=as.name(metric_name)))
    }
    
    lin_Vol_m <- lm(s_lin_Vol_model_eq, data=df_list$male)
    lin_Vol_f <- lm(s_lin_Vol_model_eq, data=df_list$female)

    quad_Vol_m <- lm(s_quad_Vol_model_eq, data=df_list$male)
    quad_Vol_f <- lm(s_quad_Vol_model_eq, data=df_list$female)
    
    lin_Vol_b <- lm(b_lin_Vol_model_eq, data=df_list$both, contrasts = list(Sex="contr.sum"))
    lin_Vol_b_cdef <- lm(b_lin_Vol_model_eq, data=df_list$both)

    quad_Vol_b <- lm(b_quad_Vol_model_eq, data=df_list$both, contrasts = list(Sex="contr.sum"))
    quad_Vol_b_cdef <- lm(b_quad_Vol_model_eq, data=df_list$both)
  }
  
  lm_models <- list("lin_noVol_m"=lin_noVol_m, "lin_noVol_f"=lin_noVol_f,
                    "quad_noVol_m"=quad_noVol_m, "quad_noVol_f"=quad_noVol_f,
                    "lin_noVol_b"=lin_noVol_b, "lin_noVol_b_cdef"=lin_noVol_b_cdef,
                    "quad_noVol_b"=quad_noVol_b, "quad_noVol_b_cdef"=quad_noVol_b_cdef)
  
  if (grepl("cerebralWM", metric_name, fixed=TRUE)) {
    maskVol_models <- list("lin_maskVol_m"=lin_Vol_m, "lin_maskVol_f"=lin_Vol_f,
                           "quad_maskVol_m"=quad_Vol_m, "quad_maskVol_f"=quad_Vol_f,
                           "lin_maskVol_b"=lin_Vol_b, "lin_maskVol_b_cdef"=lin_Vol_b_cdef,
                           "quad_maskVol_b"=quad_Vol_b, "quad_maskVol_b_cdef"=quad_Vol_b_cdef)
    lm_models <- append(lm_models, maskVol_models)
  } else if  (metric_name != "eTIV") {
    eTIV_models <- list("lin_eTIV_m"=lin_Vol_m, "lin_eTIV_f"=lin_Vol_f,
                        "quad_eTIV_m"=quad_Vol_m, "quad_eTIV_f"=quad_Vol_f,
                        "lin_eTIV_b"=lin_Vol_b, "lin_eTIV_b_cdef"=lin_Vol_b_cdef,
                        "quad_eTIV_b"=quad_Vol_b, "quad_eTIV_b_cdef"=quad_Vol_b_cdef)
    lm_models <- append(lm_models, eTIV_models)
  }
  return(lm_models)
}

## Write male, female, or combined results
writeLMresults <- function(lm_models, df_list, dname, metric_name, sex) {
  fname = sprintf("%s/%s_lm_results_%s.txt", dname, metric_name, sex)
  n = ifelse(sex=="m", nrow(df_list$male),
             ifelse(sex=="f", nrow(df_list$female), nrow(df_list$both))
             )
  models = c(paste0(c("lin_noVol_", "quad_noVol_"), sex))
  
  if (grepl("cerebralWM", metric_name, fixed=TRUE)) {
    maskVol_models = c(paste0(c("lin_maskVol_", "quad_maskVol_"), sex))
    models = append(models, maskVol_models)
  } else if (metric_name !="eTIV") {
    eTIV_models = c(paste0(c("lin_eTIV_", "quad_eTIV_"), sex))
    models = append(models, eTIV_models)
  }

  sink(fname)
  cat("Model comparison: \n")
  for (m in models) {
    cat("AIC for", m, ":", extractAIC(lm_models[[m]]), "\n")
  }
  for (m in models) {
    cat("BIC for", m, ":", extractAIC(lm_models[[m]], k=log(n)), "\n")
  }
  cat("\nResults: \n")
  for (m in models) {
    cat("\nmodel:", m, "\n")
    print(summary(lm_models[[m]]))
    print(confint(lm_models[[m]]))
    print(etasq(lm_models[[m]]))
    if (sex=="b") {
      print(summary(lm_models[[paste0(m, "_cdef")]]))
      print(confint(lm_models[[paste0(m, "_cdef")]]))
      print(etasq(lm_models[[paste0(m, "_cdef")]]))
    }
  }
  sink()
}

## Save all LM model results
saveAllLMmodels <- function(lm_models, df_list, dname, metric_name) {
  writeLMresults(lm_models, df_list, dname, metric_name, sex="m")
  writeLMresults(lm_models, df_list, dname, metric_name, sex="f")
  writeLMresults(lm_models, df_list, dname, metric_name, sex="b")
}

## Return prediction df for a given sex
createPredDF <- function(yvar, lm_models, pred_df_list, model_name, sex) {
  df_name <- ifelse(sex=="m", "male", ifelse(sex=="f", "female", "both"))
  pred_df = pred_df_list[[df_name]]
  pred_df$model_name = model_name
  conf_int <- predict(lm_models[[paste0(model_name, "_", sex)]],
                      pred_df, interval="confidence", levels=0.95)
  pred_df[yvar] = conf_int[,1]
  pred_df$LCI = conf_int[,2]
  pred_df$HCI = conf_int[,3]
  pred_int <- predict(lm_models[[paste0(model_name, "_", sex)]],
                      pred_df, interval="prediction", levels=0.95)
  pred_df$LPI = pred_int[,2]
  pred_df$HPI = pred_int[,3]

  return(pred_df)
}

## Get common Ylim for plotting male, female, and combined data
getYlim <- function(yvar, delta=0.25) {
  b_min = min(u26_dat$both[yvar], na.rm = TRUE)
  b_max = max(u26_dat$both[yvar], na.rm = TRUE)
  b_range = b_max - b_min
  b_delta = b_range * delta
  b_ylim = c((b_min - b_delta), (b_max + b_delta))
  return(b_ylim)
}

## Plot function used to generate figures in the paper.
# It shows male and female plots separately, and without too many text elements.
genSimpleAgeEffectPlots <- function(lm_models, df_list_name,
                                   main_model_name,
                                   yvar, ylab,
                                   scaling_factor=1) {
  # DF to make prediction
  pred_df_list <- get(gsub("_dat", "_pred_dat", df_list_name))
  
  # DF with predicted results
  m_model_name <- ifelse(is.character(main_model_name), main_model_name, main_model_name$male)
  f_model_name <- ifelse(is.character(main_model_name), main_model_name, main_model_name$female)

  m_pred_df <- createPredDF(yvar, lm_models, pred_df_list, m_model_name, sex="m")
  f_pred_df <- createPredDF(yvar, lm_models, pred_df_list, f_model_name, sex="f")
  
  # now plot
  yvar_lim <- getYlim(yvar, delta=0.1)*scaling_factor
  # 1) male plot
  m_p <- ggplot(data=get(df_list_name)$male, aes(x=Age, y=get(yvar)*scaling_factor, color=Sex)) + 
    geom_point() +
    theme_linedraw() +
    scale_color_manual(values=c(adjustcolor("blue", alpha.f = 0.5)), 
                       breaks=c("M"), labels=c("Male data")) +
    ylim(yvar_lim) +
    geom_ribbon(data=m_pred_df, aes(ymin=LCI*scaling_factor, ymax=HCI*scaling_factor), color=NA,
                fill=adjustcolor('grey', alpha.f = 0.5)) +
    geom_line(color="black", data=m_pred_df, aes(x=Age, y=get(yvar)*scaling_factor), size=1.1) +
    rremove("xylab") +
    theme(legend.position="none",
          legend.title=element_blank(),
          axis.text=element_text(size=10))
  
  # 2) female plot
  f_p <- ggplot(data=get(df_list_name)$female, aes(x=Age, y=get(yvar)*scaling_factor, color=Sex)) + 
    geom_point() +
    theme_linedraw() +
    scale_color_manual(values=c(adjustcolor("red", alpha.f = 0.5)), 
                       breaks=c("F"), labels=c("Female data")) +
    ylim(yvar_lim) +
    geom_ribbon(data=f_pred_df, aes(ymin=LCI*scaling_factor, ymax=HCI*scaling_factor), color=NA,
                fill=adjustcolor('grey', alpha.f = 0.5)) +
    geom_line(color="black", data=f_pred_df, aes(x=Age, y=get(yvar)*scaling_factor), size=1.1) +
    rremove("xylab") +
    theme(legend.position="none",
          legend.title=element_blank(),
          axis.text=element_text(size=10))
  
  # Combine them all
  fig <- ggarrange(m_p, f_p, nrow=1)
  fig<- annotate_figure(fig,
                        left=text_grob(ylab, size=12, rot=90),
                        bottom=text_grob("Age (years)", size=12),
                        fig.lab.size=12)
  return(fig)
}

################################################################################

#######################################################
# Basic demographic summary of the data
#######################################################
dir.create("demography")

# First save the age distribution in entire sample

# Age distribution figure with upper age line (95 percentile)
age_upper <- quantile(master_dat$Age, probs=c(0.95)) 
age_dist <- ggplot(master_dat, aes(x=Age, fill=Sex, color=Sex)) +
  geom_histogram(alpha=0.5, position="identity", bins = 100) +
  geom_vline(xintercept = age_upper, linetype="dashed", color="black", size=0.8) +
  labs(x="Age (years)", y="Number of subjects") +
  scale_color_manual(values=c("red", "blue"),
                     breaks=c("F", "M"), labels=c("Female", "Male")) +
  scale_fill_manual(values=c(adjustcolor("red", alpha.f = 0.5),
                             adjustcolor("blue", alpha.f = 0.5)), 
                     breaks=c("F", "M"), labels=c("Female", "Male")) +
  scale_x_continuous(breaks=seq(18,36,2)) +
  theme_linedraw() +
  theme(legend.position="none",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))

ggsave("demography/MRiShare_age_distribution.png",
       age_dist,
       device="png",
       width=10,
       height=7,
       units="cm",
       dpi=300)

# Save age summary
writeDStxt(u26_dat, "Age", "demography/u26_Age_descriptive_stats.txt")

## Check if there is sex difference in age
# first check if variance is the same in two groups 
sex_age_var <- var.test(Age ~ Sex, data = u26_dat$both)
sex_age_var
# No significant difference in variance
sex_age_t <- t.test(Age ~ Sex, data = u26_dat$both, var.equal=TRUE)
sex_age_t
# difference is marginal (p=0.07906)


#######################################################
###  0) eTIV
### Check whether eTIV is changing at this age range
#######################################################
dir.create("lm0_eTIV")

# First save descriptive stats
writeDStxt(u26_dat, "eTIV", "lm0_eTIV/u26_eTIV_descriptive_stats.txt")

# Create models and save the results
lm0_eTIV_models <- getLMmodels(u26_dat, metric_name = "eTIV")
saveAllLMmodels(lm0_eTIV_models, u26_dat, "lm0_eTIV", metric_name = "eTIV")

## No evidence for eTIV change at this age range

#######################################################
###  1) Cortical thickness
#######################################################
dir.create("lm1_CT")

# First save descriptive stats
writeDStxt(u26_dat, "FS6_mean_CT", "lm1_CT/u26_CT_descriptive_stats.txt")

# Check for outliers
CT_mo = numOutliers(u26_dat$both$FS6_mean_CT)
CT_eo = numOutliers(u26_dat$both$FS6_mean_CT, iqr_thr = 3.0)

# Create models and save the results
lm1_CT_models <- getLMmodels(u26_dat, metric_name = "FS6_mean_CT")
saveAllLMmodels(lm1_CT_models, u26_dat, "lm1_CT", metric_name = "FS6_mean_CT")

# Plot
CT_fig <- genSimpleAgeEffectPlots(lm1_CT_models, "u26_dat", "lin_noVol",
                                   yvar="FS6_mean_CT",
                                   ylab="CT (mm) ")
ggsave("lm1_CT/CT_vs_Age.png",
       CT_fig,
       device="png",
       width=20,
       height=7,
       units="cm",
       dpi=300)


#######################################################
###  2) White and pial surface area
#######################################################
dir.create("lm2_CSA")

# First save descriptive stats
writeDStxt(u26_dat, "FS6_total_innerCSA", "lm2_CSA/u26_innerCSA_descriptive_stats.txt")
writeDStxt(u26_dat, "FS6_total_pialCSA", "lm2_CSA/u26_pialCSA_descriptive_stats.txt")

# Check for outliers
innerCSA_mo = numOutliers(u26_dat$both$FS6_total_innerCSA)
innerCSA_eo = numOutliers(u26_dat$both$FS6_total_innerCSA, iqr_thr = 3.0)

pialCSA_mo = numOutliers(u26_dat$both$FS6_total_pialCSA)
pialCSA_eo = numOutliers(u26_dat$both$FS6_total_pialCSA, iqr_thr = 3.0)

# Check how much inner vs pial CSA correlate. 
CSA_cor <- cor(u26_dat$both$FS6_total_innerCSA, u26_dat$both$FS6_total_pialCSA)
# cor=0.98

# Plot inner vs pial CSA #
CSA_min = min(c(getYlim("FS6_total_innerCSA", 0.1)[1], getYlim("FS6_total_pialCSA", 0.1)[1]))*10^-3
CSA_max = max(c(getYlim("FS6_total_innerCSA", 0.1)[2], getYlim("FS6_total_pialCSA", 0.1)[2]))*10^-3
inner_vs_pial <- ggplot(u26_dat$both, aes(x=FS6_total_innerCSA*10^-3, y=FS6_total_pialCSA*10^-3, color=Sex)) +
  geom_point() +
  xlim(c(CSA_min, CSA_max)) +
  ylim(c(CSA_min, CSA_max)) +
  geom_abline(slope=1, intercept=0, linetype=2) +
  scale_color_manual(values=c(adjustcolor("red", alpha.f = 0.5), adjustcolor("blue", alpha.f = 0.5)),
                     breaks=c("F", "M"), labels=c("Female", "Male")) +
  labs(x=expression("inner CSA"~ ( x ~ 10^{3} ~ mm^{2})),
       y=expression("pial CSA"~ ( x ~ 10^{3} ~ mm^{2}))) +
  theme_linedraw()

ggsave("lm2_CSA/inner_vs_pialCSA.png",
       inner_vs_pial,
       device="png",
       width=10,
       height=7,
       units="cm",
       dpi=300)

# Although highly correlated, perform the same analyses and save the results
lm2_iCSA_models <- getLMmodels(u26_dat, metric_name = "FS6_total_innerCSA")
saveAllLMmodels(lm2_iCSA_models, u26_dat, "lm2_CSA", metric_name = "FS6_total_innerCSA")

lm2_pCSA_models <- getLMmodels(u26_dat, metric_name = "FS6_total_pialCSA")
saveAllLMmodels(lm2_pCSA_models, u26_dat, "lm2_CSA", metric_name = "FS6_total_pialCSA")

# Plot
innerCSA_fig <- genSimpleAgeEffectPlots(lm2_iCSA_models, "u26_dat", "lin_eTIV",
                                   scaling_factor = 10^-3,
                                   yvar="FS6_total_innerCSA",
                                   ylab=expression("inner CSA" ~ ( x ~ 10^{3} ~ mm^{2})))
ggsave("lm2_CSA/innerCSA_vs_Age.png",
       innerCSA_fig,
       device="png",
       width=20,
       height=7,
       units="cm",
       dpi=300)

pialCSA_fig <- genSimpleAgeEffectPlots(lm2_pCSA_models, "u26_dat", "lin_eTIV",
                                        scaling_factor = 10^-3,
                                        yvar="FS6_total_pialCSA",
                                        ylab=expression("pial CSA" ~ ( x ~ 10^{3} ~ mm^{2})))
ggsave("lm2_CSA/pialCSA_vs_Age.png",
       pialCSA_fig,
       device="png",
       width=20,
       height=7,
       units="cm",
       dpi=300)


#######################################################
###  3) SPM GM Volume
#######################################################
dir.create("lm3_spmGM")

# First save descriptive stats
writeDStxt(u26_dat, "SPM_GM_Volume", "lm3_spmGM/u26_spmGM_descriptive_stats.txt")

# Check for outliers
GM_mo = numOutliers(u26_dat$both$SPM_GM_Volume)
GM_eo = numOutliers(u26_dat$both$SPM_GM_Volume, iqr_thr = 3.0)

# Create models and save the results
lm3_spmGM_models <- getLMmodels(u26_dat, metric_name = "SPM_GM_Volume")
saveAllLMmodels(lm3_spmGM_models, u26_dat, "lm3_spmGM", metric_name = "SPM_GM_Volume")

# Plot
GM_fig <- genSimpleAgeEffectPlots(lm3_spmGM_models, "u26_dat", "lin_eTIV",
                                   scaling_factor = 10^-3,
                                   yvar="SPM_GM_Volume",
                                   ylab="GM volume (cc)")
ggsave("lm3_spmGM/GM_vs_Age.png",
       GM_fig,
       device="png",
       width=20,
       height=7,
       units="cm",
       dpi=300)


#######################################################
###  4) SPM WM Volume
#######################################################
dir.create("lm4_spmWM")

# First save descriptive stats
writeDStxt(u26_dat, "SPM_WM_Volume", "lm4_spmWM/u26_spmWM_descriptive_stats.txt")

# Check for outliers
WM_mo = numOutliers(u26_dat$both$SPM_WM_Volume)
WM_eo = numOutliers(u26_dat$both$SPM_WM_Volume, iqr_thr = 3.0)

# Create models and save the results
lm4_spmWM_models <- getLMmodels(u26_dat, metric_name = "SPM_WM_Volume")
saveAllLMmodels(lm4_spmWM_models, u26_dat, "lm4_spmWM", metric_name = "SPM_WM_Volume")

# Plot
WM_fig <- genSimpleAgeEffectPlots(lm4_spmWM_models, "u26_dat", "lin_eTIV",
                                   scaling_factor = 10^-3,
                                   yvar="SPM_WM_Volume",
                                   ylab="WM (cc)")
ggsave("lm4_spmWM/WM_vs_Age.png",
       WM_fig,
       device="png",
       width=20,
       height=7,
       units="cm",
       dpi=300)


# For supplement, plot WM and eTIV interaction #
wm_eTIV <- ggplot(u26_dat$both, aes(x=eTIV*10^-3, y=SPM_WM_Volume*10^-3, color=Sex, fill=Sex)) +
  geom_point(aes(fill=Sex), color="transparent", shape=21) +
  stat_smooth(aes(color=Sex), method=lm, se=FALSE) +
  scale_color_manual(values=c("red3", "blue3"),
                    breaks=c("F", "M"), labels=c("Female", "Male")) +
  scale_fill_manual(values=c(adjustcolor("red", alpha.f = 0.5), adjustcolor("blue", alpha.f = 0.5)),
                    breaks=c("F", "M"), labels=c("Female", "Male")) +
  labs(x=expression(paste("eTIV",
                          scriptstyle("(cc)"))),
       y=expression(paste("Total WM volume ",
                          scriptstyle("(cc)")))) +
  theme_linedraw()

ggsave("lm4_spmWM/WM_eTIV_interaction.png",
       wm_eTIV,
       device="png",
       width=10,
       height=7,
       units="cm",
       dpi=300)


#######################################################
###  5) FA
#######################################################
dir.create("lm5_FA")
names(master_dat)
# First save descriptive stats
writeDStxt(dtinoddi_dat, "cerebralWM_mean_FA", "lm5_FA/u26_FA_descriptive_stats.txt")

# Check for outliers
FA_mo = numOutliers(dtinoddi_dat$both$cerebralWM_mean_FA)
FA_eo = numOutliers(dtinoddi_dat$both$cerebralWM_mean_FA, iqr_thr = 3.0)

# Create models and save the results
lm5_FA_models <- getLMmodels(dtinoddi_dat, metric_name = "cerebralWM_mean_FA")
saveAllLMmodels(lm5_FA_models, dtinoddi_dat, "lm5_FA", metric_name = "cerebralWM_mean_FA")

# Plot
FA_fig <- genSimpleAgeEffectPlots(lm5_FA_models, "dtinoddi_dat", "lin_noVol",
                                   scaling_factor = 1,
                                   yvar="cerebralWM_mean_FA", ylab="FA")
ggsave("lm5_FA/FA_vs_Age.png",
       FA_fig,
       device="png",
       width=20,
       height=7,
       units="cm",
       dpi=300)

#######################################################
###  6) MD
#######################################################
dir.create("lm6_MD")

# First save descriptive stats
writeDStxt(dtinoddi_dat, "cerebralWM_mean_MD", "lm6_MD/u26_MD_descriptive_stats.txt")

# Check for outliers
MD_mo = numOutliers(dtinoddi_dat$both$cerebralWM_mean_MD)
MD_eo = numOutliers(dtinoddi_dat$both$cerebralWM_mean_MD, iqr_thr = 3.0)

# Create models and save the results
lm6_MD_models <- getLMmodels(dtinoddi_dat, metric_name = "cerebralWM_mean_MD")
saveAllLMmodels(lm6_MD_models, dtinoddi_dat, "lm6_MD", metric_name = "cerebralWM_mean_MD")

# Plot
MD_fig <- genSimpleAgeEffectPlots(lm6_MD_models, "dtinoddi_dat",
                                   main_model_name = list("male"="lin_maskVol", "female"="lin_noVol"),
                                   scaling_factor = 10^4,
                                   yvar="cerebralWM_mean_MD",
                                   ylab=expression("MD"~ ( x ~ 10^{-4} ~ mm^{2}/sec)))
ggsave("lm6_MD/MD_vs_Age.png",
       MD_fig,
       device="png",
       width=20,
       height=7,
       units="cm",
       dpi=300)

#######################################################
###  7) AD
#######################################################
dir.create("lm7_AD")

# First save descriptive stats
writeDStxt(dtinoddi_dat, "cerebralWM_mean_AD", "lm7_AD/u26_AD_descriptive_stats.txt")

# Check for outliers
AD_mo = numOutliers(dtinoddi_dat$both$cerebralWM_mean_AD)
AD_eo = numOutliers(dtinoddi_dat$both$cerebralWM_mean_AD, iqr_thr = 3.0)

# Create models and save the results
lm7_AD_models <- getLMmodels(dtinoddi_dat, metric_name = "cerebralWM_mean_AD")
saveAllLMmodels(lm7_AD_models, dtinoddi_dat, "lm7_AD", metric_name = "cerebralWM_mean_AD")

# Plot
AD_fig <- genSimpleAgeEffectPlots(lm7_AD_models, "dtinoddi_dat",
                                   main_model_name = "lin_noVol",
                                   scaling_factor = 10^4,
                                   yvar="cerebralWM_mean_AD",
                                   ylab=expression("AD"~ ( x ~ 10^{-4} ~ mm^{2}/sec)))
ggsave("lm7_AD/AD_vs_Age.png",
       AD_fig,
       device="png",
       width=20,
       height=7,
       units="cm",
       dpi=300)

#######################################################
###  8) RD
#######################################################
dir.create("lm8_RD")

# First save descriptive stats
writeDStxt(dtinoddi_dat, "cerebralWM_mean_RD", "lm8_RD/u26_RD_descriptive_stats.txt")

# Check for outliers
RD_mo = numOutliers(dtinoddi_dat$both$cerebralWM_mean_RD)
RD_eo = numOutliers(dtinoddi_dat$both$cerebralWM_mean_RD, iqr_thr = 3.0)

# Create models and save the results
lm8_RD_models <- getLMmodels(dtinoddi_dat, metric_name = "cerebralWM_mean_RD")
saveAllLMmodels(lm8_RD_models, dtinoddi_dat, "lm8_RD", metric_name = "cerebralWM_mean_RD")

# Plot
RD_fig <- genSimpleAgeEffectPlots(lm8_RD_models, "dtinoddi_dat",
                                   main_model_name = list("male"="lin_maskVol", "female"="lin_noVol"),
                                   scaling_factor = 10^4,
                                   yvar="cerebralWM_mean_RD",
                                   ylab=expression("RD"~ ( x ~ 10^{-4} ~ mm^{2}/sec)))
ggsave("lm8_RD/RD_vs_Age.png",
       RD_fig,
       device="png",
       width=20,
       height=7,
       units="cm",
       dpi=300)

# Save the fig7  source data
write.csv(dtinoddi_dat$both[, c("Age", "Sex", "FS6_masked_vol_cerebral_WM", "FS6_FA_cerebral_WM", "FS6_MD_cerebral_WM", "FS6_AD_cerebral_WM", "FS6_RD_cerebral_WM")],
          file="~/Dropbox/Ami-work/MRiShare/MRiShare_paper/MRi-Share eLIfe preps/Fig_source_data/Fig7.csv",
          row.names = FALSE)

#######################################################
###  9) NDI
#######################################################
dir.create("lm9_NDI")

# First save descriptive stats
writeDStxt(dtinoddi_dat, "cerebralWM_mean_NDI", "lm9_NDI/u26_NDI_descriptive_stats.txt")

# Check for outliers
NDI_mo = numOutliers(dtinoddi_dat$both$cerebralWM_mean_NDI)
NDI_eo = numOutliers(dtinoddi_dat$both$cerebralWM_mean_NDI, iqr_thr = 3.0)

# Create models and save the results
lm9_NDI_models <- getLMmodels(dtinoddi_dat, metric_name = "cerebralWM_mean_NDI")
saveAllLMmodels(lm9_NDI_models, dtinoddi_dat, "lm9_NDI", metric_name = "cerebralWM_mean_NDI")

# Plot
NDI_fig <- genSimpleAgeEffectPlots(lm9_NDI_models, "dtinoddi_dat",
                                   main_model_name = "lin_noVol",
                                   scaling_factor = 1,
                                   yvar="cerebralWM_mean_NDI",
                                   ylab="NDI")
ggsave("lm9_NDI/NDI_vs_Age.png",
       NDI_fig,
       device="png",
       width=20,
       height=7,
       units="cm",
       dpi=300)

#######################################################
###  10) ODI
#######################################################
dir.create("lm10_ODI")

# First save descriptive stats
writeDStxt(dtinoddi_dat, "cerebralWM_mean_ODI", "lm10_ODI/u26_ODI_descriptive_stats.txt")

# Check for outliers
ODI_mo = numOutliers(dtinoddi_dat$both$cerebralWM_mean_ODI)
ODI_eo = numOutliers(dtinoddi_dat$both$cerebralWM_mean_ODI, iqr_thr = 3.0)

# Create models and save the results
lm10_ODI_models <- getLMmodels(dtinoddi_dat, metric_name = "cerebralWM_mean_ODI")
saveAllLMmodels(lm10_ODI_models, dtinoddi_dat, "lm10_ODI", metric_name = "cerebralWM_mean_ODI")

# Plot
ODI_fig <- genSimpleAgeEffectPlots(lm10_ODI_models, "dtinoddi_dat",
                                    main_model_name = "lin_noVol",
                                    scaling_factor = 1,
                                    yvar="cerebralWM_mean_ODI",
                                    ylab="ODI")
ggsave("lm10_ODI/ODI_vs_Age.png",
       ODI_fig,
       device="png",
       width=20,
       height=7,
       units="cm",
       dpi=300)


#######################################################
###  11) IsoVF
#######################################################
dir.create("lm11_IsoVF")

# First save descriptive stats
writeDStxt(dtinoddi_dat, "cerebralWM_mean_ISOVF", "lm11_IsoVF/u26_IsoVF_descriptive_stats.txt")

# Check for outliers
IsoVF_mo = numOutliers(dtinoddi_dat$both$cerebralWM_mean_ISOVF)
IsoVF_eo = numOutliers(dtinoddi_dat$both$cerebralWM_mean_ISOVF, iqr_thr = 3.0)

# Create models and save the results
lm11_IsoVF_models <- getLMmodels(dtinoddi_dat, metric_name = "cerebralWM_mean_ISOVF")
saveAllLMmodels(lm11_IsoVF_models, dtinoddi_dat, "lm11_IsoVF", metric_name = "cerebralWM_mean_ISOVF")

# Plot
IsoVF_fig <- genSimpleAgeEffectPlots(lm11_IsoVF_models, "dtinoddi_dat",
                                    main_model_name = "lin_maskVol",
                                    scaling_factor = 1,
                                    yvar="cerebralWM_mean_ISOVF",
                                    ylab="IsoVF")
ggsave("lm11_IsoVF/IsoVF_vs_Age.png",
       IsoVF_fig,
       device="png",
       width=20,
       height=7,
       units="cm",
       dpi=300)
#####
# The End!
#####