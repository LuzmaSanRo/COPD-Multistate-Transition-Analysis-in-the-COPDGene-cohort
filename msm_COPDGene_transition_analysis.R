# Multi-State Model Analysis of the COPDGene Study Data to estimate transition probabilities between GOLD severity stages
# - This script performs cleaning, recoding, MSM fitting, and transition visualization.
# - For reproducibility and clarity, extensive comments and improved structure are provided.
# - NOTE: Data files referenced (COPDGene_P1P2P3_SM_NS_Long_Nov21.csv, COPDGene_VitalStatus_SM_NS_Oct22.csv) 
#     use the original file name as provided by the COPDGene Study. These dat files are not publicly available and 
#     require appropriate access and approvals from the COPDGene Study.

#### --- Libraries --- ####
library(readr)
library(dplyr)
library(msm)
library(reshape2)
library(ggplot2)
library (minqa)
library (tidyverse)

#### --- 1. Load COPDGene Longitudinal Data --- ####
# Path to your COPDGene longitudinal CSV file- 
# File name- original fila name as provided by the COPDGene Study
COPD <- read.csv("COPDGene_P1P2P3_SM_NS_Long_Nov21.csv")  

#### ---  1.1. Data Dictionary --- ####

# - sid : Subject ID
# - smoking_status: Smoking Status (3-state) 0=Never-smoked | 1=Former smoker | 2=Current smoker
# - finalgold_visit: GOLD stage, visit
###     0=GOLD 0 Control (FEV1 >= 80% & FEV1/FVC >= 0.7) | 1=GOLD 1 (FEV1 >= 80% & FEV1/FVC < 0.7) | 2=GOLD 2 (50% <= FEV1 < 80% & FEV1/FVC < 0.7) |
###     3=GOLD 3 (30% <= FEV1 < 50% & FEV1/FVC < 0.7) | 4=GOLD 4 (FEV1 < 30% & FEV1/FVC < 0.7) 
###    -1=PRISm (Preserved Ratio Imparied Spirometry) (FEV1/FVC >= 0.7 but FEV1 < 80%) | -2=Never-smoked Normal | -9=ineligible
# - years_from_baseline: Visit, years from baseline
# - Phase_study: COPDGene 5-year study phase (replaces visitnum)


#### --- 2. Data Cleaning --- ####

# Remove missing or 'Never' smoking status
df_COPD <- COPD %>%
  filter(!is.na(smoking_status) & smoking_status != "0")

# Remove missing GOLD stage
df_COPD <- df_COPD %>%
  filter(!is.na(finalgold_visit))


# Remove subjects with only one phase of data
df_COPD <- df_COPD %>%
  group_by(sid) %>%
  filter(n() > 1) %>%
  ungroup()

# Recode GOLD stage for analysis
# Mapping: GOLD 0->1, PRISm-1 >2, GOLD 1->3, GOLD 2->4, GOLD 3&4->5
df_COPD <- df_COPD %>%
  mutate(
    GOLD_new = case_when(
      finalgold_visit == 0 ~ 1,   # GOLD 0
      finalgold_visit == -1 ~ 2,  # PRISm
      finalgold_visit == 1 ~ 3,   # GOLD 1
      finalgold_visit == 2 ~ 4,   # GOLD 2
      finalgold_visit %in% c(3, 4) ~ 5 # GOLD 3 & 4
    )
  )

#### --- 3. Load and Clean Death Data --- ####
# File name- original fila name as provided by the COPDGene Study

df_death <- read_csv("COPDGene_VitalStatus_SM_NS_Oct22.csv")

#### ---  3.1. Data Dictionary --- ####
# - vital_status: Vital status , 1=deceased, 0=alive
# - days_followed:  Days followed: net of all contact dates

# Recode vital status: 0=Alive, 1=Death -> 1=Alive, 2=Death
df_death <- df_death %>%
  mutate(
    vital_status = recode(vital_status, `0` = 1, `1` = 2),
    Phase_Death = ifelse(vital_status == 2, 6, NA), # 6 = Death state
    newv = 1,
    year_dead = days_followed / 365.25 # Calculate year of death
  )

#### --- 4. Combine Clinical and Death Data --- ####
df_COPD_death <- bind_rows(df_COPD, df_death) %>%
  mutate(GOLD_death = GOLD_new)

df_COPD_death= df_COPD_death %>%
  mutate(GOLD_death = ifelse(is.na(GOLD_death) & Phase_Death == 6, 6, GOLD_death))

# Fill missing years_from_baseline with year_dead
df_COPD_death <- df_COPD_death %>%
  mutate(years_from_baseline = ifelse(is.na(years_from_baseline), year_dead, years_from_baseline))

#### --- 5. Prepare Data for MSM Analysis --- ####

# Select analysis variables and order by subject/phase
df_new <- df_COPD_death %>%
  select(sid, Phase_study, smoking_status, gender, finalgold_baseline, finalgold_visit,
         GOLD_new, GOLD_death, years_from_baseline) %>%
  arrange(sid, Phase_study)

df_new = df_new %>%
  group_by(sid) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  filter(!is.na(GOLD_death)) # Remove NA GOLD_death

statetable.msm(GOLD_death, sid, data=df_new)  #Checked matrix of number of observations from and to GOLD stages

#### --- 6. Descriptive Table --- ####

# GOLD stage by phase
GOLD_Stage <- table(df_new$GOLD_death, df_new$Phase_study, exclude = NULL)
addmargins(GOLD_Stage)

#### --- 7. Define Multi-State Transition Matrix --- ####

# States: 1=GOLD 0, 2=PRISm, 3=GOLD 1, 4=GOLD 2, 5=GOLD 3-4, 6=Death
Qm <- rbind(
  c(0, 1, 1, 0, 0, 1), # GOLD 0
  c(1, 0, 0, 1, 0, 1), # PRISm
  c(1, 0, 0, 1, 0, 1), # GOLD 1
  c(0, 1, 1, 0, 1, 1), # GOLD 2
  c(0, 0, 0, 1, 0, 1), # GOLD 3-4
  c(0, 0, 0, 0, 0, 1)  # Death (absorbing)
)


#### --- 8. MSM Initial Parameters --- ####

# Remove subjects appearing only once for MSM
df_new <- df_new[df_new$sid %in% df_new$sid[duplicated(df_new$sid)==1], ]
df_new$newphase <- df_new$Phase_study * 5

parinits_death <- log(
  array(
    t(crudeinits.msm(GOLD_death ~ years_from_baseline, sid, data = df_new, qmatrix = Qm)
    )
  )[which(array(t(Qm)) == 1)]
)

#### --- 9. Fit Multi-State Model (no covariates) --- ####

msm_copd <- msm(GOLD_death ~ years_from_baseline, data = df_new, qmatrix = Qm,
                    subject = sid, gen.inits = TRUE, opt.method = "bobyqa", fixedpar = FALSE)

# Output estimates and sojourn times
print(msm_copd$estimates)
print(cbind(exp(msm_copd$estimates), msm_copd$ci))
print(sojourn.msm(msm_copd, cl = 0.95))

#### --- 10. Transition Probability Matrices (1 & 5 years) --- ####

msm_copd_1y <- pmatrix.msm(msm_copd, t = 1, covariates = "mean", ci = "normal", cl = 0.95)
msm_copd_5y <- pmatrix.msm(msm_copd, t = 5, covariates = "mean", ci = "normal", cl = 0.95)

#### --- 11. Visualization of Transition Probabilities --- ####

# State labels for plot axes
gp_all_noDeath <- factor(c("GOLD 0", "PRISM", "GOLD 1", "GOLD 2", "GOLD 3-4", "Death"),
                         levels = c("GOLD 0", "PRISM", "GOLD 1", "GOLD 2", "GOLD 3-4", "Death"), ordered = TRUE)
gp_all_noDeathB <- factor(c("GOLD 3-4", "GOLD 2", "GOLD 1", "PRISM", "GOLD 0"),
                          levels = c("GOLD 3-4", "GOLD 2", "GOLD 1", "PRISM", "GOLD 0"), ordered = TRUE)

# --- 1-Year Transition Matrix Plot ---
pmsm_copd_1y <- msm_copd_1y$estimates[5:1, ] # reverse for y-axis
y_copd_msm_1y <- melt(pmsm_copd_1y)
y_copd_msm_1y$valueA <- round(y_copd_msm_1y$value * 100, 1)
base_size <- 18

fig_copd_msm_1y <- ggplot(y_copd_msm_1y, aes(y = Var1, x = Var2)) +
  geom_tile(aes(fill = valueA)) +
  scale_fill_gradient2(low = "white", mid = "orange", high = "red",
                       midpoint = 50.0, limits = c(0, 100),
                       guide = guide_colorbar(title = "Avg. 1-year\ntransition\nprobability (%)")) +
  geom_text(aes(label = sprintf("%.1f", valueA)), size = 5) +
  xlab("To State") + ylab("From State") +
  ggtitle("One-year transitions") +
  theme(
    legend.text = element_text(size = base_size),
    title = element_text(size = base_size + 2),
    axis.text.x = element_text(size = base_size, angle = 50, hjust = -0.3, colour = "grey20"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = base_size, colour = "grey20"),
    axis.title = element_text(size = base_size),
    axis.line = element_blank(), panel.grid.major = element_blank(),
    legend.title = element_text(size = base_size),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) +
  scale_x_discrete(expand = c(0, 0), labels = gp_all_noDeath, position = "top") +
  scale_y_discrete(expand = c(0, 0), labels = gp_all_noDeathB)

pdf("Transitions_COPD_1y.pdf")
print(fig_copd_msm_1y)
dev.off()

# --- 5-Year Transition Matrix Plot ---
pmsm_copd_5y <- msm_copd_5y$estimates[5:1, ]
y_copd_msm_5y <- melt(pmsm_copd_5y)
y_copd_msm_5y$valueA <- round(y_copd_msm_5y$value * 100, 1)

fig_copd_msm_5y <- ggplot(y_copd_msm_5y, aes(y = Var1, x = Var2)) +
  geom_tile(aes(fill = valueA)) +
  scale_fill_gradient2(low = "white", mid = "orange", high = "red",
                       midpoint = 50.0, limits = c(0, 100),
                       guide = guide_colorbar(title = "Avg. 5-year\ntransition\nprobability (%)")) +
  geom_text(aes(label = sprintf("%.1f", valueA)), size = 5) +
  xlab("To State") + ylab("From State") +
  ggtitle("Five-year transitions") +
  theme(
    legend.text = element_text(size = base_size),
    title = element_text(size = base_size + 2),
    axis.text.x = element_text(size = base_size, angle = 50, hjust = -0.3, colour = "grey20"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = base_size, colour = "grey20"),
    axis.title = element_text(size = base_size),
    axis.line = element_blank(), panel.grid.major = element_blank(),
    legend.title = element_text(size = base_size),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) +
  scale_x_discrete(expand = c(0, 0), labels = gp_all_noDeath, position = "top") +
  scale_y_discrete(expand = c(0, 0), labels = gp_all_noDeathB)

pdf("Transitions_COPD_5y.pdf")
print(fig_copd_msm_5y )
dev.off()

# --- Sojourn-Time by GOLD severtiy stage  Plot ---

df_soj = as.data.frame(sojourn.msm(msm_copd, cl = 0.95)) %>%
  rownames_to_column(var = "State")

print(df_soj)

# COPD- GOLD states labels
state_labels <- factor(c("GOLD 0", "PRISm", "GOLD 1", "GOLD 2", "GOLD 3-4"),
                       levels = c("GOLD 0", "PRISm", "GOLD 1", "GOLD 2", "GOLD 3-4"), ordered = TRUE)

# Plot
fig_copd_soj= ggplot(df_soj,  aes(y = estimates, x = State)) +
  geom_col(fill = "gray") +
  geom_errorbar(aes(ymin = L, ymax = U), width = 0.2) +
  geom_text(aes(label =sprintf("%.2f", estimates)), 
            vjust = 0.5, 
            hjust = 0.5,
            color = "black", 
            size = 4, 
            fontface = "bold",
            position = position_stack(vjust = 0.5)) +
  labs(
    title = "Sojourn Time by GOLD Category with 95% CI",
    x = "COPD Severity Stages",
    y = "Years"
  ) +
  theme (              # extra customization
    axis.line = element_line(color = "black", size = 0.3),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  ) +
  scale_x_discrete(expand = c(0, 0), labels = state_labels, position = "bottom") +
  scale_y_continuous(breaks = seq(0, 20, by = 2))

pdf("Sojourn_COPD.pdf")
print(fig_copd_soj)
dev.off()

#### --- Descriptive  Data --- ####
mean(df_COPD_death$age_baseline, na.rm=TRUE)
sd(df_COPD_death$age_baseline, na.rm=TRUE)

# 1. ##Descriptive Table for overall sample 
Overall = table(df_new$GOLD_death, df_new$Phase_study, exclude = NULL) # no tiene la informacion de las covariables
addmargins (Overall)
prop.table(Overall,2)
addmargins (round(100*prop.table(Overall,2), 2), 2)

###Estimate hazard for progression and regression  (1st Regression 2nd progression)

#- Hazard ratio of progression and regression--#
#A.1.) FROM GOLD 0 TO GOLD 1 OR GOLD 0 to PRISm(ref) 
qratio.msm(msm_copd,ind1=c(1,3),ind2=c(1,2))
#A.2) FROM PRISm To GOLD 2 and From PRISm to GOLD 0(ref)
qratio.msm(msm_copd,ind1=c(2,4),ind2=c(2,1))
#A.3.) FROM GOLD 1 TO GOLD 2 OR GOLD 1 to GOLD 0 (ref) 
qratio.msm(msm_copd,ind1=c(3,4),ind2=c(3,1))
#A.5.) FROM  GOLD 1 TO GOLD 2 VS PRISm TO GOLD 2 (ref)
qratio.msm(msm_copd,ind1=c(3,4),ind2=c(2,4)) 
#A.4) FROM GOLD 2 TO GOLD 3-4 OR  GOLD 2 to GOLD 1(ref) 
qratio.msm(msm_copd,ind1=c(4,5),ind2=c(4,3))
#A.4) FROM GOLD 2 TO GOLD 3-4 OR  GOLD 2 to PRISm (ref) 
qratio.msm(msm_copd,ind1=c(4,5),ind2=c(4,2))
#A.3.3) GOLD 2 to GOLD 1 or Gold 2 to PRISM (ref) 
qratio.msm(msm_copd,ind1=c(4,3),ind2=c(4,2))

#### --- End of Script --- ####
# This code is intended for sharing and reproducibility in research repositories.
# Please cite the COPDGene Study and ensure data use approvals for the referenced datasets.


