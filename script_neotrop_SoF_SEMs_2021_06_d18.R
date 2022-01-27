#'---
#' title: Estimation of the scale of effect for diversity measures in the neotropics and implementation of SEMs
#' Rproject: C:\Users\lucas\OneDrive\Documentos\01_MASTER\01_raw\00_neotropics\01_SoF_SEMs\SoF_SEMs_neotrop.Rproj
#' author: Lucas Augusto Pereira
#' date: 2021-06-d18
#' description: The first part of this script estimate the scale of effect (SoF) of three landscape metrics: % of
#' foret cover, number of patches, and edge density for 169 primate communities in the neotropical region, using
#' the multifit function (Huais, 2018). The landscape metrics was extracted from 11 different buffers (1000m to 5000m
#  with 400 m interval from each distance), so the SoF for all three metrics was assessed based on these 11 buffers.
#' The second part of this script evaluate the potential spatial autocorrelation presented in the data, by testing five
#' different spatial correlation structures (exponential, gaussian, linear, rational quadratic, and spheric). The third
#' and last part of this script we implement the landscape metrics in their respectve scales of effect in structural
#' equation models to investigate the direct and indirect effects of habitat loss and fragmentation on the species
#' richness and functional diversity of primate communities in the neotropics. 
#'---

#loading packages
library(tidyverse)
library(fitdistrplus)
library(lme4)
library(nlme)
library(MuMIn)
library(DescTools)
library(ggpubr)
library(piecewiseSEM)

# import data --------------------------------------------------------------------------------------


# importing the .csv file with the diversity measures for each site, the geographical coordinates
# and the landscape metrics

db.div.land <- readr::read_delim(here::here("FD_SR_metrics_06_d18_2021.csv"), delim = ";")
db.div.land


# finding the distribution of the data --------------------------------------------------------------

# Functional diversity (gamma distribution)
descdist(db.div.land$Functional_diversity, discrete = FALSE)
FD_fit <- fitdist(db.div.land$Functional_diversity, "gamma")
plot(FD_fit)

# Species richness (poisson distribution)
descdist(db.div.land$Species_richness, discrete = TRUE)
SR_fit <- fitdist(db.div.land$Species_richness, "nbinom")
plot(SR_fit)

# number of patches (poisson distribution)
descdist(db.div.land$np_4600, discrete = TRUE)
NP_fit <- fitdist(db.div.land$np_4600, "pois")
plot(NP_fit)

# edge density (gamma distribution)
descdist(db.div.land$ed_5000, discrete = FALSE)
ED_fit <- fitdist(db.div.land$ed_5000, method = "mme", "gamma")
plot(ED_fit)



# FIRST PART: Scale of effect for species richness (SR) and functional diversity (FD) -------------------------------------



# estimating the scale of effect of % forest cover on SR
# dimension 508 x 401
pland.SoF.SR <- multifit(mod = "lme", multief = colnames(db.div.land)[7:17], 
                         formula = SR ~ multief, data = db.div.land, 
                         args = c("random = ~ 1 | region"), criterion = "AIC", 
                         signif = TRUE, alpha = 0.05, ylab = "% Forest cover", 
                         labels = c(seq(1000, 5000, 400)))

pland.SoF.SR.DB <- as.data.frame(pland.SoF.SR$summary)
write.csv(pland.SoF.SR.DB, "00_pland_SoF_SR_06_d19_2021.csv")

# and FD
pland.SoF.FD <- multifit(mod = "lme", multief = colnames(db.div.land)[7:17], 
                         formula = FD ~ multief, data = db.div.land, 
                         args = c("random = ~ 1 | region"), criterion = "AIC", 
                         signif = TRUE, alpha = 0.05, ylab = " ", 
                         labels = c(seq(1000, 5000, 400)))

pland.SoF.FD.DB <- as.data.frame(pland.SoF.FD$summary)
write.csv(pland.SoF.FD.DB, "01_pland_SoF_FD_06_d19_2021.csv")



# estimating the scale of effect of number of patches on SR
np.SoF.SR <- multifit(mod = "lme", multief = colnames(db.div.land)[18:28], 
                      formula = SR ~ multief, data = db.div.land, 
                      args = c("random = ~ 1 | region"), criterion = "AIC", 
                      signif = TRUE, alpha = 0.05, ylab = "Number of patches", 
                      labels = c(seq(1000, 5000, 400)))

np.SoF.SR.DB <- as.data.frame(np.SoF.SR$summary)
write.csv(np.SoF.SR.DB, "02_np_SoF_SR_06_d19_2021.csv")

# and FD
np.SoF.FD <- multifit(mod = "lme", multief = colnames(db.div.land)[18:28], 
                      formula = FD ~ multief, data = db.div.land, 
                      args = c("random = ~ 1 | region"), criterion = "AIC", 
                      signif = TRUE, alpha = 0.05, ylab = " ", 
                      labels = c(seq(1000, 5000, 400)))

np.SoF.FD.DB <- as.data.frame(np.SoF.FD$summary)
write.csv(np.SoF.FD.DB, "03_np_SoF_FD_06_d19_2021.csv")



# estimating the scale of effect of edge density on SR
ed.SoF.SR <- multifit(mod = "lme", multief = colnames(db.div.land)[29:39], 
                      formula = SR ~ multief, data = db.div.land, 
                      args = c("random = ~ 1 | region"), criterion = "AIC", 
                      signif = TRUE, alpha = 0.05, xlab = "Landscape size (m)", 
                      ylab = "Edge density", labels = c(seq(1000, 5000, 400)))

ed.SoF.SR.DB <- as.data.frame(ed.SoF.SR$summary)
write.csv(ed.SoF.SR.DB, "04_ed_SoF_SR_06_d19_2021.csv")

# and FD
ed.SoF.FD <- multifit(mod = "lme", multief = colnames(db.div.land)[29:39], 
                      formula = FD ~ multief, data = db.div.land, 
                      args = c("random = ~ 1 | region"), criterion = "AIC", 
                      signif = TRUE, alpha = 0.05, xlab = "Landscape size (m)", 
                      ylab = " ", labels = c(seq(1000, 5000, 400)))

ed.SoF.FD.DB <- as.data.frame(ed.SoF.FD$summary)
write.csv(ed.SoF.FD.DB, "05_ed_SoF_FD_06_d19_2021.csv")

# END OF THE FIRST PART ------------------------------------------------------------------------------------ 

# SECOND PART: Evaluating the potential spatial autocorrelation presented in the spatially distributed data ------------

# selecting variables of interest

db.div.SoF <- db.div.land %>% 
  dplyr::select(ID, FD, SR, region, long_x, lat_y, pland_5000, np_1400, np_2200, ed_5000)
db.div.SoF


# testing spatial autocorrelation for equation 1

noCorr_m1 <- lme(np_2200 ~ pland_5000, random = ~ 1 | region, data = db.div.SoF)
summary(noCorr_m1)

corExp_m1 <- lme(np_2200 ~ pland_5000, random = ~ 1 | region, 
                 correlation = corExp(form = ~ long_x + lat_y), data = db.div.SoF)
summary(corExp_m1)

corGaus_m1 <- lme(np_2200 ~ pland_5000, random = ~ 1 | region, 
                  correlation = corGaus(form = ~ long_x + lat_y), data = db.div.SoF)
summary(corGaus_m1)

corLin_m1 <- lme(np_2200 ~ pland_5000, random = ~ 1 | region, 
                 correlation = corLin(form = ~ long_x + lat_y), data = db.div.SoF)
summary(corLin_m1)

corRatio_m1 <- lme(np_2200 ~ pland_5000, random = ~ 1 | region, 
                   correlation = corRatio(form = ~ long_x + lat_y), data = db.div.SoF)
summary(corRatio_m1)

corSpher_m1 <- lme(np_2200 ~ pland_5000, random = ~ 1 | region, 
                   correlation = corSpher(form = ~ long_x + lat_y), data = db.div.SoF)
summary(corSpher_m1)

# best model (lowest AIC = 2032.286) is corExp_m1



# testing spatial autocorrelation for equation 2

noCorr_m2 <- lme(ed_5000 ~ np_2200 + pland_5000, random = ~ 1 | region, data = db.div.SoF)
summary(noCorr_m2)

corExp_m2 <- lme(ed_5000 ~ np_2200 + pland_5000, random = ~ 1 | region, 
                 correlation = corExp(form = ~ long_x + lat_y), data = db.div.SoF)
summary(corExp_m2)

corGaus_m2 <- lme(ed_5000 ~ np_2200 + pland_5000, random = ~ 1 | region, 
                  correlation = corGaus(form = ~ long_x + lat_y), data = db.div.SoF)
summary(corGaus_m2)

corLin_m2 <- lme(ed_5000 ~ np_2200 + pland_5000, random = ~ 1 | region, 
                 correlation = corLin(form = ~ long_x + lat_y), 
                 control = lmeControl(opt = "optim"), data = db.div.SoF)
summary(corLin_m2)

corRatio_m2 <- lme(ed_5000 ~ np_2200 + pland_5000, random = ~ 1 | region, 
                   correlation = corRatio(form = ~ long_x + lat_y), data = db.div.SoF)
summary(corRatio_m2)

corSpher_m2 <- lme(ed_5000 ~ np_2200 + pland_5000, random = ~ 1 | region, 
                   correlation = corSpher(form = ~ long_x + lat_y), data = db.div.SoF)
summary(corSpher_m2)

# best model (lowest AIC = 1447.799) is corRatio_m2



# testing spatial autocorrelation for equation 3

noCorr_m3 <- lme(SR ~ pland_5000 + np_2200 + ed_5000, random = ~ 1 | region, data = db.div.SoF)
summary(noCorr_m3)

corExp_m3 <- lme(SR ~ pland_5000 + np_2200 + ed_5000, random = ~ 1 | region, 
                 correlation = corExp(form = ~ long_x + lat_y), data = db.div.SoF)
summary(corExp_m3)

corGaus_m3 <- lme(SR ~ pland_5000 + np_2200 + ed_5000, random = ~ 1 | region, 
                  correlation = corGaus(form = ~ long_x + lat_y), data = db.div.SoF)
summary(corGaus_m3)

corLin_m3 <- lme(SR ~ pland_5000 + np_2200 + ed_5000, random = ~ 1 | region, 
                 correlation = corLin(form = ~ long_x + lat_y), data = db.div.SoF)
summary(corLin_m3)

corRatio_m3 <- lme(SR ~ pland_5000 + np_2200 + ed_5000, random = ~ 1 | region, 
                   correlation = corRatio(form = ~ long_x + lat_y), data = db.div.SoF)
summary(corRatio_m3)

corSpher_m3 <- lme(SR ~ pland_5000 + np_2200 + ed_5000, random = ~ 1 | region, 
                   correlation = corSpher(form = ~ long_x + lat_y), data = db.div.SoF)
summary(corSpher_m3)

# best model (lowest AIC = 737.5832) is corExp_m3



# testing spatial autocorrelation for equation 4

noCorr_m4 <- lme(FD ~ pland_5000 + np_2200 + ed_5000 + SR, random = ~ 1 | region, data = db.div.SoF)
summary(noCorr_m4)

corExp_m4 <- lme(FD ~ pland_5000 + np_2200 + ed_5000 + SR, random = ~ 1 | region, 
                 correlation = corExp(form = ~ long_x + lat_y), data = db.div.SoF)
summary(corExp_m4)

corGaus_m4 <- lme(FD ~ pland_5000 + np_2200 + ed_5000 + SR, random = ~ 1 | region, 
                  correlation = corGaus(form = ~ long_x + lat_y), data = db.div.SoF)
summary(corGaus_m4)

corLin_m4 <- lme(FD ~ pland_5000 + np_2200 + ed_5000 + SR, random = ~ 1 | region, 
                 correlation = corLin(form = ~ long_x + lat_y), 
                 control = lmeControl(opt = "optim"), data = db.div.SoF)
summary(corLin_m4)

corRatio_m4 <- lme(FD ~ pland_5000 + np_2200 + ed_5000 + SR, random = ~ 1 | region, 
                   correlation = corRatio(form = ~ long_x + lat_y), data = db.div.SoF)
summary(corRatio_m4)

corSpher_m4 <- lme(FD ~ pland_5000 + np_2200 + ed_5000 + SR, random = ~ 1 | region, 
                   correlation = corSpher(form = ~ long_x + lat_y), data = db.div.SoF)
summary(corSpher_m4)

# best model (lowest AIC = -278.5873) is corExp_m4

# END OF THE SECOND PART ------------------------------------------------------------------------------------

# THIRD PART: Assessing the causal relationships between landscape metrics and species richness/functional diversity ------------

#We needed to remove a causal connection to test  the model since the direct+indirect model 
#had no independence relationships conceptually (saturated model, no basis set to test).
#We computed the C statistics for variations of direct+indirect model, removing one causal relation at a time. 
#Then, we chose to remove the causal relationship between forest cover and functional diversity (path 7), 
# because it generated a model with the smallest C statistics and the highest p-value (C=0.103 and p=0.95)

# saturated model

sat.sem <- psem(

	lme(np_2200 ~ pland_5000, random = ~ 1 | region, 
	    correlation = corRatio(form = ~ long_x + lat_y), data = db.div.SoF),

	lme(ed_5000 ~ np_2200 + pland_5000, random = ~ 1 | region, 
	    correlation = corRatio(form = ~ long_x + lat_y), data = db.div.SoF),

	lme(SR ~ pland_5000 + np_2200 + ed_5000, random = ~ 1 | region, 
	    correlation = corExp(form = ~ long_x + lat_y), data = db.div.SoF),

	lme(FD ~ pland_5000 + np_2200 + ed_5000 + SR, random = ~ 1 | region, 
	    correlation = corExp(form = ~ long_x + lat_y), data = db.div.SoF)

	)

summary(sat.sem)

# direct and indirect model

sem.neotrop.m1 <- psem(

	lme(np_2200 ~ pland_5000, random = ~ 1 | region, 
	    correlation = corRatio(form = ~ long_x + lat_y), data = db.div.SoF),

	lme(ed_5000 ~ np_2200 + pland_5000, random = ~ 1 | region, 
	    correlation = corRatio(form = ~ long_x + lat_y), data = db.div.SoF),

	lme(SR ~ pland_5000 + np_2200 + ed_5000, random = ~ 1 | region, 
	    correlation = corExp(form = ~ long_x + lat_y), data = db.div.SoF),

	lme(FD ~ np_2200 + ed_5000 + SR, random = ~ 1 | region, 
	    correlation = corExp(form = ~ long_x + lat_y), data = db.div.SoF)

	)

summary(sem.neotrop.m1)
sem.m1.coefs <- as.data.frame(coefs(sem.neotrop.m1))
write.csv(sem.m1.coefs, "m1_SEMs_coefs_06_d20_2021.csv")

# indirect model

sem.neotrop.m2 <- psem(

	lme(np_2200 ~ pland_5000, random = ~ 1 | region, 
	    correlation = corRatio(form = ~ long_x + lat_y), data = db.div.SoF),

	lme(ed_5000 ~ np_2200 + pland_5000, random = ~ 1 | region, 
	    correlation = corRatio(form = ~ long_x + lat_y), data = db.div.SoF),

	lme(SR ~ pland_5000 + np_2200 + ed_5000, random = ~ 1 | region, 
	    correlation = corExp(form = ~ long_x + lat_y), data = db.div.SoF),

	lme(FD ~ SR, random = ~ 1 | region, 
	    correlation = corExp(form = ~ long_x + lat_y), data = db.div.SoF)

	)

summary(sem.neotrop.m2)
sem.m2.coefs <- as.data.frame(coefs(sem.neotrop.m2))
write.csv(sem.m2.coefs, "m2_SEMs_coefs_06_d20_2021.csv")

# direct model

sem.neotrop.m3 <- psem(

	lme(np_2200 ~ pland_5000, random = ~ 1 | region, 
	    correlation = corRatio(form = ~ long_x + lat_y), data = db.div.SoF),

	lme(ed_5000 ~ np_2200 + pland_5000, random = ~ 1 | region, 
	    correlation = corRatio(form = ~ long_x + lat_y), data = db.div.SoF),

	lme(SR ~ pland_5000 + np_2200 + ed_5000, random = ~ 1 | region, 
	    correlation = corExp(form = ~ long_x + lat_y), data = db.div.SoF),

	lme(FD ~ pland_5000 + np_2200 + ed_5000, random = ~ 1 | region, 
	    correlation = corExp(form = ~ long_x + lat_y), data = db.div.SoF)

	)

summary(sem.neotrop.m3)
sem.m3.coefs <- as.data.frame(coefs(sem.neotrop.m3))
write.csv(sem.m3.coefs, "m3_SEMs_coefs_06_d20_2021.csv")

# END OF THE SCRIPT ------------------------------------------------------------------------------------