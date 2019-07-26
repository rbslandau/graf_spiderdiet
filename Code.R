# R Script Graf_et_al_2019_R: R-Code Nadin Graf, Karina P. Battes, Mirela Cimpean, Martin H.
# Entling, Katharina Frisch, Moritz Link, Andreas Scharmueller, Verena C. Schreiner, Eduard
# Szoecs, Jochen P. Zubrod, Ralf B. Schaefer
#----------------------------------------------------------------------------------------
# Relationship between agricultural pesticides and the diet of riparian spiders in a field study
#----------------------------------------------------------------------------------------
# submitted to 'Environmental Sciences Europe'

# The code has been written by: 
# Nadin Graf 
# University of Koblenz-Landau
# Fortstrasse 7
# 76829 Landau
# GERMANY
# Email: graf-nadin@uni-landau.de, schaefer-ralf@uni-landau.de

# Revised by Moritz Link & Ralf B. Schaefer

# The code has been written to reproduce the results of the related publication

####################### Structure #

# Structure of the code:   
# 1.  Diet -------------------------------- 
# 1.1 Diet Pardosa  -----------------------
# 1.2 Diet Tetragnatha -------------------- 
# 2.  Biomass ----------------------------- 
# 3.  SPCA -------------------------------- 
# 3.1 proxies for agricultural land use --- 
# 3.2 Create SPCA -------------------------
# 4.  Statistical model building ---------- 
# 4.1 Biomass ----------------------------- 
# 4.2 Tetragnatha ------------------------- 
# 4.2.1 Figure Tetragnatha ----------------
# 4.3 Pardosa -----------------------------
# 4.3.1 Figure Pardosa -------------------- 
# 4.3.2 Model for abundance of Pardosa ----
#------------------------------------------

# Put your working directory
path <- "~/Gitprojects/Publications/graf_spiderdiet"
setwd(path)

# To load the data, save the required worksheets of the provided Excel spreadsheet as .csv files

####### 1. Diet ##### 1.1 Diet Pardosa #########
library(MixSIAR)

# load consumer isotope signal for P. amentata
mix <- load_mix_data(filename = "consumer par.csv", iso_names = c("d13C", "d15N"), factors = c("site"), 
    fac_random = c(TRUE), fac_nested = c(FALSE), cont_effects = NULL)

# load isotope signals prey of P. amentata
source <- load_source_data(filename = "diet.csv", source_factors = "site", conc_dep = FALSE, data_type = "raw", mix)

# load TEF for prey
discr <- load_discr_data(filename = "TEF.csv", mix)

plot_data(filename = "isospace_plot", plot_save_pdf = FALSE, plot_save_png = FALSE, mix, source, discr)

model_filename <- "MixSIAR_model.txt"  
# Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# run = extreme: Chain length 3,000,000, Burn-in: 1,500,00, Thin: 500
jags.1 <- run_model(run = "extreme", mix, source, discr, model_filename, alpha.prior = 1, resid_err, 
    process_err)

output_options <- list(summary_save = TRUE, summary_name = "summary_statistics", sup_post = FALSE, 
    plot_post_save_pdf = FALSE, plot_post_name = "posterior_density", sup_pairs = FALSE, plot_pairs_save_pdf = FALSE, 
    plot_post_name = "pairs_plot", sup_xy = FALSE, plot_xy_save_pdf = FALSE, plot_xy_name = "xy_plot", 
    gelman = TRUE, heidel = FALSE, geweke = TRUE, diag_save = TRUE, diag_name = "diagnostics", indiv_effect = FALSE, 
    plot_post_save_png = FALSE, plot_pairs_save_png = FALSE, plot_xy_save_png = FALSE)

output_JAGS(jags.1, mix, source, output_options)

####### 1.2 Diet Tetragnatha #########

# load isotope signals for cosumer T. montana
mix <- load_mix_data(filename = "consumer tet.csv", iso_names = c("d13C", "d15N"), factors = c("site"), 
    fac_random = c(TRUE), fac_nested = c(FALSE), cont_effects = NULL)

# load isotope signals for prey of T. montana
source <- load_source_data(filename = "diet tet.csv", source_factors = "site", conc_dep = FALSE, data_type = "raw", 
    mix)

# load TEF for prey
discr <- load_discr_data(filename = "TEF.csv", mix)

plot_data(filename = "isospace_plot", plot_save_pdf = FALSE, plot_save_png = FALSE, mix, source, discr)

model_filename <- "MixSIAR_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

jags.1 <- run_model(run = "extreme", mix, source, discr, model_filename, alpha.prior = 1, resid_err, 
    process_err)

output_options <- list(summary_save = TRUE, summary_name = "summary_statistics", sup_post = FALSE, 
    plot_post_save_pdf = FALSE, plot_post_name = "posterior_density", sup_pairs = FALSE, plot_pairs_save_pdf = FALSE, 
    plot_post_name = "pairs_plot", sup_xy = FALSE, plot_xy_save_pdf = FALSE, plot_xy_name = "xy_plot", 
    gelman = TRUE, heidel = FALSE, geweke = TRUE, diag_save = TRUE, diag_name = "diagnostics", indiv_effect = FALSE, 
    plot_post_save_png = FALSE, plot_pairs_save_png = FALSE, plot_xy_save_png = FALSE)

output_JAGS(jags.1, mix, source, output_options)

####### 2. Biomass ########

setwd(path)

library(dplyr)
library(tidyr)
library(plyr)

trap_days <- read.table("emergence_trap_days.csv", header = TRUE, sep = ",", na.strings = c("NR", 
    "NA"))
aquatic_abundance <- read.table("aquatic_abundance.csv", header = TRUE, sep = ",", na.strings = c("NR", 
    "NA"))
aquatic_abundance$number <- as.numeric(aquatic_abundance$number)
biomass_estimation <- read.table("biomass_estimation.csv", header = TRUE, sep = ",", na.strings = c("NR", 
    "NA"))

# choose abundance of aquatic prey
abundance_ephemeroptera <- filter(aquatic_abundance, order == "Ephemeroptera")
abundance_trichoptera <- filter(aquatic_abundance, order == "Trichoptera")
abundance_chironomidae <- filter(aquatic_abundance, family == "Chironomidae")
abundance_simuliidae <- filter(aquatic_abundance, family == "Simuliidae")
abundance_empididae <- filter(aquatic_abundance, family == "Empididae")

# combine data sets
abundance <- bind_rows(abundance_ephemeroptera, abundance_trichoptera)
abundance <- bind_rows(abundance, abundance_chironomidae)
abundance <- bind_rows(abundance, abundance_simuliidae)
abundance <- bind_rows(abundance, abundance_empididae)

# column to join tables
abundance <- unite(abundance, "combi", "order", "suborder", "family", "genus", "species", remove = FALSE)
biomass_estimation <- unite(biomass_estimation, "combi", "order", "suborder", "family", "genus", "species", 
    remove = FALSE)

abundance <- left_join(abundance, trap_days, by = "Site")
abundance1 <- ddply(abundance, c("Site", "order"), summarise, number = sum(number))

write.csv(abundance1, file = "abundance.csv", row.names = TRUE)

abundance_biomass <- left_join(abundance, biomass_estimation[, c(1, 7:9)], by = "combi")

# add dry mass dry mass = a * Length^b
abundance_biomass$mass <- with(abundance_biomass, a * length^b)
abundance_biomass$biomass <- with(abundance_biomass, mass * number)
# biomass per day
abundance_biomass$biomass_day <- with(abundance_biomass, biomass/sum_biomass_trap_day)
# sum up biomass per site
abundance_biomass <- ddply(abundance_biomass, c("Site"), summarise, sum_biomass_day = sum(biomass_day))

abundance_biomass1 <- ddply(abundance_biomass, c("Site", "order"), summarise, sum_biomass_day = sum(biomass_day))
write.csv(abundance_biomass1, file = "abundance_biomass.csv", row.names = TRUE)

####### 3. SPCA ########
setwd(path)

# check env. vars and create spca axes

lu_data <- read.csv("environmental_data.csv", header = TRUE, sep = ",")

# Aim of the SPCA is to create a sparse principal component that reflects the intensity of
# agricultural landuse at each sampling site and puts it into relation to the other sampling sites
# Decide, which variables are a proxy of agricultural land use.

summary(lu_data)
head(lu_data)

####### 3.1 proxies for agricultural land use #### select variables that are suitable proxies for
####### agricultural land use and relevant for spiders and their aquatic prey

proxy_vars <- c("cl", "ec_ms", "ph", "no3", "po4", "so4", "o2", "shore_cover_shrubs", "shore_cover_meadow", 
    "fid_factor", "land_use_meadows", "land_use_agri", "totalarea", "ext_agri", "int_agri", "max_width", 
    "min_width", "riffles", "shading")

## only keep the proxy variables
lu_dat <- lu_data[ , proxy_vars]
str(lu_dat)

par(mfrow = c(3, 3))
for (nam in c(colnames(lu_dat))) {
    y <- lu_dat[, nam]
    plot(density(y), main = "", xlab = nam)
    plot(density(sqrt(y)), main = "", xlab = nam)
    plot(density(log10(y + 1)), main = "", xlab = nam)
}
library(pcaPP)
library(vegan)

####### 3.2 create SPCA ######

# set 2 as k.max, 2 axis
k.max <- 2

oTPO_D1 <- opt.TPO(scale(lu_dat), k.max = k.max, method = "sd")
summary(oTPO_D1$pc)

oTPO_D1$pc$load

par(mfrow = c(1, k.max))
for (i in 1:k.max) plot(oTPO_D1, k = i)
# L0 gives number of variables with zero loadings on PC lambda opt is the estimated optimal
# penalty ECV1 is the empirical cumulated variance of the PCs

## Tradeoff Curves: Explained Variance vs. lambda
par(mfrow = c(1, k.max))
for (i in 1:k.max) plot(oTPO_D1, k = i, f.x = "lambda")
# lambda represents penalty term

## Extract the spc scores
spc_D1 <- sPCAgrid(scale(lu_dat), k = k.max, lambda = oTPO_D1$pc.noord$lambda, method = "sd")
load_spcD1 <- scores(spc_D1, display = "species", scaling = 0)
pca_axes_D1 <- scores(spc_D1, display = "sites", scaling = 0)

pca_axes_D1$site <- lu_dat$site

write.csv(pca_axes_D1, file = "pca_axes_2R.csv")
# First principal component

pc1_D1 <- pca_axes_D1[order(pca_axes_D1$Comp.1), c(3, 1)]
pc1_D1

# Second principal component
pc2_D1 <- pca_axes_D1[order(pca_axes_D1$Comp.2), c(3, 2)]
pc2_D1


####### 4. Statistical model building ######

setwd(path)

dat <- read.csv("statistic.csv", header = TRUE, sep = ",")
library(qpcR)

####### 4.1 Biomass #######
mas <- dat[-c(2, 7, 8, 9, 14), ]

plot(density(mas$biomass))
qqnorm(mas$biomass)
qqline(mas$biomass, col = 2)

mas1 <- lm(biomass ~ max_sumTU_ms_corr + Comp.1 + Comp.2, mas)
summary(mas1)  #-> exclude comp2
drop1(mas1)

mas2 <- lm(biomass ~ max_sumTU_ms_corr + Comp.1, mas)
summary(mas2)  #-> exclude tox
drop1(mas2)

mas3 <- lm(biomass ~ Comp.1, mas)
summary(mas3)
plot(mas$biomass ~ mas$max_sumTU_ms_corr)

AIC(mas1, mas2, mas3)
AICc(mas1)
AICc(mas2)
AICc(mas3)

####### 4.2 Tetragnatha #####

tet <- dat[-c(2, 7, 8, 9, 14, 17), ]  # exclude sites with missing values

plot(density(tet$tet_aqu_diet))


# check if explanatory variables correlate
cor.test(tet$Comp.1, tet$biomass)
cor.test(tet$Comp.1, tet$max_sumTU_ms_corr)
cor.test(tet$biomass, tet$max_sumTU_ms_corr)
cor.test(tet$Comp.2, tet$max_sumTU_ms_corr)
cor.test(tet$Comp.2, tet$biomass)
cor.test(tet$Comp.2, tet$Comp.1)

# stepwise selection

library(betareg)
tet1 <- glm(tet_aqu_diet ~ max_sumTU_ms_corr + Comp.1 + Comp.2 + biomass, family = "binomial", tet)
summary(tet1)
drop1(tet1)  #-> exclude Comp.2
plot(tet1)

tet2 <- glm(tet_aqu_diet ~ max_sumTU_ms_corr + Comp.1 + biomass, family = "binomial", tet)
summary(tet2)
drop1(tet2)  #-> exclude max_sumTU_ms_corr

tet3 <- glm(tet_aqu_diet ~ Comp.1 + biomass, family = "binomial", tet)
summary(tet3)
plot(tet$tet_aqu_diet ~ tet$Comp.1)
drop1(tet3)  #-> final model
AIC(tet1, tet2, tet3)

# calculate explained deviance
library(modEvA)
Dsquared(tet3)

AICc(tet1)
AICc(tet2)
AICc(tet3)

library(effects)
library(trelliscopejs)
library(lattice)
# plot(predictorEffects(yourmodel))
plot(predictorEffects(tet3))

# par(mfrow=c(3,1), mar=c(4.5, 8, 1,1))

####### 4.2.1 Figure Tetragnatha ####
xlabpar <- trellis.par.get("par.xlab.text")
xlabpar$cex <- 1.7
trellis.par.set("par.xlab.text", xlabpar)
ylabpar <- trellis.par.get("par.ylab.text")
ylabpar$cex <- 1.7
trellis.par.set("par.ylab.text", ylabpar)

svg("tet_mass.svg", width = 5.7, height = 5.4)
plot(predictorEffects(tet3, "biomass"), lines = list(col = {
    "black"
}), main = paste(""), id = TRUE, points = list(cex = 4), axes = list(y = list(lim = (c(0.38, 0.62)), 
    ticks = list(at = c(0.4, 0.45, 0.5, 0.55, 0.6)), lab = {
        "Proportion of consumed aquatic diet"
    }, cex = 1.2), x = list(biomass = list(lim = (c(40, 610)), ticks = list(at = c(0, 100, 200, 300, 
    400, 500, 600)), lab = "Aquatic prey dry biomass [mg]"), rug = TRUE, cex = 1.2)))
dev.off()

range(tet$Comp.1)
svg("tet_spca1.svg", width = 5.7, height = 5.4)
plot(predictorEffects(tet3, "Comp.1"), lines = list(col = {
    "black"
}), main = paste(""), id = TRUE, points = list(cex = 4), axes = list(y = list(lim = (c(0.38, 0.62)), 
    ticks = list(at = c(0.4, 0.45, 0.5, 0.55, 0.6)), lab = {
        "Proportion of consumed aquatic diet"
    }, cex = 1.2), x = list(Comp.1 = list(lim = (c(-3.5, 5.5)), ticks = list(at = c(-3, -2, -1, 0, 
    1, 2, 3, 4, 5)), lab = "First SPCA axis"), rug = TRUE, cex = 1.2)))
dev.off()

####### 4.3 Pardosa #######

par <- dat[-c(1, 2, 7, 8, 9, 10, 14, 15, 19), ]  # exclude sites with missing values

plot(density(par$par_aqu_diet))
qqnorm(par$par_aqu_diet)
qqline(par$par_aqu_diet, col = 2)

range(par$par_aqu_diet, na.rm = TRUE)

# check if explanatory variables correlate
cor.test(par$Comp.1, par$Comp.2)
cor.test(par$max_sumTU_ms_corr, par$Comp.2)
cor.test(par$max_sumTU_ms_corr, par$Comp.1)
cor.test(par$biomass, par$max_sumTU_ms_corr)
cor.test(par$biomass, par$Comp.2)
cor.test(par$biomass, par$Comp.1)

par1 <- lm(par_aqu_diet ~ max_sumTU_ms_corr + Comp.1 + Comp.2 + biomass, par)
summary(par1)
drop1(par1)  #-> exclude biomass

par2 <- lm(par_aqu_diet ~ max_sumTU_ms_corr + Comp.1 + Comp.2, par)
summary(par2)
drop1(par2)  #-> exclude Comp.1

par3 <- lm(par_aqu_diet ~ max_sumTU_ms_corr + Comp.2, par)
summary(par3)
drop1(par3)  #-> final model

AIC(par1, par2, par3)
AICc(par1)
AICc(par2)
AICc(par3)

####### 4.3.1 Figure Pardosa ####
range(par$max_sumTU_ms_corr)
svg("par_tox.svg", width = 5.7, height = 5.4)
plot(predictorEffects(par3, "max_sumTU_ms_corr"), lines = list(col = {
    "black"
}), main = paste(""), id = TRUE, points = list(cex = 4), axes = list(y = list(lim = (c(0.3, 0.55)), 
    ticks = list(at = c(0.3, 0.35, 0.4, 0.45, 0.5, 0.55)), lab = {
        "Proportion of consumed aquatic diet"
    }, cex = 1.2), x = list(max_sumTU_ms = list(lim = (c(-1.65, 0.3)), ticks = list(at = c(-1.4, -1.2, 
    -1, -0.8, -0.6, -0.4, -0.2, 0)), lab = "Toxicity [max sumTU]"), rug = TRUE, cex = 1.2)))
dev.off()

range(par$Comp.2)
svg("par_spca2.svg", width = 5.7, height = 5.4)
plot(predictorEffects(par3, "Comp.2"), lines = list(col = {
    "black"
}), main = paste(""), id = TRUE, points = list(cex = 4), axes = list(y = list(lim = (c(0.3, 0.55)), 
    ticks = list(at = c(0.3, 0.35, 0.4, 0.45, 0.5, 0.55)), lab = {
        "Proportion of consumed aquatic diet"
    }, cex = 1.2), x = list(Comp.1 = list(lim = (c(-2.5, 3.5)), ticks = list(at = c(-3, -2, -1, 0, 
    1, 2, 3, 4, 5)), lab = "Second SPCA axis"), rug = TRUE, cex = 1.2)))
dev.off()

####### 4.3.2 model for abundance of Pardosa #####
numpar <- glm(num_par ~ max_sumTU_ms_corr, family = "quasipoisson", dat)
summary(numpar)
plot(numpar)
1 - (140.57/273.21)  #-> D2 = 0.485
Dsquared(numpar)
