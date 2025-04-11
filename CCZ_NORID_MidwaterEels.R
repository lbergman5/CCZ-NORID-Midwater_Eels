```{r}
### Vertical migration and campaign comparison
##Linear models to examine differences in vertical distribution between day/night and C5b/C5e per eel family
dat <- read.csv("Eel_distrubution_data.csv", header = TRUE, sep = ",")
head(dat)
names(dat)

#check if depth is normally distributed
shapiro.test(dat$depth) #it is not, so double check diagnostic plots if ok

# Convert factors
dat$day_night <- factor(dat$day_night)
dat$Family <- factor(dat$Family)
dat$cruise <- factor(dat$cruise)

# Subset data
dat_nemi <- subset(dat, Family == "Nemichthyidae")
dat_serri <- subset(dat, Family == "Serrivomeridae")

# Fit linear models
lm_nemi <- lm(depth ~ day_night + cruise , data = dat_nemi)
lm_serri <- lm(depth ~ day_night + cruise, data = dat_serri)

# Summarize models
summary(lm_nemi)
summary(lm_serri)

# DHARMa diagnostics
library(DHARMa)
sim_nemi <- simulateResiduals(fittedModel = lm_nemi)
sim_serri <- simulateResiduals(fittedModel = lm_serri)

# Plot diagnostics
plot(sim_nemi)
plot(sim_serri)

# Export Plot diagnostics (see Supplementary Figure S1)
pdf("lm_diagnostic_plots_nemi_cruise.pdf")
plot(sim_nemi) #KS not good
dev.off()

pdf("lm_diagnostic_plots_serri_cruise.pdf")
plot(sim_serri) #diagnostic plots ok
dev.off()

###Abiotic drivers of eel distribution
##VIF for variables to check for multicollinearity
dat <- read.csv("Eel_distrubution_data.csv", header = TRUE, sep = ",")
library(usdm)

#with only 5 variables that don't have missing variables (removing flurorescence and turbidity)
selected_vars <- c("depth", "temperature", "salinity", "dissolved_oxygen", "latitude")
dat_vif <- dat[, selected_vars]
# Calculate VIF using the usdm package
vif_result <- vif(dat_vif)
# Print the VIF results
print(vif_result)

##GAMM for examining differences between eel families
library(gamm4)
library(dplyr)
library(DHARMa)

# Prepare data
selected_vars <- c("depth", "temperature_celsius", "salinity_ppt", "dissolved_oxygen_ml_L", "latitude", "dive")
dat_clean <- dat %>%
  select(Family, all_of(selected_vars)) %>%
  na.omit()
str(dat_clean)
dat_clean$dive <- factor(dat_clean$dive)

# Convert Family to binary
dat_clean$Family <- as.numeric(dat_clean$Family == "Serrivomeridae")
head(dat_clean)

# Fit GAMM models stepwise 
gamm1 <- gamm4(Family ~ s(depth) + s(temperature_celsius) + s(salinity_ppt) + s(dissolved_oxygen_ml_L) + s(latitude),
               random = ~(1 | dive), data = dat_clean, family = binomial)
summary(gamm1$gam)
gamm2 <- gamm4(Family ~ s(depth) + s(salinity_ppt) + s(dissolved_oxygen_ml_L) + s(latitude),
               random = ~(1 | dive), data = dat_clean, family = binomial)
summary(gamm2$gam)
gamm3 <- gamm4(Family ~ s(depth) + s(dissolved_oxygen_ml_L) + s(latitude),
               random = ~(1 | dive), data = dat_clean, family = binomial)
summary(gamm3$gam)
gamm4 <- gamm4(Family ~ s(depth) + s(dissolved_oxygen_ml_L),
               random = ~(1 | dive), data = dat_clean, family = binomial)
summary(gamm4$gam)

# DHARMa diagnostics for best model
summary(gamm2$mer)
simres <- simulateResiduals(fittedModel = gamm2$mer)
#export diagnostic plots (see Supplementary Figure S2)
pdf("gamm4_diagnostic_plots.pdf")
plot(simres)
dev.off()

# R²
library(performance)
citation("performance")
packageVersion("performance")

r2(gamm2$mer)
lme4::isSingular(gamm2$mer)
r2(gamm2$mer, tolerance = 1e-6)




```{r}
###Prey and indicator taxa
#Calculating diversity for prey using Hill no
library(ggplot2)
library(vegan)
library(indicspecies)
library(RVAideMemoire)
library(rstatix)
library(agricolae)
library(iNEXT)

dat <- read.csv("Prey_transect_data.csv", header = TRUE, sep = ",", row.names=1)
head(dat)
hillcalc<-iNEXT(dat, q=1, datatype="abundance")
hillcalc

#calculating differences in prey diversity between groups

dat <-read.csv("Prey_diversity_data.csv", header = TRUE, sep = ",")
head(dat)

order_levels <- c("shallower_day", "shallower_night", "deeper_day", "deeper_night")
dat$depth_group <- factor(dat$depth_group, levels = order_levels)

pdf("Shannon_by_four_depth_and_daytime_bins.pdf")
ggplot(dat, aes(x = depth_group, y = Hill_Shannon_observed, fill = depth_group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("shallower_day" = "white", "shallower_night" = "black", "deeper_day"= "white", "deeper_night" = "black")) +
  labs(title = "Shannon Index",
       x = "Groups",
       y = "Shannon Index",
       fill = "Dataset")
dev.off()

shapiro.test(dat$Hill_Shannon_observed) #p < 0.05, not normal distribution

# Kruskal-Wallis test
kruskal_result <- kruskal.test(Hill_Shannon_observed ~ depth_group, dat)
print(kruskal_result)
# Post-hoc analysis with agricolae
kruskal_groups <- with(dat, kruskal(Hill_Shannon_observed, depth_group, console = TRUE))

group_vector_day_night <- factor(prey_observations_metadata$day.night)
dispersion_day_night <- betadisper(bray_curtis_dist, group_vector_day_night)
summary(dispersion_day_night)
dispersion_anova_day_night <- anova(dispersion_day_night)
print(dispersion_anova_day_night) #p value was not significant, meaning that the assumption of homogeneity of dispersion was met

permanova_result_day_night <- adonis2(bray_curtis_dist ~ group_vector_day_night, data = as.data.frame(group_vector_day_night), permutations = 999)
  print(permanova_result_day_night) #not significant

# Perform betadisper and PERMANOVA for depth_group
group_vector_depth_group <- factor(prey_observations_metadata$depth_group)
dispersion_depth_group <- betadisper(bray_curtis_dist, group_vector_depth_group)
summary(dispersion_depth_group)
dispersion_anova_depth_group <- anova(dispersion_depth_group)
print(dispersion_anova_depth_group) #p value was not significant, meaning that the assumption of homogeneity of dispersion was met

# Perform PERMANOVA for depth_group if dispersion is not significant
  permanova_result_depth_group <- adonis2(bray_curtis_dist ~ group_vector_depth_group, data = as.data.frame(group_vector_depth_group), permutations = 999)
  print(permanova_result_depth_group)
  # Perform pairwise comparisons for depth_group with FDR adjustment if PERMANOVA is significant
    pairwise_results_depth_group <- pairwise.perm.manova(bray_curtis_dist, group_vector_depth_group, p.adjust.method = "fdr")
    print(pairwise_results_depth_group) #still not significant
```


```{r}
#calculating indicator species for prey
#groups are day_night, with night = 1, day = 2
dat <- read.csv("Prey_transect_data_indsp_dn.csv", header = TRUE, sep = ",", row.names = 1)
indval <- multipatt(dat, dat$group, control = how(nperm=999))
summary(indval)
str(indval)

#do the same with depth groups, with shallow = 1 and deep = 2
dat <- read.csv("Prey_transect_data_indsp_sd.csv", header = TRUE, sep = ",", row.names = 1)
indval <- multipatt(dat, dat$group, control = how(nperm=999))
summary(indval)
```

```{r}
###Behavior analysis
library(tidymodels)
library(tidyverse)

#species: nemichthyidae = 0; serrivomeridae = 1
#activity: resting = 0; swimming =1
#orientation: horizontal = 0; vertical = 1
glm_behavior<-read.csv("/cloud/project/Midwater Eels/Eel_behavior.csv")
glm_n<- subset(glm_behavior, species == "0")
glm_s<- subset(glm_behavior, species == "1")

test0<-glm(species~activity + orientation, data = glm_behavior, family = binomial(link="logit"))
summary(test0)

#stepwise removal of variables for nemichthyidae
glm1 <-glm(activity ~ depth + temperature + salinity + dissolved.oxygen..ml.L., data = glm_n, family = binomial(link="logit"))
summary(glm1)

glm2 <-glm(activity ~ depth + salinity + dissolved.oxygen..ml.L., data = glm_n, family = binomial(link="logit"))
summary(glm2)

glm3 <-glm(activity ~ salinity + dissolved.oxygen..ml.L., data = glm_n, family = binomial(link="logit"))
summary(glm3)

glm4 <-glm(activity ~ + dissolved.oxygen..ml.L., data = glm_n, family = binomial(link="logit"))
summary(glm4)

#now for serrivomeridae
glm<-glm(activity ~ depth + temperature + salinity + dissolved.oxygen..ml.L., data = glm_s, family = binomial(link="logit"))
summary(glm)

glm1<-glm(activity ~ depth + salinity + dissolved.oxygen..ml.L., data = glm_s, family = binomial(link="logit"))
summary(glm1)

glm2<-glm(activity ~ depth + salinity, data = glm_s, family = binomial(link="logit"))
summary(glm2)

glm3<-glm(activity ~ depth, data = glm_s, family = binomial(link="logit"))
summary(glm3)


# DHARMa diagnostics for best model
simulateResiduals(test0, plot = TRUE)

# R²
r2(test0)
```