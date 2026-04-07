# load and wrangle data sets from USGS from 1990s, and recent samples from USGS (Flint) and USFWS (2018-2022)
# create variables for exposure over background to Pb (> 0.2 ppm) and exposure at clinically toxic levels (> 1.0 ppm) as binary variables
# fit simple logistic regression models to get site- and sampling period (arrival, incubation, brood rearing) estimates of proportion of samples above thresholds with 95% CI
# create plots showing site*sampling period bar graphs for 1990s and 2019-2022 exposure above background and clinically toxic exposure

# got packages?
library(ggplot2)
library(emmeans)

# USGS data from 1990s
# load data
nineties <- read.csv("data/sourced/usgs/waterfowl_leadExposure_alaskaRussia_petersen_1993-1999.csv", header = T)

# take a look
str(nineties)

# location codes for usgs data from metadata file
code <- seq(1,13,1) # coded site values range 1-13 for each of 13 sites
# site names in sequence of code numbers, 1-13
loc <- c("Aknerkochik", "Aphrewn", "Kashunuk", "Kigigak", "Manokinak", "Naskonat", "Opyagyarak", "BigSloug", "Tutakoke", "YKD", "Indigirka", "Colville", "Prudhoe")
# combine codes with site names for a reference table
sites <-cbind(code, loc)

# format data
# format date
nineties$date <- as.POSIXct(nineties$date, format = "%Y-%m-%d")
# create julian day
nineties$julian <- format(nineties$date, "%j")
# create month
nineties$month <- format(nineties$date, format = "%m")
# format lead values as numeric
nineties$pb <- as.numeric(nineties$pb)
# recode sites
nineties$site[nineties$location == 1] <- "Aknerkochik"
nineties$site[nineties$location == 2] <- "Aphrewn"
nineties$site[nineties$location == 3] <- "Kashunuk"
nineties$site[nineties$location == 4] <- "Kigigak"
nineties$site[nineties$location == 5] <- "Manokinak"
nineties$site[nineties$location == 6] <- "Naskonat"
nineties$site[nineties$location == 7] <- "Opyagyarak"
nineties$site[nineties$location == 8] <- "BigSlough"
nineties$site[nineties$location == 9] <- "Tutakoke"
nineties$site[nineties$location == 10] <- "YKD"
nineties$site[nineties$location == 11] <- "Indigirka"
nineties$site[nineties$location == 12] <- "Colville"
nineties$site[nineties$location == 13] <- "Prudhoe"

# subset to only Kash, Kig, and Tut which have samples from 1990s and 2018-2022
keep <- c(3,4,9) # 3 = Kash, 4 = Kig, 9 = Tut
ykd90s <- nineties[nineties$location %in% keep, ]

# create sampling period based on status codes already in the data, use fewer categories (group failed, hatch, and nest as incubation, and brood and duckling as brood)
ykd90s$period[ykd90s$status == "spring"] <- "arrival"
ykd90s$period[ykd90s$status == "failed"] <- "incubation"
ykd90s$period[ykd90s$status == "hatch"] <- "incubation"
ykd90s$period[ykd90s$status == "nest"] <- "incubation"
ykd90s$period[ykd90s$status == "brood"] <- "brood"
ykd90s$period[ykd90s$status == "duckling"] <- "brood"

# create variables for Pb levels indicating exposure, and clinical toxicity (See Paine and Green 2015 and refs therein)
ykd90s$exposed <- ifelse(ykd90s$pb > 0.2,1,0)
ykd90s$clinical <- ifelse(ykd90s$pb > 1.0,1,0)

# remove records without Pb values data (NAs)
data_hist <- ykd90s[complete.cases(ykd90s$pb), ]

# subset to SPEI data (remove other species in USGS 1990s data)
data_hist <- subset(data_hist, species == "Spectacled eider")

# fit logistic regression models to get proportion of exposured samples with 95% CI by sampling period and age
# Exposure above background with threshold > 0.20 ppm)
# adults at arrival
# subset data for just arrival
data_hist_arrv <- subset(data_hist, period == "arrival")
# fit log reg model
glm_hist_arrv_exp <- glm(formula = exposed ~ 1 , family = binomial(link = "logit"), data = data_hist_arrv) # intercept only, only one sampling site in 1990s, Kashunuk
summary(glm_hist_arrv_exp)
# get marginal means
means_hist_arrv_exp <- emmeans(glm_hist_arrv_exp, ~ 1, trans = "response")
means_hist_arrv_exp

# assemble into data frame
df_exp_prob <- data.frame(data = character(),
                                 Period = character(),
                                 period_sampling = character(),
                                 site = character(),
                                 sex = character(),
                                 age = character(),
                                 cat_pb_exp = character(),
                                 mean_exp = numeric(),
                                 lcl_exp = numeric(),
                                 ucl_exp = numeric(),
                                 n_sample = numeric())

# add estimates to table
vec_hist_arrv_exp <- c("USGS", "1990s", "arrival", "Kashunuk", "both", "ahy", "exposed", 0.25, 0.16, 0.34, 88)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_hist_arrv_exp

# adults at incubation
data_hist_inc <- subset(data_hist, period == "incubation")
# fit model
glm_hist_incub_exp <- glm(formula = exposed ~ site , family = binomial(link = "logit"), data = data_hist_inc)
summary(glm_hist_incub_exp)
# sample sizes
n_hist_inc_exp <- aggregate(data_hist_inc$exposed ~ data_hist_inc$site, FUN = length)

# get marginal means
means_hist_inc_exp <- emmeans(glm_hist_incub_exp, ~ site, trans = "response")
means_hist_inc_exp
vec_hist_inc_kash_exp <- c("USGS", "1990s", "incubation", "Kashunuk", "female", "ahy", "exposed",0.3203, 0.2657, 0.375, 281)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_hist_inc_kash_exp
vec_hist_inc_kig_exp <- c("USGS", "1990s", "incubation", "Kigigak", "female", "ahy", "exposed",0.1698, 0.1246, 0.215, 265)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_hist_inc_kig_exp
vec_hist_inc_tut_exp <- c("USGS", "1990s", "incubation", "Tutakoke", "female", "ahy", "exposed",0.0417, 0, 0.122, 24)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_hist_inc_tut_exp

# adults at brood rearing
data_hist_brood_ad <- subset(data_hist, period == "brood" & age != "loc")
# fit model
glm_hist_brood_ad_exp <- glm(formula = exposed ~ site , family = binomial(link = "logit"), data = data_hist_brood_ad)
summary(glm_hist_brood_ad_exp)
# sample sizes
n_hist_brood_ad_exp <- aggregate(data_hist_brood_ad$exposed ~ data_hist_brood_ad$site, FUN = length)
n_hist_brood_ad_exp
# get marginal means
means_hist_brood_ad_exp <- emmeans(glm_hist_brood_ad_exp, ~ site, trans = "response")
means_hist_brood_ad_exp
vec_hist_kash_brood_ad_exp <- c("USGS", "1990s", "brood-adult", "Kashunuk", "female", "ahy", "exposed",0.486, 0.3924, 0.580, 109)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_hist_kash_brood_ad_exp
vec_hist_kig_brood_ad_exp <- c("USGS", "1990s", "brood-adult", "Kigigak", "female", "ahy", "exposed",0.143, 0, 0.326, 14)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_hist_kig_brood_ad_exp

# ducklings at brood rearing
data_hist_brood_loc <- subset(data_hist, period == "brood" & age == "loc")
# fit model
glm_hist_brood_loc_exp <- glm(formula = exposed ~ 1 , family = binomial(link = "logit"), data = data_hist_brood_loc)
summary(glm_hist_brood_loc_exp)
# sample sizes
n_hist_brood_loc_exp <- aggregate(data_hist_brood_loc$exposed ~ data_hist_brood_loc$site, FUN = length)
n_hist_brood_loc_exp
# get marginal means
means_hist_brood_loc_exp <- emmeans(glm_hist_brood_loc_exp, ~ 1, trans = "response")
means_hist_brood_loc_exp
vec_hist_inc_kash_loc_exp <- c("USGS", "1990s", "brood-juv", "Kashunuk", "both", "local", "exposed",0.412, 0.277, 0.547, 51)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_hist_inc_kash_loc_exp



# prob TOX-level exposure with threshold 1.0 ppm ###
# adults at arrival
# fit model
glm_hist_arrv_tox <- glm(formula = clinical ~ 1 , family = binomial(link = "logit"), data = data_hist_arrv)
summary(glm_hist_arrv_tox)
# get marginal means
means_hist_arrv_tox <- emmeans(glm_hist_arrv_tox, ~ 1, regrid  = "response")
means_hist_arrv_tox
vec_hist_arrv_tox <- c("USGS", "1990s", "arrival", "Kashunuk", "both", "ahy", "toxic", 0.102, 0.039, 0.166, 88)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_hist_arrv_tox

# adults at incubation
# fit model
glm_hist_incub_tox <- glm(formula = clinical ~ site , family = binomial(link = "logit"), data = data_hist_inc)
summary(glm_hist_incub_tox)
# get marginal means
means_hist_inc_tox <- emmeans(glm_hist_incub_tox, ~ site, trans = "response")
means_hist_inc_tox
vec_hist_inc_kash_tox <- c("USGS", "1990s", "incubation", "Kashunuk", "female", "ahy", "toxic",0.1459, 0.1046, 0.187, 281)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_hist_inc_kash_tox
vec_hist_inc_kig_tox <- c("USGS", "1990s", "incubation", "Kigigak", "female", "ahy", "toxic",0.1019, 0.0655, 0.138, 265)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_hist_inc_kig_tox
vec_hist_inc_tut_tox <- c("USGS", "1990s", "incubation", "Tutakoke", "female", "ahy", "toxic",0.0417, 0, 0.122, 24)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_hist_inc_tut_tox

# adults at brood rearing
# fit model
glm_brood_hist_ad_tox <- glm(formula = clinical ~ site , family = binomial(link = "logit"), data = data_hist_brood_ad)
summary(glm_brood_hist_ad_tox)
# get marginal means
means_hist_brood_ad_tox <- emmeans(glm_brood_hist_ad_tox, ~ site, trans = "response")
means_hist_brood_ad_tox
vec_hist_inc_kash_ad_tox <- c("USGS", "1990s", "brood-adult", "Kashunuk", "female", "ahy", "toxic",0.147, 0.00804, 0.213, 109)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_hist_inc_kash_ad_tox
vec_hist_inc_kig_ad_tox <- c("USGS", "1990s", "brood-adult", "Kigigak", "female", "ahy", "toxic",0, 0, 0, 14)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_hist_inc_kig_ad_tox


# ducklings at brood rearing
# fit model
glm_hist_brood_loc_tox <- glm(formula = clinical ~ 1 , family = binomial(link = "logit"), data = data_hist_brood_loc)
summary(glm_hist_brood_loc_tox)
# get marginal means
means_hist_brood_loc_tox <- emmeans(glm_hist_brood_loc_tox, ~ 1, trans = "response")
means_hist_brood_loc_tox
vec_hist_inc_kash_loc_tox <- c("USGS", "1990s", "brood-juv", "Kashunuk", "both", "local", "toxic",0.137, 0.0428, 0.232, 51)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_hist_inc_kash_loc_tox

# ************Recent data ####
# Load recent dataset: USGS Kashunuk 2022, FWS 2018-2022
# Kashunuk Incubation data 2022 from Flint
recent_kash_inc <- read.csv("data/sourced/usgs/waterfowl_leadExposure_alaska_flint_2022.csv", header = T)
# create site variable, all samples from Kashunuk
recent_kash_inc$site <- "Kashunuk"
# create sampling period variable, all samples during late incubation
recent_kash_inc$capture_period <- "incubation"

# subset to relevant columns
data_kash <- cbind.data.frame(recent_kash_inc$Band_ID,
              recent_kash_inc$Species,
              recent_kash_inc$Date_found,
              recent_kash_inc$site,
              recent_kash_inc$capture_period,
              recent_kash_inc$Age,
              recent_kash_inc$Sex,
              recent_kash_inc$Lead_con)
# rename columns
colnames(data_kash) <- c("band_metal", "species", "date_capture", "site", "period", "age", "sex", "pb")

# load FWS data
# FWS data 2018-2022 from Kig, Tut, Kash
recent_fws_20182022 <- read.csv("data/eider_lead_fws_2018-2022.csv", header = TRUE)
# cbind relevant columns
data_fws <- cbind.data.frame(recent_fws_20182022$band_metal,
                recent_fws_20182022$species,
                recent_fws_20182022$date_capture,
                recent_fws_20182022$site_study,
                recent_fws_20182022$capture_period,
                recent_fws_20182022$age,
                recent_fws_20182022$sex,
                recent_fws_20182022$lead_ppm)
# rename columns
colnames(data_fws) <- c("band_metal", "species", "date_capture", "site", "period", "age", "sex", "pb")

# rbind 2022 Kashunuk data (Flint) with FWS 2018-2022 data to create data frame of recent blood samples
df_recent <- rbind.data.frame(data_kash, data_fws)

# subset to only SPEI, remove the one COEI & 2 KIEI
df_recent <- subset(df_recent, species == "Spectacled eider" | species == "SPEI")

# create Pb exposure column with 2.0 ppm threshold
df_recent$exposed <- ifelse(df_recent$pb > 0.2,1,0)
# create Pb clinical toxicity column with 1.0 threshold
df_recent$clinical <- ifelse(df_recent$pb > 1.0,1,0)

# *************reformat historical data code for recent data estimates
# fit logistic regression models to get proportion of exposed samples with 95% CI by sampling period and age
# Exposure above background with threshold > 0.20 ppm)
# adults at arrival
# subset data for just arrival
data_rcnt_arrv <- subset(df_recent, period == "arrival")
# fit log reg model
glm_rcnt_arrv_exp <- glm(formula = exposed ~ site , family = binomial(link = "logit"), data = data_rcnt_arrv) #
summary(glm_rcnt_arrv_exp)
# get sample sizes
n_recent_arrv_exp <- aggregate(data_rcnt_arrv$exposed ~ data_rcnt_arrv$site, FUN = length)
n_recent_arrv_exp
# get marginal means
means_rcnt_arrv_exp <- emmeans(glm_rcnt_arrv_exp, ~ site, trans = "response")
means_rcnt_arrv_exp

# add estimates to table
vec_rcnt_arrv_kig_exp <- c("FWS", "2018-2022", "arrival", "Kigigak", "both", "ahy", "exposed", 0.025, 0.023, 0.073, 40)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_rcnt_arrv_kig_exp
vec_rcnt_arrv_tesh_exp <- c("FWS", "2018-2022", "arrival", "Teshekpuk", "both", "ahy", "exposed", 0, 0, 0, 2)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_rcnt_arrv_tesh_exp
vec_rcnt_arrv_tut_exp <- c("FWS", "2018-2022", "arrival", "Tutakoke", "both", "ahy", "exposed", 0.022, 0, 0.064, 46)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_rcnt_arrv_tut_exp
vec_rcnt_arrv_utq_exp <- c("FWS", "2018-2022", "arrival", "Utqiagvik", "both", "ahy", "exposed", 0, 0, 0, 9)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_rcnt_arrv_utq_exp

# adults at incubation
data_rcnt_inc <- subset(df_recent, period == "incubation")
# fit model
glm_rcnt_incub_exp <- glm(formula = exposed ~ site , family = binomial(link = "logit"), data = data_rcnt_inc)
summary(glm_rcnt_incub_exp)
# get sample sizes
n_recent_inc_exp <- aggregate(data_rcnt_inc$exposed ~ data_rcnt_inc$site, FUN = length)
n_recent_inc_exp
# get marginal means
means_rcnt_inc_exp <- emmeans(glm_rcnt_incub_exp, ~ site, trans = "response")
means_rcnt_inc_exp
vec_rcnt_inc_kash_exp <- c("USGS", "2018-2022", "incubation", "Kashunuk", "female", "ahy", "exposed",0.243, 0.105, 0.381, 37)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_rcnt_inc_kash_exp
vec_rcnt_inc_kig_exp <- c("FWS", "2018-2022", "incubation", "Kigigak", "female", "ahy", "exposed",0.167, 0.084, 0.249, 78)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_rcnt_inc_kig_exp
vec_rcnt_inc_utq_exp <- c("FWS", "2018-2022", "incubation", "Utqiagvik", "female", "ahy", "exposed",0, 0, 0, 3)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_rcnt_inc_utq_exp

# adults at brood rearing
data_rcnt_brood_ad <- subset(df_recent, period == "brood" & age != "local")
# fit model
glm_rcnt_brood_ad_exp <- glm(formula = exposed ~ 1 , family = binomial(link = "logit"), data = data_rcnt_brood_ad)
summary(glm_rcnt_brood_ad_exp)
# get sample sizes
n_recent_brood_ad_exp <- aggregate(data_rcnt_brood_ad$exposed ~ data_rcnt_brood_ad$site, FUN = length)
n_recent_brood_ad_exp
# get marginal means
means_rcnt_brood_ad_exp <- emmeans(glm_rcnt_brood_ad_exp, ~ 1, trans = "response")
means_rcnt_brood_ad_exp
vec_rcnt_kash_brood_ad_exp <- c("FWS", "2018-2022", "brood-adult", "Kashunuk", "female", "ahy", "exposed",0.20, 0, 0.45, 10)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_rcnt_kash_brood_ad_exp

# ducklings at brood rearing
data_rcnt_brood_loc <- subset(df_recent, period == "brood" & age == "local")
# fit model
glm_rcnt_brood_loc_exp <- glm(formula = exposed ~ 1 , family = binomial(link = "logit"), data = data_rcnt_brood_loc)
summary(glm_rcnt_brood_loc_exp)
# get sample sizes
n_recent_brood_loc_exp <- aggregate(data_rcnt_brood_loc$exposed ~ data_rcnt_brood_loc$site, FUN = length)
n_recent_brood_loc_exp
# get marginal means
means_rcnt_brood_loc_exp <- emmeans(glm_rcnt_brood_loc_exp, ~ 1, trans = "response")
means_rcnt_brood_loc_exp
vec_rcnt_inc_kash_loc_exp <- c("FWS", "2018-2022", "brood-juv", "Kashunuk", "both", "local", "exposed",0.154, 0, 0.35, 13)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_rcnt_inc_kash_loc_exp


# prob TOX-level exposure with threshold 1.0 ppm ###
# adults at arrival
# fit model
glm_rcnt_arrv_tox <- glm(formula = clinical ~ site , family = binomial(link = "logit"), data = data_rcnt_arrv)
summary(glm_rcnt_arrv_tox)
# get sample sizes
n_recent_arrv_tox <- aggregate(data_rcnt_arrv$exposed ~ data_rcnt_arrv$site, FUN = length)
# get marginal means
means_rcnt_arrv_tox <- emmeans(glm_rcnt_arrv_tox, ~ site, regrid  = "response")
means_rcnt_arrv_tox
vec_rcnt_arrv_kig_tox <- c("FWS", "2018-2022", "arrival", "Kigigak", "both", "ahy", "toxic", 0.025, 0, 0.073, 40)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_rcnt_arrv_kig_tox
vec_rcnt_arrv_tesh_tox <- c("FWS", "2018-2022", "arrival", "Teshekpuk", "both", "ahy", "toxic", 0, 0, 0, 2)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_rcnt_arrv_tesh_tox
vec_rcnt_arrv_tut_tox <- c("FWS", "2018-2022", "arrival", "Tutakoke", "both", "ahy", "toxic", 0, 0, 0, 46)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_rcnt_arrv_tut_tox
vec_rcnt_arrv_utq_tox <- c("FWS", "2018-2022", "arrival", "Utqiagvik", "both", "ahy", "toxic", 0, 0, 0, 9)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_rcnt_arrv_utq_tox

# adults at incubation
# fit model
glm_rcnt_incub_tox <- glm(formula = clinical ~ site , family = binomial(link = "logit"), data = data_rcnt_inc)
summary(glm_rcnt_incub_tox)
# get sample sizes
n_recent_inc_tox <- aggregate(data_rcnt_inc$exposed ~ data_rcnt_inc$site, FUN = length)
# get marginal means
means_rcnt_inc_tox <- emmeans(glm_rcnt_incub_tox, ~ site, trans = "response")
means_rcnt_inc_tox
vec_rcnt_inc_kash_tox <- c("USGS", "2018-2022", "incubation", "Kashunuk", "female", "ahy", "toxic",0.081, 0, 0.169, 37)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_rcnt_inc_kash_tox
vec_rcnt_inc_kig_tox <- c("FWS", "2018-2022", "incubation", "Kigigak", "female", "ahy", "toxic",0.077, 0.018, 0.136, 78)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_rcnt_inc_kig_tox
vec_rcnt_inc_tut_tox <- c("FWS", "2018-2022", "incubation", "Utqiagvik", "female", "ahy", "toxic",0, 0, 0, 3)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_rcnt_inc_tut_tox

# adults at brood rearing
# fit model
glm_brood_rcnt_ad_tox <- glm(formula = clinical ~ 1 , family = binomial(link = "logit"), data = data_rcnt_brood_ad)
summary(glm_brood_rcnt_ad_tox)
# get sample sizes
n_recent_brood_ad_tox <- aggregate(data_rcnt_brood_ad$exposed ~ data_rcnt_brood_ad$site, FUN = length)
n_recent_brood_ad_tox
# get marginal means
means_rcnt_brood_ad_tox <- emmeans(glm_brood_rcnt_ad_tox, ~ 1, trans = "response")
means_rcnt_brood_ad_tox
vec_rcnt_inc_kash_ad_tox <- c("FWS", "2018-2022", "brood-adult", "Kashunuk", "female", "ahy", "toxic",0, 0, 0, 10)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_rcnt_inc_kash_ad_tox

# ducklings at brood rearing
# fit model
glm_rcnt_brood_loc_tox <- glm(formula = clinical ~ 1 , family = binomial(link = "logit"), data = data_rcnt_brood_loc)
summary(glm_rcnt_brood_loc_tox)
# get sample sizes
n_recent_brood_loc_tox <- aggregate(data_rcnt_brood_loc$exposed ~ data_rcnt_brood_loc$site, FUN = length)
n_recent_brood_loc_tox
# get marginal means
means_rcnt_brood_loc_tox <- emmeans(glm_rcnt_brood_loc_tox, ~ 1, trans = "response")
means_rcnt_brood_loc_tox
vec_rcnt_inc_kash_loc_tox <- c("FWS", "2018-2022", "brood-juv", "Kashunuk", "both", "local", "toxic",0, 0, 0, 13)
df_exp_prob[nrow(df_exp_prob)+1, ] <- vec_rcnt_inc_kash_loc_tox

# calculate 95% CI halfwidths for plotting
df_exp_prob$mean_exp <- as.numeric(df_exp_prob$mean_exp)
df_exp_prob$lcl_exp <- as.numeric(df_exp_prob$lcl_exp)
df_exp_prob$ucl_exp <- as.numeric(df_exp_prob$ucl_exp)

df_exp_prob$lower <- df_exp_prob$mean_exp - df_exp_prob$lcl_exp
df_exp_prob$upper <- df_exp_prob$ucl_exp - df_exp_prob$mean_exp


# create plots

# 1 1990s exposure during arrival at Kashunuk ####
# subset data to historic arrival
data_hist_kash <- subset(df_exp_prob, site == "Kashunuk" & Period == "1990s" & cat_pb_exp == "exposed")
# create var to set order of bars in plot
data_hist_kash$order <- seq(1,4,1)
# plot
p.hist_exp_kash <- ggplot(data_hist_kash) +
  geom_bar(aes(reorder(x = period_sampling, order), y = mean_exp), stat = "identity", width = 0.5, fill ="royalblue") +
  geom_errorbar(aes(x = period_sampling, ymin = lcl_exp, ymax = ucl_exp), width = 0.0) +
  xlab("Sampling Period") +
  ylab("Proportion Exposed to Lead > 0.2 ppm") +
  ylim(0, 0.75)+
  geom_text(aes(period_sampling, mean_exp, label = paste("n = ", n_sample), y = 0.01), size = 2.0, color = "white") +
  geom_text(aes(period_sampling, mean_exp, label = round(mean_exp, digits = 2.0), vjust = -8.0), size = 4.0) +
  ggtitle("Kashunuk 1990s")+
  theme_bw()
p.hist_exp_kash
# save it
ggsave("output/fig_hist_exp_kash.jpg", width = 10, height = 10, units = "cm")

# 2 1990s exposure during incubation at Kigigak ####
data_hist_kig <- subset(df_exp_prob, site == "Kigigak" & Period == "1990s" & cat_pb_exp == "exposed")
# add sampling periods with no data
vec_hist_kig_arrv_nd <- c("USGS", "1990s", "arrival", "Kigigak", "both", "AHY", "exposed", "","", "", "", "", "")
data_hist_kig[nrow(data_hist_kig)+1, ] <- vec_hist_kig_arrv_nd
vec_hist_kig_brood_loc_nd <- c("USGS", "1990s", "brood-juv", "Kigigak", "both", "exposed", "toxic", "","", "", "", "", "")
data_hist_kig[nrow(data_hist_kig)+1, ] <- vec_hist_kig_brood_loc_nd
# format
data_hist_kig$mean_exp <- as.numeric(data_hist_kig$mean_exp)
data_hist_kig$lcl_exp <- as.numeric(data_hist_kig$lcl_exp)
data_hist_kig$ucl_exp <- as.numeric(data_hist_kig$ucl_exp)
# create var to set order of bars in plot
data_hist_kig$order <- c(2,3,1,4)
# plot
p.hist_exp_kig <- ggplot(data_hist_kig) +
  geom_bar(aes(reorder(x = period_sampling, order), y = mean_exp), stat = "identity", width = 0.5, fill = "royalblue") +
  geom_errorbar(aes(x = period_sampling, ymin = lcl_exp, ymax = ucl_exp), width = 0.0) +
  scale_x_discrete(drop = FALSE) +
  xlab("Sampling Period") +
  ylab("Proportion Exposed to Lead > 0.2 ppm") +
  ylim(0, 0.75)+
  geom_text(aes(period_sampling, mean_exp, label = paste("n = ", n_sample), y =0.01), size = 2.0, color = "white") +
  geom_text(aes(period_sampling, mean_exp, label = round(mean_exp, digits = 2.0), vjust = -10.0), size = 4.0) +
  ggtitle("Kigigak 1990s")+
  theme_bw()
p.hist_exp_kig
# save it
ggsave("output/fig_hist_exp_kig.jpg", width = 10, height = 10, units = "cm")


# 3 1990s exposure during incubation at Tutakoke
data_hist_tut <- subset(df_exp_prob, site == "Tutakoke" & Period == "1990s" & cat_pb_exp == "exposed")
# add sampling periods with no data
vec_hist_tut_arrv_nd <- c("USGS", "1990s", "arrival", "Tutakoke", "both", "ahy", "exposed", "","", "", "", "", "")
data_hist_tut[nrow(data_hist_tut)+1, ] <- vec_hist_tut_arrv_nd
vec_hist_kig_brood_ad_nd <- c("USGS", "1990s", "brood-adult", "Tutakoke", "female", "ahy", "exposed", "","", "", "", "", "")
data_hist_tut[nrow(data_hist_tut)+1, ] <- vec_hist_kig_brood_ad_nd
vec_hist_tut_brood_loc_nd <- c("USGS", "1990s", "brood-juv", "Tutakoke", "both", "local", "exposed", "","", "", "", "", "")
data_hist_tut[nrow(data_hist_tut)+1, ] <- vec_hist_tut_brood_loc_nd
# format
data_hist_tut$mean_exp <- as.numeric(data_hist_tut$mean_exp)
data_hist_tut$lcl_exp <- as.numeric(data_hist_tut$lcl_exp)
data_hist_tut$ucl_exp <- as.numeric(data_hist_tut$ucl_exp)
# create var to set order of bars in plot
data_hist_tut$order <- c(2,1,3,4)
# plot
p.hist_exp_tut <- ggplot(data_hist_tut) +
  geom_bar(aes(reorder(x = period_sampling, order), y = mean_exp), stat = "identity", width = 0.5, fill = "royalblue") +
  geom_errorbar(aes(x = period_sampling, ymin = lcl_exp, ymax = ucl_exp), width = 0.0) +
  scale_x_discrete(drop = FALSE) +
  xlab("Sampling Period") +
  ylab("Proportion Exposed to Lead > 0.2 ppm") +
  ylim(0, 0.75)+
  geom_text(aes(period_sampling, mean_exp, label = paste("n = ", n_sample), y = 0.01), size = 2.0, color = "white") +
  geom_text(aes(period_sampling, mean_exp, label = round(mean_exp, digits = 2.0), vjust = -6.0), size = 4.0) +
  ggtitle("Tutakoke 1990s")+
  theme_bw()
p.hist_exp_tut
# save it
ggsave("output/fig_hist_exp_tut.jpg", width = 10, height = 10, units = "cm")

# 4 Recent and 1990s exposure at Kig color bars by "Period"
data_all_kash <- subset(df_exp_prob, site == "Kashunuk" & cat_pb_exp == "exposed")
# vector of missing data
vec_all_kash_arrv_nd <- c("USGS", "2018-2022", "arrival", "Kashunuk", "both", "ahy", "exposed", "","", "", 0, "", "")
data_all_kash[nrow(data_all_kash)+1, ] <- vec_all_kash_arrv_nd
# create var to set order of bars in plot
data_all_kash$order <- c(1,2,3,4,6,7,8,5)
# format
data_all_kash$mean_exp <- as.numeric(data_all_kash$mean_exp)
data_all_kash$lcl_exp <- as.numeric(data_all_kash$lcl_exp)
data_all_kash$ucl_exp <- as.numeric(data_all_kash$ucl_exp)
# plot
p.all_exp_kash <- ggplot(data_all_kash, aes(fill = Period, y = mean_exp, x = reorder(x = period_sampling, order), fill = Period)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(ymin = lcl_exp, ymax = ucl_exp), width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual("Period", values = c("1990s" = "royalblue", "2018-2022" = "darkorange"))+
  xlab("Sampling Period") +
  ylab("Proportion Exposed to Lead > 0.2 ppm") +
  ylim(0, 0.75)+
  geom_text(aes(period_sampling, mean_exp, label = paste("n = ", n_sample), group = Period, y = 0.01),
            size = 1.5,position = position_dodge(width = 0.9), color = "white") +
  geom_text(aes(period_sampling, mean_exp, label = round(mean_exp, digits = 2.0), group = Period, y = 0.6), # vjust = -16
            size = 3.0, position = position_dodge(width = 0.9)) +
  ggtitle("Kashunuk") +
  theme_bw()
p.all_exp_kash
# save it
ggsave("output/fig_all_exp_kash.jpg", width = 20, height = 10, units = "cm")

# 5 Recent and 1990s exposure at Kig color bars by "Period"
# subset data Kigigak
data_all_kig <- subset(df_exp_prob, site == "Kigigak" & cat_pb_exp == "exposed")
# vector of missing data
vec_all_kig_arrv_nd <- c("USGS", "1990s", "arrival", "Kigigak", "both", "ahy", "exposed", "","", "", 0, "", "")
data_all_kig[nrow(data_all_kig)+1, ] <- vec_all_kig_arrv_nd
vec_all_kig_brood_loc_nd <- c("USGS", "1990s", "brood-juv", "Kigigak", "both", "loc", "exposed", "","", "", 0, "", "")
data_all_kig[nrow(data_all_kig)+1, ] <- vec_all_kig_brood_loc_nd
vec_all_kig_rec_brood_ad_nd <- c("FWS", "2018-2022", "brood-adult", "Kigigak", "female", "ahy", "exposed", "","", "", 0, "", "")
data_all_kig[nrow(data_all_kig)+1, ] <- vec_all_kig_rec_brood_ad_nd
vec_all_kig_rec_brood_loc_nd <- c("FWS", "2018-2022", "brood-juv", "Kigigak", "both", "loc", "exposed", "","", "", 0, "", "")
data_all_kig[nrow(data_all_kig)+1, ] <- vec_all_kig_rec_brood_loc_nd

# create var to set order of bars in plot
data_all_kig$order <- c(3,5,2,4,1,7,6,8)
# format
data_all_kig$mean_exp <- as.numeric(data_all_kig$mean_exp)
data_all_kig$lcl_exp <- as.numeric(data_all_kig$lcl_exp)
data_all_kig$ucl_exp <- as.numeric(data_all_kig$ucl_exp)
# plot
p.all_exp_kig <- ggplot(data_all_kig, aes(fill = Period, y = mean_exp, x = reorder(x = period_sampling, order))) +
  geom_bar(position = position_dodge(0.9), stat = "identity") +
  geom_errorbar(aes(ymin = lcl_exp, ymax = ucl_exp), width = 0.2, position = position_dodge(0.9)) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual("Period", values = c("1990s" = "royalblue", "2018-2022" = "darkorange"))+
  xlab("Sampling Period") +
  ylab("Proportion Exposed to Lead > 0.2 ppm") +
  ylim(0, 0.75)+
  geom_text(aes(period_sampling, mean_exp, label = paste("n = ", n_sample), group = Period, y = 0.01),
            size = 1.5,position = position_dodge(width = 0.9), color = "white") +
  geom_text(aes(period_sampling, mean_exp, label = round(mean_exp, digits = 2.0), group = Period, vjust = -10),
            size = 3.0, position = position_dodge(width = 0.9)) +
  ggtitle("Kigigak") +
  theme_bw()
p.all_exp_kig
# save it
ggsave("output/fig_all_exp_kig.jpg", width = 20, height = 10, units = "cm")


# Tutakoke all ####
# subset data to historic arrival
data_all_tut <- subset(df_exp_prob, site == "Tutakoke" & cat_pb_exp == "exposed")
# vector of missing data
vec_hist_tut_arrv_nd <- c("USGS", "1990s", "arrival", "Tutakoke", "both", "ahy", "exposed", "","", "", 0, "", "")
data_all_tut[nrow(data_all_tut)+1, ] <- vec_hist_tut_arrv_nd

vec_recent_tut_inc_nd <- c("USGS", "2018-2022", "incubation", "Tutakoke", "female", "ahy", "exposed", "","", "", 0, "", "")
data_all_tut[nrow(data_all_tut)+1, ] <- vec_recent_tut_inc_nd

vec_hist_tut_rec_brood_ad_nd <- c("USGS", "1990s", "brood-adult", "Tutakoke", "female", "ahy", "exposed", "","", "", 0, "", "")
data_all_tut[nrow(data_all_tut)+1, ] <- vec_hist_tut_rec_brood_ad_nd
vec_recent_tut_rec_brood_loc_nd <- c("FWS", "2018-2022", "brood-adult", "Tutakoke", "female", "ahy", "exposed", "","", "", 0, "", "")
data_all_tut[nrow(data_all_tut)+1, ] <- vec_recent_tut_rec_brood_loc_nd

vec_hist_tut_rec_brood_loc_nd <- c("USGS", "1990s", "brood-juv", "Tutakoke", "both", "loc", "exposed", "","", "", 0, "", "")
data_all_tut[nrow(data_all_tut)+1, ] <- vec_hist_tut_rec_brood_loc_nd
vec_recent_tut_rec_brood_ad_nd <- c("FWS", "2018-2022", "brood-juv", "Tutakoke", "both", "loc", "exposed", "","", "", 0, "", "")
data_all_tut[nrow(data_all_tut)+1, ] <- vec_recent_tut_rec_brood_ad_nd

# create var to set order of bars in plot
data_all_tut$order <- c(3,5,2,4,1,7,6,8)
# format
data_all_tut$mean_exp <- as.numeric(data_all_tut$mean_exp)
data_all_tut$lcl_exp <- as.numeric(data_all_tut$lcl_exp)
data_all_tut$ucl_exp <- as.numeric(data_all_tut$ucl_exp)
# plot
p.all_exp_tut <- ggplot(data_all_tut, aes(fill = Period, y = mean_exp, x = reorder(x = period_sampling, order))) +
  geom_bar(position = position_dodge(0.9), stat = "identity") +
  geom_errorbar(aes(ymin = lcl_exp, ymax = ucl_exp), width = 0.2, position = position_dodge(0.9)) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual("Period", values = c("1990s" = "royalblue", "2018-2022" = "darkorange"))+
  xlab("Sampling Period") +
  ylab("Proportion Exposed to Lead > 0.2 ppm") +
  ylim(0, 0.75)+
  geom_text(aes(period_sampling, mean_exp, label = paste("n = ", n_sample), group = Period, y = 0.01),
            size = 1.5,position = position_dodge(width = 0.9), color = "white") +
  geom_text(aes(period_sampling, mean_exp, label = round(mean_exp, digits = 2.0), group = Period, vjust = -10),
            size = 3.0, position = position_dodge(width = 0.9)) +
  ggtitle("Tutakoke") +
  theme_bw()
p.all_exp_tut
# save it
ggsave("output/fig_all_exp_tut.jpg", width = 20, height = 10, units = "cm")
