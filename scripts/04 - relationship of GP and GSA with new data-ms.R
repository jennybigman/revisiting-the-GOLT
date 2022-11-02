	# 04 - with more species (132) what is the relationship between growth performance and gill 
	#surface area?

	# Question 3 in MS

	# Question: How does the relationship of growth performance and gill surface area from 
		# Pauly's data compare to one with more data?
	
	# to answer this question, I need to re-estimate gill area index for Pauly's data with 
	# d = 0.8 because he estimated GAI with an estimated d from a relationship of d and Wmax
	# (see text) and then later used d = 0.8 to estimate GAI, but did not report data
	
	# then, I can compare the slopes estimated just with Pauly's data between a relationship 
	# where GAI was estimated with Pauly's d from 1981 and where GAI was estimated with d = 0.8
	
	# finally, I can estimate the relationship of growth performance and gill surface area with
	# my data (132 species) and compare to that with Pauly's data
		# to do this, I will have to estimate GAI with the new data

	#### A. GAI comparison with Pauly's data ####

	# re-estimate GAI for Pauly's 42 species 

	HM_dat <- read_csv(here("./data/Hughes&Morgan_GSA_supplement.csv"))
	len(HM_dat$Species)

	# Hughes & Morgan 1973b dataset included more than the 42 species that Pauly included in
	# his dataset. Pauly only included marine species and those who had age & growth data
	
	# limit Hughes & Morgan 1973b dataset to just those species in Pauly's dataset
	len(Pauly_dat$Species)

	# check to make sure taxonomy matches
	setdiff(HM_dat$Species, Pauly_dat$Species)
	setdiff(Pauly_dat$Species, HM_dat$Species)

	# update taxonomy so it is the same in both lists
	Pauly_dat$Species[Pauly_dat$Species == "Euthynnus alliteratus"] <- "Euthynnus alletteratus"

	# make sure now that there are no mismatches
	which(unique(HM_dat$Species) != Pauly_dat$Species)

	# one species, "Engraulis encrasicholus" only has relative GSA and a body mass; Pauly still
	# incorporated this species in his dataset of total GSA by multiplying the relative GSA by
	# the body mass

	HM_dat <- HM_dat %>%
						mutate(GSAcm2 = case_when(
										TotalGSACM2 >0 ~ TotalGSACM2,
										is.na(TotalGSACM2) ~ (BodyMassG * RelativeGSACM2_G)
						)) 
							
	# If more than one measure of GSA was reported in the Hughes & Morgan 1973b dataset, Pauly
	# took a mean, so average GSA by species
	
	HM_dat <- HM_dat %>%
					  dplyr::select(Species, BodyMassG, GSAcm2) %>%
						group_by(Species) %>%
						summarise(mean_GSAcm2 = mean(GSAcm2),
											mean_BodyMassG = mean(BodyMassG))
	
	HM_dat$mean_BodyMassG <- round(HM_dat$mean_BodyMassG, 2)
	HM_dat$mean_GSAcm2 <- round(HM_dat$mean_GSAcm2, 2)


	### for reference, lm results with coelacanth
	# results with the coelacanth
	
	# pauly_fit_coel <- lm(P ~ log10(GAI), data = Pauly_dat)
	# summary(pauly_fit_coel)
	# confint(pauly_fit_coel)
		# slope = 0.43 (-0.22 to 1.09)
		# intercept = 2.23 (1.50 to 2.95)

	# estimate gill area index as GAI = G/W^d, where G = gill area in cm2, W = body weight,
	# and d = 0.8
	
	GAI_formula <- function(G, W){
	  GAI = G / (W^0.8)
	  return(GAI)
	}


	HM_dat <- HM_dat %>%
	          rowwise() %>%
	          mutate(GAI_setd = GAI_formula(G = mean_GSAcm2, W = mean_BodyMassG))

	# add in growth performance data originally estimated by Pauly following P = log (K * Winf) 
	HM_dat$P <- Pauly_dat[match(HM_dat$Species, Pauly_dat$Species), 'P']
	
	# log transform GAI for analyses
	HM_dat$LogGAI_setd <- log10(HM_dat$GAI_setd)
	
	### with the original GAI estimated and reported in Pauly 1981
	pauly_fit_lm_GAIsetd <- lm(P ~ LogGAI_setd, data = HM_dat)
	summary(pauly_fit_lm_GAIsetd)
	confint(pauly_fit_lm_GAIsetd)
	# slope = 0.83 (0.20 to 1.47)
	# intercept = 1.78 (1.05 to 2.51)

	# what about outliers (just check)

	# cook's distance
	aug_pauly_GAIsetd <- augment(pauly_fit_lm_GAIsetd)
	aug_pauly_GAIsetd %>%
	mutate(i = 1:n()) %>%
	ggplot(aes(x = i, y = .cooksd)) +
	geom_hline(yintercept = .5, color = "red", size = 1) +
	geom_point(alpha = .5, size = 3) +
	ylab("Cook's distance") +
	xlab("Observation")

	# one outlier with set d

	row_num <- which(aug_pauly_GAIsetd$.cooksd > 0.5)
	outlierP <- aug_pauly_GAIsetd$P[row_num]
	outlierGSA <- aug_pauly_GAIsetd$LogGAI_setd[row_num]
	row_num2 <- which((HM_dat$P == outlierP) & (HM_dat$LogGAI_setd == outlierGSA))
	HM_dat$Species[row_num2] 
	
	# check pareto k
	#BR_setd <- 
  # brm(data = HM_dat, family = gaussian,
  #     P ~ LogGAI_setd, 
  #     prior = c(prior(student_t(3, 2.5, 2.5), class = Intercept),
  #               prior(student_t(3, 0, 10), class = b),
  #               prior(student_t(3, 0, 2.5),  class = sigma)),
  #     chains = 4,
  #     save_all_pars = TRUE,
  #     iter = 5000, warmup = 2000)
 
  #summary(BR_setd)
	#loo_BR_setd <- loo(BR_setd, moment_match = TRUE)
 
 #save(BR_setd, file = here("output/BR_setd.rds"))
 load(file = "output/BR_setd.rds")
	
 	post_BR_setd_INT <- as.data.frame(BR_setd)
	post_BR_setd_INT_mean <- stack(apply(post_BR_setd_INT, 2, mean)) %>% rename(mean = values)
	post_BR_setd_INT_LCIs <- stack(apply(post_BR_setd_INT, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_BR_setd_INT_HCIs <- stack(apply(post_BR_setd_INT, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_BR_setd_INT_sum <- merge(post_BR_setd_INT_mean, post_BR_setd_INT_LCIs, by = "ind") %>%
														 merge(post_BR_setd_INT_HCIs) 
	post_BR_setd_INT_sum <- round_df(post_BR_setd_INT_sum, 2) %>% round_df(2)

	head(post_BR_setd_INT_sum)
	
	# proportion > 0
	BR_setd_INT_GAI_prob0 <- (length(which(post_BR_setd_INT$b_LogGAI_setd >0)))/length(post_BR_setd_INT$b_LogGAI_setd)*100 
	BR_setd_INT_GAI_prob0 <- round(BR_setd_INT_GAI_prob0, 1) # 99.4

 
	#fit robust regression and compare to gaussian 
 
  #### robust regression with x & y reversed ####
	
	# robust with weak prior on nu
	#brr_weak_setd <- brm(data = HM_dat, family = student,
  #                   P ~ LogGAI_setd, 
  #                   prior = c(prior(student_t(3, 1, 2.5), class = Intercept),
  #                             prior(student_t(3, 0, 10), class = b),
  #                             prior(student_t(3, 0, 2.5),  class = sigma),
  #                             prior(gamma(2, 0.1), class = nu)),
  #                   chains = 4,
  #                   save_all_pars = TRUE,
  #                   iter = 10000, warmup = 2000)
#
	#summary(brr_weak_setd)
	#loo(brr_weak_setd, moment_match = TRUE)
	#
	##save(brr_weak_setd, file = here("output/brr_weak_setd.rds"))
  ##load(file = "output/brr_weak_setd.rds")
#
	## strong prior on nu
	#brr_strong_setd <- brm(data = HM_dat, family = student,
	#                     P ~ LogGAI_setd, 
	#                     prior = c(prior(student_t(3, 1, 2.5), class = Intercept),
	#                               prior(student_t(3, 0, 10), class = b),
	#                               prior(student_t(3, 0, 2.5),  class = sigma),
	#                               prior(gamma(4, 1), class = nu)),
	#                     chains = 4,
	#                     save_all_pars = TRUE,
	#                     iter = 10000, warmup = 2000)
#
	#summary(brr_strong_setd)
	#loo(brr_strong_setd, moment_match = TRUE)
#
  #save(brr_strong_setd, file = here("output/brr_strong_setd.rds"))
  #load(file = "output/brr_strong_setd.rds")

	### all slopes overlapped with 95% BCIs of each other, use simplest model

	
	#### B. GAI and new data #### 

	# estimate gill area index

	raw_GSA_df <- GSA_growth_dat %>% 
							  filter(RawOrMean == "raw") %>%
							  group_by(Binomial) %>%
								mutate(MeanGSAcm2 = mean(GSAcm2),
											 MeanMassG = mean(MassG),
										   MeanLogGSAcm2 = mean(LogGSAcm2),
											 MeanLogMassG = mean(LogMassG))

	unique_raw_GSA_df <- raw_GSA_df%>%
								       distinct(., Binomial, .keep_all = TRUE)
																
	mean_GSA_df <- GSA_growth_dat %>% 
		filter(RawOrMean == "mean") %>%
		distinct(., Binomial, .keep_all = TRUE) %>%
		mutate(MeanGSAcm2 = GSAcm2,
					 MeanMassG = MassG,
					 MeanLogGSAcm2 = LogGSAcm2,
					 MeanLogMassG = LogMassG) 

	`%notin%` <- Negate(`%in%`)
	mean_GSA_df <- subset(mean_GSA_df, Binomial %notin% raw_GSA_df$Binomial)

	all_mean_GSA_df <- rbind(unique_raw_GSA_df, mean_GSA_df)

	all_mean_GSA_df <- all_mean_GSA_df %>%
	              		 rowwise() %>%
	              		 mutate(GAI_setd = GAI_formula(G = MeanGSAcm2, W = MeanMassG),
	               			  		LogGAI_setd = log10(GAI_setd))


	# fit the model
	GP_GSA_lm <- lm(mean_gperf_simple ~ LogGAI_setd, data = all_mean_GSA_df)
	summary(GP_GSA_lm)
	confint(GP_GSA_lm)
	# slope = 0.87 (0.52 to 1.23)
	# intercept = 1.92 (1.53 to 2.32)
	# so, relationship (slope) is significant 

	# outliers? 
	aug_GP_GSA_lm <- augment(GP_GSA_lm)
	aug_GP_GSA_lm %>%
	mutate(i = 1:n()) %>%
	ggplot(aes(x = i, y = .cooksd)) +
	geom_hline(yintercept = .5, color = "red", size = 1) +
	geom_point(alpha = .5, size = 3) +
	ylab("Cook's distance") +
	xlab("Observation")
	# no outliers 

	# run a model and plot
	#brms_GP_GAI <- brm(
	#				  mean_gperf_simple ~ LogGAI_setd, 
	#					data = all_mean_GSA_df, family = gaussian(),
	#					prior = c(
	#              prior(student_t(3, 0, 10), "b"),
	#              prior(student_t(3, 2.7, 2.5), "Intercept"),
	#              prior(cauchy(0, 10), "sigma")
	#              ),  					
	#  			 chains = 4,
	#				 save_pars = save_pars(all = TRUE),
	#  			 iter = 5000, warmup = 1000,
	#  			 control = list(adapt_delta = 0.90,
	#  			                max_treedepth = 18))
#
	#summary(brms_GP_GAI)
	#
	#loo(brms_GP_GAI, moment_match = TRUE)
	#bayes_R2(brms_GP_GAI)
	
	#save(brms_GP_GAI, file = here("output/brms_GP_GAI.rds"))
  load(file = "output/brms_GP_GAI.rds")

	post_brms_GP_GAI_INT <- as.data.frame(brms_GP_GAI)
	post_brms_GP_GAI_INT_mean <- stack(apply(post_brms_GP_GAI_INT, 2, mean)) %>% rename(mean = values)
	post_brms_GP_GAI_INT_LCIs <- stack(apply(post_brms_GP_GAI_INT, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_brms_GP_GAI_INT_HCIs <- stack(apply(post_brms_GP_GAI_INT, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_brms_GP_GAI_INT_sum <- merge(post_brms_GP_GAI_INT_mean, post_brms_GP_GAI_INT_LCIs, by = "ind") %>%
														 merge(post_brms_GP_GAI_INT_HCIs) 
	post_brms_GP_GAI_INT_sum <- round_df(post_brms_GP_GAI_INT_sum, 2) %>% round_df(2)

	head(post_brms_GP_GAI_INT_sum)
	
	# proportion > 0
	brms_GP_GAI_INT_INT_GAI_prob0 <- (length(which(post_brms_GP_GAI_INT$b_LogGAI_setd >0)))/length(post_brms_GP_GAI_INT$b_LogGAI_setd)*100 
	brms_GP_GAI_INT_INT_GAI_prob0 <- round(brms_GP_GAI_INT_INT_GAI_prob0, 1) # 100


	# plot
	
	newdat <- expand.grid(
	  LogGAI_setd = all_mean_GSA_df$LogGAI_setd
	)
	
	fitted_mydat <- fitted(brms_GP_GAI, newdata = newdat,
	                        interval = "confidence")
	
	all_mean_GSA_df_plot <- cbind(all_mean_GSA_df, fitted_mydat) %>%
	                        rename(fitted_lm = Estimate, fit_lwr_lm = Q2.5, 
	                        			 fit_upr_lm = Q97.5)


	### without species traditionally used in aquaculture #### START HERE ######

	# estimate gill area index
	
	GSA_growth_dat_no_ac <- GSA_growth_dat %>%
		filter(Binomial %notin% aquaculture_sp_mean)

	raw_GSA_df_no_ac <- GSA_growth_dat_no_ac %>% 
							        filter(RawOrMean == "raw") %>%
							        group_by(Binomial) %>%
								      mutate(MeanGSAcm2 = mean(GSAcm2),
								      			 MeanMassG = mean(MassG),
								      		   MeanLogGSAcm2 = mean(LogGSAcm2),
								      			 MeanLogMassG = mean(LogMassG))

	unique_raw_GSA_df_no_ac <- raw_GSA_df_no_ac %>%
								      			 distinct(., Binomial, .keep_all = TRUE)

	mean_GSA_df_no_ac <- GSA_growth_dat_no_ac %>% 
											 filter(RawOrMean == "mean") %>%
											 distinct(., Binomial, .keep_all = TRUE) %>%
											 mutate(MeanGSAcm2 = GSAcm2,
											 			 MeanMassG = MassG,
											 			 MeanLogGSAcm2 = LogGSAcm2,
											 			 MeanLogMassG = LogMassG) 
	
	mean_GSA_df_no_ac <- subset(mean_GSA_df_no_ac, Binomial %notin% raw_GSA_df_no_ac$Binomial)

	all_mean_GSA_df_no_ac <- rbind(unique_raw_GSA_df_no_ac, mean_GSA_df_no_ac)

	GAI_formula <- function(G, W){
	  GAI = G / (W^0.8)
	  return(GAI)
	}

	all_mean_GSA_df_no_ac <- all_mean_GSA_df_no_ac %>%
	              		       rowwise() %>%
	              		       mutate(GAI_setd = GAI_formula(G = MeanGSAcm2, W = MeanMassG),
	               		      	  		LogGAI_setd = log10(GAI_setd))
	
	ggplot(all_mean_GSA_df_no_ac, aes(x = LogGAI_setd, y = mean_gperf_simple)) +
		geom_point()

	len(all_mean_GSA_df_no_ac$Binomial)

	# fit the model 
	GP_GSA_lm_ac <- lm(mean_gperf_simple ~ LogGAI_setd, data = all_mean_GSA_df_no_ac)
	summary(GP_GSA_lm_ac)
	confint(GP_GSA_lm_ac) %>% round(2)

	# check outliers
	
	aug_pauly <- augment(GP_GSA_lm_ac)
	aug_pauly %>%
	mutate(i = 1:n()) %>%
	ggplot(aes(x = i, y = .cooksd)) +
	geom_hline(yintercept = .5, color = "red", size = 1) +
	geom_point(alpha = .5, size = 3) +
	ylab("Cook's distance") +
	xlab("Observation") 
	# no outliers


	brms_GP_GSA_ac <- brm(
					  mean_gperf_simple ~ LogGAI_setd, 
						data = all_mean_GSA_df_no_ac, family = gaussian(),
						prior = c(
	              prior(student_t(3, 0, 10), "b"),
	              prior(student_t(3, 0, 10), "Intercept"),
	              prior(cauchy(0, 10), "sigma")
	              ),  					
	  			 chains = 4,
	  			 iter = 5000, warmup = 1000,
	  			 control = list(adapt_delta = 0.90,
	  			                max_treedepth = 18))

	save(brms_GP_GSA_ac, file = here("./output/brms_GP_GAI_NO_AQUACULTURE.rds"))

	summary(brms_GP_GSA_ac)

	# outlier check

 loo_brms_GP_GSA_ac <- loo(brms_GP_GSA_ac)
 

	### without air-breathers
 
	## my data 

	# estimate gill area index
	GSA_growth_dat_no_abr <- GSA_growth_dat %>%
		filter(Binomial %notin% airbreathing_mean)
	
	raw_GSA_df_no_abr <- GSA_growth_dat_no_abr %>% 
							         filter(RawOrMean == "raw") %>%
							         group_by(Binomial) %>%
								       mutate(MeanGSAcm2 = mean(GSAcm2),
								       			  MeanMassG = mean(MassG),
								       		    MeanLogGSAcm2 = mean(LogGSAcm2),
								       			  MeanLogMassG = mean(LogMassG))
	
	unique_raw_GSA_df_no_abr <- raw_GSA_df_no_abr %>%
								      			 distinct(., Binomial, .keep_all = TRUE)
	
	mean_GSA_df_no_abr <- GSA_growth_dat_no_abr %>% 
											 filter(RawOrMean == "mean") %>%
											 distinct(., Binomial, .keep_all = TRUE) %>%
											 mutate(MeanGSAcm2 = GSAcm2,
											 			 MeanMassG = MassG,
											 			 MeanLogGSAcm2 = LogGSAcm2,
											 			 MeanLogMassG = LogMassG) 
	
	mean_GSA_df_no_abr <- subset(mean_GSA_df_no_abr, Binomial %notin% raw_GSA_df_no_abr$Binomial)
	
	all_mean_GSA_df_no_abr <- rbind(unique_raw_GSA_df_no_abr, mean_GSA_df_no_abr)
	
	GAI_formula <- function(G, W){
	  GAI = G / (W^0.8)
	  return(GAI)
	}

	all_mean_GSA_df_no_abr <- all_mean_GSA_df_no_abr %>%
	              		        rowwise() %>%
	              		        mutate(GAI_setd = GAI_formula(G = MeanGSAcm2, W = MeanMassG),
	               		      	  		 LogGAI_setd = log10(GAI_setd))

	ggplot(all_mean_GSA_df_no_abr, aes(x = LogGAI_setd, y = mean_gperf_simple)) +
		geom_point()

	len(all_mean_GSA_df_no_abr$Binomial)

	# fit the model 
	GP_GSA_lm <- lm(mean_gperf_simple ~ LogGAI_setd, data = all_mean_GSA_df_no_abr)
	summary(GP_GSA_lm)
	confint(GP_GSA_lm) %>% round(2)

	# check outliers
	
	aug_pauly <- augment(GP_GSA_lm)
	aug_pauly %>%
	mutate(i = 1:n()) %>%
	ggplot(aes(x = i, y = .cooksd)) +
	geom_hline(yintercept = .5, color = "red", size = 1) +
	geom_point(alpha = .5, size = 3) +
	ylab("Cook's distance") +
	xlab("Observation") 
	# no outliers


	brms_GP_GSA_no_abr <- brm(
					  mean_gperf_simple ~ LogGAI_setd, 
						data = all_mean_GSA_df_no_abr, family = gaussian(),
						prior = c(
	              prior(student_t(3, 0, 10), "b"),
	              prior(student_t(3, 0, 10), "Intercept"),
	              prior(cauchy(0, 10), "sigma")
	              ),  					
	  			 chains = 4,
	  			 iter = 5000, warmup = 1000,
	  			 control = list(adapt_delta = 0.90,
	  			                max_treedepth = 18))
	
	save(brms_GP_GSA_no_abr, file = here("./output/brms_GP_GAI_NO_AIRBREATHERS.rds"))
	
	summary(brms_GP_GSA_no_abr)
