# 05 - calculating gill area index three different ways (empirically est d, predicted d from Pauly 1981, constant d)

	# Question 4 in MS

	# here, we are working with the raw dataset, but to estimate gill area index, we need to take a mean of 
	# gill surface area and body mass 

	mean_dat <- RawGSA8_phylo %>% 
						  group_by(Binomial, mean_gperf_simple) %>%
						  summarise(meanGSAcm2 = mean(GSAcm2),
						  					meanMassG = mean(MassG),
						  					mean_Winfinity = mean(mean_Winfinity), # need to keep Winf in dataset
						  					Log_Winfinity = log10(mean_Winfinity))
	
	#### calculate gill area index ####
	
	# for the gill area index with an empirically estimated d, the gill area index will be calculated in 
	# the stan model

	#### calculate gill area index using d = 0.8 ####
	
	GAI_setd_formula <- function(G, W){
	  GAI = G / (W^0.8)
	  return(GAI)
	}
	
	mean_dat <- mean_dat %>%
	            rowwise() %>%
	            mutate(GAI_setd = GAI_setd_formula(G = meanGSAcm2, W = meanMassG),
	            			 LogGAI_setd = log10(GAI_setd))


	#### calculate gill area index using d predicted from relationship in Pauly (1981) ####
	
	predd_formula <- function(x){
	  d = 0.6742 + 0.03574 * x
	  return(d)
	}
	
	GAI_pred_formula <- function(G, W, x){
	  GAI = G / (W^x)
	  return(GAI)
	}
	
	mean_dat <- mean_dat %>%
	            rowwise() %>%
	            mutate(pred_d = predd_formula(x = Log_Winfinity),
	            			 GAI_predd = GAI_pred_formula(G = meanGSAcm2, W = meanMassG, x = pred_d),
	            			 LogGAI_predd = log10(GAI_predd))
	
	# standardize all predictors
	mean_setd <- mean(mean_dat$LogGAI_setd)
	sd_setd <- sd(mean_dat$LogGAI_setd)

	mean_predd <- mean(mean_dat$LogGAI_predd)
	sd_predd <- sd(mean_dat$LogGAI_predd)

	mean_dat <- mean_dat %>%
		mutate(LogGAI_setd_std = ((LogGAI_setd - mean_setd) / sd_setd),
					 LogGAI_predd_std = ((LogGAI_predd - mean_predd) / sd_predd))
	
	#### fit models of GP ~ GAI ####
	
	# 1. GP ~ gill area index with empirically estimated d (using stan)
	
	RawGSA8_phylo$id = as.numeric(factor(RawGSA8_phylo$Binomial))

	dat2 <- list(
							N=nrow(RawGSA8_phylo),
							J=length(mean_dat$Binomial),
							sp=RawGSA8_phylo$id,
							mean_GSA=mean_dat$meanGSAcm2,
							mean_mass =mean_dat$meanMassG,
							LogGSAcm2=RawGSA8_phylo$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo$LogCenterMassG,
							GP= mean_dat$Log_Winfinity
							)

	MLM_GAI <- stan(file = here("./stan code/MLM_GAI_CORRECT.stan"),
	              		  data = dat2,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.99,
	              		  							 "max_treedepth" = 18))

	save(MLM_GAI, file = ("./output/MLM_GAI_CORRECT.rds"))
	load(file = here("./output/MLM_GAI_CORRECT.rds"))

	post_MLM_GAI <- as.data.frame(MLM_GAI)
	names_post_MLM_GAI <- names(post_MLM_GAI)
	
	post_MLM_GAI_mean <- stack(apply(post_MLM_GAI, 2, mean)) %>% rename(mean = values)
	
	post_MLM_GAI_LCIs <- stack(apply(post_MLM_GAI, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_MLM_GAI_HCIs <- stack(apply(post_MLM_GAI, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	
	post_MLM_GAI_sum <- merge(post_MLM_GAI_mean, post_MLM_GAI_LCIs, by = "ind") %>%
											merge(post_MLM_GAI_HCIs) %>% round_df(2)
	
	post_MLM_GAI_sd <- stack(apply(post_MLM_GAI, 2, sd)) %>% rename(mean = values)

	# proportion > 0
	GP_GAI_empd_prob0 <- (length(which(post_MLM_GAI$bGP_GAI > 0))) / 
													 length(post_MLM_GAI$bGP_GAI)*100 
	GP_GAI_empd_prob0 <- round(GP_GAI_empd_prob0, 1) 

	
	# 2. with GAI from constant (set) d
	
	#get_prior(mean_gperf_simple ~ LogGAI_setd_std,
	#          data = mean_dat,
	#          family = gaussian())
	# 
	#GP_GAI_setd <-  brm(data = mean_dat, family = gaussian,
	#      						      mean_gperf_simple ~ LogGAI_setd_std,
	#      						      prior = c(prior(student_t(3, 2.9, 2.5), class = Intercept),
	#      						           	    prior(student_t(3, 5, 10), class = b),
	#      						           	    prior(student_t(3, 0, 2.5),  class = sigma)),
	#      						           	chains = 4,
	#      						           	iter = 5000, warmup = 1000)
	#
	#save(GP_GAI_setd, file = here("./output/GP_GAI_setd.rds"))
	load(file = here("./output/GP_GAI_setd.rds"))

	# summary 
	post_GP_GAI_setd <- as.data.frame(GP_GAI_setd)
	post_GP_GAI_setd_mean <- stack(apply(post_GP_GAI_setd, 2, mean)) %>% rename(mean = values)
	post_GP_GAI_setd_LCIs <- stack(apply(post_GP_GAI_setd, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_GP_GAI_setd_HCIs <- stack(apply(post_GP_GAI_setd, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_GP_GAI_setd_sum <- merge(post_GP_GAI_setd_mean, post_GP_GAI_setd_LCIs, by = "ind") %>%
														      merge(post_GP_GAI_setd_HCIs) 
	post_GP_GAI_setd_sum <- round_df(post_GP_GAI_setd_sum, 2) %>% round_df(2)

	head(post_GP_GAI_setd_sum)
	
	# proportion > 0
	GP_GP_GAI_setd_prob0 <- (length(which(post_GP_GAI_setd$b_LogGAI_setd_std >0))) / 
													 length(post_GP_GAI_setd$b_LogGAI_setd_std)*100 
	GP_GP_GAI_setd_prob0 <- round(GP_GP_GAI_setd_prob0, 1) 
	
	# 3. with GAI from predicted d
	
	#get_prior(mean_gperf_simple ~ LogGAI_predd_std,
	#          data = mean_dat,
	#          family = gaussian())
	# 
	#GP_GAI_predd <-  brm(data = mean_dat, family = gaussian,
	#      						      mean_gperf_simple ~ LogGAI_predd_std,
	#      						      prior = c(prior(student_t(3, 2.9, 2.5), class = Intercept),
	#      						           	    prior(student_t(3, 5, 10), class = b),
	#      						           	    prior(student_t(3, 0, 2.5),  class = sigma)),
	#      						           	chains = 4,
	#      						           	iter = 5000, warmup = 1000)
	#
	#save(GP_GAI_predd, file = here("./output/GP_GAI_predd.rds"))
	load(file = here("./output/GP_GAI_predd.rds"))

	# summary 
	post_GP_GAI_predd <- as.data.frame(GP_GAI_predd)
	post_GP_GAI_predd_mean <- stack(apply(post_GP_GAI_predd, 2, mean)) %>% rename(mean = values)
	post_GP_GAI_predd_LCIs <- stack(apply(post_GP_GAI_predd, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_GP_GAI_predd_HCIs <- stack(apply(post_GP_GAI_predd, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_GP_GAI_predd_sum <- merge(post_GP_GAI_predd_mean, post_GP_GAI_predd_LCIs, by = "ind") %>%
														      merge(post_GP_GAI_predd_HCIs) 
	post_GP_GAI_predd_sum <- round_df(post_GP_GAI_predd_sum, 2) %>% round_df(2)

	head(post_GP_GAI_predd_sum)
	
	# proportion > 0
	GP_GP_GAI_predd_prob0 <- (length(which(post_GP_GAI_predd$b_LogGAI_predd_std >0))) / 
														length(post_GP_GAI_predd$b_LogGAI_predd_std)*100 
	GP_GP_GAI_predd_prob0 <- round(GP_GP_GAI_predd_prob0, 1) 
	
	
	#########################################
	#### With phylogeny added ####
	#########################################
	
	# prep phylo
	
	# identity matrix for stan
	d_mat <- diag(1, 32, 32)
	
	# for brms, need to take phylo and turn it into tabular form
	A <- ape::vcv.phylo(tree_pruned, corr = FALSE)
	
	#make object with rownames as species in dataset to match tree
	sp_list <- as_tibble(unique(mean_dat$Binomial)) %>%
	           rename(., Species = value)
	
	# identify if species in dataset are on fishtree
	Binomial <- as.character(sp_list$Species)
	Species <- lapply(strsplit(as.character(Binomial), "\\ "), "[", 2)
	Genus <- lapply(strsplit(as.character(Binomial), "\\ "), "[", 1)
	
	species_list <- tibble(Binomial, Genus, Species)
	
	species_list$Binomial2 <- paste(species_list$Genus, species_list$Species,
                                sep = "_")

	phylo_list <- species_list[match(mean_dat$Binomial,
                                species_list$Binomial),
                               'Binomial2']

	mean_dat$phylo <- phylo_list$Binomial2

	# make sure they match
	setdiff(rownames(A), mean_dat$phylo)

	# 1. with GAI from empirically estimated d (using stan)
	
	RawGSA8_phylo$id = as.numeric(factor(RawGSA8_phylo$Binomial))

	dat2 <- list(
							d_mat = d_mat,
							vcov_mat = A,
							N=nrow(RawGSA8_phylo),
							J=length(mean_dat$Binomial),
							sp=RawGSA8_phylo$id,
							mean_GSA=mean_dat$meanGSAcm2,
							mean_mass =mean_dat$meanMassG,
							LogGSAcm2=RawGSA8_phylo$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo$LogCenterMassG,
							GP= mean_dat$mean_gperf_simple
							)

	MLM_GAI_phylo <- stan(file = here("./stan code/MLM_GAI_phylo.stan"),
	              		  data = dat2,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.99,
	              		  							 "max_treedepth" = 18),
	              		  pars = c("global_slope",
	              		  				 "global_int",
	              		           "beta_int",
	              		  				 "beta_slope",
	              		  				 "beta_slopes",
	              		           "sigma",
	              		  				 "log_lik",
	              		  				 "GAI",
	              		  				 "LogGAI",
	              		  				 "LogGAI_std",
	              		  				 "aGP_GAI",
	              		  				 "bGP_GAI",
	              		  				 "sigma_GP_GAI",
	              		  				 "lambda"))
	

	save(MLM_GAI_phylo, file = ("./output/MLM_GAI_phylo.rds"))
	load(file = here("./output/MLM_GAI_phylo.rds"))

	post_MLM_GAI_phylo <- as.data.frame(MLM_GAI_phylo)
	names_post_MLM_GAI_phylo <- names(post_MLM_GAI_phylo)
	
	post_MLM_GAI_mean_phylo <- stack(apply(post_MLM_GAI_phylo, 2, mean)) %>% rename(mean = values)
	
	post_MLM_GAI_LCIs_phylo <- stack(apply(post_MLM_GAI_phylo, 2, quantile, prob = (0.025))) %>% 
		rename(LCI = values)
	post_MLM_GAI_HCIs_phylo <- stack(apply(post_MLM_GAI_phylo, 2, quantile, prob = (0.975))) %>% 
		rename(HCI = values)
	
	post_MLM_GAI_phylo_sum <- merge(post_MLM_GAI_mean_phylo, 
																	 post_MLM_GAI_LCIs_phylo, by = "ind") %>%
													   merge(post_MLM_GAI_HCIs_phylo) %>% round_df(2)
	
	post_MLM_GAI_sd_phylo <- stack(apply(post_MLM_GAI_phylo, 2, sd)) %>% rename(mean = values)


	# proportion > 0
	GP_GAI_empd_prob0_phylo <- (length(which(post_MLM_GAI_phylo$bGP_GAI > 0))) / 
													    length(post_MLM_GAI_phylo$bGP_GAI) * 100 
	GP_GAI_empd_prob0_phylo <- round(GP_GAI_empd_prob0_phylo, 1) 


	# 2. with GAI from constant (set) d
	
	GP_GAI_setd_phylo <-  brm(data = mean_dat, family = gaussian,
	      						      mean_gperf_simple ~ LogGAI_setd_std + (1|gr(phylo, cov = A)),
												data2 = list(A = A),
	      						      prior = c(prior(student_t(3, 2.9, 2.5), class = Intercept),
	      						           	    prior(student_t(3, 5, 10), class = b),
	      						           	    prior(student_t(3, 0, 2.5),  class = sigma)),
	      						           	chains = 4, cores = 4,
	      						           	iter = 14000, warmup = 4000)
	
	save(GP_GAI_setd_phylo, file = here("./output/GP_GAI_setd_phylo.rds"))
	load(file = here("./output/GP_GAI_setd_phylo.rds"))

	# summary 
	post_GP_GAI_setd_phylo <- as.data.frame(GP_GAI_setd_phylo)
	
	post_GP_GAI_setd_phylo_mean <- stack(apply(post_GP_GAI_setd_phylo, 
																												2, mean)) %>%
																						rename(mean = values)
	
	post_GP_GAI_setd_phylo_LCI <- stack(apply(post_GP_GAI_setd_phylo, 
																												2, quantile, prob = (0.025))) %>% 
																						rename(LCI = values)
	
	post_GP_GAI_setd_phylo_HCI <- stack(apply(post_GP_GAI_setd_phylo, 
																												2, quantile, prob = (0.975))) %>%
																						rename(HCI = values)
	
	post_GP_GAI_setd_phylo_sum <- merge(post_GP_GAI_setd_phylo_mean, 
																		  post_GP_GAI_setd_phylo_LCI, by = "ind") %>%
														    merge(post_GP_GAI_setd_phylo_HCI) 
	
	post_GP_GAI_setd_phylo_sum <- round_df(post_GP_GAI_setd_phylo_sum, 2) %>% 
																					 round_df(2)

	head(post_GP_GAI_setd_phylo_sum)
	
	# proportion > 0
	GP_GP_GAI_setd_prob0_phylo <- (length(which(post_GP_GAI_setd_phylo$b_LogGAI_setd_std >0))) / 
													 length(post_GP_GAI_setd_phylo$b_LogGAI_setd_std)*100 
	
	GP_GP_GAI_setd_prob0_phylo <- round(GP_GP_GAI_setd_prob0_phylo, 1)
	
	# phylo signal
	hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
	(hyp <- hypothesis(GP_GAI_setd_phylo, hyp, class = NULL))
	
	# 3. with GAI from predicted d
	
	GP_GAI_predd_phylo <-  brm(data = mean_dat, family = gaussian,
	      						      mean_gperf_simple ~ LogGAI_predd_std + (1|gr(phylo, cov = A)),
												data2 = list(A = A),
	      						      prior = c(prior(student_t(3, 2.9, 2.5), class = Intercept),
	      						           	    prior(student_t(3, 5, 10), class = b),
	      						           	    prior(student_t(3, 0, 2.5),  class = sigma)),
	      						           	chains = 4,
	      						           	iter = 20000, warmup = 4000)
	
	save(GP_GAI_predd_phylo, file = here("./output/GP_GAI_predd_phylo.rds"))
	load(file = here("./output/GP_GAI_predd_phylo.rds"))

	# summary 
	post_GP_GAI_predd_phylo <- as.data.frame(GP_GAI_predd_phylo)
	
	post_GP_GAI_predd_phylo_mean <- stack(apply(post_GP_GAI_predd_phylo, 2, mean)) %>% 
																	rename(mean = values)
	
	post_GP_GAI_predd_phylo_LCIs <- stack(apply(post_GP_GAI_predd_phylo, 2, quantile, prob = (0.025))) %>% 
																	rename(LCI = values)
	
	post_GP_GAI_predd_phylo_HCIs <- stack(apply(post_GP_GAI_predd_phylo, 2, quantile, prob = (0.975))) %>% 
																	rename(HCI = values)
	
	post_GP_GAI_predd_phylo_sum <- merge(post_GP_GAI_predd_phylo_mean,
																			 post_GP_GAI_predd_phylo_LCIs, by = "ind") %>%
														     merge(post_GP_GAI_predd_phylo_HCIs) 
	
	post_GP_GAI_predd_phylo_sum <- round_df(post_GP_GAI_predd_phylo_sum, 2) %>% round_df(2)

	head(post_GP_GAI_predd_phylo_sum)
	# proportion > 0
	GP_GP_GAI_predd_phylo_prob0 <- (length(which(post_GP_GAI_predd_phylo$b_LogGAI_predd_std > 0))) / 
														length(post_GP_GAI_predd_phylo$b_LogGAI_predd_std) * 100 
	
	GP_GP_GAI_predd_phylo_prob0 <- round(GP_GP_GAI_predd_phylo_prob0, 1) 
	
	hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
	(hyp <- hypothesis(GP_GAI_predd_phylo, hyp, class = NULL))

	#########################################
	#### Without aquaculture species ####
	#########################################
	
	# trim dataset to exclude aquacultured species
	
	aq_sp <- c("Cirrhinus mrigala",
						 "Oreochromis niloticus")
		
	RawGSA8_phylo_noaq <- RawGSA8_phylo %>%
		filter(Binomial %notin% aq_sp)
	
	mean_dat_noaq <- mean_dat %>%
		filter(Binomial %notin% aq_sp)

	#### fit models of GP ~ GAI ####
	
	# 1. GP ~ gill area index with empirically estimated d (using stan)
	
	RawGSA8_phylo_noaq$id = as.numeric(factor(RawGSA8_phylo_noaq$Binomial))

	dat_noaq <- list(
							N=nrow(RawGSA8_phylo_noaq),
							J=length(mean_dat_noaq$Binomial),
							sp=RawGSA8_phylo_noaq$id,
							mean_GSA=mean_dat_noaq$meanGSAcm2,
							mean_mass =mean_dat_noaq$meanMassG,
							LogGSAcm2=RawGSA8_phylo_noaq$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo_noaq$LogCenterMassG,
							GP= mean_dat_noaq$Log_Winfinity
							)

	MLM_GAI_noaq <- stan(file = here("./stan code/MLM_GAI_CORRECT.stan"),
	              		  data = dat_noaq,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.99,
	              		  							 "max_treedepth" = 18))

	save(MLM_GAI_noaq, file = ("./output/MLM_GAI_CORRECT_noaq.rds"))
	load(file = here("./output/MLM_GAI_CORRECT_noaq.rds"))

	post_MLM_GAI_noaq <- as.data.frame(MLM_GAI_noaq)
	names_post_MLM_GAI_noaq <- names(post_MLM_GAI_noaq)
	
	post_MLM_GAI_mean_noaq <- stack(apply(post_MLM_GAI_noaq, 2, mean)) %>% rename(mean = values)
	
	post_MLM_GAI_LCIs_noaq <- stack(apply(post_MLM_GAI_noaq, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_MLM_GAI_HCIs_noaq <- stack(apply(post_MLM_GAI_noaq, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	
	post_MLM_GAI_sum_noaq <- merge(post_MLM_GAI_mean_noaq, post_MLM_GAI_LCIs_noaq, by = "ind") %>%
											merge(post_MLM_GAI_HCIs_noaq) %>% round_df(2)
	

	# proportion > 0
	GP_GAI_empd_prob0_noaq <- (length(which(post_MLM_GAI_noaq$bGP_GAI > 0))) / 
													 length(post_MLM_GAI_noaq$bGP_GAI)*100 
	GP_GAI_empd_prob0_noaq <- round(GP_GAI_empd_prob0_noaq, 1) 

	
	#2. with GAI from constant (set) d
	
	GP_GAI_setd_noaq <-  brm(data = mean_dat_noaq, family = gaussian,
	      						      mean_gperf_simple ~ LogGAI_setd_std,
	      						      prior = c(prior(student_t(3, 2.9, 2.5), class = Intercept),
	      						           	    prior(student_t(3, 5, 10), class = b),
	      						           	    prior(student_t(3, 0, 2.5),  class = sigma)),
	      						           	chains = 4,
	      						           	iter = 5000, warmup = 1000)
	
	save(GP_GAI_setd_noaq, file = here("./output/GP_GAI_setd_noaq.rds"))
	load(file = here("./output/GP_GAI_setd_noaq.rds"))

	# summary 
	post_GP_GAI_setd_noaq <- as.data.frame(GP_GAI_setd_noaq)
	post_GP_GAI_setd_mean_noaq <- stack(apply(post_GP_GAI_setd_noaq, 2, mean)) %>% rename(mean = values)
	post_GP_GAI_setd_LCIs_noaq <- stack(apply(post_GP_GAI_setd_noaq, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_GP_GAI_setd_HCIs_noaq <- stack(apply(post_GP_GAI_setd_noaq, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_GP_GAI_setd_sum_noaq <- merge(post_GP_GAI_setd_mean_noaq, post_GP_GAI_setd_LCIs_noaq, by = "ind") %>%
														      merge(post_GP_GAI_setd_HCIs_noaq) 
	post_GP_GAI_setd_sum_noaq <- round_df(post_GP_GAI_setd_sum_noaq, 2) %>% round_df(2)

	head(post_GP_GAI_setd_sum_noaq)
	
	# proportion > 0
	GP_GP_GAI_setd_prob0_noaq <- (length(which(post_GP_GAI_setd_noaq$b_LogGAI_setd_std >0))) / 
													 length(post_GP_GAI_setd_noaq$b_LogGAI_setd_std)*100 
	GP_GP_GAI_setd_prob0_noaq <- round(GP_GP_GAI_setd_prob0_noaq, 1) 
	
	# 3. with GAI from predicted d
	
	GP_GAI_predd_noaq <-  brm(data = mean_dat_noaq, family = gaussian,
	      						      mean_gperf_simple ~ LogGAI_predd_std,
	      						      prior = c(prior(student_t(3, 2.9, 2.5), class = Intercept),
	      						           	    prior(student_t(3, 5, 10), class = b),
	      						           	    prior(student_t(3, 0, 2.5),  class = sigma)),
	      						           	chains = 4,
	      						           	iter = 5000, warmup = 1000)
	
	save(GP_GAI_predd_noaq, file = here("./output/GP_GAI_predd_noaq.rds"))
	load(file = here("./output/GP_GAI_predd_noaq.rds"))

	# summary 
	post_GP_GAI_predd_noaq <- as.data.frame(GP_GAI_predd_noaq)
	post_GP_GAI_predd_mean_noaq <- stack(apply(post_GP_GAI_predd_noaq, 2, mean)) %>% rename(mean = values)
	post_GP_GAI_predd_LCIs_noaq <- stack(apply(post_GP_GAI_predd_noaq, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_GP_GAI_predd_HCIs_noaq <- stack(apply(post_GP_GAI_predd_noaq, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_GP_GAI_predd_sum_noaq <- merge(post_GP_GAI_predd_mean_noaq, post_GP_GAI_predd_LCIs_noaq, by = "ind") %>%
														      merge(post_GP_GAI_predd_HCIs_noaq) 
	post_GP_GAI_predd_sum_noaq <- round_df(post_GP_GAI_predd_sum_noaq, 2) %>% round_df(2)

	head(post_GP_GAI_predd_sum_noaq)
	
	# proportion > 0
	GP_GP_GAI_predd_prob0_noaq <- (length(which(post_GP_GAI_predd_noaq$b_LogGAI_predd_std >0))) / 
														length(post_GP_GAI_predd_noaq$b_LogGAI_predd_std)*100 
	GP_GP_GAI_predd_prob0_noaq <- round(GP_GP_GAI_predd_prob0_noaq, 1) 
	
	
	#########################################
	#### Without airbreathers ####
	#########################################
	
	# trim dataset to exclude airbreathing species
	
	ab_sp <- c("Cirrhinus mrigala",
						 "Oreochromis niloticus")
		
	RawGSA8_phylo_noab <- RawGSA8_phylo %>%
		filter(Binomial %notin% ab_sp)
	
	mean_dat_noab <- mean_dat %>%
		filter(Binomial %notin% ab_sp)

	#### fit models of GP ~ GAI ####
	
	# 1. GP ~ gill area index with empirically estimated d (using stan)
	
	RawGSA8_phylo_noab$id = as.numeric(factor(RawGSA8_phylo_noab$Binomial))

	dat_noaq <- list(
							N=nrow(RawGSA8_phylo_noab),
							J=length(mean_dat_noab$Binomial),
							sp=RawGSA8_phylo_noab$id,
							mean_GSA=mean_dat_noab$meanGSAcm2,
							mean_mass =mean_dat_noab$meanMassG,
							LogGSAcm2=RawGSA8_phylo_noab$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo_noab$LogCenterMassG,
							GP= mean_dat_noab$Log_Winfinity
							)

	MLM_GAI_noab <- stan(file = here("./stan code/MLM_GAI_CORRECT.stan"),
	              		  data = dat_noaq,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.99,
	              		  							 "max_treedepth" = 18))

	save(MLM_GAI_noab, file = ("./output/MLM_GAI_CORRECT_noab.rds"))
	load(file = here("./output/MLM_GAI_CORRECT_noab.rds"))

	post_MLM_GAI_noab <- as.data.frame(MLM_GAI_noab)
	names_post_MLM_GAI_noab <- names(post_MLM_GAI_noab)
	
	post_MLM_GAI_mean_noab <- stack(apply(post_MLM_GAI_noab, 2, mean)) %>% rename(mean = values)
	
	post_MLM_GAI_LCIs_noab <- stack(apply(post_MLM_GAI_noab, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_MLM_GAI_HCIs_noab <- stack(apply(post_MLM_GAI_noab, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	
	post_MLM_GAI_sum_noab <- merge(post_MLM_GAI_mean_noab, post_MLM_GAI_LCIs_noab, by = "ind") %>%
											merge(post_MLM_GAI_HCIs_noab) %>% round_df(2)
	

	# proportion > 0
	GP_GAI_empd_prob0_noab <- (length(which(post_MLM_GAI_noab$bGP_GAI > 0))) / 
													 length(post_MLM_GAI_noab$bGP_GAI)*100 
	GP_GAI_empd_prob0_noab <- round(GP_GAI_empd_prob0_noab, 1) 

	
	#2. with GAI from constant (set) d
	
	GP_GAI_setd_noab <-  brm(data = mean_dat_noab, family = gaussian,
	      						      mean_gperf_simple ~ LogGAI_setd_std,
	      						      prior = c(prior(student_t(3, 2.9, 2.5), class = Intercept),
	      						           	    prior(student_t(3, 5, 10), class = b),
	      						           	    prior(student_t(3, 0, 2.5),  class = sigma)),
	      						           	chains = 4,
	      						           	iter = 5000, warmup = 1000)
	
	save(GP_GAI_setd_noab, file = here("./output/GP_GAI_setd_noab.rds"))
	load(file = here("./output/GP_GAI_setd_noab.rds"))

	# summary 
	post_GP_GAI_setd_noab <- as.data.frame(GP_GAI_setd_noab)
	post_GP_GAI_setd_mean_noab <- stack(apply(post_GP_GAI_setd_noab, 2, mean)) %>% rename(mean = values)
	post_GP_GAI_setd_LCIs_noab <- stack(apply(post_GP_GAI_setd_noab, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_GP_GAI_setd_HCIs_noab <- stack(apply(post_GP_GAI_setd_noab, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_GP_GAI_setd_sum_noab <- merge(post_GP_GAI_setd_mean_noab, post_GP_GAI_setd_LCIs_noab, by = "ind") %>%
														      merge(post_GP_GAI_setd_HCIs_noab) 
	post_GP_GAI_setd_sum_noab <- round_df(post_GP_GAI_setd_sum_noab, 2) %>% round_df(2)

	head(post_GP_GAI_setd_sum_noab)
	
	# proportion > 0
	GP_GP_GAI_setd_prob0_noab <- (length(which(post_GP_GAI_setd_noab$b_LogGAI_setd_std >0))) / 
													 length(post_GP_GAI_setd_noab$b_LogGAI_setd_std)*100 
	GP_GP_GAI_setd_prob0_noab <- round(GP_GP_GAI_setd_prob0_noab, 1) 
	
	# 3. with GAI from predicted d
	
	GP_GAI_predd_noab <-  brm(data = mean_dat_noab, family = gaussian,
	      						      mean_gperf_simple ~ LogGAI_predd_std,
	      						      prior = c(prior(student_t(3, 2.9, 2.5), class = Intercept),
	      						           	    prior(student_t(3, 5, 10), class = b),
	      						           	    prior(student_t(3, 0, 2.5),  class = sigma)),
	      						           	chains = 4,
	      						           	iter = 5000, warmup = 1000)
	
	save(GP_GAI_predd_noab, file = here("./output/GP_GAI_predd_noab.rds"))
	load(file = here("./output/GP_GAI_predd_noab.rds"))

	# summary 
	post_GP_GAI_predd_noab <- as.data.frame(GP_GAI_predd_noab)
	post_GP_GAI_predd_mean_noab <- stack(apply(post_GP_GAI_predd_noab, 2, mean)) %>% rename(mean = values)
	post_GP_GAI_predd_LCIs_noab <- stack(apply(post_GP_GAI_predd_noab, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_GP_GAI_predd_HCIs_noab <- stack(apply(post_GP_GAI_predd_noab, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_GP_GAI_predd_sum_noab <- merge(post_GP_GAI_predd_mean_noab, post_GP_GAI_predd_LCIs_noab, by = "ind") %>%
														      merge(post_GP_GAI_predd_HCIs_noab) 
	post_GP_GAI_predd_sum_noab <- round_df(post_GP_GAI_predd_sum_noab, 2) %>% round_df(2)

	head(post_GP_GAI_predd_sum_noab)
	
	# proportion > 0
	GP_GP_GAI_predd_prob0_noab <- (length(which(post_GP_GAI_predd_noab$b_LogGAI_predd_std >0))) / 
														length(post_GP_GAI_predd_noab$b_LogGAI_predd_std)*100 
	GP_GP_GAI_predd_prob0_noab <- round(GP_GP_GAI_predd_prob0_noab, 1) 
	
	###########################
	#### Species with over a order of magnitude in body mass ####
	###########################
	
	sp_bm <- GSA_growth_dat %>%
						   filter(., RawOrMean == "raw") %>%
               group_by(Binomial) %>%
            	 summarise(LogMinMass = min(LogMassG),
            	           LogMaxMass = max(LogMassG)) %>% 
            	 mutate(LogMassDiff = LogMaxMass - LogMinMass) %>%
            	 filter(LogMassDiff >= 1)
	
	sp_bm_list <- unique(sp_bm$Binomial)
               
	Raw_GSA_bm <- GSA_growth_dat %>%
 		filter(Binomial %in% sp_bm_list)
 
	mean_dat_bm <- Raw_GSA_bm %>% 
					  group_by(Binomial, mean_gperf_simple) %>%
					  summarise(meanGSAcm2 = mean(GSAcm2),
					  					meanMassG = mean(MassG),
					  					mean_Winfinity = mean(mean_Winfinity), # need to keep Winf in dataset
					  					Log_Winfinity = log10(mean_Winfinity))
 
 # since new mean_dat df, need to calculate GAI
 
 mean_dat_bm <- mean_dat_bm %>%
	            	rowwise() %>%
	            	mutate(GAI_setd = GAI_setd_formula(G = meanGSAcm2, W = meanMassG),
	            				 LogGAI_setd = log10(GAI_setd))
	
 mean_dat_bm <- mean_dat_bm %>%
	            	rowwise() %>%
	            	mutate(pred_d = predd_formula(x = Log_Winfinity),
	            				 GAI_predd = GAI_pred_formula(G = meanGSAcm2, W = meanMassG, x = pred_d),
	            				 LogGAI_predd = log10(GAI_predd))

	# standardize all predictors
	mean_setd <- mean(mean_dat_bm$LogGAI_setd)
	sd_setd <- sd(mean_dat_bm$LogGAI_setd)

	mean_predd <- mean(mean_dat_bm$LogGAI_predd)
	sd_predd <- sd(mean_dat_bm$LogGAI_predd)

	mean_dat_bm <- mean_dat_bm %>%
		mutate(LogGAI_setd_std = ((LogGAI_setd - mean_setd) / sd_setd),
					 LogGAI_predd_std = ((LogGAI_predd - mean_predd) / sd_predd))
	

 #### fit models of GP ~ GAI ####
	
	# 1. GP ~ gill area index with empirically estimated d (using stan)
	
	Raw_GSA_bm$id = as.numeric(factor(Raw_GSA_bm$Binomial))

	dat_bm <- list(
							N=nrow(Raw_GSA_bm),
							J=length(mean_dat_bm$Binomial),
							sp=Raw_GSA_bm$id,
							mean_GSA=mean_dat_bm$meanGSAcm2,
							mean_mass =mean_dat_bm$meanMassG,
							LogGSAcm2=Raw_GSA_bm$LogGSAcm2,
							LogCenterMassG = Raw_GSA_bm$LogCenterMassG,
							GP= mean_dat_bm$Log_Winfinity
							)

	MLM_GAI_bm <- stan(file = here("./stan code/MLM_GAI_CORRECT.stan"),
	              		  data = dat_bm,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.99,
	              		  							 "max_treedepth" = 18))

	save(MLM_GAI_bm, file = ("./output/MLM_GAI_CORRECT_bm.rds"))
	load(file = here("./output/MLM_GAI_CORRECT_bm.rds"))

	post_MLM_GAI_bm <- as.data.frame(MLM_GAI_bm)
	names_post_MLM_GAI_bm <- names(post_MLM_GAI_bm)
	
	post_MLM_GAI_bm_mean <- stack(apply(post_MLM_GAI_bm, 2, mean)) %>% rename(mean = values)
	
	post_MLM_GAI_bm_LCI <- stack(apply(post_MLM_GAI_bm, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_MLM_GAI_bm_HCI <- stack(apply(post_MLM_GAI_bm, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	
	post_MLM_GAI_BM_sum <- merge(post_MLM_GAI_bm_mean, post_MLM_GAI_bm_LCI, by = "ind") %>%
											   merge(post_MLM_GAI_bm_HCI) %>% round_df(2)
	
	# proportion > 0
	MLM_GAI_bm_prob0 <- (length(which(post_MLM_GAI_bm$bGP_GAI > 0))) / 
													 length(post_MLM_GAI_bm$bGP_GAI)*100 
	MLM_GAI_bm_prob0 <- round(MLM_GAI_bm_prob0, 1) 
	
	#2. with GAI from constant (set) d
	
	GP_GAI_setd_bm <-  brm(data = mean_dat_bm, family = gaussian,
	      						      mean_gperf_simple ~ LogGAI_setd_std,
	      						      prior = c(prior(student_t(3, 2.9, 2.5), class = Intercept),
	      						           	    prior(student_t(3, 5, 10), class = b),
	      						           	    prior(student_t(3, 0, 2.5),  class = sigma)),
	      						           	chains = 4,
	      						           	iter = 5000, warmup = 1000)
	
	save(GP_GAI_setd_bm, file = here("./output/GP_GAI_setd_bm.rds"))
	load(file = here("./output/GP_GAI_setd_bm.rds"))

	# summary 
	post_GP_GAI_setd_bm <- as.data.frame(GP_GAI_setd_bm)
	post_GP_GAI_setd_bm_mean <- stack(apply(post_GP_GAI_setd_bm, 2, mean)) %>% rename(mean = values)
	post_GP_GAI_setd_bm_LCI <- stack(apply(post_GP_GAI_setd_bm, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_GP_GAI_setd_bm_HCI <- stack(apply(post_GP_GAI_setd_bm, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_GP_GAI_setd_bm_sum <- merge(post_GP_GAI_setd_bm_mean, post_GP_GAI_setd_bm_LCI, by = "ind") %>%
														      merge(post_GP_GAI_setd_bm_HCI) 
	post_GP_GAI_setd_bm_sum <- round_df(post_GP_GAI_setd_bm_sum, 2) %>% round_df(2)

	head(post_GP_GAI_setd_bm_sum)
	
	# proportion > 0
	GP_GAI_setd_bm_prob0 <- (length(which(post_GP_GAI_setd_bm$b_LogGAI_setd_std >0))) / 
													 length(post_GP_GAI_setd_bm$b_LogGAI_setd_std)*100 
	GP_GAI_setd_bm_prob0 <- round(GP_GAI_setd_bm_prob0, 1) 
	
	# 3. with GAI from predicted d
	
	GP_GAI_predd_bm <-  brm(data = mean_dat_bm, family = gaussian,
	      						      mean_gperf_simple ~ LogGAI_predd_std,
	      						      prior = c(prior(student_t(3, 2.9, 2.5), class = Intercept),
	      						           	    prior(student_t(3, 5, 10), class = b),
	      						           	    prior(student_t(3, 0, 2.5),  class = sigma)),
	      						           	chains = 4,
	      						           	iter = 5000, warmup = 1000)
	
	save(GP_GAI_predd_bm, file = here("./output/GP_GAI_predd_bm.rds"))
	load(file = here("./output/GP_GAI_predd_bm.rds"))

	# summary 
	post_GP_GAI_predd_bm <- as.data.frame(GP_GAI_predd_bm)
	post_GP_GAI_predd_bm_mean <- stack(apply(post_GP_GAI_predd_bm, 2, mean)) %>% rename(mean = values)
	post_GP_GAI_predd_bm_LCI <- stack(apply(post_GP_GAI_predd_bm, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_GP_GAI_predd_bm_HCI <- stack(apply(post_GP_GAI_predd_bm, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_GP_GAI_predd_bm_sum <- merge(post_GP_GAI_predd_bm_mean, post_GP_GAI_predd_bm_LCI, by = "ind") %>%
														  merge(post_GP_GAI_predd_bm_HCI) 
	post_GP_GAI_predd_bm_sum <- round_df(post_GP_GAI_predd_bm_sum, 2) %>% round_df(2)

	head(post_GP_GAI_predd_bm_sum)
	
	# proportion > 0
	GP_GAI_predd_bm_prob0 <- (length(which(post_GP_GAI_predd_bm$b_LogGAI_predd_std >0))) / 
														length(post_GP_GAI_predd_bm$b_LogGAI_predd_std)*100 
	GP_GAI_predd_bm_prob0 <- round(GP_GAI_predd_bm_prob0, 1) 

