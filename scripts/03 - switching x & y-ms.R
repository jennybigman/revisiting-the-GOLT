# 03-switching x & y

	# Question 2 in MS

	# leverage with switching x & y
	pauly_fit_lm_rev <- lm(P ~ LogGAI, data = Pauly_dat)
	summary(pauly_fit_lm_rev)
 
	# cook's distance
	aug_pauly_rev <- augment(pauly_fit_lm_rev)
	glimpse(aug_pauly_rev)
 
	aug_pauly_rev %>%  
	mutate(i = 1:n()) %>%
	ggplot(aes(x = i, y = .cooksd)) +
	geom_hline(yintercept = .5, color = "red", size = 1) +
	geom_point(alpha = .5, size = 3) +
	ylab("Cook's distance") +
	xlab("Observation") +
	ylim(0, 0.6) +
	xlim(0, 42) +
	theme_bw() +
	theme(
	  axis.text.y =element_text(size= 12, colour = "black"),
	  axis.title.y = element_text(size=12),
	  axis.title.x = element_blank(),
	  axis.text.x = element_blank(),
	  axis.ticks.x = element_blank(),
	  panel.grid.major = element_blank(),
	  panel.grid.minor = element_blank(),
	  panel.border = element_rect(linetype = 1, size= 1, fill = NA))
 
	row_num <- which(aug_pauly_rev$.cooksd > 0.5)
	outlierP <- aug_pauly_rev$P[row_num]
	outlierGSA <- aug_pauly_rev$LogGAI[row_num]
	row_num2 <- which((Pauly_dat$P == outlierP) & (Pauly_dat$LogGAI == outlierGSA))
	Pauly_dat$Species[row_num2] 
 
	# Pareto K

	 #b_pauly_normal_rev <- 
	 #  brm(data = Pauly_dat, family = gaussian,
	 #      P ~ LogGAI, 
	 #      prior = c(prior(student_t(3, 2.5, 2.5), class = Intercept),
	 #                prior(student_t(3, 0, 10), class = b),
	 #                prior(student_t(3, 0, 2.5),  class = sigma)),
	 #      chains = 4,
	 #      save_all_pars = TRUE,
	 #      iter = 10000, warmup = 2000)
	 
	#summary(b_pauly_normal_rev)
	#save(b_pauly_normal_rev, file = here("output/b_pauly_normal_rev.rds"))
	load(file = "output/b_pauly_normal_rev.rds")
 
	loo_Pauly_all_normal_rev <- loo(b_pauly_normal_rev, moment_match = TRUE)
	
	
	post_pauly_normal_rev_INT <- as.data.frame(b_pauly_normal_rev)
	post_pauly_normal_rev_INT_mean <- stack(apply(post_pauly_normal_rev_INT, 2, mean)) %>% rename(mean = values)
	post_pauly_normal_rev_INT_LCIs <- stack(apply(post_pauly_normal_rev_INT, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_pauly_normal_rev_INT_HCIs <- stack(apply(post_pauly_normal_rev_INT, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_pauly_normal_rev_INT_sum <- merge(post_pauly_normal_rev_INT_mean, 
																				 post_pauly_normal_rev_INT_LCIs, by = "ind") %>%
																	 merge(post_pauly_normal_rev_INT_HCIs) 
	post_pauly_normal_rev_INT_sum <- round_df(post_pauly_normal_rev_INT_sum, 2) %>% round_df(2)

	head(post_pauly_normal_rev_INT_sum)
	
	# proportion > 0
	pauly_normal_rev_prob0 <- (length(which(post_pauly_normal_rev_INT$b_LogGAI >0)))/length(post_pauly_normal_rev_INT$b_LogGAI)*100 
	pauly_normal_rev_prob0 <- round(pauly_normal_rev_prob0, 1) # 89.9

	
	pareto_k_table(loo_Pauly_all_normal_rev) # summary of distribution of pareto-k values
	pareto_k_ids(loo_Pauly_all_normal_rev) # identify which case(s) is bad
	
	paretok_rev <- loo_Pauly_all_normal_rev$diagnostics$pareto_k %>% as_tibble()
 
	# plot
	 paretok_plot_rev <- 
	 	loo_Pauly_all_normal_rev$diagnostics$pareto_k %>%  
  	as_tibble() %>%
  	mutate(i = 1:n()) %>%
  	rename(pareto_k = value) %>%
  	ggplot(aes(x = i, y = pareto_k)) +
  	geom_hline(yintercept = c(.7, 1), color = "red", size =1) +
  	geom_point(alpha = .5, size = 3) +
  	ylim(-0.074, 1.05) +
  	xlim(0, 42) +
  	ylab("Pareto k value") +
  	xlab("Observation") +
  	theme_bw()

	#### robust regression with x & y reversed ####

	# weak prior on nu
	#brr_weak_rev <- brm(data = Pauly_dat, family = student,
	#                    P ~ LogGAI, 
	#                    prior = c(prior(student_t(3, 1, 2.5), class = Intercept),
	#                              prior(student_t(3, 0, 10), class = b),
	#                              prior(student_t(3, 0, 2.5),  class = sigma),
	#                              prior(gamma(2, 0.1), class = nu)),
	#                    chains = 4,
	#                    save_all_pars = TRUE,
	#                    iter = 10000, warmup = 2000)
	#
	#summary(brr_weak_rev)
	#loo(brr_weak_rev, moment_match = TRUE)
	
	#save(brr_weak_rev, file = here("output/brr_weak_rev.rds"))
	load(file = "output/brr_weak_rev.rds")
	 
	post_brr_weak_INT <- as.data.frame(brr_weak_rev)
	post_brr_weak_INT_mean <- stack(apply(post_brr_weak_INT, 2, mean)) %>% rename(mean = values)
	post_brr_weak_INT_LCIs <- stack(apply(post_brr_weak_INT, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_brr_weak_INT_HCIs <- stack(apply(post_brr_weak_INT, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_brr_weak_INT_sum <- merge(post_brr_weak_INT_mean, post_brr_weak_INT_LCIs, by = "ind") %>%
												 merge(post_brr_weak_INT_HCIs) 
	post_brr_weak_INT_sum <- round_df(post_brr_weak_INT_sum, 2) %>% round_df(2)

	head(post_brr_weak_INT_sum)
	
	# proportion > 0
	brr_weak_GAI_prob0 <- (length(which(post_brr_weak_INT$b_LogGAI >0)))/length(post_brr_weak_INT$b_LogGAI)*100 
	brr_weak_GAI_prob0 <- round(brr_weak_GAI_prob0, 1) # 91.3

 
	# strong prior on nu
	#brr_strong_rev <- brm(data = Pauly_dat, family = student,
	#                     P ~ LogGAI, 
	#                     prior = c(prior(student_t(3, 1, 2.5), class = Intercept),
	#                               prior(student_t(3, 0, 10), class = b),
	#                               prior(student_t(3, 0, 2.5),  class = sigma),
	#                               prior(gamma(4, 1), class = nu)),
	#                     chains = 4,
	#                     save_all_pars = TRUE,
	#                     iter = 10000, warmup = 2000)
	#
	#summary(brr_strong_rev)
	#loo(brr_strong_rev, moment_match = TRUE)
 
	#save(brr_strong_rev, file = here("output/brr_strong_rev.rds"))
	load(file = "output/brr_strong_rev.rds")
 
	post_brr_strong_INT <- as.data.frame(brr_strong_rev)
	post_brr_strong_INT_mean <- stack(apply(post_brr_strong_INT, 2, mean)) %>% rename(mean = values)
	post_brr_strong_INT_LCIs <- stack(apply(post_brr_strong_INT, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_brr_strong_INT_HCIs <- stack(apply(post_brr_strong_INT, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_brr_strong_INT_sum <- merge(post_brr_strong_INT_mean, post_brr_strong_INT_LCIs, by = "ind") %>%
														 merge(post_brr_strong_INT_HCIs) 
	post_brr_strong_INT_sum <- round_df(post_brr_strong_INT_sum, 2) %>% round_df(2)

	head(post_brr_strong_INT_sum)
	
	# proportion > 0
	brr_strong_GAI_prob0 <- (length(which(post_brr_strong_INT$b_LogGAI >0)))/length(post_brr_strong_INT$b_LogGAI)*100 
	brr_strong_GAI_prob0 <- round(brr_strong_GAI_prob0, 1) # 92.9

	#### R squared ####
	
	summary(b_pauly_normal_rev)
  bayes_R2(b_pauly_normal_rev)
 
  summary(brr_weak_rev)
  bayes_R2(brr_weak_rev)

	summary(brr_strong_rev)
	bayes_R2(brr_strong_rev)
	
 
 
