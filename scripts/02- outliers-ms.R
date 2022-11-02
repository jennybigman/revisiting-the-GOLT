#02-outliers

	# Question 1 in MS

	#### regular lm ####
	pauly_fit_lm <- lm(log10(GAI) ~ P, data = Pauly_dat)
	summary(pauly_fit_lm)

	fit_newdat_lm <- expand.grid(
	  P = Pauly_dat$P
	)

	fitted_pauly_dat_lm <- predict(pauly_fit_lm, newdata = fit_newdat_lm,
	                            interval = "confidence")
	
	pred_pauly_dat_lm <- predict(pauly_fit_lm, newdata = fit_newdat_lm,
	                            interval = "prediction")
	
	Pauly_dat <- cbind(Pauly_dat, fitted_pauly_dat_lm) %>%
	             rename(fitted_lm = fit, fit_lwr_lm = lwr, fit_upr_lm = upr)

	Pauly_dat <- cbind(Pauly_dat, pred_pauly_dat_lm) %>%
	             rename(pred_lm = fit, pred_lwr = lwr, pred_upr = upr)

	## measures of leverage on model coefficients ##

	# cook's distance
	aug_pauly <- augment(pauly_fit_lm)
	glimpse(aug_pauly)

  cooksd_plot <-
  	aug_pauly %>%  
  	mutate(i = 1:n()) %>%
  	ggplot(aes(x = i, y = .cooksd)) +
  	geom_hline(yintercept = .5, color = "red", size = 1) +
  	geom_point(alpha = .5, size = 3) +
  	ylab("Cook's distance") +
  	xlab("Observation") +
  	ylim(0, 0.52) +
  	xlim(0, 42) +
  	theme_bw() +
  	theme(
  	  axis.text.y =element_text(size= 8, colour = "black"),
  	  axis.title.y = element_text(size=10),
  	  axis.title.x = element_blank(),
  	  axis.text.x = element_blank(),
  	  axis.ticks.x = element_blank(),
  	  panel.grid.major = element_blank(),
  	  panel.grid.minor = element_blank(),
  	  panel.border = element_rect(linetype = 1, size= 1, fill = NA)) +
  	annotate("text", x = 0, y = 0.52,
             label = "(a)", size = 4, fontface = "bold")

 
  ## Pareto k ##
	#b_pauly_normal <- 
	#  brm(data = Pauly_dat, family = gaussian,
	#      log10(GAI) ~ P, 
	#      prior = c(prior(student_t(3, 1, 2.5), class = Intercept),
	#                prior(student_t(3, 0, 10), class = b),
	#                prior(student_t(3, 0, 2.5),  class = sigma)),
	#      chains = 4,
	#      save_all_pars = TRUE,
	#      iter = 10000, warmup = 2000)
	
	#save(b_pauly_normal, file = here("output/b_pauly_normal.rds"))
	load(file = "output/b_pauly_normal.rds")
  
	loo_Pauly_all_normal <- loo(b_pauly_normal)
	#loo_Pauly_all_normal <- loo(b_pauly_normal, moment_match = TRUE)

	pareto_k_table(loo_Pauly_all_normal) # summary of distribution of pareto-k values
	pareto_k_ids(loo_Pauly_all_normal) # identify which case(s) is bad

	paretok <- loo_Pauly_all_normal$diagnostics$pareto_k %>% as_tibble()

	# plot
	pareto_k_plot <-
	loo_Pauly_all_normal$diagnostics$pareto_k %>%  
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
  	white_theme() +
		annotate("text", x = 0, y = 1.05,
              label = "(b)", size = 4, fontface = "bold")

	
	# plot together
	plot1 <- cooksd_plot + theme(plot.margin = unit(c(0.2, 0.2, 0, 0.2), "in"))
	plot2 <- pareto_k_plot + theme(plot.margin = unit(c(0, 0.2, 0.2, 0.2), "in"))
	
	outlier_plot <- plot1/plot2 +  plot_layout(heights = c(1, 1.1))
	
	ggsave(path = here("output/"), filename = "outlier_plot.png", 
	       plot = outlier_plot,
	       height = 7.5, width = 7.5, units = "in")

	#### quantreg ####

	# robust regression using quantreg package
	pauly_fit_qr <- rq(log10(GAI) ~ P, data = Pauly_dat)
	summary(pauly_fit_qr)

	fit_newdat_qr <- expand.grid(
	  P = Pauly_dat$P
	)

	fitted_pauly_dat_qr <- predict(pauly_fit_qr, newdata = fit_newdat_qr,
	                            interval = "confidence")

	Pauly_dat <- cbind(Pauly_dat, fitted_pauly_dat_qr) %>%
	             rename(fitted_qr = fit, fit_lwr_qr = lower, fit_upr_qr = higher)

	
	#### robust regression using iteratively reqeighting least squares (IRLS) ####
	pauly_fit_IRLS <- rlm(log10(GAI) ~ P, data = Pauly_dat)
	summary(pauly_fit_IRLS)
	
	fit_newdat_IRLS <- expand.grid(
	  P = Pauly_dat$P
	)

	fitted_pauly_dat_IRLS <- predict(pauly_fit_IRLS, newdata = fit_newdat_IRLS,
	                            interval = "confidence")

	Pauly_dat <- cbind(Pauly_dat, fitted_pauly_dat_IRLS) %>%
	             rename(fitted_IRLS = fit, fit_lwr_IRLS = lwr, fit_upr_IRLS = upr)

	# estimate CIs
	IRLS_output <- summary(pauly_fit_IRLS)[[4]]

	int <- IRLS_output[1, 1]
	int_se <- IRLS_output[1, 2]
	int_uprCI <- int + (1.96*int_se)
	int_lwrCI <- int - (1.96*int_se)
	
	slope <- IRLS_output[2, 1]
	slope_se <- IRLS_output[2, 2]
	slope_uprCI <- slope + (1.96*slope_se)
	slope_lwrCI <- slope - (1.96*slope_se)

	# remove the 5 outliers identified visually

	Lchal_x <- Pauly_dat$P[4]
	Lchal_y <- Pauly_dat$GAI[4]
	
	Tlept_x <- Pauly_dat$P[24]
	Tlept_y <- Pauly_dat$GAI[24]
	
	Btyr_x <- Pauly_dat$P[6]
	Btyr_y <- Pauly_dat$GAI[6]
	
	Eenc_x <- Pauly_dat$P[7]
	Eenc_y <- Pauly_dat$GAI[7]
	
	Cgob_x <- Pauly_dat$P[29]
	Cgob_y <- Pauly_dat$GAI[29]


	# create new df
	outliers1 <- data.frame(
	  x = c(Lchal_x, Tlept_x),
	  y = c(Lchal_y, Tlept_y)
	)

	# create new df
	outliers2 <- data.frame(
	  x = c(Btyr_x, Eenc_x, Cgob_x),
	  y = c(Btyr_y, Eenc_y, Cgob_y)
	)

	outliers <- rbind(outliers1, outliers2)

	Pauly_dat_outliers_all <- Pauly_dat %>%
	                      dplyr::select(Species, Winf, K, d, P, GAI) %>%
	                      filter(Species != "Latimeria chalumnae",
	                             Species != "Trichiurus lepturus",
	                             Species != "Brevoortia tyrannus",
	                             Species != "Engraulis encrasicholus",
	                             Species != "Cottus gobio")

	Pauly_dat_outliers_2 <- Pauly_dat %>%
	                      dplyr::select(Species, Winf, K, d, P, GAI) %>%
	                      filter(Species != "Latimeria chalumnae",
	                             Species != "Trichiurus lepturus"
	                             )


	# reduced major axis regression (what pauly did, but called it
	# functional regression)
	pauly_out_MA <- lmodel2(log10(GAI) ~ P, data = Pauly_dat_outliers_all, nperm = 99)

	pauly_out_2_MA <- lmodel2(log10(GAI) ~ P, data = Pauly_dat_outliers_2, nperm = 99)

	pauly_no_out_MA <- lmodel2(log10(GAI) ~ P, data = Pauly_dat, nperm = 99)


	#### bayesian regression ####

	# weak prior or nu

	brr_weak <- brm(data = Pauly_dat, family = student,
	    log10(GAI) ~ P, 
	    prior = c(prior(student_t(3, 1, 2.5), class = Intercept),
	              prior(student_t(3, 0, 10), class = b),
	              prior(student_t(3, 0, 2.5),  class = sigma),
	              prior(gamma(2, 0.1), class = nu)),
	    chains = 4,
	    iter = 10000, warmup = 2000)
	
	# strong gamma prior on nu
	brr_strong <- 
		  update(brr_weak,
		         prior = c(prior(student_t(3, 1, 2.5), class = Intercept),
		                   prior(student_t(3, 0, 10), class = b),
		                   prior(student_t(3, 0, 2.5),  class = sigma),
		                   prior(gamma(4, 1),   class = nu)))
		

 