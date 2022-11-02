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
	
	#### calculate gill area index using empirically estimated d ####
	#get_prior(LogGSAcm2 ~ LogMassG,
	#					data = RawGSA8_phylo,
	#					family = gaussian())
	#
	#GSA_mass_fit <- brm(LogGSAcm2 ~ LogMassG * Binomial,
	#										family = gaussian(),
	#										data = RawGSA8_phylo,
	#										prior = c(
	#											prior(student_t(3, 3, 3), class = Intercept),
	#											prior(student_t(3, 5, 10), class = b),
	#											prior(student_t(3, 0, 2.5), class = sigma)),
	#										control = list(max_treedepth = 15),
	#										chains = 4, cores = 4, iter = 5000, warmup = 1000)
	#
	#save(GSA_mass_fit, file = here("./output/GSA_mass_fit.rds"))
	load(file = here("./output/GSA_mass_fit.rds"))
	
	post_GSA_mass <- as.data.frame(GSA_mass_fit)
	
	intercept <- post_GSA_mass[]
	post_GSA_mass_mean <- stack(apply(post_GSA_mass, 2, mean)) %>% rename(mean = values)
	
	
	#### calculate gill area index using d = 0.8 ####
	
	GAI_formula <- function(G, W){
	  GAI = G / (W^0.8)
	  return(GAI)
	}
	
	mean_dat <- mean_dat %>%
	            rowwise() %>%
	            mutate(GAI_setd = GAI_formula(G = meanGSAcm2, W = meanMassG))


	#### calculate gill area index using d predicted from relationship in Pauly (1981) ####
	
	d_formula <- function(x){
	  d = 0.6742 + 0.03574 * x
	  return(d)
	}
	
	GAI_formula_pred <- function(G, W, x){
	  GAI = G / (W^x)
	  return(GAI)
	}
	
	mean_dat <- mean_dat %>%
	            rowwise() %>%
	            mutate(pred_d = d_formula(x = Log_Winfinity),
	            			 GAI_predd = GAI_formula_pred(G = meanGSAcm2, W = meanMassG, x = pred_d))