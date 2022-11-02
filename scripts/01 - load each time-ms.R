# 01-load all data

	# load libraries
	library(ape)
	library(beepr)
	library(broom)
	library(brms)
	library(car)
	library(cowplot)
	library(data.table)
	library(dplyr)
	library(forcats)
	library(forestmangr)
	#library(geiger)
	library(ggeffects)
	library(ggplot2)
	library(here)
	library(lme4)
#	library(lmodel2)
	library(loo)
	library(magick)
	library(MASS)
	#library(MCMCglmm)
	library(patchwork)
	library(phangorn)
	library(phytools)
	library(purrr)
	library(readr)
	library(reshape2)
#	library(rfishbase)
	library(rstan)
	library(stringr)
	library(tidyr)
	library(tidyverse)
#	library(ggdist)
	library(quantreg)
	
  
	
	`%notin%` <- Negate(`%in%`)

	# set options for stan
	rstan_options(auto_write = TRUE)
	options(mc.cores = parallel::detectCores())

	# helpful functions
	len <- function(y) {
	  z <-length(unique(y))
	  return(z)
	}
	
	# load gill surface area & growth data for 132 species
	GSA_growth_dat <- read.csv(here("data/GSA_Growth_132SP_FINAL.csv"),
                             header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)  %>%
                  	dplyr::select(., -X, -X.1, -X.2)

	# change "UNK" to mean in datasheet
	GSA_growth_dat$RawOrMean[GSA_growth_dat$RawOrMean == "UNK"] <- "mean"

	# transformations and standardizations
	GSA_growth_dat <- GSA_growth_dat %>%
                    mutate(GSAcm2 = case_when(
                              GSA_units %in% "cm2" ~ GSA,
                              GSA_units %in% "mm2" ~ GSA/100))

	GSA_growth_dat <- GSA_growth_dat %>%
	                  mutate(MassG = case_when(
	                              MassUnits %in% "grams" ~ Mass))
 
	GSA_growth_dat <- GSA_growth_dat %>%
	                  mutate(LogGSAcm2 = log10(GSAcm2),
	                         LogMassG = log10(Mass))

	# center body mass
	body_mass_avg <- GSA_growth_dat %>% group_by(Binomial) %>% summarise(mean_mass = mean(LogMassG))

	mean(body_mass_avg$mean_mass)

	GSA_growth_dat$LogCenterMassG <- GSA_growth_dat$LogMassG - (log10(300))

	# if a species has raw GSA data from more than one study, I need to only include that dataset
	# that has a larger range of body mass
	GSA_growth_dat$SpeciesStudy <- paste(GSA_growth_dat$Binomial, GSA_growth_dat$Source)

	# add a column for the number of determinations 
	GSA_growth_dat_study_sum <- GSA_growth_dat %>% 
	  count(Binomial, SpeciesStudy, RawOrMean)
	
	GSA_growth_dat_study_sum <- GSA_growth_dat_study_sum %>%
	                            arrange(Binomial, desc(n), .by_group = TRUE) %>%
	                            distinct(Binomial, .keep_all = TRUE)
                    
	to_keep <- GSA_growth_dat_study_sum$SpeciesStudy

	GSA_growth_dat <- GSA_growth_dat %>% filter(., SpeciesStudy %in% to_keep)

	# dataset for species that have at least 8 individual estimates per species
	RawGSA <- GSA_growth_dat %>%
		filter(RawOrMean == "raw")
	
	RawGSA8 <- RawGSA %>%
		group_by(Binomial) %>%
	  filter(n()>=8) 
	
	len(RawGSA8$Binomial)
	
	# load phylogeny
	fish_elasmo_supertree <- read.tree(file = here("./data/fish_elasmo_supertree.tre"))
	len(fish_elasmo_supertree$tip.label)
	# 125 species, missing 7 teleost species from the phylogeny (that's ok, just have to exclude them from )

	# which species are on the tree?
	Species <- lapply(strsplit(as.character(RawGSA8$Binomial), "\\ "), "[", 2)
	Genus <- lapply(strsplit(as.character(RawGSA8$Binomial), "\\ "), "[", 1)
	species_list <- tibble(Genus, Species)
	RawGSA8$phylo <- paste(species_list$Genus, species_list$Species,
	                                sep = "_")
	len(RawGSA8$Binomial)

	sp_drop <- setdiff(fish_elasmo_supertree$tip.label, RawGSA8$phylo)

	tree_pruned <- drop.tip(fish_elasmo_supertree, sp_drop)

	len(tree_pruned$tip.label) # 32 species in the dataset that are on the tree
	
	phylo_list <- tree_pruned$tip.label

	len(unique(RawGSA8$Binomial))
	sp_df <- unique(RawGSA8$Binomial)

	# dataset of species with raw data (at least n = 8) and on phylogeny
	RawGSA8_phylo <- filter(RawGSA8, phylo %in% phylo_list)

	setdiff(RawGSA8_phylo$phylo, tree_pruned$tip.label)
	setdiff(tree_pruned$tip.label, RawGSA8_phylo$phylo)
	
	## create an object for later use in filtering out species that are traditionally used in
	# aquaculture and those that are air-breathers

	aquaculture_sp_mean <- c("Ameiurus nebulosus",
	                          "Channa punctata",
	                          "Channa striata",
	                          "Cirrhinus mrigala",
	                          "Clarias batrachus",
	                          "Ctenopharyngodon idella",
	                          "Oncorhynchus mykiss",
	                          "Oreochromis niloticus",
	                          "Salmo trutta",
	                          "Tinca tinca")

	# species capable of air-breathing
	airbreathing_mean <- c("Anabas testudineus",
	                       "Anguilla anguilla",
	                       "Anguilla rostrata",
	                       "Boleophthalmus boddarti",
	                       "Channa punctata",
	                       "Channa striata",
	                       "Clarias batrachus",
	                       "Cobitis taenia",
	                       "Heteropneustes fossilis",
	                       "Lipophrys pholis",
	                       "Misgurnus fossilis",
	                       "Periophthalmus barbarus",
	                       "Periophthalmus chrysospilos",
	                       "Taurulus bubalis",
	                       "Zoarces viviparus")
	
	### Pauly's data
	Pauly_dat <- read.csv(here("data/Pauly_GAI_data.csv"), header = TRUE,
	                      strip.white = TRUE, stringsAsFactors = FALSE) %>% na.omit()
	
	Pauly_dat$LogGAI <- log10(Pauly_dat$GAI)

	Pauly_dat_no_coel <- Pauly_dat %>%
                    	 filter(., Species != "Latimeria chalumnae")
	# remove coelacanth

	Pauly_dat_no_coel$LogGAI <- log10(Pauly_dat_no_coel$GAI)

	Pauly_sp <- unique(Pauly_dat$BinomialResolved)

	OverlapPaulyNew <- GSA_growth_dat %>% filter(., Binomial %in% Pauly_sp)

	overlap_sp <- unique(OverlapPaulyNew$Binomial) 

	setdiff(Pauly_sp, overlap_sp)

