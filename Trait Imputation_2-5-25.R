#Script for skelevision data input and cleanup + imputation
#Written by Jacob S. Berv, jacob.berv@gmail.com

#library(devtools)
require(Rphylopars)
require(phytools)
require(ape)
require(vioplot)
library(dplyr)

sessionInfo()

#Record of working session

# > sessionInfo()
# R version 4.4.2 (2024-10-31)
# Platform: aarch64-apple-darwin24.1.0
# Running under: macOS Sequoia 15.1.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /opt/homebrew/Cellar/r/4.4.2_2/lib/R/lib/libRlapack.dylib;  LAPACK version 3.12.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/Los_Angeles
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] dplyr_1.1.4       vioplot_0.5.0     zoo_1.8-12        sm_2.2-6.0        phytools_2.3-0    maps_3.4.2.1      Rphylopars_0.3.10
# [8] ape_5.8          
# 
# loaded via a namespace (and not attached):
#   [1] fastmatch_1.1-4         gtable_0.3.6            phylolm_2.6.5           ggplot2_3.5.1           lattice_0.22-6         
# [6] numDeriv_2016.8-1.1     quadprog_1.5-8          vctrs_0.6.5             tools_4.4.2             generics_0.1.3         
# [11] parallel_4.4.2          tibble_3.2.1            pkgconfig_2.0.3         Matrix_1.7-1            scatterplot3d_0.3-44   
# [16] lifecycle_1.0.4         compiler_4.4.2          microbenchmark_1.5.0    munsell_0.5.1           mnormt_2.1.1           
# [21] combinat_0.0-8          codetools_0.2-20        pillar_1.10.1           tidyr_1.3.1             MASS_7.3-61            
# [26] clusterGeneration_1.3.8 iterators_1.0.14        boot_1.3-31             foreach_1.5.2           nlme_3.1-166           
# [31] parallelly_1.39.0       phangorn_2.12.1         Deriv_4.1.6             tidyselect_1.2.1        digest_0.6.37          
# [36] future_1.34.0           purrr_1.0.2             listenv_0.9.1           cowplot_1.1.3           grid_4.4.2             
# [41] colorspace_2.1-1        expm_1.0-0              cli_3.6.3               magrittr_2.0.3          optimParallel_1.0-2    
# [46] broom_1.0.7             future.apply_1.11.3     scales_1.3.0            backports_1.5.0         DEoptim_2.2-8          
# [51] modelr_0.1.11           globals_0.16.3          igraph_2.1.1            coda_0.19-4.1           doParallel_1.0.17      
# [56] doBy_4.6.24             rlang_1.1.5             Rcpp_1.0.13-1           glue_1.8.0              rstudioapi_0.17.1      
# [61] R6_2.5.1               


#set up working directory
setwd('~/Directory')
source('data_cleaning_functions.R')

#10_24_23 dataset from Brian Weeks
load('skelevision_vert_avo_combo_10-24-23')

#filter out anything not passeriformes
cv<-cv[cv$Order=="Passeriformes",]
length(unique(cv$phylo))
#2079 species

#read in the consensus tree from brian
jetztest<- read.tree(file="Bird_Phylogeny_8-17-21.tre")

#which tips in skelevision v1 do not match the tip names in jetz tree
nomatch<-unique(cv$phylo[!cv$phylo %in% jetztest$tip.label])
#manually fixing the names to increase matching between jetz and the dataset

{
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'brissonii')$matching_indices[1:2]] <- "Cyanocompsa_brissonii"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'brunneinucha')$matching_indices[1:13]] <- "Arremon_brunneinucha" #done
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Casiornis')$matching_indices] <- "Casiornis_rufus"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'rixosus')$matching_indices] <- "Machetornis_rixosa"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Cissopis')$matching_indices] <- "Cissopis_leverianus"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Automolus_dorsalis')$matching_indices] <- "Anabazenops_dorsalis"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Todirostrum_latirostre')$matching_indices] <- "Poecilotriccus_latirostris"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Atlapetes_torquatus')$matching_indices] <- "Arremon_torquatus"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'parellina')$matching_indices] <- "Cyanocompsa_parellina"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Certhiaxis_cinnamomea')$matching_indices] <- "Certhiaxis_cinnamomeus"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Pipra_coronata')$matching_indices] <- "Lepidothrix_coronata"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Chloropipo_holochlora')$matching_indices] <- "Xenopipo_holochlora"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Schiffornis_turdinus')$matching_indices] <- "Schiffornis_turdina"
  #cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Passerina_cyanoides')$matching_indices] <- "Cyanoloxia_glaucocaerulea"
  #Passerina_cyanoides ## need to revisit
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Todirostrum_calopterum')$matching_indices] <- "Poecilotriccus_calopterus"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Ramphotrigon_megacephala')$matching_indices] <- "Ramphotrigon_megacephalum"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Pyrrhomyias_cinnamomea')$matching_indices] <- "Pyrrhomyias_cinnamomeus"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Pipra_coeruleocapilla')$matching_indices] <- "Lepidothrix_coeruleocapilla"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Fluvicola_leucocephala')$matching_indices] <- "Arundinicola_leucocephala"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Phylloscartes_flaveolus')$matching_indices] <- "Capsiempis_flaveola"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Lipaugus_subalaris')$matching_indices] <- "Snowornis_subalaris"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Vermivora_gutturalis')$matching_indices] <- "Parula_gutturalis"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Metopothrix_aurantiacus')$matching_indices] <- "Metopothrix_aurantiaca"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Zimmerius_cinereicapillus')$matching_indices] <- "Zimmerius_cinereicapilla"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Pitylus_grossus')$matching_indices] <- "Saltator_grossus"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Muscisaxicola_alpina')$matching_indices] <- "Muscisaxicola_alpinus"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Automolus_roraimae')$matching_indices] <- "Syndactyla_roraimae"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Pipra_serena')$matching_indices] <- "Lepidothrix_serena"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Pipra_nattereri')$matching_indices] <- "Lepidothrix_nattereri"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Caryothraustes_humeralis')$matching_indices] <- "Parkerthraustes_humeralis"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Schoeniophylax_phryganophila')$matching_indices] <- "Schoeniophylax_phryganophilus"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Todirostrum_sylvia')$matching_indices] <- "Poecilotriccus_sylvia"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Todirostrum_russatum')$matching_indices] <- "Poecilotriccus_russatus"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Chloropipo_uniformis')$matching_indices] <- "Xenopipo_uniformis"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Tiaris_olivacea')$matching_indices] <- "Tiaris_olivaceus"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Tiaris_canora')$matching_indices] <- "Tiaris_canorus"
  cv$phylo[find_string_matches(string_vector = cv$phylo, search_string = 'Ammodramus_sandwichensis')$matching_indices] <- "Passerculus_sandwichensis"
}
unique(cv$phylo[!cv$phylo %in% jetztest$tip.label])
#cv[find_string_matches(string_vector = cv$phylo, search_string = 'Cyanocompsa_brissonii')$matching_indices,] #$matching_indices[1:2]
#table(cv[find_string_matches(string_vector = cv$phylo, search_string = 'Vidua')$matching_indices,]$phylo)

#human readable name change summary (dataset name changes to match jetz tree names)
{
  # Passerina_brissonii -> Cyanocompsa_brissonii
  # Atlapetes_brunneinucha -> Arremon_brunneinucha
  # Casiornis_rufa -> Casiornis_rufus
  # Machetornis_rixosus -> Machetornis_rixosa
  # Cissopis_leveriana -> Cissopis_leverianus
  # Automolus_dorsalis -> Anabazenops_dorsalis
  # Todirostrum_latirostre -> Poecilotriccus_latirostris
  # Atlapetes_torquatus -> Arremon_torquatus
  # Passerina_parellina -> Cyanocompsa_parellina
  # Certhiaxis_cinnamomea -> Certhiaxis_cinnamomeus
  # Pipra_coronata -> Lepidothrix_coronata
  # Chloropipo_holochlora -> Xenopipo_holochlora
  # Schiffornis_turdinus -> Schiffornis_turdina
  # Todirostrum_calopterum -> Poecilotriccus_calopterus
  # Ramphotrigon_megacephala -> Ramphotrigon_megacephalum
  # Pyrrhomyias_cinnamomea -> Pyrrhomyias_cinnamomeus
  # Pipra_coeruleocapilla -> Lepidothrix_coeruleocapilla
  # Fluvicola_leucocephala -> Arundinicola_leucocephala
  # Phylloscartes_flaveolus -> Capsiempis_flaveola
  # Lipaugus_subalaris -> Snowornis_subalaris
  # Vermivora_gutturalis -> Parula_gutturalis
  # Metopothrix_aurantiacus -> Metopothrix_aurantiaca
  # Zimmerius_cinereicapillus -> Zimmerius_cinereicapilla
  # Pitylus_grossus -> Saltator_grossus
  # Muscisaxicola_alpina -> Muscisaxicola_alpinus
  # Automolus_roraimae -> Syndactyla_roraimae
  # Pipra_serena -> Lepidothrix_serena
  # Pipra_nattereri -> Lepidothrix_nattereri
  # Caryothraustes_humeralis -> Parkerthraustes_humeralis
  # Schoeniophylax_phryganophila -> Schoeniophylax_phryganophilus
  # Todirostrum_sylvia -> Poecilotriccus_sylvia
  # Todirostrum_russatum -> Poecilotriccus_russatus
  # Chloropipo_uniformis -> Xenopipo_uniformis
  # Tiaris_olivacea -> Tiaris_olivaceus
  # Tiaris_canora -> Tiaris_canorus
  # Ammodramus_sandwichensis -> Passerculus_sandwichensis
  # Passerina_cyanoides is in the dataset, but there is no synonymous tip in jetz, cutting. 
}

#filter out rows that do not match the consensus tree
cv<-cv[cv$phylo %in% jetztest$tip.label,]
length(unique(cv$phylo))
#2057 species remaining

#pruning the tree to retain only the tips we can easily match
tree.pruned<-keep.tip(jetztest,tip =  cv$phylo)

#dataset filtering and re-assembly
taxonlist<-cv$phylo
ID<-cv$museum_number #cv$bone_id
translation<-as.data.frame(cbind(phylo=cv$phylo, ID=cv$museum_number))
prob_filter<-0.95

#filtering the individual datasets to retain only those with 0.95 bprob
{
  #carpometacarpus
  {
    carpometacarpus.1<-data.frame(ID, taxonlist, cv$carpometacarpus.1, cv$carpometacarpus.1.bprob)
    carpometacarpus.1$cv.carpometacarpus.1<-ifelse(carpometacarpus.1$cv.carpometacarpus.1.bprob > prob_filter, cv$carpometacarpus.1, NA) #mask those entries without prob_filter
    carpometacarpus.1$cv.carpometacarpus.1.bprob<-NULL # clean the probability column
    
    carpometacarpus.2<-data.frame(ID, taxonlist, cv$carpometacarpus.2, cv$carpometacarpus.2.bprob)
    carpometacarpus.2$cv.carpometacarpus.2<-ifelse(carpometacarpus.2$cv.carpometacarpus.2.bprob > prob_filter, cv$carpometacarpus.2, NA) #mask those entries without prob_filter
    carpometacarpus.2$cv.carpometacarpus.2.bprob<-NULL # clean the probability column
    
    #take the average of both bones when there are 2
    carpometacarpus<-merge(carpometacarpus.1, carpometacarpus.2, by='ID', all=TRUE)
    carpometacarpus$taxonlist.y<-NULL #clean redundant labels
    carpometacarpus$taxonlist.x<-NULL #clean redundant labels
    carpometacarpus$carpometacarpus <- rowMeans(cbind(carpometacarpus$cv.carpometacarpus.1, carpometacarpus$cv.carpometacarpus.2), na.rm=T)
    
    
    #generate a new dataset with standard errors and means per species
    carpometacarpus_se<-merge(translation,carpometacarpus, by="ID")[,c(1,2,5)]
    library(dplyr)
    carpometacarpus_se <- carpometacarpus_se %>%
      group_by(phylo) %>%
      summarize(
        carpometacarpus.Mean = mean(carpometacarpus),
        carpometacarpus.Mean.log = log(mean(carpometacarpus)),
        carpometacarpus.StdErr = sd(carpometacarpus) / sqrt(n()),
        carpometacarpus.StdErr.log = sd(log(carpometacarpus)) / sqrt(n())
      )
    carpometacarpus_se<-as.data.frame(carpometacarpus_se)
    #replace NAs with the mean value of SE for each column, per family
    carpometacarpus_se$carpometacarpus.StdErr.complete <- zoo::na.aggregate(carpometacarpus_se$carpometacarpus.StdErr)
    carpometacarpus_se$carpometacarpus.StdErr.log.complete <- zoo::na.aggregate(carpometacarpus_se$carpometacarpus.StdErr.log)
  }
  
  #ifelse(is.na(carpometacarpus_se$carpometacarpus.StdErr), NA, carpometacarpus_se$carpometacarpus.StdErr.complete)
  #carpometacarpus[is.na(carpometacarpus$cv.carpometacarpus.1) & !is.na(carpometacarpus$cv.carpometacarpus.2), ]
  #no rows where first column has NA and second column has a number
  
  #coracoid, no data
  {
    coracoid.1<-data.frame(ID, taxonlist, cv$coracoid.1, cv$coracoid.1.bprob)
    coracoid.1$cv.coracoid.1<-ifelse(coracoid.1$cv.coracoid.1.bprob > prob_filter, cv$coracoid.1, NA) #mask those entries without prob_filter
    #coracoid.1<-coracoid.1[coracoid.1$cv.coracoid.1.bprob > prob_filter,] #filter out rows without confidence
    #coracoid.1<-coracoid.1[complete.cases(coracoid.1),] #filter out rows that don't have data
    coracoid.1$cv.coracoid.1.bprob<-NULL # clean the probability column
    #no data for coracoid
  }
  
  #femur
  {
    femur.1<-data.frame(ID, taxonlist, cv$femur.1, cv$femur.1.bprob)
    #femur.1<-femur.1[femur.1$cv.femur.1.bprob > prob_filter,] #filter out rows without confidence
    #femur.1<-femur.1[complete.cases(femur.1),] #filter out rows that don't have data
    femur.1$cv.femur.1<-ifelse(femur.1$cv.femur.1.bprob > prob_filter, cv$femur.1, NA) #mask those entries without prob_filter
    femur.1$cv.femur.1.bprob<-NULL # clean the probability column
    
    femur.2<-data.frame(ID, taxonlist, cv$femur.2, cv$femur.2.bprob)
    #femur.2<-femur.2[femur.2$cv.femur.2.bprob > prob_filter,] #filter out rows without confidence
    #femur.2<-femur.2[complete.cases(femur.2),] #filter out rows that don't have data
    femur.2$cv.femur.2<-ifelse(femur.2$cv.femur.2.bprob > prob_filter, cv$femur.2, NA) #mask those entries without prob_filter
    femur.2$cv.femur.2.bprob<-NULL # clean the probability column
    
    #take the average of both bones when there are 2
    femur<-merge(femur.1, femur.2, by='ID', all=TRUE)
    femur$taxonlist.y<-NULL #clean redundant labels
    femur$taxonlist.x<-NULL #clean redundant labels
    femur$femur <- rowMeans(cbind(femur$cv.femur.1, femur$cv.femur.2), na.rm=T)
    
    #generate a new dataset with standard errors and means per species
    femur_se<-merge(translation,femur, by="ID")[,c(1,2,5)]
    library(dplyr)
    femur_se <- femur_se %>%
      group_by(phylo) %>%
      summarize(
        femur.Mean = mean(femur),
        femur.Mean.log = log(mean(femur)),
        femur.StdErr = sd(femur) / sqrt(n()),
        femur.StdErr.log = sd(log(femur)) / sqrt(n())
      )
    femur_se<-as.data.frame(femur_se)
    #replace NAs with the mean value of SE for each column, per family
    femur_se$femur.StdErr.complete <- zoo::na.aggregate(femur_se$femur.StdErr)
    femur_se$femur.StdErr.log.complete <- zoo::na.aggregate(femur_se$femur.StdErr.log)
  }
  
  #femur[is.na(femur$cv.femur.1) & !is.na(femur$cv.femur.2), ]
  #no rows where first column has NA and second column has a number
  
  #fibula, no data
  {
    fibula.1<-data.frame(ID, taxonlist, cv$fibula.1, cv$fibula.1.bprob)
    #fibula.1<-fibula.1[fibula.1$cv.fibula.1.bprob > prob_filter,] #filter out rows without confidence
    #fibula.1<-fibula.1[complete.cases(fibula.1),] #filter out rows that don't have data
    fibula.1$cv.fibula.1<-ifelse(fibula.1$cv.fibula.1.bprob > prob_filter, cv$fibula.1, NA) #mask those entries without prob_filter
    fibula.1$cv.fibula.1.bprob<-NULL # clean the probability column
    
    fibula.2<-data.frame(ID, taxonlist, cv$fibula.2, cv$fibula.2.bprob)
    #fibula.2<-fibula.2[fibula.2$cv.fibula.2.bprob > prob_filter,] #filter out rows without confidence
    #fibula.2<-fibula.2[complete.cases(fibula.2),] #filter out rows that don't have data
    fibula.2$cv.fibula.2<-ifelse(fibula.2$cv.fibula.2.bprob > prob_filter, cv$fibula.2, NA) #mask those entries without prob_filter
    fibula.2$cv.fibula.2.bprob<-NULL # clean the probability column
    #no data for fibula
  }
  
  #furcula
  {
    furcula.1<-data.frame(ID, taxonlist, cv$furcula.1, cv$furcula.1.bprob)
    #furcula.1<-furcula.1[furcula.1$cv.furcula.1.bprob > prob_filter,] #filter out rows without confidence
    #furcula.1<-furcula.1[complete.cases(furcula.1),] #filter out rows that don't have data
    furcula.1$cv.furcula.1<-ifelse(furcula.1$cv.furcula.1.bprob > prob_filter, cv$furcula.1, NA) #mask those entries without prob_filter
    furcula.1$cv.furcula.1.bprob<-NULL # clean the probability column
    
    #generate a new dataset with standard errors and means per species
    furcula.1_se<-merge(translation,furcula.1, by="ID")[,c(1,2,4)]
    library(dplyr)
    furcula.1_se <- furcula.1_se %>%
      group_by(phylo) %>%
      summarize(
        furcula.Mean = mean(cv.furcula.1),
        furcula.Mean.log = log(mean(cv.furcula.1)),
        furcula.StdErr = sd(cv.furcula.1) / sqrt(n()),
        furcula.StdErr.log = sd(log(cv.furcula.1)) / sqrt(n())
      )
    furcula.1_se<-as.data.frame(furcula.1_se)
    #replace NAs with the mean value of SE for each column, per family
    furcula.1_se$furcula.StdErr.complete <- zoo::na.aggregate(furcula.1_se$furcula.StdErr)
    furcula.1_se$furcula.StdErr.log.complete <- zoo::na.aggregate(furcula.1_se$furcula.StdErr.log)
  }
  
  #humerus
  {
    humerus.1<-data.frame(ID, taxonlist, cv$humerus.1, cv$humerus.1.bprob)
    #humerus.1<-humerus.1[humerus.1$cv.humerus.1.bprob > prob_filter,] #filter out rows without confidence
    #humerus.1<-humerus.1[complete.cases(humerus.1),] #filter out rows that don't have data
    humerus.1$cv.humerus.1<-ifelse(humerus.1$cv.humerus.1.bprob > prob_filter, cv$humerus.1, NA) #mask those entries without prob_filter
    humerus.1$cv.humerus.1.bprob<-NULL # clean the probability column
    
    humerus.2<-data.frame(ID, taxonlist, cv$humerus.2, cv$humerus.2.bprob)
    #humerus.2<-humerus.2[humerus.2$cv.humerus.2.bprob > prob_filter,] #filter out rows without confidence
    #humerus.2<-humerus.2[complete.cases(humerus.2),] #filter out rows that don't have data
    humerus.2$cv.humerus.2<-ifelse(humerus.2$cv.humerus.2.bprob > prob_filter, cv$humerus.2, NA) #mask those entries without prob_filter
    humerus.2$cv.humerus.2.bprob<-NULL # clean the probability column
    
    #take the average of both bones when there are 2
    humerus<-merge(humerus.1, humerus.2, by='ID', all=TRUE)
    humerus$taxonlist.y<-NULL #clean redundant labels
    humerus$taxonlist.x<-NULL #clean redundant labels
    humerus$humerus <- rowMeans(cbind(humerus$cv.humerus.1, humerus$cv.humerus.2), na.rm=T)
    
    #generate a new dataset with standard errors and means per species
    humerus_se<-merge(translation,humerus, by="ID")[,c(1,2,5)]
    library(dplyr)
    humerus_se <- humerus_se %>%
      group_by(phylo) %>%
      summarize(
        humerus.Mean = mean(humerus),
        humerus.Mean.log = log(mean(humerus)),
        humerus.StdErr = sd(humerus) / sqrt(n()),
        humerus.StdErr.log = sd(log(humerus)) / sqrt(n())
      )
    humerus_se<-as.data.frame(humerus_se)
    #replace NAs with the mean value of SE for each column, per family
    humerus_se$humerus.StdErr.complete <- zoo::na.aggregate(humerus_se$humerus.StdErr)
    humerus_se$humerus.StdErr.log.complete <- zoo::na.aggregate(humerus_se$humerus.StdErr.log)
    
  }
  
  #keel
  {
    keel.1<-data.frame(ID, taxonlist, cv$keel.1, cv$keel.1.bprob)
    #keel.1<-keel.1[keel.1$cv.keel.1.bprob > prob_filter,] #filter out rows without confidence
    #keel.1<-keel.1[complete.cases(keel.1),] #filter out rows that don't have data
    keel.1$cv.keel.1<-ifelse(keel.1$cv.keel.1.bprob > prob_filter, cv$keel.1, NA) #mask those entries without prob_filter
    keel.1$cv.keel.1.bprob<-NULL # clean the probability column
    
    #generate a new dataset with standard errors and means per species
    keel.1_se<-merge(translation,keel.1, by="ID")[,c(1,2,4)]
    library(dplyr)
    keel.1_se <- keel.1_se %>%
      group_by(phylo) %>%
      summarize(
        keel.Mean = mean(cv.keel.1),
        keel.Mean.log = log(mean(cv.keel.1)),
        keel.StdErr = sd(cv.keel.1) / sqrt(n()),
        keel.StdErr.log = sd(log(cv.keel.1)) / sqrt(n())
      )
    keel.1_se<-as.data.frame(keel.1_se)
    #replace NAs with the mean value of SE for each column, per family
    keel.1_se$keel.StdErr.complete <- zoo::na.aggregate(keel.1_se$keel.StdErr)
    keel.1_se$keel.StdErr.log.complete <- zoo::na.aggregate(keel.1_se$keel.StdErr.log)
    
  }
  
  #metatarsus
  {
    metatarsus.1<-data.frame(ID, taxonlist, cv$metatarsus.1, cv$metatarsus.1.bprob)
    #metatarsus.1<-metatarsus.1[metatarsus.1$cv.metatarsus.1.bprob > prob_filter,] #filter out rows without confidence
    #metatarsus.1<-metatarsus.1[complete.cases(metatarsus.1),] #filter out rows that don't have data
    metatarsus.1$cv.metatarsus.1<-ifelse(metatarsus.1$cv.metatarsus.1.bprob > prob_filter, cv$metatarsus.1, NA) #mask those entries without prob_filter
    metatarsus.1$cv.metatarsus.1.bprob<-NULL # clean the probability column
    
    metatarsus.2<-data.frame(ID, taxonlist, cv$metatarsus.2, cv$metatarsus.2.bprob)
    #metatarsus.2<-metatarsus.2[metatarsus.2$cv.metatarsus.2.bprob > prob_filter,] #filter out rows without confidence
    #metatarsus.2<-metatarsus.2[complete.cases(metatarsus.2),] #filter out rows that don't have data
    metatarsus.2$cv.metatarsus.2<-ifelse(metatarsus.2$cv.metatarsus.2.bprob > prob_filter, cv$metatarsus.2, NA) #mask those entries without prob_filter
    metatarsus.2$cv.metatarsus.2.bprob<-NULL # clean the probability column
    
    #take the average of both bones when there are 2
    metatarsus<-merge(metatarsus.1, metatarsus.2, by='ID', all=TRUE)
    metatarsus$taxonlist.y<-NULL #clean redundant labels
    metatarsus$taxonlist.x<-NULL #clean redundant labels
    metatarsus$metatarsus <- rowMeans(cbind(metatarsus$cv.metatarsus.1, metatarsus$cv.metatarsus.2), na.rm=T)
    
    #generate a new dataset with standard errors and means per species
    metatarsus_se<-merge(translation,metatarsus, by="ID")[,c(1,2,5)]
    library(dplyr)
    metatarsus_se <- metatarsus_se %>%
      group_by(phylo) %>%
      summarize(
        metatarsus.Mean = mean(metatarsus),
        metatarsus.Mean.log = log(mean(metatarsus)),
        metatarsus.StdErr = sd(metatarsus) / sqrt(n()),
        metatarsus.StdErr.log = sd(log(metatarsus)) / sqrt(n())
      )
    metatarsus_se<-as.data.frame(metatarsus_se)
    #replace NAs with the mean value of SE for each column, per family
    metatarsus_se$metatarsus.StdErr.complete <- zoo::na.aggregate(metatarsus_se$metatarsus.StdErr)
    metatarsus_se$metatarsus.StdErr.log.complete <- zoo::na.aggregate(metatarsus_se$metatarsus.StdErr.log)
    
    
  }
  
  #sclerotic_ring
  {
    sclerotic_ring.1<-data.frame(ID, taxonlist, cv$sclerotic_ring.1, cv$sclerotic_ring.1.bprob)
    #sclerotic_ring.1<-sclerotic_ring.1[sclerotic_ring.1$cv.sclerotic_ring.1.bprob > prob_filter,] #filter out rows without confidence
    #sclerotic_ring.1<-sclerotic_ring.1[complete.cases(sclerotic_ring.1),] #filter out rows that don't have data
    sclerotic_ring.1$cv.sclerotic_ring.1<-ifelse(sclerotic_ring.1$cv.sclerotic_ring.1.bprob > prob_filter, cv$sclerotic_ring.1, NA) #mask those entries without prob_filter
    sclerotic_ring.1$cv.sclerotic_ring.1.bprob<-NULL # clean the probability column
    
    sclerotic_ring.2<-data.frame(ID, taxonlist, cv$sclerotic_ring.2, cv$sclerotic_ring.2.bprob)
    #sclerotic_ring.2<-sclerotic_ring.2[sclerotic_ring.2$cv.sclerotic_ring.2.bprob > prob_filter,] #filter out rows without confidence
    #sclerotic_ring.2<-sclerotic_ring.2[complete.cases(sclerotic_ring.2),] #filter out rows that don't have data
    sclerotic_ring.2$cv.sclerotic_ring.2<-ifelse(sclerotic_ring.2$cv.sclerotic_ring.2.bprob > prob_filter, cv$sclerotic_ring.2, NA) #mask those entries without prob_filter
    sclerotic_ring.2$cv.sclerotic_ring.2.bprob<-NULL # clean the probability column
    
    #take the average of both bones when there are 2
    sclerotic_ring<-merge(sclerotic_ring.1, sclerotic_ring.2, by='ID', all=TRUE)
    sclerotic_ring$taxonlist.y<-NULL #clean redundant labels
    sclerotic_ring$taxonlist.x<-NULL #clean redundant labels
    sclerotic_ring$sclerotic_ring <- rowMeans(cbind(sclerotic_ring$cv.sclerotic_ring.1, sclerotic_ring$cv.sclerotic_ring.2), na.rm=T)
    
    #generate a new dataset with standard errors and means per species
    sclerotic_ring_se<-merge(translation,sclerotic_ring, by="ID")[,c(1,2,5)]
    library(dplyr)
    sclerotic_ring_se <- sclerotic_ring_se %>%
      group_by(phylo) %>%
      summarize(
        sclerotic_ring.Mean = mean(sclerotic_ring),
        sclerotic_ring.Mean.log = log(mean(sclerotic_ring)),
        sclerotic_ring.StdErr = sd(sclerotic_ring) / sqrt(n()),
        sclerotic_ring.StdErr.log = sd(log(sclerotic_ring)) / sqrt(n())
      )
    sclerotic_ring_se<-as.data.frame(sclerotic_ring_se)
    #replace NAs with the mean value of SE for each column, per family
    sclerotic_ring_se$sclerotic_ring.StdErr.complete <- zoo::na.aggregate(sclerotic_ring_se$sclerotic_ring.StdErr)
    sclerotic_ring_se$sclerotic_ring.StdErr.log.complete <- zoo::na.aggregate(sclerotic_ring_se$sclerotic_ring.StdErr.log)
    
    
  }
  
  #second_digit
  {
    second_digit_p_1.1<-data.frame(ID, taxonlist, cv$second_digit_p_1.1, cv$second_digit_p_1.1.bprob)
    #second_digit_p_1.1<-second_digit_p_1.1[second_digit_p_1.1$cv.second_digit_p_1.1.bprob > prob_filter,] #filter out rows without confidence
    #second_digit_p_1.1<-second_digit_p_1.1[complete.cases(second_digit_p_1.1),] #filter out rows that don't have data
    second_digit_p_1.1$cv.second_digit_p_1.1<-ifelse(second_digit_p_1.1$cv.second_digit_p_1.1.bprob > prob_filter, cv$second_digit_p_1.1, NA) #mask those entries without prob_filter
    second_digit_p_1.1$cv.second_digit_p_1.1.bprob<-NULL # clean the probability column
    
    second_digit_p_2.1<-data.frame(ID, taxonlist, cv$second_digit_p_2.1, cv$second_digit_p_2.1.bprob)
    #second_digit_p_2.1<-second_digit_p_2.1[second_digit_p_2.1$cv.second_digit_p_2.1.bprob > prob_filter,] #filter out rows without confidence
    #second_digit_p_2.1<-second_digit_p_2.1[complete.cases(second_digit_p_2.1),] #filter out rows that don't have data
    second_digit_p_2.1$cv.second_digit_p_2.1<-ifelse(second_digit_p_2.1$cv.second_digit_p_2.1.bprob > prob_filter, cv$second_digit_p_2.1, NA) #mask those entries without prob_filter
    second_digit_p_2.1$cv.second_digit_p_2.1.bprob<-NULL # clean the probability column
    
    #take the average of both bones when there are 2
    second_digit<-merge(second_digit_p_1.1, second_digit_p_2.1, by='ID', all=TRUE)
    second_digit$taxonlist.y<-NULL #clean redundant labels
    second_digit$taxonlist.x<-NULL #clean redundant labels
    second_digit$second_digit <- rowMeans(cbind(second_digit$cv.second_digit_p_1, sclerotic_ring$cv.second_digit_p_2), na.rm=T)
    
    #generate a new dataset with standard errors and means per species
    second_digit_se<-merge(translation,second_digit, by="ID")[,c(1,2,5)]
    library(dplyr)
    second_digit_se <- second_digit_se %>%
      group_by(phylo) %>%
      summarize(
        second_digit.Mean = mean(second_digit),
        second_digit.Mean.log = log(mean(second_digit)),
        second_digit.StdErr = sd(second_digit) / sqrt(n()),
        second_digit.StdErr.log = sd(log(second_digit)) / sqrt(n())
      )
    second_digit_se<-as.data.frame(second_digit_se)
    #replace NAs with the mean value of SE for each column, per family
    second_digit_se$second_digit.StdErr.complete <- zoo::na.aggregate(second_digit_se$second_digit.StdErr)
    second_digit_se$second_digit.StdErr.log.complete <- zoo::na.aggregate(second_digit_se$second_digit.StdErr.log)
    
  }
  
  #skull
  {
    skull.1<-data.frame(ID, taxonlist, cv$skull.1, cv$skull.1.bprob)
    #skull.1<-skull.1[skull.1$cv.skull.1.bprob > prob_filter,] #filter out rows without confidence
    #skull.1<-skull.1[complete.cases(skull.1),] #filter out rows that don't have data
    skull.1$cv.skull.1<-ifelse(skull.1$cv.skull.1.bprob > prob_filter, cv$skull.1, NA) #mask those entries without prob_filter
    skull.1$cv.skull.1.bprob<-NULL # clean the probability column
    
    #generate a new dataset with standard errors and means per species
    skull.1_se<-merge(translation,skull.1, by="ID")[,c(1,2,4)]
    library(dplyr)
    skull.1_se <- skull.1_se %>%
      group_by(phylo) %>%
      summarize(
        skull.Mean = mean(cv.skull.1),
        skull.Mean.log = log(mean(cv.skull.1)),
        skull.StdErr = sd(cv.skull.1) / sqrt(n()),
        skull.StdErr.log = sd(log(cv.skull.1)) / sqrt(n())
      )
    skull.1_se<-as.data.frame(skull.1_se)
    #replace NAs with the mean value of SE for each column, per family
    skull.1_se$skull.StdErr.complete <- zoo::na.aggregate(skull.1_se$skull.StdErr)
    skull.1_se$skull.StdErr.log.complete <- zoo::na.aggregate(skull.1_se$skull.StdErr.log)
    
    
  }
  
  #sternum, no data
  {
    sternum.1<-data.frame(ID, taxonlist, cv$sternum.1, cv$sternum.1.bprob)
    #sternum.1<-sternum.1[sternum.1$cv.sternum.1.bprob > prob_filter,] #filter out rows without confidence
    #sternum.1<-sternum.1[complete.cases(sternum.1),] #filter out rows that don't have data
    sternum.1$cv.sternum.1<-ifelse(sternum.1$cv.sternum.1.bprob > prob_filter, cv$sternum.1, NA) #mask those entries without prob_filter
    sternum.1$cv.sternum.1.bprob<-NULL # clean the probability column
    #no data for sternum
  }
  
  #tarsus
  {
    tarsus.1<-data.frame(ID, taxonlist, cv$tarsus.1, cv$tarsus.1.bprob)
    #tarsus.1<-tarsus.1[tarsus.1$cv.tarsus.1.bprob > prob_filter,] #filter out rows without confidence
    #tarsus.1<-tarsus.1[complete.cases(tarsus.1),] #filter out rows that don't have data
    tarsus.1$cv.tarsus.1<-ifelse(tarsus.1$cv.tarsus.1.bprob > prob_filter, cv$tarsus.1, NA) #mask those entries without prob_filter
    tarsus.1$cv.tarsus.1.bprob<-NULL # clean the probability column
    
    tarsus.2<-data.frame(ID, taxonlist, cv$tarsus.2, cv$tarsus.2.bprob)
    #tarsus.2<-tarsus.2[tarsus.2$cv.tarsus.2.bprob > prob_filter,] #filter out rows without confidence
    #tarsus.2<-tarsus.2[complete.cases(tarsus.2),] #filter out rows that don't have data
    tarsus.2$cv.tarsus.2<-ifelse(tarsus.2$cv.tarsus.2.bprob > prob_filter, cv$tarsus.2, NA) #mask those entries without prob_filter
    tarsus.2$cv.tarsus.2.bprob<-NULL # clean the probability column
    
    #take the average of both bones when there are 2
    tarsus<-merge(tarsus.1, tarsus.2, by='ID', all=TRUE)
    tarsus$taxonlist.y<-NULL #clean redundant labels
    tarsus$taxonlist.x<-NULL #clean redundant labels
    tarsus$tarsus <- rowMeans(cbind(tarsus$cv.tarsus.1, tarsus$cv.tarsus.2), na.rm=T)
    
    #generate a new dataset with standard errors and means per species
    tarsus_se<-merge(translation,tarsus, by="ID")[,c(1,2,5)]
    library(dplyr)
    tarsus_se <- tarsus_se %>%
      group_by(phylo) %>%
      summarize(
        tarsus.Mean = mean(tarsus),
        tarsus.Mean.log = log(mean(tarsus)),
        tarsus.StdErr = sd(tarsus) / sqrt(n()),
        tarsus.StdErr.log = sd(log(tarsus)) / sqrt(n())
      )
    tarsus_se<-as.data.frame(tarsus_se)
    #replace NAs with the mean value of SE for each column, per family
    tarsus_se$tarsus.StdErr.complete <- zoo::na.aggregate(tarsus_se$tarsus.StdErr)
    tarsus_se$tarsus.StdErr.log.complete <- zoo::na.aggregate(tarsus_se$tarsus.StdErr.log)
    
  }
  
  #ulna
  {
    ulna.1<-data.frame(ID, taxonlist, cv$ulna.1, cv$ulna.1.bprob)
    #ulna.1<-ulna.1[ulna.1$cv.ulna.1.bprob > prob_filter,] #filter out rows without confidence
    #ulna.1<-ulna.1[complete.cases(ulna.1),] #filter out rows that don't have data
    ulna.1$cv.ulna.1<-ifelse(ulna.1$cv.ulna.1.bprob > prob_filter, cv$ulna.1, NA) #mask those entries without prob_filter
    ulna.1$cv.ulna.1.bprob<-NULL # clean the probability column
    
    ulna.2<-data.frame(ID, taxonlist, cv$ulna.2, cv$ulna.2.bprob)
    #ulna.2<-ulna.2[ulna.2$cv.ulna.2.bprob > prob_filter,] #filter out rows without confidence
    #ulna.2<-ulna.2[complete.cases(ulna.2),] #filter out rows that don't have data
    ulna.2$cv.ulna.2<-ifelse(ulna.2$cv.ulna.2.bprob > prob_filter, cv$ulna.2, NA) #mask those entries without prob_filter
    ulna.2$cv.ulna.2.bprob<-NULL # clean the probability column
    
    #take the average of both bones when there are 2
    ulna<-merge(ulna.1, ulna.2, by='ID', all=TRUE)
    ulna$taxonlist.y<-NULL #clean redundant labels
    ulna$taxonlist.x<-NULL #clean redundant labels
    ulna$ulna <- rowMeans(cbind(ulna$cv.ulna.1, ulna$cv.ulna.2), na.rm=T)
    
    #generate a new dataset with standard errors and means per species
    ulna_se<-merge(translation,ulna, by="ID")[,c(1,2,5)]
    library(dplyr)
    ulna_se <- ulna_se %>%
      group_by(phylo) %>%
      summarize(
        ulna.Mean = mean(ulna),
        ulna.Mean.log = log(mean(ulna)),
        ulna.StdErr = sd(ulna) / sqrt(n()),
        ulna.StdErr.log = sd(log(ulna)) / sqrt(n())
      )
    ulna_se<-as.data.frame(ulna_se)
    #replace NAs with the mean value of SE for each column, per family
    ulna_se$ulna.StdErr.complete <- zoo::na.aggregate(ulna_se$ulna.StdErr)
    ulna_se$ulna.StdErr.log.complete <- zoo::na.aggregate(ulna_se$ulna.StdErr.log)
    
  }
  
  #radius
  {
    radius.1<-data.frame(ID, taxonlist, cv$radius.1, cv$radius.1.bprob)
    #radius.1<-radius.1[radius.1$cv.radius.1.bprob > prob_filter,] #filter out rows without confidence
    #radius.1<-radius.1[complete.cases(radius.1),] #filter out rows that don't have data
    radius.1$cv.radius.1<-ifelse(radius.1$cv.radius.1.bprob > prob_filter, cv$radius.1, NA) #mask those entries without prob_filter
    radius.1$cv.radius.1.bprob<-NULL # clean the probability column
    
    radius.2<-data.frame(ID, taxonlist, cv$radius.2, cv$radius.2.bprob)
    #radius.2<-radius.2[radius.2$cv.radius.2.bprob > prob_filter,] #filter out rows without confidence
    #radius.2<-radius.2[complete.cases(radius.2),] #filter out rows that don't have data
    radius.2$cv.radius.2<-ifelse(radius.2$cv.radius.2.bprob > prob_filter, cv$radius.2, NA) #mask those entries without prob_filter
    radius.2$cv.radius.2.bprob<-NULL # clean the probability column
    
    #take the average of both bones when there are 2
    radius<-merge(radius.1, radius.2, by='ID', all=TRUE)
    radius$taxonlist.y<-NULL #clean redundant labels
    radius$taxonlist.x<-NULL #clean redundant labels
    radius$radius <- rowMeans(cbind(radius$cv.radius.1, radius$cv.radius.2), na.rm=T)
    
    #generate a new dataset with standard errors and means per species
    radius_se<-merge(translation,radius, by="ID")[,c(1,2,5)]
    library(dplyr)
    radius_se <- radius_se %>%
      group_by(phylo) %>%
      summarize(
        radius.Mean = mean(radius),
        radius.Mean.log = log(mean(radius)),
        radius.StdErr = sd(radius) / sqrt(n()),
        radius.StdErr.log = sd(log(radius)) / sqrt(n())
      )
    radius_se<-as.data.frame(radius_se)
    #replace NAs with the mean value of SE for each column, per family
    radius_se$radius.StdErr.complete <- zoo::na.aggregate(radius_se$radius.StdErr)
    radius_se$radius.StdErr.log.complete <- zoo::na.aggregate(radius_se$radius.StdErr.log)
    
  }
  
  #mass
  {
    mass<-data.frame(ID, taxonlist, cv$vert_massing, cv$mass)
    colnames(mass)<-c("ID", "taxonlist", "vertnet_mass", "mass")
    #mass<-mass[mass$cv.mass.bprob > prob_filter,] #filter out rows without confidence
    #mass<-mass[complete.cases(mass),] #filter out rows that don't have data
    #mass$cv.mass<-ifelse(mass$cv.mass.bprob > prob_filter, cv$mass, NA) #mask those entries without prob_filter
    #mass$cv.mass.bprob<-NULL # clean the probability column
    
    #generate a new dataset with standard errors and means per species
    mass_se<-merge(translation, mass, by="ID")[,c(1,2,4,5)]
    library(dplyr)
    mass_se <- mass_se %>%
      group_by(phylo) %>%
      summarize(
        mass.vertnet.Mean = mean(vertnet_mass),
        mass.Mean = mean(mass),
        #sd = sd(cv.mass),
        mass.vertnet.log = log(mean(vertnet_mass)),
        mass.vertnet.StdErr = sd(vertnet_mass) / sqrt(n()),
        mass.vertnet.StdErr.log = sd(log(vertnet_mass)) / sqrt(n())
      )
    mass_se<-as.data.frame(mass_se)
    #replace NAs with the mean value of SE for each column, per family
    mass_se$mass.StdErr.complete <- zoo::na.aggregate(mass_se$mass.vertnet.StdErr)
    mass_se$mass.StdErr.log.complete <- zoo::na.aggregate(mass_se$mass.vertnet.StdErr.log)
    
  }
  
  #sex
  {
    sex<-data.frame(ID, taxonlist, cv$vert_sex)
    colnames(sex)<-c("ID", "taxonlist", "vertnet_sex")
    sex<-sex[c(sex$vertnet_sex == 'male' | sex$vertnet_sex == 'female'),]
    sex$vertnet_sex <- (as.factor(sex$vertnet_sex))
    
  }
  
}
#re assembling the dataset (warning is ok)
{
  #need to add individual masses
  dataframes<-list(tarsus, ulna, skull.1, second_digit, sclerotic_ring, metatarsus, keel.1, humerus, furcula.1, femur, carpometacarpus, radius, mass)
  
  #added specimen masses?
  dataframes_se<-list(tarsus_se, ulna_se, skull.1_se, second_digit_se, sclerotic_ring_se, metatarsus_se, keel.1_se, humerus_se, furcula.1_se, femur_se, carpometacarpus_se, radius_se, mass_se)
  
  # Define a function for a full outer join on a list of data frames using base R
  full_outer_join_list_baseR <- function(df_list, join_column = "ID") {
    if (length(df_list) < 2) {
      stop("At least two data frames are required for a full outer join.")
    }
    
    # Initialize the merged_df as the first data frame
    merged_df <- df_list[[1]]
    
    # Loop through the rest of the data frames and perform full outer joins
    for (i in 2:length(df_list)) {
      merged_df <- merge(merged_df, df_list[[i]], by = join_column, all = TRUE)
    }
    
    return(merged_df)
  }
  
  #merge and subset
  filtered_data <- full_outer_join_list_baseR(dataframes)
  #add back sex data for those specimens that make it this far
  filtered_data<-merge(x=sex, y=filtered_data, by = "ID", all.y=T)
  
  
  filtered_data <- filtered_data[,c("ID", "tarsus", "ulna", "cv.skull.1", "second_digit", "sclerotic_ring", "metatarsus", "cv.keel.1", "humerus", "cv.furcula.1", "femur", "carpometacarpus", "radius", "vertnet_mass", "mass", "vertnet_sex")] #"vertnet_sex"
  filtered_data <- merge(translation, filtered_data, by="ID", all=TRUE)
  
  #merge and subset (currently ignores standard error)
  filtered_data_se <- full_outer_join_list_baseR(dataframes_se, join_column="phylo")
  filtered_data_se.mass <- filtered_data_se[,c("phylo", "tarsus.Mean", "ulna.Mean", "skull.Mean", "second_digit.Mean", "sclerotic_ring.Mean", "metatarsus.Mean", "keel.Mean", "humerus.Mean", "furcula.Mean", "femur.Mean", "carpometacarpus.Mean", "radius.Mean", "mass.vertnet.Mean", "mass.Mean")]
  filtered_data_se <- filtered_data_se[,c("phylo", "tarsus.Mean", "ulna.Mean", "skull.Mean", "second_digit.Mean", "sclerotic_ring.Mean", "metatarsus.Mean", "keel.Mean", "humerus.Mean", "furcula.Mean", "femur.Mean", "carpometacarpus.Mean", "radius.Mean")]
}

#estimating descriptions of the dataset
calculate_na_percentage(filtered_data[,-c(1,2, 15, 16, 17)])
# Total non-NA values: 112769 
# Total NA values: 60259 
# Total values: 173028 
# [1] 34.82616

calculate_na_percentage(filtered_data_se[,-c(1)])
# Total non-NA values: 9768 
# Total NA values: 14916 
# Total values: 24684 
# [1] 60.42781

calculate_na_percentage(filtered_data_se.mass)
# Total non-NA values: 14383 
# Total NA values: 16472 
# Total values: 30855 
# [1] 53.38519

Skelevision_cleaned_specimens_v1.0<- filtered_data[,-c(15, 16, 17)] #excludes mass
Skelevision_cleaned_means_v1.0<- filtered_data_se

calculate_specimen_statistics(filtered_data)
# $variance
# [1] 142.6828
# $standard_deviation
# [1] 11.94499
# $range_min
# [1] 1
# $range_max
# [1] 161
# $mean
# [1] 7.009723
# $median
# [1] 3

#generate a diagnostic plot of Mass as a function of % missing data by bins
calculate_missing_percentage(filtered_data, plot=T)

#generating multivariate input data for phylopars
{
  phylopars.data.mv<-filtered_data
  #phylopars.data.mv$ID<-NULL
  phylopars.data.mv$species<-phylopars.data.mv$phylo
  phylopars.data.mv$phylo<-NULL
  
  # Move the last column to the first column
  phylopars.data.mv <- phylopars.data.mv[, c(ncol(phylopars.data.mv), 1:(ncol(phylopars.data.mv) - 1))]
  #log transform all data except for the first column
  phylopars.data.mv[, -c(1,2, 17)] <- log(phylopars.data.mv[, -c(1,2, 17)]) ### data is log transformed here
  # Count non-NA entries in each column
  non_na_count <- colSums(!is.na(phylopars.data.mv))
  
  #phylopars.data.mv.complete <- phylopars.data.mv[complete.cases(phylopars.data.mv),]
  #phylopars.data.mv.complete.nomass <- phylopars.data.mv[complete.cases(phylopars.data.mv[,-c(15,16, 17)]),]
  #str(phylopars.data.mv.complete.nomass) #up to ~600 individuals with complete data
  
  
}
{
  phylopars.data.mv.se<-filtered_data_se
  #phylopars.data.mv.se$ID<-NULL
  phylopars.data.mv.se$species<-phylopars.data.mv.se$phylo
  phylopars.data.mv.se$phylo<-NULL
  
  # Move the last column to the first column
  phylopars.data.mv.se <- phylopars.data.mv.se[, c(ncol(phylopars.data.mv.se), 1:(ncol(phylopars.data.mv.se) - 1))]
  #log transform all data except for the first column
  phylopars.data.mv.se[, -c(1)] <- log(phylopars.data.mv.se[, -c(1)]) ### data is log transformed here
  # Count non-NA entries in each column
  non_na_count <- colSums(!is.na(phylopars.data.mv.se))
  phylopars.data.mv.se.complete <- phylopars.data.mv.se[complete.cases(phylopars.data.mv.se),]

}
{
  phylopars.data.mv.se.mass<-filtered_data_se.mass[,-14] #exclude vertnet mass
  #phylopars.data.mv.se$ID<-NULL
  phylopars.data.mv.se.mass$species<-phylopars.data.mv.se.mass$phylo
  phylopars.data.mv.se.mass$phylo<-NULL
  
  # Move the last column to the first column
  phylopars.data.mv.se.mass <- phylopars.data.mv.se.mass[, c(ncol(phylopars.data.mv.se.mass), 1:(ncol(phylopars.data.mv.se.mass) - 1))]
  #log transform all data except for the first column
  phylopars.data.mv.se.mass[, -c(1)] <- log(phylopars.data.mv.se.mass[, -c(1)]) ### data is log transformed here
  # Count non-NA entries in each column
  non_na_count <- colSums(!is.na(phylopars.data.mv.se.mass))
  phylopars.data.mv.se.mass.complete <- phylopars.data.mv.se.mass[complete.cases(phylopars.data.mv.se.mass),]

  
}

##mondo imputation without mass
#phylopars.corr.mvBM <- phylopars(phylopars.data.mv[,-c(2, 15, 16, 17)], tree = tree.pruned, phylo_correlated = T, pheno_correlated = T)
#saveRDS(object = phylopars.corr.mvBM, file="phylopars.corr.mvBM.RDS") #need to correct the name after it runs
phylopars.corr.mvBM <- readRDS('phylopars.corr.mvBM.RDS')

##mondo imputation with mass (specimen masses)
#phylopars.corr.mass.mvBM <- phylopars(phylopars.data.mv[,-c(2, 16, 17)], tree = tree.pruned, phylo_correlated = T, pheno_correlated = T)
#saveRDS(object = phylopars.corr.mass.mvBM, file="phylopars.corr.mass.mvBM.RDS") #need to correct the name after it runs
phylopars.corr.mass.mvBM <- readRDS('phylopars.corr.mass.mvBM.RDS')

###creating data export objects
{
###creating objects to share with brian and charlotte
# Creating the mask
# mask <- is.na(phylopars.data.mv[,-c(2, 15, 16)])
# tmp<-phylopars.corr.mvBM$ind_recon
# tmp[mask] <- NA
# check_full_equivalence_with_tolerance(phylopars.data.mv[,-c(2, 15, 16)], tmp)
# 
# filtered_dataset<- cbind(ID = phylopars.data.mv[,c(2)], phylopars.data.mv[,-c(2, 15, 16)])
# ### DONE
# 
# tmp<-phylopars.corr.mvBM$ind_recon
# tmp[!mask] <- NA #inverse mask to mask the measured data
# imputed_dataset<- cbind(ID = phylopars.data.mv[,2], tmp)
# ### DONE
# 
# completed_dataset <- cbind(ID = phylopars.data.mv[,c(2)], species= phylopars.corr.mvBM$ind_recon[, 1], phylopars.corr.mvBM$ind_recon[,-1]) 
# ### DONE
# 
# str(imputed_dataset)
# str(filtered_dataset)
# str(completed_dataset)


#function that generates masks and the final dataset objects
individual_data.nomass<-create_datasets(ID = as.data.frame(phylopars.data.mv[,c(2)]), original_data = phylopars.data.mv[,-c(2, 15, 16, 17)], reconstructed_data = phylopars.corr.mvBM$ind_recon)
individual_data.mass<-create_datasets(ID = as.data.frame(phylopars.data.mv[,c(2)]), original_data = phylopars.data.mv[,-c(2, 16, 17)], reconstructed_data = phylopars.corr.mass.mvBM$ind_recon)
#saveRDS(individual_data.nomass, file="individual_data.mvBM.nomass.12.22.23.RDS")
#saveRDS(individual_data.mass, file="individual_data.mvBM.mass.12.22.23.RDS")
individual_data.nomass<- readRDS('individual_data.mvBM.nomass.12.22.23.RDS')

#write CSV files for brian, in linear space
#write.csv(lapply(individual_data.nomass$filtered_dataset, function(x) if(is.numeric(x)) exp(x) else x), file = "12.22.23_specimens_filtered_dataset.csv", quote=F, row.names=F)
#write.csv(lapply(individual_data.nomass$completed_dataset, function(x) if(is.numeric(x)) exp(x) else x), file = "12.22.23_specimens_completed_dataset.csv", quote=F, row.names=F)

calculate_r_squared_matching_columns(individual_data.nomass$imputed_dataset[,-c(1,2)], individual_data.mass$imputed_dataset[,-c(1,2)])
calculate_correlations_matching_columns(individual_data.nomass$imputed_dataset[,-c(1,2)], individual_data.mass$imputed_dataset[,-c(1,2)])
#including mass makes essentially no difference in the imputed datasets

## construct the species datasets
species_data.nomass<-construct_CI_df(recon_values = phylopars.corr.mvBM$anc_recon[tree.pruned$tip.label,], variances =  phylopars.corr.mvBM$anc_var[tree.pruned$tip.label,])
species_data.mass<-construct_CI_df(recon_values = phylopars.corr.mass.mvBM$anc_recon[tree.pruned$tip.label,], variances =  phylopars.corr.mass.mvBM$anc_var[tree.pruned$tip.label,])
#saveRDS(species_data.nomass, file="species_data.mvBM.nomass.12.22.23.RDS")
#saveRDS(species_data.mass, file="species_data.mvBM.mass.12.22.23.RDS")
species_data.nomass<-readRDS('species_data.mvBM.nomass.12.22.23.RDS')

#write CSV files for brian, in linear space
#write.csv(lapply(species_data.nomass$Recon, function(x) if(is.numeric(x)) exp(x) else x), file = "12.22.23_species_completed_dataset.csv", quote=F, row.names=F)
#write.csv(lapply(species_data.nomass$SE, function(x) if(is.numeric(x)) exp(x) else x), file = "12.22.23_species_completed_SE.csv", quote=F, row.names=F)

calculate_na_percentage(individual_data.nomass$filtered_dataset[grep(x=individual_data.nomass$filtered_dataset$species, pattern='Vidua'), ])
#hist(individual_data.nomass$filtered_dataset[grep(x=individual_data.nomass$filtered_dataset$species, pattern='Vidua'), ][,3])
#individual_data.nomass$filtered_dataset[grep(x=individual_data.nomass$filtered_dataset$species, pattern='Vidua'), ][71,]

}

### setting up multi-masking to estimate RSME and p-bias
{
##analyses without mass
{
    {
      ### 10%
      {
        set.seed(seed=2) #for reproducibility
        mask_set.10<-mask_and_track_percentage(data = phylopars.data.mv[,-c(2, 15, 16, 17)], mask_percentage = 10, repetitions = 10, exclude_cols = 'species')
        #setting up replicates (10%)
        {
        # Extract just the masked data from each replicate
        masked_data_list <- lapply(mask_set.10, function(x) x$masked_data)

        plan(multisession, workers=10)
        # Apply imputation in parallel
        imputation_results.10.2 <- future_lapply(masked_data_list, function(masked_data) {
          phylopars(masked_data, tree = tree.pruned, phylo_correlated = TRUE, pheno_correlated = TRUE)
        }, future.seed=T)

        # Reset future plan to sequential for safety
        plan(sequential)
        }
       
      }
      {
        # Extract masked positions from mask_set
        mask_set.10_masked_positions <- lapply(mask_set.10, function(x) x$masked_positions)
        # Extract imputed datasets from imputation_results
        imputation_results.10.imputed <- lapply(imputation_results.10, function(x) x$ind_recon)
        
        # Compute the metrics
        metrics.10 <- compute_metrics(phylopars.data.mv[,-c(2, 15, 16, 17)], imputation_results.10.imputed, mask_set.10_masked_positions)
        metrics.10.mean.rmse <- mean(metrics.10$rmse) 
        
        #plot the results
        plot_metrics(metrics.10, error_type = 'ci')
      }
      
      ### 20%
      {
        set.seed(seed=1) #for reproducibility
        mask_set.20<-mask_and_track_percentage(data = phylopars.data.mv[,-c(2, 15, 16, 17)], mask_percentage = 20, repetitions = 10, exclude_cols = 'species')
        setting up replicates (20%)
        {
          # Extract just the masked data from each replicate
          masked_data_list <- lapply(mask_set.20, function(x) x$masked_data)

          plan(multisession, workers=10)
          # Apply imputation in parallel
          imputation_results.20 <- future_lapply(masked_data_list, function(masked_data) {
            phylopars(masked_data, tree = tree.pruned, phylo_correlated = TRUE, pheno_correlated = TRUE)
          }, future.seed=T)

          # Reset future plan to sequential for safety
          plan(sequential)
        }
       
      }
      {
        # Extract masked positions from mask_set
        mask_set.20_masked_positions <- lapply(mask_set.20, function(x) x$masked_positions)
        # Extract imputed datasets from imputation_results
        imputation_results.20.imputed <- lapply(imputation_results.20, function(x) x$ind_recon)
        
        # Compute the metrics
        metrics.20 <- compute_metrics(phylopars.data.mv[,-c(2, 15, 16, 17)], imputation_results.20.imputed, mask_set.20_masked_positions)
        metrics.20.mean.rmse <- mean(metrics.20$rmse) 
        
        #plot the results
        plot_metrics(metrics.20, error_type = 'ci')
      }
      
      ### 30%
      {
        set.seed(seed=1) #for reproducibility
        mask_set.30<-mask_and_track_percentage(data = phylopars.data.mv[,-c(2, 15, 16, 17)], mask_percentage = 30, repetitions = 10, exclude_cols = 'species')
        #setting up replicates (30%)
        {
          # Extract just the masked data from each replicate
          masked_data_list <- lapply(mask_set.30, function(x) x$masked_data)

          plan(multisession, workers=10)
          # Apply imputation in parallel
          imputation_results.30 <- future_lapply(masked_data_list, function(masked_data) {
            phylopars(masked_data, tree = tree.pruned, phylo_correlated = TRUE, pheno_correlated = TRUE)
          }, future.seed=T)

          # Reset future plan to sequential for safety
          plan(sequential)
        }
        
      }
      {
        # Extract masked positions from mask_set
        mask_set.30_masked_positions <- lapply(mask_set.30, function(x) x$masked_positions)
        # Extract imputed datasets from imputation_results
        imputation_results.30.imputed <- lapply(imputation_results.30, function(x) x$ind_recon)
        
        # Compute the metrics
        metrics.30 <- compute_metrics(phylopars.data.mv[,-c(2, 15, 16, 17)], imputation_results.30.imputed, mask_set.30_masked_positions)
        metrics.30.mean.rmse <- mean(metrics.30$rmse) 
        
        #plot the results
        plot_metrics(metrics.30, error_type = 'ci')
      }
      
      ### 40%
      {
        set.seed(seed=1) #for reproducibility
        mask_set.40<-mask_and_track_percentage(data = phylopars.data.mv[,-c(2, 15, 16, 17)], mask_percentage = 40, repetitions = 10, exclude_cols = 'species')
        #setting up replicates (40%)
        {
          # Extract just the masked data from each replicate
          masked_data_list <- lapply(mask_set.40, function(x) x$masked_data)

          plan(multisession, workers=10)
          # Apply imputation in parallel
          imputation_results.40 <- future_lapply(masked_data_list, function(masked_data) {
            phylopars(masked_data, tree = tree.pruned, phylo_correlated = TRUE, pheno_correlated = TRUE)
          }, future.seed=T)

          # Reset future plan to sequential for safety
          plan(sequential)
        }
        
      }
      {
        # Extract masked positions from mask_set
        mask_set.40_masked_positions <- lapply(mask_set.40, function(x) x$masked_positions)
        # Extract imputed datasets from imputation_results
        imputation_results.40.imputed <- lapply(imputation_results.40, function(x) x$ind_recon)
        
        # Compute the metrics
        metrics.40 <- compute_metrics(phylopars.data.mv[,-c(2, 15, 16, 17)], imputation_results.40.imputed, mask_set.40_masked_positions)
        metrics.40.mean.rmse <- mean(metrics.40$rmse) 
        
        #plot the results
        plot_metrics(metrics.40, error_type = 'ci')
      }
      
      ### 50%
      {
        set.seed(seed=1) #for reproducibility
        mask_set.50<-mask_and_track_percentage(data = phylopars.data.mv[,-c(2, 15, 16, 17)], mask_percentage = 50, repetitions = 10, exclude_cols = 'species')
        #setting up replicates (50%)
        {
          # Extract just the masked data from each replicate
          masked_data_list <- lapply(mask_set.50, function(x) x$masked_data)

          plan(multisession, workers=10)
          # Apply imputation in parallel
          imputation_results.50 <- future_lapply(masked_data_list, function(masked_data) {
            phylopars(masked_data, tree = tree.pruned, phylo_correlated = TRUE, pheno_correlated = TRUE)
          }, future.seed=T)

          # Reset future plan to sequential for safety
          plan(sequential)
        }
       
      }
      {
        # Extract masked positions from mask_set
        mask_set.50_masked_positions <- lapply(mask_set.50, function(x) x$masked_positions)
        # Extract imputed datasets from imputation_results
        imputation_results.50.imputed <- lapply(imputation_results.50, function(x) x$ind_recon)
        
        # Compute the metrics
        metrics.50 <- compute_metrics(phylopars.data.mv[,-c(2, 15, 16, 17)], imputation_results.50.imputed, mask_set.50_masked_positions)
        metrics.50.mean.rmse <- mean(metrics.50$rmse) 
        
        #plot the results
        plot_metrics(metrics.50, error_type = 'ci')
      }
      
      ##60%
      {
        set.seed(seed=1) #for reproducibility
        mask_set.60<-mask_and_track_percentage(data = phylopars.data.mv[,-c(2, 15, 16, 17)], mask_percentage = 60, repetitions = 20, exclude_cols = 'species')

        {
          # Extract just the masked data from each replicate
          masked_data_list <- lapply(mask_set.60, function(x) x$masked_data)

          plan(multisession, workers=20)
          # Apply imputation in parallel with error handling
          imputation_results.60 <- future_lapply(masked_data_list, function(masked_data) {
            tryCatch({
              phylopars(masked_data, tree = tree.pruned, phylo_correlated = TRUE, pheno_correlated = TRUE)
            }, error = function(e) {
              NULL  # Return NULL if an error occurs
            })
          }, future.seed = TRUE)

          # Reset future plan to sequential for safety
          plan(sequential)
        }
       
        imputation_results.60<-imputation_results.60[!unlist(lapply(imputation_results.60, is.null))] #filter out nulls
        length(imputation_results.60) #12 replicates worked.
        imputation_results.60 <- sample(imputation_results.60, 10) #sample 10 of them
      }
      {
        # Extract masked positions from mask_set
        mask_set.60_masked_positions <- lapply(mask_set.60, function(x) x$masked_positions)
        # Extract imputed datasets from imputation_results
        imputation_results.60.imputed <- lapply(imputation_results.60, function(x) x$ind_recon)
        
        # Compute the metrics
        metrics.60 <- compute_metrics(phylopars.data.mv[,-c(2, 15, 16, 17)], imputation_results.60.imputed, mask_set.60_masked_positions)
        metrics.60.mean.rmse <- mean(metrics.60$rmse)
        #0.07946772
        
        #plot the results
        plot_metrics(metrics.60, error_type = 'ci')
      }
      
      ##70%
      {
        set.seed(seed=1) #for reproducibility
        mask_set.70<-mask_and_track_percentage(data = phylopars.data.mv[,-c(2, 15, 16, 17)], mask_percentage = 70, repetitions = 30, exclude_cols = 'species')

        {
          # Extract just the masked data from each replicate
          masked_data_list <- lapply(mask_set.70, function(x) x$masked_data)

          plan(multisession, workers=30)
          # Apply imputation in parallel with error handling
          imputation_results.70 <- future_lapply(masked_data_list, function(masked_data) {
            tryCatch({
              phylopars(masked_data, tree = tree.pruned, phylo_correlated = TRUE, pheno_correlated = TRUE)
            }, error = function(e) {
              NULL  # Return NULL if an error occurs
            })
          }, future.seed = TRUE)

          # Reset future plan to sequential for safety
          plan(sequential)
        }
       
        imputation_results.70<-imputation_results.70[!unlist(lapply(imputation_results.70, is.null))] #filter out nulls
        length(imputation_results.70) #12 replicates worked.
        imputation_results.70 <- sample(imputation_results.70, 10) #sample 10 of them
      }
      {
        # Extract masked positions from mask_set
        mask_set.70_masked_positions <- lapply(mask_set.70, function(x) x$masked_positions)
        # Extract imputed datasets from imputation_results
        imputation_results.70.imputed <- lapply(imputation_results.70, function(x) x$ind_recon)
        
        # Compute the metrics
        metrics.70 <- compute_metrics(phylopars.data.mv[,-c(2, 15, 16, 17)], imputation_results.70.imputed, mask_set.70_masked_positions)
        metrics.70.mean.rmse <- mean(metrics.70$rmse)
        #0.07946772
        
        #plot the results
        plot_metrics(metrics.70, error_type = 'ci')
      }
      
      ##80%
      {
        set.seed(seed=1) #for reproducibility
        mask_set.80<-mask_and_track_percentage(data = phylopars.data.mv[,-c(2, 15, 16, 17)], mask_percentage = 80, repetitions = 30, exclude_cols = 'species')

        {
          # Extract just the masked data from each replicate
          masked_data_list <- lapply(mask_set.80, function(x) x$masked_data)

          plan(multisession, workers=30)
          # Apply imputation in parallel with error handling
          imputation_results.80 <- future_lapply(masked_data_list, function(masked_data) {
            tryCatch({
              phylopars(masked_data, tree = tree.pruned, phylo_correlated = TRUE, pheno_correlated = TRUE)
            }, error = function(e) {
              NULL  # Return NULL if an error occurs
            })
          }, future.seed = TRUE)

          # Reset future plan to sequential for safety
          plan(sequential)
        }
       
        imputation_results.80<-imputation_results.80[!unlist(lapply(imputation_results.80, is.null))] #filter out nulls
        length(imputation_results.80) #12 replicates worked.
        imputation_results.80 <- sample(imputation_results.80, 10) #sample 10 of them
      }
      {
        # Extract masked positions from mask_set
        mask_set.80_masked_positions <- lapply(mask_set.80, function(x) x$masked_positions)
        # Extract imputed datasets from imputation_results
        imputation_results.80.imputed <- lapply(imputation_results.80, function(x) x$ind_recon)
        
        # Compute the metrics
        metrics.80 <- compute_metrics(phylopars.data.mv[,-c(2, 15, 16, 17)], imputation_results.80.imputed, mask_set.80_masked_positions)
        metrics.80.mean.rmse <- mean(metrics.80$rmse)
        #0.07946772
        
        #plot the results
        plot_metrics(metrics.80, error_type = 'ci')
      }
      
      ##90%
      {
        set.seed(seed=1) #for reproducibility
        mask_set.90<-mask_and_track_percentage(data = phylopars.data.mv[,-c(2, 15, 16, 17)], mask_percentage = 90, repetitions = 30, exclude_cols = 'species')

        {
          # Extract just the masked data from each replicate
          masked_data_list <- lapply(mask_set.90, function(x) x$masked_data)

          plan(multisession, workers=30)
          # Apply imputation in parallel with error handling
          imputation_results.90 <- future_lapply(masked_data_list, function(masked_data) {
            tryCatch({
              phylopars(masked_data, tree = tree.pruned, phylo_correlated = TRUE, pheno_correlated = TRUE)
            }, error = function(e) {
              NULL  # Return NULL if an error occurs
            })
          }, future.seed = TRUE)

          # Reset future plan to sequential for safety
          plan(sequential)
        }
       
        imputation_results.90<-imputation_results.90[!unlist(lapply(imputation_results.90, is.null))] #filter out nulls
        length(imputation_results.90) #12 replicates worked.
        imputation_results.90 <- sample(imputation_results.90, 10) #sample 10 of them
      }
      {
        # Extract masked positions from mask_set
        mask_set.90_masked_positions <- lapply(mask_set.90, function(x) x$masked_positions)
        # Extract imputed datasets from imputation_results
        imputation_results.90.imputed <- lapply(imputation_results.90, function(x) x$ind_recon)
        
        # Compute the metrics
        metrics.90 <- compute_metrics(phylopars.data.mv[,-c(2, 15, 16, 17)], imputation_results.90.imputed, mask_set.90_masked_positions)
        metrics.90.mean.rmse <- mean(metrics.90$rmse)
        #0.07946772
        
        #plot the results
        plot_metrics(metrics.90, error_type = 'ci')
      }
    }
    
    rmse.dat <- rbind(metrics.10$rmse, metrics.20$rmse, metrics.30$rmse, metrics.40$rmse, metrics.50$rmse, metrics.60$rmse, metrics.70$rmse, metrics.80$rmse, metrics.90$rmse)
    pbias.dat <- rbind(metrics.10$p_bias, metrics.20$p_bias, metrics.30$p_bias, metrics.40$p_bias, metrics.50$p_bias, metrics.60$p_bias, metrics.70$p_bias, metrics.80$p_bias, metrics.90$p_bias)
    x_values <- c(10, 20, 30, 40, 50, 60, 70, 80, 90)
    #plot_metrics(rmse.dat.list)
    #plot_rmse(rmse_dat = rmse.dat, x_values)
    plot_rmse_and_pbias(rmse_dat =  rmse.dat, pbias_dat = pbias.dat, x_values, xlim_rmse = c(10,50), xlim_pbias =  c(10,50))
    
    plot_rmse_and_pbias(rmse_dat =  rmse.dat, pbias_dat = pbias.dat, x_values)
    
    
    #export for brian
    saveRDS(rmse.dat, file='rmse.dat.RDS')
    saveRDS(pbias.dat, file='pbias.dat.RDS')
    saveRDS(x_values, file='x_values.RDS')
    
    
  }

## analyses with mass
{
    {
      ### 10%
      {
        set.seed(seed=1) #for reproducibility
        mask_set.10<-mask_and_track_percentage(data = phylopars.data.mv[,-c(2, 16, 17)], mask_percentage = 10, repetitions = 10, exclude_cols = 'species')
        #setting up replicates (10%)
        {
        # Extract just the masked data from each replicate
        masked_data_list <- lapply(mask_set.10, function(x) x$masked_data)

        plan(multisession, workers=10)
        # Apply imputation in parallel
        imputation_results.mass.10 <- future_lapply(masked_data_list, function(masked_data) {
          phylopars(masked_data, tree = tree.pruned, phylo_correlated = TRUE, pheno_correlated = TRUE)
        }, future.seed=T)

        # Reset future plan to sequential for safety
        plan(sequential)
        }
       
      }
      {
        # Extract masked positions from mask_set
        mask_set.10_masked_positions <- lapply(mask_set.10, function(x) x$masked_positions)
        # Extract imputed datasets from imputation_results
        imputation_results.mass.10.imputed <- lapply(imputation_results.mass.10, function(x) x$ind_recon)
        
        # Compute the metrics
        metrics.10 <- compute_metrics(phylopars.data.mv[,-c(2, 16, 17)], imputation_results.mass.10.imputed, mask_set.10_masked_positions)
        metrics.10.mean.rmse <- mean(metrics.10$rmse) #0.07600862
        
        #plot the results
        plot_metrics(metrics.10, error_type = 'ci')
      }
      
      ### 20%
      {
        set.seed(seed=1) #for reproducibility
        mask_set.20<-mask_and_track_percentage(data = phylopars.data.mv[,-c(2, 16, 17)], mask_percentage = 20, repetitions = 10, exclude_cols = 'species')
        #setting up replicates (20%)
        {
          # Extract just the masked data from each replicate
          masked_data_list <- lapply(mask_set.20, function(x) x$masked_data)

          plan(multisession, workers=10)
          # Apply imputation in parallel
          imputation_results.mass.20 <- future_lapply(masked_data_list, function(masked_data) {
            phylopars(masked_data, tree = tree.pruned, phylo_correlated = TRUE, pheno_correlated = TRUE)
          }, future.seed=T)

          # Reset future plan to sequential for safety
          plan(sequential)
        }
        
      }
      {
        # Extract masked positions from mask_set
        mask_set.20_masked_positions <- lapply(mask_set.20, function(x) x$masked_positions)
        # Extract imputed datasets from imputation_results
        imputation_results.mass.20.imputed <- lapply(imputation_results.mass.20, function(x) x$ind_recon)
        
        # Compute the metrics
        metrics.20 <- compute_metrics(phylopars.data.mv[,-c(2, 16, 17)], imputation_results.mass.20.imputed, mask_set.20_masked_positions)
        metrics.20.mean.rmse <- mean(metrics.20$rmse) #0.07495493
        
        #plot the results
        plot_metrics(metrics.20, error_type = 'ci')
      }
      
      ### 30%
      {
        set.seed(seed=1) #for reproducibility
        mask_set.30<-mask_and_track_percentage(data = phylopars.data.mv[,-c(2, 16, 17)], mask_percentage = 30, repetitions = 10, exclude_cols = 'species')
        # #setting up replicates (30%)
        {
          # Extract just the masked data from each replicate
          masked_data_list <- lapply(mask_set.30, function(x) x$masked_data)

          plan(multisession, workers=10)
          # Apply imputation in parallel
          imputation_results.mass.30 <- future_lapply(masked_data_list, function(masked_data) {
            phylopars(masked_data, tree = tree.pruned, phylo_correlated = TRUE, pheno_correlated = TRUE)
          }, future.seed=T)

          # Reset future plan to sequential for safety
          plan(sequential)
        }
       
      }
      {
        # Extract masked positions from mask_set
        mask_set.30_masked_positions <- lapply(mask_set.30, function(x) x$masked_positions)
        # Extract imputed datasets from imputation_results
        imputation_results.mass.30.imputed <- lapply(imputation_results.mass.30, function(x) x$ind_recon)
        
        # Compute the metrics
        metrics.30 <- compute_metrics(phylopars.data.mv[,-c(2, 16, 17)], imputation_results.mass.30.imputed, mask_set.30_masked_positions)
        metrics.30.mean.rmse <- mean(metrics.30$rmse) #0.07685427
        
        #plot the results
        plot_metrics(metrics.30, error_type = 'ci')
      }
      
      ### 40%
      {
        set.seed(seed=1) #for reproducibility
        mask_set.40<-mask_and_track_percentage(data = phylopars.data.mv[,-c(2, 16, 17)], mask_percentage = 40, repetitions = 10, exclude_cols = 'species')
        # #setting up replicates (40%)
        {
          # Extract just the masked data from each replicate
          masked_data_list <- lapply(mask_set.40, function(x) x$masked_data)

          plan(multisession, workers=10)
          # Apply imputation in parallel
          imputation_results.mass.40 <- future_lapply(masked_data_list, function(masked_data) {
            phylopars(masked_data, tree = tree.pruned, phylo_correlated = TRUE, pheno_correlated = TRUE)
          }, future.seed=T)

          # Reset future plan to sequential for safety
          plan(sequential)
        }
        
      }
      {
        # Extract masked positions from mask_set
        mask_set.40_masked_positions <- lapply(mask_set.40, function(x) x$masked_positions)
        # Extract imputed datasets from imputation_results
        imputation_results.mass.40.imputed <- lapply(imputation_results.mass.40, function(x) x$ind_recon)
        
        # Compute the metrics
        metrics.40 <- compute_metrics(phylopars.data.mv[,-c(2, 16, 17)], imputation_results.mass.40.imputed, mask_set.40_masked_positions)
        metrics.40.mean.rmse <- mean(metrics.40$rmse) #0.07752544
        
        #plot the results
        plot_metrics(metrics.40, error_type = 'ci')
      }
      
      ### 50%
      {
        set.seed(seed=1) #for reproducibility
        mask_set.50<-mask_and_track_percentage(data = phylopars.data.mv[,-c(2, 16, 17)], mask_percentage = 50, repetitions = 10, exclude_cols = 'species')
        # #setting up replicates (50%)
        {
          # Extract just the masked data from each replicate
          masked_data_list <- lapply(mask_set.50, function(x) x$masked_data)

          plan(multisession, workers=10)
          # Apply imputation in parallel
          imputation_results.mass.50 <- future_lapply(masked_data_list, function(masked_data) {
            phylopars(masked_data, tree = tree.pruned, phylo_correlated = TRUE, pheno_correlated = TRUE)
          }, future.seed=T)

          # Reset future plan to sequential for safety
          plan(sequential)
        }
       
      }
      {
        # Extract masked positions from mask_set
        mask_set.50_masked_positions <- lapply(mask_set.50, function(x) x$masked_positions)
        # Extract imputed datasets from imputation_results
        imputation_results.mass.50.imputed <- lapply(imputation_results.mass.50, function(x) x$ind_recon)
        
        # Compute the metrics
        metrics.50 <- compute_metrics(phylopars.data.mv[,-c(2, 16, 17)], imputation_results.mass.50.imputed, mask_set.50_masked_positions)
        metrics.50.mean.rmse <- mean(metrics.50$rmse) #0.07946772
        
        #plot the results
        plot_metrics(metrics.50, error_type = 'ci')
      }
      
      ##60%
      {
        set.seed(seed=1) #for reproducibility
        mask_set.60<-mask_and_track_percentage(data = phylopars.data.mv[,-c(2, 16, 17)], mask_percentage = 60, repetitions = 20, exclude_cols = 'species')

              {
                # Extract just the masked data from each replicate
                masked_data_list <- lapply(mask_set.60, function(x) x$masked_data)

                plan(multisession, workers=20)
                # Apply imputation in parallel with error handling
                imputation_results.mass.60 <- future_lapply(masked_data_list, function(masked_data) {
                  tryCatch({
                    phylopars(masked_data, tree = tree.pruned, phylo_correlated = TRUE, pheno_correlated = TRUE)
                  }, error = function(e) {
                    NULL  # Return NULL if an error occurs
                  })
                }, future.seed = TRUE)

                # Reset future plan to sequential for safety
                plan(sequential)
              }
       
        imputation_results.mass.60<-imputation_results.mass.60[!unlist(lapply(imputation_results.mass.60, is.null))] #filter out nulls
        length(imputation_results.mass.60) #12 replicates worked.
        imputation_results.mass.60 <- sample(imputation_results.mass.60, 10) #sample 10 of them
      }
      {
        # Extract masked positions from mask_set
        mask_set.60_masked_positions <- lapply(mask_set.60, function(x) x$masked_positions)
        # Extract imputed datasets from imputation_results
        imputation_results.mass.60.imputed <- lapply(imputation_results.mass.60, function(x) x$ind_recon)
        
        # Compute the metrics
        metrics.60 <- compute_metrics(phylopars.data.mv[,-c(2, 16, 17)], imputation_results.mass.60.imputed, mask_set.60_masked_positions)
        metrics.60.mean.rmse <- mean(metrics.60$rmse) #0.06357431
        
        #plot the results
        plot_metrics(metrics.60, error_type = 'ci')
      }
      
      ##70%
      {
        set.seed(seed=1) #for reproducibility
        mask_set.70<-mask_and_track_percentage(data = phylopars.data.mv[,-c(2, 16, 17)], mask_percentage = 70, repetitions = 30, exclude_cols = 'species')

              {
                # Extract just the masked data from each replicate
                masked_data_list <- lapply(mask_set.70, function(x) x$masked_data)

                plan(multisession, workers=30)
                # Apply imputation in parallel with error handling
                imputation_results.mass.70 <- future_lapply(masked_data_list, function(masked_data) {
                  tryCatch({
                    phylopars(masked_data, tree = tree.pruned, phylo_correlated = TRUE, pheno_correlated = TRUE)
                  }, error = function(e) {
                    NULL  # Return NULL if an error occurs
                  })
                }, future.seed = TRUE)

                # Reset future plan to sequential for safety
                plan(sequential)
              }
        
        imputation_results.mass.70<-imputation_results.mass.70[!unlist(lapply(imputation_results.mass.70, is.null))] #filter out nulls
        length(imputation_results.mass.70) #12 replicates worked.
        imputation_results.mass.70 <- sample(imputation_results.mass.70, 10) #sample 10 of them
      }
      {
        # Extract masked positions from mask_set
        mask_set.70_masked_positions <- lapply(mask_set.70, function(x) x$masked_positions)
        # Extract imputed datasets from imputation_results
        imputation_results.mass.70.imputed <- lapply(imputation_results.mass.70, function(x) x$ind_recon)
        
        # Compute the metrics
        metrics.70 <- compute_metrics(phylopars.data.mv[,-c(2, 16, 17)], imputation_results.mass.70.imputed, mask_set.70_masked_positions)
        metrics.70.mean.rmse <- mean(metrics.70$rmse) #0.07282172
        
        #plot the results
        plot_metrics(metrics.70, error_type = 'ci')
      }
      
      ##80%
      {
        set.seed(seed=1) #for reproducibility
        mask_set.80<-mask_and_track_percentage(data = phylopars.data.mv[,-c(2, 16, 17)], mask_percentage = 80, repetitions = 30, exclude_cols = 'species')

              {
                # Extract just the masked data from each replicate
                masked_data_list <- lapply(mask_set.80, function(x) x$masked_data)

                plan(multisession, workers=30)
                # Apply imputation in parallel with error handling
                imputation_results.mass.80 <- future_lapply(masked_data_list, function(masked_data) {
                  tryCatch({
                    phylopars(masked_data, tree = tree.pruned, phylo_correlated = TRUE, pheno_correlated = TRUE)
                  }, error = function(e) {
                    NULL  # Return NULL if an error occurs
                  })
                }, future.seed = TRUE)

                # Reset future plan to sequential for safety
                plan(sequential)
              }
        
        imputation_results.mass.80<-imputation_results.mass.80[!unlist(lapply(imputation_results.mass.80, is.null))] #filter out nulls
        length(imputation_results.mass.80) #12 replicates worked.
        imputation_results.mass.80 <- sample(imputation_results.mass.80, 10) #sample 10 of them
      }
      {
        # Extract masked positions from mask_set
        mask_set.80_masked_positions <- lapply(mask_set.80, function(x) x$masked_positions)
        # Extract imputed datasets from imputation_results
        imputation_results.mass.80.imputed <- lapply(imputation_results.mass.80, function(x) x$ind_recon)
        
        # Compute the metrics
        metrics.80 <- compute_metrics(phylopars.data.mv[,-c(2, 16, 17)], imputation_results.mass.80.imputed, mask_set.80_masked_positions)
        metrics.80.mean.rmse <- mean(metrics.80$rmse) #0.08358489
        
        #plot the results
        plot_metrics(metrics.80, error_type = 'ci')
      }
      
      ##90%
      {
        set.seed(seed=1) #for reproducibility
        mask_set.90<-mask_and_track_percentage(data = phylopars.data.mv[,-c(2, 16, 17)], mask_percentage = 90, repetitions = 30, exclude_cols = 'species')

              {
                # Extract just the masked data from each replicate
                masked_data_list <- lapply(mask_set.90, function(x) x$masked_data)

                plan(multisession, workers=30)
                # Apply imputation in parallel with error handling
                imputation_results.mass.90 <- future_lapply(masked_data_list, function(masked_data) {
                  tryCatch({
                    phylopars(masked_data, tree = tree.pruned, phylo_correlated = TRUE, pheno_correlated = TRUE)
                  }, error = function(e) {
                    NULL  # Return NULL if an error occurs
                  })
                }, future.seed = TRUE)

                # Reset future plan to sequential for safety
                plan(sequential)
              }
        
        imputation_results.mass.90<-imputation_results.mass.90[!unlist(lapply(imputation_results.mass.90, is.null))] #filter out nulls
        length(imputation_results.mass.90) #12 replicates worked.
        imputation_results.mass.90 <- sample(imputation_results.mass.90, 10) #sample 10 of them
      }
      {
        # Extract masked positions from mask_set
        mask_set.90_masked_positions <- lapply(mask_set.90, function(x) x$masked_positions)
        # Extract imputed datasets from imputation_results
        imputation_results.mass.90.imputed <- lapply(imputation_results.mass.90, function(x) x$ind_recon)
        
        # Compute the metrics
        metrics.90 <- compute_metrics(phylopars.data.mv[,-c(2, 16, 17)], imputation_results.mass.90.imputed, mask_set.90_masked_positions)
        metrics.90.mean.rmse <- mean(metrics.90$rmse) #0.1057793
        
        #plot the results
        plot_metrics(metrics.90, error_type = 'ci')
      }
    }
    
    rmse.dat <- rbind(metrics.10$rmse, metrics.20$rmse, metrics.30$rmse, metrics.40$rmse, metrics.50$rmse, metrics.60$rmse, metrics.70$rmse, metrics.80$rmse, metrics.90$rmse)
    pbias.dat <- rbind(metrics.10$p_bias, metrics.20$p_bias, metrics.30$p_bias, metrics.40$p_bias, metrics.50$p_bias, metrics.60$p_bias, metrics.70$p_bias, metrics.80$p_bias, metrics.90$p_bias)
    x_values <- c(10, 20, 30, 40, 50, 60, 70, 80, 90)
    #plot_metrics(rmse.dat.list)
    #plot_rmse(rmse_dat = rmse.dat, x_values)
    plot_rmse_and_pbias(rmse_dat =  rmse.dat, pbias_dat = pbias.dat, x_values)
    
  }

}




