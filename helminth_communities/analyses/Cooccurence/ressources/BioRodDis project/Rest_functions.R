# FUNCTIONS

## Function to create a matrix from a host-pathogen table created with for example
## dat2 %>% dplyr::select(Species, Site_Code) %>% table() -> m1

matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[,1]
  m
}

## FUNCTION TO CALCULATE NETWORK STATISTICS FOR A GIVEN SITE


net.site <- function(site){
  
  dat3[dat3$Site_Code == site,] -> gegevens
  
  gegevens %>% group_by(Species) %>% count(Bartonella) %>% filter(Bartonella == 1) %>% dplyr::select(Species, n) %>% rename(Bartonella = n) ->d1
  gegevens %>% group_by(Species) %>% count(Mycoplasma_haemo) %>% filter(Mycoplasma_haemo == 1) %>% dplyr::select(Species, n) %>% rename(Mycoplasma_haemo = n)->d2
  gegevens %>% group_by(Species) %>% count(Mycoplasma_coco) %>% filter(Mycoplasma_coco == 1) %>% dplyr::select(Species, n) %>% rename(Mycoplasma_coco = n)->d3
  gegevens %>% group_by(Species) %>% count(Candidatus_Neoehrlichia) %>% filter(Candidatus_Neoehrlichia == 1) %>% dplyr::select(Species, n) %>% rename(Candidatus_Neoehrlichia = n)->d4
  gegevens %>% group_by(Species) %>% count(Anaplasma) %>% filter(Anaplasma == 1) %>% dplyr::select(Species, n) %>% rename(Anaplasma = n)->d5
  gegevens %>% group_by(Species) %>% count(Francisella) %>% filter(Francisella == 1) %>% dplyr::select(Species, n) %>% rename(Francisella = n)->d6
  gegevens %>% group_by(Species) %>% count(Borrelia_1) %>% filter(Borrelia_1 == 1) %>% dplyr::select(Species, n) %>% rename(Borrelia_1 = n)->d7
  gegevens %>% group_by(Species) %>% count(Borrelia_2) %>% filter(Borrelia_2 == 1) %>% dplyr::select(Species, n) %>% rename(Borrelia_2 = n)->a7
  gegevens %>% group_by(Species) %>% count(Orientia) %>% filter(Orientia == 1) %>% dplyr::select(Species, n) %>% rename(Orientia = n)->d8
  gegevens %>% group_by(Species) %>% count(Sarcocystidae) %>% filter(Sarcocystidae == 1) %>% dplyr::select(Species, n) %>% rename(Sarcocystidae = n)->d9
  gegevens %>% group_by(Species) %>% count(Mycoplasma_penetrans) %>% filter(Mycoplasma_penetrans == 1) %>% dplyr::select(Species, n) %>% rename(Mycoplasma_penetrans = n)->d10
  gegevens %>% group_by(Species) %>% count(Mycoplasma_ravipulmonis_1) %>% filter(Mycoplasma_ravipulmonis_1 == 1) %>% dplyr::select(Species, n) %>% rename(Mycoplasma_ravipulmonis_1 = n)->d11
  
  gegevens %>% group_by(Species) %>% count(Mycoplasma_ravipulmonis_2) %>% filter(Mycoplasma_ravipulmonis_2 == 1) %>% dplyr::select(Species, n) %>% rename(Mycoplasma_ravipulmonis_2 = n)->a11
  
  gegevens %>% group_by(Species) %>% count(Mycoplasma_microti) %>% filter(Mycoplasma_microti == 1) %>% dplyr::select(Species, n) %>% rename(Mycoplasma_microti = n)->d12
  gegevens %>% group_by(Species) %>% count(Mycoplasma_insons) %>% filter(Mycoplasma_insons == 1) %>% dplyr::select(Species, n) %>% rename(Mycoplasma_insons = n)->d13
  gegevens %>% group_by(Species) %>% count(Ehrlichia) %>% filter(Ehrlichia == 1) %>% dplyr::select(Species, n) %>% rename(Ehrlichia = n)->d14
  gegevens %>% group_by(Species) %>% count(Chlamydia_1) %>% filter(Chlamydia_1 == 1) %>% dplyr::select(Species, n) %>% rename(Chlamydia_1 = n)->d15
  gegevens %>% group_by(Species) %>% count(Chlamydia_2) %>% filter(Chlamydia_2 == 1) %>% dplyr::select(Species, n) %>% rename(Chlamydia_2 = n)->a15
  gegevens %>% group_by(Species) %>% count(Rickettsia_1) %>% filter(Rickettsia_1 == 1) %>% dplyr::select(Species, n) %>% rename(Rickettsia_1 = n)->d16
  gegevens %>% group_by(Species) %>% count(Rickettsia_2) %>% filter(Rickettsia_2 == 1) %>% dplyr::select(Species, n) %>% rename(Rickettsia_2 = n)->a16
  gegevens %>% group_by(Species) %>% count(Streptobacillus_1) %>% filter(Streptobacillus_1 == 1) %>% dplyr::select(Species, n) %>% rename(Streptobacillus_1 = n)->d17
  gegevens %>% group_by(Species) %>% count(Streptobacillus_2) %>% filter(Streptobacillus_2 == 1) %>% dplyr::select(Species, n) %>% rename(Streptobacillus_2 = n)->a17
  gegevens %>% group_by(Species) %>% count(IFA_POX) %>% filter(IFA_POX == 1) %>% dplyr::select(Species, n) %>% rename(IFA_POX = n)->d18
  gegevens %>% group_by(Species) %>% count(IFA_PUUV) %>% filter(IFA_PUUV == 1) %>% dplyr::select(Species, n) %>% rename(IFA_PUUV = n)->d19
  gegevens %>% group_by(Species) %>% count(IFA_DOBRV) %>% filter(IFA_DOBRV == 1) %>% dplyr::select(Species, n) %>% rename(IFA_DOBRV = n)->d20
  gegevens %>% group_by(Species) %>% count(IFA_LCMV) %>% filter(IFA_LCMV == 1) %>% dplyr::select(Species, n) %>% rename(IFA_LCMV = n)->d21
  gegevens %>% group_by(Species) %>% count(Lepto) %>% filter(Lepto == 1) %>% dplyr::select(Species, n) %>% rename(Lepto = n)->d22
  gegevens %>% select(Species) %>% unique() %>% left_join(d1) %>% left_join(d2) %>% left_join(d3) %>% left_join(d4)%>% left_join(d5) %>% left_join(d6)%>% left_join(d7) %>% left_join(a7) %>% left_join(d8)%>% left_join(d9) %>% left_join(d10)%>% left_join(d11) %>% left_join(a11) %>%left_join(d12)%>% left_join(d13) %>% left_join(d14)%>% left_join(d15) %>% left_join(a15) %>% left_join(d16)%>% left_join(a16)%>% left_join(d17) %>% left_join(a17) %>%left_join(d18)%>% left_join(d19) %>% left_join(d20)%>% left_join(d21) %>% left_join(d22) %>% drop_na(Species) %>% replace(is.na(.), 0) -> m2
  
  
  #colnames(m2)
  colnames(m2)[2] <- "Bartonella sp."
  colnames(m2)[3] <- "M.haemo"
  colnames(m2)[4] <- "M.coco"
  colnames(m2)[5] <- "C.Neoehrlichia"
  colnames(m2)[6] <- "Anaplasma sp."
  colnames(m2)[7] <- "Francisella sp."
  colnames(m2)[8] <- "Borrelia_1 sp."
  colnames(m2)[9] <- "Borrelia_2 sp."
  colnames(m2)[10] <- "Orientia sp."
  colnames(m2)[11] <- "Sarcocystidae sp."
  colnames(m2)[12] <- "M.penetrans"
  colnames(m2)[13] <- "M.ravipulmonis_1 sp."
  colnames(m2)[14] <- "M.ravipulmonis_2 sp."
  colnames(m2)[15] <- "M.microti"
  colnames(m2)[16] <- "M.insons"
  colnames(m2)[17] <- "Ehrlichia sp."
  colnames(m2)[18] <- "Chlamydia_1 sp."
  colnames(m2)[19] <- "Chlamydia_2 sp."
  colnames(m2)[20] <- "Rickettsia_1 sp."
  colnames(m2)[21] <- "Rickettsia_2 sp."
  colnames(m2)[22] <- "Streptobacillus_1 sp."
  colnames(m2)[23] <- "Streptobacillus_2 sp."
  colnames(m2)[24] <- "Orthopox sp."
  colnames(m2)[25] <- "Orthohanta_PUUV sp."
  colnames(m2)[26] <- "Orthohanta_DOBRV sp."
  colnames(m2)[27] <- "Mammarena_LCMV sp."
  colnames(m2)[28] <- "Leptospira sp."
  
  m<-as.matrix(m2[,-1])
  rownames(m)<-m[,1]
  output <-networklevel(m)
  return(output)
}


## SIMILAR FUNCTION BUT WITH A HOST-PATHO MATRIX output per site

net.matrix <- function(site){
  
  dat3[dat3$Site_Code == site,] -> gegevens
  
  gegevens %>% group_by(Species) %>% count(Bartonella) %>% filter(Bartonella == 1) %>% dplyr::select(Species, n) %>% rename(Bartonella = n) ->d1
  gegevens %>% group_by(Species) %>% count(Mycoplasma_haemo) %>% filter(Mycoplasma_haemo == 1) %>% dplyr::select(Species, n) %>% rename(Mycoplasma_haemo = n)->d2
  gegevens %>% group_by(Species) %>% count(Mycoplasma_coco) %>% filter(Mycoplasma_coco == 1) %>% dplyr::select(Species, n) %>% rename(Mycoplasma_coco = n)->d3
  gegevens %>% group_by(Species) %>% count(Candidatus_Neoehrlichia) %>% filter(Candidatus_Neoehrlichia == 1) %>% dplyr::select(Species, n) %>% rename(Candidatus_Neoehrlichia = n)->d4
  gegevens %>% group_by(Species) %>% count(Anaplasma) %>% filter(Anaplasma == 1) %>% dplyr::select(Species, n) %>% rename(Anaplasma = n)->d5
  gegevens %>% group_by(Species) %>% count(Francisella) %>% filter(Francisella == 1) %>% dplyr::select(Species, n) %>% rename(Francisella = n)->d6
  gegevens %>% group_by(Species) %>% count(Borrelia_1) %>% filter(Borrelia_1 == 1) %>% dplyr::select(Species, n) %>% rename(Borrelia_1 = n)->d7
  gegevens %>% group_by(Species) %>% count(Borrelia_2) %>% filter(Borrelia_2 == 1) %>% dplyr::select(Species, n) %>% rename(Borrelia_2 = n)->a7
  gegevens %>% group_by(Species) %>% count(Orientia) %>% filter(Orientia == 1) %>% dplyr::select(Species, n) %>% rename(Orientia = n)->d8
  gegevens %>% group_by(Species) %>% count(Sarcocystidae) %>% filter(Sarcocystidae == 1) %>% dplyr::select(Species, n) %>% rename(Sarcocystidae = n)->d9
  gegevens %>% group_by(Species) %>% count(Mycoplasma_penetrans) %>% filter(Mycoplasma_penetrans == 1) %>% dplyr::select(Species, n) %>% rename(Mycoplasma_penetrans = n)->d10
  gegevens %>% group_by(Species) %>% count(Mycoplasma_ravipulmonis_1) %>% filter(Mycoplasma_ravipulmonis_1 == 1) %>% dplyr::select(Species, n) %>% rename(Mycoplasma_ravipulmonis_1 = n)->d11
  
  gegevens %>% group_by(Species) %>% count(Mycoplasma_ravipulmonis_2) %>% filter(Mycoplasma_ravipulmonis_2 == 1) %>% dplyr::select(Species, n) %>% rename(Mycoplasma_ravipulmonis_2 = n)->a11
  
  gegevens %>% group_by(Species) %>% count(Mycoplasma_microti) %>% filter(Mycoplasma_microti == 1) %>% dplyr::select(Species, n) %>% rename(Mycoplasma_microti = n)->d12
  gegevens %>% group_by(Species) %>% count(Mycoplasma_insons) %>% filter(Mycoplasma_insons == 1) %>% dplyr::select(Species, n) %>% rename(Mycoplasma_insons = n)->d13
  gegevens %>% group_by(Species) %>% count(Ehrlichia) %>% filter(Ehrlichia == 1) %>% dplyr::select(Species, n) %>% rename(Ehrlichia = n)->d14
  gegevens %>% group_by(Species) %>% count(Chlamydia_1) %>% filter(Chlamydia_1 == 1) %>% dplyr::select(Species, n) %>% rename(Chlamydia_1 = n)->d15
  gegevens %>% group_by(Species) %>% count(Chlamydia_2) %>% filter(Chlamydia_2 == 1) %>% dplyr::select(Species, n) %>% rename(Chlamydia_2 = n)->a15
  gegevens %>% group_by(Species) %>% count(Rickettsia_1) %>% filter(Rickettsia_1 == 1) %>% dplyr::select(Species, n) %>% rename(Rickettsia_1 = n)->d16
  gegevens %>% group_by(Species) %>% count(Rickettsia_2) %>% filter(Rickettsia_2 == 1) %>% dplyr::select(Species, n) %>% rename(Rickettsia_2 = n)->a16
  gegevens %>% group_by(Species) %>% count(Streptobacillus_1) %>% filter(Streptobacillus_1 == 1) %>% dplyr::select(Species, n) %>% rename(Streptobacillus_1 = n)->d17
  gegevens %>% group_by(Species) %>% count(Streptobacillus_2) %>% filter(Streptobacillus_2 == 1) %>% dplyr::select(Species, n) %>% rename(Streptobacillus_2 = n)->a17
  gegevens %>% group_by(Species) %>% count(IFA_POX) %>% filter(IFA_POX == 1) %>% dplyr::select(Species, n) %>% rename(IFA_POX = n)->d18
  gegevens %>% group_by(Species) %>% count(IFA_PUUV) %>% filter(IFA_PUUV == 1) %>% dplyr::select(Species, n) %>% rename(IFA_PUUV = n)->d19
  gegevens %>% group_by(Species) %>% count(IFA_DOBRV) %>% filter(IFA_DOBRV == 1) %>% dplyr::select(Species, n) %>% rename(IFA_DOBRV = n)->d20
  gegevens %>% group_by(Species) %>% count(IFA_LCMV) %>% filter(IFA_LCMV == 1) %>% dplyr::select(Species, n) %>% rename(IFA_LCMV = n)->d21
  gegevens %>% group_by(Species) %>% count(Lepto) %>% filter(Lepto == 1) %>% dplyr::select(Species, n) %>% rename(Lepto = n)->d22
  gegevens %>% select(Species) %>% unique() %>% left_join(d1) %>% left_join(d2) %>% left_join(d3) %>% left_join(d4)%>% left_join(d5) %>% left_join(d6)%>% left_join(d7) %>% left_join(a7) %>% left_join(d8)%>% left_join(d9) %>% left_join(d10)%>% left_join(d11) %>% left_join(a11) %>%left_join(d12)%>% left_join(d13) %>% left_join(d14)%>% left_join(d15) %>% left_join(a15) %>% left_join(d16)%>% left_join(a16)%>% left_join(d17) %>% left_join(a17) %>%left_join(d18)%>% left_join(d19) %>% left_join(d20)%>% left_join(d21) %>% left_join(d22) %>% drop_na(Species) %>% replace(is.na(.), 0) -> m2
  
  
  #colnames(m2)
  colnames(m2)[2] <- "Bartonella sp."
  colnames(m2)[3] <- "M.haemo"
  colnames(m2)[4] <- "M.coco"
  colnames(m2)[5] <- "C.Neoehrlichia"
  colnames(m2)[6] <- "Anaplasma sp."
  colnames(m2)[7] <- "Francisella sp."
  colnames(m2)[8] <- "Borrelia_1 sp."
  colnames(m2)[9] <- "Borrelia_2 sp."
  colnames(m2)[10] <- "Orientia sp."
  colnames(m2)[11] <- "Sarcocystidae sp."
  colnames(m2)[12] <- "M.penetrans"
  colnames(m2)[13] <- "M.ravipulmonis_1 sp."
  colnames(m2)[14] <- "M.ravipulmonis_2 sp."
  colnames(m2)[15] <- "M.microti"
  colnames(m2)[16] <- "M.insons"
  colnames(m2)[17] <- "Ehrlichia sp."
  colnames(m2)[18] <- "Chlamydia_1 sp."
  colnames(m2)[19] <- "Chlamydia_2 sp."
  colnames(m2)[20] <- "Rickettsia_1 sp."
  colnames(m2)[21] <- "Rickettsia_2 sp."
  colnames(m2)[22] <- "Streptobacillus_1 sp."
  colnames(m2)[23] <- "Streptobacillus_2 sp."
  colnames(m2)[24] <- "Orthopox sp."
  colnames(m2)[25] <- "Orthohanta_PUUV sp."
  colnames(m2)[26] <- "Orthohanta_DOBRV sp."
  colnames(m2)[27] <- "Mammarena_LCMV sp."
  colnames(m2)[28] <- "Leptospira sp."
  
  
  m<-matrix.please(m2)
  return(m)
}


### NULL model FUNCTION

nulmodel_r2d <- function(weblist, k){
  
  nulllist <- lapply(weblist, nullmodel, N=1, method="r2d")  #r2d
  observed <- meandiff(weblist)
  
  res <- 1:k
  
  for (i in 1:k){ # takes a few minutes !!
    nulllist <- sapply(weblist, nullmodel, N=1, method="r2d")
    res[i] <- meandiff(nulllist)
  }
  return(res)
}



