## README file - BioRodDis project

This repository contains the data and code used in the publication  
**<i>Drivers of host-pathogen community assemblies in European forests
and urban green spaces</i>**

More information on the **BioRodDis project** can be found
[here](https://biodiversa-bioroddis.hub.inrae.fr/bioroddis-project)

# Drivers of host-pathogen community assemblies in European forests and urban green spaces

Vincent Sluydts<sup>\*,1</sup>, Marie Bouilloud<sup>2,3</sup>,Maxime
Galan<sup>4</sup>, Hussein Alburkat<sup>\#,5</sup>, Anais
Bordes<sup>\#,4</sup>, Vincent Bourret<sup>\#,6</sup>, Valeria
Colombo<sup>\#,1</sup>, Luc DeBruyn<sup>\#,1</sup>, Lara
Dutra<sup>\#,5</sup>, Jana Eccard<sup>\#,7</sup>, Jasmin
Firozpoor<sup>\#,7</sup>, Romain Gallet<sup>\#,4</sup>, Maciej
Grzybek<sup>\#,8</sup>, Heikki Henttonen<sup>\#,9</sup>, Jens
Jacob<sup>\#,10</sup>, Andrew McManus<sup>\#,11</sup>, Tarja
Sironen<sup>\#,5</sup>, Peter Stuart<sup>\#,11</sup>, Caroline
Tatard<sup>\#,4</sup>, Benjamin Roche<sup>3</sup> Herwig
Leirs<sup>1</sup>, Nathalie Charbonnel<sup>\*,4</sup>

<sup>\#</sup> Equal contribution of coauthors (alphabetic order)

<b>Affiliations</b>

<sup>1</sup> University of Antwerp, Department of Biology, Evolutionary
Ecology Group, University of Antwerp, 2610 Wilrijk, Belgium;  
<sup>2</sup> CBGP, IRD, INRAE, CIRAD, Institut Agro, Univ Montpellier,
34000 Montpellier, France;  
<sup>3</sup> MIVEGEC, IRD, CNRS and Université de Montpellier,
Montpellier, France;  
<sup>4</sup> Univ.Helsinki, Veterinary Medicine, 00790 Helsinki,
Finland;  
<sup>5</sup> CBGP, INRAE, IRD, CIRAD, Institut Agro, Univ Montpellier,
34000 Montpellier, France;  
<sup>6</sup> CEFS, INRAE, Université de Toulouse, 31326
Castanet-Tolosan, France;  
<sup>7</sup> University of Potsdam, Institute of Biochemistry and
Biology, 14469 Potsdam, Germany;  
<sup>8</sup> Department of Tropical Parasitology, Institute of Maritime
and Tropical Medicine, Medical University of Gdansk, 81-519 Gdynia,
Poland;  
<sup>9</sup> Wildlife Ecology, Natural Research Institute, 00790
Helsinki, Finland;  
<sup>10</sup> Julius Kühn-Institute, Federal Research Institute for
Cultivated Plants, Institute for Epidemiology and Pathogen Diagnostics,
Rodent Research, Münster, Germany;  
<sup>11</sup> MunsterTechnological University, Department of Biological
and Pharmaceutical Sciences, Clash, Tralee, Ireland, V92 CX88.

<sup>\*</sup> Corresponding authors: <vincent.sluydts@uantwerpen.be>;
<Nathalie.Charbonnel@inrae.fr>

## Repository and reuse information

This open access repository accompanies the article “Drivers of
host-pathogen community assemblies in European forests and urban green
spaces” (Sluydts et al. 2024. BioRxiv). It provides data and code needed
to repeat the main analyses and reproduce the figures.  
When using the code and/or data available here, please refer to the
project license and cite the article using the following citation:

Sluydts V. et al. “Drivers of host-pathogen community assemblies in
European forests and urban green spaces”. bioRxiv 2024.

Bibtex record:

    @article{sluydts_drivers_2024,
    title = {Drivers of host-pathogen community assemblies in European forests and urban green spaces}, 
    journal = {bioRxiv},
    author = {Sluydts, Vincent and Bouilloud, Maria and Galan, Maxime and Alburkat, Hussein and Bordes Anais and Bourret, Vincent and Colombo, Valeria and DeBruyn, Luc and Dutra, Lara and Eccard, Jana and Firozpoor, Jasmin and Gallet, Romain and Grzybek, Maciej and Henttonen, Heiki and Jacob, jens and McManus, andrew and Sironen, Tarja and Stuart, Peter and Tatard, Caroline and Roche, Benjamin and Leirs, Herwig and Charbonnel, Nathalie}, 
    year = {2024}   
    }

ZENODO repository associated with the 16Sv4 rRNA gene data from the
spleen:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12518286.svg)](https://doi.org/10.5281/zenodo.12518286)

ZENODO repository associated with the 16Sv4 rRNA gene data from the
colon:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12527197.svg)](https://doi.org/10.5281/zenodo.12527197)

## Content description

- **[`Statistics_Bioroddis.Rmd`](Statistics_Bioroddis.Rmd)** : this script contains the code to reproduce the analysis regarding the community composition, the co-occurrences and the network statistics
  - [`Rest_functions.R`](Rest_functions.R): helper functions used in the Statistics_Bioroddis script.
  - [`lifehist.rds`](lifehist.rds): additional life-history data for each host extracted from Albery et al., 2022 Nature Ecology & Evolution [based on Plourde et al., 2017 PlosOne]

- **[`HMSC_setup.Rmd`](HMSC_setup.Rmd)**: this script sets up the
  required data and models for running the HMSC. It produces a object
  containing the different model settings, data and model names to be
  used in S1_run_models_A.R. The example provided is for the final model
  presented, but can be easily adapted to other models.
  - **Script Output**: an R object named `unfitted_BR4` to be used in
    `S1_run_models_A.R`.
- [`S1_run_models_A.R`](S1_run_models_A.R): this script takes the
  unfitted model input and runs the MCMC samples over a specified amount
  of samples and chains. It provides an Rdata file containing the MCMC
  outputs.
  - **Script Output**: an Rdata file saved in the directories models
    folder `models/models_BR4_100_samples_5000_chains_4.Rdata` sampling
    4 MCMC chains for 5000 samples with a thinning of 100.
    
- [`S2_compute_model_fit_A.R`](S2_compute_model_fit_A.R): this script
  takes the fitted model and calculates explanatory and predictive
  performances, as well as WAIC for a given model input.

Credit for these two automated scripts belongs to the
[team](https://www.helsinki.fi/en/researchgroups/statistical-ecology/software/hmsc)
behind the HMSC software.

- [`Fig2_Pie_Charts.R`](Fig2_Pie_Charts.R): R markdown containing all
  code needed to reproduce pie chart figure 2.  

- [`Fig3_bubble_plot.R`](Fig3_bubble_plot.R): R markdown containing all
  code needed to reproduce bubble plot figure 3.

- [`Fig4_Coinfection_plot.R`](Fig4_Coinfection_plot.R): R markdown
  containing all code needed to reproduce figure 4.

- **[`bioroddis.rds`](bioroddis.rds)**: R data object containing all
  data required to reproduce the output presented in the publication.
  This include the information on the specimens collected, the pathogens
  screened and the additional covariates extracted. Below a content
  descriptions of variables in the bioroddis data set.

**Individual host characteristics**  
- <i>“Id_code”</i> \[factor\]: Unique code per specimen collected.  
- <i>“Species”</i> \[factor\]: species name of each specimen.  
- <i>“Genus”</i> \[factor\]: Genus name of the specimen.  
- <i>“Sex”</i> \[factor\]: Sex of the specimen, M for Male and F for
Female.  
- <i>“Weight”</i> \[numeric\]: Weight of the specimen to the nearest
tenth of a gram.  
- <i>“Body_length”</i> \[numeric\]: Body length of the specimen to the
nearest tenth of a cm.  
- <i>“Tail_length”</i> \[numeric\]: Tail length of the specimen to the
nearest tenth of a cm.  
- <i>“Maturity”</i> \[factor\]: Sexual maturity of the specimen, ie:
Mature or Immature  
- <i>“Richness_gut_species”</i> \[numeric\]: Richness of the gut
microbiome from 16S analyses.  
- <i>“lrich”</i> \[numeric\]: natural logarithm transformed
Richness_gut_species.  
- <i>“Shannon_gut_species”</i> \[numeric\]: Shannon index of the gut
microbiome.  
- <i>“Simpson_gut_species”</i> \[numeric\]: Simpson index of the gut
microbiome.

**Trap session information**  
- <i>“Country”</i> \[factor\]: Country of trapped specimen.  
- <i>“Site_code”</i> \[factor\]: Code of site in country of trapped
specimen.  
- <i>“Longitude”</i> \[numeric\]: Longitude coordinates for locality.  
- <i>“Latitude”</i> \[numeric\]: Latitude coordinates for locality.  
- <i>“Season”</i> \[factor\]: Season in which the trapping was conducted
(Spring or Autumn)  
- <i>“Year”</i> \[factor\]: Year in which trapping was conducted
(2020-2021-2022). - <i>“Period”</i> \[factor\]: Combo of season and
year.  
- <i>“Spat_temp”</i> \[factor\]: Combo of site_code, season and year.

**Covariates**

- <i>“grass”</i> \[numeric\]: Percent grassland within a 1km circular
  buffer zone around the site location.  
- <i>“f.broad”</i> \[numeric\]: Percent of broadleaved forest within a
  1km circular buffer zone around the site location.  
- <i>“f.coni”</i> \[numeric\]: Percent of coniferous forest within a 1km
  circular buffer zone around the site location.
- <i>“w.perm”</i> \[numeric\]: Percent of permanent water bodies within
  a 1km circular buffer zone around the site location.
- <i>“w.temp”</i> \[numeric\]: Percent of temporary water bodies within
  a 1km circular buffer zone around the site location.
- <i>“hab2”</i> \[factor\]: Habitat either forested \[forest\] or urban
  green spaces \[urban\]
- <i>“hab3”</i> \[factor\]: Additional habitat classification with
  either forested \[forest\] or urban green spaces \[urban\], or sampled
  in a zoo \[zoo\].
- <i>“total_rain”</i> \[numeric\]: Accumulated rainfall the previous
  year since the data of sampling in a 1 km buffer around the sampling
  location.
- <i>“max_temp1”</i> \[numeric\]: Maximum temperature during the
  previous summer period
- <i>“min_temp1”</i> \[numeric\]: Minimum temperature over the previous
  winter period
- <i>“imper”</i> \[numeric\]: Percent imperviousness around a site
  locality
- <i>“pop”</i> \[numeric\]: Population density
- <i>“HSD”</i> \[numeric\]: Host Species Diversity, extracted from 10km
  buffer zone and taken as the minimum  
- <i>“HFI”</i> \[numeric\]: Human Footprint Index, extracted from within
  a 5km buffer zone

**Pathogen screening**  
Note that all screened pathogens have a numeric output with a 1 for a
positive, a 0 for negative test and an NA if the sample was not tested,
or unsuccessfully tested.  
- “Bartonella”  
- “Mycoplasma_haemo”  
- “Mycoplasma_coco”  
- “Candidatus_Neoehrlichia”  
- “Anaplasma”  
- “Francisella”  
- “Borrelia_1”  
- “Borrelia_2”  
- “Orientia”  
- “Sarcocystidae”  
- “Mycoplasma_penetrans”  
- “Mycoplasma_ravipulmonis_1”  
- “Mycoplasma_ravipulmonis_2”  
- “Mycoplasma_microti”  
- “Mycoplasma_insons”  
- “Ehrlichia”  
- “Chlamydia_1”  
- “Chlamydia_2”  
- “Rickettsia_1”  
- “Rickettsia_2”  
- “Leptospira”  
- “lepto.lipl32”  
- “Streptobacillus_1”  
- “Streptobacillus_2”  
- “IFA_POX”  
- “IFA_PUUV”  
- “IFA_DOBRV”  
- “IFA_LCMV”

<i>disease_traits3.csv</i>: CSV file containing the disease traits used
in the HMSC analysis.

### TODO:

- [ ] add rds datafile including covariate + check names
- [ ] add a cleaned script for coinfection figure 4
- [ ] check reference to final published version, title etc…
- [ ] add final model run and directories?
