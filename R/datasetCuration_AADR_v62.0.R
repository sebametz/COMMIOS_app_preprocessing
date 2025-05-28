library(tidyverse)
library(admixtools)
library(foreach)
library(doParallel)
library(leaflet)
library(shiny)
source("R/utils.R")

# Author: Metz Sebastian
# Date: 25/06/2024
# Mail: sebastian.metz[at]gmail.com

# ---- 
# Section 1: Variables definition, reference samples collection and dataset preparation
# ---- 

wdir <- file.path("/Users/lbm536/Documents/COMMIOS_app_preprocessing")
outdir <- file.path("/Users/lbm536/Documents/COMMIOS_app_preprocessing/data")


# Path to reference files: 
# Reference files are a list of individuals manually selected for the analysis. 
# These individuals are used as references and were previously published in other papers.
# The list is updated when new samples are published or new papers with better models come out.
list_aadr_files <- list(aadr_1240k = file.path(wdir, "raw-data", "AADR_V62", "v62.0_1240k_public"),
                      aadr_HO = file.path(wdir, "raw-data", "AADR_V62", "v62.0_HO_public"),
                      aadr_UK_filtered = file.path(wdir, "raw-data", "AADR_V62", "v62.0_HO_public.UKandIreland"),
                      aadr_modern_reference = file.path(wdir, "raw-data", "AADR_V62", "v62_Modern.tsv"),
                      aadr_ougroup_reference = file.path(wdir, "raw-data", "AADR_V62", "v62_Outgroup.tsv")
                      )


# Read annotation file and get aadr ancient UK genotypes ----

dataset <- read_tsv(str_c(list_aadr_files$aadr_UK_filtered, ".tsv"))

# Map view - check if localities are uniques!

dataset_v0.1 <- dataset |>
  mutate(label = sprintf("Individual ID: %s; Site: %s", `Genetic ID`, Locality))

source("R/InteractiveMap.R")
interactiveMap(df = dataset_v0.1, Lat = "Latitude", Long = "Longitude", label = "label")

rm(dataset_v0.1)


# Add the period of the individuals based on Patterson et al. 2022 ----
# NOTES:
# Neolithic (3950–2450 bc), 
# C/EBA (2450–1550 bc), 
# BA (1550–750 bc), 
# IA (750 bc to ad 43) , 
# Roman (ad 43 to ad 410), 
# Early Medieval / Vikings (ad 410 to 1066)

dataset_v1 <- dataset |>
  mutate(Period = case_match(`Date Mean in BP`,
                             c(13250:5900) ~ "Mesolithic",
                             c(5899:4400) ~ "Neolithic",
                             c(4399:3500) ~ "C/EBA",
                             c(3499:2700) ~ "BA",
                             c(2699:1907) ~ "IA",
                             c(1906:1540) ~ "Romans",
                             c(1539:885) ~ "Early Medieval/Vikings",
                             .default = "Medieval"))

# We calculate a quality of the samples based on SNPs number, the sex assignation, mtDNA coverage and Assessment.
# When sample was sequenced more than one time, 
# it is assigned as duplicated and only the first one with higher number of SNPs 
# and quality score is considered as Reference

dataset_v1 <- calculate_quality_index(dataset_v1, `SNPs hit on autosomal (1240k snpset)`) 
dataset_v2 <- dataset_v1 |>
  # If it is duplicated, use the one with higher Quality Score
  group_by(`Master ID`) |>
  mutate(`Usage Note` = if_else(Quality_score == max(Quality_score) & `SNPs hit on autosomal (1240k snpset)` == max(`SNPs hit on autosomal (1240k snpset)`), "Reference", "Duplicated")) |>
  ungroup() |>
  mutate(`Usage Note` = if_else(Quality_score < 8, str_c(`Usage Note`, "Quality Alert", sep = ";"), `Usage Note`))


# Save dataframe as an R data ----
aadr_uk_v62 <- dataset_v2
save(aadr_uk_v62, file = "data/aadr_uk_v62.rda")
write_tsv(aadr_uk_v62, "data/AADR_v62_curated.tsv")


####

# ---- 
# Section 2: Read AADR annotation and extract modern an outgroup samples
# ---- 
# - The individuals selected for outgroup and modern references are in /raw-data/Modern_and_Outgroup_individuals.tsv

dataset <- read_tsv(str_c(list_aadr_files$aadr_HO, ".anno"))
colnames(dataset) <- colnames(aadr_uk_v62)[1:ncol(dataset)]

references <- read_tsv(list_aadr_files$aadr_modern_reference)

references |>
  group_by(Note) |>
  count()

# Modern
aadr_modern_v62 <- dataset |>
  filter(`Genetic ID` %in% references$`Genetic ID`)
save(aadr_modern_v62, file = "data/aadr_modern_v62.rda")

# Outgroup
references <- read_tsv(list_aadr_files$aadr_ougroup_reference)

references |>
  group_by(Note) |>
  count()

aadr_outgroup_v62 <- dataset |>
  filter(`Genetic ID` %in% references$`Genetic ID`)

save(aadr_outgroup_v62, file = "data/aadr_outgroup_v62.rda")


######

# ---- 
# Section 3: Read AADR and extract References
# ---- 

dataset <- read_tsv(str_c(list_aadr_files$aadr_1240k, ".anno"))
colnames(dataset)
colnames(dataset) <- colnames(aadr_uk_v62[1:ncol(dataset)])
colnames(dataset)


# Countries to filter
political_entity <- c("Abkhazia", "Russia", "Albania", "Armenia", "Austria", #Adygei.HO <- Russia Assyrian <- Armenia, turkey, iran
                      "Bulgaria", "Channel Islands", "Cyprus", "Czech Republic", "Czechoslovakia", 
                      "Croatia", "Denmark", "Estonia", "Finland", "France", 
                      "Georgia", "Germany", "Gibraltar", "Greece", "Hungary",
                      "Iceland", "Iran", "Iraq", "Ireland", "Isle of Man", "Israel",
                      "Italy", "Libya", "Poland", "Spain", "Morocco", "Tunisia", 
                      "Turkey", "Yemen", "Jordan", "Lebanon", "Lithuania", "Malta",
                      "Moldova", "Netherlands", "Norway", "Portugal", "Romania",
                      "Saudi Arabia", "Serbia", "Slovakia", "Slovenia", "Sweden",
                      "Switzerland", "Syria", "Ukraine", "United Kingdom", "Belarus") 

# prepare data ----
dataset_v0 <- dataset |> 
  filter(between(`Date Mean in BP`, 750, 10600)) |>
  filter(`SNPs hit on autosomal (1240k snpset)` > 100000) |>
  filter(ASSESSMENT %in% c("PASS")) |>
  filter(`Political Entity` %in% political_entity) |>
  mutate(Period = case_match(`Date Mean in BP`,
                             c(13250:5900) ~ "Mesolithic",
                             c(5899:4400) ~ "Neolithic",
                             c(4399:3500) ~ "C/EBA",
                             c(3499:2700) ~ "BA",
                             c(2699:1907) ~ "IA",
                             c(1906:1540) ~ "Romans",
                             c(1539:885) ~ "Early Medieval/Vikings",
                             .default =  "Medieval"))

# remove duplicates based on number of SNPs
dataset_v0.1 <- dataset_v0 |>
  group_by(`Master ID`) |>
  arrange(desc(`SNPs hit on autosomal (1240k snpset)`)) |>
  distinct(`Master ID`, .keep_all = T) |>
  ungroup()


dataset_v0.2 <- calculate_quality_index(dataset_v0.1, `SNPs hit on autosomal (1240k snpset)`) 

# Ancient Ref ----

# Exclude relatives
dataset_ancient <- dataset_v0.2 |>
  filter(!`Political Entity` %in% c("United Kingdom", "Ireland")) |>
  filter(!str_detect(`Group ID`, "son\\.")) |>
  filter(!str_detect(`Group ID`, "sister\\.")) |>
  filter(!str_detect(`Group ID`, "brother\\."))|>
  filter(!str_detect(`Group ID`, "mother\\.")) |>
  filter(!str_detect(`Group ID`, "father\\.")) |>
  filter(!str_detect(`Group ID`, "daughter\\.")) |>
  filter(!str_detect(`Group ID`, "\\.rel\\."))

unique(dataset_ancient$`Group ID`)

dataset_ancient_v0.1 <- dataset_ancient |>
  mutate(Group = str_remove(`Group ID`, "\\.[ATWDSG]{2}")) 

unique(dataset_ancient_v0.1$Group)

groups <- dataset_ancient_v0.1 |>
  group_by(Group) |>
  count() |>
  filter(n > 2)

dataset_ancient_v0.2 <- dataset_ancient_v0.1 |>
  filter(Group %in% groups$Group)

# Manually curation?
write_tsv(dataset_ancient_v0.2, file.path(outdir, "AADR_v62_ancient_references.tsv"))

aadr_ancient_v62 <- dataset_ancient_v0.2

save(aadr_ancient_v62, file = "data/aadr_ancient_v62.rda")



#####
# -----
# # Section 4: Prepare dataset to extract Eigenstrat
# ----

data("aadr_uk_v62")
data("aadr_ancient_v62")
data("aadr_modern_v62")
data("aadr_outgroup_v62")

aadr_uk.sub <- select(aadr_uk_v62, `Genetic ID`, `Group ID`, `Master ID`, `Usage Note`) |>
  mutate(DataRef = "AADR_UK") 

aadr_ancient.sub <- select(aadr_ancient_v62, `Genetic ID`, `Group ID`, `Master ID`) |>
  mutate(`Usage Note` = "Use") |>
  mutate(DataRef = "AADR_Ancient") 

aadr_modern.sub <- select(aadr_modern_v62, `Genetic ID`, `Master ID`,`Group ID`) |>
  mutate(`Usage Note` = "Modern")  |>
  mutate(DataRef = "AADR_Modern")

aadr_outgroup.sub <- select(aadr_outgroup_v62, `Genetic ID`,`Master ID`, `Group ID`) |>
  mutate(`Usage Note` = "Outgroup")  |>
  mutate(DataRef = "AADR_Outgroup")

aadr_ref <- aadr_uk.sub |>
  add_case(aadr_ancient.sub) |>
  add_case(aadr_modern.sub) |>
  add_case(aadr_outgroup.sub) 

aadr_uniques <- aadr_ref


# -----
# # Section 4: Extract data from AADR HO
# ----

tempdir <- file.path(outdir, "temp_v62")
system2("mkdir", args = sprintf("-p %s", tempdir))
# 
write_tsv(select(aadr_uniques, 1:3), file.path(tempdir, "aadr_uniques.tsv"))
# 
aadr_ho <- update_reference(input = file.path(tempdir, "aadr_uniques.tsv"), output = file.path(tempdir, "aadr_HO_v62_sub"),
                            aadr_path = "/Users/lbm536/Documents/COMMIOS_app_preprocessing/raw-data/AADR_V62/v62.0_HO_public", by_id = T)

# Check if all the individual are in the aadr_ho
aadr_ho_ind <- read_table(str_c(aadr_ho, ".ind"), col_names = c("Genetic ID", "sex", "groupid"))

length(unique(aadr_ho_ind$`Genetic ID`))
length(unique(aadr_ref$`Genetic ID`))

# -----
# # Section 5: PCA
# ----
# NOTES:
# PCA performed with all the modern individuals and those that are annotated as References
# The rest of the samples are projected separately

# PCA all References ----

## Extract samples to use [Reference (AADR_UK), Use (AADR_Ancinet), all (AADR_Modern)]

out <- file.path(outdir, "PCA")
system2("mkdir", args = sprintf("-p %s", out))

# aadr_uniques <- read_tsv(file.path(tempdir, "aadr_uniques.tsv"))

aadr_pca <- aadr_ref |>
  filter(str_detect(`Usage Note`, "Modern|Use|Reference")) |>
  mutate(Group = if_else(`Usage Note` == "Modern", `Group ID`, `Genetic ID`))

path_aadr_pca <- file.path(out, "aadr_PCA_References_v62.tsv")

write_tsv(aadr_pca, path_aadr_pca)

geno_pca <- update_reference(input = path_aadr_pca, output = file.path(out, "aadr_pca_ref_v62"), 
                             aadr_path = aadr_ho, by_id = T)

aadr_pca_ind <- edit_indfile(input = str_c(geno_pca, ".ind"), table = aadr_pca, output = str_c(geno_pca, "_modv62.ind"),
                             by_col = "Genetic ID", remp_col = "Group")


## Prepare for PCA
poplist <- aadr_pca |>
  filter(`Usage Note` == "Modern") |>
  group_by(Group) |>
  distinct(Group) 

write_tsv(poplist, file = "data/PCA/modern.poplist", col_names = F)

run_pca(genotype = geno_pca, poplist = file.path(out, "modern.poplist"), ind = aadr_pca_ind, outdir = out)

# Add PCA to Ref
pca.eval <- read_table("data/PCA/aadr_pca_ref_v62.evec", skip = 1, col_names = c("Genetic ID", str_c("PCA", 1:10), "PCA Group ID"))

aadr_ref.n <- left_join(aadr_ref, pca.eval, by = c("Genetic ID" = "Genetic ID"))

pca_v62 <- aadr_ref.n
save(pca_v62, file = "data/pca_v62.rda")

## - Chekc plot
aadr_ref.n |>
  mutate(Group = if_else(str_detect(`Group ID`, "England|Scotland|Ireland"), "B", "A")) |>
  mutate(Group = if_else(str_detect(`Group ID`, "Spain"), "C", Group)) |>
  arrange(Group) |>
  ggplot() +
  geom_point(aes(PCA1, PCA2, fill = Group), shape = 21, size = 2) +
  scale_fill_manual(values = c("blue", "gray30", "red"))

# -----
# # Section 6: qpAdmixture with ancient individuals
# ----
# NOTES:
# qpAdmixture performed with outgroups from Patterson et al., 2022
# The admixture is performed to all AADR_UK individuals and AADR_Ancient individuals

# Prepare genotypes 1240K and outgroups ----
out <- file.path(outdir, "qpAdmixture")
system2("mkdir", args = sprintf("-p %s", out))

# Reference file of selected individuals as outgroup
outgroup_references <- read_tsv(list_aadr_files$aadr_ougroup_reference)

# Prepare aadr_ref.n UK for admixtoos
uk_references_subset <- aadr_uk_v62 |>
  select(`Genetic ID`, `Master ID`, `Group ID`,`Political Entity`, Period, `Usage Note`) |>
  rename(Note = `Usage Note`, Population = `Political Entity`) |>
  mutate(Group = `Genetic ID`)

aadr_admixture <- outgroup_references |>
  add_case(uk_references_subset)

path_aadr_qpadm <- file.path(out, "aadr_qpAdm_v62.tsv")
write_tsv(select(aadr_admixture, `Genetic ID`, `Master ID`, `Group ID`), path_aadr_qpadm)

geno_qpadm <- update_reference(input = path_aadr_qpadm, output = file.path(out, "aadr_apadm"), 
                               aadr_path = list_aadr_files$aadr_1240k, by_id = T)


# Modify Group ID of outgroup to have the same groups than in Patterson et al., 2022 for the admixture analysis ----
write_tsv(aadr_admixture, path_aadr_qpadm)
aadr_qpadm_ind <- edit_indfile(input = str_c(geno_qpadm, ".ind"), table = aadr_admixture, output = str_c(geno_qpadm, ".ind"),
                               by_col = "Genetic ID", remp_col = "Group")

# qpAdmixture ancestry composition Neolithics ----
left <- c("WHG", "EEF")
right <- c("OldAfrica", "Afanasievo", "Turkey_N", "WHGB")
targets <- filter(aadr_admixture, Group == `Genetic ID` & Period == "Neolithic")$`Genetic ID`

# Activate paralell
# Parallel configuration
parallel::detectCores()
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

qpadm_ancestry_neolithic_results <- foreach(
  i = 1:length(targets),
  .combine = 'rbind'
) %dopar% {
  library(admixtools)
  library(tidyverse)
  run_qpadm_ancestry(geno_qpadm, left, right, targets[i])
}
# Deactivate parallel
parallel::stopCluster(cl = my.cluster)
write_rds(qpadm_ancestry_neolithic_results, "data/apadam_v62_neolithics.rds")
# qpAdmixture ancestry composition ----
left <- c("OldSteppe","WHG", "EEF")
right <- c("OldAfrica", "Afanasievo", "Turkey_N", "WHGB")
targets <- filter(aadr_admixture, Group == `Genetic ID` & Period != "Neolithic")$`Genetic ID`

# Activate paralell
# Parallel configuration
parallel::detectCores()
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

qpadm_ancestry_results <- foreach(
  i = 1:length(targets),
  .combine = 'rbind'
) %dopar% {
  library(admixtools)
  library(tidyverse)
  run_qpadm_ancestry(geno_qpadm, left, right, targets[i])
}
# Deactivate parallel
parallel::stopCluster(cl = my.cluster)

write_rds(qpadm_ancestry_results, "data/apadam_v62_others.rds")

# Prepare for merge
qpadm <- qpadm_ancestry_results

qpadm.f <- NULL
for (i in 1:nrow(qpadm)) {
  aux <-  qpadm[i,1][[1]]
  aux.f <- tibble(`Genetic ID` = unique(aux$target), `qpAdm Pvalue` = qpadm[i,2][[1]]$p[1],
                  Steppe = aux$weight[1], `Steppe SE` = aux$se[1],
                  WHG = aux$weight[2], `WHG SE` = aux$se[2],
                  EEF = aux$weight[3], `EEF SE` = aux$se[3])
  qpadm.f <- if(length(qpadm.f)) add_case(qpadm.f, aux.f) else aux.f
}

qpadm <- qpadm_ancestry_neolithic_results
for (i in 1:nrow(qpadm)) {
  aux <-  qpadm[i,1][[1]]
  aux.f <- tibble(`Genetic ID` = unique(aux$target), `qpAdm Pvalue` = qpadm[i,2][[1]]$p[1],
                  Steppe = 0, `Steppe SE` = 0,
                  WHG = aux$weight[1], `WHG SE` = aux$se[1],
                  EEF = aux$weight[2], `EEF SE` = aux$se[2])
  qpadm.f <- if(length(qpadm.f)) add_case(qpadm.f, aux.f) else aux.f
}


qpadm_v62 <- qpadm.f
save(qpadm_v62, file = "data/qpadm_v62.rda")

# -----
# # Section 7: Final table aadr + pca + qpAdm
# ----
# Merge all the data into a database to share


data("aadr_uk_v62")
data("aadr_ancient_v62")
data("aadr_modern_v62")
data("aadr_outgroup_v62")
data("pca_v62")
data("qpadm_v62")

# Add group to aadr_Uk
aadr_uk_v62 <- aadr_uk_v62 |> 
  mutate(Group = str_remove_all(`Group ID`, "\\.[ATWDSG]{2}")) |>
  mutate(`Published Year` = as.character(`Published Year`))

unique(aadr_uk_v62$Group)

# Add Usage Note to ancient
map_chr(aadr_uk_v62, typeof) == map_chr(aadr_ancient_v62, typeof) 

aadr_ancient_v62 <- aadr_ancient_v62 |>
  mutate(`Usage Note` = "Ancient") |>
  mutate(Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude)) |>
  mutate(`SNPs hit on autosomal (1240k snpset)` = as.numeric(`SNPs hit on autosomal (1240k snpset)`)) |>
  mutate(`mtDNA coverage` = as.numeric(`mtDNA coverage`))

colnames(aadr_ancient_v62) %in% colnames(aadr_uk_v62)  
  
# Add Period, Group, Quality_Score and Usage Note to modern and outgroup
reference_modern <- read_tsv(list_aadr_files$aadr_modern_reference)

map_chr(aadr_uk_v62, typeof) == map_chr(aadr_modern_v62, typeof) 

aadr_modern_v62 <- aadr_modern_v62 |>
  mutate(Period = "Present", Quality_score = as.numeric(NA), `Usage Note` = "Modern") |>
  mutate(Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude)) |>
  mutate(`SNPs hit on autosomal (1240k snpset)` = as.numeric(`SNPs hit on autosomal (1240k snpset)`)) |>
  mutate(`mtDNA coverage` = as.numeric(`mtDNA coverage`)) |>
  left_join(select(reference_modern, `Genetic ID`, Group))

reference_outgroup <- read_tsv(list_aadr_files$aadr_ougroup_reference)

map_chr(aadr_uk_v62, typeof) == map_chr(aadr_outgroup_v62, typeof) 

aadr_outgroup_v62 <- aadr_outgroup_v62 |>
  mutate(Period = "-", Quality_score = as.numeric(NA), `Usage Note` = "Outgroup") |>
  mutate(Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude)) |>
  mutate(`SNPs hit on autosomal (1240k snpset)` = as.numeric(`SNPs hit on autosomal (1240k snpset)`)) |>
  mutate(`mtDNA coverage` = as.numeric(`mtDNA coverage`)) |>
  left_join(select(reference_modern, `Genetic ID`, Group))

# 
aadr_database_v62 <- aadr_uk_v62 |>
  add_case(aadr_ancient_v62) |>
  add_case(aadr_modern_v62) |>
  add_case(aadr_outgroup_v62) |>
  left_join(select(pca_v62, !c(`Group ID`, `Master ID`, `Usage Note`, DataRef, `PCA Group ID`))) |>
  left_join(qpadm_v62)

save(aadr_database_v62, file = "data/aadr_database_v62.rda")
write_tsv(aadr_database_v62, "data/aadr_database_v62.tsv")


# Prepare for app
data("aadr_database_v62")

aadr_database_v62.f <- aadr_database_v62 |>
  mutate(Country = `Political Entity`) |>
  mutate(Country = if_else(`Political Entity` == "United Kingdom", str_extract(Locality, "Scotland|England|Wales"), Country)) |>
  mutate(Country = if_else(str_detect(Locality,"Gloucestershire|Oxfordshire|Lincolnshire"), "England", Country)) |>
  mutate(Country = if_else(str_detect(Locality, "Orkney"), "Scotland", Country)) |>
  mutate(Country = if_else(str_detect(Locality, "Isle of Man"), "Isle of Man", Country)) |>
  #DataRef = "AADR_Modern"     "AADR_References" "AADR_UK"
  mutate(DataRef = if_else(str_detect(`Usage Note`, "Reference|Duplicated"), "AADR_UK", NA)) |>
  mutate(DataRef = if_else(str_detect(`Usage Note`, "Ancient"), "AADR_References", DataRef)) |>
  mutate(DataRef = if_else(str_detect(`Usage Note`, "Modern"), "AADR_Modern", DataRef)) |>
  mutate(DataRef = if_else(str_detect(`Usage Note`, "Outgroup"), "AADR_Outgroup", DataRef)) |>
  #Reference = Country;Period;Country_Period
  mutate(Period = if_else(Period == "Present", "Modern", Period)) |>
  mutate(Reference = str_replace(str_c(str_c(Country, Period, sep = ";"), str_c(Country, Period, sep = "_"), sep = ";"), " ", "_")) |>
  #other corrections
  mutate(Locality = if_else(`Genetic ID` == "I16412_d.AG", "Coneypark Cairn [Cist 1] (Scotland, Stirling)", Locality)) |>
  mutate(Locality = if_else(str_detect(Locality, "\\(JLU15\\)"), str_replace(Locality, "\\(JLU15\\)", "[JLU15]"), Locality)) |>
  mutate(Locality = if_else(str_detect(Locality, "\\(field 9\\)"), str_replace(Locality, "\\(field 9\\)", "[field 9]"), Locality)) |>
  mutate(Locality = if_else(str_detect(Locality, "\\(field 10\\)"), str_replace(Locality, "\\(field 10\\)", "[field 10]"), Locality)) |>
  mutate(Locality = if_else(str_detect(Locality, "\\(field 16\\)"), str_replace(Locality, "\\(field 16\\)", "[field 16]"), Locality)) |>
  mutate(Locality = if_else(str_detect(Locality, "\\(field 13\\)"), str_replace(Locality, "\\(field 13\\)", "[field 13]"), Locality)) |>
  mutate(Locality = if_else(str_detect(Locality, "\\(2005\\)"), str_replace(Locality, "\\(2005\\)", "[2005]"), Locality)) |>
  mutate(Locality = if_else(str_detect(Locality, "\\(Baugebiet-110\\)"), str_replace(Locality, "\\(Baugebiet-110\\)", "[Baugebiet-110]"), Locality)) |>
  mutate(Locality = if_else(Locality == "Farman (Île-de-France, Paris)", "Farman (Paris)", Locality)) |>
  mutate(Locality = if_else(Locality == "Ovilava (Upper Austria, Wels)", "Ovilava (Austria, Wels)", Locality))|>
  mutate(Locality = if_else(Locality == "Ebla / Tell Mardikh (northwestern Syria, Idlib Governorate)", "Ebla / Tell Mardikh (Syria, Idlib Governorate)", Locality))
# 

# For duplicated samples copy values from Reference
duplicated <- unique(filter(aadr_database_v62.f, str_detect(`Usage Note`, "Duplicated"))$`Master ID`)
for (i in duplicated) {
  ref <- filter(aadr_database_v62.f, `Master ID` == i, str_detect(`Usage Note`, "Reference") )
  if(nrow(ref)>0) {
    ref <- select(ref, 43:66)
    aadr_database_v62.f[grep(T, aadr_database_v62.f$`Master ID` == i), 43:66] <- ref
  } else {
    aadr_database_v62.f <- aadr_database_v62.f |>
      mutate(`Usage Note` = if_else(`Master ID` == i, str_c(`Usage Note`, "Ignore", sep = ";"), `Usage Note`))
  }
}

# Family
aadr_database_v62.ff <- aadr_database_v62.f |>
  mutate(Family = if_else(str_detect(`Group ID`, "(son|sister|brother|mother|father|daughter|1d\\.rel|possible\\.1d).*"), 
                             str_extract(`Group ID`, "(son|sister|brother|mother|father|daughter|1d\\.rel|possible\\.1d).*"),
                             NA)) |>
  mutate(Group = if_else(!is.na(Family), 
                         str_remove(Group, "_(son|sister|brother|mother|father|daughter|1d\\.rel|possible\\.1d).*"),
                         Group))

# Add locality with the format needed for app
aadr_database_v62.fff <- aadr_database_v62.ff |>
  mutate(aux = str_c(str_extract(Locality, "\\((.*)\\)", group = 1), str_extract(Locality, "(.*)\\(", group = 1), sep = ", ")) |>
  mutate(aux = if_else(is.na(aux), Locality, aux)) |>
  mutate(aux = if_else(str_detect(aux, Country), aux, str_c(Country, aux, sep = ", "))) |>
  rename(`Unique Locality Name` = aux)


unique(aadr_database_v62.fff$`Unique Locality Name`)
aadr_ref_v62 <- aadr_database_v62.fff

save(aadr_ref_v62, file = "data/aadr_ref_v62.rda")
save(aadr_ref_v62, file = "../COMMIOS_app/data/aadr_ref_v62.rda")
write_tsv(aadr_ref_v62, "data/aadr_ref_v62.tsv")


# Extract final genotypes
data("aadr_ref_v62")



tempdir <- file.path(outdir, "genotypes_v62")
system2("mkdir", args = sprintf("-p %s", tempdir))
# 
write_tsv(select(aadr_ref_v62, 1:3), file.path(tempdir, "aadr_ref_v62.tsv"))
# 
aadr_ho <- update_reference(input = file.path(tempdir, "aadr_ref_v62.tsv"), output = file.path(tempdir, "aadr_HO_v62"),
                            aadr_path = "/Users/lbm536/Documents/COMMIOS_app_preprocessing/raw-data/AADR_V62/v62.0_HO_public", by_id = T)

aadr_1240 <- update_reference(input = file.path(tempdir, "aadr_ref_v62.tsv"), output = file.path(tempdir, "aadr_1240k_v62"),
                            aadr_path = "/Users/lbm536/Documents/COMMIOS_app_preprocessing/raw-data/AADR_V62/v62.0_1240k_public", by_id = T)

# anno file
aadr_anno <- read_tsv(str_c(list_aadr_files$aadr_HO, ".anno"))
aadr_ref_v62.anno <- aadr_anno |>
  filter(`Genetic ID (suffixes: ".DG" is a high coverage shotgun genome with diploid genotype calls, ".AG" is shotgun data with each position in the genome represented by a randomly chosen sequence, ".HO" is Affymetrix Human Origins genotype data)` %in% aadr_ref_v62$`Genetic ID`)
write_tsv(aadr_ref_v62.anno, file.path(tempdir, "aadr_HO_v62.anno"))



##### Figures for paper ----

# Plot COMMIOS Paper
commios <- read_tsv("data/commios_v63.tsv")

aadr_commios <- aadr_uk_v62 |>
  filter(`Genetic ID` %in% commios$`Genetic ID`)

length(unique(aadr_uk_v62$`Master ID`))
length(unique(aadr_commios$`Master ID`))


# map / sample distribution:
# 1: The samples from commios in blue and yellow the ones from aadr

# Libraries
library(ggplot2)
library(dplyr)

# Get the world polygon and extract UK
library(giscoR)
UK <- gisco_get_countries(country = c("UK", "Ireland"), resolution = 1)


# Left chart
p <- ggplot() +
  geom_sf(data = UK, fill = "lightblue", color = "darkblue", alpha = 0.7) +
  geom_point(data = filter(aadr_uk_v62, !`Master ID` %in% commios$`Master ID`), aes(x = Longitude, y = Latitude, fill = "AADR individual"), shape = 21, size = 2) +
  geom_point(data = filter(aadr_uk_v62, `Master ID` %in% commios$`Master ID`), aes(x = Longitude, y = Latitude, fill = "COMMIOS individual"), shape = 21, size = 2) +
  scale_fill_manual(
    name = "Individuals",  # Legend title
    values = c("AADR individual" = "#FFD700",
               "COMMIOS individual" = "#FF4500")) +
  theme_bw() +
  theme(
    legend.position = "bottom",  # Moves legend below the plot
    legend.title = element_text(size = 12),  # Title font size
    legend.text = element_text(size = 10)    # Labels font size
  ) +
  ylim(50, 61)

ggsave("Figure1A_map.pdf", p,  units = "cm", device = "pdf", dpi = 600)


aux <- aadr_commios |>
  select(`Genetic ID`, `Master ID`, Period, `Date Mean in BP`) |>
  mutate(Source = "COMMIOS") |>
  add_case(select(filter(aadr_uk_v62, !`Master ID` %in% aadr_commios$`Master ID`), `Genetic ID`, `Master ID`, Period, `Date Mean in BP`) |>
             mutate(Source = "AADR"))


library(ggplot2)
library(scales)


# mutate(Period = case_match(`Date Mean in BP`,
#                            c(13250:5900) ~ "Mesolithic",
#                            c(5899:4400) ~ "Neolithic",
#                            c(4399:3500) ~ "C/EBA",
#                            c(3499:2700) ~ "BA",
#                            c(2699:1907) ~ "IA",
#                            c(1906:1540) ~ "Romans",
#                            c(1539:1885) ~ "Early Medieval/Vikings",
#                            .default =  "Medieval"))
# Define periods for annotation
periods <- data.frame(
  start = c(-5899, -4399, -3499, -2699, -1906, -1539),  # Start year BC
  end = c(-4400, -3500, -2700, -1907, -1540, -885),      # End year (BC or AD)
  label = c("Neolithic (3950–2450 BC)", 
            "C/EBA (2450–1550 BC)", 
            "BA (1550–750 BC)", 
            "IA (750 BC to AD 43)",
            "Roman (AD 43 to AD 410)",
            "Early Medieval/Vikings (AD 410 to 885)")
)

# Sample color mapping (same as map)
sample_colors <- c(
  "AADR" = "#FFD700",  # Yellow
  "COMMIOS" = "#FF4500"       # Red-Orange
)

# Create the plot
p2 <- ggplot() +
  # Jittered points
  geom_jitter(data = aux, 
              aes(x = -`Date Mean in BP`, y = "Britain", colour = Source),  # Flip BP to cal. BC
              width = 50, height = 0.1, size = 2, shape = 20) +
  # Add period annotations
  geom_rect(data = periods, 
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = label), 
            alpha = 0.2, inherit.aes = FALSE, colour = "gray80") +
  # geom_text(data = periods, 
  #           aes(x = (start + end) / 2, y = 2, label = label), 
  #           inherit.aes = FALSE, size = 3, angle = 45, hjust = 1, vjust = 1) +
  
  # Adjust colors
  scale_colour_manual(name = "Source", values = sample_colors) +
  
  # Customize x-axis for years in BC
  scale_x_continuous(name = "Time (years cal. BC)", 
                     labels = scales::comma_format(),  # Add commas for large numbers
                     breaks = seq(-6000, 100, by = 500)) +  # Adjust ticks
  
  # Custom styling
  theme_minimal() +
  labs(y = "", colour = "Source") +
  xlim(-6000, 0) +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80")
  )

ggsave("Figure1C_transect.pdf", p2,  units = "cm", device = "pdf", dpi = 600)

