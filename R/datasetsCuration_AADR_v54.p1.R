library(tidyverse)
library(admixtools)
library(foreach)
library(doParallel)
library(leaflet)
library(shiny)
# source("~/Documents/Reports/Reports_Gen/EIG_functions.R")
source("R/utils.R")

# Author: Metz Sebastian
# Date: 25/06/2024
# Mail: sebastian.metz[at]gmail.com



# ---- 
# Section 1: Read AADR annotation and data curation for UK and Ireland samples
# ---- 

# - Read AADR from the webpage - OK
# - Download files into a version control folder (AADR_Version) - OK
# - Extract all the UK samples - OK
# - Curate dataset:
#   - Sites are in the UK OK
#   - Sites <- Check each site have only one lat. long. OK
#   - Manually review the sites and update sites in the annotation table OK
#   - Annotate the Periods based on Patterson et al. 2022 OK
#   - Check the quality of the samples:
#       - N of SNPs
#       - Potential contamination
#       - Is duplicated (select the one with high number of SNPs)
#       - Family annotation (If have relatives get the one with high number of SNPs for that Family)

wdir <- file.path("/Users/lbm536/Documents/COMMIOS_app_preprocessing/raw-data","AADR_v54.1.p1")
outdir <- file.path("/Users/lbm536/Documents/COMMIOS_app_preprocessing/data")

# Download files from web page and untar ----

list_files <- list(path_1240k = "https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V54/V54.1.p1/SHARE/public.dir/v54.1.p1_1240K_public.tar",
                   path_HO = "https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V54/V54.1.p1/SHARE/public.dir/v54.1.p1_HO_public.tar")
  
# download_aadr_files(save_path = wdir, download_all = list_files)

list_aadr_file <- list(`1240k` = file.path(wdir, str_replace(basename(list_files$path_1240k), ".tar", "")),
                       HO = file.path(wdir, str_replace(basename(list_files$path_HO), ".tar", "")))

# untar(tarfile = str_c(list_aadr_file$`1240k`, ".tar"), exdir = dirname(list_aadr_file$`1240k`), verbose = TRUE)
# untar(tarfile = str_c(list_aadr_file$HO, ".tar"), exdir =  dirname(list_aadr_file$HO), verbose = TRUE)

#
# Read anno file and get aadr ancient UK genotypes ----

dataset <- read_tsv(str_c(list_aadr_file$`1240k`, ".anno"))
dataset_v0 <- dataset |>
  mutate(`Full Date One of two formats. (Format 1) 95.4% CI calibrated radiocarbon age (Conventional Radiocarbon Age BP, Lab number) e.g. 2624-2350 calBCE (3990±40 BP, Ua-35016). (Format 2) Archaeological context range, e.g. 2500-1700 BCE` = if_else(`Genetic ID` == "I21275", "384-202 calBCE (2229±24 BP, SUERC-104565)", `Full Date One of two formats. (Format 1) 95.4% CI calibrated radiocarbon age (Conventional Radiocarbon Age BP, Lab number) e.g. 2624-2350 calBCE (3990±40 BP, Ua-35016). (Format 2) Archaeological context range, e.g. 2500-1700 BCE`)) |>
  filter(`Political Entity` %in% c("United Kingdom", "Isle of Man", "Ireland") & `Method for Determining Date; unless otherwise specified, calibrations use 95.4% intervals from OxCal v4.4.2 Bronk Ramsey (2009); r5; Atmospheric data from Reimer et al (2020)` != "Modern")

rm(dataset)
#
# Use Lat. Long. to map and filter ----
unique(filter(dataset_v0, Lat. == ".." | Long. == "..")$Locality)

### Correction:
  # "West Sussex, Chichester, Apple Down"
  # This location is in Lat. = 50.916667	and Long. = -0.5 (Google Maps)

dataset_v0.1 <- dataset_v0 |>
  mutate(Latitude = if_else(Lat. == "..", "50.916667", Lat.), Longitude = if_else(Long. == "..", "-0.5", Long.)) |>
  mutate(Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude)) |>
  select(!c(Lat., Long.)) |>
  mutate(label = sprintf("Individual ID: %s; Site: %s", `Genetic ID`, Locality))

source("R/InteractiveMap.R")
interactiveMap(df = dataset_v0.1, Lat = "Latitude", Long = "Longitude", label = "label")

rm(dataset_v0)
#

# Site ID generation ----
### NOTES:
# Site names have different Lat. Long. => Check name to see if they are caves or somthing like that
# Assign a new id (Loc1:n) based on Lat and Long
# Site will be represented by ID / Name / Location / Country


sites_collection <- dataset_v0.1 |>
  select(Locality, Latitude, Longitude) |>
  mutate(LocalityF = if_else(Locality == "Anglesey, Wales, Glyn Llanbedrgoch", "Wales, Anglesey, Glyn Llanbedrgoch", Locality)) |>
  mutate(LocalityF = if_else(Locality == "Ballynahatty,  County Down", "Ireland, Ballynahatty, County Down", LocalityF)) |>
  mutate(LocalityF = if_else(Locality == "Carrowmore", "Ireland, Carrowmore", LocalityF)) |>
  mutate(LocalityF = if_else(Locality == "Dublin, Finglas", "Ireland, Dublin, Finglas", LocalityF)) |>
  mutate(LocalityF = if_else(Locality == "Dublin, Islandbridge", "Ireland, Dublin, Islandbridge", LocalityF)) |>
  mutate(LocalityF = if_else(Locality == "Dublin, Ship Street Great", "Ireland, Dublin, Ship Street Great", LocalityF)) |>
  mutate(LocalityF = if_else(Locality == "Eyrephort", "Ireland, Eyrephort", LocalityF)) |>
  mutate(LocalityF = if_else(Locality == "IsleOfMan, Balladoole", "Isle of Man, Balladoole", LocalityF)) |>
  mutate(LocalityF = if_else(Locality == "Newgrange", "Ireland, Newgrange", LocalityF)) |>
  mutate(LocalityF = if_else(Locality == "Orkney, Westray, Knowe of Skea", "Scotland, Orkney, Westray, Knowe of Skea", LocalityF)) |>
  mutate(LocalityF = if_else(Locality == "Orkney, Westray, Links of Noltland", "Scotland, Orkney, Westray, Links of Noltland", LocalityF)) |>
  mutate(LocalityF = if_else(Locality == "Primrose Grange", "Ireland, Primrose Grange", LocalityF)) |>
  mutate(LocalityF = if_else(Locality == "Rathlin Island,  County Antrim", "Ireland, Rathlin Island, County Antrim", LocalityF)) |>
  mutate(LocalityF = if_else(Locality == "Roscommon, Kilteashe", "Ireland, Roscommon, Kilteashe", LocalityF)) |>
  mutate(LocalityF = if_else(Locality == "West Sussex, Chichester, Apple Down", "England, West Sussex, Chichester, Apple Down", LocalityF)) |>
  mutate(LocalityF = if_else(Locality == "England, Co. Durham, Hartlepool, Catcote", "England, County Durham, Hartlepool, Catcote", LocalityF))

sites_collectionF <- sites_collection |>
  mutate(LocalityF = str_replace_all(LocalityF, "\\s{2,}", "\\s"))|>
  mutate(LocalityF = str_replace_all(LocalityF, "\\'", ""))|>
  arrange(LocalityF) |>
  mutate(tmp = str_c(Latitude, Longitude)) |>
  group_by(tmp) |>
  mutate(LocID = Vectorize(digest::sha1)(tmp)) |>
  ungroup() |>
  #Sites with multiple given names
  mutate(CurrentName = if_else(LocID == "4cd750251a1206c349ccbd5f5fe8ec71fdddbcb4", "England, East Riding of Yorkshire, Melton Quarry", LocalityF)) |>
  mutate(CurrentName = if_else(LocID == "0dd1635a4dfb6dcac2ebe940431a4ee5ec334cfb", "England, North Yorkshire, Scorton Quarry", CurrentName)) |>
  mutate(CurrentName = if_else(LocID == "2f1f1427b8cedc8d5b5c5332c2fe7f7df225f71f", "England, East Riding of Yorkshire, East Coast Pipeline (field 16)", CurrentName)) |>
  mutate(CurrentName = if_else(LocID == "48414d8882c330707067cc6ef77c40d27448a788", "England, Dorset, Worth Matravers", CurrentName)) |>
  mutate(CurrentName = if_else(LocID == "65a1436de463afeb11877d99f4093cb1022157da", "Scotland, Argyll and Bute, Oban, Raschoille Cave", CurrentName)) |>
  mutate(CurrentName = if_else(LocID == "98f53cf5b6254e1160a494dc1ec49692ebb13759", "England, North Yorkshire, Stockton-on-Tees, Norton, Norton Bishops Mill", CurrentName)) |>
  mutate(CurrentName = if_else(LocID == "a6b4d5894a7ea29679ecaf5da20f5b27800d7bab", "England, North Yorkshire, Stockton-on-Tees, Norton, Norton East Mill", CurrentName)) |>
  mutate(CurrentName = if_else(LocID == "c48358f7f0c993e42a24a33da6354e37bb9a668b", "England, Somerset, South Cadbury, Cadbury Castle", CurrentName)) |>
  mutate(CurrentName = if_else(LocID == "c7182524eba245f4198940a896ef9a7f7282c602", "Wales, West Glamorgan, Gower Peninsula, Port Eynon, Culver Hole Cave", CurrentName)) |>
  mutate(CurrentName = if_else(LocID == "e91c7b4306175a61456fbecd8e0266c272dd5fdd", "England, Blaydon, Tyne and Wear, Bewes Hill", CurrentName)) |>
  group_by(LocID, Locality) |>
  reframe(CorrectLocalityName = LocalityF, UniqueLocalityName = CurrentName, Latitude = Latitude, Longitude = Longitude, nSamples = n()) |>
  ungroup()

write_tsv(sites_collection, file.path(outdir, "site_collections.tsv"))


add_site <- function(GenID = str(), Loc = str(), Lat = num(), Long = num()) {
  filter <- filter(sites_collectionF, Locality == Loc, Latitude == Lat, Longitude == Long)[1,]
  filter$`Genetic ID` <- GenID
  filter
}

df_v0.2 <- dataset_v0.1 
sites_to_add <- NULL
for (i in 1:nrow(df_v0.2)) {
  aux <- add_site(df_v0.2$`Genetic ID`[i], df_v0.2$Locality[i], df_v0.2$Latitude[i], df_v0.2$Longitude[i])
  sites_to_add <- if(length(sites_to_add)) add_case(sites_to_add, aux) else aux
}

df_v0.3 <- df_v0.2 |>
  left_join(select(sites_to_add, !c("Locality", "Latitude", "Longitude")), by = c("Genetic ID"))  #



# Samples periods ----
  # NOTES:
  # Neolithic (3950–2450 bc), 
  # C/EBA (2450–1550 bc), 
  # BA (1550–750 bc), 
  # IA (750 bc to ad 43) , 
  # Roman (ad 43 to ad 410), 
  # Early Medieval / Vikings (ad 410 to 1066)

df_v0.4 <- df_v0.3 |>
  mutate(DateMeanInBP = as.numeric(`Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]`)) 

df_v0.4$Period <- case_match(df_v0.4$DateMeanInBP,
                             c(13250:5899) ~ "Mesolithic",
                             c(5900:4399) ~ "Neolithic",
                             c(4400:3499) ~ "C/EBA",
                             c(3500:2699) ~ "BA",
                             c(2700:1906) ~ "IA",
                             c(1907:1539) ~ "Romans",
                             c(1540:884) ~ "Early Medieval/Vikings")

# Samples quality control ----
  # Note:
  # Test 1: ASSESMENT = PASS | PASS_PUBLICATION = 1 / 0
  # Test 2: Group ID match Ignore | outlier = If Ignore => 0, If Outlier = 0.5
  # Test 3: nSNPs = Rank from Min to Max | Make it log? = Rank
  # Test 4: Sex was determined | Sex Ratio high confidence = If M | F => 1 else 0.5
  # Test 4: Has Y-Chromosome and mtDNA chromosome  = If M and Y-Chromosome & mtDNA => 1 If F & mtDNA => 1 else 0.5
  # Test 5: Has family

df_v0.5 <- df_v0.4 |>
  mutate(Quality = if_else(`SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)` >= 100000, "Ref;SNPs 1240k > 100k", NA)) |>
  mutate(Quality = if_else(between(`SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)`, 25000,100000), "IndOnly;SNPs 1240k > 25k", Quality)) |>
  mutate(Quality = if_else(`SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)` < 25000, "Alert;SNPs 1240k < 25k", Quality)) |>
  mutate(Quality = if_else(`Molecular Sex` == "U", "Alert;Sex not assigned", Quality)) |>
  mutate(Quality = if_else(ASSESSMENT %in% c("IGNORE"), "Alert;Poor Quality", Quality)) |>
  mutate(Quality = if_else(ASSESSMENT %in% c("QUESTIONABLE", "QUESTIONABLE_CRITICAL"), "Alert;Questionable", Quality)) |>
  # Has biological relationships
  mutate(Quality = if_else(!`Family ID and position within family` %in% c("n/a (no relatives detected)", ".."), sprintf("%s;Has_relatives", Quality), Quality)) |>
  # Is duplicated, use the one with higher SNPs and the other replace first part with usePCA
  group_by(`Master ID`) |>
  mutate(`Usage Note` = if_else(Quality == "Ref;SNPs 1240k > 100k" & `SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)` == max(`SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)`), "Reference", "Other")) |>
  ungroup() 

# Save dataframe as an R data ----
aadr_uk <- df_v0.5
save(aadr_uk, file = "data/aadr_uk.rda")

write_tsv(aadr_uk, "data/AADR_v54.1_Curated_UK.tsv")

rm(list = ls())
#


#####

# ---- 
# Section 2: Read AADR annotation and data curation for global references
# ---- 
# - Read AADR dataset
# - Extract all modern samples 
# - Extract all samples from countries of interest (WE)
# - Curate dataset:
#   - Annotate the Periods based on Patterson et al. 2022 OK
#   - Check the quality of the samples:
#       - N of SNPs
#       - Potential contamination
#       - Is duplicated (select the one with high number of SNPs)
#       - Family annotation (If have relatives get the one with high number of SNPs for that Family)
#   - Annotate usage: PCA Modern, PCA Reference, OutgroupType [Patterson, Olalde], AncientGroup

dataset <- read_tsv(str_c(list_aadr_file$HO, ".anno"))

# prepare data ----
dataset_v0 <- dataset |> 
  mutate(DateMeanInBP = as.numeric(`Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]`)) |>
  filter(DateMeanInBP == 0 | between(DateMeanInBP, 800, 14000)) |>
  filter(`SNPs hit on autosomal targets (Computed using easystats on HO snpset)` > 100000) |>
  filter(ASSESSMENT %in% c("PASS", "PASS_PUBLISHED")) |>
  filter(!str_detect(`Group ID`, "Ignore")) |>
  mutate(`Political Entity` = if_else(`Political Entity` == "Gernamy", "Germany", `Political Entity`))
  
  

political_entity <- c("Abkhazia", "Russia", "Albania", "Armenia", "Austria", #Adygei.HO <- Russia Assyrian <- Armenia, turkey, iran
                      "Bulgaria", "Channel Islands", "Cyprus", "Czech Republic", "Czechoslovakia", 
                      "Croatia", "Denmark", "Estonia", "Finland", "France", 
                      "Georgia", "Germany", "Gibraltar", "Greece", "Hungary", # Spain
                      "Iceland", "Iran", "Iraq", "Ireland", "Isle of Man", "Israel",
                      "Italy", "Libya", "Poland", "Spain", "Morocco", "Tunisia", 
                      "Turkey", "Yemen", "Jordan", "Lebanon", "Lithuania", "Malta",
                      "Moldova", "Netherlands", "Norway", "Portugal", "Romania",
                      "Saudi Arabia", "Serbia", "Slovakia", "Slovenia", "Sweden",
                      "Switzerland", "Syria", "Ukraine", "United Kingdom", "Belarus") 

dataset_v0.1 <- dataset_v0 |>
  filter(`Political Entity` %in% political_entity)

# use the same period names for all the individuals
dataset_v0.1$Period <- case_match(dataset_v0.1$DateMeanInBP,
                             c(14000:5899) ~ "Mesolithic",
                             c(5900:4399) ~ "Neolithic",
                             c(4400:3499) ~ "C/EBA",
                             c(3500:2699) ~ "BA",
                             c(2700:1906) ~ "IA",
                             c(1907:1539) ~ "Romans",
                             c(1540:700) ~ "Early Medieval/Vikings")
 
dataset_v0.2 <- dataset_v0.1 |>
  mutate(Period = if_else(DateMeanInBP == 0, "Present", Period)) 

# remove duplicates based on number of SNPs
dataset_v0.3 <- dataset_v0.2 |>
  group_by(`Master ID`) |>
  arrange(desc(`SNPs hit on autosomal targets (Computed using easystats on HO snpset)`)) |>
  distinct(`Master ID`, .keep_all = T) |>
  ungroup()

# Modern ref ----
# Save present references and manually control
dataset_modern <- dataset_v0.3 |>
  filter(Period == "Present") 

# write_tsv(dataset_modern, file.path(outdir, "AADR_v54.1_modern_references.tsv"), quote = "all")

# read Manually curated dataset and save Modern Reference table
dataset_modern_curated <- read_tsv(file.path(outdir, "ref_items/AADR_v54.1_modern_references_curated.tsv"))

dataset_modern_v0.1 <- dataset_modern |>
  left_join(select(dataset_modern_curated, `Genetic ID`, Reference), by = "Genetic ID") |>
  mutate(Reference = if_else(is.na(Reference), "Not use", Reference))

aadr_modern <- dataset_modern_v0.1
save(aadr_modern, file = file.path(outdir, "aadr_modern.rda"))

# Ancient Ref ----
dataset_ancient <- dataset_v0.3 |>
  filter(Period != "Present") |>
  filter(!`Political Entity` %in% c("United Kingdom", "Isle of Man", "Ireland")) |>
  filter(`Family ID and position within family` %in% c("n/a (no relatives detected)", "..")) |>
  filter(!str_detect(`Group ID`, "_o.{0,}")) |>
  filter(!str_detect(`Group ID`, "son\\.")) |>
  filter(!str_detect(`Group ID`, "dup\\."))  |>
  filter(!str_detect(`Group ID`, "sister\\.")) |>
  filter(!str_detect(`Group ID`, "brother\\."))|>
  filter(!str_detect(`Group ID`, "mother\\."))|>
  filter(!str_detect(`Group ID`, "father\\."))|>
  filter(!str_detect(`Group ID`, "possible"))|>
  filter(!str_detect(`Group ID`, "_noUDG"))|>
  filter(!str_detect(`Group ID`, "_rel\\."))
  
dataset_ancient_v0.1 <- dataset_ancient |>
  mutate(Group = str_remove(`Group ID`, "\\.[DSG]{2}")) 

groups <- dataset_ancient_v0.1 |>
  group_by(Group) |>
  count() |>
  filter(n > 2)

dataset_ancient_v0.2 <- dataset_ancient_v0.1 |>
  filter(Group %in% groups$Group) |>
  mutate(Reference = str_c(`Political Entity`, ";", Group, ";", Period))

# Manually curation
write_tsv(dataset_ancient_v0.2, file.path(outdir, "AADR_v54.1_ancient_references.tsv"))

dataset_curated <- read_tsv(file.path(outdir, "ref_items/AADR_v54.1_ancient_references_curated.tsv"))

aadr_ancient <- dataset_ancient_v0.1 |>
  left_join(select(dataset_curated, `Genetic ID`, Reference, `Note usage`), by = "Genetic ID")


save(aadr_ancient, file = file.path(outdir, "aadr_ancient.rda"))

# outgroups ----
aadr_outgroup <- read_tsv("data/ref_items/AADR_v54.1_outgroup.tsv")
save(aadr_outgroup, file = file.path(outdir, "aadr_outgroup.rda"))

#####
# -----
# # Section 3: Prepare dataset to extract Eigenstrat
# ----
rm(list = ls())

data("aadr_uk")
data("aadr_ancient")
data("aadr_modern")
data("aadr_outgroup")

aadr_uk.sub <- select(aadr_uk, `Genetic ID`, `Group ID`, `Master ID`, `Usage Note`) |>
  mutate(DataRef = "AADR_UK") 

aadr_ancient.sub <- select(aadr_ancient, `Genetic ID`, `Group ID`, `Master ID`,`Note usage`) |>
  rename(`Usage Note` = `Note usage`) |>
  mutate(DataRef = "AADR_Ancient") 

# Filter modern References from file "PCA_modern_references.date.tsv
aux <- read_tsv("data/PCA/PCA_modern_references.Aug_2024.tsv")

aadr_modern.sub <- select(aadr_modern, `Genetic ID`, `Master ID`,`Group ID`) |>
  mutate(`Usage Note` = "Modern")  |>
  mutate(DataRef = "AADR_Modern") |>
  filter(`Group ID` %in% aux$`Group ID`)

aadr_modern2 <- aadr_modern.sub
save(aadr_modern, file = file.path(outdir, "aadr_modern2.rda"))

aadr_outgroup.sub <- select(aadr_outgroup, `Genetic ID`,`Master ID`, `Group ID`) |>
  mutate(`Usage Note` = "Outgroup")  |>
  mutate(DataRef = "AADR_Outgroup") 

aadr_ref <- aadr_uk.sub |>
  add_case(aadr_ancient.sub) |>
  add_case(aadr_modern.sub) |>
  add_case(aadr_outgroup.sub) 
  
aadr_uniques <- aadr_ref

# aadr_uniques <- aadr_ref |>
#   group_by(`Master ID`) |>
#   distinct(`Master ID`, .keep_all = T) |>
#   ungroup()

# -----
# # Section 4: Extract data from AADR HO
# ----

tempdir <- file.path(outdir, "temp")
system2("mkdir", args = sprintf("-p %s", tempdir))
# 
write_tsv(select(aadr_uniques, 1:3), file.path(tempdir, "aadr_uniques.tsv"))
# 
aadr_ho <- update_reference(input = file.path(tempdir, "aadr_uniques.tsv"), output = file.path(tempdir, "aadr_HO_v54.1_sub"),
                            aadr_path = "/Users/lbm536/Documents/COMMIOS_app_preprocessing/raw-data/AADR_v54.1.p1/v54.1.p1_HO_public", by_id = T)

# -----
# # Section 5: PCA
# ----
  # NOTES:
  # PCA performed with all the modern individuals and those that are annotated as References
  # The rest of the samples are projected separately [For i in lowQ do PCA]

# PCA all References ----

## Extract samples to use [Reference (AADR_UK), Use (AADR_Ancinet), all (AADR_Modern)]

out <- file.path(outdir, "PCA")
system2("mkdir", args = sprintf("-p %s", out))

# aadr_uniques <- read_tsv(file.path(tempdir, "aadr_uniques.tsv"))

aadr_pca <- aadr_ref |>
  filter(`Usage Note` %in% c("Modern", "Use", "Reference")) |>
  mutate(Group = if_else(`Usage Note` == "Modern", `Group ID`, `Genetic ID`))

path_aadr_pca <- file.path(out, "aadr_PCA_References.tsv")

write_tsv(aadr_pca, path_aadr_pca)

geno_pca <- update_reference(input = path_aadr_pca, output = file.path(out, "aadr_pca_ref"), 
                             aadr_path = aadr_ho, by_id = T)

aadr_pca_ind <- edit_indfile(input = str_c(geno_pca, ".ind"), table = aadr_pca, output = str_c(geno_pca, "_mod.ind"),
                             by_col = "Genetic ID", remp_col = "Group")


## Prepare for PCA
poplist <- aadr_pca |>
  filter(`Usage Note` == "Modern") |>
  group_by(Group) |>
  distinct(Group) 

write_tsv(poplist, file = "data/PCA/modern.poplist", col_names = F)

run_pca(genotype = geno_pca, poplist = file.path(out, "modern.poplist"), ind = aadr_pca_ind, outdir = out)

# PCA other [lowQ and With Biological Relationships] ----

## Extract samples to use ["Other" (AADR_UK), all (AADR_Modern)]
aadr_pca_other <- aadr_uniques |>
  filter(`Usage Note` %in% c("Modern", "Other")) |>
  mutate(Group = if_else(`Usage Note` == "Modern", `Group ID`, `Genetic ID`))

path_aadr_pca <- file.path(out, "aadr_PCA_Other.tsv")
write_tsv(aadr_pca_other, path_aadr_pca)

geno_pca_other <- update_reference(input = path_aadr_pca, output = file.path(out, "aadr_pca_other"), 
                             aadr_path = aadr_ho, by_id = T)

aadr_pca_other_ind <- edit_indfile(input = str_c(geno_pca_other, ".ind"), table = aadr_pca_other, output = str_c(geno_pca_other, "_mod.ind"),
                             by_col = "Genetic ID", remp_col = "Group")
## Run PCA
run_pca(genotype = geno_pca_other, poplist = file.path(out, "modern.poplist"), ind = aadr_pca_other_ind, outdir = out)


# -----
# # Section 6: qpAdmixture with ancient individuals
# ----
# NOTES:
  # qpAdmixture performed with outgroups from Patterson et al., 2022
  # The admixture is performed to all AADR_UK individuals and AADR_Ancient individuals

# Prepare genotypes 1240K and outgroups ----
out <- file.path(outdir, "qpAdmixture")
system2("mkdir", args = sprintf("-p %s", out))

aadr_admixture <- select(aadr_outgroup, `Genetic ID`, `Master ID`, `Group ID`, Group, Note) |>
  rename(`Usage Note` = Note) |>
  add_case(
    select(aadr_uk, `Genetic ID`, `Master ID`, `Group ID`, `Usage Note`) |>
      mutate(Group = `Genetic ID`) |>
      filter(`Usage Note` != "Outgroup")
  ) |>
  add_case(
    select(aadr_ancient, `Genetic ID`, `Master ID`, `Group ID`, `Note usage`) |>
      rename(`Usage Note` = `Note usage`) |>
      filter(!is.na(`Usage Note`)) |>
      filter(`Usage Note` != "Outgroup") |>
      mutate(Group = `Genetic ID`)
  )
  
path_aadr_qpadm <- file.path(out, "aadr_qpAdm.tsv")
write_tsv(select(aadr_admixture, `Genetic ID`, `Master ID`, `Group ID`), path_aadr_qpadm)

geno_qpadm <- update_reference(input = path_aadr_qpadm, output = file.path(out, "aadr_apadm"), 
                             aadr_path = list_aadr_file$`1240k`, by_id = T)

# Modify Group ID of outgroup to have the same groups than in Patterson et al., 2022 for the admixture analysis ----
write_tsv(aadr_admixture, path_aadr_qpadm)
aadr_qpadm_ind <- edit_indfile(input = str_c(geno_qpadm, ".ind"), table = aadr_admixture, output = str_c(geno_qpadm, ".ind"),
                                   by_col = "Genetic ID", remp_col = "Group")


# Run qpWave ancestry ----

run_qpadm_ancestry <- function(prefix = "", left = "", right = "", target = "") {
  aux <- qpadm(prefix, left, right, target, fudge_twice = T)
  aux_r <- aux
  aux_r$popdrop$target <- target
  aux_r
}

# qpAdmixture ancestry composition of Danebury ----
left <- c("OldSteppe","WHG", "EEF")
right <- c("OldAfrica", "Afanasievo", "Turkey_N", "WHGB")

targets <- filter(aadr_admixture, Group == `Genetic ID`)$`Genetic ID`

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

qpadm_ancestry_reults <- foreach(
  i = 1:length(targets),
  .combine = 'rbind'
) %dopar% {
  library(admixtools)
  library(tidyverse)
  run_qpadm_ancestry(geno_qpadm, left, right, targets[i])
}
# Deactivate parallel
parallel::stopCluster(cl = my.cluster)

write_rds(qpadm_ancestry_reults, "data/apadam.rds")
#----

# -----
# # Section 7: Final table aadr + pca + qpAdm
# ----
# NOTE:
  # Table with unified columns
  # Column DataRef is used to extract information
  # Admixture and PCA are added to the end of the table
  # Locality of AADR_UK will be a combination of the information for that Lat and Long (curated data)

data("aadr_uk")
data("aadr_modern")
data("aadr_ancient")

aadr_uk$DataRef <- "AADR_UK"
aadr_ancient$DataRef <- "AADR_References"
aadr_modern$DataRef <- "AADR_Modern"

# Prepare aadr_uk for the merge ----
aadr_uk.f <- aadr_uk |>
  mutate(Locality = sprintf("%s;%s;%s;%s;%s;", LocID, Latitude, Longitude, CorrectLocalityName, UniqueLocalityName)) |>
  select(!c(LocID, CorrectLocalityName, UniqueLocalityName, DateMeanInBP, label, nSamples)) |>
  mutate(`Usage Note` = sprintf("%s;%s;%s;%s;", `Usage Note`, ASSESSMENT, `ASSESSMENT WARNINGS (Xcontam interval is listed if lower bound is >0.005, "QUESTIONABLE" if lower bound is 0.01-0.02, "QUESTIONABLE_CRITICAL" or "FAIL" if lower bound is >0.02) (mtcontam confidence interval is listed if coverage >2 and upper bound is <0.`,
                          Quality)) |>
  select(!c(`ASSESSMENT WARNINGS (Xcontam interval is listed if lower bound is >0.005, "QUESTIONABLE" if lower bound is 0.01-0.02, "QUESTIONABLE_CRITICAL" or "FAIL" if lower bound is >0.02) (mtcontam confidence interval is listed if coverage >2 and upper bound is <0.`,
            Quality, ASSESSMENT)) |>
  mutate(Reference = sprintf("%s;%s;%s", `Political Entity`, `Group ID`, Period)) 
  

colnames(aadr_uk.f)

# Prepare aadr_modern for the merge ----
colnames(aadr_modern)

aadr_modern.f <- aadr_modern |>
  filter(Reference != "Not use") |>
  mutate(`Usage Note` = "Refernce Modern") |>
  rename(Latitude = Lat.,
         Longitude = Long.) |>
  select(!c(`ASSESSMENT WARNINGS (Xcontam interval is listed if lower bound is >0.005, "QUESTIONABLE" if lower bound is 0.01-0.02, "QUESTIONABLE_CRITICAL" or "FAIL" if lower bound is >0.02) (mtcontam confidence interval is listed if coverage >2 and upper bound is <0.`,
            ASSESSMENT, DateMeanInBP)) |>
  mutate(`SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)` = as.numeric(`SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)`)) |>
  mutate(Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude))

colnames(aadr_modern.f) %in% colnames(aadr_uk.f)

# Prepare aadr_ancient for the merge ----
colnames(aadr_ancient)
aadr_ancient.f <- aadr_ancient |>
  filter(!is.na(`Note usage`)) |>
  rename(`Usage Note` = `Note usage`,
         Latitude = Lat.,
         Longitude = Long.) |>
  select(!c(`ASSESSMENT WARNINGS (Xcontam interval is listed if lower bound is >0.005, "QUESTIONABLE" if lower bound is 0.01-0.02, "QUESTIONABLE_CRITICAL" or "FAIL" if lower bound is >0.02) (mtcontam confidence interval is listed if coverage >2 and upper bound is <0.`,
            ASSESSMENT, DateMeanInBP, Group))  |>
  mutate(`SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)` = as.numeric(`SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)`)) |>
  mutate(Latitude = as.numeric(Latitude), Longitude = as.numeric(Longitude))

colnames(aadr_ancient.f) %in% colnames(aadr_uk.f)

# Merge the tables ----
aadr_db <- aadr_uk.f |>
  add_case(aadr_ancient.f) |>
  add_case(aadr_modern.f) 



# Read admixture, prepare and merge ----

qpadm <- read_rds("data/apadam.rds")
qpadm.f <- NULL
for (i in 1:nrow(qpadm)) {
  aux <-  qpadm[i,1][[1]]
  aux.f <- tibble(`Genetic ID` = unique(aux$target), `qpAdm Pvalue` = qpadm[i,2][[1]]$p[1],
                  Steppe = aux$weight[1], `Steppe SE` = aux$se[1],
                  WHG = aux$weight[2], `WHG SE` = aux$se[2],
                  EEF = aux$weight[3], `EEF SE` = aux$se[3])
  qpadm.f <- if(length(qpadm.f)) add_case(qpadm.f, aux.f) else aux.f
}

# Merge with aadr_db

aadr_db.v0.1 <- aadr_db |>
  left_join(qpadm.f)

# Read PCA, prepare and merge ----
pca.1 <- read_table("data/PCA/aadr_pca_ref.evec", skip = 1, col_names = c("Genetic ID", str_c("PCA", 1:10), "PCA Group ID"))
pca.2 <- read_table("data/PCA/aadr_pca_other.evec", skip = 1, col_names = c("Genetic ID", str_c("PCA", 1:10), "PCA Group ID"))

PCA <- add_case(pca.1, pca.2)

# Merge with aadr_db

aadr_db.v0.2 <- aadr_db.v0.1 |>
  left_join(PCA)


aadr_db <- aadr_db.v0.2

# Simplify colnames and save

save(aadr_db, file = "data/aadr_db.rda")
write_tsv(aadr_db, "data/aadr_db.tsv", quote = "all")
write_rds(aadr_db, file = "data/aadr_db.rds")

# -----
# # Section 8: data cleanning and preparation for app
# ----

data("aadr_db")


# Rename long colnames
df <- aadr_db |>
  rename(`Published Year` = `Year data from this individual was first published [for a present-day individuals we give the data of the data reported here; missing GreenScience 2010 (Vi33.15, Vi33.26), Olalde2018 (I2657), RasmussenNature2010 (Australian)]`) |>
  rename(`Dating method` = `Method for Determining Date; unless otherwise specified, calibrations use 95.4% intervals from OxCal v4.4.2 Bronk Ramsey (2009); r5; Atmospheric data from Reimer et al (2020)`) |>
  rename(`Date Mean in BP` = `Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]`) |>
  rename(`Date SD in BP` = `Date standard deviation in BP [OxCal sigma for a direct radiocarbon date, and standard deviation of the uniform distribution between the two bounds for a contextual date]`) |>
  rename(`Full date` = `Full Date One of two formats. (Format 1) 95.4% CI calibrated radiocarbon age (Conventional Radiocarbon Age BP, Lab number) e.g. 2624-2350 calBCE (3990±40 BP, Ua-35016). (Format 2) Archaeological context range, e.g. 2500-1700 BCE`) |>
  rename(`Age at Death` = `Age at Death from physical anthropology`) |>
  rename(`120k coverage` = `1240k coverage (taken from original pulldown where possible)`) |>
  rename(`SNPs hit on autosomal (1240k snpset)` = `SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)`) |>
  rename(`SNPs hit on autosomal (HO snpset)` = `SNPs hit on autosomal targets (Computed using easystats on HO snpset)`) |>
  rename(`Family ID` = `Family ID and position within family`) |>
  rename(`Y Haplogroup MF` = `Y haplogroup (manual curation in terminal mutation format)`) |>
  rename(`Y Haplogroup ISOGG` = `Y haplogroup (manual curation in ISOGG format)`) |>
  rename(`mtDNA coverage` = `mtDNA coverage (merged data)`) |>
  rename(`mtDNA Haplogroup` = `mtDNA haplogroup if >2x or published`) |>
  rename(`mtDNA match` = `mtDNA match to consensus if >2x (merged data)`) |>
  rename(`Damage rate` = `Damage rate in first nucleotide on sequences overlapping 1240k targets (merged data)`) |>
  rename(`Sex ratio` = `Sex ratio [Y/(Y+X) counts] (merged data)`) |>
  rename(`Library type` = `Library type (minus=no.damage.correction, half=damage.retained.at.last.position, plus=damage.fully.corrected, ds=double.stranded.library.preparation, ss=single.stranded.library.preparation)`)
  

# Separate Locality only for AADR_UK
df_aadr_uk <- df |> 
  filter(DataRef == "AADR_UK") |>
  # Local ID
  separate(Locality, sep = ";", into = c("Locality ID", "LatitudeX", "LongitudeY", "Correct Locality Name", "Unique Locality Name"),
           remove = F, extra = "drop") |>
  select(!c("LatitudeX", "LongitudeY")) |>
  mutate(`Correct Locality Name` = if_else(`Correct Locality Name` == "Roscommon, Kilteasheen", "Ireland, Roscommon, Kilteasheen", 
                                           `Correct Locality Name`)) |>
  mutate(`Unique Locality Name` = if_else(`Unique Locality Name` == "Roscommon, Kilteasheen", "Ireland, Roscommon, Kilteasheen", 
                                          `Unique Locality Name`)) |>
  mutate(Period = if_else(`Genetic ID` %in% c("I3044", "I3045"), "Early Medieval/Vikings", Period))

df_v1.0 <- df |>
  filter(DataRef != "AADR_UK") |>
  # For the merge
  mutate(`Locality ID` = NA, LatitudeX = NA, LongitudeY = NA, `Correct Locality Name` = NA, `Unique Locality Name` = NA) |>
  mutate(`Correct Locality Name` = str_c(`Political Entity`, Locality, sep = ", ")) |>
  mutate(`Correct Locality Name` = if_else(str_detect(`Correct Locality Name`, "France, Yonne, Gurgy"), "France, Yonne, Gurgy", `Correct Locality Name`)) |>
  mutate(`Correct Locality Name` = str_replace_all(`Correct Locality Name`, "_", " ")) |>
  # Join with AADR
  add_case(df_aadr_uk) |>
  # Get Country name
  separate(`Correct Locality Name`, sep = ",", into = c("Country"), remove = F, extra = "drop")  |>
  arrange(Country)


# Correction of Locality names
df_v1.1 <- df_v1.0 |>
  mutate(Country = if_else(Period == "Present" & Country == "United Kingdom", str_remove(Locality, ",.*"), Country)) |>
  mutate(Country = case_match(Country,
                              "Orkney Islands" ~ "Scotland",
                              "Kent" ~ "England",
                              "Cornwall" ~ "England",
                              "Argyll and Bute" ~ "Scotland",
                              ".." ~ "Scotland",
                              "Shetland" ~ "Scotland",
                              .default = Country)) |>
  mutate(`Correct Locality Name` = case_match(`Correct Locality Name`,
                                              "United Kingdom, Scotland, Orkney" ~ "Scotland, Orkney",
                                              "United Kingdom, England, Kent" ~ "England, Kent",
                                              "United Kingdom, England, Cornwall" ~ "England, Cornwall",
                                              "United Kingdom, Orkney Islands" ~ "Scotland, Orkney",
                                              "United Kingdom, Kent" ~ "England, Kent",
                                              "United Kingdom, Cornwall" ~ "England, Cornwall",
                                              "United Kingdom, Argyll and Bute" ~ "Scotland, Argyll and Bute",
                                              "United Kingdom, .." ~ "Scotland, Orkney",
                                              "United Kingdom, England, Cornwall and Devon" ~ "England, Cornwall",
                                              "England,sCounty Durham, Stockton-on-Tees, Windmill Fields" ~ "England, County Durham, Stockton-on-Tees, Windmill Fields" ,
                                              .default = `Correct Locality Name`)) |>
  mutate(`Correct Locality Name` = if_else(!is.na(`Unique Locality Name`) &  `Unique Locality Name` != `Correct Locality Name`, `Unique Locality Name`, `Correct Locality Name`))

# Correction of References, Country;Period;GroupID[Country+Period]
df_v1.2 <- df_v1.1 |>
  mutate(Period = if_else(Period == "Present", "Modern", Period)) |>
  mutate(Reference = str_replace(str_c(str_c(Country, Period, sep = ";"), str_c(Country, Period, sep = "_"), sep = ";"), " ", "_"))


# Simplify Groups ID
aux <- unique(df_v1.2$`Group ID`)


df_v1.3 <- df_v1.2 |>
  mutate(Group = str_remove(`Group ID`, "\\.(HO|SG|DG|WGA)")) |>
  mutate(Group = str_remove(Group, "(_o.*|_lc|_1d.*)")) |>
  mutate(Group = str_remove(Group, "(_oEEF|_high.*|_medium.*|_low.*)")) |>
  mutate(Group = str_remove(Group, "(_daughter.*|_father.*|\\.rel\\..*|_mother.*|_brother.*|_sister.*|_noUDG.*|_son.*|_possible.1d)")) |>
  mutate(Group = str_remove(Group, "\\.(HO)"))|>
  mutate(Group = str_remove(Group, "(_contam)")) |>
  mutate(Group = case_match(Group,
                            "IBS" ~ "Spanish",
                            "GBR" ~ "British",
                            "Orcadian" ~ "Scottish",
                            .default = Group))





data_ref <- df_v1.3

save(data_ref, file = "data/data_ref.rda")


# patch _ change british group

data("data_ref")
data_ref <- data_ref |>
  mutate(Group = if_else(Group == "British", 
                         if_else(Country == "England", "English", "Scottish"), Group))

save(data_ref, file = "data/data_ref.rda")


############### 
# modern <- update_reference(input = file.path(getwd(), "raw-data","modern_pca.samples"),
#                                   output = file.path(getwd(), "data", "modern_reference"),
#                                   aadr_path = file.path(getwd(), "raw-data/v54.1_HO_public"))
# 
# #############################################################
# # Section 2: Update projected dataset and merge with modern #
# #############################################################
# 
# projected <- update_reference(input = file.path(getwd(), "raw-data","projected_pca.samples"),
#                                      output = file.path(getwd(), "data", "project_reference"),
#                                      aadr_path = file.path(getwd(), "raw-data/v54.1_HO_public"), by_id = T)
# 
# #############################
# # Section 3: merge datasets #
# #############################
# 
# # Edit indfile modern
# modern_anno <- read_tsv(file.path("raw-data/modern_pca.anno"))
# modern <- file.path(getwd(), "data/modern_reference")
# 
# modern_ind_file <- edit_indfile(input = str_c(modern, ".ind"), table = modern_anno, 
#                                 output = file.path(getwd(), "data/modern_reference_mod.ind"), 
#                                 by_col = "Genetic ID", remp_col = "Name")
# 
# # Edit indfile projection
# projected_anno <- read_tsv(file.path("raw-data/projected_pca.anno"))
# projected <- file.path(getwd(), "data/project_reference")
# projected_ind_file <- edit_indfile(input = str_c(projected, ".ind"), table = projected_anno, 
#                                 output = str_c(projected, "_mod.ind"), 
#                                 by_col = "Genetic ID", remp_col = "Genetic ID")
# 
# # Merge datasets and annotation
# merge_pca_ref <- merge2genetypes(geno1 = modern, geno2 = projected, outdir = file.path(getwd(), "data"), 
#                                  outname = "modern_projected_pca", 
#                                  ind1 = modern_ind_file, ind2 = projected_ind_file)
# 
# #####################################################
# # Section 4: Update UK dataset HO and merge for pca #
# #####################################################
# 
# 
# UK_projected <- update_reference(input = file.path(getwd(), "raw-data","uk_projected_pca.samples"),
#                               output = file.path(getwd(), "data", "uk_project_reference"),
#                               aadr_path = file.path(getwd(), "raw-data/v54.1_HO_public"), by_id = T)
# 
# all_pca_ref <- merge2genetypes(geno1 = merge_pca_ref, geno2 = UK_projected, outdir = file.path(getwd(), "data"), 
#                                  outname = "all_pca")
# 
# ######################
# # Section 5: run pca #
# ######################
# 
# pca_output <- file.path(getwd(),"results/PCA_reference/")
# 
# tmp_modern_pop <- file.path(getwd(),"data/modern.population")
# write.table(unique(modern_anno$Name), tmp_modern_pop, col.names = F, quote = F, row.names = F)
# 
# run_pca(genotype = all_pca_ref,
#   poplist = tmp_modern_pop,
#   numthreads = 8, outdir = pca_output)
# 
# smartPCA <- read_table(file.path("data/all_pca.evec"), skip = 1, col_names = F)
# colnames(smartPCA) <- c("Genetic ID", 
#                         str_c("PCA", rep(1:(ncol(smartPCA)-2))), "GroupID")
# 
# write_tsv(smartPCA, file.path(pca_output, "smartPCA_output.tsv"))
# 
# #######################
# # Section 5: plot pca #
# #######################
