#####
library(httr)

download_aadr_files <- function(save_path = file.path(), download_all = list()){
  
  # Validate input
  if (missing(save_path)) {
    stop("save path must be provided.")
  }
  if (missing(download_all)) {
    stop("list with download path must be provided.")
  }
  
  
  if(!dir.exists(save_path)){
    
    dir.create(save_path)
    cat("The folder was created.\n")
    
  } else {
    
    cat(sprintf("Files to be saved into %s", save_path))
  }
  
  # Download files
  for (i in download_all) {
    
    url <- i
    
    destfile <- file.path(save_path, basename(i))
    
    # Download the file
    response <- GET(url, write_disk(destfile, overwrite = TRUE))
    
    # Check the status of the download
    if (status_code(response) == 200) {
      cat("File downloaded successfully to", destfile, "\n")
    } else {
      cat("Failed to download file. Status code:", status_code(response), "\n")
    }
  }
  
}



#####

extract_from_genotype <- function(genotype = file.path(), outname = NULL, 
                                  outdir = file.path(), ind = NULL, 
                                  poplist = file.path(), convertf = "/Users/lbm536/Software/EIG/src/convertf"){
  # Extract samples from a genotype [save poplist in a temp file]
  
 # Validate input
  if (missing(genotype)) {
    stop("genotype file path must be provided.")
  }
  if (missing(poplist)) {
    stop("poplist file path must be provided.")
  }
  
  # Construct output paths
  out <- if (is.null(outname)) {
    file.path(outdir, paste0(basename(genotype), "_substract"))
  } else {
    file.path(outdir, basename(outname))
  }
  
  indfile <- if (is.null(ind)) genotype else ind
  
  parameters <- file.path(paste0(out, ".par"))
  
  # Debug prints
  print(paste("genotype:", genotype))
  print(paste("indfile:", indfile))
  print(paste("out:", out))
  print(paste("parameters:", parameters))
  print(paste("poplist:", poplist))

  # Write parameters to file
  cat(
    sprintf("
            genotypename: %s.geno
            snpname: %s.snp
            indivname: %s.ind
            outputformat:    PACKEDANCESTRYMAP
            genotypeoutname: %s.geno
            snpoutname: %s.snp
            indivoutname: %s.ind
            poplistname: %s",
            genotype, genotype, indfile,
            out, out, out,
            poplist),
      file = parameters
  )
 # Run the convertf command
 system2(command = convertf, args = c("-p ", parameters))
  
 unlink(parameters)
 # Return path of genotype
 out
}

#####

update_reference <- function(input = file.path(), output = file.path(), 
                             aadr_path = file.path(), by_id = F){
  # Update dataset function
  
  input_file <- read_tsv(input)
  colnames(input_file) <-  c("geneticid", "masterid", "groupid")
  
  tmp_poplist <- file.path(dirname(input), "poplist.tmp")
  
  # if we want to extract only selected individuals
  if(by_id){
    tmp_indfile <- file.path(dirname(output), "tmp.ind")
    ind_file <- read.table(str_c(aadr_path, ".ind"), col.names = c("geneticid", "sex", "groupid"))
    ind_file <- mutate(ind_file, groupid = geneticid)
    
    # save individuals temp file
    write_tsv(ind_file, tmp_indfile, col_names = F)
    
    # save poplist temp file
    write.table(unique(input_file$geneticid), tmp_poplist, col.names = F, quote = F, row.names = F)
    
  } else {
    # save poplist temp file
    tmp_indfile <- str_c(aadr_path, ".ind")
    write.table(unique(input_file$groupid), tmp_poplist, col.names = F, quote = F, row.names = F)
  }
  timeout(2)
  outfile <- extract_from_genotype(genotype = aadr_path, 
                               outname = basename(output), 
                               outdir = dirname(output),
                               ind = if_else(by_id, str_replace(tmp_indfile, ".ind", ""), aadr_path), poplist = tmp_poplist)
  
  unlink(tmp_poplist)
  unlink(tmp_indfile)
  outfile
}

#######

edit_indfile <- function(input = file.path(), table = tibble(), output = file.path(), by_col = NULL, remp_col = NULL) {
  ind_file <- read_table(input, col_names = c("Genetic ID", "Sex", "GroupInd"))
  ind_new <- ind_file |>
    left_join(table, by = setNames("Genetic ID", by_col)) |>
    mutate(GroupInd = .data[[remp_col]]) |>
    select(`Genetic ID`, Sex, GroupInd)
  
  write_tsv(ind_new, output, col_names = F)
  output
}

######

merge2genetypes <- function(geno1 = file.path(), geno2 = file.path(), outdir = file.path(), outname = str(), ind1 = NULL, ind2 = NULL){
  # The inputs are the dataset 1, dataset 2, path to output directory, the name of the output files and the 
  # ind files if they are different from the original  
  
  out <- file.path(outdir, outname)
  
  ind1 <- ifelse(!is.null(ind1), str_remove(ind1, ".ind"), geno1)
  ind2 <- ifelse(!is.null(ind2), str_remove(ind2, ".ind"), geno2)
  
  tmp_par <- str_c(out, ".par")
  # Merge data with database - will create a parameter file in the directory of data
  cat(sprintf(
    "
  geno1: %s.geno
  snp1: %s.snp
  ind1: %s.ind 
  geno2: %s.geno
  snp2: %s.snp
  ind2: %s.ind
  genooutfilename: %s.geno
  snpoutfilename: %s.snp
  indoutfilename: %s.ind", 
    geno1, geno1, ind1, 
    geno2, geno2, ind2, 
    out, out, out), file = tmp_par)
  
  # Run mergeit
  system2(command = "/Users/lbm536/Software/EIG/src/mergeit", 
          args = str_c("-p ", str_c(out, ".par")))
  
  unlink(tmp_par)
  out
}

######

run_pca <- function(genotype = file.path(), poplist = file.path(), lsqproject = "YES",
                    autoshrink = "YES", numoutvec = 4, numthreads = 4,
                    outdir = file.path(), ind = NULL,
                    smartpca = "/Users/lbm536/Software/EIG/src/eigensrc/smartpca"){
  
  ind <- ifelse(!is.null(ind), str_remove(ind, ".ind"), genotype)
  
  tmp_par <- file.path(outdir, "smartpca.par")
  
  # Parameters
  cat(sprintf("
  genotypename: %s.geno
  snpname: %s.snp
  indivname: %s.ind
  evecoutname: %s.evec
  evaloutname: %s.eval
  poplistname: %s
  lsqproject: %s
  autoshrink: %s
  numoutvec: %d
  numthreads: %d", 
              genotype, genotype, ind, 
              genotype, genotype, poplist, 
              lsqproject, autoshrink, 
              numoutvec, numthreads), file = tmp_par)
  
  
  
  # # Run smartPCA
  system2(command = smartpca,
          args = c("-p ", tmp_par))
  
  unlink(tmp_par)
}

# 
# run_pca(genotype = all_pca_ref,
#         poplist = tmp_modern_pop,
#         numthreads = 8, outdir = pca_output)

######

# Quality Index for the samples
# Based on:
# - SNPs hit on autosomal (1240k snpset) scale to [0-10].
# - Sex ratio
# - Damage rate
# - ANGSD
# - hapConX
# - ASSESSMENT
# - endogenous by library

calculate_quality_index <- function(df, by) {
  
  df <- df |>
    # SNPs score
    mutate(SNPs_score = pmin({{by}} / 100000 * 10, 10)) |>
    # Sex ratio not meassured
    mutate(Sex_score = case_when(
      `Sex ratio` == ".." ~ 2,
      TRUE ~ 10
    )) |>
    # Assessment score
    mutate(Assessment_score = case_when(
      `ASSESSMENT` == "PASS" ~ 10,
      `ASSESSMENT` == "MERGE_PASS" ~ 8,
      `ASSESSMENT` == "QUESTIONABLE" ~ 3,
      `ASSESSMENT` == "CRITICAL" ~ 0,
      TRUE ~ 2
    )) |>
    # mtDNA coverage
    mutate(`mtDNA coverage` = as.numeric(`mtDNA coverage`))|>
    mutate(mtDNA_coverage_score = case_when(
      `mtDNA coverage` >= 30 ~ 10,
      `mtDNA coverage` >= 15 & `mtDNA coverage` < 30 ~ 7,
      `mtDNA coverage` >= 5 & `mtDNA coverage` < 15 ~ 5,
      `mtDNA coverage` >= 1 & `mtDNA coverage` < 5 ~ 2,
      TRUE ~ 0
    ))
  
  df <- df |>
    mutate(Quality_score = round(
      0.5 * SNPs_score +
        0.1 * Sex_score +
        0.3 * Assessment_score +
        0.1 * mtDNA_coverage_score, digits = 3)) |>
    select(!c(SNPs_score, Sex_score, Assessment_score, mtDNA_coverage_score))
  
  return(df)
}

# Run qpWave ancestry ----

run_qpadm_ancestry <- function(prefix = "", left = "", right = "", target = "") {
  aux <- qpadm(prefix, left, right, target, fudge_twice = T)
  aux_r <- aux
  aux_r$popdrop$target <- target
  aux_r
}

######
# 
# is_in_europe <- function(lat, long) {
#   # Define the bounding box for Europe
#   europe_bbox <- st_bbox(c(xmin = -31.266001, ymin = 34.5428, xmax = 39.869301, ymax = 81.008797), crs = st_crs(4326))
#   
#   # Define the polygon for Europe
#   europe_polygon <- st_as_sfc(st_bbox(europe_bbox))
#   
#   # Create a point from the given latitude and longitude
#   point <- st_sfc(st_point(c(long, lat)), crs = st_crs(4326))
#   
#   # Check if the point is within the Europe polygon
#   is_inside <- st_contains(europe_polygon, point, sparse = FALSE)
#   
#   return(is_inside)
# }
# 
# # # Example usage
# # lat <- 48.8566  # Latitude for Paris
# # long <- 2.3522  # Longitude for Paris
# # x <- is_in_europe(lat, long)