library(dplyr)
library(tidyr)
library(purrr)
# Read all .txt files in a folder and combine them into a single list
# The list will contain all the data from the .txt files
# The list will be saved as a .rds file
dir <- "/Users/nirwantandukar/Documents/Research/results/GWAS/SAP/SoLD/lowinput/annotation"
files <- list.files(dir, pattern = ".txt", full.names = TRUE)

# Save as list for each file
# Function to extract the base name (e.g., TG(60:5)) from file name
get_lipid_name <- function(path) {
  name <- basename(path)
  clean <- sub("_mod.*$", "", name)           # remove from _mod onward
  clean <- sub("_annotation\\.txt$", "", clean)  # remove trailing _annotation.txt
  return(clean)
}

# Read all files into a named list
# Read and extract only selected summary columns
# Create named list of summarized data frames
annotation_list <- setNames(
  lapply(files, function(file) {
    df <- vroom::vroom(file, show_col_types = FALSE)
    if (nrow(df) == 0) return(NULL)
    
    df %>%
      mutate(`log(p)` = -log10(MinPValue)) %>%
      ### select MinPVlaue <= 0.00001 
      dplyr::filter(MinPValue <= 0.00001) %>%
      dplyr::select(GeneID, Chromosome,`log(p)`)
      
  }),
  sapply(files, get_lipid_name)
)


annotation_list[["TG(60:5)"]]

# Optional Save to RDS
saveRDS(annotation_list, "all_annotations_lowinput.rds")




### RENAMING THE MANHATTAN PLOTS

# Set path to your folder
manhattan_dir <- paste0(getwd(),"/figures/manhattans/lowinput")
getwd()
# List all manhattan image files
files <- list.files(manhattan_dir, pattern = "_manhattan\\.jpg$", full.names = TRUE)

# Go through each file
for (file in files) {
  original_name <- basename(file)
  
  # Remove prefix and _mod suffix
  clean_name <- original_name |>
    gsub("^Rect_Manhtn\\.", "", x = _) |>
    gsub("_mod.*", "", x = _)
  
  # Add .jpg extension
  new_path <- file.path(manhattan_dir, paste0(clean_name, ".jpg"))
  
  # Rename
  file.rename(from = file, to = new_path)
}


