# download BLUEPRINT results from API
# ref: https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tutorials/API_v2/eQTL_API_tutorial.md
library("tidyverse")
library("httr")
library("glue")

library("jsonlite")
library("ggrepel")
library("qvalue")
library("optparse")

setwd("/project/nmancuso_8/data/eQTL-Catalogue/allsumstats")

option_list <- list(
  make_option(c("--i"), type="integer", default=NULL, 
              help="row index of metadata")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Change parameters
max_pulled_rows <- 1000 #All datasets will be pulled if this parameter is bigger than the actual number of datasets

URL <- glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size={max_pulled_rows}")

# Make a request
r <- GET(URL, accept_json())
# Check status
status_code(r)

# Extract content
cont <- content(r, "text", encoding = "UTF-8")
# Convert content to dataframe
datasets <- fromJSON(cont)

table(datasets$study_label)

# write a function to extract association summary stats
request_associations_from_api <- function(
    dataset_id, 
    pos="",
    variant="", 
    rsid="",
    molecular_trait_id="",
    gene_id="",
    nlog10p=""){
  
  size = 1000
  start = 0
  
  parameter_values = c(dataset_id,pos,variant,rsid,molecular_trait_id, 
                       gene_id,nlog10p)
  parameter_names = c('dataset_id','pos','variant','rsid','molecular_trait_id', 
                      'gene_id','nlog10p')
  
  while (T) {
    URL = glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/{dataset_id}/associations?size={size}&start={start}")
    
    #Adding defined parameters to the request
    for (i in 1:length(parameter_values)) {
      par = parameter_values[i]
      par_name = parameter_names[i]
      if (par != "")
        URL = glue("{URL}&{par_name}={par}")
    }
    
    r <- GET(URL, accept_json())
    cont <- content(r, "text", encoding = "UTF-8")
    
    # If the request was unsuccessful
    if (status_code(r) != 200) {
      #If we get no results at all, print error
      if (start == 0) {
        print(glue("Error {status_code(r)}"))
        print(cont)
        return ()
      }
      #else just break
      break
    }
    
    cont_df <- fromJSON(cont)
    
    if (start == 0) {
      responses <- cont_df
    }
    else{
      responses <- rbind(responses, cont_df)
    }
    start <- start + size
  }
  return(responses)
}

datasets <- datasets %>% filter(quant_method == "ge" & condition_label == "naive") %>% 
  filter(study_label %in% c("BLUEPRINT", # isolated cells
                            "GENCORD", # cord blood T cells
                            "Schmiedel_2018", # DICE: multiple CT
                            "Perez_2022", # SLE and healthy control patients
                            "Nathan_2022" # unstimulated cells
                            )) %>% 
  mutate(group = case_when(tissue_label %in% c("CD4+ T cell", "T cell", "CD8+ T cell", "CD4+ memory T cell") ~ "T_cell",
                           tissue_label %in% c("monocyte", "CD16+ monocyte") ~ "monocyte",
                           tissue_label %in% c("B cell") ~ "B_cell",
                           tissue_label %in% c("NK cell") ~ "NK",
                           tissue_label %in% c("dendritic cell") ~ "DC",
                           .default = NA)) %>% 
  filter(!is.na(group)) %>% 
  mutate(out_name = glue("{study_label}_{sample_group}"))

#datasets %>% write_tsv("CT_metadata")

datasets <- read_tsv("CT_metadata")
i <- as.integer(opt$i)

id <- datasets$dataset_id[i]
outname <- datasets$out_name[i]

assoc <- request_associations_from_api(dataset_id=id)

assoc %>% data.table::fwrite(glue("{outname}.allpairs.tsv.gz"), sep="\t")


