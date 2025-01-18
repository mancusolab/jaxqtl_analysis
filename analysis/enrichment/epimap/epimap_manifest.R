# Epimap
# annot dir: /project/gazal_569/artem/Pipeline_v0.2/data/beds/Epimap_eleanor/

library(tidyverse)

epimap <- readxl::read_excel("../data/OneK1K/annotation/epimap/Annots_final_v5.xlsx") %>% 
  filter(Source == "Epimap_v3")

mydir <- "/project/gazal_569/artem/Pipeline_v0.2/data/beds/Epimap_eleanor"

table(epimap$Name)

epimap %>% View()

table(epimap$Cell.Type)

epimap %>% filter(grepl("monocyt", Name)) %>% View

celltype_list <- c("B cells", "CD4+ T cells", 
                   "CD8+ T cells", "Dendritic cells", "NK cell", "T cells", 
                   "Mononuclear Phagocytes", "Mono_C")

manifest <- epimap %>% 
  mutate(Cell.Type = ifelse(grepl("monocyte", Name), "Mono_C", Cell.Type)) %>% 
  filter(Cell.Type %in% celltype_list) %>% 
  filter(!grepl("^Cancer", Name)) %>% 
  mutate(Path = gsub("/project/gazal_569/artem/Pipeline_v0.2/data/annotations/SLDSC/Epimap_v3",
                     mydir, Path)) %>% 
  mutate(Cell.Type.supergroup = case_when(Cell.Type == "B cells" ~ "Bcells",
                                          Cell.Type == "CD4+ T cells" ~ "CD4_Tcells",
                                          Cell.Type == "CD8+ T cells" ~ "CD8_Tcells",
                                          Cell.Type == "Dendritic cells" ~ "DC",
                                          Cell.Type == "NK cell" ~ "NK",
                                          Cell.Type == "T cells" ~ "Tcells",
                                          Cell.Type == "Mononuclear Phagocytes" ~ "PBMC",
                                          .default = as.character(Cell.Type))) %>% 
  arrange(Cell.Type.supergroup) 

manifest %>% 
  select(-c(Subgroup, Additional.Group, `Embryonic.&.Stem.cells`, X11)) %>% 
  rename(group = `Cell.Type.supergroup`) %>% 
  xlsx::write.xlsx("../data/OneK1K/annotation/epimap/epimap_manifest.xlsx")

manifest %>% 
  select(Name, Cell.Type.supergroup, Path) %>% 
  arrange(Cell.Type.supergroup) %>% pull(Cell.Type.supergroup) %>% table()
  write_tsv("./analysis/enrichment/epimap/params_bed", col_names = F)
  
epimap %>% 
    mutate(Cell.Type = ifelse(grepl("monocyte", Name), "Mono_C", Cell.Type)) %>% 
    filter(Cell.Type %in% celltype_list) %>% 
    filter(!grepl("^Cancer", Name)) %>% 
    mutate(Path = gsub("/project/gazal_569/artem/Pipeline_v0.2/data/annotations/SLDSC/Epimap_v3",
                       mydir, Path)) %>% 
    mutate(Cell.Type.supergroup = case_when(Cell.Type == "B cells" ~ "Bcells",
                                            Cell.Type == "CD4+ T cells" ~ "CD4_Tcells",
                                            Cell.Type == "CD8+ T cells" ~ "CD8_Tcells",
                                            Cell.Type == "Dendritic cells" ~ "DC",
                                            Cell.Type == "NK cell" ~ "NK",
                                            Cell.Type == "T cells" ~ "Tcells",
                                            Cell.Type == "Mononuclear Phagocytes" ~ "PBMC",
                                            .default = as.character(Cell.Type))) %>% 
    select(Name, Cell.Type.supergroup, Path) %>% 
  distinct(Cell.Type.supergroup, .keep_all = T) %>% 
  arrange(Cell.Type.supergroup) %>% 
    write_tsv("./analysis/enrichment/epimap/params_onesample_bed", col_names = F)
  