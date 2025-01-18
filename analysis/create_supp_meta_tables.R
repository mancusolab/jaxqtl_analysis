# supplemental tables

tables_dir <- "../result/cis/tables/"

# cell type library size
lib_size <- read_tsv("../data/OneK1K/pheno/celltype16_new/metadata/celltype_14_total_libsize.tsv")
summary(lib_size$libsize)
max(lib_size$libsize)

lib_size %>% 
  arrange(desc(libsize)) %>% 
  rename(total_UMI = libsize) %>% 
  xlsx::write.xlsx(file = paste0(tables_dir, "Tables.xlsx"),
                   sheetName = "TableS1", append = FALSE)

# write cell level information
genes_pass_tab <- NB_all_df %>% filter(express_percent >= 0.01) %>% 
  group_by(celltype) %>% 
  summarize(n_genes_pass = n()) %>% 
  ungroup()

table_s1 <- ind_celltype_ct %>% 
  filter(counts>0) %>% 
  group_by(celltype) %>% 
  summarize(total_cell = sum(counts),
            n = n()) %>% 
  ungroup() %>% 
  left_join(genes_pass_tab, by="celltype") %>% 
  left_join(lib_size %>% 
              arrange(desc(libsize)) %>% 
              rename(total_UMI = libsize))

table_s1 %>% 
  arrange(desc(total_cell)) %>% 
  xlsx::write.xlsx(file = paste0(tables_dir, "Tables.xlsx"),
                              sheetName = "TableS1", append = FALSE)

tables_dir
