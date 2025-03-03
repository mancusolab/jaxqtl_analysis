indir <- "../result/replication"
plots_dir
pi1_bonf <- "fdr" # or "bonf", "pi1", "fdr"

if (pi1_bonf %in% c("bonf", "fdr")){
  ylab_name <- "Replication rate"
}else {
  ylab_name <- bquote(pi[1] ~ "statistics")
}

#### Perez jaxqtl results ####

prefix <- "Perez_2022_Asian_EUR_jaxqtl_nb"
express_cutoff <- 0.01
perez_summary <- read_tsv(glue("{indir}/{prefix}_express{express_cutoff}_rep_summary.tsv")) %>% 
  mutate(bonf_rep = bonf_rep / num_snp_found)

pop_cols <- c("Asian" = "#9972AF", "European" = "#D55E00")
perez_order <- c("T4", "T8", "B", "NK", "cM", "ncM", "cDC")

if (pi1_bonf == "pi1"){
  perez_summary["rep_rate"] <- perez_summary$pi1
}else if (pi1_bonf == "bonf"){
  perez_summary["rep_rate"] <- perez_summary$bonf_rep
}else {
  perez_summary["rep_rate"] <- perez_summary$fdr
}

perez_summary %>% arrange(pop, rep_rate) %>% View

perez_summary %>% group_by(pop) %>% 
  summarize(n_eqtl = sum(num_snp),
            found_eqtl = sum(num_snp_found))

p_perez <- perez_summary %>% 
  mutate(rep_CT = fct_relevel(rep_CT, perez_order),
         yazar_CT = fct_reorder(yazar_CT, desc(rep_rate)),
         population = fct_relevel(pop, "European")) %>% 
  ggplot(aes(x = yazar_CT, y = rep_rate, fill = population)) +
  geom_bar(position = position_dodge(width=0.5), stat="identity", width = 0.5) +
  scale_fill_manual(values = pop_cols) +
  facet_grid(.~rep_CT, scales = "free_x", space = "free") +
  scale_x_discrete(
    labels = celltype_bquote
  )+
  xlab("") + ylab(ylab_name) + ylim(0, 1) +
  axis_theme+
  theme(axis.text.x = element_text(angle=35, hjust = 1),
        legend.position = c(.18, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.title = element_blank(),
        legend.text=element_text(size=7),
        legend.margin = margin(1, 1, 1, 1),
        legend.key.size = unit(0.2, units = "cm"))
p_perez

ggsave2(glue("{plots_dir}/Rep_jaxqtl_perez_Asian_EUR_{pi1_bonf}_express{express_cutoff}.png"),
        width = full_width, height = full_height*1.2, units = 'cm', dpi = 300)

perez_res <- read_tsv(glue("{indir}/{prefix}_express{express_cutoff}_rep.tsv.gz")) %>% 
  mutate(slope = ifelse(sign_flip == TRUE, slope * -1, slope)) %>% 
  mutate(ambiguous_snp = ifelse(ref == "A" & alt == "T" | ref == "T" & alt == "A" |
                                  ref == "G" & alt == "C" | ref == "C" & alt == "G", TRUE, FALSE)) %>% 
  select(-express_percent_perez) %>% 
  left_join(NB_egenes_df %>% select(phenotype_id, celltype, slope_se_wald),by=c("celltype", "phenotype_id")) %>% 
  left_join(celltype_N %>% rename("N_onek1k"="N"), by=c("celltype")) %>%
  mutate(N_perez = 88)

axis_window <- 2
ct_toplot <- c("CD4_NC", "CD8_ET", "NK", "B_IN", "Mono_C", "DC")
CT.lab3 = c(
  `CD4_NC` = bquote(CD4[NC]), 
  `CD8_ET` = bquote(CD8[ET]), 
  `NK` = bquote("NK"),
  `B_IN` = bquote(B[IN]),
  `Mono_C` = bquote(Mono[C]),
  `DC` = bquote("DC")
)

vlabeller <- function (variable, value) {
  return(CT.lab3[value])
}

p_perez_slope_EUR <- perez_res %>% 
  filter(pop == "European" & celltype %in% ct_toplot) %>%
  mutate(slope_adj_onek1k = (slope_wald / slope_se_wald) / sqrt(N_onek1k),
         slope_adj_perez = (slope / slope_se) / sqrt(N_perez)) %>% 
  mutate(paint_celltype = ifelse(qval >= 0.05, "null", celltype),
         sig = ifelse(qval >= 0.05, "no", "yes")) %>% 
  mutate(celltype = fct_relevel(celltype, celltype_order)) %>% 
  # filter(abs(slope_adj_onek1k) < axis_window & abs(slope_adj_perez) < axis_window) %>%
  arrange(sig) %>% 
  ggplot(aes(x = slope_adj_onek1k, y = slope_adj_perez)) + 
  geom_point(aes(color = paint_celltype), size = 0.5) +
  geom_smooth(method = "lm", linewidth=0.5,
              aes(color = celltype), fill = "grey") +
  # facet_wrap(.~celltype) +
  facet_grid(.~celltype, labeller = vlabeller) +
  geom_hline(yintercept = 0, color = "black", linewidth=0.1) + 
  geom_vline(xintercept = 0, color = "black", linewidth=0.1) +
  xlab("OneK1K estimates") + ylab("CLUES European estimates") +
  scale_color_manual(values = c(celltype_cols, "null"="grey")) + 
  labs(colour="Cell type") +
  xlim(-axis_window, axis_window) + ylim(-axis_window, axis_window) + 
  axis_theme; p_perez_slope_EUR + theme(legend.position = "none")

p_perez_slope_Asian <- perez_res %>% 
  filter(pop == "Asian" & celltype %in% ct_toplot) %>%
  mutate(slope_adj_onek1k = (slope_wald / slope_se_wald) / sqrt(N_onek1k),
         slope_adj_perez = (slope / slope_se) / sqrt(N_perez)) %>% 
  mutate(paint_celltype = ifelse(qval >= 0.05, "null", celltype),
         sig = ifelse(qval >= 0.05, "no", "yes")) %>% 
  mutate(celltype = fct_relevel(celltype, celltype_order)) %>% 
  # filter(abs(slope_adj_onek1k) < axis_window & abs(slope_adj_perez) < axis_window) %>%
  arrange(sig) %>% 
  ggplot(aes(x = slope_adj_onek1k, y = slope_adj_perez)) + 
  geom_point(aes(color = paint_celltype), size = 0.5) +
  geom_smooth(method = "lm", linewidth=0.5,
              aes(color = celltype), fill = "grey") +
  # facet_wrap(.~celltype) +
  facet_grid(.~celltype, labeller = vlabeller) +
  geom_hline(yintercept = 0, color = "black", linewidth=0.1) + 
  geom_vline(xintercept = 0, color = "black", linewidth=0.1) +
  xlab("OneK1K estimates") + ylab("CLUES Asian estimates") +
  scale_color_manual(values = c(celltype_cols, "null"="grey")) + 
  labs(colour="Cell type") +
  xlim(-axis_window, axis_window) + ylim(-axis_window, axis_window) + 
  axis_theme; p_perez_slope_Asian + theme(legend.position = "none")

plot_grid(
  p_perez,
  p_perez_slope_EUR + theme(legend.position = "none"), 
  p_perez_slope_Asian + theme(legend.position = "none"), 
  align = 'v',
  hjust = -1, # -1
  nrow = 3,
  rel_heights = c(1, 0.8, 0.8), rel_widths=c(1, 1), axis="l",
  labels = c("A", "B"), label_size = 8
)

ggsave2(glue("{plots_dir}/Main_Rep_jaxqtl_perez_slope_trunc_express{express_cutoff}.png"),
        width = full_width, height = full_height*2, units = 'cm', dpi = 300)

