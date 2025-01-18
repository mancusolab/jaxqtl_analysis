## eqtl sharing characteristics

#### TSS for shared vs. sp eqtl ####
# use all points to test distance ~ number of shared cell type
glm(n ~ dist, data = snp_df %>% count(eqtl, dist), 
    family = "poisson") %>% 
  broom::tidy()

p_share_dist <- snp_df %>%
  count(eqtl, dist) %>% 
  group_by(n) %>% 
  summarize(mean = mean(abs(dist))) %>%
  ggplot(aes(x = n, y = mean)) + 
  geom_point() +
  scale_x_continuous(breaks=seq(1,14,1))+
  geom_smooth(method = "lm") + 
  ylab("mean distance to TSS") + xlab("number of shared cell types") + axis_theme
p_share_dist

#### cell type proportion vs. # sp eqtl ####
n1_eqtl_cellcts <- specific_eqtl_raw %>% 
  count(celltype) %>% # pull(n) %>% sum()
  left_join(cell_counts, by="celltype") %>% 
  mutate(prop = total_cell / sum(total_cell))

p_ct_spec <- n1_eqtl_cellcts %>% 
  ggplot(aes(x = prop, y = n)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  xlab("Cell type proportion") + ylab("count of cell-type-specific eQTLs")+
  axis_theme

fit <- glm(n ~ prop, data=n1_eqtl_cellcts, family = "gaussian") # "poisson", "gaussian"
broom::tidy(fit)

# linear model R2: 91%
# poisson model R2: 87%
1-sum((fit$data$n - fit$fitted.values)^2)/sum((fit$data$n - mean(fit$data$n))^2)
fit$R


#### upset plot for 2256 eQTL with PIP > 0.5 ####
lt <- list()
for (cell in celltypes){
  lt[[cell]] <- unique(snp_df$eqtl[which(snp_df$celltype == cell)])
}

m <- make_comb_mat(lt)
m2 <- m[, comb_size(m) >= 10] # for nb_cs_pip0.9

p_upset <- grid.grabExpr(draw(UpSet(m2,
                                    comb_order = order(comb_size(m2), decreasing = TRUE),
                                    set_order = c("CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET", "CD8_S100B",
                                                  "NK", "NK_R",
                                                  "B_Mem", "B_IN", "Plasma",
                                                  "Mono_C","Mono_NC", "DC"),
                                    top_annotation = upset_top_annotation(m2,
                                                                          axis_param = list(gp = gpar(fontsize = 5)),
                                                                          annotation_name_gp = gpar(fontsize = 5),
                                                                          height = unit(1, "cm")),
                                    right_annotation = upset_right_annotation(m2,
                                                                              axis_param = list(gp = gpar(fontsize = 5)),
                                                                              annotation_name_gp = gpar(fontsize = 5),
                                                                              width = unit(1, "cm")),
                                    row_names_gp = gpar(fontsize = 5))))

plot_grid(
  p_upset,
  align = 'h',
  hjust = 0.1, # -1
  nrow = 1
)

ggsave2(paste0(plots_dir, "upset_shared_finemap_eqtl.png"),
        width = onecol_width*1.2, height = onecol_height, units = 'cm', dpi = 300)
