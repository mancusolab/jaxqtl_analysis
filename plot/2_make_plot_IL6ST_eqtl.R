
page_width <- 8.5 # 8.5
page_height <- 11
png(filename=paste0(plots_dir, gene, ".pchic.blankannot.png"), 
    width = page_width,
    height = page_height,
    units="in",
    res=300)

hl_width <- 5000

pageCreate(width = page_width,
           height = page_height,
           default.units = "inches")

p_gwas <- plotManhattan(
  data = df_gwas,
  sigVal = 0,
  chrom=chr, chromstart = chr_start, chromend = chr_end,
  fill = colorby("sig",
                 palette = colorRampPalette(c("grey50"))),
  leadSNP = list(snp = leadsnp, pch = 18, cex = 0.7, fill = "red",
                 fontsize = 0),
  assembly = assembly_name, 
  x=0.5, y=0.3, 
  width=8, 
  height=1.6,
  default.units = "inches"
)

annoHighlight(
  plot = p_gwas,
  chrom = chr,
  chromstart = lead_pos-hl_width, chromend = lead_pos+hl_width,
  x=0.5, y=0.3, height = 1.6, just = c("left", "top"),
  default.units = "inches"
)

# plotText(
#   label = paste0("GWAS"), 
#   x = 1, # 1
#   y = 1.6, # 1
#   fontsize = 10, fontface = "bold", just = "center",
#   default.units = "inches"
# )

annoYaxis(
  plot = p_gwas,
  at = c(0, 12),
  #x = 0, y = 0,
  axisLine = TRUE, fontsize = 6
)

p_cd4nc <- plotManhattan(
  data = df_eqtl %>% filter(celltype == "CD4_NC"),
  sigVal = 0,
  chrom=chr, chromstart = chr_start, chromend = chr_end,
  fill = colorby("sig",
                 palette = colorRampPalette(c(col1))),
  leadSNP = list(snp = leadsnp, pch = 18, cex = 0.7, fill = "red",
                 fontsize = 0),
  assembly = assembly_name,
  x=0.5, y=2,
  width=8, height=0.57
)

annoHighlight(
  plot = p_cd4nc,
  chrom = chr,
  chromstart = lead_pos-hl_width, chromend = lead_pos+hl_width,
  x=0.5, y=2,
  height = 0.57, just = c("left", "top"),
  default.units = "inches"
)

annoYaxis(
  plot = p_cd4nc,
  at = c(0, 12),
  axisLine = TRUE, fontsize = 6
)

# plotText(
#   label = "CD4_NC",  # expression("CD4"["NC"])
#   x = 1, # 1
#   y = 2.3, # 2
#   fontsize = 10, fontface = "bold", just = "center",
#   default.units = "inches"
# )

p_cd8nc <- plotManhattan(
  data = df_eqtl %>% filter(celltype == "CD8_NC"),
  sigVal = 0,
  chrom=chr, chromstart = chr_start, chromend = chr_end,
  fill = colorby("sig",
                 palette = colorRampPalette(c(col2))),
  leadSNP = list(snp = leadsnp, pch = 18, cex = 0.7, fill = "red",
                 fontsize = 0),
  assembly = assembly_name,
  x=0.5, y=2.9,
  width=8, height=0.38
)

annoHighlight(
  plot = p_cd8nc,
  chrom = chr,
  chromstart = lead_pos-hl_width, chromend = lead_pos+hl_width,
  x=0.5, y=2.9,
  height = 0.4, just = c("left", "top"),
  default.units = "inches"
)

annoYaxis(
  plot = p_cd8nc,
  at = c(0, 12),
  axisLine = TRUE, fontsize = 6
)

# plotText(
#   label = paste0("CD8_NC"), 
#   x = 1, # 1
#   y = 3, # 3
#   fontsize = 10, fontface = "bold", just = "center",
#   default.units = "inches"
# )

p_Bmem <- plotManhattan(
  data = df_eqtl %>% filter(celltype == "B_Mem"),
  sigVal = 0,
  pch = 1,
  chrom=chr, chromstart = chr_start, chromend = chr_end,
  fill = colorby("sig",
                 palette = colorRampPalette(c(col3))),
  leadSNP = list(snp = leadsnp, pch = 18, cex = 0.7, fill = "red",
                 fontsize = 0),
  assembly = assembly_name,
  x=0.5, y=3.8,
  width=8, height=0.18
)

annoYaxis(
  plot = p_Bmem,
  at = c(0, 12),
  axisLine = TRUE, fontsize = 6
)

# plotText(
#   label = paste0("B_Mem"), 
#   x = 1, # 1
#   y = 3.7, # 3.7
#   fontsize = 10, fontface = "bold", just = "center",
#   default.units = "inches"
# )

p_nk <- plotManhattan(
  data = df_eqtl %>% filter(celltype == "NK"),
  sigVal = 0,
  pch = 1,
  chrom=chr, chromstart = chr_start, chromend = chr_end,
  fill = colorby("sig",
                 palette = colorRampPalette(c(col4))),
  leadSNP = list(snp = leadsnp, pch = 18, cex = 0.7, fill = "red",
                 fontsize = 0),
  assembly = assembly_name,
  x=0.5, y=4.5,
  width=8, height=0.13
)

annoYaxis(
  plot = p_nk,
  at = c(0, 12),
  axisLine = TRUE, fontsize = 6
)

# plotText(
#   label = paste0("NK"), 
#   x = 1,  # 1
#   y = 4.4, # 4.4
#   fontsize = 10, fontface = "bold", just = "center",
#   default.units = "inches"
# )

p_monoc <- plotManhattan(
  data = df_eqtl %>% filter(celltype == "Mono_C"),
  sigVal = 0,
  pch = 1,
  chrom=chr, chromstart = chr_start, chromend = chr_end,
  fill = colorby("sig",
                 palette = colorRampPalette(c(col5))),
  leadSNP = list(snp = leadsnp, pch = 18, cex = 0.7, fill = "red",
                 fontsize = 0),
  assembly = assembly_name,
  x=0.5, y=5.2,
  width=8, height=0.14
)

annoYaxis(
  plot = p_monoc,
  at = c(0, 12),
  axisLine = TRUE, fontsize = 6
)

# plotText(
#   label = paste0("Mono_C"), 
#   x = 1, # 1
#   y = 5, # 5.0
#   fontsize = 10, fontface = "bold", just = "center",
#   default.units = "inches"
# )

p_bulk <- plotManhattan(
  data = df_eqtl %>% filter(celltype == "allcells"),
  sigVal = 0,
  pch = 1,
  chrom=chr, chromstart = chr_start, chromend = chr_end,
  fill = colorby("sig",
                 palette = colorRampPalette(c("black"))),
  leadSNP = list(snp = leadsnp, pch = 18, cex = 0.7, fill = "red",
                 fontsize = 0),
  assembly = assembly_name,
  x=0.5, y=5.8,
  width=8, height=0.24
)

# plotText(
#   label = paste0("Bulk"), 
#   x = 1, # 1
#   y = 5.7, # 5.5
#   fontsize = 10, fontface = "bold", just = "center",
#   default.units = "inches"
# )

annoYaxis(
  plot = p_bulk,
  at = c(0, 12),
  axisLine = TRUE, fontsize = 6
)

# # plot Arc
# HiC$length <-
#   (HiC$start2 - HiC$start1) / 1000
# 
# ## Translate lengths into heights
# HiC$h <-
#   HiC$length / max(HiC$length)

# 5_55444683

# coac curoff 0.25
55443409-55444683
55290821-55273640
# tss: 55290821, 55290821
# 55273640, 55278757

# search pchic
HiC <- HiC %>% 
  filter(chrom1 == chr & start2 <= lead_pos & end2 >= lead_pos & chrom2 == chr) %>%
  filter(start1 <= tss_start & end1 >= tss_start) %>% 
  #filter(start1 == 55293505) %>% 
  filter(celltype %in% c("nCD4", "nCD8") & score >= 5)

# HiC %>% filter(chrom1 == chr & start2 <= 55444683 & end2 >= 55444683 & chrom2 == chr) %>%
#   #filter(start1 == 55293505) %>% 
#   filter(celltype %in% c("nCD4", "nCD8"))

#hic_col <- c("#F1932D", "#F6C141", "#5289C7", "#882E72")
hic_col <- c("#5289C7", "#882E72")
my_pal <- colorRampPalette(hic_col)

archPlot <- plotPairsArches(
  data = HiC, params = params_c,
  x = 0.5, y = 6.3, height = 0.5, width = 8, alpha=1,
  fill = colorby("score", palette = my_pal),
  just = c("left", "top"), default.units = "inches",
  flip = FALSE,
  archHeight = "score"
)

plotText(
  label = leadsnp, x = 5.7, y = 6.15, fontcolor ="firebrick2",
  fontsize = 10, fontface = "bold", just = "center",
  default.units = "inches"
)

annoYaxis(
  plot = archPlot,
  at = c(0, round(HiC$score, digit=2)),
  axisLine = TRUE, fontsize = 6
)

annoHighlight(
  plot = archPlot,
  chrom = chr,
  chromstart = lead_pos-hl_width, chromend = lead_pos+hl_width,
  x=0.5, y=6.3, height = 0.5, just = c("left", "top"),
  default.units = "inches"
)

# plotText(
#   label = paste0("PCHi-C in T cells"),
#   x = 6.7, 
#   y = 6.5, #rot = 90,
#   fontsize = 10, fontface = "bold", just = "center",
#   default.units = "inches"
# )


# archPlot <- as.ggplot(grid.grabExpr(plot_connections(coac %>%
#                                                        filter(Peak2 %in% c("chr5_55443409_55444796") &
#                                                                 Peak1 %in% c("chr5_55273640_55274515")),
#                                                      chr, chr_start, chr_end,
#                                                      viewpoint = "chr5_55443409_55444796",
#                                                      # gene_model = gene_anno %>%
#                                                      #   filter(symbol %in% c("IL6ST", "ANKRD55",
#                                                      #                        "SLC38A9", "DDX4", "IL31RA", "IL6ST-DT", "RNF138P1")),
#                                                      coaccess_cutoff = 0.3,
#                                                      connection_width = .5,
#                                                      include_axis_track = FALSE,
#                                                      collapseTranscripts = FALSE)))
# plotGG(
#   plot = archPlot,
#   x = 0.5, y = 6,
#   width = 8, height = 1, just = c("left", "top")
# )

# gene annotation
plotGenes(
  assembly = assembly_name,
  chrom = chr, chromstart = chr_start, chromend = chr_end,
  x=0.5, y = 6.8, width=8, height=0.5, fontsize=8.5,
  geneHighlights = data.frame(
    "gene" = gene,
    "color" = "firebrick2"
  )
)

p_chipseq1 <- plotSignal(data = df_chipseq1,
                 params = c(params_c, ctcf_range1),
                 fill = col1, linecolor = col1,
                 x = 0.5, y = 7.35,
                 height = 0.5, width = 8,
                 default.units = "inches")

annoHighlight(
  plot = p_chipseq1,
  chrom = chr,
  chromstart = lead_pos-hl_width, chromend = lead_pos+hl_width,
  x=0.5, y=7.35, height = 0.5, just = c("left", "top"),
  default.units = "inches"
)

# plotText(
#   label = paste0("CD4+ T"), 
#   x = 1.2, # 1.2
#   y = 7.5, # 7.4
#   fontsize = 10, fontface = "bold", just = "center",
#   default.units = "inches"
# )

p_chipseq2 <- plotSignal(data = df_chipseq2, 
                         params = c(params_c, ctcf_range2),
                 fill = col2, linecolor = col2,
                 x = 0.5, y = 8,
                 height = 0.5, width = 8,
                 default.units = "inches")

annoHighlight(
  plot = p_chipseq2,
  chrom = chr,
  chromstart = lead_pos-hl_width, chromend = lead_pos+hl_width,
  x=0.5, y=8, height = 0.5, just = c("left", "top"),
  default.units = "inches"
)

# plotText(
#   label = paste0("CD8+ T"),
#   x = 1.2, # 1.2
#   y = 8.2, # 8.1
#   fontsize = 10, fontface = "bold", just = "center",
#   default.units = "inches"
# )

p_chipseq3 <- plotSignal(data = df_chipseq3, 
                         params = c(params_c, ctcf_range3),
                 fill = col3, linecolor = col3,
                 x=0.5, y=8.6, 
                 height=0.5, width=8,
                 default.units = "inches")

# plotText(
#   label = paste0("B cell"), 
#   x = 1.2, # 1.2
#   y = 8.8, # 8.7
#   fontsize = 10, fontface = "bold", just = "center",
#   default.units = "inches"
# )

p_chipseq4 <- plotSignal(data = df_chipseq4, 
                         params = c(params_c, ctcf_range4),
                         fill = col4, linecolor = col4,
                         x =0.5, y=9.3,
                         height=0.5, width=8,
                         default.units = "inches")

# plotText(
#   label = paste0("NK"), 
#   x = 1.2, # 1.2
#   y = 9.45, # 9.3
#   fontsize = 10, fontface = "bold", just = "center",
#   default.units = "inches"
# )

p_chipseq5 <- plotSignal(data = df_chipseq5, 
                         params = c(params_c, ctcf_range5),
                         fill = col5, linecolor = col5,
                         x =0.5, y=9.9, height=0.5, width=8,
                         default.units = "inches")

# plotText(
#   label = paste0("Monocyte"), 
#   x = 1.2, # 1.2
#   y = 10.1, # 10
#   fontsize = 10, fontface = "bold", just = "center",
#   default.units = "inches"
# )

plotText(
  label = paste0("H3K27ac"), 
  x = 0.2, y = 8.7, # 8.7
  rot = 90,
  fontsize = 10, fontface = "bold", just = "center",
  default.units = "inches"
)

plotText(
  label = paste0("-log10(Pval)"), 
  x = 0.2, y = 4, # 8.7
  rot = 90,
  fontsize = 10, fontface = "bold", just = "center",
  default.units = "inches"
)

plotGenomeLabel(
  chrom = chr, chromstart = chr_start, chromend = chr_end,
  assembly = "hg19",
  x = 0.5, y = 10.5, length = 8, default.units = "inches", fontsize = 10,
  scale = "Mb"
)

# remove grid lines
pageGuideHide()


dev.off()

