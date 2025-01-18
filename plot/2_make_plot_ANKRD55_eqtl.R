
png(filename=paste0(plots_dir, gene, ".png"), 
    width = 8.5, 
    height = 9,
    units="in",
    res=300)

pageCreate(width = 8.5, 
           height = 9, # 6, 11
           default.units = "inches")

hl_width <- 5000

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
  height=1,
  default.units = "inches"
)

annoHighlight(
  plot = p_gwas,
  chrom = chr,
  chromstart = lead_pos-hl_width, chromend = lead_pos+hl_width,
  x=0.5, y=0.3, height = 1, just = c("left", "top"),
  default.units = "inches"
)

plotText(
  label = paste0("GWAS"), 
  x = 1, # 1
  y = 1, # 1
  fontsize = 10, fontface = "bold", just = "center",
  default.units = "inches"
)

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
  x=0.5, y=1.5,
  width=8, height=2
)

annoHighlight(
  plot = p_cd4nc,
  chrom = chr,
  chromstart = lead_pos-hl_width, chromend = lead_pos+hl_width,
  x=0.5, y=1.5,
  height = 2, just = c("left", "top"),
  default.units = "inches"
)

annoYaxis(
  plot = p_cd4nc,
  at = c(0, 12),
  axisLine = TRUE, fontsize = 6
)

plotText(
  label = paste0("CD4_NC"), 
  x = 1, # 1
  y = 3.2, # 2
  fontsize = 10, fontface = "bold", just = "center",
  default.units = "inches"
)

p_cd8nc <- plotManhattan(
  data = df_eqtl %>% filter(celltype == "CD8_NC"),
  sigVal = 0,
  chrom=chr, chromstart = chr_start, chromend = chr_end,
  fill = colorby("sig",
                 palette = colorRampPalette(c(col2))),
  leadSNP = list(snp = leadsnp, pch = 18, cex = 0.7, fill = "red",
                 fontsize = 0),
  assembly = assembly_name,
  x=0.5, y=3.7,
  width=8, height=0.3
)

annoHighlight(
  plot = p_cd8nc,
  chrom = chr,
  chromstart = lead_pos-hl_width, chromend = lead_pos+hl_width,
  x=0.5, y=3.7,
  height = 0.3, just = c("left", "top"),
  default.units = "inches"
)

annoYaxis(
  plot = p_cd8nc,
  at = c(0, 12),
  axisLine = TRUE, fontsize = 6
)

plotText(
  label = paste0("CD8_NC"), 
  x = 1, # 1
  y = 3.7, # 3
  fontsize = 10, fontface = "bold", just = "center",
  default.units = "inches"
)

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
  x=0.5, y=4.4,
  width=8, height=0.2
)

annoYaxis(
  plot = p_Bmem,
  at = c(0, 12),
  axisLine = TRUE, fontsize = 6
)

plotText(
  label = paste0("B_Mem"),
  x = 1, # 1
  y = 4.3, # 3.7
  fontsize = 10, fontface = "bold", just = "center",
  default.units = "inches"
)

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
  x=0.5, y=5.1,
  width=8, height=0.1
)

annoYaxis(
  plot = p_nk,
  at = c(0, 12),
  axisLine = TRUE, fontsize = 6
)

plotText(
  label = paste0("NK"),
  x = 1,  # 1
  y = 5, # 4.4
  fontsize = 10, fontface = "bold", just = "center",
  default.units = "inches"
)

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
  x=0.5, y=5.7,
  width=8, height=0.14
)

annoYaxis(
  plot = p_monoc,
  at = c(0, 12),
  axisLine = TRUE, fontsize = 6
)

plotText(
  label = paste0("Mono_C"),
  x = 1, # 1
  y = 5.6, # 5.0
  fontsize = 10, fontface = "bold", just = "center",
  default.units = "inches"
)

p_bulk <- plotManhattan(
  data = df_eqtl %>% filter(celltype == "allcells"),
  sigVal = 0,
  chrom=chr, chromstart = chr_start, chromend = chr_end,
  fill = colorby("sig",
                 palette = colorRampPalette(c("black"))),
  leadSNP = list(snp = leadsnp, pch = 18, cex = 0.7, fill = "red",
                 fontsize = 0),
  assembly = assembly_name,
  x=0.5, y=6.2,
  width=8, height=1.5
)


annoHighlight(
  plot = p_bulk,
  chrom = chr,
  chromstart = lead_pos-hl_width, chromend = lead_pos+hl_width,
  x=0.5, y=6.2,
  height = 1.5, just = c("left", "top"),
  default.units = "inches"
)

plotText(
  label = paste0("Bulk"),
  x = 1, # 1
  y = 7.3, # 5.5
  fontsize = 10, fontface = "bold", just = "center",
  default.units = "inches"
)

annoYaxis(
  plot = p_bulk,
  at = c(0, 12),
  axisLine = TRUE, fontsize = 6
)

plotText(
  label = leadsnp, x = 3.9, y = 7.8, fontcolor ="firebrick2",
  fontsize = 10, fontface = "bold", just = "center",
  default.units = "inches"
)

# gene annotation
plotGenes(
  assembly = assembly_name,
  chrom = chr, chromstart = chr_start, chromend = chr_end,
  x=0.5, y = 8, width=8, height=0.5, fontsize=8.5,
  geneHighlights = data.frame(
    "gene" = gene,
    "color" = "firebrick2"
  )
)

plotText(
  label = paste0("-log10(Pval)"), 
  x = 0.1, y = 4, # 8.7
  rot = 90,
  fontsize = 10, fontface = "bold", just = "center",
  default.units = "inches"
)

plotGenomeLabel(
  chrom = chr, chromstart = chr_start, chromend = chr_end,
  assembly = "hg19",
  x = 0.5, y = 8.6, length = 8, default.units = "inches", fontsize = 10,
  scale = "Mb"
)

# remove grid lines
pageGuideHide()


dev.off()
