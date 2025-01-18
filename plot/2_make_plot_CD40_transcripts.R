
png(filename=paste0(plots_dir, gene, ".pchic.trx.png"), 
    width = 8.5, 
    height = 11,
    units="in",
    res=300)

pageCreate(width = 8.5, 
           height = 11, # 6, 11
           default.units = "inches")

hl_width <- 5000
p_gwas <- plotManhattan(
  data = df_gwas,
  sigVal = 0,
  chrom=chr, chromstart = chr_start, chromend = chr_end,
  fill = colorby("sig",
                 palette = colorRampPalette(c("grey50"))),
  leadSNP = list(snp = leadsnp, pch = 18, cex = 0.75, fill = "red",
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
  label = paste0("GWAS"), x = 1, y = 1, 
  fontsize = 8, fontface = "bold", just = "center",
  default.units = "inches"
)

annoYaxis(
  plot = p_gwas,
  at = c(0, 3),
  #x = 0, y = 0,
  axisLine = TRUE, fontsize = 6
)

p_plasma <- plotManhattan(
  data = df_eqtl %>% filter(celltype == "Plasma"),
  sigVal = 0,
  chrom=chr, chromstart = chr_start, chromend = chr_end,
  fill = colorby("sig",
                 palette = colorRampPalette(c(col_plasma))),
  leadSNP = list(snp = eqtl_snp, pch = 18, cex = 0.75, fill = "red",
                 fontsize = 0),
  assembly = assembly_name,
  x=0.5, y=1.4,
  width=8, height=1
)

annoHighlight(
  plot = p_plasma,
  chrom = chr,
  chromstart = lead_pos-hl_width, chromend = lead_pos+hl_width,
  x=0.5, y=1.4, height = 1, just = c("left", "top"),
  default.units = "inches"
)


annoYaxis(
  plot = p_plasma,
  at = c(0, 3),
  axisLine = TRUE, fontsize = 6
)

plotText(
  label = paste0("Plasma"), x = 1, y = 2, 
  fontsize = 8, fontface = "bold", just = "center",
  default.units = "inches"
)

p_B_in <- plotManhattan(
  data = df_eqtl %>% filter(celltype == "B_IN"),
  sigVal = 0,
  pch = 1,
  chrom=chr, chromstart = chr_start, chromend = chr_end,
  fill = colorby("sig",
                 palette = colorRampPalette(c(col_cd4nc))),
  leadSNP = list(snp = eqtl_snp, pch = 18, cex = 0.75, fill = "red",
                 fontsize = 0),
  assembly = assembly_name,
  x=0.5, y=2.6,
  width=8, height=0.3
)

annoYaxis(
  plot = p_B_in,
  at = c(0, 3),
  axisLine = TRUE, fontsize = 6
)

plotText(
  label = paste0("B_IN"), x = 1, y = 2.6, 
  fontsize = 8, fontface = "bold", just = "center",
  default.units = "inches"
)

p_B_Mem <- plotManhattan(
  data = df_eqtl %>% filter(celltype == "B_Mem"),
  sigVal = 0,
  pch = 1,
  chrom=chr, chromstart = chr_start, chromend = chr_end,
  fill = colorby("sig",
                 palette = colorRampPalette(c(col_cd8nc))),
  leadSNP = list(snp = eqtl_snp, pch = 18, cex = 0.75, fill = "red",
                 fontsize = 0),
  assembly = assembly_name,
  x=0.5, y=3.1,
  width=8, height=0.3
)

annoYaxis(
  plot = p_B_Mem,
  at = c(0, 3),
  axisLine = TRUE, fontsize = 6
)

plotText(
  label = paste0("B_Mem"), x = 1, y = 3.1, 
  fontsize = 8, fontface = "bold", just = "center",
  default.units = "inches"
)

p_monoc <- plotManhattan(
  data = df_eqtl %>% filter(celltype == "Mono_C"),
  sigVal = 0,
  pch = 1,
  chrom=chr, chromstart = chr_start, chromend = chr_end,
  fill = colorby("sig",
                 palette = colorRampPalette(c(col_monoc))),
  leadSNP = list(snp = eqtl_snp, pch = 18, cex = 0.75, fill = "red",
                 fontsize = 0),
  assembly = assembly_name,
  x=0.5, y=3.6,
  width=8, height=0.5
)

annoYaxis(
  plot = p_monoc,
  at = c(0, 3),
  axisLine = TRUE, fontsize = 6
)

plotText(
  label = paste0("Mono_C"), x = 1, y = 3.8, 
  fontsize = 8, fontface = "bold", just = "center",
  default.units = "inches"
)

p_mononc <- plotManhattan(
  data = df_eqtl %>% filter(celltype == "Mono_NC"),
  sigVal = 0,
  pch = 1,
  chrom=chr, chromstart = chr_start, chromend = chr_end,
  fill = colorby("sig",
                 palette = colorRampPalette(c(col_mononc))),
  leadSNP = list(snp = eqtl_snp, pch = 18, cex = 0.75, fill = "red",
                 fontsize = 0),
  assembly = assembly_name,
  x=0.5, y=4.3,
  width=8, height=0.5
)

annoYaxis(
  plot = p_mononc,
  at = c(0, 3),
  axisLine = TRUE, fontsize = 6
)

plotText(
  label = paste0("Mono_NC"), x = 1, y = 4.5, 
  fontsize = 8, fontface = "bold", just = "center",
  default.units = "inches"
)

p_bulk <- plotManhattan(
  data = df_eqtl %>% filter(celltype == "allcells"),
  sigVal = 0,
  pch = 1,
  chrom=chr, chromstart = chr_start, chromend = chr_end,
  fill = colorby("sig",
                 palette = colorRampPalette(c("black"))),
  leadSNP = list(snp = eqtl_snp, pch = 18, cex = 0.75, fill = "red",
                 fontsize = 0),
  assembly = assembly_name,
  x=0.5, y=5.2,
  width=8, height=0.2
)

plotText(
  label = paste0("Bulk"), x = 1, y = 5.2, 
  fontsize = 8, fontface = "bold", just = "center",
  default.units = "inches"
)

annoYaxis(
  plot = p_bulk,
  at = c(0, 3),
  axisLine = TRUE, fontsize = 6
)

p_lm <- plotManhattan(
  data = lm_eqtl %>% filter(celltype == "Plasma"),
  sigVal = 0,
  pch = 1,
  chrom=chr, chromstart = chr_start, chromend = chr_end,
  fill = colorby("sig",
                 palette = colorRampPalette(c("black"))),
  leadSNP = list(snp = eqtl_snp, pch = 18, cex = 0.75, fill = "red",
                 fontsize = 0),
  assembly = assembly_name,
  x=0.5, y=5.5,
  width=8, height=0.5
)

plotText(
  label = paste0("Plasma (Linear)"), x = 1, y = 5.6, 
  fontsize = 8, fontface = "bold", just = "center",
  default.units = "inches"
)

annoYaxis(
  plot = p_lm,
  at = c(0, 3),
  axisLine = TRUE, fontsize = 6
)
# coac curoff 0.25
55443409-55444683
55290821-55273640
# tss: 55290821, 55290821
# 55273640, 55278757

# search pchic
HiC <- HiC %>% 
  mutate(start_mid = (start1 + end1)/2, end_mid = (start2 + end2)/2,
         dist_1 = abs(lead_pos-start_mid), 
         dist_2 = abs(lead_pos-end_mid)) %>% 
  filter(chrom1 == chr & chrom2 == chr) %>% 
  arrange(dist_2) %>% 
  slice_head(n=1)

archPlot <- plotPairsArches(
  data = HiC, params = params_c,
  x = 0.5, y = 6.2, height = 0.5, width = 8, alpha=1,
  fill = col_plasma,
  just = c("left", "top"), default.units = "inches",
  flip = FALSE,
  archHeight = "score"
)
plotText(
  label = eqtl_snp, x = 4.7, y = 6.1, fontcolor ="firebrick2",
  fontsize = 8, fontface = "bold", just = "center",
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
  x=0.5, y=6.2, height = 0.5, just = c("left", "top"),
  default.units = "inches"
)

plotText(
  label = paste0("PCHi-C in B cells"),
  x = 5.7, 
  y = 6.4, #rot = 90,
  fontsize = 8, fontface = "bold", just = "center",
  default.units = "inches"
)


# # gene annotation
# plotGenes(
#   assembly = assembly_name,
#   chrom = chr, chromstart = chr_start, chromend = chr_end,
#   x=0.5, y = 6.7, width=8, height=0.5,
#   geneHighlights = data.frame(
#     "gene" = gene,
#     "color" = "firebrick2"
#   )
# )

p_chipseq1 <- plotSignal(data = df_chipseq1,
                         params = c(params_c, ctcf_range1),
                         fill = col_plasma, linecolor = col_plasma,
                         x = 0.5, y = 6.7,
                         height = 0.5, width = 8,
                         default.units = "inches")

plotText(
  label = paste0("B cell"), x = 1.2, y = 6.9, #rot = 90,
  fontsize = 8, fontface = "bold", just = "center",
  default.units = "inches"
)

annoHighlight(
  plot = p_chipseq1,
  chrom = chr,
  chromstart = lead_pos-hl_width, chromend = lead_pos+hl_width,
  x=1.2, y=6.7, height = 0.5, just = c("left", "top"),
  default.units = "inches"
)

## Plot genome label
plotGenomeLabel(
  chrom = chr, chromstart = chr_start, chromend = chr_end,
  assembly = "hg19",
  x = 0.5, y = 7.3, length = 8, default.units = "inches",
  scale = "Mb"
)


## add segment lines
plotSegments(
  x0 = 4.55, y0 = 7.5, x1 = 4.55, y1 = 7.8,
  default.units = "inches",
  lwd = 1, lty = 1
)

plotSegments(
  x0 = 4.55, y0 = 7.8, x1 = 1, y1 = 8,
  default.units = "inches",
  lwd = 1, lty = 1
)

plotSegments(
  x0 = 4.55, y0 = 7.8, x1 = 8, y1 = 8,
  default.units = "inches",
  lwd = 1, lty = 1
)

p_trx <- wiggleplotr::plotTranscripts(cd40_exons, #cd40_cdss, cd40_metadata, 
                                      rescale_introns = FALSE,
                                      region_coords = c(lead_pos-3000, lead_pos+3000))
layer_scales(p_trx)$y$range$range
layer_scales(p_trx)$x$range$range
plotGG(
  plot = p_trx,
    #annotate("text", x = lead_pos, y = 1, colour = "red", label="test"),
    #theme_void(), 
  x = 0.5, y = 8, width = 8, height = 2,
  just = c("left", "top"), default.units = "inches"
)

# add two SNPs
plotCircle(
  x = 1.75, y = 8.45, r = 0.05,
  default.units = "inches", fill = "orange"
)

# plotCircle(
#   x = 4.2, y = 8.45, r = 0.05, fill = "red",
#   default.units = "inches"
# )

plotRect(
  x = 4.2, y = 8.4, width = 0.1, height = 0.1,
  just = c("left", "top"), default.units = "inches",
  lwd = 1, fill = "red",
)

# read TFBS motif
tfbs_motif <- png::readPNG("../result/annotation/figures/ZNF263_TFBS.png")

plotRaster(
  image = tfbs_motif,
  x = 1.5, y = 7.5, width = 1.5, height = 0.5,
  just = c("left", "top")
)
# plotTranscripts(
#   chrom = chr, chromstart = tss_start-10000, chromend = tss_start+10000,
#   #assembly = "hg19", 
#   labels = "both",
#   assembly = ,
#   x = 0.5, y = 8, width = 8, height = 2.5,
#   just = c("left", "top"), default.units = "inches"
# )

plotText(
  label = paste0("H3K27ac"), 
  x = 0.3, y = 7, rot = 90,
  fontsize = 8, fontface = "bold", just = "center",
  default.units = "inches"
)

plotText(
  label = paste0("-log10(Pval)"), 
  x = 0.2, y = 4, # 8.7
  rot = 90,
  fontsize = 8, fontface = "bold", just = "center",
  default.units = "inches"
)

# remove grid lines
pageGuideHide()


dev.off()

