### plot funcs
library(tidyverse)
library(VennDiagram)

axis_theme <- theme_classic()+
  theme(axis.title.y = element_text(size=8, color = "black"),
        axis.title.x = element_text(size=8, color = "black"),
        axis.text.x = element_text(size=8, color = "black"),
        # title = element_text(size=8, color = "black"),
        plot.title = element_text(size=8, color = "black"),
        axis.text.y = element_text(size=8, color = "black"))

axis_theme_font6 <- theme_classic()+
  theme(axis.title.y = element_text(size=6, color = "black"),
        axis.title.x = element_text(size=6, color = "black"),
        axis.text.x = element_text(size=6, color = "black"),
        # title = element_text(size=8, color = "black"),
        plot.title = element_text(size=6, color = "black"),
        axis.text.y = element_text(size=6, color = "black"))

# AJHG style
full_width <- 17.4 # cm
onecol_width <- 8.5 # cm
full_height <- 7.6 # cm
onecol_height <- 6.5 # cm

# method_cols <- c("tqtl" = "#56B4E9", 
#                  "jaxqtl" = "#D55E00",
#                  "both" = "grey",
#                  "pois_score" = "#CC79A7")


# color palette for continuous scale
pal_red <- colorRampPalette(c("white", "red3"), space = "rgb")(200)
pal_blue <- colorRampPalette(c("white", "blue"), space = "rgb")(200)
pal_br <- colorRampPalette(c("blue", "white", "red3"), space = "rgb")(200)
pal_bry <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                              "#E0F3F8","#91BFDB","#4575B4")))(64)

pal_bry <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                                  "#E0F3F8","#91BFDB","#4575B4")))(64)

method_cols <- c("tensorqtl" = "#0072B2",
                 "tensorQTL" = "#0072B2",
                 "linear" = "#56B4E9", "Linear" = "#56B4E9",
                 "jaxQTL-Linear" = "#56B4E9",
                 "tensorQTL (Zhou et.al)"="#0072B2",
                 "NegBinom" = "#D55E00",
                 "negbinom" = "#D55E00",
                 "SAIGE-QTL" = "#F0E442",
                 "jaxQTL" = "#D55E00",
                 "jaxQTL-NegBinom" =  "#D55E00",
                 "jaxQTL-NegBinom-GPU" =  "#D55E00",#"tomato",
                 "jaxQTL-NegBinom-TPU" =  "#D55E00",#"tomato",
                 "jaxQTL-NegBinom-acat" = "#D55E00",
                 "both" = "grey",
                 "Poisson" = "#CC79A7",
                 "none" = "white")

enhancer_promoter_cols <- c("enhancer"="#009E73",
                            "promoter"="#E69F00")

sc_bulk_cols <- c("bulk-eQTL" = "black",
                  "sc-eQTL" = "#D55E00")

celltype_cols <- c("CD4_NC" = "#882E72",
                   "CD4_ET" = "#B178A6",
                   "CD4_SOX4" = "#D6C1DE",
                   "CD8_ET" = "#1965B0",
                   "CD8_NC" = "#5289C7",
                   "CD8_S100B" = "#7BAFDE",
                   "NK" = "#4EB265",
                   "NK_R" = "#90C987",
                   "Plasma" = "#CAE0AB",
                   "B_Mem" = "#F7EE55",
                   "B_IN" = "#F6C141",
                   "Mono_C" = "#F1932D",
                   "Mono_NC" = "#E8601C",
                   "DC" = "#DC050C",
                   "allcells" = "black",
                   "bulk-eQTL" = "black",
                   "non-eGenes" = "grey")

celltype_bquote <- c("CD4_NC" = bquote(CD4[NC]),
                   "CD4_ET" = bquote(CD4[ET]),
                   "CD4_SOX4" = bquote(CD4[SOX4]),
                   "CD8_ET" = bquote(CD8[ET]),
                   "CD8_NC" = bquote(CD8[NC]),
                   "CD8_S100B" = bquote(CD8[S100B]),
                   "NK" = bquote(NK),
                   "NK_R" = bquote(NK[R]),
                   "Plasma" = bquote(Plasma),
                   "B_Mem" = bquote(B[Mem]),
                   "B_IN" = bquote(B[IN]),
                   "Mono_C" = bquote(Mono[C]),
                   "Mono_NC" = bquote(Mono[NC]),
                   "DC" = bquote(DC))

sim_method_cols <- c("lm_gtex" = "#009E73",
                     "lm" = "#56B4E9", 
                     "lm_wald" = "#56B4E9",
                     "NegBinom_score_thr" = "#E69F00",
                     "pois_score" = "#CC79A7",
                     "NegBinom_score" = "#D55E00")

pal <- colorRampPalette(c("blue", "white", "red3"), space = "rgb")(200)
pal_blue_red <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                              "#E0F3F8","#91BFDB","#4575B4")))(64)


# meta analysis plot

forest_plot_theme <- function(p){
  p + theme_classic() + theme(panel.border = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.text.x.bottom = element_text(size = 8, colour = "black"),
            axis.title.x = element_text(size = 8, colour = "black"),
            legend.text=element_text(size=6),
            legend.title = element_text(size=6),
            legend.position = "right",
            title =element_text(size=8))
}

plot_ggvenn <- function(set1, set2, cell_type, gtex_1per="1per", model="nb"){
  # set1: jaxqtl cis genes
  # set2: tqtl cis genes
  res <- list()
  combined <- union(set1, set2)
  idx <- 1:length(combined)
  res['jaxqtl'] <- list(idx[combined %in% set1])
  res['tqtl'] <- list(idx[combined %in% set2])
  
  venn.diagram(res, paste0(plots_dir, "Venn/cisgene_venn_tqtl_", gtex_1per, "_", "jaxqtl_", model,"_", cell_type, ".png"), 
               fill = c("#D55E00", "#56B4E9"),
               imagetype="png", alpha=0.8, sub.fontfamily="Arial", cex=2, 
               main = cell_type, main.fontfamily="Arial",main.cex=2,
               cat.fontfamily="Arial",
               disable.logging = TRUE)
}

