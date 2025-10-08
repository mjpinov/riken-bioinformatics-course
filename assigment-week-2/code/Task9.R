################# define functions 
RESULTS_FOLDER <- "/Users/marjos/Documents/courses/ki/RIKEN-bioinformatics/riken-bioinformatics-course/assigment-week-2/results/Task9/redo"

load.packages <- function() {
  package.list = c(
    'ggplot2','dplyr','ggpubr', 'ggrepel','tidyr','ggplotify', 'scatterpie', 'grid', 'data.table'
  )
  suppressMessages(suppressWarnings(sapply(package.list, library, character.only = TRUE, quietly = TRUE)))
  return()
}
load.packages()
set.param <- function() {
  size.list <- list(
    task09.atlas_UMAP=list(width=7, height=9),
    task09.endo_lv2_UMAP =list(width=5, height=4.75),
    task09.endo_tissue_UMAP =list(width=5, height=4.5),
    task09.all =list(width=7, height=9)
  )
  font.size.list <- list(
    axis.text = 10,
    legend.text = 10,
    axis.title = 10,
    legend.title = 10,
    plot.title = 10
  )
  point.size.list <- list(
    encircle = 5,
    metacell = 1,
    white = 3
  )
  alpha.list <- list(
    encircle = 1
  )
  param.list <- list(
    point.size.list = point.size.list,
    alpha.list = alpha.list,
    size.list = size.list,
    font.size.list = font.size.list
  )
  return (param.list)
}


################################################################# defined run functions
param.list <- set.param()
task9_path <- "/Users/marjos/Documents/courses/ki/RIKEN-bioinformatics/riken-bioinformatics-course/assigment-week-2/data/Task9"
run.task09.atlas_UMAP <- function() {
  
  plot.title.size <-param.list$font.size.list$plot.title
  legend.text.size <-param.list$font.size.list$legend.text
  legend.title.size <-param.list$font.size.list$legend.title
  encircle.size <-param.list$point.size.list$encircle
  metacell.size <-param.list$point.size.list$metacell
  white.size <-param.list$point.size.list$white
  encircle.alpha <-param.list$alpha.list$encircle

  task09.atlas_UMAP.sc.df <- fread(file.path(task9_path, '00_data/03_subcluster.tsv'))
  l1_cols <- fread(file.path(task9_path, '00_data/L1_palette.tsv'), sep='\t')
  mc_map <- fread(file.path(task9_path, '00_data/metacell_v2.metacells.csv'))
  mc <- fread(file.path(task9_path, '00_data/metacell_v2_cre.tsv.gz'))

  setkey(mc_map, barcode)
  set.seed(123)
  
  task09.atlas_UMAP.sc.df[, mc :=  mc_map[barcode, SEACell]]
  task09.atlas_UMAP.mc.df <- task09.atlas_UMAP.sc.df[, .(umap_1 = median(umap_1), umap_2 = median(umap_2)), by=.(mc, Broad_Celltype)]
  task09.atlas_UMAP.sc.subsample.df <- task09.atlas_UMAP.sc.df[, .SD[sample(1:.N, .N/50)],by=.(Broad_Celltype)]

  celltype_legend <- c('Astrocytes' = '01 Astrocytes','B' = '02 B','Blood.Endothelial' = '03 Blood.Endothelial','Chondrocyte' = '04 Chondrocyte','Cumulus' = '05 Cumulus','Dendritic' = '06 Dendritic','Epi.NonSecretory' = '07 Epi.NonSecretory','Epi.Secretory' = '08 Epi.Secretory','Fibroblast' = '09 Fibroblast',  'Lymph.Endothelial' = '10 Lymph.Endothelial',  'MAST' = '11 MAST','Melanocyte' = '12 Melanocyte',  'Microglia' = '13 Microglia',  'Mono.Mac' = '14 Mono.Mac','Neuron.GABA' = '15 Neuron.GABA','Neuron.GLU' = '16 Neuron.GLU',  'OPC' = '17 OPC',  'Oligodendrocyte' = '18 Oligodendrocyte','Plasma' = '19 Plasma','Smooth.Muscle' = '20 Smooth.Muscle','T' = '21 T')
  
  min_x <- min(task09.atlas_UMAP.sc.subsample.df$umap_1)
  min_y <- min(task09.atlas_UMAP.sc.subsample.df$umap_2)
  max_x <- max(task09.atlas_UMAP.sc.subsample.df$umap_1)
  max_y <- max(task09.atlas_UMAP.sc.subsample.df$umap_2)
  len_x <- max_x - min_x
  len_y <- max_y - min_y
  min_y <- min_y-0.1*len_y
  
  task09.atlas_UMAP.plot <- ggplot() +
    geom_point(data=task09.atlas_UMAP.sc.subsample.df, mapping=aes(umap_1, umap_2, color=Broad_Celltype), shape=16, size=11, alpha=1) +
    geom_point(data=task09.atlas_UMAP.sc.subsample.df, mapping=aes(umap_1, umap_2), fill="#ffffff", shape=21, color="NA", size=7, alpha=1) +
    geom_point(data=task09.atlas_UMAP.sc.subsample.df, aes(umap_1, umap_2, color=Broad_Celltype), shape=21, size=2, alpha=0.75) +
    geom_point(data=task09.atlas_UMAP.mc.df, mapping=aes(umap_1, umap_2), shape=21, color="#ffffff", stroke=2, size=5, alpha=1) +
    geom_point(data=task09.atlas_UMAP.mc.df, mapping=aes(umap_1, umap_2, fill=Broad_Celltype), shape=21, color="#ffffff", stroke=0.75, size=4, alpha=1) +
    theme_bw() +
    guides(
      color = guide_legend(title.position = "top", title.hjust = 0.5, override.aes = list(shape=16, size=5, alpha=1), ncol=4),
      fill = guide_colourbar()
    ) +
    scale_color_manual(values = l1_cols$hex, breaks=l1_cols$celltype, labels = celltype_legend, name = "Lv1 cell types")+
    scale_fill_manual(values = l1_cols$hex, breaks=l1_cols$celltype, labels = celltype_legend, name = "Lv1 cell types")+
    coord_fixed() +
    theme(
      plot.title = element_text(size = plot.title.size+1, hjust = 0.5),
      legend.position="top", 
      legend.text=element_text(size=legend.text.size),
      legend.title=element_text(size=legend.title.size),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank()
    ) + 
    # Add arrows starting from the minimum x and y values
    annotate("segment", x = min_x, xend = min_x + len_x*0.2, y = min_y, yend = min_y, arrow = arrow(type = "closed", length = unit(0.3, "cm"))) +
    annotate("segment", x = min_x, xend = min_x, y = min_y + len_y*0.2, yend = min_y, arrow = arrow(type = "closed", length = unit(0.3, "cm"))) +
    
    # Add text labels for UMAP 1 and UMAP 2
    annotate("text", x = min_x + 0.3, y = min_y - 0.5, size = 4, label = "UMAP 1", vjust = 0.5, hjust = 0) +
    annotate("text", x = min_x - 0.5, y = min_y + 0.3, size = 4, label = "UMAP 2", vjust = 0.5, hjust = 0, angle = 90)
  
  ggsave(file.path(RESULTS_FOLDER, "task09.atlas_UMAP.screwed_up.pdf"), task09.atlas_UMAP.plot, width=param.list$size.list$task09.atlas_UMAP$width,height=param.list$size.list$task09.atlas_UMAP$height)

  return(task09.atlas_UMAP.plot)
}
run.task09.endo_lv2_UMAP <- function() {
  
  plot.title.size <-param.list$font.size.list$plot.title
  legend.text.size <-param.list$font.size.list$legend.text
  legend.title.size <-param.list$font.size.list$legend.title
  encircle.size <-param.list$point.size.list$encircle
  metacell.size <-param.list$point.size.list$metacell
  white.size <-param.list$point.size.list$white
  encircle.alpha <-param.list$alpha.list$encircle
  
  set.seed(123)
  task09.endo_lv2_UMAP.sc.df <- fread(file.path(task9_path, '00_data/03_subcluster.tsv'))
  mc_map <- fread(file.path(task9_path, '00_data/metacell_v2.metacells.csv'))
  mc <- fread(file.path(task9_path, '00_data/metacell_v2_cre.tsv.gz'))
  setkey(mc_map, barcode)
  task09.endo_lv2_UMAP.sc.df[, mc :=  mc_map[barcode, SEACell]]
  task09.endo_lv2_UMAP.sc.df <- task09.endo_lv2_UMAP.sc.df[Broad_Celltype == 'Blood.Endothelial', .SD[sample(1:.N, min(.N, 20000))]]
  task09.endo_lv2_UMAP.mc.df <- task09.endo_lv2_UMAP.sc.df[Broad_Celltype == 'Blood.Endothelial', .(sub_umap_harmony_1 = median(sub_umap_harmony_1), sub_umap_harmony_2 = median(sub_umap_harmony_2)), by=.(mc, Broad_Celltype, Narrow_Celltype)]

  min_x <- min(task09.endo_lv2_UMAP.sc.df$sub_umap_harmony_1)
  min_y <- min(task09.endo_lv2_UMAP.sc.df$sub_umap_harmony_2)
  max_x <- max(task09.endo_lv2_UMAP.sc.df$sub_umap_harmony_1)
  max_y <- max(task09.endo_lv2_UMAP.sc.df$sub_umap_harmony_2)
  len_x <- max_x - min_x
  len_y <- max_y - min_y
  min_y <- min_y-0.1*len_y
  
  celltype_color <-c(
    'BEC.Arterial' = '#e31a1c',
    'BEC.Capilliary.1' = '#bae4b3', 'BEC.Capilliary.2' = '#74c476', 'BEC.Capilliary.3' = '#31a354', 'BEC.Capilliary.4' = '#006d2c',
    'BEC.Venous.EC.1' = '#c6dbef','BEC.Venous.EC.2' = '#9ecae1','BEC.Venous.EC.3' = '#6baed6','BEC.Venous.EC.4' = '#3182bd','BEC.Venous.EC.5' = '#08519c'
  )
  
  task09.endo_lv2_UMAP.plot <- ggplot() +
    #geom_point(data=task09.endo_lv2_UMAP.sc.df, mapping=aes(sub_umap_harmony_1, sub_umap_harmony_2, color = Narrow_Celltype), shape=16, size=0.1, alpha=0.2) +
    #geom_point(data=task09.endo_lv2_UMAP.mc.df, mapping=aes(sub_umap_harmony_1, sub_umap_harmony_2, fill = Narrow_Celltype), shape=21, alpha=0.7, size=4) +
    geom_point(data=task09.endo_lv2_UMAP.sc.df, mapping=aes(sub_umap_harmony_1, sub_umap_harmony_2, color=Narrow_Celltype), shape=16, size=11, alpha=1) +
    geom_point(data=task09.endo_lv2_UMAP.sc.df, mapping=aes(sub_umap_harmony_1, sub_umap_harmony_2), fill="#ffffff", shape=21, color="NA", size=7, alpha=1) +
    geom_point(data=task09.endo_lv2_UMAP.sc.df, aes(sub_umap_harmony_1, sub_umap_harmony_2, color=Narrow_Celltype), shape=16, size=1, alpha=0.25) +
    geom_point(data=task09.endo_lv2_UMAP.mc.df, mapping=aes(sub_umap_harmony_1, sub_umap_harmony_2), shape=16, color="#ffffff", size=8, alpha=0.75) +
    geom_point(data=task09.endo_lv2_UMAP.mc.df, mapping=aes(sub_umap_harmony_1, sub_umap_harmony_2, fill=Narrow_Celltype), shape=21, color="#000000", stroke=0.75, size=5, alpha=1) +
    theme_bw() +
    #guides(color = guide_legend(override.aes = list(shape=16, size=4), ncol=4)) +
    scale_color_manual(values=celltype_color)+
    scale_fill_manual(values=celltype_color)+
    guides(color = guide_legend(override.aes = list(size=5, shape=16, alpha=1), nrow=5, title = 'Blood.Endothelial\nLv2 Cell type'), fill=guide_none()) +
    coord_fixed() +
    theme(
      plot.margin = margin(t=0.1, r=0.25, b=-0.5, l=-0.25,  "cm"), 
      plot.title = element_text(size = plot.title.size, hjust = 0.5),
      legend.position="top", 
      legend.margin=margin(0.5,0.5,0.5,0.5),
      legend.box.margin = margin(0.5, 0.5, 0.5, 0.5),
      legend.text=element_text(size=legend.text.size),
      legend.title=element_text(size=legend.title.size),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank()
    ) + 
    # Add arrows starting from the minimum x and y values
    annotate("segment", x = min_x, xend = min_x + len_x*0.15, y = min_y, yend = min_y, arrow = arrow(type = "closed", length = unit(0.3, "cm"))) +
    annotate("segment", x = min_x, xend = min_x, y = min_y, yend = min_y + len_y*0.25, arrow = arrow(type = "closed", length = unit(0.3, "cm"))) +
    
    # Add text labels for UMAP 1 and UMAP 2
    annotate("text", x = min_x + 0.3, y = min_y - 0.5, size = 4, label = "UMAP 1", vjust = 0.5, hjust = 0) +
    annotate("text", x = min_x - 0.5, y = min_y + 0.3, size = 4, label = "UMAP 2", vjust = 0.5, hjust = 0, angle = 90)

  ggsave(file.path(RESULTS_FOLDER, "task09.endo_lv2_UMAP.screwed_up.pdf"), task09.endo_lv2_UMAP.plot, width=param.list$size.list$task09.endo_lv2_UMAP$width,height=param.list$size.list$task09.endo_lv2_UMAP$height)
  
  return(task09.endo_lv2_UMAP.plot)
}
run.task09.endo_tissue_UMAP <- function() {
  
  plot.title.size <-param.list$font.size.list$plot.title
  legend.text.size <-param.list$font.size.list$legend.text
  legend.title.size <-param.list$font.size.list$legend.title

  tissue_label.df <- data.frame(
    original = c("Bile_duct", "Bladder", "Blood", "Brain", "Cervix", "Colon (immune)", 
                 "Esophagus", "Heart", "Kidney", "Liver", "Lung", "Lymph_node", 
                 "Marrow", "Muscle", "Rectum", "Salivary Glands", "Skin", 
                 "Small_intestine", "Spleen", "Stomach", "Thymus", "Trachea", "Uterus"),
    display = c("Bile duct", "Bladder", "Blood", "Brain", "Cervix", "Colon", 
                "Esophagus", "Heart", "Kidney", "Liver", "Lung", "Lymph node", 
                "Bone marrow", "Joint", "Rectum", "Salivary gland", "Skin", 
                "Small intestine", "Spleen", "Stomach", "Thymus", "Trachea", "Ovary")
  )
  
  set.seed(123)
  task09.endo_tissue_UMAP.sc.df <- fread(file.path(task9_path, '00_data/03_subcluster.tsv'))
  tissue_cols <- fread(file.path(task9_path, '00_data/tissue_palette.tsv'))
  mc_map <- fread(file.path(task9_path, '00_data/metacell_v2.metacells.csv'))
  mc <- fread(file.path(task9_path, '00_data/metacell_v2_cre.tsv.gz'))
  setkey(mc_map, barcode)
  task09.endo_tissue_UMAP.sc.df[, mc :=  mc_map[barcode, SEACell]]
  task09.endo_tissue_UMAP.sc.df <- task09.endo_tissue_UMAP.sc.df[Broad_Celltype == 'Blood.Endothelial', .SD[sample(1:.N, min(.N, 20000))]]
  task09.endo_tissue_UMAP.mc.df <- task09.endo_tissue_UMAP.sc.df[Broad_Celltype == 'Blood.Endothelial', .(sub_umap_harmony_1 = median(sub_umap_harmony_2), sub_umap_harmony_2 = median(sub_umap_harmony_1)), by=.(mc, Broad_Celltype, Narrow_Celltype)]

  task09.endo_tissue_UMAP.sc.tissue.df <- dcast(task09.endo_tissue_UMAP.sc.df, mc~Tissue)
  task09.endo_tissue_UMAP.mc.tissue.df <- merge(task09.endo_tissue_UMAP.mc.df, task09.endo_tissue_UMAP.sc.tissue.df, by='mc')

  tissue_cols <- tissue_cols %>%
    left_join(tissue_label.df, by = c("Tissue" = "original"))
    
  min_x <- min(task09.endo_tissue_UMAP.sc.df$sub_umap_harmony_1)
  min_y <- min(task09.endo_tissue_UMAP.sc.df$sub_umap_harmony_2)
  max_x <- max(task09.endo_tissue_UMAP.sc.df$sub_umap_harmony_1)
  max_y <- max(task09.endo_tissue_UMAP.sc.df$sub_umap_harmony_2)
  len_x <- max_x - min_x
  len_y <- max_y - min_y
  min_y <- min_y-0.1*len_y
  
  celltype_color <-c(
    'BEC.Arterial' = '#e31a1c',
    'BEC.Capilliary.1' = '#bae4b3', 'BEC.Capilliary.2' = '#74c476', 'BEC.Capilliary.3' = '#31a354', 'BEC.Capilliary.4' = '#006d2c',
    'BEC.Venous.EC.1' = '#c6dbef','BEC.Venous.EC.2' = '#9ecae1','BEC.Venous.EC.3' = '#6baed6','BEC.Venous.EC.4' = '#3182bd','BEC.Venous.EC.5' = '#08519c'
  )
  
  task09.endo_tissue_UMAP.plot <- ggplot() +
    #geom_point(data=task09.endo_tissue_UMAP.sc.df, mapping=aes(sub_umap_harmony_1, sub_umap_harmony_2, color = Narrow_Celltype), shape=16, size=0.1, alpha=0.2) +
    #geom_point(data=task09.endo_tissue_UMAP.mc.df, mapping=aes(sub_umap_harmony_1, sub_umap_harmony_2, fill = Narrow_Celltype), shape=21, alpha=0.7, size=4) +
    geom_point(data=task09.endo_tissue_UMAP.sc.df, mapping=aes(sub_umap_harmony_1, sub_umap_harmony_2, color=Narrow_Celltype), shape=16, size=11, alpha=1) +
    geom_point(data=task09.endo_tissue_UMAP.sc.df, mapping=aes(sub_umap_harmony_1, sub_umap_harmony_2), fill="#ffffff", shape=21, color="NA", size=7, alpha=1) +
    geom_point(data=task09.endo_tissue_UMAP.sc.df, aes(sub_umap_harmony_1, sub_umap_harmony_2, color=Narrow_Celltype), shape=16, size=1, alpha=0.25) +
    geom_point(data=task09.endo_tissue_UMAP.mc.tissue.df, aes(x=sub_umap_harmony_2, y=sub_umap_harmony_1), shape=16, color='#ffffff', size=10, alpha=0.75) +
    geom_point(data=task09.endo_tissue_UMAP.mc.tissue.df, aes(x=sub_umap_harmony_2, y=sub_umap_harmony_1), shape=16, color='#000000', size=7, alpha=1) +
    geom_scatterpie(data=task09.endo_tissue_UMAP.mc.tissue.df, aes(x=sub_umap_harmony_2, y=sub_umap_harmony_1, group=mc, r=0.3), color='NA', cols=task09.endo_tissue_UMAP.sc.df[Broad_Celltype == 'Blood.Endothelial', sort(unique(Tissue))], alpha=1) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(size=5, shape=16, alpha=1), title = 'Tissue'), color=guide_none()) +
    scale_fill_manual(values = tissue_cols$hex, breaks=tissue_cols$Tissue, labels = tissue_cols$display) +
    scale_color_manual(values=celltype_color)+
    coord_fixed() +
    theme(
      plot.margin = margin(t=0.1, r=0.15, b=-0.5, l=-0.25,  "cm"), 
      plot.title = element_text(size = plot.title.size, hjust = 0.5),
      legend.position="top", 
      legend.margin=margin(0.5,0.5,0.5,0.5),
      legend.box.margin = margin(0.5, 0.5, 0.5, 0.5),
      legend.text=element_text(size=legend.text.size),
      legend.title=element_text(size=legend.title.size),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank()
    ) + 
    # Add arrows starting from the minimum x and y values
    annotate("segment", x = min_x, xend = min_x + len_x*0.15, y = min_y, yend = min_y, arrow = arrow(type = "closed", length = unit(0.3, "cm"))) +
    annotate("segment", x = min_x, xend = min_x, y = min_y, yend = min_y + len_y*0.25, arrow = arrow(type = "closed", length = unit(0.3, "cm"))) +
    
    # Add text labels for UMAP 1 and UMAP 2
    annotate("text", x = min_x + 0.3, y = min_y - 0.5, size = 4, label = "UMAP 1", vjust = 0.5, hjust = 0) +
    annotate("text", x = min_x - 0.5, y = min_y + 0.3, size = 4, label = "UMAP 2", vjust = 0.5, hjust = 0, angle = 90)
  
  ggsave(file.path(RESULTS_FOLDER, "task09.endo_tissue_UMAP.screwed_up.pdf"), task09.endo_tissue_UMAP.plot, width=param.list$size.list$task09.endo_tissue_UMAP$width,height=param.list$size.list$task09.endo_tissue_UMAP$height)
  
  return(task09.endo_tissue_UMAP.plot)
}
################################################################# run the functions
task09.atlas_UMAP.plot <- run.task09.atlas_UMAP()
task09.endo_lv2_UMAP.plot <- run.task09.endo_lv2_UMAP()
task09.endo_tissue_UMAP.plot <- run.task09.endo_tissue_UMAP()
  



######## umap per tissue ########
## alluvial plot

install.packages("ggalluvial")
library(ggalluvial)

df <- task09.atlas_UMAP.sc.df %>%
  group_by(Tissue, Broad_Celltype) %>%
  summarise(Freq = n(), .groups = "drop") 


ggplot(df, aes(axis1 = Tissue, axis2 = Broad_Celltype, y = Freq)) +
  geom_alluvium(aes(fill = Broad_Celltype), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size=3) +
  scale_x_discrete(limits = c("Tissue", "Broad Cell Type"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Alluvial Plot of Tissue and Broad Cell Type") +
  theme_minimal() +
  theme(legend.position = "none")

# stacked plot with the celltypes proportion per tissue
# Step 1: Calculate proportion of cell types per tissue
df_prop <- task09.atlas_UMAP.sc.df %>%
  group_by(Tissue, Broad_Celltype) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Tissue) %>%
  mutate(Prop = Count / sum(Count)) %>%
  ungroup()

# Step 2: Create named rainbow color palette
celltypes <- unique(df_prop$Broad_Celltype)
n_celltypes <- length(celltypes)
rainbow_colors <- setNames(c("#96b06a",
"#9458ca",
"#61b444",
"#c856c0",
"#b3b235",
"#5970d8",
"#d99a4a",
"#5e99d4",
"#d15135",
"#4bbd80",
"#dd3a6b",
"#4dbfb7",
"#cc4b95",
"#397f4d",
"#d091d5",
"#707329",
"#7563a8",
"#a26333",
"#9c5883",
"#e58487",
"#ac495b"), celltypes)

# Step 3: Plot with custom colors
ggplot(df_prop, aes(x = Tissue, y = Prop, fill = Broad_Celltype)) +
  geom_bar(stat = "identity", position = "stack", color = "white", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = rainbow_colors) +
  ggtitle("Cell type composition per tissue") +
  ylab("Proportion of cells") +
  xlab(NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.6, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(file.path(RESULTS_FOLDER, "task09.atlas_UMAP.stacked_barplot.pdf"), width=10,height=6)


######

df_prop <- task09.atlas_UMAP.sc.df %>%
  group_by(Tissue, Broad_Celltype) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Tissue) %>%
  mutate(Prop = Count / sum(Count)) %>%
  ungroup()

ggplot(df_prop, aes(x = Tissue, y = Prop, fill = Broad_Celltype)) +
  geom_bar(stat = "identity", position = "stack", color = "white", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = rainbow_colors) +
  ggtitle("Cell type composition per tissue") +
  ylab("Proportion of cells") +
  xlab(NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.6, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank()
  )
