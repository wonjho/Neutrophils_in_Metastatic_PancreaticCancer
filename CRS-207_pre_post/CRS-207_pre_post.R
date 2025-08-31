#Won Jin Ho Laboratory
#Johns Hopkins University School of Medicine
#R code for IMC data analysis of patient liver metastatic lesions
#Platform aarch64-apple-darwin20
#R version 4.3.2 (2023-10-31)

rm(list = ls())


####READ and CLUSTER FUNCTIONS####
returnfcs <- function(FDR_cutoff=.05,
                      metaDataFile=NULL,
                      panelDataFile=NULL,
                      dataDirectory=NULL,
                      shape_conditions=NULL,
                      color_conditions=NULL){
  ## This function generates an fcs file, subtype_markers, colors and shapes for clustering 
  require(scales);require(readxl);require(dplyr);require(flowCore)
  ## Directory and metadatafile checking
  if(!dir.exists(dataDirectory)) {stop('ERR: cannot find data directory')}
  if(!file.exists(metaDataFile)) {stop('ERR: cannot find metadata.xlsx or .csv file')}
  ## Read-in metadata and clean
  ifelse(grepl(metaDataFile,pattern='.xls'),md <- read_excel(metaDataFile),md <- read.csv(metaDataFile,header = TRUE))#must be in xl format or csv
  md$File_order <- factor(md$File_order)
  md$Treatment <- factor(md$Treatment)
  rownames(md) = md$sample_id
  md$sample_id <- md$sample_id
  
  ## Make sure all files in metadata present in datadirectory
  if(!all(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])){
    print(paste('ERR: not all filenames in metadata present in data folder - missing',md$file_name[!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])],'Subsetting...'))
    md <- md[-c(!which(md$file_name %in% list.files(dataDirectory)[grep(list.files(dataDirectory),pattern = '.fcs')])),]
  }
  
  ## Read fcs into fcs_raw
  fcs_raw <- read.flowSet(paste0(dataDirectory,"/",md$file_name), transformation = FALSE, truncate_max_range = FALSE)
  panel <- read_excel(paste0(work,'/Config/panel.xlsx'))
  head(data.frame(panel))
  panel$Parameter <- gsub('-', '_', panel$Parameter)
  
  ## Export out the parameter/panel data from the flowFrame to edit
  ## use panel$Antigen to fix description in panel_fcs
  ## use metal+isotope as mapping between panel from xlsx and panel from the fcs files
  
  panel_fcs <- pData(parameters(fcs_raw[[1]]))
  panel_fcs$desc <- gsub('-', '_', panel_fcs$desc)
  panel_fcs$desc[is.na(panel_fcs$desc)] <- paste0('NA_',which(is.na(panel_fcs$desc)))  
  
  rownames(panel_fcs) = panel_fcs$name
  
  ## Replace desc with revised Name
  panel_fcs[panel$Parameter,]$desc<-panel$Name
  
  ## Replace parameter data in flowSet with edits
  pData(parameters(fcs_raw[[1]])) <- panel_fcs
  
  ## Assign objects to marker lists
  subtype_markers <- panel$Name[panel$Subtype == 1]
  functional_markers <- panel$Name[panel$Functional == 1]
  otherparameters <- panel$Name[panel$Other ==1]
  cluster_by <- panel$Name[panel$Cluster == 1]
  
  ## Check marker lists
  if(!all(subtype_markers %in% panel_fcs$desc)){stop('ERR: Not all subtype_markers in panel_fcs$desc')}
  if(!all(functional_markers %in% panel_fcs$desc)){stop('ERR: Not all functional_markers in panel_fcs$desc')}
  if(!all(otherparameters %in% panel_fcs$desc)){stop('ERR: Not all otherparameters in panel_fcs$desc')}
  if(!all(cluster_by %in% panel_fcs$desc)){stop('ERR: Not all cluster markers in panel_fcs$desc')}
  
  fcs <- fsApply(fcs_raw, function(x, cofactor = 0.8){
    colnames(x) <- panel_fcs$desc
    expr <- exprs(x)
    expr <- asinh(expr[, union(subtype_markers,functional_markers)] / cofactor)
    exprs(x) <- expr
    x
  })
  
  ## Save out the original rownames from the parameter list from the fcs flowFrame
  panellist <- rownames(pData(parameters(fcs[[1]])))
  
  ## Create dummy list to save all expression data
  exprTr_list<-c()
  
  ## Save arc sin transformed expression + spatial data for each flowFrame
  for(i in 1:length(md$file_name)){
    
    ## Expression data is combined with spatial parameter data
    
    exprRaw<-exprs(fcs_raw[[i]])
    
    colnames(exprRaw)<-panel_fcs$desc
    
    expr<-cbind(exprs(fcs[[i]])[, union(subtype_markers,functional_markers)],exprRaw[,otherparameters])
    
    ## Combine other (spatial) data with the protein data
    colnames(expr)<-c(colnames(exprs(fcs[[i]])),otherparameters)
    
    ## Add to list
    exprTr_list[[i]]<-expr
    
  }
  
  ## Create a flowSet based on the list of expression data
  fcs1<-flowSet(sapply(exprTr_list,flowFrame))
  
  ## Change parameter rownames
  panel_fcs1 <- pData(parameters(fcs1[[1]]))
  rownames(pData(parameters(fcs1[[1]]))) <- rownames(panel_fcs[panel_fcs$desc %in% pData(parameters(fcs1[[1]]))$name,])
  
  ###to scale every flowframe
  ## Save out the original rownames from the parameter list from the fcs flowFrame
  panellist <- rownames(pData(parameters(fcs[[1]])))
  
  ## Create dummy list to save all expression data
  exprTr_list<-c()
  
  ## Save arc sin transformed expression + spatial data for each flowFrame
  for(i in 1:length(md$file_name)){
    
    ## Expression data is combined with spatial parameter data
    
    expr<-exprs(fcs[[i]])
    
    expr<-t(scale(t(expr)))
    
    ## Add to list
    exprTr_list[[i]]<-expr
    
  }
  
  ## Create a flowSet based on the list of expression data
  fcs2<-flowSet(sapply(exprTr_list,flowFrame))
  
  ## Get sample ids
  sample_ids <- rep(md$sample_id, fsApply(fcs1, nrow))
  
  ## Return: 
  ## fcs (only marker expressions arcsin transformed), 
  ## fcs1 (arcsin transformed + spatial parameters), 
  ## fcs2 (scaled arcsin expr per flowframe)
  ## fcsraw (all raw data), and all marker/parameter lists
  return(list('fcs'=fcs,'fcs1'=fcs1,'fcs2'=fcs2,'fcsraw'=fcs_raw,'subtype_markers'=subtype_markers,'functional_markers'=functional_markers,'otherparameters'=otherparameters,'cluster_by'=cluster_by,'sample_ids'=sample_ids,'meta_data'=md))
}

clusterfcs <- function(fcs=output$fcs,
                       cluster_by = output$cluster_by,
                       seed=1234,plottitle='consensus_plots',
                       scaleoption=F,
                       numclusters=40){
  ## Cell population identification with FlowSOM and ConsensusClusterPlus
  require(dplyr);require(FlowSOM);require(ConsensusClusterPlus)
  set.seed(seed)
  som <- ReadInput(fcs, transform = FALSE, scale = scaleoption) %>% BuildSOM(colsToUse = cluster_by)
  
  ## Get the cell clustering into 100 SOM codes
  cell_clustering_som <- som$map$mapping[,1]
  
  ## Metaclustering into numclusters with ConsensusClusterPlus
  codes <- som$map$codes
  mc <- ConsensusClusterPlus(t(codes), maxK = numclusters, reps = 100,
                             pItem = 0.9, pFeature = 1, title = plottitle, 
                             plot = "png", clusterAlg = "hc", 
                             innerLinkage = "average", finalLinkage = "average",
                             distance = "euclidean", seed = 1234)
  
  ## Get cluster ids for each cell
  code_clustering <- mc[[numclusters]]$consensusClass#metaclusters consensus
  cell_clustering <- code_clustering[cell_clustering_som]#cell clustering from som
  return(list('code_clustering'=code_clustering,'cell_clustering'=cell_clustering,'metaclusters'=mc))
}


####DIAGNOSTIC HEATMAP FUNCTIONS ####
plot_clustering_heatmap_wrapper <- function(fcs, cell_clustering, nclusters=40,
                                            color_clusters=colorassigned, cluster_merging = NULL, 
                                            cluster_by=output$cluster_by,
                                            clusterMergeFile=NULL,
                                            fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales);require(pals);require(ComplexHeatmap)
  ## Will output the heatmap object and print it 
  
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,cluster_by]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,cluster_by]
  
  
  ## Calculate the mean expression##################################################
  
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  ## Annotation for the merged clusters
  
  if(!is.null(clusterMergeFile)){
    ifelse(grepl(clusterMergeFile,pattern='.xls'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Merged <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Merged <- color_clusters2
  }
  
  ## Colors for the heatmap
  
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  
  clusternames<-sort(unique(cluster_merging$new_cluster))
  
  colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
  
  names(colorassigned)<-clusternames
  
  color_list = list(clusters=colorassigned)
  
  color_list_byoriginal = colorassigned[match(cluster_merging$new_cluster,names(colorassigned))]
  
  cp<-rowAnnotation(clusters=cluster_merging$new_cluster,
                    col=color_list,
                    gp = gpar(col = "white", lwd = .5),
                    prop=anno_barplot(
                      clustering_prop, 
                      gp = gpar(fill=color_list_byoriginal, col=F),
                      border = F,
                      bar_width = 0.75, 
                      width = unit(2,"cm")))
  
  q <- Heatmap(expr_heat, name="scaled",
               col=rev(brewer.rdbu(100)),
               row_order = cluster_merging[order(cluster_merging$new_cluster),]$original_cluster,
               cluster_columns = T,
               cluster_rows = T,
               border = NA,
               rect_gp = gpar(col = "white", lwd = .5),
               right_annotation = cp,
               show_row_names = T,
               row_names_gp = gpar(fontsize=7),
               column_names_gp = gpar(fontsize=10),
               heatmap_legend_param = list(at=seq(from = 0, to = 1, by = 0.2)),
               width = unit(10, "cm"))
  
  print('Colors:')
  print(color_clusters)
  
  pdf(fileName, width=8, height=6) 
  
  return(q)
  
  dev.off() 
  
}

plot_clustering_heatmap_wrapper2 <- function(fcs, cell_clustering, nclusters=40,
                                             color_clusters=colorassigned,
                                             colorbar=rev(brewer.rdbu(100)),
                                             subtype_markers,
                                             fileName = 'clusteringheatmap.pdf'){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales)
  ## Will output the heatmap object and print it 
  #get expression
  expr <- fsApply(fcs, exprs);expr <-expr[,subtype_markers]
  ## Scale expression of all markers to values between 0 and 1
  rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0; expr01[expr01 > 1] <- 1;expr01 <-expr01[,subtype_markers]
  ## Calculate the mean expression##################################################
  pdf(fileName, width=8, height=11) 
  expr_mean <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  expr01_mean <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(mean))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_mean[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_mean[, colnames(expr01)])
  rownames(expr_heat) <- expr01_mean$cell_clustering
  
  ## Colors for the heatmap
  
  labels_row <- paste0(expr01_mean$cell_clustering, " ")
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_mean$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  p <- pheatmap(expr_heat, 
                color = colorbar, 
                cluster_cols = F,
                cluster_rows = F, 
                labels_row = labels_row,
                #scale="row",
                display_numbers = F, 
                number_color = "black",
                fontsize = 9, fontsize_number = 6,  
                #legend_breaks = legend_breaks,
                annotation_row = annotation_row, 
                annotation_colors = annotation_colors,
                cellwidth = 8,
                cellheight = 8,
                border_color = "white",
                annotation_legend = F
  )
  dev.off() 
  print('Colors:')
  print(color_clusters)
  print(p);return(p)
}


####REQUIRED LIBRARIES####
library(readxl)
library(stringr)
library(matrixcalc)
library(matrixStats)
library(Hmisc)
library(reshape2)
library(dplyr)
library(plotrix)
library(multcomp)
library(flowCore)
library(sf)
library(clusterSim)
library(limma)
library(corrplot)
library(packcircles)
library(cowplot)
library(autoimage)
library(ggplot2)
library(ggpubr)
library(ggiraphExtra)
library(ggridges)
library(gridExtra)
library(igraph)
library(qgraph)
library(circlize)
library(scales)
library(RColorBrewer)
library(ComplexHeatmap)
library(pheatmap)
library(pals)
library(plot3D)
library(akima)
library(basetheme)


####DATA LOADING####
setwd("~/Library/CloudStorage/Dropbox/Hopkins/Independent study/Jaffee Lab/Manuscripts/Neutrophils in metastatic disease/Codes/CRS-207_pre_post") ## Set working directory
work<-getwd()
output <- returnfcs(metaDataFile = paste0(work,"/Config/metadata.xlsx"),
                    panelDataFile = paste0(work,"/Config/panel.xlsx"),
                    dataDirectory = paste0(work,"/Data"))

## Set up levels
samplevels=c(output$meta_data$sample_id)
fileorderlevels=c(output$meta_data$File_order)
treatmentlevels=c("Pre","Post")


####DIAGNOSTICS####
## Spot check - number of cells per sample
cell_table <- table(output$sample_ids)
ggdf <- data.frame(sample_id = names(cell_table), 
                   cell_counts = as.numeric(cell_table))
ggdf$Treatment <- factor(output$meta_data$Treatment[match(ggdf$sample_id,output$meta_data$sample_id)], levels = treatmentlevels)
ggdf$File_order <- factor(output$meta_data$File_order[match(ggdf$sample_id,output$meta_data$sample_id)], levels = fileorderlevels)
ggp<-ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = Treatment)) + 
  geom_bar(stat = 'identity') + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=5)) +  
  scale_x_discrete(drop = FALSE)
pdf('Diagnostics_cellcounts.pdf',width=6, height=4);ggp; dev.off()

## Multi-dimensional scaling plot to show similarities between samples
## Get the mean marker expression per sample
expr_mean_sample_tbl <- data.frame(sample_id = output$sample_ids, fsApply(output$fcs,exprs)) %>%
  group_by(sample_id) %>%  summarize_all(funs(mean))
expr_mean_sample <- t(expr_mean_sample_tbl[, -1])
colnames(expr_mean_sample) <- expr_mean_sample_tbl$sample_id
mds <- plotMDS(expr_mean_sample, plot = FALSE)
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
                   sample_id = colnames(expr_mean_sample))
ggdf$Treatment <- factor(output$meta_data$Treatment[match(ggdf$sample_id,output$meta_data$sample_id)], levels = treatmentlevels)
ggdf$File_order <- factor(output$meta_data$File_order[match(ggdf$sample_id,output$meta_data$sample_id)], levels = fileorderlevels)
ggp<-ggplot(ggdf, aes(x = MDS1, y = MDS2, color = Treatment)) +
  geom_point(size = 2.5) +
  theme_bw()+
  theme(plot.background = element_rect(fill="white"),
        panel.background = element_rect(fill="white"),
        panel.grid = element_blank(),
        axis.line = element_line(color="black"),
        axis.text = element_text(color="black"),
        axis.title = element_text(color="black"),
        legend.background = element_rect(fill="white"),
        legend.key = element_rect(fill="white"),
        legend.text = element_text(color="white"))
pdf('Diagnostics_MDS_sample.pdf',width=6, height=6);ggp; dev.off()

## Colors for the heatmap
color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
pdf('Diagnostics_Heatmap.pdf',width=8, height=8)
pheatmap(expr_mean_sample_tbl[,output$cluster_by], color = color, display_numbers = TRUE,
         number_color = "black", fontsize_number = 3, 
         annotation_colors = annotation_colors, clustering_method = "average")
dev.off()


####CLUSTERING####
##Revised loading depending on the diagnostics if needed
output <- returnfcs(metaDataFile = paste0(work,"/Config/metadata.xlsx"),
                    panelDataFile = paste0(work,"/Config/panel.xlsx"),
                    dataDirectory = paste0(work,"/Data"))

##Clustering
output[(length(output)+1):(length(output)+3)] <- clusterfcs(fcs=output$fcs, numclusters=40, scaleoption = T) ## Set cluster number
names(output)[(length(output)-2):(length(output))] <- c('code_clustering','cell_clustering','metaclusters')


####ANNOTATIONS OF CLUSTERS####
## Load merge file
## Assign colors
clusterMergeFile = paste0(work,"/Config/merge.xlsx")
cluster_merging <- read_excel(clusterMergeFile)
clusterlevels=c(1:22) ## Set cluster number
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
clusternames<-clusterlevels
names(colorassigned)<-clusternames
mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)
cell_clustering1m <- cluster_merging$new_cluster[mm1]
output$cell_clustering1m <- cell_clustering1m

## Metacluster heatmaps
plot_clustering_heatmap_wrapper(fcs=output$fcs,
                                color_clusters = kovesi.rainbow_bgyrm_35_85_c69(50),
                                cell_clustering = output$cell_clustering, 
                                cluster_by=output$cluster_by,
                                clusterMergeFile = clusterMergeFile,
                                fileName = 'Clusteringheatmap_all.pdf'); dev.off()

clusterMergeFile = paste0(work,"/Config/merge.xlsx")
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c("Neutrophil_I", "Neutrophil_II", "Myeloid_I", "Myeloid_II", "Myeloid_III", "Myeloid_IV", "Myeloid_V", "Myeloid_VI", "Myeloid_VII", "Myeloid_VIII",
                "CD8_T_I", "CD8_T_II", "CD8_T_III", "CD4_T", "Treg", "B", "Plasma", "Tumor_I", "Tumor_II", "Tumor_III", "Stroma", "UA")

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
clusternames<-clusterlevels
names(colorassigned)<-clusternames
mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)
cell_clustering1m <- cluster_merging$new_cluster[mm1]
output$cell_clustering1m <- cell_clustering1m

plot_clustering_heatmap_wrapper2(fcs=output$fcs,
                                 colorbar = kovesi.diverging_bwr_40_95_c42(100),
                                 subtype_markers = c("CD45", "CD45RA", "CD45RO", "MPO", "CD15", "CD68", "CD163", "HLADR", "DC_SIGN", "CD206", "LAIR_1", "iNOS",
                                                     "CD8", "CD4", "FOXP3", "GZMB", "IFNy", "TNFa", "PD1", "Tox_Tox2", "CD69", "CD20", "CD138",
                                                     "CK", "Ecadherin", "pSTAT1", "pSTAT3", "CD103", "PTPN22", "KI67",
                                                     "Collagen", "SMA", "VIM", "Podoplanin"),
                                color_clusters = colorassigned,
                                cell_clustering = factor(output$cell_clustering1m, levels=clusterlevels), 
                                fileName = 'Clusteringheatmap_merged.pdf') ###FIGURE S12A###
dev.off()

## Save output list
saveRDS(output, file="CRS-207_pre_post_backup_output.rds")


####DENSITY PLOTS####
## Proportion calculations
counts_table <- table(output$cell_clustering1m, output$sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

areas <- read_xlsx(paste0(work,'/Config/areas.xlsx'))
densities <- t(t(counts)/areas$TotalArea)

write.csv(densities, 'Results_densities.csv') ###FIGURE 5C AND FIGURE S12B###

## Set up the data frame for proportional plotting
ggdf <- melt(data.frame(cluster = rownames(props),props, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "sample_id")
ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdf$cluster <- factor(ggdf$cluster, levels=clusterlevels)
ggdf$Treatment <- factor(output$meta_data$Treatment[match(ggdf$sample_id,output$meta_data$sample_id)], levels = treatmentlevels)
ggdf$File_order <- factor(output$meta_data$File_order[match(ggdf$sample_id,output$meta_data$sample_id)], levels = fileorderlevels)

ggdfd <- melt(data.frame(cluster = rownames(densities),densities, check.names = FALSE),
              id.vars = "cluster", value.name = "densities", 
              variable.name = "sample_id")
ggdfd$sample_id <- factor(ggdfd$sample_id, levels=samplevels)
ggdfd$cluster <- factor(ggdf$cluster, levels=clusterlevels)
ggdfd$Treatment <- factor(output$meta_data$Treatment[match(ggdfd$sample_id,output$meta_data$sample_id)], levels=treatmentlevels)

## DENSITY PLOTS
ggp2<-ggplot(ggdfd,aes(x=paste0(Treatment),y=densities,fill=paste0(Treatment)))+
  geom_boxplot(outlier.color=NA, lwd=0.25)+
  geom_jitter(width=0, size=.5)+
  facet_wrap(~cluster,ncol=5,scales="free")+
  ylab("Density of Cells")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.25, color="black"),
        axis.line.y = element_line(size=0.25, color="black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=7, color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white")  
  )+
  scale_fill_manual(values=cbbPalette)

pdf('Density_box.pdf',width=8,height=6) ###FIGURE 5C AND FIGURE S12B###
ggp2
dev.off()


####NETWORK ANALYSIS####
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
hex <- hue_pal()(9)
clusternames<-clusterlevels
names(colorassigned)<-clusternames

allcelltypes<-clusterlevels
legendctype<-as.data.frame(cbind(paste0("ctype",1:length(allcelltypes)),allcelltypes))
legendctype$maintype<-1
legendctype$maintype[str_detect(legendctype$allcelltypes,"Neutrophil")]<-"Neutrophil"
legendctype$maintype[str_detect(legendctype$allcelltypes,"T")]<-"Lym"
legendctype$maintype[str_detect(legendctype$allcelltypes,"B")]<-"Lym"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Plasma")]<-"Lym"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Tumor")]<-"Tumor"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Myeloid")]<-"Myeloid"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Stroma")]<-"Stroma"
legendctype$maintype[str_detect(legendctype$allcelltypes,"UA")]<-"Other"

## Calculate the number of each cell type in the dataset
ggdft <- melt(data.frame(cluster = rownames(counts), counts, check.names = FALSE),
              id.vars = "cluster", value.name = "counts", 
              variable.name = "sample_id")
ggdft$sample_id <- factor(ggdft$sample_id, levels=samplevels)
ggdft$cluster <- factor(ggdft$cluster, levels=clusterlevels)
ggdft$Treatment <- factor(output$meta_data$Treatment[match(ggdft$sample_id,output$meta_data$sample_id)], levels=treatmentlevels)

totalcounts<-ggdft %>% group_by(cluster, Treatment) %>% summarize_at(vars(counts),funs(sum))
write.csv(totalcounts,"Totalcounts.csv")

totalcounts_celltype <- ggdft %>% group_by(cluster) %>% summarize_at(vars(counts),funs(sum))
write.csv(totalcounts_celltype, "Totalcounts_celltype.csv")

totalcounts_treatment <- ggdft %>% group_by(cluster, Treatment) %>% summarize_at(vars(counts),funs(sum))
write.csv(totalcounts_treatment, "Totalcounts_treatment.csv")

## Calculate the percentage of each respective total
totalcounts_pre <- totalcounts[totalcounts$Treatment =="Pre",]
totalcounts_post <- totalcounts[totalcounts$Treatment =="Post",]

totalpre<-sum(totalcounts_pre$counts)
totalpost<-sum(totalcounts_post$counts)

pct_pre <- totalcounts_pre$counts/totalpre*100
pct_post <- totalcounts_post$counts/totalpost*100

names(pct_pre)<- totalcounts_pre$cluster
names(pct_post)<- totalcounts_post$cluster

## Identify which cell types in either of the data subsets are very rare - have less than 0.01%
exclude_pre<-which(pct_pre<.001)
exclude_pre<-legendctype$V1[match(names(exclude_pre), legendctype$allcelltypes)]

exclude_post<-which(pct_post<.001)
exclude_post<-legendctype$V1[match(names(exclude_post), legendctype$allcelltypes)]

## Get out expression levels, X, and Y coordinates
expr <- fsApply(output$fcs1, exprs) #create expression matrix

expr0<-data.frame(expr[,c(union(output$subtype_markers,output$functional_markers),"CellId","X_position","Y_position")],
                  cluster=output$cell_clustering1m,
                  sample_id=output$sample_ids)
expr0$Treatment <- factor(output$meta_data$Treatment[match(expr0$sample_id,output$meta_data$sample_id)], levels=treatmentlevels)

expr0<-expr0[expr0$cluster!="UA",]
expr0_pre<-expr0[expr0$Treatment=="Pre",]
expr0_post<-expr0[expr0$Treatment=="Post",]

## Create distance matrices
PREsamp <-unique(expr0_pre$sample_id)
POSTsamp <-unique(expr0_post$sample_id)

## PRE
expr_pre<-c()

for(k in 1:length(PREsamp)){
  expr_k<-expr0[expr0$sample_id==PREsamp[k],] 
  
  #create placer cols for shortest distance to each cell type
  dummy<-matrix(nrow=nrow(expr_k),ncol=length(allcelltypes)) 
  colnames(dummy)<-paste0("ctype",1:length(allcelltypes))
  expr_k<-data.frame(expr_k,d=dummy)
  
  #get distances from index cell (row) to all cells
  dist1<-as.matrix(dist(cbind(expr_k$X_position,expr_k$Y_position)))
  dist1<-as.data.frame(dist1)
  
  for(i in 1:nrow(expr_k)){
    #looking at cell i in core k
    celldistances1<-data.frame(dist=dist1[,i],cluster=expr_k$cluster)
    celldistances1$clustercode<-match(celldistances1$cluster,legendctype$allcelltypes)
    for(j in 1:length(allcelltypes)){
      #among all cells of cell type j, what non-zero (non-self) distance is the closest
      min_d_j<-min(celldistances1[celldistances1$cluster==allcelltypes[j],]$dist[celldistances1[celldistances1$cluster==allcelltypes[j],]$dist!=0])
      expr_k[i,paste0("d.ctype",j)]<-min_d_j
    }
  }
  
  expr_k$ctype_no<-paste0("ctype",match(expr_k$cluster,legendctype$allcelltypes))
  
  expr_k_m<-as.matrix(expr_k[,colnames(expr_k)[str_detect(colnames(expr_k),"d.c")]])
  rownames(expr_k_m)<-expr_k$ctype_no
  colnames(expr_k_m)<-legendctype$V1
  expr_k_m[is.infinite(expr_k_m)]<-NA
  
  expr_pre<-rbind(expr_pre,expr_k_m)
}

saveRDS(expr_pre,'CRS-207_pre_post_backup_dist_pre.rds')

#create expression + distance data frames
expr_pre_combined <- cbind(expr0_pre,expr_pre)

#improving robustness of distance relationships in the dataset by removing cell types without relationships or very rare cell types that would be overrepresented (sampling bias)

#remove cell type columns where there are no cells
expr_pre_rm0<-expr_pre[,colnames(expr_pre) %nin% colnames(expr_pre)[colSums(expr_pre, na.rm = T)==0]]

#remove rows/columns where there are very rare (<0.01%) cell types
expr_pre_rm<-expr_pre_rm0[rownames(expr_pre_rm0) %nin% exclude_pre, colnames(expr_pre_rm0) %nin% exclude_pre]


## POST
expr_post<-c()

for(k in 1:length(POSTsamp)){
  expr_k<-expr0[expr0$sample_id==POSTsamp[k],] 
  
  #create placer cols for shortest distance to each cell type
  dummy<-matrix(nrow=nrow(expr_k),ncol=length(allcelltypes)) 
  colnames(dummy)<-paste0("ctype",1:length(allcelltypes))
  expr_k<-data.frame(expr_k,d=dummy)
  
  #get distances from index cell (row) to all cells
  dist1<-as.matrix(dist(cbind(expr_k$X_position,expr_k$Y_position)))
  dist1<-as.data.frame(dist1)
  
  for(i in 1:nrow(expr_k)){
    #looking at cell i in core k
    celldistances1<-data.frame(dist=dist1[,i],cluster=expr_k$cluster)
    celldistances1$clustercode<-match(celldistances1$cluster,legendctype$allcelltypes)
    for(j in 1:length(allcelltypes)){
      #among all cells of cell type j, what non-zero (non-self) distance is the closest
      min_d_j<-min(celldistances1[celldistances1$cluster==allcelltypes[j],]$dist[celldistances1[celldistances1$cluster==allcelltypes[j],]$dist!=0])
      expr_k[i,paste0("d.ctype",j)]<-min_d_j
    }
  }
  
  expr_k$ctype_no<-paste0("ctype",match(expr_k$cluster,legendctype$allcelltypes))
  
  expr_k_m<-as.matrix(expr_k[,colnames(expr_k)[str_detect(colnames(expr_k),"d.c")]])
  rownames(expr_k_m)<-expr_k$ctype_no
  colnames(expr_k_m)<-legendctype$V1
  expr_k_m[is.infinite(expr_k_m)]<-NA
  
  expr_post<-rbind(expr_post,expr_k_m)
}

saveRDS(expr_post,'CRS-207_pre_post_backup_dist_post.rds')

#create expression + distance data frames
expr_post_combined <- cbind(expr0_post,expr_post)

#improving robustness of distance relationships in the dataset by removing cell types without relationships or very rare cell types that would be overrepresented (sampling bias)

#remove cell type columns where there are no cells
expr_post_rm0<-expr_post[,colnames(expr_post) %nin% colnames(expr_post)[colSums(expr_post, na.rm = T)==0]]

#remove rows/columns where there are very rare (<0.01%) cell types
expr_post_rm<-expr_post_rm0[rownames(expr_post_rm0) %nin% exclude_post, colnames(expr_post_rm0) %nin% exclude_post]


####NETWORK VISUALIZATION####
###by a distance object (of matrix) of shortest distance between cell types

## PRE
mat_pre=aggregate(x=expr_pre_rm, by=list(rownames(expr_pre_rm)), FUN=mean, na.rm=T)
groupnames<-mat_pre$Group.1
mat_pre<-as.matrix(mat_pre[,2:ncol((mat_pre))])
rownames(mat_pre)<-str_sub(groupnames,6,) # add the rownames
colnames(mat_pre)<-str_sub(colnames(mat_pre),6,) #simplify colnames
dist_pre<-mat_pre

## POST
mat_post=aggregate(x=expr_post_rm, by=list(rownames(expr_post_rm)), FUN=mean, na.rm=T)
groupnames<-mat_post$Group.1
mat_post<-as.matrix(mat_post[,2:ncol((mat_post))])
rownames(mat_post)<-str_sub(groupnames,6,) # add the rownames
colnames(mat_post)<-str_sub(colnames(mat_post),6,) #simplify colnames
dist_post<-mat_post

#color of legends revised
cols <- cbbPalette
colorlist <- cols[as.numeric(as.factor(legendctype$maintype))]

#generate visualizations of distance relationships using network qgraph
pdf("Distance_Relationships_Pre.pdf",width=8,height=6) ###FIGURE 5D###

xx<-dist_pre
yy<-xx[as.character(sort(as.numeric(rownames(xx)))),#make sure all rows and cols are ordered
       as.character(sort(as.numeric(colnames(xx))))]

colorlist_pre<-colorlist[as.numeric(colnames(mat_pre))] #set color
names(colorlist_pre)<-as.character(colnames(mat_pre)) #set color reference names

g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_pre[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy,
          title="Pre",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_pre, 1/dist_post, na.rm=T), #max distance has to be chosen properly to scale across all plots (need to have defined dist parameters first)
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlist_pre[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype)), 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)
dev.off()

g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_pre[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="Pre (to each cluster node)",
          layout='spring', repulsion=10, 
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlist_pre[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"], 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)
dev.off()

#generate visualizations of distance relationships using network qgraph
pdf("Distance_Relationships_Post.pdf",width=8,height=6) ###FIGURE 5E###

xx<-dist_post
yy<-xx[as.character(sort(as.numeric(rownames(xx)))),#make sure all rows and cols are ordered
       as.character(sort(as.numeric(colnames(xx))))]

colorlist_post<-colorlist[as.numeric(colnames(mat_post))] #set color
names(colorlist_post)<-as.character(colnames(mat_post)) #set color reference names

g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_post[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy,
          title="Post",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_post, 1/dist_post, na.rm=T), #max distance has to be chosen properly to scale across all plots (need to have defined dist parameters first)
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlist_post[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype)), 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)
dev.off()

g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_post[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="Post (to each cluster node)",
          layout='spring', repulsion=10, 
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlist_post[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"], 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)
dev.off()