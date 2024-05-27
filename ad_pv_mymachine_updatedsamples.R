#-----------------------------------------------------------------------------------------------------#
#Expanded PV v AD
#-----------------------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------------------------#
#Load the packages I need hinny
#-----------------------------------------------------------------------------------------------------#

#library(future)
#plan("multisession", workers = 4) # Enable parallelization from future package


library(ggpubr)
library(Seurat)
library(ggplot2)
library(future)
library(RColorBrewer)
library(openxlsx)
library(raster)
library(hdf5r)
library(future)
library(R.utils)
library(monocle3)
library(raster)
library(DESeq2)
library(ggnewscale)
library(ggrepel)


set.seed(123)

getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#summary function
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length(na.omit(xx[[col]])),
                     mean = mean   (na.omit(xx[[col]])),
                     sd   = sd     (na.omit(xx[[col]]))
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#extract: median IQR
format_summary = function(x) {
  output = paste(round(x[1],0)," (",round(x[2],0),", ",round(x[3],0),")",sep="")
  return(output)
}

#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
#NOW START DEG LOOP moving through cell types
#----------------------------------------------------------------------------#
#----------------------------------------------------------------------------#
#Move to my machine
# scp mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/ad_pv_project/suerat_obj/ad_pv_merged_v2.RDS /Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/data/scRNAseq/seurat_object

suerat_obj = readRDS("/Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/data/scRNAseq/seurat_object/ad_pv_merged_v2.RDS")

#Make master dir
master_dir = "/c4/home/mtaylor4/derm/ad_pv_project/results/deg_rashx/"

#Make deg output dir
deg_dir = paste(master_dir,
                sep="")
dir.create(deg_dir)

#Make text file dir
deg_textfile_fir = paste(deg_dir,
                         "ind_files/",
                         sep="")
dir.create(deg_textfile_fir)

#Make sample id - dis fiedl
suerat_obj@meta.data$sample_id_dis = paste(suerat_obj@meta.data$sample_id, suerat_obj@meta.data$dis)
unique(suerat_obj@meta.data$sample_id_dis)


#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#Write a function to recaptulate Yale's genes
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
yale_signature_method = function(input_cds, results_name, normalization_method, proportion_samples){
  
  # input_cds=suerat_obj
  # results_name="all_paper_samples"
  # normalization_method="LogNormalize"
  # proportion_samples=0.8
  
  #Get samples in a cell type subset
  all_sample_names = unique(input_cds@meta.data$sample_id_dis)
  hc_sample_names = all_sample_names[grep("Healthy",all_sample_names)]
  pv_sample_names = c( all_sample_names[grep("Psoriasis Vulgaris",all_sample_names)])
  ad_sample_names = c( all_sample_names[grep("Atopic Dermatitis",all_sample_names)] )
  
  input_cds <- NormalizeData(input_cds, normalization.method = normalization_method)
  input_cds <- FindVariableFeatures(input_cds, selection.method = "vst", nfeatures = 10000)
  
  library(presto)
  #------------------------------------------------------------------------#
  #Yale Step 1: PV DEGs 
  #------------------------------------------------------------------------#
  deg_df_pv <- FindMarkers(input_cds,
                           ident.1 = row.names(input_cds@meta.data[input_cds@meta.data$sample_id_dis %in% pv_sample_names,]),
                           ident.2 = row.names(input_cds@meta.data[input_cds@meta.data$sample_id_dis %in% hc_sample_names,])
  )
  
  deg_df_pv = deg_df_pv[order(-deg_df_pv$avg_log2FC),] #Re-order from high to low FC
  deg_df_pv = cbind(row.names(deg_df_pv),deg_df_pv) #Add col of gene symbols
  names(deg_df_pv)[1] = "gene"
  
  #------------------------------------------------------------------------#
  #Yale Step 2: AD DEGs 
  #------------------------------------------------------------------------#
  deg_df_ad <- FindMarkers(input_cds,
                           ident.1 = row.names(input_cds@meta.data[input_cds@meta.data$sample_id_dis %in% ad_sample_names,]),
                           ident.2 = row.names(input_cds@meta.data[input_cds@meta.data$sample_id_dis %in% hc_sample_names,])
  )
  
  deg_df_ad = deg_df_ad[order(-deg_df_ad$avg_log2FC),] #Re-order from high to low FC
  deg_df_ad = cbind(row.names(deg_df_ad),deg_df_ad) #Add col of gene symbols
  names(deg_df_ad)[1] = "gene"
  
  #------------------------------------------------------------------------#
  #Yale Step 3: AD v PV DEGs 
  #------------------------------------------------------------------------#
  #Only use genes that are sig for a disease vs normal
  genes_to_test = c(subset(deg_df_ad, avg_log2FC > 0.25 & p_val_adj < 0.05)$gene,
                    subset(deg_df_pv, avg_log2FC > 0.25 & p_val_adj < 0.05)$gene)
  deg_df_adpv <- FindMarkers(input_cds,
                             ident.1 = row.names(input_cds@meta.data[input_cds@meta.data$sample_id_dis %in% ad_sample_names,]),
                             ident.2 = row.names(input_cds@meta.data[input_cds@meta.data$sample_id_dis %in% pv_sample_names,]),
                             genes.use = genes_to_test
  )
  
  deg_df_adpv = deg_df_adpv[order(-deg_df_adpv$avg_log2FC),] #Re-order from high to low FC
  deg_df_adpv = cbind(row.names(deg_df_adpv),deg_df_adpv) #Add col of gene symbols
  names(deg_df_adpv)[1] = "gene"
  head(deg_df_adpv)
  
  #Now get AD and PV specific genes from this AD v PV DEGs
  ad_V_pv_genes = subset(deg_df_adpv, avg_log2FC > 0.25 & p_val_adj < 0.05)
  pv_V_ad_genes = subset(deg_df_adpv, avg_log2FC < -0.25 & p_val_adj < 0.05)
  
  #Now filter AD-v-HC and PV-v-HC by AD-PV and PV-AD respectively
  ad_specific_global = deg_df_ad[deg_df_ad$gene %in% ad_V_pv_genes$gene, ]
  pv_specific_global = deg_df_pv[deg_df_pv$gene %in% pv_V_ad_genes$gene, ]
  
  global_dis_genes = c(ad_specific_global$gene,pv_specific_global$gene)
  
  #--------------------------------------------------------------------------------------------#
  #Yale Step 4: #Get sample specific indvidiual disease vs grouped HC DEGs
  #--------------------------------------------------------------------------------------------#
  
  all_sample_names = unique(input_cds@meta.data$sample_id_dis)
  hc_sample_names = all_sample_names[grep("Healthy control",all_sample_names)]
  pv_sample_names = c( all_sample_names[grep("Psoriasis Vulgaris",all_sample_names)])
  ad_sample_names = c( all_sample_names[grep("Atopic Dermatitis",all_sample_names)] )
  dis_sample_names = c(pv_sample_names,ad_sample_names)
  
  #Loop through each individual sample and compare to grouped normals
  ind_dis_deg = do.call(rbind, lapply(1:length(dis_sample_names), function(y) {
    
    #y=1
    print(paste("y=",y,sep=""))
    #------------------------------------------------------------------------#
    #Individual disease vs grouped normals
    #------------------------------------------------------------------------#
    
    deg_df <- FindMarkers(input_cds,
                          ident.1 = row.names(input_cds@meta.data[input_cds@meta.data$sample_id_dis %in% dis_sample_names[y],]),
                          ident.2 = row.names(input_cds@meta.data[input_cds@meta.data$sample_id_dis %in% hc_sample_names,]),
                          genes.use=global_dis_genes
    )
    
    deg_df = deg_df[order(-deg_df$avg_log2FC),] #Re-order from high to low FC
    deg_df = cbind(dis_sample_names[y],row.names(deg_df),deg_df) #Add col of gene symbols
    names(deg_df)[1:2] = c("sample","gene")
    deg_df
    
  })) #End indivdiual sample DEG loop
  
  #---------------------------------------------------------#
  ##Yale Step 5: #Retain only individual DEGs that appear in at least 80% of dis samples
  #---------------------------------------------------------#
  
  #Individual sample DEG filters
  ind_dis_deg = subset(ind_dis_deg, p_val_adj<=0.05 & avg_log2FC>0.25) 
  
  #Get genes that occur at least in 80% of samples
  min_samples = round(length(ad_sample_names)*proportion_samples,0)
  ind_ad_deg = ind_dis_deg[ind_dis_deg$sample %in% ad_sample_names,]
  ad_gene_n = data.frame(table(ind_ad_deg$gene))
  ad_gene_n = ad_gene_n[order(-ad_gene_n$Freq),]
  ad_gene_n = subset(ad_gene_n, Freq>=min_samples)
  
  min_samples = round(length(pv_sample_names)*proportion_samples,0)
  ind_pv_deg = ind_dis_deg[ind_dis_deg$sample %in% pv_sample_names,]
  pv_gene_n = data.frame(table(ind_pv_deg$gene))
  pv_gene_n = pv_gene_n[order(-pv_gene_n$Freq),]
  pv_gene_n = subset(pv_gene_n, Freq>=min_samples)
  
  #---------------------------------------------------------#
  ##Yale Step 6: #Retest AD v PV using only individual DEGs that appear in at least 80% of samples
  #---------------------------------------------------------#
  genes_to_test = c(ad_gene_n$Var1,pv_gene_n$Var1)
  deg_df_adpv <- FindMarkers(input_cds,
                             ident.1 = row.names(input_cds@meta.data[input_cds@meta.data$sample_id_dis %in% ad_sample_names,]),
                             ident.2 = row.names(input_cds@meta.data[input_cds@meta.data$sample_id_dis %in% pv_sample_names,]),
                             genes.use = genes_to_test
  )
  
  deg_df_adpv = deg_df_adpv[order(-deg_df_adpv$avg_log2FC),] #Re-order from high to low FC
  deg_df_adpv = cbind(row.names(deg_df_adpv),deg_df_adpv) #Add col of gene symbols
  names(deg_df_adpv)[1] = "gene"
  head(deg_df_adpv)
  
  #Now get AD and PV specific genes from this AD v PV DEGs
  ad_V_pv_genes = subset(deg_df_adpv, avg_log2FC > 0.25 & p_val_adj < 0.05)
  pv_V_ad_genes = subset(deg_df_adpv, avg_log2FC < -0.25 & p_val_adj < 0.05)
  pv_V_ad_genes = pv_V_ad_genes[order(-pv_V_ad_genes$avg_log2FC),]  #re-oder PV from low to high FC bc NEGATIVE logFC is HIGHER IN PV
  
  #---------------------------------------------------------#
  ##Yale Step 7: #Final
  #---------------------------------------------------------#    
  deg_df_pv_consensus = pv_V_ad_genes
  deg_df_ad_consensus = ad_V_pv_genes
  
  deg_df_pv_consensus$dis = "PV"
  deg_df_ad_consensus$dis = "AD"
  
  deg_df_consensus = rbind(deg_df_pv_consensus,deg_df_ad_consensus)
  
  OUT <- createWorkbook()
  
  #Write out to excel workbook
  sheet_name = cell_types[p]
  addWorksheet(OUT, sheet_name ) # Add some sheets to the workbook
  writeData(OUT, sheet = sheet_name, x = deg_df_consensus) #Write data to a shee
  
  # Export the Excel file
  saveWorkbook(OUT,
               paste("/c4/home/mtaylor4/derm/ad_pv_project/results/yale_recapiulation/",results_name,"_",p,"_","TRM1_remapped",
                     "_",normalization_method,"_",
                     "remapped_DEGs_yale_method_v2.xlsx",sep=""),
               overwrite=T)
  
  
  #------------------------------------------------------------#
  #Compare these new Yale-derived DEGs to Yale's original DEGs
  #------------------------------------------------------------#
  
  #Read in Yale's original gene signatures
  yale_signatures = read.xlsx(
    "/c4/home/mtaylor4/derm/ad_pv_project/data/yale_signature_genes/cluster_2_expanded_degs.xlsx",
    sheet=1)
  length(intersect(deg_df_consensus$gene,yale_signatures$gene)) / nrow(deg_df_consensus)
  length(intersect(deg_df_consensus$gene,yale_signatures$gene)) / nrow(yale_signatures)
  
  
  #Compare to yale's original lists
  ad_genes_yale = c("TWIST1","LGALS1",    "IL32",   "CAPG",  "ITM2C", "MFHAS1", "ANXA1","SOS1", "CSGALNACT1","LMO4",  "IFITM2","S100A10",  "MT-ND5",  "CYSLTR1","PLA2G16", "SYNE2", "THADA",    "NEAT1","IL17RB", "RPL36A","ARHGAP21",    "NBAS",  "ACTG1","PRKX", "TGFBR3",   "TIMP1","TNFSF10", "AHNAK",    "MT-ND2",  "ISG15",  "RPL17",  "LONRF2",   "CD99","TSHZ2", "MMP25",   "IFITM1","MT-ND1",  "BIRC3",  "FAM102A", "LPCAT2","NRIP3", "CRIP1",  "CLU",   "PLP2", "ZFP36",  "ZFP36L2",  "TUBA1B","GATA3","SLC5A3",    "SFXN1", "FANK1","TAGLN2") 
  pv_genes_yale = c("CXCL13", "MTRNR2L12",  "CD7",   "MGAT4A","FTH1", "LAYN", "IL17F", "KLRB1", "GNLY","CPM", "CTSH", "GBP5", "SOX4", "CLEC2B",    "GZMB", "CD2",   "CEBPD","ODF2L","LAG3", "LRRN3", "ARHGEF12",  "PTPN13",   "TNFAIP3", "TRPS1", "SNX9", "METRNL",  "BTG1", "JUN","SPOCK2",   "GABARAPL1",  "PMEPA1", "HIST1H1E","RBPJ", "LINC01871", "MAP3K4","H1FX", "UBC",  "GALNT1",  "PNRC1","GABPB1-AS1", "RPS26", "MUC20-OT1",  "CHN1",  "NAP1L4",   "PTMS",  "F2R",   "CTLA4", "DAPK2","RAP1B","CCR6", "B3GALT2", "YPEL2",  "FYN",   "PPDPF","SLA2", "CBLB", "ADGRG1","SARAF") 
  
  length(intersect(deg_df_ad_consensus$gene,ad_genes_yale)) / length(ad_genes_yale)
  length(intersect(deg_df_pv_consensus$gene,pv_genes_yale)) / length(ad_genes_yale)
  
  #--------------------------------------------------------------------#
  #RashX plotting
  #--------------------------------------------------------------------#
  
  #reduce to unique genes that wil be plotted
  n_genes_signature = 49
  ad_genes = deg_df_ad_consensus$gene[1:n_genes_signature]
  pv_genes = deg_df_pv_consensus$gene[1:n_genes_signature]
  
  #--------------------------------------------------------------------#
  #AD gene data_matrix
  ad_genes = na.omit(row.names(input_cds)[match(ad_genes,row.names(input_cds))])
  ad_mat = t(input_cds@assays$RNA@data[ad_genes,])
  ad_mat[1:10, 1:10]
  ad_int = data.frame(
    input_cds$dis,
    input_cds$sample_id_renamed,
    rowSums(ad_mat)
  )
  names(ad_int) = c("dis", "sample", "gene_sig")
  names(ad_int) = c("dis", "sample", "gene_sig")
  ad_dfc = summarySE(ad_int,
                     measurevar = 'gene_sig',
                     groupvars = c("dis", "sample"))
  names(ad_dfc)[c(4, 6)] = c("ad_gene_sig", "ad_se")
  
  #--------------------------------------------------------------------#
  #PV gene data_matrix
  pv_genes = na.omit(row.names(input_cds)[match(pv_genes,row.names(input_cds))])
  pv_mat = t(input_cds@assays$RNA@data[pv_genes,])
  pv_mat[1:10, 1:10]
  pv_int = data.frame(
    input_cds$dis,
    input_cds$sample_id_renamed,
    rowSums(pv_mat)
  )
  names(pv_int) = c("dis", "sample", "gene_sig")
  names(pv_int) = c("dis", "sample", "gene_sig")
  pv_dfc = summarySE(pv_int,
                     measurevar = 'gene_sig',
                     groupvars = c("dis", "sample"))
  head(pv_dfc)
  names(pv_dfc)[c(4, 6)] = c("pv_gene_sig", "pv_se")
  
  #--------------------------------------------------------------------#
  re_int = cbind(ad_dfc[, c(1, 2, 4, 6)], pv_dfc[, c(4, 6)]) #combine the gene signatures
  
  #----------------------------------------------------------------#
  ### Plot: Sample-level means without individual cells
  #----------------------------------------------------------------#
  re_int$dis = gsub("Healthy control","HC",re_int$dis)
  re_int$dis = gsub("Atopic Dermatitis","AD",re_int$dis)
  re_int$dis = gsub("Psoriasis Vulgaris","PV",re_int$dis)
  
  re_int_ad = subset(re_int, dis == "AD")
  re_int_hc = subset(re_int, dis == "HC")
  re_int_pv = subset(re_int, dis == "PV")
  
  p1.5 <- ggplot(data = re_int, aes(x = ad_gene_sig, y = pv_gene_sig)) +
    stat_density_2d(data=re_int_ad, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="burlywood4") +
    new_scale_fill() +   # same as `new_scale("color")` 
    stat_density_2d(data=re_int_hc, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="darkblue") +
    new_scale_fill() +   # same as `new_scale("color")` 
    stat_density_2d(data=re_int_pv, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="olivedrab4") +
    geom_point(alpha = 1, aes(color = dis, size=dis)) +
    geom_text_repel(
      data = re_int,
      aes(x = ad_gene_sig, y = pv_gene_sig, label = sample),
      box.padding = 0.4
    ) + #add labels
    geom_errorbarh(aes(
      xmax = ad_gene_sig + ad_se,
      xmin = ad_gene_sig - ad_se,
      color = dis
    )) +
    geom_errorbar(aes(
      ymax = pv_gene_sig + pv_se,
      ymin = pv_gene_sig - pv_se,
      color = dis
    )) +
    ggtitle(paste("Trm1","Remapped genes Yale's method",normalization_method)) +
    theme_classic() +
    xlab("AD-specific genes") +
    ylab("PV-specific genes") +
    scale_color_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "Atypical Eczema" = "darkred",
        "AD" = "burlywood4",
        "PV" = "olivedrab4",
        "Hand dermatitis" = "purple"
      )
    ) +
    scale_size_manual(
      name = "Disease",
      values = c(
        "HC" = 2,
        "PRP" = 4,
        "AD" = 2,
        "PV" = 2
      )
    ) +
    theme(legend.position = "none")
  
  
  #------------------------------------------------------------------------#
  #------------------------------------------------------------------------#
  #------------------------------------------------------------------------#
  #Orginal RashX genes
  
  ad_genes = c("TWIST1","LGALS1",    "IL32",   "CAPG",  "ITM2C", "MFHAS1", "ANXA1","SOS1", "CSGALNACT1","LMO4",  "IFITM2","S100A10",  "MT-ND5",  "CYSLTR1","PLA2G16", "SYNE2", "THADA",    "NEAT1","IL17RB", "RPL36A","ARHGAP21",    "NBAS",  "ACTG1","PRKX", "TGFBR3",   "TIMP1","TNFSF10", "AHNAK",    "MT-ND2",  "ISG15",  "RPL17",  "LONRF2",   "CD99","TSHZ2", "MMP25",   "IFITM1","MT-ND1",  "BIRC3",  "FAM102A", "LPCAT2","NRIP3", "CRIP1",  "CLU",   "PLP2", "ZFP36",  "ZFP36L2",  "TUBA1B","GATA3","SLC5A3",    "SFXN1", "FANK1","TAGLN2") 
  pv_genes = c("CXCL13", "MTRNR2L12",  "CD7",   "MGAT4A","FTH1", "LAYN", "IL17F", "KLRB1", "GNLY","CPM", "CTSH", "GBP5", "SOX4", "CLEC2B",    "GZMB", "CD2",   "CEBPD","ODF2L","LAG3", "LRRN3", "ARHGEF12",  "PTPN13",   "TNFAIP3", "TRPS1", "SNX9", "METRNL",  "BTG1", "JUN","SPOCK2",   "GABARAPL1",  "PMEPA1", "HIST1H1E","RBPJ", "LINC01871", "MAP3K4","H1FX", "UBC",  "GALNT1",  "PNRC1","GABPB1-AS1", "RPS26", "MUC20-OT1",  "CHN1",  "NAP1L4",   "PTMS",  "F2R",   "CTLA4", "DAPK2","RAP1B","CCR6", "B3GALT2", "YPEL2",  "FYN",   "PPDPF","SLA2", "CBLB", "ADGRG1","SARAF") 
  
  ad_genes = ad_genes[1:n_genes_signature]
  pv_genes = pv_genes[1:n_genes_signature]
  
  input_cds@assays$RNA[1:10,1:10]
  
  input_cds@assays$RNA@data
  #--------------------------------------------------------------------#
  #Get the genes of interest
  #--------------------------------------------------------------------#
  #AD gene data_matrix
  ad_genes = na.omit(row.names(input_cds)[match(ad_genes,row.names(input_cds))])
  ad_mat = t(input_cds@assays$RNA@data[ad_genes,])
  ad_mat[1:10, 1:10]
  ad_int = data.frame(
    input_cds$dis,
    input_cds$sample_id_renamed,
    rowSums(ad_mat)
  )
  names(ad_int) = c("dis", "sample", "gene_sig")
  names(ad_int) = c("dis", "sample", "gene_sig")
  ad_dfc = summarySE(ad_int,
                     measurevar = 'gene_sig',
                     groupvars = c("dis", "sample"))
  names(ad_dfc)[c(4, 6)] = c("ad_gene_sig", "ad_se")
  
  #--------------------------------------------------------------------#
  #PV gene data_matrix
  pv_genes = na.omit(row.names(input_cds)[match(pv_genes,row.names(input_cds))])
  pv_mat = t(input_cds@assays$RNA@data[pv_genes,])
  pv_mat[1:10, 1:10]
  pv_int = data.frame(
    input_cds$dis,
    input_cds$sample_id_renamed,
    rowSums(pv_mat)
  )
  names(pv_int) = c("dis", "sample", "gene_sig")
  names(pv_int) = c("dis", "sample", "gene_sig")
  pv_dfc = summarySE(pv_int,
                     measurevar = 'gene_sig',
                     groupvars = c("dis", "sample"))
  head(pv_dfc)
  names(pv_dfc)[c(4, 6)] = c("pv_gene_sig", "pv_se")
  
  #--------------------------------------------------------------------#
  re_int = cbind(ad_dfc[, c(1, 2, 4, 6)], pv_dfc[, c(4, 6)]) #combine the gene signatures
  
  #----------------------------------------------------------------#
  ### Plot: Sample-level means without individual cells
  #----------------------------------------------------------------#
  re_int$dis = unique(re_int$dis)
  re_int$dis = gsub("Healthy control","HC",re_int$dis)
  re_int$dis = gsub("Atopic Dermatitis","AD",re_int$dis)
  re_int$dis = gsub("Psoriasis Vulgaris","PV",re_int$dis)
  
  re_int_ad = subset(re_int, dis == "AD")
  re_int_hc = subset(re_int, dis == "HC")
  re_int_pv = subset(re_int, dis == "PV")
  
  p1 <- ggplot(data = re_int, aes(x = ad_gene_sig, y = pv_gene_sig)) +
    stat_density_2d(data=re_int_ad, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="burlywood4") +
    new_scale_fill() +   # same as `new_scale("color")` 
    stat_density_2d(data=re_int_hc, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="darkblue") +
    new_scale_fill() +   # same as `new_scale("color")` 
    stat_density_2d(data=re_int_pv, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="olivedrab4") +
    #geom_point(alpha = 1, aes(color = dis, size=dis)) +
    geom_text_repel(
      data = re_int,
      aes(x = ad_gene_sig, y = pv_gene_sig, label = sample),
      box.padding = 0.4
    ) + #add labels
    geom_errorbarh(aes(
      xmax = ad_gene_sig + ad_se,
      xmin = ad_gene_sig - ad_se,
      color = dis
    )) +
    geom_errorbar(aes(
      ymax = pv_gene_sig + pv_se,
      ymin = pv_gene_sig - pv_se,
      color = dis
    )) +
    ggtitle(paste("Trm1","JID RashX genes",normalization_method)) +
    theme_classic() +
    xlab("AD-specific genes") +
    ylab("PV-specific genes") +
    scale_color_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "Atypical Eczema" = "darkred",
        "AD" = "burlywood4",
        "PV" = "olivedrab4",
        "Hand dermatitis" = "purple"
      )
    ) +
    scale_size_manual(
      name = "Disease",
      values = c(
        "HC" = 2,
        "PRP" = 4,
        "AD" = 2,
        "PV" = 2
      )
    ) +
    theme(legend.position = "none")
  
  
  #Combine both 
  layout <- cowplot::plot_grid(
    p1, p1.5, ncol = 2, nrow = 1)
  
  png(paste("/c4/home/mtaylor4/derm/ad_pv_project/results/yale_recapiulation/",results_name,"_",i,"_",cell_types[i],
            "_",normalization_method,"_",
            "hyperD_comparison_plots.png",sep=""),
      width=16,height=8,units="in",res=400)
  print(layout)
  dev.off()
  
}


#Now implement Yale's method  
yale_signature_method(suerat_obj_sub_orginal_samples,"SI_paper_samples","LogNormalize",0.8)
yale_signature_method(suerat_obj_sub,"AD-PV_paper_samples","LogNormalize",0.8)
yale_signature_method(suerat_obj_sub,"AD-PV_paper_samples","CLR",0.8)
yale_signature_method(suerat_obj_sub,"AD-PV_paper_samples","RC",0.8)

#Transfer results to my machine
scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/ad_pv_project/results/yale_recapiulation/ /Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/results


#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#Pseudo-bulk and replot https://satijalab.org/seurat/articles/de_vignette
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#


pseudobulk_analysis_fx = function(input_cds, results_name, normalization_method) {
  
  #input_cds=suerat_obj
  #results_name="AD_PV_paper_samples"
  #normalization_method="CLR"
  
  #Normalize
  # input_cds <- NormalizeData(input_cds, normalization.method = normalization_method)
  
  #Run on normalized counts to produce scaled pseudobulk data for RashX plotting
  pseudobulk_dat <- AggregateExpression(input_cds, assays = "RNA",
                                        return.seurat = T, group.by = c("dis", "sample_id_renamed"))
  
  sample_metadat_df = t(data.frame(strsplit(row.names(pseudobulk_dat@meta.data),"_")))
  pseudobulk_dat@meta.data$dis = sample_metadat_df[,1]
  pseudobulk_dat@meta.data$sample_id_renamed = sample_metadat_df[,2]
  
  
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #DEGs from pseudobulked data
  #Must re-run pseudobulking on un-normalized counts bc DESeq2 requires integers
  
  pseudobulk_dat_deg <- AggregateExpression(input_cds, assays = "RNA", slot="counts",
                                            return.seurat = T, group.by = c("dis", "sample_id_renamed"),
                                            normalization.method=NULL)
  
  sample_metadat_df = t(data.frame(strsplit(row.names(pseudobulk_dat_deg@meta.data),"_")))
  pseudobulk_dat_deg@meta.data$dis = sample_metadat_df[,1]
  pseudobulk_dat_deg@meta.data$sample_id_renamed = sample_metadat_df[,2]
  
  #Check that there is pseudobulked integer data
  pseudobulk_dat_deg@assays$RNA@data[1:10,1:10]
  pseudobulk_dat_deg@assays$RNA$counts[1:10,1:10]
  
  
  #--------------------------------------------------------#
  #--------------------------------------------------------#
  #--------------------------------------------------------#
  
  #violin plots of Ray's genes
  
  to_plot = data.frame(pseudobulk_dat_deg@meta.data$dis,
             pseudobulk_dat_deg@meta.data$sample_id_renamed,
             t(pseudobulk_dat_deg@assays$RNA@data[c("PTGER3","RIMS2","IL1RL1",
                                                  "TWIST1","NBAS","IL17RB"),]))
  names(to_plot)[1:2] = c("dis","sample")
  dfm = melt(to_plot, id=c("dis","sample")) #Melt so can be plotted
  dfm = dfm[dfm$dis %in% c("Atopic Dermatitis","Healthy control","Psoriasis Vulgaris"),]
  dfm$dis <- gsub("Atopic Dermatitis","AD",dfm$dis)
  dfm$dis <- gsub("Healthy control","PV",dfm$dis)
  dfm$dis <- gsub("Psoriasis Vulgaris","HC",dfm$dis)
  
  dis_figs =  ggplot(dfm, aes(x=dis,y=value, fill=dis, color=dis)) +
    facet_wrap(variable ~ ., nrow=2) +
    #stat_compare_means(comparisons = my_comparisons, method="wilcox.test")  +
    geom_violin(width=1, alpha=0.4) +
    geom_jitter(width=0.25, size=0.5, alpha=0.9) +
    geom_boxplot(outlier.size=0, width=0.3, color="black", size=0.7, alpha=0.4) +
    theme_classic() +
    ggtitle("Pseudobulked Trm1 expression") + 
    xlab("") + ylab("expression") + 
    theme(axis.text.x=element_text(size=20,color='black'),
          axis.text.y=element_text(size=15,color='black'),
          title=element_text(size=25),
          axis.ticks.x = element_line(size=0),
          strip.text = element_text(size = 20),
          legend.position="none") +
    scale_color_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "AD" = "burlywood4",
        "PV" = "olivedrab4"
      )
    ) +
    scale_fill_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "AD" = "burlywood4",
        "PV" = "olivedrab4"
      )
    )
  
  
  png("/Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/results/pseudobulk/ray_genes.png",
      width=12,height=12,units="in",res=400)
  print(dis_figs)
  dev.off()
  
  
  
  #--------------------------------------------------------#
  
  
  #--------------------------------------------------------#
  #Function to get DEGs for factor levels, setting factor as IDENT on the Seurat object
  #eg to get dis factor level DEGs, you do Idents(input_pseudobulk_dat) <- "dis" before calling the function
  #--------------------------------------------------------#
  pseudobulk_deg_fx = function(input_pseudobulk_dat,numerator_factor_level, denominator_factor_level) {
    
    #input_pseudobulk_dat=pseudobulk_dat_deg
    #numerator_factor_level="Atopic Dermatitis"
    #denominator_factor_level="Psoriasis Vulgaris"
    
    pseudo_bulk_deg_adpv <- FindMarkers(object = input_pseudobulk_dat, 
                                        ident.1 = numerator_factor_level, 
                                        ident.2 = denominator_factor_level,
                                        assay="RNA",
                                        slot="counts",
                                        test.use = "DESeq2",
                                        min.cells.group = 1)
    
    deg_df_adpv = pseudo_bulk_deg_adpv[order(-pseudo_bulk_deg_adpv$avg_log2FC),] #Re-order from high to low FC
    deg_df_adpv = cbind(row.names(deg_df_adpv),deg_df_adpv) #Add col of gene symbols
    names(deg_df_adpv)[1] = "gene"
    
    return(deg_df_adpv)
  }
  
  
  fold_change_cutoff = 0.5
  
  #--------------------------------------------------------------#
  #AD v PV
  #-----------------#
  Idents(pseudobulk_dat_deg) <- "dis"
  unique(pseudobulk_dat_deg@meta.data$dis)
  deg_dfadpv = pseudobulk_deg_fx(pseudobulk_dat_deg,"Atopic Dermatitis","Psoriasis Vulgaris")
  ad_V_pv_genes = subset(deg_dfadpv, avg_log2FC > fold_change_cutoff & p_val_adj < 0.05)
  pv_V_ad_genes = subset(deg_dfadpv, avg_log2FC < -fold_change_cutoff & p_val_adj < 0.05)
  pv_V_ad_genes = pv_V_ad_genes[order(-pv_V_ad_genes$avg_log2FC),]  #re-oder PV from low to high FC bc NEGATIVE logFC is HIGHER IN PV
  
  #Write out ad v pv group
  write.table(deg_dfadpv,
              file="/Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/results/pseudobulk/results/DEG/ad_v_pv_full_res_updated_samples.txt",
              row.names=F, quote=F, sep="\t")
  
  
  #--------------------------------------------------------------#
  #AD/PV v HC
  #-----------------#
  #AD v HC
  deg_dfadhc = pseudobulk_deg_fx(pseudobulk_dat_deg,"Atopic Dermatitis","Healthy control")
  ad_V_hc_genes = subset(deg_dfadhc, avg_log2FC > fold_change_cutoff & p_val_adj < 0.05)
  
  deg_dfadhc[deg_dfadhc$gene %in% c("TWIST1","NBAS","IL17RB"), ]
  
  #-----------------#
  #PV v HC
  deg_dfpvhc = pseudobulk_deg_fx(pseudobulk_dat_deg,"Psoriasis Vulgaris","Healthy control")
  pv_V_hc_genes = subset(deg_dfpvhc, avg_log2FC > fold_change_cutoff & p_val_adj < 0.05)
  
  #-----------------#
  ad_only <- ad_V_hc_genes$gene[!(ad_V_hc_genes$gene %in% pv_V_hc_genes$gene)]
  pv_only <- pv_V_hc_genes$gene[!(pv_V_hc_genes$gene %in% ad_V_hc_genes$gene)]
  
  #--------------------------------------------------------------#
  
  #--------------------------------------------------------------#
  #Now filter out the ad_only genes that ALSO occur in AD v PV
  #same thing for pv-only genes that ALSO occur in PV v AD
  #-----------------#
  
  ad_only_ad = ad_only[ad_only %in% ad_V_pv_genes$gene]
  pv_genes_pv = pv_only[pv_only %in% pv_V_ad_genes$gene]
  
  ad_genes = ad_only_ad
  pv_genes = pv_genes_pv
  
  data.frame(ad_genes)
  data.frame(pv_genes)
  
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #BATCH-CORRECDT PSEUDOBULK SAMPLES : RashX on pseudobulked samples with pseudobulk DEGs
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  bulk_dat = t(pseudobulk_dat@assays$RNA@data)
  bulk_dat[1:10,1:10]
  
  row.names(bulk_dat) = as.vector(data.frame(strsplit(row.names(bulk_dat),"_"))[2,])
  
  #Get batch effect vectors
  sample_name_vec = row.names(bulk_dat)
  batch_effect_vec = rep("nml_cell_n", nrow(bulk_dat))
  batch_effect_vec[match("HD2",sample_name_vec)] = "low_cell_n" #replace batch effect vec value w new batch call
  batch_effect_vec[match("At2",sample_name_vec)] = "low_cell_n" #replace batch effect vec value w new batch call
  batch_effect_vec[match("AD17",sample_name_vec)] = "low_cell_n" #replace batch effect vec value w new batch call
  batch_effect_vec[match("At9",sample_name_vec)] = "low_cell_n" #replace batch effect vec value w new batch call

  batch_effect_vec[match("AD6",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("AD11",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("AD12",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("AD13",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("AD14",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("AD15",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("AD16",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("AD17",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("At6a",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("At7",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("At8",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("At10",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("At11",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("At18",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("At19",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("At20",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("HC14",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("PV16",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("PV17",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("PV18",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("PV19",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("PV20",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("PV21",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("PV22",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  batch_effect_vec[match("HD2",sample_name_vec)] = "5_prime" #replace batch effect vec value w new batch call
  
  
  # Perform ComBat batch correction
  library(sva)
  combat_data <- ComBat(t(bulk_dat), batch = batch_effect_vec, par.prior = TRUE, prior.plots = FALSE)# Create a grouping variable indicating the study for each sample
  combat_data[1:10,1:5]
  
  #--------------------------------------------------------------------#
  #AD gene data_matrix
  ad_mat = t(combat_data[ad_genes,])
  ad_mat[1:10,1:5]
  ad_int = data.frame(
    row.names(ad_mat),
    row.names(ad_mat),
    rowSums(ad_mat)
  )
  names(ad_int) = c("dis", "sample", "gene_sig")
  
  #--------------------------------------------------------------------#
  #PV gene data_matrix
  pv_mat = t(combat_data[pv_genes,])
  pv_mat[1:10,1:5]
  pv_int = data.frame(
    row.names(pv_mat),
    row.names(pv_mat),
    rowSums(pv_mat)
  )
  names(pv_int) = c("dis", "sample", "gene_sig")
  
  #--------------------------------------------------------------------#
  #combine the gene signatures into single df
  re_int = cbind(ad_int[, c(1, 2, 3)], pv_int[, c(3)])
  names(re_int)[c(3,4)] = c("ad_gene_sig","pv_gene_sig")
  
  #----------------------------------------------------------------#
  ### Plot: Sample-level means without individual cells
  #----------------------------------------------------------------#
  re_int$dis[grep("AD",re_int$dis)] <- "AD"
  re_int$dis[grep("PV",re_int$dis)] <- "PV"
  re_int$dis[grep("HC",re_int$dis)] <- "HC"
  re_int$dis[grep("At",re_int$dis)] <- "At"
  
  re_int_ad = subset(re_int, dis == "AD")
  re_int_hc = subset(re_int, dis == "HC")
  re_int_pv = subset(re_int, dis == "PV")
  
  unique(re_int)
  
  #----------------------------------------------------------------#
  ### Plot: Sample-level means without individual cells
  #----------------------------------------------------------------#
  
  p1.9 <- ggplot(data = re_int, aes(x = ad_gene_sig, y = pv_gene_sig)) +
    stat_density_2d(data=re_int_ad, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="burlywood4") +
    new_scale_fill() +   # same as `new_scale("color")` 
    stat_density_2d(data=re_int_hc, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="darkblue") +
    new_scale_fill() +   # same as `new_scale("color")` 
    stat_density_2d(data=re_int_pv, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="olivedrab4") +
    geom_point(alpha = 1, aes(color = dis)) +
    geom_text_repel(
      data = re_int,
      aes(x = ad_gene_sig, y = pv_gene_sig, label = sample),
      box.padding = 0.4
    ) + #add labels
    ggtitle(paste("Trm1","Pseudobulked samples DESeq2 top AD-PV genes, batch-corrected")) +
    theme_classic() +
    xlab("AD-specific genes") +
    ylab("PV-specific genes") +
    scale_color_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "At" = "darkred",
        "AD" = "burlywood4",
        "PV" = "olivedrab4"
        
      )
    ) +
    theme(legend.position = "none")
  #p1.9
  
  
  layout <- cowplot::plot_grid(
    p1, p1.9, ncol = 2, nrow = 1)
  
  png(paste("/c4/home/mtaylor4/derm/ad_pv_project/results/pseudobulk/results/",results_name,"_",p,"_",cell_types[p],
            "_",normalization_method,"_",
            "_batch_hyperD_comparison_plots.png",sep=""),
      width=16,height=8,units="in",res=400)
  print(layout)
  dev.off()
  
  #--------------------------------------------------------------------#
  #Calculate dis centroid of pseudobulked samples
  
  re_int_centroid = re_int[re_int$dis %in% c("AD","PV","HC","At"),]
  
  ad_dfc = summarySE(re_int_centroid,
                     measurevar = 'ad_gene_sig',
                     groupvars = c("dis"))
  names(ad_dfc)[c(3, 5)] = c("ad_gene_sig", "ad_se")
  
  pv_dfc = summarySE(re_int_centroid,
                     measurevar = 'pv_gene_sig',
                     groupvars = c("dis"))
  names(pv_dfc)[c(3, 5)] = c("pv_gene_sig", "pv_se")
  
  re_int_centroid_2 = cbind(ad_dfc[, c(1, 3, 5)], pv_dfc[, c(3, 5)]) #combine the gene signatures

  p1.5.2 <- ggplot(data = re_int_centroid_2, aes(x = ad_gene_sig, y = pv_gene_sig)) +
    geom_point(alpha = 1, aes(color = dis),size=9) +
    geom_text_repel(
      data = re_int_centroid_2,
      aes(x = ad_gene_sig, y = pv_gene_sig, label = dis),
      box.padding = 4, size=10
    ) + #add labels
    geom_errorbarh(aes(
      xmax = ad_gene_sig + ad_se,
      xmin = ad_gene_sig - ad_se,
      color = dis
    ),height=10) +
    geom_errorbar(aes(
      ymax = pv_gene_sig + pv_se,
      ymin = pv_gene_sig - pv_se,
      color = dis
    ),width=3) +
    ggtitle(paste("Trm1","Pseudobulked samples DESeq2, batch-corrected, disease centroids")) +
    theme_classic() +
    xlab("AD-specific genes") +
    ylab("PV-specific genes") +
    scale_color_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "Atypical Eczema" = "darkred",
        "AD" = "burlywood4",
        "PV" = "olivedrab4",
        "Hand dermatitis" = "purple"
      )
    ) +
    theme(legend.position = "none")
  
  
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #Save atypical scores
  head(re_int)
  table(re_int$dis)
  atypic_df = subset(re_int, dis=="At")
  

  #Loop through atypicals and calculate distance from HC, PV, and AD centroids
  ad_centroid = c(subset(re_int_centroid_2, dis=="AD")$ad_gene_sig,
                  subset(re_int_centroid_2, dis=="AD")$pv_gene_sig)
  pv_centroid = c(subset(re_int_centroid_2, dis=="PV")$ad_gene_sig,
                  subset(re_int_centroid_2, dis=="PV")$pv_gene_sig)
  hc_centroid = c(subset(re_int_centroid_2, dis=="HC")$ad_gene_sig,
                  subset(re_int_centroid_2, dis=="HC")$pv_gene_sig)
  
  atypic_distances = do.call(rbind, lapply(1:nrow(atypic_df), function(i){
    
    #i=1
    ad_dist = pointDistance(p1=c(atypic_df$ad_gene_sig[i],atypic_df$pv_gene_sig[i]),
                                   p2=ad_centroid, lonlat=F) #calculate atypical sample vs centroid dist
    pv_dist = pointDistance(p1=c(atypic_df$ad_gene_sig[i],atypic_df$pv_gene_sig[i]),
                            p2=pv_centroid, lonlat=F) #calculate atypical sample vs centroid dist
    hc_dist = pointDistance(p1=c(atypic_df$ad_gene_sig[i],atypic_df$pv_gene_sig[i]),
                            p2=hc_centroid, lonlat=F) #calculate atypical sample vs centroid dist
    
    data.frame(atypic_df[i,],ad_dist,pv_dist,hc_dist)
    
  }))
  
  write.table(atypic_distances,
              file="/Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/results/pseudobulk/results/atypical_analysis/atypical_distances.txt",
              sep="\t",row.names = F, quote=F)
  
  
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  
  #--------------------------------------------------------------------#
  #RashX plot on Yale's features optimized in the JID paper
  #--------------------------------------------------------------------#
  #Compare to yale's original lists
  ad_genes_yale = c("TWIST1","LGALS1",    "IL32",   "CAPG",  "ITM2C", "MFHAS1", "ANXA1","SOS1", "CSGALNACT1","LMO4",  "IFITM2","S100A10",  "MT-ND5",  "CYSLTR1","PLA2G16", "SYNE2", "THADA",    "NEAT1","IL17RB", "RPL36A","ARHGAP21",    "NBAS",  "ACTG1","PRKX", "TGFBR3",   "TIMP1","TNFSF10", "AHNAK",    "MT-ND2",  "ISG15",  "RPL17",  "LONRF2",   "CD99","TSHZ2", "MMP25",   "IFITM1","MT-ND1",  "BIRC3",  "FAM102A", "LPCAT2","NRIP3", "CRIP1",  "CLU",   "PLP2", "ZFP36",  "ZFP36L2",  "TUBA1B","GATA3","SLC5A3",    "SFXN1", "FANK1","TAGLN2") 
  pv_genes_yale = c("CXCL13", "MTRNR2L12",  "CD7",   "MGAT4A","FTH1", "LAYN", "IL17F", "KLRB1", "GNLY","CPM", "CTSH", "GBP5", "SOX4", "CLEC2B",    "GZMB", "CD2",   "CEBPD","ODF2L","LAG3", "LRRN3", "ARHGEF12",  "PTPN13",   "TNFAIP3", "TRPS1", "SNX9", "METRNL",  "BTG1", "JUN","SPOCK2",   "GABARAPL1",  "PMEPA1", "HIST1H1E","RBPJ", "LINC01871", "MAP3K4","H1FX", "UBC",  "GALNT1",  "PNRC1","GABPB1-AS1", "RPS26", "MUC20-OT1",  "CHN1",  "NAP1L4",   "PTMS",  "F2R",   "CTLA4", "DAPK2","RAP1B","CCR6", "B3GALT2", "YPEL2",  "FYN",   "PPDPF","SLA2", "CBLB", "ADGRG1","SARAF") 
  
  ad_genes = ad_genes_yale
  pv_genes = pv_genes_yale
  
  
  #--------------------------------------------------------------------#
  #AD gene data_matrix
  ad_genes = ad_genes[ad_genes %in% row.names(combat_data)]
  ad_mat = t(combat_data[ad_genes,])
  ad_mat[1:10,1:5]
  ad_int = data.frame(
    row.names(ad_mat),
    row.names(ad_mat),
    rowSums(ad_mat)
  )
  names(ad_int) = c("dis", "sample", "gene_sig")
  
  #--------------------------------------------------------------------#
  #PV gene data_matrix
  pv_genes = pv_genes[pv_genes %in% row.names(combat_data)]
  pv_mat = t(combat_data[pv_genes,])
  pv_mat[1:10,1:5]
  pv_int = data.frame(
    row.names(pv_mat),
    row.names(pv_mat),
    rowSums(pv_mat)
  )
  names(pv_int) = c("dis", "sample", "gene_sig")
  
  #--------------------------------------------------------------------#
  #combine the gene signatures into single df
  re_int = cbind(ad_int[, c(1, 2, 3)], pv_int[, c(3)])
  names(re_int)[c(3,4)] = c("ad_gene_sig","pv_gene_sig")
  
  #----------------------------------------------------------------#
  ### Plot: Sample-level means without individual cells
  #----------------------------------------------------------------#
  re_int$dis[grep("AD",re_int$dis)] <- "AD"
  re_int$dis[grep("PV",re_int$dis)] <- "PV"
  re_int$dis[grep("HC",re_int$dis)] <- "HC"
  re_int$dis[grep("At",re_int$dis)] <- "At"
  
  re_int_ad = subset(re_int, dis == "AD")
  re_int_hc = subset(re_int, dis == "HC")
  re_int_pv = subset(re_int, dis == "PV")
  
  unique(re_int)
  
  #----------------------------------------------------------------#
  ### Plot: Sample-level means without individual cells
  #----------------------------------------------------------------#
  
  p1 <- ggplot(data = re_int, aes(x = ad_gene_sig, y = pv_gene_sig)) +
    stat_density_2d(data=re_int_ad, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="burlywood4") +
    new_scale_fill() +   # same as `new_scale("color")` 
    stat_density_2d(data=re_int_hc, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="darkblue") +
    new_scale_fill() +   # same as `new_scale("color")` 
    stat_density_2d(data=re_int_pv, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="olivedrab4") +
    geom_point(alpha = 1, aes(color = dis)) +
    geom_text_repel(
      data = re_int,
      aes(x = ad_gene_sig, y = pv_gene_sig, label = sample),
      box.padding = 0.4
    ) + #add labels
    ggtitle(paste("Trm1","JID AD-PV gene signatures, batch-corrected")) +
    theme_classic() +
    xlab("AD-specific genes") +
    ylab("PV-specific genes") +
    scale_color_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "At" = "darkred",
        "AD" = "burlywood4",
        "PV" = "olivedrab4"
        
      )
    ) +
    theme(legend.position = "none")
  #p1.9
  
  
  #--------------------------------------------------------------------#
  #Calculate dis centroid of pseudobulked samples
  
  re_int_centroid = re_int[re_int$dis %in% c("AD","PV","HC","At"),]
  
  ad_dfc = summarySE(re_int_centroid,
                     measurevar = 'ad_gene_sig',
                     groupvars = c("dis"))
  names(ad_dfc)[c(3, 5)] = c("ad_gene_sig", "ad_se")
  
  pv_dfc = summarySE(re_int_centroid,
                     measurevar = 'pv_gene_sig',
                     groupvars = c("dis"))
  names(pv_dfc)[c(3, 5)] = c("pv_gene_sig", "pv_se")
  
  re_int_centroid_2 = cbind(ad_dfc[, c(1, 3, 5)], pv_dfc[, c(3, 5)]) #combine the gene signatures
  
  p1.1.2 <- ggplot(data = re_int_centroid_2, aes(x = ad_gene_sig, y = pv_gene_sig)) +
    geom_point(alpha = 1, aes(color = dis),size=9) +
    geom_text_repel(
      data = re_int_centroid_2,
      aes(x = ad_gene_sig, y = pv_gene_sig, label = dis),
      box.padding = 4, size=10
    ) + #add labels
    geom_errorbarh(aes(
      xmax = ad_gene_sig + ad_se,
      xmin = ad_gene_sig - ad_se,
      color = dis
    ),height=10) +
    geom_errorbar(aes(
      ymax = pv_gene_sig + pv_se,
      ymin = pv_gene_sig - pv_se,
      color = dis
    ),width=3) +
    ggtitle(paste("Trm1","JID AD-PV gene signatures, batch-corrected, disease centroids")) +
    theme_classic() +
    xlab("AD-specific genes") +
    ylab("PV-specific genes") +
    scale_color_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "Atypical Eczema" = "darkred",
        "AD" = "burlywood4",
        "PV" = "olivedrab4",
        "Hand dermatitis" = "purple"
      )
    ) +
    theme(legend.position = "none")
  
  
  
  #--------------------------------------------------------------------#
  #Save individual samples as a svg: JID signature + new bulk signature
  #--------------------------------------------------------------------#
  
  layout <- cowplot::plot_grid(
    p1, p1.9, ncol = 2, nrow = 1)
  
  ggsave(file=paste("/Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/results/pseudobulk/ad_pv_updated_samples.svg",sep=""), 
         plot=layout, 
         width=16,height=8)
  
  
  #--------------------------------------------------------------------#
  #Save DISEASE CENTROIDS as a svg: JID signature + new bulk signature
  #--------------------------------------------------------------------#
  
  layout <- cowplot::plot_grid(
    p1.1.2, p1.5.2, ncol = 2, nrow = 1)
  
  ggsave(file=paste("/Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/results/pseudobulk/ad_pv_updated_samples_disease_centroids.svg",sep=""), 
         plot=layout, 
         width=10,height=5)
  
  
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#

  
  
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  
  #Lopo through and show individual genes' performance: taking the top PSO gene and the top AD genes, and moving on
  #AD gene data_matrix
  min_gene_set = min(length(ad_genes), length(pv_genes))
  
  ind_gene_covs = do.call(rbind, lapply(1:min_gene_set, function(h) {
    
    #h=2
    
    ad_mat = data.frame(combat_data[ad_genes[h],])
    ad_int = data.frame(
      row.names(ad_mat),
      row.names(ad_mat),
      ad_mat[,1]
    )
    names(ad_int) = c("dis", "sample", "gene_sig")
    
    #--------------------------------------------------------------------#
    #PV gene data_matrix
    pv_mat = data.frame(combat_data[pv_genes[h],])
    pv_int = data.frame(
      row.names(pv_mat),
      row.names(pv_mat),
      pv_mat[,1]
    )
    names(pv_int) = c("dis", "sample", "gene_sig")
    
    #--------------------------------------------------------------------#
    #combine the gene signatures into single df
    re_int = cbind(ad_int[, c(1, 2, 3)], pv_int[, c(3)])
    names(re_int)[c(3,4)] = c("ad_gene_sig","pv_gene_sig")
    
    #----------------------------------------------------------------#
    ### CENTER DATA
    center = function(x) {x - mean(x)}
    #----------------------------------------------------------------#
    re_int$ad_gene_sig = center(re_int$ad_gene_sig)
    re_int$pv_gene_sig = center(re_int$pv_gene_sig)
    
    
    #----------------------------------------------------------------#
    ### Plot: Sample-level means without individual cells
    #----------------------------------------------------------------#
    re_int$dis[grep("AD",re_int$dis)] <- "AD"
    re_int$dis[grep("PV",re_int$dis)] <- "PV"
    re_int$dis[grep("HC",re_int$dis)] <- "HC"
    
    re_int_ad = subset(re_int, dis == "AD")
    re_int_hc = subset(re_int, dis == "HC")
    re_int_pv = subset(re_int, dis == "PV")
    
    #----------------------------------------------------------------#
    ### Plot: Sample-level means without individual cells
    #----------------------------------------------------------------#
    
    p2 <- ggplot(data = re_int, aes(x = ad_gene_sig, y = pv_gene_sig)) +
      stat_ellipse(geom="polygon" ,alpha = 0.2, aes(fill=dis)) +
      geom_point(alpha = 1, aes(color = dis, size=dis)) +
      ggtitle(paste(cell_types[p],"Pseudobulked DESeq2 individual gene rank:",h)) +
      theme_classic() +
      xlab(paste("AD-specific gene: ",ad_genes[h],sep="")) +
      ylab(paste("PV-specific gene: ",pv_genes[h],sep="")) +
      scale_color_manual(
        name = "Disease",
        values = c(
          "HC" = "darkblue",
          "PRP" = "darkred",
          "AD" = "burlywood4",
          "PV" = "olivedrab4"
        )
      ) +
      scale_fill_manual(
        name = "Disease",
        values = c(
          "HC" = "darkblue",
          "PRP" = "darkred",
          "AD" = "burlywood4",
          "PV" = "olivedrab4"
          
        )
      ) +
      scale_size_manual(
        name = "Disease",
        values = c(
          "HC" = 2,
          "PRP" = 4,
          "AD" = 2,
          "PV" = 2
        )
      ) +
      theme(legend.position = "none")
    #p2
    
    png(paste("/c4/home/mtaylor4/derm/ad_pv_project/results/individual_gene_var/pseudobulk_trm1/deg_rank_",
              h,".png",sep=""),
        width=6,height=6,units="in",res=400)
    print(p2)
    dev.off()
    
    #Calculate coefficients of variation
    
    ad_cov <- sd(re_int_ad$ad_gene_sig) / mean(re_int_ad$ad_gene_sig) * 100
    pv_cov <- sd(re_int_pv$pv_gene_sig) / mean(re_int_pv$pv_gene_sig) * 100
    return( data.frame(h, ad_genes[h],pv_genes[h],ad_cov,pv_cov) )
    
  }))  #End individual gene cov loop
  
  
  scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/ad_pv_project/results/individual_gene_var/pseudobulk_trm1 /Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/results/pseudobulk/results/individual_gene_variation
  
  
  #----------------------------------------------------------------------------------------#
  #----------------------------------------------------------------------------------------#
  #----------------------------------------------------------------------------------------#
  #----------------------------------------------------------------------------------------#
  #----------------------------------------------------------------------------------------#
  #LDA on AD and PV genes
  
  #--------------------------------------------------------------------#
  #Get matrix for top AD-PV genes
  adpv_mat = t(combat_data[c(ad_genes,pv_genes),])
  adpv_mat[1:10,1:5]
  int_mat = data.frame(
    row.names(adpv_mat),
    row.names(adpv_mat),
    adpv_mat
  )
  names(int_mat)[1:2] = c("dis", "sample")
  
  int_mat$dis[grep("AD",int_mat$dis)] <- "AD"
  int_mat$dis[grep("PV",int_mat$dis)] <- "PV"
  int_mat$dis[grep("HC",int_mat$dis)] <- "HC"
  
  genes_to_consider = colnames(int_mat)[3:ncol(int_mat)]
  int_mat$sample = as.factor(int_mat$sample)
  formula_str = paste('dis~',paste(genes_to_consider,collapse = '+'),sep = '')
  linear <- lda(as.formula(formula_str), int_mat)
  
  set.seed(123)
  lda_cmb_tissuetype = ggord(linear, int_mat$dis,txt = NULL,vectyp = 0, ellipse=F, hull=T,
                             size=2) +
    theme_classic() +
    ggtitle("AD-PV gene LDA") + 
    theme(axis.text.x=element_text(size=15,color='black'),
          axis.text.y=element_text(size=15,color='black'),
          title=element_text(size=20),
          axis.ticks.x = element_line(size=0),
          strip.text = element_text(size = 20),
          legend.position="bottom",
          legend.text=element_text(size=15)) +
    scale_color_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "PRP" = "darkred",
        "AD" = "burlywood4",
        "PV" = "olivedrab4"
      )
    ) +
    scale_fill_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "PRP" = "darkred",
        "AD" = "burlywood4",
        "PV" = "olivedrab4"
      )
    ) +
    xlim(-10,15) + ylim(-6,6) + #nsk only
    xlab("LD1") + ylab("LD2")
  
  png(paste("/c4/home/mtaylor4/derm/ad_pv_project/results/pseudobulk/results/lda/figs/trm1_lda.png" ,sep=""),
      width=8,height=8,units="in",res=400)
  print(lda_cmb_tissuetype)
  dev.off()
  
  
  lda_cmb_loadings = ggord(linear, int_mat$dis,vectyp = 1, ellipse=F, hull=T,
                           size=2, alpha=0.1) +
    theme_classic() +
    ggtitle("AD-PV gene LDA loadings") + 
    theme(axis.text.x=element_text(size=15,color='black'),
          axis.text.y=element_text(size=15,color='black'),
          title=element_text(size=20),
          axis.ticks.x = element_line(size=0),
          strip.text = element_text(size = 20),
          legend.position="bottom",
          legend.text=element_text(size=15)) +
    scale_color_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "PRP" = "darkred",
        "AD" = "burlywood4",
        "PV" = "olivedrab4"
      )
    ) +
    scale_fill_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "PRP" = "darkred",
        "AD" = "burlywood4",
        "PV" = "olivedrab4"
      )
    ) 
  
  
  png(paste("/c4/home/mtaylor4/derm/ad_pv_project/results/pseudobulk/results/lda/figs/trm1_lda_loadings.png" ,sep=""),
      width=8,height=8,units="in",res=400)
  print(lda_cmb_loadings)
  dev.off()
  
  #plot top and bttm loadings
  loadings = linear$scaling
  n_loadings = 10 #number of gene loadings to plot 
  
  i=1
  print(i)
  pc_ordered = data.frame(loadings[order(-loadings[,i]),])
  pc_ordered$gene = row.names(pc_ordered)
  pc_ordered$gene = factor(pc_ordered$gene,levels=pc_ordered$gene) #re-order genes
  
  #plot top gene loadings
  top_pc = head(pc_ordered,n_loadings)
  top_pc = top_pc[order(top_pc[,i]),]
  top_pc$gene = factor(top_pc$gene,levels=top_pc$gene) #re-order genes
  
  pos_load_fig_pc1_cmb = ggplot() +
    geom_bar(data=top_pc, aes(x=top_pc[,i], y=gene), stat="identity",color="white",fill="darkolivegreen") +
    theme_classic() +
    ggtitle(paste("LD",i," POSITIVE loadings", sep="")) +
    xlab("coefficient") +
    ylab("") +
    theme(legend.position="none") +
    theme(axis.text.x=element_text(size=12,color='black'),
          axis.text.y=element_text(size=17,color='black'),
          axis.title=element_text(size=12,color='black'),
          plot.title=element_text(size=20,color='black'),
          legend.position = "none")
  
  #plot bttm gene loadings
  bttm_pc = tail(pc_ordered,n_loadings)
  
  neg_load_fig_pc1_cmb = ggplot() +
    geom_bar(data=bttm_pc, aes(x=-top_pc[,i], y=gene), stat="identity",color="white",fill="darkblue") +
    theme_classic() +
    ggtitle(paste("LD",i," NEGATIVE loadings", sep="")) +
    xlab("coefficient") +
    ylab("") +
    theme(legend.position="none") +
    theme(axis.text.x=element_text(size=12,color='black'),
          axis.text.y=element_text(size=17,color='black'),
          axis.title=element_text(size=12,color='black'),
          plot.title=element_text(size=20,color='black'),
          legend.position = "none")
  
  i=2
  print(i)
  pc_ordered = data.frame(loadings[order(-loadings[,i]),])
  pc_ordered$gene = row.names(pc_ordered)
  pc_ordered$gene = factor(pc_ordered$gene,levels=pc_ordered$gene) #re-order genes
  
  #plot top gene loadings
  top_pc = head(pc_ordered,n_loadings)
  top_pc = top_pc[order(top_pc[,i]),]
  top_pc$gene = factor(top_pc$gene,levels=top_pc$gene) #re-order genes
  
  pos_load_fig_pc2_cmb = ggplot() +
    geom_bar(data=top_pc, aes(x=top_pc[,i], y=gene), stat="identity",color="white",fill="tan") +
    theme_classic() +
    ggtitle(paste("LD",i," POSITIVE loadings", sep="")) +
    xlab("coefficient") +
    ylab("") +
    theme(legend.position="none") +
    theme(axis.text.x=element_text(size=12,color='black'),
          axis.text.y=element_text(size=17,color='black'),
          axis.title=element_text(size=12,color='black'),
          plot.title=element_text(size=20,color='black'),
          legend.position = "none")
  
  #plot bttm gene loadings
  bttm_pc = tail(pc_ordered,n_loadings)
  
  neg_load_fig_pc2_cmb = ggplot() +
    geom_bar(data=bttm_pc, aes(x=-top_pc[,i], y=gene), stat="identity",color="white",fill="darkblue") +
    theme_classic() +
    ggtitle(paste("LD",i," NEGATIVE loadings", sep="")) +
    xlab("coefficient") +
    ylab("") +
    theme(legend.position="none") +
    theme(axis.text.x=element_text(size=12,color='black'),
          axis.text.y=element_text(size=17,color='black'),
          axis.title=element_text(size=12,color='black'),
          plot.title=element_text(size=20,color='black'),
          legend.position = "none")
  
  
  layout <- cowplot::plot_grid(
    pos_load_fig_pc1_cmb, neg_load_fig_pc1_cmb,
    pos_load_fig_pc2_cmb, neg_load_fig_pc2_cmb, ncol = 2, nrow = 2)
  
  png(paste("/c4/home/mtaylor4/derm/ad_pv_project/results/pseudobulk/results/lda/figs/","lda_barchart_loadings.png" ,sep=""),
      width=16,height=16,units="in",res=400)
  print(layout)
  dev.off()
  
  scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/ad_pv_project/results/pseudobulk/results/lda/figs /Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/results/pseudobulk/results/LDA 
  
  #----------------------------------------------------------------------------------------#
  #----------------------------------------------------------------------------------------#
  #----------------------------------------------------------------------------------------#
  #----------------------------------------------------------------------------------------#
  #----------------------------------------------------------------------------------------#
  #Model selection on AD and PV genes: logistic regression with stepwise selection
  
  int_mat[1:5,1:5]
  glm_dat = int_mat[,-2] #get rid of sample column
  glm_dat[1:5,1:5]
  
  # Fit initial logistic regression model with all predictors
  glm_dat$dis = as.factor(glm_dat$dis)
  initial_model <- glm(dis ~ ., data = glm_dat, family = binomial)
  
  # Perform stepwise model selection using AIC
  final_model <- stepAIC(initial_model, direction = "both", trace = FALSE)
  
  # Print final model summary
  summary(final_model)
  
  #--------------------------------------------------------------------#
  #individual gene RashX model
  
  ad_mat = data.frame(combat_data["SMKR1",])
  ad_int = data.frame(
    row.names(ad_mat),
    row.names(ad_mat),
    ad_mat[,1]
  )
  names(ad_int) = c("dis", "sample", "gene_sig")
  
  #--------------------------------------------------------------------#
  #PV gene data_matrix
  pv_mat = data.frame(combat_data["LINC01147",])
  pv_int = data.frame(
    row.names(pv_mat),
    row.names(pv_mat),
    pv_mat[,1]
  )
  names(pv_int) = c("dis", "sample", "gene_sig")
  
  #--------------------------------------------------------------------#
  #combine the gene signatures into single df
  re_int = cbind(ad_int[, c(1, 2, 3)], pv_int[, c(3)])
  names(re_int)[c(3,4)] = c("ad_gene_sig","pv_gene_sig")
  
  
  #----------------------------------------------------------------#
  ### Plot: Sample-level means without individual cells
  #----------------------------------------------------------------#
  re_int$dis[grep("AD",re_int$dis)] <- "AD"
  re_int$dis[grep("PV",re_int$dis)] <- "PV"
  re_int$dis[grep("HC",re_int$dis)] <- "HC"
  
  re_int_ad = subset(re_int, dis == "AD")
  re_int_hc = subset(re_int, dis == "HC")
  re_int_pv = subset(re_int, dis == "PV")
  
  p3 <- ggplot(data = re_int, aes(x = ad_gene_sig, y = pv_gene_sig)) +
    stat_ellipse(geom="polygon" ,alpha = 0.2, aes(fill=dis)) +
    geom_point(alpha = 1, aes(color = dis, size=dis)) +
    ggtitle("Stepwise logistic model selection optimal genes") +
    theme_classic() +
    xlab(paste("AD-specific gene: ","SMKR1",sep="")) +
    ylab(paste("PV-specific gene: ","LINC01147",sep="")) +
    scale_color_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "PRP" = "darkred",
        "AD" = "burlywood4",
        "PV" = "olivedrab4"
      )
    ) +
    scale_fill_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "PRP" = "darkred",
        "AD" = "burlywood4",
        "PV" = "olivedrab4"
        
      )
    ) +
    scale_size_manual(
      name = "Disease",
      values = c(
        "HC" = 2,
        "PRP" = 4,
        "AD" = 2,
        "PV" = 2
      )
    ) +
    theme(legend.position = "none")
  
  
  #----------------------------------------------------------------------------------------#
  #----------------------------------------------------------------------------------------#
  #----------------------------------------------------------------------------------------#
  #LASSO model selection
  
  int_mat[1:5,1:5]
  glm_dat = int_mat[,-2] #get rid of sample column
  glm_dat[1:5,1:5]
  
  # Create matrix of predictors
  X <- as.matrix(glm_dat[, -which(names(glm_dat) == "dis")])
  
  # Convert response variable to a factor
  glm_dat$dis <- as.factor(glm_dat$dis)
  
  # Fit LASSO logistic regression model
  library(glmnet)
  
  lasso_model <- cv.glmnet(X, as.factor(glm_dat$dis), family = "multinomial", alpha = 1, type.measure = "class")
  
  # Extract the best lambda value
  best_lambda <- lasso_model$lambda.min
  
  # Get AD and PV coefficients for the best model
  best_model = coef(lasso_model, s = best_lambda)
  ad_model = data.frame(row.names(best_model$AD),best_model$AD)[-1,] #remove intercept
  ad_model = ad_model[ad_model[,2]>0,]
  
  pv_model = data.frame(row.names(best_model$PV),best_model$PV)[-1,] #remove intercept
  pv_model = pv_model[pv_model[,2]>0,]
  
  ad_genes_lasso = ad_model[,1]
  pv_genes_lasso = pv_model[,1]
  
  #--------------------------------------------------------------------#
  #AD gene data_matrix
  ad_mat = t(combat_data[ad_genes_lasso,])
  ad_int = data.frame(
    row.names(ad_mat),
    row.names(ad_mat),
    rowSums(ad_mat)
  )
  names(ad_int) = c("dis", "sample", "gene_sig")
  
  #--------------------------------------------------------------------#
  #PV gene data_matrix
  pv_mat = t(combat_data[pv_genes_lasso,])
  pv_int = data.frame(
    row.names(pv_mat),
    row.names(pv_mat),
    rowSums(pv_mat)
  )
  names(pv_int) = c("dis", "sample", "gene_sig")
  
  #--------------------------------------------------------------------#
  #combine the gene signatures into single df
  re_int = cbind(ad_int[, c(1, 2, 3)], pv_int[, c(3)])
  names(re_int)[c(3,4)] = c("ad_gene_sig","pv_gene_sig")
  
  #----------------------------------------------------------------#
  ### Plot: Sample-level means without individual cells
  #----------------------------------------------------------------#
  re_int$dis[grep("AD",re_int$dis)] <- "AD"
  re_int$dis[grep("PV",re_int$dis)] <- "PV"
  re_int$dis[grep("HC",re_int$dis)] <- "HC"
  
  re_int_ad = subset(re_int, dis == "AD")
  re_int_hc = subset(re_int, dis == "HC")
  re_int_pv = subset(re_int, dis == "PV")
  
  #----------------------------------------------------------------#
  ### Plot: Sample-level means without individual cells
  #----------------------------------------------------------------#
  
  p2.1 <- ggplot(data = re_int, aes(x = ad_gene_sig, y = pv_gene_sig)) +
    stat_density_2d(data=re_int_ad, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="burlywood4") +
    new_scale_fill() +   # same as `new_scale("color")` 
    stat_density_2d(data=re_int_hc, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="darkblue") +
    new_scale_fill() +   # same as `new_scale("color")` 
    stat_density_2d(data=re_int_pv, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="olivedrab4") +
    geom_point(alpha = 1, aes(color = dis, size=dis)) +
    geom_text_repel(
      data = re_int,
      aes(x = ad_gene_sig, y = pv_gene_sig, label = sample),
      box.padding = 0.4
    ) + #add labels
    ggtitle(paste(cell_types[p],"LASSO model selection")) +
    theme_classic() +
    xlab("AD-specific genes") +
    ylab("PV-specific genes") +
    scale_color_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "PRP" = "darkred",
        "AD" = "burlywood4",
        "PV" = "olivedrab4"
        
      )
    ) +
    scale_size_manual(
      name = "Disease",
      values = c(
        "HC" = 2,
        "PRP" = 4,
        "AD" = 2,
        "PV" = 2
      )
    ) +
    theme(legend.position = "none")
  #p2.1
  
  png(paste("/c4/home/mtaylor4/derm/ad_pv_project/results/pseudobulk/results/lasso/",results_name,"_",p,"_",cell_types[p],
            "_",normalization_method,"_",
            "_LASSO_hyperD_comparison_plots.png",sep=""),
      width=8,height=8,units="in",res=400)
  print(p2.1)
  dev.off()
  
  scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/ad_pv_project/results/pseudobulk/results/lasso/ /Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/results/pseudobulk/results/ 
    
    
    #----------------------------------------------------------------------------------------#
    #----------------------------------------------------------------------------------------#
    #----------------------------------------------------------------------------------------#
    #GBM model selection
    
    library(caret)
  library(gbm)
  library(randomForest)
  
  # Set up the training control
  control <- trainControl(method = "cv", number = 5)
  
  # Run GBRF model
  model <- train(dis ~ ., data = glm_dat, method = "rf", trControl = control)
  
  # Get importance of predictors
  importance_measures <- varImp(model)
  imporatnce_df = data.frame ( row.names(importance_measures$importance) , importance_measures$importance)
  names(imporatnce_df)[1] = "predictor_gene"
  imporatnce_df = imporatnce_df[order(-imporatnce_df$Overall),] # order from high to low importance
  
  imporatnce_df$rank = seq(1:nrow(imporatnce_df))
  
  #Examine curve of importance and find inflection point of gene number to pass to the next argument to graph: n_genes_from_gbm
  ggplot(data = imporatnce_df, aes(x = rank, y = Overall)) +
    geom_point() +
    geom_line() +
    ggtitle(paste(cell_types[p],"GBRF model selection")) +
    theme_classic() +
    xlab("Rank") +
    ylab("Importance")
  
  
  #--------------------------------------------------------------------#
  n_genes_from_gbm = 5
  #Get genes that intersect with AD and PV pseudbulk consensus genes defined above
  ad_genes_gbm = intersect(imporatnce_df$predictor_gene[1:n_genes_from_gbm],ad_genes)
  pv_genes_gbm = intersect(imporatnce_df$predictor_gene[1:n_genes_from_gbm],pv_genes)
  
  #--------------------------------------------------------------------#
  #AD gene data_matrix
  ad_mat = t(combat_data[ad_genes_gbm,])
  ad_int = data.frame(
    row.names(ad_mat),
    row.names(ad_mat),
    rowSums(ad_mat)
  )
  names(ad_int) = c("dis", "sample", "gene_sig")
  
  #--------------------------------------------------------------------#
  #PV gene data_matrix
  pv_mat = t(combat_data[pv_genes_gbm,])
  pv_int = data.frame(
    row.names(pv_mat),
    row.names(pv_mat),
    rowSums(pv_mat)
  )
  names(pv_int) = c("dis", "sample", "gene_sig")
  
  #--------------------------------------------------------------------#
  #combine the gene signatures into single df
  re_int = cbind(ad_int[, c(1, 2, 3)], pv_int[, c(3)])
  names(re_int)[c(3,4)] = c("ad_gene_sig","pv_gene_sig")
  
  #----------------------------------------------------------------#
  ### Plot: Sample-level means without individual cells
  #----------------------------------------------------------------#
  re_int$dis[grep("AD",re_int$dis)] <- "AD"
  re_int$dis[grep("PV",re_int$dis)] <- "PV"
  re_int$dis[grep("HC",re_int$dis)] <- "HC"
  
  re_int_ad = subset(re_int, dis == "AD")
  re_int_hc = subset(re_int, dis == "HC")
  re_int_pv = subset(re_int, dis == "PV")
  
  #----------------------------------------------------------------#
  ### Plot: Sample-level means without individual cells
  #----------------------------------------------------------------#
  
  p2.2 <- ggplot(data = re_int, aes(x = ad_gene_sig, y = pv_gene_sig)) +
    stat_density_2d(data=re_int_ad, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="burlywood4") +
    new_scale_fill() +   # same as `new_scale("color")` 
    stat_density_2d(data=re_int_hc, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="darkblue") +
    new_scale_fill() +   # same as `new_scale("color")` 
    stat_density_2d(data=re_int_pv, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="olivedrab4") +
    geom_point(alpha = 1, aes(color = dis, size=dis)) +
    geom_text_repel(
      data = re_int,
      aes(x = ad_gene_sig, y = pv_gene_sig, label = sample),
      box.padding = 0.4
    ) + #add labels
    ggtitle(paste(cell_types[p]," GBM model selection, top ",n_genes_from_gbm," genes",sep="")) +
    theme_classic() +
    xlab("AD-specific genes") +
    ylab("PV-specific genes") +
    scale_color_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "PRP" = "darkred",
        "AD" = "burlywood4",
        "PV" = "olivedrab4"
        
      )
    ) +
    scale_size_manual(
      name = "Disease",
      values = c(
        "HC" = 2,
        "PRP" = 4,
        "AD" = 2,
        "PV" = 2
      )
    ) +
    theme(legend.position = "none")
  #p2.2
  
  png(paste("/c4/home/mtaylor4/derm/ad_pv_project/results/pseudobulk/results/lasso/",results_name,"_",p,"_",cell_types[p],
            "_",normalization_method,"_",
            "_GBM_hyperD_comparison_plots.png",sep=""),
      width=8,height=8,units="in",res=400)
  print(p2.1)
  dev.off()
  
  scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/ad_pv_project/results/pseudobulk/results/lasso/ /Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/results/pseudobulk/results/ 
    
    
    #----------------------------------------------------------------------------------------#
    #----------------------------------------------------------------------------------------#
    #----------------------------------------------------------------------------------------#
    #----------------------------------------------------------------------------------------#
    #----------------------------------------------------------------------------------------#
    #----------------------------------------------------------------------------------------#
    #----------------------------------------------------------------------------------------#
    #----------------------------------------------------------------------------------------#
    #----------------------------------------------------------------------------------------#
  #----------------------------------------------------------------------------------------#
  #----------------------------------------------------------------------------------------#
  #----------------------------------------------------------------------------------------#
  #----------------------------------------------------------------------------------------#
  #----------------------------------------------------------------------------------------#
  #----------------------------------------------------------------------------------------#
  
  #Heatmap on batch-corrected data
  while (!is.null(dev.list()))  dev.off()
  
  #Scale matrix
  result_mat = t(as.matrix(combat_data[c(ad_genes,pv_genes),] ))
  result_mat = scale(result_mat)
  
  # Pairwise correlation between rows (samples)
  rows.cor <- cor(t(result_mat), use = "pairwise.complete.obs", method = "spearman")
  
  ## Row- and column-wise clustering using correlation 
  hclust.row <- hclust(as.dist(1-rows.cor))
  
  library(gplots)
  
  png(paste("/c4/home/mtaylor4/derm/ad_pv_project/results/pseudobulk/results/heatmaps/pseudobulk_",
            p,"_",cell_types[p],".png",sep=""),
      width=36,height=24,units="in",res=200)
  print(
    
    heatmap.2(result_mat, col = bluered(100), 
              trace = "none", density.info = "none",
              Colv = T,
              labCol=as.expression(lapply(colnames(result_mat), function(a) bquote(italic(.(a))))),
              Rowv = as.dendrogram(hclust.row),
              lhei=c(1, 10),
              cexCol = 1.4,
              cexRow = 2.5,
              cex.main = 7,
              dendrogram="row",
              main = paste("Pseudobulked ",cell_types[p],sep=""),
              margins=c(12,9),
              offsetCol=0.0,
              offsetRow=0.2
    )
  ) # end print call
  dev.off()
  
  #scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/ad_pv_project/results/pseudobulk/results/heatmaps/ /Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/results/pseudobulk/results
  
  
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  
  #--------------------------------------------------------------------#
  #Single-cell RashX from pseudobulk DEGs
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  #AD gene data_matrix
  ad_genes = na.omit(row.names(input_cds)[match(ad_genes,row.names(input_cds))])
  ad_mat = t(input_cds@assays$SCT@data[ad_genes,])
  ad_mat[1:10, 1:10]
  ad_int = data.frame(
    input_cds$dis,
    input_cds$sample_id_renamed,
    rowSums(ad_mat)
  )
  names(ad_int) = c("dis", "sample", "gene_sig")
  names(ad_int) = c("dis", "sample", "gene_sig")
  ad_dfc = summarySE(ad_int,
                     measurevar = 'gene_sig',
                     groupvars = c("dis", "sample"))
  names(ad_dfc)[c(4, 6)] = c("ad_gene_sig", "ad_se")
  
  #--------------------------------------------------------------------#
  #PV gene data_matrix
  pv_genes = na.omit(row.names(input_cds)[match(pv_genes,row.names(input_cds))])
  pv_mat = t(input_cds@assays$SCT@data[pv_genes,])
  pv_mat[1:10, 1:10]
  pv_int = data.frame(
    input_cds$dis,
    input_cds$sample_id_renamed,
    rowSums(pv_mat)
  )
  names(pv_int) = c("dis", "sample", "gene_sig")
  names(pv_int) = c("dis", "sample", "gene_sig")
  pv_dfc = summarySE(pv_int,
                     measurevar = 'gene_sig',
                     groupvars = c("dis", "sample"))
  head(pv_dfc)
  names(pv_dfc)[c(4, 6)] = c("pv_gene_sig", "pv_se")
  
  #--------------------------------------------------------------------#
  re_int = cbind(ad_dfc[, c(1, 2, 4, 6)], pv_dfc[, c(4, 6)]) #combine the gene signatures
  
  #----------------------------------------------------------------#
  ### Plot: Sample-level means without individual cells
  #----------------------------------------------------------------#
  re_int$dis = gsub("Healthy control","HC",re_int$dis)
  re_int$dis = gsub("Atopic Dermatitis","AD",re_int$dis)
  re_int$dis = gsub("Psoriasis Vulgaris","PV",re_int$dis)
  
  re_int_ad = subset(re_int, dis == "AD")
  re_int_hc = subset(re_int, dis == "HC")
  re_int_pv = subset(re_int, dis == "PV")
  
  
  p2 <- ggplot(data = re_int, aes(x = ad_gene_sig, y = pv_gene_sig)) +
    stat_density_2d(data=re_int_ad, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="burlywood4") +
    new_scale_fill() +   # same as `new_scale("color")` 
    stat_density_2d(data=re_int_hc, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="darkblue") +
    new_scale_fill() +   # same as `new_scale("color")` 
    stat_density_2d(data=re_int_pv, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                        group=dis, fill = ..level..), 
                    geom = "polygon", alpha=0.4)+
    scale_fill_gradient(low="white", high="olivedrab4") +
    geom_point(alpha = 1, aes(color = dis, size=dis)) +
    geom_text_repel(
      data = re_int,
      aes(x = ad_gene_sig, y = pv_gene_sig, label = sample),
      box.padding = 0.4
    ) + #add labels
    geom_errorbarh(aes(
      xmax = ad_gene_sig + ad_se,
      xmin = ad_gene_sig - ad_se,
      color = dis
    )) +
    geom_errorbar(aes(
      ymax = pv_gene_sig + pv_se,
      ymin = pv_gene_sig - pv_se,
      color = dis
    )) +
    ggtitle(paste(cell_types[p],"Single-cell Pseudobulked DEGs")) +
    theme_classic() +
    xlab("AD-specific genes") +
    ylab("PV-specific genes") +
    scale_color_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "PRP" = "darkred",
        "AD" = "burlywood4",
        "PV" = "olivedrab4"
        
      )
    ) +
    scale_size_manual(
      name = "Disease",
      values = c(
        "HC" = 2,
        "PRP" = 4,
        "AD" = 2,
        "PV" = 2
      )
    ) +
    theme(legend.position = "none")
  
  
  
  
  
  #--------------------------------------------------------------------#
  #--------------------------------------------------------------------#
  distance_optimiation_bulk_fx = function(ad_genes_input, pv_genes_input) {
    
    #ad_genes_input = ad_V_pv_genes
    #pv_genes_input = pv_V_ad_genes
    
    performance_met_df = do.call(rbind, lapply(2:nrow(ad_genes_input), function(i){
      #tryCatch({ #suppress error
      #i=2
      
      performance_met_df_inner = do.call(rbind, lapply(2:nrow(pv_genes_input), function(j){
        #j=2
        print(paste("i=",i, " j=",j,sep=""))
        
        #genes to be plotted
        ad_genes = ad_genes_input$gene[1:i]
        pv_genes = pv_genes_input$gene[1:j]
        
        #--------------------------------------------------------------------#
        #AD gene data_matrix
        ad_genes = na.omit(row.names(pseudobulk_dat)[match(ad_genes,row.names(pseudobulk_dat))])
        ad_mat = t(pseudobulk_dat@assays$RNA@data[ad_genes,])
        ad_int = data.frame(
          pseudobulk_dat$dis,
          pseudobulk_dat$sample_id_renamed,
          rowSums(ad_mat)
        )
        names(ad_int) = c("dis", "sample", "gene_sig")
        
        #--------------------------------------------------------------------#
        #PV gene data_matrix
        pv_genes = na.omit(row.names(pseudobulk_dat)[match(pv_genes,row.names(pseudobulk_dat))])
        pv_mat = t(pseudobulk_dat@assays$RNA@data[pv_genes,])
        pv_int = data.frame(
          pseudobulk_dat$dis,
          pseudobulk_dat$sample_id_renamed,
          rowSums(pv_mat)
        )
        names(pv_int) = c("dis", "sample", "gene_sig")
        
        #--------------------------------------------------------------------#
        re_int = cbind(ad_int[, c(1, 2, 3)], pv_int[, c(3)]) #combine the gene signatures
        
        re_int$dis = gsub("Healthy control","HC",re_int$dis)
        re_int$dis = gsub("Atopic Dermatitis","AD",re_int$dis)
        re_int$dis = gsub("Psoriasis Vulgaris","PV",re_int$dis)
        names(re_int)[c(3,4)] = c("ad_gene_sig","pv_gene_sig")
        
        
        #-----------------------------------
        #Calculate centers
        ad_centroid = data.frame( geometric.mean(subset(re_int, dis=="AD")$ad_gene_sig), geometric.mean(subset(re_int, dis=="AD")$pv_gene_sig) )
        pv_centroid = data.frame( geometric.mean(subset(re_int, dis=="PV")$ad_gene_sig), geometric.mean(subset(re_int, dis=="PV")$pv_gene_sig) )
        hc_centroid = data.frame( geometric.mean(subset(re_int, dis=="HC")$ad_gene_sig), geometric.mean(subset(re_int, dis=="HC")$pv_gene_sig) )
        names(ad_centroid) = c("xpt","ypt")
        names(pv_centroid) = c("xpt","ypt")
        names(hc_centroid) = c("xpt","ypt")
        
        centroid_mat = as.matrix(rbind(ad_centroid,pv_centroid,hc_centroid))
        row.names(centroid_mat) = c("AD","PV","HC")
        
        #-----------------------------------
        #Calculate within-disease stdevs
        ad_var = data.frame( sd(subset(re_int, dis=="AD")$ad_gene_sig), sd(subset(re_int, dis=="AD")$pv_gene_sig) )
        pv_var = data.frame( sd(subset(re_int, dis=="PV")$ad_gene_sig), sd(subset(re_int, dis=="PV")$pv_gene_sig) )
        hc_var = data.frame( sd(subset(re_int, dis=="HC")$ad_gene_sig), sd(subset(re_int, dis=="HC")$pv_gene_sig) )
        names(ad_var) = c("ad_variance","pv_variance")
        names(pv_var) = c("ad_variance","pv_variance")
        names(hc_var) = c("ad_variance","pv_variance")
        
        var_mat = as.matrix(rbind(ad_var,pv_var,hc_var))
        row.names(var_mat) = c("AD","PV","HC")
        
        #get SUMMED standard deviations across all samples
        sample_range = sum(as.vector(var_mat))
        
        #--------------------------------------------------------------------#
        #Calculate 20-80% quantile range among all samples
        
        #sample_level_range = hypotenus of AD range vs PV range
        # ad_range = quantile(re_int$ad_gene_sig,prob=0.8) - quantile(re_int$ad_gene_sig,prob=0.2)
        #  pv_range = quantile(re_int$pv_gene_sig,prob=0.8) - quantile(re_int$pv_gene_sig,prob=0.2)
        #  sample_range = (ad_range^2 + pv_range^2)^(1/2)
        
        #--------------------------------------------------------------------#
        #calculate distance between centroids
        dist_centroids = pointDistance(centroid_mat,allpairs=T, lonlat=F) #calculate all-vs-all distance
        
        #get SUMMED distances between centroids
        dist_centroids = sum(as.vector(dist_centroids)[as.vector(dist_centroids)>0])
        
        #calculate performance metric = amount of sample range taken up by distance between disease centroids
        performance_metric = dist_centroids / sample_range
        
        #Put them togehter 
        res_df = data.frame(i, j, ad_range, pv_range, sample_range, dist_centroids, performance_metric)
        names(res_df)[1:2] = c("ad_genes_n","pv_genes_n")
        res_df
        
      }))
      performance_met_df_inner
      
    }))
    
    
    performance_met_df = performance_met_df[order(-performance_met_df$performance_metric),]
    head(performance_met_df)
    
    #performance_met_df = performance_met_df[order(-performance_met_df$dist_centroids),]
    
    
    #--------------------------------------------------------------------#
    #RashX plotting of optimized gene number
    #--------------------------------------------------------------------#
    
    #reduce to unique genes that wil be plotted
    ad_genes = ad_V_pv_genes$gene[1:performance_met_df$ad_genes_n[300]]
    pv_genes = pv_V_ad_genes$gene[1:performance_met_df$pv_genes_n[300]]
    
    #--------------------------------------------------------------------#
    #AD gene data_matrix
    ad_genes = na.omit(row.names(pseudobulk_dat)[match(ad_genes,row.names(pseudobulk_dat))])
    ad_mat = t(pseudobulk_dat@assays$RNA@data[ad_genes,])
    ad_int = data.frame(
      pseudobulk_dat$dis,
      pseudobulk_dat$sample_id_renamed,
      rowSums(ad_mat)
    )
    names(ad_int) = c("dis", "sample", "gene_sig")
    
    #--------------------------------------------------------------------#
    #PV gene data_matrix
    pv_genes = na.omit(row.names(pseudobulk_dat)[match(pv_genes,row.names(pseudobulk_dat))])
    pv_mat = t(pseudobulk_dat@assays$RNA@data[pv_genes,])
    pv_int = data.frame(
      pseudobulk_dat$dis,
      pseudobulk_dat$sample_id_renamed,
      rowSums(pv_mat)
    )
    names(pv_int) = c("dis", "sample", "gene_sig")
    
    #--------------------------------------------------------------------#
    re_int = cbind(ad_int[, c(1, 2, 3)], pv_int[, c(3)]) #combine the gene signatures
    
    #----------------------------------------------------------------#
    ### Plot: Sample-level means without individual cells
    #----------------------------------------------------------------#
    re_int$dis = gsub("Healthy control","HC",re_int$dis)
    re_int$dis = gsub("Atopic Dermatitis","AD",re_int$dis)
    re_int$dis = gsub("Psoriasis Vulgaris","PV",re_int$dis)
    names(re_int)[c(3,4)] = c("ad_gene_sig","pv_gene_sig")
    
    re_int_ad = subset(re_int, dis == "AD")
    re_int_hc = subset(re_int, dis == "HC")
    re_int_pv = subset(re_int, dis == "PV")
    
    #----------------------------------------------------------------#
    ### Plot: Sample-level means without individual cells
    #----------------------------------------------------------------#
    
    p1.5 <- ggplot(data = re_int, aes(x = ad_gene_sig, y = pv_gene_sig)) +
      stat_density_2d(data=re_int_hc, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                          group=dis, fill = ..level..), 
                      geom = "polygon", alpha=0.4)+
      scale_fill_gradient(low="white", high="darkblue") +
      new_scale_fill() +   # same as `new_scale("color")` 
      stat_density_2d(data=re_int_ad, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                          group=dis, fill = ..level..), 
                      geom = "polygon", alpha=0.4)+
      scale_fill_gradient(low="white", high="burlywood4") +
      new_scale_fill() +   # same as `new_scale("color")` 
      stat_density_2d(data=re_int_pv, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                          group=dis, fill = ..level..), 
                      geom = "polygon", alpha=0.4)+
      scale_fill_gradient(low="white", high="olivedrab4") +
      geom_point(alpha = 1, aes(color = dis, size=dis)) +
      geom_text_repel(
        data = re_int,
        aes(x = ad_gene_sig, y = pv_gene_sig, label = sample),
        box.padding = 0.4
      ) + #add labels
      ggtitle(paste(cell_types[i],"Pseudobulked samples DESeq2 top 50 AD-PV genes")) +
      theme_classic() +
      xlab("AD-specific genes") +
      ylab("PV-specific genes") +
      scale_color_manual(
        name = "Disease",
        values = c(
          "HC" = "darkblue",
          "PRP" = "darkred",
          "AD" = "burlywood4",
          "PV" = "olivedrab4"
          
        )
      ) +
      scale_size_manual(
        name = "Disease",
        values = c(
          "HC" = 2,
          "PRP" = 4,
          "AD" = 2,
          "PV" = 2
        )
      ) +
      theme(legend.position = "none")
    
    p1.5
    
    
    
    #--------------------------------------------------------------------#
    #loop through DEGs from psuedobulk DEG analysis and replot!
    #--------------------------------------------------------------------#
    lowest_common_den = min(nrow(ad_V_pv_genes),nrow(pv_V_ad_genes))
    results_name="AD-PV_paper_samples"
    for(i in 2:lowest_common_den) {
      
      #i=2
      print(i)
      #reduce to unique genes that wil be plotted
      ad_genes = ad_V_pv_genes$gene[1:i]
      pv_genes = pv_V_ad_genes$gene[1:i]
      
      #--------------------------------------------------------------------#
      #AD gene data_matrix
      ad_genes = na.omit(row.names(pseudobulk_dat)[match(ad_genes,row.names(pseudobulk_dat))])
      ad_mat = t(pseudobulk_dat@assays$RNA@data[ad_genes,])
      ad_int = data.frame(
        pseudobulk_dat$dis,
        pseudobulk_dat$sample_id_renamed,
        rowSums(ad_mat)
      )
      names(ad_int) = c("dis", "sample", "gene_sig")
      
      #--------------------------------------------------------------------#
      #PV gene data_matrix
      pv_genes = na.omit(row.names(pseudobulk_dat)[match(pv_genes,row.names(pseudobulk_dat))])
      pv_mat = t(pseudobulk_dat@assays$RNA@data[pv_genes,])
      pv_int = data.frame(
        pseudobulk_dat$dis,
        pseudobulk_dat$sample_id_renamed,
        rowSums(pv_mat)
      )
      names(pv_int) = c("dis", "sample", "gene_sig")
      
      #--------------------------------------------------------------------#
      re_int = cbind(ad_int[, c(1, 2, 3)], pv_int[, c(3)]) #combine the gene signatures
      
      #----------------------------------------------------------------#
      ### Plot: Sample-level means without individual cells
      #----------------------------------------------------------------#
      re_int$dis = gsub("Healthy control","HC",re_int$dis)
      re_int$dis = gsub("Atopic Dermatitis","AD",re_int$dis)
      re_int$dis = gsub("Psoriasis Vulgaris","PV",re_int$dis)
      names(re_int)[c(3,4)] = c("ad_gene_sig","pv_gene_sig")
      
      re_int_ad = subset(re_int, dis == "AD")
      re_int_hc = subset(re_int, dis == "HC")
      re_int_pv = subset(re_int, dis == "PV")
      
      #----------------------------------------------------------------#
      ### Plot: Sample-level means without individual cells
      #----------------------------------------------------------------#
      
      p1.5 <- ggplot(data = re_int, aes(x = ad_gene_sig, y = pv_gene_sig)) +
        stat_density_2d(data=re_int_hc, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                            group=dis, fill = ..level..), 
                        geom = "polygon", alpha=0.4)+
        scale_fill_gradient(low="white", high="darkblue") +
        new_scale_fill() +   # same as `new_scale("color")` 
        stat_density_2d(data=re_int_ad, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                            group=dis, fill = ..level..), 
                        geom = "polygon", alpha=0.4)+
        scale_fill_gradient(low="white", high="burlywood4") +
        new_scale_fill() +   # same as `new_scale("color")` 
        stat_density_2d(data=re_int_pv, aes(x = ad_gene_sig, y = pv_gene_sig, 
                                            group=dis, fill = ..level..), 
                        geom = "polygon", alpha=0.4)+
        scale_fill_gradient(low="white", high="olivedrab4") +
        geom_point(alpha = 1, aes(color = dis, size=dis)) +
        geom_text_repel(
          data = re_int,
          aes(x = ad_gene_sig, y = pv_gene_sig, label = sample),
          box.padding = 0.4
        ) + #add labels
        ggtitle(paste(cell_types[8]," Pseudobulked samples DESeq2 AD-PV genes, n=",i,sep="")) +
        theme_classic() +
        xlab("AD-specific genes") +
        ylab("PV-specific genes") +
        scale_color_manual(
          name = "Disease",
          values = c(
            "HC" = "darkblue",
            "PRP" = "darkred",
            "AD" = "burlywood4",
            "PV" = "olivedrab4"
            
          )
        ) +
        scale_size_manual(
          name = "Disease",
          values = c(
            "HC" = 2,
            "PRP" = 4,
            "AD" = 2,
            "PV" = 2
          )
        ) +
        theme(legend.position = "none")
      
      png(paste("/c4/home/mtaylor4/derm/ad_pv_project/results/pseudobulk/results/hyper_d_plot_loop/",results_name,"_",i,"_",cell_types[8],
                "hyperD_comparison_plots.png",sep=""),
          width=8,height=8,units="in",res=400)
      print(p1.5)
      dev.off()
    }
    
    scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/ad_pv_project/results/pseudobulk/results/hyper_d_plot_loop/ /Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/results/pseudobulk
    
    
    
    #--------------------------------------------------------------------#
    #--------------------------------------------------------------------#
    #Batch-correct
    #Correct for sequencing platform bitch ;)
    
    #Must re-run pseudobulking on un-normalized counts bc DESeq2 requires integers
    pseudobulk_dat <- AggregateExpression(input_cds, assays = "RNA", slot="counts",
                                          return.seurat = T, group.by = c("dis", "sample_id_renamed"),
                                          normalization.method=NULL)
    
    sample_metadat_df = t(data.frame(strsplit(row.names(pseudobulk_dat@meta.data),"_")))
    pseudobulk_dat@meta.data$dis = sample_metadat_df[,1]
    pseudobulk_dat@meta.data$sample_id_renamed = sample_metadat_df[,2]
    
    #Check average and total reads
    pseudobulk_dat@assays$RNA@data[1:10,1:10]
    pseudobulk_dat@assays$RNA@counts[1:10,1:10]
    
    #get raw reads counts
    pseudobulk_count_dat = pseudobulk_dat@assays$RNA@counts
    
    #Calculate sample-wise MEAN reads across genes
    pseudobulk_count_means = data.frame(colnames(pseudobulk_count_dat),colMeans(pseudobulk_count_dat))
    names(pseudobulk_count_means) = c("sample","mean_counts")
    pseudobulk_count_means = pseudobulk_count_means[order(pseudobulk_count_means$mean_counts),]
    
    #Calculate sample-wise TOTAL reads across genes
    pseudobulk_count_sums = data.frame(colnames(pseudobulk_count_dat),colSums(pseudobulk_count_dat))
    names(pseudobulk_count_sums) = c("sample","sum_counts")
    pseudobulk_count_sums = pseudobulk_count_sums[order(pseudobulk_count_sums$sum_counts),]
    
    #Single cell counts
    cell_counts = data.frame(table(suerat_obj@meta.data$sample_id_renamed))
    cell_counts = cell_counts[order(cell_counts$Freq),]
    
    #--------------------------------------------------------------------#
    #--------------------------------------------------------------------#
    
  } #end pseudo-bulk function
  
  scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/ad_pv_project/results/pseudobulk/results/ /Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/results/pseudobulk
  
  
  
  
  #DEG recovery function definition
  
  #downsample to similar numbers of cells
  Idents(suerat_obj_sub) <- "sample_id_dis"
  suerat_obj_sub = subset(x = suerat_obj_sub, downsample = 300)
  data.frame(table(suerat_obj_sub@meta.data$sample_id_dis))
  
  #Get rid of 394
  #suerat_obj_sub = subset(x = suerat_obj_sub, subset = sample_id_dis != "Atopic Dermatitis 394")
  
  #Get samples in a cell type subset
  all_sample_names = unique(suerat_obj_sub@meta.data$sample_id_dis)
  hc_sample_names = all_sample_names[grep("Healthy control",all_sample_names)]
  pv_sample_names = c( all_sample_names[grep("Psoriasis Vulgaris",all_sample_names)])
  ad_sample_names = c( all_sample_names[grep("Atopic Dermatitis",all_sample_names)] )
  
  suerat_obj_sub <- FindVariableFeatures(suerat_obj_sub, selection.method = "vst", nfeatures = 10000)
  #all.genes <- rownames(suerat_obj_sub)
  #suerat_obj_sub <- ScaleData(suerat_obj_sub, features = all.genes)
  
  #------------------------------------------------------------------------#
  #PV DEGs 
  #------------------------------------------------------------------------#
  
  deg_df_pv <- FindMarkers(suerat_obj_sub,
                           ident.1 = row.names(suerat_obj_sub@meta.data[suerat_obj_sub@meta.data$sample_id_dis %in% pv_sample_names,]),
                           ident.2 = row.names(suerat_obj_sub@meta.data[suerat_obj_sub@meta.data$sample_id_dis %in% hc_sample_names,])
  )
  
  deg_df_pv = deg_df_pv[order(-deg_df_pv$avg_log2FC),] #Re-order from high to low FC
  deg_df_pv = cbind(row.names(deg_df_pv),deg_df_pv) #Add col of gene symbols
  names(deg_df_pv)[1] = "gene"
  
  #------------------------------------------------------------------------#
  #AD DEGs 
  #------------------------------------------------------------------------#
  deg_df_ad <- FindMarkers(suerat_obj_sub,
                           ident.1 = row.names(suerat_obj_sub@meta.data[suerat_obj_sub@meta.data$sample_id_dis %in% ad_sample_names,]),
                           ident.2 = row.names(suerat_obj_sub@meta.data[suerat_obj_sub@meta.data$sample_id_dis %in% hc_sample_names,])
  )
  
  deg_df_ad = deg_df_ad[order(-deg_df_ad$avg_log2FC),] #Re-order from high to low FC
  deg_df_ad = cbind(row.names(deg_df_ad),deg_df_ad) #Add col of gene symbols
  names(deg_df_ad)[1] = "gene"
  
  
  
  #------------------------------------------------------------------------#
  #RashX cell type
  #------------------------------------------------------------------------#
  tryCatch({ #suppress error
    
    library(ggforce)
    library("concaveman")
    
    min_percenteage_expressed = 0.4 #the amount of cells that a gene must be detected in
    gene_num = 1000 #gene number to include in hyperdimensionality plot for AD and PV
    
    pv_genes_all = subset(deg_df_pv, avg_log2FC>0 & p_val_adj<0.001 & pct.1 > min_percenteage_expressed)$gene
    ad_genes_all = subset(deg_df_ad, avg_log2FC>0 & p_val_adj<0.001 & pct.1 > min_percenteage_expressed)$gene
    
    #Now get genes that are UNIQUE to PV or AD
    pv_genes = setdiff(pv_genes_all,ad_genes_all)
    ad_genes = setdiff(ad_genes_all,pv_genes_all)
    
    #Get top 50 genes if there are more than gene_num (defined above)
    if(length(pv_genes) > gene_num ) {pv_genes = pv_genes[1:gene_num]} else {pv_genes = pv_genes}
    if(length(ad_genes) > gene_num ) {ad_genes = ad_genes[1:gene_num]} else {ad_genes = ad_genes}
    
    
    #--------------------------------------------------------------------#
    #Get the genes of interest
    #--------------------------------------------------------------------#
    #AD gene data_matrix
    ad_mat = t(suerat_obj_sub@assays$RNA[ad_genes,])
    ad_mat[1:10, 1:10]
    ad_int = data.frame(
      suerat_obj_sub$dis,
      suerat_obj_sub$sample_id_renamed,
      rowSums(ad_mat)
    )
    names(ad_int) = c("dis", "sample", "gene_sig")
    names(ad_int) = c("dis", "sample", "gene_sig")
    ad_dfc = summarySE(ad_int,
                       measurevar = 'gene_sig',
                       groupvars = c("dis", "sample"))
    names(ad_dfc)[c(4, 6)] = c("ad_gene_sig", "ad_se")
    
    #--------------------------------------------------------------------#
    #PV gene data_matrix
    pv_mat = t(suerat_obj_sub@assays$RNA[pv_genes,])
    pv_mat[1:10, 1:10]
    pv_int = data.frame(
      suerat_obj_sub$dis,
      suerat_obj_sub$sample_id_renamed,
      rowSums(pv_mat)
    )
    names(pv_int) = c("dis", "sample", "gene_sig")
    names(pv_int) = c("dis", "sample", "gene_sig")
    pv_dfc = summarySE(pv_int,
                       measurevar = 'gene_sig',
                       groupvars = c("dis", "sample"))
    head(pv_dfc)
    names(pv_dfc)[c(4, 6)] = c("pv_gene_sig", "pv_se")
    
    #--------------------------------------------------------------------#
    re_int = cbind(ad_dfc[, c(1, 2, 4, 6)], pv_dfc[, c(4, 6)]) #combine the gene signatures
    
    #----------------------------------------------------------------#
    ### Plot: Sample-level means without individual cells
    #----------------------------------------------------------------#
    re_int$dis = gsub("Healthy control","HC",re_int$dis)
    re_int$dis = gsub("Atopic Dermatitis","AD",re_int$dis)
    re_int$dis = gsub("Psoriasis Vulgaris","PV",re_int$dis)
    unique(re_int$dis)
    
    #re_int = subset(re_int, dis != "HC")
    
    p1 <- ggplot(data = re_int, aes(x = ad_gene_sig, y = pv_gene_sig)) +
      geom_point(alpha = 1, aes(color = dis, size=dis)) +
      geom_text_repel(
        data = re_int,
        aes(x = ad_gene_sig, y = pv_gene_sig, label = sample),
        box.padding = 0.4
      ) + #add labels
      geom_errorbarh(aes(
        xmax = ad_gene_sig + ad_se,
        xmin = ad_gene_sig - ad_se,
        color = dis
      )) +
      geom_errorbar(aes(
        ymax = pv_gene_sig + pv_se,
        ymin = pv_gene_sig - pv_se,
        color = dis
      )) +
      ggtitle(paste(cell_types[i],"no sample filter")) +
      theme_classic() +
      xlab("AD-specific genes") +
      ylab("PV-specific genes") +
      theme(legend.position = "right") +
      scale_color_manual(
        name = "Disease",
        values = c(
          "HC" = "darkblue",
          "PRP" = "darkred",
          "AD" = "burlywood4",
          "PV" = "olivedrab4"
          
        )
      ) +
      scale_fill_manual(
        name = "Disease",
        values = c(
          "HC" = "darkblue",
          "PRP" = "darkred",
          "AD" = "burlywood4",
          "PV" = "olivedrab4"
        )
      ) +
      geom_mark_hull(data = re_int, aes(fill = dis, label = dis), concavity = 3) +
      theme(
        axis.text = element_text(size = 15, color = 'black'),
        axis.title = element_text(size = 15, color = 'black'),
        plot.title = element_text(size = 20, color = 'black')
      ) +
      scale_size_manual(
        name = "Disease",
        values = c(
          "HC" = 2,
          "PRP" = 4,
          "AD" = 2,
          "PV" = 2
        )
      )
    
    png(paste(deg_dir,i,"_",cell_types[i],"_hyperdim_plot_nosamplefilter.png"),
        width=8,height=8,units="in",res=400)
    print(p1)
    dev.off()
  }, error=function(e){})
  
  
  #-----------------------------------------------------------------------------------------#
  #-----------------------------------------------------------------------------------------#
  #-----------------------------------------------------------------------------------------#
  #INDIVIDUAL GENE ACROSS-SAMPLE FILTERS
  tryCatch({ 
    
    min_avg_expression = 1 #this is the cutoff for a particular sample, the min average expression of a single gene for all the cells in that sample
    min_samples_percent = 0.8  #This is the cutoff for a particular gene, how many samples must have the min_avg_expression to be included as a gene in the gene set
    
    #--------------------------------------------------------------------#
    #AD gene data_matrix
    ad_mat = t(suerat_obj_sub@assays$RNA[ad_genes,])
    ad_mat[1:10, 1:10]
    ad_int = data.frame(
      suerat_obj_sub$dis,
      suerat_obj_sub$sample_id_renamed,
      ad_mat
    )
    ad_int[1:10,1:10]
    names(ad_int)[1:2] = c("dis", "sample")
    
    #get only AD samples
    ad_int = subset(ad_int, dis=="Atopic Dermatitis")
    sample_num = length(unique(ad_int$sample))
    sample_expt_cutoff = round(sample_num*min_samples_percent) #get number of samples a gene must be expressed in to be considered
    
    # Group by factor levels and calculate the mean for each group
    library(dplyr)
    averages <- data.frame(ad_int %>%
                             group_by(dis, sample) %>%
                             summarise_all(mean))
    
    #Now find the genes that are robstly expressed on average across all my lovely samples
    gene_breadth = do.call(rbind, lapply(3:ncol(averages), function(x){
      #x=3
      sample_num_above_min = data.frame(names(averages)[x], length( averages[,x][ averages[,x] > min_avg_expression ] ) )
      names(sample_num_above_min) = c("gene","sample_above")
      sample_num_above_min
    }))
    
    ad_gene_robust = subset(gene_breadth, sample_above>=sample_expt_cutoff)
    ad_genes_samplefiltered = ad_gene_robust$gene
    
    #--------------------------------------------------------------------#
    #PV gene data_matrix
    pv_mat = t(suerat_obj_sub@assays$RNA[pv_genes,])
    pv_mat[1:10, 1:10]
    pv_int = data.frame(
      suerat_obj_sub$dis,
      suerat_obj_sub$sample_id_renamed,
      pv_mat
    )
    pv_int[1:10,1:10]
    names(pv_int)[1:2] = c("dis", "sample")
    
    #get only PV samples
    pv_int = subset(pv_int, dis=="Psoriasis Vulgaris")
    sample_num = length(unique(pv_int$sample))
    sample_expt_cutoff = round(sample_num*min_samples_percent) #get number of samples a gene must be expressed in to be considered
    
    # Group by factor levels and calculate the mean for each group
    averages <- data.frame(pv_int %>%
                             group_by(dis, sample) %>%
                             summarise_all(mean))
    
    #Now find the genes that are robstly expressed on average across all my lovely samples
    gene_breadth = do.call(rbind, lapply(3:ncol(averages), function(x){
      #x=3
      sample_num_above_min = data.frame(names(averages)[x], length( averages[,x][ averages[,x] > min_avg_expression ] ) )
      names(sample_num_above_min) = c("gene","sample_above")
      sample_num_above_min
    }))
    
    pv_gene_robust = subset(gene_breadth, sample_above>=sample_expt_cutoff)
    pv_genes_samplefiltered = pv_gene_robust$gene
    
    #-----------------------------------------------------------------------------------------#
    #RashX on sample filtered genes
    
    #--------------------------------------------------------------------#
    #Get the genes of interest
    #--------------------------------------------------------------------#
    #AD gene data_matrix
    tryCatch({ unloadNamespace("dplyr") }, error=function(e){}) #unload dplyr so that summarySE will run
    ad_mat = t(suerat_obj_sub@assays$RNA[ad_genes_samplefiltered,])
    ad_int = data.frame(
      suerat_obj_sub$dis,
      suerat_obj_sub$sample_id_renamed,
      rowSums(ad_mat)
    )
    names(ad_int) = c("dis", "sample", "gene_sig")
    names(ad_int) = c("dis", "sample", "gene_sig")
    ad_dfc = summarySE(ad_int,
                       measurevar = 'gene_sig',
                       groupvars = c("dis", "sample"))
    names(ad_dfc)[c(4, 6)] = c("ad_gene_sig", "ad_se")
    
    #--------------------------------------------------------------------#
    #PV gene data_matrix
    pv_mat = t(suerat_obj_sub@assays$RNA[pv_genes_samplefiltered,])
    pv_int = data.frame(
      suerat_obj_sub$dis,
      suerat_obj_sub$sample_id_renamed,
      rowSums(pv_mat)
    )
    names(pv_int) = c("dis", "sample", "gene_sig")
    names(pv_int) = c("dis", "sample", "gene_sig")
    pv_dfc = summarySE(pv_int,
                       measurevar = 'gene_sig',
                       groupvars = c("dis", "sample"))
    head(pv_dfc)
    names(pv_dfc)[c(4, 6)] = c("pv_gene_sig", "pv_se")
    
    #--------------------------------------------------------------------#
    re_int = cbind(ad_dfc[, c(1, 2, 4, 6)], pv_dfc[, c(4, 6)]) #combine the gene signatures
    
    #----------------------------------------------------------------#
    ### Plot: Sample-level means without individual cells
    #----------------------------------------------------------------#
    re_int$dis = gsub("Healthy control","HC",re_int$dis)
    re_int$dis = gsub("Atopic Dermatitis","AD",re_int$dis)
    re_int$dis = gsub("Psoriasis Vulgaris","PV",re_int$dis)
    unique(re_int$dis)
    
    p2 <- ggplot(data = re_int, aes(x = ad_gene_sig, y = pv_gene_sig)) +
      geom_point(alpha = 1, aes(color = dis, size=dis)) +
      geom_text_repel(
        data = re_int,
        aes(x = ad_gene_sig, y = pv_gene_sig, label = sample),
        box.padding = 0.4
      ) + #add labels
      geom_errorbarh(aes(
        xmax = ad_gene_sig + ad_se,
        xmin = ad_gene_sig - ad_se,
        color = dis
      )) +
      geom_errorbar(aes(
        ymax = pv_gene_sig + pv_se,
        ymin = pv_gene_sig - pv_se,
        color = dis
      )) +
      ggtitle(paste(cell_types[i],"ind gene across-sample filter")) +
      theme_classic() +
      xlab("AD-specific genes") +
      ylab("PV-specific genes") +
      theme(legend.position = "right") +
      scale_color_manual(
        name = "Disease",
        values = c(
          "HC" = "darkblue",
          "PRP" = "darkred",
          "AD" = "burlywood4",
          "PV" = "olivedrab4"
          
        )
      ) +
      scale_fill_manual(
        name = "Disease",
        values = c(
          "HC" = "darkblue",
          "PRP" = "darkred",
          "AD" = "burlywood4",
          "PV" = "olivedrab4"
        )
      ) +
      geom_mark_hull(data = re_int, aes(fill = dis, label = dis), concavity = 3) +
      theme(
        axis.text = element_text(size = 15, color = 'black'),
        axis.title = element_text(size = 15, color = 'black'),
        plot.title = element_text(size = 20, color = 'black')
      ) +
      scale_size_manual(
        name = "Disease",
        values = c(
          "HC" = 2,
          "PRP" = 4,
          "AD" = 2,
          "PV" = 2
        )
      )
    
    png(paste(deg_dir,i,"_",cell_types[i],"_hyperdim_plot_indgenefilter.png"),
        width=8,height=8,units="in",res=400)
    print(p2)
    dev.off()
  }, error=function(e){})
  
}, error=function(e){})



pseudobulk_analysis_fx(suerat_obj_sub, "AD_PV_paper_samples","CLR")
pseudobulk_analysis_fx(suerat_obj_sub, "AD_PV_paper_samples","CLR")


#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#Get sample specific indvidiual disease vs grouped HC DEGs
tryCatch({ #suppress error
  
  suerat_obj_sub = subset(suerat_obj, reference_celltype == cell_types[i])
  suerat_obj_sub <- FindVariableFeatures(suerat_obj_sub, selection.method = "vst", nfeatures = 20000)
  
  all_sample_names = unique(suerat_obj_sub@meta.data$sample_id_dis)
  hc_sample_names = all_sample_names[grep("Healthy control",all_sample_names)]
  pv_sample_names = c( all_sample_names[grep("Psoriasis Vulgaris",all_sample_names)])
  ad_sample_names = c( all_sample_names[grep("Atopic Dermatitis",all_sample_names)] )
  dis_sample_names = c(pv_sample_names,ad_sample_names)
  
  #Loop through each individual sample and compare to grouped normals
  ind_dis_deg = do.call(rbind, lapply(1:length(dis_sample_names), function(y) {
    
    #y=1
    print(paste("y=",y,sep=""))
    #------------------------------------------------------------------------#
    #Individual disease vs grouped normals
    #------------------------------------------------------------------------#
    
    deg_df <- FindMarkers(suerat_obj_sub,
                          ident.1 = row.names(suerat_obj_sub@meta.data[suerat_obj_sub@meta.data$sample_id_dis %in% dis_sample_names[y],]),
                          ident.2 = row.names(suerat_obj_sub@meta.data[suerat_obj_sub@meta.data$sample_id_dis %in% hc_sample_names,])
    )
    
    deg_df = deg_df[order(-deg_df$avg_log2FC),] #Re-order from high to low FC
    deg_df = cbind(dis_sample_names[y],row.names(deg_df),deg_df) #Add col of gene symbols
    names(deg_df)[1:2] = c("sample","gene")
    deg_df
    
  })) #End indivdiual sample DEG loop
  
  #For each sample get the genes that are unique to AD or PV and are individual DEGs
  #Now get a gene set for each sample
  ind_dis_deg_filtered = subset(ind_dis_deg, avg_log2FC>0 & p_val_adj<0.001 & pct.1 > min_percenteage_expressed)
  
  #-------------------------------
  #-------------------------------
  #-------------------------------
  #Individual PV gene set
  pv_ind_gene_set = do.call(rbind, lapply(1:length(pv_sample_names), function(nn){
    #nn=1
    ind_pv_genes = subset(ind_dis_deg_filtered, sample==pv_sample_names[nn])$gene
    ind_pv_genes <- ind_pv_genes[!(ind_pv_genes %in% ad_genes)] #Remove any AD specific genes
    return_df = data.frame(pv_sample_names[nn],ind_pv_genes)
    names(return_df) = c("sample","pv_specific_gene")
    return_df
  }))
  
  #Get lowest common denominator for a sample gene
  df_num =  data.frame(table(pv_ind_gene_set$sample))
  lowest_den = min(df_num$Freq)
  lowest_den_pv = lowest_den
  
  #-------------------------------
  #AD programs
  ad_ind_gene_set = do.call(rbind, lapply(1:length(ad_sample_names), function(nn){
    #nn=1
    ind_ad_genes = subset(ind_dis_deg_filtered, sample==ad_sample_names[nn])$gene
    ind_ad_genes <- ind_ad_genes[!(ind_ad_genes %in% pv_genes)] #Remove any PV specific genes
    return_df = data.frame(ad_sample_names[nn],ind_ad_genes)
    names(return_df) = c("sample","ad_specific_gene")
    return_df
  }))
  
  #Get lowest common denominator for a sample gene
  df_num = data.frame(table(ad_ind_gene_set$sample))
  lowest_den = min(df_num$Freq)
  lowest_den_ad = lowest_den
  
  #-------------------------------
  if(length(ad_genes)<lowest_den_ad) {num_ad_genes_to_use = length(ad_genes)} else {num_ad_genes_to_use = lowest_den_ad}
  if(length(pv_genes)<lowest_den_pv) {num_pv_genes_to_use = length(pv_genes)} else {num_pv_genes_to_use = lowest_den_pv}
  
  #-------------------------------
  #Get specific sample's DEG's as well as global AD sample's AD genes
  
  pv_ind_gene_set_df = do.call(cbind, lapply(1:length(pv_sample_names), function(nn){
    #nn=1
    tryCatch({ 
      ind_pv_genes = subset(ind_dis_deg_filtered, sample==pv_sample_names[nn])$gene #subset to a specific sample's DEGs
      ind_pv_genes <- ind_pv_genes[!(ind_pv_genes %in% ad_genes)] #Remove any AD specific genes
      return_df = data.frame(pv_sample_names[nn],ind_pv_genes)
      names(return_df) = c("sample","pv_specific_gene")
      sample_gene_df = data.frame(return_df$pv_specific_gene[1:num_pv_genes_to_use])
      names(sample_gene_df) = pv_sample_names[nn]
      sample_gene_df
    }, error=function(e){})
  }))
  
  pv_ind_gene_set = do.call(rbind, lapply(1:length(pv_sample_names), function(nn){
    #nn=1
    tryCatch({ 
      ind_pv_genes = subset(ind_dis_deg_filtered, sample==pv_sample_names[nn])$gene #subset to a specific sample's DEGs
      ind_pv_genes <- ind_pv_genes[!(ind_pv_genes %in% ad_genes)] #Remove any AD specific genes
      return_df = data.frame(pv_sample_names[nn],ind_pv_genes)
      names(return_df) = c("sample","pv_specific_gene")
      gene_size = return_df$pv_specific_gene[1:num_pv_genes_to_use]
      
      #AD gene data_matrix
      ad_mat = t(suerat_obj_sub@assays$RNA[ad_genes[1:num_ad_genes_to_use],])
      ad_int = data.frame(
        suerat_obj_sub$dis,
        suerat_obj_sub$sample_id_renamed,
        suerat_obj_sub$sample_id_dis,
        rowSums(ad_mat)
      )
      
      names(ad_int) = c("dis", "sample","sample_name" ,"gene_sig")
      ad_int = subset(ad_int, sample_name==pv_sample_names[nn])
      ad_int = ad_int[,-3]
      
      ad_dfc = summarySE(ad_int,
                         measurevar = 'gene_sig',
                         groupvars = c("dis", "sample"))
      names(ad_dfc)[c(4, 6)] = c("ad_gene_sig", "ad_se")
      
      
      #PV gene data_matrix
      pv_mat = t(suerat_obj_sub@assays$RNA[gene_size,])
      pv_int = data.frame(
        suerat_obj_sub$dis,
        suerat_obj_sub$sample_id_renamed,
        suerat_obj_sub$sample_id_dis,
        rowSums(pv_mat)
      )
      names(pv_int) = c("dis", "sample","sample_name" ,"gene_sig")
      pv_int = subset(pv_int, sample_name==pv_sample_names[nn])
      pv_int = pv_int[,-3]
      
      pv_dfc = summarySE(pv_int,
                         measurevar = 'gene_sig',
                         groupvars = c("dis", "sample"))
      names(pv_dfc)[c(4, 6)] = c("pv_gene_sig", "pv_se")
      
      re_int = cbind(ad_dfc[, c(1, 2, 4, 6)], pv_dfc[, c(4, 6)]) #combine the gene signatures
      
      return(re_int)
    }, error=function(e){})
  }))
  
  #-----------------#
  
  ad_ind_gene_set_df = do.call(cbind, lapply(1:length(ad_sample_names), function(nn){
    #nn=1
    tryCatch({ 
      ind_ad_genes = subset(ind_dis_deg_filtered, sample==ad_sample_names[nn])$gene
      ind_ad_genes <- ind_ad_genes[!(ind_ad_genes %in% pv_genes)] #Remove any PV specific genes
      return_df = data.frame(ad_sample_names[nn],ind_ad_genes)
      names(return_df) = c("sample","ad_specific_gene")
      sample_gene_df = data.frame(return_df$ad_specific_gene[1:num_ad_genes_to_use])
      names(sample_gene_df) = ad_sample_names[nn]
      sample_gene_df
    }, error=function(e){})
  }))
  
  ad_ind_gene_set = do.call(rbind, lapply(1:length(ad_sample_names), function(nn){
    #nn=1
    tryCatch({ 
      print(nn)
      ind_ad_genes = subset(ind_dis_deg_filtered, sample==ad_sample_names[nn])$gene
      ind_ad_genes <- ind_ad_genes[!(ind_ad_genes %in% pv_genes)] #Remove any PV specific genes
      return_df = data.frame(ad_sample_names[nn],ind_ad_genes)
      names(return_df) = c("sample","ad_specific_gene")
      gene_size = return_df$ad_specific_gene[1:num_ad_genes_to_use]
      
      #AD gene data_matrix
      ad_mat = t(suerat_obj_sub@assays$RNA[gene_size,])
      ad_int = data.frame(
        suerat_obj_sub$dis,
        suerat_obj_sub$sample_id_renamed,
        suerat_obj_sub$sample_id_dis,
        rowSums(ad_mat)
      )
      
      names(ad_int) = c("dis", "sample","sample_name" ,"gene_sig")
      ad_int = subset(ad_int, sample_name==ad_sample_names[nn])
      ad_int = ad_int[,-3]
      
      ad_dfc = summarySE(ad_int,
                         measurevar = 'gene_sig',
                         groupvars = c("dis", "sample"))
      names(ad_dfc)[c(4, 6)] = c("ad_gene_sig", "ad_se")
      
      
      #PV gene data_matrix
      pv_mat = t(suerat_obj_sub@assays$RNA[pv_genes[1:num_pv_genes_to_use],])
      pv_int = data.frame(
        suerat_obj_sub$dis,
        suerat_obj_sub$sample_id_renamed,
        suerat_obj_sub$sample_id_dis,
        rowSums(pv_mat)
      )
      names(pv_int) = c("dis", "sample","sample_name" ,"gene_sig")
      pv_int = subset(pv_int, sample_name==ad_sample_names[nn])
      pv_int = pv_int[,-3]
      
      pv_dfc = summarySE(pv_int,
                         measurevar = 'gene_sig',
                         groupvars = c("dis", "sample"))
      names(pv_dfc)[c(4, 6)] = c("pv_gene_sig", "pv_se")
      
      re_int = cbind(ad_dfc[, c(1, 2, 4, 6)], pv_dfc[, c(4, 6)]) #combine the gene signatures
      
      return(re_int)
    }, error=function(e){})
  }))
  
  hc_ind_gene_set = do.call(rbind, lapply(1:length(hc_sample_names), function(nn){
    #nn=1
    tryCatch({ 
      #AD gene data_matrix
      ad_mat = t(suerat_obj_sub@assays$RNA[ad_genes[1:num_ad_genes_to_use],])
      ad_int = data.frame(
        suerat_obj_sub$dis,
        suerat_obj_sub$sample_id_renamed,
        suerat_obj_sub$sample_id_dis,
        rowSums(ad_mat)
      )
      
      names(ad_int) = c("dis", "sample","sample_name" ,"gene_sig")
      ad_int = subset(ad_int, sample_name==hc_sample_names[nn])
      ad_int = ad_int[,-3]
      
      ad_dfc = summarySE(ad_int,
                         measurevar = 'gene_sig',
                         groupvars = c("dis", "sample"))
      names(ad_dfc)[c(4, 6)] = c("ad_gene_sig", "ad_se")
      
      
      #PV gene data_matrix
      pv_mat = t(suerat_obj_sub@assays$RNA[pv_genes[1:num_pv_genes_to_use],])
      pv_int = data.frame(
        suerat_obj_sub$dis,
        suerat_obj_sub$sample_id_renamed,
        suerat_obj_sub$sample_id_dis,
        rowSums(pv_mat)
      )
      names(pv_int) = c("dis", "sample","sample_name" ,"gene_sig")
      pv_int = subset(pv_int, sample_name==hc_sample_names[nn])
      pv_int = pv_int[,-3]
      
      pv_dfc = summarySE(pv_int,
                         measurevar = 'gene_sig',
                         groupvars = c("dis", "sample"))
      names(pv_dfc)[c(4, 6)] = c("pv_gene_sig", "pv_se")
      
      re_int = cbind(ad_dfc[, c(1, 2, 4, 6)], pv_dfc[, c(4, 6)]) #combine the gene signatures
      
      return(re_int)
    }, error=function(e){})
  }))
  
  
  #----------------------------------------------------------------#
  ### Plot: Sample-level means without individual cells
  #----------------------------------------------------------------#
  re_int = rbind(ad_ind_gene_set,pv_ind_gene_set,hc_ind_gene_set)
  re_int$dis = gsub("Healthy control","HC",re_int$dis)
  re_int$dis = gsub("Atopic Dermatitis","AD",re_int$dis)
  re_int$dis = gsub("Psoriasis Vulgaris","PV",re_int$dis)
  unique(re_int$dis)
  
  #re_int = subset(re_int, dis != "HC")
  
  p3 <- ggplot(data = re_int, aes(x = ad_gene_sig, y = pv_gene_sig)) +
    geom_point(alpha = 1, aes(color = dis, size=dis)) +
    geom_text_repel(
      data = re_int,
      aes(x = ad_gene_sig, y = pv_gene_sig, label = sample),
      box.padding = 0.4
    ) + #add labels
    geom_errorbarh(aes(
      xmax = ad_gene_sig + ad_se,
      xmin = ad_gene_sig - ad_se,
      color = dis
    )) +
    geom_errorbar(aes(
      ymax = pv_gene_sig + pv_se,
      ymin = pv_gene_sig - pv_se,
      color = dis
    )) +
    ggtitle(paste(cell_types[i],"ind sample DEG filter")) +
    theme_classic() +
    xlab("AD-specific genes") +
    ylab("PV-specific genes") +
    theme(legend.position = "right") +
    scale_color_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "PRP" = "darkred",
        "AD" = "burlywood4",
        "PV" = "olivedrab4"
        
      )
    ) +
    scale_fill_manual(
      name = "Disease",
      values = c(
        "HC" = "darkblue",
        "PRP" = "darkred",
        "AD" = "burlywood4",
        "PV" = "olivedrab4"
      )
    ) +
    geom_mark_hull(data = re_int, aes(fill = dis, label = dis), concavity = 3) +
    theme(
      axis.text = element_text(size = 15, color = 'black'),
      axis.title = element_text(size = 15, color = 'black'),
      plot.title = element_text(size = 20, color = 'black')
    ) +
    scale_size_manual(
      name = "Disease",
      values = c(
        "HC" = 2,
        "PRP" = 4,
        "AD" = 2,
        "PV" = 2
      )
    )
  
  png(paste(deg_dir,i,"_",cell_types[i],"_hyperdim_plot_indsampledeg_filter.png"),
      width=8,height=8,units="in",res=400)
  print(p3)
  dev.off()
  
  #Plot all hyper-D plots together
  png(paste(deg_dir,i,"_",cell_types[i],"_all_hypderdim_plots.png"),
      width=24,height=8,units="in",res=400)
  print(multiplot(p1,p2,p3,cols=3))
  dev.off()
  
}, error=function(e){})
#-----------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------#
#Now write out sample specific DEGs to excel workbook

# OUT <- createWorkbook() #create excel workbook for a sample https://stackoverflow.com/questions/27713310/easy-way-to-export-multiple-data-frame-to-multiple-excel-worksheets
#PV no filter
sheet_name = paste(cell_types[i]," PV",sep="")
addWorksheet(OUT, sheet_name) # Add some sheets to the workbook
tryCatch({ writeData(OUT, sheet = sheet_name, x = deg_df_pv) }, error=function(e){}) #Write degs ncol(deg_df)

#AD no filter
sheet_name = paste(cell_types[i]," AD",sep="")
addWorksheet(OUT, sheet_name) # Add some sheets to the workbook
tryCatch({ writeData(OUT, sheet = sheet_name, x = deg_df_ad) }, error=function(e){}) #Write degs ncol(deg_df)

#AD-PV specific genes no filter
sheet_name = paste(cell_types[i]," AD-PV",sep="")
addWorksheet(OUT, sheet_name) # Add some sheets to the workbook
ad_genes_df = data.frame(ad_genes)
names(ad_genes_df) = "AD-specific genes"
pv_genes_df = data.frame(pv_genes)
names(pv_genes_df) = "PV-specific genes"
tryCatch({ writeData(OUT, sheet = sheet_name, x = ad_genes_df, startCol=1) }, error=function(e){})
tryCatch({ writeData(OUT, sheet = sheet_name, x = pv_genes_df, startCol=3) }, error=function(e){})

#AD-PV ind gene sample filter
sheet_name = paste(cell_types[i]," AD-PV indgene filt",sep="")
addWorksheet(OUT, sheet_name) # Add some sheets to the workbook
tryCatch({ ad_genes_df = data.frame(ad_genes_samplefiltered)}, error=function(e){})
tryCatch({ names(ad_genes_df) = "AD-specific genes"}, error=function(e){})
tryCatch({ pv_genes_df = data.frame(pv_genes_samplefiltered)}, error=function(e){})
tryCatch({ names(pv_genes_df) = "PV-specific genes"}, error=function(e){})
tryCatch({ writeData(OUT, sheet = sheet_name, x = ad_genes_df, startCol=1) }, error=function(e){})
tryCatch({ writeData(OUT, sheet = sheet_name, x = pv_genes_df, startCol=3) }, error=function(e){})

#AD ind sample DEG
sheet_name = paste(cell_types[i]," AD ind-samp DEG",sep="")
addWorksheet(OUT, sheet_name) # Add some sheets to the workbook
tryCatch({ writeData(OUT, sheet = sheet_name, x = ad_ind_gene_set_df) }, error=function(e){})

#PV ind sample DEG
sheet_name = paste(cell_types[i]," PV ind-samp DEG",sep="")
addWorksheet(OUT, sheet_name) # Add some sheets to the workbook
tryCatch({ writeData(OUT, sheet = sheet_name, x = pv_ind_gene_set_df) }, error=function(e){})

}, error=function(e){})#End cell type loop

}

#------------------------------------------------#
# Export the Excel file
#------------------------------------------------#
saveWorkbook(OUT, 
             paste(deg_dir, "cell_type_specific_deg_v3.xlsx",sep=""),
             overwrite = T)


scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/ad_pv_project/results/deg_rashx /Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/results
scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/ad_pv_project/results/deg_rashx/*all_hypderdim_plots.png /Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/results

#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#Find most common genes within each personalized DEG disease set

sheet_names = getSheetNames("/c4/home/mtaylor4/derm/ad_pv_project/results/deg_rashx/cell_type_specific_deg_v3.xlsx")

#Get personzlied gene signatures sheets
sheets_to_loop_over = sheet_names[grep("ind-samp",sheet_names)]

OUT <- createWorkbook()
do.call(rbind, lapply(1:length(sheets_to_loop_over), function(i){
  tryCatch({ 
    #i=20
    print(i)
    deg_sheet = read.xlsx(
      xlsxFile="/c4/home/mtaylor4/derm/ad_pv_project/results/deg_rashx/cell_type_specific_deg_v3.xlsx",
      sheet=sheets_to_loop_over[i])
    
    df_vector <- unlist(deg_sheet)
    
    # Use table to count occurrences of each unique element
    element_counts <- data.frame(table(df_vector))
    element_counts = element_counts[order(-element_counts$Freq),]
    names(element_counts) = c("gene","count")
    
    #Write sample to Excel work book
    sheet_name = sheets_to_loop_over[i]
    sheet_name = gsub(" ind-samp DEG","",sheet_name)
    addWorksheet(OUT, sheet_name) # Add some sheets to the workbook
    writeData(OUT, sheet = sheet_name, x = element_counts) 
    
  }, error=function(e){})
  
}))

saveWorkbook(OUT, 
             "/c4/home/mtaylor4/derm/ad_pv_project/results/deg_rashx/personalized_gene_count_summaries.xlsx",
             overwrite = T)

scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/ad_pv_project/results/deg_rashx/personalized_gene_count_summaries.xlsx /Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/results

#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#All samples DEG OF TRMS with higheset IL17A and IL17F

library(openxlsx)
library(Seurat)
library(fgsea)
library(topGO)
library(msigdbr)

#---------------------------------------------------------#
#NON-imputed: Plot genes of interest
#---------------------------------------------------------#

#Read in Seurat object that has combined all samples
suerat_obj = readRDS("/c4/home/mtaylor4/derm/PRP_case_report/all_samples_combined.RDS")

#Remove hand dermatitis
suerat_obj <- subset(x = suerat_obj, subset = dis != "Hand dermatitis")

#Within each disease find top IL17A and IL17F-producing Trms:

#Plot IL17A and IL17F for all Trm classes
unique(suerat_obj@meta.data$reference_celltype)
trms = unique(suerat_obj@meta.data$reference_celltype)[grep("Trm",unique(suerat_obj@meta.data$reference_celltype))]


#Get genes and cell-types to plot
to_plot = data.frame(suerat_obj@meta.data$dis, suerat_obj@meta.data$reference_celltype,
                     suerat_obj@assays$RNA@data["IL17A",],
                     suerat_obj@assays$RNA@data["IL17F",],
                     suerat_obj@assays$RNA@data["IFNG",],
                     suerat_obj@assays$RNA@data["IL26",])

names(to_plot) = c("dis","celltype","IL17A","IL17F","IFNG","IL26")


#get only trms
to_plot = to_plot[to_plot$celltype %in% trms,]

library(reshape2)
dfm = melt(to_plot, id=c("dis","celltype")) #Melt so can be plotted

dfm$dis = gsub("Healthy control","HC",dfm$dis)
dfm$dis = gsub("Psoriasis Vulgaris","PV",dfm$dis)

unimputed_trm_ILs_fig =  ggplot(dfm, aes(x=dis,y=value, fill=dis, color=dis)) +
  facet_grid(variable ~ celltype) +
  #stat_compare_means(comparisons = my_comparisons, method="wilcox.test")  +
  geom_violin(width=1, alpha=0.4) +
  geom_jitter(width=0.25, size=0.5, alpha=0.5) +
  geom_boxplot(outlier.size=0, width=0.3, color="black", size=0.7, alpha=0.4) +
  theme_classic() +
  ggtitle("Non-imputed") + 
  xlab("") + ylab("expression") + 
  theme(axis.text.x=element_text(size=20,color='black'),
        axis.text.y=element_text(size=15,color='black'),
        title=element_text(size=25),
        axis.ticks.x = element_line(size=0),
        strip.text = element_text(size = 20),
        legend.position="none") +
  scale_fill_manual(name="",values = c("HC" = "darkblue",
                                       "PRP" = "darkred",
                                       "PV" = "darkorange3")) +
  scale_color_manual(name="",values = c("HC" = "darkblue",
                                        "PRP" = "darkred",
                                        "PV" = "darkorange3"))

png("/Users/marktaylor/Desktop/Projects/Cho\ derm/Case\ Report\ PRP/results/IL_expression_figs/unimputed_ILs.png",
    width=15,height=13.5,units="in",res=400)
print(unimputed_trm_ILs_fig)
dev.off()


#---------------------------------------------------------#
#IMPUTATION: Impute on only top 100 genes bc otherwise runs out of memory: ARLA impmlemented by SeuratWrappers https://github.com/satijalab/seurat-wrappers/blob/master/docs/alra.md
#---------------------------------------------------------#
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(dplyr)

#Read in Seurat object that has combined all samples
suerat_obj = readRDS("/c4/home/mtaylor4/derm/PRP_case_report/all_samples_combined.RDS")

#Get variable genes
suerat_obj <- FindVariableFeatures(suerat_obj, selection.method = "vst", nfeatures = 1000)
var_genes = VariableFeatures(suerat_obj)

# Initial processing and visualization
suerat_obj <- SCTransform(suerat_obj) %>% RunPCA() %>% RunUMAP(dims = 1:30)

#Cobine variable genes with genes of interest
to_impute_genes = unique(c(var_genes, "IL17A","IL17F","IFNG","IL26"))

# run ALRA, creates alra assay of imputed values

suerat_obj <- RunALRA(suerat_obj,genes.use = to_impute_genes)

saveRDS(suerat_obj,file="/c4/home/mtaylor4/derm/PRP_case_report/all_samples_combined_imputed.RDS")


#---------------------------------------------------------#
#Imputed: Plot genes of interest
#---------------------------------------------------------#

#Read in Seurat object that has combined all samples
suerat_obj = readRDS("/c4/home/mtaylor4/derm/PRP_case_report/all_samples_combined_imputed.RDS")

#Remove hand dermatitis
suerat_obj <- subset(x = suerat_obj, subset = dis != "Hand dermatitis")

#Within each disease find top IL17A and IL17F-producing Trms:

#Plot IL17A and IL17F for all Trm classes
unique(suerat_obj@meta.data$reference_celltype)
trms = unique(suerat_obj@meta.data$reference_celltype)[grep("Trm",unique(suerat_obj@meta.data$reference_celltype))]

suerat_obj@assays$alra
#Get genes and cell-types to plot
to_plot = data.frame(suerat_obj@meta.data$dis, suerat_obj@meta.data$reference_celltype,
                     suerat_obj@assays$alra@data["IL17A",],
                     suerat_obj@assays$alra@data["IL17F",],
                     suerat_obj@assays$alra@data["IFNG",],
                     suerat_obj@assays$alra@data["IL26",])

names(to_plot) = c("dis","celltype","IL17A","IL17F","IFNG","IL26")


#get only trms
to_plot = to_plot[to_plot$celltype %in% trms,]

library(reshape2)
dfm = melt(to_plot, id=c("dis","celltype")) #Melt so can be plotted

dfm$dis = gsub("Healthy control","HC",dfm$dis)
dfm$dis = gsub("Psoriasis Vulgaris","PV",dfm$dis)

arla_mputed_trm_ILs_fig =  ggplot(dfm, aes(x=dis,y=value, fill=dis, color=dis)) +
  facet_grid(variable ~ celltype) +
  #stat_compare_means(comparisons = my_comparisons, method="wilcox.test")  +
  geom_violin(width=1, alpha=0.4) +
  geom_jitter(width=0.25, size=0.5, alpha=0.5) +
  geom_boxplot(outlier.size=0, width=0.3, color="black", size=0.7, alpha=0.4) +
  theme_classic() +
  ggtitle("ARLA-imputed") + 
  xlab("") + ylab("expression") + 
  theme(axis.text.x=element_text(size=20,color='black'),
        axis.text.y=element_text(size=15,color='black'),
        title=element_text(size=25),
        axis.ticks.x = element_line(size=0),
        strip.text = element_text(size = 20),
        legend.position="none") +
  scale_fill_manual(name="",values = c("HC" = "darkblue",
                                       "PRP" = "darkred",
                                       "PV" = "darkorange3")) +
  scale_color_manual(name="",values = c("HC" = "darkblue",
                                        "PRP" = "darkred",
                                        "PV" = "darkorange3"))

png("/c4/home/mtaylor4/derm/PRP_case_report/results/IL_expression_figs/arla_imputed_ILs.png",
    width=15,height=13.5,units="in",res=400)
print(arla_mputed_trm_ILs_fig)
dev.off()

scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/PRP_case_report/results/IL_expression_figs/* /Users/marktaylor/Desktop/Projects/Cho\ derm/Case\ Report\ PRP/results/IL_expression_figs

#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#DEG on ARLA imputed Trms
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#

#Set up GSEA stuff
all_gene_sets = msigdbr(species = "Homo sapiens")
pathwaysH = split(x = all_gene_sets$gene_symbol, f = all_gene_sets$gs_name) # fixing format to work with fgsea

#Function to return a list of GO and GSEA on a ranked vector of genes
go_gsea_fun = function(gene_rank_vector, geneUniverse) {
  
  #gene_rank_vector = as.vector(abs(deg_df_pv$avg_log2FC)) #make ranked and named vector to put into fgsea
  #names(gene_rank_vector) = row.names(deg_df_pv)
  #geneUniverse <- row.names(suerat_obj)
  
  #Get sig upregulated DEGs
  #--------------------------------------------#
  #GSEA 
  #--------------------------------------------#
  fgseaRes <- fgsea(pathways = pathwaysH, 
                    stats    = gene_rank_vector,
                    minSize  = 10,
                    maxSize  = 500,
                    scoreType="pos")
  ordered_gsea = fgseaRes[order(fgseaRes$pval),]
  
  #--------------------------------------------#
  #GO analysis on top 100 pos and neg DEGs
  #--------------------------------------------#
  #Top 100 POSITIVE DEGs
  genesOfInterest = names(gene_rank_vector)
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes=geneList,
                #geneSel = selection,
                annot = annFUN.org, 
                mapping = "org.Hs.eg.db",
                ID = "symbol")
  
  resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "fisher") #Fisher exact test
  tab <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 120)
  tab = subset(tab, raw.p.value <=0.01) #keep only sig results
  go_candidates = genesOfInterest
  selcGenes <- genesInTerm(GOdata, whichGO=tab$GO.ID)
  genes_in_go = do.call(rbind, lapply(1:length(selcGenes), function(s){
    #i=1
    data.frame(names(selcGenes)[s],  paste( c(na.omit(selcGenes[[s]][match(go_candidates,  selcGenes[[s]])])) , collapse=", ") )
  }))
  names(genes_in_go) = c("GO_ID","sig_genes_in_term")
  tab = cbind(tab, genes_in_go)
  
  return(list(ordered_gsea,tab))
}


#Make a new pan-Trm class: eg Trm1 --> Trm
suerat_obj@meta.data$reference_celltype_pantrm = suerat_obj@meta.data$reference_celltype
suerat_obj@meta.data$reference_celltype_pantrm[grep("Trm",suerat_obj@meta.data$reference_celltype_pantrm)] <- "Trm"
unique(suerat_obj@meta.data$reference_celltype_pantrm)

#Set up factor levels to loop through
quantile_thresholds = c(0.5,0.65,0.8,0.95)
celltypes_of_interest = c("Trm")
genes_of_interest = c("IL17A","IL17F","IFNG","IL26")
prp_samples = c("PRP 330_cd45","PRP 331_cd45")

#Reset default assay to RNA
DefaultAssay(suerat_obj) <- 'RNA'  #Make default assay RNA

#Variable genes, scale, PCA
suerat_obj <- FindVariableFeatures(suerat_obj, selection.method = "vst", nfeatures = 10000)

#Make master dir
master_dir = "/c4/home/mtaylor4/derm/PRP_case_report/results/arla_imputed_topcells/"
dir.create(master_dir)

#Make deg output dir
deg_dir = paste(master_dir,"deg/",sep="")
dir.create(deg_dir)

#Make text file dir
deg_textfile_fir = paste(deg_dir,"ind_files/",sep="")
dir.create(deg_textfile_fir)

#Get vectors of sample names
all_sample_names = unique(suerat_obj@meta.data$sample_id_dis)
hc_sample_names = all_sample_names[grep("Healthy control",all_sample_names)]
pv_sample_names = c( all_sample_names[grep("Psoriasis Vulgaris",all_sample_names)])

meta_data = suerat_obj@meta.data

OUT <- createWorkbook() #create excel workbook for a sample https://stackoverflow.com/questions/27713310/easy-way-to-export-multiple-data-frame-to-multiple-excel-worksheets

# for(i in 1:length(dis_sample_names)) {
for(i in 1:length(quantile_thresholds)) { #quantile loop
  tryCatch({ #suppress error
    #i=1
    
    for(k in 1:length(celltypes_of_interest)) { #cell type loop
      tryCatch({ #suppress error
        #k=1
        
        for(t in 1:length(genes_of_interest)) { #genes loop
          tryCatch({ #suppress error
            #t=1
            
            for(y in 1:length(prp_samples)) { #PRP samples loop
              tryCatch({ #suppress error
                #y=1
                
                print(paste("i=",i," k=",k," t=",t," y=",y,sep=""))
                
                #Now it's all about getting the right cell IDs
                
                #HC
                hc_cells = row.names(suerat_obj@meta.data[((suerat_obj@meta.data$sample_id_dis %in% hc_sample_names) & 
                                                             (suerat_obj@meta.data$reference_celltype_pantrm %in% celltypes_of_interest[k])),])
                
                #PRP
                prp_cells = row.names(suerat_obj@meta.data[((suerat_obj@meta.data$sample_id_dis %in% prp_samples[y]) & 
                                                              (suerat_obj@meta.data$reference_celltype_pantrm %in% celltypes_of_interest[k])),])
                
                #For a gene of interest get the cells expressing more than the quantile amount of ARLA-imputed value
                expression_threshold = quantile(suerat_obj@assays$alra[genes_of_interest[t],prp_cells], probs = quantile_thresholds[i], na.rm = T)
                expression_vector = t(data.frame(suerat_obj@assays$alra[genes_of_interest[t],prp_cells])) #get df of a gene's imputed values for cell subset
                expression_cut_off_cellnames = names(expression_vector[expression_vector[,1] >= expression_threshold, ]) #get only cells above a certain expression level
                
                #Replace decimals with dashes
                expression_cut_off_cellnames = gsub("\\.", "-", expression_cut_off_cellnames)
                
                deg_df <- FindMarkers(suerat_obj,
                                      ident.1 = expression_cut_off_cellnames,
                                      ident.2 = hc_cells)
                
                deg_df = deg_df[order(-deg_df$avg_log2FC),] #Re-order from high to low FC
                deg_df = cbind(row.names(deg_df),deg_df) #Add col of gene symbols
                names(deg_df)[1] = "gene"
                
                #GO and GSEA: go_gsea_fun(gene_rank_vector, geneUniverse)
                
                #Get  NEGATIVE SIGNIFICANT DEGs
                tryCatch({ #suppress error
                  deg_df_neg = subset(deg_df, avg_log2FC<0 & p_val_adj<0.05)
                  neg_vec = as.vector(abs(deg_df_neg$avg_log2FC)) #make ranked and named vector to put into fgsea
                  names(neg_vec) = row.names(deg_df_neg)
                  neg_enrichment = go_gsea_fun(neg_vec,row.names(suerat_obj))
                }, error=function(e){})
                
                #Get POSITIVE SIGNIFICANT DEGs
                tryCatch({ #suppress error
                  deg_df_neg = subset(deg_df, avg_log2FC>0 & p_val_adj<0.05)
                  neg_vec = as.vector(abs(deg_df_neg$avg_log2FC)) #make ranked and named vector to put into fgsea
                  names(neg_vec) = row.names(deg_df_neg)
                  pos_enrichment = go_gsea_fun(neg_vec,row.names(suerat_obj))
                }, error=function(e){})
                
                #Write out text files
                
                factor_level_name_combo = paste(prp_samples[y],celltypes_of_interest[k],genes_of_interest[t],quantile_thresholds[i])
                
                #Have to convert GSEA results to df bc it contains a list that won't write to text file like a goddam pain in the ass
                pos_gsea_df <- as.data.frame(pos_enrichment[[1]])
                neg_gsea_df <- as.data.frame(neg_enrichment[[1]])
                
                write.table(deg_df, file=paste(deg_textfile_fir,factor_level_name_combo," DEGs.txt",sep=""),sep="\t",row.names=F,quote=F) 
                tryCatch({ write.table(pos_enrichment[[2]], file=paste(deg_textfile_fir,factor_level_name_combo," POS SIG GO.txt",sep=""),sep="\t",row.names=F,quote=F) }, error=function(e){})
                tryCatch({ write.table(pos_gsea_df[,1:7], file=paste(deg_textfile_fir,factor_level_name_combo," POS SIG GSEA.txt",sep=""),sep="\t",row.names=F,quote=F) }, error=function(e){})
                tryCatch({ write.table(neg_enrichment[[2]], file=paste(deg_textfile_fir,factor_level_name_combo," NEG SIG GO.txt",sep=""),sep="\t",row.names=F,quote=F) }, error=function(e){})
                tryCatch({ write.table(neg_gsea_df[,1:7], file=paste(deg_textfile_fir,factor_level_name_combo," NEG SIG GSEA.txt",sep=""),sep="\t",row.names=F,quote=F) }, error=function(e){})
                
                #Now write out sample specific DEGs to excel workbook
                sheet_name = paste(factor_level_name_combo, sep="")
                addWorksheet(OUT, sheet_name) # Add some sheets to the workbook
                writeData(OUT, sheet = sheet_name, x = deg_df) #Write degs ncol(deg_df)
                tryCatch({ writeData(OUT, sheet = sheet_name, x = pos_enrichment[[2]], startCol=8) }, error=function(e){}) #Write GO res ncol(tab)
                tryCatch({ writeData(OUT, sheet = sheet_name, x = pos_gsea_df[,1:7], startCol=18) }, error=function(e){}) #Write GSEA res ncol(tab)
                tryCatch({ writeData(OUT, sheet = sheet_name, x = neg_enrichment[[2]], startCol=27) }, error=function(e){}) #Write GO res ncol(tab)
                tryCatch({ writeData(OUT, sheet = sheet_name, x = neg_gsea_df[,1:7], startCol=37) }, error=function(e){}) #Write GSEA res ncol(tab)
                
                
              }, error=function(e){}) }
          }, error=function(e){}) }
      }, error=function(e){}) }
  }, error=function(e){}) }


#------------------------------------------------#
# Export the Excel file
#------------------------------------------------#
saveWorkbook(OUT, 
             paste(deg_dir, "imputed_celltype_genecutoff_deg.xlsx",sep=""),
             overwrite = T)

scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/PRP_case_report/results/arla_imputed_topcells /Users/marktaylor/Desktop/Projects/Cho\ derm/Case\ Report\ PRP/results

#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#

#Now compare AD and PV signatures in Trm cells of the 2 PRP samples

pv_genes = c("CXCL13","CD7","MGAT4A","FTH1","LAYN","IL17F","KLRB1","GNLY","CPM","CTSH","GBP5","SOX4","CLEC2B","GZMB","CD2","CEBPD","ODF2L","LAG3","LRRN3","ARHGEF12","PTPN13","TNFAIP3","TRPS1","SNX9","METRNL","BTG1","JUN","SPOCK2","GABARAPL1","PMEPA1","RBPJ","LINC01871","MAP3K4","UBC","GALNT1","PNRC1","GABPB1-AS1","RPS26","MUC20-OT1","CHN1","NAP1L4","PTMS","F2R","CTLA4","DAPK2","RAP1B","CCR6","B3GALT2","YPEL2","FYN","PPDPF","SLA2","CBLB","ADGRG1","SARAF","DUSP1")
ad_genes = c("TWIST1","LGALS1","IL32","CAPG","ITM2C","MFHAS1","ANXA1","SOS1","CSGALNACT1","LMO4","IFITM2","S100A10","SYNE2","THADA","NEAT1","IL17RB","ARHGAP21","NBAS","ACTG1","TGFBR3","TNFSF10","AHNAK","ISG15","RPL17","LONRF2","TSHZ2","MMP25","IFITM1","BIRC3","FAM102A","LPCAT2","NRIP3","CRIP1","CLU","ZFP36","ZFP36L2","TUBA1B","GATA3","SLC5A3","SFXN1","FANK1","TAGLN2","RPS17","SMKR1","TC2N","MYL12A","LINC02195","CRNDE","CHDH","HLA-DRB1","SLC4A7","MIB1","SEC61G","RPS29","TRERF1","S100A4")

#Get only genes taht are in the object
pv_genes = pv_genes[pv_genes %in% row.names(suerat_obj)]

#Subset to PRP sample Trms
suerat_obj_sub <- subset(x = suerat_obj, subset = reference_celltype_pantrm == "Trm")
Idents(object = suerat_obj_sub) <- "sample_id_dis"
suerat_obj_sub <- subset(x = suerat_obj_sub, idents = c("PRP 330_cd45","PRP 331_cd45"))
unique(suerat_obj_sub$sample_id_dis)
unique(suerat_obj_sub$reference_celltype_pantrm)

signature_plot = cbind(data.frame(suerat_obj_sub$sample_id_dis), 
                       data.frame(rowSums(t(suerat_obj_sub@assays$RNA[ad_genes,]))),
                       data.frame(rowSums(t(suerat_obj_sub@assays$RNA[pv_genes,])))
)
names(signature_plot) = c("sample","AD_signature","PV_signature")

library(reshape2)
dfm = melt(signature_plot, id=c("sample")) #Melt so can be plotted

unique(dfm$sample)
dfm$sample = gsub("PRP 330_cd45","330",dfm$sample)
dfm$sample = gsub("PRP 331_cd45","331",dfm$sample)
dfm$variable = gsub("AD_signature","AD signature",dfm$variable)
dfm$variable = gsub("PV_signature","PV signature",dfm$variable)

my_comparisons <- list( c("330", "331"))

PRP_sample_AD_PV_signatures =  ggplot(dfm, aes(x=sample,y=value, fill=sample, color=sample)) +
  facet_grid(. ~ variable) +
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test")  +
  geom_violin(width=1, alpha=0.4) +
  geom_jitter(width=0.25, size=0.5, alpha=0.5) +
  geom_boxplot(outlier.size=0, width=0.3, color="black", size=0.7, alpha=0.4) +
  theme_classic() +
  ggtitle("PRP sample signature expression") + 
  xlab("sample") + ylab("integrated expression") + 
  theme(axis.text.x=element_text(size=20,color='black'),
        axis.text.y=element_text(size=15,color='black'),
        title=element_text(size=25),
        axis.ticks.x = element_line(size=0),
        strip.text = element_text(size = 20),
        legend.position="none") +
  scale_fill_manual(name="",values = c("330" = "goldenrod4",
                                       "331" = "orangered3")) +
  scale_color_manual(name="",values = c("330" = "goldenrod4",
                                        "331" = "orangered3"))

png("/c4/home/mtaylor4/derm/PRP_case_report/results/IL_expression_figs/arla_imputed_ILs.png",
    width=15,height=13.5,units="in",res=400)
print(arla_mputed_trm_ILs_fig)
dev.off()

#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#

#Now show that PRP cells cluster better with PV than HC

#Read in data
suerat_obj = readRDS("/c4/home/mtaylor4/derm/PRP_case_report/all_samples_combined_imputed.RDS")

#Remove hand dermatitis
suerat_obj <- subset(x = suerat_obj, subset = dis != "Hand dermatitis")

#Make a new pan-Trm class: eg Trm1 --> Trm
suerat_obj@meta.data$reference_celltype_pantrm = suerat_obj@meta.data$reference_celltype
suerat_obj@meta.data$reference_celltype_pantrm[grep("Trm",suerat_obj@meta.data$reference_celltype_pantrm)] <- "Trm"
unique(suerat_obj@meta.data$reference_celltype_pantrm)

#Get cell types to loop over
celltypes_of_interest = unique(suerat_obj$reference_celltype_pantrm)

clustering_proportion_results = do.call(rbind, lapply(1:length(celltypes_of_interest), function(k) {
  
  tryCatch({ #suppress error
    #k=26
    
    suerat_obj_sub <- subset(x = suerat_obj, subset = reference_celltype_pantrm == celltypes_of_interest[k])
    
    #Set active assay
    DefaultAssay(suerat_obj_sub) <- "RNA"
    
    suerat_obj_sub <- FindVariableFeatures(suerat_obj_sub, selection.method = "vst", nfeatures = 5000)
    all.genes <- rownames(suerat_obj_sub)
    suerat_obj_sub <- ScaleData(suerat_obj_sub, features = all.genes)
    suerat_obj_sub <- RunPCA(suerat_obj_sub, features = VariableFeatures(object = suerat_obj_sub))
    
    suerat_obj_sub <- FindNeighbors(suerat_obj_sub, dims = 1:10)
    suerat_obj_sub <- FindClusters(suerat_obj_sub, resolution = 0.5)
    suerat_obj_sub@meta.data$cell_cluster = suerat_obj_sub@active.ident #Add cell cluster id as metadata column
    
    #Run UMAP
    suerat_obj_sub <- RunUMAP(suerat_obj_sub, dims=1:10)
    
    #Now find proportion of PRP-PV cells in the same cluster
    n_prp = nrow( subset(suerat_obj_sub@meta.data, dis=="PRP") )
    n_pv = nrow( subset(suerat_obj_sub@meta.data, dis=="Psoriasis Vulgaris") )
    n_hc = nrow( subset(suerat_obj_sub@meta.data, dis=="Healthy control") )
    
    clusts = unique(suerat_obj_sub@meta.data$cell_cluster)
    
    co_occuring_clust_df = do.call(rbind, lapply(1:length(clusts), function(o){
      tryCatch({ #suppress error
        #o=1
        print(o)
        suerat_obj_sub_sub <- subset(x = suerat_obj_sub, subset = cell_cluster == clusts[o])
        prp_pv = round ( (nrow( suerat_obj_sub_sub@meta.data[suerat_obj_sub_sub@meta.data$dis %in% c("PRP","Psoriasis Vulgaris"),] ) / (nrow(suerat_obj_sub_sub@meta.data)) )*100, 2)
        prp_hc = round ( (nrow( suerat_obj_sub_sub@meta.data[suerat_obj_sub_sub@meta.data$dis %in% c("PRP","Healthy control"),] ) / (nrow(suerat_obj_sub_sub@meta.data)) )*100, 2)
        
        co_occuring_clust = data.frame(clusts[o], prp_pv, prp_hc)
        names(co_occuring_clust)[1] = "cluster"
        co_occuring_clust
      }, error=function(e){})
    })) #end cluster loop
    
    #combine with cell type and return
    celltype_dis_clust = data.frame(celltypes_of_interest[k],co_occuring_clust_df)
    names(celltype_dis_clust)[1] = "celltype"
    return(celltype_dis_clust)
    
  }, error=function(e){})
}))    

clustering_proportion_results   

#Now summarize and plot

#summary function
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length(na.omit(xx[[col]])),
                     mean = mean   (na.omit(xx[[col]])),
                     sd   = sd     (na.omit(xx[[col]]))
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

clustering_proportion_results = clustering_proportion_results[,-2]
dfm = melt(clustering_proportion_results, id=c("celltype")) #Melt so can be plotted

head(dfm)
dfc = summarySE(dfm, measurevar='value', groupvars=c('celltype','variable'))

dfc$variable = gsub("prp_pv","PRP-PV",dfc$variable)
dfc$variable = gsub("prp_hc","PRP-HC",dfc$variable)

clust_togeth_fig =  ggplot(dfc, aes(x=celltype,y=value, fill=variable, color=variable, alpha=variable)) +
  #facet_grid(. ~ variable) +
  #stat_compare_means(comparisons = my_comparisons, method="wilcox.test")  +
  geom_bar(stat = "identity", position="dodge") +
  theme_classic() +
  ggtitle("PRP unsupervised clustering with PV vs HC") + 
  xlab("") + ylab("% single cells clustered together") + 
  theme(axis.text.x=element_text(size=15,color='black',angle=90),
        axis.text.y=element_text(size=15,color='black'),
        title=element_text(size=25),
        axis.ticks.x = element_line(size=2),
        strip.text = element_text(size = 20),
        legend.position="bottom") +
  theme(legend.text=element_text(size=rel(3))) +
  scale_fill_manual(name="",values = c("PRP-PV" = "darkred",
                                       "PRP-HC" = "darkblue")) +
  scale_color_manual(name="",values = c("PRP-PV" = "darkred",
                                        "PRP-HC" = "darkblue")) +
  scale_alpha_manual(name="",values = c("PRP-PV" = 1,
                                        "PRP-HC" = 0.5))


#Trm alone
dfc_sub = dfc[dfc$celltype %in% c("Trm"),]

clust_togeth_fig =  ggplot(dfc_sub, aes(x=celltype,y=value, fill=variable, color=variable, alpha=variable)) +
  #facet_grid(. ~ variable) +
  #stat_compare_means(comparisons = my_comparisons, method="wilcox.test")  +
  geom_bar(stat = "identity", position="dodge") +
  theme_classic() +
  geom_errorbar(aes(x=celltype, ymin=value-se, ymax=value+se), width=0.2, colour="black", alpha=1, size=1.3, position=position_dodge(width=0.9)) +
  ggtitle("PRP unsupervised clustering") + 
  xlab("") + ylab("% single cells co-clustering") + 
  theme(axis.text.x=element_text(size=25,color='black',angle=0),
        axis.text.y=element_text(size=15,color='black'),
        title=element_text(size=25),
        axis.ticks.x = element_line(size=2),
        strip.text = element_text(size = 20),
        legend.position="bottom") +
  theme(legend.text=element_text(size=rel(3))) +
  scale_fill_manual(name="",values = c("PRP-PV" = "darkred",
                                       "PRP-HC" = "darkblue")) +
  scale_color_manual(name="",values = c("PRP-PV" = "darkred",
                                        "PRP-HC" = "darkblue")) +
  scale_alpha_manual(name="",values = c("PRP-PV" = 1,
                                        "PRP-HC" = 0.5))

#----------------------------------------------------------------------------------------------#
#Plot out specific cell type's unsupervised clustering + disease
#----------------------------------------------------------------------------------------------#
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

clustering_proportion_results = do.call(rbind, lapply(1:length(celltypes_of_interest), function(k) {
  
  tryCatch({ #suppress error
    #k=26
    
    suerat_obj_sub <- subset(x = suerat_obj, subset = reference_celltype_pantrm == celltypes_of_interest[k])
    
    #Set active assay
    DefaultAssay(suerat_obj_sub) <- "RNA"
    
    suerat_obj_sub <- FindVariableFeatures(suerat_obj_sub, selection.method = "vst", nfeatures = 5000)
    all.genes <- rownames(suerat_obj_sub)
    suerat_obj_sub <- ScaleData(suerat_obj_sub, features = all.genes)
    suerat_obj_sub <- RunPCA(suerat_obj_sub, features = VariableFeatures(object = suerat_obj_sub))
    
    suerat_obj_sub <- FindNeighbors(suerat_obj_sub, dims = 1:10)
    suerat_obj_sub <- FindClusters(suerat_obj_sub, resolution = 0.5)
    suerat_obj_sub@meta.data$cell_cluster = suerat_obj_sub@active.ident #Add cell cluster id as metadata column
    
    #Run UMAP
    suerat_obj_sub <- RunUMAP(suerat_obj_sub, dims=1:10)
    
    #Unsupervised cluster
    colourCount = length(unique(suerat_obj_sub@meta.data$cell_cluster))
    p1 <- DimPlot(suerat_obj_sub, reduction = "umap", group.by = "cell_cluster",
                  cols=getPalette(colourCount),
                  raster=FALSE) +
      ggtitle(paste(celltypes_of_interest[k],"unsupervised cell cluster")) + xlab("") + ylab("") +
      theme(axis.text=element_text(size=0,color='black'),
            axis.title=element_text(size=0),
            title=element_text(size=15),
            axis.ticks = element_line(size=0))
    
    #Disease
    unique(suerat_obj_sub@meta.data$dis)
    suerat_obj_sub@meta.data$dis = gsub("Healthy control","HC",suerat_obj_sub@meta.data$dis)
    suerat_obj_sub@meta.data$dis = gsub("Psoriasis Vulgaris","PV",suerat_obj_sub@meta.data$dis)
    
    p2 <- DimPlot(suerat_obj_sub, reduction = "umap", group.by = "dis",
                  #cols=getPalette(colourCount),
                  raster=FALSE) +
      ggtitle(paste(celltypes_of_interest[k],"disease")) + xlab("") + ylab("") +
      theme(axis.text=element_text(size=0,color='black'),
            axis.title=element_text(size=0),
            title=element_text(size=15),
            axis.ticks = element_line(size=0))  +
      scale_color_manual(name="",values = c("HC" = "darkblue",
                                            "PRP" = "darkred",
                                            "PV" = "darkorange3"))
    
    
    
    png(paste("/c4/home/mtaylor4/derm/PRP_case_report/results/cluster_pca_analysis/cluster_dig_fig/",
              k,"_",celltypes_of_interest[k],".png",sep=""),
        width=12,height=6,units="in",res=500)
    print(multiplot(p1,p2,cols=2))
    dev.off()
    
  }, error=function(e){})
}))    

#transfer results to my machine
scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/PRP_case_report/results/cluster_pca_analysis/cluster_dig_fig/ /Users/marktaylor/Desktop/Projects/Cho\ derm/Case\ Report\ PRP/results/cluster_pca_analysis


#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#Loop to plot through metadata metrics to a subset of seurat_object
#------------------------------------------------------------------------#
names(suerat_obj_sub@meta.data)
metadat_col_index = match(c("sample_id_dis"),names(suerat_obj_sub@meta.data))

metadat_df = suerat_obj_sub@meta.data #pull out all the cell metadata

#Find metadata column index of metrics that I want to color cells by 
umap_dir = paste("/c4/home/mtaylor4/derm/PRP_case_report/results/cluster_pca_analysis/cell_type_umap/",
                 "umap/",
                 sep="")
dir.create(umap_dir)

for(k in 1:length(metadat_col_index)) {
  tryCatch({ #suppress error
    
    #k=1
    colourCount = length(unique(suerat_obj@meta.data[,metadat_col_index[k]]))
    
    #Fig: all factor levels together
    p1 <- DimPlot(suerat_obj_sub, reduction = "umap", group.by = names(suerat_obj@meta.data)[metadat_col_index[k]],
                  cols=getPalette(colourCount),
                  raster=FALSE) +
      ggtitle( paste("Trm", names(suerat_obj_sub@meta.data)[metadat_col_index[k]],sep=" ") ) + xlab("") + ylab("") +
      theme(axis.text=element_text(size=0,color='black'),
            axis.title=element_text(size=0),
            title=element_text(size=15),
            axis.ticks = element_line(size=0))
    
    png(paste(umap_dir,
              k,".0_",names(suerat_obj_sub@meta.data)[metadat_col_index[k]] ,".png",sep=""),
        width=9,height=8,units="in",res=300)
    print(p1)
    dev.off()
    
    factor_levels = na.omit(unique(suerat_obj_sub@meta.data[,metadat_col_index[k]]))
    
    for(z in 1:length(factor_levels)){
      tryCatch({ #suppress error
        #z=1
        print(paste("k=",k," z=",z,sep=""))
        
        #Subset to factor level
        metadat_sub = metadat_df[metadat_df[,metadat_col_index[k]] %in% factor_levels[z],]
        
        p2 <- DimPlot(suerat_obj_sub, reduction = "umap",
                      cols="lightgrey",
                      cells.highlight = row.names(metadat_sub),
                      cols.highlight = "darkred",
                      sizes.highlight=0.3,
                      raster=FALSE) +
          ggtitle( paste("Trm"," " ,names(suerat_obj_sub@meta.data)[metadat_col_index[k]],": ", factor_levels[z],sep="")) + 
          xlab("") + ylab("") +
          theme(axis.text=element_text(size=0,color='black'),
                axis.title=element_text(size=0),
                title=element_text(size=15),
                axis.ticks = element_line(size=0),
                legend.position = "none")
        
        png(paste(umap_dir,
                  k,".",z,"_",names(suerat_obj_sub@meta.data)[metadat_col_index[k]],"_",factor_levels[z],".png",sep=""),
            width=7,height=7,units="in",res=300)
        print(p2)
        dev.off()
        
        
      }, error=function(e){})
      
    } #End factor level loop
    
  }, error=function(e){})
  
} # end metadata factor loop

#transfer results to my machine
scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/PRP_case_report/results/cluster_pca_analysis/cell_type_umap/ /Users/marktaylor/Desktop/Projects/Cho\ derm/Case\ Report\ PRP/results/cluster_pca_analysis


#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#GET SUPERCLUSTERS

suerat_obj_sub@meta.data$dis_supcluster = suerat_obj_sub@meta.data$cell_cluster
suerat_obj_sub@meta.data$dis_supcluster = gsub( "\\b0\\b", "PV-HC supercluster", suerat_obj_sub@meta.data$dis_supcluster,perl=TRUE)
suerat_obj_sub@meta.data$dis_supcluster = gsub( "\\b1\\b", "PV-HC supercluster", suerat_obj_sub@meta.data$dis_supcluster,perl=TRUE)
suerat_obj_sub@meta.data$dis_supcluster = gsub( "\\b2\\b", "PV-HC supercluster", suerat_obj_sub@meta.data$dis_supcluster,perl=TRUE)
suerat_obj_sub@meta.data$dis_supcluster = gsub( "\\b3\\b", "PRP-PV supercluster", suerat_obj_sub@meta.data$dis_supcluster,perl=TRUE)
suerat_obj_sub@meta.data$dis_supcluster = gsub( "\\b4\\b", "PV-HC supercluster", suerat_obj_sub@meta.data$dis_supcluster,perl=TRUE)
suerat_obj_sub@meta.data$dis_supcluster = gsub( "\\b5\\b", "PV-HC supercluster", suerat_obj_sub@meta.data$dis_supcluster,perl=TRUE)
suerat_obj_sub@meta.data$dis_supcluster = gsub( "\\b6\\b", "PV-HC supercluster", suerat_obj_sub@meta.data$dis_supcluster,perl=TRUE)
suerat_obj_sub@meta.data$dis_supcluster = gsub( "\\b7\\b", "PRP-PV supercluster", suerat_obj_sub@meta.data$dis_supcluster,perl=TRUE)
suerat_obj_sub@meta.data$dis_supcluster = gsub( "\\b8\\b", "PRP-PV supercluster", suerat_obj_sub@meta.data$dis_supcluster,perl=TRUE)
suerat_obj_sub@meta.data$dis_supcluster = gsub( "\\b9\\b", "PRP-PV supercluster", suerat_obj_sub@meta.data$dis_supcluster,perl=TRUE)
suerat_obj_sub@meta.data$dis_supcluster = gsub( "\\b10\\b", "PV-HC supercluster", suerat_obj_sub@meta.data$dis_supcluster,perl=TRUE)
suerat_obj_sub@meta.data$dis_supcluster = gsub( "\\b11\\b", "PRP-PV supercluster", suerat_obj_sub@meta.data$dis_supcluster,perl=TRUE)
suerat_obj_sub@meta.data$dis_supcluster = gsub( "\\b12\\b", "PV-HC supercluster", suerat_obj_sub@meta.data$dis_supcluster,perl=TRUE)
suerat_obj_sub@meta.data$dis_supcluster = gsub( "\\b13\\b", "PV-HC supercluster", suerat_obj_sub@meta.data$dis_supcluster,perl=TRUE)
suerat_obj_sub@meta.data$dis_supcluster = gsub( "\\b14\\b", "PRP-PV supercluster", suerat_obj_sub@meta.data$dis_supcluster,perl=TRUE)
suerat_obj_sub@meta.data$dis_supcluster = gsub( "\\b15\\b", "Disjunct", suerat_obj_sub@meta.data$dis_supcluster,perl=TRUE)
suerat_obj_sub@meta.data$dis_supcluster = gsub( "\\b16\\b", "PRP-PV supercluster", suerat_obj_sub@meta.data$dis_supcluster,perl=TRUE)

#Unsupervised cluster
colourCount = length(unique(suerat_obj_sub@meta.data$dis_supcluster))
p1 <- DimPlot(suerat_obj_sub, reduction = "umap", group.by = "dis_supcluster",
              cols=getPalette(colourCount),
              raster=FALSE) +
  ggtitle("Supercluster") + xlab("") + ylab("") +
  theme(axis.text=element_text(size=0,color='black'),
        axis.title=element_text(size=0),
        title=element_text(size=15),
        axis.ticks = element_line(size=0))


#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#DEG between PVs in supercluster 
#Set up GSEA stuff
all_gene_sets = msigdbr(species = "Homo sapiens")
pathwaysH = split(x = all_gene_sets$gene_symbol, f = all_gene_sets$gs_name) # fixing format to work with fgsea

#Function to return a list of GO and GSEA on a ranked vector of genes
go_gsea_fun = function(gene_rank_vector, geneUniverse) {
  
  #gene_rank_vector = as.vector(abs(deg_df_pv$avg_log2FC)) #make ranked and named vector to put into fgsea
  #names(gene_rank_vector) = row.names(deg_df_pv)
  #geneUniverse <- row.names(suerat_obj)
  
  #Get sig upregulated DEGs
  #--------------------------------------------#
  #GSEA 
  #--------------------------------------------#
  fgseaRes <- fgsea(pathways = pathwaysH, 
                    stats    = gene_rank_vector,
                    minSize  = 10,
                    maxSize  = 500,
                    scoreType="pos")
  ordered_gsea = fgseaRes[order(fgseaRes$pval),]
  
  #--------------------------------------------#
  #GO analysis on top 100 pos and neg DEGs
  #--------------------------------------------#
  #Top 100 POSITIVE DEGs
  genesOfInterest = names(gene_rank_vector)
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes=geneList,
                #geneSel = selection,
                annot = annFUN.org, 
                mapping = "org.Hs.eg.db",
                ID = "symbol")
  
  resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "fisher") #Fisher exact test
  tab <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 120)
  tab = subset(tab, raw.p.value <=0.01) #keep only sig results
  go_candidates = genesOfInterest
  selcGenes <- genesInTerm(GOdata, whichGO=tab$GO.ID)
  genes_in_go = do.call(rbind, lapply(1:length(selcGenes), function(s){
    #i=1
    data.frame(names(selcGenes)[s],  paste( c(na.omit(selcGenes[[s]][match(go_candidates,  selcGenes[[s]])])) , collapse=", ") )
  }))
  names(genes_in_go) = c("GO_ID","sig_genes_in_term")
  tab = cbind(tab, genes_in_go)
  
  return(list(ordered_gsea,tab))
}

#Reset default assay to RNA
DefaultAssay(suerat_obj_sub) <- 'RNA'  #Make default assay RNA

#Variable genes, scale, PCA
suerat_obj_sub <- FindVariableFeatures(suerat_obj_sub, selection.method = "vst", nfeatures = 10000)

#Make master dir
master_dir = "/c4/home/mtaylor4/derm/PRP_case_report/results/cluster_pca_analysis/supercluster_deg/"
dir.create(master_dir)

#Make deg output dir
deg_dir = paste(master_dir,"deg/",sep="")
dir.create(deg_dir)

#Make text file dir
deg_textfile_fir = paste(deg_dir,"ind_files/",sep="")
dir.create(deg_textfile_fir)


suerat_obj_sub@meta.data

pvcells_supercluster_prp_pv = row.names(suerat_obj_sub@meta.data[suerat_obj_sub@meta.data$dis_supcluster %in% "PRP-PV supercluster" & 
                                                                   suerat_obj_sub@meta.data$dis %in% "PV",])
pvcells_supercluster_pv_hc = row.names(suerat_obj_sub@meta.data[suerat_obj_sub@meta.data$dis_supcluster %in% "PV-HC supercluster" & 
                                                                  suerat_obj_sub@meta.data$dis %in% "PV",])

deg_df <- FindMarkers(suerat_obj_sub,
                      ident.1 = pvcells_supercluster_prp_pv,
                      ident.2 = pvcells_supercluster_pv_hc)

deg_df = deg_df[order(-deg_df$avg_log2FC),] #Re-order from high to low FC
deg_df = cbind(row.names(deg_df),deg_df) #Add col of gene symbols
names(deg_df)[1] = "gene"

#GO and GSEA: go_gsea_fun(gene_rank_vector, geneUniverse)

#Get  NEGATIVE SIGNIFICANT DEGs
tryCatch({ #suppress error
  deg_df_neg = subset(deg_df, avg_log2FC<0 & p_val_adj<0.05)
  neg_vec = as.vector(abs(deg_df_neg$avg_log2FC)) #make ranked and named vector to put into fgsea
  names(neg_vec) = row.names(deg_df_neg)
  neg_enrichment = go_gsea_fun(neg_vec,row.names(suerat_obj))
}, error=function(e){})

#Get POSITIVE SIGNIFICANT DEGs
tryCatch({ #suppress error
  deg_df_neg = subset(deg_df, avg_log2FC>0 & p_val_adj<0.05)
  neg_vec = as.vector(abs(deg_df_neg$avg_log2FC)) #make ranked and named vector to put into fgsea
  names(neg_vec) = row.names(deg_df_neg)
  pos_enrichment = go_gsea_fun(neg_vec,row.names(suerat_obj))
}, error=function(e){})

#Write out text files

factor_level_name_combo = "Trm PRP-PV PV-HC superclust"

#Have to convert GSEA results to df bc it contains a list that won't write to text file like a goddam pain in the ass
pos_gsea_df <- as.data.frame(pos_enrichment[[1]])
neg_gsea_df <- as.data.frame(neg_enrichment[[1]])

write.table(deg_df, file=paste(deg_textfile_fir,factor_level_name_combo," DEGs.txt",sep=""),sep="\t",row.names=F,quote=F) 
tryCatch({ write.table(pos_enrichment[[2]], file=paste(deg_textfile_fir,factor_level_name_combo," POS SIG GO.txt",sep=""),sep="\t",row.names=F,quote=F) }, error=function(e){})
tryCatch({ write.table(pos_gsea_df[,1:7], file=paste(deg_textfile_fir,factor_level_name_combo," POS SIG GSEA.txt",sep=""),sep="\t",row.names=F,quote=F) }, error=function(e){})
tryCatch({ write.table(neg_enrichment[[2]], file=paste(deg_textfile_fir,factor_level_name_combo," NEG SIG GO.txt",sep=""),sep="\t",row.names=F,quote=F) }, error=function(e){})
tryCatch({ write.table(neg_gsea_df[,1:7], file=paste(deg_textfile_fir,factor_level_name_combo," NEG SIG GSEA.txt",sep=""),sep="\t",row.names=F,quote=F) }, error=function(e){})

OUT <- createWorkbook() #create excel workbook for a sample https://stackoverflow.com/questions/27713310/easy-way-to-export-multiple-data-frame-to-multiple-excel-worksheets

#Now write out sample specific DEGs to excel workbook
sheet_name = paste(factor_level_name_combo, sep="")
addWorksheet(OUT, sheet_name) # Add some sheets to the workbook
writeData(OUT, sheet = sheet_name, x = deg_df) #Write degs ncol(deg_df)
tryCatch({ writeData(OUT, sheet = sheet_name, x = pos_enrichment[[2]], startCol=8) }, error=function(e){}) #Write GO res ncol(tab)
tryCatch({ writeData(OUT, sheet = sheet_name, x = pos_gsea_df[,1:7], startCol=18) }, error=function(e){}) #Write GSEA res ncol(tab)
tryCatch({ writeData(OUT, sheet = sheet_name, x = neg_enrichment[[2]], startCol=27) }, error=function(e){}) #Write GO res ncol(tab)
tryCatch({ writeData(OUT, sheet = sheet_name, x = neg_gsea_df[,1:7], startCol=37) }, error=function(e){}) #Write GSEA res ncol(tab)


#------------------------------------------------#
# Export the Excel file
#------------------------------------------------#
saveWorkbook(OUT, 
             paste(deg_dir, "trm_supercluster_deg.xlsx",sep=""),
             overwrite = T)

scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/PRP_case_report/results/cluster_pca_analysis/supercluster_deg/ /Users/marktaylor/Desktop/Projects/Cho\ derm/Case\ Report\ PRP/results

#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#
#------------------------------------------------------------------------#


png("/c4/home/mtaylor4/derm/PRP_case_report/results/IL_expression_figs/arla_imputed_ILs.png",
    width=15,height=13.5,units="in",res=400)
print(arla_mputed_trm_ILs_fig)
dev.off()

#------------------------------------------------------------------------#
#Loop to plot through metadata metrics
#------------------------------------------------------------------------#
names(suerat_obj@meta.data)
metadat_col_index = match(c("cell_cluster","dis_level_1","dis_level_2",
                            "sample.id","project","data_ext_int"),names(suerat_obj@meta.data))

metadat_df = suerat_obj@meta.data #pull out all the cell metadata

#Find metadata column index of metrics that I want to color cells by 
umap_dir = paste(master_dir,
                 "umap/",
                 sep="")
dir.create(umap_dir)

for(k in 1:length(metadat_col_index)) {
  tryCatch({ #suppress error
    
    #k=2
    colourCount = length(unique(suerat_obj@meta.data[,metadat_col_index[k]]))
    
    #Fig: all factor levels together
    p1 <- DimPlot(suerat_obj, reduction = "umap", group.by = names(suerat_obj@meta.data)[metadat_col_index[k]],
                  cols=getPalette(colourCount),
                  raster=FALSE) +
      ggtitle( paste(celltype_file_shortname[i], names(suerat_obj@meta.data)[metadat_col_index[k]],sep=" ") ) + xlab("") + ylab("") +
      theme(axis.text=element_text(size=0,color='black'),
            axis.title=element_text(size=0),
            title=element_text(size=15),
            axis.ticks = element_line(size=0))
    
    png(paste(umap_dir,
              k,".0_",names(suerat_obj@meta.data)[metadat_col_index[k]] ,".png",sep=""),
        width=9,height=8,units="in",res=300)
    print(p1)
    dev.off()
    
    factor_levels = na.omit(unique(suerat_obj@meta.data[,metadat_col_index[k]]))
    
    for(z in 1:length(factor_levels)){
      tryCatch({ #suppress error
        #z=1
        print(paste("k=",k," z=",z,sep=""))
        
        #Subset to factor level
        metadat_sub = metadat_df[metadat_df[,metadat_col_index[k]] %in% factor_levels[z],]
        
        p2 <- DimPlot(suerat_obj, reduction = "umap",
                      cols="lightgrey",
                      cells.highlight = row.names(metadat_sub),
                      cols.highlight = "darkred",
                      sizes.highlight=0.3,
                      raster=FALSE) +
          ggtitle( paste(celltype_file_shortname[i]," " ,names(suerat_obj@meta.data)[metadat_col_index[k]],": ", factor_levels[z],sep="")) + 
          xlab("") + ylab("") +
          theme(axis.text=element_text(size=0,color='black'),
                axis.title=element_text(size=0),
                title=element_text(size=15),
                axis.ticks = element_line(size=0),
                legend.position = "none")
        
        png(paste(umap_dir,
                  k,".",z,"_",names(suerat_obj@meta.data)[metadat_col_index[k]],"_",factor_levels[z],".png",sep=""),
            width=7,height=7,units="in",res=300)
        print(p2)
        dev.off()
        
        
      }, error=function(e){})
      
    } #End factor level loop
    
  }, error=function(e){})
  
} # end metadata factor loop







#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#All samples within CELL TYPE DEGs

library(openxlsx)
library(Seurat)
library(fgsea)
library(topGO)
library(msigdbr)

#Set up GSEA stuff
all_gene_sets = msigdbr(species = "Homo sapiens")
pathwaysH = split(x = all_gene_sets$gene_symbol, f = all_gene_sets$gs_name) # fixing format to work with fgsea

#Function to return a list of GO and GSEA on a ranked vector of genes
go_gsea_fun = function(gene_rank_vector, geneUniverse) {
  
  #gene_rank_vector = as.vector(abs(deg_df_pv$avg_log2FC)) #make ranked and named vector to put into fgsea
  #names(gene_rank_vector) = row.names(deg_df_pv)
  #geneUniverse <- row.names(suerat_obj)
  
  #Get sig upregulated DEGs
  #--------------------------------------------#
  #GSEA 
  #--------------------------------------------#
  fgseaRes <- fgsea(pathways = pathwaysH, 
                    stats    = gene_rank_vector,
                    minSize  = 10,
                    maxSize  = 500,
                    scoreType="pos")
  ordered_gsea = fgseaRes[order(fgseaRes$pval),]
  
  #--------------------------------------------#
  #GO analysis on top 100 pos and neg DEGs
  #--------------------------------------------#
  #Top 100 POSITIVE DEGs
  genesOfInterest = names(gene_rank_vector)
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes=geneList,
                #geneSel = selection,
                annot = annFUN.org, 
                mapping = "org.Hs.eg.db",
                ID = "symbol")
  
  resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "fisher") #Fisher exact test
  tab <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 120)
  tab = subset(tab, raw.p.value <=0.01) #keep only sig results
  go_candidates = genesOfInterest
  selcGenes <- genesInTerm(GOdata, whichGO=tab$GO.ID)
  genes_in_go = do.call(rbind, lapply(1:length(selcGenes), function(s){
    #i=1
    data.frame(names(selcGenes)[s],  paste( c(na.omit(selcGenes[[s]][match(go_candidates,  selcGenes[[s]])])) , collapse=", ") )
  }))
  names(genes_in_go) = c("GO_ID","sig_genes_in_term")
  tab = cbind(tab, genes_in_go)
  
  return(list(ordered_gsea,tab))
}


#Read in Seurat object that has combined all samples
suerat_obj = readRDS("/c4/home/mtaylor4/derm/PRP_case_report/all_samples_combined.RDS")

#Remove hand dermatitis
suerat_obj <- subset(x = suerat_obj, subset = dis != "Hand dermatitis")

#Get table of cell counts vs disease
table(suerat_obj@meta.data$reference_celltype,suerat_obj@meta.data$dis)

#Variable genes, scale, PCA
suerat_obj <- FindVariableFeatures(suerat_obj, selection.method = "vst", nfeatures = 10000)
all.genes <- rownames(suerat_obj)
suerat_obj <- ScaleData(suerat_obj, features = all.genes)

#Make master dir
master_dir = "/c4/home/mtaylor4/derm/PRP_case_report/results/all_sample_merged/deg/"

#Make deg output dir
deg_dir = paste(master_dir,
                "celltype_specific/",
                sep="")
dir.create(deg_dir)

#Make text file dir
deg_textfile_fir = paste(deg_dir,
                         "ind_files/",
                         sep="")
dir.create(deg_textfile_fir)

#Make vector of all cell types
cell_types = unique(suerat_obj@meta.data$reference_celltype)

#Initiate excel workbook
OUT <- createWorkbook() #create excel workbook for a sample https://stackoverflow.com/questions/27713310/easy-way-to-export-multiple-data-frame-to-multiple-excel-worksheets

#Start loop through cell types
for(i in 1:length(cell_types)) {
  # for(i in 1:2) {
  
  tryCatch({ #suppress error
    
    #i=1
    
    suerat_obj_sub <- subset(x = suerat_obj, subset = reference_celltype == cell_types[i])
    diseases = unique(suerat_obj_sub@meta.data$dis)
    diseases =  na.omit(c( diseases[grep("Psoriasis Vulgaris",diseases)],diseases[grep("PRP",diseases)] ))
    
    #Start loop through diseases within each cell type
    for(j in 1:length(diseases)) {
      
      tryCatch({ #suppress error
        
        #j=1    
        deg_df <- FindMarkers(suerat_obj_sub,
                              ident.1 = row.names(suerat_obj_sub@meta.data[suerat_obj_sub@meta.data$dis %in% diseases[j],]),
                              ident.2 = row.names(suerat_obj_sub@meta.data[suerat_obj_sub@meta.data$dis %in% "Healthy control",])
        )
        
        deg_df = deg_df[order(-deg_df$avg_log2FC),] #Re-order from high to low FC
        deg_df = cbind(row.names(deg_df),deg_df) #Add col of gene symbols
        names(deg_df)[1] = "gene"
        
        #GO and GSEA: go_gsea_fun(gene_rank_vector, geneUniverse)
        
        #Get  NEGATIVE SIGNIFICANT DEGs
        tryCatch({ #suppress error
          deg_df_neg = subset(deg_df, avg_log2FC<0 & p_val_adj<0.05)
          neg_vec = as.vector(abs(deg_df_neg$avg_log2FC)) #make ranked and named vector to put into fgsea
          names(neg_vec) = row.names(deg_df_neg)
          neg_enrichment = go_gsea_fun(neg_vec,row.names(suerat_obj))
        }, error=function(e){})
        
        #Get POSITIVE SIGNIFICANT DEGs
        tryCatch({ #suppress error
          deg_df_neg = subset(deg_df, avg_log2FC>0 & p_val_adj<0.05)
          neg_vec = as.vector(abs(deg_df_neg$avg_log2FC)) #make ranked and named vector to put into fgsea
          names(neg_vec) = row.names(deg_df_neg)
          pos_enrichment = go_gsea_fun(neg_vec,row.names(suerat_obj))
        }, error=function(e){})
        
        #Write out text files
        
        write.table(deg_df, file=paste(cell_types[i],diseases[j]," DEGs.txt",sep=""),sep="\t",row.names=F,quote=F) 
        tryCatch({ write.table(pos_enrichment[[2]], file=paste(deg_textfile_fir,cell_types[i],diseases[j],"POS SIG GO.txt",sep=" "),sep="\t",row.names=F,quote=F) }, error=function(e){})
        tryCatch({ write.table(pos_enrichment[[1]], file=paste(deg_textfile_fir,cell_types[i],diseases[j],"POS SIG GSEA.txt",sep=" "),sep="\t",row.names=F,quote=F) }, error=function(e){})
        tryCatch({ write.table(neg_enrichment[[2]], file=paste(deg_textfile_fir,cell_types[i],diseases[j],"NEG SIG GO.txt",sep=" "),sep="\t",row.names=F,quote=F) }, error=function(e){})
        tryCatch({ write.table(neg_enrichment[[1]], file=paste(deg_textfile_fir,cell_types[i],diseases[j],"NEG SIG GSEA.txt",sep=" "),sep="\t",row.names=F,quote=F) }, error=function(e){})
        
        
        #Now write out sample specific DEGs to excel workbook
        sheet_name = paste(cell_types[i],diseases[j], sep=" ")
        addWorksheet(OUT, sheet_name) # Add some sheets to the workbook
        writeData(OUT, sheet = sheet_name, x = deg_df) #Write degs ncol(deg_df)
        tryCatch({ writeData(OUT, sheet = sheet_name, x = pos_enrichment[[2]], startCol=8) }, error=function(e){}) #Write GO res ncol(tab)
        tryCatch({ writeData(OUT, sheet = sheet_name, x = pos_enrichment[[1]], startCol=18) }, error=function(e){}) #Write GSEA res ncol(tab)
        tryCatch({ writeData(OUT, sheet = sheet_name, x = neg_enrichment[[2]], startCol=27) }, error=function(e){}) #Write GO res ncol(tab)
        tryCatch({ writeData(OUT, sheet = sheet_name, x = neg_enrichment[[1]], startCol=37) }, error=function(e){}) #Write GSEA res ncol(tab)
        
      }, error=function(e){})
    } #End disease-celltype loop
    
  }, error=function(e){}) 
} #End cell type loop
#------------------------------------------------#
# Export the Excel file
#------------------------------------------------#
saveWorkbook(OUT, 
             paste(deg_dir, "celltype_dis_deg.xlsx",sep=""),
             overwrite = T)

scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/PRP_case_report/results/all_sample_merged /Users/marktaylor/Desktop/Projects/Cho\ derm/Case\ Report\ PRP/results

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################

out_dirs = list.dirs("/c4/home/mtaylor4/derm/internal_data/cellranger_output",recursive=F)
out_dir_names = list.dirs("/c4/home/mtaylor4/derm/internal_data/cellranger_output",recursive=F,full.names = F)

mapped_files = do.call(rbind, lapply(1:length(out_dirs), function(i){
  #i=1
  print(i)
  files_present = list.files(paste(out_dirs[i],"/outs/filtered_feature_bc_matrix/",sep=""),full.names = F,pattern=".gz")
  res_file = data.frame(out_dir_names[i], length(files_present), paste(files_present,collapse=" "))
  names(res_file) = c("sample","filtered_files","filtered_files_names")
  res_file
}))


#Combine cell types from across
library(Seurat)
library(future)

rds_files = list.files("/c4/home/mtaylor4/derm/PRP_case_report/seurat_objects_cell_typed_combined_reference",full.names = T) #get all RDS objects in this directory
sampe_names = as.data.frame(t(data.frame(strsplit(rds_files,"/c4/home/mtaylor4/derm/PRP_case_report/seurat_objects_cell_typed_combined_reference"))))
sampe_names = gsub(".RDS","", sampe_names$V2)
sampe_names = gsub("/","", sampe_names)

#Specify cell types
cell_types = c("CTLem","eTreg1","CTLex","Tmm1","ILC/NK","NK","Trm3","Trm1","Tcm","Tmm2","Tet","Trm2","cmTreg","Tn",
               "Tmm3","eTreg2","ILC2","moDC1","DC1","Mac2","Mono","moDC2","B","InfMono","LC1","moDC3",
               "LC2","Mac1","LC3","DC2","DC3","Plasma","Mac3","Mac4","migDC","Mast","Mast-c","Trm-c","Treg-c",
               "CTL-c","ILC/NK-c",
               "Keratinocyte diff","Keratinocyte undiff","Fibroblast","Pericyte","Vascular EC","Melanocyte","Lymphatic EC","Erythrocyte")

#specify the minimum number cells that have to be present in a sample
min_sample_cells = 30

#Loop through all cell types in internal and external data and integrate into a single object
for(j in 1:length(cell_types)){
  
  tryCatch({ #suppress error
    #j=1
    rds_list = do.call(list, lapply(1:length(rds_files), function(i){
      # rds_list = do.call(list, lapply(c(1,3:10), function(i){
      
      tryCatch({ #suppress error
        #i=1
        print(i)
        suerat_obj = readRDS(rds_files[i])
        suerat_obj = subset(x = suerat_obj, subset = reference_celltype == cell_types[j]) #Subset to only a specific cell type
        DefaultAssay(suerat_obj) <- 'RNA'  #Make default assay RNA
        
        #assign each cds object to a new name
        name <- paste("rds.", i, sep = "")
        return(assign(name, suerat_obj)) #Create a new object name that points to the Suerat object that was just read in
      }, error=function(e){})
    }))
    
    #Remove null list elements (samples with no cells)
    remove_rds_vec = which(sapply(rds_list, is.null))
    if(length(remove_rds_vec)==0) {rds_list = rds_list} else {rds_list = rds_list[-which(sapply(rds_list, is.null))]}
    
    
    #-----------------------------------------#
    #remove samples with too few cells
    #-----------------------------------------#
    #Find indices of samples in list with too few cells
    small_samples_indices = do.call(c, lapply(1:length(rds_list), function(k) {
      #k=1
      if ( dim(rds_list[[k]])[2] <= min_sample_cells ) {small_samples_index = k } else {small_samples_index = ""}
      small_samples_index
    }))
    samples_to_drop = na.omit(as.numeric(small_samples_indices))
    rds_list = rds_list[-c(samples_to_drop)] #Remove Seurat objects with fewer cells than min_sample_cells
    
    #Seurat integration: https://satijalab.org/seurat/articles/integration_introduction.html
    features <- SelectIntegrationFeatures(object.list = rds_list)
    
    #Integration anchors
    #Speed up FindInetegrationAnchors: https://github.com/satijalab/seurat/discussions/3999
    plan("multisession", workers = 4) # Enable parallelization from future package
    options(future.globals.maxSize = 20000 * 1024^2)
    cell.anchors <- FindIntegrationAnchors(object.list = rds_list, 
                                           anchor.features = features,
                                           #reference = c(1,2),
                                           reduction = "rpca",
                                           dims = 1:30
    )
    
    #Find smallest number of anchor pairs in order to set k.weight smaller downstream
    cell_anch_pairs = paste(cell.anchors@anchors$dataset1,cell.anchors@anchors$dataset2)
    cell_anch_pairs_freq = data.frame(table(cell_anch_pairs))
    cell_anch_pairs_freq = cell_anch_pairs_freq[order(cell_anch_pairs_freq$Freq),]
    min_anchors = cell_anch_pairs_freq$Freq[1]
    
    #Now actually integrate the data
    plan("sequential") # Stop paralellization
    immune.combined <- IntegrateData(anchorset = cell.anchors,
                                     k.weight = 10) # https://github.com/satijalab/seurat/issues/3930
    
    saveRDS(immune.combined,
            file=paste("/c4/home/mtaylor4/derm/PRP_case_report/seurat_objects_celltype_integrated/",
                       j,"_",cell_types[j],
                       ".RDS",
                       sep=""))
    
  }, error=function(e){})
} #end cell type loop incrementing on j

#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------#

#Transfer dis markers here

scp /Users/marktaylor/Desktop/Projects/Cho\ derm/Case\ Report\ PRP/dis_marker_sets/dis_markers.xlsx mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/PRP_case_report/dis_markers

#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#Unsupervised cell cluster, DEG, and PCA analysis
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#

# BiocManager::install("org.Hs.eg.db")

#Loop through all cell type RDS's and plot out 
library(openxlsx)
library(Seurat)
library(ggplot2)
library(fgsea)
library(topGO)
library(msigdbr)
library(RColorBrewer)
library(ggrepel)

set.seed(123)

getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#summary function
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length(na.omit(xx[[col]])),
                     mean = mean   (na.omit(xx[[col]])),
                     sd   = sd     (na.omit(xx[[col]]))
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#---------------------------------------------------------#
#---------------------------------------------------------#
#---------------------------------------------------------#
#---------------------------------------------------------#

version_dir = "v1"

version_dir_tomake = paste("/c4/home/mtaylor4/derm/PRP_case_report/results/cluster_pca_analysis/",
                           version_dir,"/",
                           sep="")
dir.create(version_dir_tomake)

celltype_files = list.files("/c4/home/mtaylor4/derm/PRP_case_report/seurat_objects_celltype_integrated/",full.names = T)
celltype_file_shortname = list.files("/c4/home/mtaylor4/derm/PRP_case_report/seurat_objects_celltype_integrated/",full.names = F)
celltype_file_shortname = unlist(strsplit(celltype_file_shortname,".RDS"))

#Sample ID and disease matches
dis_markers = read.xlsx(
  "/c4/home/mtaylor4/derm/PRP_case_report/dis_markers/dis_markers.xlsx"
)

#Set up GSEA stuff
all_gene_sets = msigdbr(species = "Homo sapiens")
pathwaysH = split(x = all_gene_sets$gene_symbol, f = all_gene_sets$gs_name) # fixing format to work with fgsea

# Set up integrated a priori gene sets to test DE in aggregate gene set expression
gene_sets <- msigdbr(species = "Homo sapiens")
gene_sets <- gene_sets[gene_sets$gs_cat %in% "C7" & gene_sets$gs_subcat %in% "IMMUNESIGDB",]

#Collect all gene sets
gene_set_cats = unique(gene_sets$gs_name)

#Function to return a list of GO and GSEA on a ranked vector of genes
go_gsea_fun = function(gene_rank_vector, geneUniverse) {
  
  #gene_rank_vector = as.vector(abs(deg_df_pv$avg_log2FC)) #make ranked and named vector to put into fgsea
  #names(gene_rank_vector) = row.names(deg_df_pv)
  #geneUniverse <- row.names(suerat_obj)
  
  #Get sig upregulated DEGs
  #--------------------------------------------#
  #GSEA 
  #--------------------------------------------#
  fgseaRes <- fgsea(pathways = pathwaysH, 
                    stats    = gene_rank_vector,
                    minSize  = 10,
                    maxSize  = 500,
                    scoreType="pos")
  ordered_gsea = fgseaRes[order(fgseaRes$pval),]
  
  #--------------------------------------------#
  #GO analysis on top 100 pos and neg DEGs
  #--------------------------------------------#
  #Top 100 POSITIVE DEGs
  genesOfInterest = names(gene_rank_vector)
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes=geneList,
                #geneSel = selection,
                annot = annFUN.org, 
                mapping = "org.Hs.eg.db",
                ID = "symbol")
  
  resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "fisher") #Fisher exact test
  tab <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 120)
  tab = subset(tab, raw.p.value <=0.01) #keep only sig results
  go_candidates = genesOfInterest
  selcGenes <- genesInTerm(GOdata, whichGO=tab$GO.ID)
  genes_in_go = do.call(rbind, lapply(1:length(selcGenes), function(s){
    #i=1
    data.frame(names(selcGenes)[s],  paste( c(na.omit(selcGenes[[s]][match(go_candidates,  selcGenes[[s]])])) , collapse=", ") )
  }))
  names(genes_in_go) = c("GO_ID","sig_genes_in_term")
  tab = cbind(tab, genes_in_go)
  
  return(list(ordered_gsea,tab))
}

#---------------------------------------------------------#
#---------------------------------------------------------#
#---------------------------------------------------------#

#Start loop!!!!!!!
for(i in 1:length(celltype_files)){
  tryCatch({ #suppress error
    
    #i=1
    print(i)
    
    #Output dir
    #Make output directory
    master_dir = paste("/c4/home/mtaylor4/derm/PRP_case_report/results/cluster_pca_analysis/",
                       version_dir,"/",
                       celltype_file_shortname[i],"/",sep="")
    dir.create(master_dir)
    
    #Read in object
    suerat_obj = readRDS(celltype_files[i])
    
    #make integrated the default object
    DefaultAssay(suerat_obj) <- "integrated"
    
    #Remove hand dermatitis
    suerat_obj <- subset(x = suerat_obj, subset = dis != "Hand dermatitis")
    
    #Variable genes, scale, PCA
    suerat_obj <- FindVariableFeatures(suerat_obj, selection.method = "vst", nfeatures = 5000)
    all.genes <- rownames(suerat_obj)
    suerat_obj <- ScaleData(suerat_obj, features = all.genes)
    suerat_obj <- RunPCA(suerat_obj, features = VariableFeatures(object = suerat_obj))
    
    suerat_obj <- FindNeighbors(suerat_obj, dims = 1:10)
    suerat_obj <- FindClusters(suerat_obj, resolution = 0.5)
    suerat_obj@meta.data$cell_cluster = suerat_obj@active.ident #Add cell cluster id as metadata column
    
    #Run UMAP
    suerat_obj <- RunUMAP(suerat_obj, dims=1:10)
    
    #------------------------------------------------------------------------#
    #Loop to plot through metadata metrics
    #------------------------------------------------------------------------#
    names(suerat_obj@meta.data)
    table(suerat_obj@meta.data$sample_id)
    table(suerat_obj@meta.data$dis)
    table(suerat_obj@meta.data$cell_cluster)
    
    metadat_col_index = match(c("dis","cell_cluster","sample_id"),names(suerat_obj@meta.data))
    
    metadat_df = suerat_obj@meta.data #pull out all the cell metadata
    
    #Find metadata column index of metrics that I want to color cells by 
    umap_dir = paste(master_dir,
                     "umap/",
                     sep="")
    dir.create(umap_dir)
    
    for(k in 1:length(metadat_col_index)) {
      tryCatch({ #suppress error
        
        #k=1
        colourCount = length(unique(suerat_obj@meta.data[,metadat_col_index[k]]))
        
        #Fig: all factor levels together
        p1 <- DimPlot(suerat_obj, reduction = "umap", group.by = names(suerat_obj@meta.data)[metadat_col_index[k]],
                      cols=getPalette(colourCount),
                      raster=FALSE) +
          ggtitle( paste(celltype_file_shortname[i], names(suerat_obj@meta.data)[metadat_col_index[k]],sep=" ") ) + xlab("") + ylab("") +
          theme(axis.text=element_text(size=0,color='black'),
                axis.title=element_text(size=0),
                title=element_text(size=15),
                axis.ticks = element_line(size=0))
        
        png(paste(umap_dir,
                  k,".0_",names(suerat_obj@meta.data)[metadat_col_index[k]] ,".png",sep=""),
            width=9,height=8,units="in",res=300)
        print(p1)
        dev.off()
        
        factor_levels = na.omit(unique(suerat_obj@meta.data[,metadat_col_index[k]]))
        
        for(z in 1:length(factor_levels)){
          tryCatch({ #suppress error
            #z=1
            print(paste("k=",k," z=",z,sep=""))
            
            #Subset to factor level
            metadat_sub = metadat_df[metadat_df[,metadat_col_index[k]] %in% factor_levels[z],]
            
            p2 <- DimPlot(suerat_obj, reduction = "umap",
                          cols="lightgrey",
                          cells.highlight = row.names(metadat_sub),
                          cols.highlight = "darkred",
                          sizes.highlight=0.3,
                          raster=FALSE) +
              ggtitle( paste(celltype_file_shortname[i]," " ,names(suerat_obj@meta.data)[metadat_col_index[k]],": ", factor_levels[z],sep="")) + 
              xlab("") + ylab("") +
              theme(axis.text=element_text(size=0,color='black'),
                    axis.title=element_text(size=0),
                    title=element_text(size=15),
                    axis.ticks = element_line(size=0),
                    legend.position = "none")
            
            png(paste(umap_dir,
                      k,".",z,"_",names(suerat_obj@meta.data)[metadat_col_index[k]],"_",factor_levels[z],".png",sep=""),
                width=7,height=7,units="in",res=300)
            print(p2)
            dev.off()
            
            
          }, error=function(e){})
          
        } #End factor level loop
        
      }, error=function(e){})
      
    } # end metadata factor loop
    
    
    ###############################################################################################
    #DEGs among dis
    ###############################################################################################
    
    #Make deg output dir
    deg_dir = paste(master_dir,
                    "deg/",
                    sep="")
    dir.create(deg_dir)
    
    #10k most variable genes
    suerat_obj <- FindVariableFeatures(suerat_obj, selection.method = "vst", nfeatures = 10000)
    
    #PV
    deg_df_pv <- FindMarkers(suerat_obj,
                             ident.1 = row.names(suerat_obj@meta.data[suerat_obj@meta.data$dis %in% "Psoriasis Vulgaris",]),
                             ident.2 = row.names(suerat_obj@meta.data[suerat_obj@meta.data$dis %in% "Healthy control",]))
    deg_df_pv = deg_df_pv[order(-deg_df_pv$avg_log2FC),] #Re-order from high to low FC
    deg_df_pv = cbind(row.names(deg_df_pv),deg_df_pv) #Add col of gene symbols
    names(deg_df_pv)[1] = "gene"
    
    #PRP
    deg_df_prp <- FindMarkers(suerat_obj,
                              ident.1 = row.names(suerat_obj@meta.data[suerat_obj@meta.data$dis %in% "PRP",]),
                              ident.2 = row.names(suerat_obj@meta.data[suerat_obj@meta.data$dis %in% "Healthy control",]))
    deg_df_prp = deg_df_prp[order(-deg_df_prp$avg_log2FC),] #Re-order from high to low FC
    deg_df_prp = cbind(row.names(deg_df_prp),deg_df_prp) #Add col of gene symbols
    names(deg_df_prp)[1] = "gene"
    
    
    #GO and GSEA: go_gsea_fun(gene_rank_vector, geneUniverse)
    pv_vec = as.vector(abs(deg_df_pv$avg_log2FC)) #make ranked and named vector to put into fgsea
    names(pv_vec) = row.names(deg_df_pv)
    pv_enrichment = go_gsea_fun(pv_vec,row.names(suerat_obj))
    
    prp_vec = as.vector(abs(deg_df_prp$avg_log2FC)) #make ranked and named vector to put into fgsea
    names(prp_vec) = row.names(deg_df_prp)
    prp_enrichment = go_gsea_fun(prp_vec,row.names(suerat_obj))
    
    
    #Write out results
    OUT <- createWorkbook() #create excel workbook for a sample https://stackoverflow.com/questions/27713310/easy-way-to-export-multiple-data-frame-to-multiple-excel-worksheets
    
    sheet_name = paste(celltype_file_shortname[i]," PV", sep="")
    addWorksheet(OUT, sheet_name) # Add some sheets to the workbook
    writeData(OUT, sheet = sheet_name, x = deg_df_pv) #Write degs ncol(deg_df)
    writeData(OUT, sheet = sheet_name, x = pv_enrichment[[2]], startCol=8) #Write GO res ncol(tab)
    writeData(OUT, sheet = sheet_name, x = pv_enrichment[[1]], startCol=18) #Write GSEA res ncol(tab)
    
    sheet_name = paste(celltype_file_shortname[i]," PRP", sep="")
    addWorksheet(OUT, sheet_name) # Add some sheets to the workbook
    writeData(OUT, sheet = sheet_name, x = deg_df_prp) #Write degs ncol(deg_df)
    writeData(OUT, sheet = sheet_name, x = prp_enrichment[[2]], startCol=8) #Write GO res ncol(tab)
    writeData(OUT, sheet = sheet_name, x = prp_enrichment[[1]], startCol=18) #Write GSEA res ncol(tab)
    
    
    #---------------------------------------------------#
    #Now test differential expression in integrated a priori immune pathways
    #---------------------------------------------------#
    
    gene_set_comparisons = do.call(rbind, lapply(1:length(gene_set_cats), function(p) {
      
      tryCatch({ #suppress error
        #p=1
        print(p)
        gene_set_sub = subset(gene_sets, gs_name==gene_set_cats[p])
        
        metagene =  data.frame(suerat_obj@meta.data$sample_id, suerat_obj@meta.data$dis,
                               colSums(suerat_obj@assays$RNA@data[na.omit(match(gene_set_sub$gene_symbol,row.names(suerat_obj@assays$RNA@data))),]))
        names(metagene) = c("sample","dis","integrated_program")
        
        #get summaries of each dis
        dfc = summarySE(metagene, measurevar='integrated_program', groupvars=c('dis'))
        
        hc_n = subset(dfc, dis=="Healthy control")$N
        hc_mean = subset(dfc, dis=="Healthy control")$integrated_program
        hc_se = subset(dfc, dis=="Healthy control")$se
        
        pv_n = subset(dfc, dis=="Psoriasis Vulgaris")$N
        pv_mean = subset(dfc, dis=="Psoriasis Vulgaris")$integrated_program
        pv_se = subset(dfc, dis=="Psoriasis Vulgaris")$se
        
        prp_n = subset(dfc, dis=="PRP")$N
        prp_mean = subset(dfc, dis=="PRP")$integrated_program
        prp_se = subset(dfc, dis=="PRP")$se
        
        #Calculate FC
        prp_fc = prp_mean/hc_mean
        pv_fc = pv_mean/hc_mean
        prppv_fc = prp_mean/pv_mean
        
        #wilcox rank sum
        pv_wilcox = wilcox.test(subset(metagene, dis=="Healthy control")$integrated_program, subset(metagene, dis=="Psoriasis Vulgaris")$integrated_program, alternative = "two.sided")
        pv_W = pv_wilcox$statistic
        pv_p = pv_wilcox$p.value
        
        prp_wilcox = wilcox.test(subset(metagene, dis=="Healthy control")$integrated_program, subset(metagene, dis=="PRP")$integrated_program, alternative = "two.sided")
        prp_W = prp_wilcox$statistic
        prp_p = prp_wilcox$p.value
        
        prppv_wilcox = wilcox.test(subset(metagene, dis=="PRP")$integrated_program, subset(metagene, dis=="Psoriasis Vulgaris")$integrated_program, alternative = "two.sided")
        prppv_W = prp_wilcox$statistic
        prppv_p = prp_wilcox$p.value
        
        res_df = data.frame(gene_set_sub$gs_name[1],gene_set_sub$gs_id[1],gene_set_sub$gs_pmid[1],gene_set_sub$gs_geoid[1],gene_set_sub$gs_description[1],
                            hc_n, hc_mean, hc_se,
                            pv_n, pv_mean, pv_se,
                            prp_n, prp_mean, prp_se,
                            pv_fc,pv_W,pv_p,
                            prp_fc,prp_W,prp_p,
                            prppv_fc,prppv_W,prppv_p)
        
        names(res_df)[1:4] = c("gs_name","gs_id","gs_pmid","gs_description")
        return(res_df)
      }, error=function(e){})
    }))
    
    #Re-order from high to low PRP vs PV differential
    gene_set_comparisons = gene_set_comparisons[order(-gene_set_comparisons$prppv_fc),] #order from high to low fc: PRP  / PV
    sheet_name = paste("a priori pathway","", sep="")
    addWorksheet(OUT, sheet_name) # Add some sheets to the workbook
    writeData(OUT, sheet = sheet_name, x = gene_set_comparisons) #Write degs ncol(deg_df)
    
    #------------------------------------------------#
    # Export the Excel file
    #------------------------------------------------#
    saveWorkbook(OUT, 
                 paste(deg_dir, "dis_deg.xlsx",sep=""),
                 overwrite = T)
    
    
    #------------------------------------------------#
    #Now plot out a prior pathway expression differences
    #------------------------------------------------#
    #Top and bttm plots
    num_plots = 50
    
    #Make fig directory
    pathway_deg_fig_dir = paste(deg_dir, "a_priori_pathway_fig/",sep="")
    dir.create(pathway_deg_fig_dir)
    
    #-------------#
    #Top PRP vs PV
    #-------------#
    my_comparisons <- list( c("HC", "PRP"), c("HC", "PV"), c("PV", "PRP"))
    
    for(v in 1:num_plots){
      
      #v=1
      print(v)
      gene_set_sub = subset(gene_sets, gs_name==gene_set_comparisons$gs_name[v])
      
      metagene =  data.frame(suerat_obj@meta.data$sample_id, suerat_obj@meta.data$dis,
                             colSums(suerat_obj@assays$RNA@data[na.omit(match(gene_set_sub$gene_symbol,row.names(suerat_obj@assays$RNA@data))),]))
      names(metagene) = c("sample","dis","integrated_program")
      
      metagene$dis = gsub("Healthy control","HC", metagene$dis)
      metagene$dis = gsub("Psoriasis Vulgaris","PV", metagene$dis)
      
      my_comparisons <- list( c("HC", "PRP"), c("HC", "PV"),c("PRP", "PV"))
      
      
      pathway_expression_plot = ggplot(metagene, aes(x=dis,y=integrated_program, fill=dis, color=dis)) +
        stat_compare_means(comparisons = my_comparisons, method="wilcox.test")  +
        geom_violin(width=1, alpha=0.4) +
        geom_jitter(width=0.25, size=0.8, alpha=0.1) +
        geom_boxplot(outlier.size=0, width=0.3, color="black", size=0.7, alpha=0.4) +
        theme_classic() +
        ggtitle(paste(v,". ",gene_set_comparisons$gs_name[v],sep="")) + 
        xlab("") + ylab("gene set integrated expression") + 
        theme(axis.text.x=element_text(size=12,color='black'),
              axis.text.y=element_text(size=10,color='black'),
              title=element_text(size=12),
              axis.ticks.x = element_line(size=0),
              strip.text.x = element_text(size = 0),
              legend.position="none") +
        scale_fill_manual(name="",values = c("HC" = "darkblue",
                                             "PRP" = "darkred",
                                             "PV" = "darkorange3")) +
        scale_color_manual(name="",values = c("HC" = "darkblue",
                                              "PRP" = "darkred",
                                              "PV" = "darkorange3"))
      
      png(paste(pathway_deg_fig_dir,
                v,".",gene_set_comparisons$gs_name[v],".png",sep=""),
          width=5,height=7,units="in",res=200)
      print(pathway_expression_plot)
      dev.off()
    }
    
    
  }, error=function(e){})
  
} #End celltype integrated across sample DEG loop!!


scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/PRP_case_report/results/cluster_pca_analysis/v1 /Users/marktaylor/Desktop/Projects/Cho\ derm/Case\ Report\ PRP/results/testing_dir


###############################################################################################
###############################################################################################
###############################################################################################
#DEGs between PRP samples WITHIN celltype

#i=35
print(i)

#Output dir
#Make output directory
master_dir = paste("/c4/home/mtaylor4/derm/PRP_case_report/results/prp_sample_deg",
                   "/",
                   celltype_file_shortname[i],"/",sep="")
dir.create(master_dir)

#Read in object
suerat_obj = readRDS(celltype_files[i])

#make integrated the default object
DefaultAssay(suerat_obj) <- "integrated"

#Remove hand dermatitis
suerat_obj <- subset(x = suerat_obj, subset = dis != "Hand dermatitis")

#Variable genes, scale, PCA
suerat_obj <- FindVariableFeatures(suerat_obj, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(suerat_obj)
suerat_obj <- ScaleData(suerat_obj, features = all.genes)
suerat_obj <- RunPCA(suerat_obj, features = VariableFeatures(object = suerat_obj))

suerat_obj <- FindNeighbors(suerat_obj, dims = 1:10)
suerat_obj <- FindClusters(suerat_obj, resolution = 0.5)
suerat_obj@meta.data$cell_cluster = suerat_obj@active.ident #Add cell cluster id as metadata column

#Run UMAP
suerat_obj <- RunUMAP(suerat_obj, dims=1:10)


#Make deg output dir
deg_dir = paste(master_dir,
                "deg/",
                sep="")
dir.create(deg_dir)

#10k most variable genes
suerat_obj <- FindVariableFeatures(suerat_obj, selection.method = "vst", nfeatures = 10000)

#PV
unique(suerat_obj@meta.data$sample_id)
deg_df_pv <- FindMarkers(suerat_obj,
                         ident.1 = row.names(suerat_obj@meta.data[suerat_obj@meta.data$sample_id %in% "330_cd45",]),
                         ident.2 = row.names(suerat_obj@meta.data[suerat_obj@meta.data$sample_id %in% "331_cd45",]))
deg_df_pv = deg_df_pv[order(-deg_df_pv$avg_log2FC),] #Re-order from high to low FC
deg_df_pv = cbind(row.names(deg_df_pv),deg_df_pv) #Add col of gene symbols
names(deg_df_pv)[1] = "gene"

#GO and GSEA: go_gsea_fun(gene_rank_vector, geneUniverse)
pv_vec = as.vector(abs(deg_df_pv$avg_log2FC)) #make ranked and named vector to put into fgsea
names(pv_vec) = row.names(deg_df_pv)
pv_enrichment = go_gsea_fun(pv_vec,row.names(suerat_obj))

#Write out results
OUT <- createWorkbook() #create excel workbook for a sample https://stackoverflow.com/questions/27713310/easy-way-to-export-multiple-data-frame-to-multiple-excel-worksheets

sheet_name = paste(celltype_file_shortname[i],"PRP samples", sep="")
addWorksheet(OUT, sheet_name) # Add some sheets to the workbook
writeData(OUT, sheet = sheet_name, x = deg_df_pv) #Write degs ncol(deg_df)
writeData(OUT, sheet = sheet_name, x = pv_enrichment[[2]], startCol=8) #Write GO res ncol(tab)
writeData(OUT, sheet = sheet_name, x = pv_enrichment[[1]], startCol=18) #Write GSEA res ncol(tab)

#------------------------------------------------#
# Export the Excel file
#------------------------------------------------#
saveWorkbook(OUT, 
             paste(deg_dir, "prp_sample_deg.xlsx",sep=""),
             overwrite = T)



nrow(deg_df_pv[deg_df_pv$p_val_adj<0.05,])

mean(deg_df_pv[deg_df_pv$p_val_adj<0.05 & deg_df_pv$avg_log2FC>0,]$avg_log2FC)
mean(deg_df_pv[deg_df_pv$p_val_adj<0.05 & deg_df_pv$avg_log2FC<0,]$avg_log2FC)


scp mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/PRP_case_report/results/prp_sample_deg/8_Trm1/deg/prp_sample_deg.xlsx /Users/marktaylor/Desktop/Projects/Cho\ derm/Case\ Report\ PRP/results

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
#DEGs between PRP samples ACROSS celltype

#-----------------------------------------------------------------------------------#
#Combine PRP samples
#-----------------------------------------------------------------------------------#
rds_files = list.files("/c4/home/mtaylor4/derm/PRP_case_report/seurat_objects_cell_typed_combined_reference",full.names = T) #get all RDS objects in this directory
sampe_names = as.data.frame(t(data.frame(strsplit(rds_files,"/c4/home/mtaylor4/derm/PRP_case_report/seurat_objects_cell_typed_combined_reference"))))
sampe_names = gsub(".RDS","", sampe_names$V2)
sampe_names = gsub("/","", sampe_names)

#Get only PRP samples 
prp_samples = c("330_cd45","331_cd45")
prp_indices = match(prp_samples, sampe_names)

rds_list = do.call(list, lapply(prp_indices, function(i){
  # rds_list = do.call(list, lapply(c(1,3:10), function(i){
  
  tryCatch({ #suppress error
    #i=1
    print(i)
    suerat_obj = readRDS(rds_files[i])
    DefaultAssay(suerat_obj) <- 'RNA'  #Make default assay RNA
    
    #assign each cds object to a new name
    name <- paste("rds.", i, sep = "")
    return(assign(name, suerat_obj)) #Create a new object name that points to the Suerat object that was just read in
  }, error=function(e){})
}))


#Seurat integration: https://satijalab.org/seurat/articles/integration_introduction.html
features <- SelectIntegrationFeatures(object.list = rds_list)

#Integration anchors
#Speed up FindInetegrationAnchors: https://github.com/satijalab/seurat/discussions/3999
library(future)
plan("multisession", workers = 4) # Enable parallelization from future package
options(future.globals.maxSize = 20000 * 1024^2)
cell.anchors <- FindIntegrationAnchors(object.list = rds_list, 
                                       anchor.features = features,
                                       #reference = c(1,2,3,4,5),
                                       reduction = "rpca",
                                       dims = 1:30
)


#Now actually integrate the data
plan("sequential") # Stop paralellization
suerat_obj <- IntegrateData(anchorset = cell.anchors,
                            k.weight = 30) # https://github.com/satijalab/seurat/issues/3930

saveRDS(suerat_obj,
        file=paste("/c4/home/mtaylor4/derm/PRP_case_report/",
                   "prp_samples_combined.RDS",
                   sep=""))

#-----------------------------------------------------------------------------------#
#DEGs between PRP samples
#-----------------------------------------------------------------------------------#
#Read in object
suerat_obj = readRDS("/c4/home/mtaylor4/derm/PRP_case_report/prp_samples_combined.RDS")

#make integrated the default object
DefaultAssay(suerat_obj) <- "integrated"

#Variable genes, scale, PCA
suerat_obj <- FindVariableFeatures(suerat_obj, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(suerat_obj)
suerat_obj <- ScaleData(suerat_obj, features = all.genes)
suerat_obj <- RunPCA(suerat_obj, features = VariableFeatures(object = suerat_obj))

suerat_obj <- FindNeighbors(suerat_obj, dims = 1:10)
suerat_obj <- FindClusters(suerat_obj, resolution = 0.5)
suerat_obj@meta.data$cell_cluster = suerat_obj@active.ident #Add cell cluster id as metadata column

#Make deg output dir
deg_dir = paste("/c4/home/mtaylor4/derm/PRP_case_report/results/prp_intersample_deg/")
dir.create(deg_dir)

#PV
deg_df_pv <- FindMarkers(suerat_obj,
                         ident.1 = row.names(suerat_obj@meta.data[suerat_obj@meta.data$sample_id %in% "330_cd45",]),
                         ident.2 = row.names(suerat_obj@meta.data[suerat_obj@meta.data$sample_id %in% "331_cd45",]))
deg_df_pv = deg_df_pv[order(-deg_df_pv$avg_log2FC),] #Re-order from high to low FC
deg_df_pv = cbind(row.names(deg_df_pv),deg_df_pv) #Add col of gene symbols
names(deg_df_pv)[1] = "gene"

#GO and GSEA: go_gsea_fun(gene_rank_vector, geneUniverse)
pv_vec = as.vector(abs(deg_df_pv$avg_log2FC)) #make ranked and named vector to put into fgsea
names(pv_vec) = row.names(deg_df_pv)
pv_enrichment = go_gsea_fun(pv_vec,row.names(suerat_obj))

#Write out results
OUT <- createWorkbook() #create excel workbook for a sample https://stackoverflow.com/questions/27713310/easy-way-to-export-multiple-data-frame-to-multiple-excel-worksheets

sheet_name = paste("PRP samples", sep="")
addWorksheet(OUT, sheet_name) # Add some sheets to the workbook
writeData(OUT, sheet = sheet_name, x = deg_df_pv) #Write degs ncol(deg_df)
writeData(OUT, sheet = sheet_name, x = pv_enrichment[[2]], startCol=8) #Write GO res ncol(tab)
writeData(OUT, sheet = sheet_name, x = pv_enrichment[[1]], startCol=18) #Write GSEA res ncol(tab)

#------------------------------------------------#
# Export the Excel file
#------------------------------------------------#
saveWorkbook(OUT, 
             paste(deg_dir, "prp_intersample_deg.xlsx",sep=""),
             overwrite = T)


scp mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/PRP_case_report/results/prp_intersample_deg/prp_intersample_deg.xlsx /Users/marktaylor/Desktop/Projects/Cho\ derm/Case\ Report\ PRP/results

###############################################################################################
###############################################################################################

#Merge / Combine individual healthy control samples into a single combined HC object

library(Seurat)
library(future)

rds_files = list.files("/c4/home/mtaylor4/derm/PRP_case_report/seurat_objects_cell_typed_combined_reference",full.names = T) #get all RDS objects in this directory
sampe_names = as.data.frame(t(data.frame(strsplit(rds_files,"/c4/home/mtaylor4/derm/PRP_case_report/seurat_objects_cell_typed_combined_reference"))))
sampe_names = gsub(".RDS","", sampe_names$V2)
sampe_names = gsub("/","", sampe_names)

#----------------------------------------------------------------------------------#
#Make combined normal
#----------------------------------------------------------------------------------#
hc_samples = c("149","150","154","155","169","195","204","207","239","241","242","247","248","267","286")
hc_indices = match(hc_samples, sampe_names)

rds_list = do.call(list, lapply(hc_indices, function(i){
  # rds_list = do.call(list, lapply(c(1,3:10), function(i){
  
  tryCatch({ #suppress error
    #i=1
    print(i)
    suerat_obj = readRDS(rds_files[i])
    DefaultAssay(suerat_obj) <- 'RNA'  #Make default assay RNA
    
    #assign each cds object to a new name
    name <- paste("rds.", i, sep = "")
    return(assign(name, suerat_obj)) #Create a new object name that points to the Suerat object that was just read in
  }, error=function(e){})
}))

hc.merged <- merge(rds_list[[1]], y = c(rds_list[2:3]), project = "healthy_controls",
                   merge.data = TRUE)


hc.merged <- merge(rds_list[[1]], y = c(rds_list[2:length(rds_list)]), project = "healthy_controls",
                   merge.data = TRUE)

saveRDS(hc.merged,
        file=paste("/c4/home/mtaylor4/derm/PRP_case_report/",
                   "healthy_controls_combined.RDS",
                   sep=""))


###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
#PCA
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################

#Find PRP-specific 
deg_res_dirs = list.dirs(version_dir_tomake, recursive=F)
deg_res_dir_names = list.dirs(version_dir_tomake, recursive=F,full.names = F)

OUT <- createWorkbook()
sheet_name = paste("PRP-sepcific DEGs",sep="")
addWorksheet(OUT, sheet_name) # Add some sheets to the workbook

prp_unique = do.call(list, lapply(1:length(deg_res_dirs), function(d){
  #prp_unique = do.call(list, lapply(1:3, function(d){
  
  tryCatch({ #suppress error
    #d=3
    print(d)
    
    #Read in DEG files
    pv_df = read.xlsx( paste(deg_res_dirs[d],"/deg/dis_deg.xlsx",sep="") ,sheet=1)
    pv_df = subset(pv_df, p_val_adj<0.05)
    
    prp_df = read.xlsx( paste(deg_res_dirs[d],"/deg/dis_deg.xlsx",sep="") ,sheet=2)
    prp_df = subset(prp_df, p_val_adj<0.05)
    
    #Get unique PRP genes that do not appear in PV list
    # Create two character vectors
    prp_specific_sig_degs = data.frame(setdiff(prp_df$gene, pv_df$gene))
    names(prp_specific_sig_degs) = deg_res_dir_names[d]
    
    #Write out to excel workbook
    writeData(OUT, sheet = sheet_name, x = prp_specific_sig_degs, startCol=d ) #Write data to a sheet
    
    return(prp_specific_sig_degs)
  }, error=function(e){})
}))

# Export the Excel file
saveWorkbook(OUT, 
             paste(version_dir_tomake,"/0_prp_specific_degs.xlsx",sep=""),
             overwrite=T)

#Transfer just PRP_unique genes
scp mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/PRP_case_report/results/cluster_pca_analysis/v1/1_prp_specific_degs.xlsx /Users/marktaylor/Desktop/Projects/Cho\ derm/Case\ Report\ PRP/results/
  
  
  #Transfer all DEG results
  scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/PRP_case_report/results/cluster_pca_analysis/ /Users/marktaylor/Desktop/Projects/Cho\ derm/Case\ Report\ PRP/results/
  
  
  
  
  ###############################################################################################
###############################################################################################
###############################################################################################

#Find metadata column index of metrics that I want to color cells by 
pca_dir = paste(master_dir,
                "pca/",
                sep="")
dir.create(pca_dir)

#Add UMAP and PCA cell embeddings to the Seurat object
metadat_df = suerat_obj@meta.data
metadat_df = cbind(metadat_df,suerat_obj@reductions$umap@cell.embeddings,suerat_obj@reductions$pca@cell.embeddings)

###############################################################################################
#Umap means for comparison
###############################################################################################

umap_factor_function = function(metadat_df_input, factor_1) {
  
  # metadat_df_input=metadat_df
  # factor_1 = "dis_level_1"
  
  dfc1 = summarySE(metadat_df_input, measurevar="UMAP_1", groupvars=c(factor_1))
  dfc2 = summarySE(metadat_df_input, measurevar="UMAP_2", groupvars=c(factor_1))
  
  #combine the gene signatures
  re_int = cbind(dfc1[,c(1,3,5)], dfc2[,c(3,5)])
  names(re_int) = c("factor_level","umap1","umap1_se","umap2","umap2_se")
  
  p3 = ggplot(data=re_int, aes(x=umap1, y=umap2)) +
    geom_point(alpha=1,aes(color=factor_level), size=4) +
    geom_linerange(aes(xmax = umap1 + umap1_se, xmin = umap1 - umap1_se, color=factor_level),size=1) +
    geom_linerange(aes(ymax = umap2 + umap2_se, ymin = umap2 - umap2_se, color=factor_level),size=1) +
    geom_text_repel(data=re_int, aes(x=umap1, y=umap2,label=factor_level),box.padding = 0.4) + #add labels
    #stat_ellipse(aes(color=dis)) +
    theme_classic() +
    ggtitle(paste("UMAP centroids for ",factor_1,sep="")) +
    xlab("Axis 1 mean") +
    ylab("Axis 2 mean") +
    theme(legend.position="right") +
    scale_color_manual(values=getPalette(length(re_int$factor_level))) +
    theme(axis.text=element_text(size=15,color='black'),
          axis.title=element_text(size=15,color='black'),
          plot.title=element_text(size=15,color='black'),
          legend.position = "none") 
  
  png(paste(pca_dir,
            factor_1,"UMAP_means.png",sep=""),
      width=10,height=10,units="in",res=300)
  print(p3)
  dev.off()
  
}

factors_to_analyze = c("cell_cluster","dis_level_1","dis_level_2")

for(g in 1:length(factors_to_analyze)){
  tryCatch({ #suppress error
    print(g)
    umap_factor_function(metadat_df, factors_to_analyze[g])
  }, error=function(e){})
}


#----------------------------------------------------------------------------------#
#Now look for PCs
#----------------------------------------------------------------------------------#

#Map all PCs
pcafig_dir = paste(pca_dir,
                   "pc_fig/",
                   sep="")
dir.create(pcafig_dir)

for(h in 1:30) {
  
  #h=1
  print(h)
  pc_plot = ggplot(data=metadat_df, aes_string(x=paste("PC_",h,sep=""), y=paste("PC_",(h+1),sep=""))) +
    geom_point(alpha=0.6,aes(color=cell_cluster), size=0.4) +
    theme_classic() +
    ggtitle(paste("PC ",h," & ", (h+1), sep="")) +
    xlab(paste("PC ",h,sep="")) +
    ylab(paste("PC ",(h+1),sep="")) +
    theme(legend.position="right") +
    scale_color_manual(values=getPalette(length(unique(metadat_df$cell_cluster)))) +
    theme(axis.text=element_text(size=15,color='black'),
          axis.title=element_text(size=15,color='black'),
          plot.title=element_text(size=15,color='black'),
          legend.position = "right") +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    geom_hline(yintercept = 0, linetype="dotted",size=1) +
    geom_vline(xintercept = 0, linetype="dotted",size=1)
  
  
  png(paste(pcafig_dir,
            "PC_",h,"_", (h+1),".png",sep=""),
      width=11,height=10,units="in",res=300)
  print(pc_plot)
  dev.off()
  
} #end PC loop

#Now factor centroids
pcafactorfig_dir = paste(pca_dir,
                         "pcfactor_fig/",
                         sep="")
dir.create(pcafactorfig_dir)

pca_factor_function = function(metadat_df_input, factor_1) {
  
  # metadat_df_input=metadat_df
  # factor_1 = "predicted.id"
  
  for(h in 1:30) {
    tryCatch({ #suppress error
      
      #h=1
      dfc1 = summarySE(metadat_df_input, measurevar=paste("PC_",h,sep=""), groupvars=c(factor_1))
      dfc2 = summarySE(metadat_df_input, measurevar=paste("PC_",(h+1),sep=""), groupvars=c(factor_1))
      
      #combine the gene signatures
      re_int = cbind(dfc1[,c(1,3,5)], dfc2[,c(3,5)])
      names(re_int) = c("factor_level","umap1","umap1_se","umap2","umap2_se")
      
      p1 = ggplot(data=re_int, aes(x=umap1, y=umap2)) +
        geom_point(alpha=1,aes(color=factor_level), size=4) +
        geom_linerange(aes(xmax = umap1 + umap1_se, xmin = umap1 - umap1_se, color=factor_level),size=1) +
        geom_linerange(aes(ymax = umap2 + umap2_se, ymin = umap2 - umap2_se, color=factor_level),size=1) +
        geom_text_repel(data=re_int, aes(x=umap1, y=umap2,label=factor_level),box.padding = 0.4,size=7) + #add labels
        #stat_ellipse(aes(color=dis)) +
        theme_classic() +
        ggtitle(paste("PCA centroids for ",factor_1,": PC ",h," & ", (h+1), sep="")) +
        xlab(paste("PC ",h,sep="")) +
        ylab(paste("PC ",(h+1),sep="")) +
        theme(legend.position="right") +
        scale_color_manual(values=getPalette(length(re_int$factor_level))) +
        theme(axis.text=element_text(size=20,color='black'),
              axis.title=element_text(size=20,color='black'),
              plot.title=element_text(size=20,color='black'),
              legend.position = "none") +
        geom_hline(yintercept = 0, linetype="dotted",size=1) +
        geom_vline(xintercept = 0, linetype="dotted",size=1)
      
      png(paste(pcafactorfig_dir,
                factor_1,"_","PC_",h,"_", (h+1),".png",sep=""),
          width=10,height=10,units="in",res=300)
      print(p1)
      dev.off()
    }, error=function(e){})
    
  } #end PC loop
} #end PC factor function

#Call on a bunch of metadata fields
factors_to_analyze = c("cell_cluster","dis_level_1","dis_level_2")

for(g in 1:length(factors_to_analyze)){
  tryCatch({ #suppress error
    print(g)
    pca_factor_function(metadat_df, factors_to_analyze[g])
  }, error=function(e){})
}

#----------------------------------------------------------------------------------#
#Find gene loadings onto the PCs

#Now plot loadings centroids
pca_loadings_fig = paste(pca_dir,
                         "pc_loading_fig/",
                         sep="")
dir.create(pca_loadings_fig)


loadings = suerat_obj@reductions$pca@feature.loadings
n_loadings = 50 #number of gene loadings to plot 

for(e in 1:ncol(loadings)) {
  tryCatch({ #suppress error
    #e=1
    print(e)
    pc_ordered = data.frame(loadings[order(-loadings[,e]),])
    pc_ordered$gene = row.names(pc_ordered)
    pc_ordered$gene = factor(pc_ordered$gene,levels=pc_ordered$gene) #re-order genes
    
    #plot top gene loadings
    top_pc = head(pc_ordered,n_loadings)
    top_pc = top_pc[order(top_pc[,e]),]
    top_pc$gene = factor(top_pc$gene,levels=top_pc$gene) #re-order genes
    
    pos_load_fig = ggplot() +
      geom_bar(data=top_pc, aes(x=top_pc[,e], y=gene), stat="identity",color="white",fill="darkred") +
      theme_classic() +
      ggtitle(paste("PC",e," POSITIVE loadings", sep="")) +
      xlab("loading") +
      ylab("") +
      theme(legend.position="none") +
      scale_color_manual(values=getPalette(length(re_int$factor_level))) +
      theme(axis.text.x=element_text(size=12,color='black'),
            axis.text.y=element_text(size=17,color='black'),
            axis.title=element_text(size=12,color='black'),
            plot.title=element_text(size=20,color='black'),
            legend.position = "none")
    
    #plot bttm gene loadings
    bttm_pc = tail(pc_ordered,n_loadings)
    
    neg_load_fig = ggplot() +
      geom_bar(data=bttm_pc, aes(x=-top_pc[,e], y=gene), stat="identity",color="white",fill="darkblue") +
      theme_classic() +
      ggtitle(paste("PC",e," NEGATIVE loadings", sep="")) +
      xlab("loading") +
      ylab("") +
      theme(legend.position="none") +
      scale_color_manual(values=getPalette(length(re_int$factor_level))) +
      theme(axis.text.x=element_text(size=12,color='black'),
            axis.text.y=element_text(size=17,color='black'),
            axis.title=element_text(size=12,color='black'),
            plot.title=element_text(size=20,color='black'),
            legend.position = "none")
    
    png(paste(pca_loadings_fig,
              "PC_",e,"_","loadings.png",sep=""),
        width=10,height=15,units="in",res=300)
    print(multiplot(pos_load_fig,neg_load_fig,cols=2))
    dev.off()
  }, error=function(e){})
}

OUT <- createWorkbook() #create excel workbook for a sample https://stackoverflow.com/questions/27713310/easy-way-to-export-multiple-data-frame-to-multiple-excel-worksheets

for(r in 1:20) {
  tryCatch({ #suppress error
    #r=1
    print(r)
    pc_ordered = data.frame(loadings[order(-loadings[,r]),])
    pc_ordered$gene = row.names(pc_ordered)
    pc_ordered$gene = factor(pc_ordered$gene,levels=pc_ordered$gene) #re-order genes
    
    #plot top gene loadings
    top_pc = head(pc_ordered,n_loadings)
    
    #plot bttm gene loadings
    bttm_pc = tail(pc_ordered,n_loadings)
    
    
    #--------------------------------------------#
    #GSEA: top
    #--------------------------------------------#
    gene_ranks = as.vector(abs(top_pc[,i])) #make ranked and named vector to put into fgsea
    names(gene_ranks) = row.names(top_pc)
    
    fgseaRes <- fgsea(pathways = pathwaysH, 
                      stats    = gene_ranks,
                      minSize  = 10,
                      maxSize  = 500,
                      scoreType="pos")
    ordered_gsea = fgseaRes[order(fgseaRes$pval),]
    
    #--------------------------------------------#
    #GO: top
    #--------------------------------------------#
    geneUniverse <- row.names(suerat_obj)
    genesOfInterest = row.names(top_pc)
    geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
    names(geneList) <- geneUniverse
    GOdata <- new("topGOdata",
                  ontology = "BP",
                  allGenes=geneList,
                  #geneSel = selection,
                  annot = annFUN.org, 
                  mapping = "org.Hs.eg.db",
                  ID = "symbol")
    
    resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "fisher") #Fisher exact test
    tab <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 120)
    tab = subset(tab, raw.p.value <=0.01) #keep only sig results
    go_candidates = genesOfInterest
    selcGenes <- genesInTerm(GOdata, whichGO=tab$GO.ID)
    genes_in_go = do.call(rbind, lapply(1:length(selcGenes), function(p){
      #i=1
      data.frame(names(selcGenes)[p],  paste( c(na.omit(selcGenes[[p]][match(go_candidates,  selcGenes[[p]])])) , collapse=", ") )
    }))
    names(genes_in_go) = c("GO_ID","sig_genes_in_term")
    tab = cbind(tab, genes_in_go)
    
    #Write out to excel workbook
    topload_genes = data.frame(row.names(top_pc))
    names(topload_genes) = paste("PC",r," pos loading genes",sep="")
    addWorksheet(OUT, paste("PC",r," pos",sep="")) # Add some sheets to the workbook
    writeData(OUT, sheet = paste("PC",r," pos",sep=""), x = topload_genes) #Write degs ncol(deg_df)
    writeData(OUT, sheet = paste("PC",r," pos",sep=""), x = tab, startCol=3) #Write GO res ncol(tab)
    writeData(OUT, sheet = paste("PC",r," pos",sep=""), x = ordered_gsea, startCol=13) #Write GO res ncol(tab)
    
    #--------------------------------------------#
    #GSEA: bttm
    #--------------------------------------------#
    gene_ranks = as.vector(abs(bttm_pc[,r])) #make ranked and named vector to put into fgsea
    names(gene_ranks) = row.names(bttm_pc)
    
    fgseaRes <- fgsea(pathways = pathwaysH, 
                      stats    = gene_ranks,
                      minSize  = 10,
                      maxSize  = 500,
                      scoreType="pos")
    ordered_gsea = fgseaRes[order(fgseaRes$pval),]
    
    #--------------------------------------------#
    #GO: bttm
    #--------------------------------------------#
    geneUniverse <- row.names(suerat_obj)
    genesOfInterest = row.names(bttm_pc)
    geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
    names(geneList) <- geneUniverse
    GOdata <- new("topGOdata",
                  ontology = "BP",
                  allGenes=geneList,
                  #geneSel = selection,
                  annot = annFUN.org, 
                  mapping = "org.Hs.eg.db",
                  ID = "symbol")
    
    resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "fisher") #Fisher exact test
    tab <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 120)
    tab = subset(tab, raw.p.value <=0.01) #keep only sig results
    go_candidates = genesOfInterest
    selcGenes <- genesInTerm(GOdata, whichGO=tab$GO.ID)
    genes_in_go = do.call(rbind, lapply(1:length(selcGenes), function(p){
      #i=1
      data.frame(names(selcGenes)[p],  paste( c(na.omit(selcGenes[[p]][match(go_candidates,  selcGenes[[p]])])) , collapse=", ") )
    }))
    names(genes_in_go) = c("GO_ID","sig_genes_in_term")
    tab = cbind(tab, genes_in_go)
    
    #Write out to excel workbook
    topload_genes = data.frame(row.names(bttm_pc))
    names(topload_genes) = paste("PC",r," neg loading genes",sep="")
    addWorksheet(OUT, paste("PC",r," neg",sep="")) # Add some sheets to the workbook
    writeData(OUT, sheet = paste("PC",r," neg",sep=""), x = topload_genes) #Write degs ncol(deg_df)
    writeData(OUT, sheet = paste("PC",r," neg",sep=""), x = tab, startCol=3) #Write GO res ncol(tab)
    writeData(OUT, sheet = paste("PC",r," neg",sep=""), x = ordered_gsea, startCol=13) #Write GO res ncol(tab)
    
    
  }, error=function(e){})
  
}

# Export the Excel file
saveWorkbook(OUT, 
             paste(pca_dir,"PC_gene_loadings.xlsx",sep=""),
             overwrite = T)



#}, error=function(e){})
#} #End cell type integration loop


#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#

#Now permute through UMAP settings and accentuate cell difference

celltype_files = list.files("/c4/home/mtaylor4/derm/combined_data/data/cell_types_integrated/",full.names = T)
celltype_file_shortname = list.files("/c4/home/mtaylor4/derm/combined_data/data/cell_types_integrated/",full.names = F)
celltype_file_shortname = unlist(strsplit(celltype_file_shortname,".RDS"))

sample_metadat_match = read.xlsx(
  "/c4/home/mtaylor4/derm/combined_data/data/sample_dis_match.xlsx"
)
sample_metadat_match = subset(sample_metadat_match, Dis != "Sarcoid")


version_dir = "v1"
for(i in 1:length(celltype_files)){
  tryCatch({ #suppress error
    
    #i=40
    
    #Output dir
    #Concatenate string for location directory
    results_dir = paste("/c4/home/mtaylor4/derm/combined_data/results/cluster_pca_analysis/umap_solutions/",
                        version_dir,"/",sep="")
    dir.create(results_dir)
    
    results_dir = paste("/c4/home/mtaylor4/derm/combined_data/results/cluster_pca_analysis/umap_solutions/",
                        version_dir,"/",
                        celltype_file_shortname[i],"/",sep="")
    dir.create(results_dir)
    
    #Read in object
    suerat_obj = readRDS(celltype_files[i])
    DefaultAssay(suerat_obj) <- "integrated"
    
    #Add metadata columns with sample and then revamp in loop below
    suerat_obj@meta.data$dis_level_1 = suerat_obj@meta.data$sample.id
    suerat_obj@meta.data$dis_level_2 = suerat_obj@meta.data$sample.id
    
    #Run through all the sample
    for(j in 1:nrow(sample_metadat_match)){
      #j=1
      print(j)
      # suerat_obj@meta.data$dis_level_1 = gsub(paste("\\b",sample_metadat_match$sample.id[j],"\\b",sep=""), sample_metadat_match$Dis[j],suerat_obj@meta.data$dis_level_1)
      # suerat_obj@meta.data$dis_level_2 = gsub(sample_metadat_match$sample.id[j], sample_metadat_match$Dis_2[j],suerat_obj@meta.data$dis_level_2)
      suerat_obj@meta.data$dis_level_1[ grep(sample_metadat_match$sample.id[j],suerat_obj@meta.data$dis_level_1) ] <- sample_metadat_match$Dis[j]
      suerat_obj@meta.data$dis_level_2[ grep(sample_metadat_match$sample.id[j],suerat_obj@meta.data$dis_level_2) ] <- sample_metadat_match$Dis_2[j]
      
    }
    
    #Remove NA's in the Suerat object
    suerat_obj@meta.data[is.na(suerat_obj@meta.data)] <- "No_match" #replace NA's in metadata with No_match
    suerat_obj <- subset(x = suerat_obj, subset = dis_level_1 != "No_match")
    
    
    suerat_obj <- FindVariableFeatures(suerat_obj, selection.method = "vst", nfeatures = 5000)
    all.genes <- rownames(suerat_obj)
    suerat_obj <- ScaleData(suerat_obj, features = all.genes)
    suerat_obj <- RunPCA(suerat_obj, features = VariableFeatures(object = suerat_obj))
    
    suerat_obj <- FindNeighbors(suerat_obj, dims = 1:10)
    suerat_obj <- FindClusters(suerat_obj, resolution = 0.5)
    suerat_obj@meta.data$cell_cluster = suerat_obj@active.ident #Add cell cluster id as metadata column
    
    #Find metadata column index of metrics that I want to color cells by 
    metadat_col_index = match(c("dis_level_1"),names(suerat_obj@meta.data))
    
    #umap dir
    umap_dir = paste(results_dir,
                     "umap/",sep="")
    dir.create(umap_dir)
    
    for(k in 1:length(metadat_col_index)) {
      tryCatch({ #suppress error
        #k=1
        
        min_distances = seq(1,0.1,by=-0.1)
        n.neighbor_vec = seq(30,300,by=30)
        
        mindist_nn_df = do.call(rbind, lapply(1:length(min_distances), function(g) { #NEAREST NEIGHBOR LOOP
          #for(g in 1:length(min_distances)){ #MIN DISTANCE LOOPS
          tryCatch({ #suppress error
            #g=1
            
            nn_df =  do.call(rbind, lapply(1:length(n.neighbor_vec), function(h) { #NEAREST NEIGHBOR LOOP
              #for(h in 1:length(n.neighbor_vec)){ #NEAREST NEIGHBOR LOOP
              tryCatch({ #suppress error
                #h=1
                
                print(paste("g=",g," h=",h,sep=""))
                suerat_obj <- RunUMAP(suerat_obj, dims=1:10,
                                      n.neighbors = n.neighbor_vec[h],
                                      min.dist=min_distances[g])
                
                colourCount = length(unique(suerat_obj@meta.data[,metadat_col_index[k]]))
                
                #Fig: all factor levels together
                p1 <- DimPlot(suerat_obj, reduction = "umap", group.by = names(suerat_obj@meta.data)[metadat_col_index[k]],
                              cols=getPalette(colourCount),
                              raster=FALSE) +
                  ggtitle( paste(celltype_file_shortname[i], names(suerat_obj@meta.data)[metadat_col_index[k]],
                                 "Min dist=",min_distances[g], " Nn=",n.neighbor_vec[h],sep=" ") ) + xlab("") + ylab("") +
                  theme(axis.text=element_text(size=15,color='black'),
                        axis.title=element_text(size=0),
                        title=element_text(size=15),
                        axis.ticks = element_line(size=0))
                
                png(paste(umap_dir,
                          g,".",h,".MinDist_",min_distances[g],"_Nn_",n.neighbor_vec[h],".png",sep=""),
                    width=9,height=8,units="in",res=300)
                print(p1)
                dev.off()
                
                #Get all umaps and diseases
                umap_coords = data.frame(suerat_obj$dis_level_1,suerat_obj[["umap"]]@cell.embeddings)
                names(umap_coords)[1] = "dis"
                
                umap1_dfc = summarySE(umap_coords, measurevar="UMAP_1", groupvars=c("dis"))
                names(umap1_dfc)[c(5)] = c("umap1_se")
                umap_1_variance_mean = mean(umap1_dfc$umap1_se)
                
                umap2_dfc= summarySE(umap_coords, measurevar="UMAP_2", groupvars=c("dis"))
                names(umap2_dfc)[c(5)] = c("umap2_se")
                umap_2_variance_mean = mean(umap2_dfc$umap2_se)
                
                #Get all umaps and diseases
                all_coords_mat = as.matrix(data.frame(umap1_dfc$UMAP_1,umap2_dfc$UMAP_2))
                row.names(all_coords_mat) = umap1_dfc$dis
                dist_mat = pointDistance(all_coords_mat, allpairs=T, lonlat=F) #calculate all-vs-all distance
                rownames(dist_mat) = umap1_dfc$dis
                colnames(dist_mat) = umap1_dfc$dis
                distance_vec = as.vector(dist_mat) #convert all the distances to a vector
                distance_vec = distance_vec[distance_vec!=0]
                distance_mean = mean(distance_vec)
                variance_normalized_mean = distance_mean/umap_1_variance_mean/umap_2_variance_mean #normalize mean by variance
                
                distance_matrix = data.frame(celltype_file_shortname[i],min_distances[g],n.neighbor_vec[h],distance_mean,umap_1_variance_mean,umap_2_variance_mean,variance_normalized_mean)
                names(distance_matrix) = c("celltype","min_dist","NearestNeighb","interdisease_mean_dist","umap_1_variance_mean","umap_2_variance_mean","variance_normalized_mean")
                distance_matrix   
              }, error=function(e){})
              
            })) #End nearest neighbor loop
            
            nn_df #Return nearest neighbor results
            
          }, error=function(e){})
          
        })) #End min dist loop
        
        write.table(mindist_nn_df,
                    file=paste(results_dir,"umap_clustering_analysis_table.txt"),
                    quote=F,row.names=F,sep="\t")
        
      }, error=function(e){})
    } #end metadat col for loop
    
  }, error=function(e){})
} #End cell type loop


#scp results to my machine
scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/combined_data/results/cluster_pca_analysis/umap_solutions /Users/marktaylor/Desktop/Projects/Cho\ derm/pan-rash/Results/cluster_pca_analysis

#-----------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------#
#Pathologic axis analysis to find most constributory cells
#-----------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------#

celltype_files = list.files("/c4/home/mtaylor4/derm/combined_data/data/cell_types_integrated/",full.names = T)
celltype_file_shortname = list.files("/c4/home/mtaylor4/derm/combined_data/data/cell_types_integrated/",full.names = F)
celltype_file_shortname = unlist(strsplit(celltype_file_shortname,".RDS"))

sample_metadat_match = read.xlsx(
  "/c4/home/mtaylor4/derm/combined_data/data/sample_dis_match.xlsx"
)
sample_metadat_match = subset(sample_metadat_match, Dis != "Sarcoid")


version_dir = "v2"
for(i in 1:length(celltype_files)){
  tryCatch({ #suppress error
    
    #i=1
    
    #Read in Seurat object
    suerat_obj = readRDS(celltype_files[i])
    DefaultAssay(suerat_obj) <- "integrated"
    
    #Add metadata columns with sample and then revamp in loop below
    suerat_obj@meta.data$dis_level_1 = suerat_obj@meta.data$sample.id
    suerat_obj@meta.data$dis_level_2 = suerat_obj@meta.data$sample.id
    
    #Run through all the sample
    for(j in 1:nrow(sample_metadat_match)){
      #j=1
      print(j)
      # suerat_obj@meta.data$dis_level_1 = gsub(paste("\\b",sample_metadat_match$sample.id[j],"\\b",sep=""), sample_metadat_match$Dis[j],suerat_obj@meta.data$dis_level_1)
      # suerat_obj@meta.data$dis_level_2 = gsub(sample_metadat_match$sample.id[j], sample_metadat_match$Dis_2[j],suerat_obj@meta.data$dis_level_2)
      suerat_obj@meta.data$dis_level_1[ grep(sample_metadat_match$sample.id[j],suerat_obj@meta.data$dis_level_1) ] <- sample_metadat_match$Dis[j]
      suerat_obj@meta.data$dis_level_2[ grep(sample_metadat_match$sample.id[j],suerat_obj@meta.data$dis_level_2) ] <- sample_metadat_match$Dis_2[j]
      
    }
    
    #Remove NA's in the Suerat object
    suerat_obj@meta.data[is.na(suerat_obj@meta.data)] <- "No_match" #replace NA's in metadata with No_match
    suerat_obj <- subset(x = suerat_obj, subset = dis_level_1 != "No_match")
    
    
    suerat_obj <- FindVariableFeatures(suerat_obj, selection.method = "vst", nfeatures = 5000)
    all.genes <- rownames(suerat_obj)
    suerat_obj <- ScaleData(suerat_obj, features = all.genes)
    
    #Read in DEG file
    deg_res_file = paste("/c4/home/mtaylor4/derm/combined_data/results/cluster_pca_analysis/v",
                         version_dir,"/",
                         celltype_file_shortname[i],"/deg/cluster_specific_deg.xlsx",sep="")
    
    #Get DEG sheets
    sheet_names = getSheetNames(file=deg_res_file)
    dis_sheets = sheet_names[grep("dis_level_1",sheet_names)]
    
    go_thru_dis_markers = do.call(rbind, lapply(1:length(dis_sheets), function(j){
      tryCatch({ #suppress error
        
        #j=1
        print(paste("i=",i," j=",j,sep=""))
        
        marker_file = read.xlsx(
          xlsxFile=deg_res_file,
          sheet=dis_sheets[j],
          startRow = 1,
          colNames = TRUE,
          rowNames = FALSE,
          detectDates = FALSE,
          skipEmptyRows = TRUE,
          skipEmptyCols = TRUE,
          rows = NULL,
          cols = c(1:6),
          check.names = FALSE,
          sep.names = ".",
          namedRegion = NULL,
          na.strings = "NA",
          fillMergedCells = FALSE
        )
        
        marker_file_pos = subset(marker_file, avg_log2FC>0)
        
        dis_name = gsub("dis_level_1","",dis_sheets[j])
        
        #Calculate marker score for all cells
        cell_marker_score = data.frame(rowSums(t(suerat_obj@assays$integrated@data[marker_file_pos$gene,])))
        marker_score_df = data.frame(suerat_obj$dis_level_1,cell_marker_score)
        names(marker_score_df)[2] = "cell_marker_score"
        marker_score_df$disease_target = marker_score_df$suerat_obj.dis_level_1
        marker_score_df$disease_target[marker_score_df$disease_target!=dis_name] <- "Non-target phenotype"
        marker_score_df$disease_target[marker_score_df$disease_target==dis_name] <- "Target phenotype"
        
        target_pheno_fig = ggplot(marker_score_df, aes(x=disease_target,y=cell_marker_score, fill=disease_target)) +
          #stat_compare_means(comparisons = my_comparisons) +
          geom_violin(aes(fill=disease_target), width=1) +
          geom_jitter(width = 0.15, size=0.1) + 
          geom_boxplot(outlier.size=0, width=0.2, color="black") +
          ggtitle(paste(celltype_file_shortname[i],dis_name,"vs all others marker score")) + 
          xlab("") + ylab("cell marker score") + 
          theme(axis.text.x=element_text(size=25,color='black'),
                axis.text.y=element_text(size=25,color='black'),
                axis.title=element_text(size=25),
                axis.ticks.x = element_line(size=0),
                strip.text.x = element_text(size = 0),
                legend.position="none",
                title = element_text(size = 25, face="bold")) +
          theme(
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            element_line(colour = 'black', size = 0.3),
            # Change axis line
            axis.line = element_line(colour = "black"),
            strip.background = element_rect(colour="white",
                                            fill="white")
          ) +
          scale_fill_manual(values = c("Non-target phenotype" = "blue",
                                       "Target phenotype" = "red")) + theme(legend.position="none")
        
        png(paste("/c4/home/mtaylor4/derm/combined_data/results/cluster_pca_analysis/marker_score_dim_analysis/",
                  celltype_file_shortname[i],"_",dis_name,"_marker_score_dist.png",sep=""),
            width=8,height=8,units="in",res=300)
        print(target_pheno_fig)
        dev.off()
        
        #-------------------#
        #Now array all cells along this dimensions and drop and recalculate mean differences
        target_markerscore_df = marker_score_df[marker_score_df$disease_target %in% "Target phenotype", ]
        target_markerscore_df = target_markerscore_df[order(target_markerscore_df$cell_marker_score),]
        
        non_target_markerscore_df = marker_score_df[marker_score_df$disease_target %in% "Non-target phenotype", ]
        non_target_markerscore_mean = mean(non_target_markerscore_df$cell_marker_score)
        
        #Now sequentially drop cells until I get to the end and calculate mean difference
        celldrop_df_mean = do.call(rbind, lapply(1:nrow(target_markerscore_df), function(u){
          #u=1
          diff_mean = mean(target_markerscore_df$cell_marker_score[u:nrow(target_markerscore_df)]) - non_target_markerscore_mean
          data.frame(u,row.names(target_markerscore_df)[u],diff_mean)
        }))
        names(celldrop_df_mean)[2] = "cell_id"
        
        #Now calculate the 1st derivative of these monotonic durves
        first_deriv_df = do.call(rbind, lapply(1:(nrow(celldrop_df_mean)-1), function(u){
          #u=1
          slope = celldrop_df_mean$diff_mean[u+1] - celldrop_df_mean$diff_mean[u]
          data.frame(u, slope)
        }))
        
        # first_deriv_df = subset(second_deriv_df, slope>0.01)
        
        #Do no consider the last 10% of cells
        ten_percent_cell_count = round(nrow(first_deriv_df) *0.1)
        first_deriv_df_sub = first_deriv_df [1:(nrow(first_deriv_df)-ten_percent_cell_count), ]
        
        #Now calculate the 2nd derivative of these monotonic durves
        second_deriv_df = do.call(rbind, lapply(1:(nrow(first_deriv_df_sub)-1), function(u){
          #u=2
          slope = (first_deriv_df_sub$slope[u+1] - first_deriv_df_sub$slope[u])
          data.frame(u, slope)
        }))
        
        second_deriv_df = second_deriv_df[order(-second_deriv_df$slope),]
        
        marker_score_fig = ggplot(data=celldrop_df_mean,aes(x=u,y=diff_mean)) +
          geom_point() +
          geom_line() +
          theme_classic() +
          geom_vline(xintercept = second_deriv_df$u[1], 
                     color = "red", size=1) +
          ggtitle(paste(celltype_file_shortname[i],dis_name,"marker score-sorted\nsignal optimization")) + 
          xlab("Cells dropped") + ylab("Marker score difference") + 
          theme(axis.text.x=element_text(size=25,color='black'),
                axis.text.y=element_text(size=25,color='black'),
                axis.title=element_text(size=25),
                axis.ticks.x = element_line(size=0),
                strip.text.x = element_text(size = 0),
                legend.position="none",
                title = element_text(size = 25, face="bold"))
        
        png(paste("/c4/home/mtaylor4/derm/combined_data/results/cluster_pca_analysis/marker_score_dim_analysis/",
                  celltype_file_shortname[i],"_",dis_name,"_signal_cells.png",sep=""),
            width=8,height=8,units="in",res=300)
        print(marker_score_fig)
        dev.off()
        
        #Get cells with biggest 2nd derivative not considering the last 10% of cells   
        contributory_cells = subset(celldrop_df_mean, u>second_deriv_df$u[1])
        
        contributory_cells
        
      }, error=function(e){})
      
    }))
    
    #Now subset to these that are top pathalogic cells
    #Make cell ID a column in the metadata
    suerat_obj@meta.data$cell_id = row.names(suerat_obj@meta.data)
    
    # Set identity classes to an existing column in meta data
    Idents(object = suerat_obj) <- "cell_id"
    suerat_obj_sub <- subset(suerat_obj, subset = cell_id %in% go_thru_dis_markers$cell_id ) #subset to top pathologic cells
    
    suerat_obj_sub <- FindVariableFeatures(suerat_obj_sub, selection.method = "vst", nfeatures = 5000)
    all.genes <- rownames(suerat_obj_sub)
    suerat_obj_sub <- ScaleData(suerat_obj_sub, features = all.genes)
    suerat_obj_sub <- RunPCA(suerat_obj_sub, features = VariableFeatures(object = suerat_obj_sub))
    suerat_obj_sub <- RunUMAP(suerat_obj_sub, dims=1:10)
    
    #Plot top pathologic cells
    metadat_col_index = match(c("dis_level_1"),names(suerat_obj_sub@meta.data))
    colourCount = length(unique(suerat_obj@meta.data[,metadat_col_index[k]]))
    
    #Fig: top pathologic cells
    p1 <- DimPlot(suerat_obj_sub, reduction = "umap", group.by = names(suerat_obj@meta.data)[metadat_col_index[k]],
                  cols=getPalette(colourCount),
                  raster=FALSE) +
      ggtitle( paste(celltype_file_shortname[i], names(suerat_obj@meta.data)[metadat_col_index[k]],
                     " top pathologic cells",sep=" ") ) + xlab("") + ylab("") +
      theme(axis.text=element_text(size=15,color='black'),
            axis.title=element_text(size=0),
            title=element_text(size=15),
            axis.ticks = element_line(size=0))
    
    
    
    #Fig: all pathologic cells
    suerat_obj <- RunPCA(suerat_obj, features = VariableFeatures(object = suerat_obj_sub))
    suerat_obj <- RunUMAP(suerat_obj, dims=1:10)
    p2 <- DimPlot(suerat_obj, reduction = "umap", group.by = names(suerat_obj@meta.data)[metadat_col_index[k]],
                  cols=getPalette(colourCount),
                  raster=FALSE) +
      ggtitle( paste(celltype_file_shortname[i], names(suerat_obj@meta.data)[metadat_col_index[k]],
                     " all pathologic cells",sep=" ") ) + xlab("") + ylab("") +
      theme(axis.text=element_text(size=15,color='black'),
            axis.title=element_text(size=0),
            title=element_text(size=15),
            axis.ticks = element_line(size=0))
    
    
    png(paste("/c4/home/mtaylor4/derm/combined_data/results/cluster_pca_analysis/marker_score_dim_analysis/",
              celltype_file_shortname[i],"_","_refined_cell_umap.png",sep=""),
        width=16,height=8,units="in",res=300)
    print(multiplot(p2,p1,cols=2))
    dev.off()
    
  }, error=function(e){})
  
}# end cell type loop  

#scp results to my machine
scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/combined_data/results/cluster_pca_analysis/marker_score_dim_analysis/ /Users/marktaylor/Desktop/Projects/Cho\ derm/pan-rash/Results/cluster_pca_analysis


scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/combined_data/results/cluster_pca_analysis/marker_score_dim_analysis/9_* /Users/marktaylor/Desktop/Projects/Cho\ derm/pan-rash/Results/cluster_pca_analysis/marker_score_dim_analysis



#Scrape summary f

library(rvest)

sample_paths = list.dirs("/c4/home/mtaylor4/derm/internal_data/cellranger_output/", recursive=F, full.names = T)
sample_path_names = list.dirs("/c4/home/mtaylor4/derm/internal_data/cellranger_output/", recursive=F, full.names = F)

#Copyt all summaries into the same directory

for(i in 1:length(sample_paths)) {
  tryCatch({ #suppress error
    
    #i=1
    print(i)
    file.copy(from = paste(sample_paths[i],"/outs/web_summary.html",sep=""), 
              to = paste("/c4/home/mtaylor4/derm/internal_data/cellranger_websummaries/", sample_path_names[i], "__web_summary.html",sep="" ),
              overwrite=T
    )
  }, error=function(e){})
}

scp -r mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/internal_data/cellranger_websummaries/ /Users/marktaylor/Desktop/Projects/Cho\ derm/Case\ Report\ PRP/cellranger/web_summaries


#------------------------------------------------------------------------------------------------------------#


#Web scraping code
do.call(rbind, lapply(1:length(sample_paths), function(i){
  
  #i=1
  page <- read_html(paste(sample_paths[i],"/outs/web_summary.html",sep=""))
  
  estimated_cells <- page %>% html_node("div:contains('estimatedNumberOfCellsValue') + div")
  
  
  sewd("")
  
  #Convert html to text
  writeLines(page, paste(sample_paths[i],"/outs/web_summary.txt",sep=""))  # Replace "output.txt" with the desired file name
  
  
  # Extract text from specific HTML elements using CSS selectors
  text_elements <- page %>% html_nodes("p")  # Adjust the CSS selector as needed
  
  target_element <- as.vector(page %>% html_nodes(":contains('Estimated Number of Cells')"))
  
  
  text_content <- html_text(page)
  
  str(text_content)
  
  html_node(page, "#Sequencing")
  
  
  
}))



# Read the HTML file
page <- read_html(file_path)

# Extract text from specific HTML elements using CSS selectors
text_elements <- page %>% html_nodes("p")  # Adjust the CSS selector as needed



#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#Re-analysis of Shao et al: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183134
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#

# BiocManager::install("oligo")
# BiocManager::install("pd.hugene.2.1.st")

library(oligo)
library(pd.hugene.2.1.st)
library(limma)

# Replace "your_data_folder" with the path to the folder containing your .CEL files
data_folder <- "/c4/home/mtaylor4/derm/PRP_case_report/bulk_data/shao_et_al/GSE183134_RAW/"

# Read .CEL files
celFiles <- list.files(data_folder, pattern = ".CEL", full.names = TRUE)
affyBatch <- read.celfiles(celFiles)

# Background correction (RMA background correction)
affyBatch <- rma(affyBatch)

# Find matching probes
probe_id = read.csv("/c4/home/mtaylor4/derm/PRP_case_report/bulk_data/shao_et_al/HuGene-2_1-st-v1.na36.hg19.probeset.csv",
                    skip=22)

#Get new column of gene IDs

head(probe_id)

gene_symbols_df = do.call(rbind, lapply(1:nrow(probe_id), function(i){
  tryCatch({ #suppress error
    
    #i=1
    print(i)
    gene_symbols = data.frame(strsplit(as.character(data.frame(strsplit(probe_id$gene_assignment[i]," // "))[2,]) , " /// "))[1,]
    data.frame(gene_symbols, probe_id$probeset_id[i])
  }, error=function(e){})
  
}))


gene_probe_match = gene_symbols_df[grep("IL17",gene_symbols_df$gene_symbols),]

gene_probe_match$probe_id.probeset_id.i. %in% row.names(affyBatch@assayData$exprs)



#-------------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#
#Tissue specificity

#----------------------------------------------#
#1. Download Reynolds et al AD, PV, HC data: https://developmental.cellatlas.io/diseased-skin

#----------------------------------------------#
#2. Convert Renolds et al h5ad to Seurat object
# remotes::install_github("pmbio/MuDataSeurat")
library(MuDataSeurat)

bm <- ReadH5AD("/c4/home/mtaylor4/derm/external_data/processed_data/h5ad/human_pso_ad_healthy.h5ad")
saveRDS(bm,"/c4/home/mtaylor4/derm/external_data/processed_data/seurat_objects/reynold_2021.rds")

#Transfer to my machine
scp mtaylor4@c4-log1.ucsf.edu:/c4/home/mtaylor4/derm/external_data/processed_data/seurat_objects/reynold_2021.rds /Users/marktaylor/Desktop/Projects/Cho\ derm/ad_pv_project/data/punch_bx/scRNAseq/Reynolds_et_al

#----------------------------------------------#
#3. Normalize data
bm = readRDS("/c4/home/mtaylor4/derm/external_data/processed_data/seurat_objects/reynold_2021.rds")

library(Seurat)
bm <- NormalizeData(bm)

bm@assays$RNA$data[1:10,1:10]

data.frame(table(bm@meta.data$final_clustering))


suerat_obj_sub = subset(bm, Status == c("Eczema","Psoriasis") & )

#----------------------------------------------#
#4.For any given signature find which cell class it's most highly expressed in

head(bm@meta.data)
ad_bulk_markers = read.table("/c4/home/mtaylor4/derm/ad_pv_project/data/ad_pv_markers/ad_genes_he_et_al.txt",header=T)
pv_bulk_markers = read.table("/c4/home/mtaylor4/derm/ad_pv_project/data/ad_pv_markers/pv_genes_he_et_al.txt",header=T)

head(ad_bulk_markers)
head(pv_bulk_markers)

#Loop through diseases
unique(head(bm@meta.data$Status))
unique(head(bm@meta.data$Site))
unique(head(bm@meta.data$Tissue))

sc_dis =  c("Eczema","Psoriasis")


library(dplyr)

for(i in 1:length(sc_dis)) {
  
  #i=1
  
  #Subset to a disease
  suerat_obj_sub = subset(bm, Status == sc_dis[i])
  
  #Set dis markers as either AD or PV
  if(sc_dis[i]=="Eczema") {dis_markers = ad_bulk_markers} else {dis_markers = pv_bulk_markers}
  
  #Loop through dis markers and calculate means
  #I don't understand thow this
  do.call(rbind, lapply(1:nrow(dis_markers), function(j){
    #j=1
    
    #Get dataframe of celltype and normalized expression
    df = data.frame(suerat_obj_sub@meta.data$final_clustering, suerat_obj_sub@assays$RNA@data[dis_markers[j,1],] )
    names(df) = c("celltype","expression") 
    
    #Calculate geometric means for each celltype
    geometric_means <-t(data.frame( df %>%
                                      group_by(celltype) %>%
                                      summarize(geometric_mean = exp(mean(log(expression))))
    ))
    data.frame(geometric_means)
    #Add gene name
    names(geometric_means) = geometric_means[1,]
    geometric_means[2,] = as.numeric(geometric_means[2,])
    geometric_means$gene = dis_markers[j,1]
    data.frame(geometric_means)
    
    geometric_means$undiffkc_th = geometric_means$
      
      
      data.frame(dis_markers[j,1],geometric_means)
    
    #Get
    
  }))
  
}



