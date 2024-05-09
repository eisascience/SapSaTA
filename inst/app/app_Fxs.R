


ImputeSDA2SerV2 <- function(SerObj, sda_loadings, keepComps=NULL, sdaObjID="", 
         plot=F, doAsinh = T, MakeReduc=T, assay = "RNA"){
  
  genes_overlap = intersect(rownames(SerObj), colnames(sda_loadings))
  
  if(is.null(keepComps)) keepComps = 1:nrow(sda_loadings)
  
  if(!assay %in% names(SerObj@assays)){
    print(paste0("Available Assays: ", paste0(names(SerObj@assays), collapse = ", ") )  )
    stop("Assay not in object ")
  }
  
  sda_scores = Matrix::as.matrix(Matrix::t(sda_loadings[keepComps, genes_overlap] %*% SerObj@assays[[assay]]$data[genes_overlap, ]))
  colnames(sda_scores) = paste0("sda.", sdaObjID, ".V", keepComps)
  
  if(doAsinh) {
    SerObj = AddMetaData(SerObj, as.data.frame(asinh(sda_scores)))
  } else {
    SerObj = AddMetaData(SerObj, as.data.frame(sda_scores))
  }
  
  if(plot){
    print(Seurat::FeaturePlot(SerObj, features = colnames(sda_scores), order = T) & 
            ggplot2::scale_colour_gradientn(colours = c('navy', 'dodgerblue', 'white', 'gold', 'red')))
  }
  
  
  
  if(MakeReduc){
    CommonLoadingMat = t(as.matrix(sda_loadings)[keepComps, genes_overlap])
    colnames(CommonLoadingMat) = paste0("sda.", sdaObjID, ".V", keepComps)
    
    
    
    SerObj = SDAScoreMeta2Reduction(SerObj = SerObj,
                                    sdaComps = colnames(CommonLoadingMat), 
                                    loadingMat = as.matrix(CommonLoadingMat),
                                    reduction.key =  paste0("sda", sdaObjID, "_"), 
                                    includeLoading = T,
                                    assayName = "RNA", reduction.name = paste0("sda_", sdaObjID) )
    
  }
  
  
  
  
  
  return(SerObj)
  
}




plot_loadings_coordinates <- function(SDARedDataLS,
         reduction,
         # mart = NULL,
         genes = NULL,
         dimN,
         highlight_genes = NULL, 
         TopNpos = 10, TopNneg=10,
         data_set = "hsapiens_gene_ensembl", #mmulatta_gene_ensembl
         invertWeights=F, includeUnMapped = T, geneLocPath=NULL ) {
  
  # library(biomaRt)
  # library(ggplot2)
  # library(ggrepel)
  
  # Use the Ensembl Mart for human genes
  # if(is.null(mart)) mart <- useMart("ensembl", data_set)
  
  # Get gene coordinates
  # if(is.null(genes)) genes <- getBM(attributes = c("chromosome_name", "start_position", "end_position"), mart = mart)
  
  # Calculate chromosome lengths
  chromosomes <- unique(genes$chromosome_name)
  
  chromosome_length_red <- sapply(chromosomes, function(chr) {
    chr_genes <- subset(genes, chromosome_name == chr)
    max_end <- max(chr_genes$end_position)
    return(max_end)
  })
  
  chromosome_length_red = chromosome_length_red[names(chromosome_length_red) %in% c(1:22, "X", "Y", "MT")]
  chromosome_lengthsDF = data.frame(chromosome = names(chromosome_length_red),
                                    length = as.numeric(chromosome_length_red))
  
  
  
  Genes2Map = colnames(SDARedDataLS$loadings[[reduction]]$loadings)
  
  
  
  if(!is.null(geneLocPath)) {
    
    if(file.exists(geneLocPath)){
      gene_locations = readRDS(geneLocPath)
      
    } 
    # if(!file.exists(geneLocPath)){
    #   print("file does not exist, downloading new results")
    #   
    #   gene_locations <- SDAtools:::get.location(
    #     gene.symbols = Genes2Map,
    #     data_set = data_set,
    #     gene_name = "external_gene_name"
    #   )
    #   
    #   saveRDS(gene_locations, geneLocPath)
    # }
  } else{
    # if(is.null(geneLocPath)){
    # gene_locations <- get.location(
    #   gene.symbols = Genes2Map,
    #   data_set = data_set,
    #   gene_name = "external_gene_name"
    # )
    # } 
  }
  
  
  gene_locations[is.na(gene_locations$start_position), ]$chromosome_name = "?"
  gene_locations[is.na(gene_locations$start_position), ]$start_position = 1
  
  cl = lapply(unique(gene_locations$chromosome_name), function(xN){
    data.frame(xN, round( (max(subset(gene_locations, chromosome_name %in% xN)$start_position) - min(subset(gene_locations, chromosome_name %in% xN)$start_position))/2))
  }) %>% data.table::rbindlist() %>% as.data.frame()
  colnames(cl) = colnames=c("chr", "center")
  
  if(includeUnMapped){
    if(sum(grepl("^LOC", Genes2Map)) > 0){
      LOCgenes = Genes2Map[grepl("^LOC", Genes2Map)]
      
      geneLocPath_fix = gsub(".rds", "_LOCfix.rds", geneLocPath)
      
      
      
      gene_locations = rbind(gene_locations,
                             data.frame(gene_symbol = LOCgenes, 
                                        chromosome_name = rep("?", length(LOCgenes)),
                                        start_position = rep(1, length(LOCgenes))))
      
    }
  }
  
  
  temp <- merge(gene_locations, chromosome_lengthsDF, by.x = "chromosome_name", by.y = "chromosome", all.x = TRUE) %>% 
    as.data.frame()
  temp$genomic_position <- temp$start_position + temp$length / 2
  
  
  component = SDARedDataLS$loadings[[reduction]]$loadings[dimN,]
  
  if(invertWeights){
    component = component * -1
  }
  
  # Subset genes with weights
  label_data <- data.frame(
    gene_symbol = names(component),
    loading = component
  )
  
  # Merge with genomic positions
  label_data <- merge(label_data, temp, by = "gene_symbol", all.x = TRUE)
  
  
  label_data[is.na(label_data$start_position), ]$chromosome_name = "?"
  label_data[is.na(label_data$start_position), ]$length = 3
  if(includeUnMapped & sum(is.na(label_data$start_position))>0) label_data[is.na(label_data$start_position), ]$genomic_position = max(label_data$genomic_position, na.rm = T) + 2
  label_data[is.na(label_data$start_position), ]$start_position = 1
  
  if(sum(is.na(label_data$length))>0) label_data[is.na(label_data$length), ]$length = 3
  if(includeUnMapped & sum(is.na(label_data$genomic_position))>0) label_data[is.na(label_data$genomic_position), ]$genomic_position = max(label_data$genomic_position, na.rm = T) + 2
  
  
  label_data$gene_symbol_show = ""
  label_data$gene_symbol_show[order(label_data$loading, decreasing = TRUE)[1:TopNpos]] = label_data$gene_symbol[order(label_data$loading, decreasing = TRUE)[1:TopNpos]] 
  label_data$gene_symbol_show[order(label_data$loading, decreasing = FALSE)[1:TopNneg]] = label_data$gene_symbol[order(label_data$loading, decreasing = FALSE)[1:TopNneg]] 
  
  
  
  
  P <- ggplot(label_data, aes(genomic_position, loading, size = abs(loading)^2)) + 
    geom_point(stroke = 0, aes(alpha = (abs(loading))^0.7, 
                               color = ifelse(abs(loading)>.05, "cornflowerblue", "grey") ) ) + 
    scale_color_manual(values =  c("cornflowerblue", "grey"))  + 
    # scale_x_continuous(breaks = cl$center, 
    #                    labels = cl$chr, minor_breaks = NULL) +
    # scale_color_manual(values = rep_len("black", length(unique(label_data$chromosome_name))+1))+
    # scale_colour_manual(values = c(rep_len(c("black", "cornflowerblue"), 
    #                                        length(unique(label_data$chromosome_name))), "grey")) +
    xlab("Genomic Coordinate") + ylab("Weight")  +
    geom_label_repel(aes(label = gene_symbol_show), max.overlaps = max(c(TopNneg,TopNpos ))*2+1,
                     size = 3, box.padding = unit(0.5, "lines"), 
                     point.padding = unit(0.1, "lines"), force = 10, segment.size = 0.2, segment.color = "blue") +
    theme_minimal() + theme(legend.position = "none"); P
  
  # Highlight genes if necessary
  if (!is.null(highlight_genes)) {
    P <- P + geom_point(data = label_data[label_data$gene_symbol %in% highlight_genes, ], color = "red")
  }
  
  return(P)
}



getTopLoadedSDAgenes <- function(envv, topNfeats){
  
  if(is.null(envv$Top_loaded_genes)){
    
    Top_Pos = lapply(names(envv$SDARedDataLS$loadings), function(xN){
      # xN =1
      tempDF = lapply(1:envv$SDARedDataLS$loadings[[xN]]$ncomp, function(cN){
        x = envv$SDARedDataLS$loadings[[xN]]$loadings[cN,]
        c(sort(x, decreasing = T)[1:topNfeats] %>% names())
      }) %>% data.frame()
      
      colnames(tempDF) = paste0("sda.",xN,".V", 1:envv$SDARedDataLS$loadings[[xN]]$ncomp)
      tempDF
    }) %>% data.frame()
    
    
    Top_Neg =  lapply(names(envv$SDARedDataLS$loadings), function(xN){
      # xN =1
      tempDF = lapply(1:envv$SDARedDataLS$loadings[[xN]]$ncomp, function(cN){
        x = envv$SDARedDataLS$loadings[[xN]]$loadings[cN,]
        c(sort(x, decreasing = F)[1:topNfeats] %>% names())
      }) %>% data.frame()
      
      colnames(tempDF) = paste0("sda.",xN,".V", 1:envv$SDARedDataLS$loadings[[xN]]$ncomp)
      tempDF
    }) %>% data.frame()
    
    
    
    envv$Top_loaded_genes = list(Top_Pos = Top_Pos,
                                 Top_Neg = Top_Neg)
    
  }
  
  return(envv)

  
}


