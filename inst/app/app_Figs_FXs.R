clean_GeneSet <- function(GeneSet){
  
  if(length(grep(",", GeneSet)) == 0){
    
    if(length(grep('"', GeneSet)) + length(grep("'", GeneSet))>0) {
      GeneSet <- unlist(strsplit(gsub("'", '', gsub('"', '', GeneSet)), " "))
    } else {
      GeneSet <- unlist(strsplit(GeneSet, " "))
    }
    
    
  } else {
    GeneSet <- (unlist(strsplit(gsub(" ", "", gsub("'", '', gsub('"', '', GeneSet))), ",")))
    
  }
  
  GeneSet
}

# Render_GeneExpr_MoDSTA_UMAPintg <- function(envv, input, GeneSet){
#   
#   GeneSet = clean_GeneSet(GeneSet)
#   
#   GeneSet = GeneSet[GeneSet %in% rownames(envv$MoDSTA)]
#   
#   print(GeneSet)
#   
#   if(is.null(envv$MoDSTA) | length(GeneSet) < 1){
#     plot(x=0, y=0, main="Load MoDSTA dataset")
#   } else {
#     
#     if(length(GeneSet)==1) {
#       
#       FeaturePlot(envv$MoDSTA , reduction = "umapSup.harmony", features = GeneSet[1], 
#                   # raster = T, raster.dpi = c(800, 800),
#                   order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
#         theme_classic(base_size = 14) +
#         theme(axis.line = element_blank(),
#               axis.text.x = element_blank(),
#               axis.text.y = element_blank(),
#               axis.ticks = element_blank(),
#               axis.title = element_blank() #,plot.title = element_blank()
#         )  &
#         ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
#       
#     } else{
#       
#       if(input$GeneExprOpr == "sum"){
#         envv$MoDSTA$MultiGeneExpr = colSums(GetAssayData(envv$MoDSTA, assay = "RNA", layer = "data")[GeneSet, ])
#         
#       } else {
#         envv$MoDSTA$MultiGeneExpr = colMeans(GetAssayData(envv$MoDSTA, assay = "RNA", layer = "data")[GeneSet, ])
#         
#       }
#       
#       
#       FeaturePlot(envv$MoDSTA , reduction = "umapSup.harmony", features = "MultiGeneExpr",
#                   # raster = T, raster.dpi = c(800, 800),
#                   order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
#         theme_classic(base_size = 14) +
#         theme(axis.line = element_blank(),
#               axis.text.x = element_blank(),
#               axis.text.y = element_blank(),
#               axis.ticks = element_blank(),
#               axis.title = element_blank() #,plot.title = element_blank()
#         ) + ggtitle(ifelse(input$GeneExprOpr == "sum", "Summed Expression", "Mean of Expression Set"))  &
#         ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
#       
#       
#     }
#     
#     
#     
#   }
#   
#   
#   
# }
# 
# Threshold_GeneExpr_MoDSTA_UMAPintg <- function(envv, input, GeneSet){
#   
#   GeneSet = clean_GeneSet(GeneSet)
#   
#   
#   GeneSet = GeneSet[GeneSet %in% rownames(envv$MoDSTA)]
#   
#   print(GeneSet)
#   
#   
#   
#   if(is.null(envv$MoDSTA) | length(GeneSet) < 1){
#     plot(x=0, y=0, main="Load STseq H1 dataset")
#   } else {
#     
#     
#     if(length(GeneSet)==1) {
#       
#       envv$MoDSTA$GeneExpr = colSums(GetAssayData(envv$MoDSTA, assay = "RNA", layer = "data")[GeneSet, ])
#       
#       
#     } else{
#       
#       if(input$GeneExprOpr == "sum"){
#         envv$MoDSTA$GeneExpr = colSums(GetAssayData(envv$MoDSTA, assay = "RNA", layer = "data")[GeneSet, ])
#         
#       } else {
#         envv$MoDSTA$GeneExpr = colMeans(GetAssayData(envv$MoDSTA, assay = "RNA", layer = "data")[GeneSet, ])
#         
#       }
#       
#       
#       
#     }
#     
#     
#     tempExprVec = envv$MoDSTA$GeneExpr
#     
#     envv$MoDSTA$ExprBL = (ifelse( tempExprVec > max( tempExprVec ) * input$ExprThresh,
#                                                "Above", "Below"))
#     
#     
#     DimPlot(envv$MoDSTA, reduction = "umapSup.harmony",
#                    group.by = "ExprBL",
#     ) +facet_wrap(~ExprBL, ncol=4) +
#       theme_classic(base_size = 14) +
#       guides(colour = guide_legend(override.aes = list(size = 1.5, alpha=1), ncol=1)) +
#       theme(legend.position = "none",
#             axis.line = element_blank(),
#             axis.text.x = element_blank(),
#             axis.text.y = element_blank(),
#             axis.ticks = element_blank(),
#             axis.title = element_blank() #,plot.title = element_blank()
#       ) + ggtitle(paste0("Threshold cut above: ", round(max(tempExprVec) * input$ExprThresh,3)))
#     
#   }
#   
#   
# }
# 
# Render_GeneExpr_STseqH1_UMAP <- function(envv, input, GeneSet){
#   
#   GeneSet = clean_GeneSet(GeneSet)
#   
#   
#   GeneSet = GeneSet[GeneSet %in% rownames(envv$MSTseqCells1)]
#   
#   print(GeneSet)
#   
#   if(is.null(envv$MSTseqCells1) | length(GeneSet) < 1){
#     plot(x=0, y=0, main="Load STseq H1 dataset")
#   } else {
#     
#     if(length(GeneSet)==1) {
#       
#       FeaturePlot(envv$MSTseqCells1 , reduction = "umap", features = GeneSet[1], 
#                   # raster = T, raster.dpi = c(800, 800),
#                   order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
#         theme_classic(base_size = 14) +
#         theme(axis.line = element_blank(),
#               axis.text.x = element_blank(),
#               axis.text.y = element_blank(),
#               axis.ticks = element_blank(),
#               axis.title = element_blank() #,plot.title = element_blank()
#         )  &
#         ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
#       
#     } else{
#       
#       if(input$GeneExprOpr == "sum"){
#         envv$MSTseqCells1$MultiGeneExpr = colSums(GetAssayData(envv$MSTseqCells1, assay = "SCT", layer = "data")[GeneSet, ])
#         
#       } else {
#         envv$MSTseqCells1$MultiGeneExpr = colMeans(GetAssayData(envv$MSTseqCells1, assay = "SCT", layer = "data")[GeneSet, ])
#         
#       }
#       
#       
#       FeaturePlot(envv$MSTseqCells1 , reduction = "umap", features = "MultiGeneExpr",
#                   # raster = T, raster.dpi = c(800, 800),
#                   order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
#         theme_classic(base_size = 14) +
#         theme(axis.line = element_blank(),
#               axis.text.x = element_blank(),
#               axis.text.y = element_blank(),
#               axis.ticks = element_blank(),
#               axis.title = element_blank() #,plot.title = element_blank()
#         ) + ggtitle(ifelse(input$GeneExprOpr == "sum", "Summed Expression", "Mean of Expression Set"))  &
#         ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
#       
#     }
#     
#     
#     
#   }
# }

Render_GeneExpr_STseqH1_Spatial <- function(envv, input, GeneSet){
  
  GeneSet = clean_GeneSet(GeneSet)
  
  GeneSet = GeneSet[GeneSet %in% rownames(envv$MSTseqCells1)]
  
  print(GeneSet)
  
  if(is.null(envv$MSlideSeqV1CellsH1) | length(GeneSet) < 1){
    plot(x=0, y=0, main="Load STseq H1 dataset")
  } else {
    
    if(length(GeneSet)==1) {
      
      gg1 = SpatialFeaturePlot(envv$MSTseqCells1 , features = GeneSet[1], 
                               # raster = T, raster.dpi = c(800, 800),
                               max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    } else{
      
      if(input$GeneExprOpr == "sum"){
        envv$MSTseqCells1$MultiGeneExpr = colSums(GetAssayData(envv$MSTseqCells1, assay = "SCT", layer = "data")[GeneSet, ])
        
      } else {
        envv$MSTseqCells1$MultiGeneExpr = colMeans(GetAssayData(envv$MSTseqCells1, assay = "SCT", layer = "data")[GeneSet, ])
        
      }
      
      
      gg1 = SpatialFeaturePlot(envv$MSTseqCells1 , features = "MultiGeneExpr",
                               # raster = T, raster.dpi = c(800, 800),
                               max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    }
    if(is.null(ranges$GeneExpr_STseqH1_Spatial_x)){
      gg1
    } else {
      xrange = ggplot_build(gg1)$layout$panel_params[[1]]$x.range
      yrange = ggplot_build(gg1)$layout$panel_params[[1]]$y.range
      
      
      gg1 + coord_cartesian(xlim = c((xrange[1] + (abs(xrange[1] - xrange[2])*(min(ranges$GeneExpr_STseqH1_Spatial_x) - 0.2*min(ranges$GeneExpr_STseqH1_Spatial_x)) )), 
                                     (xrange[1] + (abs(xrange[1] - xrange[2])*(max(ranges$GeneExpr_STseqH1_Spatial_x) + 0.2*max(ranges$GeneExpr_STseqH1_Spatial_x)) ))), 
                            ylim = c((yrange[1] + (abs(yrange[1] - yrange[2])*(min(ranges$GeneExpr_STseqH1_Spatial_x) - 0.2*min(ranges$GeneExpr_STseqH1_Spatial_x)) )), 
                                     (yrange[1] + (abs(yrange[1] - yrange[2])*(max(ranges$GeneExpr_STseqH1_Spatial_x) + 0.2*max(ranges$GeneExpr_STseqH1_Spatial_x)) )) ), expand = T)
      
    }
    
    
  }
  
}

Render_GeneExpr_SlideSeqV1H1_UMAP <- function(envv, input, GeneSet){
  
  GeneSet = clean_GeneSet(GeneSet)
  
  
  GeneSet = GeneSet[GeneSet %in% rownames(envv$MSlideSeqV1CellsH1)]
  
  print(GeneSet)
  
  if(is.null(envv$MSlideSeqV1CellsH1) | length(GeneSet) < 1){
    plot(x=0, y=0, main="Load STseq H1 dataset")
  } else {
    
    if(length(GeneSet)==1) {
      
      FeaturePlot(envv$MSlideSeqV1CellsH1 , reduction = "umap", features = GeneSet[1], 
                  # raster = T, raster.dpi = c(800, 800),
                  order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    } else{
      
      if(input$GeneExprOpr == "sum"){
        envv$MSlideSeqV1CellsH1$MultiGeneExpr = colSums(GetAssayData(envv$MSlideSeqV1CellsH1, assay = "SCT", layer = "data")[GeneSet, ])
        
      } else {
        envv$MSlideSeqV1CellsH1$MultiGeneExpr = colMeans(GetAssayData(envv$MSlideSeqV1CellsH1, assay = "SCT", layer = "data")[GeneSet, ])
        
      }
      
      
      FeaturePlot(envv$MSlideSeqV1CellsH1 , reduction = "umap", features = "MultiGeneExpr",
                  # raster = T, raster.dpi = c(800, 800),
                  order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        ) + ggtitle(ifelse(input$GeneExprOpr == "sum", "Summed Expression", "Mean of Expression Set"))  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    }
    
    
    
  }
}

Render_GeneExpr_SlideSeqV1H1_Spatial <- function(envv, input, GeneSet){
  
  GeneSet = clean_GeneSet(GeneSet)
  
  
  GeneSet = GeneSet[GeneSet %in% rownames(envv$MSlideSeqV1CellsH1)]
  
  print(GeneSet)
  
  
  
  if(is.null(envv$MSlideSeqV1CellsH1) | length(GeneSet) < 1){
    plot(x=0, y=0, main="Load STseq H1 dataset")
  } else {
    
    
    if(length(GeneSet)==1) {
      
    gg1=  SpatialFeaturePlot(envv$MSlideSeqV1CellsH1 , features = GeneSet[1], 
                         # raster = T, raster.dpi = c(800, 800),
                         max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    } else{
      
      if(input$GeneExprOpr == "sum"){
        envv$MSlideSeqV1CellsH1$MultiGeneExpr = colSums(GetAssayData(envv$MSlideSeqV1CellsH1, assay = "SCT", layer = "data")[GeneSet, ])
        
      } else {
        envv$MSlideSeqV1CellsH1$MultiGeneExpr = colMeans(GetAssayData(envv$MSlideSeqV1CellsH1, assay = "SCT", layer = "data")[GeneSet, ])
        
      } 
      
      
     gg1= SpatialFeaturePlot(envv$MSlideSeqV1CellsH1 , features = "MultiGeneExpr",
                         # raster = T, raster.dpi = c(800, 800),
                         max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    }
    
    
    # print(paste0("x: ", ranges$x*10000, " y: ", ranges$y*10000))
    
    if(is.null(ranges$GeneExpr_SlideSeqV1H1_Spatial_x)){
      gg1
    } else {
      xrange = ggplot_build(gg1)$layout$panel_params[[1]]$x.range
      yrange = ggplot_build(gg1)$layout$panel_params[[1]]$y.range
      
      
      gg1 + coord_cartesian(xlim = c((xrange[1] + (abs(xrange[1] - xrange[2])*(min(ranges$GeneExpr_SlideSeqV1H1_Spatial_x) - 0.2*min(ranges$GeneExpr_SlideSeqV1H1_Spatial_x)) )), 
                                     (xrange[1] + (abs(xrange[1] - xrange[2])*(max(ranges$GeneExpr_SlideSeqV1H1_Spatial_x) + 0.2*max(ranges$GeneExpr_SlideSeqV1H1_Spatial_x)) ))), 
                            ylim = c((yrange[1] + (abs(yrange[1] - yrange[2])*(min(ranges$GeneExpr_SlideSeqV1H1_Spatial_x) - 0.2*min(ranges$GeneExpr_SlideSeqV1H1_Spatial_x)) )), 
                                     (yrange[1] + (abs(yrange[1] - yrange[2])*(max(ranges$GeneExpr_SlideSeqV1H1_Spatial_x) + 0.2*max(ranges$GeneExpr_SlideSeqV1H1_Spatial_x)) )) ), expand = T)
      
    }
 
  }
    
}

Threshold_GeneExpr_SlideSeqV1H1_Spatial <- function(envv, input, GeneSet){
  
  GeneSet = clean_GeneSet(GeneSet)
  
  
  GeneSet = GeneSet[GeneSet %in% rownames(envv$MSlideSeqV1CellsH1)]
  
  print(GeneSet)
  
  
  
  if(is.null(envv$MSlideSeqV1CellsH1) | length(GeneSet) < 1){
    plot(x=0, y=0, main="Load STseq H1 dataset")
  } else {
    
    
    if(length(GeneSet)==1) {
      
      envv$MSlideSeqV1CellsH1$GeneExpr = colSums(GetAssayData(envv$MSlideSeqV1CellsH1, assay = "SCT", layer = "data")[GeneSet, ])
      
      
    } else{
      
      if(input$GeneExprOpr == "sum"){
        envv$MSlideSeqV1CellsH1$GeneExpr = colSums(GetAssayData(envv$MSlideSeqV1CellsH1, assay = "SCT", layer = "data")[GeneSet, ])
        
      } else {
        envv$MSlideSeqV1CellsH1$GeneExpr = colMeans(GetAssayData(envv$MSlideSeqV1CellsH1, assay = "SCT", layer = "data")[GeneSet, ])
        
      }
      
      
      
    }
    
    
    tempExprVec = envv$MSlideSeqV1CellsH1$GeneExpr
    
    envv$MSlideSeqV1CellsH1$ExprBL = (ifelse( tempExprVec > max( tempExprVec ) * input$ExprThresh,
                                               "Above", "Below"))
    
    
    
  SpatialDimPlot(envv$MSlideSeqV1CellsH1 ,
                 group.by = "ExprBL",
  ) +facet_wrap(~ExprBL, ncol=4) +
    theme_classic(base_size = 14) +
    guides(colour = guide_legend(override.aes = list(size = 1.5, alpha=1), ncol=1)) +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank() #,plot.title = element_blank()
    ) + ggtitle(paste0("Threshold cut above: ", round(max(tempExprVec) * input$ExprThresh,3)))
   
  }
 
  
}

Threshold_GeneExpr_STseqH1_Spatial<- function(envv, input, GeneSet){
  
  GeneSet = clean_GeneSet(GeneSet)
  
  
  GeneSet = GeneSet[GeneSet %in% rownames(envv$MSTseqCells1)]
  
  print(GeneSet)
  
  
  
  if(is.null(envv$MSTseqCells1) | length(GeneSet) < 1){
    plot(x=0, y=0, main="Load STseq H1 dataset")
  } else {
    
    
    if(length(GeneSet)==1) {
      
      envv$MSTseqCells1$GeneExpr = colSums(GetAssayData(envv$MSTseqCells1, assay = "SCT", layer = "data")[GeneSet, ])
      
      
    } else{
      
      if(input$GeneExprOpr == "sum"){
        envv$MSTseqCells1$GeneExpr = colSums(GetAssayData(envv$MSTseqCells1, assay = "SCT", layer = "data")[GeneSet, ])
        
      } else {
        envv$MSTseqCells1$GeneExpr = colMeans(GetAssayData(envv$MSTseqCells1, assay = "SCT", layer = "data")[GeneSet, ])
        
      }
      
      
      
    }
    
    
    tempExprVec = envv$MSTseqCells1$GeneExpr
    
    envv$MSTseqCells1$ExprBL = (ifelse( tempExprVec > max( tempExprVec ) * input$ExprThresh,
                                               "Above", "Below"))
    
    
    
    SpatialDimPlot(envv$MSTseqCells1 ,
                   group.by = "ExprBL",
    ) +facet_wrap(~ExprBL, ncol=4) +
      theme_classic(base_size = 14) +
      guides(colour = guide_legend(override.aes = list(size = 1.5, alpha=1), ncol=1)) +
      theme(legend.position = "none",
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank() #,plot.title = element_blank()
      ) + ggtitle(paste0("Threshold cut above: ", round(max(tempExprVec) * input$ExprThresh,3)))
    
  }
  
  
}

Render_GeneExpr_SlideSeqV1H2_UMAP <- function(envv, input, GeneSet){
  
  GeneSet = clean_GeneSet(GeneSet)
  
  
  GeneSet = GeneSet[GeneSet %in% rownames(envv$MSlideSeqV1CellsH2)]
  
  print(GeneSet)
  
  if(is.null(envv$MSlideSeqV1CellsH2) | length(GeneSet) < 1){
    plot(x=0, y=0, main="Load STseq H2 dataset")
  } else {
    
    if(length(GeneSet)==1) {
      
      FeaturePlot(envv$MSlideSeqV1CellsH2 , reduction = "umap", features = GeneSet[1], 
                  # raster = T, raster.dpi = c(800, 800),
                  order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    } else{
      
      if(input$GeneExprOpr == "sum"){
        envv$MSlideSeqV1CellsH2$MultiGeneExpr = colSums(GetAssayData(envv$MSlideSeqV1CellsH2, assay = "SCT", layer = "data")[GeneSet, ])
        
      } else {
        envv$MSlideSeqV1CellsH2$MultiGeneExpr = colMeans(GetAssayData(envv$MSlideSeqV1CellsH2, assay = "SCT", layer = "data")[GeneSet, ])
        
      }
      
      
      FeaturePlot(envv$MSlideSeqV1CellsH2 , reduction = "umap", features = "MultiGeneExpr",
                  # raster = T, raster.dpi = c(800, 800),
                  order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        ) + ggtitle(ifelse(input$GeneExprOpr == "sum", "Summed Expression", "Mean of Expression Set"))  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    }
    
    
    
  }
}

Render_GeneExpr_SlideSeqV1H2_Spatial <- function(envv, input, GeneSet){
  
  GeneSet = clean_GeneSet(GeneSet)
  
  
  GeneSet = GeneSet[GeneSet %in% rownames(envv$MSlideSeqV1CellsH2)]
  
  print(GeneSet)
  
  
  
  if(is.null(envv$MSlideSeqV1CellsH2) | length(GeneSet) < 1){
    plot(x=0, y=0, main="Load STseq H2 dataset")
  } else {
    
    
    if(length(GeneSet)==1) {
      
      gg1=  SpatialFeaturePlot(envv$MSlideSeqV1CellsH2 , features = GeneSet[1], 
                               # raster = T, raster.dpi = c(800, 800),
                               max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    } else{
      
      if(input$GeneExprOpr == "sum"){
        envv$MSlideSeqV1CellsH2$MultiGeneExpr = colSums(GetAssayData(envv$MSlideSeqV1CellsH2, assay = "SCT", layer = "data")[GeneSet, ])
        
      } else {
        envv$MSlideSeqV1CellsH2$MultiGeneExpr = colMeans(GetAssayData(envv$MSlideSeqV1CellsH2, assay = "SCT", layer = "data")[GeneSet, ])
        
      } 
      
      
      gg1= SpatialFeaturePlot(envv$MSlideSeqV1CellsH2 , features = "MultiGeneExpr",
                              # raster = T, raster.dpi = c(800, 800),
                              max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    }
    
    
    # print(paste0("x: ", ranges$x*10000, " y: ", ranges$y*10000))
    
    if(is.null(ranges$GeneExpr_SlideSeqV1H2_Spatial_x)){
      gg1
    } else {
      xrange = ggplot_build(gg1)$layout$panel_params[[1]]$x.range
      yrange = ggplot_build(gg1)$layout$panel_params[[1]]$y.range
      
      
      gg1 + coord_cartesian(xlim = c((xrange[1] + (abs(xrange[1] - xrange[2])*(min(ranges$GeneExpr_SlideSeqV1H2_Spatial_x) - 0.2*min(ranges$GeneExpr_SlideSeqV1H2_Spatial_x)) )), 
                                     (xrange[1] + (abs(xrange[1] - xrange[2])*(max(ranges$GeneExpr_SlideSeqV1H2_Spatial_x) + 0.2*max(ranges$GeneExpr_SlideSeqV1H2_Spatial_x)) ))), 
                            ylim = c((yrange[1] + (abs(yrange[1] - yrange[2])*(min(ranges$GeneExpr_SlideSeqV1H2_Spatial_x) - 0.2*min(ranges$GeneExpr_SlideSeqV1H2_Spatial_x)) )), 
                                     (yrange[1] + (abs(yrange[1] - yrange[2])*(max(ranges$GeneExpr_SlideSeqV1H2_Spatial_x) + 0.2*max(ranges$GeneExpr_SlideSeqV1H2_Spatial_x)) )) ), expand = T)
      
    }
    
  }
  
}

Threshold_GeneExpr_SlideSeqV1H2_Spatial <- function(envv, input, GeneSet){
  
  GeneSet = clean_GeneSet(GeneSet)
  
  
  GeneSet = GeneSet[GeneSet %in% rownames(envv$MSlideSeqV1CellsH2)]
  
  print(GeneSet)
  
  
  
  if(is.null(envv$MSlideSeqV1CellsH2) | length(GeneSet) < 1){
    plot(x=0, y=0, main="Load STseq H2 dataset")
  } else {
    
    
    if(length(GeneSet)==1) {
      
      envv$MSlideSeqV1CellsH2$GeneExpr = colSums(GetAssayData(envv$MSlideSeqV1CellsH2, assay = "SCT", layer = "data")[GeneSet, ])
      
      
    } else{
      
      if(input$GeneExprOpr == "sum"){
        envv$MSlideSeqV1CellsH2$GeneExpr = colSums(GetAssayData(envv$MSlideSeqV1CellsH2, assay = "SCT", layer = "data")[GeneSet, ])
        
      } else {
        envv$MSlideSeqV1CellsH2$GeneExpr = colMeans(GetAssayData(envv$MSlideSeqV1CellsH2, assay = "SCT", layer = "data")[GeneSet, ])
        
      }
      
      
      
    }
    
    
    tempExprVec = envv$MSlideSeqV1CellsH2$GeneExpr
    
    envv$MSlideSeqV1CellsH2$ExprBL = (ifelse( tempExprVec > max( tempExprVec ) * input$ExprThresh,
                                               "Above", "Below"))
    
    
    
    SpatialDimPlot(envv$MSlideSeqV1CellsH2 ,
                   group.by = "ExprBL",
    ) +facet_wrap(~ExprBL, ncol=4) +
      theme_classic(base_size = 14) +
      guides(colour = guide_legend(override.aes = list(size = 1.5, alpha=1), ncol=1)) +
      theme(legend.position = "none",
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank() #,plot.title = element_blank()
      ) + ggtitle(paste0("Threshold cut above: ", round(max(tempExprVec) * input$ExprThresh,3)))
    
  }
  
  
}

Threshold_GeneExpr_STseqH2_Spatial<- function(envv, input, GeneSet){
  
  GeneSet = clean_GeneSet(GeneSet)
  
  
  GeneSet = GeneSet[GeneSet %in% rownames(envv$MSTseqCells1)]
  
  print(GeneSet)
  
  
  
  if(is.null(envv$MSTseqCells1) | length(GeneSet) < 1){
    plot(x=0, y=0, main="Load STseq H2 dataset")
  } else {
    
    
    if(length(GeneSet)==1) {
      
      envv$MSTseqCells1$GeneExpr = colSums(GetAssayData(envv$MSTseqCells1, assay = "SCT", layer = "data")[GeneSet, ])
      
      
    } else{
      
      if(input$GeneExprOpr == "sum"){
        envv$MSTseqCells1$GeneExpr = colSums(GetAssayData(envv$MSTseqCells1, assay = "SCT", layer = "data")[GeneSet, ])
        
      } else {
        envv$MSTseqCells1$GeneExpr = colMeans(GetAssayData(envv$MSTseqCells1, assay = "SCT", layer = "data")[GeneSet, ])
        
      }
      
      
      
    }
    
    
    tempExprVec = envv$MSTseqCells1$GeneExpr
    
    envv$MSTseqCells1$ExprBL = (ifelse( tempExprVec > max( tempExprVec ) * input$ExprThresh,
                                        "Above", "Below"))
    
    
    
    SpatialDimPlot(envv$MSTseqCells1 ,
                   group.by = "ExprBL",
    ) +facet_wrap(~ExprBL, ncol=4) +
      theme_classic(base_size = 14) +
      guides(colour = guide_legend(override.aes = list(size = 1.5, alpha=1), ncol=1)) +
      theme(legend.position = "none",
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank() #,plot.title = element_blank()
      ) + ggtitle(paste0("Threshold cut above: ", round(max(tempExprVec) * input$ExprThresh,3)))
    
  }
  
  
}

Render_GeneExpr_SlideSeqV1MT3_UMAP <- function(envv, input, GeneSet){
  
  GeneSet = clean_GeneSet(GeneSet)
  
  
  GeneSet = GeneSet[GeneSet %in% rownames(envv$MSlideSeqV1CellsWT3)]
  
  print(GeneSet)
  
  if(is.null(envv$MSlideSeqV1CellsWT3) | length(GeneSet) < 1){
    plot(x=0, y=0, main="Load STseq MT3 dataset")
  } else {
    
    if(length(GeneSet)==1) {
      
      FeaturePlot(envv$MSlideSeqV1CellsWT3 , reduction = "umap", features = GeneSet[1], 
                  # raster = T, raster.dpi = c(800, 800),
                  order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    } else{
      
      if(input$GeneExprOpr == "sum"){
        envv$MSlideSeqV1CellsWT3$MultiGeneExpr = colSums(GetAssayData(envv$MSlideSeqV1CellsWT3, assay = "SCT", layer = "data")[GeneSet, ])
        
      } else {
        envv$MSlideSeqV1CellsWT3$MultiGeneExpr = colMeans(GetAssayData(envv$MSlideSeqV1CellsWT3, assay = "SCT", layer = "data")[GeneSet, ])
        
      }
      
      
      FeaturePlot(envv$MSlideSeqV1CellsWT3 , reduction = "umap", features = "MultiGeneExpr",
                  # raster = T, raster.dpi = c(800, 800),
                  order = T, max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        ) + ggtitle(ifelse(input$GeneExprOpr == "sum", "Summed Expression", "Mean of Expression Set"))  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    }
    
    
    
  }
}

Render_GeneExpr_SlideSeqV1MT3_Spatial <- function(envv, input, GeneSet){
  
  GeneSet = clean_GeneSet(GeneSet)
  
  
  GeneSet = GeneSet[GeneSet %in% rownames(envv$MSlideSeqV1CellsWT3)]
  
  print(GeneSet)
  
  
  
  if(is.null(envv$MSlideSeqV1CellsWT3) | length(GeneSet) < 1){
    plot(x=0, y=0, main="Load STseq MT3 dataset")
  } else {
    
    
    if(length(GeneSet)==1) {
      
    gg1=  SpatialFeaturePlot(envv$MSlideSeqV1CellsWT3 , features = GeneSet[1], 
                         # raster = T, raster.dpi = c(800, 800),
                         max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    } else{
      
      if(input$GeneExprOpr == "sum"){
        envv$MSlideSeqV1CellsWT3$MultiGeneExpr = colSums(GetAssayData(envv$MSlideSeqV1CellsWT3, assay = "SCT", layer = "data")[GeneSet, ])
        
      } else {
        envv$MSlideSeqV1CellsWT3$MultiGeneExpr = colMeans(GetAssayData(envv$MSlideSeqV1CellsWT3, assay = "SCT", layer = "data")[GeneSet, ])
        
      } 
      
      
     gg1= SpatialFeaturePlot(envv$MSlideSeqV1CellsWT3 , features = "MultiGeneExpr",
                         # raster = T, raster.dpi = c(800, 800),
                         max.cutoff = 'q99', min.cutoff = 'q01')  +
        theme_classic(base_size = 14) +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank() #,plot.title = element_blank()
        )  &
        ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))
      
    }
    
    
    # print(paste0("x: ", ranges$x*10000, " y: ", ranges$y*10000))
    
    if(is.null(ranges$GeneExpr_SlideSeqV1MT3_Spatial_x)){
      gg1
    } else {
      xrange = ggplot_build(gg1)$layout$panel_params[[1]]$x.range
      yrange = ggplot_build(gg1)$layout$panel_params[[1]]$y.range
      
      
      gg1 + coord_cartesian(xlim = c((xrange[1] + (abs(xrange[1] - xrange[2])*(min(ranges$GeneExpr_SlideSeqV1MT3_Spatial_x) - 0.2*min(ranges$GeneExpr_SlideSeqV1MT3_Spatial_x)) )), 
                                     (xrange[1] + (abs(xrange[1] - xrange[2])*(max(ranges$GeneExpr_SlideSeqV1MT3_Spatial_x) + 0.2*max(ranges$GeneExpr_SlideSeqV1MT3_Spatial_x)) ))), 
                            ylim = c((yrange[1] + (abs(yrange[1] - yrange[2])*(min(ranges$GeneExpr_SlideSeqV1MT3_Spatial_x) - 0.2*min(ranges$GeneExpr_SlideSeqV1MT3_Spatial_x)) )), 
                                     (yrange[1] + (abs(yrange[1] - yrange[2])*(max(ranges$GeneExpr_SlideSeqV1MT3_Spatial_x) + 0.2*max(ranges$GeneExpr_SlideSeqV1MT3_Spatial_x)) )) ), expand = T)
      
    }
 
  }
    
}

Threshold_GeneExpr_SlideSeqV1MT3_Spatial <- function(envv, input, GeneSet){
  
  GeneSet = clean_GeneSet(GeneSet)
  
  
  GeneSet = GeneSet[GeneSet %in% rownames(envv$MSlideSeqV1CellsWT3)]
  
  print(GeneSet)
  
  
  
  if(is.null(envv$MSlideSeqV1CellsWT3) | length(GeneSet) < 1){
    plot(x=0, y=0, main="Load STseq MT3 dataset")
  } else {
    
    
    if(length(GeneSet)==1) {
      
      envv$MSlideSeqV1CellsWT3$GeneExpr = colSums(GetAssayData(envv$MSlideSeqV1CellsWT3, assay = "SCT", layer = "data")[GeneSet, ])
      
      
    } else{
      
      if(input$GeneExprOpr == "sum"){
        envv$MSlideSeqV1CellsWT3$GeneExpr = colSums(GetAssayData(envv$MSlideSeqV1CellsWT3, assay = "SCT", layer = "data")[GeneSet, ])
        
      } else {
        envv$MSlideSeqV1CellsWT3$GeneExpr = colMeans(GetAssayData(envv$MSlideSeqV1CellsWT3, assay = "SCT", layer = "data")[GeneSet, ])
        
      }
      
      
      
    }
    
    
    tempExprVec = envv$MSlideSeqV1CellsWT3$GeneExpr
    
    envv$MSlideSeqV1CellsWT3$ExprBL = (ifelse( tempExprVec > max( tempExprVec ) * input$ExprThresh,
                                               "Above", "Below"))
    
    
    
  SpatialDimPlot(envv$MSlideSeqV1CellsWT3 ,
                 group.by = "ExprBL",
  ) +facet_wrap(~ExprBL, ncol=4) +
    theme_classic(base_size = 14) +
    guides(colour = guide_legend(override.aes = list(size = 1.5, alpha=1), ncol=1)) +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank() #,plot.title = element_blank()
    ) + ggtitle(paste0("Threshold cut above: ", round(max(tempExprVec) * input$ExprThresh,3)))
   
  }
 
  
}

Threshold_GeneExpr_STseqMT3_Spatial<- function(envv, input, GeneSet){
  
  GeneSet = clean_GeneSet(GeneSet)
  
  
  GeneSet = GeneSet[GeneSet %in% rownames(envv$MSTseqCells1)]
  
  print(GeneSet)
  
  
  
  if(is.null(envv$MSTseqCells1) | length(GeneSet) < 1){
    plot(x=0, y=0, main="Load STseq MT3 dataset")
  } else {
    
    
    if(length(GeneSet)==1) {
      
      envv$MSTseqCells1$GeneExpr = colSums(GetAssayData(envv$MSTseqCells1, assay = "SCT", layer = "data")[GeneSet, ])
      
      
    } else{
      
      if(input$GeneExprOpr == "sum"){
        envv$MSTseqCells1$GeneExpr = colSums(GetAssayData(envv$MSTseqCells1, assay = "SCT", layer = "data")[GeneSet, ])
        
      } else {
        envv$MSTseqCells1$GeneExpr = colMeans(GetAssayData(envv$MSTseqCells1, assay = "SCT", layer = "data")[GeneSet, ])
        
      }
      
      
      
    }
    
    
    tempExprVec = envv$MSTseqCells1$GeneExpr
    
    envv$MSTseqCells1$ExprBL = (ifelse( tempExprVec > max( tempExprVec ) * input$ExprThresh,
                                               "Above", "Below"))
    
    
    
    SpatialDimPlot(envv$MSTseqCells1 ,
                   group.by = "ExprBL",
    ) +facet_wrap(~ExprBL, ncol=4) +
      theme_classic(base_size = 14) +
      guides(colour = guide_legend(override.aes = list(size = 1.5, alpha=1), ncol=1)) +
      theme(legend.position = "none",
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank() #,plot.title = element_blank()
      ) + ggtitle(paste0("Threshold cut above: ", round(max(tempExprVec) * input$ExprThresh,3)))
    
  }
  
  
}


