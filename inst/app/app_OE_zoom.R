


# When a double-click happens, check if there's a brush on the plot.
# If so, zoom to the brush bounds; if not, reset the zoom.
observeEvent(input$GeneExpr_SlideSeqV1MT1_Spatial_dblclick, {
  brush <- input$GeneExpr_SlideSeqV1MT1_Spatial_brush
  
  if (!is.null(brush)) {
    ranges$GeneExpr_SlideSeqV1MT1_Spatial_x <- c(brush$xmin, brush$xmax)
    ranges$GeneExpr_SlideSeqV1MT1_Spatial_y <- c(brush$ymin, brush$ymax)
    # print(paste0("x: ", ranges$x, " y: ", ranges$y))
    
  } else {
    ranges$GeneExpr_SlideSeqV1MT1_Spatial_x <- NULL
    ranges$GeneExpr_SlideSeqV1MT1_Spatial_y <- NULL
  }
})


observeEvent(input$GeneExpr_SlideSeqV1MT2_Spatial_dblclick, {
  brush <- input$GeneExpr_SlideSeqV1MT2_Spatial_brush
  
  if (!is.null(brush)) {
    ranges$GeneExpr_SlideSeqV1MT2_Spatial_x <- c(brush$xmin, brush$xmax)
    ranges$GeneExpr_SlideSeqV1MT2_Spatial_y <- c(brush$ymin, brush$ymax)
    # print(paste0("x: ", ranges$x, " y: ", ranges$y))
    
  } else {
    ranges$GeneExpr_SlideSeqV1MT2_Spatial_x <- NULL
    ranges$GeneExpr_SlideSeqV1MT2_Spatial_y <- NULL
  }
})



observeEvent(input$GeneExpr_SlideSeqV1MT3_Spatial_dblclick, {
  brush <- input$GeneExpr_SlideSeqV1MT3_Spatial_brush
  
  if (!is.null(brush)) {
    ranges$GeneExpr_SlideSeqV1MT3_Spatial_x <- c(brush$xmin, brush$xmax)
    ranges$GeneExpr_SlideSeqV1MT3_Spatial_y <- c(brush$ymin, brush$ymax)
    # print(paste0("x: ", ranges$x, " y: ", ranges$y))
    
  } else {
    ranges$GeneExpr_SlideSeqV1MT3_Spatial_x <- NULL
    ranges$GeneExpr_SlideSeqV1MT3_Spatial_y <- NULL
  }
})


# When a double-click happens, check if there's a brush on the plot.
# If so, zoom to the brush bounds; if not, reset the zoom.
observeEvent(input$GeneExpr_STseqMT1_Spatial_dblclick, {
  brush <- input$GeneExpr_STseqMT1_Spatial_brush
  
  if (!is.null(brush)) {
    ranges$GeneExpr_STseqMT1_Spatial_x <- c(brush$xmin, brush$xmax)
    ranges$GeneExpr_STseqMT1_Spatial_y <- c(brush$ymin, brush$ymax)
    print(paste0("x: ", ranges$GeneExpr_STseqMT1_Spatial_x, " y: ", ranges$GeneExpr_STseqMT1_Spatial_y))
    
  } else {
    ranges$GeneExpr_STseqMT1_Spatial_x <- NULL
    ranges$GeneExpr_STseqMT1_Spatial_y <- NULL
  }
})