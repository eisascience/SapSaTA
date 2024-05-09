

output$GeneExpr_SlideSeqV1H1_Spatial <- renderPlot({
  
  GeneSet <- input$GeneSet
  

  
  if(input$ExprThresh %in% c(0, 1) ){
    
    Render_GeneExpr_SlideSeqV1H1_Spatial(envv, input, GeneSet)
    
  } else {
    
    # gg1
    
    Threshold_GeneExpr_SlideSeqV1H1_Spatial(envv, input, GeneSet)
      
    }
    
 
  
})


output$GeneExpr_SlideSeqV1H2_Spatial <- renderPlot({
  
  GeneSet <- input$GeneSet
  
  
  
  if(input$ExprThresh %in% c(0, 1) ){
    
    Render_GeneExpr_SlideSeqV1H2_Spatial(envv, input, GeneSet)
    
  } else {
    
    # gg1
    
    Threshold_GeneExpr_SlideSeqV1H2_Spatial(envv, input, GeneSet)
    
  }
  
  
  
})



output$GeneExpr_SlideSeqV1MT3_Spatial <- renderPlot({
  
  GeneSet <- input$GeneSet
  
  
  
  if(input$ExprThresh %in% c(0, 1) ){
    
    Render_GeneExpr_SlideSeqV1MT3_Spatial(envv, input, GeneSet)
    
  } else {
    
    # gg1
    
    Threshold_GeneExpr_SlideSeqV1MT3_Spatial(envv, input, GeneSet)
    
  }
  
  
  
})



output$GeneExpr_STseqH1_Spatial <- renderPlot({
  
  GeneSet <- input$GeneSet
  

  if(input$ExprThresh %in% c(0, 1) ){
    
    Render_GeneExpr_STseqH1_Spatial(envv, input, GeneSet)
    
  } else {
    
    
    Threshold_GeneExpr_STseqH1_Spatial(envv, input, GeneSet)
    
   
    
  }
  
  
})
