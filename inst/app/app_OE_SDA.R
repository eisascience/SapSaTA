observeEvent(input$Apply2MoDSTA, {
  
  envv$InfoBox_SDABrowser = "Projected SDA to the Mouse scRNA-Seq Testis Atlas (MoSTA)"
  
  selected_SDArun = input$sda.run
  # selected_compN = input$sda.comp.N
  
  if(!is.null(selected_SDArun)){
    print(selected_SDArun)
    
    envv$MoDSTA = ImputeSDA2SerV2(SerObj = envv$MoDSTA,
                                  sda_loadings = envv$SDARedDataLS$loadings[[selected_SDArun]]$loadings,
                                  # keepComps = unique(c(1, 2, CompN)),
                                  sdaObjID = selected_SDArun, plot=F, MakeReduc = F)
    print("Projection complete")
    envv$SDAcomps = colnames(envv$MoDSTA@meta.data)[grep("sda.", colnames(envv$MoDSTA@meta.data))]
    envv$SDAcomps = naturalsort::naturalsort(envv$SDAcomps)
    envv$commands$SDAproj_MoDSTA = T
    # envv$CompN = as.numeric( strsplit(selected_compN, "_")[[1]][2])
    # envv$CompNname = paste0("sda.", selected_SDArun, ".V", envv$CompN)
    
  } else{ 
    print("selected_SDArun is NULL")
    }
  
  
  
})

observeEvent(input$Apply2STseqMT1, {
  
  envv$InfoBox_SDABrowser = "Projected SDA results to the Mouse StereoSeq Spatial data"
  
  selected_SDArun = input$sda.run
  
  if(!is.null(selected_SDArun)){
    
  envv$MSTseqCells1 = ImputeSDA2SerV2(SerObj = envv$MSTseqCells1 ,
                                      sda_loadings = envv$SDARedDataLS$loadings[[selected_SDArun]]$loadings,
                                      # keepComps = unique(c(1, 2, CompN)),
                                      sdaObjID = selected_SDArun, plot=F, MakeReduc = F, assay = "SCT")
  print("Projection complete")
  envv$commands$SDAproj_MSTseqCells1 = T
  
  # envv$SDAcomps = colnames(envv$MSTseqCells1@meta.data)[grep("sda.", colnames(envv$MSTseqCells1@meta.data))]
  # envv$SDAcomps = naturalsort::naturalsort(envv$SDAcomps)
  } else{ 
    print("selected_SDArun is NULL")
  }
  
})


observeEvent(input$Apply2SlideSeqV1MT1, {
  
  envv$InfoBox_SDABrowser = "Projected SDA results to the Mouse SlideSeqV1 Spatial data"
  
  selected_SDArun = input$sda.run
  
  if(!is.null(selected_SDArun)){
    
    envv$MSlideSeqV1CellsWT1 = ImputeSDA2SerV2(SerObj = envv$MSlideSeqV1CellsWT1 ,
                                        sda_loadings = envv$SDARedDataLS$loadings[[selected_SDArun]]$loadings,
                                        # keepComps = unique(c(1, 2, CompN)),
                                        sdaObjID = selected_SDArun, plot=F, MakeReduc = F, assay = "SCT")
    print("Projection complete")
    envv$commands$SDAproj_MSlideSeqV1CellsWT1 = T
    
    # envv$SDAcomps = colnames(envv$MSlideSeqV1CellsWT1@meta.data)[grep("sda.", colnames(envv$MSlideSeqV1CellsWT1@meta.data))]
    # envv$SDAcomps = naturalsort::naturalsort(envv$SDAcomps)
  } else{ 
    print("selected_SDArun is NULL")
  }
  
})

observeEvent(input$NextComp, {
  selected_compN = input$sda.comp.N
  CompN = as.numeric( strsplit(selected_compN, "_")[[1]][2])
 
  n_comps <- nrow(envv$SDARedDataLS$loadings[[selected_SDArun()]]$loadings)
  
  
  if (CompN < n_comps) {
    new_comp <- CompN + 1
    updateSelectInput(session, "sda.comp.N", selected = paste0("Comp_", new_comp))
  }
  
})


observeEvent(input$PreviousComp, {
  selected_compN = input$sda.comp.N
  CompN = as.numeric( strsplit(selected_compN, "_")[[1]][2])
  
  n_comps <- nrow(envv$SDARedDataLS$loadings[[selected_SDArun()]]$loadings)
  
  
  if (CompN > 1) {
    new_comp <- CompN - 1
    updateSelectInput(session, "sda.comp.N", selected = paste0("Comp_", new_comp))
  }
  
})

# observeEvent(input$Gene2SearchSDA,{
#   
#   req(input$Gene2SearchSDA)
#   req(envv$SDARedDataLS)
# 
#    
#   
# 
#   PosHits = lapply(envv$Top_loaded_genes$Top_Pos, function(x){
#     input$Gene2SearchSDA %in% x
#   }) %>% unlist() %>% which() %>% names() %>% naturalsort::naturalsort()
#   
#   PosHits = gsub("sda.", "", PosHits)
#   
#   NegHits = lapply(envv$Top_loaded_genes$Top_Neg, function(x){
#     input$Gene2SearchSDA %in% x
#   }) %>% unlist() %>% which() %>% names() %>% naturalsort::naturalsort()
#   
#   NegHits = gsub("sda.", "", NegHits)
#   
#   
#   if(length(PosHits)==0) PosHits = "Gene not found top pos loaded"
#   if(length(NegHits)==0) NegHits = "Gene not found top neg loaded"
#   
#   # Update UI with results
#   output$displayPosHits <- renderText({ PosHits })
#   output$displayNegHits <- renderText({ NegHits })
#   
# })


PosHits <- reactive({
  req(input$Gene2SearchSDA)
  req(envv$SDARedDataLS)
  hits = lapply(envv$Top_loaded_genes$Top_Pos, function(x) {
    input$Gene2SearchSDA %in% x
  }) %>% unlist() %>% which() %>% names() %>% naturalsort::naturalsort()
  hits = gsub("sda.", "", hits)
  if(length(hits) == 0) "Gene not found in top positive loaded" else hits
})

NegHits <- reactive({
  req(input$Gene2SearchSDA)
  req(envv$SDARedDataLS)
  hits = lapply(envv$Top_loaded_genes$Top_Neg, function(x) {
    input$Gene2SearchSDA %in% x
  }) %>% unlist() %>% which() %>% names() %>% naturalsort::naturalsort()
  hits = gsub("sda.", "", hits)
  if(length(hits) == 0) "Gene not found in top negative loaded" else hits
})

# Display the results
output$displayPosHits <- renderText({ PosHits() })
output$displayNegHits <- renderText({ NegHits() })

# Ensure the reactive is evaluated at startup
observe({
  input$Gene2SearchSDA
  PosHits()
  NegHits()
})