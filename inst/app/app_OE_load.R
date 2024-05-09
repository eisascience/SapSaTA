# observeEvent(input$MoDSTA, {
#   
#   envv$InfoBox_Load = "Loaded Mouse scRNA-Seq Testis Atlas (MoSTA)"
#   
#   envv = Load_MoDSTA(envv)
#   
# 
# })

observeEvent(input$SlideSeqV1Sapien1, {
  
  envv$InfoBox_Load = "Loaded SlideSeq V1 Spatial omics homo sapien testis 1"
  
  envv = Load_SlideSeqH1(envv)
  
})

observeEvent(input$SlideSeqV1Sapien2, {
  
  envv$InfoBox_Load = "Loaded SlideSeq V1 Spatial omics homo sapien testis 2"
  
  envv = Load_SlideSeqH2(envv)
  
})
# 
# observeEvent(input$SlideSeqV1Mouse3, {
#   
#   envv$InfoBox_Load = "Loaded SlideSeq V1 Spatial omics mouse testis 3"
#   
#   envv = Load_SlideSeqMT3(envv)
#   
# })
# 
# 
# observeEvent(input$STseqMouse1, {
#   
#   envv$InfoBox_Load = "Loaded Stereo-seq Spatial omics mouse testis 1"
#   
#   envv = Load_STseqH1(envv)
#   
# })
# 
# observeEvent(input$SDAres, {
#   
#   envv$InfoBox_Load = "Loaded SDA results trained on Mouse scRNA-Seq Testis Atlas (MoSTA)"
#   
#   envv = Load_SDA(envv)
#   
# })



# observeEvent(input$BiomaRtLoad, {
#   
#   # loadStatus("Loading data...")
#   # Disable the button while loading to prevent repeated clicks
#   # shinyjs::disable("loadData")
#   
#   # BioMart related operations
#   # library(biomaRt)
#   # envv$mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#   # envv$genes <- getBM(attributes = c("chromosome_name", "start_position", "end_position"), mart = mart)
#   envv = Load_genes(envv)
#   # envv$genes = readRDS("./data/genes.rds")
#   # Assuming 'genes' and 'mart' are to be used in global environment
#   # assign("genes", genes, envir = .GlobalEnv)
#   # assign("mart", mart, envir = .GlobalEnv)
#   
#   # Re-enable button and update status after load
#   # shinyjs::enable("loadData")
#   # loadStatus("Data loaded successfully.")
# })

observeEvent(input$AllInputs, {
  
  envv$InfoBox_Load = "Loaded all inputs: Spatial H1 & H2"
  
  
  # envv = Load_MoDSTA(envv)
  # envv = Load_STseqH1(envv)
  envv = Load_SlideSeqH1(envv)
  envv = Load_SlideSeqH2(envv)
  # envv = Load_SlideSeqMT3(envv)
  # envv = Load_SDA(envv)
  # envv = Load_genes(envv)
})