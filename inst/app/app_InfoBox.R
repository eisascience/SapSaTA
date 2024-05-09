
output$InfoBox_Load <- renderValueBox({
  valueBox(
    value = "Info Bar", #format(Sys.time(), "%a %b %d %X %Y %Z"),
    subtitle = envv$InfoBox_Load,
    icon = icon("area-chart"),
    color = "yellow" #if (downloadRate >= input$rateThreshold) "yellow" else "aqua"
  )
})

output$InfoBox_SDABrowser <- renderValueBox({
  valueBox(
    value = "Info Bar", #format(Sys.time(), "%a %b %d %X %Y %Z"),
    subtitle = envv$InfoBox_SDABrowser,
    icon = icon("area-chart"),
    color = "yellow" #if (downloadRate >= input$rateThreshold) "yellow" else "aqua"
  )
})
