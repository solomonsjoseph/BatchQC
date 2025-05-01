### KBET TAB ###

## run kBET
observeEvent(input$run_kBET, {
  req(reactivevalue$se, input$kbet_assay, input$kbet_batch)
  msg <- sprintf('Generating plot for: %s...', input$kbet_assay)
  withProgress(message = msg, {
    setProgress(.5, 'Running  kBET...')
    results <- run_kBET(se = reactivevalue$se,
                        assay_to_normalize = input$kbet_assay,
                        batch = input$kbet_batch,
                        k0 = input$knet_k0)
    setProgress(.8, 'Displaying figure...')
    output$kBET_plot <- renderPlot({
      plot_kBET(results)})
    setProgress(1, 'Complete.')
  })
})
