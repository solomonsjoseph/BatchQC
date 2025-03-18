### UMAP TAB ###

## Plot umap
observeEvent(input$umap_plot, {
    req(reactivevalue$se)
    validate(need(input$num_neighbors < dim(reactivevalue$se)[2], "Please
            select a value for the number of neighbors that is less than the
            size of your dataset"),
        need(input$distance > 0 && input$distance < input$spread, "Please
            select a distance that is greater than 0 and less than spread"),
        need(!input$exploratory ||
                (input$exploratory && (dim(reactivevalue$se)[2] > 15)),
        "Exploratory plot is only possible for datasets with more than 15
        samples."))
    assay <- input$umap_assay
    msg <- sprintf('Generating plot for: %s...', paste(assay, collapse = ', '))
    withProgress(message = msg, {
        results <- umap(reactivevalue$se,
            input$umap_assay,
            input$batch_variable,
            input$num_neighbors,
            input$distance,
            input$spread,
            input$exploratory)
        setProgress(.8, 'Displaying figure...')
        output$umap <- renderPlot({
            results})
        setProgress(1, 'Complete.')
    })
})
