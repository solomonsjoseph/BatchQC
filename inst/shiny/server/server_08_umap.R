### UMAP TAB ###

## Plot umap
observeEvent(input$umap_plot, {
    req(reactivevalue$se)
    validate(need(input$num_neighbors < dim(reactivevalue$se)[1] &&
            input$distance > 0, "Please select a value for the number of
        neighbors that is less than the size of your dataset and a distance
        greater than 0"))
    assay <- input$umap_assay
    msg <- sprintf('Generating plot for: %s...', paste(assay, collapse = ', '))
    withProgress(message = msg, {
        results <- umap(reactivevalue$se,
            input$umap_assay,
            input$batch_variable,
            input$num_neighbors,
            input$distance)
        setProgress(.8, 'Displaying figure...')
        output$umap <- renderPlot({
            validate(need(input$num_neighbors < dim(reactivevalue$se)[1] &&
                    input$distance > 0, "Please select a value for the number of
                    neighbors that is less than the size of your dataset and a
                    distance greater than 0"))
            results})
        setProgress(1, 'Complete.')
    })
})
