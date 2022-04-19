init <- function() {
  
  type <- 'plot'
  box_title <- 'Features Identified across Gradient'
  help_text <- 'The frequency of precursor identifications based on the Dinosaur search is plotted across the chromatographic gradient.'
  source_file <- 'features'
  
  .validate <- function(data, input) {
    validate(need(data()[['features']], paste0('Please provide a features.tsv file')))
  }
  
  .plotdata <- function(data, input) {
    plotdata <- data()[['features']][,c('Raw.file', 'mz', 'rtStart','charge','rtApex','rtEnd')]
    plotdata$Category = 'z = 1'
    plotdata$Category[plotdata$charge > 1] <- 'z > 1'

    return(plotdata)
  }
  
  .plot <- function(data, input) {
    .validate(data, input)
    plotdata <- .plotdata(data, input)
    validate(need((nrow(plotdata) > 1), paste0('No Rows selected')))
    
    maxRT <- max(plotdata$rtApex)
    
    ggplot(plotdata, aes(x=rtApex, color = Category)) + 
      facet_wrap(~Raw.file, nrow = 1, scales = "free_x") + 
      
      stat_bin(aes(y=..count..), size = 0.8, bins=100,position = "identity",geom="step")+
      coord_flip() + 
      labs(x='Retention Time (min)', y='Number of Precursors') +
      theme_base(input=input) +
      
      xlim(0, maxRT) +
      scale_color_manual(values=c(custom_colors[[1]], custom_colors[[6]]))+
      custom_theme
    
  }

  
  return(list(
    type=type,
    box_title=box_title,
    help_text=help_text,
    source_file=source_file,
    validate_func=.validate,
    plotdata_func=.plotdata,
    plot_func=.plot,
    dynamic_width=150,
    dynamic_width_base=150
  ))
}
