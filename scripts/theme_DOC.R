theme_DOC <- function(){
  font="sans" 
  
  theme_minimal() %+replace%  #Use classic theme as base, replace elements i'd like to change
    
    theme(
      axis.text = element_text(
        size=14,
        family=font),
      
      axis.title = element_text(
        size=18,
        family=font),
      
      axis.line = element_line(
        size=0),
      
      panel.grid.major = element_line(
        size=0.5),
      
      panel.grid.minor.x = element_line(
        size=0.2),
      
      plot.title = element_text(
        family=font,
        size=18,
        face="bold",
        hjust=0, #Align left
        vjust=2), #Raise slightly)
      
      strip.text = element_text(
        family=font,
        size=14,
        face="italic",
        margin=margin(t=7)),
      
      #plot.background = element_rect(fill=NA,colour=NA),
  
      #strip.background = element_rect(fill=NA,colour=NA),
      
      legend.title = element_text(size = 14),  # Increase size of legend title
  
      legend.text = element_text(size = 12))    # Increase size of legend text
}

