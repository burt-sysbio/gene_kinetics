theme_R <- function(){ 
  font <- "Georgia"   #assign font family up front
  
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 8),                #font size
  
      axis.title = element_text(              #axis text
        family = font,            #axis famuly
        size = 9),   
      
      panel.background = element_rect(size = 1, color = "black")
      #since the legend often requires manual tweaking 
      #based on plot content, don't define it here
    )
}