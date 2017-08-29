#Load the libraries we will need
library(shiny)

#Page single page Shiny app
fluidPage(
  
  #What to call app?
  titlePanel(""),
  
  #Create a fluid row to create an overall plot title
  fluidRow(column(9, offset = 3, h4(textOutput('title'), align = "center"))),
  
  #Create fluid row for well panel and plots
  fluidRow(
    
    #Add column for side panel
    column(3,
           
      #Create side panel
      wellPanel(
      
        #Select measurment
        selectInput("meas",label = h4("Measurement"),
                    choices = list("OGI" = 1,"OCI"  =2),selected = 1),

        #Length of interval
        sliderInput("intLength",label = h4("CI Length"),min = 90,max = 99,value = 95,step = 1,post = "%"),
      
        #Add button to download data
        fluidRow(
          column(4,actionButton("reset", "Reset")),
          column(8,downloadButton('download', 'Download'))
        ),
      
        #Add in some simple help text.
        h6(helpText('Note: You can exclude a study by clicking on it.'))
      
      )
      
    ),
    
    #Add column for spagetti plot
    column(9,
      plotOutput('metaPlot',width = "100%",height="600px",click = "metaClick")
    )
  
  )
      
)
