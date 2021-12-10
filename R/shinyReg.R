#' Visualize two-predictor regression using shiny
#'
#' Opens shiny window to allow basic interactive exploration of a dataset.
#'
#' @param dat data frame that you would like to explore
#' @param overallTitle What to you want your window's big title to be?
#' 
#' @export
#' @examples
#' shinyReg(dat_cochraneEtAl_2019_PLOSOne)
#' 

shinyReg <- function(dat,overallTitle = 'Data Exploration'){
  
  # # TO DO: 
  # include more instructions
  # include warnings about what's problematic
  # include filtering if possible (e.g., looking at subgroups)
  # more description and example

  options(warn = -1)

require(shiny)  ;  require(lmSupport) ; require(effects) ; library(psych)

ui <- fluidPage(

  # Application title
  titlePanel(overallTitle)

  # , sliderInput(inputId = 'input1name',label = 'input1label',min=1,max=100,value=10)
  , selectInput(inputId = 'x1Ax',label ='x1Ax : ', choices=colnames(dat))
  , selectInput(inputId = 'x2Ax',label ='x2Ax : ', choices=colnames(dat))
  , selectInput(inputId = 'yAx',label ='yAx : ', choices=colnames(dat))
  ,selectInput(inputId='interactionOrNot',label='Interaction? * = "yes"; + = "no"',choices=c('*','+'))
  # numericInput(inputId = 'input2name',label = 'input2label')


  , plotOutput(outputId = 'scatter')

  , verbatimTextOutput(outputId='modExpl')
  , verbatimTextOutput(outputId='modTxt')
  , verbatimTextOutput(outputId='effTxt')
  , verbatimTextOutput(outputId='pairPanelTxt')
  , plotOutput(outputId = 'pairPanel')
  # , plotOutput(outputId = 'histX2')
  # , plotOutput(outputId = 'histY')
  # tableOutput(outputId = 'output2name')


)


server <- function(input, output) {

  curDat <- reactive({
    cbind(dat[input$x1Ax],dat[input$x2Ax],dat[input$yAx])
  })

  regForm <- reactive({
    as.formula(paste(as.character(input$yAx),'~',as.character(input$x1Ax),input$interactionOrNot,as.character(input$x2Ax)))
  })

  regMod <- reactive({
    lm(regForm(),data=curDat())
  })

  output$modExpl <- renderPrint(cat(paste('Regression using ',input$x1Ax,' and ',input$x2Ax,' to predict ',input$yAx,': ',sep='')))
  output$modTxt <- renderPrint({modelSummary(regMod(),t=F)})
  output$effTxt <- renderPrint({modelEffectSizes(regMod())})

  output$scatter <- renderPlot({
    plot(allEffects(regMod(),partial.residuals=T))
    # plot(curDat()[,1],curDat()[,2],main='Scatterplot',xlab=input$x1Ax,ylab=input$yAx)
    # abline(regMod()$coefficients[1],regMod()$coefficients[2])
    # curPlot <- ggplot(dat) + geom_point(input$x1Ax,input$yAx) ; print(curPlot)
  })

  output$pairPanel <- renderPlot({
    pairs.panels(curDat())# ,freq=floor(length(curDat()[,1])/2))
  })

  output$pairPanelTxt <- renderPrint(cat('Histograms and correlations of the selected variables:'))

  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"
  )
}

# Run the application
capture.output({suppressWarnings({
shinyApp(ui = ui, server = server)
})}, file='NUL')
# options(warn = 1)
}
