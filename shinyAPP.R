library(shiny)
library('affy')

ui = navbarPage(title = 'Aplikacja GFP',
  
  tabPanel('Wczytaj plik',

    
      
      fileInput("file", "Wybierz pliki",multiple = TRUE, accept=c( ".CEL",'.txt')),
      
      htmlOutput('output_name')
    
  
  ),
    
    tabPanel('Obrazy',
    
      verbatimTextOutput('output', placeholder = FALSE),
      plotOutput('images',width = 1500,height = 1500)
    
    ),
  
  tabPanel('Histogramy(dane nieznormalizowane)',
           
           h2('Histogramy dla danych nieznormalizowanych'),
           plotOutput('histogram1',height = 700),
           plotOutput('histogram2',height = 700),
           plotOutput('histogram3',height = 700)

  ),
  tabPanel('Wykres degradacji RNA',
           
           plotOutput('degradation',height = 700),
           
  ),
  tabPanel('Wykresy po znormalizowaniu',
       
           plotOutput('normalizacja1',height = 700),
           
  ),
  tabPanel('Wykresy po znormalizowaniu osobne',
         
           plotOutput('normalizacja2',height = 700)
           
  ),
  )



server = function(input, output) {
  
  output$title = renderPrint({
    HTML(paste0("<center><h2><b>",'Aplikacja',"</b></h2></center>"))
  })
  
  #pokazuje jaki plik wczytany
  output$output_name = renderPrint({
    if(length(input$file$name)>0){
      a = 'Wczytano:'
      b = input$file$name
      c = paste(a,b)
      d = '\n,'

      opis = input$file$datapath[5]
      
      Dane = ReadAffy(phenoData = opis)
      
    #wczytanie
    Dane = ReadAffy(phenoData = 'opis.txt')
    if (length(b)>0){
      HTML(paste0("<h7><b>",c,d,"</b></h7>"))
    }
    }
  })
  
  output$images = renderPlot({
    if(length(input$file$name)>0){

      opis = input$file$datapath[5]
      Dane = ReadAffy(phenoData = opis)
      sample = sampleNames(Dane)
      par(mfrow=c(2,2))
      image(Dane[,1])
      image(Dane[,2])
      image(Dane[,3])
      image(Dane[,4])
      #wczytanie
      
      
    }
  })
  
  
  output$histogram1 = renderPlot({
    if(length(input$file$name)>0){
      
      opis = input$file$datapath[5]
      Dane = ReadAffy(phenoData = opis)
      sample = sampleNames(Dane)
  
      hist(Dane, col=c('red','blue','green','black'), main='Histogram w skali log',lwd = 2)
      legend('right',legend = c(sample),lty = 5,col = c('red','blue','green','black'), cex=0.9,lwd = 3)

      #wczytanie
      
      
    }
  })
  
  output$histogram2 = renderPlot({
    if(length(input$file$name)>0){
      
      opis = input$file$datapath[5]
      Dane = ReadAffy(phenoData = opis)
      sample = sampleNames(Dane)
      

      hist(Dane,log=FALSE, col=c('red','blue','green','black'), main='Histogram w skali liniowej',lwd = 2)
      legend('right', legend = c(sample),lty = 2,col = c('red','blue','green','black'),cex=0.9,lwd = 2)
      #wczytanie
      
      
    }
  })

  output$histogram3 = renderPlot({
    if (length(input$file$name)>0){
      opis = input$file$datapath[5]
      Dane = ReadAffy(phenoData = opis)
      sample = sampleNames(Dane)
      boxplot(Dane)
      
      
    }
  })
  
  output$degradation = renderPlot({
    if (length(input$file$name)>0){
      opis = input$file$datapath[5]
      Dane = ReadAffy(phenoData = opis)
      sample = sampleNames(Dane)
      rnaplot=AffyRNAdeg(Dane)
      plotAffyRNAdeg(rnaplot, col=c('red','blue','green','black'),lwd = 2)
      legend('topleft',legend = c(sample),lty = 5,col = c('red','blue','green','black'), cex=0.9,lwd = 3)
      
      
    }
  })
  
  output$normalizacja1 = renderPlot({
    if (length(input$file$name)>0){
      opis = input$file$datapath[5]
      Dane = ReadAffy(phenoData = opis)
      sample = sampleNames(Dane)
      normMAS= expresso(Dane, bgcorrect.method = 'mas',
                        normalize.method = 'constant',
                        pmcorrect.method = 'mas',
                        summary.method = 'mas')
      normRMA= expresso(Dane, bgcorrect.method = 'rma',
                        normalize.method = 'quantiles',
                        pmcorrect.method = 'pmonly',
                        summary.method = 'medianpolish')
      
      
      colors = c('red','blue','green','black')
      par(mfrow=c(1,2))
      plotDensity(exprs(normMAS), col = colors,
                  main = "histogram MAS", xlab = "intensity",lwd = 2)
      sampleMAS = sampleNames(normMAS)
      legend('right',legend = c(sampleMAS),lty = 5,col = colors)
      plotDensity(exprs(normRMA), col = rainbow(4),
                  main = "Histogram RMA", xlab = "intensity")
      sampleRMA = sampleNames(normRMA)
      legend('topright',legend = c(sampleRMA),lty = 5,col = colors)
      
      
    }else{
      h2('Laduje')
    }
  })
  
  output$normalizacja2 = renderPlot({
    if (length(input$file$name)>0){
      
      opis = input$file$datapath[5]
      samplesM=c('HCT_p53MK1.CEL','HCT_p53M0h.CEL')
      colors_M = c('black','red')
      data_M=ReadAffy(filenames = samplesM)
      MASdlaMM=expresso(data_M, bgcorrect.method = "mas", normalize.method = "constant",
                        pmcorrect.method = "mas", summary.method = "mas")
      RMAdlaM=expresso(data_M, bgcorrect.method = "rma", normalize.method = "quantiles",
                       pmcorrect.method = "pmonly", summary.method = "medianpolish")
      

      
      
      
      samplesP=c('HCT_p53PK1.CEL','HCT_p53P0h.CEL')
      colors_P = c('black','red')
      data_P=ReadAffy(filenames = samplesP)
      MASdlaP=expresso(data_P, bgcorrect.method = "mas", normalize.method = "constant",
                       pmcorrect.method = "mas", summary.method = "mas")
      RMAdlaP=expresso(data_P, bgcorrect.method = "rma", normalize.method = "quantiles",
                       pmcorrect.method = "pmonly", summary.method = "medianpolish")
      
      par(mfrow=c(2,2))
      plotDensity(exprs(MASdlaMM), col =colors_M,
                  main = "Histogram MAS HCT_p53M ", xlab="intensity")
      legend('topright', legend = samplesM, lty = 2, cex = 0.75, col = colors_M)
      plotDensity(exprs(RMAdlaM), col = colors_M,
                  main = "Histogram RMA HCT_p53M", xlab="intensity")
      legend('topright', legend = samplesM, lty = 2, cex = 0.75, col = colors_M)
      
      plotDensity(exprs(MASdlaP), col =colors_M,
                  main = "Histogram MAS HCT_p53P ", xlab="intensity")
      legend('topright', legend = samplesP, lty = 2, col = colors_M)
      plotDensity(exprs(RMAdlaP), col = colors_M,
                  main = "Histogram RMA HCT_p53P", xlab="intensity")
      legend('topright', legend = samplesP, lty = 2, col = colors_M)

      
      
    }
  })
  
  

}


shinyApp(ui = ui, server = server)
