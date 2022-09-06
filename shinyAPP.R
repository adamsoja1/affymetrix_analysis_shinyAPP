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
           shinycssloaders::withSpinner(plotOutput('histogram1',height = 700)),
           shinycssloaders::withSpinner(plotOutput('histogram2',height = 700)),
           shinycssloaders::withSpinner(plotOutput('histogram3',height = 700)),

  ),
  tabPanel('Wykres degradacji RNA',
           
           plotOutput('degradation',height = 700),
           
  ),
  tabPanel('Wykresy po znormalizowaniu',
       
           shinycssloaders::withSpinner(plotOutput('normalizacja1',height = 2000)),
           
  ),
  tabPanel('Wykresy po znormalizowaniu osobne',
         
           shinycssloaders::withSpinner(plotOutput('normalizacja2',height = 700)),
           
  ),
  tabPanel('FOLD CHANGE',
           
           shinycssloaders::withSpinner(verbatimTextOutput('normalizacja3')),
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
      
      MAS_skala = exprs(normMAS)
      RMA_skala = 2^exprs(normRMA)
      
      colors = c('red','blue','green','black')
      par(mfrow=c(2,1))
      plotDensity(MAS_skala, col = colors,
                  main = "histogram MAS", xlab = "intensity",lwd = 2)
      sampleMAS = sampleNames(normMAS)
      legend('right',legend = c(sampleMAS),lty = 5,col = colors)
      plotDensity(RMA_skala, col = rainbow(4),
                  main = "Histogram RMA", xlab = "intensity")
      sampleRMA = sampleNames(normRMA)
      legend('topright',legend = c(sampleRMA),lty = 5,col = colors)
      
      
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
      
      MAS_skala_m = exprs(MASdlaMM)
      RMA_skala_m = 2^exprs(RMAdlaM)
      
      
      samplesP=c('HCT_p53PK1.CEL','HCT_p53P0h.CEL')
      colors_P = c('black','red')
      data_P=ReadAffy(filenames = samplesP)
      MASdlaP=expresso(data_P, bgcorrect.method = "mas", normalize.method = "constant",
                       pmcorrect.method = "mas", summary.method = "mas")
      RMAdlaP=expresso(data_P, bgcorrect.method = "rma", normalize.method = "quantiles",
                       pmcorrect.method = "pmonly", summary.method = "medianpolish")
      
      MAS_skala_p = exprs(MASdlaMM)
      RMA_skala_p = 2^exprs(RMAdlaM)
      
      par(mfrow=c(2,2))
      plotDensity(MAS_skala_m, col =colors_M,
                  main = "Histogram MAS HCT_p53M ", xlab="intensity")
      legend('topright', legend = samplesM, lty = 2, cex = 0.75, col = colors_M)
      
      plotDensity(RMA_skala_m, col = colors_M,
                  main = "Histogram RMA HCT_p53M", xlab="intensity")
      legend('topright', legend = samplesM, lty = 2, cex = 0.75, col = colors_M)
      
      plotDensity(MAS_skala_p, col =colors_M,
                  main = "Histogram MAS HCT_p53P ", xlab="intensity")
      legend('topright', legend = samplesP, lty = 2, col = colors_M)
      
      plotDensity(RMA_skala_p, col = colors_M,
                  main = "Histogram RMA HCT_p53P", xlab="intensity")
      
      legend('topright', legend = samplesP, lty = 2, col = colors_M)

      
      
    }
  })
  
  
  
  
  
  
  
  
  
  
  
  
  
  


  
  
output$normalizacja3 = renderPrint({
  
    

    
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
  
  MAS_skala = exprs(normMAS)
  RMA_skala = 2^exprs(normRMA)
  
  
  
  MASdlaMzad6 =MAS_skala[,4]/MAS_skala[,2]
  iloscMASMzad6=length(MASdlaMzad6[which(MASdlaMzad6 >= 1.2)])

  
  MASdlaPzad6=MAS_skala[,3]/MAS_skala[,1]
  iloscMASPzad6=length(MASdlaPzad6[which(MASdlaPzad6 >= 1.2)])

  
  RMAdlaMzad6=RMA_skala[,4]/RMA_skala[,2]
  iloscRMAMzad6=length(RMAdlaMzad6[which(RMAdlaMzad6 >= 1.2)])

  
  RMAdlaPzad6=RMA_skala[,3]/RMA_skala[,1]
  iloscRMAPzad6=length(RMAdlaPzad6[which(RMAdlaPzad6 >= 1.2)])

  
  
  samplesM=c('HCT_p53MK1.CEL','HCT_p53M0h.CEL')
  colors_M = c('black','red')
  data_M=ReadAffy(filenames = samplesM)
  MASdlaMM=expresso(data_M, bgcorrect.method = "mas", normalize.method = "constant",
                    pmcorrect.method = "mas", summary.method = "mas")
  RMAdlaM=expresso(data_M, bgcorrect.method = "rma", normalize.method = "quantiles",
                   pmcorrect.method = "pmonly", summary.method = "medianpolish")
  
  MAS_skala_m = exprs(MASdlaMM)
  RMA_skala_m = 2^exprs(RMAdlaM)
  
  
  samplesP=c('HCT_p53PK1.CEL','HCT_p53P0h.CEL')
  colors_P = c('black','red')
  data_P=ReadAffy(filenames = samplesP)
  MASdlaP=expresso(data_P, bgcorrect.method = "mas", normalize.method = "constant",
                   pmcorrect.method = "mas", summary.method = "mas")
  RMAdlaP=expresso(data_P, bgcorrect.method = "rma", normalize.method = "quantiles",
                   pmcorrect.method = "pmonly", summary.method = "medianpolish")
  
  MAS_skala_p = exprs(MASdlaMM)
  RMA_skala_p = 2^exprs(RMAdlaM)
 
  
  MASdlaMzad9=MAS_skala_m[,2]/MAS_skala_m[,1]
  iloscMASMzad9=length(MASdlaMzad9[which(MASdlaMzad9 >= 1.2)])

  print(paste('Geny ktore wzrosly przynajmniej o 20% po napromienieniu wzgledem
kontroli dla MAS5 dla normalizacji  tak jak w punkcie 9 zmienione o 20% (p53M)',iloscMASMzad9))
  
  
  MASdlaPzad9=MAS_skala_p[,2]/MAS_skala_p[,1]
  iloscMASPzad9=length(MASdlaPzad9[which(MASdlaPzad9 >= 1.2)])
  
  print(paste('Geny ktore wzrosly przynajmniej o 20% po napromienieniu wzgledem
kontroli dla MAS5 dla normalizacji  tak jak w punkcie 9 (p53P)',iloscMASPzad9))
  #RMA zad 9
  RMAdlaMzad9=RMA_skala_m[,2]/RMA_skala_m[,1]
  iloscRMAMzad9=length(RMAdlaMzad9[which(RMAdlaMzad9 >= 1.2)])

  print(paste('Geny ktore wzrosly przynajmniej o 20% po napromienieniu wzgledem
kontroli dla RMA dla normalizacji  tak jak w punkcie 9 (p53M)',iloscRMAMzad9))
  
  RMAdlaPzad9=RMA_skala_p[,2]/RMA_skala_p[,1]
  iloscRMAPzad9=length(RMAdlaPzad9[which(RMAdlaPzad9 >= 1.2)])

  print(paste('Geny ktore wzrosly przynajmniej o 20% po napromienieniu wzgledem
kontroli dla  RMA dla normalizacji  tak jak w punkcie 9 (p53P)',iloscRMAPzad9))

  
    

      print(paste(' Geny ktore wzrosly przynajmniej o 20% po napromienieniu wzgledem
kontroli dla  MAS  dla normalizacji  tak jak w punkcie 6 (p53M)',iloscMASMzad6))
      print(paste(' Geny ktore wzrosly przynajmniej o 20% po napromienieniu wzgledem
kontroli dla  MASdla normalizacji  tak jak w punkcie 6 (p53P)',iloscMASPzad6))
      print(paste(' Geny ktore wzrosly przynajmniej o 20% po napromienieniu wzgledem
kontroli dla  RMA dla normalizacji  tak jak w punkcie 6 (p53M)',iloscRMAMzad6))
      print(paste(' Geny ktore wzrosly przynajmniej o 20% po napromienieniu wzgledem
kontroli dla  RMA dla normalizacji   tak jak w punkcie 6 (p53P)',iloscRMAPzad6))
      
      
  }
      
      
      
      
      
    })
    
  
  



}


shinyApp(ui = ui, server = server)
