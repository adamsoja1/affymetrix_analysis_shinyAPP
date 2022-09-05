if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("affy")
library('affy')


opis=read.table('opis.txt')
opis
#wczytanie
Dane = ReadAffy(phenoData = 'opis.txt')
sample = sampleNames(Dane)
#obrazy
image(Dane[,1])
image(Dane[,2])
image(Dane[,3])
image(Dane[,4])



#histogramy

hist(Dane, col=c('red','blue','green','black'), main='Histogram w skali log',lwd = 2)
legend('right',legend = c(sample),lty = 5,col = c('red','blue','green','black'), cex=0.9,lwd = 3)
hist(Dane,log=FALSE, col=c('red','blue','green','black'), main='Histogram w skali liniowej',lwd = 2)
legend('right', legend = c(sample),lty = 2,col = c('red','blue','green','black'),cex=0.9,lwd = 2)




boxplot(Dane, main="Wykresy pude³kowe")



rnaplot=AffyRNAdeg(Dane)
plotAffyRNAdeg(rnaplot, col=c('red','blue','green','black'),lwd = 2)
legend('topleft',legend = c(sample),lty = 5,col = c('red','blue','green','black'), cex=0.9,lwd = 3)


normMAS= expresso(Dane, bgcorrect.method = 'mas',
              normalize.method = 'constant',
              pmcorrect.method = 'mas',
              summary.method = 'mas')
normRMA= expresso(Dane, bgcorrect.method = 'rma',
              normalize.method = 'quantiles',
              pmcorrect.method = 'pmonly',
              summary.method = 'medianpolish')


colors = c('red','blue','green','black')
plotDensity(exprs(MAS), col = colors,
            main = "histogram MAS", xlab = "intensity",lwd = 2)
sampleMAS = sampleNames(MAS)
legend('right',legend = c(sampleMAS),lty = 5,col = colors)
plotDensity(exprs(RMA), col = rainbow(4),
            main = "Histogram RMA", xlab = "intensity")
sampleRMA = sampleNames(RMA)
legend('topright',legend = c(sampleRMA),lty = 5,col = colors)








 

samplesM=c('HCT_p53MK1.CEL','HCT_p53M0h.CEL')
colors_M = c('black','red')
data_M=ReadAffy(filenames = samplesM)
MASdlaMM=expresso(data_M, bgcorrect.method = "mas", normalize.method = "constant",
               pmcorrect.method = "mas", summary.method = "mas")
RMAdlaM=expresso(data_M, bgcorrect.method = "rma", normalize.method = "quantiles",
               pmcorrect.method = "pmonly", summary.method = "medianpolish")

plotDensity(exprs(MASdlaMM), col =colors_M,
            main = "Histogram MAS HCT_p53M ", xlab="intensity")
legend('topright', legend = samplesM, lty = 2, cex = 0.75, col = colors_M)
plotDensity(exprs(RMAdlaM), col = colors_M,
            main = "Histogram RMA HCT_p53M", xlab="intensity")
legend('topright', legend = samplesM, lty = 2, cex = 0.75, col = colors_M)







samplesP=c('HCT_p53PK1.CEL','HCT_p53P0h.CEL')
colors_P = c('black','red')
data_P=ReadAffy(filenames = samplesP)
MASdlaP=expresso(data_P, bgcorrect.method = "mas", normalize.method = "constant",
                  pmcorrect.method = "mas", summary.method = "mas")
RMAdlaP=expresso(data_P, bgcorrect.method = "rma", normalize.method = "quantiles",
                 pmcorrect.method = "pmonly", summary.method = "medianpolish")

plotDensity(exprs(MASdlaP), col =colors_M,
            main = "Histogram MAS HCT_p53P ", xlab="intensity")
legend('topright', legend = samplesP, lty = 2, col = colors_M)
plotDensity(exprs(RMAdlaP), col = colors_M,
            main = "Histogram RMA HCT_p53P", xlab="intensity")
legend('topright', legend = samplesP, lty = 2, col = colors_M)







MASdlaMzad6 =(exprs(normMAS[,4]))/(exprs(normMAS[,2]))

iloscMASMzad6=length(MASdlaMzad6[which(MASdlaMzad6 >= 1.2),1])
samplesMASMzad6=names(MASdlaMzad6[which(MASdlaMzad6 >= 1.2),1])
MASdlaPzad6=(exprs(normMAS[,3]))/(exprs(normMAS[,1]))
iloscMASPzad6=length(MASdlaPzad6[which(MASdlaPzad6 >= 1.2),1])
samplesMASPzad6=names(MASdlaPzad6[which(MASdlaPzad6 >= 1.2),1])

#RMA zad 6
RMAdlaMzad6=(exprs(normRMA[,4]))/(exprs(normRMA[,2]))
iloscRMAMzad6=length(RMAdlaMzad6[which(RMAdlaMzad6 >= 1.2),1])
samplesRMAMzad6=names(RMAdlaMzad6[which(RMAdlaMzad6 >= 1.2),1])
RMAdlaPzad6=(exprs(normRMA[,3]))/(exprs(normRMA[,1]))
iloscRMAPzad6=length(RMAdlaPzad6[which(RMAdlaPzad6 >= 1.2),1])
samplesMAPzad6=names(RMAdlaPzad6[which(RMAdlaPzad6 >= 1.2),1])
#MAS zad 9

MASdlaMzad9=(exprs(MASdlaMM[,2]))/(exprs(MASdlaMM[,1]))
iloscMASMzad9=length(MASdlaMzad9[which(MASdlaMzad9 >= 1.2),1])
samplesMASMzad9=names(MASdlaMzad9[which(MASdlaMzad9 >= 1.2),1])

MASdlaPzad9=(exprs(MASdlaP[,2]))/(exprs(MASdlaP[,1]))
iloscMASPzad9=length(MASdlaPzad9[which(MASdlaPzad9 >= 1.2),1])
samplesMASPzad9=names(MASdlaPzad9[which(MASdlaPzad9 >= 1.2),1])
#RMA zad 9
RMAdlaMzad9=(exprs(RMAdlaM[,2]))/(exprs(RMAdlaM[,1]))
iloscRMAMzad9=length(RMAdlaMzad9[which(RMAdlaMzad9 >= 1.2),1])
samplesRMAMzad9=names(RMAdlaMzad9[which(RMAdlaMzad9 >= 1.2),1])

RMAdlaPzad9=(exprs(RMAdlaP[,2]))/(exprs(RMAdlaP[,1]))
iloscRMAPzad9=length(RMAdlaPzad9[which(RMAdlaPzad9 >= 1.2),1])
samplesRMAPzad9=names(RMAdlaPzad9[which(RMAdlaPzad9 >= 1.2),1])




sampleszad6=c(samplesMASMzad6,samplesMASPzad6,samplesRMAMzad6,samplesMAPzad6)
sampleszad9=c(samplesMASMzad9,samplesMASPzad9,samplesRMAMzad9,samplesRMAPzad9)

a1=match(zad6,samplesMASMzad6,nomatch="0")
powtarz1=length(a1[which(a1!=0)])
specyf1 = length(sampleszad6)-powtarz1

a2=match(zad6,samplesMASPzad6,nomatch="0")
powtarz2=length(a2[which(a2!=0)])
specyf2 = length(sampleszad6)-powtarz2

a3=match(zad6,samplesRMAMzad6,nomatch="0")
powtarz3=length(a3[which(a3!=0)])
specyf3 = length(sampleszad6)-powtarz3

a4=match(zad6,samplesMAPzad6,nomatch="0")
powtarz4=length(a4[which(a4!=0)])
specyf4 = length(sampleszad6)-powtarz4

a5 = match(zad9,samplesMASMzad9,nomatch = "0")
powtarz5 = length(a5[which(a5!= 0)])
specyf5 = length(sampleszad9)-powtarz5

a6=match(zad9,samplesMASPzad9,nomatch="0")
powtarz6=length(a6[which(a6!=0)])
specyf6 = length(sampleszad9)-powtarz6

a7=match(zad9,samplesRMAMzad9,nomatch="0")
powtarz7=length(a7[which(a7!=0)])
specyf7 = length(sampleszad9)-powtarz7

a8=match(zad9,samplesRMAPzad9,nomatch="0")
powtarz8=length(a8[which(a8!=0)])
specyf8 = length(sampleszad9)-powtarz8

print('Liczba genów powtarzaj¹cych siê po normalizacji MAS dla zad 6 dla M:')
powtarz1

print('Liczba genów specyficznych po normalizacji MAS dla zad 6 dla M:')
specyf1

print('Liczba genów powtarzaj¹cych siê po normalizacji MAS dla zad 6 dla P:')
powtarz2

print('Liczba genów specyficznych po normalizacji MAS dla zad 6 dla P:')
specyf2

print('Liczba genów powtarzaj¹cych siê po normalizacji RMA dla zad 6 dla M:')
powtarz3

print('Liczba genów specyficznych po normalizacji RMA dla zad 6 dla M:')
specyf3

print('Liczba genów powtarzaj¹cych siê po normalizacji RMA dla zad 6 dla P:')
powtarz4

print('Liczba genów specyficznych po normalizacji RMA dla zad 6 dla P:')
specyf4

print('Liczba genów powtarzaj¹cych siê po normalizacji MAS dla zad 9 dla M:')
powtarz5

print('Liczba genów specyficznych po normalizacji MAS dla zad 9 dla M:')
specyf5

print('Liczba genów powtarzaj¹cych siê po normalizacji MAS dla zad 9 dla P:')
powtarz6

print('Liczba genów specyficznych po normalizacji MAS dla zad 9 dla P:')
specyf6

print('Liczba genów powtarzaj¹cych siê po normalizacji RMA dla zad 9 dla M:')
powtarz7

print('Liczba genów specyficznych po normalizacji RMA dla zad 9 dla M:')
specyf7

print('Liczba genów powtarzaj¹cych siê po normalizacji RMA dla zad 9 dla P:')
powtarz8

print('Liczba genów specyficznych po normalizacji RMA dla zad 9 dla P:')
specyf8








