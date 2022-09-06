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

MAS_skala = exprs(normMAS)
RMA_skala = 2^exprs(normRMA)

colors = c('red','blue','green','black')
par(mfrow=c(1,2))
plotDensity(MAS_skala, col = colors,
            main = "histogram MAS", xlab = "intensity",lwd = 2)
sampleMAS = sampleNames(normMAS)
legend('right',legend = c(sampleMAS),lty = 5,col = colors)
plotDensity(RMA_skala, col = rainbow(4),
            main = "Histogram RMA", xlab = "intensity")
sampleRMA = sampleNames(normRMA)
legend('topright',legend = c(sampleRMA),lty = 5,col = colors)







 

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







MASdlaMzad6 =MAS_skala[,4]/MAS_skala[,2]
iloscMASMzad6=length(MASdlaMzad6[which(MASdlaMzad6 >= 1.2)])
samplesMASMzad6=names(MASdlaMzad6[which(MASdlaMzad6 >= 1.2)])

print(paste('Geny które wzros³y przynajmniej o 20% po napromienieniu wzglêdem
kontroli dla MAS5  dla normalizacji p53M tak jak w punkcie 6',iloscMASMzad6))

MASdlaPzad6=MAS_skala[,3]/MAS_skala[,1]
iloscMASPzad6=length(MASdlaPzad6[which(MASdlaPzad6 >= 1.2)])
samplesMASPzad6=names(MASdlaPzad6[which(MASdlaPzad6 >= 1.2)])

print(paste('Geny które wzros³y przynajmniej o 20% po napromienieniu wzglêdem
kontroli dla MAS5  dla normalizacji p53P tak jak w punkcie 6',iloscMASpzad6))

#RMA zad 6
RMAdlaMzad6=RMA_skala[,4]/RMA_skala[,2]
iloscRMAMzad6=length(RMAdlaMzad6[which(RMAdlaMzad6 >= 1.2)])
samplesRMAMzad6=names(RMAdlaMzad6[which(RMAdlaMzad6 >= 1.2)])

print(paste('Geny które wzros³y przynajmniej o 20% po napromienieniu wzglêdem
kontroli dla RMA  dla normalizacji p53M tak jak w punkcie 6',iloscRMAMzad6))

RMAdlaPzad6=RMA_skala[,3]/RMA_skala[,1]
iloscRMAPzad6=length(RMAdlaPzad6[which(RMAdlaPzad6 >= 1.2)])
samplesMAPzad6=names(RMAdlaPzad6[which(RMAdlaPzad6 >= 1.2)])

print(paste('Geny które wzros³y przynajmniej o 20% po napromienieniu wzglêdem
kontroli dla RMA  dla normalizacji p53P tak jak w punkcie 6',samplesMAPzad6))




#MAS zad 9



MAS_skala_m = exprs(MASdlaMM)
RMA_skala_m = 2^exprs(RMAdlaM)


MAS_skala_p = exprs(MASdlaMM)
RMA_skala_p = 2^exprs(RMAdlaM)

MASdlaMzad9=MAS_skala_m[,2]/MAS_skala_m[,1]
iloscMASMzad9=length(MASdlaMzad9[which(MASdlaMzad9 >= 1.2)])
samplesMASMzad9=names(MASdlaMzad9[which(MASdlaMzad9 >= 1.2)])

print(paste('Geny które wzros³y przynajmniej o 20% po napromienieniu wzglêdem
kontroli dla MAS5  dla normalizacji p53M tak jak w punkcie 9',iloscMASMzad9))


MASdlaPzad9=MAS_skala_p[,2]/MAS_skala_p[,1]
iloscMASPzad9=length(MASdlaPzad9[which(MASdlaPzad9 >= 1.2)])
samplesMASPzad9=names(MASdlaPzad9[which(MASdlaPzad9 >= 1.2)])
print(paste('Geny które wzros³y przynajmniej o 20% po napromienieniu wzglêdem
kontroli dla MAS5  dla normalizacji p53P tak jak w punkcie 9',iloscMASPzad9))
#RMA zad 9
RMAdlaMzad9=RMA_skala_m[,2]/RMA_skala_m[,1]
iloscRMAMzad9=length(RMAdlaMzad9[which(RMAdlaMzad9 >= 1.2)])
samplesRMAMzad9=names(RMAdlaMzad9[which(RMAdlaMzad9 >= 1.2)])

print(paste('Geny które wzros³y przynajmniej o 20% po napromienieniu wzglêdem
kontroli dla RMA  dla normalizacji p53M tak jak w punkcie 9',iloscRMAMzad9))

RMAdlaPzad9=RMA_skala_p[,2]/RMA_skala_p[,1]
iloscRMAPzad9=length(RMAdlaPzad9[which(RMAdlaPzad9 >= 1.2)])
samplesRMAPzad9=names(RMAdlaPzad9[which(RMAdlaPzad9 >= 1.2)])

print(paste('Geny które wzros³y przynajmniej o 20% po napromienieniu wzglêdem
kontroli dla RMA  dla normalizacji p53P tak jak w punkcie 9',iloscRMAPzad9))



sampleszad6=c(samplesMASMzad6,samplesMASPzad6,samplesRMAMzad6,samplesMAPzad6)
sampleszad9=c(samplesMASMzad9,samplesMASPzad9,samplesRMAMzad9,samplesRMAPzad9)

a1=match(sampleszad6,samplesMASMzad6,nomatch="0")
powtarz1=length(a1[which(a1!=0)])
specyf1 = length(sampleszad6)-powtarz1

a2=match(sampleszad6,samplesMASPzad6,nomatch="0")
powtarz2=length(a2[which(a2!=0)])
specyf2 = length(sampleszad6)-powtarz2

a3=match(sampleszad6,samplesRMAMzad6,nomatch="0")
powtarz3=length(a3[which(a3!=0)])
specyf3 = length(sampleszad6)-powtarz3

a4=match(sampleszad6,samplesMAPzad6,nomatch="0")
powtarz4=length(a4[which(a4!=0)])
specyf4 = length(sampleszad6)-powtarz4

a5 = match(sampleszad9,samplesMASMzad9,nomatch = "0")
powtarz5 = length(a5[which(a5!= 0)])
specyf5 = length(sampleszad9)-powtarz5

a6=match(sampleszad9,samplesMASPzad9,nomatch="0")
powtarz6=length(a6[which(a6!=0)])
specyf6 = length(sampleszad9)-powtarz6

a7=match(sampleszad9,samplesRMAMzad9,nomatch="0")
powtarz7=length(a7[which(a7!=0)])
specyf7 = length(sampleszad9)-powtarz7

a8=match(sampleszad9,samplesRMAPzad9,nomatch="0")
powtarz8=length(a8[which(a8!=0)])
specyf8 = length(sampleszad9)-powtarz8

print(paste('Liczba genów powtarzaj¹cych siê po normalizacji MAS dla zad 6 dla M: ',powtarz1))


print(paste('Liczba genów specyficznych po normalizacji MAS dla zad 6 dla M: ',specyf1))


print(paste('Liczba genów powtarzaj¹cych siê po normalizacji MAS dla zad 6 dla P: ',powtarz2))

print(paste('Liczba genów specyficznych po normalizacji MAS dla zad 6 dla P: ',specyf2))


print(paste('Liczba genów powtarzaj¹cych siê po normalizacji RMA dla zad 6 dla M: ',powtarz3))


print(paste('Liczba genów specyficznych po normalizacji RMA dla zad 6 dla M: ',specyf3))


print(paste('Liczba genów powtarzaj¹cych siê po normalizacji RMA dla zad 6 dla P: ',powtarz4))


print(paste('Liczba genów specyficznych po normalizacji RMA dla zad 6 dla P: ',specyf4))


print(paste('Liczba genów powtarzaj¹cych siê po normalizacji MAS dla zad 9 dla M: ',powtarz5))


print(paste('Liczba genów specyficznych po normalizacji MAS dla zad 9 dla M:',specyf5))


print(paste('Liczba genów powtarzaj¹cych siê po normalizacji MAS dla zad 9 dla P: ',powtarz6))


print('Liczba genów specyficznych po normalizacji MAS dla zad 9 dla P: ',specyf6))


print(paste('Liczba genów powtarzaj¹cych siê po normalizacji RMA dla zad 9 dla M: ',powtarz7))


print(paste('Liczba genów specyficznych po normalizacji RMA dla zad 9 dla M: ',specyf7))


print(paste('Liczba genów powtarzaj¹cych siê po normalizacji RMA dla zad 9 dla P:',powtarz8))


print(paste('Liczba genów specyficznych po normalizacji RMA dla zad 9 dla P:',specyf8))









