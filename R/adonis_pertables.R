`adonis_pertables` <-
function(formula=X~., data, permutations = 5, method = "bray", by=NULL){

#require(vegan)

fmla<-fmla0<- formula
fmla[[2]]<-substitute(X)
Xsim<-eval(fmla0[[2]])[[2]] # selecciona el segundo elemento de un objeto pertables, que es la lista de tablas permutadas
Xraw<-eval(fmla0[[2]])[[3]] # selecciona el tercer elemento, que es la raw data table
    
# calculo de adonis sobre los raw data
X<-Xraw
adonis.raw<-adonis2(formula=fmla, data=data, permutations = permutations, method = method, by = by)
namestab<-rownames(adonis.raw)

# cálculo de adonis sobre los datos simulados
X<-Xsim
results<- lapply(X, function(X)adonis2(formula=fmla, data=data, permutations=permutations,method=method, by=by))

F<- sapply(results, function(x) x$F)
rownames(F) <- namestab

R2<- sapply(results, function(x) x$R2)
rownames(R2) <- namestab

p<- sapply(results, function(x) x$`Pr(>F)`)
rownames(p) <- namestab
R2.quant<- apply(R2, 1, quantile,c(0,0.005,0.025,0.5,0.975,0.995,1),na.rm=TRUE)[,1:(length(namestab)-2)] #quitamos las columnas de total y residuales
p.quant<- apply(p, 1, quantile,c(0,0.005,0.025,0.5,0.975,0.995,1),na.rm=TRUE)[,1:(length(namestab)-2)] #quitamos las columnas de total y residuales

F.raw<-adonis.raw$F
Fls <- data.frame(cbind(F.raw, F))
Fls <- Fls[-c((length(F.raw)-1):length(F.raw)), ]
ptax<- apply(Fls, 1, function(x) (rank(x)/(length(x)))[1])
ptax<-ifelse(ptax<=0.5,ptax,1-ptax)
adonis.raw$Prtax <- c(ptax, NA, NA)
#adonis.raw$call <- match.call()
adonis.output<- list(raw=adonis.raw, simulation=list(F=F, R2 = R2, pvalue = p, R2.quant = R2.quant, p.quant = p.quant))
#class(adonis.output) <- c("adonis.pertables",class(adonis.output))
class(adonis.output) <- c("adonis_pertables",class(adonis.output))
return(adonis.output)
}
