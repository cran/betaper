`rda_pertables` <-
function (fml, data, scale = FALSE, ...) 
{
    #require(vegan)
    dataname <- deparse(substitute(data))

#rda.var <- function(X, formula, dataname, scale) {
#        rda.t <- do.call(what = rda.formula, list(formula = formula,  data = eval(parse(text = dataname)), scale = scale))
    rda.var <-function(cosa) {                 # 14/02/2020 
        rda.t <-  eval(parse(text =cosa))  # 14/02/2020    
	tot.chi <- rda.t$CCA$tot.chi/rda.t$tot.chi
        pseudoF <- with(rda.t, (CCA$tot.chi/CCA$rank)/(CA$tot.chi/CA$rank))
        pvalor <- anova(rda.t)$P[[1]]
        rda.scores <- scores(rda.t, display = c("sites", "bp"))
        return(list(result = c(tot.chi, pseudoF, pvalor), rda.scores = rda.scores))
    }
 
  

    fmla <- fmla0 <- fml
    fmla[[2]] <- substitute(X)
    #X <- eval(fmla0[[2]])[[2]]
    Xall <- eval(fmla0[[2]])[[2]] # 14/02/2020
    Xraw <- eval(fmla0[[2]])[[3]]
    
    
   # expression del RDA que se va a evaluar para cada tabla simulada # 14/02/2020
    cosa <-paste(c( "rda(", fmla, ", data=", dataname,", scale=", scale,")"), collapse="")
    
    #var.results <- lapply(X, FUN = rda.var, formula = fmla, dataname = dataname,     scale = scale)

   # loop incluido en vez del lapply # 14/02/2020
    var.results <- list()
    for(i in 1:length(Xall)) {
       # cat(i,"\n") # 14/02/2020
        X <- Xall[[i]]
        var.results[[i]] <- rda.var(cosa)
   }
	
    sites <- lapply(var.results, function(x) x$rda.scores$sites)
    env <- lapply(var.results, function(x) x$rda.scores$biplot * 
        attr(x$rda.scores, "const"))
    results <- sapply(var.results, function(x) x$result)
    row.names(results) <- c("Rsquared", "pseudoF", "p-value")
    rda.quant <- apply(results, 1, quantile, c(0, 0.005, 0.025, 
        0.5, 0.975, 0.995, 1))
    X <- Xraw
#    rda.raw <- do.call(what = rda.formula, list(formula = fmla,  data = eval(parse(text = dataname)), scale = scale))
    rda.raw <- eval(parse(text =cosa))  # 14/02/2020 
	# OJO: Ni los pseudo F ni los p(tax) se salvan, ni se exportan, ni se imprimen ni nada de nada)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   pseudoF.raw <- with(rda.raw, (CCA$tot.chi/CCA$rank)/(CA$tot.chi/CA$rank))
    ptax <- ((rank(c(pseudoF.raw, results[2, ])))/(length(results[2, 
        ]) + 1))[1]
    ptax <- ifelse(ptax <= 0.5, ptax, 1 - ptax)
    raw.anova <- anova(rda.raw)
    #raw.anova$Pr <- c(ptax, NA)
    #names(raw.anova)[6] <- "Pr(tax)"
    raw.anova$"Pr(>tax)" <- c(ptax, NA) # 14/02/2020
        
    rda.output <- list(raw = list(rda.raw = rda.raw), simulation = list(results = results, 
        rda.quant = rda.quant, sites = sites, biplot = env))
    #class(rda.output) <- c("rda.pertables", class(rda.output))
        class(rda.output) <- c("rda_pertables", class(rda.output))  # 14/02/2020
    return(rda.output)
}

