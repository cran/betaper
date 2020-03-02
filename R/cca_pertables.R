cca_pertables <-
function (fml, data,  scale = FALSE,...) 
{
    #require(vegan) # 14/02/2020
    dataname <- deparse(substitute(data))
    
    
    #cca.var <- function(X, formula, data, ...) { # 14/02/2020 
    #     cca.t <- do.call(what = cca, list(formula = fml,   data = eval(parse(text = dataname)), scale = scale))
   cca.var <- function(cosa) {                 # 14/02/2020 
        cca.t <-  eval(parse(text =cosa))  # 14/02/2020    
        tot.chi <- cca.t$CCA$tot.chi/cca.t$tot.chi
        pseudoF <- with(cca.t, (CCA$tot.chi/CCA$rank)/(CA$tot.chi/CA$rank))
        pvalor <- anova(cca.t)$P[[1]]
        cca.scores <- scores(cca.t, display = c("sites", "bp"))
        return(list(result = c(tot.chi, pseudoF, pvalor), cca.scores = cca.scores))
    }
    
     

    fmla <- fmla0 <- fml
    fmla[[2]] <- substitute(X)
    #X <- eval(fmla0[[2]])[[2]]
    Xall <- eval(fmla0[[2]])[[2]] # 14/02/2020
    Xraw <- eval(fmla0[[2]])[[3]]
   
   # expression del CCA que se va a evaluar para cada tabla simulada # 14/02/2020
    cosa <-paste(c( "cca(", fmla, ", data=", dataname,", scale=", scale,")"), collapse="")
    
    #var.results <- lapply(X, FUN = cca.var, formula = fmla, dataname = dataname,     scale = scale)
    # loop incluido en vez del lapply # 14/02/2020
    var.results <- list()
    for(i in 1:length(Xall)) {
       # cat(i,"\n") # 14/02/2020
        X <- Xall[[i]]
        var.results[[i]] <- cca.var(cosa)
   }
	
    sites <- lapply(var.results, function(x) x$cca.scores$sites)
    env <- lapply(var.results, function(x) x$cca.scores$biplot)
    results <- sapply(var.results, function(x) x$result)
    row.names(results) <- c("Rsquared", "pseudoF", "p-value")
    cca.quant <- apply(results, 1, quantile, c(0, 0.005, 0.025, 
        0.5, 0.975, 0.995, 1))
    X <- Xraw
    
    #cca.raw <- do.call(what = cca, list(formula = fmla,   data = eval(parse(text = dataname)), scale = scale))
    cca.raw <- eval(parse(text =cosa))  # 14/02/2020 
    pseudoF.raw <- with(cca.raw, (CCA$tot.chi/CCA$rank)/(CA$tot.chi/CA$rank))
    ptax <- ((rank(c(pseudoF.raw, results[2, ])))/(length(results[2, 
        ]) + 1))[1]
    ptax <- ifelse(ptax <= 0.5, ptax, 1 - ptax)
    raw.anova <- anova(cca.raw)
    # raw.anova$Pr <- c(ptax, NA)
     #names(raw.anova)[6] <- "Pr(tax)" 
     raw.anova$"Pr(tax)" <- c(ptax, NA) # 14/02/2020
     
     # OJO: Ni los pseudoF ni los Pr(tax) se devuelven en ningun sitio !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
    cca.output <- list(raw = list(cca.raw = cca.raw), simulation = list(results = results, 
        cca.quant = cca.quant, sites = sites, biplot = env))
    #class(cca.output) <- c("cca.pertables", class(cca.output))
    class(cca.output) <- c("cca_pertables", class(cca.output)) # 14/02/2020
    return(cca.output)
}

