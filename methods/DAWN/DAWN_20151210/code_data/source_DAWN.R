library(glmnet)

## Preprocess data and call functions to perform complete DAWN analysis (PNS and HMRF)
identify_risk_genes <- function(GLA_pvalue, GLA_genenames, expdata, covGene, fixedGene, pthres_pns, corthres, lambda, 
    pthres_hmrf, iter, trimthres, file_out_pns = "pns.network", file_out_hmrf = "DAWN_result.csv") {
    addCovGene <- hasArg(covGene)
    addFixedGene <- hasArg(fixedGene)
    cat("Include additional covariate? ", addCovGene, "\nFix gene hidden states? ", addFixedGene, "\n")
    t_start <- proc.time()
    
    cat("Preprocessing data... ", "\n")
    missing_data <- which(is.na(GLA_genenames))
    if (length(missing_data)>0) {
      GLA_pvalue <- GLA_pvalue[-missing_data]
      GLA_genenames <- GLA_genenames[-missing_data]
    }   
    dup_in_GLA <- duplicated(GLA_genenames)
    if (sum(dup_in_GLA) > 0) {
        ## remove duplicated genesz
        GLA_pvalue <- GLA_pvalue[-which(dup_in_GLA)]
        GLA_genenames <- GLA_genenames[-which(dup_in_GLA)]
    }
    common_genes <- intersect(row.names(expdata), GLA_genenames)  ##keep genes have both expression data and GLA p-values
    if (length(common_genes) == 0) {
        cat("\n\nERROR: no common genes shared by GLA data and expression data, please check gene names again...\n")
        return()
    }
    expdata <- expdata[common_genes, ]
    GLA_pvalue <- as.data.frame(GLA_pvalue)
    rownames(GLA_pvalue) <- GLA_genenames
    GLA_pvalue <- GLA_pvalue[common_genes, , drop = FALSE]
    GLA_pvalue <- GLA_pvalue
    cat("Preprocessing finished. ", "\n\n")
    
    ####################################################################################################################################
    ## Estimate the network
    ## pthres is the threshold for p-value. When pthres gets larger, the number of nodes in the network will increase. 
    ## corthres is the threshold for pairwise-correlation. When corthres is smaller, more neighbors of risk genes are kept. 
    ## lambda is the tuning parameter for lasso. When lambda gets larger, the estiamted network will have fewer edges. 
    ####################################################################################################################################
    cat("Start PNS... ", "\n")
    t_pns_start <- proc.time()
    graphres_PNS <- PNS_algorithm(expdata = expdata, GLA_pvalue = GLA_pvalue, pthres_pns = pthres_pns, 
        corthres = corthres, lambda = lambda)
    t_pns_end <- proc.time()
    if (!is.list(graphres_PNS)) {
        cat("\tnumber of genes after PNS: 0\n")
        cat("\nERROR: cannot perform HMRF on empty PNS graph, please adjust your parameters.\n")
        return(0)
    }
    cat("PNS finished; writing to ", file_out_pns, "...\n")
    graphfinal_PNS <- graphres_PNS$graphfinal  ##estimated network
    genelist_PNS <-  graphres_PNS$Gene ##a list of genes in the estimated network
    rownames(graphfinal_PNS) <- genelist_PNS
    colnames(graphfinal_PNS) <- genelist_PNS
    write.csv(graphfinal_PNS, file_out_pns)
    cat("Writing to ", file_out_pns, " finished. ", "\n\n")
    ####################################################################################################################################
    ##DAWN main algorithm
    ####################################################################################################################################
    cat("Start HMRF analysis... ", "\n")
    t_hmrf_start <- proc.time()
    if (addCovGene) {
        if (addFixedGene) {
          result_DAWN <- DAWN_HMRF(pv = graphres_PNS$pvfinal, graph = graphres_PNS$graphfinal, covGene = covGene, fixedGene = fixedGene, 
                                   pthres_hmrf = pthres_hmrf, iter = iter, trimthres = trimthres)
        } else {
          result_DAWN <- DAWN_HMRF(pv = graphres_PNS$pvfinal, graph = graphres_PNS$graphfinal, covGene = covGene, pthres_hmrf = pthres_hmrf, 
                                   iter = iter, trimthres = trimthres)
        }
        cat("\nEstimated parameters:\n")
        cat("b0 =", result_DAWN$b0, ", b1 =", result_DAWN$b1, ", b2 =", result_DAWN$b2, ", mu1 =", result_DAWN$mu1, ", sigmas =", result_DAWN$sigmas, "\n")
    } else {
        if (addFixedGene) {
          result_DAWN <- DAWN_HMRF(pv = graphres_PNS$pvfinal, graph = graphres_PNS$graphfinal, fixedGene = fixedGene, 
                                   pthres_hmrf = pthres_hmrf, iter = iter, trimthres = trimthres)
        } else {
          result_DAWN <- DAWN_HMRF(pv = graphres_PNS$pvfinal, graph = graphres_PNS$graphfinal, pthres_hmrf = pthres_hmrf, 
                                   iter = iter, trimthres = trimthres)
        }
        cat("\nEstimated parameters:\n")
        cat("b0 = ", result_DAWN$b0, ", b1 = ", result_DAWN$b1, ", mu1 =", result_DAWN$mu1, ", sigmas =", result_DAWN$sigmas, "\n")
    }
    t_hmrf_end <- proc.time()
    cat("HMRF analysis finished; writing to ", file_out_hmrf, "...\n")
    report <- result_report(finalposter = (1 - result_DAWN$post), genes = graphres_PNS$Gene, pv = graphres_PNS$pvfinal)  ## generate report
    write.csv(report, file_out_hmrf)
    cat("Writing to ", file_out_hmrf, " finished. ", "\n\n")
    t_end <- proc.time()
    
    cat("Total time:\t", (t_end - t_start)[3], "s\n")
    cat("\tPNS:\t", (t_pns_end - t_pns_start)[3], "s\n")
    cat("\tHMRF:\t", (t_hmrf_end - t_hmrf_start)[3], "s\n")
}

## Main function for PNS algorithm
## Input:
##      expdata - dataframe, expression data
##      GLA_pvalue - dataframe, GLA p-values
##      pthres_pns - threshold for p-value screening
##      corthres - threshold for pairwise-correlation
##      lambda - tuning parameter for lasso
## Output:
##      graphres - a list object representing the estimated PNS network
##          graphres$Gene - genes in the network
##          graphres$graphfinal - estimated unweighted, undirected network in matrix format
##          graphres$pvfinal - p-values of genes in the network
PNS_algorithm <- function(expdata,GLA_pvalue, pthres_pns, corthres, lambda) {  
    genenames <- rownames(GLA_pvalue)
    pv <- GLA_pvalue[,1]
    exp_gene <- as.character(row.names(expdata))
    pvkeep <- rep(0, length(exp_gene))
    
    ## permutate rows in GLA data to match the gene ordering in expression data
    pvi <- matrix(c(1:length(genenames)), length(genenames), 1)
    rownames(pvi) <- genenames
    pvkeep <- pv[pvi[exp_gene, ]]
    cat("\tnumber of genes before PNS: ", length(pvkeep), "\n")
    
    ## screen nodes (step 1, 2, and 3 in Algorithm 1)
    cleaned_graph <- cleanup_f(pv = pvkeep, expdata = expdata, pthres_pns = pthres_pns, corthres = corthres)
    keepi <- c(cleaned_graph$kscreeni, cleaned_graph$allnreal)
    if (length(cleaned_graph$yexp) > 0) {
        ## construct graph by applying Meinshausen and Bühlmann (2006)’s regression based approach
        graphres <- meishausen_f(pv = pvkeep, cleaned_graph = cleaned_graph, lambda = lambda)
        graphres$Gene <- exp_gene[keepi]
        cat("\tnumber of genes after PNS: ", length(keepi), "\n")
        return(graphres)
    } else {
        return(0)
    }
}

## Main function for HMRF analysis: perform HMRF analysis to find risk genes from the network estimated by PNS algorithm 
## Input:
##      pv - vector, p-values of genes in the network
##      graph - matrix, estimated network in matrix format
##      covGene - vector, names of genes in the additional covariate of HMRF analysis
##      fixedGene - vector, names of genes to have fixed hidden states 1 in HMRF analysis
##      pthres_hmrf - threshold for p-values in HMRF initialization, default=0.05
##      iter - number of iteration, default=100
##      trimthres - threshold for Z score trimming, default=5
## Output:
##      result_DAWN - a list object representing the estimated statistics of genes
##          result_DAWN$Iupdate - updated I
##          result_DAWN$post - posterior distribution of I
##          result_DAWN$b0 - estimated b
##          result_DAWN$b1 - estimated c
##          result_DAWN$b2 - estimated coefficient for additional covariate (default=0 if additional covariate is not included)
##          result_DAWN$mu1 - estimated mean
##          result_DAWN$sigmas - estimated variance
DAWN_HMRF <- function(pv, graph, covGene, fixedGene, pthres_hmrf, iter, trimthres) {
    ## initialization
    zscore <- qnorm(1 - pv)
    if (missing(pthres_hmrf)) {
        pthres_hmrf <- 0.05
    }
    if (missing(iter)) {
        iter <- 100
    }
    if (missing(trimthres)) {
        trimthres <- 5
    }
    addCovGene <- hasArg(covindex)
    addFixedGene <- hasArg(fixindex)
    if (addCovGene) {
      covindex <- rep(0, length(genelist_PNS))
      covindex[(genelist_PNS) %in% covGene] <- 1
    }
    if (addFixedGene) {
      fixindex <- rep(0, length(genelist_PNS))
      fixindex[(genelist_PNS) %in% fixedGene] <- 1
      fixindex <- (zscore > trimthres) | fixindex + 0
    } else {
        fixindex <- (zscore > trimthres) + 0
    }
    
    Istart <- (pv < pthres_hmrf) + 0
    Istart[which(fixindex==1)] <- 1
    b1_1 <- 0
    b0_1 <- 0
    b2_1 <- 0
    Iupdate <- Istart
    mu1 <- mean(zscore[Istart == 1 & fixindex == 0])
    sigmas1 <- (sd(zscore[Istart == 0]))^2
    sigmas2 <- (sd(zscore[Istart == 1 & fixindex == 0]))^2
    posterior <- rep(0, length(zscore))
    
    for (iteri in 1:iter) {
        ## Algorithm2, step2(a) optimize b, c (and d)
        res <- optimize_params(graph, Iupdate, addCovGene, covindex, 20)
        b0 <- res[1]
        b1 <- res[2]
        b2 <- res[3]
        if (addCovGene) {
          cat("iteration", iteri-1, ": b0 =", b0_1,", b1 =", b1_1, ", b2 =", b2_1, ", mu1 =", mu1, ", sigmas1 =", sigmas1, "\n")
        } else {
          cat("iteration", iteri-1, ": b0 =", b0_1, ", b1 =", b1_1, ", mu1 =", mu1, ", sigmas1 =", sigmas1, "\n")
        }
        if (abs(b1_1 - b1) < 0.001 & abs(b0_1 - b0) < 0.001 & abs(b2_1 - b2) < 1e-05) {
          break
        }
        b1_1 <- b1
        b0_1 <- b0
        b2_1 <- b2
        ## Algorithm2, step2(b) update I
        Ibefore = Iupdate
        for (i in 1:length(pv)) {            
            new_vec <- compute_new(i, b0, b1, b2, Iupdate, graph[i, ], addCovGene, covindex[i])       
            new1 <- new_vec[1]
            new2 <- new_vec[2] 
            p1 <- dnorm(zscore[i], mu1 * Iupdate[i], sqrt(sigmas2 * Iupdate[i] + sigmas1 * (1 - Iupdate[i]))) /(1 + 
                exp(new2-new1))  ##sigmas2 is the alternative, sigmas1 is the null var
            p2 <- dnorm(zscore[i], mu1 * (1 - Iupdate[i]), sqrt(sigmas2 * (1 - Iupdate[i]) + sigmas1 * Iupdate[i])) /
                (1 + exp(new1-new2))
          
            if (Iupdate[i] == 1) 
                posterior[i] <- p1/(p1 + p2)
            if (Iupdate[i] == 0) 
                posterior[i] <- p2/(p1 + p2)
            if (p2 > p1) {
                Iupdate[i] <- 1 - Iupdate[i]
            }
            if (fixindex[i] != 0) {
              Iupdate[i] <- 1
            }
        }
        ## Algorithm2, step2(c) update parameters
        mu1 <- sum(posterior[fixindex == 0] * zscore[fixindex == 0])/sum(posterior[fixindex == 0])
        sigmas2 <- sum(posterior[fixindex == 0] * (zscore[fixindex == 0] - mu1)^2)/sum(posterior[fixindex == 0])
        sigmas1 <- sum((1 - posterior[fixindex == 0]) * (zscore[fixindex == 0])^2)/sum(1 - posterior[fixindex == 0])
        sigmas <- (sigmas1 * sum(posterior[fixindex == 0]) + sigmas2 * sum(1 - posterior[fixindex == 0]))/length(posterior)
        sigmas2 <- sigmas
        sigmas1 <- sigmas
    }
    result_DAWN <- list()
    result_DAWN$Iupdate <- Iupdate
    result_DAWN$post <- posterior
    result_DAWN$b0 <- b0
    result_DAWN$b1 <- b1
    result_DAWN$b2 <- b2
    result_DAWN$mu1 <- mu1
    result_DAWN$sigmas <- sigmas
    if (!check_b0(I = Iupdate, b0 = b0)) {
      cat('\n\nWARNING: DAWN identified a large number of risk genes. 
          Assumptions of the model may be false. 
          The set of risk genes likely contains many false positives.\n')
    }
    if (!check_b1(G = graph, I = Iupdate, b0 = b0, b1 = b1)) {
      cat('\n\nWARNING: Weak connectivity among risk genes in the input graph. 
          Assumptions of the model appear to be false. 
          The set of risk genes likely contains many false positives.\n')
    }
    return(result_DAWN)
}


##  Formating final results
##  Input:
##      finalposter - final posterior of hidden states I
##      kgene - genes in the graph
##      pv - p-values of genes in the graph
## Output:
##      report - a data frame containing information for genes and their corresponding posteriors, GLA p-values, and FDR
result_report <- function(finalposter, genes, pv) {
    usepost <- finalposter
    rankpost <- sort(usepost)
    localfdr <- rep(0, length(usepost))
    
    for (i in 1:length(localfdr)) {
        localfdr[i] <- mean(rankpost[1:i])
    }
    
    flocalfdr <- rep(0, length(localfdr))
    rankp <- rank(usepost, ties.method = "random")
    flocalfdr <- localfdr[rankp]
    
    report <- data.frame(genes, usepost, pv, flocalfdr)
    names(report) <- c("gene", "posterior", "GLA.pvalue", "FDR")
    return(report)
}

####### Supporting functions for PNS algorithm #######

## Generate a plot of square of correlation (R^2) vs lambda to help user choose a reasonable lambda
## Input:
##      expdata - dataframe, expression data
##      GLA_pvalue - dataframe, GLA p-values
##      pthres_pns - threshold for p-value screening
##      corthres - threshold for pairwise-correlation
## Output:
##      lambda_R2 - matrix, every row represents a lambda and R^2 value pair
##      a plot visualizing lambda_R2
#choose_lambda(expdata=expdata, GLA_pvalue=GLA_pvalue, pthres_pns=pthres_pns, corthres=corthres)
choose_lambda <- function(expdata, GLA_pvalue, pthres_pns, corthres) {
  cat('Running choose_lambda()...\n')
  lower_bound <- 0.05
  upper_bound <- 0.6
  step_size <- 0.02
  lambda_vec <- seq(lower_bound, upper_bound, step_size)
  cat('lower_bound =', lower_bound, '\n')
  cat('upper_bound =', upper_bound, '\n')
  cat('step_size =', step_size, '\n')
  lambda_R2 <- c()
  total_rounds <- length(lambda_vec)
  nround <- 0
  for (lambda in lambda_vec) {
    nround <- nround + 1
    cat('Running', nround, '/', total_rounds, 'round...\n')
    graphres_PNS <- PNS_algorithm(expdata = expdata, GLA_pvalue = GLA_pvalue, pthres_pns = pthres_pns, 
                                  corthres = corthres, lambda = lambda)
    R2 <- getR2(graph = graphres_PNS$graphfinal)
    lambda_R2 <- rbind(lambda_R2, c(lambda, R2))
    cat('Estimated: lambda =', lambda, ', R2 =', R2, '\n')
  }
  plot(lambda_R2[,1], lambda_R2[,2], xlab = 'lambda', ylab = 'R^2', type = 'l', 
       main = 'R^2 vs lambda')
}

## Compute the square of correlation (R^2) for PNS graph
## Input:
##      graph - estimated graph by PNS
## Output:
##      R2 - square of correlation for input graph
getR2 <- function(graph) {
  pk <- table(colSums(graph))
  if (names(pk)[1]=='0') {
    pk <- pk[-1]
  }
  stopifnot(length(pk) > 0)
  k <- as.numeric(names(pk))
  R2 <- cor(log(pk), log(k)) ^ 2
  return(R2)
}

## Generate a graph based on meinshausen's neighborhood selection
## Input:
##      pv - GLA p-values
##      cleaned_graph - a list object representing nodes in subnetwork
##      lambda - tuning parameter for lasso
## Output:
##      graphres - an list object representing the PNS network
##          graphres$graphfinal - estimated network in matrix format
##          graphres$pvfinal - p-values of genes in the network
##          graphres$zscoref - Z scores of genes in the network
meishausen_f <- function(pv, cleaned_graph, lambda) {
    if (missing(lambda)) {
        lambda <- 0.24
    }
    yexp <- cleaned_graph$yexp  ##expression data of core genes
    xexp <- cleaned_graph$xexp  ##expression data of real neighbors
    allnreal <- cleaned_graph$allnreal  ##indexes of real neighbors
    kscreeni <- cleaned_graph$kscreeni  ##indexes of core genes
    corm <- cleaned_graph$corm  ##correlation matrix for the whole expression data
    numcg <- length(kscreeni)
    numng <- length(allnreal)
    if (numng > 0) {
        xyexp <- rbind(yexp, xexp)
    }
    
    
    ## perform lasso regression for each subnetwork gene on each core gene
    edgef2 <- function(i) {
        testr1 <- glmnet(t(xyexp[-i, ]), t(yexp[i, ]), family = "gaussian", lambda = lambda)
        vectorr <- rep(0, numcg + numng)
        vectorr[-i] <- (c(abs(testr1$beta[, 1])) > 0) + 0
        return(vectorr)
    }
    fullm1 <- rbind(sapply(c(1:numcg), edgef2))
    
    coregenem1 <- 1 - (1 - fullm1[1:numcg, ]) * (1 - t(fullm1[1:numcg, ]))  ##edges among core genes
    if (numng > 0) {
        ngenem <- fullm1[c((numcg + 1):(numcg + numng)), ]  ##edges between core genes and real neighbors    
        pvfinal <- pv[c(kscreeni, allnreal)]
        graphfinal <- matrix(0, numcg + numng, numcg + numng)
        graphfinal[c(1:numcg), c(1:numcg)] <- coregenem1
        graphfinal[c(1:numcg), numcg + c(1:numng)] <- t(ngenem)
        graphfinal[numcg + c(1:numng), c(1:numcg)] <- ngenem
    } else {
        pvfinal <- pv[kscreeni]
        graphfinal <- matrix(0, numcg, numcg)
        graphfinal[c(1:numcg), c(1:numcg)] <- coregenem1
    }
    
    
    graphres <- list()
    graphres$graphfinal <- graphfinal  ##estimated network in matrix format
    graphres$pvfinal <- pvfinal  ##p-values of genes in the network
    return(graphres)
}

## Perform p-value screening, correlation screening, and neighbor retrieving (step 1, 2, and 3 in Algorithm 1)
## Input:
##      pv - GLA p-values
##      expdata - expression data
##      pthres_pns - threshold for p-value, default=0.1
##      corthres - threshold for pairwise-correlation, default=0.7
## Output:
##      cleanres - a list object representing nodes in subnetwork
##          cleanres$yexp - expression data of core genes
##          cleanres$xexp - expression data of real neighbors
##          cleanres$allnreal - indexes of real neighbors
##          cleanres$kscreeni - indexes of core genes
##          cleanres$corm - correlation matrix for the whole expression data
cleanup_f <- function(pv, expdata, pthres_pns, corthres) {
    corm <- abs(cor(t(expdata)))  ##correlation matrix
    
    ## p-value screening
    if (missing(pthres_pns)) {
        pthres_pns <- 0.1
    }
    if (missing(corthres)) {
        corthres <- 0.7
    }
    screeni <- c(1:length(pv))[pv < pthres_pns]  ##indexes of key genes (in all genes)
    screenpv <- pv[screeni]
    screencorm <- corm[screeni, screeni]
    
    ## correlation screening
    screengraph <- (screencorm > corthres) + 0 - diag(1, dim(screencorm)[1])
    screeni2 <- c(1:(dim(screencorm)[1]))[colSums(screengraph) > 0]  ##indexes of core genes(in key genes)
    screengraph2 <- screengraph[screeni2, screeni2]
    
    ## retrieving neighbors
    kscreeni <- screeni[screeni2]  ##indexes of core genes (in all genes)
    ograph <- (corm > corthres) + 0 - diag(1, dim(corm)[1])  ##graph only screened by correlation
    first.order.find <- function(i) {
        ki <- c(1:(dim(corm)[1]))[ograph[i, ] == 1]
        return(ki)
    }
    ab <- sapply(kscreeni, first.order.find)  ##find indexes of neighbors for each core gene
    if (length(ab) == 0) {
        allnreal <- c()
        xexp <- data.frame()
    } else {
        allneighbor <- c()
        for (i in 1:length(kscreeni)) {
            allneighbor <- c(allneighbor, ab[[i]])
        }
        alln <- unique(allneighbor)
        allnreal <- alln[-which(alln %in% kscreeni)]  ##indexes of real neighbors that are not core genes
        xexp <- expdata[allnreal, ]
    }   
    yexp <- expdata[kscreeni, ]
    
    cleanres <- list()
    cleanres$yexp <- yexp
    cleanres$xexp <- xexp
    cleanres$allnreal <- allnreal
    cleanres$kscreeni <- kscreeni
    cleanres$corm <- corm
    return(cleanres)
}


####### supporting functions for HMRF analysis #######

## Compute intermediate results for Algorithm2, step2(a)
## Input:
##      i - index of current gene
##      b0 - parameter b in Ising model
##      b1 - parameter c in Ising model
##      b2 - parameter d in Ising model (coefficient for additional covariate, default=0 if additional covariate is not used)
##      Iupdate - hidden states I
##      graph_i - connectivity of current gene to other genes in PNS (i'th row in PNS graph)
##      addCovGene - whether to use covGene information or not
##      covindex_i - if current gene is in covindex
## Output:
##      a vector of intermediate results (new1, new2)
compute_new <- function(i, b0, b1, b2, Iupdate, graph_i, addCovGene, covindex_i) {
  if (addCovGene) {
    new1 <- (b0 * Iupdate[i] + b1 * Iupdate[i] * t(graph_i) %*% Iupdate + b2 * Iupdate[i] * covindex_i)
    new2 <- (b0 * (1 - Iupdate[i]) + b1 * (1 - Iupdate[i]) * t(graph_i) %*% Iupdate + b2 * (1 - Iupdate[i]) * 
                  covindex_i)
  } else {
    new1 <- (b0 * Iupdate[i] + b1 * Iupdate[i] * t(graph_i) %*% Iupdate)
    new2 <- (b0 * (1 - Iupdate[i]) + b1 * (1 - Iupdate[i]) * t(graph_i) %*% Iupdate)
  }
  return(c(new1, new2))
}

## Estimate parameters in Ising model by maximizing the pseudo likelihood
## Input:
##      graph - graph estimated by PNS algorithm
##      Iupdate - hidden indicator vector I
##      addCovGene - whether to use covGene information or not
##      covindex - vector, indicates whether a gene is in the addtional covariate
##      times - number of iterations
## Output:
##      a vector of estimated paramter values - c(b0, b1, b2)
##        b0 - parameter b in Ising model
##        b1 - parameter c in Ising model
##        b2 - parameter d in Ising model (coefficient for additional covariate, default=0 if additional covariate is not used)
optimize_params <- function(graph1, Iupdate, addCovGene, covindex, times) {
    ## initialization
    b0_iter <- 0
    b1_iter <- 0
    b2_iter <- 0
    risknode.num <- Iupdate %*% graph1
    if (addCovGene) {
        for (k in 1:times) {
            b0_iter_1 <- optimize(plf_covGene, c(-20, 0), b1 = b1_iter, b2 = b2_iter, risknode.num = risknode.num, zupdate = Iupdate, 
                covindex = covindex, maximum = T)$maximum
            b1_iter_1 <- optimize(plf_covGene, c(0, 10), b0 = b0_iter_1, b2 = b2_iter, risknode.num = risknode.num, zupdate = Iupdate, 
                covindex = covindex, maximum = T)$maximum
            b2_iter_1 <- optimize(plf_covGene, c(0, 10), b1 = b1_iter_1, b0 = b0_iter_1, risknode.num = risknode.num, zupdate = Iupdate, 
                covindex = covindex, maximum = T)$maximum
            if (abs(b0_iter_1 - b0_iter) < 10^-5 & abs(b1_iter_1 - b1_iter) < 10^-5 & abs(b2_iter_1 - b2_iter) < 10^-5) {
                break
            }
            b0_iter <- b0_iter_1
            b1_iter <- b1_iter_1
            b2_iter <- b2_iter_1
        }
    } else {
        for (k in 1:times) {
            b0_iter_1 <- optimize(plf, c(-20, 0), b1 = b1_iter, risknode.num = risknode.num, zupdate = Iupdate, maximum = T)$maximum
            b1_iter_1 <- optimize(plf, c(0, 10), b0 = b0_iter_1, risknode.num = risknode.num, zupdate = Iupdate, maximum = T)$maximum
            if (abs(b0_iter_1 - b0_iter) < 10^-5 & abs(b1_iter_1 - b1_iter) < 10^-5) {
                break
            }
            b0_iter <- b0_iter_1
            b1_iter <- b1_iter_1
        }
    }
    
    return(c(b0_iter, b1_iter, b2_iter))
}

## Compute pseudo likelihood
## Input:
##      b0 - parameter b in Ising model
##      b1 - parameter c in Ising model
##      risknode.num - number of risk nodes connected to each node in the graph estimated by PNS algorithm
##      zupdate - hidden indicators I
## Output:
##      fvalue - pseudo likelihood value
plf <- function(b0, b1, risknode.num, zupdate) {
    fv.fun <- function(i) {
        new1 <- exp(b0 * zupdate[i] + b1 * zupdate[i] * risknode.num[i])
        new2 <- exp(b0 * (1 - zupdate[i]) + b1 * (1 - zupdate[i]) * risknode.num[i])
        fvalue <- log(new1/(new1 + new2))
    }
    fvalue <- sum(sapply(c(1:length(zupdate)), fv.fun))
    return(fvalue)
}

## Compute pseudo likelihood with additional covariate added
## Input:
##      b0 - parameter b in Ising model
##      b1 - parameter c in Ising model
##      b2 - parameter d in Ising model (coefficient for covGenes)
##      risknode.num - number of risk nodes connected to each node in the graph estimated by PNS algorithm
##      covindex - vector, indicates whether a gene is in the addtional covariate
##      zupdate - hidden indicator vector I
## Output:
##      fvalue - pseudo likelihood value
plf_covGene <- function(b0, b1, b2, risknode.num, covindex, zupdate) {
    fv.fun <- function(i) {
        new1 <- exp(b0 * zupdate[i] + b1 * zupdate[i] * risknode.num[i] + b2 * covindex[i] * zupdate[i])
        new2 <- exp(b0 * (1 - zupdate[i]) + b1 * (1 - zupdate[i]) * risknode.num[i] + b2 * covindex[i] * (1 - zupdate[i]))
        fvalue <- log(new1/(new1 + new2))
    }
    fvalue <- sum(sapply(c(1:length(zupdate)), fv.fun))
    return(fvalue)
}

## Check if b0 value is reasonable
## Input:
##      I - estimated hidden indicator vector
##      b0 - estimated b0
##      thres_b0 - threshold for evaluating b0
## Output:
##      pass - a boolean variable, True if b0 passes the test, False otherwise
check_b0 <- function(I, b0, thres_b0 = 0.05){
  ## Check if P(sum(I)>thres) < thres under the distribution that:
  ## P(sum(I)) ~ binom(n,exp(b0)/(1+exp(b0)))
  n <- length(I)
  q <- floor(n * thres_b0);
  p <- exp(b0) / (1 + exp(b0))
  tail_prob <- pbinom(q = q, size=n, prob = p, lower.tail = F) ##P(sum(I)>thres)
  pass <- (tail_prob < thres_b0)
  return(pass)
}

## Check if b1 value is reasonable
## Input:
##      G - adjacency matrix for graph
##      I - estimated hidden indicator vector
##      b0 - estimated b0
##      b1 - estimated b1
##      thres_b1 - threshold for evaluating b1
## Output:
##      pass - a boolean variable, True if b0 passes the test, False otherwise
check_b1 <- function(G, I, b0, b1, thres_b1 = 0.1) {
  ## Check if b1*I'GI is significant comparing to b0*I
  pass <- abs(b1*I%*%G%*%I / (b0*sum(I)))[1,1] > thres_b1
  return(pass)
}
  
 
