require(hash)
require(assertthat)
require(biomaRt)

#load in whatever gene tools is available
load("gene-tools-kevinl1-2016-03-29.RData")

.control <- setClass("control", representation(in.place = "logical",
  verbose = "logical", place.holder = "logical", enforce.uniq = "logical"),
  prototype(in.place = T, verbose = T, place.holder = T, enforce.uniq = F))

.convert.list2control <- function(lis){
  con = .control()
  
  if(!is.null(lis$in.place)) con@in.place = lis$in.place
  if(!is.null(lis$verbose)) con@verbose = lis$verbose
  if(!is.null(lis$place.holder)) con@place.holder = lis$place.holder
  if(!is.null(lis$enforce.uniq)) con@enforce.uniq = lis$enforce.uniq

  con
}

convert.naming <- function(gene.vec, from, to, control = list(
  in.place = T, verbose = T, place.holder = T, enforce.uniq = F)){

  vec = c(from, to)
  assert_that(all(vec %in% c("symbol", "entrez", "ensembl")))
  vec = sapply(vec, function(x){
    if(x == "symbol") return("hgnc_symbol")
    if(x == "entrez") return("entrezgene")
    if(x == "ensembl") return("ensembl_gene_id")
  })

  con = .convert.list2control(control)

  ensembl = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

  #need to account for different symbol names
  if(from == "symbol"){
    gene.idx = numeric(length(gene.vec))
    for(i in 1:length(gene.idx)){
      gene.idx[i] = eval(parse(text = paste0("gene.tools$hash$", gene.vec[i])))
    }
  
    #include ALL genes
    gene.vec = c(names(gene.tools$syn.list)[gene.idx], unlist(gene.tools$syn.list[gene.idx]))
    gene.vec = as.character(gene.vec)
  }

  gene.info = getBM(attributes = vec, filters = vec[1], values = gene.vec, mart = ensembl)

  #in.place means that the ordering of the outputs matters
  #TO FINISH
  #question: how to handle multiplicity. how to handle non-uniquness?
  if(con@in.place){
  
  } else {
    idx = which(colnames(gene.info) == vec[2])
    unique(gene.info[,idx])
  }

}
