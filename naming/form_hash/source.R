#the goal is to form a hash table and list of all the gene synonyms

library(hash)

setwd("~/kathryn.git.all/central.git/naming/form_hash")
files = list.files()
idx = grep("^hgnc.*.txt", files)

dat = vector("list", length(idx))
for(i in 1:length(idx)){
  dat[[i]] = read.csv(files[idx[i]], sep = "\t")
}

dat = do.call(rbind, dat)

#determine unique "approved symbols"
tmp.gene = as.character(dat[,1])
uniq.gene = sort(unique(tmp.gene))
gene.list = vector("list", length(uniq.gene))
names(gene.list) = uniq.gene

#create hash table and start filling it in
has = hash()

for(i in 1:length(uniq.gene)){
  idx = which(dat[,1] == uniq.gene[i])
  
  tmp.dat = dat[idx, c(2:3)]
  tmp.dat = apply(tmp.dat, 2, as.character)
  alt.names = sort(unique(as.vector(tmp.dat)))

  has[[uniq.gene[i]]] = i

  #remove blanks  
  idx.alt = which(sapply(alt.names, nchar) == 0)
  if(length(idx.alt)) alt.names = alt.names[-idx.alt]

  if(length(alt.names)){ 
    for(j in 1:length(alt.names)) has[[alt.names[j]]] = i
  }

  gene.list[[i]] = alt.names

  if(i %% floor(length(uniq.gene)/10) == 0) cat('*')
}

gene.tools = list(hash = has, syn.list = gene.list)
save(gene.tools, file = paste0("~/kathryn.git.all/central.git/naming/gene-tools-kevinl1-", Sys.Date(), ".RData"))
