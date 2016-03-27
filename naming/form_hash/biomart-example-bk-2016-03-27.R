rm(list=ls())

#download biomaRt from Bioconductor
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")

# load biomaRt 
library(biomaRt)
ensembl=useMart(biomart = "ensembl",dataset="hsapiens_gene_ensembl")

# this lists all fields you can filter your data on
availble.filters = listFilters(ensembl)
# this lists all fields you can put out
availble.attributes = listAttributes(ensembl)

# here I am looking for information on a gene for which I have the ensembl gene id
gene.id=c("ENSG00000128833")
# I want information on all transcripts and which parts of the exons these transcripts are using. My filter is based on ensembl_gene_id, and the value I am 
# selecting is the gene id above
gene.info=getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id","ensembl_exon_id","exon_chrom_start","exon_chrom_end"),
              filters="ensembl_gene_id",values=gene.id,mart=ensembl)

# in addition I want the hgnc gene name
gene.info=getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id","ensembl_exon_id","exon_chrom_start","exon_chrom_end","hgnc_symbol"),
                filters="ensembl_gene_id",values=gene.id,mart=ensembl)
# you name get an error stating that you are selecting from multiple attribute parts
# to solve this you do a two step query and then merge the results
gene.info.s1=getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id","ensembl_exon_id","exon_chrom_start","exon_chrom_end"),
                filters="ensembl_gene_id",values=gene.id,mart=ensembl)
gene.info.s2=getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                   filters="ensembl_gene_id",values=gene.id,mart=ensembl)

## finally if you need alias and previous names you should go to
# http://biomart.genenames.org/martform/#!/default/HGNC?datasets=hgnc_gene_mart (ensembl currently does not have this ability)
# here you can bulk upload gene names and then select the filter alias and previous names (2 runs) to get all possible known names for a gene
# select attributes "approved symbol" "alias symbol" and "previous symbol"

gene.id=c("CHD8","DUPLIN", "AUTS18", "KIAA1564", "HELSNF1")
gene.info=getBM(attributes=c("ensembl_gene_id","hgnc_id","hgnc_symbol"),
                   filters="hgnc_symbol",values=gene.id,mart=ensembl)
