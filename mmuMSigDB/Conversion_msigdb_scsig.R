rm(list = ls())


# ===================== Read gene set files =======================

#Author: Minghui Wang
read.gmt=function(filename){
  if(! file.exists(filename)) stop('File ',filename,' not available\n')
  dat=readLines(filename)
  n=length(dat)
  res=list(genesets=vector(mode = "list", length = n),geneset.names=vector(mode = "character", length = n),geneset.descriptions=vector(mode = "character", length = n))
  for(i in 1:n){
    s=strsplit(dat[i],'\t')[[1]]
    res$genesets[[i]]=s[-c(1:2)]
    res$geneset.names[i]=s[1]
    res$geneset.descriptions[i]=s[2]
  }
  names(res$genesets)=res$geneset.names
  res
}
#
write.gmt=function(obj,filename){
  conn=file(filename,'w')
  for(i in 1:length(obj$genesets)){
    cat(obj$geneset.names[i],obj$geneset.descriptions[i],obj$genesets[[i]],file=conn,sep='\t')
    cat('\n',file=conn)
  }
  close(conn)
  return(invisible())
}

SCSig <- read.gmt(filename = "./msigdb_v7/scsig.all.v1.0.symbols.gmt")
SCSig


#====================== Build Homolog Conversion Table =====================

library("biomaRt")

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

human2mouse = getLDS(attributes = c("hgnc_symbol","ensembl_gene_id","entrezgene_id"), mart = human, 
                     attributesL = c("mgi_symbol","ensembl_gene_id","entrezgene_id"), martL = mouse, uniqueRows=T)



#########################################################################
# ====================== Conversion for C5.BP ==========================
gs <-SCSig

library(progress)
pb <- progress_bar$new(total = length(gs$genesets))
for (gsid in 1:length(gs$genesets)){
  
  gvec<-gs$genesets[[gsid]]
  # go through all genes in one gene set
  for (gene in gvec){
    gene.rowsel      <-  human2mouse$HGNC.symbol==gene
    symbol.converted <-( human2mouse[gene.rowsel, ]$MGI.symbol ) [1]  # if multiple hits exist, use first one. (just in case)
    if ( symbol.converted!="" & !is.na(symbol.converted) ){ # if conversion found
      gvec[gvec==gene] <- symbol.converted
    }else{ # or removed the unmatched item
      gvec<-gvec[gvec!=gene]
    }
  }
  gs$genesets[[gsid]] <- gvec
  pb$tick()
}
gs
write.gmt(obj = gs, filename = "SCSig.v1.0.symbols.gmt")
rm(gs)
