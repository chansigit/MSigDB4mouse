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

c2.cp <- read.gmt(filename = "./msigdb_v7/c2.cp.v7.0.symbols.gmt")
c5.bp <- read.gmt(filename = "./msigdb_v7/c5.bp.v7.0.symbols.gmt")
c7.all<- read.gmt(filename = "./msigdb_v7/c7.all.v7.0.symbols.gmt")
#length(c2.cp$genesets)


#====================== Build Homolog Conversion Table =====================

library("biomaRt")

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

human2mouse = getLDS(attributes = c("hgnc_symbol","ensembl_gene_id","entrezgene_id"), mart = human, 
                     attributesL = c("mgi_symbol","ensembl_gene_id","entrezgene_id"), martL = mouse, uniqueRows=T)



# ====================== Conversion Quality Check ==========================
#msigdb.gene.symbols <- unique(
#  c(unique(unlist(c2.cp$genesets)),
#    unique(unlist(c5.bp$genesets)),
#    unique(unlist(c7.all$genesets))) )
#msigdb.gene.symbols

#gene.symbol.missed <-msigdb.gene.symbols[! msigdb.gene.symbols %in% human2mouse$HGNC.symbol]
#which(human2mouse$MGI.symbol, human2mouse$MGI.symbol==NA)



# ====================== Conversion for C2.CP ==========================
gs <-c2.cp

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

write.gmt(obj = gs, filename = "c2.cp.mmu.v7.0.symbols.gmt")
rm(gs)


# ====================== Conversion for C5.BP ==========================
gs <-c5.bp

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

write.gmt(obj = gs, filename = "c5.bp.mmu.v7.0.symbols.gmt")
rm(gs)

# ====================== Conversion for C7.ALL ==========================
gs <-c7.all

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

write.gmt(obj = gs, filename = "c7.all.mmu.v7.0.symbols.gmt")
rm(gs)

