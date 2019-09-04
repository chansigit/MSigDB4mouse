library("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV2 = getLDS(attributes = c("mgi_symbol"), 
                 filters = "mgi_symbol", 
                 values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

mart = useMart('ensembl') 
listDatasets(mart)

listAttributes(human)->hsAttr


human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

human2mouse = getLDS(attributes = c("hgnc_symbol","ensembl_gene_id","entrezgene_id"), mart = human, 
                     attributesL = c("mgi_symbol","ensembl_gene_id","entrezgene_id"), martL = mouse, uniqueRows=T)

mouse2human = getLDS(attributes = c("mgi_symbol","ensembl_gene_id","entrezgene_id"),  mart = mouse,
                     attributesL= c("hgnc_symbol","ensembl_gene_id","entrezgene_id"), martL = human, uniqueRows=T)
