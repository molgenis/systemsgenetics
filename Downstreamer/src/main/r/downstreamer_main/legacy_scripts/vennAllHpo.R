setwd("D:\\UMCG\\Genetica\\Projects\\Depict2Pgs")



load("keyCisTransRare/key_gene_overlap_files_for_patrick.RData")

row.names(hpo.gwas.link)



str(key.genes)

gwasPthreshold <- 0.05/nrow(gene.pvalues)

gene <- rownames(gene.pvalues)[2] 
gene <- "ENSG00000179455"
gene <- "ENSG00000173452"
trait <- row.names(hpo.gwas.link)[1] 




keyEnrichedGwas <- read.delim("keyCisTransRare/enrichedKeyGenesGwas.txt")[,1]
str(keyEnrichedGwas)


geneBox <- sapply(rownames(gene.pvalues), function(gene){
  
  anyGwas <- FALSE
  anyKey <- FALSE
  anyHpo <- FALSE
  anyGwasKeyPair <- FALSE
  anyGwasHpoPair <- FALSE
  anyKeyHpoPair <- FALSE
  anyGWasKeyHpoTriple <- FALSE
  
  #row.names(hpo.gwas.link)
  for(trait in keyEnrichedGwas){
    
    relatedHpo <- hpo.gwas.link[trait,"Related HPO ID"]
    
    gwasSig <- gene.pvalues[gene,trait] <= gwasPthreshold
    isKey <- gene %in% key.genes[[trait]]
    
    isHpo <- hpo[gene,relatedHpo] > 0
    if(is.na(isHpo)){
      isHpo <- FALSE
    }
    
    if(gwasSig){
      anyGwas <- TRUE
    }
    if(isKey){
      anyKey <- TRUE
    }
    if(isHpo){
      anyHpo <- TRUE
    }
    if(isHpo & isKey){
      anyKeyHpoPair <- TRUE
    }
    if(gwasSig & isKey){
      anyGwasKeyPair <- TRUE
    }
    if(gwasSig & isHpo){
      anyGwasHpoPair <- TRUE
    }
    if(gwasSig & isHpo & isKey){
      anyGWasKeyHpoTriple <- TRUE
    }
    
    
  }
  
  if(anyGWasKeyHpoTriple){
    return("GWAS&KEY&HPO")
  }
  if(anyGwasKeyPair){
    return("GWAS&KEY")
  }
  if(anyGwasHpoPair){
    return("GWAS&HPO")
  }
  if(anyKeyHpoPair){
    return("KEY&HPO")
  }
  if(anyHpo){
    return("HPO")
  }
  if(anyKey){
    return("KEY")
  }
  if(anyGwas){
    return("GWAS")
  }
  
})

geneBox2 <- geneBox[!sapply(geneBox, is.null)]

geneBox3 <- do.call(c,geneBox2)

eulerData <- table(geneBox3)
eulerData2 <- as.vector(eulerData)
names(eulerData2) <- names(eulerData)
plot(euler(combinations = eulerData2, shape = "ellipse"))

eulerData3 <- eulerData2
names(eulerData3) <- gsub("_", "&", names(eulerData2))
pdf("keyCisTransRare/significantGwasKeyGenesHpo.pdf")
plot(euler(combinations = eulerData3, shape = "ellipse"), quantities = TRUE)
dev.off()

(604+84)/7465


(272+84)/(2308+628+84+272)


(276+81)/(2194+567+81+276)


(513+81)/(5570+513+81+276)


#sapply(row.names(hpo.gwas.link), function(trait){

i <- 1;

trait <- row.names(hpo.gwas.link)[i] 
print(trait)
pdf(paste0("keyCisTransRare/eulerKeyCisRare", trait, ".pdf"))  
relatedHpo <- hpo.gwas.link[trait,"Related HPO ID"]

gwasGenes <- row.names(gene.pvalues)[gene.pvalues[,trait] <= gwasPthreshold]
keyGenes <- key.genes[[trait]]
hpoGenes <- row.names(hpo)[hpo[,relatedHpo] > 0]

geneLists <- list(`GWAS genes` = gwasGenes, hpoGenes = hpoGenes, keyGenes = keyGenes)

plot(euler(geneLists, shape = "ellipse"), quantities = TRUE, main = hpo.gwas.link[trait,"GWAS Trait"])


dev.off()

i <- i + 1

})

str(hpo.gwas.link)
