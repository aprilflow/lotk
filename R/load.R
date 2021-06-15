#### Load script for most RNA based analyses

## Tidyverse packages, could be replaced by library(tidyverse)
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(scales)
library(stringr)

## Enhances to tidyverse
library(tidylog)
library(cowplot)
library(ggrepel)
library(ggExtra)
library(ggdendro)
library(ggforce)

## For differential gene transcript expression analysis
library(limma)
library(edgeR)
library(DESeq2)
library(goseq)
library(GO.db)

## Other visualization tools
library(Rtsne)
library(umap)
library(dendsort)
library(UpSetR)
library(venn)
library(patchwork)
library(pheatmap)

## Color palettes
library(RColorBrewer)
library(colorspace)
library(viridis)


library(data.table)

#####
options(tibble.width = Inf)

##### For do_DE scripts
# To get metadata for samples
getSampleInfo <- function(dataDir) {
  return(read.table(file=paste0(dataDir, "/sampleInfo.txt"), sep="\t", header=TRUE))
}
# To get list of genes (with annotations)
getGenes <- function(dataDir) {
  data_df <- read.table(paste0(dataDir, "/featureCounts/hg19/counts_all.txt"), header=TRUE)
  data_df <- data_df[1:(dim(data_df)[1]-92),]
  genes_df <- data_df[,1:6]
  return(genes_df)
}
# To get transcript quantification in count form
# Transcript quants were done with featureCounts from the subread package
getCounts <- function(dataDir) {
  sampleInfo <- getSampleInfo(dataDir)
  data_df <- read.table(paste0(dataDir, "/featureCounts/hg19/counts_all.txt"), header=TRUE)
  # Exclude ERCCs
  data_df <- data_df[1:(dim(data_df)[1]-92),]
  counts_df <- data_df[7:dim(data_df)[2]]
  names(counts_df) <- as.vector(sampleInfo$SampleId)
  return(counts_df)
}
# To get perturbation groups info
getGrouping <- function(dataDir, exclude) {
  sampleInfo <- getSampleInfo(dataDir)
  sampleInfo <- sampleInfo[-which(sampleInfo$SampleId %in% exclude),]
  N <- dim(sampleInfo)[1]
  grouping <- vector(length=N)
  for (i in 1:N) {
    sample_name <- as.character(sampleInfo$Notes[i])
    l <- as.vector(strsplit2(sample_name, "_"))
    if ("V5" %in% l) {
      l <- l[-which(l=="V5")]
    }
    grouping[i] <- paste(l[1:(length(l)-3)], collapse = "_")
  }
  return(grouping)
}


#####
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  g <- ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_smooth(method = "lm", col = "red") +
    # geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) +
    labs(subtitle = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                          "\nn =", dim(fit$model)[1],
                          "\nIntercept =",signif(fit$coef[[1]],5 ),
                          "\nSlope =",signif(fit$coef[[2]], 5),
                          "\np-value =",signif(summary(fit)$coef[2,4], 5)))
  
  return(g)
  
}

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

## y is matrix of rpkm/fpkm values (rows = genes, columns = samples)

fpkm_to_tpm <- function ( fpkm ) {
    tpm <- exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
    return(tpm)
}

log_fpkm_to_tpm <- function ( fpkm ) {
    tpm <- log2(fpkm) - log2(sum(fpkm)) + log2(1e6)
    return(tpm)
}


#### GO Analysis

do_goseq_DE <- function (data, perturbation, control, direction = NA, test = "GO") {
    foo <- data %>%
        filter(Perturbation == perturbation) %>%
        filter(Control == control)
    if (is.na(direction)) {
        genes = as.integer(p.adjust(foo$PValue[foo$logFC!=0],
                                method="BH")<.05)
        names(genes)=foo[foo$logFC!=0,]$Geneid        
    } else if (direction=="up") {
        genes = as.integer(p.adjust(foo$PValue[foo$logFC>0],
                                method="BH")<.05)
        names(genes)=foo[foo$logFC>0,]$Geneid
    } else if (direction=="down") {
        genes = as.integer(p.adjust(foo$PValue[foo$logFC<0],
                                method="BH")<.05)
        names(genes)=foo[foo$logFC<0,]$Geneid
    }
    ## Fit the Probability Weighting Function (PWF)
    pwf = nullp(genes, "hg19", "geneSymbol")
    print(head(pwf))
    ## Fetch relationship between gene symbols and GO categories
    ## Then use Wallenius approximation
    if (test == "GO") {
        GO.wall=goseq(pwf, "hg19", "geneSymbol")
        return(GO.wall)
    } else if (test == "MF") {
        GO.MF = goseq(pwf, "hg19", "geneSymbol", test.cats=c("GO:MF"))
        return(GO.MF)
    } else if (test == "BP") {
        GO.BP = goseq(pwf, "hg19", "geneSymbol", test.cats=c("GO:BP"))
        return(GO.BP)
    } else if (test == "hallmarks") {
        hallmarks <- read.table("/fh/fast/berger_a/grp/bergerlab_shared/Projects/RNA_eVIP/mSigDB/h.all.v6.2.symbols.gmt", sep = "\t", fill = TRUE)
        names(hallmarks) <- c("Category", "Gene", seq(200))
        hallmarks <- data.frame(hallmarks) %>%
            gather("foo", "Gene", 3:202) %>%
            dplyr::select("Category", "Gene") %>%
            filter(Gene != "")
        GO.hallmarks <- goseq(pwf, gene2cat = hallmarks)
        GO.hallmarks$ontology <- "mSigH"
        return(GO.hallmarks)
        }
    }

