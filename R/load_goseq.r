#### GO Analysis
library(goseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GO.db)

do_goseq_DE <- function(data, perturbation, control,
                        direction = NA, test = "GO") {
    foo <- data %>%
        filter(Perturbation == perturbation) %>%
        filter(Control == control)
    if (is.na(direction)) {
        genes <- as.integer(p.adjust(foo$PValue[foo$logFC != 0],
            method = "BH"
        ) < .05)
        names(genes) <- foo[foo$logFC != 0, ]$Geneid
    } else if (direction == "up") {
        genes <- as.integer(p.adjust(foo$PValue[foo$logFC > 0],
            method = "BH"
        ) < .05)
        names(genes) <- foo[foo$logFC > 0, ]$Geneid
    } else if (direction == "down") {
        genes <- as.integer(p.adjust(foo$PValue[foo$logFC < 0],
            method = "BH"
        ) < .05)
        names(genes) <- foo[foo$logFC < 0, ]$Geneid
    }
    ## Fit the Probability Weighting Function (PWF)
    pwf <- nullp(genes, "hg19", "geneSymbol")
    print(head(pwf))
    ## Fetch relationship between gene symbols and GO categories
    ## Then use Wallenius approximation
    if (test == "GO") {
        go_wall <- goseq(pwf, "hg19", "geneSymbol")
        return(go_wall)
    } else if (test == "MF") {
        go_mf <- goseq(pwf, "hg19", "geneSymbol", test.cats = c("GO:MF"))
        return(go_mf)
    } else if (test == "BP") {
        go_bp <- goseq(pwf, "hg19", "geneSymbol", test.cats = c("GO:BP"))
        return(go_bp)
    } else if (test == "hallmarks") {
        hallmarks <- read.table(
            paste0(
                "/fh/fast/berger_a/grp/bergerlab_shared/Projects/RNA_eVIP/",
                "mSigDB/h.all.v6.2.symbols.gmt"
            ),
            sep = "\t", fill = TRUE
        )
        names(hallmarks) <- c("Category", "Gene", seq(200))
        hallmarks <- data.frame(hallmarks) %>%
            gather("foo", "Gene", 3:202) %>%
            dplyr::select("Category", "Gene") %>%
            filter(Gene != "")
        go_hallmarks <- goseq(pwf, gene2cat = hallmarks)
        go_hallmarks$ontology <- "mSigH"
        return(go_hallmarks)
    }
}