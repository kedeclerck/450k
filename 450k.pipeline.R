

Rscript.dir <- c("D:/Dropbox/Dropbox/Phd/Rscripts/450k calcua/")


# import 450k methylation data --------------------------------------------------------------

source(file.path(Rscript.dir, "data.import.R"))

rnb.set <- data.import(idat.dir="D:/Dropbox/Dropbox/Phd/R analysis/Idatfiles",
                       sample.annotation="D:/Dropbox/Dropbox/Phd/R analysis/Sample_annotation.txt",
                       identifiers.column = "Sample_ID",
                       export=F,
                       output.dir = "D:/Dropbox/Dropbox/Phd/R analysis/output")


# preprocessing (normalization + filtering) --------------------------------------------------------

source(file.path(Rscript.dir, "preprocessing.R"))

rnb.set <- preprocess(rnb.set,
                      prefiltering=list(SNP="3", detPVal=0.01, missing=0),
                      normalization=list(method="bmiq", bgcorr.method="none"),
                      postfiltering=list(sex=TRUE, context=TRUE),
                      export=TRUE,
                      output.dir="D:/Dropbox/Dropbox/Phd/R analysis/output")


# exploration --------------------------------------------------------------------------------

source(file.path(Rscript.dir, "explore.R"))

explore(rnb.set,
        traits=c("Sample_group", "Sentrix_ID", "Sentrix_Position"),
        mycols=c("steelblue", "orange2", "limegreen", "pink2", "red2", "black", "gray", "darkmagenta", "burlywood4","aquamarine3","deeppink2"),
        density.plot=list(main="",
                          xlab="Beta-values",
                          lwd=3),
        PCA.plot=list(screeplot=TRUE, 
                      xlab="Principal component 1",
                      ylab="Principal component 2", 
                      main="", 
                      cex=1,
                      pch=16,
                      xlim=NA,
                      ylim=NA,
                      dots=TRUE,
                      labels=c(TRUE, "Sample_ID"),
                      PC.assoc=TRUE),
        cluster=list(distance="euclidean",
                     linkage="average",
                     xlab="",
                     ylab="Height",
                     main=""),
        export=T,
        output.dir="D:/Dropbox/Dropbox/Phd/R analysis/output")


# batch correction -----------------------------------------------------------------------

source(file.path(Rscript.dir, "batch.correction.R"))

combat.methyl <- batch(rnb.set,
                       GOI=list(GOI.column="Sample_group", type="factor"),
                       covariates.column=list(factor="Sentrix_Position", cont=NULL),
                       batch.column="Sentrix_ID",
                       export=F,
                       output.dir="D:/Dropbox/Dropbox/Phd/R analysis/output")


# differential methylation -----------------------------------------------------------------

source(file.path(Rscript.dir, "diff.meth.R"))

load("D:/Dropbox/Dropbox/Phd/Rscripts/450kshort.RData")
source("D:/Dropbox/Dropbox/Phd/Rscripts/HM450.annotate.R")
source("D:/Dropbox/Dropbox/Phd/Rscripts/annotate.450kDMR.R")

# create design matrix
group.limma <- factor(pheno(rnb.set)$Sample_group)

design <- model.matrix(~0+group.limma)
colnames(design) <- c(levels(group.limma))

# create contrast matrix
contrast.matrix <- makeContrasts(Atherosclerosis-Flavanol,
                                 levels=design)

diff.methyl <- diff.meth(rnb.set,
                         design=design,
                         contrast.matrix=contrast.matrix,
                         DMRcate.arg=list(lambda=1000, C=2,
                                          pcutoff=0.0001,
                                          betacutoff=NULL),
                         comb.p=F,
                         export=F,
                         output.dir="D:/Dropbox/Dropbox/Phd/R analysis/output")
                      

# genomic enrichment -------------------------------------------------------------------

source(file.path(Rscript.dir, "gen.enrich.R"))

source("D:/Dropbox/Dropbox/Phd/Rscripts/FisherBarplot.R")

gen.enrich.methyl <- gen.enrich(limma.out=diff.meth.output$limma$`Atherosclerosis - Flavanol`,
                                p.cutoff=list(P.Value=0.05 ,adj.P.Val=0.05),
                                beta.cutoff=0,
                                dir=FALSE,
                                enrichment=list(gene=TRUE,
                                                CGI=TRUE,
                                                states=TRUE,
                                                TFBS=TRUE),
                                bars=list(ylab="Relative overlap",
                                          col=c("gray47", "limegreen"),
                                          asterisk=TRUE),
                                export=TRUE,
                                output.dir="D:/Dropbox/Dropbox/Phd/R analysis/output")


                       




         