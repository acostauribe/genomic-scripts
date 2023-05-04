#title: "PC-AiR"
#author: "Juliana Acosta-Uribe"
#date: '2022-JUL-08'

#Bioconductor: <http://www.bioconductor.org/packages/release/bioc/html/GENESIS.html>
#Tutorial: <http://bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/pcair.html>
#Citation: Gogarten, S.M., Sofer, T., Chen, H., Yu, C., Brody, J.A., Thornton, T.A., Rice, K.M., and Conomos, M.P. (2019). Genetic association testing using the GENESIS R/Bioconductor package. Bioinformatics. <doi:10.1093/bioinformatics/btz567>.


# 0. Install packages 

#install.packages('dplyr')
library(dplyr)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.15")
#BiocManager::install("GENESIS")
library(GENESIS)
#BiocManager::install("GWASTools")
library(GWASTools)
#https://rdrr.io/bioc/GWASTools/man/GWASTools-package.html
#BiocManager::install("SNPRelate")
library(SNPRelate)


# 1. Set up variables

# Set the prefix of your plink files
prefix="CLM-Esc-HiQual"
# Set the number of minor allelic frequencies you want to test
mafs = c(0.05, 0.1) 
# Set up the standard deviation limit and you want to use for outlier detection
sd.limit = 6
# Set the number of iterations for outlier removal
iterations = 5
# Set the number of principal components you want to use for outlier filtering
pcs = 6

# 2. Create a GDS file using SNPRelate

snpgdsBED2GDS(bed.fn = paste0(prefix,".bed"), 
              bim.fn = paste0(prefix,".bim"), 
              fam.fn = paste0(prefix,".fam"), 
              out.gdsfn = paste0(prefix,".gds"))

geno = GdsGenotypeReader(filename = paste0(prefix,".gds"))
genoData = GenotypeData(geno)
sample.set = getScanID(geno)


# 3. Run KING using SNP relate

# LD pruning is not recommended in KING, so run it with the original files
gds = snpgdsOpen(paste0(prefix,".gds"), 
                 allow.duplicate=TRUE, 
                 allow.fork=TRUE)
king = snpgdsIBDKING(gds,
                     type=c("KING-robust"))
KINGmat = kingToMatrix(king)


# 3.  Run PC-AiR on pruned SNPs

# Minor allele frequency (MAF) filtering of rare variants is recommended for PC-AiR (original paper recomends 5%)
for (m in mafs) { 
  treshold = m
  sample.maf = sample.set
  
  for (i in 1:iterations) {
    sample.iteration = sample.maf
  
    # Perform LD prunning to select a set of independent SNPs for analysis. 
    # We use the snpgdsLDpruning in the SNPRelate package, which returns a list of snp IDs
    snpset = snpgdsLDpruning(gds, 
                             method = "corr", 
                             slide.max.bp = 10e6, 
                             ld.threshold=  sqrt(0.1), 
                             verbose = FALSE,
                             sample.id = sample.iteration,
                             maf = treshold)
    pruned = unlist(snpset, 
                    use.names=FALSE)
    # Perform the PCA
    mypcair = pcair(genoData, 
                    kinobj = KINGmat, 
                    divobj = KINGmat,
                    snp.include = pruned,
                    sample.include = sample.iteration)
    # Plot it
    plot(mypcair, 
         main = paste0("MAF of ", m, " iteration # ", i))
    # The default is to plot PC values as black dots and blue pluses for individuals in the “unrelated subset” and “related subsets” respectively
    
    # Highlight outliers (shows which individuasl are beyond 6 standard deviations) 
    eigenvec = data.frame(mypcair[["vectors"]])
    pc1.mean = mean(eigenvec[,1])
    pc1.sd = sd(eigenvec[,1])
    pc2.mean = mean(eigenvec[,2])
    pc2.sd = sd(eigenvec[,2])
    plot(mypcair,  main = paste0("MAF of ", m, " iteration # ", i), cex = 1)
    abline(v = c(pc1.mean, (pc1.mean + sd.limit*pc1.sd),(pc1.mean - sd.limit*pc1.sd)), col=c("gray", "red", "red"))
    abline(h = c(pc2.mean, (pc2.mean + sd.limit*pc2.sd),(pc2.mean - sd.limit*pc2.sd)),col=c("gray", "red", "red"))
    legend("topleft", inset=.05, c("PC mean", paste(">", sd.limit, "SD", sep= " ")), fill=c("gray", "red"), horiz=TRUE, cex =0.8)
    
    # Generate a list of outliers of the top 10 PCs
    eigenvec = data.frame(mypcair[["vectors"]])
    for (pc in 1:pcs) {
      # Organize your dataframe
      pc.col = as.vector(paste0("PC",pc))
      names(eigenvec)[pc] = "test.column"
      # Do some math
      pc.mean = mean(eigenvec[,pc])
      pc.sd = sd(eigenvec[,pc])
      upper.limit = as.vector(pc.mean + sd.limit*pc.sd)
      lower.limit = as.vector(pc.mean - sd.limit*pc.sd)
      # Filter Rows
      outliers = eigenvec %>% filter(test.column>upper.limit | test.column<lower.limit)
      outlier.ids = data.frame(rownames(outliers))
      names(outlier.ids) = pc.col
      assign(paste(pc.col, "out", sep= "."), as.list(outlier.ids))
      
      # Restore the column name
      names(eigenvec)[pc] = pc.col
    }
    
    # Merge the outliers for each iteration
    all.outliers = mget(ls(pattern = 'PC')) 
    ###all.outliers = c(PC1.out, PC2.out, PC3.out, PC4.out, PC5.out, PC6.out, PC7.out, PC8.out, PC9.out, PC10.out)
    sample.remove = unique(as.integer(unlist(all.outliers)))
    assign(paste('removed.outliers.maf', m , "iteration", i, sep = "."), sample.remove) 
    
    # Select sample for next iteration
    sample.maf = sample.iteration[! sample.iteration %in% sample.remove]
    
    # Perform a scree plot to analyse variance explained by each PC
    eigenval = data.frame(mypcair[["values"]])
    colnames(eigenval) = c('values')
    eigenval_sum = sum(eigenval$values)
    eigenval = transform(eigenval, proportion = ((100*eigenval$values)/eigenval_sum))
    eigenval_vector = mypcair[["values"]]
    plot(eigenval$proportion,
         xlab = "Principal Component",
         ylab = "% of variance explained",
         main = paste0("MAF of ", m, " iteration # ", i))
    
    # Save the pcair object with a new name
    pcair.maf.it = paste('pcair.maf', m, "iteration", i, sep = ".") 
    assign(pcair.maf.it, mypcair)
  }
}
#close the gds
snpgdsClose(gds)

  
  
  
