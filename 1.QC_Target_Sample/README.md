# 1. QC of Target Data

# 1.0 Obtaining the target data

Target data consist of individual-level genotype-phenotype data, no sum stats!

!!! note
    We will need `plink` for this setp.
    
  `module load plink`

# 1.01 Running the GWAS?

*TBC*

*unsure about what sort of data MoBa provides, if GWAS needs to ran, see PhD scripts*

# 1.1 QC of Target data

Below are the QC steps that comprise the QC checklist for the Target data.

# 1.11 Sample size

Disc and Target samples should be of a similar n.

# 1.12 Genome build

As stated in the Discovery data section, check that the genome build for our Disc and Target data is the same, as it should be.

# 1.13 Standard GWAS QC

The Target data must be quality controlled to at least the standards 
implemented in GWAS studies, e.g. removing SNPs with low genotyping rate, 
low minor allele frequency, out of Hardy-Weinberg Equilibrium, removing
individuals with low genotyping rate 
(see [Marees et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/)).

The following `plink` command applies some of these QC metrics to the target data:


```bash
plink \
    --bfile EUR \
    --maf 0.01 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
    --write-snplist \
    --make-just-fam \
    --out EUR.QC
```

Each of the parameters corresponds to the following

| Parameter | Value | Description|
|:-:|:-:|:-|
| bfile | EUR | Informs `plink` that the input genotype files should have a prefix of `EUR` |
| maf | 0.01 | Removes all SNPs with minor allele frequency less than 0.01. Genotyping errors typically have a larger influence on SNPs with low MAF. Studies with large sample sizes could apply a lower MAF threshold|
| hwe | 1e-6 | Removes SNPs with low P-value from the Hardy-Weinberg Equilibrium Fisher's exact or chi-squared test. SNPs with significant P-values from the HWE test are more likely affected by genotyping error or the effects of natural selection. Filtering should be performed on the control samples to avoid filtering SNPs that are causal (under selection in cases). When phenotype information is included, plink will automatically perform the filtering in the controls. |
| geno | 0.01 | Excludes SNPs that are missing in a high fraction of subjects. A two-stage filtering process is usually performed (see [Marees et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/)).|
| mind | 0.01 | Excludes individuals who have a high rate of genotype missingness, since this may indicate problems in the DNA sample or processing. (see [Marees et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/) for more details).|
| make-just-fam | - | Informs `plink` to only generate the QC'ed sample name to avoid generating the .bed file.  |
| write-snplist | - | Informs `plink` to only generate the QC'ed SNP list to avoid generating the .bed file. |
| out | EUR.QC | Informs `plink` that all output should have a prefix of `EUR.QC` |


!!! note
    Normally, we can generate a new genotype file using the new sample list.
    However,  this will use up a lot of storage space. Using `plink`'s
    `--extract`, `--exclude`, `--keep`, `--remove`, `--make-just-fam` and `--write-snplist` functions, we can work 
    solely on the list of samples and SNPs without duplicating the 
    genotype file, reducing the storage space usage.  
    
Very high or low heterozygosity rates in individuals could be due to DNA contamination or to high levels of inbreeding. Therefore, samples with extreme heterozygosity are typically removed prior to downstream analyses. 

First, we perform pruning to remove highly correlated SNPs:

```bash
plink \
    --bfile EUR \
    --keep EUR.QC.fam \
    --extract EUR.QC.snplist \
    --indep-pairwise 200 50 0.25 \
    --out EUR.QC
```

Each of the parameters corresponds to the following

| Parameter | Value | Description|
|:-:|:-:|:-|
| bfile | EUR | Informs `plink` that the input genotype files should have a prefix of `EUR` |
| keep | EUR.QC.fam | Informs `plink` that we only want to use samples in `EUR.QC.fam` in the analysis |
| extract | EUR.QC.snplist | Informs `plink` that we only want to use SNPs in `EUR.QC.snplist` in the analysis |
|indep-pairwise| 200 50 0.25 | Informs `plink` that we wish to perform pruning with a window size of 200 variants, sliding across the genome with step size of 50 variants at a time, and filter out any SNPs with LD $r^2$ higher than 0.25|
| out | EUR.QC | Informs `plink` that all output should have a prefix of `EUR.QC` |


This will generate two files 1) **EUR.QC.prune.in** and 2) **EUR.QC.prune.out**. All SNPs within **EUR.QC.prune.in** have a pairwise $r^2 < 0.25$. 


Heterozygosity rates can then be computed using `plink`:

```bash
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.fam \
    --het \
    --out EUR.QC
```

This will generate the **EUR.QC.het** file, which contains F coefficient estimates for assessing heterozygosity.
We will remove individuals with F coefficients that are more than 3 standard deviation (SD) units from the mean, which can be performed using the following `R` command (assuming that you have R downloaded, then you can open an `R` session by typing `R` in your terminal):


  ```R
  dat <- read.table("EUR.QC.het", header=T) # Read in the EUR.het file, specify it has header
  m <- mean(dat$F) # Calculate the mean  
  s <- sd(dat$F) # Calculate the SD
  valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
  write.table(valid[,c(1,2)], "EUR.valid.sample", quote=F, row.names=F) # print FID and IID for valid samples
  q() # exit R
  ```




# 1.14 Ambiguous SNPs

These were removed during the Discovery data QC.

# 1.15 Mismatching SNPs

SNPs that have mismatching alleles reported in the Disc and Target data may be resolvable by strand-flipping the alleles to their complementary alleles in e.g. the Target data, such as for a SNP with A/C in the Disc data and G/T in the Target. This can be achieved with the following steps:


!!! note
    Most PGS software will perform strand-flipping automatically, thus this step is usually not required.


1\. Load the bim file, the summary statistic and the QC SNP list into R

```R
# Read in bim file
bim <- read.table("EUR.bim")
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
# Read in QCed SNPs
qc <- read.table("EUR.QC.snplist", header = F, stringsAsFactors = F)
# Read in the GWAS data
height <-
  read.table(gzfile("Height.QC.gz"),
    header = T,
    stringsAsFactors = F, 
    sep="\t")
# Change all alleles to upper case for easy comparison
height$A1 <- toupper(height$A1)
height$A2 <- toupper(height$A2)
bim$B.A1 <- toupper(bim$B.A1)
bim$B.A2 <- toupper(bim$B.A2)
```

2\. Identify SNPs that require strand flipping 

```R
# Merge summary statistic with target
info <- merge(bim, height, by = c("SNP", "CHR", "BP"))
# Filter QCed SNPs
info <- info[info$SNP %in% qc$V1,]
# Function for finding the complementary allele
complement <- function(x) {
  switch (
    x,
    "A" = "T",
    "C" = "G",
    "T" = "A",
    "G" = "C",
    return(NA)
  )
}
# Get SNPs that have the same alleles across base and target
info.match <- subset(info, A1 == B.A1 & A2 == B.A2)
# Identify SNPs that are complementary between base and target
info$C.A1 <- sapply(info$B.A1, complement)
info$C.A2 <- sapply(info$B.A2, complement)
info.complement <- subset(info, A1 == C.A1 & A2 == C.A2)
# Update the complementary alleles in the bim file
# This allow us to match the allele in subsequent analysis
complement.snps <- bim$SNP %in% info.complement$SNP
bim[complement.snps,]$B.A1 <-
  sapply(bim[complement.snps,]$B.A1, complement)
bim[complement.snps,]$B.A2 <-
  sapply(bim[complement.snps,]$B.A2, complement)
```



3\. Identify SNPs that require recoding in the target (to ensure the coding allele in the target data is the effective allele in the base summary statistic)

```R
    # identify SNPs that need recoding
    info.recode <- subset(info, A1 == B.A2 & A2 == B.A1)
    # Update the recode SNPs
    recode.snps <- bim$SNP %in% info.recode$SNP
    tmp <- bim[recode.snps,]$B.A1
    bim[recode.snps,]$B.A1 <- bim[recode.snps,]$B.A2
    bim[recode.snps,]$B.A2 <- tmp

    # identify SNPs that need recoding & complement
    info.crecode <- subset(info, A1 == C.A2 & A2 == C.A1)
    # Update the recode + strand flip SNPs
    com.snps <- bim$SNP %in% info.crecode$SNP
    tmp <- bim[com.snps,]$B.A1
    bim[com.snps,]$B.A1 <- as.character(sapply(bim[com.snps,]$B.A2, complement))
    bim[com.snps,]$B.A2 <- as.character(sapply(tmp, complement))
    
    # Output updated bim file
    write.table(
        bim[,c("SNP", "B.A1")],
        "EUR.a1",
        quote = F,
        row.names = F,
        col.names = F,
        sep="\t"
    )
    ```

=== "With data.table and magrittr"

    ```R
    # identify SNPs that need recoding
    recode.snps <- info[B.A1==A2 & B.A2==A1, SNP]
    # Update the bim file
    bim[SNP %in% recode.snps, c("B.A1", "B.A2") :=
            list(B.A2, B.A1)]

    # identify SNPs that need recoding & complement
    com.recode <- info[sapply(B.A1, complement) == A2 &
                        sapply(B.A2, complement) == A1, SNP]
    # Now update the bim file
    bim[SNP %in% com.recode, c("B.A1", "B.A2") :=
            list(sapply(B.A2, complement),
                sapply(B.A1, complement))]
    # Write the updated bim file
    fwrite(bim[,c("SNP", "B.A1")], "EUR.a1", col.names=F, sep="\t")
    ```


4\. Identify SNPs that have different allele in base and target (usually due to difference in genome build or Indel)

=== "Without data.table"

    ```R
    mismatch <-
        bim$SNP[!(bim$SNP %in% info.match$SNP |
                    bim$SNP %in% info.complement$SNP | 
                    bim$SNP %in% info.recode$SNP |
                    bim$SNP %in% info.crecode$SNP)]
    write.table(
        mismatch,
        "EUR.mismatch",
        quote = F,
        row.names = F,
        col.names = F
    )
    q() # exit R
    ```

=== "With data.table"

    ``` R
    mismatch <- bim[!(SNP %in% info.match |
                        SNP %in% com.snps |
                        SNP %in% recode.snps |
                        SNP %in% com.recode), SNP]
    write.table(mismatch, "EUR.mismatch", quote=F, row.names=F, col.names=F)
    q() # exit R
    ```

We can then use the **EUR.a1** file to update the A1 alleles

## \# Duplicate SNPs
Make sure to remove any duplicate SNPs in your target data (these target data were simulated and so include no duplicated SNPs).

## \# Sex chromosomes 
Sometimes sample mislabelling can occur, which may lead to invalid results. One indication of a mislabelled sample is a difference between reported sex and that indicated by the sex chromosomes. While this may be due to a difference in sex and gender identity, it could also reflect mislabeling of samples or misreporting and, thus, individuals in which there is a mismatch between biological and reported sex are typically removed. A sex check can be performed in PLINK, in which individuals are called as females if their X chromosome homozygosity estimate (F statistic) is < 0.2 and as males if the estimate is > 0.8.

Before performing a sex check, pruning should be performed (see [here](target.md#35-standard-gwas-qc)).
A sex check can then easily be conducted using `plink`
```bash
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.valid.sample \
    --check-sex \
    --out EUR.QC
```

This will generate a file called **EUR.QC.sexcheck** containing the F-statistics for each individual. Individuals are typically called as being biologically male if the F-statistic is > 0.8 and biologically female if F < 0.2.

=== "Without library"

    ```R
    # Read in file
    valid <- read.table("EUR.valid.sample", header=T)
    dat <- read.table("EUR.QC.sexcheck", header=T)
    valid <- subset(dat, STATUS=="OK" & FID %in% valid$FID)
    write.table(valid[,c("FID", "IID")], "EUR.QC.valid", row.names=F, col.names=F, sep="\t", quote=F) 
    q() # exit R
    ```

=== "With data.table"

    ```R
    library(data.table)
    # Read in file
    valid <- fread("EUR.valid.sample")
    dat <- fread("EUR.QC.sexcheck")[FID%in%valid$FID]
    fwrite(dat[STATUS=="OK",c("FID","IID")], "EUR.QC.valid", sep="\t") 
    q() # exit R
    ```

??? note "How many samples were excluded due mismatched Sex information?"
    - `4` samples were excluded

## \# Sample overlap
Since the target data were simulated there are no overlapping samples between the base and target data here (see the relevant section of [the paper](https://www.nature.com/articles/s41596-020-0353-1) for discussion of the importance of avoiding sample overlap). 

## \# Relatedness
Closely related individuals in the target data may lead to overfitted results, limiting the generalisability of the results. 

Before calculating the relatedness, pruning should be performed (see [here](target.md#35-standard-gwas-qc)).
Individuals that have a first or second degree relative in the sample ($\hat{\pi} > 0.125$) can be removed with the following command:

```bash
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.valid \
    --rel-cutoff 0.125 \
    --out EUR.QC
```

??? note "How many related samples were excluded?"
    - `0` samples were excluded

!!! note
    A greedy algorithm is used to remove closely related individuals in a way that optimizes the size of the sample retained.                However, the algorithm is dependent on the random seed used, which can generate different results. Therefore, to reproduce
    the same result, you will need to specify the same random seed. 
    
    PLINK's algorithm for removing related individuals does not account for the phenotype under study. 
    To minimize the removal of cases of a disease, the following algorithm can be used instead: 
    [GreedyRelated](https://gitlab.com/choishingwan/GreedyRelated).

## Generate final QC'ed target data file
After performing the full analysis, you can generate a QC'ed data set with the following command:

```bash
plink \
    --bfile EUR \
    --make-bed \
    --keep EUR.QC.rel.id \
    --out EUR.QC \
    --extract EUR.QC.snplist \
    --exclude EUR.mismatch \
    --a1-allele EUR.a1
```

Each of the parameters corresponds to the following

| Parameter | Value | Description|
|:-:|:-:|:-|
| bfile | EUR | Informs `plink` that the input genotype files should have a prefix of `EUR` |
| keep | EUR.QC.rel.id | Informs `plink` that we only want to keep samples in `EUR.QC.rel.id` |
| extract | EUR.QC.snplist | Informs `plink` that we only want to use SNPs in `EUR.QC.snplist` in the analysis |
| exclude | EUR.mismatch | Informs `plink` that we wish to remove any SNPs in `EUR.mismatch`|
| a1-allele |  EUR.a1 | Fix all A1 alleles to those specified in `EUR.a1` |
| out | EUR.QC | Informs `plink` that all output should have a prefix of `EUR.QC` |

