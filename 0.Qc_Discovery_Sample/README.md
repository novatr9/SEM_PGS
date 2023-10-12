# 0. Gwas on Discovery sample

The first run will be done with simulated data to test the pipeline. TBD if we bypass the GWAS step and use sumstats directly.

Empirical Discovery sample will be UKB data. Will start with height as a trait to proof-of-concept the pipeline.

# 0.1 QC of Discovery sample already-ran GWAS data

The following abbreviations that can be found in column headers of a sum stat file:

- **CHR**: The chromosome in which the SNP resides
- **BP**: Chromosomal co-ordinate of the SNP
- **SNP**: SNP ID, usually in the form of rs-ID
- **A1**: The effect allele of the SNP
- **A2**: The non-effect allele of the SNP
- **N**: Number of samples used to obtain the effect size estimate
- **SE**: The standard error (SE) of the effect size esimate
- **P**: The P-value of association between the SNP genotypes and the base phenotype
- **OR**: The effect size estimate of the SNP, if the outcome is binary/case-control. If the outcome is continuous or treated as continuous then this will usually be BETA
- **INFO**: The imputation information score
- **MAF**: The minor allele frequency (MAF) of the SNP

# 0.11 Heritability check

PGS analyses are performed on base data with a chip-heritability estimate $h_{snp}^{2} > 0.05$. The chip-heritability of a GWAS can be estimated using e.g. LD Score Regression (LDSC). 

# 0.12 Effect allele

It is important to know which allele is the effect allele and which is the non-effect allele for PGS association results to be in the correct direction.

!!! Important

Some GWAS results files do not make clear which allele is the effect allele and which is the non-effect allele. If the incorrect assumption is made in computing the PRS, then the effect of the PRS in the target data will be in the wrong direction.

To avoid misleading conclusions the effect allele from the base (GWAS) data must be known.

# 0.13 Check file transfer

A common problem is that the downloaded base data file can be corrupted during download, which can cause PRS software to crash or to produce errors in results. 

However, a `md5sum` hash is generally included in files so that file integrity can be checked. The following command performs this ```md5sum``` check:

 ```bash
 md5sum trait.gwas.txt.gz
 ```

if the file is intact, then ```md5sum``` generates a string of characters. If a different string is generated than the one included in the file, then it is corrupted.

### 0.14 Genome build

You must check that your Discovery and Target samples are on the same genome build. If they are not then use a tool such as [LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver)  to make the builds consistent across the data sets.

# 0.15 Standard GWAS QC

Both the Discovery and Target data should be subjected to the standard stringent QC steps performed in GWAS. If the Discovery data have been obtained as summary statistics from a public source, then the typical QC steps we are able to perform are to filter the SNPs according to INFO score and MAF. SNPs with low minor allele frequency (MAF) or imputation information score (INFO) are more likely to generate false positive results due to their lower statistical power (and higher probability of genotyping errors in the case of low MAF). Therefore, SNPs with low MAF and INFO are typically removed before performing downstream analyses. Best practice to remove SNPs with MAF < 1% and INFO < 0.8 (with very large base sample sizes these thresholds could be reduced if sensitivity checks indicate reliable results). These SNP filters can be achieved using the following code:


``` bash, 
gunzip -c trait.gwas.txt.gz |\
awk 'NR==1 || ($11 > 0.01) && ($10 > 0.8) {print}' |\
gzip  > trait.gz
```

The bash code above does the following:

1. Decompresses and reads the **trait.gwas.txt.gz** file
2. Prints the header line (```NR==1```)
3. Prints any line with MAF above 0.01 (```$11``` because the eleventh column of the file contains the MAF information)
4. Prints any line with INFO above 0.8 (```$10``` because the tenth column of the file contains the INFO information)
5. Compresses and writes the results to **trait.gz**

# 0.16 Mismatching SNPs

SNPs that have mismatching alleles reported in the Disc and Target samples are either resolvable by "strand-flipping" the alleles to their complementary alleles in e.g. the Target sample, such as for a SNP with A/C in the Disc data and G/T in the Target, or non-resolvable, such as for a SNP with C/G in the Disc and C/T in the Target. Most polygenic score software perform strand-flipping automatically for SNPs that are resolvable, and remove non-resolvable mismatching SNPs.

--> Talk to **Emmanuel**

Since we need the Target sample to know which SNPs have mismatching alleles, we will perform this strand-flipping in the Target sample.

# 0.17 Duplicate SNPs

If an error has occurred in the generation of the Disc sample, then there may be duplicated SNPs in the Disc sample file. Most PGS software do not allow duplicated SNPs in the Disc sample input and thus they should be removed, using a command such as the one below:

```bash
gunzip -c trait.gz |\
awk '{seen[$3]++; if(seen[$3]==1){ print}}' |\
gzip - > trait.nodup.gz
```

The above command does the following:

1. Decompresses and reads the **trait.gz** file
2. Count number of time SNP ID was observed, assuming the third column contian the SNP ID (`seen[$3]++`). If this is the first time seeing this SNP ID, print it. 
3. Compresses and writes the results to **trait.nodup.gz**


# 0.18  Ambiguous SNPs

If the Disc and Target data were generated using different genotyping chips and the chromosome strand (+/-) that was used for either is unknown, then it is not possible to pair-up the alleles of ambiguous SNPs (i.e. those with complementary alleles, either C/G or A/T SNPs) across the data sets, because it will be unknown whether the Disc and Target data are referring to the same allele or not. While allele frequencies could be used to infer which alleles are on the same strand, the accuracy of this could be low for SNPs with MAF close to 50% or when the base and target data are from different populations. Therefore, we recommend removing all ambiguous SNPs to avoid introducing this potential source of systematic error.

Ambiguous SNPs can be removed in the Disc data and then there will be no such SNPs in the subsequent analyses, since analyses are performed only on SNPs that overlap between the Disc and Target data.

Nonambiguous SNPs can be retained using the following:

```bash
gunzip -c trait.nodup.gz |\
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' |\
    gzip > trait.QC.gz
```

# 0.19 Sex chromosomes

Sometimes sample mislabelling can occur, which may lead to invalid results. One indication of a mislabelled sample is a difference between reported sex and that indicated by the sex chromosomes. While this may be due to a difference in sex and gender identity, it could also reflect mislabeling of samples or misreporting and, thus, individuals in which there is a mismatch between biological and reported sex are typically removed. See the Target Data section in which a sex-check is performed.

# 0.20 Sample overlap

In simulated data there are no overlapping samples between the Disc and Target data. 

### 1. GWAS on Target sample

The first run will be done with simulated data to test the pipeline. TBD if we bypass the GWAS step and use sumstats directly.

Empirical Target sample will be Moba data. Also starting with height as a proof-of-concept.

### 2. 

*TO BE DETERMINED HOW TO LOOK FOR OVERLAP*

# 0.21 Relatedness

Closely related individuals within and between the Disc and the Target data may lead to overfitted results, limiting the generalizability of the results. Relatedness within the Target data is tested in the Target Data section.

The **trait.QC.gz** base data are now ready for using in downstream analyses.
