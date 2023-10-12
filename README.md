# SEM_PGS

### This is a repository containing the workflow for multivariate SEM_PG

This README file provides information re background and versions of the pipeline. 

The current version of directory will host all scripts of the pipeline being used for quality control, partitioning of transmitted and non transmitted nucleotides and calculation of polygenic scores for simulated and then empirical data.
Briefly, the pipeline takes care of *insert here steps, in shape of a table?bullet points?* etc. between Discovery and Target samples, and calculation of final scores. 



This directory is temporary and being improved, and a more generic version of pipeline / container for customizable genetic data and non-cloud based environment will be added.

no.valenza@icloud.com

Read this [tutorial](https://choishingwan.github.io/PRS-Tutorial/)

## Workflow

### 0. Gwas on Discovery sample

The first run will be done with simulated data to test the pipeline. TBD if we bypass the GWAS step and use sumstats directly.

Empirical Discovery sample will be UKB data. Will start with height as a trait to proof-of-concept the pipeline.

#### 0.1 QC of Discovery sample already-ran GWAS data

##### 0.11 Heritability check

PGS analyses are performed on base data with a chip-heritability estimate h2 snp > 0.05. The chip-heritability of a GWAS can be estimated using e.g. LD Score Regression (LDSC). 

#### 0.12 Effect allele

It is important to know which allele is the effect allele and which is the non-effect allele for PGS association results to be in the correct direction.

/!\ Some GWAS results files do not make clear which allele is the effect allele and which is the non-effect allele. If the incorrect assumption is made in computing the PRS, then the effect of the PRS in the target data will be in the wrong direction.

To avoid misleading conclusions the effect allele from the base (GWAS) data must be known.

##### 0.13 Check file transfer

A common problem is that the downloaded base data file can be corrupted during download, which can cause PRS software to crash or to produce errors in results. However, a md5sum hash is generally included in files so that file integrity can be checked. The following command performs this md5sum check:


Linux
md5sum Height.gwas.txt.gz

OS X
if the file is intact, then md5sum generates a string of characters, which in this case should be: a2b15fb6a2bbbe7ef49f67959b43b160. If a different string is generated, then the file is corrupted.

### 1. GWAS on Target sample

The first run will be done with simulated data to test the pipeline. TBD if we bypass the GWAS step and use sumstats directly.

Empirical Target sample will be Moba data. Also starting with height as a proof-of-concept.

### 2. 
