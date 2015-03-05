# HLAassoc: Tests for association between disease and HLA alleles.

## News

* v1.0 (23 Jan 2015): initial release. 
* v1.1 (16 Feb 2015): linear and logistic regression were added

## Introduction

## Requirement
* [Python 2.7](https://www.python.org/)
* [pandas](http://pandas.pydata.org/)
* [SciPy](http://www.scipy.org/)
* [StatsModels](http://statsmodels.sourceforge.net/)

## Installation:

- **Install Python 2.7**     
Python may already installed on your system. Else, you can follow this [guild](https://wiki.python.org/moin/BeginnersGuide/Download) to install it.

- **Install Python modules**   

```
sudo pip install pandas
sudo pip install git+http://github.com/scipy/scipy/
sudo pip install statsmodels
```

**Note:** There are several free scientific python distributions such as [Anaconda](http://continuum.io/downloads) and [Enthought Canopy](https://www.enthought.com/products/canopy/) which are already integrated the core scientific analytic and scientific Python packages such as `SciPy`, `pandas` and `StatsModels`.

- **Download HLAassoc**   
HLAassoc is available [here](https://github.com/felixfan/HLAassoc/archive/v1.1.tar.gz).

## Optionals and Usage

```
python HLAassoc.py -h
usage: HLAassoc.py [-h] [-v] -i FILE [-d {2,4}] [-m {allelic,dom,rec}]
                   [-t {chisq,fisher,logistic,linear}] [-c COVAR]
                   [-n COVARNAME] [-f FREQ] [-a {FDR,Bonferroni,Holm}]
                   [-o OUT] [-V {False,True}] [-p PERM]

HLA Association Analysis

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -i FILE, --file FILE  input file
  -d {2,4}, --digits {2,4}
                        digits to test, default 4
  -m {allelic,dom,rec}, --model {allelic,dom,rec}
                        genetic model, default allelic
  -t {chisq,fisher,logistic,linear}, --test {chisq,fisher,logistic,linear}
                        statistical test method, default chisq
  -c COVAR, --covar COVAR
                        covariants file
  -n COVARNAME, --covarname COVARNAME
                        select a particular subset of covariates
  -f FREQ, --freq FREQ  minimal frequency, default 0.05
  -a {FDR,Bonferroni,Holm}, --adjust {FDR,Bonferroni,Holm}
                        p value correction, default FDR
  -o OUT, --out OUT     output file
  -V {False,True}, --print {False,True}
                        print output to screen
  -p PERM, --perm PERM  number of permutation
```

### 1) Genotype Input (-i or --file)  

The input file is a white-space (space or tab) delimited file. The first two columns are mandatory: Individual ID and Phenotype. The Individual IDs are alphanumeric and should uniquely identify a person. The second column is phenotype which can be either a quantitative trait or an affection status. Affection status should be coded as 1 and 2 for unaffected and affected, respectively.

Genotypes (column 3 onwards) should also be white-space delimited. Every gene must have two alleles specified. All alleles (see [Nomenclature of HLA Alleles](http://hla.alleles.org/nomenclature/naming.html)) do not need to have the same digits. However, if you want to test association at 4 digits, all alleles should have at least 4 digits resolution. Missing genotype is denoted as `NA`.

No header row should be given. For example, here are two individuals typed for 6 genes (one row = one person):  

```
0001 2 A*02:07:01 A*11:01:01 B*51:01:01 B*51:01:01 C*14:02:01 C*14:02:01 DQA1*01:04:01 DQA1*01:04:01 DQB1*03:03:02 DQB1*05:02:01 DRB1*07:01:01 DRB1*14:54:01
0002 1 A*24:02:01 A*33:03:01 B*15:25:01 B*58:01:01 C*03:02:02 C*04:03 NA NA DQB1*03:01:01 DQB1*03:01:01 DRB1*03:01:01 DRB1*12:02:01
```

**Note:** There are one case and one control. The six genes are: HLA-A, HLA-B, HLA-C, HLA-DQA1, HLA-DQB1 and HLA-DRB1. Each gene has two columns.

**Note:** Individual `0002` does not have genotype for HLA-DQA1 (two `NA`). All alleles have six digits resolution except that one allele of HLA-C of individual `0002` only has four digits resolution. It is fine if we only want to test association at two or four digits resolution.  

### 2) Digits resolution (-d or --digits)   

Test of association using two digits or four digits. When two was used, alleles such as `A*02:01` and `A*02:06` will be combined as `A*02`. Default value is 4.

### 3) Genetic model to test (-m or --model)   


## Output



