# GAMBIT

A C++ tool for Gene-based Analysis with oMniBus, Integrative Tests

- Implements SKAT, burden, and ACAT gene-based test methods using variant- or region-based functional annotations
- Calculates annotation-stratified gene-based tests (e.g., TWAS/PrediXcan tests using eSNPs, gene-based tests using only coding variants, and gene-based tests using enhancer-to-target-gene maps)
- Calculates omnibus gene-based tests by aggregating across annotation classes
- Inputs: GWAS association summary statistics file (chromosome, position, ref/alt allele, and z-score or beta-hat + se), annotation files, and LD reference panel


## GWAS Summary Statistics

- GWAS summary statistics files can be specified via `--gwas my_summary_stats.txt.gz`. Input files must be ordered by chromosome and genomic position, with input fields as shown below:

```
#CHR  POS     REF  ALT  SNP_ID      N         ZSCORE   ANNO
1     721290  C    G    rs12565286  58663.62  0.86661  Intergenic
1     752566  G    A    rs3094315   57135     0.5521   Intergenic
1     775659  A    G    rs2905035   54570     1.12098  Intron:LOC643837
1     777122  A    T    rs2980319   54570     1.11906  Exon:LOC643837
```

- The first four fields and `ZSCORE` are required, while `SNP_ID`, `ANNO` and `N` (effective sample size) are optional. 
- See `format_gwas_summary_stats.sh` for annotating GWAS summary statistics files using `EPACTS/TabAnno`.

## Annotation-Stratified Gene-Based Tests

#### Gene-Based Analysis with Regulatory Elements
- To compute gene-based tests using regulatory element annotations, specify an annotation bed file with regulatory-element-to-target-gene weights via `--anno-bed my_reg_elems.txt.gz`, formatted

```
#CHR  START   END     CLASS     ELEMENT_ID          TARGET_GENES                     ANNO
chr1  567400  567600  Enhancer  chr1:567400:567600  MIB2:4.12|CPTP:2.53|GLTPD1:2.53  .
chr1  568000  568200  Enhancer  chr1:568000:568200  ATAD3A:2.75                      .
chr1  758600  758800  Enhancer  chr1:758600:758800  C1orf170:2.57|PERM1:2.57         .
chr1  769200  769400  Enhancer  chr1:769200:769400  C1orf170:3.36|PERM1:3.36         .
```

- Association tests for individual regulatory elements is reported in `*.stratified_out.txt` files, and gene-based p-values (aggregating across regulatory elements for each gene) in `*.summary_out.txt` files.

- **Aggregation Methods for Regulatory Elements.** By default, GAMBIT aggregates test statistics across variants in regulatory elements using a weighted sum of single-variant chi-squared statistics (SKAT gene-based test).  To instead use weighted ACAT to combine single-variant p-values, specify `--acat`.

#### Gene-Based Analysis with Coding and Other Annotated Variants 

- To compute gene-based tests using coding and other variants, GAMBIT relies on the `ANNO` field in GWAS summary statistics and an annotation hierarchy definitions file specified via `--anno-defs my_defs.txt`, formatted as below:

```
#CLASS    SUBCLASS          ANNO_TERMS
Coding    Protein_Altering  Nonsynonymous,Start_Loss,Stop_Gain,Stop_Loss,CodonGain,CodonLoss,Frameshift
Coding    Splice_Site       Essential_Splice_Site,Normal_Splice_Site
Coding    Exon_Other        Exon,Synonymous
UTR       UTR3              Utr3
UTR       UTR5              Utr5
```

- The `ANNO_TERMS` field specifies a comma-separated list of annotation terms (matching terms from the GWAS summary statistics file's `ANNO` field), and `CLASS` and `SUBCLASS` determine the annotation hierarchy and classes reported in output files. 

- **Gene-Based Test Output.** Test statistics stratified by gene and annotation subclass are provided in `*.stratified_out.txt` files, and gene-based p-values (aggregating across annotation classes for each gene) in `*.summary_out.txt` files.

- **Variant Aggregation Methods.** By default, GAMBIT aggregates test statistics across variants using a weighted sum of single-variant chi-squared statistics (SKAT gene-based test).  To instead use weighted ACAT to combine single-variant p-values, specify `--acat`.

#### TWAS Analysis
- To compute TWAS/PrediXcan gene-based tests using GAMBIT, specify an eWeight file via `--eweights my_eWeights.txt.gz`, formatted

```
##TISSUE_IDS=0:Adipose_Subcutaneous,1:Adipose_Visceral_Omentum,2:Adrenal_Gland,3:Artery_Aorta
#CHR  POS     RSID       REF  ALT  BETAS
1     752566  rs3094315  G    A    C1orf159=3.92e-02@0|UBE2J2=-1.49e-02@0|FAM87B=2.75e-01@1;1.25e-01@2;1.17e-01@3
1     752721  rs3131972  A    G    LINC00115=1.15e-01@0;1.75e-02@3;4.90e-02@4|RP11-206L10.8=3.21e-02@1
1     754182  rs3131969  A    G    LINC00115=-2.1e-02@1|RP5-857K21.2=-8.27e-02@2|RP11-206L10.9=-1.11e-01@2
1     760912  rs1048488  C    T    C1orf159=3.35e-04@0|TTLL10=-1.4e-02@3|FAM87B=1.75e-01@1;1.12e-01@2;9.51e-02@3|SAMD11=-1.27e-02@2
```

- The `BETAS` field format is `eGene_A=Weight_A1@Tissue_A1;Weight_A2@Tissue_A2|eGene_B=Weight_B1@Tissue_B1`, and labels for tissue IDs can be specified in the header. 
- **Subsetting tissues.** To restrict analysis to a subset of tissues/cell-types, specify a comma-separated list of tissues following the `--tissues` flag. By default, GAMBIT includes all tissues/cell-types present in the eWeight file. 
- **Tissue Aggregation for Omnibus tests.** GAMBIT reports both single-tissue TWAS/PrediXcan analysis results, and omnibus tests results aggregating across all specified tissues/cell-types for each eGene. Omnibus p-values for multi-tissue TWAS/PrediXcan analysis can be calculated in GAMBIT using either 1) the maximum single-tissue test statistic based on the joint distribution of single-tissue statistics, 2) the sum of squared single-tissue z-scores (analogous to SKAT), or 3) ACAT [default]. Omnibus test method for multi-tissue analysis can be specified via `--tissue-aggreg` (`ACAT`, `MinP`, `SKAT`, or `ALL`).
- **Single-tissue and omnibus test output.** Gene-based tests and p-values for each eGene-tissue pair are reported in `*.stratified_out.txt` files, and omnibus p-values (aggregating across all tissues for each eGene) in `*.summary_out.txt` files.

#### dTSS-Weighted Gene-Based Tests
- To incorporate un-annotated regulatory variants in gene-based analysis, GAMBIT implements a dTSS (distance to Transcription Start Site) weighted gene-based test, which aggregates all single-variant p-values within a specified window from each gene's TSS using ACAT and assigns higher weight to variants nearer the TSS using an exponential decay function. 
- To compute dTSS-weighted gene-based tests, specify a TSS bed file via `--tss-bed my_tss_bed.bed.gz`, fomatted 

```
#CHR  START   END     SYMBOL      GENE             GENE_ANNO
1     11868   11869   DDX11L1     ENSG00000223972  transcribed_unprocessed_pseudogene
1     62947   62948   OR4G11P     ENSG00000240361  transcribed_unprocessed_pseudogene
1     69090   69091   OR4F5       ENSG00000186092  protein_coding
1     131024  131025  CICP27      ENSG00000233750  processed_pseudogene
```

- **Window size.** The window size for dTSS-weighted gene-based tests can be modified by specifying `--tss-window BASEPAIRS` (500 Kbp by default).
- **dTSS decay function.** The relative weight assigned to variants nearer/farther from the TSS can be modified by specifying `--tss-alpha ALPHA`, where alpha=0 implies all variants receive equal weight, and larger values confer more weight to variants nearer the TSS. `--tss-alpha` also accepts comma-separated lists of alpha values, in which case GAMBIT computes global test p-values across all specified values (individual p-values are reported in `INFO` output field). By default, GAMBIT uses dTSS alpha values `1e-4,5e-5,1e-5,5e-6`. 

## Methods References

Statistical methods implemented in GAMBIT:
- Sequence Kernel Association Test (SKAT): [Wu et al. (2011), *AJHG*](https://doi.org/10.1016/j.ajhg.2011.05.029)
- TWAS: [Gusev et al. (2016), *Nat Genet*](https://www.nature.com/articles/ng.3506)
- PrediXcan: [Gamazon et al. (2015), *Nat Genet*](https://www.nature.com/articles/ng.3367) and [Barbeira et al. (2018), *Nat Comm*](https://www.nature.com/articles/s41467-018-03621-1)
- Aggregated Cauchy Association Test (ACAT): [Liu et al. (2019), *AJHG*](https://doi.org/10.1016/j.ajhg.2019.01.002) and [Liu and Xie (2018), *arXiv*](https://arxiv.org/abs/1808.09011)
- Asymptotically exact Harmonic Mean P-value (HMP): [Wilson (2019), *PNAS*](https://doi.org/10.1073/pnas.1814092116)

## Software References

Libraries and resources used or adapted in GAMBIT:<br/>
** PDF, CDF, and quantile functions**<br/>
- [ROOT system library, CERN](https://root.cern.ch/)
- [CDFLIB, Brown et al.](https://people.sc.fsu.edu/~jburkardt/cpp_src/cdflib/cdflib.html)
- [eigenmvn, Benazera et al.](https://github.com/beniz/eigenmvn)
- [libMvtnorm, Zhan et al.](https://github.com/zhanxw/libMvtnorm)
** Hopscotch hashing**<br/>
- [tsl library, Tessil et al.](https://github.com/Tessil/hopscotch-map)
** Tabix and HTSLIB**<br/>
- [tabixpp, ekg et al.](https://github.com/ekg/tabixpp)
- [htslib, samtools team](https://github.com/samtools/htslib)
** Matrix libraries**<br/>
- [Eigen matrix library](http://eigen.tuxfamily.org)

## Feedback and bug reports

- Feel free to contact Corbin Quick (corbinq@gmail.com) with questions, bug reports, or feedback

