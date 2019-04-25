
#nanotatoR: structural variant annotation and classification
#author: "Surajit Bhattacharya,Hayk Barsheghyan, Emmanuele C Delot and Eric Vilain"

# Introduction
Short-read sequencing (**SRS**) is the predominant technique of DNA sequencing used for clinical diagnosis. It utilizes flowcells covered with millions of surface-bound oligonucleotides that allow parallel sequencing of hundreds of millions of independent short reads. However, as the average sequencing read length is approximately 150 bp, large structural variants (**SVs**) and copy number variants (**CNVs**) are challenging to observe. This creates a diagnostic gap between the clinical phenotypes and the underlying genetic mechanisms in the field of biomedical sciences. Novel technologies such as optical genome mapping and long-read sequencing have partially addressed the issues of SV and CNV detection; however, the identification of pathogenic variants among thousands of called SVs/CNVs throughout the genome has proven to be challenging as the analytical pipelines available for single nucleotide variants are not applicable to SV analysis. Thus, we have built an R-based annotation package “nanotatoR” specifically for structural variants to provide a multitude of critical functional annotations for large genomic variations.

#Package Installation

**nanotatoR** is currently available from the GitHub repository. Installation method is as follows:
```{r eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("nanotatoR", version = "3.8")
```

```{r eval=TRUE}
library("nanotatoR")
```

The **nanotatoR** package is compatible with R versions ≥ 3.6.

#Package Functionalities

Given a list of structural variants with chromosome number, SV type, and variant start/end positions, nanotatoR can perform the following functions (both GRCh37/38 can be used):

##Structural Variant Frequency

Genomic variant frequencies are one of the most important filtration characteristics for identification of rare, possibly pathogenic, variants. NanotatoR uses the Database of Genomic Variants (**DGV**) and DECIPHER (**DECIPHER**), publicly available reference control database, to estimate structural variant frequencies in the general population. Compared with the traditional single nucleotide variant frequency calculations, frequency estimates for structural variants pose larger difficulty due to the much higher breakpoint variability between “same” structural variants. In order to provide accurate estimates of population frequencies, nanotatoR recognizes five categories of SVs: insertions, deletions, duplications, inversions and translocations (e.g. “gains” of DNA material are considered as “insertions/duplications” and “losses” of DNA material as “deletions”). In order for two SVs to be considered as “same”, nanotatoR, by default, checks whether they belong to the same category (e.g. deletion), have greater than 50% size similarity and the SV breakpoints (start and end positions) are within 10 kbp from each other. Currently the 50% size similarity cutoff is not applied for duplications, inversions and translocations as the sizes for these structural variants are not computed by most SV callers; however, a breakpoint cutoff of 50 kbp is applied in order to identify matching variants (note: duplication breakpoint cutoff is 10 kbp). Default breakpoint cutoffs and percent similarity determined by Bionano Genomics (see the references). 
In order for the natotatoR to estimate SV frequencies it requires the following input files: DECIPHER reference file(**decipherpath**) and “smap” file containing structural variant information generated via default Bionano Genomics Solve/Tools Pipelines for optical mapping-based SV calling (**smappath**). The default input parameters for SV breakpoints and percent similarities are as follows: insertions/duplications/deletions (**win_indel**) of 10,000 bases, inversions and translocations (**win_inv_trans**) of 50,000 bases and percentage similarity (**perc_similarity**) of 0.5 or 50%. The output from this function can be of 2 types: an R object (**dataFrame**) or a text file (**Text**).

```{r eval=TRUE}
decipherpath <- system.file("extdata", "population_cnv.txt", package = "nanotatoR")
smapName <- "F1.1_TestSample1_solo_hg19.smap"
smappath <- system.file("extdata", package = "nanotatoR")
win_indel=10000;win_inv_trans=50000;perc_similarity=0.5
decipherext<-Decipher_extraction(decipherpath, input_fmt = "Text", smappath, 
smap= smapName,
win_indel = 10000, perc_similarity = 0.5, returnMethod="dataFrame")
```
Other than **DGV** and **DECIPHER**, nanotatoR can also take as input a data set of variants created by Bionano Genomics comprising of 204 samples. Input parameters are similar to that of **DECIPHER**, with an addition of following parameters: confidence score threshold for insertion and deletion (**indelconf**; Default is 0.5), inversion (**invconf**; Default is 0.01) and translocation (**transconf**; Default is 0.1 (determined by Bionano Genomics).
Note: The Bionano genomic reference database can be downloaded using the following command wget http://bnxinstall.com/solve/Solve3.3_10252018.tar.gz , followed by decompressing it using the tar -xvzf Solve3.3_10252018.tar.gz. The folder containing the database is in the config file which can be found in the folowing directory $PWD/Solve3.3_10252018/VariantAnnotation/10252018/config/. The reference files are named  based on the variant type and reference genome. For example: current_ctrl_dup_hg19_anonymize.txt, would be the duplication reference variant file for hg19 reference genome.

```{r eval=TRUE}
mergedFiles <- system.file("extdata", "BNSOLO2_merged.txt", package = "nanotatoR")
smapName <- "F1.1_TestSample1_solo_hg19.smap"
smappath <- system.file("extdata",  package = "nanotatoR")
win_indel = 10000; win_inv_trans = 50000; perc_similarity = 0.5;
indelconf = 0.5; invconf = 0.01;transconf = 0.1
cohortFreq<-cohortFrequency(internalBNDB = mergedFiles , smappath , 
smap=smapName, input_fmt ="Text", buildBNInternalDB=FALSE, win_indel, win_inv_trans, 
perc_similarity , indelconf, invconf, transconf,returnMethod="dataFrame")
```

##Cohort Analysis

The cohort analysis is designed to provide internal variant frequency, and parental zygosity within the cohort. The function will first merge all of the available individual smaps to form an “internal reference” (**buildSVInternalDB**=*TRUE*) or use an existing file if the internal SV database is already available (**buildSVInternalDB**=*FALSE*). The function requires the paths of the query smap file, the merged internal database file, as well as the following parameters: confidence score threshold for insertion and deletion (**indelconf**; Default is 0.5), inversion (**invconf**; Default is 0.01) and translocation (**transconf**; Default is 0.1 (determined by Bionano Genomics). The output can be a (**dataFrame**) or a text file (**Text**).

```{r eval=TRUE}
mergedFiles <- system.file("extdata", "nanotatoR_merged.txt", package = "nanotatoR")
smapName <- "F1.1_TestSample1_solo_hg19.smap"
smappath <- system.file("extdata",  package = "nanotatoR")
intFreq <- internalFrequency(mergedFiles = mergedFiles, smappath = smappath ,
smap = smapName, buildSVInternalDB=FALSE, win_indel=10000, 
win_inv_trans=50000, 
perc_similarity=0.5, indelconf=0.5, invconf=0.01, 
transconf=0.1, limsize=1000, win_indel_parents=5000,input_fmt="Text",
win_inv_trans_parents=40000,
returnMethod="dataFrame")
intFreq [1,]
```

##Gene Overlap 

Incorporates known gene and non-coding RNA genomic locations that overlap with or are near the identified structural variants. NanotatoR automatically determines the number of overlapping genes for a given structural variant, provides gene strandedness (+ or -) as well as the percent overlap of SVs with the gene. The function also provides gene names and corresponding distances from SVs for the nearest genes (user selected, default 3 on each side) that are upstream and downstream.
The function requires an input BED file (**inputfmt**=*BED*) or Bionano compliant BED file (**inputfmt**=*BNBED*), which re-codes the X and Y chromosome notations to 23 and 24, respectively. The BED files are used to overlap structural variants of the query smap file (**smap**) with overlapping and non-overlapping upstream and downstream genes (**n**;default is 3). The output from this function can be of 2 types: an R object (**dataFrame**) or a text file (**Text**).



```{r eval=TRUE}
smapName="F1.1_TestSample1_solo_hg19.smap"
smap = system.file("extdata", smapName, package="nanotatoR")
bedFile <- system.file("extdata", "Homo_sapiens.Hg19_BN_bed.txt", package="nanotatoR")
outpath <- system.file("extdata",  package="nanotatoR")
datcomp <- compSmapbed(smap, bed=bedFile, inputfmtBed =  "BNBED", n = 3, returnMethod_bedcomp = c("dataFrame"))
datcomp[1,]
```

##Entrez Extract

Generates a list of genes involved with the patient phenotype and overlaps it with gene names that span structural variants. The user provided, phenotypic key words are used to generate gene lists from selected databases such as ClinVar, OMIM, GTR and Gene Registry. The generated lists are used to prioritize structural variants that occur in genes known to be associated with patient’s phenotype.
The input to the function is a term, which can be provided as a single character input (**method**="Single"), or a vector of terms (**method**="Multiple") or as a text file (**method**="Text"). The output can be selected as a (**dataFrame**) or a text file (**Text**).

```{r eval=TRUE}
 terms="Muscle Weakness"
 gene <- gene_list_generation(method="Single", term=terms, returnMethod="dataFrame")
 gene[1:10,]
```


##Filter Variants 

Provides easy to use, user-selected, filtration criteria to segregate variants into corresponding groups (such as de novo, inherited or occurring in the gene list). In this step the generated or available gene lists are appended to the smap file.
The input file for this function can either be an (**smap**) or a dataframe. Both raw and nanotatoR-annotated smaps can serve as inputs. It has an option to take the input of the smap (**input_fmt_svMap**) and genelist (**input_fmt_geneList**) as a dataFrame or a text file, and produces an excel as the output.


```{r eval=FALSE}
terms <- "Muscle Weakness"
gene <- gene_list_generation(
method = "Single", term = terms,
returnMethod_GeneList = "dataFrame"
)
RzipFile = "zip.exe"
RZIPpath <- system.file("extdata", RzipFile, package = "nanotatoR")
smapName <- "F1.1_TestSample1_solo_hg19.smap"
smappath <- system.file("extdata", smapName, package = "nanotatoR")
smappath1 <- system.file("extdata", package = "nanotatoR")
run_bionano_filter(input_fmt_geneList = c("dataframe"), 
input_fmt_svMap = c("Text"),
SVFile = smappath, dat_geneList = gene, 
outpath = smappath1, outputFilename = "test", 
RZIPpath = RZIPpath)
```

##Main

The main function consecutively runs the available nanotatoR sub-functions by appending the outputs from each function. 
It takes as an input the smap file, DGV file, BED file, internal database file, phenotype term list as an input. It also takes in the output path and the filename for the final excel file. Individual, sub-function, input parameters are also available for user selections.


```{r eval=FALSE}
terms <- "Muscle Weakness"
gene <- gene_list_generation(
method = "Single", term = terms,
returnMethod = "dataFrame"
)
mergedFiles <- system.file ("extdata", "BNSOLO2_merged.txt", 
package = "nanotatoR")
RzipFile = "zip.exe"
RZIPpath <- system.file("extdata", RzipFile, package = "nanotatoR")
smapName <- "F1.1_TestSample1_solo_hg19.smap"
smappath <- system.file("extdata", smapName, package = "nanotatoR")
path <- system.file("extdata", "SoloFile/", package = "nanotatoR")
hgpath <- system.file ("extdata", "GRCh37_hg19_variants_2016-05-15.txt", package = "nanotatoR")
decipherpath <- system.file("extdata", "population_cnv.txt", package = "nanotatoR")
bedFile <- system.file("extdata", "Homo_sapiens.Hg19.bed", package="nanotatoR")
pattern <- "_hg19.smap"
nanotatoR_main(smap = smappath, bed = bedFile,
inputfmtBed = c("BNBED"),
n = 3,  buildSVInternalDB = TRUE, soloPath = path, solopattern = pattern,
input_fmt_INF = c("dataFrame"), buildBNInternalDB = FALSE,
returnMethod_bedcomp = c("dataFrame"), returnMethod_DGV = c("dataFrame"),
returnMethod_Internal = c("dataFrame"), input_fmt_DGV = c("dataFrame"),
hgpath = hgpath, smapName = smapName, limsize=1000, win_indel_parents=5000,
decipherpath = decipherpath, dbOutput_Int = "dataframe",
win_inv_trans_parents=40000, win_indel_DGV = 10000,
input_fmt_geneList = c("dataFrame"), input_fmt_svMap = c("dataFrame"),
input_fmt_decipher = "dataFrame",input_fmt_BN = "dataFrame",
returnMethod_GeneList = c("dataFrame"),returnMethod_BNCohort =  c("dataFrame"),
returnMethod_decipher = c("dataFrame"), mergedFiles_BN = mergedFiles,
dat_geneList = gene , method_entrez = "", 
outpath = smappath, outputFilename = "test",
RZIPpath = RZIPpath)
```

#References

1. Hayk Barseghyan, Wilson Tang, Richard T. Wang, Miguel Almalvez, Eva Segura, Matthew S. Bramble, Allen Lipson, Emilie D. Douine, Hane Lee, Emmanuèle C. Délot, Stanley F. Nelson and Eric Vilain.Next-generation mapping: a novel approach for detection of pathogenic structural variants with a potential utility in clinical diagnosis. Genome Medicine 2017 9:90. https://doi.org/10.1186/s13073-017-0479-0

2. Winter, D. J.  rentrez: an R package for the NCBI eUtils API The R Journal 2017 9(2):520-526

3. Christopher Brown.  hash: Full feature implementation of hash/associated arrays/dictionaries.2013. R package version 2.2.6. https://CRAN.R-project.org/package=hash

4. Hadley Wickham. stringr: Simple, Consistent Wrappers for Common String Operations. 2018. R package version 1.3.1. https://CRAN.R-project.org/package=stringr

5. Bionano Genomics. Theory Of Operation – Structural Variant Calling. https://bionanogenomics.com/wp-content/uploads/2018/04/30110-Bionano-Solve-Theory-of-Operation-Structural-Variant-Calling.pdf

6. Bionano Genomics. Theory of Operation - Variant Annotation Pipeline. https://bionanogenomics.com/wp-content/uploads/2018/04/30190-Bionano-Solve-Theory-of-Operation-Variant-Annotation-Pipeline.pdf
