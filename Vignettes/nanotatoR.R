## ----eval=FALSE----------------------------------------------------------
#  library("devtools")
#  devtools::install_github("VilainLab/Nanotator")

## ----eval=TRUE-----------------------------------------------------------
library("nanotatoR")

## ----eval=TRUE-----------------------------------------------------------
smap="F1.1_GM24385_DLE-1_P_trio_hg19.smap"
smappath = system.file("extdata", smap, package="nanotatoR")
hgpath=system.file("extdata", "GRCh37_hg19_variants_2016-05-15.txt", package="nanotatoR")
win_indel_DGV=10000;win_inv_trans_DGV=50000;perc_similarity_DGV=0.5
datDGV<- DGV_extraction (hgpath, smappath,input_fmt_DGV = "Text", smap, 
win_indel_DGV = 10000, win_inv_trans_DGV = 50000,perc_similarity_DGV = 0.5,returnMethod_DGV="dataFrame")


## ----eval=FALSE----------------------------------------------------------
#  
#  path <- system.file("extdata", "SoloFile", package="nanotatoR")
#  pattern="_hg19.smap"
#  smapName="F1.1_GM24385_DLE-1_P_trio_hg19.smap"
#  smappath = system.file("extdata", smapName, package="nanotatoR")
#  win_indel = 10000
#  perc_similarity = 0.5
#  indelconf = 0.5
#  invconf = 0.01
#  transconf =0.1
#  limsize=1000
#  intFreq<-internalFrequency(smappath=smappath , buildSVInternalDB=TRUE, soloPath=path, solopattern=pattern,outpath=path,input_fmt_INF="Text",win_indel,limsize =limsize,
#  win_inv_trans=50000, perc_similarity ,indelconf, invconf ,transconf,returnMethod_Internal="dataFrame")
#  intFreq[1:2,]

## ----eval=TRUE-----------------------------------------------------------

smapName="F1.1_GM24385_DLE-1_P_trio_hg19.smap"
smap = system.file("extdata", smapName, package="nanotatoR")
bedFile <- system.file("extdata", "Homo_sapiens.Hg19.bed", package="nanotatoR")
outpath <- system.file("extdata",  package="nanotatoR")
datcomp<-compSmapbed(smap, bed=bedFile, inputfmtBed =  "BED", outpath, 
n = 3, returnMethod_bedcomp = c("Text", "dataFrame"))

## ----eval=TRUE-----------------------------------------------------------
terms="Muscle Weakness"
gene<-gene_list_generation(method="Single", term=terms, 
   	returnMethod_GeneList="dataFrame")
gene[1:10,]

## ----eval=FALSE----------------------------------------------------------
#  smapName="F1.1_GM24385_DLE-1_P_trio_hg19.smap"
#  smappath = system.file("extdata", smapName, package="nanotatoR")
#  terms="Muscle Weakness"
#  gene<-gene_list_generation(method="Single", term=terms,
#  returnMethod_GeneList="dataFrame")
#   run_bionano_filter(SVFile=smappath,fileName,input_fmt_geneList="dataFrame",
#      input_fmt_svMap="Text",RtoolsZIPpath="")

## ----eval=FALSE----------------------------------------------------------
#  terms="Muscle Weakness"
#  gene<-gene_list_generation(method="Single", term=terms,
#    	returnMethod_GeneList="dataFrame")
#  
#  smapName="F1.1_GM24385_DLE-1_P_trio_hg19.smap"
#  smappath = system.file("extdata", smapName, package="nanotatoR")
#  nanotatoR_main(smap, bed, inputfmtBed = c("BNBED"),
#   n=3,mergedFiles , buildSVInternalDB=TRUE, soloPath, solopattern,
#  input_fmt_INF=c("dataframe"),returnMethod_GeneList=c("dataframe"),
#  returnMethod_bedcomp=c("dataframe"),returnMethod_DGV=c("dataframe"),
#  returnMethod_Internal=c("dataframe"),input_fmt_DGV=c("dataframe"),	
#  hgpath, smapName,method=c("Single"), term, thresh=5,
#  input_fmt_geneList=c("dataframe"),input_fmt_svMap=c("dataframe"),
#  svData,dat_geneList,outpath="",outputFilename="",RZIPpath="")

