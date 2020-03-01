source("X:/Suro/Annotator/R_02152019/R_BioconductoR_release_oct_2019/InternalFrequency_DLE_Solo.r")
source("X:/Suro/Annotator/R_02152019/R_BioconductoR_release_oct_2019/Bed_SV_Comp.r")
source("X:/Suro/Annotator/R_02152019/R_BioconductoR_release_oct_2019/DGV_extraction_WIP.r")
source("X:/Suro/Annotator/R_02152019/R_BioconductoR_release_oct_2019/entrez_extract.r")
#source("Y:/Suro/Annotator/R_01182019/entrez_extract.r")
source("X:/Suro/Annotator/R_02152019/R_BioconductoR_release_oct_2019/filterread_DLE_Solo.r")
source("X:/Suro/Annotator/R_02152019/R_BioconductoR_release_oct_2019/RNASEQ_Analysis_solo.r")
source("X:/Suro/Annotator/R_02152019/R_BioconductoR_release_oct_2019/DecipherExtraction.r")
source("X:/Suro/Annotator/R_02152019/R_BioconductoR_release_oct_2019/cohortDB_comp_WIP.r")
source("X:/Suro/Annotator/R_02152019/R_BioconductoR_release_oct_2019/mergesmap.r")
mergedFiles="X:/Hayks_Materials/BNG/Projects/Control/Merged/nano_merged_09042019.txt"
smap="X:/Hayks_Materials/BNG/Projects/Control/NA12878_Sample/NA12878_DLE1_VAP_solo4.smap"
bed="X:/Suro/Annotator/Data/Homo_sapiens_GRCH19_BN.bed"
hgpath="X:/Suro/Annotator/Data/GRCh37_hg19_variants_2016-05-15.txt"
termPath = "X:/Hayks_Materials/BNG/Projects/Control/NA12878_Terms.csv"
#termPath="X:/Hayks_Materials/BNG/Projects/DSD/DSD_Terms/DSD_Draft_terms_cleaned_GeneList1.txt"
outpath="X:/Hayks_Materials/BNG/Projects/VSVM_nanotatoR_EXCELs/"
outputFilename="NA12878_DLEEnzyme_NR_07222019"
RNASeqDir="X:/Hayks_Materials/UDN_RNA_Seq_Counts/Counts/NA12878/"
decipherpath="X:/Suro/Annotator/Data/population_cnv.txt"
win_indel_INF = 10000; win_inv_trans_INF = 50000; 
perc_similarity_INF= 0.5; indelconf = 0.5; invconf = 0.01; 
transconf = 0.1;win_indel_DGV = 10000; win_inv_trans_DGV = 50000; 
perc_similarity_DGV = 0.5;limsize=1000
perc_similarity_DGV = 0.5;pattern_Proband="*_P_*"
pattern_Sibling= "*_AS_*"
buildBNInternalDB=TRUE;buildSVInternalDB=FALSE;RNAseqcombo=TRUE
path="X:/Suro/bionano/BNG_tools/variant_annotation/10252018/config"
pattern="_hg19_anonymize.txt"
RZIPpath="X:/Suro/Rtools/bin/zip.exe"
indexfile = "X:/Hayks_Materials/BNG/Projects/Sample_NIDKeys_nanotatoR.csv"
clinvar = "X:/Suro/nanotatoRDatabases/clinvar_10292019.txt"
gtr = "X:/Suro/nanotatoRDatabases/gtr_10292019.txt"
omim = "X:/Suro/nanotatoRDatabases/mim2gene.txt"

library(openxlsx)
library(hash)
#library(biomaRt)
library(rentrez)
library(stringr)
library(EnsDb.Hsapiens.v79)
library(org.Hs.eg.db)
library(tidyverse)
library(httr)
library(AnnotationDbi)
library(tidyverse)

termListPresent = TRUE; RNASeqDatasetPresent = TRUE

print("####PipeLine Starts####")
	if(termListPresent == TRUE){
	start_time <- Sys.time()
	
	dat_geneList <- tryCatch(
	    gene_list_generation(
	        method_entrez = c("Text"), 
            termPath = termPath, 
            thresh = 5, 
            returnMethod = c("dataFrame"), omim = omim, 
			clinvar = clinvar, 
			gtr = gtr,
            removeClinvar = FALSE, removeGTR = FALSE, 
			downloadClinvar = TRUE, downloadGTR = TRUE,
			url_clinvar = "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
			url_gtr = "ftp://ftp.ncbi.nlm.nih.gov/pub/GTR/data/test_condition_gene.txt"
			),
			error = function(e) {
				    print(paste("gene_list_generation fails"))
					return (NA)
			}
			)
	} else{ dat_geneList <- NULL}
	#start_time <- Sys.time()
	end_time <- Sys.time()
    print(paste("Time taken to run gene_list_generation is:" , start_time-end_time))
	start_time <- Sys.time()
	datcompSmap <- tryCatch(compSmapbed(smap=smap,
	        bed=bed, 
			EnzymeType = "DLE",
			inputfmt = "BNBED",
			n = 3,
			returnMethod = "dataFrame", 
			input_fmt_smap = "Text" ),
			error = function(e) {
				    stop(paste("nanotatoR Pipeline cannot work"))
					return (NA)
			})
	end_time <- Sys.time()
    print(paste("Time taken to run compSmapbed is:" , start_time-end_time))
	start_time <- Sys.time()
	datDGV <- tryCatch(
	    DGV_extraction(
		    hgpath = hgpath, 
			smap_data = datcompSmap, 
			input_fmt_DGV = c("DataFrame"), 
            win_indel_DGV = win_indel_DGV, 
			win_inv_trans_DGV = win_inv_trans_DGV, 
	        perc_similarity_DGV = perc_similarity_DGV,
			returnMethod = c("dataFrame")
			),
			error = function(e) {
				    print(paste("DGV_extraction cannot work"))
					return (datcompSmap)
			}
			)
	end_time <- Sys.time()
    print(paste("Time taken to run DGV_extraction is:" , start_time-end_time))
	dim(datDGV)
	start_time <- Sys.time()
	if(buildSVInternalDB==FALSE){
	datInf <- tryCatch(
	    internalFrequency_DLE_Solo(
	    mergedFiles = mergedFiles , 
		buildSVInternalDB = FALSE, 
		smapdata = datDGV, 
		input_fmt = "dataFrame",
		win_indel = win_indel_INF , 
		win_inv_trans = win_inv_trans_INF, 
		smap=smap, 
        perc_similarity = perc_similarity_INF, 
		indelconf = indelconf , 
		invconf = invconf, 
		transconf = transconf,
		limsize = limsize, 
		returnMethod="dataFrame", 
		indexfile = indexfile),
		error = function(e) {
				    print(paste("internalFrequency cannot work"))
					return (datDGV)
			}
			)
	} else{
    datInf <- tryCatch(
	        internalFrequency_DLE_Solo(
	        buildSVInternalDB = TRUE, 
			path = path, 
			pattern = pattern, 
			outpath = outpath, 
			smap=smap,
            smapdata = datDGV, 
			smapName = smapName, 
			input_fmt = "dataFrame",
			limsize = limsize, 
			win_indel = win_indel_INF, 
            win_inv_trans = win_inv_trans_INF, 
			perc_similarity = perc_similarity_INF, 
			indelconf = indelconf, 
			invconf = invconf, 
			transconf = transconf, 
            returnMethod=c("dataFrame"), 
			indexfile = indexfile),
			error = function(e) {
				    print(paste("internalFrequency cannot work"))
					return (datDGV)
			}
			)
    
}
     end_time <- Sys.time()
     print(paste("Time taken to run internalFrequency is:" , start_time-end_time))
	 dim(datInf)
	start_time <- Sys.time()
if(buildBNInternalDB==FALSE){
        datchort <- tryCatch(
		cohortFrequency(
		internalBNDB = internalBNDB, 
		buildBNInternalDB = FALSE, 
		smapdata = datInf, 
		input_fmt = c("dataFrame"), 
		win_indel = win_indel_INF, 
		win_inv_trans = win_inv_trans_INF, 
		perc_similarity = perc_similarity_INF, 
		indelconf = indelconf, 
		invconf = invconf, 
		limsize = limsize,
		transconf = transconf,
		returnMethod = c("dataFrame")),
		error = function(e) {
				    print(paste("cohortFrequency cannot work"))
					return (datInf)
			}
		)
} else{
    datchort <- tryCatch(
	cohortFrequency(
	    buildBNInternalDB = TRUE, 
		dbOutput = c("dataframe"),
		smapdata = datInf, 
		input_fmt = c("dataFrame"), 
		InternalDBpath = path, 
		InternalDBpattern = pattern, 
		win_indel = win_indel_INF, 
		win_inv_trans = win_inv_trans_INF, 
		perc_similarity = perc_similarity_INF, 
		indelconf = indelconf, 
		invconf = invconf, 
		limsize=limsize, 
		transconf = transconf,
		returnMethod = c("dataFrame")),
		error = function(e) {
				    print(paste("cohortFrequency cannot work"))
					return (datInf)
			}
		)
    
}
    
	end_time <- Sys.time()
     print(paste("Time taken to run cohortFrequency is:" , start_time-end_time))
	 dim(datchort)
start_time <- Sys.time()
datdecipher <- tryCatch(
        Decipher_extraction(
		decipherpath = decipherpath, 
		smap_data = datchort,
		input_fmt = c("dataFrame"),
        win_indel = win_indel_INF, 
		perc_similarity = perc_similarity_INF, 
		returnMethod = c("dataFrame")),
		error = function(e) {
				    print(paste("Decipher_extraction cannot work"))
					return (datchort)
			}
		)
		
end_time <- Sys.time()
print(paste("Time taken to run Decipher_extraction is:" , start_time-end_time))
dim(datdecipher)

start_time <- Sys.time()
if(RNASeqDatasetPresent == TRUE){
if(RNAseqcombo==TRUE){
    RNASeqData <- tryCatch(RNAseqcombine_solo(RNASeqDir = RNASeqDir, 
	            returnMethod="dataFrame"),
	error = function(e) {
	    print(paste("RNAseqcombine cannot work"))
		return (NA)
		}
	)
    datRNASeq <- tryCatch(
        SmapRNAseqquery_solo(
		input_fmt_SV = "dataFrame", 
		smapdata = datdecipher, 
		smappath = smappath, 
		input_fmt_RNASeq = "dataFrame", 
		RNASeqData = RNASeqData, 
		outputfmt = "datFrame",
		EnzymeType = "DLE",
		pattern_Proband = pattern_Proband),
	    error = function(e) {
	        print(paste("RNAseqcombine cannot work"))
		    return (datdecipher)
		}
	)
}else{
    datRNASeq <- tryCatch(
        SmapRNAseqquery_solo(
		input_fmt_SV = "dataFrame", 
		smapdata = datdecipher, 
		smappath = smappath, 
		input_fmt_RNASeq = "Text",
		RNASeqPATH = RNASeqDir,
		EnzymeType = "DLE",
		outputfmt = "datFrame",
		pattern_Proband = pattern_Proband),
	error = function(e) {
	    print(paste("RNAseqcombine cannot work"))
		return (datdecipher)
		}
	)
}
end_time <- Sys.time()
print(paste("Time taken to run SmapRNAseqquery is:" , start_time-end_time))
start_time <- Sys.time()
'r3 <-read.csv("C:/Annotator/Data/RNASeq_F1.1.csv")
datRNASeq <- cbind(datdecipher, r3[, 2:ncol(r3)])
datRNASeq[is.na(datRNASeq)]<-"-"
datRNASeq[is.na(datRNASeq$Type2)]<-"-\"
#dat_geneList<-read.table("C:/Annotator/Data/F3.1_UDN992683_GeneList.txt",header=TRUE,sep=" ")'
tryCatch( 
run_bionano_filter_DLE_solo(
input_fmt_geneList = "dataFrame", 
input_fmt_svMap = "dataFrame", 
SVFile = NULL,
svData = datRNASeq, 
dat_geneList = dat_geneList,
outpath = outpath,
outputFilename = outputFilename,
RZIPpath = RZIPpath, primaryGenesPresent = TRUE),
error = function(e) {
	    print(paste("run_bionano_filter cannot work"))
		return (datRNASeq)
		}
	)
} else{
tryCatch( 
run_bionano_filter_DLE_solo(
input_fmt_geneList = "dataFrame", 
input_fmt_svMap = "dataFrame", 
SVFile = NULL,
svData = datdecipher, 
dat_geneList = dat_geneList,
outpath = outpath,
outputFilename = outputFilename,
RZIPpath = RZIPpath, primaryGenesPresent = TRUE),
error = function(e) {
	    print(paste("run_bionano_filter cannot work"))
		return (datdecipher)
		}
	)}
end_time <- Sys.time()
print(paste("Time taken to run run_bionano_filter is:" , start_time-end_time))

