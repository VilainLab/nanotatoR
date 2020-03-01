#' Annotation and visualisation of Bionano SV, of SVmerge Dyad samples.
#'
#' @param method  character. Input Method for terms. Choices are 
#'                "Single","Multiple" and "Text".
#' @param termPath  character. Path and file name for textfile.
#' @param terms  character. Single or Multiple Terms.
#' @param outpath character. Path where gene lists are saved.
#' @param omim character. omim2gene file name and location.
#' @param clinvar character. clinvar file name and location.
#' @param gtr character. gtr file name and location.
#' @param removeClinvar logical. Deletes the Clinvar database if TRUE.
#' @param removeGTR logical. Deletes the GTR database if TRUE.
#' @param downloadClinvar logical. Downloads the Clinvar database if TRUE.
#' @param downloadGTR logical. Downloads the GTR database if TRUE.
#' @param url_gtr character. url for GTR.
#' @param url_clinvar character. url for clinvar.
#' @param thresh integer. Threshold for the number of terms sent to entrez.
#'                Note if large lists are sent to ncbi, it might fail to get
#'                processed. Default is 5.
#' @param smap  character. Path to SMAP file.
#' @param bed  Text Bionano Bed file.
#' @param EnzymeType Character. Type of enzyme. Options Dual and DLE.
#' @param inputfmt character Whether the bed input is UCSC bed or Bionano bed.
#' @param outpath  character Path for the output files.
#' @param n  numeric Number of genes to report which are nearest to the breakpoint.
#' Default is  	3.
#' @param termListPresent  logical Checks whether term list is provided by the user.
#' @param returnMethod Character. Type of output Dataframe or in Text format.
#' @param mergedFiles  character. Path to the merged SV files.
#' @param smappath  character. path to the query smap file.
#' @param smap  character. File name for the smap
#' @param buildSVInternalDB  boolean. Checking whether the merged solo 
#' file database exist.
#' @param inputfmt character. Choice between Text and DataFrame.
#' @param path  character. Path to the solo file database.
#' @param pattern  character. pattern of the file names to merge.
#' @param outpath  character. Path to merged SV solo datasets.
#' @param win_indel  Numeric. Insertion and deletion error window.
#' @param win_inv_trans  Numeric. Inversion and translocation error window.
#' @param perc_similarity  Numeric . ThresholdPercentage similarity 
#' of the query SV and reference SV.
#' @param url_gtr character. url for GTR.
#' @param url_clinvar character. url for clinvar.
#' @param indelconf  Numeric. Threshold for insertion and deletion confidence.
#' @param indelconf  Numeric. Threshold for inversion confidence.
#' @param indelconf  Numeric. Threshold for translocation confidence.
#' @param returnMethod character. Choice between Text and DataFrame.
#' @param hgpath  character. Path to Database of Genomic Variants (DGV)
#'                Text file. 
#' @param smappath  character. Path and file name for textfile.
#' @param terms  character. Single or Multiple Terms.
#' @param outpath character. Path where gene lists are saved.
#' @param input_fmt character. Choice between text or data frame as 
#' an input to the DGV frequency calculator.
#' @param outputType Variants in excel tabs or in different csv files.
#'        Options Excel or csv.
#' @param fileprefix Prefix to use for each of the files in the directory.
#' @param directoryName Directory name where individual SV files will be stored.
#' @param smap_data Dataset containing smap data.
#' @param thresh integer. Threshold for the number of terms sent to entrez.
#'                Note if large lists are sent to ncbi, it might fail to get
#'                processed. Default is 5.
#' @param returnMethod character. Choice between text or data frame as the output.
#' @param input_fmt_geneList character. Choice of gene list input 
#'        Text or Dataframe.
#' @param input_fmt_svMap  character. Choice of gene list input 
#'        Text or Dataframe.
#' @param SVFile  character. SV file name.
#' @param svData Dataframe Input data containing SV data.
#' @param dat_geneList Dataframe Input data containing geneList data.
#' @param fileName Character Name of file containing Gene List data.
#' @param outpath Character Directory to the output file.
#' @param outFileName Character Output filename.
#' @return Excel file containing the annotated SV map, tabs divided based on
#' type of SVs.
#' @return Text files containg gene list and terms associated with them
#'         are stored as text files.
#' @examples
#' terms="Muscle Weakness"
#' gene_list_generation(method="Single", terms, outpath=systemDir())
#' @importFrom stats na.omit 
#' @export
nanotatoR_SVmerge_Dyad<-function(
	smap, bed, inputfmt = c("bed", "BNBed"), 
    n=3,EnzymeType = c("Dual", "DLE"),
    mergedFiles , smappath ,  
	buildSVInternalDB=FALSE, path, pattern, 
	win_indel_INF = 10000, win_inv_trans_INF = 50000, 
	perc_similarity_INF= 0.5, indelconf = 0.5, invconf = 0.01, 
	transconf = 0.1,
	hgpath, win_indel_DGV = 10000, win_inv_trans_DGV = 50000, 
    perc_similarity_DGV = 0.5,
	method_entrez=c("Single","Multiple","Text"), smapName,termPath, 
    term, thresh=5, omim, clinvar, gtr, 
	removeClinvar = FALSE, removeGTR = FALSE, 
	downloadClinvar = FALSE, downloadGTR = FALSE,
	url_gtr = "ftp://ftp.ncbi.nlm.nih.gov/pub/GTR/data/test_condition_gene.txt",
	RNAseqcombine = TRUE,RNASeqDir,returnMethod = "dataFrame",
	input_fmt_SV = "dataFrame",
    input_fmt_RNASeq = "dataFrame",
    RNASeqData, RNASeqPATH, outputfmt = "datFrame",
    pattern_Proband = NA,pattern_Mother = NA,pattern_Father = NA,
	pattern_Sibling = NA,
	outpath,outputFilename="", 
	termListPresent = TRUE,
	InternaldatabasePresent = TRUE,
	RNASeqDatasetPresent = TRUE,
	outputType = c("Excel", "csv"), 
	directoryName, fileprefix)
	{
print("####PipeLine Starts####")
	if(termListPresent == TRUE){
	start_time <- Sys.time()
	
	dat_geneList <- tryCatch(
	    gene_list_generation(
	        method_entrez = c("Single"), 
            term = term, 
            thresh = 5, 
            returnMethod = c("dataFrame"), omim = omim, 
			gtr = gtr,
            removeGTR = FALSE, 
			downloadClinvar = TRUE, downloadGTR = FALSE,
			),
			error = function(e) {
				    print(paste("gene_list_generation fails"))
					return (NA)
			}
		)
	}else if (datgeneList_Present == TRUE)
	{ 
	    dat_geneList <- read.table(datGeneListPath, header = TRUE)
	}else{ dat_geneList <- NULL }
	start_time <- Sys.time()
	end_time <- Sys.time()
    print(paste("Time taken to run gene_list_generation is:" , start_time-end_time))
	start_time <- Sys.time()
	datcompSmap <- tryCatch(compSmapbed(smap=smap,
	        EnzymeType = "Dual",
	        bed=bed, 
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
	datInf <- tryCatch(internalFrequencyTrio_Dyad(
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
	        internalFrequencyTrio_Dyad(
	        buildSVInternalDB = TRUE, 
			path = path, 
			pattern = pattern, 
			outpath = outpath, 
			smap=smap,
            smapdata = datDGV, 
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
run_bionano_filter_Dual_duo(
input_fmt_geneList = "dataFrame", 
input_fmt_svMap = "dataFrame", 
SVFile = NULL,
svData = datRNASeq, 
dat_geneList = dat_geneList,
outpath = outpath,
outputFilename = outputFilename,
RZIPpath = RZIPpath,
primaryGenesPresent = TRUE,
outputType = c("Excel")),
error = function(e) {
	    print(paste("run_bionano_filter cannot work"))
		return (datRNASeq)
		}
	)
} else{
tryCatch( 
run_bionano_filter_Dual_duo(
input_fmt_geneList = "dataFrame", 
input_fmt_svMap = "dataFrame", 
SVFile = NULL,
svData = datdecipher, 
dat_geneList = dat_geneList,
outpath = outpath,
outputFilename = outputFilename,
RZIPpath = RZIPpath,
primaryGenesPresent = TRUE,
outputType = c("Excel")),
error = function(e) {
	    print(paste("run_bionano_filter cannot work"))
		return (datdecipher)
		}
	)}
end_time <- Sys.time()
print(paste("Time taken to run run_bionano_filter is:" , start_time-end_time))
}
