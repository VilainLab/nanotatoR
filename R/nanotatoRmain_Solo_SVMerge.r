#' Annotation and visualisation of Bionano SV, of Solo SVMerge samples.
#'
#' @param method_entrez  character. Input Method for terms. Choices are 
#'                "Single","Multiple" and "Text".
#' @param decipherpath  character. Decipher database path.
#' @param pattern_Proband  character. Pattern for proband. 
#' @param termPath  character. Path and file name for textfile.
#' @param term  character. Single or Multiple Terms.
#' @param outpath character. Path where gene lists are saved.
#' @param omim character. omim2gene file name and location.
#' @param omimID character. Omim ID.
#' @param RNASeqData  dataFrame. RNAseq data with gene names.
#' @param RNASeqPATH  character. RNAseq dataset path . 
#' @param clinvar character. clinvar file name and location.
#' @param gtr character. gtr file name and location.
#' @param removeClinvar logical. Deletes the Clinvar database if TRUE.
#' @param removeGTR logical. Deletes the GTR database if TRUE.
#' @param downloadClinvar logical. Downloads the Clinvar database if TRUE.
#' @param downloadGTR logical. Downloads the GTR database if TRUE.
#' @param url_gtr character. url for GTR.
#' @param thresh integer. Threshold for the number of terms sent to entrez.
#'                Note if large lists are sent to ncbi, it might fail to get
#'                processed. Default is 5.
#' @param internalBNDB character. internak Bionano merged databse.
#' @param smap  character. Path to SMAP file.
#' @param labelType character. Type of labels used for mapping. 
#'        Choices are Dual, DLE and Both.
#' @param SVMerge_path    character. Path for the Dual labelled cmap
#' @param SVMerge_pattern    character. pattern of the dual files.
#' @param SE_path    character. Path for the Dual labelled cmap
#' @param SE_pattern    character. pattern of the dual files.
#' @param Samplecodes character. File containing relations and IDs 
#' associated to them.
#' @param mergeKey character. File containing sample ID and relation.
#' @param mergedKeyoutpath character. File path storing sample name and nanoID
#'        key information.
#' @param mergedKeyFname character. File name storing sample name and nanoID
#'        key information.
#' @param bed  Text Bionano Bed file.
#' @param EnzymeType Character. Type of enzyme. Options Dual and DLE.
#' @param inputfmtBed character Whether the bed input is UCSC bed or Bionano bed.
#' @param outpath  character Path for the output files.
#' @param n  numeric Number of genes to report which are nearest to the breakpoint.
#' Default is  	3.
#' @param termListPresent  logical Checks whether term list is provided by the user.
#' @param geneListPresent  logical Checks whether gene list is provided by the user.
#' @param primaryGenesPresent  logical Checks whether 
#' primarygene list is provided by the user.
#' @param returnMethod Character. Type of output Dataframe or in Text format.
#' @param mergedFiles  character. Path to the merged SV files.
#' @param smappath  character. path to the query smap file.
#' @param smap  character. File name for the smap
#' @param buildSVInternalDB  boolean. Checking whether the merged solo 
#' file database exist.
#' @param InternaldatabasePresent  boolean. Checking whether
#' internal DB present.
#' @param RNASeqDatasetPresent  boolean. Checking whether
#' RNASeq database present or not.
#' @param buildBNInternalDB  boolean. Checking whether the merged BNDB 
#' file database exist.
#' @param path  character. Path to the solo file database.
#' @param pattern  character. pattern of the file names to merge.
#' @param outpath  character. Path to merged SV solo datasets.
#' @param win_indel_INF  Numeric. Insertion and deletion error window.
#' @param win_inv_trans_INF  Numeric. Inversion and translocation error window.
#' @param perc_similarity_INF  Numeric . ThresholdPercentage similarity 
#' of the query SV and reference SV.
#' @param win_indel_DGV  Numeric. Insertion and deletion error window for DGV.
#' @param win_inv_trans_DGV  Numeric. Inversion and translocation error window
#' for DGV.
#' @param perc_similarity_DGV  Numeric . ThresholdPercentage similarity 
#' of the query SV and reference SV, for DGV..
#' @param indelconf  Numeric. Threshold for insertion and deletion confidence.
#' @param invconf  Numeric. Threshold for inversion confidence.
#' @param transconf  Numeric. Threshold for translocation confidence.
#' @param returnMethod character. Choice between Text and DataFrame.
#' @param hgpath  character. Path to Database of Genomic Variants (DGV)
#'                Text file. 
#' @param RZIPpath  character. Path to RZippath.
#' @param labelType character. Type of labels used for mapping. 
#'        Choices are Dual, DLE and Both.
#' @param SVMerge_path    character. Path for the Dual labelled cmap
#' @param SVMerge_pattern    character. pattern of the dual files.
#' @param SE_path    character. Path for the Dual labelled cmap
#' @param SE_pattern    character. pattern of the dual files.
#' @param Samplecodes character. File containing relations and IDs 
#' associated to them.
#' @param mergeKey character. File containing sample ID and relation.
#' @param mergedKeyoutpath character. File path storing sample name and nanoID
#'        key information.
#' @param smappath  character. Path and file name for textfile.
#' @param outpath character. Path where gene lists are saved.
#' @param indexfile character. indexfile containing nano ID 
#' and sample relation.
#' @param thresh integer. Threshold for the number of terms sent to entrez.
#'                Note if large lists are sent to ncbi, it might fail to get
#'                processed. Default is 5.
#' @param returnMethod character. Choice between text or data frame as the output.
#' @param RNAseqcombo boolean whether RNASeq datasets are combined or not.
#' @param RNASeqDir boolean Directory for RNASeq.
#' @param datGeneListPath Character Path for genelist.
#' @param outpath Character Directory to the output file.
#' @param limsize  Numeric. Minimum size for SV. Default 1000.
#' @param outputType Variants in excel tabs or in different csv files.
#'        Options Excel or csv.
#' @param fileprefix  character Prefix to use for each of the files in the directory.
#' @param directoryName Directory name where individual SV files will be stored.
#' @param outputFilename Character Output filename.
#' @return Excel file containing the annotated SV map, tabs divided based on
#' type of SVs.
#' @return Text files containg gene list and terms associated with them
#'         are stored as text files.
#' @examples
#' smapName="NA12878_Q.S_VAP_SVmerge_solo5.txt"
#' smap = system.file("extdata", smapName, package="nanotatoR")
#' bedFile <- system.file("extdata", "HomoSapienGRCH19_lift37.bed", package="nanotatoR")
#' hgpath=system.file("extdata", "GRCh37_hg19_variants_2016-05-15.txt", package="nanotatoR")
#' decipherpath = system.file("extdata", "population_cnv.txt", package="nanotatoR")
#' omim = system.file("extdata", "mim2gene.txt", package="nanotatoR") 
#' clinvar = system.file("extdata", "localPDB/", package="nanotatoR") 
#' gtr = system.file("extdata", "gtrDatabase.txt", package="nanotatoR")
#' labelType = c("SE")
#' SE_path = system.file("extdata", "SoloFile/", package="nanotatoR")
#' SE_pattern = "*_DLE1_*"
#' Samplecodes = system.file("extdata", "nanotatoR_sample_codes.csv", package="nanotatoR")
#' mergeKey = system.file("extdata", "nanotatoR_control_sample_codes.csv", package="nanotatoR")
#' mergedKeyoutpath = system.file("extdata", package="nanotatoR")
#' mergedKeyFname = "Sample_index.csv" 
#' RNASeqDir = system.file("extdata", "NA12878_P_Blood_S1.genes.results", package="nanotatoR")
#' path = system.file("extdata", "Bionano_config/", package = "nanotatoR")
#' pattern = "_hg19.txt"
#' outputFilename <- "NA12878_Q.S_VAP_SVmerge_solo5_out"
#' outpath <- system.file("extdata", smapName, package = "nanotatoR")
#' RZIPpath <- system.file("extdata", "zip.exe", package = "nanotatoR")
#' nanotatoR_main_Solo_SVmerge(
#' smap = smap, bed = bedFile, inputfmtBed = c("bed"), 
#' n=3,EnzymeType = c("SVMerge"),
#' buildBNInternalDB=TRUE,
#'  path = path , pattern = pattern, 
#' buildSVInternalDB = TRUE,
#' labelType = c("SE"), decipherpath = decipherpath,
#' SE_path = SE_path, SE_pattern = SE_pattern,
#' win_indel_INF = 10000, win_inv_trans_INF = 50000, 
#' perc_similarity_INF= 0.5, indelconf = 0.5, invconf = 0.01, 
#' transconf = 0.1, 
#' hgpath = hgpath, win_indel_DGV = 10000, 
#' win_inv_trans_DGV = 50000, 
#' perc_similarity_DGV = 0.5, limsize = 1000,
#' method_entrez=c("Single"), 
#' term = "Liver cirrhosis", RZIPpath = RZIPpath,
#' omim = omim, clinvar = clinvar, gtr = gtr, 
#' removeClinvar = TRUE, removeGTR = TRUE, 
#' downloadClinvar = FALSE, downloadGTR = FALSE,
#' RNASeqDatasetPresent = TRUE,
#' RNAseqcombo = TRUE, geneListPresent = FALSE,
#' RNASeqDir = RNASeqDir, returnMethod = "dataFrame",
#' pattern_Proband = "*_P_*",
#' outpath = outpath,
#' indexfile = system.file("extdata", "Sample_index.csv",package="nanotatoR"),
#' primaryGenesPresent = FALSE,
#' outputFilename = outputFilename, 
#' termListPresent = FALSE,
#' InternaldatabasePresent = TRUE,
#' outputType = c("Excel"))
#' @importFrom stats na.omit 
#' @export
nanotatoR_main_Solo_SVmerge <-function(
    smap, bed, inputfmtBed = c("bed", "BNBed"), 
    n=3,
    buildBNInternalDB = TRUE,
    mergedFiles , smappath ,  
    buildSVInternalDB = FALSE, 
    path, pattern, 
    win_indel_INF = 10000, win_inv_trans_INF = 50000, 
    perc_similarity_INF= 0.5, indelconf = 0.5, invconf = 0.01, 
    transconf = 0.1,
    hgpath, win_indel_DGV = 10000, win_inv_trans_DGV = 50000, 
    perc_similarity_DGV = 0.5, 
    method_entrez = c("Single","Multiple","Text"), termPath, 
    term, thresh=5, limsize = 1000,
    EnzymeType = c("SVmerge", "SE"),
    labelType = c("SVMerge", "SE", "Both"),
    SVMerge_path ,SVMerge_pattern , 
    SE_path , SE_pattern,
    Samplecodes ,mergeKey,
    mergedKeyoutpath, mergedKeyFname,
    RNAseqcombo=TRUE,RNASeqDir,returnMethod="dataFrame",
    RNASeqData,RNASeqPATH,
    pattern_Proband=NA,
    outpath,outputFilename="", 
    termListPresent = TRUE,
    internalBNDB, clinvar,
    InternaldatabasePresent = TRUE,
    RNASeqDatasetPresent = TRUE,
    geneListPresent = TRUE,
    omim, gtr, 
    removeClinvar = FALSE, removeGTR = FALSE, 
    downloadClinvar = FALSE, downloadGTR = FALSE,
    url_gtr, omimID,
    RZIPpath, directoryName, fileprefix,
    datGeneListPath, decipherpath,
    indexfile, 
    primaryGenesPresent = TRUE, 
    outputType = c("Excel", "csv"))
    {
    print("####PipeLine Starts####")
	 start_time <- Sys.time()
    termListPresent = termListPresent
    if(termListPresent == TRUE){
   
    dat_geneList <- tryCatch(
        gene_list_generation(
            method_entrez = c("Text"), 
            termPath = termPath, 
            thresh = 5, 
            omimID = omimID,
            returnMethod = c("dataFrame"), omim = omim, 
            clinvar = clinvar, 
            gtr = gtr,
            removeGTR = removeGTR, 
            downloadClinvar = downloadClinvar, 
            downloadGTR = downloadGTR,
            removeClinvar = removeClinvar,
            url_gtr = url_gtr
            ),
            error = function(e) {
                    print(paste("gene_list_generation fails"))
                    return (NA)
            }
        )
    }else if (geneListPresent == TRUE)
    { 
        dat_geneList <- read.table(datGeneListPath, header = TRUE)
    }else{ dat_geneList <- NULL } 
    #start_time <- Sys.time()
    end_time <- Sys.time()
    print(paste("Time taken to run gene_list_generation is:" , start_time-end_time))
    start_time <- Sys.time()
    datcompSmap <- tryCatch(overlapnearestgeneSearch(smap=smap,
            bed=bed, 
            EnzymeType = EnzymeType,
            inputfmtBed = inputfmtBed,
            n = 3,
            returnMethod_bedcomp = "dataFrame", 
            input_fmt_SV = "Text" ),
            error = function(e) {
                    stop(paste("nanotatoR Pipeline cannot work"))
                    return (NA)
            })
    end_time <- Sys.time()
    print(paste("Time taken to run compSmapbed is:" , start_time-end_time))
    start_time <- Sys.time()
    datDGV <- tryCatch(
        DGVfrequency(
            hgpath = hgpath, 
            smap_data = datcompSmap, 
            win_indel_DGV = win_indel_DGV, 
            input_fmt_SV = "dataFrame",
            EnzymeType  = EnzymeType,
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
	buildSVInternalDB = buildSVInternalDB
    if(buildSVInternalDB==FALSE){
    datInf <- tryCatch(
        internalFrequency_solo(
            mergedFiles = mergedFiles , 
            buildSVInternalDB = FALSE, 
            smapdata = datDGV, 
            input_fmt = "dataFrame",
            win_indel = win_indel_INF , 
            win_inv_trans = win_inv_trans_INF, 
            smap=smap, 
            input_fmt_SV = "dataFrame",
            EnzymeType  = EnzymeType,
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
        internalFrequency_solo(
            buildSVInternalDB = TRUE, 
            smapdata = datDGV, 
            labelType = labelType,
            SVMerge_path = SVMerge_path, 
            SVMerge_pattern = SVMerge_pattern, 
            SE_path = SE_path, 
            SE_pattern = SE_pattern, 
            input_fmt_SV = "dataFrame",
            EnzymeType  = EnzymeType,
            Samplecodes = Samplecodes, 
            mergeKey = mergeKey,
            outpath = outpath,
            mergedKeyoutpath = mergedKeyoutpath,
            mergedKeyFname = mergedKeyFname,
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
	buildBNInternalDB = buildBNInternalDB
if(buildBNInternalDB==FALSE){
        datchort <- tryCatch(
        BNDBfrequency(
        internalBNDB = internalBNDB, 
        buildBNInternalDB = FALSE, 
        smapdata = datInf, 
        win_indel = win_indel_INF, 
        win_inv_trans = win_inv_trans_INF, 
        input_fmt_SV = "dataFrame",
        EnzymeType  = EnzymeType, 
        perc_similarity = perc_similarity_INF, 
        indelconf = indelconf, 
        invconf = invconf, 
        limsize = limsize,
        transconf = transconf,
        returnMethod = c("dataFrame")),
        error = function(e) {
                    print(paste("BNDBfrequency cannot work"))
                    return (datInf)
            }
        )
} else{
    datchort <- tryCatch(
    BNDBfrequency(
        buildBNInternalDB = TRUE, 
        dbOutput = c("dataframe"),
        smapdata = datInf, 
        BNDBpath = path, 
        BNDBpattern = pattern, 
        win_indel = win_indel_INF, 
        input_fmt_SV = "dataFrame",
        EnzymeType  = EnzymeType,
        win_inv_trans = win_inv_trans_INF, 
        perc_similarity = perc_similarity_INF, 
        indelconf = indelconf, 
        invconf = invconf, 
        limsize=limsize, 
        transconf = transconf,
        returnMethod = c("dataFrame")),
        error = function(e) {
                    print(paste("BNDBfrequency cannot work"))
                    return (datInf)
            }
        )
    
}
    
    end_time <- Sys.time()
     print(paste("Time taken to run BNDBfrequency is:" , start_time-end_time))
     dim(datchort)
start_time <- Sys.time()
datdecipher <- tryCatch(
        Decipherfrequency(
        decipherpath = decipherpath, 
        smap_data = datchort,
        win_indel = win_indel_INF, 
        input_fmt_SV = "dataFrame",
        EnzymeType  = EnzymeType,
        perc_similarity = perc_similarity_INF, 
        returnMethod = c("dataFrame")),
        error = function(e) {
                    print(paste("Decipherfrequency cannot work"))
                    return (datchort)
            }
        )
        
end_time <- Sys.time()
print(paste("Time taken to run Decipherfrequency is:" , start_time-end_time))
dim(datdecipher)

start_time <- Sys.time()
RNASeqDatasetPresent = RNASeqDatasetPresent
RNAseqcombo = RNAseqcombo
if(RNASeqDatasetPresent == TRUE){
if(RNAseqcombo==TRUE){
    RNASeqData <- tryCatch(RNAseqcombine_solo(RNASeqDir = RNASeqDir, 
                returnMethod="dataFrame"),
    error = function(e) {
        print(paste("RNAseqcombo cannot work"))
        return (NA)
        }
    )
    datRNASeq <- tryCatch(
        SVexpression_solo(
        input_fmt_SV = "dataFrame", 
        smapdata = datdecipher, 
        smappath = smappath, 
        input_fmt_RNASeq = "dataFrame", 
        RNASeqData = RNASeqData, 
        outputfmt = "datFrame",
        EnzymeType  = EnzymeType,
        pattern_Proband = pattern_Proband),
        error = function(e) {
            print(paste("RNAseqcombo cannot work"))
            return (datdecipher)
        }
    )
}else{
    datRNASeq <- tryCatch(
        SVexpression_solo(
        input_fmt_SV = "dataFrame", 
        smapdata = datdecipher, 
        smappath = smappath, 
        input_fmt_RNASeq = "Text",
        RNASeqPATH = RNASeqDir,
        EnzymeType  = EnzymeType,
        outputfmt = "datFrame",
        pattern_Proband = pattern_Proband),
        error = function(e) {
        print(paste("RNAseqcombo cannot work"))
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
run_bionano_filter_SVMerge_solo(
input_fmt_geneList = "dataFrame", 
input_fmt_SV = "dataFrame", 
svData = datRNASeq, 
dat_geneList = dat_geneList,
outpath = outpath,
EnzymeType  = EnzymeType,
outputType = outputType,
outputFilename = outputFilename,
RZIPpath = RZIPpath, 
primaryGenesPresent = primaryGenesPresent),
error = function(e) {
        print(paste("run_bionano_filter_SE_solo cannot work"))
        return (datRNASeq)
        }
    )
} else{
tryCatch( 
run_bionano_filter_SVMerge_solo(
input_fmt_geneList = "dataFrame", 
input_fmt_SV = "dataFrame", 
svData = datdecipher, 
dat_geneList = dat_geneList,
outpath = outpath,
EnzymeType  = EnzymeType,
outputFilename = outputFilename,
outputType = outputType,
RZIPpath = RZIPpath, 
primaryGenesPresent = primaryGenesPresent),
error = function(e) {
        print(paste("run_bionano_filter_SE_solo cannot work"))
        return (datdecipher)
        }
    )}
end_time <- Sys.time()
print(paste("Time taken to run run_bionano_filter_SE_solo is:" , start_time-end_time))

}
