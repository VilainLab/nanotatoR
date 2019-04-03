#' Annotation of Bionano SV.
#'
#' @param smap  character. Path to SMAP file.
#' @param inputfmtBed character. Choice between Text and DataFrame
#' as input for bed file.
#' @param bed  Text Choice between UCSC bed or Bionano bed.
#' @param n  numeric Number of genes to report which are nearest to 
#' the breakpoint.
#' Default is   3.
#' @param mergedFiles_BN  character. Path to the merged BN SV files.
#' @param mergedFiles_INF  character. Path to the merged BN SV files.
#' @param buildSVInternalDB  boolean. Checking whether the merged solo
#' file database exist or you need to build it. Default= TRUE.
#' @param buildBNInternalDB  boolean. Checking whether the merged Bionano
#' file database exist or you need to build it. Default= TRUE.
#' @param soloPath  character. Path to the solo file database.
#' @param solopattern  character. pattern of the file names to merge.
#' @param InternalDBpath  character. Path to the BNFile file database.
#' @param InternalDBpattern  character. pattern of the BNFile names to merge.
#' @param smapName character. Name of the smap file.
#' @param win_indel_INF  Numeric. Insertion and deletion error window.
#' @param win_inv_trans_INF  Numeric. Inversion and translocation error window.
#' @param perc_similarity_INF  Numeric . ThresholdPercentage similarity
#' of the query SV and reference SV.
#' @param input_fmt_DGV character. Choice between Text and DataFrame
#' for input to DGV.
#' @param input_fmt_INF character. Choice between Text and DataFrame
#' for input to INF.
#' @param dbOutput_BN  character. Output of merged bionano data.
#' @param fname_BN  character. Filename in case dbOutput_BN = Text.
#' @param dbOutput_Int  character. Output of solo bionano data.
#' @param fname_Int  character. Filename in case dbOutput_Int = Text.
#' @param returnMethod_GeneList character. Return Methods from the 
#' gene_list_generation
#' function, choice between Text and Dataframe.
#' @param returnMethod_bedcomp character. Return Methods from the compSmapbed
#' function, choice between Text and Dataframe.
#' @param returnMethod_DGV character. Return Methods from 
#' the DGV_extraction function, choice between
#' Text and Dataframe.
#' @param returnMethod_Internal character. Return Methods 
#' from the internalFrequency function, choice between
#' Text and Dataframe.
#' @param returnMethod_BNCohort character. Return 
#' Methods from the Bionano function, choice between
#' Text and Dataframe.
#' @param returnMethod_decipher character.
#' Return Methods from the decipher Frequency function, choice between
#' Text and Dataframe.
#' @param limsize Numeric. Minimum size of SV that can be determined 
#' accurately by the Bionano SV caller. Default 1000.
#' @param win_indel_parents  Numeric. Insertion and deletion error window to 
#' determine zygosity in case of parents. Default 5000.
#' @param win_inv_trans_parents  Numeric. Inversion and translocation error 
#' window to determine zygosity in case of parents. Default 40000.
#' @param indelconf Numeric. Threshold for insertion and deletion confidence.
#' @param invconf  Numeric. Threshold for inversion confidence.
#' @param transconf  Numeric. Threshold for translocation confidence.
#' @param hgpath  character. Path to Database of Genomic Variants (DGV)
#'                Text file.
#' @param decipherpath  character. Path to DECIPHER. Text file.
#' @param win_indel_DGV  Numeric. Insertion and deletion error window.
#' @param win_inv_trans_DGV  Numeric. Inversion and translocation error window.
#' @param perc_similarity_DGV  Numeric . ThresholdPercentage similarity
#' of the query SV and reference SV.
#' @param method_entrez  character. Input Method for terms for entrez. Choices 
#'                are "Single","Multiple" and "Text".
#' @param termPath  character. Path and file name for textfile for terms.
#' @param term  character. Single or Multiple Terms as vectord.
#' @param thresh integer. Threshold for the number of terms sent to entrez.
#'                Note if large lists are sent to ncbi, it might fail to get
#'                processed. Default is 5.
#' @param input_fmt_geneList character. Choice of gene list input
#'        Text or Dataframe.
#' @param input_fmt_BN character. Choice of Bionano dataset input
#'        Text or Dataframe.
#' @param input_fmt_decipher character. Choice of gene list input
#'        Text or Dataframe.
#' @param input_fmt_svMap  character. Choice of SVmap input for final step
#'        Text or Dataframe.
#' @param dat_geneList Dataframe Input data containing geneList data.
#' @param outpath Character Directory to the output file.
#' @param outputFilename Character Output filename for the annotated smap.
#' @param RZIPpath Character. Path for the Rtools zip.exe
#' @return Excel file containing the annotated SV map, tabs divided based on
#' type of SVs.
#' @return Text files containg gene list 
#' and terms associated with them are stored as text files.
#' \dontrun{
#' @examples
#' terms <- "Muscle Weakness"
#' gene <- gene_list_generation(
#'   method = "Single", term = terms,
#'   returnMethod = "dataFrame"
#' )
#' mergedFiles <- system.file ("extdata", "BNSOLO2_merged.txt", 
#' package = "nanotatoR")
#' RzipFile = "zip.exe"
#' RZIPpath <- system.file("extdata", RzipFile, package = "nanotatoR")
#' smapName <- "F1.1_TestSample1_solo_hg19.smap"
#' smappath <- system.file("extdata", smapName, package = "nanotatoR")
#' path <- system.file("extdata", "SoloFile/", package = "nanotatoR")
#' hgpath <- system.file ("extdata", 
#' "GRCh37_hg19_variants_2016-05-15.txt", package = "nanotatoR")
#' decipherpath <- system.file("extdata", "population_cnv.txt", package = 
#' "nanotatoR")
#' bedFile <- system.file("extdata", "Homo_sapiens.Hg19.bed", 
#' package="nanotatoR")
#' pattern <- "_hg19.smap"
#' nanotatoR_main(smap = smappath, bed = bedFile,
#'  inputfmtBed = c("BNBED"),
#'  n = 3,  buildSVInternalDB = TRUE, soloPath = path, solopattern = pattern,
#'  input_fmt_INF = c("dataFrame"), buildBNInternalDB = FALSE,
#'  returnMethod_bedcomp = c("dataFrame"), returnMethod_DGV = c("dataFrame"),
#'  returnMethod_Internal = c("dataFrame"), input_fmt_DGV = c("dataFrame"),
#'  hgpath = hgpath, smapName = smapName, limsize=1000, win_indel_parents=5000,
#'  decipherpath = decipherpath, dbOutput_Int = "dataframe",
#'  win_inv_trans_parents=40000, win_indel_DGV = 10000,
#'  input_fmt_geneList = c("dataFrame"), input_fmt_svMap = c("dataFrame"),
#'  input_fmt_decipher = "dataFrame",input_fmt_BN = "dataFrame",
#'  returnMethod_GeneList = c("dataFrame"),returnMethod_BNCohort = 
#'  c("dataFrame"),
#'  returnMethod_decipher = c("dataFrame"), mergedFiles_BN = mergedFiles,
#'  dat_geneList = gene , method_entrez = "", 
#'  outpath = smappath, outputFilename = "test",
#'  RZIPpath = RZIPpath
#'  )
#' @importFrom stats na.omit 
#' @export
nanotatoR_main<-function(
    smap, bed, inputfmtBed = c("BED", "BNBED"), 
    n=3, InternalDBpath, InternalDBpattern, dbOutput_Int,
    fname_Int, dbOutput_BN, fname_BN,
    buildSVInternalDB=FALSE, soloPath, solopattern, 
    returnMethod_bedcomp = c("Text","dataFrame"), mergedFiles_BN,
    win_indel_INF = 10000, win_inv_trans_INF = 50000, 
    perc_similarity_INF = 0.5, indelconf = 0.5, invconf = 0.01, 
    transconf = 0.1,returnMethod_DGV = c("Text","dataFrame"),
    hgpath, win_indel_DGV = 10000, win_inv_trans_DGV = 50000, 
    perc_similarity_DGV = 0.5,returnMethod_Internal = c("Text","dataFrame"),
    input_fmt_DGV = c("Text","dataFrame"), input_fmt_BN = c("Text","dataFrame"),
    input_fmt_INF = c("Text","dataFrame"), 
    input_fmt_decipher = c("Text","dataFrame"),
    input_fmt_svMap = c("Text","dataFrame"), dat_geneList,
    decipherpath, input_fmt_geneList = c("Text","dataFrame"),
    returnMethod_GeneList = c("Text","dataFrame"), buildBNInternalDB = FALSE,
    returnMethod_BNCohort = c("Text","dataFrame"), mergedFiles_INF,
    returnMethod_decipher = c("Text","dataFrame"),
    method_entrez = c("Single","Multiple","Text"), smapName, termPath, 
    term, thresh = 5, limsize = 1000, win_indel_parents = 5000,
    win_inv_trans_parents = 40000, 
    outpath, outputFilename = ""  , RZIPpath)
    {
        if (method_entrez =="Text"){
            dat_geneList <- gene_list_generation(method_entrez = "Text",
                termPath=termPath,
                returnMethod = returnMethod_GeneList,thresh = 20)
        }
        else if (method_entrez == "Multiple"){
            dat_geneList <- gene_list_generation(method_entrez = "Multiple",
                term = term, returnMethod = returnMethod_GeneList,thresh = 20)
        }
        else if (method_entrez == "Single"){
            dat_geneList <- gene_list_generation(method_entrez = "Single",
                term = term, returnMethod = returnMethod_GeneList)
        }
        else if (method_entrez == "" && length(dat_geneList)> 0){
            dat_geneList <- dat_geneList
        }
        else{
            stop("GeneList or term required for analysis!!!")
        }
        datcompSmap<-compSmapbed(smap = smap,bed = bed, 
            inputfmtBed = inputfmtBed, n = n,
            returnMethod_bedcomp = returnMethod_bedcomp) 
        datDGV<- DGV_extraction(hgpath = hgpath, smap_data 
            = datcompSmap,input_fmt_DGV = input_fmt_DGV,
            win_indel_DGV = win_indel_DGV, 
            win_inv_trans_DGV = win_inv_trans_DGV, 
            perc_similarity_DGV = perc_similarity_DGV,returnMethod 
            = returnMethod_DGV)
        if(buildSVInternalDB==FALSE){
            datInf<-internalFrequency(mergedFiles 
                = mergedFiles_INF , buildSVInternalDB=FALSE, 
                smapdata=datDGV, input_fmt=input_fmt_INF, win_indel 
                = win_indel_INF , 
                win_inv_trans = win_inv_trans_INF, smap = smapName,
                perc_similarity =perc_similarity_INF, indelconf, invconf, 
                transconf, 
                limsize = limsize,returnMethod = returnMethod_Internal)
        }
        else{
            datInf<-internalFrequency(buildSVInternalDB = TRUE, 
                path = soloPath, pattern 
                = solopattern, outpath, dbOutput = dbOutput_Int,
                smapdata = datDGV, smap = smapName,input_fmt = input_fmt_INF,
                limsize = limsize, win_indel = win_indel_INF, 
                win_inv_trans = win_inv_trans_INF, 
                perc_similarity 
                = perc_similarity_INF, indelconf , invconf, transconf, 
                returnMethod= returnMethod_Internal)
          
    }
     
    
    if(buildBNInternalDB==FALSE){
        datchort<-cohortFrequency(internalBNDB 
            = mergedFiles_BN, buildBNInternalDB = FALSE, 
            smapdata = datInf, input_fmt 
            = input_fmt_BN, win_indel = win_indel_INF, 
            win_inv_trans = win_inv_trans_INF, 
            perc_similarity = perc_similarity_INF, 
            indelconf = indelconf, invconf = invconf, limsize=limsize,  
            transconf = transconf,returnMethod = returnMethod_BNCohort)
    }
    else{
        datchort<-cohortFrequency(buildBNInternalDB = TRUE, 
            smapdata = datInf, input_fmt = input_fmt_BN, 
            BNDBPath = InternalDBpath, 
            BNDBPattern = InternalDBpattern, win_indel = 
            win_indel_INF, dbOutput = dbOutput_BN,
            win_inv_trans = win_inv_trans_INF,
            perc_similarity = perc_similarity_INF, indelconf = indelconf, 
            invconf = invconf, limsize=limsize, 
            transconf = transconf,returnMethod = returnMethod_BNCohort)
          
    }
    datdecipher<-Decipher_extraction(decipherpath = 
        decipherpath, smap_data = datchort,
        input_fmt = input_fmt_decipher,win_indel = win_indel_INF, 
        perc_similarity = perc_similarity_INF, 
        returnMethod = returnMethod_decipher)
    
    
    run_bionano_filter(input_fmt_geneList = "dataFrame",
        input_fmt_svMap = input_fmt_svMap,
        SVFile = NULL,svData = datdecipher,dat_geneList = dat_geneList,
        outpath = outpath,outputFilename = outputFilename,RZIPpath = RZIPpath)
    }
    