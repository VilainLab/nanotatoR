#' Annotation and visualisation of Bionano SV.
#'
#' @param smap  character. Path to SMAP file.
#' @param inputfmtBed character. Choice between Text and DataFrame
#' as input for bed file.
#' @param bed  Text Choice between UCSC bed or Bionano bed.
#' @param n  numeric Number of genes to report which are nearest to the breakpoint.
#' Default is  	3.
#' @param mergedFiles  character. Path to the merged SV solo files.
#' @param buildSVInternalDB  boolean. Checking whether the merged solo 
#' file database exist or you need to build it. Default= TRUE.
#' @param soloPath  character. Path to the solo file database.
#' @param solopattern  character. pattern of the file names to merge.
#' @param smapName character. Name of the smap file.
#' @param win_indel_INF  Numeric. Insertion and deletion error window.
#' @param win_inv_trans_INF  Numeric. Inversion and translocation error window.
#' @param perc_similarity_INF  Numeric . ThresholdPercentage similarity 
#' of the query SV and reference SV.
#' @param input_fmt_DGV character. Choice between Text and DataFrame 
#' for input to DGV.
#' @param input_fmt_INF character. Choice between Text and DataFrame
#' for input to INF.
#' @param returnMethod_GeneList character. Return Methods from the gene_list_generation 
#' function, choice between Text and Dataframe.
#' @param returnMethod_bedcomp character. Return Methods from the compSmapbed 
#' function, choice between Text and Dataframe.
#' @param returnMethod_DGV character. Return Methods from the DGV_extraction function, choice between
#' Text and Dataframe.
#' @param returnMethod_Internal character. Return Methods from the internalFrequency function, choice between
#' Text and Dataframe.
#' @param limsize  Numeric. Threshold for size limit for the breakpoint, 
#' for checking that the breakpint size is valid or not. Default is 1000 bases.
#' @param indelconf Numeric. Threshold for insertion and deletion confidence.
#' @param invconf  Numeric. Threshold for inversion confidence.
#' @param transconf  Numeric. Threshold for translocation confidence.
#' @param hgpath  character. Path to Database of Genomic Variants (DGV)
#'                Text file. 
#' @param win_indel_DGV  Numeric. Insertion and deletion error window.
#' @param win_inv_trans_DGV  Numeric. Inversion and translocation error window.
#' @param perc_similarity_DGV  Numeric . ThresholdPercentage similarity 
#' of the query SV and reference SV.
#' @param method  character. Input Method for terms for entrez. Choices are 
#'                "Single","Multiple" and "Text".
#' @param termPath  character. Path and file name for textfile for terms.
#' @param term  character. Single or Multiple Terms as vectord.
#' @param thresh integer. Threshold for the number of terms sent to entrez.
#'                Note if large lists are sent to ncbi, it might fail to get
#'                processed. Default is 5.
#' @param input_fmt_geneList character. Choice of gene list input 
#'        Text or Dataframe.
#' @param input_fmt_svMap  character. Choice of gene list input 
#'        Text or Dataframe.
#' @param svData Dataframe Input data containing SV data.
#' @param dat_geneList Dataframe Input data containing geneList data.
#' @param outpath Character Directory to the output file.
#' @param outputFilename Character Output filename for the annotated smap.
#' @param RZIPpath Character. Path for the Rtools zip.exe
#' @return Excel file containing the annotated SV map, tabs divided based on
#' type of SVs.
#' @return Text files containg gene list and terms associated with them are stored as text files.
#' \dontrun{
#' @examples
#' terms="Muscle Weakness"
#' gene<-gene_list_generation(method="Single", term=terms, 
#'   	returnMethod_GeneList="dataFrame")
#' 
#' smapName="F1.1_GM24385_DLE-1_P_trio_hg19.smap"
#' smappath = system.file("extdata", smapName, package="nanotatoR")
#' nanotatoR_main(smap, bed, inputfmtBed = c("BNBED"), 
#'  n=3,mergedFiles , buildSVInternalDB=TRUE, soloPath, solopattern, 
#' input_fmt_INF=c("dataframe"),returnMethod_GeneList=c("dataframe"),
#' returnMethod_bedcomp=c("dataframe"),returnMethod_DGV=c("dataframe"),
#' returnMethod_Internal=c("dataframe"),input_fmt_DGV=c("dataframe"),	
#' hgpath, smapName,method=c("Single"), term, thresh=5,
#' input_fmt_geneList=c("dataframe"),input_fmt_svMap=c("dataframe"),
#' svData,dat_geneList,outpath="",outputFilename="",RZIPpath="")
#'	}
#' @export
nanotatoR_main<-function(smap, bed, inputfmtBed = c("BED", "BNBED"), 
    n=3,
   	mergedFiles , buildSVInternalDB=FALSE, soloPath, solopattern, 
	win_indel_INF = 10000, win_inv_trans_INF = 50000, 
	perc_similarity_INF= 0.5, indelconf = 0.5, invconf = 0.01, 
	transconf = 0.1,limsize=1000,input_fmt_INF=c("Text","dataframe"),
	returnMethod_GeneList=c("Text","dataframe"),returnMethod_bedcomp=c("Text","dataframe"),
	returnMethod_DGV=c("Text","dataframe"),returnMethod_Internal=c("Text","dataframe"),
	input_fmt_DGV=c("Text","dataframe"),
	hgpath, win_indel_DGV = 10000, win_inv_trans_DGV = 50000, 
    perc_similarity_DGV = 0.5,smapName,
	method=c("Single","Multiple","Text"), termPath="", 
    term, thresh=5,input_fmt_geneList=c("Text","dataframe"),
	input_fmt_svMap=c("Text","dataframe"),svData,dat_geneList,
	outpath="",outputFilename="",RZIPpath="")
	{
	
	 dat_geneList<-gene_list_generation(method="Text",termPath=termPath,returnMethod_GeneList="dataFrame")
	 datcompSmap<-compSmapbed(smap=smap,bed=bed, inputfmtBed="BED",n=3,returnMethod_bedcomp="dataFrame") 
	 datDGV<- DGV_extraction(hgpath=hgpath, smap_data=datcompSmap,input_fmt_DGV=c("dataFrame"),
     win_indel_DGV, win_inv_trans_DGV, perc_similarity_DGV,
	 returnMethod_DGV=c("dataFrame"))
	 if(buildSVInternalDB==FALSE){
	    datInf<-internalFrequency(mergedFiles=mergedFiles , buildSVInternalDB=FALSE, smapdata=datDGV, input_fmt_INF="dataFrame", win_indel_INF , win_inv_trans_INF,smapName=smapName,
		perc_similarity_INF, indelconf , invconf , transconf ,limsize=limsize,returnMethod_Internal="dataFrame")
		}
	else{
	    datInf<-internalFrequency(buildSVInternalDB=TRUE, soloPath=soloPath, solopattern=solopattern, outpath,  smapdata=datDGV, smapName=smapName,input_fmt_INF="dataFrame",limsize=limsize, win_indel_INF, win_inv_trans_INF, perc_similarity_INF, indelconf , invconf, transconf, returnMethod_Internal=c("dataFrame"))
	      
	}
	 
	
	run_bionano_filter(input_fmt_geneList="dataFrame",input_fmt_svMap="dataFrame",
                    svData=datInf,dat_geneList=dat_geneList,outpath="C://Annotator//Data//",outputFilename="F1.1_UDN287643_P",RZIPpath="")
	}
        	 
	