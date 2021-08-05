#' Merging DLE labelled smaps
#'
#' @param path character. Path to the solo files directory.
#' @param pattern    character. Pattern for the solo files.
#' @param outMode    character. The ouput mode. Choices, dataframe or Text.
#' @param outpath character. Path where the dual labelled 
#'        merged samples are kept. Is mandatory if outMode is Text.
#' @param pipeline Analysis pieline used. Options Rare variant Pipeline (RVP)
#' or Denovo Variant pipeline(DVP).
#' @return Text files containg merged smaps from different samples
#' @examples
#' mergedSmap <- mergingSMAP_SE (
#' path = system.file("extdata", "SoloFile/", package="nanotatoR"),
#' pattern = "*_DLE1_*", outMode = "dataframe", 
#' outpath = system.file("extdata", "Merged/", package="nanotatoR"))
#' @import hash
#' @import tidyverse
#' @export

mergingSMAP_internal <- function(
path ,pattern,outMode = c("Text", "dataframe"),
outpath, pipeline = c("DVP", "RVP"))
{
l <- list.files(path = path, pattern = pattern, full.names = TRUE)
nam <- c()
datfinal <- data.frame()
pipeline = pipeline
###Merging
for (ii in seq_along((l)))
{
    print(ii)
    smap = l[ii]
    con <- file(smap, "r")
    r10 <- readLines(con, n = -1)
    close(con)
    # datfinal<-data.frame()
    g1 <- grep("RefEndPos", r10)
    g2 <- grep("RefStartPos", r10)
    'gg1 <- grep("# BSPQI Sample", r10)
    stt <- strsplit(r10[gg1], split = ":")
    fname_temp <- stt[[1]][2]
        
    if(length(grep("UDN*", fname_temp)) ==1){
    ###UDN
    stt1 <- strsplit(fname_temp, split = "_P_BspQI_assembly*")
    fname <- stt1[[1]][1]
    } else{
    ###DSD
    stt1 <- strsplit(fname_temp, split = "_BspQI_assembly*")
    fname <- stt1[[1]][1]
    }
    
    stt1 <- strsplit(fname_temp, split = "_BspQI_assembly*")
    fname <- stt1[[1]][1]'
    
    
    #print (paste0("SampleName:", fname))
    if (g1 == g2) {
        dat <- gsub("#h ", "", as.character(r10))
        # dat<-gsub('\t',' ',r10)
        dat4 <- textConnection(dat[g1:length(dat)])
        r1 <- read.table(dat4, sep = "\t", header = TRUE)
        close(dat4)
    } else {
    stop("column names doesnot Match")
    }
	if(length(unique(r1$Sample)) > 1){
	    r1$Sample <- gsub("^-$", "ExperimentLabel", r1$Sample)
	}else{r1$Sample <- r1$Sample}
	if(any(unique(r1$Sample) == "ExperimentLabel") == TRUE){
	    g1 <- strsplit(smap, split = "/")
		g2 <- strsplit(g1[[1]][length(g1[[1]])], split = ".smap")
		#g3 <- strsplit(g1[[1]][length(g1[[1]])], split = "_")
	    #Samp <- as.character(g2[[1]][1])
        Samp <- as.character(g2[[1]][1])
	}else{
	   	Samp <- as.character(unique(r1$Sample))
	}
    'st1 <- strsplit(Samp, split = "*_DLE")
    SampleID <- st1[[1]][1]'
    if(length(grep("*_BspQI_*", Samp)) >= 1){
        stt1 <- strsplit(Samp, split = "*_BspQI")
        fname_temp <- stt1[[1]][1]
    }else if(length(grep("*_BssSI_*", Samp)) >= 1){
        stt1 <- strsplit(Samp, split = "*_BssSI")
        fname_temp <- stt1[[1]][1]
    }else if(length(grep("*_DLE_*", Samp)) >= 1){
        stt1 <- strsplit(Samp, split = "*_DLE")
        fname_temp <- stt1[[1]][1]
    }else{print("Sample pattern not Found")}
    SampleID <- fname_temp
    r1 <- cbind(SampleID = rep(str_squish(as.character(SampleID)), times = nrow(r1)), r1)
    #### print(dim(r1))
    datfinal <- data.frame(rbind(datfinal, r1))
    'if(length(grep(names(r1), "Sample")) >= 1){
        #fname2 <- unique(r1$Sample)
        if(length(grep(fname2, 
        stt1 <- strsplit(fname_temp, split = "_BspQI_assembly*")
        fname_temp <- stt1[[1]][1]
    } else{fname2 <- fname}'
    ''
    #datfinal <- datfinal
    #### print(dim(datfinal))
    
}
        Type1 <- as.character(c(rep ("-", times = nrow(datfinal))))
        Type2 <- as.character(c(rep ("-", times = nrow(datfinal))))
        Fail_BSPQI_assembly_chimeric_score = as.character(
            c(rep ("-", times = nrow(datfinal)))
            )
        Fail_BSSSI_assembly_chimeric_score = as.character(
            c(rep ("-", times = nrow(datfinal)))
            )
        Found_in_self_BSPQI_molecules = as.character(
            c(rep ("-", times = nrow(datfinal)))
            )
        Found_in_self_BSSSI_molecules = as.character(
            c(rep ("-", times = nrow(datfinal)))
            )
        Method <- as.character(c(rep ("SE", times = nrow(datfinal))))
        SampleID = datfinal$SampleID
		if(pipeline == "RVP"){
		    Zygosity = as.character(c(rep ("-", times = nrow(datfinal))))
		}
		else{
		    Zygosity = as.character(datfinal$Zygosity)
		}	
        'strsample <- strsplit(as.character(datfinal$Sample), 
            split = "*_DLE_*")
        SampleID = sapply(strsample, function(x) x[1])'
        datFinal_DLE <- data.frame(
		    SampleID = datfinal$SampleID,
		    datfinal[,1:10],
			datfinal[,27:29],
            Zygosity = Zygosity,
            datfinal[,31:35],			
            Type1 = Type1, Type2 = Type2, 
            Fail_BSPQI_assembly_chimeric_score = Fail_BSPQI_assembly_chimeric_score, 
            Fail_BSSSI_assembly_chimeric_score =Fail_BSSSI_assembly_chimeric_score,
            Fail_assembly_chimeric_score_SE = as.character(
                datfinal$Fail_assembly_chimeric_score
            ),
			datfinal[,37],
            Found_in_self_BSPQI_molecules = Found_in_self_BSPQI_molecules, 
            Found_in_self_BSSSI_molecules = Found_in_self_BSPQI_molecules,
            Found_in_SE_self_molecules = as.character(
                datfinal$Found_in_self_molecules
                ), 
            Method = Method,
			datfinal[,38:39]
			
        )
        ##Writing output dataframe or text
        if (outMode == "Text"){
            len <- length(l)
            fileOut <- paste0("DLE_1_",len,".txt")
            write.table(datFinal_DLE, file.path(outpath, fileOut), 
                row.names = FALSE
                )
        }
        else if (outMode == "dataframe"){
            return(datFinal_DLE)
        }
        else {
            stop ("Your Outmode is absent !!!")
        }
    }
