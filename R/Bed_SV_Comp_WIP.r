#' Reads Bionano Bedfiles
#'
#' @param BNFile  character. Path to Bionano Bed File.
#' @return Data Frame Contains the gene information.
#' @examples
#' BNFile <- system.file("extdata", "Homo_sapiens.GRCH19_BN.bed", package="nanotatoR")
#' bed<-readBNBedFiles(BNFile)
#' @import utils
#' @export
readBNBedFiles <- function(BNFile) {
    ## Reading the Bed File con<-file(BNFile,'r') r10<-readLines(con,n=-1)
    ## close(con) Converting the data to data frame
    ## dat4<-textConnection(r10)
    ## r12<-read.table(dat4,sep='\t',header=FALSE) Extracting data
    r12 <- read.table(BNFile, header = TRUE)
    chrom <- r12[, 1]
    chromstart <- r12[, 2]
    chromend <- r12[, 3]
    gene <- as.character(r12[, 4])
    strand <- as.character(r12[, 5])
    ## Combining the data to form a data frame
    dat1 <- data.frame(
        Chromosome = chrom, Chromosome_Start = chromstart,
        Chromosome_End = chromend, Gene = gene, Strand = strand
    )
    return(dat1)
}


#' Reads BED files to produce bionano Bed files
#'
#' @param bedFile    character. Path to UCSC Bed File.
#' @param outdir    character. Path to output directory.
#' @param returnMethod    character. Path to output directory.
#' @return Data Frame or text file. Contains the gene information.
#' @examples
#' bedFile <- system.file("extdata", "Homo_sapiens.Hg19.bed", 
#'       package="nanotatoR")
#' bed<-buildrunBNBedFiles(bedFile,returnMethod="dataFrame")
#' @import utils
#' @import stringr
#' @export
buildrunBNBedFiles <- function(bedFile, returnMethod = c("Text", "dataFrame"), 
    outdir) {
    ## Reading the Bed File
    con <- file(bedFile, "r")
    r10 <- readLines(con, n = -1)
    close(con)
    ## Converting the data to data frame
    dat4 <- textConnection(r10)
    r12 <- read.table(dat4, sep = "\t", header = FALSE)
    ## Extracting data
    # print(dim(r12))
    chrom <- stringr::str_trim(r12[, 1])
    chromstart <- stringr::str_trim(r12[, 2])
    chromend <- stringr::str_trim(r12[, 3])
    gene <- stringr::str_trim(as.character(r12[, 4]))
    strand <- stringr::str_trim(as.character(r12[, 5]))
    ## Changing the chromosome Start Name
    chrom1 <- gsub("chr", "", x = chrom)
    ## Male and Female chromosome association
    if (length(grep("X", chrom1)) > 1 & length(grep("Y", chrom1)) > 1) {
        chrom1 <- gsub("X", 23, chrom1)
        chrom1 <- gsub("Y", 24, chrom1)
    } else if (length(grep("X", chrom1)) > 1) {
        chrom1 <- gsub("X", 23, chrom1)
    } else {
        print("Genome doesnot have any Sex Chromosome")
    }
    ## Combining the data to form a data frame
    num <- seq_len(length(strand))
    dat1 <- data.frame(
        Chromosome = as.character(chrom1),
        Chromosome_Start = as.numeric(chromstart),
        Chromosome_End = as.numeric(chromend), 
        Gene = as.character(gene),
        num, Strand = as.character(strand),
        chromStart = as.numeric(chromstart), chromEnd = as.numeric(chromend),
        colors = rep("128,0,128", length(strand)), row.names = NULL
    )
    ## writing Bionano Bedfiles
    if (returnMethod == "Text") {
        if (length(grep("\\\\", bedFile)) >= 1) {
            st <- strsplit(bedFile, split = "\\\\")
            fname <- st[[1]][4]
            st1 <- strsplit(fname, split = ".bed")
            fname1 <- paste(st1[[1]][1], "_BN.bed", sep = "")
        } else {
            fname <- bedFile
            st1 <- strsplit(fname, split = ".bed")
            fname1 <- paste(st1[[1]][1], "_BN.bed", sep = "")
        }
        write.table(
            dat1, paste(outdir, "/", fname1, sep = ""),
            row.names = FALSE, col.names = FALSE
        )
    } else if (returnMethod == "dataFrame") {
        return(dat1)
    } else {
        stop("Method of Return improper")
    }
}
#' Reads SMAP files to extract information
#'
#' @param smap    character. Path to SMAP file.
#' @return Data Frame or text file. Contains the SMAP information.
#' @examples
#' smapName="F1.1_GM24385_DLE-1_P_trio_hg19.smap"
#' smap = system.file("extdata", smapName, package="nanotatoR")
#' readSMap(smap)
#' @import utils
#' @export

readSMap <- function(smap, input_fmt_smap = c("Text","dataFrame")) {
    ## reading the smap text file
        if(input_fmt_smap == "Text"){
            con <- file(smap, "r")
            r10 <- readLines(con, n = -1)
            close(con)
            # datfinal<-data.frame()
            g1 <- grep("RawConfidence", r10)
            g2 <- grep("RefStartPos", r10)
            
            gg1 <- grep("# BSPQI Sample", r10)
			gg2 <- grep("# BSSSI Sample", r10)
	        if(length(gg1) > 0 & length(gg2) > 0){
                stt <- strsplit(r10[gg1], split = ":")
				stt35 <- strsplit(r10[gg2], split = ":")
				fname_temp <- stt[[1]][2]
				fname_temp1 <- stt35[[1]][2]
                stt1 <- strsplit(fname_temp, split = "_BspQI_assembly*")
				fname1 <- stt1[[1]][1]
                stt2 <- strsplit(fname_temp1, split = "_BssSI_assembly*")
                fname2<- stt2[[1]][1]
	        	if(fname1 == fname2){
                    fname <- str_squish(fname2)
	        	} else{stop("Mismatch in File Names")}
                g1 <- grep("RefEndPos", r10)
                g2 <- grep("RefStartPos", r10)
                
                # r10<-as.character(r10)
                if (g1 == g2){
                  #g3 <- grep("# ",r10)
                    dat <- gsub("# ", "", as.character(r10))
                    
                    # dat<-gsub('\t',' ',as.character(dat))
                    dat4 <- textConnection(dat[g1:length(dat)])
                    ##### print(dim(dat4))
                    r1 <- read.table(dat4, sep = "\t", header = TRUE)
                    close(dat4)
	        	}else{stop("File Incorrect!!!")}
	        }else{
	            g1 <- grep("RefEndPos", r10)
                g2 <- grep("RefStartPos", r10)
                
                # r10<-as.character(r10)
                if (g1 == g2){
                  #g3 <- grep("#h ",r10)
                    dat <- gsub("#h ", "", as.character(r10))
                    
                    # dat<-gsub('\t',' ',as.character(dat))
                    dat4 <- textConnection(dat[g1:length(dat)])
                    ##### print(dim(dat4))
                    r1 <- read.table(dat4, sep = "\t", header = TRUE)
                    close(dat4)
	        	}else{stop("File Incorrect!!!")}
	            fname_temp <- as.character(unique(r1$Sample))
	        	g2 <- grep("*_BspQI*", fname_temp)
	        	g3 <-  grep( "*_BssSI*", fname_temp)
	        	if(length(g2)> 0 & length(g3)>0){
	        	    stt1 <- strsplit(fname_temp, split = "*_BspQI*")
	        	    stt2 <- strsplit(fname_temp, split = "*_BssSI*")
	        	    fname1 <- stt1[[1]][1]
	        	    fname2<- stt2[[1]][1]
	        		if(fname1 == fname2){
                        fname <- str_squish(fname1)
	        	    }else{stop("Mismatch in File Names")}
	        	}else if (length(g2)> 0 & length(g3) == 0){
	        	    stt1 <- strsplit(fname_temp, split = "*_BspQI*")
	        	    fname1 <- stt1[[1]][1]
	        	    fname <- str_squish(fname1)
	        	}else if (length(g2) == 0 & length(g3) > 0){
	        	    stt2 <- strsplit(fname_temp, split = "*_BssSI*")
	        	    fname2<- stt2[[1]][1]
	        	    fname <- str_squish(fname2)
	        	}else{stop("No Sample ID")}
	        }
		    sampnam <- grep("SampleID", names(r1))
		    if(length(sampnam) != 1){
                r1 <- cbind(SampleID = rep(
                    str_squish(as.character(fname)), 
                    times = nrow(r1)), r1
                    )
            } else { r1 <- r1 }
			
        }else if(input_fmt_smap == "dataFrame"){
                r1<-smapdata
        }
        else{
                stop("Input Format incorrect")
        }
        #sampID <- str_squish(as.character(unique(r1$SampleID)))
    return(r1)
}
#' Reads DLE SMAP files to extract information
#'
#' @param smap    character. Path to SMAP file.
#' @return Data Frame or text file. Contains the SMAP information.
#' @examples
#' smapName="F1.1_GM24385_DLE-1_P_trio_hg19.smap"
#' smap = system.file("extdata", smapName, package="nanotatoR")
#' readSMap(smap)
#' @import utils
#' @export

readSMap_DLE <- function(smap, input_fmt_smap = "Text") {
  ## reading the smap text file
        if(input_fmt_smap == "Text"){
            con <- file(smap, "r")
            r10 <- readLines(con, n = -1)
            close(con)
            # datfinal<-data.frame()
            g1 <- grep("RawConfidence", r10)
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
        Samp <- as.character(unique(r1$Sample))
        st1 <- strsplit(Samp, split = "*_DLE")
        SampleID <- st1[[1]][1]
        r1 <- cbind(SampleID = rep(str_squish(as.character(SampleID)),
            times = nrow(r1)), r1)
        }
        else if(input_fmt_smap == "dataFrame"){
            r1<-smapdata
        }
        else{
            stop("Input Format incorrect")
        }
    #sampID <- str_squish(as.character(unique(r1$SampleID)))
  return(r1)
}

#' Calculates Genes that overlap the SV region
#'
#' @param bed    Text Bionano Bed file.
#' @param chrom    character SVmap chromosome.
#' @param startpos    numeric starting position of the breakpoints.
#' @param endpos    numeric end position of the breakpoints.
#' @param svid    numeric Structural variant identifier (Bionano generated).
#' @return Data Frame. Contains the SVID,Gene name,strand information and
#' percentage of SV covered.
#' @examples
#' smapName="F1.1_GM24385_DLE-1_P_trio_hg19.smap"
#' smap = system.file("extdata", smapName, package="nanotatoR")
#' bedFile <- system.file("extdata", "Homo_sapiens.Hg19.bed", 
#'        package="nanotatoR")
#' bed<-buildrunBNBedFiles(bedFile,returnMethod="dataFrame")
#' smap<-readSMap(smap)
#' chrom<-smap$RefcontigID1
#' startpos<-smap$RefStartPos
#' endpos<-smap$RefEndPos
#' if (length(grep("SVIndex",names(smap)))>0){
#'        svid <- smap$SVIndex
#'    }else{
#'     svid <- smap$SmapEntryID
#'     }
#' overlapGenes(bed, chrom, startpos, endpos, svid)
#' @import utils
#' @export
overlapGenes <- function(bed, chrom, startpos, endpos, svid, chrom2, SVTyp,
    bperrorindel = 3000, bperrorinvtrans = 10000) {
    ## Initialising data
    print("***Overlap Genes***")
    data1 <- data.frame()
    gnsInf <- c()
    SVID <- c()
	chrom2 = chrom2
	chrom = chrom
    ## Getting the midpoint of the geme bed$geneMid <-
    ## ceiling(bed$Chromosome_Start + ((bed$Chromosome_End -
    ## bed$Chromosome_Start)/2)) Detecting Genes present between the SV
    ## breakpoints and Calculating the percentage coverage of genes across
    ## the breakpoints per chromosome.
    for (ii in seq_len(length(chrom)))
    #for (ii in 1:5)
    { #print(paste("OverLap:",ii))
        # Checking for genes in the breakpoint
		
		if(startpos[ii] == -1 | endpos[ii] == -1){
		    #print(ii)
			SVID <- c(SVID, svid[ii])
		    gnsInf = c(gnsInf,"-")
        }
        else{		
            dat10 <- bed[which(bed$Chromosome == chrom[ii]), ]
			dat15 <- bed[which(bed$Chromosome == chrom2[ii]), ]
            if(SVTyp[ii] == "insertion" 
			    | SVTyp[ii] == "deletion"
				| SVTyp[ii] == "duplication"
				| SVTyp[ii] == "duplication_split"
				| SVTyp[ii] == "duplication_inverted"){
				start_adj <- startpos[ii] - bperrorindel
				start_adj1 <- startpos[ii] + bperrorindel
				end_adj <- endpos[ii] + bperrorindel
				end_adj1 <- endpos[ii] - bperrorindel
			    dat11 <- dat10[which(((
					    (dat10$Chromosome_Start >= start_adj 
					    & dat10$Chromosome_Start <= start_adj1)
					    & (dat10$Chromosome_End <= end_adj
						& dat10$Chromosome_End >= end_adj1))
						| (dat10$Chromosome_Start <= start_adj 
						& dat10$Chromosome_End >= end_adj)
						| (dat10$Chromosome_Start >= start_adj1 
						& dat10$Chromosome_End <= end_adj1)
						| ((dat10$Chromosome_Start >= start_adj1
						& dat10$Chromosome_Start <= end_adj1) 
                        & dat10$Chromosome_End > end_adj)
						| (dat10$Chromosome_Start <= start_adj
						& (dat10$Chromosome_End >= start_adj1
						& dat10$Chromosome_End <= end_adj1))
						| (dat10$Chromosome_End >= end_adj
						& (dat10$Chromosome_End >= end_adj1
						& dat10$Chromosome_End <= end_adj))
						| (dat10$Chromosome_Start <= start_adj
						& (dat10$Chromosome_End >= start_adj
						& dat10$Chromosome_End <= start_adj1)))), ]
				if (nrow(dat11) > 1) {
                    ## Extracting strand and chromosome start information
                    chromos <- dat11$Chromosome
                    chromStart <- dat11$Chromosome_Start
                    chromEnd <- dat11$Chromosome_End
                    gen1 <- dat11$Gene
                    strnd <- dat11$Strand
                    genMid <- dat11$geneMid
                    geneInfo <- c()
                    ## Calculating the coverage of the gene
                    for (k in seq_len(length(gen1)))
                    {
                        ##Checking if orientation of the gene and calculating 
                        ##length based on that
                        if(chromStart[k]>=chromEnd[k]){
                            lengthgene <- chromStart[k] - chromEnd[k]
                        }else{
                            lengthgene <- chromEnd[k]-chromStart[k]
                        }
                        lengthbp <- end_adj-start_adj
                        # print(paste('distfromStart:',distfromStart),sep='')
                        # print(paste('distfromEnd:',distfromEnd),sep='') Checking for
                        # partial/complete inclusion in the SV
                        if ((chromStart[k] >= start_adj & chromEnd[k] <= end_adj)) {
                            percentage <- 100
                        } 
                        else if ((chromStart[k] < start_adj 
                            & chromEnd[k] > end_adj) 
                            & (lengthgene > lengthbp)) {
                            percentage <- round(abs(((start_adj - end_adj) / (chromStart[k] - chromEnd[k])) * 100), digits = 2)
                        } 
                        else if (((chromStart[k] < start_adj) 
                                    & (chromEnd[k] > start_adj) 
                                    & (chromEnd[k] <= end_adj))) {
                                        distfromStart <- chromEnd[k] - start_adj
                                        percentage <- round(abs((
                                        distfromStart / (chromEnd[k] - chromStart[k])) * 100), 
                                        digits = 2)
                        } 
                        else if (((chromStart[k] >= start_adj) 
                                    & (chromStart[k] < end_adj) 
                                    & (chromEnd[k] > end_adj))) {
                                        distfromEnd <- end_adj - chromStart[k]
                                        percentage <- round(abs((
                                            distfromEnd / (chromEnd[k] - chromStart[k])) * 100), digits = 2)
                        } 
                        else {
                            percentage <- "NA"
                        }
                        if (is.na(percentage)) {
                            geneInfo <- "-"
                        }
                        else {
                            geneInfo <- c(geneInfo, paste(
                                gen1[k], "(", strnd[k], ":",
                                percentage, ")", sep = ""
                            ))
                        }
                    }
                    geneInfo1 <- paste(geneInfo, collapse = ";")
                    #print(geneInfo1)
                    gnsInf <- c(gnsInf, str_squish(str_trim((geneInfo1),side="both")))
                    # print(paste("1:",ii,":",svid[ii],sep=""))
                    SVID <- c(SVID, svid[ii])
                } else if (nrow(dat11) == 1) {
                    chromos <- dat11$Chromosome
                    chromStart <- dat11$Chromosome_Start
                    chromEnd <- dat11$Chromosome_End
                    gen1 <- dat11$Gene
                    strnd <- dat11$Strand
                    genMid <- dat11$geneMid
                    geneInfo <- c()
                    distfromStart <- start_adj - chromStart
                    distfromEnd <- end_adj - chromEnd
                    lengthgene <- abs(chromStart - chromEnd)
                    lengthbp <- abs(start_adj - end_adj)
		        
                    if ((chromStart >= start_adj & chromEnd <= end_adj)) {
                            percentage <- 100
                        } else if ((chromStart < start_adj & chromEnd >
                            end_adj) & (lengthgene > lengthbp)) {
                            percentage <- round(abs(((
                                    start_adj - end_adj) / (chromStart - chromEnd)) * 100
                                    ), digits = 2)
                        } else if (((chromStart < start_adj 
                                    & chromEnd > start_adj)
                                    & (chromEnd <= end_adj))) {
                                    distfromStart <- chromEnd - start_adj
                                    percentage <- round(abs((
                                        distfromStart / (chromEnd - chromStart)) * 100
                                        ), digits = 2)
                        } else if (((chromStart >= start_adj 
                                    & chromStart < end_adj) 
                                    & (chromEnd > end_adj))) {
                                        distfromEnd <- end_adj - chromStart
                                        percentage <- round(
                                            abs((distfromEnd / (chromEnd - chromStart)) * 100
                                            ), digits = 2)
                        } else {
                            percentage <- "NA"
                        }
                        if (is.na(percentage)) {
                            geneInfo <- "-"
                        }
                        else {
                            geneInfo <- c(geneInfo, paste(
                                gen1, "(", strnd, ":",
                                percentage, ")", sep = ""
                            ))
                        }
                #print(geneInfo)
                    gnsInf <- c(gnsInf, str_squish(str_trim((geneInfo),side="both")))
                    SVID <- c(SVID, svid[ii])
                } else {
                    SVID <- c(SVID, svid[ii])
                    gnsInf <- c(gnsInf, "-")
                    # print(paste("1:",ii,":",svid[ii],sep=""))
                }
		    }else if((length(grep("inversion", SVTyp[ii])) >= 1)){
					start_adj <- startpos[ii] - bperrorinvtrans
					start_adj1 <- startpos[ii] + bperrorinvtrans
				    end_adj <- endpos[ii] + bperrorinvtrans
					end_adj1 <- endpos[ii] - bperrorinvtrans
			        dat9 <- dat10[which(
					    (dat10$Chromosome_Start <= start_adj 
					    & dat10$Chromosome_End >= start_adj1)
					    | (dat10$Chromosome_Start <= start_adj 
					    & (dat10$Chromosome_End >= start_adj 
						& dat10$Chromosome_End <= start_adj1))
						| (dat10$Chromosome_End >= start_adj1
					    & (dat10$Chromosome_Start >= start_adj 
						& dat10$Chromosome_Start <= start_adj1))), ]
					dat8 <- dat10[which(
					    (dat10$Chromosome_Start <= end_adj1 
					    & dat10$Chromosome_End >= end_adj)
					    | (dat10$Chromosome_Start <= end_adj1 
					    & (dat10$Chromosome_End >= end_adj1 
						& dat10$Chromosome_End <= end_adj))
						| (dat10$Chromosome_End >= end_adj
					    & (dat10$Chromosome_Start >= end_adj1 
						& dat10$Chromosome_Start <= end_adj))), ]
					dat11 <- rbind(dat8, dat9)
				    if (nrow(dat11) > 1) {
                        ## Extracting strand and chromosome start information
                        chromos <- dat11$Chromosome
                        chromStart <- dat11$Chromosome_Start
                        chromEnd <- dat11$Chromosome_End
                        gen1 <- dat11$Gene
                        strnd <- dat11$Strand
                        genMid <- dat11$geneMid
                        geneInfo <- c()
                        ## Calculating the coverage of the gene
                        for (k in seq_len(length(gen1)))
                        {
                            ##Checking if orientation of the gene and calculating 
                            ##length based on that
                            if(chromStart[k]>=chromEnd[k]){
                                lengthgene <- chromStart[k] - chromEnd[k]
                            }else{
                                lengthgene <- chromEnd[k]-chromStart[k]
                            }
                            lengthbp1 <- start_adj1 - start_adj
							lengthbp2 <- end_adj - end_adj1
                            # print(paste('distfromStart:',distfromStart),sep='')
                            # print(paste('distfromEnd:',distfromEnd),sep='') Checking for
                            # partial/complete inclusion in the SV
                            if ((chromStart[k] >= start_adj & chromEnd[k] <= start_adj1) 
							    | (chromStart[k] >= end_adj1 & chromEnd[k] <= end_adj)) {
                                percentage <- 100
                            } 
                            else if ((chromStart[k] < start_adj 
                                & chromEnd[k] > start_adj1) 
                                & (lengthgene > lengthbp1)) {
                                percentage <- round(abs(((start_adj1 - start_adj) / (chromStart[k] - chromEnd[k])) * 100), digits = 2)
                            } 
							else if ((chromStart[k] < end_adj1
                                & chromEnd[k] > end_adj) 
                                & (lengthgene > lengthbp2)) {
                                percentage <- round(abs(
								    ((end_adj - end_adj1) / (chromStart[k] - chromEnd[k])) * 100), digits = 2)
                            } 
                            else if (((chromStart[k] < start_adj) 
                                        & (chromEnd[k] > start_adj) 
                                        & (chromEnd[k] <= start_adj1))) {
                                            distfromStart <- chromEnd[k] - start_adj
                                            percentage <- round(abs((
                                            distfromStart / (chromEnd[k] - chromStart[k])) * 100), 
                                            digits = 2)
                            } 
							else if (((chromStart[k] < end_adj1) 
                                        & (chromEnd[k] > end_adj1) 
                                        & (chromEnd[k] <= end_adj))) {
                                            distfromStart <- chromEnd[k] - end_adj1
                                            percentage <- round(abs((
                                            distfromStart / (chromEnd[k] - chromStart[k])) * 100), 
                                            digits = 2)
                            } 
                            else if (((chromStart[k] >= start_adj) 
                                        & (chromStart[k] < start_adj1
                                        & chromEnd[k] > start_adj1))) {
                                            distfromEnd <- start_adj1 - chromStart[k]
                                            percentage <- round(abs((
                                                distfromEnd / (chromEnd[k] - chromStart[k])) * 100), digits = 2)
                            } 
							else if (((chromStart[k] >= end_adj1) 
                                        & (chromStart[k] < end_adj
                                        & chromEnd[k] > end_adj))) {
                                            distfromEnd <- end_adj - chromStart[k]
                                            percentage <- round(abs((
                                                distfromEnd / (chromEnd[k] - chromStart[k])) * 100), digits = 2)
                            }
                            else {
                                percentage <- "NA"
                            }
                            if (is.na(percentage)) {
                                geneInfo <- "-"
                            }
                            else {
                                geneInfo <- c(geneInfo, paste(
                                    gen1[k], "(", strnd[k], ":",
                                    percentage, ")", sep = ""
                                ))
                            }
                        }
                        geneInfo1 <- paste(geneInfo, collapse = ";")
                        #print(geneInfo1)
                        gnsInf <- c(gnsInf, str_squish(str_trim((geneInfo1),side="both")))
                        # print(paste("1:",ii,":",svid[ii],sep=""))
                        SVID <- c(SVID, svid[ii])
                    } else if (nrow(dat11) == 1) {
                        chromos <- dat11$Chromosome
                        chromStart <- dat11$Chromosome_Start
                        chromEnd <- dat11$Chromosome_End
                        gen1 <- dat11$Gene
                        strnd <- dat11$Strand
                        genMid <- dat11$geneMid
                        geneInfo <- c()
                        distfromStart <- start_adj - chromStart
                        distfromEnd <- end_adj - chromEnd
                        'lengthgene <- abs(chromStart - chromEnd)
                        lengthbp <- abs(start_adj - end_adj)'
		                
                        if(chromStart >= chromEnd){
                                lengthgene <- chromStart - chromEnd
                            }else{
                                lengthgene <- chromEnd - chromStart
                            }
                            lengthbp1 <- start_adj1 - start_adj
							lengthbp2 <- end_adj - end_adj1
                            # print(paste('distfromStart:',distfromStart),sep='')
                            # print(paste('distfromEnd:',distfromEnd),sep='') Checking for
                            # partial/complete inclusion in the SV
                            if ((chromStart >= start_adj & chromEnd <= start_adj1) 
							    | (chromStart >= end_adj1 & chromEnd <= end_adj)) {
                                percentage <- 100
                            } 
                            else if ((chromStart < start_adj 
                                & chromEnd > start_adj1) 
                                & (lengthgene > lengthbp1)) {
                                percentage <- round(abs(((start_adj1 - start_adj) / (chromStart - chromEnd)) * 100), digits = 2)
                            } 
							else if ((chromStart < end_adj1
                                & chromEnd > end_adj) 
                                & (lengthgene > lengthbp2)) {
                                percentage <- round(abs(
								    ((end_adj - end_adj1) / (chromStart - chromEnd)) * 100), digits = 2)
                            } 
                            else if (((chromStart < start_adj) 
                                        & (chromEnd > start_adj) 
                                        & (chromEnd <= start_adj1))) {
                                            distfromStart <- chromEnd - start_adj
                                            percentage <- round(abs((
                                            distfromStart / (chromEnd - chromStart)) * 100), 
                                            digits = 2)
                            } 
							else if (((chromStart < end_adj1) 
                                        & (chromEnd > end_adj1) 
                                        & (chromEnd <= end_adj))) {
                                            distfromStart <- chromEnd - end_adj1
                                            percentage <- round(abs((
                                            distfromStart / (chromEnd - chromStart)) * 100), 
                                            digits = 2)
                            } 
                            else if (((chromStart >= start_adj) 
                                        & (chromStart < start_adj1
                                        & chromEnd > start_adj1))) {
                                            distfromEnd <- start_adj1 - chromStart
                                            percentage <- round(abs((
                                                distfromEnd / (chromEnd - chromStart)) * 100), digits = 2)
                            } 
							else if (((chromStart >= end_adj1) 
                                        & (chromStart < end_adj
                                        & chromEnd > end_adj))) {
                                            distfromEnd <- end_adj - chromStart
                                            percentage <- round(abs((
                                                distfromEnd / (chromEnd - chromStart)) * 100), digits = 2)
                            }
                            else {
                                percentage <- "NA"
                            }
                            if (is.na(percentage)) {
                                geneInfo <- "-"
                            }
                            else {
                                geneInfo <- c(geneInfo, paste(
                                    gen1[k], "(", strnd[k], ":",
                                    percentage, ")", sep = ""
                                ))
                            }
                    #print(geneInfo)
                        gnsInf <- c(gnsInf, str_squish(str_trim((geneInfo),side="both")))
                        SVID <- c(SVID, svid[ii])
                    } else {
                        SVID <- c(SVID, svid[ii])
                        gnsInf <- c(gnsInf, "-")
                        # print(paste("1:",ii,":",svid[ii],sep=""))
                    }
				    
					
				}else if(SVTyp[ii] == "translocation_interchr"
				    | (SVTyp[ii] == "translocation"
					| SVTyp[ii] == "translocation_intrachr")){
					start_adj <- startpos[ii] - bperrorinvtrans
					start_adj1 <- startpos[ii] + bperrorinvtrans
				    end_adj <- endpos[ii] + bperrorinvtrans
					end_adj1 <- endpos[ii] - bperrorinvtrans
			        dat9 <- dat10[which(
					    (dat10$Chromosome_Start <= start_adj 
					    & dat10$Chromosome_End >= start_adj1)
					    | (dat10$Chromosome_Start <= start_adj 
					    & (dat10$Chromosome_End >= start_adj 
						& dat10$Chromosome_End <= start_adj1))
						| (dat10$Chromosome_End >= start_adj1
					    & (dat10$Chromosome_Start >= start_adj 
						& dat10$Chromosome_Start <= start_adj1))), ]
					dat8 <- dat15[which(
					    (dat15$Chromosome_Start <= end_adj1 
					    & dat15$Chromosome_End >= end_adj)
					    | (dat15$Chromosome_Start <= end_adj1 
					    & (dat15$Chromosome_End >= end_adj1 
						& dat15$Chromosome_End <= end_adj))
						| (dat15$Chromosome_End >= end_adj
					    & (dat15$Chromosome_Start >= end_adj1 
						& dat15$Chromosome_Start <= end_adj))), ]
					#dat11 <- rbind(dat8, dat9)
					if (nrow(dat8) > 1 & nrow(dat9) > 1) {
                        ## Extracting strand and chromosome start information
						print("1st chromosome")
                        chromos1 <- dat9$Chromosome
                        chromStart1 <- dat9$Chromosome_Start
                        chromEnd1 <- dat9$Chromosome_End
                        gen1 <- dat9$Gene
                        strnd1 <- dat9$Strand
                        genMid1 <- dat9$geneMid
						print("2nd chromosome")
						chromos2 <- dat8$Chromosome
                        chromStart2 <- dat8$Chromosome_Start
                        chromEnd2 <- dat8$Chromosome_End
                        gen2 <- dat8$Gene
                        strnd2 <- dat8$Strand
                        genMid2 <- dat8$geneMid
                        geneInfo <- c()
						## Calculating the coverage of the gene
                        for (k in seq_len(length(gen1)))
                        {
                            ##Checking if orientation of the gene and calculating 
                            ##length based on that
                            if(chromStart1[k]>=chromEnd1[k]){
                                lengthgene1 <- chromStart1[k] - chromEnd1[k]
                            }else{
                                lengthgene1 <- chromEnd1[k]-chromStart1[k]
                            }
							
                            lengthbp1 <- start_adj1 - start_adj
							lengthbp2 <- end_adj - end_adj1
                            # print(paste('distfromStart:',distfromStart),sep='')
                            # print(paste('distfromEnd:',distfromEnd),sep='') Checking for
                            # partial/complete inclusion in the SV
                            if ((chromStart1[k] >= start_adj & chromEnd1[k] <= start_adj1)) {
                                percentage <- 100
                            } 
                            else if ((chromStart1[k] < start_adj 
                                & chromEnd1[k] > start_adj1) 
                                & (lengthgene1 > lengthbp1)) {
                                percentage <- round(abs(((start_adj1 - start_adj) / (chromEnd1[k]-chromStart1[k])) * 100), digits = 2)
                            } 
							
                            else if (((chromStart1[k] < start_adj) 
                                        & (chromEnd1[k] > start_adj) 
                                        & (chromEnd1[k] <= start_adj1))) {
                                            distfromStart <- chromEnd1[k] - start_adj
                                            percentage <- round(abs((
                                            distfromStart / (chromEnd1[k] - chromStart1[k])) * 100), 
                                            digits = 2)
                            } 
							else if (((chromStart2[k] < end_adj1) 
                                        & (chromEnd2[k] > end_adj1) 
                                        & (chromEnd2[k] <= end_adj))) {
                                            distfromStart <- chromEnd2[k] - end_adj1
                                            percentage <- round(abs((
                                            distfromStart / (chromEnd2[k] - chromStart2[k])) * 100), 
                                            digits = 2)
                            } 
                            else if (((chromStart1[k] >= start_adj) 
                                        & (chromStart1[k] < start_adj1
                                        & chromEnd1[k] > start_adj1))) {
                                            distfromEnd <- start_adj1 - chromStart1[k]
                                            percentage <- round(abs((
                                                distfromEnd / (chromEnd1[k] - chromStart1[k])) * 100), digits = 2)
                            } 
							else if (((chromStart2[k] >= end_adj1) 
                                        & (chromStart2[k] < end_adj
                                        & chromEnd2[k] > end_adj))) {
                                            distfromEnd <- end_adj - chromStart2[k]
                                            percentage <- round(abs((
                                                distfromEnd / (chromEnd2[k] - chromStart2[k])) * 100), digits = 2)
                            }
                            else {
                                percentage <- "NA"
                            }
                            if (is.na(percentage)) {
                                geneInfo <- "-"
                            }
                            else {
                                geneInfo <- c(geneInfo, paste(
                                    gen1[k], "(", strnd1[k], ":",
                                    percentage, ")", sep = ""
                                ))
                            }
                        }
						
						for (k in seq_len(length(gen2)))
                        {
                            ##Checking if orientation of the gene and calculating 
                            ##length based on that
                            if(chromStart2[k]>=chromEnd2[k]){
                                lengthgene2 <- chromStart2[k] - chromEnd2[k]
                            }else{
                                lengthgene2 <- chromEnd2[k]-chromStart2[k]
                            }
                            #lengthbp1 <- start_adj1 - start_adj
							lengthbp2 <- end_adj - end_adj1
                            # print(paste('distfromStart:',distfromStart),sep='')
                            # print(paste('distfromEnd:',distfromEnd),sep='') Checking for
                            # partial/complete inclusion in the SV
                            if ((chromStart2[k] >= end_adj1 & chromEnd2[k] <= end_adj)) {
                                percentage <- 100
                            } 
                            else if ((chromStart2[k] < end_adj1
                                & chromEnd2[k] > end_adj) 
                                & (lengthgene > lengthbp2)) {
                                percentage <- round(abs(
								    ((end_adj - end_adj1) / (chromEnd2[k] - chromStart2[k])) * 100), digits = 2)
                            } 
                            else if (((chromStart2[k] < end_adj1) 
                                        & (chromEnd2[k] > end_adj1) 
                                        & (chromEnd2[k] <= end_adj))) {
                                            distfromStart <- chromEnd2[k] - end_adj1
                                            percentage <- round(abs((
                                            distfromStart / (chromEnd2[k] - chromStart2[k])) * 100), 
                                            digits = 2)
                            } 
                            else if (((chromStart2[k] >= end_adj1) 
                                        & (chromStart2[k] < end_adj
                                        & chromEnd2[k] > end_adj))) {
                                            distfromEnd <- end_adj - chromStart2[k]
                                            percentage <- round(abs((
                                                distfromEnd / (chromEnd2[k] - chromStart2[k])) * 100), digits = 2)
                            }
                            else {
                                percentage <- "NA"
                            }
                            if (is.na(percentage)) {
                                geneInfo <- "-"
                            }
                            else {
                                geneInfo <- c(geneInfo, paste(
                                    gen2[k], "(", strnd2[k], ":",
                                    percentage, ")", sep = ""
                                ))
                            }
                        }
                        geneInfo1 <- paste(geneInfo, collapse = ";")
                        #print(geneInfo1)
                        gnsInf <- c(gnsInf, str_squish(str_trim((geneInfo1),side="both")))
                        # print(paste("1:",ii,":",svid[ii],sep=""))
                        SVID <- c(SVID, svid[ii])
                    } else if (nrow(dat8) > 1 & nrow(dat9) > 1) {
                        ## Extracting strand and chromosome start information
						print("1st chromosome")
                        chromos1 <- dat9$Chromosome
                        chromStart1 <- dat9$Chromosome_Start
                        chromEnd1 <- dat9$Chromosome_End
                        gen1 <- dat9$Gene
                        strnd1 <- dat9$Strand
                        genMid1 <- dat9$geneMid
						print("2nd chromosome")
						chromos2 <- dat8$Chromosome
                        chromStart2 <- dat8$Chromosome_Start
                        chromEnd2 <- dat8$Chromosome_End
                        gen2 <- dat8$Gene
                        strnd2 <- dat8$Strand
                        genMid2 <- dat8$geneMid
                        geneInfo <- c()
						## Calculating the coverage of the gene
                        for (k in seq_len(length(gen1)))
                        {
                            ##Checking if orientation of the gene and calculating 
                            ##length based on that
                            if(chromStart1[k]>=chromEnd1[k]){
                                lengthgene1 <- chromStart1[k] - chromEnd1[k]
                            }else{
                                lengthgene1 <- chromEnd1[k]-chromStart1[k]
                            }
							
                            lengthbp1 <- start_adj1 - start_adj
							lengthbp2 <- end_adj - end_adj1
                            # print(paste('distfromStart:',distfromStart),sep='')
                            # print(paste('distfromEnd:',distfromEnd),sep='') Checking for
                            # partial/complete inclusion in the SV
                            if ((chromStart1[k] >= start_adj & chromEnd1[k] <= start_adj1)) {
                                percentage <- 100
                            } 
                            else if ((chromStart1[k] < start_adj 
                                & chromEnd1[k] > start_adj1) 
                                & (lengthgene1 > lengthbp1)) {
                                percentage <- round(abs(((start_adj1 - start_adj) / (chromEnd1[k]-chromStart1[k])) * 100), digits = 2)
                            } 
							
                            else if (((chromStart1[k] < start_adj) 
                                        & (chromEnd1[k] > start_adj) 
                                        & (chromEnd1[k] <= start_adj1))) {
                                            distfromStart <- chromEnd1[k] - start_adj
                                            percentage <- round(abs((
                                            distfromStart / (chromEnd1[k] - chromStart1[k])) * 100), 
                                            digits = 2)
                            } 
							else if (((chromStart2[k] < end_adj1) 
                                        & (chromEnd2[k] > end_adj1) 
                                        & (chromEnd2[k] <= end_adj))) {
                                            distfromStart <- chromEnd2[k] - end_adj1
                                            percentage <- round(abs((
                                            distfromStart / (chromEnd2[k] - chromStart2[k])) * 100), 
                                            digits = 2)
                            } 
                            else if (((chromStart1[k] >= start_adj) 
                                        & (chromStart1[k] < start_adj1
                                        & chromEnd1[k] > start_adj1))) {
                                            distfromEnd <- start_adj1 - chromStart1[k]
                                            percentage <- round(abs((
                                                distfromEnd / (chromEnd1[k] - chromStart1[k])) * 100), digits = 2)
                            } 
							else if (((chromStart2[k] >= end_adj1) 
                                        & (chromStart2[k] < end_adj
                                        & chromEnd2[k] > end_adj))) {
                                            distfromEnd <- end_adj - chromStart2[k]
                                            percentage <- round(abs((
                                                distfromEnd / (chromEnd2[k] - chromStart2[k])) * 100), digits = 2)
                            }
                            else {
                                percentage <- "NA"
                            }
                            if (is.na(percentage)) {
                                geneInfo <- "-"
                            }
                            else {
                                geneInfo <- c(geneInfo, paste(
                                    gen1[k], "(", strnd1[k], ":",
                                    percentage, ")", sep = ""
                                ))
                            }
                        }
						
						for (k in seq_len(length(gen2)))
                        {
                            ##Checking if orientation of the gene and calculating 
                            ##length based on that
                            if(chromStart2[k]>=chromEnd2[k]){
                                lengthgene2 <- chromStart2[k] - chromEnd2[k]
                            }else{
                                lengthgene2 <- chromEnd2[k]-chromStart2[k]
                            }
                            #lengthbp1 <- start_adj1 - start_adj
							lengthbp2 <- end_adj - end_adj1
                            # print(paste('distfromStart:',distfromStart),sep='')
                            # print(paste('distfromEnd:',distfromEnd),sep='') Checking for
                            # partial/complete inclusion in the SV
                            if ((chromStart2[k] >= end_adj1 & chromEnd2[k] <= end_adj)) {
                                percentage <- 100
                            } 
                            else if ((chromStart2[k] < end_adj1
                                & chromEnd2[k] > end_adj) 
                                & (lengthgene > lengthbp2)) {
                                percentage <- round(abs(
								    ((end_adj - end_adj1) / (chromEnd2[k] - chromStart2[k])) * 100), digits = 2)
                            } 
                            else if (((chromStart2[k] < end_adj1) 
                                        & (chromEnd2[k] > end_adj1) 
                                        & (chromEnd2[k] <= end_adj))) {
                                            distfromStart <- chromEnd2[k] - end_adj1
                                            percentage <- round(abs((
                                            distfromStart / (chromEnd2[k] - chromStart2[k])) * 100), 
                                            digits = 2)
                            } 
                            else if (((chromStart2[k] >= end_adj1) 
                                        & (chromStart2[k] < end_adj
                                        & chromEnd2[k] > end_adj))) {
                                            distfromEnd <- end_adj - chromStart2[k]
                                            percentage <- round(abs((
                                                distfromEnd / (chromEnd2[k] - chromStart2[k])) * 100), digits = 2)
                            }
                            else {
                                percentage <- "NA"
                            }
                            if (is.na(percentage)) {
                                geneInfo <- "-"
                            }
                            else {
                                geneInfo <- c(geneInfo, paste(
                                    gen2[k], "(", strnd2[k], ":",
                                    percentage, ")", sep = ""
                                ))
                            }
                        }
                        geneInfo1 <- paste(geneInfo, collapse = ";")
                        #print(geneInfo1)
                        gnsInf <- c(gnsInf, str_squish(str_trim((geneInfo1),side="both")))
                        # print(paste("1:",ii,":",svid[ii],sep=""))
                        SVID <- c(SVID, svid[ii])
                    } else if (nrow(dat8) >= 1 & nrow(dat9)== 0) {
                        ## Extracting strand and chromosome start information
						
						print("2nd chromosome")
						chromos2 <- dat8$Chromosome
                        chromStart2 <- dat8$Chromosome_Start
                        chromEnd2 <- dat8$Chromosome_End
                        gen2 <- dat8$Gene
                        strnd2 <- dat8$Strand
                        genMid2 <- dat8$geneMid
                        geneInfo <- c()
						
						for (k in seq_len(length(gen2)))
                        {
                            ##Checking if orientation of the gene and calculating 
                            ##length based on that
                            if(chromStart2[k]>=chromEnd2[k]){
                                lengthgene2 <- chromStart2[k] - chromEnd2[k]
                            }else{
                                lengthgene2 <- chromEnd2[k]-chromStart2[k]
                            }
                            #lengthbp1 <- start_adj1 - start_adj
							lengthbp2 <- end_adj - end_adj1
                            # print(paste('distfromStart:',distfromStart),sep='')
                            # print(paste('distfromEnd:',distfromEnd),sep='') Checking for
                            # partial/complete inclusion in the SV
                            if ((chromStart2[k] >= end_adj1 & chromEnd2[k] <= end_adj)) {
                                percentage <- 100
                            } 
                            else if ((chromStart2[k] < end_adj1
                                & chromEnd2[k] > end_adj) 
                                & (lengthgene > lengthbp2)) {
                                percentage <- round(abs(
								    ((end_adj - end_adj1) / (chromEnd2[k] - chromStart2[k])) * 100), digits = 2)
                            } 
                            else if (((chromStart2[k] < end_adj1) 
                                        & (chromEnd2[k] > end_adj1) 
                                        & (chromEnd2[k] <= end_adj))) {
                                            distfromStart <- chromEnd2[k] - end_adj1
                                            percentage <- round(abs((
                                            distfromStart / (chromEnd2[k] - chromStart2[k])) * 100), 
                                            digits = 2)
                            } 
                            else if (((chromStart2[k] >= end_adj1) 
                                        & (chromStart2[k] < end_adj
                                        & chromEnd2[k] > end_adj))) {
                                            distfromEnd <- end_adj - chromStart2[k]
                                            percentage <- round(abs((
                                                distfromEnd / (chromEnd2[k] - chromStart2[k])) * 100), digits = 2)
                            }
                            else {
                                percentage <- "NA"
                            }
                            if (is.na(percentage)) {
                                geneInfo <- "-"
                            }
                            else {
                                geneInfo <- c(geneInfo, paste(
                                    gen2[k], "(", strnd2[k], ":",
                                    percentage, ")", sep = ""
                                ))
                            }
                        }
                        geneInfo1 <- paste(geneInfo, collapse = ";")
                        #print(geneInfo1)
                        gnsInf <- c(gnsInf, str_squish(str_trim((geneInfo1),side="both")))
                        # print(paste("1:",ii,":",svid[ii],sep=""))
                        SVID <- c(SVID, svid[ii])
                    }else if (nrow(dat8) == 0 & nrow(dat9) >= 1) {
                        print("1st chromosome")
                        chromos1 <- dat9$Chromosome
                        chromStart1 <- dat9$Chromosome_Start
                        chromEnd1 <- dat9$Chromosome_End
                        gen1 <- dat9$Gene
                        strnd1 <- dat9$Strand
                        genMid1 <- dat9$geneMid
						
						## Calculating the coverage of the gene
                        for (k in seq_len(length(gen1)))
                        {
                            ##Checking if orientation of the gene and calculating 
                            ##length based on that
                            if(chromStart1[k]>=chromEnd1[k]){
                                lengthgene1 <- chromStart1[k] - chromEnd1[k]
                            }else{
                                lengthgene1 <- chromEnd1[k]-chromStart1[k]
                            }
							
                            lengthbp1 <- start_adj1 - start_adj
							lengthbp2 <- end_adj - end_adj1
                            # print(paste('distfromStart:',distfromStart),sep='')
                            # print(paste('distfromEnd:',distfromEnd),sep='') Checking for
                            # partial/complete inclusion in the SV
                            if ((chromStart1[k] >= start_adj & chromEnd1[k] <= start_adj1)) {
                                percentage <- 100
                            } 
                            else if ((chromStart1[k] < start_adj 
                                & chromEnd1[k] > start_adj1) 
                                & (lengthgene1 > lengthbp1)) {
                                percentage <- round(abs(((start_adj1 - start_adj) / (chromEnd1[k]-chromStart1[k])) * 100), digits = 2)
                            } 
							
                            else if (((chromStart1[k] < start_adj) 
                                        & (chromEnd1[k] > start_adj) 
                                        & (chromEnd1[k] <= start_adj1))) {
                                            distfromStart <- chromEnd1[k] - start_adj
                                            percentage <- round(abs((
                                            distfromStart / (chromEnd1[k] - chromStart1[k])) * 100), 
                                            digits = 2)
                            } 
							else if (((chromStart2[k] < end_adj1) 
                                        & (chromEnd2[k] > end_adj1) 
                                        & (chromEnd2[k] <= end_adj))) {
                                            distfromStart <- chromEnd2[k] - end_adj1
                                            percentage <- round(abs((
                                            distfromStart / (chromEnd2[k] - chromStart2[k])) * 100), 
                                            digits = 2)
                            } 
                            else if (((chromStart1[k] >= start_adj) 
                                        & (chromStart1[k] < start_adj1
                                        & chromEnd1[k] > start_adj1))) {
                                            distfromEnd <- start_adj1 - chromStart1[k]
                                            percentage <- round(abs((
                                                distfromEnd / (chromEnd1[k] - chromStart1[k])) * 100), digits = 2)
                            } 
							else if (((chromStart2[k] >= end_adj1) 
                                        & (chromStart2[k] < end_adj
                                        & chromEnd2[k] > end_adj))) {
                                            distfromEnd <- end_adj - chromStart2[k]
                                            percentage <- round(abs((
                                                distfromEnd / (chromEnd2[k] - chromStart2[k])) * 100), digits = 2)
                            }
                            else {
                                percentage <- "NA"
                            }
                            if (is.na(percentage)) {
                                geneInfo <- "-"
                            }
                            else {
                                geneInfo <- c(geneInfo, paste(
                                    gen1[k], "(", strnd1[k], ":",
                                    percentage, ")", sep = ""
                                ))
                            }
                        }
                        geneInfo1 <- paste(geneInfo, collapse = ";")
                        #print(geneInfo1)
                        gnsInf <- c(gnsInf, str_squish(str_trim((geneInfo1),side="both")))
                        # print(paste("1:",ii,":",svid[ii],sep=""))
                        SVID <- c(SVID, svid[ii])
                    }else {
                        SVID <- c(SVID, svid[ii])
                        gnsInf <- c(gnsInf, "-")
                        # print(paste("1:",ii,":",svid[ii],sep=""))
                    }
					
					
			    }else { SVID <- c(SVID, svid[ii])
                        gnsInf <- c(gnsInf, "-")}
					
			}
				
					
					
            ## Extracting strand and calculating percentage coverage information for
            ## a gene covering the breakpoints if multiple genes are found in the
            ## breakpoint region.
            
        }       
	
    ## Writing a returning a data frame with gene information and SVID
    dat3 <- data.frame(SVID, OverlapGenes_strand_perc = gnsInf)
    return(dat3)
}
#' Calculates Genes that are near to the SV region
#'
#' @param bed    Text Bionano Bed file.
#' @param chrom    character SVmap chromosome.
#' @param chrom2   character SVmap 2nd chromosome.
#' @param startpos    numeric starting position of the breakpoints.
#' @param endpos    numeric end position of the breakpoints.
#' @param svid    numeric Structural variant identifier (Bionano generated).
#' @param n    numeric Number of genes to report which are nearest to the breakpoint.
#' Default is 3.
#' @param bperrorindel Numeric. base pair error indel.
#' @param bperrorinvtrans Numeric. base pair error invtranslocation.
#' @return Data Frame. Contains the SVID,Gene name,strand information and
#' Distance from the SV covered.
#' @examples
#' smapName="F1.1_GM24385_DLE-1_P_trio_hg19.smap"
#' smap = system.file("extdata", smapName, package="nanotatoR")
#' bedFile <- system.file("extdata", "Homo_sapiens.Hg19.bed", package="nanotatoR")
#' bed<-buildrunBNBedFiles(bedFile,returnMethod="dataFrame")
#' smap<-readSMap(smap)
#' chrom<-smap$RefcontigID1
#' startpos<-smap$RefStartPos
#' endpos<-smap$RefEndPos
#' if (length(grep("SVIndex",names(smap)))>0){
#'        svid <- smap$SVIndex
#'    }else{
#'     svid <- smap$SmapEntryID
#'     }
#' n<-3
#' nonOverlapGenes(bed, chrom, startpos, endpos, svid,n)
#' @import utils
#' @export
nonOverlapGenes <- function(bed, chrom, startpos, 
    chrom2, endpos, svid, n = 3, SVTyp,
    bperrorindel = 3000, bperrorinvtrans = 10000) {
    ## Initializing data
    print("***NonOverlap Genes***")
    data1 <- data.frame()
    gnsInf_UP <- c()
    gnsInf_DN <- c()
    SVID <- c()
	SVTyp = SVTyp
    ## Extracting genes near the breakpoints and calculating distance from
    ## it.
	chrom <- chrom
	chrom2 <- chrom2
    for (ii in seq_len(length(chrom))){   
	    if(startpos[ii] == -1| startpos[ii] == -1){
	        SVID <- c(SVID, svid[ii])
	        gnsInf_UP <- c(gnsInf_UP, "-")
			gnsInf_DN <- c(gnsInf_DN, "-")
		}
		else{
		    if(SVTyp[ii] == "insertion" 
			    | SVTyp[ii] == "deletion"
				| SVTyp[ii] == "duplication"
				| SVTyp[ii] == "duplication_split"
				| SVTyp[ii] == "duplication_inverted"){
				start_adj <- startpos[ii] - bperrorindel
				end_adj <- endpos[ii] + bperrorindel
				dat12 <- bed[which(as.character(bed$Chromosome) == chrom[ii]), ]
                datup <- dat12[which((dat12$Chromosome_Start < start_adj &
                dat12$Chromosome_End < start_adj) &    (dat12$Strand=="-")), ]
		    
                # datup_minus <- dat12[which(as.character(bed$Chromosome) == chrom[ii]
                # & ((bed$Chromosome_End < start_adj & bed$Chromosome_Start <
                # end_adj)) & bed$Strand == '-'), ]
                datdn <- dat12[which((dat12$Chromosome_Start > end_adj & dat12$Chromosome_End >
                    end_adj) & (dat12$Strand=="+")), ]
                # datdn_minus <- bed[which(as.character(bed$Chromosome) == chrom[ii] &
                # ((bed$Chromosome_End > end_adj & bed$Chromosome_Start >
                # start_adj)) & bed$Strand == '-'), ] print(datup);print(datdn)
                if (nrow(datup) > 0 & nrow(datdn) > 0) {
                    ## Upstream Genes
                #print("1")
                    datup$diff_up <- abs(round((start_adj - datup$Chromosome_Start)/1000,digits=3))
                    
                    dat_up <- datup[order(datup$diff_up), ]
                    dat_up <- dat_up[which(dat_up$diff_up > 0), ]
                    genes_up <- dat_up$Gene[1:n]
                    strands_up <- dat_up$Strand[1:n]
                    diff_up <- dat_up$diff_up[1:n]
                    genesUP <- c()
                    ### Storing the gene, strand and distance information in vectors
                    for (k in 1:n)
                    {
                        pas <- paste(
                            genes_up[k], "(", strands_up[k], ":", diff_up[k],
                            ")", sep = ""
                        )
                        genesUP <- c(genesUP, pas)
                    }
                    genesUP1 <- paste(genesUP, collapse = ";")
                    #pri    nt(genesUP1))
                    gnsInf_UP <- c(gnsInf_UP, str_squish(str_trim((genesUP1),side="both")))
                    # Downstram Genes
                    datdn$diff_dn<-abs(round((
                            end_adj - datdn$Chromosome_Start)/1000,digits=3)
                            )
                    
                    dat_dn <- datdn[order(datdn$diff_dn), ]
                    dat_dn <- dat_dn[which(dat_dn$diff_dn > 0), ]
                    genes_dn <- dat_dn$Gene[1:n]
                    strands_dn <- dat_dn$Strand[1:n]
                    diff_dn <- dat_dn$diff_dn[1:n]
                    genesDN <- c()
                    ### Storing the gene, strand and distance information in vectors
                    for (k in 1:n)
                    {
                        pas <- paste(
                            genes_dn[k], "(", strands_dn[k], ":", diff_dn[k],
                            ")", sep = ""
                        )
                        genesDN <- c(genesDN, pas)
                    }
                    genesDN1 <- paste(genesDN, collapse = ";")
                    #pri    nt(genesDN1)
                    gnsInf_DN <- c(
                            gnsInf_DN,str_squish(str_trim((genesDN1), side="both"))
                            )
                }       
                else if (nrow(datup) == 0 & nrow(datdn) > 0) {
                        gnsInf_UP <- c(gnsInf_UP, "-")
                    #print("2")
                     # Downstram Genes
                        datdn$diff_dn<-abs(round((
                                end_adj - datdn$Chromosome_Start)/1000,digits=3)
                                )
                        'diffStartGene_startbp <- abs(round((start_adj - datdn$Chromosome_Start)/1000,digits=3))
                            diffEndGene_startbp <-abs(round((end_adj - datdn$Chromosome_Start)/1000,digits=3))
                            diffStartGene_endbp <- abs(round((start_adj - datdn$Chromosome_End)/1000,digits=3))
                            diffEndGene_endbp <- abs(round((end_adj - datdn$Chromosome_End)/1000,digits=3))'
                        'for (kk in 1:length(diffStartGene_startbp))
                        {
                            
                            if ((datdn$Chromosome_Start[kk] > start_adj & datdn$Chromosome_End[kk] >
                                end_adj) & (diffStartGene_startbp[kk] < 0 & diffEndGene_startbp[kk] <
                                0) & (diffEndGene_startbp[kk] > diffEndGene_endbp[kk])) {
                                datdn$diff_dn[kk] <- as.numeric(diffEndGene_startbp[kk])
                            } else if ((datdn$Chromosome_Start[kk] > start_adj & datdn$Chromosome_End[kk] >
                                end_adj) & (diffStartGene_startbp < 0 & diffEndGene_startbp <
                                0) & (diffEndGene_startbp < diffEndGene_endbp)) {
                                datdn$diff_dn <- as.numeric(diffEndGene_endbp[kk])
                            } else {
                                datdn$diff_dn[kk] <- 0
                            }
                        }'
                        dat_dn <- datdn[order(datdn$diff_dn), ]
                        dat_dn <- dat_dn[which(dat_dn$diff_dn > 0), ]
                        genes_dn <- dat_dn$Gene[1:n]
                        strands_dn <- dat_dn$Strand[1:n]
                        diff_dn <- dat_dn$diff_dn[1:n]
                        genesDN <- c()
                        ### Storing the gene, strand and distance information in vectors
                        for (k in 1:n)
                        {
                            pas <- paste(
                                genes_dn[k], "(", strands_dn[k], ":", diff_dn[k],
                                ")", sep = ""
                            )
                            genesDN <- c(genesDN, pas)
                        }
                        genesDN1 <- paste(genesDN, collapse = ";")
                        #print(genesDN1)
                        gnsInf_DN <- c(
                            gnsInf_DN, str_squish(str_trim((genesDN1),side = "both"))
                            )
                        
                } 
                else if (nrow(datup) > 0 & nrow(datdn) == 0) {
                    #print("3")
                        
                    ## Upstream Genes
                    #print("1")
                    datup$diff_up <- abs(
                             round((start_adj - datup$Chromosome_Start)/1000,digits=3
                            ))
                        dat_up <- datup[order(datup$diff_up), ]
                        dat_up <- dat_up[which(dat_up$diff_up > 0), ]
                        genes_up <- dat_up$Gene[1:n]
                        strands_up <- dat_up$Strand[1:n]
                        diff_up <- dat_up$diff_up[1:n]
                        genesUP <- c()
                        ### Storing the gene, strand and distance information in vectors
                        for (k in 1:n)
                        {
                            pas <- paste(
                                genes_up[k], "(", strands_up[k], ":", diff_up[k],
                                ")", sep = ""
                            )
                            genesUP <- c(genesUP, pas)
                        }
                        genesUP1 <- paste(genesUP, collapse = ";")
                        #print(genesUP1)
                        gnsInf_UP <- c(gnsInf_UP, str_squish(str_trim((genesUP1),side="both")))
                        gnsInf_DN <- c(gnsInf_DN, "-")
                    } 
                else{
                        gnsInf_UP <- c(gnsInf_UP, "-")
                        gnsInf_DN <- c(gnsInf_DN, "-")
                    }   
			    
			    
                SVID <- c(SVID, svid[ii])
                
		    }else if((length(grep("inversion", SVTyp[ii])) >= 1)){
		        start_adj <- startpos[ii] - bperrorinvtrans
		        start_adj1 <- startpos[ii] + bperrorinvtrans
		        end_adj <- endpos[ii] + bperrorinvtrans
		        end_adj1 <- endpos[ii] - bperrorinvtrans
				dat12 <- bed[which(as.character(bed$Chromosome) == chrom[ii]), ]
                    datup <- dat12[which((dat12$Chromosome_Start < start_adj &
                    dat12$Chromosome_End < start_adj) &    (dat12$Strand=="-")), ]
		        
                # datup_minus <- dat12[which(as.character(bed$Chromosome) == chrom[ii]
                # & ((bed$Chromosome_End < start_adj & bed$Chromosome_Start <
                # end_adj)) & bed$Strand == '-'), ]
                datdn <- dat12[which((dat12$Chromosome_Start > end_adj & dat12$Chromosome_End >
                    end_adj) & (dat12$Strand=="+")), ]
                # datdn_minus <- bed[which(as.character(bed$Chromosome) == chrom[ii] &
                # ((bed$Chromosome_End > end_adj & bed$Chromosome_Start >
                # start_adj)) & bed$Strand == '-'), ] print(datup);print(datdn)
                if (nrow(datup) > 0 & nrow(datdn) > 0) {
                    ## Upstream Genes
                #print("1")
                    datup$diff_up <- abs(round((start_adj - datup$Chromosome_Start)/1000,digits=3))
                    
                    dat_up <- datup[order(datup$diff_up), ]
                    dat_up <- dat_up[which(dat_up$diff_up > 0), ]
                    genes_up <- dat_up$Gene[1:n]
                    strands_up <- dat_up$Strand[1:n]
                    diff_up <- dat_up$diff_up[1:n]
                    genesUP <- c()
                    ### Storing the gene, strand and distance information in vectors
                    for (k in 1:n)
                    {
                        pas <- paste(
                            genes_up[k], "(", strands_up[k], ":", diff_up[k],
                            ")", sep = ""
                        )
                        genesUP <- c(genesUP, pas)
                    }
                    genesUP1 <- paste(genesUP, collapse = ";")
                #print(genesUP1))
                    gnsInf_UP <- c(gnsInf_UP, str_squish(str_trim((genesUP1),side="both")))
                    # Downstram Genes
                    datdn$diff_dn<-abs(round((
                            end_adj - datdn$Chromosome_Start)/1000,digits=3)
                            )
                    
                    dat_dn <- datdn[order(datdn$diff_dn), ]
                    dat_dn <- dat_dn[which(dat_dn$diff_dn > 0), ]
                    genes_dn <- dat_dn$Gene[1:n]
                    strands_dn <- dat_dn$Strand[1:n]
                    diff_dn <- dat_dn$diff_dn[1:n]
                    genesDN <- c()
                    ### Storing the gene, strand and distance information in vectors
                    for (k in 1:n)
                    {
                        pas <- paste(
                            genes_dn[k], "(", strands_dn[k], ":", diff_dn[k],
                            ")", sep = ""
                        )
                        genesDN <- c(genesDN, pas)
                    }
                    genesDN1 <- paste(genesDN, collapse = ";")
                #print(genesDN1)
                    gnsInf_DN <- c(
                            gnsInf_DN,str_squish(str_trim((genesDN1), side="both"))
                            )
                } 
                else if (nrow(datup) == 0 & nrow(datdn) > 0) {
                        gnsInf_UP <- c(gnsInf_UP, "-")
                    #print("2")
                    # Downstram Genes
                    datdn$diff_dn<-abs(round((
                            end_adj - datdn$Chromosome_Start)/1000,digits=3)
                            )
                    
                    dat_dn <- datdn[order(datdn$diff_dn), ]
                    dat_dn <- dat_dn[which(dat_dn$diff_dn > 0), ]
                    genes_dn <- dat_dn$Gene[1:n]
                    strands_dn <- dat_dn$Strand[1:n]
                    diff_dn <- dat_dn$diff_dn[1:n]
                    genesDN <- c()
                    ### Storing the gene, strand and distance information in vectors
                    for (k in 1:n)
                    {
                        pas <- paste(
                            genes_dn[k], "(", strands_dn[k], ":", diff_dn[k],
                            ")", sep = ""
                        )
                        genesDN <- c(genesDN, pas)
                    }
                    genesDN1 <- paste(genesDN, collapse = ";")
                    #print(genesDN1)
                    gnsInf_DN <- c(
                        gnsInf_DN, str_squish(str_trim((genesDN1),side = "both"))
                        )
                    
                } 
                else if (nrow(datup) > 0 & nrow(datdn) == 0) {
                #print("3")
                   
                ## Upstream Genes
                #print("1")
                    datup$diff_up <- abs(
                         round((start_adj - datup$Chromosome_Start)/1000,digits=3
                        ))
                    
                    dat_up <- datup[order(datup$diff_up), ]
                    dat_up <- dat_up[which(dat_up$diff_up > 0), ]
                    genes_up <- dat_up$Gene[1:n]
                    strands_up <- dat_up$Strand[1:n]
                    diff_up <- dat_up$diff_up[1:n]
                    genesUP <- c()
                    ### Storing the gene, strand and distance information in vectors
                    for (k in 1:n)
                    {
                        pas <- paste(
                            genes_up[k], "(", strands_up[k], ":", diff_up[k],
                            ")", sep = ""
                        )
                        genesUP <- c(genesUP, pas)
                    }
                    genesUP1 <- paste(genesUP, collapse = ";")
                    #print(genesUP1)
                    gnsInf_UP <- c(gnsInf_UP, str_squish(str_trim((genesUP1),side="both")))
                    gnsInf_DN <- c(gnsInf_DN, "-")
                } 
                else {
                    gnsInf_UP <- c(gnsInf_UP, "-")
                    gnsInf_DN <- c(gnsInf_DN, "-")
                }
	            
                    ## Upstream genes
	            
                    ## Downstream Genes
	            
                    ## Getting the genes and strand data
	            
                    ### Genes_Up_minusEnd Getting Down streamplus genes
	            
                    ## Getting Genes minus Downstream
                    SVID <- c(SVID, svid[ii])
		    }else if((length(grep("translocation", SVTyp[ii])) >= 1)){
		        		
		        start_adj1 <- startpos[ii] + bperrorinvtrans
		        end_adj <- endpos[ii] + bperrorinvtrans
		        end_adj1 <- endpos[ii] - bperrorinvtrans
				dat12 <- bed[which(as.character(bed$Chromosome) == chrom[ii]), ]
                datup <- dat12[which((dat12$Chromosome_Start < start_adj &
                dat12$Chromosome_End < start_adj) &    (dat12$Strand=="-")), ]
		        
                # datup_minus <- dat12[which(as.character(bed$Chromosome) == chrom[ii]
                # & ((bed$Chromosome_End < start_adj & bed$Chromosome_Start <
                # end_adj)) & bed$Strand == '-'), ]
				dat13 <- bed[which(as.character(bed$Chromosome) == chrom2[ii]), ]
                datdn <- dat13[which((dat13$Chromosome_Start > end_adj & dat13$Chromosome_End >
                    end_adj) & (dat13$Strand=="+")), ]
                # datdn_minus <- bed[which(as.character(bed$Chromosome) == chrom[ii] &
                # ((bed$Chromosome_End > end_adj & bed$Chromosome_Start >
                # start_adj)) & bed$Strand == '-'), ] print(datup);print(datdn)
                if (nrow(datup) > 0 & nrow(datdn) > 0) {
                    ## Upstream Genes
                #print("1")
                    datup$diff_up <- abs(round((start_adj - datup$Chromosome_Start)/1000,digits=3))
                    
                    dat_up <- datup[order(datup$diff_up), ]
                    dat_up <- dat_up[which(dat_up$diff_up > 0), ]
                    genes_up <- dat_up$Gene[1:n]
                    strands_up <- dat_up$Strand[1:n]
                    diff_up <- dat_up$diff_up[1:n]
                    genesUP <- c()
                    ### Storing the gene, strand and distance information in vectors
                    for (k in 1:n)
                    {
                        pas <- paste(
                            genes_up[k], "(", strands_up[k], ":", diff_up[k],
                            ")", sep = ""
                        )
                        genesUP <- c(genesUP, pas)
                    }
                    genesUP1 <- paste(genesUP, collapse = ";")
                #print(genesUP1))
                    gnsInf_UP <- c(gnsInf_UP, str_squish(str_trim((genesUP1),side="both")))
                    # Downstram Genes
                    datdn$diff_dn<-abs(round((
                            end_adj - datdn$Chromosome_Start)/1000,digits=3)
                            )
                    
                    dat_dn <- datdn[order(datdn$diff_dn), ]
                    dat_dn <- dat_dn[which(dat_dn$diff_dn > 0), ]
                    genes_dn <- dat_dn$Gene[1:n]
                    strands_dn <- dat_dn$Strand[1:n]
                    diff_dn <- dat_dn$diff_dn[1:n]
                    genesDN <- c()
                    ### Storing the gene, strand and distance information in vectors
                    for (k in 1:n)
                    {
                        pas <- paste(
                            genes_dn[k], "(", strands_dn[k], ":", diff_dn[k],
                            ")", sep = ""
                        )
                        genesDN <- c(genesDN, pas)
                    }
                    genesDN1 <- paste(genesDN, collapse = ";")
                #print(genesDN1)
                    gnsInf_DN <- c(
                            gnsInf_DN,str_squish(str_trim((genesDN1), side="both"))
                            )
                } 
                else if (nrow(datup) == 0 & nrow(datdn) > 0) {
                        gnsInf_UP <- c(gnsInf_UP, "-")
                    #print("2")
                    # Downstram Genes
                    datdn$diff_dn<-abs(round((
                            end_adj - datdn$Chromosome_Start)/1000,digits=3)
                            )
                    
                    dat_dn <- datdn[order(datdn$diff_dn), ]
                    dat_dn <- dat_dn[which(dat_dn$diff_dn > 0), ]
                    genes_dn <- dat_dn$Gene[1:n]
                    strands_dn <- dat_dn$Strand[1:n]
                    diff_dn <- dat_dn$diff_dn[1:n]
                    genesDN <- c()
                    ### Storing the gene, strand and distance information in vectors
                    for (k in 1:n)
                    {
                        pas <- paste(
                            genes_dn[k], "(", strands_dn[k], ":", diff_dn[k],
                            ")", sep = ""
                        )
                        genesDN <- c(genesDN, pas)
                    }
                    genesDN1 <- paste(genesDN, collapse = ";")
                    #print(genesDN1)
                    gnsInf_DN <- c(
                        gnsInf_DN, str_squish(str_trim((genesDN1),side = "both"))
                        )
                    
                } 
                else if (nrow(datup) > 0 & nrow(datdn) == 0) {
                #print("3")
                   
                ## Upstream Genes
                #print("1")
                    datup$diff_up <- abs(
                         round((start_adj - datup$Chromosome_Start)/1000,digits=3
                        ))
                    
                    dat_up <- datup[order(datup$diff_up), ]
                    dat_up <- dat_up[which(dat_up$diff_up > 0), ]
                    genes_up <- dat_up$Gene[1:n]
                    strands_up <- dat_up$Strand[1:n]
                    diff_up <- dat_up$diff_up[1:n]
                    genesUP <- c()
                    ### Storing the gene, strand and distance information in vectors
                    for (k in 1:n)
                    {
                        pas <- paste(
                            genes_up[k], "(", strands_up[k], ":", diff_up[k],
                            ")", sep = ""
                        )
                        genesUP <- c(genesUP, pas)
                    }
                    genesUP1 <- paste(genesUP, collapse = ";")
                    #print(genesUP1)
                    gnsInf_UP <- c(gnsInf_UP, str_squish(str_trim((genesUP1),side="both")))
                    gnsInf_DN <- c(gnsInf_DN, "-")
                } 
                else {
                    gnsInf_UP <- c(gnsInf_UP, "-")
                    gnsInf_DN <- c(gnsInf_DN, "-")
                }
	            
                    ## Upstream genes
	            
                    ## Downstream Genes
	            
                    ## Getting the genes and strand data
	            
                    ### Genes_Up_minusEnd Getting Down streamplus genes
	            
                    ## Getting Genes minus Downstream
                    SVID <- c(SVID, svid[ii])
		    }else{
                SVID <- c(SVID, svid[ii]) 
				gnsInf_UP <- c(gnsInf_UP, "-")
                gnsInf_DN <- c(gnsInf_DN, "-")
				                    
            }
	    }
	}
    ## Writing and storing the data
    dat4 <- data.frame(
        SVID, Upstream_nonOverlapGenes_dist_kb = gnsInf_UP,
        Downstream_nonOverlapGenes_dist_kb = gnsInf_DN
    )
    return(dat4)
}
#' Extracts gene information from bed files
#'
#' @param smap    character. Path to SMAP file.
#' @param bed    Text. Normal Bed files or Bionano Bed file.
#' @param inputfmtBed character Whether the bed input is UCSC bed or Bionano bed.
#' Note: extract in bed format to be read by bedsv:
#' awk '{print $1,$4,$5,$18,$7}' gencode.v19.annotation.gtf >HomoSapienGRCH19.bed
#' @param outpath    character Path for the output files.
#' @param n    numeric Number of genes to report which are nearest to the
#' breakpoint. Default is        3.
#' @param returnMethod_bedcomp Character. Type of output Dataframe or in Text format.
#' @param EnzymeType Character. Type of enzyme. Options Dual and DLE.
#' @param bperrorindel Numeric. base pair error indel.
#' @param bperrorinvtrans Numeric. base pair error invtranslocation.
#' @return Data Frame and Text file. Contains the smap with additional Gene Information.
#' @examples
#' smapName="F1.1_GM24385_DLE-1_P_trio_hg19.smap"
#' smap = system.file("extdata", smapName, package="nanotatoR")
#' bedFile <- system.file("extdata", "Homo_sapiens.Hg19.bed", package="nanotatoR")
#' outpath <- system.file("extdata",    package="nanotatoR")
#' datcomp<-compSmapbed(smap, bed=bedFile, inputfmtBed =    "BED", outpath,
#' n = 3, returnMethod_bedcomp = c("dataFrame"), EnzymeType = "Dual")
#' datcomp[1,]
#' @import utils
#' @export
compSmapbed <- function(smap, bed, inputfmtBed = c("BED", "BNBED"), 
                        EnzymeType = c("Dual", "DLE"), outpath,
                        n = 3, returnMethod_bedcomp = c("Text", "dataFrame"), 
                        input_fmt_smap = c("Text","dataFrame"),
						bperrorindel = 3000, 
						bperrorinvtrans = 10000) {
    print("###Comparing SVs and Beds###")
    ## Checking if bed file is from UCSC or BN bed files.
    if (inputfmtBed == "BNBED") {
        r1 <- readBNBedFiles(BNFile = bed)
    } else if (inputfmtBed == "BED") {
        r1 <- buildrunBNBedFiles(bed, returnMethod = "dataFrame")
    } else {
        stop("Incorrect File format")
    }
    ## reading SMAPs and extracting data
    if(EnzymeType == "Dual"){
        r2 <- readSMap(smap, input_fmt_smap)
    }
    else{
        r2 <- readSMap_DLE(smap, input_fmt_smap)
    }
    chrom <- r2$RefcontigID1
	chrom2 <- r2$RefcontigID2
    startpos <- r2$RefStartPos
    endpos <- r2$RefEndPos
    if (length(grep("SVIndex", names(r2))) > 0) {
        svid <- r2$SVIndex
    } else {
        svid <- r2$SmapEntryID
    }
    SVTyp <- r2$Type
    ## Calls for the functions to extract overlap/non-overlap genes and
    ## calculate informations for the same.
    dat3 <- overlapGenes(r1, chrom, startpos, endpos, svid, chrom2, SVTyp = SVTyp, 
	    bperrorindel = bperrorindel, bperrorinvtrans = bperrorinvtrans)
    # print(names(dat3))
    dat4 <- nonOverlapGenes(r1, chrom, startpos, endpos, svid, chrom2, SVTyp = SVTyp,
	     bperrorindel = bperrorindel, bperrorinvtrans = bperrorinvtrans)
    # print(names(dat4))
    ## Writing results in data frame and text files
    data1 <- data.frame(r2, 
            OverlapGenes_strand_perc = dat3$OverlapGenes_strand_perc,
            dat4[, 2:ncol(dat4)])
    # print(names(data1))
    if (returnMethod_bedcomp == "Text") {
        st1 <- strsplit(smap, split = "\\\\")
        st2 <- strsplit(st1[[1]][4], split = ".txt")
        fname <- paste(st2[[1]][1], "_Final_SVMAP.txt", sep = "")
        write.table(data1, file.path(outpath, fname))
    } else if (returnMethod_bedcomp == "dataFrame") {
        return(data1)
    } else {
        stop("Method of return incorrect")
    }
}
