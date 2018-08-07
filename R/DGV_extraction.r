#' Frequency calculation of variants compared to DGV.
#'
#' @param hgpath  character. Path to Database of Genomic Variants (DGV)
#'                Text file. 
#' @param smappath  character. Path and file name for textfile.
#' @param terms  character. Single or Multiple Terms.
#' @param outpath character. Path where gene lists are saved.
#' @param input_fmt character. Choice between text or data frame as 
#' an input to the DGV frequency calculator.
#' @param smap_data Dataset containing smap data.
#' @param thresh integer. Threshold for the number of terms sent to entrez.
#'                Note if large lists are sent to ncbi, it might fail to get
#'                processed. Default is 5.
#' @param returnMethod character. Choice between text or data frame as the output.
#' @return Text and character vector containg gene list and terms associated with them
#'         are stored as text files.
#' @examples
#' \dontrun{
#' hgpath="C:\\Annotator\\Data\\GRCh37_hg19_variants_2016-05-15.txt";
#' smappath="Z:/Hayk's_Materials/Bionano/Projects/UDN/F1_UDN287643_Benic.Aria/";
#' smap="F1_UDN287643_P_Q.S_VAP_SVmerge_trio_access.txt";
#' win_indel=10000;win_inv_trans=50000;perc_similarity=0.5
#' DGV_extraction (hgpath, smappath, win_indel = 10000, win_inv_trans = 50000, 
#' perc_similarity = 0.5,returnMethod="dataFrame")
#' }
#' @export
DGV_extraction <- function(hgpath, smappath, smap_data,input_fmt=c("Text","DataFrame"),
    win_indel = 10000, win_inv_trans = 50000, perc_similarity = 0.5,
	returnMethod=c("Text","dataFrame"))
    {
    # S='F' Change the window for Inversion/translocation 50000
    print("###Calculating DGV Frequency###")
    ## Reading the DGV file and the smap file
    r <- read.table(paste(hgpath), header = TRUE, sep = "\t")
    ## Unique variant length extraction for frquency Calculation
    samp <- as.character(unique(r$samples))
    #usamp <- UniqueSample(samp)
    usamp <- 7965
    # varaccl<-length(unique(varacc)) close(con) 
	##Checking if the input format is dataframe or Text
	if(input_fmt=="Text"){
	##Pattern matching needs to be done to remove #
	    con <- file(paste(smappath, smap, sep = ""), "r")
    r10 <- readLines(con, n = -1)
    close(con)
    datfinal <- data.frame()
    g1 <- grep("RawConfidence", r10)
    g2 <- grep("RefStartPos", r10)
    if (g1 == g2)
    {
        dat <- gsub("# ", "", r10)
        # dat<-gsub('\t',' ',r1)
        dat4 <- textConnection(dat[g1:length(dat)])
        r1 <- read.table(dat4, sep = "\t", header = TRUE)
    } else
    {
        stop("column names doesnot Match")
    }
	}
	else if (input_fmt=="DataFrame"){
	r1<-smap_data
	}
	else{
	stop("Incorrect format")
	}
    ## Checking Sex male/female and assigning chromosome number accordingly
    chro <- length(unique(r1$RefcontigID1))
    chro1 <- c(1:chro)
    dataFinal <- c()
    ## Extracting Data for 1 chromosome at a time and comparing Make change
    ## in XY
    for (ii in 1:length(chro1))
    {
        #print(paste("Chromosome:", chro1[ii]))
        ## Extracting data from DGV dataset
        if (ii == 23)
        {
            kk <- "X"
        } else if (ii == 24)
        {
            kk <- "Y"
        } else
        {
            kk <- chro1[ii]
        }
        dat <- r[which(r$chr == kk), ]
        # variantType1<-dat$variantsubtype Changing the variant terms in DGV to
        # match svmap
        dat$variantsubtype <- gsub("loss", "deletion", dat$variantsubtype)
        dat$variantsubtype <- gsub("gain", "insertion", dat$variantsubtype)
        dat$variantsubtype <- gsub("gain+loss", "insertion+deletion", dat$variantsubtype)
        variantType1 <- dat$variantsubtype
        ## Extracting data from SVmap
        dat1 <- r1[which(r1$RefcontigID1 == chro1[ii]), ]
        ## Adding the windows to the breakpoints
        rf <- dat1$RefStartPos
        rf_wb_ind <- rf - win_indel
        rf_fb_ind <- rf + win_indel
        rf_wb_int <- rf - win_inv_trans
        rf_fb_int <- rf + win_inv_trans
        re <- dat1$RefEndPos
        re_wf_ind <- re + win_indel
        re_wb_ind <- re - win_indel
        re_wb_int <- re - win_inv_trans
        re_fb_int <- re + win_inv_trans
        ## Calculating size
        size_bn <- dat1$Size
        variantType2 <- as.character(dat1$Type)
        
        # countfre<-0 percn<-c()
        datf <- c()
        for (nn in 1:length(rf))
        {
            ## Comparing the conditions
            if (variantType2[nn] == "deletion" | variantType2[nn] == "insertion")
            {
                dat2 <- dat[which((dat$start <= rf[nn] & dat$start >= rf_wb_ind[nn] | 
                  dat$start >= rf[nn] & dat$start <= rf_fb_ind[nn]) & (dat$end >= 
                  re[nn] & dat$end <= re_wf_ind[nn] | dat$end <= re[nn] & 
                  dat$end >= re_wb_ind[nn])), ]
                size1 <- size_bn[nn]
                #print(dim(dat2))
                ## Writing if the dgv_match is TRUE or not
                
                
                if (nrow(dat2) > 1)
                {
                  countfre <- 0
                  type <- dat2$variantsubtype
                  for (ll in 1:nrow(dat2))
                  {
                    size_dgv <- dat2$end[ll] - dat2$start[ll]
                    perc <- (size1/size_dgv)
                    if (perc >= perc_similarity & (identical(type[ll], 
                      variantType2[nn])))
                      {
                      fre <- as.character(dat2$samples[ll])
                      if (length(fre) >= 1)
                      {
                        g <- grep(",", fre)
                        if (length(g) > 0)
                        {
                          fre1 <- as.character(strsplit(as.character(fre), 
                            split = ",")[[1]])
                        } else
                        {
                          fre1 <- fre
                        }
                      } else
                      {
                        fre1 <- NULL
                      }
                      countfre <- countfre + length(fre1)
                    } else
                    {
                      countfre <- countfre + 0
                    }
                    # ctr=ctr+1
                  }
                  data1 <- data.frame(dat1[nn, ], DGV_Frequency = format(round(((countfre/usamp) * 100),digits=2), scientific = FALSE), stringsAsFactors = FALSE)
                  datf <- rbind(datf, data1)
                } else if (nrow(dat2) == 1)
                {
                  # dgv_match=TRUE Calculating percentage similarity
                  type <- dat2$variantsubtype
                  size_dgv <- dat2$end - dat2$start
                  perc <- (size1/size_dgv)
                  if (perc >= perc_similarity & (identical(type, variantType2[nn])))
                  {
                    fre <- as.character(dat2$samples)
                    if (length(fre) >= 1)
                    {
                      g <- grep(",", fre)
                      countfre <- 0
                      if (length(g) > 0)
                      {
                        countfre <- length(as.character(strsplit(as.character(fre), 
                          split = ",")[[1]]))
                      } else
                      {
                        countfre <- length(fre)
                      }
                    } else
                    {
                      countfre <- NULL
                    }
                    # data1<-data.frame(dat1[nn,],DGV_Frequency=format(ct/usamp,scientific=FALSE),PercentageSimilarity=percn,stringsAsFactors
                    # = FALSE)
                  } else
                  {
                    countfre <- 0
                  }
                  # ctr=ctr+1
                  data1 <- data.frame(dat1[nn, ], DGV_Frequency = format(round(((countfre/usamp) * 100),digits=2), scientific = FALSE), stringsAsFactors = FALSE)
                  datf <- rbind(datf, data1)
                } else
                {
                  # dgv_match=FALSE
                  data1 <- data.frame(dat1[nn, ], DGV_Frequency = 0, stringsAsFactors = FALSE)
                  datf <- rbind(datf, data1)
                  # next;
                }
                
            } else if (length(grep("inversion", variantType2[nn])) >= 1 | 
                length(grep("translocation", variantType2[nn])) >= 1)
                {
                dat2 <- dat[which((dat$start <= rf[nn] & dat$start >= rf_wb_int[nn] | 
                  dat$start >= rf[nn] & dat$start <= rf_fb_int[nn]) & (dat$end >= 
                  re[nn] & dat$end <= re_fb_int[nn] | dat$end <= re[nn] & 
                  dat$end >= re_wb_int[nn])), ]
                size1 <- size_bn[nn]
                #print(dim(dat2))
                ## Writing if the dgv_match is TRUE or not
                countfre <- c()
                if (nrow(dat2) > 1)
                {
                  countfre <- 0
                  type <- dat2$variantsubtype
                  for (ll in 1:nrow(dat2))
                  {
                    size_dgv <- dat2$end[ll] - dat2$start[ll]
                    perc <- (size1/size_dgv)
                    if (perc >= perc_similarity & (identical(type[ll], 
                      variantType2[nn])))
                      {
                      fre <- as.character(dat2$samples[ll])
                      if (length(fre) >= 1)
                      {
                        g <- grep(",", fre)
                        if (length(g) > 0)
                        {
                          fre1 <- as.character(strsplit(as.character(fre), 
                            split = ",")[[1]])
                        } else
                        {
                          fre1 <- fre
                        }
                      } else
                      {
                        fre1 <- NULL
                      }
                      countfre <- countfre + length(fre1)
                    } else
                    {
                      countfre <- countfre + 0
                    }
                    # ctr=ctr+1
                  }
                  data1 <- data.frame(dat1[nn, ], DGV_Frequency = format(round(((countfre/usamp) * 100),digits=2), scientific = FALSE), stringsAsFactors = FALSE)
                  datf <- rbind(datf, data1)
                } else if (nrow(dat2) == 1)
                {
                  # dgv_match=TRUE Calculating percentage similarity
                  type <- dat2$variantsubtype
                  size_dgv <- dat2$end - dat2$start
                  perc <- (size1/size_dgv)
                  if (perc >= perc_similarity & (identical(type, variantType2[nn])))
                  {
                    fre <- as.character(dat2$samples)
                    if (length(fre) >= 1)
                    {
                      g <- grep(",", fre)
                      countfre <- 0
                      if (length(g) > 0)
                      {
                        countfre <- length(as.character(strsplit(as.character(fre), 
                          split = ",")[[1]]))
                      } else
                      {
                        countfre <- length(fre)
                      }
                    } else
                    {
                      countfre <- NULL
                    }
                    # data1<-data.frame(dat1[nn,],DGV_Frequency=format(ct/usamp,scientific=FALSE),PercentageSimilarity=percn,stringsAsFactors
                    # = FALSE)
                  } else
                  {
                    countfre <- 0
                  }
                  # ctr=ctr+1
                  data1 <- data.frame(dat1[nn, ], DGV_Frequency = format(round(((countfre/usamp) * 100),digits=2), scientific = FALSE), stringsAsFactors = FALSE)
                  datf <- rbind(datf, data1)
                } else
                {
                  # dgv_match=FALSE
                  data1 <- data.frame(dat1[nn, ], DGV_Frequency = 0, stringsAsFactors = FALSE)
                  datf <- rbind(datf, data1)
                  # next;
                }
                
            } else
            {
                data1 <- data.frame(dat1[nn, ], DGV_Frequency = 0, stringsAsFactors = FALSE)
                datf <- rbind(datf, data1)
            }
            
        }
        dataFinal <- rbind(dataFinal, datf)
    }
	##Return Mode Dataframe or Text
	if (returnMethod=="Text"){
    st1 <- strsplit(smap, ".txt")
    fname <- st1[[1]][1]
    row.names(dataFinal) <- c()
    write.table(dataFinal, paste(smappath, fname, "_DGV.txt", sep = ""), 
        sep = "\t", row.names = FALSE)
	}
	else if (returnMethod=="dataFrame"){
	return (dataFinal)
    }
	else{
	stop ("ReturnMethod Incorrect")
	}
}

#' Calculating number of samples.
#'
#' @param samp  character. Unique samples
#' @examples
#' \dontrun{
#' UniqueSample(samp)
#' }
#' @export
UniqueSample <- function(samp)
{
    samp1 <- c()
    for (ii in 1:length(samp))
    {
        g1 <- grep(",", samp[ii])
        if (length(g1 >= 1))
        {
            st <- strsplit(samp[ii], split = ",")
            gen <- as.character(st[[1]])
            samp1 <- c(samp1, as.character(unique(gen)))
        } else
        {
            if (samp[ii] != "")
            {
                samp1 <- c(samp1, as.character(samp[ii]))
            } else
            {
                next
            }
        }
    }
    return(samp1)
}
  
