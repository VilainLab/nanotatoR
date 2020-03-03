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
#' decipherpath="Z:/Suro/Annotator/Data/population_cnv.txt"
#' smappath="Z:/Suro/Annotator/Data/";
#' smap="F1.1_UDN287643_P_Q.S_VAP_SVmerge_trio_original2.txt";
#' win_indel=10000;win_inv_trans=50000;perc_similarity=0.5
#' DGV_extraction (hgpath, smappath, win_indel = 10000, win_inv_trans = 50000,
#' perc_similarity = 0.5,returnMethod="dataFrame")
#' }
#' @export
Decipher_extraction <- function(decipherpath, smappath, smap, smap_data, 
    input_fmt = c("Text", "dataFrame"), win_indel = 10000, perc_similarity = 0.5, 
    returnMethod = c("Text", "dataFrame"))
    {
    # S='F' Change the window for Inversion/translocation 50000
    print("###Calculating Decipher Frequency###")
    ## Reading the DGV file and the smap file
    r <- read.table(decipherpath, header = TRUE, sep = "\t")
    ## Unique variant length extraction for frquency Calculation samp <-
    ## as.character(unique(r$samples)) usamp <- UniqueSample(samp)
    #usamp <- 7965
    datfinal <- data.frame()
    # varaccl<-length(unique(varacc)) close(con) Checking if the input
    # format is dataframe or Text
     if(input_fmt == "Text"){
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
		r1 <- cbind(SampleID = rep(str_squish(as.character(SampleID)), times = nrow(r1)), r1)
		}
		else if(input_fmt == "dataFrame"){
			r1<-smap_data
		}
		else{
			stop("Input Format incorrect")
		}
    ## Checking Sex male/female and assigning chromosome number accordingly
    chro <- length(unique(r1$RefcontigID1))
    chro1 <- c(1:chro)
    dataFinal <- c()
    ## Extracting Data for 1 chromosome at a time and comparing Make change
    ## in XY
    
    for (ii in 1:length(chro1)) # for (ii in 1:20) for (ii in 1)
    {
        # print(paste('Chromosome:', chro1[ii])) Extracting data from DGV
        # dataset print(paste('ii:',ii))
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
        # match svmap dat$variantsubtype <- gsub('loss', 'deletion',
        # dat$variantsubtype) dat$variantsubtype <- gsub('gain', 'insertion',
        # dat$variantsubtype) dat$variantsubtype <- gsub('gain+loss',
        # 'insertion+deletion', dat$variantsubtype) variantType1 <-
        # dat$variantsubtype Extracting data from SVmap
        dat1 <- r1[which(r1$RefcontigID1 == chro1[ii]), ]
        ## Adding the windows to the breakpoints
        rf <- dat1$RefStartPos
        rf_wb_ind <- rf - win_indel
        rf_fb_ind <- rf + win_indel
        re <- dat1$RefEndPos
        re_wf_ind <- re + win_indel
        re_wb_ind <- re - win_indel
        ## Calculating size
        size_bn <- dat1$Size
        variantType2 <- as.character(dat1$Type)
        
        # countfre<-0 percn<-c()
        datf <- c()
        for (nn in 1:length(rf)) # for (nn in 1:20)
        {
            # print(paste('nn:',nn)) Comparing the conditions
            if (variantType2[nn] == "deletion")
            {
                dat2 <- dat[which((dat$start >= rf_wb_ind[nn]  & dat$end <= re_wf_ind[nn])), ]
                size1 <- size_bn[nn]
                # print(dim(dat2)) Writing if the dgv_match is TRUE or not
                
                
                if (nrow(dat2) > 1)
                {
                  countfre <- 0
                  # type <- dat2$variantsubtype
                  del_freq <- as.numeric(dat2$deletion_frequency)
                  del_samp <- as.numeric(dat2$sample_size)
				  dat2$size_dec <- dat2$end - dat2$start
				  dat2$perc_ref_query <- as.numeric(dat2$size_dec)/size1
				  dat2$perc_query_ref <- size1/as.numeric(dat2$size_dec)
                  type = c()
                  freq = c()
                  samp_size = c()
                  # print(paste('nrow(dat2):',nrow(dat2),sep=''))
                  for (ll in 1:nrow(dat2))
                  {
				     
                    size_dec <- dat2$end[ll] - dat2$start[ll]
                    #perc <- (size1/size_dec)
                    
                    if ((dat2$perc_ref_query[ll] >= perc_similarity & dat2$perc_query_ref[ll] >= perc_similarity) & (del_freq[ll] > 0))
                    {
                      # print(paste('del_freq:',del_freq[ll],sep=''))
                      freq <- c(freq, as.numeric(del_freq[ll]))
                      # ctr=ctr+1
                    } else
                    {
                      freq <- c(freq, 0)
                    }
                  }
                  # freq=paste(freq,collapse=',')
                  if (mean(freq) == 0)
                  {
                    freq = 0
                  } else
                  {
                    freq = round(mean(freq), 2)
                  }
                  # samp_size=paste(samp_size,collapse=',') type=paste(type,collapse=',')
                  data1 <- data.frame(dat1[nn, ], DECIPHER_Frequency = freq, 
                    stringsAsFactors = FALSE)
                  datf <- rbind(datf, data1)
                } else if (nrow(dat2) == 1)
                {
                  # dgv_match=TRUE Calculating percentage similarity
                  del_freq <- as.numeric(dat2$deletion_frequency)
                  del_samp <- as.numeric(dat2$sample_size)
                  size_dec <- dat2$end - dat2$start
				  perc_ref_query <- as.numeric(size_dec)/size1
				  perc_query_ref <- size1/as.numeric(size_dec)
                  # ctr=ctr+1
                  if ((perc_ref_query >= perc_similarity & perc_query_ref >= perc_similarity)
 				  & (del_freq > 0))
                  {
                    # print(paste('del_freq:',del_freq,sep=''))
                    freq <- round(del_freq, 2)
                    # ctr=ctr+1
                  } else
                  {
                    freq <- 0
                    
                  }
                  # freq=paste(as.numeric(freq),collapse=',')
                  # samp_size=paste(as.numeric(samp_size),collapse=',')
                  
                  data1 <- data.frame(dat1[nn, ], DECIPHER_Frequency = freq, 
                    stringsAsFactors = FALSE)
                  datf <- rbind(datf, data1)
                } else
                {
                  # dgv_match=FALSE
                  data1 <- data.frame(dat1[nn, ], DECIPHER_Frequency = 0, 
                    stringsAsFactors = FALSE)
                  datf <- rbind(datf, data1)
                  # next;
                }
                
            } else if (variantType2[nn] == "insertion" | variantType2[nn] == 
                "duplication")
                {
                dat2 <- dat[which((dat$start >= rf_wb_ind[nn]  & dat$end <= re_wf_ind[nn])), ]
                size1 <- size_bn[nn]
                # print(dim(dat2)) Writing if the dgv_match is TRUE or not
                
                
                if (nrow(dat2) > 1)
                {
                  # print(paste('nrow(dat2):',nrow(dat2),sep=''))
                  countfre <- 0
                  # type <- dat2$variantsubtype
                  ins_freq <- as.numeric(dat2$duplication_frequency)
                  ins_samp <- as.numeric(dat2$sample_size)
                  dat2$size_dec <- dat2$end - dat2$start
				  dat2$perc_ref_query <- as.numeric(dat2$size_dec)/size1
				  dat2$perc_query_ref <- size1/as.numeric(dat2$size_dec)
				  type = c()
                  freq = c()
                  samp_size = c()
                  for (ll in 1:nrow(dat2))
                  {
                    'size_dec <- dat2$end[ll] - dat2$start[ll]
                    perc <- (size1/size_dec)'
                    
                    if ((dat2$perc_ref_query[ll] >= perc_similarity & dat2$perc_query_ref[ll] >= perc_similarity) & (ins_freq[ll] > 0))
                    {
                      # print(paste('ins_freq:',ins_freq[ll],sep=''))
                      freq <- c(freq, as.numeric(ins_freq[ll]))
                      # type<-c(type,'Insertion')
                      # samp_size<-c(samp_size,as.numeric(ins_samp[ll])) ctr=ctr+1
                    } else
                    {
                      freq <- c(freq, 0)
                      # type<-c(type,'-') samp_size<-c(samp_size,'-')
                    }
                  }
                  if (mean(freq) == 0)
                  {
                    freq = 0
                  } else
                  {
                    freq = round(mean(freq), 2)
                  }
                  
                  # samp_size=paste(samp_size,collapse=',') type=paste(type,collapse=',')
                  
                  data1 <- data.frame(dat1[nn, ], DECIPHER_Frequency = freq, 
                    stringsAsFactors = FALSE)
                  datf <- rbind(datf, data1)
                } else if (nrow(dat2) == 1)
                {
                  # dgv_match=TRUE Calculating percentage similarity
                  ins_freq <- as.numeric(dat2$duplication_frequency)
                  ins_samp <- as.numeric(dat2$sample_size)
                  size_dec <- dat2$end - dat2$start
				  perc_ref_query <- as.numeric(size_dec)/size1
				  perc_query_ref <- size1/as.numeric(size_dec)
                  # ctr=ctr+1
                  if ((perc_ref_query >= perc_similarity & perc_query_ref >= perc_similarity)
     			  & (ins_freq > 0))	
                  {
                    # print(paste('ins_freq:',ins_freq,sep=''))
                    freq <- round(ins_freq, 2)
                    
                    # ctr=ctr+1
                  } else
                  {
                    freq <- 0
                    
                  }
                  
                  data1 <- data.frame(dat1[nn, ], DECIPHER_Frequency = freq, 
                    stringsAsFactors = FALSE)
                  datf <- rbind(datf, data1)
                } else
                {
                  # dgv_match=FALSE
                  data1 <- data.frame(dat1[nn, ], DECIPHER_Frequency = 0, 
                    stringsAsFactors = FALSE)
                  datf <- rbind(datf, data1)
                  # next;
                }
            } else
            {
                # dgv_match=FALSE
                data1 <- data.frame(dat1[nn, ], DECIPHER_Frequency = 0, 
                  stringsAsFactors = FALSE)
                datf <- rbind(datf, data1)
                # next;
            }
            
        }
        dataFinal <- rbind(dataFinal, datf)
    }
    
    
    ## Return Mode Dataframe or Text
    if (returnMethod == "Text")
    {
        st1 <- strsplit(smap, ".txt")
        fname <- st1[[1]][1]
        row.names(dataFinal) <- c()
        write.table(dataFinal, paste(smappath, fname, "_Decipher.txt", 
            sep = ""), sep = "\t", row.names = FALSE)
    } else if (returnMethod == "dataFrame")
    {
        return(dataFinal)
    } else
    {
        stop("ReturnMethod Incorrect")
    }
}
