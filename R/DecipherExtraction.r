#' Frequency calculation of variants compared to Decipher.
#'
#' @param decipherpath  character. Decipher Text file.
#' @param outpath character. Path where gene lists are saved.
#' @param smap_data Dataset containing smap data.
#' @param smap character Filepath for smap.
#' @param win_indel character indel window. Default 10000.
#' @param perc_similarity  Numeric . ThresholdPercentage similarity 
#' of the query SV and reference SV.
#' @param EnzymeType  boolean . Options SE and SVMerge.
#' @param input_fmt_SV  boolean . Options SE and SVMerge.
#' @param returnMethod character. Choice between 
#' text or data frame as the output.
#' @return dataframe containing decipher data.
#' are stored as text files.
#' @examples
#' decipherpath = system.file("extdata", "population_cnv.txt",
#' package="nanotatoR")
#' smappath=system.file("extdata", "GM24385_Ason_DLE1_VAP_trio5.smap", 
#' package="nanotatoR")
#' datdecipher <- Decipherfrequency (decipherpath = decipherpath, 
#' smap = smappath, win_indel = 10000,
#' EnzymeType= "SE",
#' perc_similarity = 0.5,returnMethod="dataFrame", 
#' input_fmt_SV = "Text")
#' datdecipher[1,]
#' @export
Decipherfrequency  <- function(decipherpath, smap, 
    smap_data, 
    win_indel = 10000, perc_similarity = 0.5, 
    returnMethod = c("Text", "dataFrame"),
    input_fmt_SV = c("Text", "dataFrame"),
    EnzymeType = c("SVMerge", "SE"),
    outpath)
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
    if(input_fmt_SV=="dataFrame"){
        smapdata = smap_data
        if(EnzymeType == "SVMerge"){
            #smapdata <- readSMap(smap, input_fmt_smap = "Text")
            SVID<-smapdata$SVIndex
        }
        else{
            #smapdata <- readSMap_DLE(smap, input_fmt_smap)
            SVID<-smapdata$SmapEntryID
        }
    }
    else if(input_fmt_SV=="Text"){
        if(EnzymeType == "SVMerge"){
            smapdata <- readSMap(smap, input_fmt_smap = "Text")
            SVID<-smapdata$SVIndex
        }
        else{
            smapdata <- readSMap_DLE(smap, input_fmt_smap = "Text")
            SVID<-smapdata$SmapEntryID
        }
    }
    else{
        stop("Input format for SMAP Incorrect")
    }
    r1 <- smapdata
    ## Checking Sex male/female and assigning chromosome number accordingly
    chro <- (unique(r1$RefcontigID1))
    #chro1 <- c(1:chro)
    dataFinal <- c()
    ## Extracting Data for 1 chromosome at a time and comparing Make change
    ## in XY
    
    for (ii in seq_along(chro)) # for (ii in 1:20) for (ii in 1)
    {
        # print(paste('Chromosome:', chro[ii])) Extracting data from DGV
        # dataset print(paste('ii:',ii))
        if (ii == 23)
        {
            kk <- "X"
        } else if (ii == 24)
        {
            kk <- "Y"
        } else
        {
            kk <- chro[ii]
        }
        dat <- r[which(r$chr == kk), ]
        # variantType1<-dat$variantsubtype Changing the variant terms in DGV to
        # match svmap dat$variantsubtype <- gsub('loss', 'deletion',
        # dat$variantsubtype) dat$variantsubtype <- gsub('gain', 'insertion',
        # dat$variantsubtype) dat$variantsubtype <- gsub('gain+loss',
        # 'insertion+deletion', dat$variantsubtype) variantType1 <-
        # dat$variantsubtype Extracting data from SVmap
        dat1 <- r1[which(r1$RefcontigID1 == chro[ii]), ]
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
                dat2 <- dat[which((dat$start >= rf_wb_ind[nn]  
                    & dat$end <= re_wf_ind[nn])), ]
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
                    if ((perc_ref_query >= perc_similarity 
				    & perc_query_ref >= perc_similarity)
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
                
            } else if (variantType2[nn] == "insertion" 
                | variantType2[nn] == "duplication")
                {
                dat2 <- dat[which((dat$start >= rf_wb_ind[nn]  
                    & dat$end <= re_wf_ind[nn])), ]
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
                    
                    if ((dat2$perc_ref_query[ll] >= perc_similarity 
                        & dat2$perc_query_ref[ll] >= perc_similarity)
                        & (ins_freq[ll] > 0))
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
        filename <- paste(fname, "_Decipher.txt", sep = "") 
        write.table(dataFinal, file.path(outpath, filename), 
        sep = "\t", row.names = FALSE)
    } else if (returnMethod == "dataFrame")
    {
        return(dataFinal)
    } else
    {
        stop("ReturnMethod Incorrect")
    }
}
