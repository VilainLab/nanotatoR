#' Frequency calculation of variants compared to DGV.
#'
#' @param decipherpath  character. Path to DECIPHER
#'                Text file.
#' @param smappath  character. path to the query smap file.
#' @param smap  character. File name for the smap. 
#' @param smap_data  character. Dataframe if input type chosen as dataframe.
#' @param input_fmt character. Choice between text or data frame as
#' an input to the DGV frequency calculator.
#' @param win_indel  Numeric. Insertion and deletion error window.
#' Default 10000.
#' @param perc_similarity  Numeric . ThresholdPercentage similarity 
#' of the query SV and reference SV. Default 0.5.
#' @param returnMethod character. Choice between text or data frame as the 
#' output.
#' @return Text and character vector containg gene list and terms associated 
#' with them
#'         are stored as text files.
#' @examples
#' decipherpath <- system.file("extdata", "population_cnv.txt", 
#' package = "nanotatoR")
#' smapName <- "F1.1_TestSample1_solo_hg19.smap"
#' smappath <- system.file("extdata", package = "nanotatoR")
#' win_indel=10000;win_inv_trans=50000;perc_similarity=0.5
#' decipherext<-Decipher_extraction (decipherpath, input_fmt = "Text", smappath, 
#' smap= smapName,
#' win_indel = 10000, perc_similarity = 0.5, returnMethod="dataFrame")
#' @export
Decipher_extraction <- function(decipherpath, smappath, smap, smap_data, 
    input_fmt = c("Text", "dataFrame"), win_indel = 10000, 
    perc_similarity = 0.5, 
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
    if (input_fmt == "Text")
    {
        ## Pattern matching needs to be done to remove #
        con <- file(file.path(smappath, smap), "r")
        r10 <- readLines(con, n = -1)
        close(con)
        
        g1 <- grep("RawConfidence", r10)
        g2 <- grep("RefStartPos", r10)
        if (g1 == g2)
        {
            dat <- gsub("# ", "", r10)
            # dat<-gsub('\t',' ',r1)
            dat4 <- textConnection(dat[g1:length(dat)])
            r1 <- read.table(dat4, sep = "\t", header = TRUE)
            close(dat4)
        } else
        {
            stop("column names doesnot Match")
        }
    } else if (input_fmt == "dataFrame")
    {
        r1 <- smap_data
    } else
    {
        stop("Incorrect format")
    }
    ## Checking Sex male/female and assigning chromosome number accordingly
    chro1 <- unique(r1$RefcontigID1)
    #chro1 <- c(1:chro)
    dataFinal <- c()
    ## Extracting Data for 1 chromosome at a time and comparing Make change
    ## in XY
    
    for (ii in seq_along(chro1)){
        # print(paste('Chromosome:', chro1[ii])) Extracting data from DGV
        print(paste('ii:',ii))
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
        for (nn in seq_len(length(rf))) # for (nn in 1:20)
        {
            print(paste('nn:',nn)) 
            if (variantType2[nn] == "deletion")
            {
                dat2 <- dat[which((dat$start <= rf[nn] & 
                    dat$start >= rf_wb_ind[nn] | 
                    dat$start >= rf[nn] & dat$start <= rf_fb_ind[nn]) & 
                    (dat$end >= re[nn] & dat$end <= re_wf_ind[nn] | 
                    dat$end <= re[nn] & dat$end >= re_wb_ind[nn])), ]
                size1 <- size_bn[nn]
                # print(dim(dat2)) Writing if the dgv_match is TRUE or not
                
                
                if (nrow(dat2) > 1)
                {
                    countfre <- 0
                    # type <- dat2$variantsubtype
                    del_freq <- as.numeric(dat2$deletion_frequency)
                    del_samp <- as.numeric(dat2$sample_size)
                    type = c()
                    freq = c()
                    samp_size = c()
                    # print(paste('nrow(dat2):',nrow(dat2),sep=''))
                    for (ll in seq_len(nrow(dat2))){
                        size_dec <- dat2$end[ll] - dat2$start[ll]
                        perc <- (size1/size_dec)
                    
                        if (perc >= perc_similarity & (del_freq[ll] > 0)){
                            # print(paste('del_freq:',del_freq[ll],sep=''))
                            freq <- c(freq, as.numeric(del_freq[ll]))
                            # ctr=ctr+1
                        }
                        else {
                            freq <- c(freq, 0)
                        }
                    }
                    # freq=paste(freq,collapse=',')
                    if (mean(freq) == 0){
                        freq = 0
                    } 
                    else {
                        freq = round(mean(freq), 2)
                    }
                 
                    data1 <- data.frame(dat1[nn, ], DECIPHER_Frequency = freq, 
                        stringsAsFactors = FALSE)
                    datf <- rbind(datf, data1)
                } 
                else if (nrow(dat2) == 1){
                    # dgv_match=TRUE Calculating percentage similarity
                    del_freq <- as.numeric(dat2$deletion_frequency)
                    del_samp <- as.numeric(dat2$sample_size)
                    size_dec <- dat2$end - dat2$start
                    perc <- (size1/size_dec)
                  # ctr=ctr+1
                  if (perc >= perc_similarity & (del_freq > 0)){
                    # print(paste('del_freq:',del_freq,sep=''))
                    freq <- round(del_freq, 2)
                    # ctr=ctr+1
                  } else {
                    freq <- 0
                    
                  }
                  # freq=paste(as.numeric(freq),collapse=',')
                  # samp_size=paste(as.numeric(samp_size),collapse=',')
                  
                    data1 <- data.frame(dat1[nn, ], DECIPHER_Frequency = freq, 
                        stringsAsFactors = FALSE)
                    datf <- rbind(datf, data1)
                } 
                else{
                    # dgv_match=FALSE
                    data1 <- data.frame(dat1[nn, ], DECIPHER_Frequency = 0, 
                    stringsAsFactors = FALSE)
                    datf <- rbind(datf, data1)
                    # next;
                }
                
            } 
            else if (variantType2[nn] == "insertion" | variantType2[nn] == 
                "duplication"){
                dat2 <- dat[which((dat$start <= rf[nn] & 
                    dat$start >= rf_wb_ind[nn] | 
                    dat$start >= rf[nn] & dat$start <= rf_fb_ind[nn]) & 
                    (dat$end >= re[nn] & dat$end <= re_wf_ind[nn] | 
                    dat$end <= re[nn] & dat$end >= re_wb_ind[nn])), ]
                size1 <- size_bn[nn]
                # print(dim(dat2)) Writing if the dgv_match is TRUE or not
                
                
                if (nrow(dat2) > 1) {
                    # print(paste('nrow(dat2):',nrow(dat2),sep=''))
                    countfre <- 0
                    # type <- dat2$variantsubtype
                    ins_freq <- as.numeric(dat2$duplication_frequency)
                    ins_samp <- as.numeric(dat2$sample_size)
                    type = c()
                    freq = c()
                    samp_size = c()
                    for (ll in seq_len(nrow(dat2))) {
                        size_dec <- dat2$end[ll] - dat2$start[ll]
                        perc <- (size1/size_dec)
                    
                        if (perc >= perc_similarity & (ins_freq[ll] > 0)){
                            # print(paste('ins_freq:',ins_freq[ll],sep=''))
                            freq <- c(freq, as.numeric(ins_freq[ll]))
                            # type<-c(type,'Insertion')
                      
                        } 
                        else {
                            freq <- c(freq, 0)
                            # type<-c(type,'-') samp_size<-c(samp_size,'-')
                        }
                    }
                    if (mean(freq) == 0) {
                        freq = 0
                    } 
                    else {
                        freq = round(mean(freq), 2)
                    }
                  
                 
                  
                data1 <- data.frame(dat1[nn, ], DECIPHER_Frequency = freq, 
                    stringsAsFactors = FALSE)
                datf <- rbind(datf, data1)
                } 
                else if (nrow(dat2) == 1) {
                    # dgv_match=TRUE Calculating percentage similarity
                    ins_freq <- as.numeric(dat2$duplication_frequency)
                    ins_samp <- as.numeric(dat2$sample_size)
                    size_dec <- dat2$end - dat2$start
                    perc <- (size1/size_dec)
                    # ctr=ctr+1
                    if (perc >= perc_similarity & (ins_freq > 0)) {
                        # print(paste('ins_freq:',ins_freq,sep=''))
                        freq <- round(ins_freq, 2)
                    
                    # ctr=ctr+1
                    } 
                    else {
                        freq <- 0
                    
                    }
                  
                    data1 <- data.frame(dat1[nn, ], DECIPHER_Frequency = freq, 
                        stringsAsFactors = FALSE)
                    datf <- rbind(datf, data1)
                } 
                else {
                    # dgv_match=FALSE
                    data1 <- data.frame(dat1[nn, ], DECIPHER_Frequency = 0, 
                        stringsAsFactors = FALSE)
                    datf <- rbind(datf, data1)
                    # next;
                }
            } 
            else {
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
