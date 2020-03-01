#' Frequency calculation of variants compared to DGV.
#'
#' @param hgpath  character. Path to Database of Genomic Variants (DGV)
#'                Text file.
#' @param smappath  character. Path for smap textfile.
#' @param smap character. File name for smap textfile.
#' @param input_fmt_DGV character. Choice between text or data frame as
#' an input to the DGV frequency calculator.
#' @param smap_data dataframe. Dataset containing smap data.
#' @param win_indel_DGV  Numeric. Insertion and deletion error window.Default 10000
#' bases.
#' @param win_inv_trans_DGV  Numeric. Inversion and translocation error window.
#' Default 50000 bases.
#' @param perc_similarity_DGV  Numeric . ThresholdPercentage similarity 
#' of the query SV and reference SV. Default 0.5.
#' @param returnMethod character. Choice between text or data frame as the output.
#' @param outpath character. Path where gene lists are saved.
#' @return Text and character vector containg gene list and terms associated with them
#'         are stored as text files.
#' @examples
#' \dontrun{
#' hgpath="C:\\nanotatoR\\Data\\GRCh37_hg19_variants_2016-05-15.txt";
#' smappath="Z:/Hayk's_Materials/Bionano/Projects/UDN/F1_UDN287643_Benic.Aria/";
#' smap="F1_UDN287643_P_Q.S_VAP_SVmerge_trio_access.txt";
#' win_indel_DGV=10000;win_inv_trans_DGV=50000;perc_similarity_DGV=0.5
#' DGV_extraction (hgpath, smappath, win_indel_DGV = 10000, win_inv_trans_DGV = 50000,
#' perc_similarity_DGV = 0.5,returnMethod="dataFrame")
#' }
#' @import utils
#' @export
DGV_extraction <-
  function(hgpath,
           smappath,
           smap,
           smap_data,
           input_fmt_DGV = c("Text", "DataFrame"),
           win_indel_DGV = 10000,
           win_inv_trans_DGV = 50000,
           perc_similarity_DGV = 0.5,
           returnMethod = c("Text", "dataFrame"),outpath){
    # S='F' Change the window for Inversion/translocation 50000
    
    ## Reading the DGV file and the smap file
    r <- read.table(paste(hgpath), header = TRUE, sep = "\t")
    ## Unique variant length extraction for frquency Calculation
    #samp <- as.character(unique(r$samples))
    #countsamp <- UniqueSample(samp)
    countsamp <- 54946
    # varaccl<-length(unique(varacc)) close(con)
    ##Checking if the input format is dataframe or Text
    if(input_fmt_DGV == "Text"){
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
        r1 <- cbind(
            SampleID = rep(str_squish(as.character(SampleID)), 
                times = nrow(r1)), r1)
        }
        else if(input_fmt_DGV == "DataFrame"){
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
    for (ii in 1:length(chro1))
    {
        print(paste("Chromosome:", chro1[ii]))
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
        dat <- r[which(r$chr == kk),]
        # variantType1<-dat$variantsubtype Changing the variant terms in DGV to
        # match svmap
        dat$variantsubtype <-
        gsub("loss", "deletion", dat$variantsubtype)
        dat$variantsubtype <- gsub("gain", "insertion", dat$variantsubtype)
        dat$variantsubtype <- gsub("mobile element insertion", "insertion",as.character(dat$variantsubtype))
        dat$variantsubtype <- gsub("novel sequence insertion", "insertion", as.character(dat$variantsubtype))
       
        variantType1 <- dat$variantsubtype
        ## Extracting data from SVmap
        dat1 <- r1[which(r1$RefcontigID1 == chro1[ii]),]
        ## Adding the windows to the breakpoints
        rf <- dat1$RefStartPos
        rf_wb_ind <- rf - win_indel_DGV
        rf_fb_ind <- rf + win_indel_DGV
        rf_wb_int <- rf - win_inv_trans_DGV
        rf_fb_int <- rf + win_inv_trans_DGV
        re <- dat1$RefEndPos
        re_wf_ind <- re + win_indel_DGV
        re_wb_ind <- re - win_indel_DGV
        re_wb_int <- re - win_inv_trans_DGV
        re_fb_int <- re + win_inv_trans_DGV
        ## Calculating size
        size_bn <- dat1$Size
        variantType2 <- as.character(dat1$Type)
      
        # countfre<-0 percn<-c()
        datf <- c()
        for (nn in 1:length(rf))
        {#print(nn)
            ## Comparing the conditions
            if (variantType2[nn] == "deletion" |
                variantType2[nn] == "insertion")
            {
                dat2 <-dat[which((dat$start >= rf_wb_ind[nn] & dat$end <= re_wf_ind[nn])),]
                size1 <- size_bn[nn]
                #print(dim(dat2))
                ## Writing if the dgv_match is TRUE or not
                  
                if (nrow(dat2) > 1)
                {
                    countfre <- 0
                    ##countsamp<- 0
                    countfreunfilt <- 0
                    type <- dat2$variantsubtype
                    if(variantType2[nn] == "insertion"){
                        freq <- as.numeric(dat2$observedgains)
                    }
                else{
                    freq <- as.numeric(dat2$observedlosses)
                }
                freq[is.na(freq)]<-0
                dat2$size_dgv <- dat2$end - dat2$start
                #perc <- (size1 / size_dgv)
                dat2$perc_ref_query <- as.numeric(dat2$size_dgv)/size1
                dat2$perc_query_ref <- size1/as.numeric(dat2$size_dgv)            
                samp <- as.numeric(dat2$samplesize)
                for (ll in 1:nrow(dat2)){
              
                    if ((dat2$perc_ref_query[ll] >= perc_similarity_DGV 
					    & dat2$perc_query_ref[ll] >= perc_similarity_DGV)
						& (identical(type[ll], variantType2[nn]))){
                        countfre <- countfre + as.numeric(freq[ll])
                        ##countsamp <- countsamp + as.numeric(samp[ll])
                    } else
                    {
                        countfre <- countfre + 0
                    }
                # ctr=ctr+1
                }
                if (countfre==0){
                    DGV_Freq_Perc=0; DGV_Count=0
                }else{
                    DGV_Freq_Perc = format(
                        ((countfre / countsamp) * 100), scientific = FALSE
                        )
                    DGV_Count=countfre
                }    
                data1 <- data.frame(dat1[nn,], 
				    DGV_Count=as.numeric(DGV_Count), 
					DGV_Freq_Perc = as.numeric(DGV_Freq_Perc) , 
					stringsAsFactors = FALSE)
                datf <- rbind(datf, data1)
                } else if (nrow(dat2) == 1){
                    # dgv_match=TRUE Calculating percentage similarity
                    type <- dat2$variantsubtype
                    size_dgv <- dat2$end - dat2$start
                    perc_ref_query <- as.numeric(size_dgv)/size1
                    perc_query_ref <- size1/as.numeric(size_dgv)
                    if(variantType2[nn] == "insertion"){
                        freq <- as.numeric(dat2$observedgains)
                    }
                    else{
                        freq <- as.numeric(dat2$observedlosses)
                    }
                freq[is.na(freq)]<-0        
                samp <- as.numeric(dat2$samplesize)
                if ((perc_ref_query >= perc_similarity_DGV 
				    & perc_query_ref >= perc_similarity_DGV) 
					& (identical(type, variantType2[nn]))){
              
                    countfre <- countfre + as.numeric(freq)
                } else{
                    countfre <- 0
                }
                if (countfre==0){
                    DGV_Freq_Perc=0;DGV_Count=0
                }else{
                    DGV_Freq_Perc = format(((countfre / countsamp) * 100), 
					        scientific = FALSE)
                    DGV_Count = countfre
                }    
                data1 <- data.frame(dat1[nn,], 
				    DGV_Count = DGV_Count, 
					DGV_Freq_Perc =as.numeric(DGV_Freq_Perc) , 
					stringsAsFactors = FALSE)
                datf <- rbind(datf, data1)
            } else
            {
            # dgv_match=FALSE
                data1 <- data.frame(dat1[nn,],DGV_Count = 0, 
				    DGV_Freq_Perc = 0, stringsAsFactors = FALSE)
                datf <- rbind(datf, data1)
            # next;
            }
          
        }
        else if ((
            variantType2[nn] == "duplication" 
            | variantType2[nn] == "duplication_split" 
            | variantType2[nn] == "duplication_inverted"))
        {
            dat2 <-dat[which((dat$start >= rf_wb_ind[nn] 
                & dat$end <= re_wf_ind[nn])),]
            size1 <- size_bn[nn]
            ## Writing if the dgv_match is TRUE or not
          
          
            if (nrow(dat2) > 1){
                countfre <- 0
                dat2$size_dgv <- dat2$end - dat2$start
                #perc <- (size1 / size_dgv)
                dat2$perc_ref_query <- as.numeric(dat2$size_dgv)/size1
                dat2$perc_query_ref <- size1/as.numeric(dat2$size_dgv)
            
                type <- dat2$variantsubtype
                freq <- as.numeric(dat2$observedlosses)
                samp <- as.numeric(dat2$samplesize)
                freq[is.na(freq)]<-0
                for (ll in 1:nrow(dat2)){
              
                #perc <- (size1 / size_dgv)
                if ( (identical(type[ll], "duplication") 
				    | identical(type[ll], "duplication_split")
					| identical(type[ll], "duplication_inverted"))){
                    countfre <- countfre + as.numeric(freq[ll])
                    ##countsamp <- countsamp + as.numeric(samp[ll])
                    'countfre <- countfre + length(fre1)'
                }
                else if ( (dat2$perc_ref_query[ll] >= perc_similarity_DGV 
				    & dat2$perc_query_ref[ll] >= perc_similarity_DGV)
					& (identical(type[ll], "insertion"))){
                    countfre <- countfre + as.numeric(freq[ll])
                }
                else
                {
                    countfre <- countfre + 0
                }
              # ctr=ctr+1
                }
                if (countfre==0){
                    DGV_Freq_Perc=0;DGV_Count = 0
                }else{
                    DGV_Freq_Perc = format(
                        ((countfre / countsamp) * 100), scientific = FALSE)
                    DGV_Count = countfre
                }    
                data1 <- data.frame( dat1[nn,], 
				    DGV_Count = as.numeric(countfre), 
					DGV_Freq_Perc = as.numeric(DGV_Freq_Perc), 
					stringsAsFactors = FALSE)
                datf <- rbind(datf, data1)
            } else if (nrow(dat2) == 1){
                # dgv_match=TRUE Calculating percentage similarity
                type <- dat2$variantsubtype
                size_dgv <- dat2$end - dat2$start
                #perc <- (size1 / size_dgv)
                perc_ref_query <- as.numeric(size_dgv)/size1
                perc_query_ref <- size1/as.numeric(size_dgv)
                freq <- as.numeric(dat2$observedlosses)
                freq[is.na(freq)]<-0
                samp <- as.numeric(dat2$samplesize)
                if ((identical(type, variantType2[nn]))){
                    countfre <- countfre + as.numeric(freq)
                } else if ( (perc_ref_query >= perc_similarity_DGV 
				    & perc_query_ref >= perc_similarity_DGV)
					& (identical(type[ll], "insertion"))){
                    countfre <- countfre + as.numeric(freq)
                } 
            else
            {
              countfre <- 0
            }
            
            if (countfre==0){
               DGV_Freq_Perc=0;DGV_Count = 0
            }else{
               DGV_Freq_Perc = format(((countfre / countsamp) * 100), 
                    scientific = FALSE)
               DGV_Count = as.numeric(countfre)
            }
            # ctr=ctr+1
                data1 <- data.frame(
				    dat1[nn,],DGV_Count = as.numeric(countfre), 
					DGV_Freq_Perc = as.numeric(DGV_Freq_Perc) , 
					stringsAsFactors = FALSE)
                datf <- rbind(datf, data1)
            } else
            {
                # dgv_match=FALSE
                data1 <- data.frame(dat1[nn,], 
                    DGV_Count = 0, 
                    DGV_Freq_Perc = 0, 
                    stringsAsFactors = FALSE)
                datf <- rbind(datf, data1)
                # next;
            }
          
        }
        else if (length(grep("inversion", variantType2[nn])) >= 1){
           dat2 <- dat[which(( dat$start <= rf[nn] & 
                        dat$start >= rf_wb_int[nn] | dat$start >= rf[nn] &
                        dat$start <= rf_fb_int[nn]) & 
                        (dat$end >= re[nn] & dat$end <= re_fb_int[nn] | 
                        dat$end <= re[nn] & dat$end >= re_wb_int[nn])),]
            size1 <- size_bn[nn]
            #print(dim(dat2))
            ## Writing if the dgv_match is TRUE or not
            countfre <- 0
            #countsamp<- 0
            countfreunfilt <- 0
            if (nrow(dat2) > 1){
                countfre <- 0
                countfreunfilt <- 0
                type <- dat2$variantsubtype
                samp <- as.numeric(dat2$samplesize)
                for (ll in 1:nrow(dat2))
                {#print(ll)
              
                    if ((identical(type[ll],variantType2[nn]))){#print("1")
                
                        countfre <- countfre + 1
                        ##countsamp <- countsamp + as.numeric(samp[ll])
                    } else{#print("2")
                        countfre <- countfre + 0
                    }
                    # ctr=ctr+1
                }
            if (countfre==0){
               DGV_Freq_Perc=0; DGV_Count = 0
            }else{
                DGV_Freq_Perc = format(((countfre / countsamp) * 100), 
                    scientific = FALSE)
                DGV_Count = as.numeric(countfre)
            }
            data1 <-  data.frame(dat1[nn,], 
                DGV_Count = as.numeric(countfre), 
                DGV_Freq_Perc = as.numeric(DGV_Freq_Perc) , 
                stringsAsFactors = FALSE )
            datf <- rbind(datf, data1)
            } else if (nrow(dat2) == 1){
                # dgv_match=TRUE Calculating percentage similarity
                type <- dat2$variantsubtype
                samp <- as.numeric(dat2$samplesize)
                countfre <- 0
                #countsamp<- 0
                if ((identical(type, variantType2[nn]))){   
                    countfre <- countfre + 1
                    ##countsamp <- countsamp + as.numeric(samp[ll])
              
                } else{
                    countfre <- 0
                }
                # data1<-data.frame(dat1[nn,],DGV_Freq_Perc=format(ct/countsamp,scientific=FALSE),PercentageSimilarity=percn,stringsAsFactors
                # = FALSE)
            if (countfre==0){
                DGV_Freq_Perc=0; DGV_Count = 0
            } else{
                DGV_Freq_Perc = format(((countfre / countsamp) * 100), 
                    scientific = FALSE)
                DGV_Count = as.numeric(countfre)
            }
            # ctr=ctr+1
            data1 <- data.frame( dat1[nn,], 
                DGV_Count = as.numeric(countfre), 
                DGV_Freq_Perc = as.numeric(DGV_Freq_Perc), 
                stringsAsFactors = FALSE)
            datf <- rbind(datf, data1)
            } else{
            # dgv_match=FALSE
                data1 <- data.frame(dat1[nn,], DGV_Count = 0, 
                      DGV_Freq_Perc = 0, stringsAsFactors = FALSE)
                datf <- rbind(datf, data1)
            # next;
            }
          
        } else{
            # dgv_match=FALSE
            data1 <- data.frame(dat1[nn,], DGV_Count = 0, DGV_Freq_Perc = 0, stringsAsFactors = FALSE)
            datf <- rbind(datf, data1)
            # next;
        }
        
    }
    dataFinal <- rbind(dataFinal, datf)
}
##Return Mode Dataframe or Text
if (returnMethod == "Text") {
    st1 <- strsplit(smap, ".txt")
    fname <- st1[[1]][1]
    row.names(dataFinal) <- c()
    write.table(
        dataFinal,
        paste(outpath, fname, "_DGV.txt", sep = ""),
        sep = "\t",
        row.names = FALSE
      )
    }
else if (returnMethod == "dataFrame") {
    return (dataFinal)
    }
else{
    stop ("ReturnMethod Incorrect")
    }
}


