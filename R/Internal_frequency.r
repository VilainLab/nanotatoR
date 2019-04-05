#' Merges Solo SV files to one common SV file
#'
#' @param path  character. Path to the solo files.
#' @param pattern  character. file name pattern for solo files.
#' @param dbOutput  character. Output type. Default text and dataframe.
#' @param fname  character. filename if dbOutput = text.
#' @param outpath  character. path for output file if dbOutput = text.
#' @return Text file containing all the solo SMAP files.
#' @examples
#' path <- system.file("extdata", "SoloFile/", package = "nanotatoR")
#' pattern <- "_hg19.smap"
#' mergedFiles <- makeMergedSVData(path, pattern, dbOutput = "dataframe")
#' mergedFiles[1,]
#' @importFrom stats na.omit  
#' @import utils
#' @export
makeMergedSVData <- function(path, pattern, outpath, fname,
                        dbOutput=c("dataframe","text")){
    setwd(path)
    l <- list.files(".", pattern)
    nam <- c()
    datfinal <- data.frame()
    for (ii in seq_along((l)))
    {
        print(l[ii])
        con <- file(l[ii], "r")
        r10 <- readLines(con, n = -1)
        close(con)
        g1 <- grep("RefEndPos", r10)
        g2 <- grep("RefStartPos", r10)
        # r10<-as.character(r10)
        if (g1 == g2)
        {
            dat <- gsub("# ", "", as.character(r10))
            # dat<-gsub('\t',' ',as.character(dat))
            dat4 <- textConnection(dat[g1:length(dat)])
            ##### print(dim(dat4))
            r1 <- read.table(dat4, sep = "\t", header = TRUE)
            nam <- names(r1)
            r1 <- data.frame(r1)
            st <- strsplit(l[ii], split = ".txt")
            fname <- st[[1]][1]
            fname1 <- c(rep(fname, length(r1[, 1])))
            r1 <- cbind(fname1, r1)
            #### print(dim(r1))
            datfinal <- data.frame(rbind(datfinal, r1))
            close(dat4)
            #### print(dim(datfinal))
        } else
        {
            stop("column names doesnot Match")
        }
    }
    names(datfinal)[1] <- "SVIdentifier"
    if(dbOutput == "text"){
        filename <- paste(fname, "_merged.txt", sep = "")
        write.table(datfinal, file.path(outpath, filename), sep = "\t")
    }
    else if(dbOutput == "dataframe"){
        return(datfinal)
    }
    else{stop("method not found")}
}

#' Calculates the internal frequencies of SV in internal cohorts
#'
#' @param mergedFiles  character. Path to the merged SV files.
#' @param smappath  character. path to the query smap file.
#' @param smap  character. File name for the smap
#' @param smapdata  character. Dataframe if input type chosen as dataframe.
#' @param buildSVInternalDB  boolean. Checking whether the merged solo 
#' file database exist.
#' @param input_fmt Format in which data is provided as an input to the 
#' function.
#' @param path  character. Path to the solo file database.
#' @param pattern  character. pattern of the file names to merge.
#' @param outpath  character. Path to merged SV solo datasets.
#' @param dbOutput  character. Output of merged bionano data.
#' @param fname  character. Filename in case dbOutput = Text.
#' @param win_indel  Numeric. Insertion and deletion error window.
#' Default 10000.
#' @param win_inv_trans  Numeric. Inversion and translocation error window.
#' Default 50000.
#' @param perc_similarity  Numeric . ThresholdPercentage similarity 
#' of the query SV and reference SV. Default 0.5.
#' @param indelconf  Numeric. Threshold for insertion and deletion confidence.
#' Default 0.5
#' @param invconf  Numeric. Threshold for inversion confidence.Default 0.01.
#' @param transconf  Numeric. Threshold for translocation confidence. 
#' Default 0.1.
#' @param limsize Numeric. Minimum size of SV that can be determined 
#' acurately by the Bionano SV caller. Default 1000.
#' @param win_indel_parents  Numeric. Insertion and deletion error window to 
#' determine zygosity in case of parents. Default 5000.
#' @param win_inv_trans_parents  Numeric. Inversion and translocation error 
#' window to determine zygosity in case of parents. Default 40000.
#' @param returnMethod character. Choice between Text and DataFrame.
#' @return Text file or data frames containing internalFrequency data.
#' @examples
#' path <- system.file("extdata", "SoloFile", package = "nanotatoR")
#' pattern <- "_hg19.smap"
#' smapName <- "F1.1_TestSample1_solo_hg19.smap"
#' smappath <- system.file("extdata",  package = "nanotatoR")
#' indelconf = 0.5; invconf = 0.01;transconf = 0.1;input_fmt="Text";
#' internalFrequency(path = path, pattern = pattern,
#' dbOutput =c("dataframe"), smappath = smappath ,
#' smap = smapName, buildSVInternalDB=TRUE, win_indel=10000, 
#' win_inv_trans=50000, 
#' perc_similarity=0.5, indelconf=0.5, invconf=0.01, fname= "Solo",
#' transconf=0.1, limsize=1000, win_indel_parents=5000,input_fmt="Text",
#' win_inv_trans_parents=40000,
#' returnMethod="dataFrame")
#' @importFrom stats na.omit 
#' @import hash
#' @export
internalFrequency <- function(mergedFiles, smappath , smap , 
    buildSVInternalDB=FALSE, smapdata, input_fmt=c("Text","dataFrame"), 
    path, pattern, outpath, win_indel=10000, win_inv_trans=50000, 
    perc_similarity=0.5, indelconf=0.5, invconf=0.01, fname,
    limsize=1000, win_indel_parents=5000,win_inv_trans_parents=40000, 
    transconf=0.1, dbOutput =c("dataframe","text"),
    returnMethod=c("Text","dataFrame"))
    {
    #library(hash)
    print("###Calculating the Internal Frequency###")
    if(buildSVInternalDB == TRUE){
        r <- makeMergedSVData(path, pattern, outpath, dbOutput = dbOutput)
    } 
    else{
        r <- read.table(mergedFiles, sep = "\t", header = TRUE)
    }
    usamp <- length(unique(r$SVIdentifier))
    if(input_fmt == "Text"){
        con <- file(file.path(smappath, smap), "r")
        r10 <- readLines(con, n = -1)
        close(con)
        # datfinal<-data.frame()
        g1 <- grep("RawConfidence", r10)
        g2 <- grep("RefStartPos", r10)
        if (g1 == g2)
        {
            dat <- gsub("# ", "", as.character(r10))
            # dat<-gsub('\t',' ',r10)
            dat4 <- textConnection(dat[g1:length(dat)])
            r1 <- read.table(dat4, sep = "\t", header = TRUE)
            close(dat4)
        } 
        else{
            stop("column names doesnot Match")
        }
    }
    else if(input_fmt=="dataFrame"){
        r1<-smapdata
    }
    else{
        stop("Input Format incorrect")
    }
    ufam <- as.character(unique(r$SVIdentifier))
    famid <- c()
    for (ii in seq_len(length(ufam)))
    {
        stt <- strsplit(ufam[[ii]][1], split = "[.]")
        famid <- c(famid, as.character(stt[[1]][1]))
    }
    datf1 <- data.frame(table(famid))
    ha <- hash()
    .set(ha, keys = as.character(datf1$famid), 
        values = as.character(datf1$Freq))
    ## Checking Sex male/female and assigning chromosome number accordingly
    chro1 <- (unique(r1$RefcontigID1))
    #chro1 <- c(1:chro)
    dataFinal <- c()
    for (ii in seq_along(chro1))
    {
        print(paste('Chrom:',ii,sep=''))
        dat <- r[which(r$RefcontigID1 == ii), ]
        
        # variantType1<-dat$variantsubtype Changing the variant terms in DGV to
        # match svmap
        variantType1 <- dat$Type
        BSPQI_status_DB <- as.character(dat$Found_in_self_BSPQI_molecules)
        BSSSI_status_DB <- as.character(dat$Found_in_self_BSSSI_molecules)
        ## Extracting data from SVmap
        dat1 <- r1[which(r1$RefcontigID1 == chro1[ii]), ]
        chromo2<-dat1$RefcontigID2
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
        conf1 <- dat1$Confidence
        variantType2 <- as.character(dat1$Type)
        ###Parents Zygosity check
        'win_indel_parents<-as.numeric(win_indel_parents[oo])
        win_inv_trans_parents<- as.numeric(win_inv_trans_parents[oo])
        perc_similarity_parents<-as.numeric(perc_similarity_parents[oo])'
        #win_indel_parents<-as.numeric(win_indel_parents)
        #win_inv_trans_parents<- as.numeric(win_inv_trans_parents)
        'win_inv_trans_parents<- as.numeric(win_inv_trans_parents[oo])'
        #perc_similarity_parents<-as.numeric(perc_similarity_parents)
        rf_wb_ind_parents <- rf - win_indel_parents
        rf_fb_ind_parents <- rf + win_indel_parents
        rf_wb_int_parents <- rf - win_inv_trans_parents
        rf_fb_int_parents <- rf + win_inv_trans_parents
        #re <- dat1$RefEndPos
        re_wf_ind_parents <- re + win_indel_parents
        re_wb_ind_parents <- re - win_indel_parents
        re_wb_int_parents <- re - win_inv_trans_parents
        re_fb_int_parents <- re + win_inv_trans_parents
        # svfam<-as.character(dat1$SVIdentifier)
        if((length(grep("\\\\",smap))>=1)){
            spl<-strsplit(as.character(smap), split = "\\\\")
            lenspll1<-length(spl[[1]])
            spl1 <- strsplit(as.character(spl[[1]][lenspll1]), split = "_")
            # svfamid<-as.character(spl1[[1]][2])
            spl2 <- strsplit(as.character(spl1[[1]][1]), split = "[.]")
            patID <- spl2[[1]][2]
            svfamid <- spl2[[1]][1]
        }
        else if((length(grep("/",smap))>=1)){
            spl<-strsplit(as.character(smap), split = "/")
            lenspll1<-length(spl[[1]])
            spl1 <- strsplit(as.character(spl[[1]][lenspll1]), split = "_")
            # svfamid<-as.character(spl1[[1]][2])
            spl2 <- strsplit(as.character(spl1[[1]][1]), split = "[.]")
            patID <- spl2[[1]][2]
            svfamid <- spl2[[1]][1]
        }
        else{
            spl1 <- strsplit(as.character(smap), split = "_")
            # svfamid<-as.character(spl1[[1]][2])
            spl2 <- strsplit(as.character(spl1[[1]][1]), split = "[.]")
            patID <- spl2[[1]][2]
            svfamid <- spl2[[1]][1]
        }
        # svind<-dat1$SVIndex
        BSPQI_status_Query <- as.character(dat1$Found_in_self_BSPQI_molecules)
        BSSSI_status_Query <- as.character(dat1$Found_in_self_BSSSI_molecules)
        bsssi<-c()
        bspqi<-c()
        rest<-c()
        # conf<-dat$Confidence countfre<-0 percn<-c()
        datf <- c()
        for (nn in seq_along(rf))
        {
            
            #print(paste("nn:",nn)) 
            
            if ((variantType2[nn] == "deletion" | 
                variantType2[nn] == "insertion")){       
                        
                dat2 <- dat[which(((dat$RefStartPos <= rf[nn] & 
                    dat$RefStartPos >= rf_wb_ind[nn]) | 
                    (dat$RefStartPos >= rf[nn] & dat$RefStartPos <= 
                    rf_fb_ind[nn])) & ((dat$RefEndPos >= re[nn] & 
                    dat$RefEndPos <= re_wf_ind[nn]) |
                    (dat$RefEndPos <= re[nn] & dat$RefEndPos >= re_wb_ind[nn]))
                    ), ]
                    
                size1 <- size_bn[nn]
                ## Calculating Internal Frequency
                if (nrow(dat2) > 1)
                {
                    countfre <- c();countfreunfilt<-c()
                    svfam1 <- dat2$SVIdentifier
                    # svind1<-dat2$SVIndex
                    sv1 <- strsplit(as.character(svfam1), split = "_")
                    size_internal <- dat2$Size
                    zygo <- as.character(dat2$Zygosity)
                    type <- as.character(dat2$Type)
                    typ1 <- as.character(dat2$Type1)
                    typ2 <- as.character(dat2$Type2)
                    #svid <- as.numeric(dat2$SVIndex)
                    if (length(grep("SVIndex",names(smap)))>0){
                        svid <- dat2$SVIndex
                    }
                    else{
                        svid <- dat2$SmapEntryID
                    }
                    BSPQI_status_DB <- 
                        as.character(dat2$Found_in_self_BSPQI_molecules)
                    BSSSI_status_DB <- 
                        as.character(dat2$Found_in_self_BSSSI_molecules)
                    svStrt<-as.numeric(dat2$RefStartPos)
                    svEnd<-as.numeric(dat2$RefEndPos)
                    type1<-as.character(dat2$Type1)
                    type2<-as.character(dat2$Type2)
                    conf <- dat2$Confidence
                    motherZygosity <- c()
                    fatherZygosity <- c()
                    motherZygosity_exact <- c()
                    fatherZygosity_exact <- c()
                    svfamid3 <- c()
                    svfamid3_unfilt <- c()
                    bsssi<-c()
                    bs_bq<-c()
                    bspqi<-c()
                    svSAMP<-c()
                    svSAMP_unfiltered<-c()
                    homozygo<-c()
                    for (ll in seq_len(nrow(dat2))){
                        perc <- (size1/size_internal[ll])
                        #### print(perc) print(type[ll])
                        stt <- strsplit(sv1[[ll]][1], split = "[.]")
                        patID1 <- stt[[1]][2]
                        svfamid1 <- stt[[1]][1]
                        #### print(svfamid1)
                        ###Filtered INDEL
                        ##Check whether the parents have same breakpoints 
                        ##as the proband
                    
                        
                        if (((svStrt[ll]==rf[nn]) & (svEnd[ll]==re[nn])) & 
                            (perc >= perc_similarity) & 
                            identical(type[ll], variantType2[nn]) & 
                            (identical(svfamid1, svfamid))){
                            # Family Extra column zygosity
                            if (patID == 1 & patID1 == 2){
                          
                                # svfamid3<-c(svfamid3,svfamid1)
                                motherZygosity_exact <- c(motherZygosity_exact,
                                    as.character(zygo[ll]))
                                fatherZygosity_exact <- c(fatherZygosity_exact, 
                                    NULL)
                            } 
                            else if (patID == 1 & patID1 == 3){
                                
                                motherZygosity_exact <- c(motherZygosity_exact,
                                    NULL)
                                fatherZygosity_exact <- c(fatherZygosity_exact,
                                    as.character(zygo[ll]))
                            }   
                            else if (patID == patID1){
                                motherZygosity_exact <- c(motherZygosity_exact,
                                    NULL)
                                fatherZygosity_exact <- c(fatherZygosity_exact,
                                    NULL)
                            } 
                            else{
                                # svfamid3<-c(svfamid3,svfamid1)
                          
                                motherZygosity_exact <- c(motherZygosity_exact,
                                    NULL)
                                fatherZygosity_exact <- c(fatherZygosity_exact,
                                    NULL)
                            }
                        } 
                        else if ((perc >= perc_similarity) & 
                            identical(type[ll], variantType2[nn]) & 
                            (identical(svfamid1, svfamid)) & 
                            ((( svStrt[ll] <= rf[nn] & 
                            svStrt[ll] >= rf_wb_ind_parents[nn]) | 
                            (svStrt[ll] >= rf[nn] & 
                            svStrt[ll] <= rf_fb_ind_parents[nn])) & 
                            ((svEnd[ll] >= re[nn] & 
                            svEnd[ll] <= re_wf_ind_parents[nn]) |
                            (svEnd[ll] <= re[nn] & 
                            svEnd[ll] >= re_wb_ind_parents[nn])))
                            ){
                            #Family Extra column zygosity
                                if (patID == 1 & patID1 == 2){
                        
                                    # svfamid3<-c(svfamid3,svfamid1)
                                    motherZygosity <- c(motherZygosity, 
                                        as.character(zygo[ll]))
                                    fatherZygosity <- c(fatherZygosity, NULL)
                                } 
                                else if (patID == 1 & patID1 == 3){
                                    
                                    motherZygosity <- c(motherZygosity, NULL)
                                    fatherZygosity <- c(fatherZygosity, 
                                        as.character(zygo[ll]))
                                } 
                                else if (patID == patID1){
                         
                                    motherZygosity <- c(motherZygosity, NULL)
                                    fatherZygosity <- c(fatherZygosity, NULL)
                                } 
                                else {
                                    # svfamid3<-c(svfamid3,svfamid1)
                        
                                    motherZygosity <- c(motherZygosity, NULL)
                                    fatherZygosity <- c(fatherZygosity, NULL)
                                }
                            }
                            else if ((perc >= perc_similarity) & 
                                (identical(type[ll], variantType2[nn])) &
                                ((size_internal[ll]>=limsize)) & 
                                !(identical(svfamid1, svfamid)) & 
                                (conf[ll] >= indelconf) & 
                                ((BSPQI_status_DB[ll] == "yes" & 
                                BSSSI_status_DB[ll] == "yes") 
                                | (BSPQI_status_DB[ll] == "no" & 
                                BSSSI_status_DB[ll] == "yes") | 
                                (BSPQI_status_DB[ll] == "yes" & 
                                BSSSI_status_DB[ll] == "no") | 
                                (BSPQI_status_DB[ll] == "yes" & 
                                BSSSI_status_DB[ll] == "-") | 
                                (BSPQI_status_DB[ll] == "-" & 
                                BSSSI_status_DB[ll] == "yes")) & 
                                ((typ1[ll]=="insertion" & 
                                typ2[ll]=="insertion") | 
                                (typ1[ll]=="insertion" & is.na(typ2[ll])) |
                                (is.na(typ1[ll]) & typ2[ll]=="insertion") | 
                                (typ1[ll]=="deletion" & typ2[ll]=="deletion") |
                                (typ1[ll]=="deletion" & is.na(typ2[ll])) | 
                                (is.na(typ1[ll]) & typ2[ll]=="deletion"))
                                ) {
                                #print ("Not in the same family")
                                countfre <- c(countfre,1)
                                svfamid3 <- c(svfamid3, as.character(
                                    sv1[[ll]][1]))
                                motherZygosity <- c(motherZygosity, "-")
                                fatherZygosity <- c(fatherZygosity, "-")
                                svSAMP<-c(svSAMP,ll)
                                if(as.character(zygo[ll])=="homozygous"){
                                    homozygo<-c(homozygo, 
                                    as.character((sv1[[ll]][1])))
                                }
                                else{
                                    homozygo<-c(homozygo, NULL)
                                }
                            }
                            else{
                                
                                motherZygosity <- c(motherZygosity, "-")
                                fatherZygosity <- c(fatherZygosity, "-")     
                            }
                            ###Unfiltered INDEL
                           if (perc >= perc_similarity & 
                                identical(type[ll], variantType2[nn]) &
                                !(identical(svfamid1, svfamid))){
                                countfreunfilt <- c(countfreunfilt,1)
                                svfamid3_unfilt <- c(svfamid3_unfilt, 
                                    as.character(sv1[[ll]][1]))
                                svSAMP_unfiltered<-c(svSAMP_unfiltered,ll)
                            }
                            else{
                               print ("SVs does not meet the unfiltered 
                                    criteria!!!")
                            }
                      
                     
                    
                        }
                  
                    ##Calculating filtered Frequency INDEL
                    ##Filtration
                    if ((length(svfamid3) >= 1) & length(countfre) >= 1){ 
                        zyg<-as.character(dat2$Zygosity[svSAMP])
                        dat_temp <- data.frame(svfamid3, countfre,zyg)
                        countfre1 <- 0
                        if (nrow(dat_temp) > 0){
                            #### print(paste('dat_temp',dim(dat_temp)))
                            usampid <- as.character(unique(svfamid3))
                      
                        for (u in seq_len(length(usampid))){
                            dat_temp1 <- dat_temp[which(
                                dat_temp$svfamid3 %in% usampid[u]), ]
                            zygo<-as.character(dat_temp1$zyg)
                            ct <- sum(dat_temp1$countfre)
                            if (ct >= 1 & nrow (dat_temp1)>= 1) {
                                if(length(unique(as.character(zygo)))==1){
                                    if(unique(as.character(zygo) )== 
                                        "homozygous")
                                    {
                                        countfre1 <- countfre1 + 2
                                    }
                                    else if(unique(as.character(zygo))==
                                        "Unknown"){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else {
                                        countfre1 <- countfre1 + 1
                                    }
                                }
                                else{
                                    g1<-grep("homozygous",as.character(zygo))
                                    g2<-grep("heterozygous",as.character(zygo))
                                    g3<-grep("unknown",as.character(zygo))  
                                    if(length(g1)>=1 & length(g2)>=1 
                                        & length(g3)>=1){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else if(length(g1)==0 & length(g2)>=1 & 
                                        length(g3)>=1){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else if(length(g1)>=1 & length(g2)>=1 & 
                                        length(g3)==0){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else if(length(g1)>=1 & length(g2)==0 & 
                                        length(g3)>=1){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else{
                                        countfre1 <- countfre1 + 0
                                    }
                                }
                            }
                            else{
                                countfre1 <- countfre1 + 0
                            }
                        }
                    }
                    else{
                          countfre1 <- 0
                        }
                    }
                    else{
                        countfre1 <- 0
                    }
                    #### Unfiltered
                   if ((length(svfamid3_unfilt) >= 1) & 
                        length(countfreunfilt) >= 1)
                    {
                        zyg_unfiltered<-as.character(dat2$Zygosity[
                        svSAMP_unfiltered])
                        dat_temp <- data.frame(svfamid3_unfilt, 
                            countfreunfilt,zyg_unfiltered)
                        countfreunfilt1 <- 0
                        if (nrow(dat_temp) > 0){
                            #### print(paste('dat_temp',dim(dat_temp)))
                            usampid <- as.character(unique(svfamid3_unfilt))
                      
                            for (u in seq_len(length(usampid)))
                            {
                                dat_temp1 <- dat_temp[
                                    which
                                    (dat_temp$svfamid3_unfilt %in% usampid[u]),
                                    ]
                                zygo<-as.character(dat_temp1$zyg_unfiltered)
                                ct <- sum(dat_temp1$countfreunfilt)
                                if (ct >= 1 & nrow (dat_temp1)>= 1)
                                {
                                    if(length(unique(as.character(zygo))) == 1)
                                    {
                                        if(unique(as.character(zygo)) == 
                                            "homozygous"){
                                            countfreunfilt1 <- countfreunfilt1  
                                            + 2
                                        }
                                        else if(unique(as.character(zygo)) == 
                                        "Unknown"){
                                            countfreunfilt1 <- countfreunfilt1 
                                            + 2
                                        }
                                        else {
                                            countfreunfilt1 <- countfreunfilt1 
                                            + 1
                                        }
                                    }
                                    else{
                                        g1<-grep("homozygous",as.character(zygo)
                                            )
                                        g2<-grep("heterozygous",as.character(
                                            zygo))   
                                        g3<-grep("unknown",as.character(zygo))  
                                        if(length(g1)>=1 & length(g2)>=1 
                                            & length(g3)>=1){
                                            countfreunfilt1 <- countfreunfilt1 
                                            + 2
                                        }
                                        else if(length(g1)==0 & 
                                            length(g2)>=1 & 
                                            length(g3)>=1){
                                            countfreunfilt1 <- countfreunfilt1 
                                            + 2
                                        }
                                        else if(length(g1)>=1 & length(g2)>=1 & 
                                            length(g3)==0){
                                            countfreunfilt1 <- countfreunfilt1 
                                            + 2
                                        }
                                        else if(length(g1)>=1 & length(g2)==0 &
                                            length(g3)>=1){
                                            countfreunfilt1 <- countfreunfilt1  
                                            + 2
                                        }
                                        else{
                                            countfreunfilt1 <- countfreunfilt1 
                                            + 0
                                        }
                                    }
                                }
                                else{
                                   countfreunfilt1 <- countfreunfilt1 + 0
                                }
                            }
                        }   
                        else{
                            countfreunfilt1 <- 0
                            }
                        }
                        else{
                            countfreunfilt1 <- 0
                        }
                  ###Father Zygosity
                        fatherZygosity1 <- gsub("-", NA, fatherZygosity)
                        fatherZygosity1 <- as.character(na.omit(fatherZygosity1)
                            )
                        if(length(fatherZygosity_exact)>0){
                            fatherZygosity<-as.character(unique
                            (fatherZygosity_exact))
                        } 
                        else {
                            fatherZygosity<-as.character(unique
                            (fatherZygosity1))
                        }
                        if (length(fatherZygosity) > 1){
                            fatherZygosity <- paste(as.character
                            (fatherZygosity), collapse = ",")
                        } 
                        else if (length(fatherZygosity) == 0){
                            fatherZygosity <- "-"
                        } 
                        else{
                            fatherZygosity <- as.character(fatherZygosity)
                        }
                        ###Mother Zygosity
                        motherZygosity1 <- gsub("-", NA, motherZygosity)
                        motherZygosity1 <- as.character(na.omit(motherZygosity1
                            )
                            )
                        if(length(motherZygosity_exact)>0){
                            motherZygosity<-as.character(unique(
                                motherZygosity_exact))
                        } 
                        else {
                            motherZygosity<-as.character(unique(
                                motherZygosity1))
                        }
                 
                        if (length(motherZygosity) > 1) { 
                            motherZygosity <- paste(as.character(
                                motherZygosity), collapse = ",")
                        } 
                        else if (length(motherZygosity) == 0){
                            motherZygosity <- "-"
                        } 
                        else {
                             motherZygosity <- as.character(motherZygosity)
                        }
                  
                        if(length(homozygo)>=1){
                            homozygo <- length(unique(homozygo))
                        }
                        else if(length(homozygo)==1){
                            homozygo <- length(homozygo)
                        }
                        else{
                            homozygo <- 0
                        }
                                        
                  
                       
                        famno <- as.numeric(unname(hash::values(ha, 
                            keys = svfamid)))
                        data1 <- data.frame(dat1[nn, ], 
                            Internal_Freq_Perc_Filtered = 
                            as.numeric(format(round(
                            (countfre1/(2*(usamp - famno))) * 100,digits=3)
                            , scientific = FALSE)),
                            Internal_Freq_Perc_Unfiltered = 
                            as.numeric
                            (format(round((countfreunfilt1/
                            (2*(usamp - famno))) * 100,digits=3), 
                            scientific = FALSE)),
                            Internal_Homozygotes=as.numeric(homozygo), 
                            MotherZygosity = as.character(motherZygosity), 
                            FatherZygosity = as.character(fatherZygosity),  
                            stringsAsFactors = FALSE)
                        
                        datf <- rbind(datf, data1)
                }  
                else if (nrow(dat2) == 1){
                    # dgv_match=TRUE Calculating percentage similarity
                    countfre <- 0;countfreunfilt<-0
                    conf <- dat2$Confidence
                    size_internal <- dat2$Size
                    BSPQI_status_DB <- as.character(
                        dat2$Found_in_self_BSPQI_molecules)
                    BSSSI_status_DB <- as.character(
                        dat2$Found_in_self_BSSSI_molecules)
                    perc <- (size1/size_internal)
                    svv <- strsplit(as.character(svfam1), split = "_")
                    zygo <- as.character(dat2$Zygosity)
                    stt <- strsplit(svv[[1]][1], split = "[.]")
                    patID1 <- stt[[1]][2]
                    svfamid1 <- stt[[1]][1]
                    motherZygosity <- ""
                    fatherZygosity <- ""
                    type <- as.character(dat2$Type)
                    typ1 <- as.character(dat2$Type1)
                    typ2 <- as.character(dat2$Type2)
                    homozygo<-c()
                    ##Check for parents
                    if (perc >= perc_similarity & 
                        identical(type, variantType2[nn]) & 
                        (identical(svfamid1, svfamid)) &
                        (((dat2$RefStartPos <= rf[nn] & 
                        dat2$RefStartPos >= rf_wb_ind_parents[nn]) | 
                        (dat2$  RefStartPos >= rf[nn] & 
                        dat2$RefStartPos <= rf_fb_ind_parents[nn])) &
                        ((dat2$RefEndPos >= re[nn] & 
                        dat2$RefEndPos <= re_wf_ind_parents[nn]) | 
                        (dat2$RefEndPos <= re[nn] & 
                        dat2$RefEndPos >= re_wb_ind_parents[nn])))){
                        # Family Extra column zygosity
                        if (patID == 1 & patID1 == 2){
                            motherZygosity <- as.character(zygo)
                            fatherZygosity <- "-"
                        } 
                        else if (patID == 1 & patID1 == 3){
                            motherZygosity <- "-"
                            fatherZygosity <- as.character(zygo)
                        } 
                        else if (patID == patID1){
                            motherZygosity <- "-"
                            fatherZygosity <- "-"
                        } 
                        else{
                            motherZygosity <- "-"
                            fatherZygosity <- "-"
                        }
                    }
                    else if (perc >= perc_similarity & 
                        identical(type, variantType2[nn]) & 
                        !(identical(svfamid1, svfamid)) & 
                        ((size_internal>=limsize)) & 
                        ((BSPQI_status_DB == "yes" & BSSSI_status_DB == "yes") 
                        | (BSPQI_status_DB == "no" & BSSSI_status_DB == "yes") 
                        | (BSPQI_status_DB == "yes" & BSSSI_status_DB == "no")
                        | (BSPQI_status_DB == "yes" & BSSSI_status_DB == "-") |
                        (BSPQI_status_DB == "-" & BSSSI_status_DB == "yes")) & 
                        (conf > indelconf) & 
                        ((typ1=="insertion" & typ2=="insertion") | 
                        (typ1=="insertion" & is.na(typ2)) | 
                        (is.na(typ1) & typ2=="insertion") | 
                        (typ1=="deletion" & typ2=="deletion") |
                        (typ1=="deletion" & is.na(typ2)) | 
                        (is.na(typ1) & typ2=="deletion")) ){
                        #print("not in a family")
                        if(as.character(zygo)=="homozygous"){
                            countfre <- 2
                            if(as.character(zygo)=="homozygous"){
                                homozygo<- as.character(svfamid1)
                            }
                            else{
                                homozygo<- NULL
                            }
                        }
                        else if(zygo=="unknown"){
                            countfre <- 2
                        }
                        else { 
                            countfre <- 1
                        }
                        motherZygosity <- "-"
                        fatherZygosity <- "-"
                    } 
                    else{ 
                        countfre <- 0
                        motherZygosity <- "-"
                        fatherZygosity <- "-"
                    }
                        ###Unfiltered
                    if (perc >= perc_similarity & 
                        identical(type, variantType2[nn]) & 
                        !(identical(svfamid1, svfamid))){                 
                        if(zygo=="homozygous"){
                            countfreunfilt <- 2
                        }
                        else if(zygo=="unknown"){
                            countfreunfilt <- 2
                        }
                        else{ 
                            countfreunfilt <- 1
                        }
                    } 
                    else{
                        countfreunfilt <- 0
                    }
                    
                    if(length(homozygo)>0){
                        homozygo <- length(unique(homozygo))
                    }
                    else{
                        homozygo=0
                    }
                  
                famno <- as.numeric(unname(hash::values(ha, keys = svfamid)))
                data1 <- data.frame(dat1[nn, ], 
                    Internal_Freq_Perc_Filtered = 
                    as.numeric(format(round((countfre/(2*(usamp - famno))) 
                    * 100,digits=3), scientific = FALSE)),
                    Internal_Freq_Perc_Unfiltered = 
                    as.numeric(format
                    (round((countfreunfilt/(2*(usamp - famno))) 
                    * 100,digits=3), scientific = FALSE)), 
                    Internal_Homozygotes=as.numeric(homozygo), 
                    MotherZygosity = as.character(motherZygosity), 
                    FatherZygosity = as.character(fatherZygosity),  
                    stringsAsFactors = FALSE)
                  
                datf <- rbind(datf, data1)
                } 
                else
                {
                    # dgv_match=FALSE
                    data1 <- data.frame(dat1[nn, ], 
                    Internal_Freq_Perc_Filtered = 0, 
                    Internal_Freq_Perc_Unfiltered = 0, 
                    Internal_Homozygotes= 0, 
                    MotherZygosity = "-", 
                    FatherZygosity = "-", 
                    stringsAsFactors = FALSE)
                    #### print(dim(data1)) print(names(data1)[56:58])
                    #### print(names(datf)[56:58])
                    #### print(identical((names(data1)[56]),(names(datf)[56])))
                    #### print(identical((names(data1)[57]),(names(datf)[57])))
                    #### print(identical((names(data1)[58]),(names(datf)[58])))
                    datf <- rbind(datf, data1)
                    # next;
                }
                
            }
            else if ((variantType2[nn] == "duplication" |
                    variantType2[nn] == "duplication_split" | 
                    variantType2[nn] == "duplication_inverted" |
                    variantType2[nn] == "insertion"))
                {   
                    dat2 <- dat[which((dat$RefStartPos <= rf[nn]&
                        dat$RefStartPos >= rf_wb_ind[nn] | 
                        dat$RefStartPos >= rf[nn] & 
                        dat$RefStartPos <= rf_fb_ind[nn]) & 
                        (dat$RefEndPos >= re[nn] & 
                        dat$RefEndPos <= re_wf_ind[nn] | 
                        dat$RefEndPos <= re[nn] & 
                        dat$RefEndPos >= re_wb_ind[nn])), ]
                    size1 <- size_bn[nn]
                    ## Calculating Internal Frequency
                    if (nrow(dat2) > 1){ 
                        #print("1")
                        countfre <- c();countfreunfilt<-c()
                        svfam1 <- dat2$SVIdentifier
                        # svind1<-dat2$SVIndex
                        sv1 <- strsplit(as.character(svfam1), split = "_")
                        size_internal <- dat2$Size
                        zygo <- as.character(dat2$Zygosity)
                        type <- as.character(dat2$Type)
                        typ1 <- as.character(dat2$Type1)
                        typ2 <- as.character(dat2$Type2)
                        BSPQI_status_DB <- 
                        as.character(dat2$Found_in_self_BSPQI_molecules)
                        BSSSI_status_DB <- 
                        as.character(dat2$Found_in_self_BSSSI_molecules)
                        BSPQI_chimeric_score_DB <- 
                        as.character(dat2$Fail_BSPQI_assembly_chimeric_score)
                        BSSSI_chimeric_score_DB <- 
                        as.character(dat2$Fail_BSSSI_assembly_chimeric_score)
                        motherZygosity_exact <- c()
                        fatherZygosity_exact <- c()
                        conf <- dat2$Confidence
                        svStrt<-as.numeric(dat2$RefStartPos)
                        svEnd<-as.numeric(dat2$RefEndPos)
                        motherZygosity <- c()
                        fatherZygosity <- c()
                        svfamid3 <- c()
                        svfamid3_unfilt <- c()
                        svSAMP_unfiltered<-c()
                        svSAMP<-c()
                        homozygo<-c()
                        for (ll in seq_len(nrow(dat2))){
                            #print(ll)
                            #### print(perc) print(type[ll])
                            perc <- (size1/size_internal[ll])
                            stt <- strsplit(sv1[[ll]][1], split = "[.]")
                            patID1 <- stt[[1]][2]
                            svfamid1 <- stt[[1]][1]
                            #### print(svfamid1)
                            ####Checking for Duplication
                            if (((svStrt[ll]==rf[nn]) & 
                                (svEnd[ll]==re[nn])) &  
                                (identical(type[ll], "duplication") |
                                identical(type[ll], "duplication_split")| 
                                identical(type[ll], "duplication_inverted")) & 
                                (identical(svfamid1, svfamid))) {
                                # Family Extra column zygosity
                                if (patID == 1 & patID1 == 2) {                        
                                # svfamid3<-c(svfamid3,svfamid1)
                                    motherZygosity_exact <- 
                                    c(motherZygosity_exact, as.character(
                                        zygo[ll]))
                                    fatherZygosity_exact <- 
                                    c(fatherZygosity_exact, NULL)
                                } 
                                else if (patID == 1 & patID1 == 3){
                                
                                    motherZygosity_exact <- 
                                    c(motherZygosity_exact, NULL)
                                    fatherZygosity_exact <- 
                                    c(fatherZygosity_exact, 
                                    as.character(zygo[ll]))
                                }
                                else if (patID == patID1) {                   
                                    motherZygosity_exact <- 
                                    c(motherZygosity_exact, NULL)
                                    fatherZygosity_exact <- 
                                    c(fatherZygosity_exact, NULL)
                                } 
                                else {
                                    # svfamid3<-c(svfamid3,svfamid1)
                        
                                    motherZygosity_exact <- 
                                    c(motherZygosity_exact, NULL)
                                    fatherZygosity_exact <- 
                                    c(fatherZygosity_exact, NULL)
                                }
                            } 
                            else if (((( svStrt[ll] <= rf[nn] & 
                                svStrt[ll] >= rf_wb_ind_parents[nn]) |
                                (svStrt[ll] >= rf[nn] & 
                                svStrt[ll] <= rf_fb_ind_parents[nn])) & 
                                ((svEnd[ll] >= re[nn] & 
                                svEnd[ll] <= re_wf_ind_parents[nn]) | 
                                (svEnd[ll] <= re[nn] & 
                                svEnd[ll] >= re_wb_ind_parents[nn]))) &
                                (identical(type[ll], "duplication") |
                                identical(type[ll], "duplication_split")| 
                                identical(type[ll], "duplication_inverted")) & 
                               (identical(svfamid1, svfamid))){
                                # Family Extra column zygosity
                                if (patID == 1 & patID1 == 2){
                        
                                   # svfamid3<-c(svfamid3,svfamid1)
                                    motherZygosity <- 
                                    c(motherZygosity, as.character(zygo[ll]))
                                    fatherZygosity <- 
                                    c(fatherZygosity, NULL)
                                } 
                                else if (patID == 1 & patID1 == 3){
                                   
                                    motherZygosity <- c(motherZygosity, NULL)
                                    fatherZygosity <- 
                                    c(fatherZygosity, as.character(zygo[ll]))
                                } 
                                else if (patID == patID1){
                   
                                    motherZygosity <- c(motherZygosity, NULL)
                                    fatherZygosity <- c(fatherZygosity, NULL)
                                }
                                else {
                                    # svfamid3<-c(svfamid3,svfamid1)
                        
                                    motherZygosity <- c(motherZygosity, NULL)
                                    fatherZygosity <- c(fatherZygosity, NULL)
                                }
                            }
                            else if (((identical(type[ll], "duplication") |
                                identical(type[ll], "duplication_split") | 
                                identical(type[ll], "duplication_inverted"))) & 
                                !(identical(svfamid1, svfamid)) & 
                                ((BSPQI_status_DB[ll] == "yes" & 
                                BSSSI_status_DB[ll] == "yes") | 
                                (BSPQI_status_DB[ll] == "no" & 
                                BSSSI_status_DB[ll] == "yes") |
                                (BSPQI_status_DB[ll] == "yes" & 
                                BSSSI_status_DB[ll] == "no") | 
                                (BSPQI_status_DB[ll] == "yes" & 
                                BSSSI_status_DB[ll] == "-") | 
                                (BSPQI_status_DB[ll] == "-" & 
                                BSSSI_status_DB[ll] == "yes")) & 
                                ((BSPQI_chimeric_score_DB[ll] == "pass" &
                                BSSSI_chimeric_score_DB[ll] == "pass") |
                                (BSPQI_chimeric_score_DB[ll] == "fail" & 
                                BSSSI_chimeric_score_DB[ll] == "pass") | 
                                (BSPQI_chimeric_score_DB[ll] == "pass" & 
                                BSSSI_chimeric_score_DB[ll] == "fail") | 
                                (BSPQI_chimeric_score_DB[ll] == "pass" & 
                                BSSSI_chimeric_score_DB[ll] == "-") | 
                                (BSPQI_chimeric_score_DB[ll] == "-" & 
                                BSSSI_chimeric_score_DB[ll] == "pass")) &
                                ((typ1[ll]=="duplication" 
                                & typ2[ll]=="duplication") | 
                                (typ1[ll]=="duplication" & is.na(typ2[ll])) |
                                (is.na(typ1[ll]) & typ2[ll]=="duplication") | 
                                (typ1[ll]=="duplication_split" & 
                                typ2[ll]=="duplication_split") | 
                                (typ1[ll]=="duplication_split" & 
                                is.na(typ2[ll])) | 
                                (is.na(typ1[ll]) & 
                                typ2[ll]=="duplication_split") | 
                                (typ1[ll]=="duplication_inverted" & 
                                typ2[ll]=="duplication_inverted") | 
                                (typ1[ll]=="duplication_inverted" & 
                                is.na(typ2[ll])) | 
                                (is.na(typ1[ll]) & 
                                typ2[ll]=="duplication_inverted"))){
                                    #print ("HI!")
                                    countfre <- c(countfre,1)
                                    svfamid3 <- c(svfamid3,
                                    as.character(sv1[[ll]][1]))
                                    motherZygosity <- c(motherZygosity, "-")
                                    fatherZygosity <- c(fatherZygosity, "-")
                                    svSAMP<-c(svSAMP,ll)
                                    if(as.character(zygo[ll])=="homozygous"){
                                        homozygo<-c(homozygo, 
                                        as.character(sv1[[ll]][1]))
                                    }
                                    else{
                                        homozygo<-c(homozygo, NULL)
                                    }
                                }
                                else if (((svStrt[ll]==rf[nn]) & 
                                    (svEnd[ll]==re[nn]))&
                                    (identical(type[ll], "insertion")) & 
                                    (identical(svfamid1, svfamid))){
                                        # Family Extra column zygosity
                                        if (patID == 1 & patID1 == 2){
                        
                                            # svfamid3<-c(svfamid3,svfamid1)
                                            motherZygosity_exact <- 
                                            c(motherZygosity_exact,
                                            as.character(zygo[ll]))
                                            fatherZygosity_exact <- 
                                            c(fatherZygosity_exact, NULL)
                                        } 
                                        else if (patID == 1 & patID1 == 3) {
                                            
                                            motherZygosity_exact <- 
                                            c(motherZygosity_exact,NULL)
                                            fatherZygosity_exact <- 
                                            c(fatherZygosity_exact, 
                                            as.character(zygo[ll]))
                                        } 
                                        else if (patID == patID1){
                                            motherZygosity_exact <- 
                                            c(motherZygosity_exact,NULL)
                                            fatherZygosity_exact <- 
                                            c(fatherZygosity_exact, NULL)
                                        } 
                                        else {
                                            # svfamid3<-c(svfamid3,svfamid1)
                        
                                            motherZygosity_exact <- c(
                                                motherZygosity_exact,NULL)
                                            fatherZygosity_exact <- c(
                                                fatherZygosity_exact, NULL)
                                        }
                                    } 
                                    else if ((identical(type[ll], "insertion")) 
                                        & ((( svStrt[ll] <= rf[nn] & 
                                        svStrt[ll] >= rf_wb_ind_parents[nn]) | 
                                        (svStrt[ll] >= rf[nn] & 
                                        svStrt[ll] <= rf_fb_ind_parents[nn])) &
                                        ((svEnd[ll] >= re[nn] & 
                                        svEnd[ll] <= re_wf_ind_parents[nn]) |
                                        (svEnd[ll] <= re[nn] & 
                                        svEnd[ll] >= re_wb_ind_parents[nn]))) & 
                                        (perc >= perc_similarity) & 
                                        (identical(svfamid1, svfamid)) ) {
                                            # Family Extra column zygosity
                                            if (patID == 1 & patID1 == 2){
                        
                                                # svfamid3<-c(svfamid3,svfamid1)
                                                motherZygosity <- 
                                                c(motherZygosity, 
                                                as.character(zygo[ll]))
                                                fatherZygosity <- c(
                                                fatherZygosity, NULL)
                                            } 
                                            else if (patID == 1 & patID1 == 3){
                                                
                                                motherZygosity <- c(
                                                motherZygosity, NULL)
                                                fatherZygosity <- c(
                                                fatherZygosity, 
                                                as.character(zygo[ll]))
                                            } 
                                            else if (patID == patID1){
                       
                                            motherZygosity <- c(
                                                motherZygosity, 
                                                NULL
                                                )
                                            fatherZygosity <- c(
                                                fatherZygosity, 
                                                NULL)
                                            }
                                            else{
                                                # svfamid3<-c(svfamid3,svfamid1)
                        
                                                motherZygosity <- c(
                                                    motherZygosity,
                                                    NULL)
                                                fatherZygosity <- c(
                                                    fatherZygosity, 
                                                    NULL)
                                            }
                                        } 
                                        else if ((
                                            identical(type[ll], "insertion")) & 
                                            !(identical(svfamid1, svfamid)) & 
                                            ((size_internal[ll]>=limsize)) & 
                                            ((BSPQI_status_DB[ll] == "yes" & 
                                            BSSSI_status_DB[ll] == "yes") |
                                            (BSPQI_status_DB[ll] == "no" & 
                                            BSSSI_status_DB[ll] == "yes") |
                                            (BSPQI_status_DB[ll] == "yes" & 
                                            BSSSI_status_DB[ll] == "no") |
                                            (BSPQI_status_DB[ll] == "yes" & 
                                            BSSSI_status_DB[ll] == "-") |
                                            (BSPQI_status_DB[ll] == "-" & 
                                            BSSSI_status_DB[ll] == "yes")) & 
                                            (conf[ll] >= indelconf) & 
                                            (typ1[ll]=="insertion" & 
                                            typ2[ll]=="insertion") | 
                                            (typ1[ll]=="insertion" & 
                                            is.na(typ2[ll])) | 
                                            (is.na(typ1[ll]) & 
                                            typ2[ll]=="insertion")) {
                                                #print ("HI!")
                                                countfre <- c(countfre,1)
                                                svfamid3 <- c(svfamid3, 
                                                as.character(sv1[[ll]][1]))
                                                if(as.character(zygo[ll])==
                                                    "homozygous"){
                                                    homozygo<-c(homozygo,
                                                    as.character(sv1[[ll]][1]))
                                                }
                                                else{
                                                    homozygo<-c(homozygo, NULL)
                                                }
                                            motherZygosity <- c(motherZygosity,
                                            "-")
                                            fatherZygosity <- c(fatherZygosity, 
                                            "-")
                                            svSAMP<-c(svSAMP,ll)
                                        }
                                        else {
                                             
                                                motherZygosity <- c(
                                                    motherZygosity, "-")
                                                fatherZygosity <- c(
                                                    fatherZygosity, "-")
                                        } 
                                    if (((identical(type[ll], "duplication") |
                                        identical(type[ll], 
                                        "duplication_split") | 
                                        identical(type[ll], 
                                        "duplication_inverted") |
                                        identical(type[ll], "insertion")))
                                        & !(identical(svfamid1, svfamid))) {
                                            countfreunfilt <- c(
                                            countfreunfilt,1
                                            )
                                            svfamid3_unfilt <- c(
                                            svfamid3_unfilt, 
                                            as.character(sv1[[ll]][1]))
                                            svSAMP_unfiltered<-c(
                                            svSAMP_unfiltered,ll)
                                        } 
                                        else {
                                            print("Unfiltered criterion not 
                                                met!")
                                        }
                      
                     
                        }
                        ##Calculating filtered Frequency INDEL
                        ##Filtration
                 
                        if ((length(svfamid3) >= 1) & 
                            length(countfre) >= 1){ 
                            zyg<-as.character(dat2$Zygosity[svSAMP])
                            dat_temp <- data.frame(svfamid3, countfre,zyg)
                            countfre1 <- 0
                            if (nrow(dat_temp) > 0) {
                            #### print(paste('dat_temp',dim(dat_temp)))
                            usampid <- as.character(unique(svfamid3))
                      
                        for (u in seq_len(length(usampid))){
                                dat_temp1 <- dat_temp [which(
                                    dat_temp$svfamid3 %in% usampid[u]), 
                                    ]
                                zygo <- as.character(dat_temp1$zyg)
                                ct <- sum(dat_temp1$countfre)
                                if (ct >= 1 & nrow (dat_temp1)>= 1) {
                                    if(length(unique(as.character(zygo)))==1){
                                        if(unique(as.character(zygo))==
                                            "homozygous"){
                                        countfre1 <- countfre1 + 2
                                        }
                                        else if(unique(as.character(zygo))==
                                        "Unknown"){
                                            countfre1 <- countfre1 + 2
                                    }
                                    else {
                                        countfre1 <- countfre1 + 1
                                    }
                                }
                                else{
                                    g1<-grep("homozygous",as.character(zygo))
                                    g2<-grep("heterozygous",as.character(zygo))   
                                    g3<-grep("unknown",as.character(zygo))  
                                    if(length(g1)>=1 & length(g2)>=1 & 
                                        length(g3)>=1){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else if(length(g1)==0 & length(g2)>=1 & 
                                        length(g3)>=1){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else if(length(g1)>=1 & length(g2)>=1 & 
                                        length(g3)==0){
                                            countfre1 <- countfre1 + 2
                                    }
                                    else if(length(g1)>=1 & length(g2)==0 & 
                                        length(g3)>=1){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else{
                                        countfre1 <- countfre1 + 0
                                    }
                                }
                            }
                            else{
                                countfre1 <- countfre1 + 0
                            }
                        }
                    }
                    else
                        {
                          countfre1 <- 0
                        }
                    }
                    else{
                        countfre1 <- 0
                        }
                    #### Unfiltered
                   if ((length(svfamid3_unfilt) >= 1) & 
                        length(countfreunfilt) >= 1) {
                            zyg_unfiltered<-as.character(
                                dat2$Zygosity[svSAMP_unfiltered])
                             dat_temp <- data.frame(svfamid3_unfilt, 
                                 countfreunfilt,zyg_unfiltered)
                    countfreunfilt1 <- 0
                    if (nrow(dat_temp) > 0) {
                      #### print(paste('dat_temp',dim(dat_temp)))
                      usampid <- as.character(unique(svfamid3_unfilt))
                      
                        for (u in seq_len(length(usampid))) {
                            dat_temp1 <- dat_temp[which(
                                dat_temp$svfamid3_unfilt %in% usampid[u]), ]
                            zygo<-as.character(dat_temp1$zyg_unfiltered)
                            ct <- sum(dat_temp1$countfreunfilt)
                            if (ct >= 1 & nrow (dat_temp1)>= 1)
                            {
                                if(length(unique(as.character(zygo)))==1){
                                    if(unique(as.character(zygo)) == 
                                        "homozygous")
                                    {
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }
                                    else if(unique(as.character(zygo))==
                                        "Unknown"){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }
                                    else {
                                        countfreunfilt1 <- countfreunfilt1 + 1
                                    }
                                }
                                else{
                                    g1<-grep("homozygous",as.character(zygo))
                                    g2<-grep("heterozygous",as.character(zygo))   
                                    g3<-grep("unknown",as.character(zygo))  
                                    if(length(g1)>=1 & length(g2)>=1 &
                                        length(g3)>=1){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }
                                    else if(length(g1)==0 & length(g2)>=1 &
                                        length(g3)>=1){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }
                                    else if(length(g1)>=1 & length(g2)>=1 & 
                                        length(g3)==0){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }
                                    else if(length(g1)>=1 & length(g2)==0 & 
                                        length(g3)>=1){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }
                                    else{
                                        countfreunfilt1 <- countfreunfilt1 + 0
                                    }
                                }
                            }
                            else{
                                countfreunfilt1 <- countfreunfilt1 + 0
                            }
                        }
                    }
                    else{
                            countfreunfilt1 <- 0
                        }
                    }
                    else{
                        countfreunfilt1 <- 0
                    }
                  ###Father Zygosity
                    fatherZygosity1 <- gsub("-", NA, fatherZygosity)
                    fatherZygosity1 <- as.character(na.omit(fatherZygosity1))
                    if(length(fatherZygosity_exact)>0){
                        fatherZygosity<-as.character
                        (unique(fatherZygosity_exact))
                    } 
                    else {
                        fatherZygosity<-as.character(unique(fatherZygosity1))
                    }
                  
                    if (length(fatherZygosity) > 1){
                        fatherZygosity <- paste(as.character(fatherZygosity), 
                        collapse = ",")
                    } 
                    else if (length(fatherZygosity) == 0){
                        fatherZygosity <- "-"
                    } 
                    else {
                        fatherZygosity <- as.character(fatherZygosity)
                    }
                    ###Mother Zygosity
                    motherZygosity1 <- gsub("-", NA, motherZygosity)
                    motherZygosity1 <- as.character(na.omit(motherZygosity1))
                    if(length(motherZygosity_exact)>0){
                        motherZygosity<-as.character 
                        (unique(motherZygosity_exact))
                    } 
                    else {
                        motherZygosity<-as.character(unique(motherZygosity1))
                    }
                 
                    if (length(motherZygosity) > 1){
                        motherZygosity <- paste(
                        as.character(motherZygosity), collapse = ",")
                    } 
                    else if (length(motherZygosity) == 0) {
                        motherZygosity <- "-"
                    } 
                    else{
                        motherZygosity <- as.character(motherZygosity)
                    }
                    ###Proband Zygosity
                    if(length(homozygo)>=1){
                        homozygo <- length(unique(homozygo))
                    }
                    else if(length(homozygo)==1){
                        homozygo <- length(homozygo)
                    }
                    else{
                        homozygo <- 0
                    }
                  
                   
                    famno <- as.numeric(unname(hash::values(ha, 
                    keys = svfamid)))
                    data1 <- data.frame(dat1[nn, ], 
                    Internal_Freq_Perc_Filtered = as.numeric(
                    format(round((countfre1/(2*(usamp - famno))) 
                    * 100,digits=3), scientific = FALSE)),
                    Internal_Freq_Perc_Unfiltered = 
                    as.numeric(format(round
                    ((countfreunfilt1/(2*(usamp - famno))) * 100,digits=3),
                    scientific = FALSE)), 
                    Internal_Homozygotes=as.numeric(homozygo), 
                    MotherZygosity = as.character(motherZygosity), 
                    FatherZygosity = as.character(fatherZygosity), 
                    stringsAsFactors = FALSE)
                   
                    datf <- rbind(datf, data1)
                    }
                    else if (nrow(dat2) == 1){ 
                        # dgv_match=TRUE Calculating percentage similarity
                        countfre <- 0;countfreunfilt<-0
                        size_internal <- dat2$Size
                        conf <- dat2$Confidence
                        BSPQI_status_DB <- 
                        as.character(dat2$Found_in_self_BSPQI_molecules)
                        BSSSI_status_DB <- 
                        as.character(dat2$Found_in_self_BSSSI_molecules)
                        BSPQI_chimeric_score_DB <- 
                        as.character(dat2$Fail_BSPQI_assembly_chimeric_score)
                        BSSSI_chimeric_score_DB <- 
                        as.character(dat2$Fail_BSSSI_assembly_chimeric_score)
                        perc <- (size1/size_internal)
                        svv <- strsplit(as.character(svfam1), split = "_")
                        zygo <- as.character(dat2$Zygosity)
                        stt <- strsplit(svv[[1]][1], split = "[.]")
                        patID1 <- stt[[1]][2]
                        homozygo<-c()
                        svfamid1 <- stt[[1]][1]
                        motherZygosity <- ""
                        fatherZygosity <- ""
                        type <- as.character(dat2$Type)
                        typ1 <- as.character(dat2$Type1)
                        typ2 <- as.character(dat2$Type2)
                  #print(type)
                        if ((((dat2$RefStartPos <= rf[nn] & 
                        dat2$RefStartPos >= rf_wb_ind_parents[nn]) | 
                        (dat2$  RefStartPos >= rf[nn] & 
                        dat2$RefStartPos <= rf_fb_ind_parents[nn])) & 
                        ((dat2$RefEndPos >= re[nn] & 
                        dat2$RefEndPos <= re_wf_ind_parents[nn]) | 
                        (dat2$RefEndPos <= re[nn] & 
                        dat2$RefEndPos >= re_wb_ind_parents[nn]))) & 
                        ((identical(type, "duplication") |
                        identical(type, "duplication_split")| 
                        identical(type, "duplication_inverted"))) & 
                        (identical(svfamid1, svfamid)) )
                        {
                            'print(type)
                            print(svfamid1)
                            print(svfamid)
                            print(BSPQI_status_DB)
                            print(BSSSI_status_DB)
                            print(paste(ii,":",nn))'
                            # Family Extra column zygosity
                            if (patID == 1 & patID1 == 2){
                                motherZygosity <- as.character(zygo)
                                fatherZygosity <- "-"
                            } 
                            else if (patID == 1 & patID1 == 3){
                                motherZygosity <- "-"
                                fatherZygosity <- as.character(zygo)
                            } 
                            else if (patID == patID1){
                                motherZygosity <- "-"
                                fatherZygosity <- "-"
                            } 
                            else{
                                motherZygosity <- "-"
                                fatherZygosity <- "-"
                            }
                        } 
                        else if (((identical(type, "duplication") |
                        identical(type, "duplication_split")| 
                        identical(type, "duplication_inverted"))) & 
                        !(identical(svfamid1, svfamid)) & 
                        ((BSPQI_status_DB == "yes" & 
                        BSSSI_status_DB == "yes") | 
                        (BSPQI_status_DB == "no" & 
                        BSSSI_status_DB == "yes") | 
                        (BSPQI_status_DB == "yes" & 
                        BSSSI_status_DB == "no") | 
                        (BSPQI_status_DB == "yes" & 
                        BSSSI_status_DB == "-") | 
                        (BSPQI_status_DB == "-" & 
                        BSSSI_status_DB == "yes")) & 
                        ((BSPQI_chimeric_score_DB == "pass" & 
                        BSSSI_chimeric_score_DB == "pass") | 
                        (BSPQI_chimeric_score_DB == "fail" & 
                        BSSSI_chimeric_score_DB == "pass") | 
                        (BSPQI_chimeric_score_DB == "pass" & 
                        BSSSI_chimeric_score_DB == "fail") | 
                        (BSPQI_chimeric_score_DB == "pass" & 
                        BSSSI_chimeric_score_DB == "-") | 
                        (BSPQI_chimeric_score_DB == "-" & 
                        BSSSI_chimeric_score_DB == "pass")) & 
                        ((typ1=="duplication" & typ2=="duplication") |
                        (typ1=="duplication" & is.na(typ2)) | 
                        (is.na(typ1) & typ2=="duplication") | 
                        (typ1=="duplication_split" & 
                        typ2=="duplication_split") | 
                        (typ1=="duplication_split" & is.na(typ2)) |
                        (is.na(typ1) & typ2=="duplication_split") |
                        (typ1=="duplication_inverted" & 
                        typ2=="duplication_inverted") | 
                        (typ1=="duplication_inverted" & is.na(typ2)) |
                        (is.na(typ1) & typ2=="duplication_inverted"))){
                            if(zygo=="homozygous"){
                                countfre <- 2
                            }
                            else if(zygo=="unknown"){
                                countfre <- 2
                            }
                            else{ 
                                countfre <- 1
                            }
                            if(as.character(zygo)=="homozygous"){
                                homozygo <- as.character(svfamid1)
                            }
                            else{
                                homozygo<- NULL
                            }
                            motherZygosity <- "-"
                            fatherZygosity <- "-"
                            'motherZygosity <- "-"
                            fatherZygosity <- "-"'
                        } 
                        else if (((identical(type, "insertion"))) & 
                        (identical(svfamid1, svfamid)) & 
                        (((dat2$RefStartPos <= rf[nn] & 
                        dat2$RefStartPos >= rf_wb_ind_parents[nn]) | 
                        (dat2$  RefStartPos >= rf[nn] & 
                        dat2$RefStartPos <= rf_fb_ind_parents[nn])) &
                        ((dat2$RefEndPos >= re[nn] & 
                        dat2$RefEndPos <= re_wf_ind_parents[nn]) |
                        (dat2$RefEndPos <= re[nn] &
                        dat2$RefEndPos >= re_wb_ind_parents[nn]))) &
                        (perc >= perc_similarity)){
                            # Family Extra column zygosity
                            if (patID == 1 & patID1 == 2){
                                motherZygosity <- as.character(zygo)
                                fatherZygosity <- "-"
                            }
                            else if (patID == 1 & patID1 == 3){
                                motherZygosity <- "-"
                                fatherZygosity <- as.character(zygo)
                            } 
                            else if (patID == patID1){
                                motherZygosity <- "-"
                                fatherZygosity <- "-"
                            }
                            else {
                                motherZygosity <- "-"
                                fatherZygosity <- "-"
                            }
                        } 
                        else if (((identical(type, "insertion"))) & 
                        !(identical(svfamid1, svfamid)) & 
                        perc >= perc_similarity &
                        ((BSPQI_status_DB == "yes" &
                        BSSSI_status_DB == "yes") |
                        (BSPQI_status_DB == "no" & 
                        BSSSI_status_DB == "yes") |
                        (BSPQI_status_DB == "yes" & 
                        BSSSI_status_DB == "no") |
                        (BSPQI_status_DB == "yes" & 
                        BSSSI_status_DB == "-") |
                        (BSPQI_status_DB == "-" & 
                        BSSSI_status_DB == "yes")) & 
                        ((size_internal>=limsize)) & 
                        (conf > indelconf) & 
                        ((typ1=="insertion" & typ2=="insertion") |
                        (typ1=="insertion" & is.na(typ2)) | 
                        (is.na(typ1) & typ2=="insertion"))){
                            if(zygo=="homozygous"){
                                countfre <- 2
                            }
                            else if(zygo=="unknown"){
                                countfre <- 2
                            }
                            else{ 
                                countfre <- 1
                            }
                            if(as.character(zygo)=="homozygous"){
                                homozygo <- as.character(svfamid1)
                            }
                            else{
                                homozygo <- NULL
                            }
                            motherZygosity <- "-"
                            fatherZygosity <- "-"
                        } 
                        else {
                            countfre <- 0
                    
                            motherZygosity <- "-"
                            fatherZygosity <- "-"
                   
                        }
                        if (((identical(type, "duplication") |
                        identical(type, "duplication_split")| 
                        identical(type, "duplication_inverted")|
                        identical(type, "insertion"))) & 
                        !(identical(svfamid1, svfamid))){
                            if(zygo=="homozygous"){
                                countfreunfilt <- 2
                            }
                            else if(zygo=="unknown"){
                                countfreunfilt <- 2
                            }
                            else{ 
                                countfreunfilt <- 1
                            }
                    
                        }
                        else{
                            countfreunfilt <- 0
                        }
                    
                    if(length(homozygo)>0){
                        homozygo <- length(homozygo)
                    }
                    else{
                        homozygo= 0
                    }
                  
                    famno <- as.numeric(unname(hash::values(ha, keys = svfamid)
                        ))
                    data1 <- data.frame(dat1[nn, ], 
                    Internal_Freq_Perc_Filtered = 
                    as.numeric(format(round((countfre/(2*(usamp - famno))) 
                    * 100,digits=3), scientific = FALSE)),
                    Internal_Freq_Perc_Unfiltered = 
                    as.numeric(format(round((countfreunfilt/(2*(usamp - famno)))
                    * 100,digits=3), scientific = FALSE)),
                    Internal_Homozygotes=as.numeric(homozygo), 
                    MotherZygosity = as.character(motherZygosity), 
                    FatherZygosity = as.character(fatherZygosity), 
                    stringsAsFactors = FALSE)
                    
                    datf <- rbind(datf, data1)
                   } 
                   else
                   {
                    # dgv_match=FALSE
                    data1 <- data.frame(dat1[nn, ], 
                    Internal_Freq_Perc_Filtered = 0, 
                    Internal_Freq_Perc_Unfiltered = 0,
                    Internal_Homozygotes= 0, 
                    MotherZygosity = "-", 
                    FatherZygosity = "-", 
                    stringsAsFactors = FALSE)
                  #### print(dim(data1)) print(names(data1)[56:58])
                  #### print(names(datf)[56:58])
                  #### print(identical((names(data1)[56]),(names(datf)[56])))
                  #### print(identical((names(data1)[57]),(names(datf)[57])))
                  #### print(identical((names(data1)[58]),(names(datf)[58])))
                     datf <- rbind(datf, data1)
                  # next;
                }
                
            }
            else if ((length(grep("inversion", variantType2[nn])) >= 1) | 
            (length(grep("translocation", variantType2[nn])) >= 1)) {
                dat2 <- dat[which((dat$RefStartPos <= rf[nn] & 
                dat$RefStartPos >= rf_wb_int[nn] | 
                dat$RefStartPos >= rf[nn] & 
                dat$RefStartPos <= rf_fb_int[nn]) & 
                (dat$RefEndPos >= re[nn] &
                dat$RefEndPos <= re_fb_int[nn] |
                dat$RefEndPos <= re[nn] & 
                dat$RefEndPos >= re_wb_int[nn])), ]
                size1 <- size_bn[nn]
                
                ## Writing if the dgv_match is TRUE or not
                
                if (nrow(dat2) > 1){
                    countfre <- c();countfreunfilt<-c()
                    svfam1 <- dat2$SVIdentifier
                    # svind1<-dat2$SVIndex
                    sv1 <- strsplit(as.character(svfam1), split = "_")
                    #size_internal <- dat2$Size
                    zygo <- as.character(dat2$Zygosity)
                    type <- as.character(dat2$Type)
                    typ1 <- as.character(dat2$Type1)
                    typ2 <- as.character(dat2$Type2)
                    svStrt<-as.numeric(dat2$RefStartPos)
                    svEnd<-as.numeric(dat2$RefEndPos)
                    BSPQI_status_DB <- as.character(
                    dat2$Found_in_self_BSPQI_molecules)
                    BSSSI_status_DB <- as.character(
                    dat2$Found_in_self_BSSSI_molecules)
                    BSPQI_chimeric_score_DB <- as.character(
                        dat2$Fail_BSPQI_assembly_chimeric_score)
                    BSSSI_chimeric_score_DB <- as.character(
                        dat2$Fail_BSSSI_assembly_chimeric_score)
                    conf <- dat2$Confidence
                    motherZygosity <- c()
                    fatherZygosity <- c()
                    homozygo<-c()
                    motherZygosity_exact <- c()
                    fatherZygosity_exact <- c()
                    svfamid3 <- c()
                  
                    svfamid3_unfilt <- c()
                    chrom2<-dat2$RefcontigID2
                    svSAMP<-c()
                    svSAMP_unfiltered<-c()
                    for (ll in seq_len(nrow(dat2))){
                    #perc <- (size1/size_internal[ll])
                    
                        stt <- strsplit(sv1[[ll]][1], 
                        split = "[.]")
                        patID1 <- stt[[1]][2]
                        svfamid1 <- stt[[1]][1]
                        if (((svStrt[ll]==rf[nn]) & 
                        (svEnd[ll]==re[nn])) & 
                        (chrom2[ll]==chromo2[nn]) & 
                        identical(type[ll], variantType2[nn]) &
                        (identical(svfamid1, svfamid))){
                            # Family Extra column zygosity
                            if (patID == 1 & patID1 == 2){
                                motherZygosity_exact <- 
                                c(motherZygosity_exact, as.character(zygo[ll]))
                                fatherZygosity_exact <- 
                                c(fatherZygosity_exact, NULL)
                            } 
                            else if (patID == 1 & patID1 == 3){
                        
                                motherZygosity_exact <- 
                                c(motherZygosity_exact, NULL)
                                fatherZygosity_exact <- 
                                c(fatherZygosity_exact, as.character(zygo[ll]))
                            } 
                            else if (patID == patID1){                        
                                motherZygosity_exact <- 
                                c(motherZygosity_exact, NULL)
                                fatherZygosity_exact <- 
                                c(fatherZygosity_exact, NULL)
                            } 
                            else{
                        
                                motherZygosity_exact <- 
                                c(motherZygosity_exact, NULL)
                                fatherZygosity_exact <- 
                                c(fatherZygosity_exact, NULL)
                            }
                        } 
                        else if ( (((svStrt[ll] <= rf[nn] & 
                        svStrt[ll] >= rf_wb_int_parents[nn] |
                        svStrt[ll] >= rf[nn] & 
                        svStrt[ll] <= rf_fb_int_parents[nn]) &
                        (svEnd[ll] >= re[nn] & 
                        svEnd[ll] <= re_fb_int_parents[nn] |
                        svEnd[ll] <= re[nn] & 
                        svEnd[ll] >= re_wb_int_parents[nn]))) & 
                        identical(chrom2[ll],chromo2[nn]) & 
                        identical(type[ll], variantType2[nn]) & 
                        (identical(svfamid1, svfamid))){
                            # Family Extra column zygosity
                            if (patID == 1 & patID1 == 2) {
                                motherZygosity <- c(motherZygosity, 
                                    as.character(zygo[ll]))
                                fatherZygosity <- c(fatherZygosity, NULL)
                            } 
                            else if (patID == 1 & patID1 == 3){
                        
                                motherZygosity <- c(motherZygosity, NULL)
                                fatherZygosity <- c(fatherZygosity, 
                                as.character(zygo[ll]))
                            } 
                            else if (patID == patID1){
                        
                                motherZygosity <- c(motherZygosity, NULL)
                                fatherZygosity <- c(fatherZygosity, NULL)
                            }
                            else {
                                motherZygosity <- c(motherZygosity, NULL)
                                fatherZygosity <- c(fatherZygosity, NULL)
                            }
                        } 
                        else if (identical(chrom2[ll],chromo2[nn]) & 
                        identical(type[ll], variantType2[nn]) & 
                        !(identical(svfamid1, svfamid)) & 
                        ((BSPQI_status_DB[ll] == "yes" & 
                        BSSSI_status_DB[ll] == "yes") | 
                        (BSPQI_status_DB[ll] == "no" &
                        BSSSI_status_DB[ll] == "yes") |
                        (BSPQI_status_DB[ll] == "yes" &
                        BSSSI_status_DB[ll] == "no") |
                        (BSPQI_status_DB[ll] == "yes" &
                        BSSSI_status_DB[ll] == "-") |
                        (BSPQI_status_DB[ll] == "-" &
                        BSSSI_status_DB[ll] == "yes")) & 
                        ((identical(type[ll], "inversion") & 
                        conf[ll] >= invconf) |
                        (identical(type[ll], "translocation") & 
                        conf[ll] >= transconf)) &
                        ((BSPQI_chimeric_score_DB[ll] == "pass" & 
                        BSSSI_chimeric_score_DB[ll] == "pass") | 
                        (BSPQI_chimeric_score_DB[ll] == "fail" & 
                        BSSSI_chimeric_score_DB[ll] == "pass") |
                        (BSPQI_chimeric_score_DB[ll] == "pass" & 
                        BSSSI_chimeric_score_DB[ll] == "fail") |
                        (BSPQI_chimeric_score_DB[ll] == "pass" & 
                        BSSSI_chimeric_score_DB[ll] == "-") |
                        (BSPQI_chimeric_score_DB[ll] == "-" & 
                        BSSSI_chimeric_score_DB[ll] == "pass")) &  
                        ((typ1[ll]=="inversion" & typ2[ll]=="inversion") | 
                        (typ1[ll]=="inversion" & is.na(typ2[ll])) | 
                        (is.na(typ1[ll]) & typ2[ll]=="inversion") | 
                        (typ1[ll]=="translocation" & 
                        typ2[ll]=="translocation") | 
                        (typ1[ll]=="translocation" & is.na(typ2[ll])) |
                        (is.na(typ1[ll]) & typ2[ll]=="translocation"))){
                            countfre <- c(countfre,1)
                            svfamid3 <- c(svfamid3, as.character(sv1[[ll]][1]))
                            motherZygosity <- c(motherZygosity, "-")
                            fatherZygosity <- c(fatherZygosity, "-")
                            if(as.character(zygo[ll])=="homozygous"){
                                homozygo<-c(homozygo, as.character(sv1[[ll]][1])
                                )
                            }
                            else{
                                homozygo<-c(homozygo, NULL)
                            }
                            svSAMP<-c(svSAMP,as.numeric(ll))
                        } 
                        else{
                            motherZygosity <- c(motherZygosity, "-")
                            fatherZygosity <- c(fatherZygosity, "-")
                      
                        }
                         ##Unfiltered calculation
                    if (identical(chrom2[ll],chromo2[nn]) & 
                    identical(type[ll], variantType2[nn]) &
                    !(identical(svfamid1, svfamid)) ) {
                        countfreunfilt <- c(countfreunfilt,1)
                        svfamid3_unfilt <- c(svfamid3_unfilt, as.character(
                            sv1[[ll]][1]))
                        svSAMP_unfiltered<-c(svSAMP_unfiltered,as.numeric(ll))
                    }
                    else{
                        print ("SV does not meet unfiltered criteria!!!")
                        
                      }
                      
                    
                    }
                  ##Calculating filtered Frequency INDEL
                  ##Filtration
                    if ((length(svfamid3) >= 1) & 
                        length(countfre) >= 1){ 
                        zyg<-as.character(dat2$Zygosity[svSAMP])
                        dat_temp <- data.frame(svfamid3, countfre,zyg)
                        countfre1 <- 0
                        if (nrow(dat_temp) > 0){
                            #### print(paste('dat_temp',dim(dat_temp)))
                            usampid <- as.character(unique(svfamid3))
                      
                            for (u in seq_len(length(usampid))){
                            dat_temp1 <- dat_temp[
                                which(
                                dat_temp$svfamid3 %in% usampid[u]
                                ), ]
                            zygo<-as.character(dat_temp1$zyg)
                            ct <- sum(dat_temp1$countfre)
                            if (ct >= 1 & nrow (dat_temp1)>= 1){
                                if(length(unique(as.character(zygo))) == 1){
                                    if(unique(as.character(zygo)) 
                                    == "homozygous"){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else if(unique(as.character(zygo)) == 
                                        "Unknown"){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else {
                                        countfre1 <- countfre1 + 1
                                    }
                                }
                                else{
                                    g1<-grep("homozygous",as.character(zygo))
                                    g2<-grep("heterozygous",as.character(zygo))   
                                    g3<-grep("unknown",as.character(zygo))  
                                    if(length(g1)>=1 & length(g2)>=1 & 
                                        length(g3)>=1){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else if(length(g1)==0 & 
                                    length(g2)>=1 & 
                                    length(g3)>=1){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else if(length(g1)>=1 & 
                                    length(g2)>=1 & 
                                    length(g3)==0){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else if(length(g1)>=1 & length(g2)==0 &
                                    length(g3)>=1){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else{
                                        countfre1 <- countfre1 + 0
                                    }
                                }
                            }
                            else{
                                countfre1 <- countfre1 + 0
                            }
                        }
                        }else
                        {
                          countfre1 <- 0
                        }
                    }else{
                        countfre1 <- 0
                        }
                    #### Unfiltered
                   if ((length(svfamid3_unfilt) >= 1) & length(countfreunfilt)
                        >= 1)
                  {
                    zyg_unfiltered<-as.character(dat2$Zygosity[svSAMP_unfiltered])
                    dat_temp <- data.frame(svfamid3_unfilt, countfreunfilt,zyg_unfiltered)
                    countfreunfilt1 <- 0
                    if (nrow(dat_temp) > 0)
                    {
                      #### print(paste('dat_temp',dim(dat_temp)))
                      usampid <- as.character(unique(svfamid3_unfilt))
                      
                        for (u in seq_len(length(usampid))) {
                            dat_temp1 <- dat_temp[which(
                                dat_temp$svfamid3_unfilt %in% usampid[u]), ]
                            zygo<-as.character(dat_temp1$zyg_unfiltered)
                            ct <- sum(dat_temp1$countfreunfilt)
                            if (ct >= 1 & nrow (dat_temp1)>= 1)
                            {
                                if(length(unique(as.character(zygo)))==1){
                                    if(unique(as.character(zygo)) ==
                                        "homozygous"){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }
                                    else if(unique(as.character(zygo)) == 
                                        "Unknown"){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }
                                    else {
                                        countfreunfilt1 <- countfreunfilt1 + 1
                                    }
                                }
                                else{
                                    g1<-grep("homozygous",as.character(zygo))
                                    g2<-grep("heterozygous",as.character(zygo))   
                                    g3<-grep("unknown",as.character(zygo))  
                                    if(length(g1)>=1 & length(g2)>=1 & 
                                        length(g3)>=1){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }
                                    else if(length(g1)==0 & length(g2)>=1 & 
                                    length(g3)>=1){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }
                                    else if(length(g1)>=1 & length(g2)>=1 & 
                                        length(g3)==0){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }
                                    else if(length(g1)>=1 & length(g2)==0 & 
                                        length(g3)>=1){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }
                                    else{
                                        countfreunfilt1 <- countfreunfilt1 + 0
                                    }
                                }
                            }else{
                                countfreunfilt1 <- countfreunfilt1 + 0
                            }
                        }
                    }
                    else
                        {
                            countfreunfilt1 <- 0
                        }
                    }
                    else{
                        countfreunfilt1 <- 0
                    }
                 ###Father Zygosity
                  fatherZygosity1 <- gsub("-", NA, fatherZygosity)
                  fatherZygosity1 <- as.character(na.omit(fatherZygosity1))
                  if(length(fatherZygosity_exact)>0){
                    fatherZygosity<-as.character(unique(fatherZygosity_exact))
                    } else {
                        fatherZygosity<-as.character(unique(fatherZygosity1))
                    }
                  
                  if (length(fatherZygosity) > 1)
                  {
                    fatherZygosity <- paste(as.character(fatherZygosity), 
                      collapse = ",")
                  } else if (length(fatherZygosity) == 0)
                  {
                    fatherZygosity <- "-"
                  } else
                  {
                    fatherZygosity <- as.character(fatherZygosity)
                  }
                  ###Mother Zygosity
                  motherZygosity1 <- gsub("-", NA, motherZygosity)
                  motherZygosity1 <- as.character(na.omit(motherZygosity1))
                  if(length(motherZygosity_exact)>0){
                    motherZygosity<-as.character(unique(motherZygosity_exact))
                    } else 
                  {
                  motherZygosity<-as.character(unique(motherZygosity1))
                  }
                 
                  if (length(motherZygosity) > 1)
                  {
                    motherZygosity <- paste(as.character(motherZygosity), 
                      collapse = ",")
                  } else if (length(motherZygosity) == 0)
                  {
                    motherZygosity <- "-"
                  } else
                  {
                    motherZygosity <- as.character(motherZygosity)
                  }
                    if(length(homozygo)>=1){
                        homozygo <- length(unique(homozygo))
                    }else if(length(homozygo)==1){
                        homozygo <- length(homozygo)
                    }else{
                        homozygo <- 0
                    }
                 
                  famno <- as.numeric(unname(hash::values(ha, keys = svfamid)))
                  data1 <- data.frame(dat1[nn, ], 
                  Internal_Freq_Perc_Filtered = as.numeric(
                  format(round((countfre1/(2*(usamp - famno))) * 100,digits=3), 
                  scientific = FALSE)),
                  Internal_Freq_Perc_Unfiltered = as.numeric(
                  format(round((countfreunfilt1/(2*(usamp - famno))) * 
                  100,digits=3), scientific = FALSE)), 
                  Internal_Homozygotes=as.numeric(homozygo), 
                  MotherZygosity = as.character(motherZygosity), 
                  FatherZygosity = as.character(fatherZygosity), 
                  stringsAsFactors = FALSE)
                  #### print(dim(data1))
                  datf <- rbind(datf, data1)
                } 
                else if (nrow(dat2) == 1)
                {
                    ## dgv_match=TRUE Calculating percentage similarity
                    countfre<-0;countfreunfilt<-0
                    chrom2<-dat2$RefcontigID2
                    svv <- strsplit(as.character(svfam1), split = "_")
                    zygo <- as.character(dat2$Zygosity)
                    stt <- strsplit(svv[[1]][1], split = "[.]")
                    type <- as.character(dat2$Type)
                    typ1 <- as.character(dat2$Type1)
                    typ2 <- as.character(dat2$Type2)
                    patID1 <- stt[[1]][2]
                    svfamid1 <- stt[[1]][1]
                    conf <- dat2$Confidence
                    homozygo<-c()
                    BSPQI_status_DB <- as.character(
                    dat2$Found_in_self_BSPQI_molecules)
                    BSSSI_status_DB <- as.character(
                    dat2$Found_in_self_BSSSI_molecules)
                    BSPQI_chimeric_score_DB <- as.character(
                    dat2$Fail_BSPQI_assembly_chimeric_score)
                    BSSSI_chimeric_score_DB <- as.character(
                    dat2$Fail_BSSSI_assembly_chimeric_score)
                    if ((((dat2$RefStartPos <= rf[nn] & 
                        dat2$RefStartPos >= rf_wb_int_parents[nn] | 
                        dat2$RefStartPos >= rf[nn] & 
                        dat2$RefStartPos <= rf_fb_int_parents[nn]) & 
                        (dat2$RefEndPos >=  re[nn] & 
                        dat2$RefEndPos <= re_fb_int_parents[nn] | 
                        dat2$RefEndPos <= re[nn] & 
                        dat2$RefEndPos >= re_wb_int_parents[nn]))) & 
                        identical(chrom2,chromo2[nn])& 
                        identical(type, variantType2[nn]) & 
                        (identical(svfamid1, svfamid))){
                    # Family Extra column zygosity
                        if (patID == 1 & patID1 == 2){
                            countfre <- countfre + 0
                            motherZygosity <- as.character(zygo)
                            fatherZygosity <- "-"
                        } 
                        else if (patID == 1 & patID1 == 3){
                            countfre <- countfre + 0
                            motherZygosity <- "-"
                            fatherZygosity <- as.character(zygo)
                        } 
                        else if (patID == patID1){
                            countfre <- countfre + 0
                            motherZygosity <- "-"
                            fatherZygosity <- "-"
                        } 
                        else{
                      
                            motherZygosity <- "-"
                            fatherZygosity <- "-"
                        }
                    }  
                    else if (identical(chrom2,chromo2[nn]) &
                    identical(type,variantType2[nn]) & 
                    !(identical(svfamid1, svfamid)) & 
                    ((BSPQI_status_DB == "yes" & 
                    BSSSI_status_DB == "yes") | 
                    (BSPQI_status_DB == "no" &
                    BSSSI_status_DB == "yes") | 
                    (BSPQI_status_DB == "yes" & 
                    BSSSI_status_DB == "no") | 
                    (BSPQI_status_DB == "yes" &
                    BSSSI_status_DB == "-") | 
                    (BSPQI_status_DB == "-" & 
                    BSSSI_status_DB == "yes")) & 
                    ((identical(type, "inversion") & conf > invconf) | 
                    (identical(type, "translocation") & conf > transconf)) & 
                    ((BSPQI_chimeric_score_DB == "pass" & 
                    BSSSI_chimeric_score_DB == "pass") | 
                    (BSPQI_chimeric_score_DB == "fail" & 
                    BSSSI_chimeric_score_DB == "pass") | 
                    (BSPQI_chimeric_score_DB == "pass" & 
                    BSSSI_chimeric_score_DB == "fail") |
                    (BSPQI_chimeric_score_DB == "pass" & 
                    BSSSI_chimeric_score_DB == "-") | 
                    (BSPQI_chimeric_score_DB == "-" & 
                    BSSSI_chimeric_score_DB == "pass")) & 
                    ((typ1=="inversion" & typ2=="inversion") | 
                    (typ1=="inversion" & is.na(typ2)) | 
                    (is.na(typ1) & typ2=="inversion") | 
                    (typ1=="translocation" & typ2=="translocation") |
                    (typ1=="translocation" & is.na(typ2)) | 
                    (is.na(typ1) & typ2=="translocation"))){
                        if(zygo=="homozygous"){
                            countfre <- 2
                        }
                        else if(zygo=="unknown"){
                            countfre <- 2
                        }
                        else{ 
                            countfre <- 1
                        }
                        if(as.character(zygo)=="homozygous"){
                            homozygo<- as.character(svfamid1)
                        }
                        else{
                            homozygo<-  NULL
                        }
                        motherZygosity <- "-"
                        fatherZygosity <- "-"
                    } 
                    else {
                        countfre <- countfre + 0
                        motherZygosity <- "-"
                        fatherZygosity <- "-"
                    }
                    if (identical(chrom2,chromo2[nn]) & 
                        identical(type[ll], variantType2[nn]) & 
                        !(identical(svfamid1, svfamid)) ){
                      
                        if(zygo=="homozygous"){
                            countfreunfilt <- 2
                        }
                        else if(zygo=="unknown"){
                            countfreunfilt <- 2
                        }
                        else{ 
                            countfreunfilt <- 1
                        }
                        }
                        else{
                            countfreunfilt <-  0
                        }
                     
                     
                    ##Checking for whether there are any homozygous variants or 
                    ##not
                    if(length(homozygo)>0){
                        homozygo <- length(homozygo)
                    }
                    else{
                        homozygo= 0
                    }
                  ##### print(paste('INVTRANS:',countfre,sep=''))
                    famno <- as.numeric(unname(hash::values(ha, 
                        keys = svfamid)))
                    data1 <- data.frame(dat1[nn, ], 
                    Internal_Freq_Perc_Filtered = 
                    as.numeric(format(round((countfre/(2*(usamp - famno))) 
                    * 100,digits=3), scientific = FALSE)),
                    Internal_Freq_Perc_Unfiltered = 
                    as.numeric(format(round((countfreunfilt/(2*(usamp - famno))) 
                    * 100,digits=3), scientific = FALSE)), 
                    Internal_Homozygotes=as.numeric(homozygo), 
                    MotherZygosity = as.character(motherZygosity), 
                    FatherZygosity = as.character(fatherZygosity),  
                    stringsAsFactors = FALSE)
                  
                  #### print(dim(data1))
                    datf <- rbind(datf, data1)
                } 
                else{
                 
                    data1 <- data.frame(dat1[nn, ], 
                    Internal_Freq_Perc_Filtered= 0, 
                    Internal_Freq_Perc_Unfiltered =0,
                    Internal_Homozygotes= 0, MotherZygosity = "-", 
                    FatherZygosity = "-",  stringsAsFactors = FALSE)
                  #### print(dim(data1))
                    datf <- rbind(datf, data1)
                   # next;
                }
                
            }
            else{
                #### print(paste('QueryData:',dim(dat1[nn,]),sep=''))
                
                data1 <- data.frame(dat1[nn, ], Internal_Freq_Perc_Filtered= 0,
                Internal_Freq_Perc_Unfiltered = 0, 
                Internal_Homozygotes= 0, MotherZygosity = "-", 
                FatherZygosity = "-", stringsAsFactors = FALSE)
                #### print(dim(data1))
                datf <- rbind(datf, data1)
            }
            
        }
        dataFinal <- rbind(dataFinal, datf)
    }
    if(returnMethod=="Text"){
        st1 <- strsplit(smap, ".txt")
        fname <- st1[[1]][1]
        row.names(dataFinal) <- c()
        filenam<-paste(fname,"_Int.txt",sep="")
        write.table(dataFinal, file.path(smappath, fname), sep = "\t", 
            row.names = FALSE)
    }
    else if (returnMethod=="dataFrame"){
        return(dataFinal)
    }
    else{
        stop("returnMethod Incorrect")
    }
}



