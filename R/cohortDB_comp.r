#' Merges Solo SV files to one common SV file
#'
#' @param path  character. Path to the solo files.
#' @param pattern  character. file name pattern for solo files.
#' @param outpath  character. file path for the output file.
#' @param fname  character. file name for the output file.
#' @param dbOutput  character. Output option database or text.
#' @return Text file containing all the solo SMAP files.
#' @examples
#' path <- system.file("extdata", "Bionano_config/", package = "nanotatoR")
#' pattern <- "_hg19.txt"
#' mergedSmap <- makeInternalBNDatabase(path = path, 
#'    pattern = pattern, dbOutput = "dataframe")
#' mergedSmap[1,]
#' @importFrom stats na.omit 
#' @export
makeInternalBNDatabase <- function(path, pattern, outpath,fname,
                            dbOutput=c("dataframe","text")) {
    #setwd(path)
    allfiles <- list.files(path, pattern)
    nam <- c()
    datfinal <- data.frame()
    for (files in allfiles) {
        print(files)
        ##### print(dim(dat4))
        r1 <- read.table(file.path(path, files), header = TRUE)
        print(dim(r1))
        samp<-as.character(r1$sample)
        typ<-as.character(r1$type)
        zygo<-as.character(r1$zygosity)
        if(typ[1]=="trans"){
        chro1<-as.numeric(r1$chrA)
        chro2<-as.numeric(r1$chrB)
        bpstrt<-as.numeric(r1$bkptA)
        bpend<-as.numeric(r1$bkptB)
        siz<-as.numeric(rep(-1,length(bpend)))
        }
        else{
        chro1<-as.numeric(r1$chr)
        chro2<-as.numeric(r1$chr)
        bpstrt<-as.numeric(r1$start)
        bpend<-as.numeric(r1$end)
        siz<-as.numeric(r1$size)
        }
        score<-as.numeric(r1$score)
        smapID<-as.numeric(r1$smapId)
        typ<-gsub("^ins$","insertion",typ)
        typ<-gsub("^del$","deletion",typ)
        typ<-gsub("^inv$","inversion",typ)
        typ<-gsub("^trans$","translocation",typ)
        dattemp<-data.frame(SmapID=smapID,Sample=samp,Type=typ,Chromosome1=chro1,
                Chromosome2=chro2,BreakPntStrt=bpstrt,BreakPntEnd=bpend,
                Size=siz,Zygosity=zygo,Score=score)
        datfinal<-rbind(datfinal,dattemp)
        

    }
    if(dbOutput=="text"){
    filename <- paste(fname, "_merged.txt", sep = "")
    write.table(datfinal, file.path(outpath, filename), sep = "\t")
    }else if(dbOutput=="dataframe"){
    return(datfinal)
    }else{stop("method not found")}
}

#' Calculates the internal frequencies of BNDB cohorts
#'
#' @param internalBNDB  character. Path to the merged SV files.
#' @param smappath  character. path to the query smap file.
#' @param smap  character. File name for the smap 
#' @param smapdata  dataframe. smapdata in the form of dataframe. 
#' @param buildBNInternalDB  boolean. Checking whether the merged BNDB 
#' file database exist.
#' @param EnzymeType Character. Type of enzyme. Options SVmerge and SE.
#' @param BNDBpath  character. Path to the BNDB file database. 
#' @param dbOutput  character. database output type.
#' Options dataframe or text.
#' @param BNDBpattern  character. pattern of the file names to merge.
#' @param outpath  character. Path to merged SV solo datasets.
#' @param win_indel  Numeric. Insertion and deletion error window.
#' @param win_inv_trans  Numeric. Inversion and translocation error window.
#' @param perc_similarity  Numeric . ThresholdPercentage similarity 
#' of the query SV and reference SV.
#' @param indelconf  Numeric. Threshold for insertion and deletion Score.
#' @param invconf  Numeric. Threshold for inversion Score.
#' @param transconf  Numeric. Threshold for translocation Score.
#' @param limsize  Numeric. SV size limit.
#' @param returnMethod character. Choice between Text and DataFrame.
#' @param input_fmt_SV character. Choice between Text and DataFrame.
#' @return Text file or data frames containing internalFrequency data.
#' @examples
#' path <- system.file("extdata", "Bionano_config/", package = "nanotatoR")
#' pattern <- "*_hg19_*"
#' smapName="GM24385_Ason_DLE1_VAP_trio5.smap"
#' smap = system.file("extdata", smapName, package="nanotatoR")
#' BNDBfrequency(smap = smap, 
#'  buildBNInternalDB=TRUE, 
#'    input_fmt_SV = "Text",
#'    dbOutput="dataframe",
#'    BNDBpath = path, 
#'    BNDBpattern = pattern, 
#'    outpath, 
#'    win_indel = 10000,
#'    win_inv_trans = 50000, 
#'  perc_similarity = 0.5, 
#'    indelconf = 0.5, 
#'    invconf = 0.01, 
#'    limsize = 1000,
#'  transconf = 0.1,
#'    returnMethod=c("dataFrame"), 
#'    EnzymeType = c("SE"))
#' @importFrom stats na.omit 
#' @export
BNDBfrequency  <- function(internalBNDB, smappath , smap , 
    buildBNInternalDB = FALSE, smapdata, 
    input_fmt_SV = c("Text","dataFrame"),
    dbOutput = c("dataframe","text"),
    BNDBpath, 
    BNDBpattern, 
    outpath, 
    win_indel = 10000,
    win_inv_trans = 50000, 
    perc_similarity = 0.5, 
    indelconf = 0.5, 
    invconf = 0.01, 
    limsize = 1000,
    transconf = 0.1,
    returnMethod=c("Text","dataFrame"), 
    EnzymeType = c("SVmerge", "SE"))
    {
    #library(hash)
    print("###Cohort Frequency Calculation###")
    if(buildBNInternalDB==TRUE){
        print("###Building Database###")
        r<-makeInternalBNDatabase(path = BNDBpath, 
        pattern = BNDBpattern,dbOutput=dbOutput)
    } else{
        print("###Reading Database###")
        r <- read.table(internalBNDB, sep = "\t", header = TRUE)
    }
    usampt<-as.character(unique(r$Sample))
    strst<-strsplit(usampt,split="_")
    usampp<-c()
    for(o in 1:length(strst)){
    usampp<-c(usampp,as.character(strst[[o]][1]))
    }
    
    usamp<- length((usampp))
    print("###Calculating BN Frequency###")
    if(input_fmt_SV=="dataFrame"){
        smapdata = smapdata
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
    sampID <- str_squish(as.character(unique(r1$SampleID)))
    ## Checking Sex male/female and assigning chromosome number accordingly
    chro <-unique(r1$RefcontigID1)
    dataFinal <- c()
    for (ii in seq_along(chro))
    {
        print(paste('Chrom:',chro[ii],sep=''))
        dat <- r[which(r$Chromosome1 == chro[ii]), ]
        
        # variantType1<-dat$variantsubtype Changing the variant terms in DGV to
        # match svmap
        variantType1 <- as.character(dat$Type)
        #BSPQI_status_DB <- as.character(dat$Found_in_self_BSPQI_molecules)
        #BSSSI_status_DB <- as.character(dat$Found_in_self_BSSSI_molecules)
        ## Extracting data from SVmap
        dat1 <- r1[which(r1$RefcontigID1 == chro[ii]), ]
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
                
        # conf<-dat$Confidence countfre<-0 percn<-c()
        datf <- c()
        for (nn in 1:length(rf))
        {
            #print(paste("nn:",nn)) 
            
            if ((variantType2[nn] == "deletion" 
                | variantType2[nn] == "insertion"))
            {
                
                dat2 <-dat[which((dat$BreakPntStrt >= rf_wb_ind[nn] 
                    & dat$BreakPntEnd <= re_wf_ind[nn])),]
                size1 <- size_bn[nn]
                ## Calculating Internal Frequency
                if (nrow(dat2) > 1)
                {
                    countfre <- c();countfreunfilt<-c()
                    svStrt<-as.numeric(dat2$BreakPntStrt)
                    svEnd<-as.numeric(dat2$BreakPntEnd)
                    size_internal <- as.numeric(dat2$Size)
                  
                    type <- as.character(dat2$Type)
                    conf <- as.numeric(dat2$Score)
                    zygo <- as.character(dat2$Zygosity)
                    samp <- as.character(dat2$Sample)
                  
                    sampid<-c();sampid_unfiltered<-c();
                    zyg<-c();zyg_unfiltered<-c()
                    homozygo<-c(); unmapped<-c()
                for (ll in 1:nrow(dat2))
                  
                {
                    perc_ref_query <- as.numeric(dat2$Size)/size1
                    perc_query_ref <- size1/as.numeric(dat2$Size)
                    
                    #perc <- (size1/size_internal[ll])
                    #### print(perc) print(type[ll])
                    #print(perc)
                    
                    #### print(svfamid1)
                    if ((perc_ref_query[ll] >= perc_similarity 
                        & perc_query_ref[ll] >= perc_similarity) 
                        & identical(type[ll], variantType2[nn]) 
                        & (conf[ll] >= indelconf) 
                        & size_internal[ll] >=limsize)
                      {
                      
                        countfre <- c(countfre,1)
                        sampid<-c(sampid,samp[ll])                       
                        zyg<-c(zyg,zygo[ll])
                        if(as.character(zygo[ll])=="homozygous"){
                            homozygo<-c(homozygo, as.character(samp[ll]))
                       }
                        else{
                            homozygo<-c(homozygo, NULL)
                       }
                       } else {
                            unmapped<-c(unmapped, ll)
                        }
                        ##Unfiltered
                        if ((perc_ref_query[ll] >= perc_similarity 
                            & perc_query_ref[ll] >= perc_similarity) 
                            & identical(type[ll], variantType2[nn]))
                        {
                        
                       
                            countfreunfilt <- c(countfreunfilt, 1)
                            sampid_unfiltered<-c(sampid_unfiltered,as.character(samp[ll]))
                            zyg_unfiltered<-c(zyg_unfiltered,as.character(zygo[ll]))
                    
                       } 
                  }
                  ##Calculating filtered Frequency
                  ##Filtration
                   if ((length(sampid) >= 1) & length(countfre) >= 1)
                  { 
                    
                    dat_temp <- data.frame(sampid, countfre,zyg)
                    countfre1 <- 0
                    if (nrow(dat_temp) > 0)
                    {
                      #### print(paste('dat_temp',dim(dat_temp)))
                      usampid <- as.character(unique(sampid))
                      
                        for (u in 1:length(usampid))
                        {
                            dat_temp1 <- dat_temp[which(dat_temp$sampid %in% usampid[u]), ]
                            zygo<-as.character(dat_temp1$zyg)
                            ct <- sum(dat_temp1$countfre)
                            if (ct >= 1 & nrow (dat_temp1)>= 1)
                            {
                                if(length(unique(as.character(zygo)))==1){
                                    if(
                                        unique(as.character(zygo)) == "homozygous"
                                    ){
                                        countfre1 <- countfre1 + 2
                                    }else if(
                                        unique(as.character(zygo)) == "Unknown"
                                        ){
                                        countfre1 <- countfre1 + 2
                                    }else {
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
                                    }else if(length(g1)==0 
                                        & length(g2)>=1 & length(g3)>=1){
                                            countfre1 <- countfre1 + 2
                                    }else if(length(g1)>=1 & length(g2)>=1 
                                        & length(g3)==0){
                                            countfre1 <- countfre1 + 2
                                    }
                                    else if(length(g1)>=1 
                                        & length(g2)==0 
                                        & length(g3)>=1){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else{
                                        countfre1 <- countfre1 + 0
                                    }
                                }
                            }else{
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
                    if ((length(sampid_unfiltered) >= 1) 
                        & length(countfreunfilt) >= 1)
                    {
                    
                        dat_temp <- data.frame(sampid_unfiltered, 
                            countfreunfilt,zyg_unfiltered)
                        countfreunfilt1 <- 0
                        if (nrow(dat_temp) > 0){
                            #### print(paste('dat_temp',dim(dat_temp)))
                            usampid <- as.character(unique(sampid_unfiltered))
                      
                                for (u in 1:length(usampid)){
                                    dat_temp1 <- dat_temp[
                                        which(
                                        dat_temp$sampid_unfiltered %in% usampid[u]
                                        ), ]
                                zygo<-as.character(dat_temp1$zyg_unfiltered)
                                ct <- sum(dat_temp1$countfreunfilt)
                                    if (ct >= 1 & nrow (dat_temp1)>= 1){
                                        if(length(unique(
                                            as.character(zygo)))==1){
                                            if(unique(
                                                as.character(zygo))=="homozygous"
                                            ){
                                                countfreunfilt1 <- countfreunfilt1 + 2
                                            }else if(unique(
                                                as.character(zygo))=="Unknown"){
                                                countfreunfilt1 <- countfreunfilt1 + 2
                                            }else {
                                                countfreunfilt1 <- countfreunfilt1 + 1
                                            }
                                        }
                                        else{
                                            g1<-grep(
                                                "homozygous",as.character(zygo)
                                            )
                                            g2<-grep(
                                                "heterozygous",as.character(zygo)
                                            )   
                                            g3<-grep(
                                                "unknown",as.character(zygo))    
                                            if(length(g1)>=1 & length(g2)>=1 
                                                & length(g3)>=1){
                                                countfreunfilt1 <- countfreunfilt1 + 2
                                            }else if(length(g1)==0 & length(g2)>=1 
                                                & length(g3)>=1){
                                                countfreunfilt1 <- countfreunfilt1 + 2
                                            }else if(length(g1)>=1 & length(g2)>=1 
                                                & length(g3)==0){
                                                countfreunfilt1 <- countfreunfilt1 + 2
                                            }else if(length(g1)>=1 & length(g2)==0 
                                                & length(g3)>=1){
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
                            }else{
                                countfreunfilt1 <- 0
                            }
                    }else{
                        countfreunfilt1 <- 0
                    }
                    ##Extracting the number of Homozygotes
                    if(length(homozygo)>=1){
                        homozygo <- length(unique(homozygo))
                    }else if(length(homozygo)==1){
                        homozygo <- length(homozygo)
                    }else{
                        homozygo <- 0
                    }
                  ##### print(names(datf)[56:58]) print(paste('INSDEL:',countfre1,sep=''))
                  
                    data1 <- data.frame(
                            dat1[nn, ], 
                            BNG_Freq_Perc_Filtered = as.numeric(format(round((countfre1/(2*(usamp))) * 100,digits=3), scientific = FALSE)),
                            BNG_Freq_Perc_UnFiltered = as.numeric(format(round((countfreunfilt1/(2*(usamp))) * 100,digits=3), scientific = FALSE)), 
                            BNG_Homozygotes=as.numeric(homozygo), stringsAsFactors = FALSE)
                  #### print(dim(data1)) print(dim(data1))
                  #### print(identical(names(data1),names(datf))) print(names(data1)[56:58])
                   datf <- rbind(datf, data1)
                } else if (nrow(dat2) == 1)
                {
                    # dgv_match=TRUE Calculating percentage similarity
                    countfre <- 0;countfreunfilt<-0
                    conf <- dat2$Score
                    size_internal <- dat2$Size
                    perc_ref_query <- as.numeric(dat2$Size)/size1
                    perc_query_ref <- size1/as.numeric(dat2$Size)
                    zygo <- as.character(dat2$Zygosity)
                    type <- as.character(dat2$Type)
                    samp <- as.character(dat2$Sample)
                    homozygo<-c()
                    if ((perc_ref_query >= perc_similarity 
                        & perc_query_ref >= perc_similarity) 
                        & identical(type, variantType2[nn])
                        & conf >= indelconf 
                        & size_internal>=limsize)
                        {
                            if(zygo =="homozygous"){
                                countfre <- countfre + 2
                                homozygo <- as.character(samp)
                            }else if(zygo =="unknown"){
                                countfre <- countfre + 2
                            }else{ 
                                countfre <- countfre+1
                            }
                        } else{
                           countfre <- 0
                        }
                    if ((perc_ref_query >= perc_similarity 
                        & perc_query_ref >= perc_similarity) 
                        & identical(type, variantType2[nn]))
                        {
                        if(zygo =="homozygous"){
                            countfreunfilt <- countfreunfilt + 2
                        }else if(zygo =="unknown"){
                            countfreunfilt <- countfreunfilt + 2
                        }else{ 
                            countfreunfilt <- countfreunfilt + 1
                        }
                    }
                    ##### print(paste('INSDEL:',countfre,sep='')) print(names(datf)[56:58])
                  ##Calculating the length
                    if(length(homozygo)>0){
                    homozygo <- length(homozygo)
                    }else{
                        homozygo=0
                    }
                  #print(paste("countfre:",countfre,"HomoZygo:",homozygo,sep=""))
                    data1 <- data.frame(dat1[nn, ], 
                        BNG_Freq_Perc_Filtered = as.numeric(
                            format(round((countfre/(2*(usamp))) * 100,digits=3
                            ), 
                            scientific = FALSE)),
                            BNG_Freq_Perc_UnFiltered = as.numeric(
                                format(round((countfreunfilt/(2*(usamp))) * 100,digits=3),
                                scientific = FALSE)
                                ), BNG_Homozygotes=as.numeric(homozygo), stringsAsFactors = FALSE)
                  #### print(dim(data1)) print(names(data1)[56:58])
                  #### print(identical(names(data1),names(datf))) print(names(data1))
                  datf <- rbind(datf, data1)
                } else
                {
                  # dgv_match=FALSE
                  data1 <- data.frame(dat1[nn, ], BNG_Freq_Perc_Filtered = 0, BNG_Freq_Perc_UnFiltered = 0, BNG_Homozygotes = 0,stringsAsFactors = FALSE)
                  #### print(dim(data1)) print(names(data1)[56:58])
                  #### print(names(datf)[56:58])
                  #### print(identical((names(data1)[56]),(names(datf)[56])))
                  #### print(identical((names(data1)[57]),(names(datf)[57])))
                  #### print(identical((names(data1)[58]),(names(datf)[58])))
                  datf <- rbind(datf, data1)
                  # next;
                }
                
            }else if ((variantType2[nn] == "duplication" 
                | variantType2[nn] == "duplication_split"
                |variantType2[nn] == "duplication_inverted"
                | variantType2[nn] == "insertion"))
            {
                dat2 <-dat[which((dat$BreakPntStrt >= rf_wb_int[nn] 
                    & dat$BreakPntEnd <= re_fb_int[nn])),]
                size1 <- size_bn[nn]
                ## Calculating Internal Frequency
                if (nrow(dat2) > 1)
                {
                    countfre <- c();countfreunfilt<- c()
                    #print (paste(nn, ":", ll))
                    # svind1<-dat2$SVIndex
                    
                    size_internal <- dat2$Size
                    type <- as.character(dat2$Type)
                    conf <- dat2$Score
                    zygo <- as.character(dat2$Zygosity)
                    samp <- as.character(dat2$Sample)
                    sampid<-c();sampid_unfiltered<-c();
                    zyg<-c();zyg_unfiltered<-c()
                    unmapped<-c()
                    homozygo<-c()
                    for (ll in 1:nrow(dat2))
                    {
                        #print (ll)
                        #print(type[ll])
                        perc_ref_query <- as.numeric(dat2$Size)/size1
                        perc_query_ref <- size1/as.numeric(dat2$Size)
                        #### print(perc) 
                        #### print(svfamid1)
                        if ((identical(type[ll], "duplication") |
                            identical(type[ll], "duplication_split")| 
                            identical(type[ll], "duplication_inverted")
                            & (perc_ref_query[ll] >= perc_similarity 
                            & perc_query_ref[ll] >= perc_similarity)))
                        {
                            countfre <- c(countfre,1)
                            sampid<-c(sampid,samp[ll])                       
                            zyg<-c(zyg,zygo[ll])
                            if(as.character(zygo[ll])=="homozygous"){
                                homozygo<-c(homozygo, as.character(samp[ll]))
                           }
                            else{
                                homozygo<-c(homozygo, NULL)
                           }
                        } else if ((identical(type[ll], "insertion")) 
                            & (perc_ref_query[ll] >= perc_similarity 
                            & perc_query_ref[ll] >= perc_similarity) 
                            & (size_internal[ll]>=limsize) 
                            & (conf[ll] >= indelconf))
                        {
                            countfre <- c(countfre,1)
                            sampid<-c(sampid,samp[ll])                       
                            zyg<-c(zyg,zygo[ll])
                            if(as.character(zygo[ll])=="homozygous"){
                                homozygo<-c(homozygo, as.character(samp[ll]))
                            }
                            else{
                                homozygo<-c(homozygo, NULL)
                            }
                        }                    
                        else {
                               unmapped<-c(unmapped, ll)
                            }
                            ##Unfiltered
                            if ((identical(type[ll], "duplication") |
                            identical(type[ll], "duplication_split")| 
                            identical(type[ll], "duplication_inverted")|
                            identical(type[ll], "insertion")
                            & (perc_ref_query[ll] >= perc_similarity 
                            & perc_query_ref[ll] >= perc_similarity)))
                            {   
                           
                            countfreunfilt <- c(countfreunfilt, 1)
                            sampid_unfiltered<-c(
                                sampid_unfiltered,as.character(samp[ll])
                                )
                            zyg_unfiltered<-c(zyg_unfiltered,as.character(zygo[ll])
                            )
                        
                            } 
                        
                    }
                  ##Calculating filtered Frequency
                 ##Filtration
                  ##Filtration
                   if ((length(sampid) >= 1) & length(countfre) >= 1)
                  { 
                    
                    dat_temp <- data.frame(sampid, countfre,zyg)
                    countfre1 <- 0
                    if (nrow(dat_temp) > 0)
                    {
                      #### print(paste('dat_temp',dim(dat_temp)))
                      usampid <- as.character(unique(sampid))
                      
                        for (u in 1:length(usampid))
                        {
                            dat_temp1 <- dat_temp[which(dat_temp$sampid %in% usampid[u]), ]
                            zygo<-as.character(dat_temp1$zyg)
                            ct <- sum(dat_temp1$countfre)
                            if (ct >= 1 & nrow (dat_temp1)>= 1)
                            {
                                if(length(unique(as.character(zygo)))==1){
                                    if(unique(as.character(zygo))=="homozygous"){
                                        countfre1 <- countfre1 + 2
                                    }else if(unique(as.character(zygo))=="Unknown"){
                                        countfre1 <- countfre1 + 2
                                    }else {
                                        countfre1 <- countfre1 + 1
                                    }
                                }
                                else{
                                    g1<-grep("homozygous",as.character(zygo))
                                    g2<-grep("heterozygous",as.character(zygo))   
                                    g3<-grep("unknown",as.character(zygo))    
                                    if(length(g1)>=1 & length(g2)>=1 & length(g3)>=1){
                                        countfre1 <- countfre1 + 2
                                    }else if(length(g1)==0 & length(g2)>=1 & length(g3)>=1){
                                        countfre1 <- countfre1 + 2
                                    }else if(length(g1)>=1 & length(g2)>=1 & length(g3)==0){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else if(length(g1)>=1 & length(g2)==0 & length(g3)>=1){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else{
                                        countfre1 <- countfre1 + 0
                                    }
                                }
                            }else{
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
                   if ((length(sampid_unfiltered) >= 1) & length(countfreunfilt) >= 1)
                  {
                    
                    dat_temp <- data.frame(sampid_unfiltered, countfreunfilt,zyg_unfiltered)
                    countfreunfilt1 <- 0
                    if (nrow(dat_temp) > 0)
                    {
                      #### print(paste('dat_temp',dim(dat_temp)))
                      usampid <- as.character(unique(sampid_unfiltered))
                      
                        for (u in 1:length(usampid))
                        {
                            dat_temp1 <- dat_temp[which(
                                dat_temp$sampid_unfiltered %in% usampid[u]), ]
                            zygo<-as.character(dat_temp1$zyg_unfiltered)
                            ct <- sum(dat_temp1$countfreunfilt)
                            if (ct >= 1 & nrow (dat_temp1)>= 1)
                            {
                                if(length(unique(as.character(zygo)))==1){
                                    if(unique(as.character(zygo))=="homozygous"){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }else if(unique(as.character(zygo))=="Unknown"){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }else {
                                        countfreunfilt1 <- countfreunfilt1 + 1
                                    }
                                }
                                else{
                                    g1<-grep("homozygous",as.character(zygo))
                                    g2<-grep("heterozygous",as.character(zygo))   
                                    g3<-grep("unknown",as.character(zygo))    
                                    if(length(g1)>=1 
                                        & length(g2)>=1 & length(g3)>=1){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }else if(length(g1)==0 
                                        & length(g2)>=1 & length(g3)>=1){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }else if(length(g1)>=1 
                                        & length(g2)>=1 & length(g3)==0){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }
                                    else if(length(g1)>=1 
                                        & length(g2)==0 & length(g3)>=1){
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
                    }else
                        {
                          countfreunfilt1 <- 0
                        }
                    }else{
                        countfreunfilt1 <- 0
                        }
                    ##Extracting the number of Homozygotes
                    if(length(homozygo)>=1){
                        homozygo <- length(unique(homozygo))
                    }else if(length(homozygo)==1){
                        homozygo <- length(homozygo)
                    }else{
                        homozygo <- 0
                    }
                  ##### print(names(datf)[56:58]) print(paste('INSDEL:',countfre1,sep=''))
                  
                    data1 <- data.frame(dat1[nn, ], 
                        BNG_Freq_Perc_Filtered = as.numeric(
                            format(round((countfre1/(2*(usamp))) * 100,digits=3), scientific = FALSE)),
                            BNG_Freq_Perc_UnFiltered = as.numeric(format(round((countfreunfilt1/(2*(usamp))) * 100,digits=3), scientific = FALSE)),
                            BNG_Homozygotes=as.numeric(homozygo), stringsAsFactors = FALSE)
                  #### print(dim(data1)) print(dim(data1))
                  #### print(identical(names(data1),names(datf))) print(names(data1)[56:58])
                  datf <- rbind(datf, data1)
                } else if (nrow(dat2) == 1)
                {
                  # dgv_match=TRUE Calculating percentage similarity
                  countfre <- 0;countfreunfilt<-0
                    conf <- dat2$Score
                    size_internal <- dat2$Size
                    samp <- as.character(dat2$Sample)
                    type <- as.character(dat2$Type)
                    zygo <- as.character(dat2$Zygosity)
                    perc_ref_query <- as.numeric(dat2$Size)/size1
                    perc_query_ref <- size1/as.numeric(dat2$Size)
                    homozygo<-c()
                  if ((identical(type, "duplication") |
                        identical(type, "duplication_split")| 
                        identical(type, "duplication_inverted")
                        & (perc_ref_query >= perc_similarity 
                        & perc_query_ref >= perc_similarity)))
                    {
                        if(zygo =="homozygous"){
                        countfre <- countfre+2
                        homozygo<-as.character(samp)
                       }else if(zygo =="unknown"){
                        countfre <- countfre+2
                       }else{ 
                        countfre <- countfre+1
                       }
                    } else if ((identical(type, "insertion")) & (size_internal >=limsize)
                                & (conf >= indelconf) 
                                & (perc_ref_query >= perc_similarity 
                                & perc_query_ref >= perc_similarity))
                    {
                        if(zygo =="homozygous"){
                        countfre <- countfre+2
                        homozygo<-as.character(samp)
                       }else if(zygo =="unknown"){
                        countfre <- countfre+2
                       }else{ 
                        countfre <- countfre+1
                       }
                    }                    
                    else {
                           countfre <- 0
                           
                        }
                        ##Unfiltered
                    if ((identical(type, "duplication") |
                        identical(type, "duplication_split")| 
                        identical(type, "duplication_inverted")|
                        identical(type, "insertion"))
                        & (perc_ref_query >= perc_similarity 
                        & perc_query_ref >= perc_similarity))
                        {   
                       
                        if(zygo =="homozygous"){
                            countfreunfilt <- countfreunfilt + 2
                        }else if(zygo =="unknown"){
                            countfreunfilt <- countfreunfilt + 2
                        }else{ 
                            countfreunfilt <- countfreunfilt + 1
                        }
                        } else{
                            countfreunfilt <- 0
                            } 
                    
                    'if ((identical(type, "duplication") |
                        identical(type, "duplication_split")| 
                        identical(type, "duplication_inverted")|
                        identical(type, "insertion")))
                    { ##Filtered
                      if(zygo =="homozygous"){
                        countfre <- countfre+2
                        homozygo<-as.character(samp)
                       }else if(zygo =="unknown"){
                        countfre <- countfre+2
                       }else{ 
                        countfre <- countfre+1
                       }
                      ##Un-Filtered
                      if(zygo =="homozygous"){
                        countfreunfilt <- countfreunfilt + 2
                       }else if(zygo =="unknown"){
                        countfreunfilt <- countfreunfilt + 2
                       }else{ 
                        countfreunfilt <- countfreunfilt + 1
                        }
                    
                    } else
                    {
                    
                      
                    }'
                    ##Calculating the length
                    if(length(homozygo)>0){
                        homozygo <- length(homozygo)
                    }else{
                        homozygo=0
                    }
                  ##### print(paste('INSDEL:',countfre,sep='')) print(names(datf)[56:58])
                  #countfre1=countfre
                  #print(paste("countfre:",countfre,"HomoZygo:",homozygo,sep=""))
                data1 <- data.frame(dat1[nn, ],
                     BNG_Freq_Perc_Filtered = as.numeric(format(round((countfre/(2*(usamp))) * 100,digits=3), scientific = FALSE)),
                    BNG_Freq_Perc_UnFiltered = as.numeric(format(round((countfreunfilt/(2*(usamp))) * 100,digits=3), scientific = FALSE)), 
                    BNG_Homozygotes=as.numeric(homozygo), stringsAsFactors = FALSE)
                  #### print(dim(data1)) print(names(data1)[56:58])
                  #### print(identical(names(data1),names(datf))) print(names(data1))
                  datf <- rbind(datf, data1)
                } else
                {
                  # dgv_match=FALSE
                  data1 <- data.frame(dat1[nn, ], BNG_Freq_Perc_Filtered = 0, BNG_Freq_Perc_UnFiltered = 0,BNG_Homozygotes= 0,stringsAsFactors = FALSE)
                  #### print(dim(data1)) print(names(data1)[56:58])
                  #### print(names(datf)[56:58])
                  #### print(identical((names(data1)[56]),(names(datf)[56])))
                  #### print(identical((names(data1)[57]),(names(datf)[57])))
                  #### print(identical((names(data1)[58]),(names(datf)[58])))
                  datf <- rbind(datf, data1)
                  # next;
                }
                
            } else if ((length(grep("inversion", variantType2[nn])) >= 1) | 
                (length(grep("translocation", variantType2[nn])) >= 1))
                {
                #print(nn)
                dat2 <- dat[which((dat$BreakPntStrt <= rf[nn] & 
                    dat$BreakPntStrt >= rf_wb_int[nn] | 
                    dat$BreakPntStrt >= rf[nn] & 
                    dat$BreakPntStrt <= rf_fb_int[nn]) & 
                    (dat$BreakPntEnd >= re[nn] & 
                    dat$BreakPntEnd <= re_fb_int[nn] | 
                    dat$BreakPntEnd <= re[nn] & 
                    dat$BreakPntEnd >= re_wb_int[nn])), ]
                
                size1 <- size_bn[nn]
                
                ## Writing if the dgv_match is TRUE or not
                countfre <- c()
                if (nrow(dat2) > 1)
                {
                  countfre <- c();countfreunfilt<-c()
                  #size_internal <- dat2$Size
                  type <- as.character(dat2$Type)
                  conf <- as.numeric(dat2$Score)
                  zygo <- as.character(dat2$Zygosity)
                  samp <- as.character(dat2$Sample)
                  
                  chrom2<-dat2$Chromosome2
                  sampid<-c();sampid_unfiltered<-c();zyg<-c();zyg_unfiltered<-c()
                  unmapped<-c()
                  homozygo<-c()
                  for (ll in 1:nrow(dat2))
                  {
                    #perc <- (size1/size_internal[ll])
                    
                   
                    if ((chrom2[ll]==chromo2[nn]) & identical(type[ll], 
                      variantType2[nn]) & ((identical(type[ll], "inversion") & 
                      conf[ll] >= invconf) | (identical(type[ll], "translocation") & 
                      conf[ll] >= transconf)))
                      {
                        countfre <- c(countfre,1)
                        sampid<-c(sampid,samp[ll])
                        zyg<-c(zyg,zygo[ll])
                        if(as.character(zygo[ll])=="homozygous"){
                            homozygo<-c(homozygo, as.character(samp[ll]))
                        }
                        else{
                            homozygo<-c(homozygo, NULL)
                        }
                      
                      
                    } 
                    else{ unmapped<-c(unmapped,ll)} 
                    ###Un filtered
                    if (identical(type[ll], variantType2[nn]) & (chrom2[ll]==chromo2[nn]))
                    {
                        countfreunfilt <- c(countfreunfilt,1)
                        sampid_unfiltered<-c(sampid_unfiltered,as.character(samp[ll]))
                        zyg_unfiltered<-c(zyg_unfiltered,as.character(zygo[ll]))
                       
                    }
                    
                  }
                  ##Filtration
                   if ((length(sampid) >= 1) & length(countfre) >= 1)
                  { 
                    
                    dat_temp <- data.frame(sampid, countfre,zyg)
                    countfre1 <- 0
                    if (nrow(dat_temp) > 0)
                    {
                      #### print(paste('dat_temp',dim(dat_temp)))
                      usampid <- as.character(unique(sampid))
                      
                        for (u in 1:length(usampid))
                        {
                            dat_temp1 <- dat_temp[which(dat_temp$sampid %in% usampid[u]), ]
                            zygo<-as.character(dat_temp1$zyg)
                            ct <- sum(dat_temp1$countfre)
                            if (ct >= 1 & nrow (dat_temp1)>= 1)
                            {
                                if(length(unique(as.character(zygo)))==1){
                                    if(unique(as.character(zygo))=="homozygous"){
                                        countfre1 <- countfre1 + 2
                                    }else if(unique(as.character(zygo))=="Unknown"){
                                        countfre1 <- countfre1 + 2
                                    }else {
                                        countfre1 <- countfre1 + 1
                                    }
                                }
                                else{
                                    g1<-grep("homozygous",as.character(zygo))
                                    g2<-grep("heterozygous",as.character(zygo))   
                                    g3<-grep("unknown",as.character(zygo))    
                                    if(length(g1)>=1 
                                        & length(g2)>=1 
                                        & length(g3)>=1){
                                        countfre1 <- countfre1 + 2
                                    }else if(length(g1)==0 
                                        & length(g2)>=1 
                                        & length(g3)>=1){
                                        countfre1 <- countfre1 + 2
                                    }else if(length(g1)>=1 
                                        & length(g2)>=1 
                                        & length(g3)==0){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else if(length(g1)>=1 
                                        & length(g2)==0 
                                        & length(g3)>=1){
                                        countfre1 <- countfre1 + 2
                                    }
                                    else{
                                        countfre1 <- countfre1 + 0
                                    }
                                }
                            }else{
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
                   if ((length(sampid_unfiltered) >= 1) & length(countfreunfilt) >= 1)
                  {
                    
                    dat_temp <- data.frame(sampid_unfiltered, countfreunfilt,zyg_unfiltered)
                    countfreunfilt1 <- 0
                    if (nrow(dat_temp) > 0)
                    {
                      #### print(paste('dat_temp',dim(dat_temp)))
                      usampid <- as.character(unique(sampid_unfiltered))
                      
                        for (u in 1:length(usampid))
                        {
                            dat_temp1 <- dat_temp[which(dat_temp$sampid_unfiltered %in% usampid[u]), ]
                            zygo<-as.character(dat_temp1$zyg_unfiltered)
                            ct <- sum(dat_temp1$countfreunfilt)
                            if (ct >= 1 & nrow (dat_temp1)>= 1)
                            {
                                if(length(unique(as.character(zygo)))==1){
                                    if(unique(as.character(zygo))=="homozygous"){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }else if(unique(as.character(zygo))=="Unknown"){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }else {
                                        countfreunfilt1 <- countfreunfilt1 + 1
                                    }
                                }
                                else{
                                    g1<-grep("homozygous",as.character(zygo))
                                    g2<-grep("heterozygous",as.character(zygo))   
                                    g3<-grep("unknown",as.character(zygo))    
                                    if(length(g1)>=1 & length(g2)>=1 & length(g3)>=1){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }else if(length(g1)==0 & length(g2)>=1 & length(g3)>=1){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }else if(length(g1)>=1 & length(g2)>=1 & length(g3)==0){
                                        countfreunfilt1 <- countfreunfilt1 + 2
                                    }
                                    else if(length(g1)>=1 & length(g2)==0 & length(g3)>=1){
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
                    }else
                        {
                          countfreunfilt1 <- 0
                        }
                    }else{
                        countfreunfilt1 <- 0
                        }
                 
                  ##### print(paste('INVTRANS:',countfre,sep=''))
                  if(length(homozygo)>=1){
                        homozygo <- length(unique(homozygo))
                    }else if(length(homozygo)==1){
                        homozygo <- length(homozygo)
                    }else{
                        homozygo <- 0
                    }
                  
                data1 <- data.frame(dat1[nn, ], 
                    BNG_Freq_Perc_Filtered = as.numeric(format(round((countfre1/(2*(usamp))) * 100,digits=3), scientific = FALSE)),
                    BNG_Freq_Perc_UnFiltered = as.numeric(format(round((countfreunfilt1/(2*(usamp)) * 100),digits=3), scientific = FALSE)), 
                    BNG_Homozygotes=as.numeric(homozygo), stringsAsFactors = FALSE)
                  #### print(dim(data1))
                  datf <- rbind(datf, data1)
                } else if (nrow(dat2) == 1)
                {
                  ## dgv_match=TRUE Calculating percentage similarity
                  countfre<-0;countfreunfilt<-0
                  chrom2<-dat2$Chromosome2
                  type <- as.character(dat2$Type)
                  conf <- dat2$Score
                  samp <- as.character(dat2$Sample)
                  zygo <- as.character(dat2$Zygosity)
                  homozygo<-c()
                if ((chrom2==chromo2[nn]) & identical(type,variantType2[nn]) & 
                    ((identical(type, "inversion") & conf >= invconf) | 
                      (identical(type, "translocation") & conf >= transconf)))
                      {
                       if(zygo =="homozygous"){
                        countfre <- countfre+2
                        homozygo<-as.character(samp)
                       }else if(zygo =="unknown"){
                        countfre <- countfre+2
                       }else{ 
                        countfre <- countfre+1
                       }
                    
                    } else
                    {
                        countfre <-  0
                    
                    }
                  ##Unfiletered Frequency
                if (identical(type, variantType2[nn]) & (chrom2==chromo2[nn])){
                    if(zygo =="homozygous"){
                        countfreunfilt <- countfreunfilt + 2
                        
                    }else if(zygo =="unknown"){
                        countfreunfilt <- countfreunfilt + 2
                    }else{ 
                        countfreunfilt <- countfreunfilt + 1
                    }
                }
                  ##### print(paste('INVTRANS:',countfre,sep=''))
                  ##Calculating the length
                if(length(homozygo)>0){
                    homozygo <- length(homozygo)
                }else{
                    homozygo = 0
                }
                  #countfre1=countfre
                  #print(paste("countfre:",countfre,"HomoZygo:",homozygo,sep=""))
                data1 <- data.frame(dat1[nn, ],
                    BNG_Freq_Perc_Filtered = as.numeric(format(round((countfre/(2*(usamp))) * 100,digits=3), scientific = FALSE)),
                    BNG_Freq_Perc_UnFiltered = as.numeric(format(round((countfreunfilt/(2*(usamp))) * 100,digits=3), scientific = FALSE)), 
                    BNG_Homozygotes=as.numeric(homozygo), stringsAsFactors = FALSE)
                  
                  #### print(dim(data1))
                  datf <- rbind(datf, data1)
                } else
                {
                  # dgv_match=FALSE print(paste('QueryData:',dim(dat1[nn,]),sep=''))
                data1 <- data.frame(dat1[nn, ], 
                    BNG_Freq_Perc_Filtered= 0, BNG_Freq_Perc_UnFiltered =0, 
                    BNG_Homozygotes= 0, stringsAsFactors = FALSE)
                  #### print(dim(data1))
                  datf <- rbind(datf, data1)
                  # next;
                }
                
            } else
            {
                #### print(paste('QueryData:',dim(dat1[nn,]),sep=''))
                
                data1 <- data.frame(dat1[nn, ], BNG_Freq_Perc_Filtered= 0, BNG_Freq_Perc_UnFiltered =0, BNG_Homozygotes= 0,
                    stringsAsFactors = FALSE)
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



