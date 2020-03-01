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
#' @param EnzymeType Character. Type of enzyme. Options Dual and DLE.
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
#' @param indexfile File containing connection between sample and nanoIDs
#' @param returnMethod character. Choice between Text and DataFrame.
#' @return Text file or data frames containing internalFrequency data.
#' Required if you want to calculate internal frequency.
#' @examples
#' mergedFiles <- system.file("extdata", "nanotatoR_merged.txt", package = "nanotatoR")
#' smapName <- "F1.1_TestSample1_solo_hg19.smap"
#' smappath <- system.file("extdata", package = "nanotatoR")
#' indelconf = 0.5; invconf = 0.01;transconf = 0.1;input_fmt="Text";
#' internalFrequency(mergedFiles = mergedFiles, smappath = smappath , smap = smapName, 
#' buildSVInternalDB=FALSE, win_indel=10000, 
#' win_inv_trans=50000, 
#' perc_similarity=0.5, indelconf=0.5, invconf=0.01, 
#' transconf=0.1, limsize=1000, win_indel_parents=5000,input_fmt="Text",
#' win_inv_trans_parents=40000,
#' returnMethod="dataFrame")
#' @importFrom stats na.omit 
#' @import hash
#' @export
internalFrequency_DLE_Solo <- function(mergedFiles, smappath , smap , 
    buildSVInternalDB=FALSE, smapdata, input_fmt = c("Text","dataFrame"), 
    path, pattern, outpath, win_indel=10000, win_inv_trans=50000, 
    perc_similarity=0.5, indelconf=0.5, invconf=0.01, fname, 
    limsize=1000, win_indel_parents=5000,win_inv_trans_parents=40000, 
    transconf=0.1, dbOutput =c("dataframe","text"),
    returnMethod=c("Text","dataFrame"),
    indexfile = "X:/Hayks_Materials/BNG/Projects/Sample_NIDKeys.csv",
	EnzymeType = c("Dual", "DLE"))
    {
    #library(hash)
    print("###Calculating the Internal Frequency###")
    if(buildSVInternalDB == TRUE){
        r <- Merging_DLE_Dual(
            labelType = labelType, 
            Dual_path = Dual_path, 
            Dual_pattern = Dual_pattern, 
            DLE_path = DLE_path, 
            DLE_pattern = DLE_pattern, 
            Samplecodes = Samplecodes, 
            mergeKey = mergeKey, 
            outpath = outpath, 
            filename = filename,
	        mergedKeyoutpath = mergedKeyoutpath,
	        mergedKeyFname = mergedKeyFname)
        }
        else{
            r <- read.table(mergedFiles, sep = "\t", header = TRUE)
        }
    usamp <- length(unique(r$nanoID))
    if(input_fmt == "Text"){
		if(EnzymeType == "Dual"){
            r1 <- readSMap(smap, input_fmt_smap)
        }
        else{
            r1 <- readSMap_DLE(smap, input_fmt_smap)
        }
    }
	else if(input_fmt == "dataFrame"){
		r1<-smapdata
	}
	else{
		stop("Input Format incorrect")
	}
    sampID <- str_squish(as.character(unique(r1$SampleID)))
    
    
    ufam <- as.character(unique(r$nanoID))
    famid <- c()
    for (qp in seq_len(length(ufam)))
    {
        stt <- strsplit(ufam[[qp]][1], split = "[.]")
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
       
        dat <- r[which(r$RefcontigID1 == ii), ]
        #perc_similarity_parents = perc_similarity_parents
		perc_similarity = perc_similarity
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
        #print(indexfile)
        if(buildSVInternalDB == FALSE){
            famindexfile <-read.csv(indexfile)
		}
		else{
		    famindexfile <- read.csv(file.path(
			    mergedKeyoutpath = mergedKeyoutpath, 
				mergedKeyFname = mergedKeyFname
				))
			}
        ha1 <- hash()
        .set(ha1, keys = famindexfile$SampleID, values =  famindexfile$NID)
        val1 <-as.character(hash::values(ha1, keys = sampID))
        stt1 <- strsplit(val1, split = "\\.")
        patID <- stt1[[1]][2]
        svfamid <- stt1[[1]][1]
		ha2 <- hash()
        .set(ha2, keys = famindexfile$SampleID, values =  famindexfile$Tag)
        'if((length(grep("\\\\",smap))>=1)){
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
        }'
        # svind<-dat1$SVIndex
        bsssi<-c()
        bspqi<-c()
        rest<-c()
        # conf<-dat$Confidence countfre<-0 percn<-c()
        datf <- c()
		print(paste('Chrom:',ii,sep=''))
		#for (nn in 1:20)
        for (nn in seq_along(rf))
		{
            
            #print(paste("nn:",nn)) 
           
            if ((variantType2[nn] == "deletion" | 
                variantType2[nn] == "insertion")){ 
                        
                dat2 <-dat[which((dat$RefStartPos >= rf_wb_ind[nn] & dat$RefEndPos <= re_wf_ind[nn])),]
                    
                size1 <- size_bn[nn]
                ## Calculating Internal Frequency
                if (nrow(dat2) > 1)
                {
                    dat2$perc_ref_query <- as.numeric(dat2$Size)/size1
					dat2$perc_query_ref <- size1/as.numeric(dat2$Size)
					stt12 <- str_split(as.character(dat2$nanoID), pattern = "\\.")
					famid <- c(); relnid <- c()
					
					'for(ee in seq_along(stt12)){
					    famid <- c(famid, as.character(stt12[[ee]][1]))
					    relnid <- c(relnid, as.character(stt12[[ee]][2]))
					}
					'
					dat2$FamilyID <- sapply(stt12, function(x) x[1])
					dat2$RelationID <- sapply(stt12, function(x) x[2])
					'dat2MomEqual <- subset(dat2,(((RefStartPos ==rf[nn]) & (RefEndPos==re[nn])) & 
                            (perc_ref_query >= perc_similarity & perc_query_ref >= perc_similarity) 
							& (as.character(Type) == variantType2[nn])
                            & (as.character(FamilyID) == svfamid)
							& (patID == 1 & RelationID == 2)))
					dat2DadEqual <- subset(dat2,(((RefStartPos ==rf[nn]) & (RefEndPos==re[nn])) & 
                            (perc_ref_query >= perc_similarity & perc_query_ref >= perc_similarity)
							& (as.character(Type) == variantType2[nn])
                            & (as.character(FamilyID) == svfamid)
							& (patID == 1 & RelationID == 3)))
					dat2MomRange <- subset(dat2,((( RefStartPos >= rf_wb_int_parents[nn]) & (RefEndPos <= re_fb_int_parents[nn])) & (perc_ref_query >= perc_similarity & perc_query_ref >= perc_similarity) 
						& (as.character(Type) == variantType2[nn])
                            & (as.character(FamilyID) == svfamid)
							& (patID == 1 & RelationID == 2)))
					dat2MomRange <- subset(dat2,((((dat2$RefStartPos <= rf[nn] & 
					    dat2$RefStartPos >= rf_wb_ind_parents[nn]) | 
						(dat2$RefStartPos >= rf[nn] & dat2$RefStartPos <= rf_wb_ind_parents[nn])) 
						& ((dat2$RefEndPos >= re[nn] & dat2$RefStartPos <= re_wf_ind_parents[nn]) 
						| (dat2$RefEndPos <= re[nn] & dat2$RefEndPos >= re_wf_ind_parents[nn]))) 
						& (dat2$perc >= perc_similarity) 
						& (as.character(dat2$Type) == variantType2[nn])
                            & (as.character(dat2$FamilyID) == svfamid)
							& (patID == 1 & dat2$RelationID == 2)))
					dat2DadRange <- subset(dat2,((( RefStartPos >= rf_wb_int_parents[nn]) & (RefEndPos <= re_fb_int_parents[nn])) 
						& (perc_ref_query >= perc_similarity & perc_query_ref >= perc_similarity)
						& (as.character(Type) == variantType2[nn])
                            & (as.character(FamilyID) == svfamid)
							& (patID == 1 & RelationID == 3)))'
					dat2Frequency <- subset(dat2,((
					    perc_ref_query >= perc_similarity 
						& perc_query_ref >= perc_similarity) 
					    & (Size >= limsize)
						& (as.character(Type) == variantType2[nn])
                        & (as.character(FamilyID) != svfamid)
						& (Confidence >= indelconf)
						& ((Found_in_self_BSPQI_molecules == "yes" & 
                        Found_in_self_BSSSI_molecules == "yes") 
                        | (Found_in_self_BSPQI_molecules == "no" 
						& Found_in_self_BSSSI_molecules == "yes") 
						| (Found_in_self_BSPQI_molecules == "yes" 
						& Found_in_self_BSSSI_molecules == "no") 
						| (Found_in_self_BSPQI_molecules == "yes" 
						& Found_in_self_BSSSI_molecules == "-") 
						| (Found_in_self_BSPQI_molecules == "-" 
						& Found_in_self_BSSSI_molecules == "yes")
						| (Found_in_DLE_self_molecules == "yes"))
						& ((Type1 =="insertion" & Type2 =="insertion") 
						| (Type1 =="insertion" & is.na(Type2)) 
						| (is.na(Type1) & Type2 =="insertion") 
						| (Type1 == "deletion" & Type2 == "deletion") 
						| (Type1 == "deletion" & is.na(Type2)) 
						| (is.na(Type1) & Type2 =="deletion")
						| (Type1 == "-" & Type2 == "-" & Type == "deletion")
						| (Type1 == "-" & Type2 == "-" & Type == "insertion")
						| (Type == "insertion" | Type == "deletion"))))
					dat2UnfilteredFrequency <- subset(dat2,((
					    perc_ref_query >= perc_similarity & 
						perc_query_ref >= perc_similarity) & 
						(as.character(Type) == variantType2[nn])
                        & (as.character(FamilyID) != svfamid)))
						 
					
					
					'countfre <- c();countfreunfilt<-c()
                    sv1 <- dat2$nanoID
					
                    # svind1<-dat2$SVIndex
                    #sv1 <- strsplit(as.character(svfam1), split = "_")
                    size_internal <- dat2$Size
                    zygo <- as.character(dat2$Zygosity)
                    type <- as.character(dat2$Type)
                    typ1 <- as.character(dat2$Type1)
                    typ2 <- as.character(dat2$Type2)
                    #svid <- as.numeric(dat2$SVIndex)
                    \'if (length(grep("SVIndex",names(smap)))>0){
                        svid <- dat2$SVIndex
                    }
                    else{
                        svid <- dat2$SmapEntryID
                    }\'
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
						print(perc)
						print (sv1[ll])
                        #### print(perc) print(type[ll])
                        stt <- strsplit(as.character(sv1[ll]), split = "[.]")
                        patID1 <- stt[[1]][2]
                        svfamid1 <- stt[[1]][1]
						#### print(svfamid1)
                        ###Filtered INDEL
                        ##Check whether the parents have same breakpoints 
                        ##as the proband
                    
                        
                        if (((svStrt[ll]==rf[nn]) & (svEnd[ll]==re[nn])) & 
                            (perc_ref_query >= perc_similarity & perc_query_ref >= perc_similarity) & 
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
                        else if ((perc_ref_query >= perc_similarity & perc_query_ref >= perc_similarity) & 
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
                            else if ((perc_ref_query >= perc_similarity & perc_query_ref >= perc_similarity) & 
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
								print(countfre)
                                svfamid3 <- c(svfamid3, as.character(sv1[ll]))
                                motherZygosity <- c(motherZygosity, "-")
                                fatherZygosity <- c(fatherZygosity, "-")
                                svSAMP<-c(svSAMP,ll)
                                if(as.character(zygo[ll])=="homozygous"){
                                    homozygo<-c(homozygo, 
                                    as.character((sv1[ll])))
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
                                    as.character(sv1[ll]))
                                svSAMP_unfiltered<-c(svSAMP_unfiltered,ll)
                            }
                            else{
                               \'print ("SVs does not meet the unfiltered 
                                    criteria!!!")\'
                            }
                      
                     
                    
                        }'
					##Filtered
                    if (nrow(dat2Frequency) > 0){
					    dat2Homozygotes <- data.frame(dat2Frequency %>% 
						mutate(Homozygo = Zygosity=="homozygous") %>%
						group_by(nanoID) %>%
						summarise (Homozygotes = sum(Homozygo)))
						cntHomozygotes <- sum(dat2Homozygotes$Homozygo)
						dat2Homozygotes [,2] <- ifelse(dat2Homozygotes[,2] > 0, 2, 0)
						dat2Heterozygotes <- data.frame(dat2Frequency %>% 
						mutate(Heterozygo = Zygosity == "heterozygous") %>%
						group_by(nanoID) %>%
						summarise (Heterozygotes = sum(Heterozygo)))
						dat2Heterozygotes [,2] <- ifelse(dat2Heterozygotes[,2] > 0, 1, 0)
						dat2Unknown <- data.frame(dat2Frequency %>% 
						mutate(Unknown = Zygosity == "unknown") %>%
						group_by(nanoID) %>%
						summarise (Unknown = sum(Unknown)))
						dat2Unknown [,2] <- ifelse(dat2Unknown[,2] > 0, 2, 0)
						dat2filtFreq <- data.frame(NR = dat2Homozygotes[,1], Homozygote = dat2Homozygotes[,2],
						Heterozygotes = dat2Heterozygotes[,2], Unknown = dat2Unknown[,2])
						row.names(dat2filtFreq) <- as.character(dat2Homozygotes[,1])
						dat2filtFreq$RowSum <- ifelse((dat2filtFreq$Heterozygotes ==1 & (dat2filtFreq$Homozygote == 0 & dat2filtFreq$Unknown == 0)), 1, 2)
						filtFreq <- sum(dat2filtFreq$RowSum)
						} else{
						    filtFreq  = 0
							cntHomozygotes = 0
                        }
						
					
					###Unfiltered
					
					if (nrow(dat2UnfilteredFrequency) > 0){
					    dat2Homozygotes <- data.frame(dat2UnfilteredFrequency %>% 
						mutate(Homozygo=Zygosity=="homozygous") %>%
						group_by(nanoID) %>%
						summarise (Homozygotes = sum(Homozygo)))
						#cntHomozygotes <- sum(dat2Homozygotes$Homozygotes)
						dat2Homozygotes [,2] <- ifelse(dat2Homozygotes[,2] > 0, 2, 0)
						dat2Heterozygotes <- data.frame(dat2UnfilteredFrequency %>% 
						mutate(Heterozygo = Zygosity == "heterozygous") %>%
						group_by(nanoID) %>%
						summarise (Heterozygotes = sum(Heterozygo)))
						dat2Heterozygotes [,2] <- ifelse(dat2Heterozygotes[,2] > 0, 1, 0)
						dat2Unknown <- data.frame(dat2UnfilteredFrequency %>% 
						mutate(Unknown = Zygosity == "unknown") %>%
						group_by(nanoID) %>%
						summarise (Unknown = sum(Unknown)))
						dat2Unknown [,2] <- ifelse(dat2Unknown[, 2] > 0, 2, 0)
						dat2unfiltFreq <- data.frame(Homozygote = dat2Homozygotes[,2],
						Heterozygotes = dat2Heterozygotes[,2], Unknown = dat2Unknown[,2])
						row.names(dat2unfiltFreq) <- as.character(dat2Homozygotes[,1])
						dat2unfiltFreq$RowSum <- ifelse((dat2unfiltFreq$Heterozygotes ==1 & (dat2unfiltFreq$Homozygote == 0 & dat2unfiltFreq$Unknown == 0)), 1, 2)
						UnfilteredFrequency  = sum(dat2unfiltFreq$RowSum)
						} else{
						    UnfilteredFrequency = 0
						}
						
                  
                    ##Calculating filtered Frequency INDEL
                    ##Zygosity Parent
					 '   if(length(dat2DadEqual$Zygosity) > 0){
						    if(length(as.character(unique
                            (dat2DadEqual$Zygosity))) == 1){
							    fatherZygosity<-as.character(unique(dat2DadEqual$Zygosity))
							}else {
							    fatherZygosity <- paste(
								as.character(dat2DadEqual$Zygosity), 
								collapse = ","
								)
							}
						} else if (length(dat2DadRange$Zygosity) > 0){
					    fatherzygosity <- as.character(dat2DadRange$Zygosity)
                        fatherZygosity1 <- gsub("-", NA, fatherzygosity)
                        fatherZygosity1 <- as.character(
						        na.omit(fatherZygosity1)
                                )
                        if(length(fatherZygosity1) == 1){
                            fatherZygosity<-as.character(
							unique(fatherZygosity1)
							)
                        } 
                        else {
                            fatherZygosity <- paste(
								as.character(fatherZygosity1), 
								collapse = ","
								)
                        }
						} else {
						    fatherZygosity <- "-"
							}
                       
                        ###Mother Zygosity
                        if(length(dat2MomEqual$Zygosity) > 0){
						    if(length(as.character(unique
                            (dat2MomEqual$Zygosity))) == 1){
							    motherZygosity<-as.character(unique(dat2MomEqual$Zygosity))
							}else {
							    motherZygosity <- paste(
								as.character(dat2MomEqual$Zygosity), 
								collapse = ","
								)
							}
						} else if (length(dat2MomRange$Zygosity) > 0){
					    motherZygosity <- as.character(dat2MomRange$Zygosity)
                        motherZygosity1 <- gsub("-", NA, motherZygosity)
                        motherZygosity1 <- as.character(
						        na.omit(motherZygosity1)
                                )
                        if(length(motherZygosity1) == 1){
                            motherZygosity<-as.character(
							unique(motherZygosity1)
							)
                        } 
                        else {
                            motherZygosity <- paste(
								as.character(motherZygosity1), 
								collapse = ","
								)
                        }
						} else {
						    motherZygosity <- "-"
							}
                       
                  
                        if(length(homozygo)>=1){
                            homozygo <- length(unique(homozygo))
                        }
                        else if(length(homozygo)==1){
                            homozygo <- length(homozygo)
                        }
                        else{
                            homozygo <- 0
                        }'
                                        
                  
                        famno <- unname(hash::values(ha2, keys = sampID))
						fampas <- paste("^",as.character(famindexfile$Tag),"$",sep="")
						famnopas <- paste("^",as.character(famno),"$",sep="")
                        famnumb <- length(grep(famnopas, fampas, fixed=TRUE))
                        data1 <- data.frame(dat1[nn, ], 
                            Internal_Freq_Perc_Filtered = 
                            as.numeric(format(round(
                            (filtFreq/(2*(usamp - famnumb))) * 100,digits=3)
                            , scientific = FALSE)),
                            Internal_Freq_Perc_Unfiltered = 
                            as.numeric
                            (format(round((UnfilteredFrequency/
                            (2*(usamp - famnumb))) * 100,digits=3), 
                            scientific = FALSE)),
                            Internal_Homozygotes=as.numeric(cntHomozygotes), 
                            MotherZygosity = "-", 
                            FatherZygosity = "-",  
                            stringsAsFactors = FALSE)
                        
                        datf <- rbind(datf, data1)
                }  
                else if (nrow(dat2) == 1){
                    # dgv_match=TRUE Calculating percentage similarity
                    countfre <- 0;countfreunfilt<-0
                    conf <- dat2$Confidence
                    size_internal <- dat2$Size
                    BSPQI_status_DB <- as.character(dat2$Found_in_self_BSPQI_molecules)
                    BSSSI_status_DB <- as.character(dat2$Found_in_self_BSSSI_molecules)
                    perc_query_ref <- (size1/size_internal)
					perc_ref_query <- (size_internal/size1)
                    sv1 <- dat2$nanoID
					Found_in_DLE_self_molecules <- dat2$Found_in_DLE_self_molecules
                    zygo <- as.character(dat2$Zygosity)
                    stt <- strsplit(as.character(sv1), split = "[.]")
                    patID1 <- stt[[1]][2]
                    svfamid1 <- stt[[1]][1]
                    motherZygosity <- ""
                    fatherZygosity <- ""
                    type <- as.character(dat2$Type)
                    typ1 <- as.character(dat2$Type1)
                    typ2 <- as.character(dat2$Type2)
                    homozygo<-c()
                    ##Check for parents
                    'if ((perc_ref_query >= perc_similarity & perc_query_ref >= perc_similarity) & 
                        identical(type, variantType2[nn]) & 
                        (identical(svfamid1, svfamid)) &
                        ((dat2$RefStartPos >= rf_wb_int_parents[nn]) & (dat2$RefEndPos <= re_fb_int_parents[nn]))){
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
                    else'
					if ((perc_ref_query >= perc_similarity 
					    & perc_query_ref >= perc_similarity) & 
                        identical(type, variantType2[nn]) & 
                        !(identical(svfamid1, svfamid)) & 
                        ((size_internal>=limsize)) & 
                        ((BSPQI_status_DB == "yes" & BSSSI_status_DB == "yes") 
                        | (BSPQI_status_DB == "no" & BSSSI_status_DB == "yes") 
                        | (BSPQI_status_DB == "yes" & BSSSI_status_DB == "no")
                        | (BSPQI_status_DB == "yes" & BSSSI_status_DB == "-") |
                        (BSPQI_status_DB == "-" & BSSSI_status_DB == "yes") | 
                        (Found_in_DLE_self_molecules == "yes")) & 
                        (conf > indelconf) & 
                        ((typ1=="insertion" & typ2=="insertion") | 
                        (typ1=="insertion" & is.na(typ2)) | 
                        (is.na(typ1) & typ2=="insertion") | 
                        (typ1=="deletion" & typ2=="deletion") |
                        (typ1=="deletion" & is.na(typ2)) | 
                        (is.na(typ1) & typ2=="deletion")
						| (typ1=="-" & typ2=="-" & type == "deletion")
						| (typ1=="-" & typ2=="-" & type == "insertion")
						| (type == "insertion" | type == "deletion"))){
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
                    if ((perc_ref_query >= perc_similarity & 
					        perc_query_ref >= perc_similarity) & 
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
                    Internal_Freq_Perc_Filtered = as.numeric(format(round((
					countfre/(2*(usamp - famno))) * 100,digits=3), 
					scientific = FALSE)),
                    Internal_Freq_Perc_Unfiltered = as.numeric(format(
					round((countfreunfilt/(2*(usamp - famno))) * 100,digits=3), scientific = FALSE)), 
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
                
            } else if ((variantType2[nn] == "duplication" |
                    variantType2[nn] == "duplication_split" | 
                    variantType2[nn] == "duplication_inverted" |
                    variantType2[nn] == "insertion"))
                {   
                   dat2 <-dat[which((dat$RefStartPos >= rf_wb_ind[nn] & dat$RefEndPos <= re_wf_ind[nn])),]
                    size1 <- size_bn[nn]
                    ## Calculating Internal Frequency
                    if (nrow(dat2) > 1){ 
					dat2$perc_ref_query <- as.numeric(dat2$Size)/size1
					dat2$perc_query_ref <- size1/as.numeric(dat2$Size)
					stt12 <- str_split(as.character(dat2$nanoID), pattern = "\\.")
					famid <- c(); relnid <- c()
					dat2$FamilyID <- sapply(stt12, function(x) x[1])
					dat2$RelationID <- sapply(stt12, function(x) x[2])
					'dat2MomEqual_duplication <- subset(dat2,(((RefStartPos ==rf[nn]) & (RefEndPos==re[nn]))  
                            & (as.character(Type) == variantType2[nn])
                            & (as.character(FamilyID) == svfamid)
							& (patID == 1 & RelationID == 2) & ((as.character(Type) == "duplication") |
						    (as.character(Type) == "duplication_split") |
							(as.character(Type) == "duplication_inverted")
						)))
					dat2DadEqual_duplication <- subset(dat2,(((RefStartPos ==rf[nn]) & (RefEndPos==re[nn]))  
                            & (as.character(Type) == variantType2[nn])
                            & (as.character(FamilyID) == svfamid)
							& (patID == 1 & RelationID == 3) & ((as.character(Type) == "duplication") |
						    (as.character(Type) == "duplication_split") |
							(as.character(Type) == "duplication_inverted")
						)))
					dat2MomRange_duplication <- subset(dat2,((( RefStartPos >= rf_wb_ind_parents[nn]) & (RefEndPos <= re_wf_ind_parents[nn]))
						& (as.character(Type) == variantType2[nn])
                            & (as.character(FamilyID) == svfamid)
							& (patID == 1 & RelationID == 2)& ((as.character(Type) == "duplication") |
						    (as.character(Type) == "duplication_split") |
							(as.character(Type) == "duplication_inverted")
						)))
					dat2DadRange_duplication <- subset(dat2,((( RefStartPos >= rf_wb_ind_parents[nn]) & (RefEndPos <= re_wf_ind_parents[nn])) 
						& (as.character(Type) == variantType2[nn])
                            & (as.character(FamilyID) == svfamid)
							& (patID == 1 & RelationID == 3) & ((as.character(Type) == "duplication") |
						    (as.character(Type) == "duplication_split") |
							(as.character(Type) == "duplication_inverted")
						)))
					dat2MomEqual_insertion <- subset(dat2,(((RefStartPos ==rf[nn]) & (RefEndPos==re[nn])) & 
                            (as.character(Type) == variantType2[nn])
							& (perc_ref_query >= perc_similarity & perc_query_ref >= perc_similarity) 
                            & (as.character(FamilyID) == svfamid)
							& (patID == 1 & RelationID == 2) & ((as.character(Type) == "insertion")
						)))
					dat2DadEqual_insertion <- subset(dat2,(((RefStartPos ==rf[nn]) & (RefEndPos==re[nn])) & 
                            (as.character(Type) == variantType2[nn])
							& (perc_ref_query >= perc_similarity & perc_query_ref >= perc_similarity) 
                            & (as.character(FamilyID) == svfamid)
							& (patID == 1 & RelationID == 3) & ((as.character(Type) == "insertion")
						)))
					dat2MomRange_insertion <- subset(dat2,((( RefStartPos >= rf_wb_ind_parents[nn]) & (RefEndPos <= re_wf_ind_parents[nn])) 
						& (perc_ref_query >= perc_similarity & perc_query_ref >= perc_similarity) 
						& (as.character(Type) == variantType2[nn])
                            & (as.character(FamilyID) == svfamid)
							& (patID == 1 & RelationID == 2)& ((as.character(Type) == "insertion")
						)))
					dat2DadRange_insertion <- subset(dat2,(((RefStartPos >= rf_wb_ind_parents[nn]) & (RefEndPos <= re_wf_ind_parents[nn]))  
						& (perc_ref_query >= perc_similarity & perc_query_ref >= perc_similarity) 
						& (as.character(Type) == variantType2[nn])
                            & (as.character(FamilyID) == svfamid)
							& (patID == 1 & RelationID == 3) & ((as.character(Type) == "insertion")
						)))'
					dat2FrequencyDuplication <- subset(dat2,(((as.character(Type) == "duplication") |
						    (as.character(Type) == "duplication_split") |
							(as.character(Type) == "duplication_inverted"))
						    & (perc_ref_query >= perc_similarity 
						    & perc_query_ref >= perc_similarity)
						    & (as.character(FamilyID) != svfamid)
						    &((Found_in_self_BSPQI_molecules == "yes" & 
                            Found_in_self_BSSSI_molecules == "yes") 
                            | (Found_in_self_BSPQI_molecules == "no" 
							& Found_in_self_BSSSI_molecules == "yes") 
							| (Found_in_self_BSPQI_molecules == "yes" 
							& Found_in_self_BSSSI_molecules == "no") 
							| (Found_in_self_BSPQI_molecules == "yes" 
							& Found_in_self_BSSSI_molecules == "-") 
							| (Found_in_self_BSPQI_molecules == "-" 
							& Found_in_self_BSSSI_molecules == "yes")
							| (Found_in_DLE_self_molecules == "yes"))
							& ((Type1 =="duplication" & 
                                Type2 =="duplication") | 
                                (Type1 =="duplication" & is.na(Type2)) |
                                (is.na(Type1) & Type2 =="duplication") 
								| (Type1 == "duplication_split" & Type2 == "duplication_split") |
                                (Type1 == "duplication_split" & is.na(Type2)) | 
                                (is.na(Type1) & Type2 =="duplication_split")
								| (Type1 == "duplication_inverted" & Type2 == "duplication_inverted") |
                                (Type1 == "duplication_inverted" & is.na(Type2)) | 
                                (is.na(Type1) & Type2 =="duplication_inverted")
								| (Type1 == "-" & Type2 == "-" & Type == "duplication_inverted")
								| (Type1 == "-" & Type2 == "-" & Type == "duplication")
								| (Type1 == "-" & Type2 == "-" & Type == "duplication_split")
								|((as.character(Type) == "duplication") |
						         (as.character(Type) == "duplication_split") |
							     (as.character(Type) == "duplication_inverted")
						    ))
							& ((Fail_BSPQI_assembly_chimeric_score == "pass" & 
                            Fail_BSSSI_assembly_chimeric_score == "pass") 
                            | (Fail_BSPQI_assembly_chimeric_score == "fail" 
							& Fail_BSSSI_assembly_chimeric_score == "pass") 
							| (Fail_BSPQI_assembly_chimeric_score == "pass" 
							& Fail_BSSSI_assembly_chimeric_score == "fail") 
							| (Fail_BSPQI_assembly_chimeric_score == "pass" 
							& Fail_BSSSI_assembly_chimeric_score == "-") 
							| (Fail_BSPQI_assembly_chimeric_score == "-" 
							& Fail_BSSSI_assembly_chimeric_score == "pass")
							|(Fail_assembly_chimeric_score_DLE == "pass"))))
					dat2Insertion <- subset(dat2,(
					    (perc_ref_query >= perc_similarity 
						& perc_query_ref >= perc_similarity) 
					    & (Size >= limsize)
						& (as.character(Type) == variantType2[nn])
                        & (as.character(FamilyID) != svfamid)
						& (Confidence >= indelconf)
						& ((Found_in_self_BSPQI_molecules == "yes" & 
                            Found_in_self_BSSSI_molecules == "yes") 
                            | (Found_in_self_BSPQI_molecules == "no" 
							& Found_in_self_BSSSI_molecules == "yes") 
							| (Found_in_self_BSPQI_molecules == "yes" 
							& Found_in_self_BSSSI_molecules == "no") 
							| (Found_in_self_BSPQI_molecules == "yes" 
							& Found_in_self_BSSSI_molecules == "-") 
							| (Found_in_self_BSPQI_molecules == "-" 
							& Found_in_self_BSSSI_molecules == "yes")
							| (Found_in_DLE_self_molecules == "yes") )
							& ((Type1 =="insertion" & 
                                Type2 =="insertion") | 
                                (Type1 =="insertion" & is.na(Type2)) |
                                (is.na(Type1) & Type2 =="insertion")
							| (Type1 == "-" & Type2 == "-" & Type == "insertion")
							| (Type == "insertion"))))
					dat2UnfilteredFrequency <- subset(dat2,(
					        ((as.character(Type) == "duplication") |
						    (as.character(Type) == "duplication_split") |
							(as.character(Type) == "duplication_inverted") |
							(as.character(Type) == "insertion"))
                            & (perc_ref_query >= perc_similarity 
						    & perc_query_ref >= perc_similarity)							
							& (as.character(FamilyID) != svfamid)))
					if (nrow(dat2FrequencyDuplication) > 0 & 
					        nrow(dat2Insertion) > 0){
							dat2Frequency <- rbind(
							    dat2FrequencyDuplication,
								dat2Insertion
								)
							}
					else if(nrow(dat2FrequencyDuplication) > 0 & 
					        nrow(dat2Insertion) == 0){
							dat2Frequency <- dat2FrequencyDuplication
							}
					else if(nrow(dat2FrequencyDuplication) == 0 & 
					        nrow(dat2Insertion) > 0){
							dat2Frequency <- dat2Insertion
							}
							
							
							
							
                        ##Calculating filtered Frequency INDEL
                        ##Filtration
                 
                        ##Filtered
                    if (nrow(dat2Frequency) > 0){
					    
					    countfre1 <- 0
						dat2Homozygotes <- data.frame(dat2Frequency %>% 
						mutate(Homozygo = Zygosity=="homozygous") %>%
						group_by(nanoID) %>%
						summarise (Homozygotes = sum(Homozygo)))
						cntHomozygotes <- sum(dat2Homozygotes$Homozygo)
						dat2Homozygotes [,2] <- ifelse(dat2Homozygotes[,2] > 0, 2, 0)
						dat2Heterozygotes <- data.frame(dat2Frequency %>% 
						mutate(Heterozygo = Zygosity == "heterozygous") %>%
						group_by(nanoID) %>%
						summarise (Heterozygotes = sum(Heterozygo)))
						dat2Heterozygotes [,2] <- ifelse(dat2Heterozygotes[,2] > 0, 1, 0)
						dat2Unknown <- data.frame(dat2Frequency %>% 
						mutate(Unknown = Zygosity == "unknown") %>%
						group_by(nanoID) %>%
						summarise (Unknown = sum(Unknown)))
						dat2Unknown [,2] <- ifelse(dat2Unknown[,2] > 0, 2, 0)
						dat2filtFreq <- data.frame(NR = dat2Homozygotes[,1], Homozygote = dat2Homozygotes[,2],
						Heterozygotes = dat2Heterozygotes[,2], Unknown = dat2Unknown[,2])
						row.names(dat2filtFreq) <- as.character(dat2Homozygotes[,1])
						dat2filtFreq$RowSum <- ifelse((dat2filtFreq$Heterozygotes ==1 & (dat2filtFreq$Homozygote == 0 & dat2filtFreq$Unknown == 0)), 1, 2)
						filtFreq <- sum(dat2filtFreq$RowSum)
						} else{
						    filtFreq  = 0
							cntHomozygotes = 0
                        }
						
					
					###Unfiltered
					
					if (nrow(dat2UnfilteredFrequency) > 0){
					    dat2Homozygotes <- data.frame(dat2UnfilteredFrequency %>% 
						mutate(Homozygo=Zygosity=="homozygous") %>%
						group_by(nanoID) %>%
						summarise (Homozygotes = sum(Homozygo)))
						#cntHomozygotes <- sum(dat2Homozygotes$Homozygotes)
						dat2Homozygotes [,2] <- ifelse(dat2Homozygotes[,2] > 0, 2, 0)
						dat2Heterozygotes <- data.frame(dat2UnfilteredFrequency %>% 
						mutate(Heterozygo = Zygosity == "heterozygous") %>%
						group_by(nanoID) %>%
						summarise (Heterozygotes = sum(Heterozygo)))
						dat2Heterozygotes [,2] <- ifelse(dat2Heterozygotes[,2] > 0, 1, 0)
						dat2Unknown <- data.frame(dat2UnfilteredFrequency %>% 
						mutate(Unknown = Zygosity == "unknown") %>%
						group_by(nanoID) %>%
						summarise (Unknown = sum(Unknown)))
						dat2Unknown [,2] <- ifelse(dat2Unknown[, 2] > 0, 2, 0)
						dat2unfiltFreq <- data.frame(Homozygote = dat2Homozygotes[,2],
						Heterozygotes = dat2Heterozygotes[,2], Unknown = dat2Unknown[,2])
						row.names(dat2unfiltFreq) <- as.character(dat2Homozygotes[,1])
						dat2unfiltFreq$RowSum <- ifelse((dat2unfiltFreq$Heterozygotes ==1 & (dat2unfiltFreq$Homozygote == 0 & dat2unfiltFreq$Unknown == 0)), 1, 2)
						UnfilteredFrequency  = sum(dat2unfiltFreq$RowSum)
						} else{
						    UnfilteredFrequency = 0
						}
						
                 
                    ##Calculating filtered Frequency INDEL
                    ##Zygosity Parent
					    'if(length(dat2DadEqual$Zygosity) > 0){
						    if(length(as.character(unique
                            (dat2DadEqual$Zygosity))) == 1){
							    fatherZygosity<-as.character(unique(dat2DadEqual$Zygosity))
							}else {
							    fatherZygosity <- paste(
								as.character(dat2DadEqual$Zygosity), 
								collapse = ","
								)
							}
						} else if (length(dat2DadRange$Zygosity) > 0){
					    fatherzygosity <- as.character(dat2DadRange$Zygosity)
                        fatherZygosity1 <- gsub("-", NA, fatherzygosity)
                        fatherZygosity1 <- as.character(
						        na.omit(fatherZygosity1)
                                )
                        if(length(fatherZygosity1) == 1){
                            fatherZygosity<-as.character(
							unique(fatherZygosity1)
							)
                        } 
                        else {
                            fatherZygosity <- paste(
								as.character(fatherZygosity1), 
								collapse = ","
								)
                        }
						} else {
						    fatherZygosity <- "-"
							}
                       
                        ###Mother Zygosity
                        if(length(dat2MomEqual$Zygosity) > 0){
						    if(length(as.character(unique
                            (dat2MomEqual$Zygosity))) == 1){
							    motherZygosity<-as.character(unique(dat2MomEqual$Zygosity))
							}else {
							    motherZygosity <- paste(
								as.character(dat2MomEqual$Zygosity), 
								collapse = ","
								)
							}
						} else if (length(dat2MomRange$Zygosity) > 0){
					    motherZygosity <- as.character(dat2MomRange$Zygosity)
                        motherZygosity1 <- gsub("-", NA, motherZygosity)
                        motherZygosity1 <- as.character(
						        na.omit(motherZygosity1)
                                )
                        if(length(motherZygosity1) == 1){
                            motherZygosity<-as.character(
							unique(motherZygosity1)
							)
                        } 
                        else {
                            motherZygosity <- paste(
								as.character(motherZygosity1), 
								collapse = ","
								)
                        }
						} else {
						    motherZygosity <- "-"
							}'
                       
                  
                        motherZygosity <- "-"
                        fatherZygosity <- "-"
                   
                        famno <- unname(hash::values(ha2, keys = sampID))
                        fampas <- paste("^",as.character(famindexfile$Tag),"$",sep="")
						famnopas <- paste("^",as.character(famno),"$",sep="")
                        famnumb <- length(grep(famnopas, fampas, fixed=TRUE))
                        data1 <- data.frame(dat1[nn, ], 
                            Internal_Freq_Perc_Filtered = 
                            as.numeric(format(round(
                            (filtFreq/(2*(usamp - famnumb))) * 100,digits=3)
                            , scientific = FALSE)),
                            Internal_Freq_Perc_Unfiltered = 
                            as.numeric (format(round(
							(UnfilteredFrequency/(2*(usamp - famnumb))) * 100,digits=3), 
                            scientific = FALSE)),
                            Internal_Homozygotes=as.numeric(cntHomozygotes), 
                            MotherZygosity = as.character(motherZygosity), 
                            FatherZygosity = as.character(fatherZygosity),  
                            stringsAsFactors = FALSE)
                   
                    datf <- rbind(datf, data1)
                    }
                    else if (nrow(dat2) == 1){ 
                        # dgv_match=TRUE Calculating percentage similarity
                        countfre <- 0;countfreunfilt<-0
                        svv <- dat2$nanoID
                        #size_internal <- dat2$Size
                        conf <- dat2$Confidence
						Found_in_self_molecules <- dat2$Found_in_self_molecules
						Fail_assembly_chimeric_score_DLE <- dat2$Fail_assembly_chimeric_score_DLE
                        BSPQI_status_DB <- as.character(dat2$Found_in_self_BSPQI_molecules)
                        BSSSI_status_DB <- as.character(dat2$Found_in_self_BSSSI_molecules)
                        BSPQI_chimeric_score_DB <- 
                        as.character(dat2$Fail_BSPQI_assembly_chimeric_score)
                        BSSSI_chimeric_score_DB <- 
                        as.character(dat2$Fail_BSSSI_assembly_chimeric_score)
                        perc_ref_query <- as.numeric(dat2$Size)/size1
					    perc_query_ref <- size1/as.numeric(dat2$Size)
                        #svv <- strsplit(as.character(svfam1), split = "_")
                        zygo <- as.character(dat2$Zygosity)
                        stt <- strsplit(as.character(svv), split = "[.]")
                        patID1 <- stt[[1]][2]
                        homozygo<-c()
                        svfamid1 <- stt[[1]][1]
                        motherZygosity <- ""
                        fatherZygosity <- ""
                        type <- as.character(dat2$Type)
                        typ1 <- as.character(dat2$Type1)
                        typ2 <- as.character(dat2$Type2)
                  #print(type)
                        'if (((dat2$RefStartPos >= rf_wb_int_parents[nn]) & 
						(dat2$RefEndPos <= re_fb_int_parents[nn])) & 
                        ((identical(type, "duplication") |
                        identical(type, "duplication_split")| 
                        identical(type, "duplication_inverted"))) & 
                        (identical(svfamid1, svfamid)) )
                        {
                            print(type)
                            print(svfamid1)
                            print(svfamid)
                            print(BSPQI_status_DB)
                            print(BSSSI_status_DB)
                            print(paste(ii,":",nn))
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
                        else'
						if (((identical(type, "duplication") |
                        identical(type, "duplication_split")| 
                        identical(type, "duplication_inverted"))) 
						& (perc_ref_query >= perc_similarity 
						& perc_query_ref >= perc_similarity)
						& !(identical(svfamid1, svfamid)) & 
                        ((BSPQI_status_DB == "yes" & 
                        BSSSI_status_DB == "yes") | 
                        (BSPQI_status_DB == "no" & 
                        BSSSI_status_DB == "yes") | 
                        (BSPQI_status_DB == "yes" & 
                        BSSSI_status_DB == "no") | 
                        (BSPQI_status_DB == "yes" & 
                        BSSSI_status_DB == "-") | 
                        (BSPQI_status_DB == "-" & 
                        BSSSI_status_DB == "yes")
						| (Found_in_DLE_self_molecules == "yes")) & 
                        ((BSPQI_chimeric_score_DB == "pass" & 
                        BSSSI_chimeric_score_DB == "pass") | 
                        (BSPQI_chimeric_score_DB == "fail" & 
                        BSSSI_chimeric_score_DB == "pass") | 
                        (BSPQI_chimeric_score_DB == "pass" & 
                        BSSSI_chimeric_score_DB == "fail") | 
                        (BSPQI_chimeric_score_DB == "pass" & 
                        BSSSI_chimeric_score_DB == "-") | 
                        (BSPQI_chimeric_score_DB == "-" & 
                        BSSSI_chimeric_score_DB == "pass")
						| (Fail_assembly_chimeric_score_DLE == "pass")) & 
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
                        (is.na(typ1) & typ2=="duplication_inverted")
						| (typ1 == "-" & typ2 == "-" & type == "duplication")
						| (typ1 == "-" & typ2 == "-" & type == "duplication_inverted")
						| (typ1 == "-" & typ2 == "-" & type == "duplication_split")
						| (identical(type, "duplication") |
                        identical(type, "duplication_split")| 
                        identical(type, "duplication_inverted")))){
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
                        !(identical(svfamid1, svfamid)) & 
                        (perc_ref_query >= perc_similarity & perc_query_ref >= perc_similarity) &
                        ((BSPQI_status_DB == "yes" &
                        BSSSI_status_DB == "yes") |
                        (BSPQI_status_DB == "no" & 
                        BSSSI_status_DB == "yes") |
                        (BSPQI_status_DB == "yes" & 
                        BSSSI_status_DB == "no") |
                        (BSPQI_status_DB == "yes" & 
                        BSSSI_status_DB == "-") |
                        (BSPQI_status_DB == "-" & 
                        BSSSI_status_DB == "yes")
						| (Found_in_DLE_self_molecules == "yes")) & 
                        ((size_internal>=limsize)) & 
                        (conf > indelconf) & 
                        ((typ1=="insertion" & typ2=="insertion") |
                        (typ1=="insertion" & is.na(typ2)) | 
                        (is.na(typ1) & typ2=="insertion")
						| (typ1 == "-" & typ2 == "-" & type == "insertion")
						|(type == "insertion"))){
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
                        identical(type, "insertion"))) 
						& (perc_ref_query >= perc_similarity 
						& perc_query_ref >= perc_similarity) 
						& !(identical(svfamid1, svfamid))){
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
                    Internal_Freq_Perc_Filtered = as.numeric(format(round((countfre/(2*(usamp - famno))) * 100,digits=3), scientific = FALSE)),
                    Internal_Freq_Perc_Unfiltered = as.numeric(format(round((countfreunfilt/(2*(usamp - famno))) * 100,digits=3), scientific = FALSE)),
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
                dat2 <- dat[which(((dat$RefStartPos <= rf[nn] & 
                    dat$RefStartPos >= rf_wb_ind[nn]) | 
                    (dat$RefStartPos >= rf[nn] & dat$RefStartPos <= 
                    rf_fb_ind[nn])) & ((dat$RefEndPos >= re[nn] & 
                    dat$RefEndPos <= re_wf_ind[nn]) |
                    (dat$RefEndPos <= re[nn] & dat$RefEndPos >= re_wb_ind[nn]))
                    ), ]
                size1 <- size_bn[nn]
                
                ## Writing if the dgv_match is TRUE or not
                
                if (nrow(dat2) > 1){
				    #dat2$perc <- as.numeric(dat2$Size)/size1
					stt12 <- str_split(as.character(dat2$nanoID), pattern = "\\.")
					famid <- c(); relnid <- c()
					dat2$FamilyID <- sapply(stt12, function(x) x[1])
					dat2$RelationID <- sapply(stt12, function(x) x[2])
					'dat2MomEqual <- subset(dat2,(((RefStartPos ==rf[nn]) & (RefEndPos==re[nn])) 
                            & (RefcontigID2 == chromo2[nn]) 
							& (as.character(Type) == variantType2[nn])
                            & (as.character(FamilyID) == svfamid)
							& (patID == 1 & RelationID == 2)))
					dat2DadEqual <- subset(dat2,(((RefStartPos ==rf[nn]) & (RefEndPos==re[nn]))           
							& (RefcontigID2 == chromo2[nn]) 
							& (as.character(Type) == variantType2[nn])
                            & (as.character(FamilyID) == svfamid)
							& (patID == 1 & RelationID == 3)))
					dat2MomRange <- subset(dat2,((( RefStartPos >= rf_wb_int_parents[nn]) & (RefEndPos <= re_fb_int_parents[nn]))  
						& (RefcontigID2 == chromo2[nn]) 
						& (as.character(Type) == variantType2[nn])
                            & (as.character(FamilyID) == svfamid)
							& (patID == 1 & RelationID == 2)))
					dat2DadRange <- subset(dat2,((( RefStartPos >= rf_wb_int_parents[nn]) & (RefEndPos <= re_fb_int_parents[nn])) 
						& (RefcontigID2 == chromo2[nn]) 
						& (as.character(Type) == variantType2[nn])
                            & (as.character(FamilyID) == svfamid)
							& (patID == 1 & RelationID == 3)))'
					dat2Frequencyinvtrans <- subset(dat2,(((
                        as.character(Type) == "inversion" &
                          as.numeric(Confidence) >= invconf)|
                           ((as.character(Type) == "translocation_intrachr"
                            | as.character(Type) == "translocation_interchr"
                            | as.character(Type) == "translocation")
 							& as.numeric(Confidence) >= transconf))
						& (as.character(FamilyID) != svfamid)
						& (RefcontigID2 == chromo2[nn]) 
						& ((Found_in_self_BSPQI_molecules == "yes" & 
                            Found_in_self_BSSSI_molecules == "yes") 
                            | (Found_in_self_BSPQI_molecules == "no" 
							& Found_in_self_BSSSI_molecules == "yes") 
							| (Found_in_self_BSPQI_molecules == "yes" 
							& Found_in_self_BSSSI_molecules == "no") 
							| (Found_in_self_BSPQI_molecules == "yes" 
							& Found_in_self_BSSSI_molecules == "-") 
							| (Found_in_self_BSPQI_molecules == "-" 
							& Found_in_self_BSSSI_molecules == "yes")
							| (Found_in_DLE_self_molecules == "yes"))
							& ((Fail_BSPQI_assembly_chimeric_score == "pass" & 
                            Fail_BSSSI_assembly_chimeric_score == "pass") 
                            | (Fail_BSPQI_assembly_chimeric_score == "fail" 
							& Fail_BSSSI_assembly_chimeric_score == "pass") 
							| (Fail_BSPQI_assembly_chimeric_score == "pass" 
							& Fail_BSSSI_assembly_chimeric_score == "fail") 
							| (Fail_BSPQI_assembly_chimeric_score == "pass" 
							& Fail_BSSSI_assembly_chimeric_score == "-") 
							| (Fail_BSPQI_assembly_chimeric_score == "-" 
							& Fail_BSSSI_assembly_chimeric_score == "pass")
							| (Fail_assembly_chimeric_score_DLE == "pass"))))
						dat2UnfilteredFrequency <- subset(dat2,(
						    (RefcontigID2 == chromo2[nn]) 
					        & ((as.character(Type) == variantType2[nn])
							& (as.character(FamilyID) != svfamid))))
					dat2Frequency <- dat2Frequencyinvtrans
                  ##Filtered
                    if (nrow(dat2Frequency) > 0){
					    print(0)
					    countfre1 <- 0
						dat2Homozygotes <- data.frame(dat2Frequency %>% 
						mutate(Homozygo = Zygosity=="homozygous") %>%
						group_by(nanoID) %>%
						summarise (Homozygotes = sum(Homozygo)))
						cntHomozygotes <- sum(dat2Homozygotes$Homozygo)
						dat2Homozygotes [,2] <- ifelse(dat2Homozygotes[,2] > 0, 2, 0)
						dat2Heterozygotes <- data.frame(dat2Frequency %>% 
						mutate(Heterozygo = Zygosity == "heterozygous") %>%
						group_by(nanoID) %>%
						summarise (Heterozygotes = sum(Heterozygo)))
						dat2Heterozygotes [,2] <- ifelse(dat2Heterozygotes[,2] > 0, 1, 0)
						dat2Unknown <- data.frame(dat2Frequency %>% 
						mutate(Unknown = Zygosity == "unknown") %>%
						group_by(nanoID) %>%
						summarise (Unknown = sum(Unknown)))
						dat2Unknown [,2] <- ifelse(dat2Unknown[,2] > 0, 2, 0)
						dat2filtFreq <- data.frame(NR = dat2Homozygotes[,1], Homozygote = dat2Homozygotes[,2],
						Heterozygotes = dat2Heterozygotes[,2], Unknown = dat2Unknown[,2])
						row.names(dat2filtFreq) <- as.character(dat2Homozygotes[,1])
						dat2filtFreq$RowSum <- ifelse((dat2filtFreq$Heterozygotes ==1 & (dat2filtFreq$Homozygote == 0 & dat2filtFreq$Unknown == 0)), 1, 2)
						filtFreq <- sum(dat2filtFreq$RowSum)
						} else{
						    filtFreq  = 0
							cntHomozygotes = 0
                        }
						
					
					###Unfiltered
					
					if (nrow(dat2UnfilteredFrequency) > 0){
					    dat2Homozygotes <- data.frame(dat2UnfilteredFrequency %>% 
						mutate(Homozygo=Zygosity=="homozygous") %>%
						group_by(nanoID) %>%
						summarise (Homozygotes = sum(Homozygo)))
						#cntHomozygotes <- sum(dat2Homozygotes$Homozygotes)
						dat2Homozygotes [,2] <- ifelse(dat2Homozygotes[,2] > 0, 2, 0)
						dat2Heterozygotes <- data.frame(dat2UnfilteredFrequency %>% 
						mutate(Heterozygo = Zygosity == "heterozygous") %>%
						group_by(nanoID) %>%
						summarise (Heterozygotes = sum(Heterozygo)))
						dat2Heterozygotes [,2] <- ifelse(dat2Heterozygotes[,2] > 0, 1, 0)
						dat2Unknown <- data.frame(dat2UnfilteredFrequency %>% 
						mutate(Unknown = Zygosity == "unknown") %>%
						group_by(nanoID) %>%
						summarise (Unknown = sum(Unknown)))
						dat2Unknown [,2] <- ifelse(dat2Unknown[, 2] > 0, 2, 0)
						dat2unfiltFreq <- data.frame(Homozygote = dat2Homozygotes[,2],
						Heterozygotes = dat2Heterozygotes[,2], Unknown = dat2Unknown[,2])
						row.names(dat2unfiltFreq) <- as.character(dat2Homozygotes[,1])
						dat2unfiltFreq$RowSum <- ifelse((dat2unfiltFreq$Heterozygotes ==1 
						& (dat2unfiltFreq$Homozygote == 0 
						& dat2unfiltFreq$Unknown == 0)), 1, 2)
						UnfilteredFrequency  = sum(dat2unfiltFreq$RowSum)
						} else{
						    UnfilteredFrequency = 0
						}
						
                 
                    ##Calculating filtered Frequency INDEL
                    ##Zygosity Parent
					    'if(length(dat2DadEqual$Zygosity) > 0){
						    if(length(as.character(unique
                            (dat2DadEqual$Zygosity))) == 1){
							    fatherZygosity<-as.character(unique(dat2DadEqual$Zygosity))
							}else {
							    fatherZygosity <- paste(
								as.character(dat2DadEqual$Zygosity), 
								collapse = ","
								)
							}
						} else if (length(dat2DadRange$Zygosity) > 0){
					    fatherzygosity <- as.character(dat2DadRange$Zygosity)
                        fatherZygosity1 <- gsub("-", NA, fatherzygosity)
                        fatherZygosity1 <- as.character(
						        na.omit(fatherZygosity1)
                                )
                        if(length(fatherZygosity1) == 1){
                            fatherZygosity<-as.character(
							unique(fatherZygosity1)
							)
                        } 
                        else {
                            fatherZygosity <- paste(
								as.character(fatherZygosity1), 
								collapse = ","
								)
                        }
						} else {
						    fatherZygosity <- "-"
							}
                       
                        ###Mother Zygosity
                        if(length(dat2MomEqual$Zygosity) > 0){
						    if(length(as.character(unique
                            (dat2MomEqual$Zygosity))) == 1){
							    motherZygosity<-as.character(unique(dat2MomEqual$Zygosity))
							}else {
							    motherZygosity <- paste(
								as.character(dat2MomEqual$Zygosity), 
								collapse = ","
								)
							}
						} else if (length(dat2MomRange$Zygosity) > 0){
					    motherZygosity <- as.character(dat2MomRange$Zygosity)
                        motherZygosity1 <- gsub("-", NA, motherZygosity)
                        motherZygosity1 <- as.character(
						        na.omit(motherZygosity1)
                                )
                        if(length(motherZygosity1) == 1){
                            motherZygosity<-as.character(
							unique(motherZygosity1)
							)
                        } 
                        else {
                            motherZygosity <- paste(
								as.character(motherZygosity1), 
								collapse = ","
								)
                        }
						} else {
						    motherZygosity <- "-"
							}'
                       
                  
                  
                        motherZygosity <- "-"
                        fatherZygosity <- "-"
                        famno <- unname(hash::values(ha2, keys = sampID))
                        fampas <- paste("^",as.character(famindexfile$Tag),"$",sep="")
						famnopas <- paste("^",as.character(famno),"$",sep="")
                        famnumb <- length(grep(famnopas, fampas, fixed=TRUE))
                        data1 <- data.frame(dat1[nn, ], 
                            Internal_Freq_Perc_Filtered = 
                            as.numeric(format(round(
                            (filtFreq/(2*(usamp - famnumb))) * 100,digits=3)
                            , scientific = FALSE)),
                            Internal_Freq_Perc_Unfiltered = 
                            as.numeric
                            (format(round((UnfilteredFrequency/(2*(usamp - famnumb))) * 100,digits=3), scientific = FALSE)),
                            Internal_Homozygotes=as.numeric(cntHomozygotes), 
                            MotherZygosity = as.character(motherZygosity), 
                            FatherZygosity = as.character(fatherZygosity),  
                            stringsAsFactors = FALSE)
						datf <- rbind(datf, data1)
                } 
                else if (nrow(dat2) == 1)
                {
                    ## dgv_match=TRUE Calculating percentage similarity
                    countfre<-0;countfreunfilt<-0
                    chrom2<-dat2$RefcontigID2
                    svv <- dat2$nanoID
                    zygo <- as.character(dat2$Zygosity)
                    stt <- strsplit(as.character(svv), split = "[.]")
                    type <- as.character(dat2$Type)
                    typ1 <- as.character(dat2$Type1)
                    typ2 <- as.character(dat2$Type2)
                    patID1 <- stt[[1]][2]
                    svfamid1 <- stt[[1]][1]
                    conf <- dat2$Confidence
                    homozygo<-c()
					Fail_assembly_chimeric_score_DLE <- dat2$Fail_assembly_chimeric_score_DLE
					Found_in_DLE_self_molecules <- dat2$Found_in_DLE_self_molecules
                    BSPQI_status_DB <- as.character(dat2$Found_in_self_BSPQI_molecules)
                    BSSSI_status_DB <- as.character(dat2$Found_in_self_BSSSI_molecules)
                    BSPQI_chimeric_score_DB <- as.character(
                    dat2$Fail_BSPQI_assembly_chimeric_score)
                    BSSSI_chimeric_score_DB <- as.character(
                    dat2$Fail_BSSSI_assembly_chimeric_score)
                    'if (((dat2$RefStartPos >= rf_wb_int_parents[nn]) & 
					    (dat2$RefEndPos <= re_fb_int_parents[nn])) & 
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
                    else'
					if (identical(chrom2,chromo2[nn]) &
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
                    BSSSI_status_DB == "yes")
					| (Found_in_DLE_self_molecules == "yes"))
					|((identical(type, "inversion") & conf > invconf) | 
                    ((identical(type, "translocation_intrachr")
                      |identical(type, "translocation_interchr")
                      | identical(type, "translocation"))
 					  & (conf > transconf))) & 
                    ((BSPQI_chimeric_score_DB == "pass" & 
                    BSSSI_chimeric_score_DB == "pass") | 
                    (BSPQI_chimeric_score_DB == "fail" & 
                    BSSSI_chimeric_score_DB == "pass") | 
                    (BSPQI_chimeric_score_DB == "pass" & 
                    BSSSI_chimeric_score_DB == "fail") |
                    (BSPQI_chimeric_score_DB == "pass" & 
                    BSSSI_chimeric_score_DB == "-") | 
                    (BSPQI_chimeric_score_DB == "-" & 
                    BSSSI_chimeric_score_DB == "pass")
					|(Fail_assembly_chimeric_score_DLE == "pass")) ){
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
                        identical(type, variantType2[nn]) & 
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
                    famno <- as.numeric(unname(hash::values(ha, keys = svfamid)))
                    data1 <- data.frame(dat1[nn, ], 
                    Internal_Freq_Perc_Filtered = 
                    as.numeric(format(round((
					countfre/(2*(usamp - famno))) * 100,digits=3), scientific = FALSE)),
                    Internal_Freq_Perc_Unfiltered = as.numeric(format(round(
					(countfreunfilt/(2*(usamp - famno))) * 100,digits=3), scientific = FALSE)), 
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
