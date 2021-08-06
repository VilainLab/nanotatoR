#' Merging dual labelled smaps
#'
#' @param path character. Path to the solo files directory.
#' @param pattern  character. Patternna for the solo files.
#' @param outMode  character. The ouput mode. Choices, dataframe or Text.
#' @param outpath character. Path where the dual labelled merged samples are kept.
#'            Is mandatory if outMode is Text.
#' @return Text files containg merged smaps from different samples
#' @examples
#' a <- mergingSMAP_SVMerge(
#' path = system.file("extdata", "SoloFile/", package="nanotatoR"),
#' pattern = "*.txt", outMode = "dataframe", 
#' outpath = system.file("extdata", "SoloFile/", package="nanotatoR"))
#' @export
mergingSMAP_SVMerge <- function(path ,
    pattern , outMode = c("dataframe","Text"), outpath)
{
# Selecting all files with given pattern
l <- list.files(path = path, pattern = pattern)
nam <- c()
datfinal <- data.frame()
###Merging
# Reading all the smap files individually and merging them into a dataframe
for (ii in seq_along((l))){
    print(l[ii])
    con <- file(file.path(path, l[ii]), "r")
    r10 <- readLines(con, n = -1)
    close(con)
    gg1 <- grep("# BSPQI Sample", r10)
    stt <- strsplit(r10[gg1], split = ":")
    fname_temp <- stt[[1]][2]
    stt1 <- strsplit(fname_temp, split = "_BspQI_assembly*")
    fname_temp <- stt1[[1]][1]
    stt2 <- strsplit(fname_temp, split = "_BssSI_assembly*")
    fname<- stt2[[1]][1]
    fname <- str_squish(fname)
    #gg7 <- grep("UDN*", fname_temp1)
    'if (length(gg7) > 0){
    stt3 <- strsplit(fname_temp1, split ="_")
    fname <- stt3[[1]][1]
    fname <- str_squish(fname)
    } else{
        fname <- str_squish(fname_temp1)
    }'
    g1 <- grep("RefEndPos", r10)
    g2 <- grep("RefStartPos", r10)
    
    # r10<-as.character(r10)
    if (g1 == g2)
    {
    #g3 <- grep("#h ",r10)
    dat <- gsub("# ", "", as.character(r10))
    
    # dat<-gsub('\t',' ',as.character(dat))
    dat4 <- textConnection(dat[g1:length(dat)])
    ##### print(dim(dat4))
    r1 <- read.table(dat4, sep = "\t", header = TRUE)
    close(dat4)
    nam <- names(r1)
    r1 <- data.frame(r1)
    g5 <- grep("PutativeGeneFusion", nam)
    if (length(g5) > 0){
        r1 <- r1[ , 1:(g5-1)]
    }
    else{
        r1 <- r1
    }
    'st <- strsplit(l[ii], split = ".txt")
    fname <- st[[1]][1]'
    
    fname1 <- c(rep(fname, length(r1[, 1])))
    r1 <- cbind(fname1, r1)
    #### print(dim(r1))
    datfinal <- data.frame(rbind(datfinal, r1))
    
    #### print(dim(datfinal))
    } else
    {
    stop("column names doesnot Match")
    }
    
}
        Found_in_SE_molecules <- as.character(
            c(rep ("-", times = nrow(datfinal)))
            )
        Fail_assembly_chimeric_score_SE <- as.character(
            c(rep ("-", times = nrow(datfinal)))
            )
        Method <- as.character(
            c(rep ("BSPQI_BSSSI", times = nrow(datfinal)))
        )
        datFinal_Dual <- data.frame(
        SVIndex = as.character(datfinal$SVIndex), 
        SampleID = as.character(datfinal$fname1), 
        RefcontigID1 = as.numeric(datfinal$RefcontigID1), 
        RefcontigID2 = as.numeric(datfinal$RefcontigID2), 
        RefStartPos = as.numeric(datfinal$RefStartPos), 
        RefEndPos = as.numeric(datfinal$RefEndPos), 
        Confidence = as.numeric(datfinal$Confidence), 
        Size = as.numeric(datfinal$Size), 
        Zygosity = as.character(datfinal$Zygosity), 
        Type = as.character(datfinal$Type), 
        Type1 = as.character(datfinal$Type1), 
        Type2 = as.character(datfinal$Type2), 
        Fail_BSPQI_assembly_chimeric_score = as.character(
        datfinal$Fail_BSPQI_assembly_chimeric_score
        ), 
        Fail_BSSSI_assembly_chimeric_score = as.character(
        datfinal$Fail_BSSSI_assembly_chimeric_score
        ),
        Fail_assembly_chimeric_score_SE = Fail_assembly_chimeric_score_SE, 
        Found_in_self_BSPQI_molecules = as.character(
        datfinal$Found_in_self_BSPQI_molecules
        ), 
        Found_in_self_BSSSI_molecules = as.character(
        datfinal$Found_in_self_BSSSI_molecules
        ), 
        Found_in_SE_self_molecules = Found_in_SE_molecules , 
        Method = Method
        )

        if (outMode == "Text"){
            len <- length(l)
            fileOut <- paste0("Dual_1_",len,".txt")
            write.table(datFinal_Dual, file.path(
                outpath, fileOut), 
                row.names = FALSE)
        }
        else if (outMode == "dataframe"){
            return(datFinal_Dual)
        }
        else {
            stop ("Your Outmode is absent !!!")
        }
}
#' Merging DLE labelled smaps
#'
#' @param path character. Path to the solo files directory.
#' @param pattern    character. Pattern for the solo files.
#' @param outMode    character. The ouput mode. Choices, dataframe or Text.
#' @param outpath character. Path where the dual labelled 
#'        merged samples are kept. Is mandatory if outMode is Text.
#' @return Text files containg merged smaps from different samples
#' @examples
#' mergedSmap <- mergingSMAP_SE (
#' path = system.file("extdata", "SoloFile/", package="nanotatoR"),
#' pattern = "*_DLE1_*", outMode = "dataframe", 
#' outpath = system.file("extdata", "Merged/", package="nanotatoR"))
#' @import hash
#' @import tidyverse
#' @export

mergingSMAP_SE <- function(
path ,pattern,outMode = c("Text", "dataframe"),
outpath)
{
l <- list.files(path = path, pattern = pattern, full.names = TRUE)
nam <- c()
datfinal <- data.frame()
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
		'if(pipeline == "RVP"){
		    Zygosity = as.character(c(rep ("-", times = nrow(datfinal))))
		}
		else{
		    Zygosity = as.character(datfinal$Zygosity)
		}	
        strsample <- strsplit(as.character(datfinal$Sample), 
            split = "*_DLE_*")
        SampleID = sapply(strsample, function(x) x[1])'
        datFinal_DLE <- data.frame(
            SVIndex = as.character(datfinal$SmapEntryID), 
            SampleID = datfinal$SampleID,
            RefcontigID1 = as.numeric(datfinal$RefcontigID1), 
            RefcontigID2 = as.numeric(datfinal$RefcontigID2),
            RefStartPos = as.numeric(datfinal$RefStartPos), 
            RefEndPos = as.numeric(datfinal$RefEndPos), 
            Confidence = as.numeric(datfinal$Confidence), 
            Size = as.numeric(datfinal$Size), 
            Zygosity = Zygosity, 
            Type = as.character(datfinal$Type),
            Type1 = Type1, Type2 = Type2, 
            Fail_BSPQI_assembly_chimeric_score = Fail_BSPQI_assembly_chimeric_score, 
            Fail_BSSSI_assembly_chimeric_score =Fail_BSSSI_assembly_chimeric_score,
            Fail_assembly_chimeric_score_SE = as.character(
                datfinal$Fail_assembly_chimeric_score
            ),
            Found_in_self_BSPQI_molecules = Found_in_self_BSPQI_molecules, 
            Found_in_self_BSSSI_molecules = Found_in_self_BSPQI_molecules,
            Found_in_SE_self_molecules = as.character(
                datfinal$Found_in_self_molecules
                ), 
            Method = Method
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
#' Mapping Realtionship to unique nanoIDs
#'
#' @param Samplecodes character. File containing relations and IDs associated to them.
#' @param mergeKey character. File containing sample ID and relation.
#' @param outMode    character. The ouput mode. Choices, dataframe or Text.
#' @param outpath character. Path where the dual labelled merged samples are kept.
#'                        Is mandatory if outMode is Text.
#' @return Text files containg merged smaps from different samples
#' @examples
#' FamilyInfoPrep(
#' Samplecodes = system.file("extdata", "nanotatoR_sample_codes.csv", package="nanotatoR"),
#' mergeKey = system.file("extdata", "nanotatoR_control_sample_codes.csv", package="nanotatoR"),
#' outMode = c("dataframe"))
#' @import hash
#' @import tidyverse
#' @export

FamilyInfoPrep <- function(
    Samplecodes = "X:/Hayks_Materials/BNG/Projects/nanotatoR_sample_codes.csv",
    mergeKey = "X:/Hayks_Materials/BNG/Projects/MergeKey.csv",
    outMode = c("Text", "dataframe"),
    outpath = "X:/Hayks_Materials/BNG/Projects/VAP_DLE1_solo_SMAPs/Merged"){
                    
        samplecodes <- read.csv(Samplecodes)
        mergecodes <- read.csv(mergeKey)
        usampIDs <- unique(mergecodes$SampleID)
        famtag <- unique(mergecodes$Tag)
        famtagp <- paste0("^",as.character(mergecodes$Tag),"$")
        #tag <- c()
        nid <- c()
        datFinal <- data.frame()
    for(ki in seq_along(famtag)){
        ftag <- gsub("^\\w","NR",famtag[ki])
        dat_Proband <- mergecodes[which(
        mergecodes$Tag == famtag[ki] & mergecodes$Relationship == "proband"
        ),]
        codeProband <- samplecodes[which(samplecodes$Relation == "proband"),]
        codeProbandID <- unique(codeProband$nanoID)
        if (nrow(dat_Proband) >= 1){
        
        dat_Proband$NID <- paste0(ftag, ".", codeProbandID)
        datFinal <- rbind(datFinal,dat_Proband)
    } else{
        print("Proband Absent")
    }
    dat_Mother <- mergecodes[which(mergecodes$Tag == famtag[ki] 
        & mergecodes$Relationship == "mother"),]
    codeMother <- samplecodes[which(samplecodes$Relation == "mother"),]
    codeMotherID <- unique(codeMother$nanoID)
    if (nrow(dat_Mother) >= 1){
        dat_Mother$NID <- paste0(ftag,".", codeMotherID)
        datFinal <- rbind(datFinal,dat_Mother)
    } 
    else{
        print("Mother Absent")
    }
    dat_Father <- mergecodes[which(mergecodes$Tag == famtag[ki] &
            mergecodes$Relationship == "father"),]
    codeFather <- samplecodes[which(samplecodes$Relation == "father"),]
    codeFatherID <- unique(codeFather$nanoID)
    
    if (nrow(dat_Father) >= 1){
        dat_Father$NID <- paste0(ftag,".", codeFatherID)
        datFinal <- rbind(datFinal,dat_Father)
    } else{
        print("Father Absent")
    }
    
    dat_Sister <- mergecodes[which(mergecodes$Tag == famtag[ki] &
                mergecodes$Relationship == "sister"),]
    codeSister <- samplecodes[which(samplecodes$Relation == "sister"),]
    codeSisterID <- unique(codeSister$nanoID)
    if (nrow(dat_Sister) > 1){
        for (i    in 1:nrow(dat_Sister)){
        dat_Sister$NID[i] <- paste0(ftag,".",codeSisterID,".",i)
    }
    datFinal <- rbind(datFinal,dat_Sister)
    }else if (nrow(dat_Sister) == 1){
        dat_Sister$NID <- paste0(ftag, ".", codeSisterID)
        datFinal <- rbind(datFinal,dat_Sister)
    }else {
        print ("Sister not present")
    }
    
    dat_Brother <- mergecodes[which(mergecodes$Tag == famtag[ki] &
                mergecodes$Relationship == "brother"),]
    codebrother <- samplecodes[which(samplecodes$Relation == "brother"),]
    codebrotherID <- unique(codebrother$nanoID)
    if (nrow(dat_Brother) > 1){
        for (i    in 1:nrow(dat_Brother)){
            dat_Brother$NID[i] <- paste0(ftag,".",codebrotherID,".",i)
        }
    datFinal <- rbind(datFinal,dat_Brother)
    } else if (nrow(dat_Brother) == 1){
            dat_Brother$NID <- paste0(ftag,".",codebrotherID)
            datFinal <- rbind(datFinal,dat_Brother)
    }else {
        print ("Brother not present")
    }
    
    dat_Sibling <- mergecodes[which(mergecodes$Tag == famtag[ki] &
                    mergecodes$Relationship == "sibling"),]
    codesibling <- samplecodes[which(samplecodes$Relation == "sibling"),]
    codesiblingID <- unique(codesibling$nanoID)
    if (nrow(dat_Sibling) > 1){
        for (i    in 1:nrow(dat_Sibling)){
            dat_Sibling$NID[i] <- paste0(ftag, ".", codesiblingID,".",i)
        }
        datFinal <- rbind(datFinal,dat_Sibling)
    }else if (nrow(dat_Sibling) == 1){
        dat_Sibling$NID <- paste0(ftag, ".", codesiblingID)
        datFinal <- rbind(datFinal,dat_Sibling)
    }else {
        print ("Sibling not present")
    }
    
    dat_Cousin <- mergecodes[which(mergecodes$Tag == famtag[ki] &
                    mergecodes$Relationship == "cousin"),]
    codecousin <- samplecodes[which(samplecodes$Relation == "cousin"),]
    codecousinID <- unique(codecousin$nanoID)
    if (nrow(dat_Cousin) > 1){
        for (i    in 1:nrow(dat_Cousin)){
            dat_Cousin$NID[i] <- paste0(ftag, ".",codecousinID, ".", i)
    }
    datFinal <- rbind(datFinal,dat_Cousin)
    }else if (nrow(dat_Cousin) == 1){
        dat_Cousin$NID <- paste0(ftag, ".", codecousinID)
        datFinal <- rbind(datFinal,dat_Cousin)
    }else {
    print ("Cousin not present")
}
    
    dat_Husband <- mergecodes[which(mergecodes$Tag == famtag[ki] &
        mergecodes$Relationship == "husband"),
        ]
    codehusband <- samplecodes[which(samplecodes$Relation == "husband"),]
    codehusbandID <- unique(codehusband$nanoID)
    if (nrow(dat_Husband) > 1){
        for (i    in 1:nrow(dat_Husband)){
            dat_Husband$NID[i] <- paste0(ftag, ".",codehusbandID, ".", i)
        }
    datFinal <- rbind(datFinal,dat_Husband)
    } else if (nrow(dat_Husband) == 1){
        dat_Husband$NID <- paste0(ftag, ".", codehusbandID)
        datFinal <- rbind(datFinal,dat_Husband)
    } else {
        print ("Husband not present")
    }
    
    dat_Wife <- mergecodes[which(mergecodes$Tag == famtag[ki] &
                                    mergecodes$Relationship == "wife"),]
    codewife <- samplecodes[which(samplecodes$Relation == "wife"),]
    codewifeID <- unique(codewife$nanoID)
    if (nrow(dat_Wife) > 1){
        for (i    in 1:nrow(dat_Wife)){
            dat_Wife$NID[i] <- paste0(ftag, ".",codewifeID, ".", i)
        }
    datFinal <- rbind(datFinal,dat_Wife)
    } else if (nrow(dat_Wife) == 1){
        dat_Wife$NID <- paste0(ftag, ".",codewifeID)
        datFinal <- rbind(datFinal,dat_Wife)
    } else {
        print ("Wife not present")
    }
    
    dat_Son <- mergecodes[which(mergecodes$Tag == famtag[ki] &
                    mergecodes$Relationship == "son"),]
    codeson <- samplecodes[which(samplecodes$Relation == "son"),]
    codesonID <- unique(codeson$nanoID)
    if (nrow(dat_Son) > 1){
        for (i    in 1:nrow(dat_Son)){
        dat_Son$NID[i] <- paste0(ftag, ".",codesonID, ".", i)
    }
    datFinal <- rbind(datFinal,dat_Son)
    } else if (nrow(dat_Son) == 1){
        dat_Son$NID <- paste0(ftag, ".",codesonID)
        datFinal <- rbind(datFinal,dat_Son)
    }else {
        print ("Son not present")
    }
    
    dat_Daughter <- mergecodes[which(mergecodes$Tag == famtag[ki] &
                    mergecodes$Relationship == "daughter"),]
    codedaughter <- samplecodes[which(samplecodes$Relation == "daughter"),]
    codedaughterID <- unique(codedaughter$nanoID)
    if (nrow(dat_Daughter) > 1){
        for (i    in 1:nrow(dat_Daughter)){
        dat_Daughter$NID[i] <- paste0(ftag, ".",codedaughterID, ".", i)
    }
    datFinal <- rbind(datFinal,dat_Daughter)
    } else if (nrow(dat_Daughter) == 1){
        dat_Daughter$NID <- paste0(ftag, ".",codedaughterID)
        datFinal <- rbind(datFinal,dat_Daughter)
    }else {
        print ("Daughter not present")
    }
    
    
    dat_maternalGrandmother <- mergecodes[which(mergecodes$Tag == famtag[ki] &
                            mergecodes$Relationship == "maternalgrandmother"),]
    codematernalGrandmother <- samplecodes[which(samplecodes$Relation == "maternalgrandmother"),]
    codematernalGrandmotherID <- unique(codematernalGrandmother$nanoID)
    if (nrow(dat_maternalGrandmother) > 1){
        for (i    in 1:nrow(dat_maternalGrandmother)){
        dat_maternalGrandmother$NID[i] <- paste0(ftag, ".", codematernalGrandmotherID, ".", i)
    }
        datFinal <- rbind(datFinal,dat_maternalGrandmother)
    }else if (nrow(dat_maternalGrandmother) == 1){
        dat_maternalGrandmother$NID <- paste0(ftag, ".", codematernalGrandmotherID)
        datFinal <- rbind(datFinal,dat_maternalGrandmother)
    }else {
        print ("Maternal Grand Mother not present")
    }
    
    dat_maternalgrandfather <- mergecodes[which(mergecodes$Tag == famtag[ki] &
                        mergecodes$Relationship == "maternalgrandfather"),]
    codegrandmaternalgrandfather <- samplecodes[which(samplecodes$Relation == "maternalgrandfather"),]
    codegrandmaternalgrandfatherID <- unique(codegrandmaternalgrandfather$nanoID)
    if (nrow(dat_maternalgrandfather) > 1){
        for (i    in 1:nrow(dat_maternalgrandfather)){
            dat_maternalgrandfather$NID[i] <- paste0(ftag, ".", codegrandmaternalgrandfatherID, ".", i)
        }
    datFinal <- rbind(datFinal,dat_maternalgrandfather)
    } else if (nrow(dat_maternalgrandfather) == 1){
        dat_maternalgrandfather$NID <- paste0(ftag, ".", codegrandmaternalgrandfatherID)
        datFinal <- rbind(datFinal,dat_maternalgrandfather)
    } else {
        print ("Maternal Grand Father not present")
    }
    dat_paternalgrandmother <- mergecodes[which(mergecodes$Tag == famtag[ki] &
                            mergecodes$Relationship == "paternalgrandmother"),]
    codegrandpaternalgrandmother <- samplecodes[which(samplecodes$Relation == "paternalgrandmother"),]
    codegrandpaternalgrandmotherID <- unique(codegrandpaternalgrandmother$nanoID)
    if (nrow(dat_paternalgrandmother) > 1){
        for (i    in 1:nrow(dat_paternalgrandmother)){
        dat_paternalgrandmother$NID[i] <- paste0(ftag, ".", codegrandpaternalgrandmotherID, ".", i)
    }
        datFinal <- rbind(datFinal,dat_paternalgrandmother)
    }else if (nrow(dat_paternalgrandmother) == 1){
        dat_paternalgrandmother$NID <- paste0(ftag, ".", codegrandpaternalgrandmotherID)
        datFinal <- rbind(datFinal,dat_paternalgrandmother)
    }else {
        print ("Paternal Grand Mother not present")
    }
    dat_paternalgrandfather <- mergecodes[which(mergecodes$Tag == famtag[ki] &
                        mergecodes$Relationship == "paternalgrandfather"),]
    codepaternalgrandfather <- samplecodes[which(samplecodes$Relation == "paternalgrandfather"),]
    codepaternalgrandfatherID <- unique(codepaternalgrandfather$nanoID)
    if (nrow(dat_paternalgrandfather) > 1){
        for (i    in 1:nrow(dat_paternalgrandfather)){
            dat_paternalgrandfather$NID[i] <- paste0(ftag, ".", codepaternalgrandfatherID, ".", i)
        }
    datFinal <- rbind(datFinal,dat_paternalgrandfather)
    } else if (nrow(dat_paternalgrandfather) == 1){
        dat_paternalgrandfather$NID <- paste0(ftag, ".", codepaternalgrandfatherID)
        datFinal <- rbind(datFinal,dat_paternalgrandfather)
    } else {
        print ("Paternal Grand Father not present")
    }
    
    dat_maternaluncle <- mergecodes[which(mergecodes$Tag == famtag[ki] &
                        mergecodes$Relationship == "maternaluncle"),]
    codematernaluncle <- samplecodes[which(samplecodes$Relation == "maternaluncle"),]
    codematernaluncleID <- unique(codematernaluncle$nanoID)
    if (nrow(dat_maternaluncle) > 1){
        for (i    in 1:nrow(dat_maternaluncle)){
            dat_maternaluncle$NID[i] <- paste0(ftag, ".", codematernaluncleID, ".", i)
        }
    datFinal <- rbind(datFinal,dat_maternaluncle)
    } else if (nrow(dat_maternaluncle) == 1){
        dat_maternaluncle$NID <- paste0(ftag, ".", codematernaluncleID)
        datFinal <- rbind(datFinal,dat_maternaluncle)
    } else {
        print ("Maternal Uncle not present")
    }
    
    dat_maternalaunt <- mergecodes[which(mergecodes$Tag == famtag[ki] &
                        mergecodes$Relationship == "maternalaunt"),]
    codematernalaunt <- samplecodes[which(samplecodes$Relation == "maternalaunt"),]
    codematernalauntID <- unique(codematernalaunt$nanoID)
    if (nrow(dat_maternalaunt) > 1){
        for (i    in 1:nrow(dat_maternalaunt)){
            dat_maternalaunt$NID[i] <- paste0(ftag, ".", codematernalauntID, ".", i)
        }
    datFinal <- rbind(datFinal,dat_maternalaunt)
    } else if (nrow(dat_maternalaunt) == 1){
        dat_maternalaunt$NID <- paste0(ftag, ".", codematernalauntID)
        datFinal <- rbind(datFinal,dat_maternalaunt)
    } else {
        print ("Maternal Aunt not present")
    }
    
    dat_paternaluncle <- mergecodes[which(mergecodes$Tag == famtag[ki] &
                        mergecodes$Relationship == "paternaluncle"),]
    codepaternaluncle <- samplecodes[which(samplecodes$Relation == "paternaluncle"),]
    codepaternaluncleID <- unique(codepaternaluncle$nanoID)
    if (nrow(dat_paternaluncle) > 1){
        for (i    in 1:nrow(dat_paternaluncle)){
            dat_paternaluncle$NID[i] <- paste0(ftag, ".", codepaternaluncleID, ".", i)
        }
    datFinal <- rbind(datFinal,dat_paternaluncle)
    } else if (nrow(dat_paternaluncle) == 1){
        dat_paternaluncle$NID <- paste0(ftag, ".", codepaternaluncleID)
        datFinal <- rbind(datFinal,dat_paternaluncle)
    } else {
        print ("Paternal Uncle not present")
    }
    
    dat_paternalaunt <- mergecodes[which(mergecodes$Tag == famtag[ki] &
                        mergecodes$Relationship == "paternalaunt"),]
    codepaternalaunt <- samplecodes[which(samplecodes$Relation == "paternalaunt"),]
    codepaternalauntID <- unique(codepaternalaunt$nanoID)
    if (nrow(dat_paternalaunt) > 1){
        for (i    in 1:nrow(dat_paternalaunt)){
            dat_paternalaunt$NID[i] <- paste0(ftag, ".", codepaternalauntID, ".", i)
        }
    datFinal <- rbind(datFinal,dat_paternalaunt)
    } else if (nrow(dat_paternalaunt) == 1){
        dat_paternalaunt$NID <- paste0(ftag, ".", codepaternalauntID)
        datFinal <- rbind(datFinal,dat_paternalaunt)
    } else {
        print ("Paternal Grand Father not present")
    }
    
}
datFinal$NID <- as.character(paste(datFinal$ProjectID,"_",datFinal$NID, sep = ""))
##Writing output dataframe or text
        if (outMode == "Text"){
            #len <- length(l)
            fileOut <- paste0("Sample_NIDKeys_nanotatoR.csv")
            write.csv(datFinal, file.path(outpath = outpath, fileOut), 
                row.names = FALSE)
        }
        else if (outMode == "dataframe"){
            return(datFinal)
        }
        else {
            stop ("Your Outmode is absent !!!")
        }
}
#' Merging Dual and DLE, and adding nanotatoR relation ID
#'
#' @param labelType character. Type of labels used for mapping. 
#'        Choices are Dual, DLE and Both.
#' @param SVMerge_path    character. Path for the Dual labelled cmap
#' @param SVMerge_pattern    character. pattern of the dual files.
#' @param SE_path    character. Path for the Dual labelled cmap
#' @param SE_pattern    character. pattern of the dual files.
#' @param Samplecodes character. File containing relations and IDs 
#' associated to them.
#' @param mergeKey character. File containing sample ID and relation.
#' @param outpath character. Path where the merged samples are kept.
#' @param filename character. Output file name.
#' @param mergedKeyoutpath character. File path storing sample name and nanoID
#'        key information.
#' @param mergedKeyFname character. File name storing sample name and nanoID
#'        key information.
#' @param outputMode character. Mode of databse output. Text or dataframe.
#' @return Text files containg merged smaps from different samples
#' @examples
#' dat1 <- merging_SE_SVMerge (
#' labelType = c("SE"),
#' SE_path = system.file("extdata", "SoloFile/", package="nanotatoR"),
#' SE_pattern = "*_DLE1_*",
#' Samplecodes = system.file("extdata", "nanotatoR_sample_codes.csv", package="nanotatoR"),
#' mergeKey = system.file("extdata", "nanotatoR_control_sample_codes.csv", package="nanotatoR"),
#' outpath = system.file("extdata", package="nanotatoR"),
#' mergedKeyoutpath = system.file("extdata", package="nanotatoR"), 
#' mergedKeyFname = "Sample_index.csv",
#' filename= "nanotatoRControl.txt", 
#' outputMode = "dataframe")
#' @import hash
#' @import tidyverse
#' @export

merging_SE_SVMerge <- function(
    labelType = c("SVMerge", "SE", "Both", "SE_Cancer"),
    SVMerge_path ,SVMerge_pattern , 
    SE_path , SE_pattern,
    Samplecodes ,mergeKey,
    outpath , mergedKeyoutpath ,
    mergedKeyFname , filename,
    outputMode = c("dataframe", "Text")){
    ###Creating the sample ID relation connection
   if(labelType == "SE_Cancer"){
	    re <- read.csv(mergeKey)
	    re$sampleIds <- paste(re$ProjectID, "_", re$Tag, "_", re$SampleID, sep = "")
	    re$NID <- paste("NR_", re$Tag, sep ="")
	    datFinal <- re
		write.csv(datFinal, file.path(mergedKeyoutpath, mergedKeyFname), 
                row.names = FALSE)
	}else{
	    datFinal <- FamilyInfoPrep(Samplecodes, mergeKey, outMode = "dataframe")
        sampleIds <- as.character(unique(datFinal$SampleID))
        write.csv(datFinal, file.path(mergedKeyoutpath, mergedKeyFname), 
                row.names = FALSE)
	}
    ###Checking for labeltype and doing the calculation accordingly
    if(labelType == "Both"){
	    datFinal <- FamilyInfoPrep(Samplecodes, mergeKey, outMode = "dataframe")
        sampleIds <- as.character(unique(datFinal$SampleID))
        write.csv(datFinal, file.path(mergedKeyoutpath, mergedKeyFname), 
            row.names = FALSE)
        datDual <- mergingSMAP_SVMerge(path = SVMerge_path, 
            pattern = SVMerge_pattern, 
            outMode = "dataframe")
        datDLE <- mergingSMAP_SE(path = SE_path, 
            pattern = SE_pattern, 
            outMode = "dataframe")
        'datFinal <- FamilyInfoPrep(Samplecodes, mergeKey, outMode = "dataframe")
        sampleIds <- as.character(unique(datFinal$SampleID))'
        datfinal <- rbind(datDual,datDLE)
    }else if(labelType == "SE"){
        datDLE <- mergingSMAP_SE(path = SE_path, 
            pattern = SE_pattern, outMode = "dataframe")
            datfinal <- datDLE
    }else if(labelType == "SVMerge"){
        datDual <- mergingSMAP_SVMerge(path = SVMerge_path, pattern = SVMerge_pattern, outMode = "dataframe")
        datfinal <- datDual
    }else if(labelType == "SE_Cancer"){
	    datDLE <- mergingSMAP_SE(path = SE_path, 
            pattern = SE_pattern, outMode = "dataframe")
        datfinal <- datDLE
	}
	else{stop("Label Type Empty !!!!")}    
    ###Adding relation and NanoId information to the merged solo files
    dataFinal <- c()
    for (n in seq_along(sampleIds)){
        print(sampleIds[n])
        dataFinal_temp <- datfinal[which(datfinal$SampleID == sampleIds[n]), ]
        dataFinal_temp$nanoID <- rep(datFinal$NID[n], 
		     times = length(dataFinal_temp$SampleID))
        dataFinal_temp$ProjectID <- rep(datFinal$ProjectID[n], 
		times = length(dataFinal_temp$SampleID))
        dataFinal_temp$TagID <- rep(
            paste(datFinal$ProjectID[n], "_" ,datFinal$Tag[n], sep =""), 
            times = length(dataFinal_temp$SampleID))
        dataFinal <- rbind(dataFinal,dataFinal_temp)
    }
     #filename <- "UDN_DSD_merged_06232019.txt"
    ###writing the file
    if(outputMode == "Text"){
        write.table(dataFinal, file.path(outpath, filename), 
            row.names =FALSE, sep ="\t")
    }
    else if(outputMode == "dataframe"){
        return(dataFinal)
    }
    else{stop("Outmode incorrect !!!")}
}
    



