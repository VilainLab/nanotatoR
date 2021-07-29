#' Getting the data from annotated smaps to extract SV
#' information based on type of variants.
#'
#' @param input_fmt_geneList character. Choice of gene list input
#'        Text or Dataframe.
#' @param input_fmt_SV  character. Choice of gene list input
#'        Text or Dataframe.
#' @param smap  character. SV file name.
#' @param svData Dataframe Input data containing SV data.
#' @param dat_geneList Dataframe Input data containing geneList data.
#' @param fileName Character Name of file containing Gene List data.
#' @param outpath Character Directory to the output file.
#' @param outputFilename Character Output filename.
#' @param RZIPpath Character Path for the Rtools Zip package.
#' @param outputType Character. Variants in excel tabs or in different csv files.
#'        Options Excel or csv.
#' @param EnzymeType Character. Enzyme type used. Options SVMerge or SE.
#' @param fileprefix Character. fileprefix to use for each of the files in the directory.
#' @param directoryName Character. Directory name where individual SV files will be stored.
#' @param primaryGenesPresent boolean Checks whether the primary gene List is present or not.
#' @return Excel file containing the annotated SV map, tabs divided based on
#' type of SVs.
#' @examples
#' smapName <- "NA12878_Q.S_VAP_SVmerge_solo5.txt"
#' outputFilename <- "NA12878_Q.S_VAP_SVmerge_solo5_out"
#' smappath <- system.file("extdata", smapName, package = "nanotatoR")
#' outpath <- system.file("extdata", smapName, package = "nanotatoR")
#' RZIPpath <- system.file("extdata", "zip.exe", package = "nanotatoR")
#' smap = system.file("extdata", smapName, package="nanotatoR")
#' bedFile <- system.file("extdata", "HomoSapienGRCH19_lift37.bed", package="nanotatoR")
#' outpath <- system.file("extdata", package="nanotatoR")
#' directoryName <- system.file("extdata", package="nanotatoR")
#' datcomp<-overlapnearestgeneSearch(smap = smap, 
#'     bed=bedFile, inputfmtBed = "bed", outpath, 
#'     n = 3, returnMethod_bedcomp = c("dataFrame"), 
#'     input_fmt_SV = "Text",
#'     EnzymeType = "SVMerge", 
#'     bperrorindel = 3000, bperrorinvtrans = 50000)
#' hgpath=system.file("extdata", "GRCh37_hg19_variants_2016-05-15.txt", package="nanotatoR")
#' datDGV <- DGVfrequency (hgpath = hgpath, 
#'     smap_data = datcomp,
#'     win_indel_DGV = 10000, 
#'     EnzymeType = "SVMerge",
#'     input_fmt_SV = "dataFrame",
#'     perc_similarity_DGV = 0.5,returnMethod="dataFrame")
#'     indelconf = 0.5; invconf = 0.01;transconf = 0.1;
#' datInf <- internalFrequency_solo(smapdata = datDGV, 
#'     buildSVInternalDB=TRUE, win_indel=10000, 
#'     win_inv_trans=50000, 
#'     labelType = c("SE"),
#'     EnzymeType = "SVMerge",
#'     SE_path = system.file("extdata", "SoloFile/", package="nanotatoR"),
#'     SE_pattern = "*_DLE1_*", perc_similarity_parents =0.9,
#'     Samplecodes = system.file("extdata", "nanotatoR_sample_codes.csv", package="nanotatoR"),
#'     mergeKey = system.file("extdata", "nanotatoR_control_sample_codes.csv", package="nanotatoR"),
#'     mergedKeyoutpath = system.file("extdata", package="nanotatoR"), 
#'     mergedKeyFname = "Sample_index.csv",
#'     indexfile = system.file("extdata", mergedKeyFname, package="nanotatoR"),
#'     perc_similarity = 0.5, indelconf = 0.5, invconf = 0.01, 
#'     transconf = 0.1, limsize = 1000, win_indel_parents = 5000,
#'     win_inv_trans_parents=40000,
#'     returnMethod="dataFrame", input_fmt_SV = "dataFrame")
#' path <- system.file("extdata", "Bionano_config/", package = "nanotatoR")
#' pattern <- "*_hg19_*"
#' datBNDB <- BNDBfrequency(smapdata = datInf, 
#'     buildBNInternalDB=TRUE, 
#'     input_fmt_SV = "dataFrame",
#'     dbOutput="dataframe",
#'     EnzymeType = "SVMerge",
#'     BNDBpath = path, 
#'     BNDBpattern = pattern, 
#'     fname, outpath, 
#'     win_indel = 10000,
#'     win_inv_trans = 50000, 
#'     perc_similarity = 0.5, 
#'     indelconf = 0.5, 
#'     invconf = 0.01, 
#'     limsize = 1000,
#'     transconf = 0.1,
#'     returnMethod=c("dataFrame"))
#' decipherpath = system.file("extdata", "population_cnv.txt", package="nanotatoR")
#' datdecipher <- Decipherfrequency (decipherpath = decipherpath, 
#'     smap_data = datBNDB, win_indel = 10000, 
#'     perc_similarity = 0.5,returnMethod="dataFrame", 
#'     input_fmt_SV = "dataFrame", EnzymeType = c("SVMerge"))
#' run_bionano_filter_SVMerge_solo (input_fmt_geneList = c("Text"),
#'     input_fmt_SV = c("dataFrame"),
#'     svData = datdecipher, 
#'     dat_geneList = dat_geneList,
#'     RZIPpath = RZIPpath, EnzymeType = c("SVMerge"),
#'     outputType = c("csv"),
#'     primaryGenesPresent = FALSE, 
#'     directoryName = directoryName,
#'     fileprefix = "AnnotatedSamplesNA12878_SVMerge")
#' @import openxlsx
#' @import hash
#' @importFrom stats na.omit
#' @export
run_bionano_filter_SVMerge_solo <- function(
        input_fmt_geneList = c("Text", "dataFrame"),
        input_fmt_SV = c("Text", "dataFrame"),
        EnzymeType = c("SE", "SVMerge"),
        smap = NULL, svData, dat_geneList, fileName, outpath,
        outputFilename = "", RZIPpath, primaryGenesPresent = TRUE,
        outputType = c("Excel", "csv"), directoryName, fileprefix) {
    # library(openxlsx)
    # library(hash)
    # setwd(path)##change the directory to your working directory
    Sys.setenv("R_ZIPCMD" = RZIPpath)
    ## GeneList Input Format
    
    ha <- hash()
    
    # ll<-list.files(pattern="txt")
    if(input_fmt_SV=="dataFrame"){
        smapdata = svData
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
    r <- smapdata
    primaryGenesPresent = primaryGenesPresent
    if (primaryGenesPresent == TRUE){
        if (input_fmt_geneList == "Text") {
            rr <- read.table(fileName, header = TRUE)
        }
        else if (input_fmt_geneList == "dataFrame") {
            rr <- dat_geneList
        }
        else {
            stop("input_fmt_geneList Incorrect!!")
        }
    
        
        ogene <- as.character(r$OverlapGenes_strand_perc)
        upgene <- as.character(r$Upstream_nonOverlapGenes_dist_kb)
        dngene <- as.character(r$Downstream_nonOverlapGenes_dist_kb)
    # nogene<-c(upgene,dngene)
    
    ### OverlapGene
        dataPGOV <- overlappingGenes (rr, ogene)
    ### Non-Overlap Up-stream Gene
        dataPGUP <- nonOverlappingUPGenes (rr, upgene)
    ### Non-Overlap Down-stream Gene
        dataPGDN <- nonOverlappingDNGenes (rr, dngene)
    ### Non-OverlapDnGene
    
    # len<-length(pagene)-length(pg)
    # genesPG<-c(as.character(pg),rep("-",len))
        data <- data.frame(cbind(
        r, Overlap_PG = as.character(dataPGOV$pagene),
        Overlap_Terms = as.character(dataPGOV$pagene_term),
        Overlap_ClinicalSig = as.character(dataPGOV$pagene_clinSig),
        Non_Overlap_UP_PG = as.character(dataPGUP$nopageneup),
        Non_Overlap_UP_Terms = as.character(dataPGUP$nopageneup_term),
        Non_Overlap_UP_ClinicalSig = as.character(dataPGUP$nopageneup_clinSig),
        Non_Overlap_DN_PG = as.character(dataPGDN$nopagenedn),
        Non_Overlap_DN_Terms = as.character(dataPGDN$nopagenedn_term),
        Non_Overlap_DN_ClinicalSig = as.character(dataPGDN$nopagenedn_clinSig)
    ))
    }else if (primaryGenesPresent == FALSE){
        data <- smapdata
    }else {stop("primaryGenesPresent Incorrect!!")}
    ##
    data$BNG_Freq_Perc_Filtered<-gsub("-", 0, as.character(
    data$BNG_Freq_Perc_Filtered))
    data$BNG_Freq_Perc_UnFiltered<-gsub("-",0,as.character(
    data$BNG_Freq_Perc_UnFiltered))
    data$DGV_Freq_Perc<-as.numeric(data$DGV_Freq_Perc)
    data$Internal_Freq_Perc_Filtered<-as.numeric(
    data$Internal_Freq_Perc_Filtered)
    data$Internal_Freq_Perc_Unfiltered<-as.numeric(
        data$Internal_Freq_Perc_Unfiltered)
    data$BNG_Freq_Perc_Filtered<-as.numeric(data$BNG_Freq_Perc_Filtered)
    data$BNG_Freq_Perc_UnFiltered<-as.numeric(data$BNG_Freq_Perc_UnFiltered)
    data$DECIPHER_Frequency<-as.numeric(data$DECIPHER_Frequency)
    
    
    
    dat <- data[which(data$Type %in% "insertion"), ]
    dat1 <- data[which(data$Type %in% "deletion"), ]
    in_del<-rbind(dat,dat1)
    dat44 <- data[which(data$Type %in% "duplication"), ]
    dat45 <- data[which(data$Type %in% "duplication_split"), ]
    dat46 <- data[which(data$Type %in% "duplication_inverted"), ]
    dat_dup_rem<-rbind(dat44,dat45,dat46)
    #in_del_dup <- rbind(in_del, dat_dup_rem)
    dat_dup_rem_FINAL<-dat_dup_rem[which(((dat_dup_rem$Type %in% "duplication") |
         (dat_dup_rem$Type %in% "duplication_split") | (dat_dup_rem$Type %in% "duplication_inverted"))&
        (((dat_dup_rem$Fail_BSPQI_assembly_chimeric_score == "pass"
        & dat_dup_rem$Fail_BSSSI_assembly_chimeric_score == "pass") |
          (dat_dup_rem$Fail_BSPQI_assembly_chimeric_score == "fail" &
            dat_dup_rem$Fail_BSSSI_assembly_chimeric_score == "pass") |
          (dat_dup_rem$Fail_BSPQI_assembly_chimeric_score == "pass" &
            dat_dup_rem$Fail_BSSSI_assembly_chimeric_score == "fail")
        | (dat_dup_rem$Fail_BSPQI_assembly_chimeric_score == "pass"
          & dat_dup_rem$Fail_BSSSI_assembly_chimeric_score == "-")
        | (dat_dup_rem$Fail_BSPQI_assembly_chimeric_score == "-" &
            dat_dup_rem$Fail_BSSSI_assembly_chimeric_score == "pass")))), ]
        
    dat3 <- rbind(dat1, dat, dat_dup_rem_FINAL)
    dat10 <- dat3[which(((dat3$Found_in_self_BSPQI_molecules == "yes"
        & dat3$Found_in_self_BSSSI_molecules == "yes") |
        (dat3$Found_in_self_BSPQI_molecules == "no"
        & dat3$Found_in_self_BSSSI_molecules == "yes") |
        (dat3$Found_in_self_BSPQI_molecules == "yes"
        & dat3$Found_in_self_BSSSI_molecules == "no") |
        (dat3$Found_in_self_BSPQI_molecules == "yes"
        & dat3$Found_in_self_BSSSI_molecules == "-") |
        (dat3$Found_in_self_BSPQI_molecules == "-" &
          dat3$Found_in_self_BSSSI_molecules == "yes"))), ]
      
    'dat11 <- dat3[which(((dat3$Found_in_self_BSPQI_molecules == "yes" &
      dat3$Found_in_self_BSSSI_molecules == "yes") |
      (dat3$Found_in_self_BSPQI_molecules == "no" &
        dat3$Found_in_self_BSSSI_molecules == "yes") |
      (dat3$Found_in_self_BSPQI_molecules == "yes" &
        dat3$Found_in_self_BSSSI_molecules == "no") |
      (dat3$Found_in_self_BSPQI_molecules == "yes" &
        dat3$Found_in_self_BSSSI_molecules == "-") |
      (dat3$Found_in_self_BSPQI_molecules == "-" &
        dat3$Found_in_self_BSSSI_molecules == "yes")) &
      ((dat3$Found_in_parents_BSPQI_molecules == "both" &
        dat3$Found_in_parents_BSSSI_molecules == "both") |
        (dat3$Found_in_parents_BSPQI_molecules == "-" &
          dat3$Found_in_parents_BSSSI_molecules == "both") |
        (dat3$Found_in_parents_BSPQI_molecules == "both" &
          dat3$Found_in_parents_BSSSI_molecules == "-") | 
    (dat3$Found_in_parents_BSPQI_molecules == "both" &
          dat3$Found_in_parents_BSSSI_molecules == "mother") |
    (dat3$Found_in_parents_BSPQI_molecules == "mother" &
          dat3$Found_in_parents_BSSSI_molecules == "both") | 
    (dat3$Found_in_parents_BSPQI_molecules == "both" &
          dat3$Found_in_parents_BSSSI_molecules == "father") |
    (dat3$Found_in_parents_BSPQI_molecules == "father" &
          dat3$Found_in_parents_BSSSI_molecules == "both"))),]
    at12 <- dat3[which(((dat3$Found_in_self_BSPQI_molecules == "yes" &
      dat3$Found_in_self_BSSSI_molecules == "yes") |
      (dat3$Found_in_self_BSPQI_molecules == "no" &
        dat3$Found_in_self_BSSSI_molecules == "yes") |
      (dat3$Found_in_self_BSPQI_molecules == "yes" &
        dat3$Found_in_self_BSSSI_molecules == "no") |
      (dat3$Found_in_self_BSPQI_molecules == "yes" &
        dat3$Found_in_self_BSSSI_molecules == "-") |
      (dat3$Found_in_self_BSPQI_molecules == "-" &
        dat3$Found_in_self_BSSSI_molecules == "yes")) &
      ((dat3$Found_in_parents_BSPQI_molecules == "mother" &
        dat3$Found_in_parents_BSSSI_molecules == "mother") |
        (dat3$Found_in_parents_BSPQI_molecules == "-" &
          dat3$Found_in_parents_BSSSI_molecules == "mother") |
        (dat3$Found_in_parents_BSPQI_molecules == "mother" &
          dat3$Found_in_parents_BSSSI_molecules == "-"))), ]
    dat13 <- dat3[which(((dat3$Found_in_self_BSPQI_molecules == "yes" &
      dat3$Found_in_self_BSSSI_molecules == "yes") |
      (dat3$Found_in_self_BSPQI_molecules == "no" &
        dat3$Found_in_self_BSSSI_molecules == "yes") |
      (dat3$Found_in_self_BSPQI_molecules == "yes" &
        dat3$Found_in_self_BSSSI_molecules == "no") |
      (dat3$Found_in_self_BSPQI_molecules == "yes" &
        dat3$Found_in_self_BSSSI_molecules == "-") |
      (dat3$Found_in_self_BSPQI_molecules == "-" &
        dat3$Found_in_self_BSSSI_molecules == "yes")) &
      ((dat3$Found_in_parents_BSPQI_molecules == "father" &
        dat3$Found_in_parents_BSSSI_molecules == "father") |
        (dat3$Found_in_parents_BSPQI_molecules == "-" &
          dat3$Found_in_parents_BSSSI_molecules == "father") |
        (dat3$Found_in_parents_BSPQI_molecules == "father" &
          dat3$Found_in_parents_BSSSI_molecules == "-"))), ]
    dat14 <- rbind(dat11, dat12, dat13)'
    
    #gg <- grep("inversion", as.character(data$Type))
    dat8 <- data[which(data$Type %in% "inversion"
            | data$Type %in% "inversion_partial"
            | data$Type %in% "inversion_paired"
            | data$Type %in% "inversion_repeat"), ]
    dat8 <- dat8[which(((dat8$Found_in_self_BSPQI_molecules == "yes"
        & dat8$Found_in_self_BSSSI_molecules == "yes") |
        (dat8$Found_in_self_BSPQI_molecules == "no" &
        dat8$Found_in_self_BSSSI_molecules == "yes") |
        (dat8$Found_in_self_BSPQI_molecules == "yes" &
        dat8$Found_in_self_BSSSI_molecules == "no") |
        (dat8$Found_in_self_BSPQI_molecules == "yes" &
        dat8$Found_in_self_BSSSI_molecules == "-") |
        (dat8$Found_in_self_BSPQI_molecules == "-" &
        dat8$Found_in_self_BSSSI_molecules == "yes")) &
        ((dat8$Fail_BSPQI_assembly_chimeric_score == "pass" &
        dat8$Fail_BSSSI_assembly_chimeric_score == "pass") |
        (dat8$Fail_BSPQI_assembly_chimeric_score == "fail"
        & dat8$Fail_BSSSI_assembly_chimeric_score == "pass") |
        (dat8$Fail_BSPQI_assembly_chimeric_score == "pass" &
        dat8$Fail_BSSSI_assembly_chimeric_score == "fail")
        | (dat8$Fail_BSPQI_assembly_chimeric_score == "pass"
        & dat8$Fail_BSSSI_assembly_chimeric_score == "-")
        | (dat8$Fail_BSPQI_assembly_chimeric_score == "-" &
        dat8$Fail_BSSSI_assembly_chimeric_score == "pass"))), ]
    gg1 <- grep("translocation", as.character(data$Type))
    dat7 <- data[gg1, ]
    dat7 <- dat7[which(((dat7$Found_in_self_BSPQI_molecules == "yes" &
        dat7$Found_in_self_BSSSI_molecules == "yes") |
        (dat7$Found_in_self_BSPQI_molecules == "no" &
        dat7$Found_in_self_BSSSI_molecules == "yes") |
        (dat7$Found_in_self_BSPQI_molecules == "yes" &
        dat7$Found_in_self_BSSSI_molecules == "no") |
        (dat7$Found_in_self_BSPQI_molecules == "yes" &
        dat7$Found_in_self_BSSSI_molecules == "-") |
        (dat7$Found_in_self_BSPQI_molecules == "-" &
        dat7$Found_in_self_BSSSI_molecules == "yes")) &
        ((dat7$Fail_BSPQI_assembly_chimeric_score == "pass"
        & dat7$Fail_BSSSI_assembly_chimeric_score == "pass") |
        (dat7$Fail_BSPQI_assembly_chimeric_score == "fail" &
        dat7$Fail_BSSSI_assembly_chimeric_score == "pass") |
        (dat7$Fail_BSPQI_assembly_chimeric_score == "pass" &
        dat7$Fail_BSSSI_assembly_chimeric_score == "fail")
        | (dat7$Fail_BSPQI_assembly_chimeric_score == "pass"
        & dat7$Fail_BSSSI_assembly_chimeric_score == "-")
        | (dat7$Fail_BSPQI_assembly_chimeric_score == "-" &
        dat7$Fail_BSSSI_assembly_chimeric_score == "pass"))), ]
    dat6 <- data[which(data$Type %in% "MisMatch"), ]
    #dat14 <- rbind(dat11, dat12, dat13)
    ovrlapPG<-which(!((as.character(data$Overlap_PG)=="-")))
    datovrlapPG<-data[ovrlapPG,]
    nonovrlapupPG<-which(!((as.character(data$Non_Overlap_UP_PG)=="-")))
    datnonovrlapUPPG<-data[nonovrlapupPG,]
    nonovrlapdnPG<-which(!((as.character(data$Non_Overlap_DN_PG)=="-")))
    datnonovrlapDNPG<-data[nonovrlapdnPG,]
    datOvrLap <- rbind(datovrlapPG, datnonovrlapUPPG, datnonovrlapDNPG)
      
      'list_of_datasets <- list(
        "all" = data, "indel_dup" = dat10,
        "indel_dup_both" = dat11,
        "indel_dup_mother" = dat12, "indel_dup_father" = dat13,
        "indel_dup_cmpdHET" = dat14, "inv" = dat8, "trans" = dat7, "mismatch" = dat6,
        "all_PG_OV" = datovrlapPG, "all_PG_Non_OV_UP" = datnonovrlapUPPG,"all_PG_Non_OV_DN" = datnonovrlapDNPG
      )'
    if (outputType == "Excel"){
        list_of_datasets <- list(
            "indel_dup" = dat10, "inv" = dat8, "trans" = dat7, 
            "all_PG_OV" = datOvrLap, "mismatch" = dat6, "all" = data
        )
        fname <- paste(outputFilename, ".xlsx", sep = "")
        write.xlsx(list_of_datasets, file = file.path(
            outpath, fname), keepNA = TRUE)
    } else if (outputType == "csv"){
        write.csv(dat10, file.path(directoryName, paste(
            fileprefix,"_indel_dup.csv",sep = "")), row.names = FALSE)
        write.csv(dat8, file.path(directoryName, paste(
            fileprefix,"_inv.csv",sep = "")), row.names = FALSE)
        write.csv(dat7, file.path(directoryName, paste(
            fileprefix,"_trans.csv",sep = "")), row.names = FALSE)
        write.csv(datOvrLap, file.path(directoryName, paste(
            fileprefix,"_all_PG_OV.csv",sep = "")), row.names = FALSE)
        write.csv(dat6, file.path(directoryName, paste(
            fileprefix,"_mismatch.csv",sep = "")), row.names = FALSE)
        write.csv(data, file.path(directoryName, paste(
            fileprefix,"_all.csv",sep = "")), row.names = FALSE)
    } else {stop(" outputType incorrect !!")}
    
      # gg<-grep("inversion",as.character(data$Type))
      # dat8<-data[gg,]
      # gg1<-grep("translocation",as.character(data$Type))
      # dat7<-data[gg1,]
      # dat6<-data[which(data$Type %in% "MisMatch"),]
      # list_of_datasets <- list("all" = data, "indel" = dat3, "indel_denovo"= dat10, "indel_both"= dat11,
      #                        "indel_mother"= dat12, "indel father"= dat13, "indel_cmpdHET"=dat14, "inv"= dat8, "trans"=dat7, "mismatch"=dat6,
      #                       "all_PG_OV"=data, "all_PG_Non_OV"= data)
      # st1<-strsplit(smap,".txt")[[1]][1]
                   
    
    
    # gg<-grep("inversion",as.character(data$Type))
    # dat8<-data[gg,]
    # gg1<-grep("translocation",as.character(data$Type))
    # dat7<-data[gg1,]
    # dat6<-data[which(data$Type %in% "MisMatch"),]
    # list_of_datasets <- list("all" = data, "indel" = dat3, "indel_denovo"= dat10, "indel_both"= dat11,
    #"indel_mother"= dat12, "indel father"= dat13, "indel_cmpdHET"=dat14, "inv"= dat8, "trans"=dat7, "mismatch"=dat6,
    # "all_PG_OV"=data, "all_PG_Non_OV"= data)
    # st1<-strsplit(smap,".txt")[[1]][1]
    
    }
    
    
    ### Running the code
    ### Copy this 1 line at a time without the hash (#)
    ### Warning1: Don't remove the # from before the command from filterread.r
    ### Warning2: give the fullpath with filename for theR code
    # source("Z:\\bionano\\Test_R_Codes\\filterread.r")##path of where the Rcode is
    #
    # run_bionano_filter(path,smap,fileName)
    