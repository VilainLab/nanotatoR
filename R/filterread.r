#' Getting the data from annotated smaps to extract SV
#' information based on type of variants.
#'
#' @param input_fmt_geneList character. Choice of gene list input
#'        Text or Dataframe.
#' @param input_fmt_svMap  character. Choice of gene list input
#'        Text or Dataframe.
#' @param SVFile  character. SV file name.
#' @param svData Dataframe Input data containing SV data.
#' @param dat_geneList Dataframe Input data containing geneList data.
#' @param fileName Character Name of file containing Gene List data.
#' @param outpath Character Directory to the output file.
#' @param outputFilename Character Output filename.
#' @param RZIPpath Character Path for the Rtools Zip package.
#' @return Excel file containing the annotated SV map, tabs divided based on
#' type of SVs.
#' @examples
#' \dontrun{
#' terms <- "Muscle Weakness"
#' gene <- gene_list_generation(
#'   method = "Single", term = terms,
#'   returnMethod_GeneList = "dataFrame"
#' )
#' RzipFile = "zip.exe"
#' RZIPpath <- system.file("extdata", RzipFile, package = "nanotatoR")
#' smapName <- "F1.1_TestSample1_solo_hg19.smap"
#' smappath <- system.file("extdata", smapName, package = "nanotatoR")
#' smappath1 <- system.file("extdata", package = "nanotatoR")
#' run_bionano_filter(input_fmt_geneList = c("dataframe"), input_fmt_svMap = c("Text"),
#'   SVFile = smappath, dat_geneList = gene, outpath = smappath1, outputFilename = "test", 
#'   RZIPpath = RZIPpath)
#' }
#' @import openxlsx
#' @import hash
#' @importFrom stats na.omit
#' @export


run_bionano_filter <- function(input_fmt_geneList = c("Text", "dataFrame"),
        input_fmt_svMap = c("Text", "dataFrame"),
        SVFile = NULL, svData, dat_geneList, fileName, outpath,
        outputFilename = "", RZIPpath) {
    # library(openxlsx)
    # library(hash)
    # setwd(path)##change the directory to your working directory
    Sys.setenv("R_ZIPCMD" = RZIPpath)
    ## GeneList Input Format
    if (input_fmt_geneList == "Text") {
        rr <- read.table(fileName, header = TRUE)
    }
    else if (input_fmt_geneList == "dataFrame") {
        rr <- dat_geneList
    }
    else {
        stop("Dataframe Incorrect!!")
    }
    ha <- hash()
    .set(ha, keys = rr$Genes, values = rr$Terms)

    pg <- as.character(rr$Genes)
    gen <- paste("^", pg, "$", sep = "")
    # ll<-list.files(pattern="txt")
    if (input_fmt_svMap == "Text") {
        r <- read.table(SVFile, sep = "\t", header = TRUE)
    }
    else if (input_fmt_svMap == "dataFrame") {
        r <- svData
    }
    else {
        stop("Input Formats Incorrect")
    }
    ogene <- as.character(r$OverlapGenes_strand_perc)
    upgene <- as.character(r$Upstream_nonOverlapGenes_dist_kb)
    dngene <- as.character(r$Downstream_nonOverlapGenes_dist_kb)
    # nogene<-c(upgene,dngene)
    pagene <- c()
    pagene_term <- c()
    nopageneup <- c()
    nopageneup_term <- c()
    nopagenedn <- c()
    nopagenedn_term <- c()
    ### OverlapGene
    for (ii in seq_len(length(ogene))) {
        dd <- grep(";", ogene[ii])
        gene_Terms <- c()
        if (length(dd) > 0) {
            st <- strsplit(ogene[ii], ";")
            gen2 <- as.character(unlist(st))
            opag <- c()

            for (l in seq_len(length(gen2))) {
                st1 <- st <- strsplit(gen2[l], "\\(")
                gen2[l] <- st1[[1]][1]
                pa <- paste("^", gen2[l], "$", sep = "")
                gg <- grep(pa, gen, fixed = TRUE)
                if (length(gg) > 0) {
                    opag <- c(opag, as.character(gen2[l]))
                    val1 <- hash::values(ha, keys = gen2[l])
                    gene_Terms <- c(gene_Terms, as.character(val1))
                }
                else {
                    opag <- c(opag, "-")
                }
            }
            opag <- as.character(unique(opag))
            opagt <- as.character(unique(gene_Terms))
            opagt <- paste(opagt, collapse = ";")
        if (length(opag) > 1) {
            opagpa <- paste(opag, collapse = ";")
            ff <- grep("-", opagpa)
            if (length(ff) > 0) {
                opagpa <- gsub("-", "", opagpa)
                opagpa <- gsub(";;", ";", opagpa)
                opagpa <- gsub("^;", "", opagpa)
                pagene <- c(pagene, opagpa)
                pagene_term <- c(pagene_term, opagt)
            }
            else {
                pagene <- c(pagene, opagpa)
                pagene_term <- c(pagene_term, opagt)
            }
        }
        else {
            pagene <- c(pagene, opag)
            pagene_term <- c(pagene_term, opagt)
        }
    }
    else {
        st1 <- st <- strsplit(ogene[ii], "\\(")
        ogene[ii] <- st1[[1]][1]
        pa <- paste("^", ogene[ii], "$", sep = "")
        gg <- grep(pa, gen, fixed = TRUE)
        if (length(gg) > 0) {
            pagene <- c(pagene, as.character(ogene[ii]))
            val1 <- hash::values(ha, keys = ogene[ii])
            pagene_term <- c(pagene_term, as.character(val1))
        }
        else {
            pagene <- c(pagene, "-")
            pagene_term <- c(pagene_term, "-")
        }
    }
}
    ### Non-OverlapUpGene
    for (jj in seq_len(length(upgene))) {
        dd1 <- grep(";", upgene[jj])
        gene_Terms_no <- c()
        if (length(dd1) > 0) {
            st1 <- strsplit(upgene[jj], ";")
            gen3 <- as.character(unlist(st1))
            nopag <- c()

            for (l in seq_len(length(gen3))) {
                st1 <- st <- strsplit(gen3[l], "\\(")
                gen3[l] <- st1[[1]][1]
                pa1 <- paste("^", gen3[l], "$", sep = "")
                gg2 <- grep(pa1, gen, fixed = TRUE)
                if (length(gg2) > 0) {
                    nopag <- c(nopag, as.character(gen3[l]))
                    val2 <- hash::values(ha, keys = gen3[l])
                    gene_Terms_no <- c(gene_Terms_no, as.character(val2))
                }
                else {
                    nopag <- c(nopag, "-")
                }
            }
        nopag <- as.character(unique(nopag))
        nopagt <- as.character(unique(gene_Terms_no))
        nopagt <- paste(nopagt, collapse = ";")
        if (length(nopag) > 1) {
            nopagpa <- paste(nopag, collapse = ";")
            ff <- grep("-", nopagpa)
            if (length(ff) > 0) {
                nopagpa <- gsub("-", "", nopagpa)
                nopagpa <- gsub(";;", ";", nopagpa)
                nopagpa <- gsub("^;", ";", nopagpa)
                nopageneup <- c(nopageneup, nopagpa)
                nopageneup_term <- c(nopageneup_term, nopagt)
            }
            else {
                nopageneup <- c(nopageneup, nopagpa)
                nopageneup_term <- c(nopageneup_term, nopagt)
            }
        }
        else {
            nopageneup <- c(nopageneup, nopag)
            nopageneup_term <- c(nopageneup_term, nopagt)
        }
    }
    else {
        st1 <- st <- strsplit(upgene[ii], "\\(")
        upgene[ii] <- st1[[1]][1]
        pa <- paste("^", upgene[jj], "$", sep = "")
        gg <- grep(pa, gen, fixed = TRUE)
        if (length(gg) > 0) {
            nopageneup <- c(nopageneup, as.character(upgene[jj]))
            val2 <- hash::values(ha, keys = upgene[jj])
            nopageneup_term <- c(nopageneup_term, as.character(val2))
        }
        else {
            nopageneup <- c(nopageneup, "-")
            nopageneup_term <- c(nopageneup_term, "-")
        }
    }
}
    ### Non-OverlapDnGene
    for (jj in seq_len(length(dngene))) {
        dd1 <- grep(";", dngene[jj])
        gene_Terms_no <- c()
        if (length(dd1) > 0) {
            st1 <- strsplit(dngene[jj], ";")
            gen3 <- as.character(unlist(st1))
            nopag <- c()

            for (l in seq_len(length(gen3))) {
                st1 <- st <- strsplit(gen3[l], "\\(")
                gen3[l] <- st1[[1]][1]
                pa1 <- paste("^", gen3[l], "$", sep = "")
                gg2 <- grep(pa1, gen, fixed = TRUE)
                if (length(gg2) > 0) {
                    nopag <- c(nopag, as.character(gen3[l]))
                    val2 <- hash::values(ha, keys = gen3[l])
                    gene_Terms_no <- c(gene_Terms_no, as.character(val2))
                }
                else {
                    nopag <- c(nopag, "-")
                }
            }
        nopag <- as.character(unique(nopag))
        nopagt <- as.character(unique(gene_Terms_no))
        nopagt <- paste(nopagt, collapse = ";")
        if (length(nopag) > 1) {
            nopagpa <- paste(nopag, collapse = ";")
            ff <- grep("-", nopagpa)
            if (length(ff) > 0) {
                nopagpa <- gsub("-", "", nopagpa)
                nopagpa <- gsub(";;", ";", nopagpa)
                nopagpa <- gsub("^;", ";", nopagpa)
                nopagenedn <- c(nopagenedn, nopagpa)
                nopagenedn_term <- c(nopagenedn_term, nopagt)
            }
            else {
                nopagenedn <- c(nopagenedn, nopagpa)
                nopagenedn_term <- c(nopagenedn_term, nopagt)
            }
        }
        else {
            nopagenedn <- c(nopagenedn, nopag)
            nopagenedn_term <- c(nopagenedn_term, nopagt)
        }
    }
    else {
        st1 <- st <- strsplit(dngene[ii], "\\(")
        dngene[ii] <- st1[[1]][1]
        pa <- paste("^", dngene[jj], "$", sep = "")
        gg <- grep(pa, gen, fixed = TRUE)
        if (length(gg) > 0) {
            nopagenedn <- c(nopagenedn, as.character(dngene[jj]))
            val2 <- hash::values(ha, keys = dngene[jj])
            nopagenedn_term <- c(nopagenedn_term, as.character(val2))
        }
        else {
            nopagenedn <- c(nopagenedn, "-")
            nopagenedn_term <- c(nopagenedn_term, "-")
        }
    }
}
    # len<-length(pagene)-length(pg)
    # genesPG<-c(as.character(pg),rep("-",len))
    data <- data.frame(cbind(
        r, Overlap_PG = as.character(pagene),
        Overlap_Terms = as.character(pagene_term),
        Non_Overlap_UP_PG = as.character(nopageneup),
        Non_Overlap_UP_Terms = as.character(nopageneup_term),
        Non_Overlap_DN_PG = as.character(nopagenedn),
        Non_Overlap_DN_Terms = as.character(nopagenedn_term)
    ))
    data$BNG_Freq_Perc_Filtered<-gsub("-",0,as.character(data$BNG_Freq_Perc_Filtered))
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
    data$Type1<-as.character(data$Type1)
    data$Type2<-as.character(data$Type2)
  
  
    dat <- data[which(data$Type %in% "insertion"), ]
    dat1 <- data[which(data$Type %in% "deletion"), ]
    dat44 <- data[which(data$Type %in% "duplication"), ]
    dat45 <- data[which(data$Type %in% "duplication_split"), ]
    dat46 <- data[which(data$Type %in% "duplication_inverted"), ]
    dat_dup_rem<-rbind(dat44,dat45,dat46)
    dat_dup_rem_FINAL<-dat_dup_rem[which(((dat_dup_rem$Type %in% "duplication") |
        (
        dat_dup_rem$Type %in% "duplication_split") | 
        (dat_dup_rem$Type %in% "duplication_inverted"))&
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
        dat3$Found_in_self_BSSSI_molecules == "yes")) &
        ((dat3$Found_in_parents_BSPQI_molecules == "none"
        & dat3$Found_in_parents_BSSSI_molecules == "none") |
        (dat3$Found_in_parents_BSPQI_molecules == "none"
        & dat3$Found_in_parents_BSSSI_molecules == "-") |
        (dat3$Found_in_parents_BSPQI_molecules == "-" &
        dat3$Found_in_parents_BSSSI_molecules == "none"))), ]
    
    dat11 <- dat3[which(((dat3$Found_in_self_BSPQI_molecules == "yes" &
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
    dat12 <- dat3[which(((dat3$Found_in_self_BSPQI_molecules == "yes" &
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
    dat14 <- rbind(dat11, dat12, dat13)

    gg <- grep("inversion", as.character(data$Type))
    dat8 <- data[gg, ]
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
    dat14 <- rbind(dat11, dat12, dat13)
    ovrlapPG<-which(!((as.character(data$Overlap_PG)=="-")))
    datovrlapPG<-data[ovrlapPG,]
    nonovrlapupPG<-which(!((as.character(data$Non_Overlap_UP_PG)=="-")))
    datnonovrlapUPPG<-data[nonovrlapupPG,]
    nonovrlapdnPG<-which(!((as.character(data$Non_Overlap_DN_PG)=="-")))
    datnonovrlapDNPG<-data[nonovrlapdnPG,]
  
  
  
    list_of_datasets <- list(
        "all" = data, "indel_dup_denovo" = dat10,
        "indel_dup_both" = dat11,
        "indel_dup_mother" = dat12, "indel_dup_father" = dat13,
        "indel_dup_cmpdHET" = dat14, "inv" = dat8, "trans" = dat7, 
        "mismatch" = dat6,
        "all_PG_OV" = datovrlapPG)


  # gg<-grep("inversion",as.character(data$Type))
  # dat8<-data[gg,]
  # gg1<-grep("translocation",as.character(data$Type))
  # dat7<-data[gg1,]
  # dat6<-data[which(data$Type %in% "MisMatch"),]
  # list_of_datasets <- list("all" = data, "indel" = dat3, "indel_denovo"= dat10, "indel_both"= dat11,
  #                        "indel_mother"= dat12, "indel father"= dat13, "indel_cmpdHET"=dat14, "inv"= dat8, "trans"=dat7, "mismatch"=dat6,
  #                       "all_PG_OV"=data, "all_PG_Non_OV"= data)
  # st1<-strsplit(SVFile,".txt")[[1]][1]
  fname <- paste(outputFilename, ".xlsx", sep = "")
  write.xlsx(list_of_datasets, file = file.path(outpath, fname),keepNA=TRUE)
}



