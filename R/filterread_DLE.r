#' Extracting terms for genes that overlap SVs
#'
#' @param rr character. dataframe with primary genes and terms
#'      associated.
#' @param ogene  character. genes that overlap the SV.
#' @return Dataframe with overlapping genes and terms.
#' @examples
#' \dontrun{
#' terms="Muscle Weakness"
#' dat_geneList <- gene_list_generation(
#'      method_entrez = c("Single"), 
#'      term = terms,
#'      returnMethod=c("dataFrame"), 
#'		omim = "Y:/Suro/nanotatoRDatabases/mim2gene.txt", 
#'		clinvar = "Y:/Suro/nanotatoRDatabases/clinvar_07172019.txt",
#		gtr = "Y:/Suro/nanotatoRDatabases/gtr_07162019.txt")
#' rr <- dat_geneList
#' smapName="F1.1_GM24385_DLE-1_P_trio_hg19.smap"
#' smap = system.file("extdata", smapName, package="nanotatoR")
#' bedFile <- system.file("extdata", "Homo_sapiens.Hg19.bed", package="nanotatoR")
#' outpath <- system.file("extdata",  package="nanotatoR")
#' datcomp<-compSmapbed(smap, bed=bedFile, inputfmtBed =  "BED", outpath,
#' n = 3, returnMethod_bedcomp = c("dataFrame"))
#' ogene <- as.character(datcomp$OverlapGenes_strand_perc)
#' datogenes <- overlappingGenes (rr, ogene)
#' )
#' }
#' @import hash
#' @importFrom stats na.omit
#' @export
overlappingGenes <- function(rr, ogene){
ha <- hash()
ha1 <- hash()
.set(ha, keys = rr$Genes, values = rr$Terms)
.set(ha1, keys = rr$Genes, values = rr$ClinicalSignificance)
pg <- as.character(rr$Genes)
pagene <- c()
pagene_term <- c()
pagene_clinSig <- c()
gen <- paste("^", pg, "$", sep = "")
for (ii in seq_len(length(ogene))) {
    dd <- grep(";", ogene[ii])
    gene_Terms <- c()
	clinSig <- c()
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
			val2 <- hash::values(ha1, keys = gen2[l])
            clinSig <- c(clinSig, as.character(val2))
        }
        else {
            opag <- c(opag, "-")
			gene_Terms <- c(gene_Terms, "-")
			clinSig <- c(clinSig, "-")
        }
      }
    opag <- as.character(unique(opag))
    opagt <- as.character(unique(gene_Terms))
    opagt <- paste(opagt, collapse = ";")
	opagcs <- as.character(unique(clinSig))
    opagcs <- paste(opagcs, collapse = ";")
        if (length(opag) > 1) {
            opagpa <- paste(opag, collapse = ";")
            ff <- grep("-", opagpa)
            if (length(ff) > 0) {
                opagpa <- gsub("-", "", opagpa)
                opagpa <- gsub(";;", ";", opagpa)
                opagpa <- gsub("^;", "", opagpa)
                pagene <- c(pagene, opagpa)
                pagene_term <- c(pagene_term, opagt)
				pagene_clinSig <- c(pagene_clinSig, opagcs)
            }
            else {
                pagene <- c(pagene, opagpa)
                pagene_term <- c(pagene_term, opagt)
				pagene_clinSig <- c(pagene_clinSig, opagcs)
            }
        }
        else {
            pagene <- c(pagene, opag)
            pagene_term <- c(pagene_term, opagt)
			pagene_clinSig <- c(pagene_clinSig, opagcs)
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
			val2 <- hash::values(ha1, keys = ogene[ii])
            pagene_clinSig <- c(pagene_clinSig, as.character(val2))
        }
        else {
            pagene <- c(pagene, "-")
            pagene_term <- c(pagene_term, "-")
			pagene_clinSig <- c(pagene_clinSig, "-")
        }
    }
}
    dataPGOV <- data.frame(pagene, pagene_term, pagene_clinSig)
    return(dataPGOV)
  
}

#' Extracting terms for genes that overlap SVs
#'
#' @param rr character. dataframe with primary genes and terms
#'      associated.
#' @param ogene  character. genes that overlap the SV.
#' @return Dataframe with overlapping genes and terms.
#' @examples
#' \dontrun{
#' terms="Muscle Weakness"
#' dat_geneList <- gene_list_generation(
#'      method_entrez = c("Single"), 
#'      term = terms,
#'      returnMethod=c("dataFrame"), 
#'		omim = "Y:/Suro/nanotatoRDatabases/mim2gene.txt", 
#'		clinvar = "Y:/Suro/nanotatoRDatabases/clinvar_07172019.txt",
#		gtr = "Y:/Suro/nanotatoRDatabases/gtr_07162019.txt")
#' rr <- dat_geneList
#' smapName="F1.1_GM24385_DLE-1_P_trio_hg19.smap"
#' smap = system.file("extdata", smapName, package="nanotatoR")
#' bedFile <- system.file("extdata", "Homo_sapiens.Hg19.bed", package="nanotatoR")
#' outpath <- system.file("extdata",  package="nanotatoR")
#' datcomp <- compSmapbed(smap, bed=bedFile, inputfmtBed =  "BED", outpath,
#' n = 3, returnMethod_bedcomp = c("dataFrame"))
#' upgene <- as.character(datcomp$Upstream_nonOverlapGenes_dist_kb)
#' dataPGUP <- nonOverlappingUPGenes (rr, upgene)
#'}
#' @import hash
#' @importFrom stats na.omit
#' @export
nonOverlappingUPGenes <- function(rr, upgene){
ha <- hash()
ha1 <- hash()
.set(ha, keys = rr$Genes, values = rr$Terms)
.set(ha1, keys = rr$Genes, values = rr$ClinicalSignificance)
pg <- as.character(rr$Genes)
nopageneup <- c()
nopageneup_term <- c()
nopageneup_clinSig <- c()
gen <- paste("^", pg, "$", sep = "")
for (jj in seq_len(length(upgene))) {
    dd1 <- grep(";", upgene[jj])
    gene_Terms_no <- c()
	up_clinSig <- c()
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
				val4 <- hash::values(ha1, keys = gen3[l])
                up_clinSig <- c(up_clinSig, as.character(val4))
            }
            else {
                nopag <- c(nopag, "-")
				gene_Terms_no <- c(gene_Terms_no, "-")
			    up_clinSig <- c(up_clinSig, "-")
            }
        }
    nopag <- as.character(unique(nopag))
    nopagt <- as.character(unique(gene_Terms_no))
    nopagt <- paste(nopagt, collapse = ";")
	nopagupcs <- as.character(unique(up_clinSig))
    nopagupcs <- paste(nopagupcs, collapse = ";")
      if (length(nopag) > 1) {
        nopagpa <- paste(nopag, collapse = ";")
        ff <- grep("-", nopagpa)
        if (length(ff) > 0) {
            nopagpa <- gsub("-", "", nopagpa)
            nopagpa <- gsub(";;", ";", nopagpa)
            nopagpa <- gsub("^;", ";", nopagpa)
            nopageneup <- c(nopageneup, nopagpa)
            nopageneup_term <- c(nopageneup_term, nopagt)
			nopageneup_clinSig <- c(nopageneup_clinSig, nopagupcs)
        }
        else {
            nopageneup <- c(nopageneup, nopagpa)
            nopageneup_term <- c(nopageneup_term, nopagt)
			nopageneup_clinSig <- c(nopageneup_clinSig, nopagupcs)
        }
      }
      else {
            nopageneup <- c(nopageneup, nopag)
            nopageneup_term <- c(nopageneup_term, nopagt)
			nopageneup_clinSig <- c(nopageneup_clinSig, nopagupcs)
        }
    }
    else {
        st1 <- st <- strsplit(upgene[jj], "\\(")
        upgene[jj] <- st1[[1]][1]
        pa <- paste("^", upgene[jj], "$", sep = "")
        gg <- grep(pa, gen, fixed = TRUE)
        if (length(gg) > 0) {
            nopageneup <- c(nopageneup, as.character(upgene[jj]))
            val2 <- hash::values(ha, keys = upgene[jj])
            nopageneup_term <- c(nopageneup_term, as.character(val2))
			val4 <- hash::values(ha1, keys = upgene[jj])
            nopageneup_clinSig <- c(nopageneup_clinSig, as.character(val4))
        }
        else {
            nopageneup <- c(nopageneup, "-")
            nopageneup_term <- c(nopageneup_term, "-")
			nopageneup_clinSig <- c(nopageneup_clinSig, "-")
        }
    }
  }
    dataPGUP <- data.frame(nopageneup, nopageneup_term, nopageneup_clinSig)
    return(dataPGUP)
}
#' Extracting terms for genes that overlap SVs
#'
#' @param rr character. dataframe with primary genes and terms
#'      associated.
#' @param ogene  character. genes that overlap the SV.
#' @return Dataframe with overlapping genes and terms.
#' @examples
#' \dontrun{
#' terms="Muscle Weakness"
#' dat_geneList <- gene_list_generation(
#'      method_entrez = c("Single"), 
#'      term = terms,
#'      returnMethod=c("dataFrame"), 
#'		omim = "Y:/Suro/nanotatoRDatabases/mim2gene.txt", 
#'		clinvar = "Y:/Suro/nanotatoRDatabases/clinvar_07172019.txt",
#		gtr = "Y:/Suro/nanotatoRDatabases/gtr_07162019.txt")
#' rr <- dat_geneList
#' smapName="F1.1_GM24385_DLE-1_P_trio_hg19.smap"
#' smap = system.file("extdata", smapName, package="nanotatoR")
#' bedFile <- system.file("extdata", "Homo_sapiens.Hg19.bed", package="nanotatoR")
#' outpath <- system.file("extdata",  package="nanotatoR")
#' datcomp <- compSmapbed(smap, bed=bedFile, inputfmtBed =  "BED", outpath,
#' n = 3, returnMethod_bedcomp = c("dataFrame"))
#' dngene <- as.character(r$Downstream_nonOverlapGenes_dist_kb)
#' dataPGDN <- nonOverlappingDNGenes (rr, dngene)
#'}
#' @import hash
#' @importFrom stats na.omit
#' @export

nonOverlappingDNGenes <- function(rr, dngene){
nopagenedn <- c()
nopagenedn_term <- c()
nopagenedn_clinSig <- c()
ha <- hash()
ha1 <- hash()
.set(ha, keys = rr$Genes, values = rr$Terms)
.set(ha1, keys = rr$Genes, values = rr$ClinicalSignificance)
pg <- as.character(rr$Genes)
gen <- paste("^", pg, "$", sep = "")
for (jj in seq_len(length(dngene))) {
    dd1 <- grep(";", dngene[jj])
    gene_Terms_no <- c()
	dn_clinSig <- c()
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
				val4 <- hash::values(ha1, keys = gen3[l])
                dn_clinSig <- c(dn_clinSig, as.character(val4))
            }
            else {
                nopag <- c(nopag, "-")
				gene_Terms_no <- c(gene_Terms_no, "-")
			    dn_clinSig <- c(dn_clinSig, "-")
            }
        }
        nopag <- as.character(unique(nopag))
        nopagt <- as.character(unique(gene_Terms_no))
        nopagt <- paste(nopagt, collapse = ";")
		nopagdncs <- as.character(unique(dn_clinSig))
        nopagdncs <- paste(nopagdncs, collapse = ";")
        if (length(nopag) > 1) {
        nopagpa <- paste(nopag, collapse = ";")
        ff <- grep("-", nopagpa)
        if (length(ff) > 0) {
          nopagpa <- gsub("-", "", nopagpa)
          nopagpa <- gsub(";;", ";", nopagpa)
          nopagpa <- gsub("^;", ";", nopagpa)
          nopagenedn <- c(nopagenedn, nopagpa)
          nopagenedn_term <- c(nopagenedn_term, nopagt)
		  nopagenedn_clinSig <- c(nopagenedn_clinSig, nopagdncs)
        }
        else {
          nopagenedn <- c(nopagenedn, nopagpa)
          nopagenedn_term <- c(nopagenedn_term, nopagt)
		  nopagenedn_clinSig <- c(nopagenedn_clinSig, nopagdncs)
        }
      }
      else {
        nopagenedn <- c(nopagenedn, nopag)
        nopagenedn_term <- c(nopagenedn_term, nopagt)
		nopagenedn_clinSig <- c(nopagenedn_clinSig, nopagdncs)
      }
    }
    else {
      st1 <- st <- strsplit(dngene[jj], "\\(")
      dngene[jj] <- st1[[1]][1]
      pa <- paste("^", dngene[jj], "$", sep = "")
      gg <- grep(pa, gen, fixed = TRUE)
      if (length(gg) > 0) {
        nopagenedn <- c(nopagenedn, as.character(dngene[jj]))
        val2 <- hash::values(ha, keys = dngene[jj])
        nopagenedn_term <- c(nopagenedn_term, as.character(val2))
		val4 <- hash::values(ha1, keys = dngene[jj])
        nopagenedn_clinSig <- c(nopagenedn_clinSig, as.character(val4))
      }
      else {
        nopagenedn <- c(nopagenedn, "-")
        nopagenedn_term <- c(nopagenedn_term, "-")
		nopagenedn_clinSig <- c(nopagenedn_clinSig, "-")
      }
    }
  }
    dataPGDN <- data.frame(nopagenedn, nopagenedn_term, nopagenedn_clinSig)
    return(dataPGDN)
 }
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
#' @param outputType Variants in excel tabs or in different csv files.
#'        Options Excel or csv.
#' @param fileprefix Prefix to use for each of the files in the directory.
#' @param directoryName Directory name where individual SV files will be stored.
#' @param RZIPpath Character Path for the Rtools Zip package.
#' @return Excel file containing the annotated SV map, tabs divided based on
#' type of SVs.
#' @return individual csv files for each variant type if csv option is chosen.
#' @examples
#' \dontrun{
#' terms <- "Muscle Weakness"
#' gene <- gene_list_generation(
#'   method = "Single", term = terms,
#'   returnMethod_GeneList = "dataFrame"
#' )
#' smapName <- "F1.1_GM24385_DLE-1_P_trio_hg19.smap"
#' smappath <- system.file("extdata", smapName, package = "nanotatoR")
#' nanotatoR_main(smap, bed,
#'   inputfmtBed = c("BNBED"),
#'   n = 3, mergedFiles, buildSVInternalDB = TRUE, soloPath, solopattern,
#'   input_fmt_INF = c("dataframe"), returnMethod_GeneList = c("dataframe"),
#'   returnMethod_bedcomp = c("dataframe"), returnMethod_DGV = c("dataframe"),
#'   returnMethod_Internal = c("dataframe"), input_fmt_DGV = c("dataframe"),
#'   hgpath, smapName, method = c("Single"), term, thresh = 5,
#'   input_fmt_geneList = c("dataframe"), input_fmt_svMap = c("dataframe"),
#'   svData, dat_geneList, outpath = "", outputFilename = "", RZIPpath = "",
#'   outputType = c("Excel", "csv"), directoryName, fileprefix)
#' }
#' @import openxlsx
#' @import hash
#' @importFrom stats na.omit
#' @export


run_bionano_filter_Trio_DLE <- function(primaryGenesPresent = TRUE, 
                                input_fmt_geneList = c("Text", "dataFrame"),
                                input_fmt_svMap = c("Text", "dataFrame"),
                                SVFile = NULL, svData, dat_geneList, fileName, outpath,
                                outputFilename = "", RZIPpath,
		                        outputType = c("Excel", "csv"), 
								directoryName, fileprefix) {
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
        Non_Overlap_DN_Terms = as.character(dataPGDN$nopagenedn_term)),
		Non_Overlap_DN_ClinicalSig = as.character(dataPGDN$nopagenedn_clinSig))
    }else if (primaryGenesPresent == FALSE){
	    data <- r
	}else {stop("primaryGenesPresent Incorrect!!")}
  data$BNG_Freq_Perc_Filtered<-gsub("-",0,as.character(data$BNG_Freq_Perc_Filtered))
  data$BNG_Freq_Perc_UnFiltered<-gsub("-",0,as.character(data$BNG_Freq_Perc_UnFiltered))
  data$DGV_Freq_Perc<-as.numeric(data$DGV_Freq_Perc)
  data$Internal_Freq_Perc_Filtered<-as.numeric(data$Internal_Freq_Perc_Filtered)
  data$Internal_Freq_Perc_Unfiltered<-as.numeric(data$Internal_Freq_Perc_Unfiltered)
  data$BNG_Freq_Perc_Filtered<-as.numeric(data$BNG_Freq_Perc_Filtered)
  data$BNG_Freq_Perc_UnFiltered<-as.numeric(data$BNG_Freq_Perc_UnFiltered)
  data$DECIPHER_Frequency<-as.numeric(data$DECIPHER_Frequency)
 
  
  dat <- data[which(data$Type %in% "insertion"), ]
  dat1 <- data[which(data$Type %in% "deletion"), ]
  dat44 <- data[which(data$Type %in% "duplication"), ]
  dat45 <- data[which(data$Type %in% "duplication_split"), ]
  dat46 <- data[which(data$Type %in% "duplication_inverted"), ]
  dat_dup_rem<-rbind(dat44,dat45,dat46)
  dat_dup_rem_FINAL<-dat_dup_rem[which(((dat_dup_rem$Type %in% "duplication") |
     (dat_dup_rem$Type %in% "duplication_split") | (dat_dup_rem$Type %in% "duplication_inverted"))&
    (dat_dup_rem$Fail_assembly_chimeric_score == "pass")), ]
   
  dat3 <- rbind(dat1, dat, dat_dup_rem_FINAL)
  dat10 <- dat3[which((dat3$Found_in_self_molecules == "yes") &
    ((dat3$Found_in_parents_molecules == "none" | dat3$Found_in_parents_molecules == "-"))), ]
    
  dat11 <- dat3[which((dat3$Found_in_self_molecules == "yes") &
    ((dat3$Found_in_parents_molecules == "both"))),]
  dat12 <- dat3[which((dat3$Found_in_self_molecules == "yes") &
    ((dat3$Found_in_parents_molecules == "mother"))),]
  dat13 <- dat3[which((dat3$Found_in_self_molecules == "yes") &
    ((dat3$Found_in_parents_molecules == "father"))),]
  dat14 <- rbind(dat12, dat13)

  'gg <- grep("inversion", as.character(data$Type))
  dat8 <- data[gg, ]'
  dat8 <- data[which(data$Type %in% "inversion"
	        | data$Type %in% "inversion_partial"
			| data$Type %in% "inversion_paired"
			| data$Type %in% "inversion_repeat"), ]
  dat8 <-dat8[which((dat8$Found_in_self_molecules == "yes") 
        & (dat8$Fail_assembly_chimeric_score == "pass")), ]
  gg1 <- grep("translocation", as.character(data$Type))
  dat7 <- data[gg1, ]
  dat7 <-dat7[which((dat7$Found_in_self_molecules == "yes") 
        & (dat7$Fail_assembly_chimeric_score == "pass")), ]
  dat6 <- data[which(data$Type %in% "MisMatch"), ]
  ovrlapPG<-which(!((as.character(data$Overlap_PG)=="-")))
  datovrlapPG<-data[ovrlapPG,]
  nonovrlapupPG<-which(!((as.character(data$Non_Overlap_UP_PG)=="-")))
  datnonovrlapUPPG<-data[nonovrlapupPG,]
  nonovrlapdnPG<-which(!((as.character(data$Non_Overlap_DN_PG)=="-")))
  datnonovrlapDNPG<-data[nonovrlapdnPG,]
  datOvrLap <- rbind(datovrlapPG,datnonovrlapUPPG,datnonovrlapDNPG)
  
  'list_of_datasets <- list(
    "all" = data, "indel_dup_denovo" = dat10,
    "indel_dup_both" = dat11,
    "indel_dup_mother" = dat12, "indel_dup_father" = dat13,
    "indel_dup_cmpdHET" = dat14, "inv" = dat8, "trans" = dat7, "mismatch" = dat6,
    "all_PG_OV" = datovrlapPG, "all_PG_Non_OV_UP" = datnonovrlapUPPG,"all_PG_Non_OV_DN" = datnonovrlapDNPG
  )'
  if (outputType == "Excel"){
        list_of_datasets <- list(
        "indel_dup_denovo" = dat10,
        "indel_dup_both" = dat11, "indel_dup_cmpdHET" = dat14,
	    "inv" = dat8, "trans" = dat7,
        "indel_dup_mother" = dat12, "indel_dup_father" = dat13,
      	"all_PG_OV" = datOvrLap, "all" = data)
		fname <- paste(outputFilename, ".xlsx", sep = "")
        write.xlsx(list_of_datasets, file = file.path(outpath, fname), keepNA = TRUE)
	} else if (outputType == "csv"){
	    write.csv(dat10, file.path(directoryName, paste(prefix,"_indel_dup.csv",sep = ""), row.names = FALSE))
		write.csv(dat11, file.path(directoryName, paste(prefix,"_indel_dup_both.csv",sep = ""), row.names = FALSE))
		write.csv(dat14, file.path(directoryName, paste(prefix,"_indel_dup_cmpdHET.csv",sep = ""), row.names = FALSE))
		write.csv(dat8, file.path(directoryName, paste(prefix,"_inv.csv",sep = ""), row.names = FALSE))
		write.csv(dat7, file.path(directoryName, paste(prefix,"_trans.csv",sep = ""), row.names = FALSE))
		write.csv(dat12, file.path(directoryName, paste(prefix,"_indel_dup_mother.csv",sep = ""), row.names = FALSE))
		write.csv(dat13, file.path(directoryName, paste(prefix,"_indel_dup_father.csv",sep = ""), row.names = FALSE))
		write.csv(overlapPG, file.path(directoryName, paste(prefix,"_all_PG_OV.csv",sep = ""), row.names = FALSE))
		write.csv(data, file.path(directoryName, paste(prefix,"_all.csv",sep = ""), row.names = FALSE))
	} else {stop(" outputType incorrect !!")}
  


  # gg<-grep("inversion",as.character(data$Type))
  # dat8<-data[gg,]
  # gg1<-grep("translocation",as.character(data$Type))
  # dat7<-data[gg1,]
  # dat6<-data[which(data$Type %in% "MisMatch"),]
  # list_of_datasets <- list("all" = data, "indel" = dat3, "indel_denovo"= dat10, "indel_both"= dat11,
  #                        "indel_mother"= dat12, "indel father"= dat13, "indel_cmpdHET"=dat14, "inv"= dat8, "trans"=dat7, "mismatch"=dat6,
  #                       "all_PG_OV"=data, "all_PG_Non_OV"= data)
  # st1<-strsplit(SVFile,".txt")[[1]][1]
  'fname <- paste(outputFilename, "_DLE.xlsx", sep = "")
  write.xlsx(list_of_datasets, file = file.path(outpath, fname),keepNA=TRUE)'
}


### Running the code
### Copy this 1 line at a time without the hash (#)
### Warning1: Don't remove the # from before the command from filterread.r
### Warning2: give the fullpath with filename for the  R code
# source("Z:\\bionano\\Test_R_Codes\\filterread.r")##path of where the Rcode is
#
# run_bionano_filter(path,SVFile,fileName)
