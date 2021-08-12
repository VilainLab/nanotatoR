#' Extracting terms for genes that overlap SVs
#'
#' @param rr character. dataframe with primary genes and terms
#'      associated.
#' @param ogene  character. genes that overlap the SV.
#' @return Dataframe with overlapping genes and terms.
#' @examples
#' terms= c("steroid_Gene","steroid synthesis_Gene")
#' genes <- c("NR1H3", "ABCC4")
#' rr <- data.frame(Genes = genes, Terms = terms)
#' smapName="GM24385_Ason_DLE1_VAP_trio5.smap"
#' smap = system.file("extdata", smapName, package="nanotatoR")
#' bedFile <- system.file("extdata", "HomoSapienGRCH19_lift37.bed", package="nanotatoR")
#' outpath <- system.file("extdata", package="nanotatoR")
#' datcomp<-overlapnearestgeneSearch(smap = smap, 
#' bed=bedFile, inputfmtBed = "bed", outpath, 
#' n = 3, returnMethod_bedcomp = c("dataFrame"), 
#' input_fmt_SV = "Text",
#' EnzymeType = "SE", 
#' bperrorindel = 3000, bperrorinvtrans = 50000)
#' ogene <- as.character(datcomp$OverlapGenes_strand_perc)
#' datogenes <- overlappingGenes (rr, ogene)
#' @import hash
#' @importFrom stats na.omit
#' @export
overlappingGenes <- function(rr, ogene){
ha <- hash()
ha1 <- hash()
    if(ncol(rr) > 1){
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
                opagpa <- gsub(";$", "", opagpa)
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
    }else{
        #.set(ha, keys = rr$Genes, values = rr$Terms)
        #.set(ha1, keys = rr$Genes, values = rr$ClinicalSignificance)
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
                        #val1 <- hash::values(ha, keys = gen2[l])
                        gene_Terms <- c(gene_Terms, "-")
                      #val2 <- hash::values(ha1, keys = gen2[l])
                        clinSig <- c(clinSig, "-")
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
                opagpa <- gsub(";$", "", opagpa)
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
                    #val1 <- hash::values(ha, keys = ogene[ii])
                    pagene_term <- c(pagene_term, "-")
                #val2 <- hash::values(ha1, keys = ogene[ii])
                    pagene_clinSig <- c(pagene_clinSig, "-")
                }
                else {
                    pagene <- c(pagene, "-")
                    pagene_term <- c(pagene_term, "-")
                pagene_clinSig <- c(pagene_clinSig, "-")
                }
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
#' @param upgene  character. genes that overlap the SV.
#' @return Dataframe with overlapping genes and terms.
#' @examples
#' terms= c("steroid_Gene","steroid synthesis_Gene")
#' genes <- c("NR1H3", "ABCC4")
#' rr <- data.frame(Genes = genes, Terms = terms)
#' smapName="GM24385_Ason_DLE1_VAP_trio5.smap"
#' smap = system.file("extdata", smapName, package="nanotatoR")
#' bedFile <- system.file("extdata", "HomoSapienGRCH19_lift37.bed", package="nanotatoR")
#' outpath <- system.file("extdata", package="nanotatoR")
#' datcomp<-overlapnearestgeneSearch(smap = smap, 
#' bed=bedFile, inputfmtBed = "bed", outpath, 
#' n = 3, returnMethod_bedcomp = c("dataFrame"), 
#' input_fmt_SV = "Text",
#' EnzymeType = "SE", 
#' bperrorindel = 3000, bperrorinvtrans = 50000)
#' upgene <- as.character(datcomp$Upstream_nonOverlapGenes_dist_kb)
#' dataPGUP <- nonOverlappingUPGenes (rr, upgene)
#' @import hash
#' @importFrom stats na.omit
#' @export
nonOverlappingUPGenes <- function(rr, upgene){
    ha <- hash()
    ha1 <- hash()
    nopageneup <- c()
    nopageneup_term <- c()
    nopageneup_clinSig <- c()
    if(ncol(rr) > 1){
        .set(ha, keys = rr$Genes, values = rr$Terms)
        .set(ha1, keys = rr$Genes, values = rr$ClinicalSignificance)
        pg <- as.character(rr$Genes)
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
                         nopagpa <- gsub(";$", ";", nopagpa)
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
    }
    else{
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
                        #val2 <- hash::values(ha, keys = gen3[l])
                        gene_Terms_no <- c(gene_Terms_no, "-")
                #val4 <- hash::values(ha1, keys = gen3[l])
                        up_clinSig <- c(up_clinSig, "-")
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
                         nopagpa <- gsub(";$", ";", nopagpa)
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
                            #val2 <- hash::values(ha, keys = upgene[jj])
                            nopageneup_term <- c(nopageneup_term, "-")
                        val4 <- hash::values(ha1, keys = upgene[jj])
                            nopageneup_clinSig <- c(nopageneup_clinSig, "-")
                        }
                        else {
                            nopageneup <- c(nopageneup, "-")
                            nopageneup_term <- c(nopageneup_term, "-")
                        nopageneup_clinSig <- c(nopageneup_clinSig, "-")
                        }
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
#' @param dngene  character. genes that overlap the SV.
#' @return Dataframe with overlapping genes and terms.
#' @examples
#' terms= c("steroid_Gene","steroid synthesis_Gene")
#' genes <- c("NR1H3", "ABCC4")
#' rr <- data.frame(Genes = genes, Terms = terms)
#' smapName="GM24385_Ason_DLE1_VAP_trio5.smap"
#' smap = system.file("extdata", smapName, package="nanotatoR")
#' bedFile <- system.file("extdata", "HomoSapienGRCH19_lift37.bed", package="nanotatoR")
#' outpath <- system.file("extdata", package="nanotatoR")
#' datcomp<-overlapnearestgeneSearch(smap = smap, 
#' bed=bedFile, inputfmtBed = "bed", outpath, 
#' n = 3, returnMethod_bedcomp = c("dataFrame"), 
#' input_fmt_SV = "Text",
#' EnzymeType = "SE", 
#' bperrorindel = 3000, bperrorinvtrans = 50000)
#' dngene <- as.character(datcomp$Downstream_nonOverlapGenes_dist_kb)
#' dataPGDN <- nonOverlappingDNGenes (rr, dngene)
#' @import hash
#' @importFrom stats na.omit
#' @export

nonOverlappingDNGenes <- function(rr, dngene){
    nopagenedn <- c()
    nopagenedn_term <- c()
    nopagenedn_clinSig <- c()
    ha <- hash()
    ha1 <- hash()
    if(ncol(rr) > 1){
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
              nopagpa <- gsub(";$", ";", nopagpa)
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
    }else{
        #.set(ha, keys = rr$Genes, values = rr$Terms)
        #.set(ha1, keys = rr$Genes, values = rr$ClinicalSignificance)
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
                        #val2 <- hash::values(ha, keys = gen3[l])
                        gene_Terms_no <- c(gene_Terms_no, "-")
                #val4 <- hash::values(ha1, keys = gen3[l])
                        dn_clinSig <- c(dn_clinSig, "-")
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
              nopagpa <- gsub(";$", ";", nopagpa)
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
                #val2 <- hash::values(ha, keys = dngene[jj])
                nopagenedn_term <- c(nopagenedn_term, "-")
            #val4 <- hash::values(ha1, keys = dngene[jj])
                nopagenedn_clinSig <- c(nopagenedn_clinSig, "-")
              }
              else {
                nopagenedn <- c(nopagenedn, "-")
                nopagenedn_term <- c(nopagenedn_term, "-")
            nopagenedn_clinSig <- c(nopagenedn_clinSig, "-")
              }
            }
        }
    }
    dataPGDN <- data.frame(nopagenedn, nopagenedn_term, nopagenedn_clinSig)
    return(dataPGDN)
 }
