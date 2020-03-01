#' Extracting genes for phenotype/diseases from NCBI.
#'
#' @param method_entrez  character. Input Method for terms. Choices are 
#'                "Single","Multiple" and "Text".
#' @param termPath  character. Path and file name for textfile.
#' @param term  character. Single or Multiple Terms.
#' @param outpath character. Path where gene lists are saved.
<<<<<<< HEAD
#' @param omim character. omim2gene file name and location.
#' @param clinvar character. clinvar file name and location.
#' @param gtr character. gtr file name and location.
#' @param removeClinvar logical. Deletes the Clinvar database if TRUE.
#' @param removeGTR logical. Deletes the GTR database if TRUE.
#' @param downloadClinvar logical. Downloads the Clinvar database if TRUE.
#' @param downloadGTR logical. Downloads the GTR database if TRUE.
#' @param url_gtr character. url for GTR.
#' @param url_clinvar character. url for clinvar.
=======
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
#' @param thresh integer. Threshold for the number of terms sent to entrez.
#'                Note if large lists are sent to ncbi, it might fail to get
#'                processed. Default is 5.
#' @param returnMethod Method of returning output. Options, Text or data.frame.
#' @return Text files containg gene list and terms associated with them
#'         are stored as text files.
#' @examples
#' terms="Muscle Weakness"
<<<<<<< HEAD
#' genes <- gene_list_generation(
#'      method_entrez = c("Single"), 
#'      term = terms,
#'      returnMethod=c("dataFrame"), 
#'		omim = "Y:/Suro/nanotatoRDatabases/mim2gene.txt", 
#'		downloadClinvar = TRUE, 
#'		gtr = "Y:/Suro/nanotatoRDatabases/gtr_07162019.txt")
=======
#' genes <- gene_list_generation(method="Single", term=terms, returnMethod="dataFrame")
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
#' @import stats
#' @import rentrez 
#' @import utils
#' @import httr
<<<<<<< HEAD
#' @import tidyverse
=======
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
#' @export
gene_list_generation<-function(method_entrez = c("Single","Multiple","Text"), 
        termPath, 
        term, outpath, thresh=5, 
<<<<<<< HEAD
        returnMethod=c("Text","dataFrame"), omim, clinvar, gtr,
        removeClinvar = FALSE, removeGTR = FALSE, 
		downloadClinvar = FALSE, downloadGTR = FALSE,
		url_gtr){
    
=======
        returnMethod=c("Text","dataFrame")){
        
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
    #setwd(path)
    thresh=5
    ##Checking for whether the input method is commandline (Single or Multiple)
    ## or from Text
    if (method_entrez == "Single"){
        terms<-term
    }
    else if (method_entrez == "Multiple"){
        terms<- as.character(term)
<<<<<<< HEAD
       
    }
    else if (method_entrez == "Text"){
        r<-read.csv(termPath)
        print(termPath)   
        terms<-as.character(r$Terms)
=======
        terms<- as.character(term) 
    }
    else if (method_entrez == "Text"){
        r<-read.csv(termPath)
        terms<-as.character(r$Term)
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
    }
    else{
        stop("method_entrez of Input Incorrect!!!")
    }
    #ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
    ##Checking whether the number of terms is above the threshold
    # and then calls the appropriate functions for getting gene lists.
    ##Note if the term number is more than 40, supply terms in batches 
    # of 15 or 20
    if(length(terms)>=thresh){
        ##Diving the number of genes to be given as an
        ##input based on the threshold.
<<<<<<< HEAD
        terms_1<-as.character(terms[1:thresh])
        terms_2<-as.character(terms[thresh:length(terms)])
        ##Extracting genes from the databases
        g<-gene_extraction(terms_1)
        g_1<-gene_extraction(terms_2)
        print(paste0("main",omim = omim))
        g_o<-omim_gene(terms_1, omim =omim)
        g_o_1<-omim_gene(terms_2, omim =omim)
        print(paste0("main",gtr = gtr))
        g_g<-gtr_gene(terms_1, 
		        gtr = gtr, 
				url_gtr = url_gtr, 
				downloadGTR = downloadGTR)
        g_g_1<-gtr_gene(terms_2, 
		    gtr = gtr,
			url_gtr = url_gtr,
			downloadGTR = downloadGTR)
        #print(paste0("main",clinvar = clinvar))
		if(downloadClinvar == FALSE){
        g_c<-clinvar_gene(terms_1, 
		    clinvar = clinvar, 
			downloadClinvar = downloadClinvar)
        g_c_1<-clinvar_gene(terms_2, 
		    clinvar = clinvar, 
			downloadClinvar = downloadClinvar)
		}else{
		print(terms_1)
		g_c<-clinvar_gene(terms_1, 
		    downloadClinvar = downloadClinvar)
        g_c_1<-clinvar_gene(terms_2, 
		    downloadClinvar = downloadClinvar)
		}
=======
        terms_1<-as.character(terms[seq(thresh)])
        terms_2<-as.character(terms[thresh:seq_len(length(terms))])
        ##Extracting genes from the databases
        g<-gene_extraction(terms_1)
        g_1<-gene_extraction(terms_2)
        g_o<-omim_gene(terms_1)
        g_o_1<-omim_gene(terms_2)
        g_g<-gtr_gene(terms_1)
        g_g_1<-gtr_gene(terms_2)
        g_c<-clinvar_gene(terms_1)
        g_c_1<-clinvar_gene(terms_2)
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
        ##Collating output from all the datasets
        Genes<-c(as.character(g$geneName),as.character(g_1$geneName),
            as.character(g_o$omimGenes),as.character(g_o_1$omimGenes),
            as.character(g_g$gtrGenes),as.character(g_g_1$gtrGenes),
            as.character(g_c$clinvarGenes),
            as.character(g_c_1$clinvarGenes))
        Terms<-c(as.character(g$Final_terms),as.character(g_1$Final_terms),
            as.character(g_o$Final_terms_OMIM),
            as.character(g_o_1$Final_terms_OMIM),
            as.character(g_g$Final_terms_GTR),
            as.character(g_g_1$Final_terms_GTR),
            as.character(g_c$Final_terms_Clinvar),
            as.character(g_c_1$Final_terms_Clinvar))
<<<<<<< HEAD
		ClinicalSig <- c(as.character(rep("-",length(g$Final_terms))),
		    as.character(rep("-",length(g_1$Final_terms))),
			as.character(rep("-",length(g_o$Final_terms_OMIM))),
			as.character(rep("-",length(g_o_1$Final_terms_OMIM))),
			as.character(rep("-",length(g_g$Final_terms_GTR))),
			as.character(rep("-",length(g_g_1$Final_terms_GTR))),
			as.character(g_c$clinicalSig),
			as.character(g_c_1$clinicalSig))
    ##Writing genes and terms in a dataset
    dat<-data.frame(Genes, Terms, ClinicalSig)
    }
    else{
        ##Extracting genes from the databases
        g<-gene_extraction(terms = terms)
        #print(dim(g))
        g_o<-omim_gene(terms = terms, omim = omim)
		if(downloadGTR == FALSE){
		    g_g<-gtr_gene(terms = terms, 
		        gtr = gtr,
			    downloadGTR = FALSE)
		}
		else{
		    g_g<-gtr_gene(terms = terms, 
		        url_gtr = url_gtr,
			    downloadGTR = TRUE)
		}
		if(downloadClinvar == FALSE){
		    g_c<-clinvar_gene(terms = terms, 
		    clinvar = clinvar,
			downloadClinvar = FALSE)
		}
		else{
		    g_c<-clinvar_gene(terms = terms, 
		    downloadClinvar = downloadClinvar)
		}
        
        
=======
    ##Writing genes and terms in a dataset
    dat<-data.frame(Genes,Terms)
    }
    else{
        ##Extracting genes from the databases
        g<-gene_extraction(terms)
        #print(dim(g))
        g_o<-omim_gene(terms)
        g_g<-gtr_gene(terms)
        g_c<-clinvar_gene(terms)
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
        ##Collating output from all the datasets
        Genes<-c(as.character(g$geneName),as.character(g_g$gtrGenes),
            as.character(g_o$omimGenes),as.character(g_c$clinvarGenes))
        Terms<-c(as.character(g$Final_terms),
            as.character(g_g$Final_terms_GTR),
            as.character(g_o$Final_terms_OMIM),
            as.character(g_c$Final_terms_Clinvar))
<<<<<<< HEAD
		ClinicalSig <- c(as.character(rep("-",length(g$Final_terms))),
		    as.character(rep("-",length(g_o$Final_terms_OMIM))),
			as.character(rep("-",length(g_g$Final_terms_GTR))),
			as.character(g_c$clinicalSig))
		##Writing genes and terms in a dataset
        dat<-data.frame(Genes, Terms, ClinicalSig)
=======
        ##Writing genes and terms in a dataset
        dat<-data.frame(Genes,Terms)
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
    }
    ##We write out all the genes and terms from the 
    ##different databases into a$description
    ##single dataframe
    gene1<-c()
    term1<-c()
<<<<<<< HEAD
	csig <- c()
	if(nrow(dat) > 0){
	    Genes <- as.character(dat$Genes)
	    Terms <- as.character(dat$Terms)
	    clinsg <- as.character(dat$ClinicalSig)
        uqGenes<-as.character(unique(na.exclude(Genes)))
        for(ii in 1:length(uqGenes)){
            idx<-which(dat$Genes %in% uqGenes[ii])
            if(length(idx)>1){
                gene1<-c(gene1,as.character(uqGenes[ii]))
                tt<-as.character(Terms[idx])
	    		pa<-paste(tt,collapse=",")
	    		paF <- paste(uqGenes[ii], "(", pa, ")", sep ="")
                term1<-c(term1,as.character(paF))
	    		nn <- as.character(clinsg[idx])
	    		if(length(unique(nn)) == 1 
	    		    & length(grep("-",unique(nn)) >= 1)){
	    		    csig <- c(csig,"-")
	    		}else{
	    		    
	    			pa1 <-paste(nn,collapse=",")
	    		    pa1F <- paste(uqGenes[ii], "(", pa1, ")", sep ="")
	    		    csig <- c(csig,as.character(pa1F))
	    		}
            }
            else{
                gene1<-c(gene1,as.character(uqGenes[ii]))
                paF <- paste(uqGenes[ii], "(", Terms[idx], ")", sep ="")
	    		term1<-c(term1,as.character(paF))
	    		nn <- as.character(clinsg[idx])
	    		if(length(unique(nn)) == 1 
	    		    & length(grep("-",unique(nn)) >= 1)){
	    		    csig <- c(csig,"-")
	    		}else{
	    		pa1F <- paste(uqGenes[ii], "(", csig[idx], ")", sep ="")
	    		csig <- c(csig,as.character(pa1F))
	    		}
	    
            }
        }
	}else { 
	    gene1 <- c(); 
	    Terms <- c(); 
	    ClinicalSignificance <- c()
	}
    ##Write the final dataframe
    dat_Final<-data.frame(Genes=gene1,Terms=term1, ClinicalSignificance = csig)
=======
    uqGenes<-as.character(unique(na.exclude(Genes)))
    for(ii in seq_along(uqGenes)){
        idx<-which(dat$Genes %in% uqGenes[ii])
        if(length(idx)>1){
            gene1<-c(gene1,as.character(uqGenes[ii]))
            tt<-as.character(Terms[idx])
            pa<-paste(tt,collapse=",")
            term1<-c(term1,as.character(pa))
        }
        else{
            gene1<-c(gene1,as.character(uqGenes[ii]))
            term1<-c(term1,as.character(Terms[idx]))
 
        }
    }
    ##Write the final dataframe
    dat_Final<-data.frame(Genes=gene1,Terms=term1)
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06

    #final_genes<-unique(final_genes)
    ##Make the filename and write the data into a text file
    if(returnMethod=="Text"){
        st<-strsplit(term,".csv")
        fname<-st[[1]][1]
        rownames(dat_Final)<-NULL
        fileName<-paste(fname,"_GeneList.txt",sep="")
        write.table(dat_Final,file.path(outpath,fileName),
        col.names = c("Genes","Terms"),row.names=FALSE,sep="\t")
    }
    else if (returnMethod=="dataFrame"){
        dat_Final<-dat_Final
    }
    else{
        stop ("Return Method Incorrect")
    }
<<<<<<< HEAD
    if (removeClinvar == TRUE){
        file.remove(clinvar)
    }else{
        print("keeping the Clinvar files")
    }
    if (removeGTR == TRUE){
        file.remove(gtr)
    }else{
        print("keeping the GTR files")
    }
=======
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
    return(dat_Final)
}
#' Extracting genes from gene database NCBI.
#'
#' @param terms  Single or Multiple Terms.
#' @return Dataframe returned containing gene lists in entrezid and Gene 
#' Symbols, and terms associated with it
#' @examples
#' terms="Muscle Weakness"
#' ge <- gene_extraction(terms)
#' @import rentrez utils
#' @importFrom stats na.omit 
#' @import httr
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#' @export
gene_extraction<-function(terms){
<<<<<<< HEAD
    ##Initializing values
    geneName<-c()
    Final_terms <-c()
    httr::set_config(httr::config(http_version = 0))
    ##Initialising data to be extracted from ensembl
    'ensembl = useMart("ensembl",host = "www.ensembl.org", 
    ensemblRedirect = FALSE,dataset="hsapiens_gene_ensembl")'
    
    ##Checking for term size and extracting gene list accordingly
    if (length(terms)>1){ #checking for terms greater than 1
        for(ll in seq_len(length(terms))){
            print(terms[ll])
            ##Searching Entrez with the terms
            print(paste0("Gene:", ll))
            patl <- paste(terms[ll])
            b<-entrez_search(db="gene",term=patl,retmax=99999999)
            ##Extracting gene IDs, and checking if length of the gene ids
            ##greater than 0.
            geneID<-b$ids
            gg1<-data.frame()
            ##Checking if geneids greater than 0, if not moved to the next term
            if(length(geneID)>0){
                #print(paste("GeneID length:",length(geneID),sep=""))
                ##uid<-c();geneNam<-c();desc<-c();status<-c()
                #geneNam<-c()
                
                #for(ii in 1:length(geneID)){
                # a<-entrez_summary(db="gene", id=geneID[ii])
                # ##uid<-c(uid,as.character(a$uid))
                # geneNam<-c(geneNam,toupper(a$name))
                # 
                # ##desc<-c(desc,a$description)
                # ##status<-c(status,a$status)
                # }
                ##extractingthe gene Symbols associated with the gene id from 
                ##Biomart.
                gn1 <- tryCatch(
                    mapIds(org.Hs.eg.db, geneID, "SYMBOL", "ENTREZID"),
                    warning = function(w) {
                    print("warning")
                    return(NA)}, 
                    error = function(e) {
                    print("Error")
                    return(NA)}
                )
                gn2 = tryCatch(    
                mapIds(org.Hs.eg.db, geneID, "ENSEMBL", "ENTREZID"),
                warning = function(w) {
                    print("warning")
                    return(NA)}, 
                    error = function(e) {
                    print("Error")
                    return(NA)}
                )
                if (length(unique(as.character(na.omit(gn1)))) ==0) {next} else {gn1 = gn1}
                if (length(unique(as.character(na.omit(gn2)))) ==0) {next} else {gn2 = gn2}
                rn <- row.names(data.frame(gn1))
                rn1 <- row.names(data.frame(gn2))
                gene1<-data.frame(
                    entrezID = as.numeric(rn),
                    ensemblID = as.character(data.frame(gn2)[,1]),
                    hgnc_symbol = as.character(data.frame(gn1)[,1])
                    )
                gn<-as.character(unique(gene1$hgnc_symbol))
                gn_1_temp <- gn[gn != ""]
                gn_1 <- as.character(na.omit(gn_1_temp))
                geneName<-c(geneName,as.character(gn_1))
                terms_list<-as.character(rep(paste(terms[ll],"_Gene",sep=""),
                    length(gn_1)))
                Final_terms<-c(Final_terms,terms_list)
            }
            else{next}
        }   
    }
    else if (length(terms)==1){ 
        #checking for terms equal to 1
        ##Searching in the entrez database
        patl <- paste(terms)
        b<-entrez_search(db="gene", term=patl, retmax=99999999)
=======
##Initializing values
geneName<-c()
Final_terms <-c()
httr::set_config(httr::config(http_version = 0))
##Initialising data to be extracted from ensembl
'ensembl = useMart("ensembl",host = "www.ensembl.org", 
ensemblRedirect = FALSE,dataset="hsapiens_gene_ensembl")'

##Checking for term size and extracting gene list accordingly
if (length(terms)>1){ #checking for terms greater than 1
    for(ll in seq_len(length(terms))){
        ##Searching Entrez with the terms
        b<-entrez_search(db="gene",term=terms[ll],retmax=99999999)
        ##Extracting gene IDs, and checking if length of the gene ids
        ##greater than 0.
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
        geneID<-b$ids
        gg1<-data.frame()
        ##Checking if geneids greater than 0, if not moved to the next term
        if(length(geneID)>0){
<<<<<<< HEAD
            print(paste("GeneID length:",length(geneID),sep=""))
            ##uid<-c();geneNam<-c();desc<-c();status<-c()
            #geneNam<-c()
=======
            #print(paste("GeneID length:",length(geneID),sep=""))
            ##uid<-c();geneNam<-c();desc<-c();status<-c()
            #geneNam<-c()
            
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
            #for(ii in 1:length(geneID)){
            # a<-entrez_summary(db="gene", id=geneID[ii])
            # ##uid<-c(uid,as.character(a$uid))
            # geneNam<-c(geneNam,toupper(a$name))
            # 
            # ##desc<-c(desc,a$description)
            # ##status<-c(status,a$status)
            # }
<<<<<<< HEAD
            ##Extracting gene symbols using BiMart
            gn1 <- mapIds(org.Hs.eg.db, geneID, "SYMBOL", "ENTREZID")
            gn2 <- mapIds(org.Hs.eg.db, geneID, "ENSEMBL", "ENTREZID")
            rn <- row.names(data.frame(gn2))
=======
            ##extractingthe gene Symbols associated with the gene id from 
            ##Biomart.
            gn1 <- mapIds(org.Hs.eg.db, geneID, "SYMBOL", "ENTREZID")
            gn2 <- mapIds(org.Hs.eg.db, geneID, "ENSEMBL", "ENTREZID")
            rn <- row.names(data.frame(gn1))
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
            rn1 <- row.names(data.frame(gn2))
            gene1<-data.frame(
                entrezID = as.numeric(rn),
                ensemblID = as.character(data.frame(gn2)[,1]),
                hgnc_symbol = as.character(data.frame(gn1)[,1])
                )
            gn<-as.character(unique(gene1$hgnc_symbol))
            gn_1_temp <- gn[gn != ""]
            gn_1 <- as.character(na.omit(gn_1_temp))
            geneName<-c(geneName,as.character(gn_1))
<<<<<<< HEAD
            terms_list<-as.character(rep(paste(terms,"_Gene",sep=""),
                length(gn_1)))
            Final_terms<-c(Final_terms,terms_list)
        }
        else{
            print ("No genes for term !!")
        }
    }
    else{
        stop("Term length is zero")
    }
    
    ##Data written to the dataframe and returned
    dat1<-data.frame(geneName,Final_terms)
    #Finaldata<-c()
    #print(paste("GeneName length:",length(geneName),sep=""))
    return(dat1)
}
#' Extracting genes from OMIM database NCBI.
#'
#' @param terms character Single or Multiple Terms.
#' @param omim character omim database location.
=======
            terms_list<-as.character(rep(paste(terms[ll],"_Gene",sep=""),
                length(gn_1)))
            Final_terms<-c(Final_terms,terms_list)
        }
        else{next}
    }   
}
else if (length(terms)==1){ 
    #checking for terms equal to 1
    ##Searching in the entrez database
    b<-entrez_search(db="gene",term=terms,retmax=99999999)
    geneID<-b$ids
    gg1<-data.frame()
    ##Checking if geneids greater than 0, if not moved to the next term
    if(length(geneID)>0){
        print(paste("GeneID length:",length(geneID),sep=""))
        ##uid<-c();geneNam<-c();desc<-c();status<-c()
        #geneNam<-c()
        #for(ii in 1:length(geneID)){
        # a<-entrez_summary(db="gene", id=geneID[ii])
        # ##uid<-c(uid,as.character(a$uid))
        # geneNam<-c(geneNam,toupper(a$name))
        # 
        # ##desc<-c(desc,a$description)
        # ##status<-c(status,a$status)
        # }
        ##Extracting gene symbols using BiMart
        gn1 <- mapIds(org.Hs.eg.db, geneID, "SYMBOL", "ENTREZID")
            gn2 <- mapIds(org.Hs.eg.db, geneID, "ENSEMBL", "ENTREZID")
            rn <- row.names(data.frame(gn1))
            rn1 <- row.names(data.frame(gn2))
            gene1<-data.frame(
                entrezID = as.numeric(rn),
                ensemblID = as.character(data.frame(gn2)[,1]),
                hgnc_symbol = as.character(data.frame(gn1)[,1])
                )
            gn<-as.character(unique(gene1$hgnc_symbol))
            gn_1_temp <- gn[gn != ""]
            gn_1 <- as.character(na.omit(gn_1_temp))
            geneName<-c(geneName,as.character(gn_1))
            terms_list<-as.character(rep(paste(terms,"_Gene",sep=""),
                length(gn_1)))
            Final_terms<-c(Final_terms,terms_list)
    }
    else{
        print ("No genes for term !!")
    }
}
else{
stop("Term length is zero")
}

##Data written to the dataframe and returned
dat1<-data.frame(geneName,Final_terms)
#Finaldata<-c()
#print(paste("GeneName length:",length(geneName),sep=""))
return(dat1)
}
#' Extracting genes from OMIM database NCBI.
#'
#' @param terms  Single or Multiple Terms.
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
#' @return Dataframe returned containing gene lists in entrezid and Gene 
#' Symbols, and terms associated with it
#' @examples
#' terms="Muscle Weakness"
<<<<<<< HEAD
#' ge <- omim_gene(terms, omim = "Y:/Suro/nanotatoRDatabases/mim2gene.txt")
#' @import rentrez 
#' @import utils
=======
#' ge <- omim_gene(terms)
#' @import rentrez utils
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#' @importFrom stats na.omit 
#' @import httr
#' @export
<<<<<<< HEAD
omim_gene<-function(terms, omim){
=======
omim_gene<-function(terms){
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
    ##Initialising variables
    omimGenes<-c()
    Final_terms_OMIM<-c()
    httr::set_config(httr::config(http_version = 0))
    ##Initialising data to be extracted from ensembl
    ####Checking for term size and extracting gene list accordingly
<<<<<<< HEAD
    print(paste0("omim_gene",omim))
    datOmim <- reading_mim2gene(omim = omim)
    if (length(terms)>1){
        ##Length of the terms input greater than 1.
        for(ll in seq_len(length(terms))){
            print(paste0("OMIM:", ll))
            ## Searching the OMIM database for terms and Gene IDs extracted
            #patl <- paste(terms[ll],"[DIS]")
            'c<-entrez_search(db="omim",term=paste(terms[ll],"[WORD]",sep=""),
                    retmax=99999999)'
            d<-entrez_search(db="omim",term=paste(terms[ll],"[DIS]",sep=""),
                    retmax=99999999)
            geneName <- c()
            #omimgeneID_word<-c$ids
            omimgeneID_dis<-d$ids
            final_omim<-unique(as.character(omimgeneID_dis))
=======
    if (length(terms)>1){
        ##Length of the terms input greater than 1.
        for(ll in seq_len(length(terms))){
            ## Searching the OMIM database for terms and Gene IDs extracted
            c<-entrez_search(db="omim",term=paste(terms[ll],"[WORD]",sep=""),
                    retmax=99999999)
            d<-entrez_search(db="omim",term=paste(terms[ll],"[DIS]",sep=""),
                    retmax=99999999)
            omimgeneID_word<-c$ids
            omimgeneID_dis<-d$ids
            final_omim<-unique(c(as.character(omimgeneID_word),
                            as.character(omimgeneID_dis)))
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
            ##Checking if geneids greater than 0, if not moved to the next term
            if(length(final_omim)>0){
                #omimgen<-c()
                #for (kk in 1:length(final_omim)){
                #k<-entrez_summary(db="omim", id=final_omim[kk])
                #gene<-as.character(k$title)
                #st<-strsplit(gene,split=";")
                #omimgen<-c(omimgen,as.character(toupper(st[[1]][2])))
                #}
<<<<<<< HEAD
                ###Extracting OMIM data
                omimid <- datOmim[,1]
                geneOmim <- datOmim[,2]
                pas <- paste0("^",omimid,"$")
                pasq <- paste0("^",final_omim,"$")
                geneID <- c()
                for(ii in 1:length(pasq)){
                    g1 <- grep(pasq[ii], pas, fixed = TRUE)
                    geneID <- c(geneID, as.character(geneOmim[g1]))
                }
                geneID <- na.omit(as.character(geneID))
                ##Get the Human gene symbols for the extracted Gene IDs
            gn1 <- tryCatch(
                mapIds(org.Hs.eg.db, geneID, "SYMBOL", "ENTREZID"),
                warning = function(w) {
                print("warning")
                return(NA)}, 
                error = function(e) {
                print("Error")
                return(NA)}
            )
            gn2 = tryCatch(    
            mapIds(org.Hs.eg.db, geneID, "ENSEMBL", "ENTREZID"),
            warning = function(w) {
                print("warning")
                return(NA)}, 
                error = function(e) {
                print("Error")
                return(NA)}
            )
            if (length(unique(as.character(na.omit(gn1)))) ==0) {next} else {gn1 = gn1}
            if (length(unique(as.character(na.omit(gn2)))) ==0) {next} else {gn2 = gn2}
                'gn1 <- mapIds(org.Hs.eg.db, geneID, "SYMBOL", "ENTREZID")
                gn2 <- mapIds(org.Hs.eg.db, geneID, "ENSEMBL", "ENTREZID")'
=======
                ##Get the Human gene symbols for the extracted Gene IDs
                gn1 <- mapIds(org.Hs.eg.db, final_omim, "SYMBOL", "ENTREZID")
                gn2 <- mapIds(org.Hs.eg.db, final_omim, "ENSEMBL", "ENTREZID")
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
                rn <- row.names(data.frame(gn1))
                rn1 <- row.names(data.frame(gn2))
                gene1<-data.frame(
                entrezID = as.numeric(rn),
                ensemblID = as.character(data.frame(gn2)[,1]),
                hgnc_symbol = as.character(data.frame(gn1)[,1])
                )
                gn<-as.character(unique(gene1$hgnc_symbol))
                gn_1_temp <- gn[gn != ""]
                gn_1 <- as.character(na.omit(gn_1_temp))
                omimGenes<-c(omimGenes,as.character(gn_1))
<<<<<<< HEAD
                terms_list<-as.character(rep(paste(terms[ll],"_OMIM",sep=""),
=======
                terms_list<-as.character(rep(paste(terms[ll],"_Gene",sep=""),
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
                    length(gn_1)))
                Final_terms_OMIM<-c(Final_terms_OMIM,terms_list)
            }
            else{
                next
            }        
        }
    }
    else if (length(terms)==1){
        ##Extracting data from OMIM
<<<<<<< HEAD
        'c<-entrez_search(db="omim",term=paste(terms,"[WORD]",sep=""),
            retmax=99999999)'
        d <- entrez_search(db="omim",term=paste(terms,"[DIS]",sep=""),
                retmax=99999999)
        #omimgeneID_word<-c$ids
        omimgeneID_dis<-d$ids
        final_omim<-unique(c(as.character(omimgeneID_dis)))
=======
        c<-entrez_search(db="omim",term=paste(terms,"[WORD]",sep=""),
            retmax=99999999)
        d<-entrez_search(db="omim",term=paste(terms,"[DIS]",sep=""),
                retmax=99999999)
        omimgeneID_word<-c$ids
        omimgeneID_dis<-d$ids
        final_omim<-unique(c(as.character(omimgeneID_word),as.character(
            omimgeneID_dis)))
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
        ##Checking if geneids greater than 0, if not move to the next term
        if(length(final_omim)>0){
            #omimgen<-c()
            #for (kk in 1:length(final_omim)){
            #k<-entrez_summary(db="omim", id=final_omim[kk])
            #gene<-as.character(k$title)
            #st<-strsplit(gene,split=";")
            #omimgen<-c(omimgen,as.character(toupper(st[[1]][2])))
            #}
            ##Extract Gene Symols throug Bionano 
<<<<<<< HEAD
            omimid <- datOmim[,1]
            geneOmim <- datOmim[,2]
            pas <- paste0("^",omimid,"$")
            pasq <- paste0("^",final_omim,"$")
            geneID <- c()
            for(ii in 1:length(pasq)){
                g1 <- grep(pasq[ii], pas, fixed = TRUE)
                geneID <- c(geneID, as.character(geneOmim[g1]))
            }
            geneID <- na.omit(as.character(geneID))
            gn1 <- tryCatch(
                mapIds(org.Hs.eg.db, geneID, "SYMBOL", "ENTREZID"),
                warning = function(w) {
                print("warning")
                return(NA)}, 
                error = function(e) {
                print("Error")
                return(NA)}
            )
            gn2 = tryCatch(    
            mapIds(org.Hs.eg.db, geneID, "ENSEMBL", "ENTREZID"),
            warning = function(w) {
                print("warning")
                return(NA)}, 
                error = function(e) {
                print("Error")
                return(NA)}
            )
            if (length(unique(as.character(na.omit(gn1)))) ==0) {gn1 <- c()} else {gn1 = gn1}
            if (length(unique(as.character(na.omit(gn2)))) ==0) {gn2 <- c()} else {gn2 = gn2}
            if(length(gn2) > 0){
			    rn <- row.names(data.frame(gn1))
                rn1 <- row.names(data.frame(gn2))
                gene1<-data.frame(
                    entrezID = as.numeric(rn),
                    ensemblID = as.character(data.frame(gn2)[,1]),
                    hgnc_symbol = as.character(data.frame(gn1)[,1])
                )
                gn<-as.character(unique(gene1$hgnc_symbol))
                gn_1_temp <- gn[gn != ""]
                gn_1 <- as.character(na.omit(gn_1_temp))
                omimGenes<-c(omimGenes,as.character(gn_1))
                    terms_list<-as.character(rep(paste(terms,"_OMIM",sep=""),
                    length(gn_1)))
                Final_terms_OMIM<-c(Final_terms_OMIM,terms_list)
			}else{print ("No genes for term !!")}
=======
            gn1 <- mapIds(org.Hs.eg.db, final_omim, "SYMBOL", "ENTREZID")
            gn2 <- mapIds(org.Hs.eg.db, final_omim, "ENSEMBL", "ENTREZID")
            rn <- row.names(data.frame(gn1))
            rn1 <- row.names(data.frame(gn2))
            gene1<-data.frame(
                entrezID = as.numeric(rn),
                ensemblID = as.character(data.frame(gn2)[,1]),
                hgnc_symbol = as.character(data.frame(gn1)[,1])
            )
            gn<-as.character(unique(gene1$hgnc_symbol))
            gn_1_temp <- gn[gn != ""]
            gn_1 <- as.character(na.omit(gn_1_temp))
            omimGenes<-c(omimGenes,as.character(gn_1))
                terms_list<-as.character(rep(paste(terms,"_Gene",sep=""),
                length(gn_1)))
            Final_terms_OMIM<-c(Final_terms_OMIM,terms_list)
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
        }
        else{
            print ("No genes for term !!")
        }
    }
    else{
        stop("Term length is zero")
    }
    ##Writing data into data frame and returning 
    dat2<-data.frame(omimGenes,Final_terms_OMIM)
    print(paste("OMIM GeneName length:",length(omimGenes),sep=""))
    return(dat2)
}
#' Extracting genes from gtr database NCBI.
#'
#' @param terms  Single or Multiple Terms.
<<<<<<< HEAD
#' @param gtr character gtr database location.
#' @param downloadGTR boolean If TRUE, download the gtr database.
#'  Default FALSE.
=======
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
#' @return Dataframe returned containing gene lists in entrezid and Gene 
#' Symbols, and terms associated with it
#' @examples
#' terms="Muscle Weakness"
<<<<<<< HEAD
#' ge <- gtr_gene(terms,gtr = "Y:/Suro/nanotatoRDatabases/gtr_07162019.txt")
#' @importFrom stats na.omit 
#' @export
gtr_gene<-function(terms, gtr, url_gtr, downloadGTR){
    ##Initialising variables
	print("###GTR data Analysis###")
    gtrGenes<-c()
    Final_terms_GTR<-c()
    print(downloadGTR)
    ##Initialising data to be extracted from ensembl
    httr::set_config(httr::config(http_version = 0))
    datGtr <- reading_GTR(gtr = gtr, 
	    url_gtr = url_gtr, 
		downloadGTR = downloadGTR)
=======
#' ge <- gtr_gene(terms)
#' @import rentrez utils
#' @import httr
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#' @importFrom stats na.omit 
#' @export
gtr_gene<-function(terms){
    ##Initialising variables
    gtrGenes<-c()
    Final_terms_GTR<-c()
    ##Initialising data to be extracted from ensembl
    httr::set_config(httr::config(http_version = 0))
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
    ####Checking for term size and extracting gene list accordingly
    if (length(terms)>1){
        for(ll in seq_len(length(terms))){##if term length greater than 1
            ##Extracting data from GTR database
<<<<<<< HEAD
            print(paste0("GTR:", ll))
            'e<-entrez_search(db="gtr",term=terms[ll],retmax=99999999)
            gtrgeneID<-e$ids'
            g1 <- grep(as.character(terms[ll]), 
                as.character(datGtr$DiseaseName), ignore.case = TRUE)
            geneID <- as.character(unique(na.omit(datGtr$GeneID[g1])))
            geneID <- gsub("N/A", NA, geneID)
            ##Checking if gene ID extracted is greater than 0
            if(length(geneID) > length(is.na(geneID))){
                gn<-as.character(unique(geneID))
                gn_1_temp <- gn[gn != ""]
                gn_1 <- as.character(na.omit(gn_1_temp))
                gtrGenes<-c(gtrGenes,as.character(gn_1))
                    terms_list<-as.character(
                        rep(paste(terms[ll],"_GTR",sep=""),
                        length(gn_1)))
                Final_terms_GTR<-c(Final_terms_GTR,terms_list)
                
            }
            else{
                print ("No genes for term !!")
=======
            e<-entrez_search(db="gtr",term=terms[ll],retmax=99999999)
            gtrgeneID<-e$ids
            ##Checking if gene ID extracted is greater than 0
            if(length(gtrgeneID)>0){
                #for (nn in 1:length(gtrgeneID)){
                #n<-entrez_summary(db="gtr", id=gtrgeneID[nn])
                #gtrGenes<-c(gtrGenes,n$testtargetlist)
                #}
                ##Extracting gene symbols through Biomart
                gn1 <- mapIds(org.Hs.eg.db, gtrgeneID, "SYMBOL", "ENTREZID")
                gn2 <- mapIds(org.Hs.eg.db, gtrgeneID, "ENSEMBL", "ENTREZID")
                rn <- row.names(data.frame(gn1))
                rn1 <- row.names(data.frame(gn2))
                gene1<-data.frame(
                    entrezID = as.numeric(rn),
                    ensemblID = as.character(data.frame(gn2)[,1]),
                    hgnc_symbol = as.character(data.frame(gn1)[,1])
                )
                gn<-as.character(unique(gene1$hgnc_symbol))
                gn_1_temp <- gn[gn != ""]
                gn_1 <- as.character(na.omit(gn_1_temp))
                gtrGenes<-c(gtrGenes,as.character(gn_1))
                    terms_list<-as.character(rep(paste(terms[ll],"_Gene",sep=""),
                    length(gn_1)))
                Final_terms_GTR <- c(Final_terms_GTR,terms_list)
                
            }
            else{
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
                next
            }
        }
    }
    else if (length(terms)==1){##For single terms
        ##Extracting data from GTR
<<<<<<< HEAD
        g1 <- grep(as.character(terms), as.character(datGtr$DiseaseName), 
                ignore.case = TRUE)
            geneID <- as.character(unique(na.omit(datGtr$GeneID[g1])))
            geneID <- gsub("N/A", NA, geneID)
            ##Checking if gene ID extracted is greater than 0
            if(length(geneID) > length(is.na(geneID))){
                gn<-as.character(unique(geneID))
                gn_1_temp <- gn[gn != ""]
                gn_1 <- as.character(na.omit(gn_1_temp))
                gtrGenes<-c(gtrGenes,as.character(gn_1))
                    terms_list<-as.character(rep(paste(terms,"_GTR",sep=""),
                    length(gn_1)))
                Final_terms_GTR<-c(Final_terms_GTR,terms_list)
                
            }
               else{
            print ("No genes for term !!")
        }        
        }
        
=======
        e<-entrez_search(db="gtr",term=terms,retmax=99999999)
        gtrgeneID<-e$ids
        ##If extracted gene length greater than 0
        if(length(gtrgeneID)>0){
            #for (nn in 1:length(gtrgeneID)){
            #n<-entrez_summary(db="gtr", id=gtrgeneID[nn])
            #gtrGenes<-c(gtrGenes,n$testtargetlist)
            #}
            ##Extracting gene symbols from Biomart
            gn1 <- mapIds(org.Hs.eg.db, gtrgeneID, "SYMBOL", "ENTREZID")
            gn2 <- mapIds(org.Hs.eg.db, gtrgeneID, "ENSEMBL", "ENTREZID")
            rn <- row.names(data.frame(gn1))
            rn1 <- row.names(data.frame(gn2))
            gene1<-data.frame(
                entrezID = as.numeric(rn),
                ensemblID = as.character(data.frame(gn2)[,1]),
                hgnc_symbol = as.character(data.frame(gn1)[,1])
            )
            gn<-as.character(unique(gene1$hgnc_symbol))
            gn_1_temp <- gn[gn != ""]
            gn_1 <- as.character(na.omit(gn_1_temp))
            gtrGenes <- c(gtrGenes,as.character(gn_1))
            terms_list<-as.character(rep(paste(terms , "_Gene",sep=""),
                length(gn_1)))
            Final_terms_GTR <- c(Final_terms_GTR,terms_list)
            
        }
        else{
            print ("No genes for term !!")
        }
    }
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
    else {##If term length is zero
        stop("Term length is zero")
    }
    ##Writing data into data frame and returning
    dat3<-data.frame(gtrGenes,Final_terms_GTR)
    print(paste("gtrGenes GeneName length:",length(gtrGenes),sep=""))
    return(dat3)
}
<<<<<<< HEAD
#' Extracting genes from gtr database NCBI.
#'
#' @param terms  Single or Multiple Terms.
#' @param gtr character gtr database location.
#' @param downloadGTR boolean If TRUE, download the gtr database.
#'  Default FALSE.
#' @param clinsig. character. clinical significance chosen. Choices are
#' "Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic",
#' "Benign", "Likely Benign", "Benign/Likely Benign","VUS", 
#' "Benign/Likely Benign/VUS", "All".
=======
#' Extracting genes from clinvar database NCBI.
#'
#' @param terms  Single or Multiple Terms.
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
#' @return Dataframe returned containing gene lists in entrezid and Gene 
#' Symbols, and terms associated with it
#' @examples
#' terms="Muscle Weakness"
<<<<<<< HEAD
#' ge <- gtr_gene(terms,gtr = "Y:/Suro/nanotatoRDatabases/gtr_07162019.txt")
#' @importFrom stats na.omit 
#' @import VarfromPDB
#' @export
#' @export
clinvar_gene<-function(terms, clinvar,downloadClinvar){
    print("###Clinvar data Analysis###")
	##Initializing variables
    clinvarGenes<-c()
    Final_terms_Clinvar<-c()
	clinicalSig <- c()
	#clinsig = clinsig
	#print(downloadClinvar)
    ##Initialising the database
    'httr::set_config(httr::config(http_version = 0))
    datclinvar <- reading_clinvar(clinvar = clinvar, 
        downloadClinvar = downloadClinvar, url_clinvar = url_clinvar)'
    if (length(terms)>1){
	    genes <- c()
		clinsigfig <- c()
        for(ll in seq_len(length(terms))){
            print(paste0("Clinvar:", ll))
            'e<-entrez_search(db="gtr",term=terms[ll],retmax=99999999)
            gtrgeneID<-e$ids'
			if(downloadClinvar == TRUE){
			    clinvar.phenotype = extract_clinvar(keyword= terms[ll])
			}else{
			    clinvar.phenotype = extract_clinvar(keyword= terms[ll], localPDB.path = clinvar)
			}
            'g1 <- grep(as.character(terms[ll]), as.character(datclinvar$DiseaseName), 
                ignore.case = TRUE)'
			dis2genes <- clinvar.phenotype$gene2dis
			disterms <- clinvar.phenotype$gene2dis$DiseaseName
			pterms <- terms[ll]
			'pterms1 <- paste0(";", terms[ll], ";")
			pterms2 <- paste0(";", terms[ll])
			pterms3 <- paste0(terms[ll], ";")'
		    g1 <- grep(pterms, as.character(disterms), ignore.case=TRUE)  
			'g2 <- grep(pterms1, as.character(disterms), ignore.case=TRUE)  
			g3 <- grep(pterms2, as.character(disterms), ignore.case=TRUE)  
			g4 <- grep(pterms3, as.character(disterms), ignore.case=TRUE)  '
			geneID <- unique(as.numeric(dis2genes$X.GeneID[g1]))
			omim <- unique(as.numeric(dis2genes$DiseaseMIM[g1]))
			varn <- clinvar.phenotype$variants
			st1 <- strsplit(as.character(varn$PhenotypeIDS), split = ",")
			varn$omimID <- sapply(st1, function(x) strsplit(x[2],split="OMIM:")[[1]][2])
			if (length(geneID) > 1){
			    
			    for(ii in 1:length(geneID)){
				    da1 <- varn[which(varn$GeneID == geneID[ii] & 
					    (varn$ClinicalSignificance == "Benign" 
						| varn$ClinicalSignificance == "Likely Benign"
						| varn$ClinicalSignificance == "Benign/Likely benign"
					    | varn$ClinicalSignificance == "Uncertain significance"
						| varn$ClinicalSignificance == "Pathogenic/Likely pathogenic"
						| varn$ClinicalSignificance == "Pathogenic"
						| varn$ClinicalSignificance == "Likely pathogenic")),]
					if(nrow(da1) > 0){
					    genes <- c(genes, as.character(unique(da1$GeneSymbol)))
					    g1 <- grep("Pathogenic/Likely pathogenic", da1$ClinicalSignificance, ignore.case = TRUE)
					    g2 <- grep("Pathogenic", da1$ClinicalSignificance, ignore.case = TRUE)
					    g3 <- grep("Likely pathogenic", da1$ClinicalSignificance, ignore.case = TRUE)
     				    if((length(g1) >= 1 & length(g2) >= 1 & length(g3) >= 1)
					       |(length(g1) >= 1 | length(g2) >= 1 | length(g3) >= 1)){
					        cs <-  "Pathogenic/Likely pathogenic"
					    }else{cs <-  as.character("-")}
				        'if(clinsig == "Pathogenic"){
			                da1 <- varn[which(varn$GeneID == geneID[ii] & 
					            (varn$ClinicalSignificance == "Pathogenic")
					        	& varn$omimID == omim[ii]),]
					    } else if(clinsig == "Likely Pathogenic"){
					         da1 <- varn[which(varn$GeneID == geneID[ii] & 
					            (varn$ClinicalSignificance == "Likely pathogenic") 
					        	& varn$omimID == omim[ii]),]
					     
					    }
					    else if(clinsig == "Pathogenic/Likely pathogenic"){
					         da1 <- varn[which(varn$GeneID == geneID[ii] & 
					            (varn$ClinicalSignificance == "Pathogenic/Likely pathogenic") 
					        	& varn$omimID == omim[ii]),]
					    }
					    else if(clinsig == "Benign"){
					         da1 <- varn[which(varn$GeneID == geneID[ii] & 
					            (varn$ClinicalSignificance == "Benign") 
					        	& varn$omimID == omim[ii]),]
					    }
					    else if(clinsig == "Benign/Likely benign"){
					         da1 <- varn[which(varn$GeneID == geneID[ii] & 
					            (varn$ClinicalSignificance == "Benign/Likely benign") 
					        	& varn$omimID == omim[ii]),]
					    	
					    }
					    else if(clinsig == "Likely benign"){
					         da1 <- varn[which(varn$GeneID == geneID[ii] & 
					            (varn$ClinicalSignificance == "Likely benign") 
					        	& varn$omimID == omim[ii]),]
					    	
					    	
					    }
					    else if(clinsig == "VUS"){
					         da1 <- varn[which(varn$GeneID == geneID[ii] & 
					            (varn$ClinicalSignificance == "Uncertain significance") 
					        	& varn$omimID == omim[ii]),]
					    	
					    }
					    else if(clinsig == "Benign/Likely Benign/VUS"){
					        da1 <- varn[which(varn$GeneID == geneID[ii] & 
					            (varn$ClinicalSignificance == "Benign" 
                                | varn$ClinicalSignificance == "Likely Benign"
					            | varn$ClinicalSignificance == "Uncertain significance")
					        	| varn$omimID == omim[ii]),]
					    	
					    	
					    }else if(clinsig == "All"){
					        
					    	
					    } else {
					        stop ("Clinicical Significance input parameter incorrect")
                        }'					
				        if(length(cs) > 1){
					        csp <- paste(cs, collapse = ",")
					    }else {csp <- cs}
				        clinsigfig <- c(clinsigfig,csp)
					
				    }else{ genes <- c(genes, NA)}
				}
			}
			else if (length(geneID) == 1){
			    da1 <- varn[which(varn$GeneID == geneID & 
					    (varn$ClinicalSignificance == "Benign" 
						| varn$ClinicalSignificance == "Likely Benign"
						| varn$ClinicalSignificance == "Benign/Likely benign"
					    | varn$ClinicalSignificance == "Uncertain significance"
						| varn$ClinicalSignificance == "Pathogenic/Likely pathogenic"
						| varn$ClinicalSignificance == "Pathogenic"
						| varn$ClinicalSignificance == "Likely pathogenic")),]
					if(nrow(da1) > 0){
				    #genes <- c(genes, as.character(unique(da1$GeneSymbol)))
					    g1 <- grep("Pathogenic/Likely pathogenic", da1$ClinicalSignificance, ignore.case = TRUE)
					    g2 <- grep("Pathogenic", da1$ClinicalSignificance, ignore.case = TRUE)
					    g3 <- grep("Likely pathogenic", da1$ClinicalSignificance, ignore.case = TRUE)
     				    if((length(g1) >= 1 & length(g2) >= 1 & length(g3) >= 1)
					       |(length(g1) >= 1 | length(g2) >= 1 | length(g3) >= 1)){
					        cs <-  "Pathogenic/Likely pathogenic"
					    }else{cs <-  as.character("-")}
			            'if(clinsig == "Pathogenic"){
			                da1 <- varn[which(varn$GeneID == geneID & 
					            (varn$ClinicalSignificance == "Pathogenic")
					        	& varn$omimID == omim),]
					    } else if(clinsig == "Likely Pathogenic"){
					         da1 <- varn[which(varn$GeneID == geneID & 
					            (varn$ClinicalSignificance == "Likely pathogenic") 
					        	& varn$omimID == omim),]
					     
					    }
					    else if(clinsig == "Pathogenic/Likely pathogenic"){
					         da1 <- varn[which(varn$GeneID == geneID & 
					            (varn$ClinicalSignificance == "Pathogenic/Likely pathogenic") 
					        	& varn$omimID == omim),]
					    }
					    else if(clinsig == "Benign"){
					         da1 <- varn[which(varn$GeneID == geneID & 
					            (varn$ClinicalSignificance == "Benign") 
					        	& varn$omimID == omim),]
					    }
					    else if(clinsig == "Benign/Likely benign"){
					         da1 <- varn[which(varn$GeneID == geneID & 
					            (varn$ClinicalSignificance == "Benign/Likely benign") 
					        	& varn$omimID == omim),]
					    }
					    else if(clinsig == "Likely benign"){
					         da1 <- varn[which(varn$GeneID == geneID & 
					            (varn$ClinicalSignificance == "Likely benign") 
					        	& varn$omimID == omim),]
					    }
					    else if(clinsig == "VUS"){
					         da1 <- varn[which(varn$GeneID == geneID & 
					            (varn$ClinicalSignificance == "Uncertain significance") 
					        	& varn$omimID == omim),]
					    }
					    else if(clinsig == "Benign/Likely Benign/VUS"){
					        da1 <- varn[which(varn$GeneID == geneID & 
					            (varn$ClinicalSignificance == "Benign" 
					    		| varn$ClinicalSignificance == "Likely Benign"
					            | varn$ClinicalSignificance == "Uncertain significance")
					        	& varn$omimID == omim),]
					    }else if(clinsig == "All"){
					        da1 <- varn[which(varn$GeneID == geneID & 
					            (varn$ClinicalSignificance == "Benign" 
					    		| varn$ClinicalSignificance == "Likely Benign"
					    		| varn$ClinicalSignificance == "Benign/Likely benign"
					            | varn$ClinicalSignificance == "Uncertain significance"
					    		| varn$ClinicalSignificance == "Pathogenic/Likely pathogenic"
					    		| varn$ClinicalSignificance == "Pathogenic"
					    		| varn$ClinicalSignificance == "Likely pathogenic")
					        	& varn$omimID == omim),]
					    } else {
					        stop ("Clinicical Significance input parameter incorrect")
                        }'
				        genes <- c(genes, as.character(unique(da1$GeneSymbol)))
				        #cs <-  as.character(unique(da1$ClinicalSignificance))
				        	if(length(cs) > 1){
				        	    csp <- paste(cs, collapse = ",")
				        	}
				        	else {csp <- cs}
				        clinsigfig <- c(clinsigfig,csp)
					}else{genes <- c(genes, NA)}
				
			}
			else{
			    genes <- c(genes, NA)		
			}
			
			
			
			
		}
		##Checking if gene ID extracted is greater than 0
        if(length(genes) > 0){
            gn<-as.character(unique(genes))
            gn_1_temp <- gn[gn != ""]
            gn_1 <- as.character(na.omit(gn_1_temp))
            clinvarGenes<-c(clinvarGenes,as.character(gn_1))
                terms_list<-as.character(rep(paste(terms[ll],"_ClinVar",sep=""),
                length(gn_1)))
				
            Final_terms_Clinvar<-c(Final_terms_Clinvar,terms_list)
            clinicalSig <- c(clinicalSig, as.character(clinsigfig))
        }else{
                print ("No genes for term !!")
                clinvarGenes <- "-"
			    Final_terms_Clinvar = "-"
			    clinicalSig = "-"
            }
    }
    else if (length(terms)==1){##Checking if term is single
            ##Extracting data from Clinvar
		    genes <- c()
			clinsigfig <- c()
			#print(clinsig)
            if(downloadClinvar == TRUE){
			    clinvar.phenotype = extract_clinvar(keyword= terms)
			}else{
			    clinvar.phenotype = extract_clinvar(keyword= terms, localPDB.path = clinvar)
			}
            'g1 <- grep(as.character(terms[ll]), as.character(datclinvar$DiseaseName), 
                ignore.case = TRUE)'
			dis2genes <- clinvar.phenotype$gene2dis
			disterms <- clinvar.phenotype$gene2dis$DiseaseName
			pterms <- terms
			'pterms1 <- paste0(";", terms[ll], ";")
			pterms2 <- paste0(";", terms[ll])
			pterms3 <- paste0(terms[ll], ";")'
		    g1 <- grep(pterms, as.character(disterms), ignore.case=TRUE)  
			'g2 <- grep(pterms1, as.character(disterms), ignore.case=TRUE)  
			g3 <- grep(pterms2, as.character(disterms), ignore.case=TRUE)  
			g4 <- grep(pterms3, as.character(disterms), ignore.case=TRUE)  '
			geneID <- unique(as.numeric(dis2genes$X.GeneID[g1]))
			omim <- unique(as.numeric(dis2genes$DiseaseMIM[g1]))
			varn <- clinvar.phenotype$variants
			st1 <- strsplit(as.character(varn$PhenotypeIDS), split = ",")
			varn$omimID <- sapply(st1, function(x) strsplit(x[2],split="OMIM:")[[1]][2])
			if (length(geneID) > 1){
			    
			    for(ii in 1:length(geneID)){
				    da1 <- varn[which(varn$GeneID == geneID[ii] & 
					    (varn$ClinicalSignificance == "Benign" 
						| varn$ClinicalSignificance == "Likely Benign"
						| varn$ClinicalSignificance == "Benign/Likely benign"
					    | varn$ClinicalSignificance == "Uncertain significance"
						| varn$ClinicalSignificance == "Pathogenic/Likely pathogenic"
						| varn$ClinicalSignificance == "Pathogenic"
						| varn$ClinicalSignificance == "Likely pathogenic")),]
					if(nrow(da1) > 0){
					    genes <- c(genes, as.character(unique(da1$GeneSymbol)))
					    g1 <- grep("Pathogenic/Likely pathogenic", da1$ClinicalSignificance, ignore.case = TRUE)
					    g2 <- grep("Pathogenic", da1$ClinicalSignificance, ignore.case = TRUE)
					    g3 <- grep("Likely pathogenic", da1$ClinicalSignificance, ignore.case = TRUE)
     				    if(((length(g1) >= 1 & length(g2) >= 1 & length(g3) >= 1)
					       |(length(g1) >= 1 | length(g2) >= 1 | length(g3) >= 1))){
					        cs <-  "Pathogenic/Likely pathogenic"
					    }else{cs <-  as.character("-")}
				        'if(clinsig == "Pathogenic"){
			                da1 <- varn[which(varn$GeneID == geneID[ii] & 
					            (varn$ClinicalSignificance == "Pathogenic")
					        	& varn$omimID == omim[ii]),]
					    } else if(clinsig == "Likely Pathogenic"){
					         da1 <- varn[which(varn$GeneID == geneID[ii] & 
					            (varn$ClinicalSignificance == "Likely pathogenic") 
					        	& varn$omimID == omim[ii]),]
					     
					    }
					    else if(clinsig == "Pathogenic/Likely pathogenic"){
					         da1 <- varn[which(varn$GeneID == geneID[ii] & 
					            (varn$ClinicalSignificance == "Pathogenic/Likely pathogenic") 
					        	& varn$omimID == omim[ii]),]
					    }
					    else if(clinsig == "Benign"){
					         da1 <- varn[which(varn$GeneID == geneID[ii] & 
					            (varn$ClinicalSignificance == "Benign") 
					        	& varn$omimID == omim[ii]),]
					    }
					    else if(clinsig == "Benign/Likely benign"){
					         da1 <- varn[which(varn$GeneID == geneID[ii] & 
					            (varn$ClinicalSignificance == "Benign/Likely benign") 
					        	& varn$omimID == omim[ii]),]
					    }
					    else if(clinsig == "Likely benign"){
					         da1 <- varn[which(varn$GeneID == geneID[ii] & 
					            (varn$ClinicalSignificance == "Likely benign") 
					        	& varn$omimID == omim[ii]),]
					    }
					    else if(clinsig == "VUS"){
					         da1 <- varn[which(varn$GeneID == geneID[ii] & 
					            (varn$ClinicalSignificance == "Uncertain significance") 
					        	& varn$omimID == omim[ii]),]
					    }
					    else if(clinsig == "Benign/Likely Benign/VUS"){
					        da1 <- varn[which(varn$GeneID == geneID[ii] & 
					            (varn$ClinicalSignificance == "Benign" 
                                | varn$ClinicalSignificance == "Likely Benign"
					            | varn$ClinicalSignificance == "Uncertain significance")
					        	& varn$omimID == omim[ii]),]
					    }else if(clinsig == "All"){
					        da1 <- varn[which(varn$GeneID == geneID[ii] & 
					            (varn$ClinicalSignificance == "Benign" 
					    		| varn$ClinicalSignificance == "Likely Benign"
					    		| varn$ClinicalSignificance == "Benign/Likely benign"
					            | varn$ClinicalSignificance == "Uncertain significance"
					    		| varn$ClinicalSignificance == "Pathogenic/Likely pathogenic"
					    		| varn$ClinicalSignificance == "Pathogenic"
					    		| varn$ClinicalSignificance == "Likely pathogenic")
					        	& varn$omimID == omim[ii]),]
					    } else {
					        stop ("Clinicical Significance input parameter incorrect")
                        }					
				        genes <- c(genes, as.character(unique(da1$GeneSymbol)))
					    cs <-  as.character(unique(da1$ClinicalSignificance))'
					    if(length(cs) > 1){
					        csp <- paste(cs, collapse = ",")
					    }
					    else {csp <- cs}
				        clinsigfig <- c(clinsigfig,csp)
					}else {genes <- c(genes, NA)}
				}
			}else if (length(geneID) == 1){
			    da1 <- varn[which(varn$GeneID == geneID & 
					    (varn$ClinicalSignificance == "Benign" 
						| varn$ClinicalSignificance == "Likely Benign"
						| varn$ClinicalSignificance == "Benign/Likely benign"
					    | varn$ClinicalSignificance == "Uncertain significance"
						| varn$ClinicalSignificance == "Pathogenic/Likely pathogenic"
						| varn$ClinicalSignificance == "Pathogenic"
						| varn$ClinicalSignificance == "Likely pathogenic")),]
					if(nrow(da1)>0){
					    genes <- c(genes, as.character(unique(da1$GeneSymbol)))
					    g1 <- grep("Pathogenic/Likely pathogenic", da1$ClinicalSignificance, ignore.case = TRUE)
					    g2 <- grep("Pathogenic", da1$ClinicalSignificance, ignore.case = TRUE)
					    g3 <- grep("Likely pathogenic", da1$ClinicalSignificance, ignore.case = TRUE)
     				    if(((length(g1) >= 1 & length(g2) >= 1 & length(g3) >= 1)
					       |(length(g1) >= 1 | length(g2) >= 1 | length(g3) >= 1))){
					        cs <-  "Pathogenic/Likely pathogenic"
					    }else{cs <-  as.character("-")}
			            'if(clinsig == "Pathogenic"){
			                da1 <- varn[which(varn$GeneID == geneID & 
					            (varn$ClinicalSignificance == "Pathogenic")
					        	& varn$omimID == omim),]
					    } else if(clinsig == "Likely Pathogenic"){
					         da1 <- varn[which(varn$GeneID == geneID & 
					            (varn$ClinicalSignificance == "Likely pathogenic") 
					        	& varn$omimID == omim),]
					     
					    }
					    else if(clinsig == "Pathogenic/Likely pathogenic"){
					         da1 <- varn[which(varn$GeneID == geneID & 
					            (varn$ClinicalSignificance == "Pathogenic/Likely pathogenic") 
					        	& varn$omimID == omim),]
					    }
					    else if(clinsig == "Benign"){
					         da1 <- varn[which(varn$GeneID == geneID & 
					            (varn$ClinicalSignificance == "Benign") 
					        	& varn$omimID == omim),]
					    }
					    else if(clinsig == "Benign/Likely benign"){
					         da1 <- varn[which(varn$GeneID == geneID & 
					            (varn$ClinicalSignificance == "Benign/Likely benign") 
					        	& varn$omimID == omim),]
					    }
					    else if(clinsig == "Likely benign"){
					         da1 <- varn[which(varn$GeneID == geneID & 
					            (varn$ClinicalSignificance == "Likely benign") 
					        	& varn$omimID == omim),]
					    }
					    else if(clinsig == "VUS"){
					         da1 <- varn[which(varn$GeneID == geneID & 
					            (varn$ClinicalSignificance == "Uncertain significance") 
					        	& varn$omimID == omim),]
					    }
					    else if(clinsig == "Benign/Likely Benign/VUS"){
					        da1 <- varn[which(varn$GeneID == geneID & 
					            (varn$ClinicalSignificance == "Benign" 
					    		| varn$ClinicalSignificance == "Likely Benign"
					            | varn$ClinicalSignificance == "Uncertain significance")
					        	& varn$omimID == omim),]
					    }else if(clinsig == "All"){
					        da1 <- varn[which(varn$GeneID == geneID & 
					            (varn$ClinicalSignificance == "Benign" 
					    		| varn$ClinicalSignificance == "Likely Benign"
					    		| varn$ClinicalSignificance == "Benign/Likely benign"
					            | varn$ClinicalSignificance == "Uncertain significance"
					    		| varn$ClinicalSignificance == "Pathogenic/Likely pathogenic"
					    		| varn$ClinicalSignificance == "Pathogenic"
					    		| varn$ClinicalSignificance == "Likely pathogenic")
					        	& varn$omimID == omim),]
					    } else {
					        stop ("Clinicical Significance input parameter incorrect")
                        }
				gene    s <- c(genes, as.character(da1$GeneSymbol))
				cs <    -  as.character(unique(da1$ClinicalSignificance))'
					    if(length(cs) > 1){
					        csp <- paste(cs, collapse = ",")
					    }
					    else {csp <- cs}
				        clinsigfig <- c(clinsigfig,csp)
					}else{genes <-c(genes, NA)}
			}else{
			    genes <- c(genes, NA)
			}
			
			
			
			
		
		##Checking if gene ID extracted is greater than 0
        if(length(genes) > 0){
            gn<-as.character(unique(genes))
            gn_1_temp <- gn[gn != ""]
            gn_1 <- as.character(na.omit(gn_1_temp))
            clinvarGenes<-c(clinvarGenes,as.character(gn_1))
                terms_list<-as.character(rep(paste(terms,"_ClinVar",sep=""),
                length(gn_1)))
            Final_terms_Clinvar<-c(Final_terms_Clinvar,terms_list)
            clinicalSig <- c(clinicalSig, as.character(clinsigfig))
        } else{
            print ("No genes for term !!")
            clinvarGenes <- "-"
			Final_terms_Clinvar = "-"
			clinicalSig = "-"
        }
    }else{
=======
#' ge <- clinvar_gene(terms)
#' @import rentrez utils
#' @import httr
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#' @importFrom stats na.omit 
#' @export
clinvar_gene<-function(terms){
    ##Initializing variables
    clinvarGenes<-c()
    Final_terms_Clinvar<-c()
    ##Initialising the database
    httr::set_config(httr::config(http_version = 0))
    if (length(terms)>1){
        for(ll in seq_len(length(terms))){
            ##Extracting data from Clinvar
            f<-entrez_search(db="clinvar",term=paste(terms[ll],"[WORD]",sep=""),
                retmax=99999999)
            g<-entrez_search(db="clinvar",term=paste(terms[ll],"[DIS]",sep=""),
                retmax=99999999)
            h<-entrez_search(db="clinvar",term=terms[ll],retmax=99999999)
            clinvargeneID_word<-f$ids
            clinvargeneID_dis<-g$ids
            clinvargeneID_all<-h$ids
            final_clinvar<-unique(c(as.character(clinvargeneID_word),
                as.character(clinvargeneID_dis),
                        as.character(clinvargeneID_all)))
            ##Checking whether extracted gene lengths greater than 0
            if(length(final_clinvar)>0){
                #clinvargen<-c()
                #for (tt in 1:length(final_clinvar)){
                #t<-entrez_summary(db="clinvar", id=final_clinvar[tt])
                #clinvargen<-c(clinvargen,as.character(t$genes$symbol))
                #}
                ##Extracting gene symbols through biomart
                gn1 <- mapIds(org.Hs.eg.db, final_clinvar, "SYMBOL", "ENTREZID")
                gn2 <- mapIds(org.Hs.eg.db, final_clinvar, "ENSEMBL", "ENTREZID")
                rn <- row.names(data.frame(gn1))
                rn1 <- row.names(data.frame(gn2))
                gene1<-data.frame(
                    entrezID = as.numeric(rn),
                    ensemblID = as.character(data.frame(gn2)[,1]),
                    hgnc_symbol = as.character(data.frame(gn1)[,1])
                )
                gn<-as.character(unique(gene1$hgnc_symbol))
                gn_1_temp <- gn[gn != ""]
                gn_1 <- as.character(na.omit(gn_1_temp))
                clinvarGenes <- c(clinvarGenes,as.character(gn_1))
                terms_list<-as.character(rep(paste(terms [ll] , 
                    "_Gene",sep=""), length(gn_1)))
                Final_terms_Clinvar <- c(Final_terms_Clinvar,terms_list)
            }
            else{
                next
            }
        }
    }
    else if (length(terms)==1){##Checking if term is single
        ##Extracting data from Clinvar
        f<-entrez_search(db="clinvar",term=paste(terms,"[WORD]",sep=""),
                retmax=99999999)
        g<-entrez_search(db="clinvar",term=paste(terms,"[DIS]",sep="")
                ,retmax=99999999)
        h<-entrez_search(db="clinvar",term=terms,retmax=99999999)
        clinvargeneID_word<-f$ids
        clinvargeneID_dis<-g$ids
        clinvargeneID_all<-h$ids
        final_clinvar<-unique(c(as.character(clinvargeneID_word),
                        as.character(clinvargeneID_dis),
                        as.character(clinvargeneID_all)))
        ##Checking whether extracted gene lengths greater than 0
        if(length(final_clinvar)>0){
            #clinvargen<-c()
            #for (tt in 1:length(final_clinvar)){
            #t<-entrez_summary(db="clinvar", id=final_clinvar[tt])
            #clinvargen<-c(clinvargen,as.character(t$genes$symbol))
            #}
            ##Extracting gene symbol through biomart
            gn1 <- mapIds(org.Hs.eg.db, final_clinvar, "SYMBOL", "ENTREZID")
            gn2 <- mapIds(org.Hs.eg.db, final_clinvar, "ENSEMBL", "ENTREZID")
            rn <- row.names(data.frame(gn1))
            rn1 <- row.names(data.frame(gn2))
            gene1<-data.frame(
                entrezID = as.numeric(rn),
                ensemblID = as.character(data.frame(gn2)[,1]),
                hgnc_symbol = as.character(data.frame(gn1)[,1])
            )
            gn<-as.character(unique(gene1$hgnc_symbol))
            gn_1_temp <- gn[gn != ""]
            gn_1 <- as.character(na.omit(gn_1_temp))
            clinvarGenes<-c(clinvarGenes,as.character(gn_1))
            terms_list<-as.character(rep(paste(terms , "_Gene",sep=""),
                length(gn_1)))
            Final_terms_Clinvar<-c(Final_terms_Clinvar,terms_list)
        }
        else{
            print ("No genes for term !!")
        }
    }
    else{
>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06
        ##If term length is zero
        stop("Term length is zero")
    }
    ##Writing gene and terms into dataframe and returning
<<<<<<< HEAD
    dat4 <- data.frame(clinvarGenes,Final_terms_Clinvar, clinicalSig)
    #print(paste("clinvarGenes GeneName length:",length(clinvarGenes),sep=""))
    return(dat4)
}
#' Reading and parsing OMIM database.
#'
#' @param omim character omim database location.
#' @return Dataframe returned containing Omim ID and gene IDs. 
#' @examples
#' reading_mim2gene(omim  = "Y:/Suro/nanotatoRDatabases/mim2gene.txt")
#' @importFrom stats na.omit 
#' @import utils
#' @export
reading_mim2gene <- function(omim){
            omim = omim
            print(omim)
            print(paste0("omim_gene",omim))
            con <- file(omim, "r")
            r10 <- readLines(con, n = -1)
            close(con)
            # datfinal<-data.frame()
            # datfinal<-data.frame()
            g1 <- grep("MIM Number", r10)
            g2 <- grep("Entrez Gene ID \\(NCBI\\)", r10)
            print(g1)
            print(g2)
            'gg1 <- grep("# BSPQI Sample", r10)
            stt <- strsplit(r10[gg1], split = ":")
            fname_temp <- stt[[1]][2]'
            #print (paste0("SampleName:", fname))
            if (g1 == g2) {
                dat <- gsub("#", "", as.character(r10))
                # dat<-gsub('\t',' ',r10)
                dat4 <- textConnection(dat[g1:length(dat)])
                r1 <- read.table(dat4, sep = "\t", header = TRUE)
                close(dat4)
            } else {
                stop("column names doesnot Match")
            }
    dat <- data.frame(OMIMID = as.character(r1[,1]), GeneID = as.character(r1[,3]))
    return(dat)

}

#' Reading and parsing gtr database.
#'
#' @param gtr character gtr database location.
#' @param downloadGTR logical if true, downloads gtr database, .
#' and store data in the gtr location, else reads dataset from
#' gtr location.
#' @return Dataframe representation of gtr database.
#' @examples
#' a <- reading_GTR(gtr  = "Y:/Suro/nanotatoRDatabases/gtr_07162019.txt",
#'        downloadGTR = FALSE)
#' @importFrom stats na.omit 
#' @import utils
#' @export
reading_GTR <- function(gtr, downloadGTR, 
        url_gtr = "ftp://ftp.ncbi.nlm.nih.gov/pub/GTR/data/test_condition_gene.txt"){
        downloadGTR = downloadGTR
		if(downloadGTR == TRUE){
        download.file(url =  url_gtr, 
        destfile = gtr)
        tab <- read.table(gtr, sep="\t", header=TRUE, comment.char="",
            na.strings=".", stringsAsFactors=FALSE,
            quote="", fill=FALSE)
        dat <- data.frame(GeneID = as.character(tab$gene_symbol), 
            DiseaseName = as.character(tab$object_name))
        }
        else{
        tab <- read.table(gtr, sep="\t", header=TRUE, comment.char="",
            na.strings=".", stringsAsFactors=FALSE,
            quote="", fill=FALSE)
        dat <- data.frame(GeneID = as.character(tab$gene_symbol), 
            DiseaseName = as.character(tab$object_name))
        }
return(dat)
}
=======
    dat4<-data.frame(clinvarGenes,Final_terms_Clinvar)
    #print(paste("clinvarGenes GeneName length:",length(clinvarGenes),sep=""))
    return(dat4)
}

>>>>>>> 04c1b43bc3aded9fbd47aace30093a5b16389a06


