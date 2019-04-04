#' Extracting genes for phenotype/diseases from NCBI.
#'
#' @param method_entrez  character. Input Method for terms. Choices are 
#'                "Single","Multiple" and "Text".
#' @param termPath  character. Path and file name for textfile.
#' @param term  character. Single or Multiple Terms.
#' @param outpath character. Path where gene lists are saved.
#' @param thresh integer. Threshold for the number of terms sent to entrez.
#'                Note if large lists are sent to ncbi, it might fail to get
#'                processed. Default is 5.
#' @param returnMethod Method of returning output. Options, Text or data.frame.
#' @return Text files containg gene list and terms associated with them
#'         are stored as text files.
#' @examples
#' terms="Muscle Weakness"
#' gene_list_generation(method="Single", term=terms, returnMethod="dataFrame")
#' @import stats
#' @import rentrez 
#' @import biomaRt
#' @import utils
#' @export
gene_list_generation<-function(method_entrez = c("Single","Multiple","Text"), 
        termPath, 
        term, outpath, thresh=5, 
        returnMethod=c("Text","dataFrame")){
        
    #setwd(path)
    thresh=5
    ##Checking for whether the input method is commandline (Single or Multiple)
    ## or from Text
    if (method_entrez == "Single"){
        terms<-term
    }
    else if (method_entrez == "Multiple"){
        terms<- as.character(term)
        terms<- as.character(term) 
    }
    else if (method_entrez == "Text"){
        r<-read.csv(termPath)
        terms<-as.character(r$Term)
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
        terms_1<-as.character(terms[1:thresh])
        terms_2<-as.character(terms[thresh:length(terms)])
        ##Extracting genes from the databases
        g<-gene_extraction(terms_1)
        g_1<-gene_extraction(terms_2)
        g_o<-omim_gene(terms_1)
        g_o_1<-omim_gene(terms_2)
        g_g<-gtr_gene(terms_1)
        g_g_1<-gtr_gene(terms_2)
        g_c<-clinvar_gene(terms_1)
        g_c_1<-clinvar_gene(terms_2)
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
        ##Collating output from all the datasets
        Genes<-c(as.character(g$geneName),as.character(g_g$gtrGenes),
            as.character(g_o$omimGenes),as.character(g_c$clinvarGenes))
        Terms<-c(as.character(g$Final_terms),
            as.character(g_g$Final_terms_GTR),
            as.character(g_o$Final_terms_OMIM),
            as.character(g_c$Final_terms_Clinvar))
        ##Writing genes and terms in a dataset
        dat<-data.frame(Genes,Terms)
    }
    ##We write out all the genes and terms from the 
    ##different databases into a$description
    ##single dataframe
    gene1<-c()
    term1<-c()
    uqGenes<-as.character(unique(na.exclude(Genes)))
    for(ii in 1:length(uqGenes)){
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
#' @import rentrez biomaRt utils
#' @importFrom stats na.omit 
#' @export
gene_extraction<-function(terms){
##Initializing values
geneName<-c()
Final_terms <-c()
##Initialising data to be extracted from ensembl
ensembl = useMart("ensembl",host = "www.ensembl.org", 
ensemblRedirect = FALSE,dataset="hsapiens_gene_ensembl")
##Checking for term size and extracting gene list accordingly
if (length(terms)>1){ #checking for terms greater than 1
    for(ll in seq_len(length(terms))){
        ##Searching Entrez with the terms
        b<-entrez_search(db="gene",term=terms[ll],retmax=99999999)
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
            gene1 = getBM(attributes = c("entrezgene", "ensembl_gene_id", 
            "hgnc_symbol"), filters = "entrezgene", 
            values = geneID, mart = ensembl)
            gn<-as.character(unique(gene1$hgnc_symbol))
            gn_1<-gn[gn != ""]
            geneName<-c(geneName,as.character(gn_1))
            terms_list<-as.character(rep(paste(terms,"_Gene",sep=""),
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
        gene1 = getBM(attributes = c("entrezgene", "ensembl_gene_id", 
        "hgnc_symbol"), filters = "entrezgene", 
        values = geneID, mart = ensembl)
        gn<-as.character(unique(gene1$hgnc_symbol))
        gn_1<-gn[gn != ""]
        geneName<-c(geneName,as.character(gn_1))
        terms_list<-as.character(rep(paste(terms,"_Gene",sep=""),length(gn_1)))
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
#' @return Dataframe returned containing gene lists in entrezid and Gene 
#' Symbols, and terms associated with it
#' @examples
#' terms="Muscle Weakness"
#' omim_gene(terms)
#' @import rentrez biomaRt utils
#' @importFrom stats na.omit 
#' @export
omim_gene<-function(terms){
    ##Initialising variables
    omimGenes<-c()
    Final_terms_OMIM<-c()
    ##Initialising data to be extracted from ensembl
    ensembl = useMart("ensembl",host = "www.ensembl.org", 
            ensemblRedirect = FALSE,dataset="hsapiens_gene_ensembl")
    ####Checking for term size and extracting gene list accordingly
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
            ##Checking if geneids greater than 0, if not moved to the next term
            if(length(final_omim)>0){
                #omimgen<-c()
                #for (kk in 1:length(final_omim)){
                #k<-entrez_summary(db="omim", id=final_omim[kk])
                #gene<-as.character(k$title)
                #st<-strsplit(gene,split=";")
                #omimgen<-c(omimgen,as.character(toupper(st[[1]][2])))
                #}
                ##Get the Human gene symbols for the extracted Gene IDs
                gene2 = getBM(attributes = c("entrezgene", "ensembl_gene_id", 
                                "hgnc_symbol"), filters = "entrezgene", 
                        values = final_omim, mart = ensembl)
                gn2<-as.character(unique(gene2$hgnc_symbol))
                gn_2<-gn2[gn2 != ""]
                omimGenes<-c(omimGenes,as.character(gn_2))
                terms_list<-as.character(rep(paste(terms[ll],
                "_OMIMGene",sep=""),length(gn_2)))
                ##Writing it into a data frame
                Final_terms_OMIM<-c(Final_terms_OMIM,terms_list)
            }else{next}        
        }
    }
    else if (length(terms)==1){
        ##Extracting data from OMIM
        c<-entrez_search(db="omim",term=paste(terms,"[WORD]",sep=""),
            retmax=99999999)
        d<-entrez_search(db="omim",term=paste(terms,"[DIS]",sep=""),
                retmax=99999999)
        omimgeneID_word<-c$ids
        omimgeneID_dis<-d$ids
        final_omim<-unique(c(as.character(omimgeneID_word),as.character(
            omimgeneID_dis)))
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
            gene2 = getBM(attributes = c("entrezgene", "ensembl_gene_id", 
                "hgnc_symbol"), 
                filters = "entrezgene", values = final_omim, mart = ensembl)
            gn2<-as.character(unique(gene2$hgnc_symbol))
            gn_2<-gn2[gn2 != ""]
            omimGenes<-c(omimGenes,as.character(gn_2))
            terms_list<-as.character(rep(paste(terms,"_OMIMGene",sep=""),
                length(gn_2)))
            Final_terms_OMIM<-c(Final_terms_OMIM,terms_list)
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
#' @return Dataframe returned containing gene lists in entrezid and Gene 
#' Symbols, and terms associated with it
#' @examples
#' terms="Muscle Weakness"
#' gtr_gene(terms)
#' @import rentrez biomaRt utils
#' @importFrom stats na.omit 
#' @export
gtr_gene<-function(terms){
    ##Initialising variables
    gtrGenes<-c()
    Final_terms_GTR<-c()
    ##Initialising data to be extracted from ensembl
    ensembl = useMart("ensembl",host = "www.ensembl.org", 
            ensemblRedirect = FALSE,dataset="hsapiens_gene_ensembl")
    ####Checking for term size and extracting gene list accordingly
    if (length(terms)>1){
        for(ll in seq_len(length(terms))){##if term length greater than 1
            ##Extracting data from GTR database
            e<-entrez_search(db="gtr",term=terms[ll],retmax=99999999)
            gtrgeneID<-e$ids
            ##Checking if gene ID extracted is greater than 0
            if(length(gtrgeneID)>0){
                #for (nn in 1:length(gtrgeneID)){
                #n<-entrez_summary(db="gtr", id=gtrgeneID[nn])
                #gtrGenes<-c(gtrGenes,n$testtargetlist)
                #}
                ##Extracting gene symbols through Biomart
                gene3 = getBM(attributes = c("entrezgene", "ensembl_gene_id", 
                    "hgnc_symbol"), 
                    filters = "entrezgene", values = gtrgeneID, mart = ensembl)
                gn3<-as.character(unique(gene3$hgnc_symbol))
                gn_3<-gn3[gn3 != ""]
                gtrGenes<-c(gtrGenes,as.character(gn_3))
                terms_list<-as.character(rep(paste(terms[ll],
                    "_GTRGene",sep=""),length(gn_3)))
                Final_terms_GTR<-c(Final_terms_GTR,terms_list)
                
            }else{next}
        }
    }
    else if (length(terms)==1){##For single terms
        ##Extracting data from GTR
        e<-entrez_search(db="gtr",term=terms,retmax=99999999)
        gtrgeneID<-e$ids
        ##If extracted gene length greater than 0
        if(length(gtrgeneID)>0){
            #for (nn in 1:length(gtrgeneID)){
            #n<-entrez_summary(db="gtr", id=gtrgeneID[nn])
            #gtrGenes<-c(gtrGenes,n$testtargetlist)
            #}
            ##Extracting gene symbols from Biomart
            gene3 = getBM(attributes = c("entrezgene", "ensembl_gene_id", 
                    "hgnc_symbol"), 
                    filters = "entrezgene", values = gtrgeneID, mart = ensembl)
            gn3<-as.character(unique(gene3$hgnc_symbol))
            gn_3<-gn3[gn3 != ""]
            gtrGenes<-c(gtrGenes,as.character(gn_3))
            terms_list<-as.character(rep(paste(terms,"_GTRGene",sep=""),
                length(gn_3)))
            Final_terms_GTR<-c(Final_terms_GTR,terms_list)
            
        }
        else{
            print ("No genes for term !!")
        }
    }
    else {##If term length is zero
        stop("Term length is zero")
    }
    ##Writing data into data frame and returning
    dat3<-data.frame(gtrGenes,Final_terms_GTR)
    print(paste("gtrGenes GeneName length:",length(gtrGenes),sep=""))
    return(dat3)
}
#' Extracting genes from clinvar database NCBI.
#'
#' @param terms  Single or Multiple Terms.
#' @return Dataframe returned containing gene lists in entrezid and Gene 
#' Symbols, and terms associated with it
#' @examples
#' terms="Muscle Weakness"
#' clinvar_gene(terms)
#' @import rentrez biomaRt utils
#' @importFrom stats na.omit 
#' @export
clinvar_gene<-function(terms){
    ##Initializing variables
    clinvarGenes<-c()
    Final_terms_Clinvar<-c()
    ##Initialising the database
    ensembl = useMart("ensembl",host = "www.ensembl.org", 
            ensemblRedirect = FALSE,dataset="hsapiens_gene_ensembl")
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
                gene4 = getBM(attributes = c("entrezgene", "ensembl_gene_id",
                                "hgnc_symbol"), filters = "entrezgene", 
                        values = final_clinvar, mart = ensembl)
                gn4<-as.character(unique(gene4$hgnc_symbol))
                gn_4<-gn4[gn4 != ""]
                clinvarGenes<-c(clinvarGenes,as.character(gn_4))
                terms_list<-as.character(rep(paste(terms[ll],"_ClinVarGene",
                    sep=""),
                                length(gn_4)))
                Final_terms_Clinvar<-c(Final_terms_Clinvar,terms_list)
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
            gene4 = getBM(attributes = c("entrezgene", "ensembl_gene_id", 
                "hgnc_symbol"), 
                filters = "entrezgene", values = final_clinvar, mart = ensembl)
            gn4<-as.character(unique(gene4$hgnc_symbol))
            gn_4<-gn4[gn4 != ""]
            clinvarGenes<-c(clinvarGenes,as.character(gn_4))
            terms_list<-as.character(rep(paste(terms,"_ClinVarGene",sep=""),
                length(gn_4)))
            Final_terms_Clinvar<-c(Final_terms_Clinvar,terms_list)
        }
        else{
            print ("No genes for term !!")
        }
    }
    else{
        ##If term length is zero
        stop("Term length is zero")
    }
    ##Writing gene and terms into dataframe and returning
    dat4<-data.frame(clinvarGenes,Final_terms_Clinvar)
    #print(paste("clinvarGenes GeneName length:",length(clinvarGenes),sep=""))
    return(dat4)
}





