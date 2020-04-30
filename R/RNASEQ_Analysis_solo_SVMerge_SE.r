#' Combining the RNAseq reads of family members in a 
#' single file.
#'
#' @param RNASeqDir  character. Directory containing RNAseq reads.
#' @param returnMethod  character. Method of returning Data.
#' @param outpath  character. Contains file path if Method of return is chosen as 
#' Text.
#' @param outFileName  character. Output file name. 
#' @return Text or Dataframe containing TPM read counts of genes in the family.
#' @examples
#' RNASeqDir = system.file("extdata", package="nanotatoR")
#' returnMethod="dataFrame"
#' datRNASeq <- RNAseqcombine_solo(RNASeqDir = RNASeqDir,
#' returnMethod = returnMethod)
#' @importFrom AnnotationDbi mapIds
#' @importFrom stats na.omit 
#' @import org.Hs.eg.db
#' @export
RNAseqcombine_solo<-function(RNASeqDir,returnMethod=c("Text","dataFrame"),
                        outpath="",outFileName=""){ 
    #library()
    #setwd(RNASeqDir)
    l <- list.files(path = RNASeqDir, 
       pattern="*.genes.results", full.names = TRUE)
    len<-length(l)
    #dat<-listDatasets(ensembl)
    #g1<-grep("sscrofa",listDatasets(ensembl)$dataset)
    'grch37 = useMart(="ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", 
                    dataset="hsapiens_gene_ensembl")'
    gen<-c(); 
    ##Need to make this function dynamic
    dat<-data.frame(matrix(ncol=len,nrow=nrow(r<-read.table(l[1],sep="\t",header=TRUE))))
    cnam<-c()
    for (ii in 1:length(l)){
        r<-read.table(l[ii],sep="\t",header=TRUE)
        gen<-c(gen,as.character(r$gene_id))
        dat[,ii]<-as.numeric(r$TPM)
        str<-strsplit(l[ii],split=".genes.results")
        #print(str[[1]][1])
        cnam<-c(cnam,str[[1]][1])
    }
    #datf<-data.frame(dat)
    gen<-unique(as.character(gen))
    st<-strsplit(gen,split="[.]")
    genes<-c()
    for(k in 1:length(gen)){
        genes<-c(genes,as.character(st[[k]][1]))
    }
    data1 <- data.frame(dat[,1])
    names(data1) <- cnam
    genesym <- c()
    ensemblid <- c()
    'gene1 = getBM(attributes = c("external_gene_name", "ensembl_gene_id"), 
                filters = "ensembl_gene_id", values = genes, mart = grch37)'
    gn1 <- mapIds(org.Hs.eg.db, genes, "SYMBOL", "ENSEMBL")
    'rn <- row.names(data.frame(gn1))
    rn1 <- row.names(data.frame(gn2))'
    gene1<-data.frame(
        ensembl_gene_id = as.character(names((gn1))),
        external_gene_name = as.character(data.frame(gn1)[,1])
        )
    genesym<-as.character(gene1$external_gene_name)
    ensemblid<-as.character(gene1$ensembl_gene_id)
    gene3<-c()
    ens<-c()
    ensemblid<-paste("^",ensemblid,"$",sep="")
    for(kk in 1:length(genes)){
        pag<-paste("^",genes[kk],"$",sep="")
        val<-grep(pag,ensemblid,fixed=TRUE)
        if(length(val)>0){
            gene3<-c(gene3,as.character(unique(genesym[val])))
            ens<-c(ens,as.character(genes[kk]))
        }
        else{
            gene3<-c(gene3,as.character("-"))
            ens<-c(ens,as.character(genes[kk]))
        }
    }
    RNASeqDat<-data.frame(GeneName = as.character(gene3),
            GeneID = as.character(ens),data1)
    if(returnMethod=="Text"){
        fname = file.path(outFileName,".csv",sep = "")
        write.csv(RNASeqDat,file.path(outpath,fname), row.names = FALSE )
    }
    else if (returnMethod=="dataFrame"){
        return (RNASeqDat)
    }
    else{
        stop("Invalid ReturnMethod")
    }
}
#' Annotating the Overlapping genes with RNAseq expression
#'
#' @param gnsOverlap  character. Vector containing overlapping genes.
#' @param SVID  character. SV Index ID.
#' @param RNASeqData  dataFrame. RNAseq data with gene names. 
#' @param pattern_Proband  character. Pattern for proband. 
#' @return Dataframe containing TPM read counts of overlapping genes.
#' @examples
#' RNASeqDir = system.file("extdata", package="nanotatoR")
#' returnMethod="dataFrame"
#' datRNASeq <- RNAseqcombine_solo(RNASeqDir = RNASeqDir,
#' returnMethod = returnMethod)
#' gnsOverlap <- c("AGL")
#' SVID = 397
#' datgnovrlap <- OverlapRNAseq_solo(gnsOverlap = gnsOverlap, 
#' SVID = SVID, RNASeqData = datRNASeq,
#' pattern_Proband = "*_P_*")
#' @importFrom stats na.omit  
#' @export
OverlapRNAseq_solo<-function(gnsOverlap, SVID, RNASeqData,
                pattern_Proband = NA){
  
  ###Finding the column names
   #print(gnsOverlap);print(length(gnsOverlap))
  'if(is.na(pattern_Father)==FALSE){
    fatherInd<-grep(pattern_Father,names(RNASeqData))
    } else{fatherInd<- NA}
    if(is.na(pattern_Mother)==FALSE){
    motherInd<-grep(pattern_Mother,names(RNASeqData))
    } else{motherInd<- NA}'
    if(is.na(pattern_Proband)==FALSE){
        probandInd <- grep(pattern_Proband,names(RNASeqData))
    } else{probandInd<- NA}
    'if(is.na(pattern_Sibling)==FALSE){
      siblingInd<-grep(pattern_Sibling,names(RNASeqData))
    } 
    else{
      siblingInd<-NA
    }'
    sv<-c()
    gene<-c()
    gnsname<-as.character(RNASeqData$GeneName)
    pasgnsname<-pasgnovlap<-paste("^",gnsname,"$",sep="")
    'overlap_ensemblgenes = select(EnsDb.Hsapiens.v79, gnsOverlap, 
                c("GENEID","GENENAME"), "SYMBOL")
    gnsOverlapID<-as.character(overlap_ensemblgenes$GENEID)'
    #print(paste("gnsOverlap:",gnsOverlap))
    #print(paste("overlap_ensemblgenes:",overlap_ensemblgenes))
    #print(paste("gnsOverlapID:",gnsOverlapID))
    #genes<-as.character(overlap_ensemblgenes$SYMBOL)
    ###Extracting Reads
    ###
    ###Genes Names Extraction
    #print(paste("fatherInd :",fatherInd))
    #print(paste("motherInd :",motherInd))
    #print(paste("probandInd :",probandInd))
    #print(paste("siblingInd :",siblingInd))
    gnsOverlapID <- as.character(gnsOverlap)
    #print(gnsOverlapID)
    if(length(gnsOverlapID)>1){
    
        datGeneInfoTemp<-data.frame()
        fatherReads<-c()
        motherReads<-c()
        probandReads<-c()
        siblingReads<-c()
        for (ki in 1:length(gnsOverlapID)){
            pasgnovlap <- paste("^", gnsOverlapID[ki],"$", sep = "")
            #print(ki)
            gg<-grep(pasgnovlap, pasgnsname, fixed = TRUE)
            dat_temp<-RNASeqData[gg, ]
        
            if(nrow(dat_temp)>1){
                #print("FALSE")
                #print(ki)
                dat_temp_1 <- mean(dat_temp[,probandInd])
                #dat_temp_1<-data.frame(dat_temp_1)
                dat_temp1 <- dat_temp[1,]
                dat_temp1[,probandInd] <- dat_temp_1
                #print(dim(dat_temp1))
            
                if(is.na(probandInd[1])==FALSE){
                    if(length(probandInd)>1){
                        for(j in probandInd){
                            probandcount<-c(probandcount,dat_temp1[,j])
                        }
                        probandReads<-c(
                            probandReads,paste(probandcount,collapse = ":"))
                    }
                    else if(length(probandInd)==1){
                        probandReads<-c(probandReads,dat_temp1[,probandInd])
                    }
                    else{
                        probandReads<-c(probandReads,0)
                    }
                }
                else{
                    probandReads<-c(probandReads,"-")
                }
            }
            else if (nrow(dat_temp)==1) {
                #print("TRUE")
                #print(dim(dat_temp1))
                probandcount<-c()
                if(is.na(probandInd [1])==FALSE){ 
                    if(length(probandInd)>1){
                        for(j in probandInd){
                            probandcount<-c(probandcount,dat_temp[,j])
                        }
                        probandReads<-c(
                        probandReads, paste(probandcount,collapse = ":"))
                    }
                    else if(length(probandInd)==1){
                        probandReads<-c(
                            probandReads, mean(dat_temp[,probandInd]))
                    }
                    else{
                        probandReads<-c(probandReads,0)
                    }
                }
                else{
                    probandReads<-c(probandReads, "-")
                }
            }
            else{
                probandReads<-c(probandReads,"-")
            }   
            gene<-c(gene,as.character(gnsOverlapID[ki]))
        }
        if(is.na(probandInd[1])==FALSE){
            ProbandGenes<-c()
            for(ii in 1:length(gene)){
                pasgene<-paste(gene[ii],"(", probandReads[ii],")", sep = "")
                ProbandGenes<-c(ProbandGenes,pasgene)
            }
            ProbandTPM<-paste(ProbandGenes,collapse=";")
            } else{
                ProbandTPM <- "-"
            }
        datGeneInfo<-data.frame(SVID=SVID, ProbandTPM=ProbandTPM)     
                
    }
    else if(length(gnsOverlapID)==1){
        pasgnovlap<-paste("^",as.character(gnsOverlapID),"$",sep="")
        #print(ki)
        gg<-grep(pasgnovlap,pasgnsname,fixed=TRUE)
      
        dat_temp<-RNASeqData[gg,]
        if(nrow(dat_temp)>1){
            dat_temp_1 <- mean(dat_temp[,probandInd])
            #dat_temp_1<-data.frame(dat_temp_1)
            dat_temp1<-dat_temp[1,]
            dat_temp1[,probandInd]<-dat_temp_1
            #dat_temp1<-cbind(dat_temp[1,1:3],dat_temp1)
            #print(dim(dat_temp1))
            probandcount<-c()
            if(is.na(probandInd[1])==FALSE){ 
                if(length(probandInd)>1){
                    for(j in probandInd){
                        probandcount <- c(probandcount,mean(dat_temp1[,j]))
                    }
                probandReads <- paste(probandcount,collapse = ":")
                }
                else if(length(probandInd)==1){
                    probandReads<-mean(dat_temp1[,probandInd])
                }
                else{
                    probandReads <- 0
                }
            } 
            else{
                probandReads <- "-"
            }
          
          
        #motherReads<-dat_temp1[,motherInd]
        #probandReads<-dat_temp1[,probandInd]
        #if(is.na(siblingInd[1])==TRUE){
        #siblingReads<-"-"
        #}
        #else{
        #siblingReads<-dat_temp1[,siblingInd]
        #}
        }
        else if (nrow(dat_temp)==1){
       
            #print(dim(dat_temp1))
            if(is.na(probandInd[1])==FALSE){ 
                if(length(probandInd)>1){
                    for(j in probandInd){
                        probandcount <- c(probandcount,mean(dat_temp[,j]))
                    }
                    probandReads<-paste(probandcount,collapse=":")
                }
                else if(length(probandInd)==1){
                    probandReads<-mean(dat_temp[,probandInd])
                }
                else{
                    probandReads<-0
                }
            } else{
                probandReads<-"-"
            }
        }
        else{
            probandReads<-"-"
        }
        #gene<-overlap_ensemblgenes$SYMBOL
        if(is.na(probandInd[1])==FALSE){
            gene<- as.character(gnsOverlapID)
        
            ProbandGenes<-paste(gene,"(",probandReads,")",sep="")
            #ProbandGenes<-c(ProbandGenes,pasgene)
            ProbandTPM<-as.character(ProbandGenes)
        }else{
            ProbandTPM <- "-"
        }
        datGeneInfo<-data.frame(SVID = SVID,ProbandTPM = ProbandTPM)
            
    }
    else{
        datGeneInfo <- data.frame(SVID = SVID,ProbandTPM = "-")
    }
  
  #print(warnings())
  
    return(datGeneInfo)
}
#' Annotating the Non-Overlapping genes with RNAseq expression
#'
#' @param gnsNonOverlap  character. Vector containing non-overlapping genes.
#' @param SVID  character. SV Index ID.
#' @param RNASeqData  dataFrame. RNAseq data with gene names. 
#' @param pattern_Proband  character. Pattern for proband. 
#' @return Dataframe containing TPM read counts of overlapping genes.
#' @examples
#' RNASeqDir = system.file("extdata", package="nanotatoR")
#' returnMethod="dataFrame"
#' datRNASeq <- RNAseqcombine_solo(RNASeqDir = RNASeqDir,
#' returnMethod = returnMethod)
#' gnsNonOverlap <- c("DDX11L1", "MIR1302-2HG", "OR4G4P")
#' SVID = 397
#' datgnnonovrlap <- nonOverlapRNAseq_solo(gnsNonOverlap = gnsNonOverlap, 
#' SVID = SVID, RNASeqData = datRNASeq,
#' pattern_Proband = "*_P_*")
#' @importFrom stats na.omit  
#' @export

nonOverlapRNAseq_solo<-function(gnsNonOverlap,SVID,RNASeqData,
                pattern_Proband = NA){
  ## annotation
  ###Checking if the input is empty; else if not empty add 
  ###expression values for each genes
    datGeneInfo<-data.frame()
    SVID=SVID
    ###Extracting the index for the the parents
    'if(is.na(pattern_Father)==FALSE){
    fatherInd<-grep(pattern_Father,names(RNASeqData))
    } else{fatherInd <- NA}
    if(is.na(pattern_Mother)==FALSE){
    motherInd<-grep(pattern_Mother,names(RNASeqData))
    } else{motherInd <- NA}'
    if(is.na(pattern_Proband)==FALSE){
         probandInd<-grep(pattern_Proband,names(RNASeqData))
    } else{probandInd <- NA}
    ##Checking for sibling
    'if(is.na(pattern_Sibling)==FALSE){
      siblingInd<-grep(pattern_Sibling,names(RNASeqData))
    } 
    else{
      siblingInd<-NA
    }'
    
    gene<-c()
    gnsname<-as.character(RNASeqData$GeneName)
    pasgnsname<-pasgnovlap<-paste("^",as.character(gnsname),"$",sep="")
    'nonoverlap_ensemblgenes = select(EnsDb.Hsapiens.v79, gnsNonOverlap, 
                c("GENEID","GENENAME"), "SYMBOL")
    gnsnonOverlapID<-as.character(nonoverlap_ensemblgenes$GENEID)'    
    ###Extracting Reads
    ###
    ###Genes Names Extraction
    gnsnonOverlapID<- as.character(gnsNonOverlap)
    if(length(gnsnonOverlapID)>1){
        #datGeneInfoTemp<-data.frame()
        probandReads<-c()
      
      
        for (ki in 1:length(gnsnonOverlapID)){
            
            pasgnnonovlap <- paste("^",as.character(
                gnsnonOverlapID[ki]),"$",sep = "")
            gg<-grep(pasgnnonovlap,pasgnsname,fixed=TRUE)
            dat_temp<-RNASeqData[gg,]
            if(nrow(dat_temp)>1){
                dat_temp_1 <- mean(dat_temp[,probandInd])
                #dat_temp_1<-data.frame(dat_temp_1)
                dat_temp1<-dat_temp[1,]
                dat_temp1[,probandInd]<-dat_temp_1
                #print(dim(dat_temp1))
                probandcount<-c()
                if(is.na(probandInd[1]) == FALSE){
                    if(length(probandInd)>1){
                        for(j in probandInd){
                            probandcount <- c(probandcount,dat_temp1[,j])
                        }
                    probandReads<-c(
                        probandReads,paste(probandcount,collapse=":")
                        )
                    }
                    else if(length(probandInd)==1){
                        probandReads<-c(probandReads,dat_temp1[,probandInd])
                    }
                    else{
                        probandReads<-c(probandReads,0)
                    }
                } else{
                    probandReads<-c(probandReads, "-")
                }
            }
            else if (nrow(dat_temp)==1) {
                #print(dim(dat_temp1))
                probandcount<-c()
                if(is.na(probandInd[1]) == FALSE){ 
                    if(length(probandInd)>1){
                        for(j in probandInd){
                            probandcount <- c(probandcount, dat_temp[,j])
                        }
                    probandReads<-c(probandReads,paste(probandcount,collapse=":"))
                }
                else if(length(probandInd)==1){
                    probandReads<-c(probandReads,dat_temp[,probandInd])
                }
                else{
                    probandReads<-c(probandReads,0)
                }
            }
            else{
                probandReads <- c(probandReads,"-")
            }
        }
        else{
            probandReads<-c(probandReads,"-")
        }
        gene<-c(gene,as.character(gnsnonOverlapID[ki]))
    }
        if(is.na(probandInd[1])==FALSE){
            ProbandGenes<-c()
        
            for(ii in 1:length(gene)){
                pasgene<-paste(gene[ii],"(",probandReads[ii],")",sep="")
                ProbandGenes<-c(ProbandGenes,pasgene)
            }
            ProbandTPM<-paste(ProbandGenes,collapse=";")
            } else{
                ProbandTPM<-"-"
            }
        datGeneInfo<-data.frame(SVID=SVID,ProbandTPM=ProbandTPM)
    } 
    else if(length(gnsnonOverlapID)==1){
        pasgnnonovlap<-paste("^",as.character(gnsnonOverlapID),"$",sep="")
        gg<-grep(pasgnnonovlap,pasgnsname,fixed=TRUE)
        #gg<-grep(pasgnnonovlap,gnsname,fixed=TRUE)
        dat_temp<-RNASeqData[gg,]
      
        if(nrow(dat_temp)>1){
            dat_temp_1 <- mean(dat_temp[,probandInd])
            #dat_temp_1<-data.frame(dat_temp_1)
            dat_temp1<-dat_temp[1,]
            dat_temp1[, probandInd]<-dat_temp_1
           #print(dim(dat_temp1))
            probandcount<-c()
            if(is.na(probandInd[1])==FALSE){
                if(length(probandInd)>1){
                    for(j in probandInd){
                        probandcount<-c(probandcount,dat_temp1[,j])
                    }
                probandReads<-paste(probandcount,collapse=":")
            }
            else if(length(probandInd)==1){
                probandReads<-dat_temp1[,probandInd]
            }
            else{
                probandReads<-0
            }
        }
        else{
            probandReads<- "-"
        }
          
          
        }
        else if (nrow(dat_temp)==1){
       
         #print(dim(dat_temp1))
            probandcount<-c()
            if(is.na(probandInd)==FALSE){
                if(length(probandInd)>1){
                    for(j in probandInd){
                        probandcount<-c(probandcount,dat_temp[,j])
                    }
                probandReads <- paste(probandcount,collapse=":")
            }
            else if(length(probandInd)==1){
                probandReads <- dat_temp[,probandInd]
            }
            else{
                probandReads<-0
            }
        }
        else{
            probandReads<- "-"
        }
        }
        else{
            probandReads<-"-"
        }
        gene<-c(gene,as.character(gnsnonOverlapID))
        if(is.na(probandInd[1])==FALSE){
            ProbandGenes <- c()
            ProbandGenes <- paste(gene,"(",probandReads,")",sep="")
            #ProbandGenes<-c(ProbandGenes,pasgene)
        #}
            ProbandTPM<-as.character(ProbandGenes)
        }else{
            ProbandTPM<-"-"
        }
        datGeneInfo <- data.frame(SVID = SVID,ProbandTPM = ProbandTPM)
    }
    else{
        datGeneInfo < -data.frame(SVID = SVID,ProbandTPM = "-")
    
    }
    return(datGeneInfo)
}
#' Annotating the Overlapping and Non-Overlapping genes with RNAseq expression
#'
#' @param input_fmt_SV  character. Input format of the SV data.Options 
#' "Text" or "DataFrame".
#' @param smapdata  dataframe. SV data dataframe.
#' @param smappath  character. smap path.
#' @param input_fmt_RNASeq  character. Input format of  
#' the RNASeq data. Options "Text" or "DataFrame"..
#' @param RNASeqData  dataFrame. RNAseq data with gene names.
#' @param RNASeqPATH  character. RNAseq dataset path . 
#' @param outputfmt  character. Output format of  
#' the result. Options "Text" or "DataFrame"..
#' @param pattern_Proband  character. Pattern for proband. 
#' @param EnzymeType  character. Enzyme used. option "SVMerge" or "SE".
#' @return Dataframe Annotated datafreme with RNASeq data.
#' @examples
#' RNASeqDir = system.file("extdata", package="nanotatoR")
#' returnMethod="dataFrame"
#' datRNASeq <- RNAseqcombine_solo(RNASeqDir = RNASeqDir,
#' returnMethod = returnMethod)
#' smapName="NA12878_DLE1_VAP_solo5.smap"
#' smap = system.file("extdata", smapName, package="nanotatoR")
#' smap = system.file("extdata", smapName, package="nanotatoR")
#' bedFile <- system.file("extdata", "HomoSapienGRCH19_lift37.bed", package="nanotatoR")
#' outpath <- system.file("extdata", package="nanotatoR")
#' datcomp<-overlapnearestgeneSearch(smap = smap, 
#'     bed=bedFile, inputfmtBed = "bed", outpath, 
#'     n = 3, returnMethod_bedcomp = c("dataFrame"), 
#'     input_fmt_SV = "Text",
#'     EnzymeType = "SE", 
#'     bperrorindel = 3000, 
#'     bperrorinvtrans = 50000)
#' datRNASeq1 <- SVexpression_solo (input_fmt_SV=c("dataFrame"),
#'     smapdata = datcomp,
#'     input_fmt_RNASeq=c("dataFrame"),
#'     RNASeqData = datRNASeq,
#'     outputfmt=c("datFrame"),
#'     pattern_Proband = "*_P_*", EnzymeType = c("SE"))
#' datRNASeq1[1,]
#' @importFrom stats na.omit  
#' @export
SVexpression_solo <- function(input_fmt_SV=c("Text","dataFrame"),
        smapdata,smappath, input_fmt_RNASeq=c("Text","dataFrame"),
        RNASeqData,RNASeqPATH,outputfmt=c("Text","datFrame"),
        pattern_Proband = NA, EnzymeType = c("SVMerge", "SE")){
  ###RNASEQ Analysis data
    if(input_fmt_RNASeq=="dataFrame"){
        RNASeqData = RNASeqData
    }
    else if(input_fmt_RNASeq=="Text"){
        RNASeqData=read.csv(RNASeqPATH)
    }
    else{
        stop("Input format for RNASeq Data Incorrect")
    }
  
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
            smapdata <- readSMap(smappath, input_fmt_smap = "Text")
            SVID<-smapdata$SVIndex
        }
        else{
            smapdata <- readSMap_DLE(smappath, input_fmt_smap = "Text")
            SVID<-smapdata$SmapEntryID
        }
    }
    else{
        stop("Input format for SMAP Incorrect")
    }
    ##Extracting Data
    overlapgenes <- stringr::str_trim(smapdata$OverlapGenes_strand_perc)
  
    
    dataOverLap <- data.frame(matrix(nrow = nrow(smapdata), ncol = 2))
    ##Extracting Overlapped Genes
    #dataOverLap<-data.frame(matrix(nrow=10,ncol=5))
    names(dataOverLap)<-c("SVID", "OverlapProbandTPM")
    print("###OverlapGenes###")
    for(kk in 1:length(overlapgenes)){
        #print(paste("kk:",kk))
        #for(kk in 1:10){
        #print(kk)
        datOverLap<-data.frame()
        #print(paste("kk:",kk,sep=""))
        svID<-as.character(SVID[kk])
        if(length(grep(";",overlapgenes[kk]))>=1){
            st1<-strsplit(as.character(overlapgenes[kk]), split = ";")
            sttemp<-as.character(st1[[1]])
            #print("1")
            gns_overlap<-c()
            for (tt in 1:length(sttemp)){
                gn_temp<-strsplit(sttemp[tt], split = "\\(")
                gns_overlap<-c(gns_overlap,as.character(gn_temp[[1]][1]))
            }
            
                datOverLap<-OverlapRNAseq_solo(
                        gnsOverlap = as.character(gns_overlap),
                        SVID = svID, 
                        RNASeqData = RNASeqData,
                        pattern_Proband = pattern_Proband)
                
        }
        else if (length(grep("\\(", as.character(overlapgenes[kk]))) >= 1){
            #print("2")
            gnsOverlap<-strsplit(as.character(overlapgenes[kk]),split="\\(")[[1]][1]
                datOverLap<-OverlapRNAseq_solo(gnsOverlap = as.character(gnsOverlap),
                    SVID = svID,
                    RNASeqData = RNASeqData,
                    pattern_Proband = pattern_Proband)
            }
            else{
                #print(paste("OverLapDNSVID:",svID))
                datOverLap<-data.frame(
                    SVID = svID, ProbandTPM = "-")
                }
                dataOverLap[kk,]<-c(
                    as.character(datOverLap$SVID),
                    Proband_OverlapGeneExpression_TPM = as.character(datOverLap$ProbandTPM))
        }
    
            ##Extracting NonOverlapped Genes
            nearestUPGenes<-smapdata$Upstream_nonOverlapGenes_dist_kb
            #datanonOverLapUP<-data.frame(matrix(nrow=nrow(smapdata),ncol=5))
            datanonOverLapUP<-data.frame(matrix(nrow = nrow(smapdata),ncol = 2))
            names(datanonOverLapUP) <- c("SVID", "NonOverlapUPProbandTPM")
            print("###NonOverlapUPStreamGenes###") 
            for(ll in 1:length(nearestUPGenes))
            {
                #for(ll in 1:10){
                #print(ll)
                datNonOverLapUP <- data.frame()
                #print(paste("llUP:",ll,sep=""))
                svID<-as.character(SVID[ll])
                if(length(grep(";",nearestUPGenes[ll]))>=1){
                    st1<-strsplit(
                        as.character(nearestUPGenes[ll]),
                        split = ";")
                        sttemp<-as.character(st1[[1]])
                        #print("1")
                        gns_nonoverlap_up<-c()
                    for (mm in 1:length(sttemp)){
                        gn_temp<-strsplit(sttemp[mm],split="\\(")
                        gns_nonoverlap_up<-c(gns_nonoverlap_up, 
                            as.character(gn_temp[[1]][1]))
                    }
                
                    datNonOverLapUP<-nonOverlapRNAseq_solo(
                        gnsNonOverlap = as.character(gns_nonoverlap_up),
                        SVID = svID,
                        RNASeqData = RNASeqData,
                        pattern_Proband = pattern_Proband)
                    
                }
                else if (length(grep(
                "\\(",as.character(nearestUPGenes[ll]))) >=1
                ){
            #print("2")
                        gnsNonOverlapUP<-strsplit(as.character(nearestUPGenes[ll]),split="\\(")[[1]][1]
                        
                        datNonOverLapUP<-nonOverlapRNAseq_solo(
                            gnsNonOverlap = as.character(gnsNonOverlapUP),
                            SVID = svID,RNASeqData = RNASeqData,
                            pattern_Proband=pattern_Proband)
                        
                    }
                else{
                    #print(paste("NonOverLapUPSVID:",svID))
                        datNonOverLapUP<-data.frame(
                            SVID = svID, ProbandTPM="-"
                        )
        }
        datanonOverLapUP[ll,]<-c(
            as.character(datNonOverLapUP$SVID), 
            Proband_Upstream_nonOverlapGeneExpression_TPM = as.character(datNonOverLapUP$ProbandTPM))
                
    }
    
  ##Extracting NonOverlapped Down Stream Genes
    nearestDNGenes<-smapdata$Downstream_nonOverlapGenes_dist_kb
    datanonOverLapDN<-data.frame(matrix(nrow = nrow(smapdata), ncol = 2))
    names(datanonOverLapDN)<-c("SVID","NonOverlapDNProbandTPM")
    print("###NonOverlapDNStreamGenes###") 
    for(nn in 1:length(nearestDNGenes))
    {
        #for(nn in 1:10){
        datNonOverLapDN<-data.frame()
        # print(paste("llDN:",ll,sep=""))
        svID<-as.character(SVID[nn])
        if(length(grep(";",nearestDNGenes[nn]))>=1){
            st1 <- strsplit(as.character(nearestDNGenes[nn]), split = ";")
            sttemp<-as.character(st1[[1]])
            #print("1")
            gns_nonoverlap_dn<-c()
            for (mm in 1:length(sttemp)){
                gn_temp <- strsplit(sttemp[mm],split="\\(")
                gns_nonoverlap_dn <- c(
                        gns_nonoverlap_dn,as.character(gn_temp[[1]][1])
                        )
            }
            datNonOverLapDN<-nonOverlapRNAseq_solo(
                    gnsNonOverlap = as.character(gns_nonoverlap_dn),
                    SVID = svID,RNASeqData = RNASeqData,
                    pattern_Proband = pattern_Proband)
            
        }
        else if (length(grep("\\(",as.character(nearestDNGenes[nn]))) >= 1){
            # print("2")
            gnsNonOverlapDN<-strsplit(as.character(
                nearestDNGenes[nn]),split="\\(")[[1]][1]
                datNonOverLapDN<-nonOverlapRNAseq_solo(
                        gnsNonOverlap = as.character(gnsNonOverlapDN),
                        SVID = svID,
                        RNASeqData = RNASeqData,
                        pattern_Proband = pattern_Proband)
                
        }
        else{
            #print(paste("NonOverLapDNSVID:",svID))
            #print ("SVID")
            datNonOverLapDN<-data.frame(
            SVID = svID,ProbandTPM="-")
        }
            datanonOverLapDN[nn,]<-c(as.character(datNonOverLapDN$SVID),
                Proband_Downstream_nonOverlapGeneExpression_TPM = as.character(datNonOverLapDN$ProbandTPM))
            }
  
            dataFinal<-data.frame(smapdata, 
                OverlapProbandEXP = dataOverLap[,2],
                NonOverlapUPprobandEXP = datanonOverLapUP[,2],
                NonOverlapDNprobandEXP = datanonOverLapDN[,2])
    return(dataFinal)

}