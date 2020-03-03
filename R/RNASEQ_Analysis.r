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
#' RNASeqDir="Z:/Suro/Annotator/Data/RNAseq"
#' returnMethod="dataFrame"
#' RNAseqcombine(RNASeqDir,returnMethod)
#' @importFrom stats na.omit
#' @import AnnotationDbi
#' @import org.Hs.eg.db
#' @import EnsDb.Hsapiens.v79
#' @export
RNAseqcombine<-function(RNASeqDir,returnMethod=c("Text","dataFrame"),
                        outpath="",outFileName=""){ 
    #library(biomaRt)
    setwd(RNASeqDir)
    l <- list.files(path = RNASeqDir, pattern="*.genes.results", full.name = TRUE)
    len<-length(l)
    #dat<-listDatasets(ensembl)
    #g1<-grep("sscrofa",listDatasets(ensembl)$dataset)
    'grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", 
                    dataset="hsapiens_gene_ensembl")'
    gen<-c(); 
    ##Need to make this function dynamic
    dat<-data.frame(matrix(ncol = len,
	    nrow = nrow(r<-read.table(l[1],sep="\t",header=TRUE))))
    cnam<-c()
    for (ii in 1:length(l)){
        r <- read.table(l[ii],sep="\t",header=TRUE)
        gen <- c(gen,as.character(r$gene_id))
        dat[,ii] <- as.numeric(r$TPM)
        str <- strsplit(l[ii],split=".genes.results")
	    #print(str[[1]][1])
        cnam <- c(cnam,str[[1]][1])
    }
    #datf<-data.frame(dat)
    gen <- unique(as.character(gen))
    st <- strsplit(gen,split="[.]")
    genes <- c()
    for(k in 1:length(gen)){
        genes<-c(genes,as.character(st[[k]][1]))
    }
    data1<-dat[,1:ncol(dat)]
    names(data1)<-cnam
    genesym<-c()
    ensemblid<-c()
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
            gene3<-c(gene3,as.character(genesym[val]))
            ens<-c(ens,as.character(genes[kk]))
        }
        else{
            gene3<-c(gene3,as.character("-"))
            ens<-c(ens,as.character(genes[kk]))
        }
    }
    RNASeqDat<-data.frame(GeneName = as.character(gene3),
	    GeneID = as.character(ens),data1)
    if(returnMethod == "Text"){
        fname = file.path(outFileName, ".csv" ,sep = "")
        write.csv(RNASeqDat,file.path(outpath,fname), 
		    row.names = FALSE)
    }
    else if (returnMethod == "dataFrame"){
        return (RNASeqDat)
    }
    else{
        stop("Invalid ReturnMethod")
    }
}

#' Extract Read counts for genes that overlap SVs.
#'
#' @param gnsOverlap  character. genes that overlap SV.
#' @param SVID  character. ID of the SVs.
#' @param RNASeqData  character. Expression of the genes.
#' @param pattern_Proband  character. Pattern to identify the proband reads.
#' @param pattern_Father  character. Pattern to identify the father reads.
#' @param pattern_Mother  character. Pattern to identify the mother reads.
#' @param pattern_Sibling  character. Pattern to identify the sibling reads.
#' @param pattern_Sibling  character. Pattern to identify the sibling reads.
#' @return Text or Dataframe containing TPM read counts of genes in the family.
#' @examples
#' RNASeqDir="Z:/Suro/Annotator/Data/RNAseq"
#' returnMethod="dataFrame"
#' OverlapRNAseq(RNASeqDir,returnMethod)
#' @importFrom stats na.omit 
#' @export
OverlapRNAseq<-function(gnsOverlap, SVID, RNASeqData,
				pattern_Proband = NA, pattern_Mother = NA,
				pattern_Father = NA, pattern_Sibling = NA){
  
    ###Finding the column names
    #print(gnsOverlap);print(length(gnsOverlap))
    if(is.na(pattern_Father)==FALSE){
        fatherInd<-grep(pattern_Father,names(RNASeqData))
	} else{fatherInd<- NA}
	if(is.na(pattern_Mother)==FALSE){
        motherInd<-grep(pattern_Mother,names(RNASeqData))
	} else{motherInd<- NA}
	if(is.na(pattern_Proband)==FALSE){
        probandInd<-grep(pattern_Proband,names(RNASeqData))
	} else{probandInd<- NA}
    if(is.na(pattern_Sibling)==FALSE){
        siblingInd<-grep(pattern_Sibling,names(RNASeqData))
    } 
	else{
        siblingInd<-NA
    }
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
	        pasgnovlap<-paste("^",gnsOverlapID[ki],"$",sep="")
		    #print(ki)
            gg<-grep(pasgnovlap,pasgnsname,fixed=TRUE)
            dat_temp<-RNASeqData[gg,]
		
            if(nrow(dat_temp)>1){
                dat_temp_1<-apply(dat_temp[,4:ncol(dat_temp)],2,mean)
		        #dat_temp_1<-data.frame(dat_temp_1)
		        dat_temp1<-dat_temp[1,]
		        dat_temp1[,4:ncol(dat_temp)]<-dat_temp_1
		         #print(dim(dat_temp1))
		        fathercount<-c();mothercount<-c();
				probandcount<-c();siblingcount<-c()
                if(is.na(fatherInd[1])== FALSE){
		            if(length(fatherInd)>1){
		                for(j in fatherInd){
                            fathercount<-c(fathercount,dat_temp1[,j])
			            }
			            fatherReads<-c(fatherReads, 
					    paste(fathercount,collapse = ":"))
                    }
		            else if(length(fatherInd)==1){
                        fatherReads<-c(fatherReads,dat_temp1[,fatherInd])
			        }
                    else{
                        fatherReads<-c(fatherReads,0)
			        }
			    }
			    else{
                    fatherReads<-c(fatherReads,"-")
                }
                if(is.na(motherInd[1])==FALSE){  
                    if(length(motherInd)>1){
		                for(j in motherInd){
                            mothercount<-c(mothercount,dat_temp1[,j])
			            }
			            motherReads<-c(motherReads,paste(mothercount, 
					        collapse = ":"))
                    }
		            else if(length(motherInd)==1){
                        motherReads<-c(motherReads,dat_temp1[,motherInd])
			        }
                    else{
                        motherReads<-c(motherReads,0)
			        }
			    }
			    else{
                    motherReads<-c(motherReads,"-")
                }
		        if(is.na(probandInd[1])==FALSE){  
                    if(length(probandInd)>1){
		                for(j in probandInd){
                            probandcount<-c(probandcount,dat_temp1[,j])
			            }
			        probandReads<-c(probandReads,paste(probandcount,collapse=":"))
                }
		        else if(length(probandInd)==1){
                    probandReads<-c(probandReads,dat_temp1[,probandInd])
			    }
                else{
                    probandReads<-c(probandReads,0)
			    }
			    }else{
                    probandReads<-c(probandReads, "-")
                }
          
            if(is.na(siblingInd[1])==FALSE){
                if(length(siblingInd)>1){
		            for(j in siblingInd){
                        siblingcount<-c(siblingcount,dat_temp1[,j])
			        }
		            siblingReads<-c(siblingReads,paste(siblingcount,collapse=":"))
                }		    
                else if(length(siblingInd)==1){
              siblingReads<-c(siblingReads,dat_temp1[,siblingInd])
            }
            else{
              siblingReads<-c(siblingReads,0)
            }
        }
        else{
            siblingReads<-c(siblingReads,"-")
          }
        }
        else if (nrow(dat_temp)==1) {
          #print(dim(dat_temp1))
		  fathercount<-c();mothercount<-c();probandcount<-c();siblingcount<-c()
		  if(is.na(fatherInd)==FALSE){
          if(length(fatherInd [1])>1){
		  for(j in fatherInd){
            fathercount<-c(fathercount,mean(dat_temp[,j]))
			}
			fatherReads<-c(fatherReads,paste(fathercount,collapse=":"))
            }
		    else if(length(fatherInd)==1){
            fatherReads<-c(fatherReads,mean(dat_temp[,fatherInd]))
			}
          else{
            fatherReads<-c(fatherReads,0)
			}
			}
			else{
            fatherReads<-c(fatherReads,"-")
          }
          if(is.na(motherInd [1])==FALSE){ 
          if(length(motherInd)>1){
		  for(j in motherInd){
            mothercount<-c(mothercount,dat_temp[,j])
			}
			motherReads<-c(motherReads,paste(mothercount,collapse=":"))
            }
		    else if(length(motherInd)==1){
            motherReads<-c(motherReads,mean(dat_temp[,motherInd]))
			}
          else{
            motherReads<-c(motherReads,0)
			}
			}
			else{
            motherReads<-c(motherReads, "-")
          }
		   if(is.na(probandInd [1])==FALSE){ 
          if(length(probandInd)>1){
		  for(j in probandInd){
            probandcount<-c(probandcount,dat_temp[,j])
			}
			probandReads<-c(probandReads,paste(probandcount,collapse=":"))
            }
		    else if(length(probandInd)==1){
            probandReads<-c(probandReads,mean(dat_temp[,probandInd]))
			}
          else{
            probandReads<-c(probandReads,0)
			}
            }
			else{
            probandReads<-c(probandReads, "-")
          }
			
          if(is.na(siblingInd[1])==FALSE){
            if(length(siblingInd)>1){
		    for(j in siblingInd){
            siblingcount<-c(siblingcount,dat_temp[,j])
			}
		    siblingReads<-c(siblingReads,paste(siblingcount,collapse=":"))
            }		    
            else if(length(siblingInd)==1){
              siblingReads<-c(siblingReads,dat_temp[,siblingInd])
            }
            else{
              siblingReads<-c(siblingReads,0)
            }
          }
          else{
            siblingReads<-c(siblingReads,"-")
          }
        }
		else{
		  fatherReads<-c(fatherReads,"-")
          motherReads<-c(motherReads,"-")
          probandReads<-c(probandReads,"-")
          if(is.na(siblingInd[1])==FALSE){
            siblingReads<-c(siblingReads,"-")
          }
          else{
            siblingReads<-c(siblingReads,"-")
          }
		}
        gene<-c(gene,as.character(gnsOverlapID[ki]))
    }
	    if(is.na(probandInd[1])==FALSE){
		ProbandGenes<-c()
		for(ii in 1:length(gene)){
		pasgene<-paste(gene[ii],"(",probandReads[ii],")",sep="")
		ProbandGenes<-c(ProbandGenes,pasgene)
		}
		ProbandTPM<-paste(ProbandGenes,collapse=";")
		} else{
		       ProbandTPM <- "-"
			}
		
		if(is.na(fatherInd[1])==FALSE){
		FatherGenes<-c()
		for(ii in 1:length(gene)){
		pasgene<-paste(gene[ii],"(",fatherReads[ii],")",sep="")
		FatherGenes<-c(FatherGenes,pasgene)
		}
		FatherTPM<-paste(FatherGenes,collapse=";")
		} else{
		       FatherTPM <- "-"
			}
		if(is.na(motherInd[1])==FALSE){
		MotherGenes<-c()
		for(ii in 1:length(gene)){
		pasgene<-paste(gene[ii],"(",motherReads[ii],")",sep="")
		MotherGenes<-c(MotherGenes,pasgene)
		}
		MotherTPM<-paste(MotherGenes,collapse=";")
		} else{
		       MotherTPM <- "-"
			}
		
		
		if(is.na(siblingInd[1])==FALSE){
		siblingGenes<-c()
		for(ii in 1:length(gene)){
		
		pasgene<-paste(gene[ii],"(",siblingReads[ii],")",sep="")
		siblingGenes<-c(siblingGenes,pasgene)
		}
		SiblingTPM<-paste(siblingGenes,collapse=";")
		}
        else{
		       SiblingTPM<-"-"
			}
		       
        datGeneInfo<-data.frame(SVID=SVID,ProbandTPM=ProbandTPM,
		        FatherTPM=FatherTPM,MotherTPM=MotherTPM,
				SiblingTPM=SiblingTPM)     
        		
    }
    else if(length(gnsOverlapID)==1){
	   pasgnovlap<-paste("^",as.character(gnsOverlapID),"$",sep="")
		#print(ki)
        gg<-grep(pasgnovlap,pasgnsname,fixed=TRUE)
      #gg<-grep(gnsOverlapID,gnsname,fixed=TRUE)
      dat_temp<-RNASeqData[gg,]
     if(nrow(dat_temp)>1){
          dat_temp_1<-apply(dat_temp[,4:ncol(dat_temp)],2,mean)
		  #dat_temp_1<-data.frame(dat_temp_1)
		  dat_temp1<-dat_temp[1,]
		  dat_temp1[,4:ncol(dat_temp)]<-dat_temp_1
		  #dat_temp1<-cbind(dat_temp[1,1:3],dat_temp1)
		   #print(dim(dat_temp1))
		  fathercount<-c();mothercount<-c();probandcount<-c();siblingcount<-c()
		  if(is.na(fatherInd[1])==FALSE){
          if(length(fatherInd)>=1){
		  for(j in fatherInd){
            fathercount<-c(fathercount,mean(dat_temp1[,j]))
			}
			fatherReads<-paste(fathercount,collapse=":")
            }
		    else if(length(fatherInd)==1){
            fatherReads<-mean(dat_temp1[,fatherInd])
			}
          else{
            fatherReads<-0
			}
			} else{
            fatherReads<- "-"
            }
          if(is.na(motherInd[1])==FALSE){ 
          if(length(motherInd)>1){
		  for(j in motherInd){
            mothercount<-c(mothercount,mean(dat_temp1[,j]))
			}
			motherReads<-paste(mothercount,collapse=":")
            }
		    else if(length(motherInd)==1){
            motherReads<-mean(dat_temp1[,motherInd])
			}
          else{
            motherReads<-0
			}
			} else{
            motherReads <- "-"
            }
		  if(is.na(probandInd[1])==FALSE){ 
          if(length(probandInd)>1){
		  for(j in probandInd){
            probandcount<-c(probandcount,mean(dat_temp1[,j]))
			}
			probandReads<-paste(probandcount,collapse=":")
            }
		    else if(length(probandInd)==1){
            probandReads<-mean(dat_temp1[,probandInd])
			}
          else{
            probandReads<-0
			}
			} else{
            probandReads <- "-"
            }
          
          if(is.na(siblingInd[1])==FALSE){
            if(length(siblingInd)>1){
		    for(j in siblingInd){
            siblingcount<-c(siblingcount,mean(dat_temp1[,j]))
			}
		    siblingReads<-paste(siblingcount,collapse=":")
            }		    
            else if(length(siblingInd)==1){
              siblingReads<-mean(dat_temp1[,siblingInd])
            }
            else{
              siblingReads <- 0
            }
          }
          else{
                siblingReads<- "-"
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
		  fathercount<-c();mothercount<-c();probandcount<-c();siblingcount<-c()
		  if(is.na(fatherInd[1])==FALSE){
          if(length(fatherInd)>1){
		  for(j in fatherInd){
            fathercount<-c(fathercount,mean(dat_temp[,j]))
			}
			fatherReads<-paste(fathercount,collapse=":")
            }
		    else if(length(fatherInd)==1){
            fatherReads<-mean(dat_temp[,fatherInd])
			}
          else{
            fatherReads<-0
			}
			}else{
		        fatherReads<-"-"
		    }
          if(is.na(motherInd[1])==FALSE){  
          if(length(motherInd)>1){
		  for(j in motherInd){
            mothercount<-c(mothercount,mean(dat_temp[,j]))
			}
			motherReads<-paste(mothercount,collapse=":")
            }
		    else if(length(motherInd)==1){
            motherReads<-mean(dat_temp[,motherInd])
			}
          else{
            motherReads<-0
			}
			} else{
		        motherReads<-"-"
		    }
			if(is.na(probandInd [1])==FALSE){  
          if(length(probandInd)>1){
		  for(j in probandInd){
            probandcount<-c(probandcount,mean(dat_temp[,j]))
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
			
          
          if(is.na(siblingInd[1])==FALSE){
            if(length(siblingInd)>1){
		    for(j in siblingInd){
            siblingcount<-c(siblingcount,mean(dat_temp[,j]))
			}
		    siblingReads<-paste(siblingcount,collapse=":")
            }		    
            else if(length(siblingInd)==1){
              siblingReads<-mean(dat_temp[,siblingInd])
            }
            else{
              siblingReads<-0
            }
          }
          else{
            siblingReads<-"-"
          }
      }
	  else{
	    fatherReads<-"-"
        motherReads<-"-"
        probandReads<-"-"
		if(is.na(siblingInd[1])==TRUE){
		siblingReads<-"-"
		}
		else{
		siblingReads<-"-"
		}
	  }
      #gene<-overlap_ensemblgenes$SYMBOL
	  if(is.na(probandInd[1])==FALSE){
	    gene<- as.character(gnsOverlapID)
	    #ProbandGenes<-c()
	    ProbandGenes<-paste(gene,"(",probandReads,")",sep="")
		#ProbandGenes<-c(ProbandGenes,pasgene)
		ProbandTPM<-as.character(ProbandGenes)
		}else{
            ProbandTPM <- "-"
        }
		if(is.na(fatherInd[1])==FALSE){
		FatherGenes<-paste(gene,"(",fatherReads,")",sep="")
		FatherTPM<-as.character(FatherGenes)
		}
		else{
            FatherTPM <- "-"
        }
		if(is.na(motherInd[1])==FALSE){
		MotherGenes<-paste(gene,"(",motherReads,")",sep="")
		#MotherGenes<-c(MotherGenes,pasgene)
		#}
		MotherTPM<-as.character(MotherGenes)
		}
		else{
            MotherTPM <- "-"
        }
		if(is.na(siblingInd[1])==FALSE){
		siblingGenes<-paste(gene,"(",siblingReads,")",sep="")
		#siblingGenes<-c(siblingGenes,pasgene)
		#}
		SiblingTPM<-as.character(siblingGenes)
		}
        else{
            	SiblingTPM<-"-"
        }				
        datGeneInfo<-data.frame(SVID=SVID,ProbandTPM=ProbandTPM,
		        FatherTPM=FatherTPM,MotherTPM=MotherTPM,
				SiblingTPM=SiblingTPM)
			
    }
	else{
	    datGeneInfo<-data.frame(SVID=SVID,ProbandTPM="-",
		        FatherTPM="-",MotherTPM="-",
				SiblingTPM="-")
	
	
	
	}
  
  #print(warnings())
  
  return(datGeneInfo)
}


#' Extract Read counts for genes that are near SVs.
#'
#' @param gnsNonOverlap  character. genes that are upstream 
#' and/or downstream of SV.
#' @param SVID  character. ID of the SVs.
#' @param RNASeqData  character. Expression of the genes.
#' @param pattern_Proband  character. Pattern to identify the proband reads.
#' @param pattern_Father  character. Pattern to identify the father reads.
#' @param pattern_Mother  character. Pattern to identify the mother reads.
#' @param pattern_Sibling  character. Pattern to identify the sibling reads.
#' @param pattern_Sibling  character. Pattern to identify the sibling reads.
#' @return Text or Dataframe containing TPM read counts of genes in the family.
#' @examples
#' RNASeqDir="Z:/Suro/Annotator/Data/RNAseq"
#' returnMethod="dataFrame"
#' RNAseqcombine(RNASeqDir,returnMethod)
#' @importFrom stats na.omit 
#' @export
nonOverlapRNAseq<-function(gnsNonOverlap,SVID,RNASeqData,
				pattern_Proband=NA,pattern_Mother=NA,
				pattern_Father=NA,pattern_Sibling=NA){
  ##Biomart annotation
  ###Checking if the input is empty; else if not empty add 
  ###expression values for each genes
	datGeneInfo<-data.frame()
    SVID=SVID
    ###Extracting the index for the the parents
	if(is.na(pattern_Father)==FALSE){
	fatherInd<-grep(pattern_Father,names(RNASeqData))
	} else{fatherInd <- NA}
	if(is.na(pattern_Mother)==FALSE){
    motherInd<-grep(pattern_Mother,names(RNASeqData))
	} else{motherInd <- NA}
	if(is.na(pattern_Proband)==FALSE){
    probandInd<-grep(pattern_Proband,names(RNASeqData))
	} else{probandInd <- NA}
	##Checking for sibling
    if(is.na(pattern_Sibling)==FALSE){
      siblingInd<-grep(pattern_Sibling,names(RNASeqData))
    } 
	else{
      siblingInd<-NA
    }
	
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
      fatherReads<-c()
      motherReads<-c()
      probandReads<-c()
      siblingReads<-c()
	  
      for (ki in 1:length(gnsnonOverlapID)){
	   pasgnnonovlap<-paste("^",as.character(gnsnonOverlapID[ki]),"$",sep="")
        gg<-grep(pasgnnonovlap,pasgnsname,fixed=TRUE)
        dat_temp<-RNASeqData[gg,]
        if(nrow(dat_temp)>1){
          dat_temp_1<-apply(dat_temp[,4:ncol(dat_temp)],2,mean)
		  #dat_temp_1<-data.frame(dat_temp_1)
		  dat_temp1<-dat_temp[1,]
		  dat_temp1[,4:ncol(dat_temp)]<-dat_temp_1
		   #print(dim(dat_temp1))
		  fathercount<-c();mothercount<-c();probandcount<-c();siblingcount<-c()
          if(is.na(fatherInd[1])==FALSE){
		  if(length(fatherInd)>1){
		  for(j in fatherInd){
            fathercount<-c(fathercount,dat_temp1[,j])
			}
			fatherReads<-c(fatherReads,paste(fathercount,collapse=":"))
            }
		    else if(length(fatherInd)==1){
            fatherReads<-c(fatherReads,dat_temp1[,fatherInd])
			}
          else{
            fatherReads<-c(fatherReads,0)
			}
			}
            else{
            fatherReads<-c(fatherReads,"-")
          }
		  if(is.na(motherInd[1])==FALSE){
          if(length(motherInd)>1){
		  for(j in motherInd){
            mothercount<-c(mothercount,dat_temp1[,j])
			}
			motherReads<-c(motherReads,paste(mothercount,collapse=":"))
            }
		    else if(length(motherInd)==1){
            motherReads<-c(motherReads,dat_temp1[,motherInd])
			}
          else{
            motherReads<-c(motherReads,0)
			}
			} else{
                motherReads<-c(motherReads, "-")
            }
		  if(is.na(probandInd[1])==FALSE){
          if(length(probandInd)>1){
		  for(j in probandInd){
            probandcount<-c(probandcount,dat_temp1[,j])
			}
			probandReads<-c(probandReads,paste(probandcount,collapse=":"))
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
          if(is.na(siblingInd[1])==FALSE){
            if(length(siblingInd)>1){
		    for(j in siblingInd){
            siblingcount<-c(siblingcount,dat_temp1[,j])
			}
		    siblingReads<-c(siblingReads,paste(siblingcount,collapse=":"))
            }		    
            else if(length(siblingInd)==1){
              siblingReads<-c(siblingReads,dat_temp1[,siblingInd])
            }
            else{
              siblingReads<-c(siblingReads,0)
            }
          }
          else{
            siblingReads<-c(siblingReads,"-")
          }
        }
        else if (nrow(dat_temp)==1) {
          #print(dim(dat_temp1))
		  fathercount<-c();mothercount<-c();probandcount<-c();siblingcount<-c()
		  if(is.na(fatherInd[1])==FALSE){
          if(length(fatherInd)>1){
		  for(j in fatherInd){
            fathercount<-c(fathercount,dat_temp[,j])
			}
			fatherReads<-c(fatherReads,paste(fathercount,collapse=":"))
            }
		    else if(length(fatherInd)==1){
            fatherReads<-c(fatherReads,dat_temp[,fatherInd])
			}
          else{
            fatherReads<-c(fatherReads,0)
			}
			}
			else{
            fatherReads<-c(fatherReads,"-")
          }
			
          if(is.na(motherInd[1])==FALSE){ 
          if(length(motherInd)>1){
		  for(j in motherInd){
            mothercount<-c(mothercount,dat_temp[,j])
			}
			motherReads<-c(motherReads,paste(mothercount,collapse=":"))
            }
		    else if(length(motherInd)==1){
            motherReads<-c(motherReads,dat_temp[,motherInd])
			}
          else{
            motherReads<-c(motherReads,0)
			}
		 }
			else{
            motherReads<-c(motherReads,"-")
          }
		   if(is.na(probandInd[1])==FALSE){ 
          if(length(probandInd)>1){
		  for(j in probandInd){
            probandcount<-c(probandcount,dat_temp[,j])
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
          if(is.na(siblingInd[1])==FALSE){
            if(length(siblingInd)>1){
		    for(j in siblingInd){
            siblingcount<-c(siblingcount,dat_temp[,j])
			}
		    siblingReads<-c(siblingReads,paste(siblingcount,collapse=":"))
            }		    
            else if(length(siblingInd)==1){
              siblingReads<-c(siblingReads,dat_temp[,siblingInd])
            }
            else{
              siblingReads<-c(siblingReads,0)
            }
          }
          else{
            siblingReads<-c(siblingReads,"-")
          }
        }
		else{
		  fatherReads<-c(fatherReads,"-")
          motherReads<-c(motherReads,"-")
          probandReads<-c(probandReads,"-")
          if(is.na(siblingInd[1])==FALSE){
            siblingReads<-c(siblingReads,"-")
          }
          else{
            siblingReads<-c(siblingReads,"-")
          }
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
		 if(is.na(fatherInd[1])==FALSE){
		FatherGenes<-c()
		for(ii in 1:length(gene)){
		pasgene<-paste(gene[ii],"(",fatherReads[ii],")",sep="")
		FatherGenes<-c(FatherGenes,pasgene)
		}
		FatherTPM<-paste(FatherGenes,collapse=";")
		}else{
		       FatherTPM<-"-"
			}
		if(is.na(motherInd[1])==FALSE){
		MotherGenes<-c()
		for(ii in 1:length(gene)){
		pasgene<-paste(gene[ii],"(",motherReads[ii],")",sep="")
		MotherGenes<-c(MotherGenes,pasgene)
		}
		MotherTPM<-paste(MotherGenes,collapse=";")
		}else{
		       MotherTPM<-"-"
			}
		
		if(is.na(siblingInd[1])==FALSE){
		siblingGenes<-c()
		for(ii in 1:length(gene)){
		pasgene<-paste(gene[ii],"(",siblingReads[ii],")",sep="")
		siblingGenes<-c(siblingGenes,pasgene)
		}
		SiblingTPM<-paste(siblingGenes,collapse=";")
		}
        else{
		       SiblingTPM<-"-"
			}
        datGeneInfo<-data.frame(SVID=SVID,ProbandTPM=ProbandTPM,
		        FatherTPM=FatherTPM,MotherTPM=MotherTPM,
				SiblingTPM=SiblingTPM)      
    }	
    else if(length(gnsnonOverlapID)==1){
	pasgnnonovlap<-paste("^",as.character(gnsnonOverlapID),"$",sep="")
	gg<-grep(pasgnnonovlap,pasgnsname,fixed=TRUE)
      #gg<-grep(pasgnnonovlap,gnsname,fixed=TRUE)
      dat_temp<-RNASeqData[gg,]
	  
	  if(nrow(dat_temp)>1){
          dat_temp_1<-apply(dat_temp[,4:ncol(dat_temp)],2,mean)
		  #dat_temp_1<-data.frame(dat_temp_1)
		  dat_temp1<-dat_temp[1,]
		  dat_temp1[,4:ncol(dat_temp)]<-dat_temp_1
		   #print(dim(dat_temp1))
		  fathercount<-c();mothercount<-c();probandcount<-c();siblingcount<-c()
          if(is.na(fatherInd[1])==FALSE){
		  if(length(fatherInd)>1){
		  for(j in fatherInd){
            fathercount<-c(fathercount,dat_temp1[,j])
			}
			fatherReads<-paste(fathercount,collapse=":")
            }
		    else if(length(fatherInd)==1){
            fatherReads<-dat_temp1[,fatherInd]
			}
          else{
            fatherReads<-0
			}
             }
          else{
            fatherReads<- "-"
          }
		  if(is.na(motherInd[1])==FALSE){
          if(length(motherInd)>1){
		  for(j in motherInd){
            mothercount<-c(mothercount,dat_temp1[,j])
			}
			motherReads<-paste(mothercount,collapse=":")
            }
		    else if(length(motherInd)==1){
            motherReads<-dat_temp1[,motherInd]
			}
          else{
            motherReads<-0
			}
			 }
          else{
            motherReads<- "-"
          }
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
          
          if(is.na(siblingInd[1])==FALSE){
            if(length(siblingInd)>1){
		    for(j in siblingInd){
            siblingcount<-c(siblingcount,dat_temp1[,j])
			}
		    siblingReads<-paste(siblingcount,collapse=":")
            }		    
            else if(length(siblingInd)==1){
              siblingReads<-dat_temp1[,siblingInd]
            }
            else{
              siblingReads <- 0
            }
          }
          else{
            siblingReads<- "-"
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
		  fathercount<-c();mothercount<-c();probandcount<-c();siblingcount<-c()
		  if(is.na(fatherInd[1])==FALSE){
          if(length(fatherInd)>1){
		  for(j in fatherInd){
            fathercount<-c(fathercount,dat_temp[,j])
			}
			fatherReads<-paste(fathercount,collapse=":")
            }
		    else if(length(fatherInd)==1){
            fatherReads<-dat_temp[,fatherInd]
			}
          else{
            fatherReads<-0
			}
			}
          else{
            fatherReads<- "-"
          }
            if(is.na(motherInd[1])==FALSE){
          if(length(motherInd)>1){
		  for(j in motherInd){
            mothercount<-c(mothercount,dat_temp[,j])
			}
			motherReads<-paste(mothercount,collapse=":")
            }
		    else if(length(motherInd)==1){
            motherReads<-dat_temp[,motherInd]
			}
          else{
            motherReads<-0
			}
			}
          else{
            motherReads<- "-"
          }
		  if(is.na(probandInd)==FALSE){
          if(length(probandInd)>1){
		  for(j in probandInd){
            probandcount<-c(probandcount,dat_temp[,j])
			}
			probandReads<-paste(probandcount,collapse=":")
            }
		    else if(length(probandInd)==1){
            probandReads<-dat_temp[,probandInd]
			}
          else{
            probandReads<-0
			}
            }
          else{
            probandReads<- "-"
          }
          if(is.na(siblingInd[1])==FALSE){
            if(length(siblingInd)>1){
		    for(j in siblingInd){
            siblingcount<-c(siblingcount,dat_temp[,j])
			}
		    siblingReads<-paste(siblingcount,collapse=":")
            }		    
            else if(length(siblingInd)==1){
              siblingReads<-dat_temp[,siblingInd]
            }
            else{
              siblingReads<-0
            }
          }
          else{
            siblingReads<-"-"
          }
      }
	  else{
	    fatherReads<-"-"
        motherReads<-"-"
        probandReads<-"-"
		if(is.na(siblingInd[1])==TRUE){
		siblingReads<-"-"
		}
		else{
		siblingReads<-"-"
		}
	  }
		gene<-c(gene,as.character(gnsnonOverlapID))
		if(is.na(probandInd[1])==FALSE){
	    ProbandGenes<-c()
		ProbandGenes<-paste(gene,"(",probandReads,")",sep="")
		#ProbandGenes<-c(ProbandGenes,pasgene)
		#}
		ProbandTPM<-as.character(ProbandGenes)
		}else{
            	ProbandTPM<-"-"
        }
		if(is.na(fatherInd[1])==FALSE){
		FatherGenes<-c()
		FatherGenes<-paste(gene,"(",fatherReads,")",sep="")
		FatherTPM<-as.character(FatherGenes)
		}else{
            	FatherTPM <- "-"
        }
		if(is.na(motherInd[1])==FALSE){
		MotherGenes<-c()
		MotherGenes<-paste(gene,"(",motherReads,")",sep="")
		MotherTPM<-as.character(MotherGenes)
		}else{
            	MotherTPM <- "-"
        }
		
		if(is.na(siblingInd[1])==FALSE){
		siblingGenes<-paste(gene,"(",siblingReads,")",sep="")
		#siblingGenes<-c(siblingGenes,pasgene)
		#}
		SiblingTPM<-as.character(siblingGenes)
		}
        else{
            	SiblingTPM<-"-"
        }				
		       
        datGeneInfo<-data.frame(SVID=SVID,ProbandTPM=ProbandTPM,
		        FatherTPM=FatherTPM,MotherTPM=MotherTPM,
				SiblingTPM=SiblingTPM)
    }
	else{
	    datGeneInfo<-data.frame(SVID=SVID,ProbandTPM="-",
		        FatherTPM="-",MotherTPM="-",
				SiblingTPM="-")
	
	}
  
  return(datGeneInfo)
}

#' Extract Read counts for genes that are near 
#' or overalapping SVs.
#'
#' @param input_fmt_SV  character. genes that are upstream 
#' and/or downstream of SV.
#' @param SVID  character. ID of the SVs.
#' @param RNASeqData  character. Expression of the genes.
#' @param pattern_Proband  character. Pattern to identify the proband reads.
#' @param pattern_Father  character. Pattern to identify the father reads.
#' @param pattern_Mother  character. Pattern to identify the mother reads.
#' @param pattern_Sibling  character. Pattern to identify the sibling reads.
#' @param pattern_Sibling  character. Pattern to identify the sibling reads.
#' @param EnzymeType  character. Enzyme used. option "Dual" or "DLE".
#' @return Text or Dataframe containing TPM read counts of genes in the family.
#' @examples
#' RNASeqDir="Z:/Suro/Annotator/Data/RNAseq"
#' returnMethod="dataFrame"
#' RNAseqcombine(RNASeqDir,returnMethod)
#' @importFrom stats na.omit  
#' @export

SmapRNAseqquery<-function(input_fmt_SV=c("Text","dataFrame"),smapdata,smappath,
                          input_fmt_RNASeq=c("Text","dataFrame"),
                          RNASeqData,RNASeqPATH,outputfmt=c("Text","datFrame"),
                          pattern_Proband=NA,pattern_Mother=NA,pattern_Father=NA,
						  pattern_Sibling=NA, EnzymeType = c("Dual", "DLE")){
  ###RNASEQ Analysis data
  if(input_fmt_RNASeq=="dataFrame"){
    RNASeqData=RNASeqData
  }
  else if(input_fmt_RNASeq=="Text"){
    RNASeqData=read.csv(RNASeqPATH)
  }
  else{
    stop("Input format for RNASeq Data Incorrect")
  }
  
  if(input_fmt_SV=="dataFrame"){
    if(EnzymeType == "Dual"){
            smapdata <- readSMap(smap, input_fmt_smap)
			SVID<-smapdata$SVIndex
        }
        else{
            smapdata <- readSMap_DLE(smap, input_fmt_smap)
			SVID<-smapdata$SmapEntryID
        }
  }
  else if(input_fmt_SV=="Text"){
    smapdata=read.csv(smappath)
  }
  else{
    stop("Input format for SMAP Incorrect")
  }
  ##Extracting Data
  overlapgenes<-str_trim(smapdata$OverlapGenes_strand_perc)
  
  
  SVID<-smapdata$SVIndex
  dataOverLap<-data.frame(matrix(nrow=nrow(smapdata),ncol=5))
  ##Extracting Overlapped Genes
  #dataOverLap<-data.frame(matrix(nrow=10,ncol=5))
  names(dataOverLap)<-c("SVID","OverlapProbandTPM",
                "OverlapFatherTPM","OverlapMotherTPM",
				"OverlapSiblingTPM")
  print("###OverlapGenes###")
  for(kk in 1:length(overlapgenes)){
  #print(kk)
  #for(kk in 1:10){
    #print(kk)
    datOverLap<-data.frame()
    #print(paste("kk:",kk,sep=""))
    svID<-as.character(SVID[kk])
	if(length(grep(";",overlapgenes[kk]))>=1){
    st1<-strsplit(as.character(overlapgenes[kk]),split=";")
    sttemp<-as.character(st1[[1]])
    #print("1")
    gns_overlap<-c()
    for (tt in 1:length(sttemp)){
        gn_temp<-strsplit(sttemp[tt],split="\\(")
        gns_overlap<-c(gns_overlap,as.character(gn_temp[[1]][1]))
      }
		if(is.na(pattern_Sibling)==TRUE){
		datOverLap<-OverlapRNAseq(gnsOverlap=as.character(gns_overlap),SVID=svID,RNASeqData=RNASeqData,
	            pattern_Proband=pattern_Proband,pattern_Mother=pattern_Mother,
				pattern_Father=pattern_Father,pattern_Sibling=NA)
			}
		else{
		datOverLap<-OverlapRNAseq(gnsOverlap=as.character(gns_overlap),SVID=svID,RNASeqData=RNASeqData,
	            pattern_Proband=pattern_Proband,pattern_Mother=pattern_Mother,
				pattern_Father=pattern_Father,pattern_Sibling=pattern_Sibling)
    }		
    }
    else if (length(grep("\\(",as.character(overlapgenes[kk])))>=1){
	  #print("2")
      gnsOverlap<-strsplit(as.character(overlapgenes[kk]),split="\\(")[[1]][1]
      if(is.na(pattern_Sibling)==TRUE){
		datOverLap<-OverlapRNAseq(gnsOverlap=as.character(gnsOverlap),SVID=svID,RNASeqData=RNASeqData,
	            pattern_Proband=pattern_Proband,pattern_Mother=pattern_Mother,
				pattern_Father=pattern_Father,pattern_Sibling=NA)
			}
		else{
		datOverLap<-OverlapRNAseq(gnsOverlap=as.character(gnsOverlap),SVID=svID,RNASeqData=RNASeqData,
	            pattern_Proband=pattern_Proband,pattern_Mother=pattern_Mother,
				pattern_Father=pattern_Father,pattern_Sibling=pattern_Sibling)
    }
    }
	else{
	#print(paste("OverLapDNSVID:",svID))
	datOverLap<-data.frame(SVID=svID,ProbandTPM="-",FatherTPM="-",MotherTPM="-",SiblingTPM="-")
	    }
	dataOverLap[kk,]<-c(as.character(datOverLap$SVID),Proband_OverlapGeneExpression_TPM=as.character(datOverLap$ProbandTPM),Father_OverlapGeneExpression_TPM=as.character(datOverLap$FatherTPM),Mother_OverlapGeneExpression_TPM=as.character(datOverLap$MotherTPM),Sibling_OverlapGeneExpression_TPM=as.character(datOverLap$SiblingTPM))
	}
    
##Extracting NonOverlapped Genes
nearestUPGenes<-smapdata$Upstream_nonOverlapGenes_dist_kb
#datanonOverLapUP<-data.frame(matrix(nrow=nrow(smapdata),ncol=5))
datanonOverLapUP<-data.frame(matrix(nrow=10,ncol=5))
names(datanonOverLapUP)<-c("SVID","NonOverlapUPProbandTPM",
                        "NonOverlapUPFatherTPM","NonOverlapUPMotherTPM",
						"NonOverlapUPSiblingTPM")
print("###NonOverlapUPStreamGenes###") 
for(ll in 1:length(nearestUPGenes)){
#for(ll in 1:10){
     
    datNonOverLapUP<-data.frame()
    #print(paste("llUP:",ll,sep=""))
    svID<-as.character(SVID[ll])
	if(length(grep(";",nearestUPGenes[ll]))>=1){
    st1<-strsplit(as.character(nearestUPGenes[ll]),split=";")
    sttemp<-as.character(st1[[1]])
    #print("1")
    gns_nonoverlap_up<-c()
    for (mm in 1:length(sttemp)){
        gn_temp<-strsplit(sttemp[mm],split="\\(")
        gns_nonoverlap_up<-c(gns_nonoverlap_up,as.character(gn_temp[[1]][1]))
      }
		if(is.na(pattern_Sibling)==TRUE){
		datNonOverLapUP<-nonOverlapRNAseq(gnsNonOverlap=as.character(gns_nonoverlap_up),SVID=svID,RNASeqData=RNASeqData,
	            pattern_Proband=pattern_Proband,pattern_Mother=pattern_Mother,
				pattern_Father=pattern_Father,pattern_Sibling=NA)
			}
		else{
		datNonOverLapUP<-nonOverlapRNAseq(gnsNonOverlap=as.character(gns_nonoverlap_up),SVID=svID,RNASeqData=RNASeqData,
	            pattern_Proband=pattern_Proband,pattern_Mother=pattern_Mother,
				pattern_Father=pattern_Father,pattern_Sibling=pattern_Sibling)
    }		
    }
    else if (length(grep("\\(",as.character(nearestUPGenes[ll])))>=1){
	  #print("2")
      gnsNonOverlapUP<-strsplit(as.character(nearestUPGenes[ll]),split="\\(")[[1]][1]
      if(is.na(pattern_Sibling)==TRUE){
		datNonOverLapUP<-nonOverlapRNAseq(gnsNonOverlap=as.character(gnsNonOverlapUP),SVID=svID,RNASeqData=RNASeqData,
	            pattern_Proband=pattern_Proband,pattern_Mother=pattern_Mother,
				pattern_Father=pattern_Father,pattern_Sibling=NA)
			}
		else{
		datNonOverLapUP<-nonOverlapRNAseq(gnsNonOverlap=as.character(gnsNonOverlapUP),SVID=svID,RNASeqData=		    RNASeqData,pattern_Proband=pattern_Proband,pattern_Mother=pattern_Mother,
				 pattern_Father=pattern_Father,pattern_Sibling=pattern_Sibling)
        }
    }
	else{
	#print(paste("NonOverLapUPSVID:",svID))
	datNonOverLapUP<-data.frame(SVID=svID,ProbandTPM="-",FatherTPM="-",MotherTPM="-",SiblingTPM="-")
	    }
	datanonOverLapUP[ll,]<-c(as.character(datNonOverLapUP$SVID),Proband_Upstream_nonOverlapGeneExpression_TPM=as.character(datNonOverLapUP$ProbandTPM),Father_Upstream_nonOverlapGeneExpression_TPM=as.character(datNonOverLapUP$FatherTPM),Mother_Upstream_nonOverlapGeneExpression_TPM=as.character(datNonOverLapUP$MotherTPM),Sibling_Upstream_nonOverlapGeneExpression_TPM=as.character(datNonOverLapUP$SiblingTPM))
	}
	
  ##Extracting NonOverlapped Down Stream Genes
nearestDNGenes<-smapdata$Downstream_nonOverlapGenes_dist_kb
datanonOverLapDN<-data.frame(matrix(nrow=nrow(smapdata),ncol=5))
names(datanonOverLapDN)<-c("SVID","NonOverlapDNProbandTPM",
                        "NonOverlapDNFatherTPM","NonOverlapDNMotherTPM",
						"NonOverlapDNSiblingTPM")
  print("###NonOverlapDNStreamGenes###") 
  for(nn in 1:length(nearestDNGenes)){
  #for(nn in 1:10){
    datNonOverLapDN<-data.frame()
    # print(paste("llDN:",ll,sep=""))
    svID<-as.character(SVID[nn])
	if(length(grep(";",nearestDNGenes[nn]))>=1){
    st1<-strsplit(as.character(nearestDNGenes[nn]),split=";")
    sttemp<-as.character(st1[[1]])
    #print("1")
    gns_nonoverlap_dn<-c()
    for (mm in 1:length(sttemp)){
        gn_temp<-strsplit(sttemp[mm],split="\\(")
        gns_nonoverlap_dn<-c(gns_nonoverlap_dn,as.character(gn_temp[[1]][1]))
      }
		if(is.na(pattern_Sibling[1])==TRUE){
		datNonOverLapDN<-nonOverlapRNAseq(gnsNonOverlap=as.character(gns_nonoverlap_dn),SVID=svID,RNASeqData=RNASeqData,
	            pattern_Proband=pattern_Proband,pattern_Mother=pattern_Mother,
				pattern_Father=pattern_Father,pattern_Sibling=NA)
			}
		else{
		datNonOverLapDN<-nonOverlapRNAseq(gnsNonOverlap=as.character(gns_nonoverlap_dn),SVID=svID,RNASeqData=RNASeqData,
	            pattern_Proband=pattern_Proband,pattern_Mother=pattern_Mother,
				pattern_Father=pattern_Father,pattern_Sibling=pattern_Sibling)
    }		
    }
    else if (length(grep("\\(",as.character(nearestDNGenes[nn])))>=1){
	 # print("2")
      gnsNonOverlapDN<-strsplit(as.character(nearestDNGenes[nn]),split="\\(")[[1]][1]
      if(is.na(pattern_Sibling[1])==TRUE){
		datNonOverLapDN<-nonOverlapRNAseq(gnsNonOverlap=as.character(gnsNonOverlapDN),SVID=svID,RNASeqData=RNASeqData,
	            pattern_Proband=pattern_Proband,pattern_Mother=pattern_Mother,
				pattern_Father=pattern_Father,pattern_Sibling=NA)
			}
		else{
		datNonOverLapDN<-nonOverlapRNAseq(gnsNonOverlap=as.character(gnsNonOverlapDN),SVID=svID,RNASeqData=		    RNASeqData,pattern_Proband=pattern_Proband,pattern_Mother=pattern_Mother,
				 pattern_Father=pattern_Father,pattern_Sibling=pattern_Sibling)
        }
    }
	else{
	#print(paste("NonOverLapDNSVID:",svID))
	#print ("SVID")
	datNonOverLapDN<-data.frame(SVID=svID,ProbandTPM="-",FatherTPM="-",MotherTPM="-",SiblingTPM="-")
	    }
	datanonOverLapDN[nn,]<-c(as.character(datNonOverLapDN$SVID),Proband_Downstream_nonOverlapGeneExpression_TPM=as.character(datNonOverLapDN$ProbandTPM),Father_Downstream_nonOverlapGeneExpression_TPM=as.character(datNonOverLapDN$FatherTPM),Mother_Downstream_nonOverlapGeneExpression_TPM=as.character(datNonOverLapDN$MotherTPM),Sibling_Downstream_nonOverlapGeneExpression_TPM=as.character(datNonOverLapDN$SiblingTPM))
	}
  
  dataFinal<-data.frame(smapdata,dataOverLap[,2:ncol(dataOverLap)],
            datanonOverLapUP[,2:ncol(datanonOverLapUP)],
            datanonOverLapDN[,2:ncol(datanonOverLapDN)])
return(dataFinal)

}

	
  
  
  
    
    
    
    