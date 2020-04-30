#' Extract the genes related to a disease or disease alias from HPO database.
#'
#' @param keyword  character. character string: keyword,
#' to search a disease, a clinical feature, or a phenotype.
#' @param localPDB.path  character. the path of 
#' localized public data bases. The default value is set in the working 
#' directory. 
#' @return subset of HPO and 
#' extract the genes and alias for a disease(phenotype), 
#' or a clinical feature. Function modified from 
#' pheno_extract_HPO function VarFromPDB.
#' @examples
#' HPO.phenotype = phenoextractHPO_mod("retinoblastoma",
#' localPDB.path = system.file("extdata", "localPDB", package="nanotatoR"))
#' @importFrom stats na.omit 
#' @importFrom curl curl_download
#' @import XML2R
#' @import XML
#' @import stringr
#' @export
phenoextractHPO_mod <- function(keyword, 
    localPDB.path){
    download.path = localPDB.path
    print(localPDB.path)
    localPDB = localPDB.path
    if(!is.null(localPDB)){
        print("***localpdb exists***")
        print(localPDB.path)
        if(file.exists(paste(localPDB.path,"phenotype_annotation.tab",sep="/"))){
            print("***HPO exists***")
            HPO <- paste(localPDB.path,"phenotype_annotation.tab",sep="/")
            }else{
                HPO <- NULL
        }        
            if(file.exists(paste(
                localPDB.path,"phenotype_to_genes.txt",sep="/"))) {  
                print("***phenotypetogene exists***")
                diseases_to_genes <- paste(localPDB.path,
				    "phenotype_to_genes.txt",sep="/")
                }else{
                      diseases_to_genes <- NULL
            }         
        }else{
            HPO <- NULL; diseases_to_genes <- NULL  
        }     
         
    #check HPO database
    if(is.null(HPO)){
        HPO <- "http://compbio.charite.de/jenkins/job/hpo.annotations/lastSuccessfulBuild/artifact/misc/phenotype_annotation.tab"
        if(is.null(download.path)){
            download.path <- paste(getwd(),"localPDB",sep="/")
        }else{
            download.path <- paste(download.path , "localPDB",sep="/")
        }
        if(!file.exists(download.path))
            dir.create(download.path )
       options(timeout = 300)
        if( !file.exists(paste(download.path,"phenotype_annotation.tab",sep="/"))){
           curl_download(HPO,paste(download.path,"phenotype_annotation.tab",sep="/"))
		}else{download.path = download.path}
       HPO <- paste(download.path,"phenotype_annotation.tab",sep="/")
    }
     HPO <- read.delim(HPO,header= FALSE)
     
     #input disease2gene dataset
    if(is.null(diseases_to_genes)){
        diseases_to_genes <- "http://compbio.charite.de/jenkins/job/hpo.annotations/lastSuccessfulBuild/artifact/util/annotation/phenotype_to_genes.txt"
        if(is.null(download.path)){
            download.path <- paste(getwd(),"localPDB",sep="/")
        }else{
            download.path <- paste(download.path , "localPDB",sep="/")
        }
        if(!file.exists(download.path))
          dir.create(download.path )
        options(timeout = 300)
        if( !file.exists(paste(download.path,"phenotype_to_genes.txt",sep="/"))){
           curl_download(diseases_to_genes,paste(download.path,"phenotype_to_genes.txt",sep="/"))
        }else{download.path = download.path}
       diseases_to_genes <- paste(download.path,"phenotype_to_genes.txt",sep="/"
            )
    }
    #diseases_to_genes <- read.delim(diseases_to_genes,fill= TRUE , flush= TRUE,header=FALSE,col.names=c("diseaseId","geneID","GeneSymbol"),comment.char = "#")
    diseases_to_genes_temp <- read.delim(diseases_to_genes,
        fill= TRUE , flush= TRUE,header=FALSE,
        col.names=c("HPOID","HPO_label","entrezgeneid", 
        "entrezgenesymbol", "addInfo", "Source", 
        "diseaseId"),comment.char = "#")
    diseaseId <- diseases_to_genes_temp$diseaseId
    geneID <- diseases_to_genes_temp$entrezgeneid
    GeneSymbol <- diseases_to_genes_temp$entrezgenesymbol
    diseases_to_genes <- data.frame(diseaseId, geneID, GeneSymbol)
    if(!is.null(keyword)){ 
        HPO.merge <- unlist(apply(HPO,1,function(x) paste(as.character(x),collapse="_")))
        HPO.j <- HPO[grep_split(keyword,HPO.merge),]
        
        HPO.j.sim <-  unique(HPO.j[,c(1:3,6,12)])
        colnames(HPO.j.sim) <- c("Database","ID","DiseaseName","reference","Synonym")
            
        db.id <- paste(HPO.j.sim[,1],HPO.j.sim[,2],sep=":")
        diseases_to_genes.j <- diseases_to_genes[is.element(diseases_to_genes[,1],db.id),]
        colnames(diseases_to_genes.j) <- c("DiseaseID","GeneID","GeneName")
        if( nrow(diseases_to_genes.j) >0 ){
            diseases_to_genes.j$Synonym <- diseases_to_genes.j$DiseaseName <-  ""
            for(i in unique(diseases_to_genes.j[,1])){
                diseases_to_genes.j[diseases_to_genes.j[,1] == i,"DiseaseName"] <- as.character(unique(HPO.j.sim[db.id==i,3]))
                diseases_to_genes.j[diseases_to_genes.j[,1] == i,"Synonym"] <- as.character(unique(HPO.j.sim[db.id==i,5]))         
            }
        }    
      }else{
        #return(list(HPO,diseases_to_genes))
        diseases_to_genes.j = NULL
    } 
      return(diseases_to_genes.j)  
}
