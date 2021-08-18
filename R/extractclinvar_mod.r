#' Extract the genes and variants related to a genetic disorder from ClinVar
#'
#' @param keyword  character. character string: keyword,
#' to search a disease, a clinical feature, or a phenotype.
#' @param localPDB.path  character. the path of
#' localized public data bases. The default value is set in the working
#' directory.
#' @param type  character. the type of the information to extract,
#' must be one of "gene", "variant","both"(default).
#' @param HPO.disease  character. MIM number of the disease.
#' The default value is NULL, which means that all the
#' OMIM number of the disease in HPO are added.
#' localized public data bases. The default value is set in the working
#' directory.
#' @param genelist  character. the gene(s) associated
#' to the disease, or the genes you are interested.
#' @param OMIM  character. whether use the information from
#' OMIM database. The default value is NULL.
#' It can be set 'yes' when you make sue you have a OMIM API key.
#' @return subset of the file gene_condition_source_id, which include
#' all the information about genes and phenotypes in ClinVar and
#' subset of the file variant_summary.txt,
#' but added sevetal colomns which describe the phenotype
#' from GeneReview, MedGen, and OMIM databases. Function modified from
#' extract_clinvar function VarFromPDB.
#' @examples
#' keyword = "retinoblastoma"
#' extract_clinvar_mod(keyword,
#' localPDB.path = system.file("extdata", "localPDB", package="nanotatoR"),
#' type = "both", HPO.disease = NULL,
#' genelist = NULL, OMIM = NULL)
#' @importFrom stats na.omit
#' @importFrom curl curl_download
#' @import XML2R
#' @import XML
#' @import stringr
#' @export
extract_clinvar_mod <- function(keyword,
    localPDB.path,
    type="both",
    HPO.disease = NULL, genelist = NULL, OMIM = NULL){
    if( !is.null(OMIM)) {
        morbidmap=paste(localPDB.path,"morbidmap.txt",sep="/")
        morbidmap <- read.delim(morbidmap,comment.char = "#")
        colnames(morbidmap) <- c("disease","gene","gene.mim.no","location")
    }
    print(localPDB.path)
    localPDB = localPDB.path
    if(!is.null(localPDB)){
        print("***localpdb exists***")
        print(localPDB.path)
        if(file.exists(paste(localPDB.path,"variant_summary.txt.gz",sep="/"))){
            print("***Clinvar File Present***")
            clinvar <- paste(localPDB.path,"variant_summary.txt.gz",sep="/")
            }else{
                clinvar <- NULL
        }

        if(file.exists(paste(localPDB.path,"gene_condition_source_id",sep="/"))){
            print("***gene2dis File Present***")
            gene2dis <- paste(localPDB.path,"gene_condition_source_id",sep="/")
            }else{
                gene2dis <- NULL
        }


    }else{
        print("***localpdb doesn't exists***")
        print(localPDB.path)
        clinvar <- NULL
        gene2dis <- NULL
    }

##check clinvar database, download the files if missing
    if(is.null(clinvar)){
        clinvar <- "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
        download.path <- localPDB.path
        clinvar.local <- paste(download.path,"variant_summary.txt.gz",sep="/")
        if(!file.exists(download.path))
            dir.create(download.path )
            options(timeout = 5000)
            if( !file.exists(clinvar.local))
                curl_download(clinvar, clinvar.local)
                clinvar <- paste(download.path,
                "variant_summary.txt.gz",sep="/")
    }

    if(is.null(gene2dis)){
        gene2dis <- "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id"
        #download.path <- paste(getwd(),"localPDB",sep="/")
        download.path <- localPDB.path
        gene2dis.local <- paste(download.path,"gene_condition_source_id",sep="/")
        if(!file.exists(download.path))
            dir.create(download.path )
        options(timeout = 5000)
        if(!file.exists(gene2dis.local) )
            curl_download(gene2dis,gene2dis.local)
        gene2dis <- paste(download.path,"gene_condition_source_id",sep="/")
    }

## input the summary file
    if(substr(clinvar,nchar(clinvar)-1,nchar(clinvar)) == "gz"){
        clinvar <- read.delim(gzfile(clinvar))
        }else{
            clinvar <- read.delim(clinvar)
    }

    gene2dis <- read.delim(gene2dis)


## HPO, the phynotypes maybe have diffrent names in HPO
        if(is.null(HPO.disease)){
            if(file.exists(paste(
                localPDB.path,"phenotype_to_genes.txt",sep="/"))){
                    HPO.disease.check <- phenoextractHPO_mod(keyword= keyword,
                    localPDB.path = localPDB.path)
            }else{
                HPO.disease.check <- phenoextractHPO_mod(keyword= keyword)
            }
            HPO.disease <- as.character(
                unique(HPO.disease.check[grep("OMIM",HPO.disease.check[,1]),1])
                )
       }

##----------------------------------------------
##begin to search
        if(!is.null(keyword)){
            gene2dis.d <- gene2dis[grep_split(keyword,gene2dis[,"DiseaseName"]),]
            pheno.yes <- as.character(gene2dis.d[,"DiseaseName"])
           }else if((is.null(keyword)&
                !is.null(HPO.disease))
                | (is.null(keyword)
                & !is.null(genelist))){
                gene2dis.d <- c()
                pheno.yes <- c()
            }else{
                gene2dis.d <- gene2dis
                pheno.yes <- c()
            }

        if(!is.null(HPO.disease)){
            HPO.disease.no <- unlist(lapply(
                HPO.disease,function(x) unlist(strsplit(x,"OMIM:"))[2]))
            gene2dis.d2 <- gene2dis[is.element(gene2dis[,"DiseaseMIM"],HPO.disease.no),]
            pheno.yes2 <- as.character(gene2dis.d2[,"DiseaseName"])
            gene2dis.d <- rbind(gene2dis.d,gene2dis.d2)
            pheno.yes <- union(pheno.yes,pheno.yes2)
       }

   #for a given genelist
        if(!is.null(genelist)){
            gene2dis.d3 <- gene2dis[is.element(gene2dis[,"AssociatedGenes"],genelist),]
            gene2dis.d <- rbind(gene2dis.d,gene2dis.d3)
       }

        gene2dis.extr <- unique(gene2dis.d)
        if(nrow(gene2dis.extr) > 0){
            gene2dis.extr$pheno.check <- "no"
            gene2dis.extr[is.element(gene2dis.extr$DiseaseName,pheno.yes),"pheno.check"] <-  "yes"

            genes <- unique(as.character(gene2dis.extr[,2]))

    ##extract the variants in the genes
         clinvar.extr <- clinvar[is.element(clinvar[,"GeneSymbol"],genes),]

    ## extract the variants from summary file directly
          clinvar.d <- clinvar[grep_split(keyword,clinvar[,"PhenotypeList"]),]

    ## merge the variants from ClinVar and other databases
        clinvar.extr <- unique(rbind(clinvar.extr,clinvar.d))
        }else{
            clinvar.extr = NULL
      }

    extract <- list(gene2dis.extr,clinvar.extr)
    names(extract) <- c("gene2dis","variants")
    return(extract)
}
