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
#' @param outFileName Character Output filename.
#' @return Excel file containing the annotated SV map, tabs divided based on
#' type of SVs.
#' @examples
#' \dontrun{
#' path="Z:/Hayk's_Materials/Bionano/Projects/UDN/F1_UDN287643_Benic.Aria"
#' SVFile="F1.1_UDN287643_P_Q.S_VAP_SVmerge_trio_DGV_Int.txt"
#' fileName="F1.1_UDN287643_P_Q.S_VAP_SVmerge_trio_DGV.txt"
#' run_bionano_filter(path,SVFile,fileName)
#' }
#' @import openxlsx
#' @import hash
#' @export
run_bionano_filter<-function(input_fmt_geneList=c("Text","dataFrame"),input_fmt_svMap=c("Text","dataFrame"),
                    SVFile=NULL,svData,dat_geneList,fileName,outpath,outputFilename){
  print("###Printing in Excel###")
  Sys.setenv("R_ZIPCMD"="C:/Rtools/bin/zip.exe")
  ##GeneList Input Format
  if(input_fmt_geneList=="Text"){
  rr<-read.table(fileName,header=TRUE)
  }
  else if (input_fmt_geneList=="dataFrame"){
  rr=dat_geneList
  }
  else{
  stop("Dataframe Incorrect!!")
  }
  ha<-hash()
  .set(ha,keys=rr$Genes,values=rr$Terms)
     
              pg<-as.character(rr$Genes)
              gen<-paste("^",pg,"$",sep="")
              if(input_fmt_svMap=="Text"){
              r<-read.table(SVFile,sep="\t",header=T)
			  }
			  else if (input_fmt_svMap=="dataFrame"){
			  r<-svData
			  }
			  else{
			  stop("Input Formats Incorrect")
			  }
              ogene<-as.character(r$overlapGenes)
              nogene<-as.character(r$nearestNonOverlapGene)
              pagene<-c()
              pagene_term<-c()
              nopagene<-c()
              nopagene_term<-c()
              for( ii in 1:length(ogene)){
                dd<-grep(";",ogene[ii])
                gene_Terms<-c()
                if(length(dd)>0){
                  st<-strsplit(ogene[ii],";")
                  gen2<-as.character(unlist(st))
                  opag<-c()
                  
                  for (l in 1:length(gen2)){
                    pa<-paste("^",gen2[l],"$",sep="")
                    gg<-grep(pa,gen,fixed=TRUE)
                    if(length(gg)>0){
                      opag<-c(opag,as.character(gen2[l]))
                      val1<-hash::values(ha,keys=gen2[l])
                      gene_Terms<-c(gene_Terms,as.character(val1))
                    }
                    else{
                      opag<-c(opag,"-")
                    }
                  }
                  
                  opag<-as.character(unique(opag))
                  opagt<-as.character(unique(gene_Terms))
                  opagt<-paste(opagt,collapse=";")
                  if(length(opag)>1){
                    opagpa<-paste(opag,collapse=";")
                    ff<-grep("-",opagpa)
                    if (length(ff)>0){
                      opagpa<-gsub("-","",opagpa)
                      opagpa<-gsub(";;",";",opagpa)
                      opagpa<-gsub("^;","",opagpa)
                      pagene<-c(pagene,opagpa)
                      pagene_term<-c(pagene_term,opagt)
                    }
                    else{
                      pagene<-c(pagene,opagpa)
                      pagene_term<-c(pagene_term,opagt)
                    }
                  }
                  else{ 
                    pagene<-c(pagene,opag)
                    pagene_term<-c(pagene_term,opagt)
                  }
                  
                }
                else{
                  pa<-paste("^",ogene[ii],"$",sep="")
                  gg<-grep(pa,gen,fixed=TRUE)
                  if(length(gg)>0){
                    pagene<-c(pagene,as.character(ogene[ii]))
                    val1<-hash::values(ha,keys=ogene[ii])
                    pagene_term<-c(pagene_term,as.character(val1))
                  }
                  else{
                    pagene<-c(pagene,"-")
                    pagene_term<-c(pagene_term,"-")
                  }
                }
              }
              for( jj in 1:length(nogene)){
                dd1<-grep(";",nogene[jj])
                gene_Terms_no<-c()
                if(length(dd1)>0){
                  st1<-strsplit(nogene[jj],";")
                  gen3<-as.character(unlist(st1))
                  nopag<-c()
                  
                  for (l in 1:length(gen3)){
                    pa1<-paste("^",gen3[l],"$",sep="")
                    gg2<-grep(pa1,gen,fixed=TRUE)
                    if(length(gg2)>0){
                      nopag<-c(nopag,as.character(gen3[l]))
                      val2<-hash::values(ha,keys=gen3[l])
                      gene_Terms_no<-c(gene_Terms_no,as.character(val2))
                    }
                    else{
                      nopag<-c(nopag,"-")
                    }
                  }
                  nopag<-as.character(unique(nopag))
                  nopagt<-as.character(unique(gene_Terms_no))
                  nopagt<-paste(nopagt,collapse=";")
                  if(length(nopag)>1){
                    nopagpa<-paste(nopag,collapse=";")
                    ff<-grep("-",nopagpa)
                    if (length(ff)>0){
                      nopagpa<-gsub("-","",nopagpa)
                      nopagpa<-gsub(";;",";",nopagpa)
                      nopagpa<-gsub("^;",";",nopagpa)
                      nopagene<-c(nopagene,nopagpa)
                      nopagene_term<-c(nopagene_term,nopagt)
                    }
                    else{
                      nopagene<-c(nopagene,nopagpa)
                      nopagene_term<-c(nopagene_term,nopagt)
                    }
                  }
                  else{ 
                    nopagene<-c(nopagene,nopag)
                    nopagene_term<-c(nopagene_term,nopagt)
                  }
                  
                }
                else{
                  pa<-paste("^",nogene[jj],"$",sep="")
                  gg<-grep(pa,gen,fixed=TRUE)
                  if(length(gg)>0){
                    nopagene<-c(nopagene,as.character(nogene[jj]))
                    val2<-hash::values(ha,keys=nogene[jj])
                    nopagene_term<-c(nopagene_term,as.character(val2))
                  }
                  else{
                    nopagene<-c(nopagene,"-")
                    nopagene_term<-c(nopagene_term,"-")
                  }
                }
              }
              data<-data.frame(cbind(r,Overlap_PG=as.character(pagene),
			  Overlap_Terms=as.character(pagene_term),
			  Non_Overlap_PG=as.character(nopagene),
			  Non_Overlap_Terms=as.character(nopagene_term)))
              dat<-data[which(data$Type %in% "insertion"),]
			  dat1<-data[which(data$Type %in% "deletion"),]
			  dat3<-rbind(dat1,dat)
			  dat10<-dat3[which(((dat3$Found_in_self_BSPQI_molecules=="yes" 
			             & dat3$Found_in_self_BSSSI_molecules=="yes") | 
						(dat3$Found_in_self_BSPQI_molecules=="no" 
						& dat3$Found_in_self_BSSSI_molecules=="yes") |
						(dat3$Found_in_self_BSPQI_molecules=="yes" 
						& dat3$Found_in_self_BSSSI_molecules=="no")| 
						(dat3$Found_in_self_BSPQI_molecules=="yes"  
						&dat3$Found_in_self_BSSSI_molecules=="-") |
						(dat3$Found_in_self_BSPQI_molecules=="-" & 
						dat3$Found_in_self_BSSSI_molecules=="yes"))&
						((dat3$Found_in_parents_BSPQI_molecules=="-"
						& dat3$Found_in_parents_BSSSI_molecules=="-") |
						(dat3$Found_in_parents_BSPQI_molecules=="none" 
						& dat3$Found_in_parents_BSSSI_molecules=="none") |
						(dat3$Found_in_parents_BSPQI_molecules=="none" 
						& dat3$Found_in_parents_BSSSI_molecules=="-")|
						(dat3$Found_in_parents_BSPQI_molecules=="-" & dat3$Found_in_parents_BSSSI_molecules=="none"))),] 
			  dat11<-dat3[which(((dat3$Found_in_self_BSPQI_molecules=="yes" & dat3$Found_in_self_BSSSI_molecules=="yes") | 
						(dat3$Found_in_self_BSPQI_molecules=="no" &dat3$Found_in_self_BSSSI_molecules=="yes") |
						(dat3$Found_in_self_BSPQI_molecules=="yes" &dat3$Found_in_self_BSSSI_molecules=="no")| 
						(dat3$Found_in_self_BSPQI_molecules=="yes" &dat3$Found_in_self_BSSSI_molecules=="-") |
						(dat3$Found_in_self_BSPQI_molecules=="-" &dat3$Found_in_self_BSSSI_molecules=="yes"))&
						((dat3$Found_in_parents_BSPQI_molecules=="both" & dat3$Found_in_parents_BSSSI_molecules=="both") |
						(dat3$Found_in_parents_BSPQI_molecules=="-" & dat3$Found_in_parents_BSSSI_molecules=="both") |
						(dat3$Found_in_parents_BSPQI_molecules=="both" & dat3$Found_in_parents_BSSSI_molecules=="-"))),]
			  dat12<-dat3[which(((dat3$Found_in_self_BSPQI_molecules=="yes" & dat3$Found_in_self_BSSSI_molecules=="yes") | 
						(dat3$Found_in_self_BSPQI_molecules=="no" &dat3$Found_in_self_BSSSI_molecules=="yes") |
						(dat3$Found_in_self_BSPQI_molecules=="yes" &dat3$Found_in_self_BSSSI_molecules=="no")| 
						(dat3$Found_in_self_BSPQI_molecules=="yes" &dat3$Found_in_self_BSSSI_molecules=="-") |
						(dat3$Found_in_self_BSPQI_molecules=="-" &dat3$Found_in_self_BSSSI_molecules=="yes"))&
						((dat3$Found_in_parents_BSPQI_molecules=="mother" & dat3$Found_in_parents_BSSSI_molecules=="mother") |
						(dat3$Found_in_parents_BSPQI_molecules=="-" & dat3$Found_in_parents_BSSSI_molecules=="mother") |
						(dat3$Found_in_parents_BSPQI_molecules=="mother" & dat3$Found_in_parents_BSSSI_molecules=="-"))),]
			 dat13<-dat3[which(((dat3$Found_in_self_BSPQI_molecules=="yes" & dat3$Found_in_self_BSSSI_molecules=="yes") | 
						(dat3$Found_in_self_BSPQI_molecules=="no" &dat3$Found_in_self_BSSSI_molecules=="yes") |
						(dat3$Found_in_self_BSPQI_molecules=="yes" &dat3$Found_in_self_BSSSI_molecules=="no")| 
						(dat3$Found_in_self_BSPQI_molecules=="yes" &dat3$Found_in_self_BSSSI_molecules=="-")|
						(dat3$Found_in_self_BSPQI_molecules=="-" &dat3$Found_in_self_BSSSI_molecules=="yes"))&
						((dat3$Found_in_parents_BSPQI_molecules=="father" & dat3$Found_in_parents_BSSSI_molecules=="father") |
						(dat3$Found_in_parents_BSPQI_molecules=="-" & dat3$Found_in_parents_BSSSI_molecules=="father") |
						(dat3$Found_in_parents_BSPQI_molecules=="father" & dat3$Found_in_parents_BSSSI_molecules=="-"))),]
			dat14<-rbind(dat11,dat12,dat13)			  

gg<-grep("inversion",as.character(data$Type))
dat8<-data[gg,]
dat8<-dat8[which(((dat8$Found_in_self_BSPQI_molecules=="yes" & dat8$Found_in_self_BSSSI_molecules=="yes") | 
						(dat8$Found_in_self_BSPQI_molecules=="no" &dat8$Found_in_self_BSSSI_molecules=="yes") |
						(dat8$Found_in_self_BSPQI_molecules=="yes" &dat8$Found_in_self_BSSSI_molecules=="no")| 
						(dat8$Found_in_self_BSPQI_molecules=="yes" &dat8$Found_in_self_BSSSI_molecules=="-") |
						(dat8$Found_in_self_BSPQI_molecules=="-" & dat8$Found_in_self_BSSSI_molecules=="yes")) &
						((dat8$Fail_BSPQI_assembly_chimeric_score=="pass" & dat8$Fail_BSSSI_assembly_chimeric_score=="pass") | 
                        (dat8$Fail_BSPQI_assembly_chimeric_score=="fail" & dat8$Fail_BSSSI_assembly_chimeric_score=="pass") |
                        (dat8$Fail_BSPQI_assembly_chimeric_score=="pass" & dat8$Fail_BSSSI_assembly_chimeric_score=="fail")
                        | (dat8$Fail_BSPQI_assembly_chimeric_score=="pass" & dat8$Fail_BSSSI_assembly_chimeric_score=="-")
                        | (dat8$Fail_BSPQI_assembly_chimeric_score=="-" & dat8$Fail_BSSSI_assembly_chimeric_score=="pass"))),]
gg1<-grep("translocation",as.character(data$Type))
dat7<-data[gg1,]
dat7<-dat7[which(((dat7$Found_in_self_BSPQI_molecules=="yes" & dat7$Found_in_self_BSSSI_molecules=="yes") | 
						(dat7$Found_in_self_BSPQI_molecules=="no" &dat7$Found_in_self_BSSSI_molecules=="yes") |
						(dat7$Found_in_self_BSPQI_molecules=="yes" &dat7$Found_in_self_BSSSI_molecules=="no")| 
						(dat7$Found_in_self_BSPQI_molecules=="yes" &dat7$Found_in_self_BSSSI_molecules=="-") |
						(dat7$Found_in_self_BSPQI_molecules=="-" & dat7$Found_in_self_BSSSI_molecules=="yes")) &
						((dat7$Fail_BSPQI_assembly_chimeric_score=="pass" & dat7$Fail_BSSSI_assembly_chimeric_score=="pass") | 
                        (dat7$Fail_BSPQI_assembly_chimeric_score=="fail" & dat7$Fail_BSSSI_assembly_chimeric_score=="pass") |
                        (dat7$Fail_BSPQI_assembly_chimeric_score=="pass" & dat7$Fail_BSSSI_assembly_chimeric_score=="fail")
                        | (dat7$Fail_BSPQI_assembly_chimeric_score=="pass" & dat7$Fail_BSSSI_assembly_chimeric_score=="-")
                        | (dat7$Fail_BSPQI_assembly_chimeric_score=="-" & dat7$Fail_BSSSI_assembly_chimeric_score=="pass"))),]
dat6<-data[which(data$Type %in% "MisMatch"),]
list_of_datasets <- list("all" = data,  "indel_denovo"= dat10, "indel_both"= dat11,
                         "indel_mother"= dat12, "indel father"= dat13, "indel_cmpdHET"=dat14, "inv"= dat8, "trans"=dat7, "mismatch"=dat6,
                         "all_PG_OV"=data, "all_PG_Non_OV"= data)
dat14<-rbind(dat11,dat12,dat13)			  
fname<-paste(outputFilename,".xlsx",sep="")
write.xlsx(list_of_datasets, file = file.path(outpath,fname))
              
}





