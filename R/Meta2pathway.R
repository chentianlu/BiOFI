#' @title Meta2pathway
#' @description
#' A function is used to predict the metabolic functional pathways according to the abundance of metabolites based on the algorithm Pathway Activity Profiling (PAPi).
#' Note: Please use compound IDs in KEGG database to replace the names of metabolites. For example, C00002 means ATP or Adenosine 5'-triphosphate in KEGG.
#' Users also can use Metaboanalyst (https://www.metaboanalyst.ca/MetaboAnalyst/upload/ConvertView.xhtml)(Pang, et al., 2022) to convert the metabolites names to compound IDs.
#' If metabolites do not have compound IDs, keep the original names.
#' See data example "Meta2pathway.eg" for details.
#'
#' @param metadf A data frame including metabolites with compound ID and abundance.
#' @return meta2path_res A data frame containing metabolites with compound ID and predicted KEGG pathways.
#' @export Meta2pathway
#' @import tidyverse
#' @importFrom  dplyr summarise_all
#' @importFrom  dplyr summarise
#' @importFrom  dplyr n
#' @importFrom  dplyr group_by
#' @examples
#' Meta2pathway_res <- Meta2pathway(metadf = Meta2pathway.eg)

Meta2pathway<-function(metadf){
  A=metadf
  metaApathdb=get0("sysdata", envir = asNamespace("BiOFI"))$metabo4path
  KEGG_MetabolitesAndPathways=metaApathdb[,c(1,3)]
  A2<-t(A)
  A3<-data.frame()
  for(i in rownames(A2)){
    m1<-KEGG_MetabolitesAndPathways[which(KEGG_MetabolitesAndPathways$Compound==i),]
    A3<-rbind(A3,m1)
  }

  A4<-as.data.frame(matrix(nrow=length(rownames(A3)),ncol = length(colnames(A2))))
  colnames(A4)=colnames(A2)
  A3<-cbind(A3,A4)
  for(i in rownames(A2)){
    for(j in 1:length(colnames(A2))){
      A3[which(A3$Compound==i),j+2]=A2[i,j]
    }
  }

  A3=A3[,-1]
  Pathway_Name <- A3[,1]
  temp_df <- A3
  # by_pathway=group_by(A3,Pathway_Name)
  # A3<-as.data.frame(summarise_all(by_pathway,sum))
  A3 <- A3 %>%
    group_by(Pathway_Name) %>%
    summarise_all(sum)
  A3<-as.data.frame(A3)
  rownames(A3)=A3[,1]
  A3=A3[,-1]

  # rank=as.data.frame(summarise(by_pathway,rank=n()))
  rank<- temp_df %>%
    group_by(Pathway_Name) %>%
    summarise(rank=n())
  rank=as.data.frame(rank)

  # by_pathway2=group_by(KEGG_MetabolitesAndPathways,Pathway_Name)
  # rank2=as.data.frame(summarise(by_pathway2,rank=n()))
  rank2<- KEGG_MetabolitesAndPathways %>%
    group_by(Pathway_Name) %>%
    summarise(rank=n())
  rank2=as.data.frame(rank2)

  value=c()
  a=c()
  for(i in 1:length(rank$Pathway_Name)){
    a[i]=rank$rank[i]/rank2[which(rank2$Pathway_Name==rank$Pathway_Name[i]),]$rank*100
    value=append(value,a[i])
  }
  rankper<-data.frame(rank$Pathway_Name,value)


  for(i in 1:length(rownames(A3))){

    A3[i,]=A3[i,]/rankper[which(rankper$rank.Pathway_Name==rownames(A3)[i]),]$value
  }

  meta2path=t(A3)
  A4=rbind(A2,A3)
  meta2path_res=t(A4)
  return(meta2path_res)
}

