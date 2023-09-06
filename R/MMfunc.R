#' @title MMfunc
#' @description A function is used to group microbes and metabolites in shared functional pathways based on different threshold setting according to functional pathway IFS,
#' microbes IFS and metabolites IFS.
#'
#' @param IFS The result returned by IFS() function.
#' @param path_IFSthre A IFS threshold value to identify the shared pathways.
#' @param mm_IFSthre A IFS threshold value to identify the microbes/metabolites affiliation with the shared pathways.
#' @return MMfunc_res A list including shared functional pathways.
#' @export MMfunc
#' @importFrom networkD3 saveNetwork
#' @importFrom networkD3 sankeyNetwork
#' @import tidyverse
#' @examples
#' IFScore <- IFS(microApath = micro.eg, metaApath = metabo.eg,
#'                conf = NULL, groupInfo = groupInfo.eg)
#' MMfunc_res <- MMfunc(IFS = IFScore)

MMfunc<-function(IFS,path_IFSthre=2.5,mm_IFSthre=2){
  mic_met_path=get0("sysdata", envir = asNamespace("BiOFI"))$mm4path
  mic_path=get0("sysdata", envir = asNamespace("BiOFI"))$micro4path
  A=as.data.frame(IFS[["IFS_completed"]])
  colnames(A)=c("nodes","ns","group")
  Metabolite=A[which(A$group=="Metabolite"& A$ns>mm_IFSthre),]$nodes
  Microbe=A[which(A$group=="Microbe"& A$ns>mm_IFSthre),]$nodes
  Metabolome_pathway=A[A$group=="Metabolome_pathway"& A$ns>path_IFSthre,]
  Microbiome_pathway=A[A$group=="Microbiome_pathway" & A$ns>path_IFSthre,]


  intersect_path <- unique(intersect(Metabolome_pathway$nodes,Microbiome_pathway$nodes))

  pathway<-rbind(subset(Metabolome_pathway,nodes %in% intersect_path),subset(Microbiome_pathway,nodes %in% intersect_path))


  metabolite2pathway<-data.frame()
  for ( i in Metabolite){
    m1<-mic_met_path[which(mic_met_path$metabolite==i),]
    m1<-m1[,c(3,4)]
    metabolite2pathway<-rbind(metabolite2pathway,m1)
  }
  colnames(metabolite2pathway)[2]='Pathway'


  meta2func=data.frame()

  for(i in pathway$nodes){
    a1=metabolite2pathway[which(metabolite2pathway$Pathway==i),]
    meta2func=rbind(meta2func,a1)
  }
  #

  # microbe and pathways ------------------------------------------------------------------

  microbename=Microbe
  microbelist=c()
  B=data.frame()
  for(i in 1:length(microbename)){
    microbe=strsplit(microbename[i],'; ')
    microbenew=vector()
    for(i in 1:length(microbe[[1]])){
      if(nchar(microbe[[1]][i])==3){
        numb=i-1
        break}else{
          numb=length(microbe[[1]])
        }}
    for(i in 1:numb){
      a=microbe[[1]][i]
      microbenew=append(microbenew,a)
    }
    microbenew=paste(microbenew,collapse ='; ')
    microbelist=append(microbelist,microbenew)
  }

  for(i in 1:length(microbelist)){
    micname=strsplit(microbelist[i],'; ')
    num=length(micname[[1]])
    if(num==3){
      if(micname[[1]][2]=='p__unclassified Bacteria'){
        A=mic_path[mic_path$Class=='c__Abyssogena phaseoliformis symbiont [TAX:596095]',]

      }else{
        A=mic_path[mic_path$Class==micname[[1]][3],]
      }
      if(length(rownames(A)!=0)){
        A$name=microbelist[i]}
    }
    if(num==4){
      A=mic_path[mic_path$Order==micname[[1]][4],]
      A$name=microbelist[i]
    }
    if(num==5){
      A=mic_path[mic_path$Family==micname[[1]][5],]
      if(length(rownames(A)!=0)){
        A$name=microbelist[i]}
    }
    if(num==6){
      micname6=gsub('g__','',micname[[1]][6])
      A41=mic_path[mic_path$Family==micname[[1]][5],]

      A=A41[grepl(micname6,A41$Species),]
      if(length(rownames(A)!=0)){
        A$name=microbelist[i]}

    }
    if(num==7){
      micname7=gsub('s__','',micname[[1]][7])
      A51=mic_path[mic_path$Order==micname[[1]][4],]
      A52=A51[A51$Family==micname[[1]][5],]
      A53=mic_path[mic_path$Genus==micname[[1]][6],]
      A=mic_path[grepl(micname7,A53$Species),]
      if(length(rownames(A)!=0)){
        A$name=microbelist[i]}
    }

    if(is.null(A)){
      A=data.frame()
    }

    B=rbind(B,A)
  }
  Micropath=unique(B[,c(11,10)])

  # revise ------------------------------------------------------------------

  micro2pathway<-data.frame()
  for(i in Microbe){
    m3<-Micropath[which(Micropath$name==i),]
    m3<-m3[,c(1,2)]
    micro2pathway<-rbind(micro2pathway,m3)
  }
  colnames(micro2pathway)[2]='Pathway'

  micro2func=data.frame()
  for(i in pathway$nodes){
    a1=micro2pathway[which(micro2pathway$Pathway==i),]
    micro2func=rbind(micro2func,a1)
  }

  pathway2=intersect(unique(meta2func$Pathway),unique(micro2func$Pathway))

  meta2func1=data.frame()
  micro2func1=data.frame()
  for(i in pathway2){
    a1=meta2func[which(meta2func$Pathway==i),]
    meta2func1=rbind(meta2func1,a1)
    a2=micro2func[which(micro2func$Pathway==i),]
    micro2func1=rbind(micro2func1,a2)
  }

  # save --------------------------------------------------------------------

  # if(!file.exists("./results/MMfunc")){
  #   dir.create("./results/MMfunc",recursive = T)
  # }

  temp_Node_score <- IFS[["IFS_completed"]][,1:2]
  colnames(meta2func1)[1] <- "nodes"
  meta2func1= merge(meta2func1,temp_Node_score,all.x = T,by="nodes")
  colnames(micro2func1)[1] <- "nodes"
  micro2func1=merge(micro2func1,temp_Node_score,all.x = T,by="nodes")

  # write.csv(meta2func1,file='./results/MMfunc/Metabolite in same pathway.csv',row.names = FALSE)
  # write.csv(micro2func1,file='./results/MMfunc/Microbe in same pathway.csv',row.names = FALSE)


  # Sankey ------------------------------------------------------------------

  rownames(meta2func1)<-1:length(rownames(meta2func1))
  colnames(meta2func1)=c("source","target","ns")

  rownames(micro2func1)<-1:length(rownames(micro2func1))
  colnames(micro2func1)=c("target","source","ns")

  MM2func4<-rbind(meta2func1,micro2func1)
  MM2func4$value<-rep(1,times=length(rownames(MM2func4)))
  nodes<- data.frame(name=c(as.character(MM2func4$source), as.character(MM2func4$target)) %>% unique())
  MM2func4$IDsource=match(MM2func4$source, nodes$name)-1
  MM2func4$IDtarget=match(MM2func4$target, nodes$name)-1

  sn <- sankeyNetwork(Links = MM2func4, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      NodeID = "name",Value='value',
                      units = "TWh", fontSize = 12, nodeWidth = 30)

  # saveNetwork(sn,file = "./results/MMfunc/sankeyNetwork.html")
  meta_in_func <- unique(meta2func1$source)
  mic_in_func <- unique(micro2func1$target)

  MMfunc_res <- list(meta2func = meta_in_func, micro2func = mic_in_func,snfig = sn)
  return(MMfunc_res)
}

