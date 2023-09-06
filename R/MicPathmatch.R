#'@title MicPathmatch
#'@description
#' A function is used to trace the affiliation between microbes (16S or Metagenome data) and functional pathways based on Kyoto Encyclopedia of Gene and Genomes (KEGG).
#' You can get two tables using this function.
#' One is a table including predicted detail next level strains and pathways of your own microbes, the other is a table including only your own microbes and corresponding pathways.
#' The format of microbe is like this:
#' k__XXX; p__XXX; c__XXX; o__XXX; f__XXX; g__XXX; s__XXX;
#' k means kingdom;
#' p means phylum;
#' c means class;
#' o means order;
#' f means family;
#' g means genus;
#' s means species.
#'
#' @param microbes A character vector containing the names of microbes.
#' @return micmatch_res A list with two data frames, one contains strains, sub-strains and affiliated pathways and the other data frame does not contain strains or sub-strains.
#' @export MicPathmatch
#' @examples
#' MicPathmatch_res <- MicPathmatch(microbes = MicPathmatch.eg)

MicPathmatch <- function(microbes){

  microbes=microbes
  microApathdb=get0("sysdata", envir = asNamespace("BiOFI"))$micro4path
  if(!length(microbes[grep("; ",microbes)])==length(microbes)){

    stop("The list of microbes may contain inappropriate format,please check the format of microbes ('The correct format is k__XXX; p__XXX; c__XXX; o__XXX; f__XXX; g__XXX; s__XXX;').")
  }
  microbelist=c()
  B=data.frame()
  for(i in 1:length(microbes)){
    microbe=strsplit(microbes[i],'; ')
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
        A=microApathdb[microApathdb$Class=='c__Abyssogena phaseoliformis symbiont [TAX:596095]',]

      }else{
        A=microApathdb[microApathdb$Class==micname[[1]][3],]
      }
      if(length(rownames(A)!=0)){
        A$name=microbelist[i]}
    }
    if(num==4){
      A=microApathdb[microApathdb$Order==micname[[1]][4],]
      if(length(rownames(A)!=0)){
        A$name=microbelist[i]}
    }
    if(num==5){
      A=microApathdb[microApathdb$Family==micname[[1]][5],]
      if(length(rownames(A)!=0)){
        A$name=microbelist[i]}
    }
    if(num==6){
      micname6=gsub('g__','',micname[[1]][6])
      A41=microApathdb[microApathdb$Family==micname[[1]][5],]

      A=A41[grepl(micname6,A41$Species),]
      if(length(rownames(A)!=0)){
        A$name=microbelist[i]}

    }
    if(num==7){
      micname7=gsub('s__','',micname[[1]][7])
      A51=microApathdb[microApathdb$Order==micname[[1]][4],]
      A52=A51[A51$Family==micname[[1]][5],]
      A53=microApathdb[microApathdb$Genus==micname[[1]][6],]
      A=microApathdb[grepl(micname7,A53$Species),]
      if(length(rownames(A)!=0)){
        A$name=microbelist[i]}
    }

    if(is.null(A)){
      A=data.frame()
    }

    B=rbind(B,A)
  }
  microbetotle=B
  B1=unique(B[,c(11,10)])
  microbepredict=B1
  micmatch_res=list('Pathways of microbes(detail)'=microbetotle,'Pathways of microbes'=microbepredict)
  return(micmatch_res)
}
