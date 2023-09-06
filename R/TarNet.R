#' @title TarNet
#' @description
#' A function is used to plot a network that consists of microbes and metabolites in the shared pathway from the result of MMfunc() function.
#'
#' @param IFS The result returned by IFS() function.
#' @param mm2path The results returned by MMfunc() function.
#' @param r A threshold value of correlation coefficient to construct the network. Default value is 0.1.
#' @param p_adjust A threshold value of correlation adjust p value to construct the network. Default  value is 0.05.
#' @return TarnetPairs A data frame including target network pairs.
#'
#' @export TarNet
#' @import igraph
#' @import ggrepel
#' @import visNetwork
#' @importFrom htmlwidgets JS
#' @import ppcor
#' @examples
#' IFScore <- IFS(microApath = micro.eg, metaApath = metabo.eg,
#'                conf = NULL, groupInfo = groupInfo.eg)
#' MMfunc_res <- MMfunc(IFS = IFScore)
#' TarNet_res <- TarNet(IFS = IFScore, mm2path = MMfunc_res)


TarNet<-function(IFS,mm2path,r = 0.1,p_adjust = 0.05){
  mic_meta_cor=IFS[["mic_meta_cor"]]
  mic_meta_cor_p.adjust=IFS[["mic_meta_cor_p.adjust"]]

  meta_in_func=mm2path[[1]]
  mic_in_func=mm2path[[2]]

  mic_meta_cor=mic_meta_cor[mic_in_func,meta_in_func]
  mic_meta_cor_p.adjust=mic_meta_cor_p.adjust[mic_in_func,meta_in_func]

  CorrDF <- function(cormat, pmat) {
    ut = matrix (TRUE, nrow = as.numeric(nrow (cormat)), ncol = as.numeric(ncol (pmat)))
    data.frame(
      from = rownames(cormat)[row(cormat)[ut]],
      to = colnames(cormat)[col(cormat)[ut]],
      r =(cormat)[ut],
      p.adjust = pmat[ut]
    )
  }

  mic_meta_cor_df <- CorrDF(mic_meta_cor,mic_meta_cor_p.adjust)
  # if(!file.exists("./results/TarNet")){
  #   dir.create("./results/TarNet",recursive = T)
  # }

  # write.csv(mic_meta_cor_df,file='./results/TarNet/network_cor_results.csv',row.names = FALSE)

  mic_meta_cor_filtered <- mic_meta_cor_df[which(abs(mic_meta_cor_df$r) > r & mic_meta_cor_df$p.adjust < p_adjust),]

  if(length(mic_meta_cor_filtered$from)==0){
    stop("Too large a value of R results in no edges that satisfy the condition!" )
  }

  centrality = "degree"
  nodeattrib <- data.frame(nodes = union(mic_meta_cor_filtered$from,mic_meta_cor_filtered$to))

  nodeattrib$group <- 0
  for (i in as.character(nodeattrib$nodes)){
    if (i %in% mic_in_func == TRUE){
      nodeattrib[nodeattrib$nodes == i,"group"] <- "Microbe"
    }else if(i %in% meta_in_func == TRUE){
      nodeattrib[nodeattrib$nodes == i,"group"] <- "Metabolite"
    }

  }

  # save file ---------------------------------------------------------------


  # write.csv(nodeattrib,file='./results/TarNet/network_nodes_results.csv',row.names = FALSE)

  co_net <- graph_from_data_frame(mic_meta_cor_filtered,directed = F, vertices = nodeattrib)
  co_net=simplify(co_net, remove.multiple=T)
  E(co_net)[ r< 0 ]$color <- "green"
  E(co_net)[ r>0 ]$color <- "red"

  network=visIgraph(co_net)%>%
    visIgraphLayout(layout = "layout_with_fr", smooth = TRUE)%>%
    visNodes(color = list(hover = "green"))%>%visInteraction(hover = TRUE)%>%visLegend(useGroups = T,stepX = 70,stepY=70)%>%
    visGroups(groupname = 'Microbe',color="purple",shape="triangle")%>%
    visGroups(groupname = 'Metabolite',color="tomato",shape="dot")%>%
    visOptions(selectedBy = 'group')
  # saveNetwork(network,file = "./results/TarNet/network.html")
  TarnetPairs <- list(tarnetpairs = mic_meta_cor_filtered,networkfig = network)
  return(TarnetPairs)
}
