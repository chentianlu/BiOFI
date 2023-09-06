#' @title IPS
#' @description
#' A function is used to compute the importance pair score (IPS) of the network returned by the TarNet() function.
#'
#' @param TarnetPairs The result returned by TarNet() function.
#' @param IFS The result returned by IFS() function.
#' @return IPS_res A data frame including the importance pair score (IPS) of the target network.
#'
#' @export IPS
#' @examples
#' IFScore <- IFS(microApath = micro.eg, metaApath = metabo.eg,
#'                conf = NULL, groupInfo = groupInfo.eg)
#' MMfunc_res <- MMfunc(IFS = IFScore)
#' TarNet_res <- TarNet(IFS = IFScore, mm2path = MMfunc_res)
#' IPScore <- IPS(TarnetPairs = TarNet_res,IFS = IFScore)

IPS<-function(TarnetPairs,IFS){
  IFSdf1 = IFS[["IFS_completed"]][,1:2]
  colnames(IFSdf1) <- c("to","to_IFS")
  IFSdf2 = IFS[["IFS_completed"]][,1:2]
  colnames(IFSdf2) <- c("from","from_IFS")
  IPS_res <- merge(TarnetPairs[["tarnetpairs"]],IFSdf1,all.x = T)
  IPS_res <- merge(IPS_res,IFSdf2,all.x = T)
  IPS_res$r_score <- abs(IPS_res$r)
  IPS_res$r_score <-  seq(nrow(IPS_res),1,-1)
  IPS_res$r_score <- (IPS_res$r_score - 1)/(nrow(IPS_res) - 1) * 2 + 1
  IPS_res$ips <- IPS_res$to_IFS + IPS_res$from_IFS + IPS_res$r_score
  IPS_res <- IPS_res[order(IPS_res$ips,decreasing = T),]
  return(IPS_res)
}
