set.overlap = function(tftargs,targs,net){
  tftargs = unique(unlist(tftargs))
  targs=intersect(targs,names(net$bygene))
  
  inter = length(intersect(tftargs,targs))  
  if(inter <3){return(NA)}
  
  fmat =matrix(c(inter,(length(tftargs)-inter),(length(targs)-inter),((length(net$bygene))-length(tftargs ) - length(targs) +inter)),nrow=2)
  return(fisher.test(fmat,alternative="greater")$p.value)   
}
merge.pvalues = function(tftargs,targs,net){
  if(length(intersect(names(targs),unlist(tftargs))) <3){return(NA)}     
  return(fishersMethod(targs[intersect(names(targs),unlist(tftargs))]))  
}

list.enriched = function(tftargs,targs,net){
  inside=which(names(targs) %in% unlist(tftargs))
  outside = setdiff(1:length(targs),inside)
  if(length(inside)<3){return(NA)}
  return(wilcox.test(targs[inside],targs[outside])$p.value)
}


fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower.tail=FALSE)