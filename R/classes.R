#defining class
setClass("coregnet",
         representation(GRN="data.frame", 
                        # all Gene Regulatory Networks. Target gene (col1), co-activators (col2 : "a b c ..."), co-repressors (col3 : "z x w ..")
                        adjacencyList="list",
                        # Same network than GRN. 2 entries, bytf and bygene. Each entry is a list, with as many entries as regulators of genes
                        # for each gene or regulators, 2 entries act and rep containing a vector of the activated/repressed targets or regulators
                        # ex1 : adjacencyList$bytf$ELF3$act
                        # ex2 : adjacencyList$bygene$FGFR3$act
                        coRegulators="data.frame",# a description of the frequent coregulators. Unused yet
                        evidences="list", #
                        evidenceDescription="data.frame",
                        inferenceParameters="list"),
         #evidenceDescription
         #"evidenceType","originalGene","originalReg","originalEvidence"
         #"commonGene","commonReg","evidences","commonEvidences", "enrichment","p.value")
         validity=function(object){
           if(nrow(object@GRN)==0){
             return("Need a set of Gene Regulatory Networks (GRN) describing the regulation of each gene")
           }else if( sum(names(object@adjacencyList)  %in% c("bytf","bygene"))<2 | length(object@adjacencyList) != 2){
             return("Need a representation of the network in an adjacency list (automatic if using the coregnet() function)")
           }else{
             return(TRUE)
           }
         }
)


#initialization
coregnet <- function(GRN,expressionDATA=NULL) {
  if(is.data.frame(GRN)){
    if(ncol(GRN) == 2 & !is.null(expressionDATA)){
      stop("not correctly implemented yet")
      
      if((sum(unique(GRN[,1]) %in% rownames(expressionDATA) )>0 ) &
           (sum(unique(GRN[,2]) %in% rownames(expressionDATA) ) > 20  ) ){
        if(length(unique(GRN[,2])) > length(unique(unique(GRN[,1])))){
          warnings("The provided GRN looks like target genes are in the second column when they actually should be in the first")}
        corels=apply(GRN,1,function(gg){return(cor(as.numeric(expressionDATA[gg[1],]),as.numeric(expressionDATA[gg[2],])))})
      }else{
        stop("To few genes of the GRN in the expression dataset.")
      }
    }
    if(ncol(GRN)>=3  ){
      if(!is.character(GRN[,1])| !is.character(GRN[,2])| !is.character(GRN[,3])){
        stop(paste("Does not look like a good input network. All columns should be of type character,"
                   ,"the first describing the target gene and columns 2 and 3 each describing" ,
                   "the set of regulators with a space seperating each regulator in the set."))
      }else{
        colnames(GRN)[1:3] =  c("Target", "coact","corep")
      }
    }
    
  }else if(is.matrix(GRN)){
    if(ncol(GRN) > nrow(GRN) ){
      warning(paste("The provided GRN looks like target genes are columns when they actually should be in line." ,
                    "Running anyway but might get a reverse target -> regulator network."))
    }
    if(is.null(colnames(GRN)) | is.null(rownames(GRN))){
      stop("Need column names and row names to identify regulators and target genes respectively.")
    }
    
    values <- as.vector(GRN)
    simpleTR <- expand.grid(dimnames(GRN),stringsAsFactors =FALSE)
    colnames(simpleTR) <- c("Target","Reg")
    simpleTR=simpleTR[which(values!=0),]
    simpleTR=data.frame(simpleTR, "score" = values[which(values!=0)])
    simpleTR$coact=NA
    simpleTR$corep=NA
    # if their are some scores below 0, we beleive that these mean inhibitor regulation
    if(min(simpleTR$score) < 0){
      message("Negative values are considered as inhibitory regulation.")
      simpleTR$R2 = abs(simpleTR$score)
    }else if(is.null(expressionDATA)){
      stop(paste("It appears that the provided score only has positive scores.",
                 "This is fine as long as a gene expression data is provided in order to diffeentiate positive and negative regulation." ,
                 "The expressionData argument is therefore missing."))
    }else if(sum(c(rownames(GRN),colnames(GRN)) %in% rownames(expressionDATA)) != (ncol(GRN)+nrow(GRN))){
      stop("The genes in the network must all be in the rownames of the input expression data.")
    }else{
      print("computing corrrelation to sign the network")
      simpleTR$R2 = simpleTR$score
      ALLCOREL = cor(t(expressionDATA[unique(c(as.character(unlist(simpleTR[,1])),as.character(unlist(simpleTR[,2])))),]))
      print("corel done")
      simpleTR=do.call(rbind,mclapply(unique(simpleTR$Reg),function(r){
        is=which(simpleTR$Reg ==r )
        tmp=simpleTR[is,]
        cr=ALLCOREL[r,]
        tmp[,"score"]=tmp$score * sign(cr[tmp$Target])
        return(tmp)
      }))
      print("assigning sign done")
    }
    
    simpleTR[which(simpleTR$score>0),"coact"]=simpleTR[which(simpleTR$score>0),"Reg"]
    simpleTR[which(simpleTR$score<0),"corep"]=simpleTR[which(simpleTR$score<0),"Reg"]
    sigrns=simpleTR[,c("Target","coact","corep","score","R2")]
    
    act = sigrns$coact
    rep = sigrns$corep
    is = 1:length(rep)
    
    sigrns=unique(sigrns)
    SIGRNS = list()
    
    bygene=apply(GRN,1,function(l,tfs){
      act=tfs[which(l>0)]
      act =unname( act[which(!is.na(act) & act !="")])
      rep=tfs[which(l<0)]
      rep =unname( rep[which(!is.na(rep) & rep !="")])
      return(list("act"=act,"rep"=rep))
    },tfs=colnames(GRN))
    
    bytf=apply(GRN,2,function(l,genes){
      act=genes[which(l>0)]
      act =unname( act[which(!is.na(act) & act !="")])
      rep=genes[which(l<0)]
      rep =unname( rep[which(!is.na(rep) & rep !="")])
      return(list("act"=act,"rep"=rep))
    },genes=rownames(GRN))
    
    coregnetwork= new("coregnet",GRN=sigrns  ,"adjacencyList"=list("bygene"=bygene,"bytf"=bytf))
    coregnetwork@coRegulators = coregulators(coregnetwork,verbose=FALSE,alpha=1)
    return(  coregnetwork )
  }else{
    stop("Wrong input format. See help.")
  }
  reshapedNet = .quicknonuniqgrnsTOSIGRNS(GRN)
  coregnetwork= new("coregnet",GRN=reshapedNet$sigrns,"adjacencyList"=reshapedNet$adjList)
  coregnetwork@coRegulators = coregulators(coregnetwork,verbose=FALSE,alpha=1)
  return(  coregnetwork )
}








setGeneric("coregnetToList", function(network) {
    standardGeneric("coregnetToList")
})

setMethod("coregnetToList", signature(network = "coregnet"),
          function(network){  
              return(lapply(lapply(lapply(network@adjacencyList$bytf,unlist),unname),unique))
              })
   

setGeneric("coregnetToDataframe", function(network) {
    standardGeneric("coregnetToDataframe")
})

setMethod("coregnetToDataframe", signature(network = "coregnet"),
          function(network){  
              adj=coregnetToList(network)
              return(data.frame("Regulator"=unlist(lapply(names(adj),function(y){rep.int(y,length(adj[[y]]))})),
                                "Target"=unlist(adj,use.names=FALSE),stringsAsFactors=FALSE))
          })
    




setGeneric("targets", function(object,regulator=NULL,type=c("regulating","activating","repressing")) {
  standardGeneric("targets")
})
setGeneric("activators", function(object,target,type=c("single","coregulators")) {
  standardGeneric("activators")
})
setGeneric("repressors", function(object,target,type=c("single","coregulators")) {
  standardGeneric("repressors")
})
setGeneric("regulators", function(object,target=NULL,type=c("single","coregulators")) {
  standardGeneric("regulators")
})

setMethod("targets",    signature(object = "coregnet"), function(object,regulator=NULL,type=c("regulating","activating","repressing")){
  type <- match.arg(type)
  if(length(regulator)>=1){
    if(length(intersect( regulator , names(regulators(object))))==0){
      return(NA)
    }
    if(type=="regulating"){
      return(c(na.omit(unique(unlist(object@adjacencyList$bytf[regulator])))))
    }else{
      type=c("repressing"="rep","activating"="act")[type]
      return(c(na.omit(unique(unlist(lapply(object@adjacencyList$bytf[regulator],function(x){return(x[[type]])}))))))
    }
    
  }else{
    return(c(na.omit(unique(unlist(object@adjacencyList$bytf)))))
  }
})
setMethod("activators", signature(object = "coregnet"), function(object,target,type=c("single","coregulators")){
  type <- match.arg(type)
  if(length(target )==1){
    if(! target %in% targets(object)){
      return(NA)
    }
    if(type=="single"){
      return(c(na.omit(unique(object@adjacencyList$bygene[[target]]$act))))
    }else{
      coreg=unique(object@GRN[which(object@GRN$Target == target ),"coact"])
      return(c(na.omit(coreg)))
    }
  }else if(length(target)>1){
    if(length(intersect( target , targets(object)))==0){
      return(NA)
    }
    if(type=="single"){
      return(c(na.omit(unique(unlist(lapply(object@adjacencyList$bygene[target],function(x){return(x$act)}))))))
    }else{
      coreg=unique(object@GRN[which(object@GRN$Target %in% target ),"coact"])
      return(c(na.omit(coreg)))
    }
    
  }else{
    return(c(na.omit(unique(unlist(lapply(object@adjacencyList$bygene,function(x){return(x$act)}))))))
  }
})
setMethod("repressors", signature(object = "coregnet"), function(object,target,type=c("single","coregulators")){
  type <- match.arg(type)
  if(length(target )==1){
    if(! target %in% targets(object)){
      return(NA)
    }
    if(type=="single"){
      return(c(na.omit(unique(object@adjacencyList$bygene[[target]]$rep))))
    }else{
      coreg=unique(object@GRN[which(object@GRN$Target == target ),"corep"])
      return(c(na.omit(coreg)))
    }
  }else if(length(target)>1){
    if(length(intersect( target , targets(object)))==0){
      return(NA)
    }
    if(type=="single"){
      return(c(na.omit(unique(unlist(lapply(object@adjacencyList$bygene[target],function(x){return(x$rep)}))))))
    }else{
      coreg=unique(object@GRN[which(object@GRN$Target %in% target ),"corep"])
      return(c(na.omit(coreg)))
    }
    
  }else{
    return(c(na.omit(unique(unlist(lapply(object@adjacencyList$bygene,function(x){return(x$act)}))))))
  }
})
setMethod("regulators", signature(object = "coregnet"), function(object,target=NULL,type=c("single","coregulators")){
  type <- match.arg(type)
  if(length(target )==1){
    if(! target %in% targets(object)){
      return(NA)
    }
    if(type=="single"){
      return(c(na.omit(unique(unlist(object@adjacencyList$bygene[[target]])))))
    }else{
      coreg=unique(object@GRN[which(object@GRN$Target == target ),c("coact","corep")])
      return(c(na.omit(coreg)))
    }
  }else if(length(target)>1){
    if(length(intersect( target , targets(object)))==0){
      return(NA)
    }
    if(type=="single"){
      return(c(na.omit(unique(unlist(object@adjacencyList$bygene[target])))))
    }else{
      coreg=unique(object@GRN[which(object@GRN$Target %in% target ),c("coact","corep")])
      return(c(na.omit(coreg)))
    }
    
  }else{
    
    regs=sort(sapply(lapply(object@adjacencyList$bytf,unlist),length),decreasing=TRUE)
    regs=regs[which(names(regs) !="NA" & names(regs) != "EMPTY" & !is.na(names(regs)))]
    return(regs)
  }
})



setGeneric("coregulators", function(object,maxcoreg=2,verbose=TRUE,minCommonGenes=ifelse(maxcoreg==2,1,10),adjustMethod="fdr",alpha=0.01) {
  standardGeneric("coregulators")
})



setMethod("coregulators", signature(object = "coregnet"), function(object,maxcoreg=2,verbose=FALSE,
                      minCommonGenes=ifelse(maxcoreg==2,1,10),adjustMethod="fdr",alpha=0.01) {
  grn=object@GRN
  if(length(grep(" ",grn$coact)) ==0 &length(grep(" ",grn$corep))==0){
    warning("No natural co-regulators found in the network. Only pairs will be returned.")
    maxcoreg=2
  }
  
  if(maxcoreg==2 ){
      if(minCommonGenes==1 & length(object@coRegulators)>0 ){
        if("fisherTest" %in% colnames(object@coRegulators)){
          object@coRegulators$adjustedPvalue = p.adjust(object@coRegulators$fisherTest,method=adjustMethod)
          return(object@coRegulators[which(object@coRegulators$adjustedPvalue <= alpha),])
        }
      }
      if(verbose){  message("Transforming co-regulation network.")}   
      # if there are some GRN with sets of regulators (several TF)
      if((sum(grepl(" ",grn$coact))+sum(grepl(" ",grn$corep)))>0){                      
        if(verbose){message("Mining co-occurences.")}
        trans=as(c(strsplit(grn$coact," "),strsplit(grn$corep," ")),"transactions")
        frcoreg=suppressWarnings(apriori(trans,parameter=list(support =(minCommonGenes)/nrow(grn),
                          maxlen=2,minlen=2,target="frequent itemsets"),control=list(verbose=FALSE)))
        coregs = data.frame(do.call(rbind,as((slot(frcoreg,"items")),"list")),"support" = as((slot(frcoreg,"quality")),"vector"))
        colnames(coregs) = c("Reg1","Reg2","Support")
        coregs[,1]=as.character( coregs[,1])
        coregs[,2]=as.character( coregs[,2])
        coregs = coregs[order(coregs[,3],decreasing=TRUE),]
        coregs$nGRN = coregs[,3] * length(trans)      
    }else{
      adjlist = object@adjacencyList
      if(verbose){  message("Transforming co-regulation network.")}
      trans=as(c(lapply(adjlist$bygene,function(x){return(x$act)}),lapply(adjlist$bygene,function(x){return(x$rep)})),"transactions")
      if(verbose){message("Mining co-occurences.")}
      frcoreg=suppressWarnings(apriori(trans,parameter=list(support =minCommonGenes/(length(adjlist$bygene)*2),
                        maxlen=2,minlen=2,target="frequent itemsets"),control=list(verbose=FALSE)))     
      if(length(frcoreg)==0){
        warning("No coregulators found with this threshold.")
        return( data.frame())
      }
      coregs = data.frame(do.call(rbind,as((slot(frcoreg,"items")),"list")),"support" = as((slot(frcoreg,"quality")),"vector"))
      colnames(coregs) = c("Reg1","Reg2","Support")
      coregs[,1]=as.character( coregs[,1])
      coregs[,2]=as.character( coregs[,2])
      coregs = coregs[order(coregs[,3],decreasing=TRUE),]      
      coregs$nGRN = coregs[,3] * length(trans)  
    }
    universe = unique(names(object@adjacencyList$bygene))    
    coregs$fisherTest=unlist(mclapply(data.frame(t(coregs[,1:2]),stringsAsFactors=FALSE),function(co){
      r1=object@adjacencyList$bytf[[co[1]]]
      r1 = list("act" = setdiff(r1$act,r1$rep),"rep"=setdiff(r1$rep,r1$act))
      r2=object@adjacencyList$bytf[[co[2]]]
      r2 = list("act" = setdiff(r2$act,r2$rep),"rep"=setdiff(r2$rep,r2$act))
      coregulatedgenes = unique(c(intersect(r1$act,r2$act),intersect(r1$rep,r2$rep)))
      anticoregulatedgenes = unique(c(intersect(r1$rep,r2$act),intersect(r1$rep,r2$act)))
      r1Only =unique(c( setdiff(r1$act,r2$act)  ,setdiff(r1$rep,r2$rep)   ))
      r2Only =unique(c( setdiff(r2$act,r1$act)  ,setdiff(r2$rep,r1$rep)   ))
      reste  = setdiff(universe,unique(unlist(c(r1,r2))))
      return( fisher.test(matrix(c(length(coregulatedgenes),length(r1Only),length(r2Only),length(reste)),nrow=2),alternative="greater")$p.value)
    }))    
    rownames(coregs)=NULL
    coregs$adjustedPvalue = p.adjust(coregs$fisherTest,method=adjustMethod)
    return(coregs[which(coregs$adjustedPvalue <= alpha),])
    
  }else if(maxcoreg > 2){
    adjlist = object@adjacencyList
    if(verbose){  message("Transforming co-regulation network.")}
    listGRNact=as(lapply(adjlist$bygene,function(x){return(x$act)}),"transactions")
    listGRNrep=as(lapply(adjlist$bygene,function(x){return(x$rep)}),"transactions")
    if(verbose){message("Mining co-occurences.")}
    frcoact=suppressWarnings(apriori(listGRNact,parameter=list(support =minCommonGenes/length(listGRNact),
                                                            maxlen=maxcoreg,minlen=2,target="maximally frequent itemsets"),control=list(verbose=FALSE)))
    frcorep=suppressWarnings(apriori(listGRNrep, parameter=list(support =minCommonGenes/length(listGRNrep),
                                                            maxlen=maxcoreg,minlen=2,target="maximally frequent itemsets"),control=list(verbose=FALSE)))
    frcoreg= union(frcorep,frcoact)
    print(length(frcoreg))
    if(length(frcoreg)==0){
      warning("No coregulators found with this threshold.")
      return(data.frame())
    }
    coregs = data.frame("CoRegulators"=sapply(as((slot(frcoreg,"items")),"list"),paste,collapse=" "),
                        "Support" = as((slot(frcoreg,"quality")),"vector"),stringsAsFactors = FALSE)
    coregs = coregs[order(coregs[,2],decreasing=TRUE),]
    colnames(coregs) = c("CoRegulators","Support")
    rownames(coregs)=NULL
  }else{stop("Wrong number of maximum coregulators. Must be 2 or more.")}
  return(coregs)
})


#new printing method
setMethod("print","coregnet",
          function(x){
            print(paste( length(x@adjacencyList$bytf), "Transcription Factors. ",length(x@adjacencyList$bygene),
                         "Target Genes. ",length(unlist(x@adjacencyList$bygene)),"Regulatory interactions."))
            if(nrow(x@evidenceDescription) >=1){
              print(x@evidenceDescription[,5:10])
            }else{
              print("No added evidences.")
            }
            
          }
)

#new method to printout
setMethod("show",signature(object="coregnet"),
          function(object){
            print(paste( length(object@adjacencyList$bytf), "Transcription Factors. ",length(object@adjacencyList$bygene),
                         "Target Genes. ",length(unlist(object@adjacencyList$bygene)),"Regulatory interactions."))
            if(nrow(object@evidenceDescription) >=1){
              print(object@evidenceDescription[,5:10])
            }else{
              print("No added evidences.")
            }
            
          }
)

# new length method
setMethod("length","coregnet",
          function(x){
            return(c("TF"=length(x@adjacencyList$bytf),"Gene"=length(x@adjacencyList$bygene),"Interactions"=length(unlist(x@adjacencyList$bygene))))
          }
)

# new dimension method
setMethod("dim","coregnet",
          function(x){
            return(c("TF"=length(x@adjacencyList$bytf),"Gene"=length(x@adjacencyList$bygene),"Interactions"=length(unlist(x@adjacencyList$bygene))))
          }
)


# new summary method
setMethod("summary","coregnet",
          function(object,...){
            print(object,...)
          }
)




.val2col<-function(z, zlim, col = NULL, breaks){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  CUT <- cut(z, breaks=breaks)
  colorlevels <- col[match(CUT, levels(CUT))] # assign colors to heights for each point
  return(colorlevels)
}


.traits <- function(m) {
  nc <- ncol(m)
  nr <- nrow(m)
  
  qui <- which(!is.na(m[1,]))
  d <- 0.5
  x <- (0:nr-d)/(nr-1)
  
  for (y in qui)
    lapply(x, function(u)
      lines(c(u,u), c((y-1-d)/(nc-1), (y-d)/(nc-1))))
  abline(h=(0:nc-d)/(nc-1))
}
