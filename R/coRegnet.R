


masterRegulator = function(coregnet,targetGenes,method=c("set.overlap","merge.pvalues","list.enriched"))
{
  method=match.arg(method)  
  ## testing the input
  
  if(class(coregnet) !="coregnet"){
    stop("At this function can only be used with a network of type CoRegNet.")
  }
  coRegNetwork = coregnet@adjacencyList
  if(method == "set.overlap" & !is.character(targetGenes)){
    stop(paste("When inferring Master Regulators using a set of genes to overlap the targetGenes", 
               "input needs to be a vector of charcter containing the set of target genes")  )
  }else if(method != "set.overlap" & (is.data.frame(targetGenes) | is.matrix(targetGenes))){
    x=targetGenes[,1]
    names(x)=rownames(targetGenes)
    targetGenes =x
  }else if((method != "set.overlap" & (is.numeric(targetGenes)|is.integer(targetGenes)) &
            length(intersect(names(targetGenes),names(coRegNetwork$bygene)))==0)){    
    stop(paste("The input of the target genes of interest is not usable by this." ,
               "Please see the manual or contact the maintainer to add a new possibility to the package."))
  }
  
  MR=switch(method,
            set.overlap = sapply(coRegNetwork$bytf,  set.overlap  ,net=coRegNetwork,targs=targetGenes),
            merge.pvalues = sapply(coRegNetwork$bytf,  merge.pvalues  ,net=coRegNetwork,targs=targetGenes),
            list.enriched = sapply(coRegNetwork$bytf,  list.enriched  ,net=coRegNetwork,targs=targetGenes))
  
  return(sort(MR))
}



setGeneric("fitCoregnet", function(network,expData,permutation=0) {
  standardGeneric("fitCoregnet")
})

setMethod("fitCoregnet", signature(network = "coregnet"), function(network,expData,permutation=0){  
  
  # making sure that genes in the network and genes in expression dataset are the same
  allgenesandregs = unique(c(unlist(network@adjacencyList$bytf),names(network@adjacencyList$bygene)))    
  if(length(intersect(allgenesandregs,rownames(expData))) < 0.1*length(allgenesandregs)){
    stop(paste("More than 90% of the network genes are not in the expression data." ,
               "The influence of the regulators cannot be computed with these settings.\nNote" ,
               ": The expression data must be given with genes in line."))}
  netregs=unique(names(network@adjacencyList$bytf))
  netgenes=unique(names(network@adjacencyList$bygene))
  expgenes<-rownames(expData)
  samples<-colnames(expData)
  expData = t(expData)
  expgenes -> colnames(expData)
  samples->rownames(expData)
  # selecting only regulators in the dataset, only genes in the data set and only genes with regulators in the dataset
  netregs=intersect(netregs,expgenes)
  netgenes=intersect(netgenes,expgenes)
  isuniqgenes=(nrow(network@GRN) == length(unique(network@GRN[,1])))
   netgenes=netgenes[which(unlist(mclapply(network@adjacencyList$bygene[netgenes],function(regulators){
    regulators=unique(unlist(regulators))
    return( length(regulators)== length(intersect(regulators,expgenes)) )
  })    ))]
  grnToTest = network@GRN[which(network@GRN[,1] %in% netgenes),]
  genexp=expData[,netgenes]
  regexp = expData[,netregs]
  
  listedgrn = data.frame(t(grnToTest[,1:3]),stringsAsFactors=FALSE)
  result=mclapply(listedgrn,.fitGRN, genexp=genexp,regexp=regexp)
  
  fittedExp =t(.getEntry(result,"fitted"))  
  errors =t(.getEntry(result,"residuals"))  
  measures = .getEntry(result,"numscores")
  R2 = as.numeric(measures["R2",]  )
  RMSE = as.numeric(measures["RMSE",])
  
  
  if(isuniqgenes){
    rownames(errors) = netgenes
    rownames(fittedExp) = netgenes
    names(R2)=netgenes
    names(RMSE)=netgenes
  }
  if(permutation >0){
    permutFit= mclapply(listedgrn,function(grn){
      
      results=sapply(1:permutation  ,function(i){   return( 
        .fitGRN(grn, genexp=genexp,regexp=regexp,permut=TRUE) 
      )})
      return(c(mean(R2[grn[1]] <= results[1,]),mean(RMSE[grn[1]] >= results[2,])))
    })
    permutFit=do.call(rbind,permutFit)
    permutR2 = permutFit[,1]    
    permutRMSE = permutFit[,2]    
    if(isuniqgenes){
      names(permutR2)=netgenes
      names(permutRMSE)=netgenes
    }
    return(list("fitted.values"=fittedExp,"fitted.residuals"=errors,"R2"=R2,"RMSE"=RMSE,"quantile.R2"=permutR2,"quantile.RMSE"=permutRMSE))
  }else{
    return(list("fitted.values"=fittedExp,"fitted.residuals"=errors,"R2"=R2,"RMSE"=RMSE))
  }
  
  
})









setGeneric("regulatorInfluence", function(object,expData,minTarg = 10,withEvidences=FALSE,addCoregulators=FALSE,is.scaled=FALSE) {
  standardGeneric("regulatorInfluence")
})

setMethod("regulatorInfluence", signature(object = "coregnet"), function(object,expData,minTarg = 10,withEvidences=FALSE,addCoregulators=FALSE,is.scaled=FALSE)
{
  adjlist = object@adjacencyList
  allgenes = unique(unlist(adjlist$bytf))  
  
  if(length(intersect(allgenes,rownames(expData))) < 0.1*length(allgenes)){
    stop(paste("More than 90% of the network genes are not in the expression data." ,
               "The influence of the regulators cannot be computed with these settings.\nNote" ,
               ": The expression data must be given with genes in line."))
  }
  
  if( sum(!allgenes %in% rownames(expData))>0){
    adjlist=.subsigrns(adjlist,rownames(expData))
  }    
  sampgenes = dimnames(expData)
  if(! is.scaled){
    expData = t(scale(t(expData),scale=FALSE))
    dimnames(expData)=sampgenes
  }
  
  tfs = names(adjlist$bytf)  
  
  if(addCoregulators & !withEvidences){
    object@adjacencyList=adjlist
    freqcotfs=coregulators(object,maxcoreg=length(object@adjacencyList$bytf),minCommonGenes=minTarg,verbose=FALSE)  
    cotfs = freqcotfs[,1]
    tfs=c(tfs ,cotfs)
  }
  if( withEvidences){
    if(is.null(object@evidenceDescription) ){stop("No evidences added.")
    }else if(sum(object@evidenceDescription$evidenceType == "regulatory")==0){
      stop("No regulatory evidence added.")
    }
    RegEv=rownames(object@evidenceDescription)[which(object@evidenceDescription[,1]=="regulatory")]  
  }
  
  
  tfScore = mclapply(tfs,function(tf){
    
    if(length(unlist(strsplit(tf," "))) > 1){
      cotf=unlist(strsplit(tf," "))
      acti = names(which(table(unlist(lapply(adjlist$bytf[cotf],function(x){return(x$act)})))== length(cotf)))
      repr =  names(which(table(unlist(lapply(adjlist$bytf[cotf],function(x){return(x$rep)})))== length(cotf)))
    }else{
      
      acti = adjlist$bytf[[tf]]$act
      repr =  adjlist$bytf[[tf]]$rep
      
      if(withEvidences& length(repr )>minTarg &(length(acti) > minTarg)){
        validTargets =unique(unlist(lapply(RegEv,function(regevname){
          ev=object@evidences[[regevname]]
          return(ev[which(ev[,1] == tf),2])
        })))
        repr=intersect(repr,validTargets)
        acti=intersect(acti,validTargets)
      }      
    }
    if(length(repr )>minTarg &(length(acti) > minTarg )) {                        
      return(
        unlist( lapply(1:ncol(expData),function(tumor){                    
          g = (t.test( (expData[acti,tumor]) , (expData[repr,tumor]),alternative="two.sided"))
          return(g$statistic )
        }))
      )
    }else{ return(NULL)  }        
  })        
  names(tfScore) = unlist( tfs)
  tfscore = do.call(rbind,tfScore)
  colnames(tfscore) = colnames(expData)
  return(tfscore)  
})







setGeneric("addCooperativeEvidences", function(object,...) {
  standardGeneric("addCooperativeEvidences")
})
setMethod("addCooperativeEvidences", signature(object = "coregnet"), function(object,...){
  
  #get names in ... and each of the data.frames or lists in ...
  evnames <-as.character(unlist( as.list(substitute(list(...)))[-1L]))
  allevidence <- list(...)
  badNames=grep("\\(|\\)|\\[|\\]",evnames)
  if(length(badNames)>0){
    evnames[badNames] = paste("cooperative", (1:length(badNames))+nrow(object@evidenceDescription) ,sep="")
  }
  
  if(ncol(object@coRegulators )== 1){
    object@coRegulators
  }
  
  oneworked=FALSE  
  for( i in 1:length(allevidence)){
    print(evnames[i])
    netaddedev =.addOneCoRegulatoryEvidence(object,allevidence[[i]],evnames[i])
    if(!is.null(netaddedev)){
      print(paste(evnames[i],"was integrated into the network."))
      oneworked=TRUE
      object =netaddedev
    }
  }
  if(!oneworked){
    message("None of the additional evidences was integrated.")
  }
  return(object)
})




setGeneric("addEvidences", function(object,...) {
  standardGeneric("addEvidences")
})

setMethod("addEvidences", signature(object = "coregnet"), function(object,...){  
  #get names in ... and each of the data.frames or lists in ...
  evnames <-as.character(unlist( as.list(substitute(list(...)))[-1L]))
  allevidence <- list(...)
  badNames=grep("\\(|\\)|\\[|\\]",evnames)
  if(length(badNames)>0){
    evnames[badNames] = paste("regulation", (1:length(badNames))+nrow(object@evidenceDescription) ,sep="")
  }  
  oneworked=FALSE  
  for( i in 1:length(allevidence)){
    netaddedev =.addOneRegulatoryEvidence(object,allevidence[[i]],evnames[i])
    if(!is.null(netaddedev)){
      message(paste(evnames[i],"was integrated into the network."))
      oneworked=TRUE
      object =netaddedev
    }        
  }
  if(!oneworked){
    message("None of the additional evidences was integrated.")
  }
  return(object)
})



setGeneric("refine", function(object,GRNselection=c("best","maximize","threshold"),integration=c("unsupervised","supervised"),
                              referenceEvidence=NULL,evidenceToMaximize="R2",threshold=NULL,verbose=TRUE) {
  standardGeneric("refine")
})

setMethod("refine", signature(object = "coregnet"), function(object,GRNselection=c("best","maximize","threshold"),
                                                             integration=c("unsupervised","supervised"), referenceEvidence=NULL,
                                                             evidenceToMaximize="R2",threshold=NULL,verbose=TRUE){
  
  # Step 1
  # input verif (kinda) and preprocessing.
  
  GRNselection <- match.arg(GRNselection)
  integration<- match.arg(integration)
  

  grns=object@GRN
  acts = strsplit(grns[,2]," ")
  reps =strsplit(grns[,3]," ")
  nRepressors=sapply((reps),function(r){sum(!is.na(r))})
  nActivators =sapply((acts),function(r){sum(!is.na(r))})
  evi=c(names(object@evidences),"R2")
  famille="gaussian"
  
  
  if(!is.null(referenceEvidence)){
    if(!referenceEvidence %in% evi){
      stop(paste("The reference column does not exist in the coregulatory network. Here are the type of evidences that can be used as reference evidence for supervied refinement:",
                 paste(evi,collapse=" ")))    
    }}  
  
  # Step 2
  # Scoring each of the GRN.
  # Here, what you get is A score per GRN (usually several potential GRN per gene)
  # In unsupervised mode it's the mean of all evidences (R2 and other available evidences)
  # In supervised mode, you need a reference evidence that will be used to learn weights to maximize these ref evidences (using lm on the proportion of ref evidence in a GRN)
  
  if(integration=="supervised" & !is.null(referenceEvidence)) {
    form=paste(referenceEvidence,"~",paste(setdiff(evi,referenceEvidence),collapse="+"))
    print(form)
    fit=glm(form,family=famille,data=grns)
    grns$MergeScore=predict(fit,type="response")
    print(summary(fit)$coefficients)  
  }else{
    if(length(evi)==1){
      grns$MergeScore=grns[,evi]
    }else{
      grns$MergeScore=apply(grns[,evi],1,mean)    
    }
  }
  
  
  # Now we have a scorefor each GRN, either a mean or a weighted sum of all evidences (including R2)
  
  # Step 3 
  # score based GRN selection
  # Three scenarios :
  # -  Maximization. This is kinda of tricky and may not be very stable but can give Excelent results.
  #                   Basically it order the GRNs per score and computes the ratio of ReferenceInteraction to predicted interactions.
  #                   Then, the score at which this  ratio is maximal but not in the extreme number of selected GRN is selected as the threshold.
  # - best. Just get one GRN per gene, the one with the best score.
  # - threshold. Select all the GRN with a score above the threshold.
  
  sigrns=NULL
  #################
  if(GRNselection == "maximize"  ){
    nreg = nRepressors+nActivators
    sc = grns$MergeScore
    if(is.null(evidenceToMaximize) ){
      evidenceToMaximize=referenceEvidence
    }
    ref=grns[,evidenceToMaximize]
    ordre = order(sc)
    nintpred = cumsum(nreg[ordre])
    nintval = cumsum(ref[ordre])
    
    maxratio = max( (nintval/nintpred)[as.integer(0.5*length(unique(grns$Target))):as.integer(0.9*nrow(grns))]  )
    index = which((nintval/nintpred)==maxratio)
    if(length(index)> 1){
      index = index[which(as.integer(0.8*length(unique(grns$Target))) &  index <as.integer(0.9*nrow(grns)) )]
      if(length(index)> 1){
        index = index[which.min(index)]
      }
    }
    
    
    thresh=sc[ordre[index]]
    
    if(verbose){
      plot(sc[ordre],(nintval/nintpred),type="l",xlab="Merged Score",ylab="Ratio of validated interactions over predicted interaction",
           main="GRN score and ratio of valid predicted interaction",
           sub="The red line maximises the ratio and is the choosen threshold")
               
      abline(v=thresh,col="red")
    }
    sigrns = grns[which(grns$MergeScore >= thresh),]
    
    #################
    
  }else if(GRNselection == "threshold" & !is.null(threshold)){
    sigrns=grns[which(grns$MergeScore >= threshold),]
    
    #################
  }else{
    sigrns=do.call(rbind,lapply(unique(grns$Target),function(ta){
      tmp=grns[which(grns$Target == ta),]
      return(tmp[which(tmp$MergeScore == max(tmp$MergeScore) ),])}) )
    
  }
  reshapedNet = .quicknonuniqgrnsTOSIGRNS(sigrns)
  object@GRN=reshapedNet$sigrns
  object@adjacencyList =reshapedNet$adjList
  object@coRegulators =  data.frame()
  object@coRegulators = coregulators(object,verbose=FALSE,alpha=1)
  for( evname in names(object@evidences)){
    object@evidenceDescription[evname,5:10]= .descriptionUpdate(object,evname,
                                                                object@evidenceDescription[evname,"evidenceType"])
  }
  return(object)
})




discretizeExpressionData = function(numericalExpression,threshold=NULL,refSamples=NULL,standardDeviationThreshold=1){
  
  numericalExpression=as.matrix(numericalExpression)
  
  
  if(!is.null(refSamples) ){
    refmeans = apply(numericalExpression[,refSamples],1,mean)
    centered  =t(scale(t(numericalExpression[,setdiff(colnames(numericalExpression),refSamples)]),scale=FALSE,center=refmeans))  
    rownames(centered) = rownames(numericalExpression)
    colnames(centered) = setdiff(colnames(numericalExpression),refSamples)
  }else if(min(matrix(numericalExpression) ) >= 0){#  means that it's raw (log or not) data
    centered  =t(scale(t(numericalExpression),scale=FALSE))  
    dimnames(centered) = dimnames(numericalExpression) 
  }else{
    centered=numericalExpression
    centered[which(is.nan(centered))]=0
  }
  if(is.null(threshold)){
    threshold=sd(centered)*standardDeviationThreshold
  }
  nco = ncol(centered)
  nro = nrow(centered)
  discreteExpression =matrix(  as.integer( centered >= threshold  ) + (-as.integer(centered <= (- threshold) )),nrow=nro,ncol=nco)
  dimnames(discreteExpression) = dimnames(centered)
  return(discreteExpression)
}








hLICORN=function( numericalExpression,discreteExpression=discretizeExpressionData(numericalExpression)
, TFlist, GeneList=setdiff(rownames(numericalExpression),TFlist),parallel = c("multicore","no", "snow"),cluster=NULL,
minGeneSupport=0.1,minCoregSupport = 0.2,maxCoreg=floor(log(ncol(numericalExpression),base=3)),
searchThresh=0.5,nGRN=NA,verbose=FALSE){
    
    
    
    #######  #######  #######  #######  #######  #######
    # INPUT VERIFICATION
    if(  sum(! unique(discreteExpression) %in% -1:1) > 0  ){
        stop("Discrete expression data should only have values in {-1, 0, 1}")}
    
    if(length(rownames(numericalExpression)) > length(unique(rownames(numericalExpression)))){
        stop("No gene duplicates are allowed in the row.names.")
    }
    
    if(nrow(numericalExpression) != nrow(discreteExpression) |
    sum(rownames(discreteExpression) != rownames(numericalExpression))>0 ){
        stop("Discrete expression and continuous expression should have the same dimensions and the same rownames (gene/tf names)") }
    
    if(length(intersect(TFlist,rownames(numericalExpression)))<=1 ){
        stop("At least 2 of the provided regulators/transcription factor (TFlist) should be in the rownames in the gene expression matrix")    }
    if(ncol(numericalExpression) > nrow(numericalExpression)){
        warning("Expression data should be in a matrix or data frame with genes in rows and samples in column.")
    }
    if(length(intersect(GeneList,rownames(numericalExpression)))==0 ){
        stop("The list of genes (GeneList) should be in the rownames in the gene expression matrix")    }
    
    
    if(verbose){
        message("Pre-process.")
    }
    
    
    # select only genes and TF with ones or minus ones in at least minGeneSupport portion of samples
    genesupport = which(apply(abs(discreteExpression), 1 , sum) > (ncol(numericalExpression)*(minGeneSupport)))
    discreteExpression=discreteExpression[genesupport,]
    numericalExpression=numericalExpression[genesupport,]
    TFlist = intersect(rownames(numericalExpression),TFlist)
    GeneList= intersect(rownames(numericalExpression),GeneList)
    
    if(length(TFlist)<5){
        stop("Less than 5 of the provided TF are suitable to infer a network. Either provide more TF, more variations in the discrete dataset (more 1 or -1) or decrease the minGeneSupport parameter to select more but less variant TFs.")
    }
    if(length(GeneList)==0){
        stop("No genes were suitable to infer regulators. Either provide more variations in the discrete dataset (more 1 or -1) or decrease the minGeneSupport parameter to allow the selection of more but less variant Genes.")
    }
    
    
    # Get all the matrices and datasets needed (gene and tf expression, numerical or discrete)
    #If only one gene is given, R will automatically make a vector. The following make sure this does not happen.
    if(length(GeneList)==1){
        geneNumExp= matrix(numericalExpression[GeneList,],nrow=1)
        geneDiscExp= matrix(discreteExpression[GeneList,],nrow=1)
        rownames(geneNumExp)=GeneList
        rownames(geneDiscExp)=GeneList
    }else{
        geneNumExp= numericalExpression[GeneList,]
        geneDiscExp= discreteExpression[GeneList,]
    }
    regNumExp= numericalExpression[TFlist,]
    regDiscExp= discreteExpression[TFlist,]
    
    ##    ##    ##    ##    ##    ##    ##    ##    ##
    ## TRANSFORMING ALL DISCRETE DATA INTO TRANSACTIONs
    # To run apriori, the discrete data must be binary. So, the discrete data is simply becoming two concatenated binary matrix
    # first n samples are positive expression values, then all negative values.
    posSamples = 1:ncol(discreteExpression)
    negSamples= (ncol(discreteExpression) +1):(ncol(discreteExpression) *2)
    regBitData =cbind(regDiscExp==+1 , regDiscExp== -1)
    transRegBitData= as(t(regBitData),"transactions")
    
    if(verbose){
        message("Mining coregulator ...")
    }
    
    ##    ##    ##    ##    ##    ##    ##    ##    ##
    ## MINING FOR FREQUENT COREGULATORS
    # using apriori instead of eclat. testing may be required for possible speed improvement.
    miningFunction=apriori
    transitemfreq =suppressWarnings(miningFunction(transRegBitData,parameter=list(support = minGeneSupport/2,maxlen=1,target="frequent itemsets")
    ,control=list(verbose=FALSE)))
    if(maxCoreg > 1){
        transitemfreq=c(transitemfreq,suppressWarnings(miningFunction(transRegBitData,parameter=list(support =minCoregSupport/2,minlen=2,maxlen=maxCoreg,target="closed frequent itemsets")
        ,control=list(verbose=FALSE))))
    }
    coregs =as(slot(transitemfreq,"items"),"list")
    
    
    if(verbose){
        message(paste("Learning a Co-Regulatory network for:\n",  length(GeneList)," target genes, ",length(TFlist)," regulators and a ",
        "total of coregulator sets ",length(coregs),"sets of potential co-regulators.\nSearch parameters :\n",
        "Maximum size of co-regulator sets : ",maxCoreg,"\nNumber of putative GRN per gene : ",
        nGRN,"\nMinimum number of differentially expressed samples to select a single gene : "
        ,minGeneSupport,"\nMinimum number of differentially expressed samples to select a set of co-regulator : "
        ,minCoregSupport,collapse=""))
        
        message("Mining GRN ...")
    }
    
    
    
    result=data.frame()
    gotNet=FALSE
    #just because it's easier toadd here 5% and remove it at the first line in the while loop, where it needs to be decrementale in case no GRNs are found
    searchThresh=  1/((1/searchThresh)-1)
    
    # In very large datasets of very heterogeneous samples (such as the large collection of unrelated cell lines ...)
    # It is possible that no GRN can be fitted with stringent threshold (usually 50%) and that no GRN is found.
    # In case this happens, the threshold is decremented step by step and if no network is found at 10%, then none can be found ...
    while(searchThresh >= 0.05 & !gotNet )
    {
        # decrements the search threshold in case nothing is found
        #(can be the case for VERY large datasets for which it can be hard to find regulators with 50% of matching +1 and -1)
        searchThresh =  1/((1/searchThresh)+1)
        #running hlicorn for each gene in a multithread way if needed.
        if(parallel =="multicore" & length(GeneList)>1 & getOption("mc.cores", 2L) > 1)
        {
            result =mclapply(GeneList,oneGeneHLICORN,geneDiscExp=geneDiscExp,regDiscExp=regDiscExp,
            coregs=coregs,transitemfreq=transitemfreq,transRegBitData=transRegBitData,searchThresh=searchThresh,
            genexp=geneNumExp,regnexp=regNumExp,nresult=nGRN)
            gotNet=TRUE
        }else if(parallel =="snow" & !is.null(cluster) & length(GeneList)>1){
            result =parLapply(cluster,GeneList,oneGeneHLICORN,geneDiscExp=geneDiscExp,regDiscExp=regDiscExp,
            coregs=coregs,transitemfreq=transitemfreq,transRegBitData=transRegBitData,searchThresh=searchThresh,
            genexp=geneNumExp,regnexp=regNumExp,nresult=nGRN)
            gotNet=TRUE
        }else if( length(GeneList)>1){
            result =lapply(GeneList,oneGeneHLICORN,geneDiscExp=geneDiscExp,regDiscExp=regDiscExp,
            coregs=coregs,transitemfreq=transitemfreq,transRegBitData=transRegBitData,searchThresh=searchThresh,
            genexp=geneNumExp,regnexp=regNumExp,nresult=nGRN)
            gotNet=TRUE
        }else{
            result=oneGeneHLICORN(GeneList,geneDiscExp,regDiscExp,coregs,transitemfreq,transRegBitData,searchThresh ,
            genexp=geneNumExp,regnexp=regNumExp,nresult=nGRN)
            gotNet=TRUE
        }
        
    }
    #if one of the above call to LICORN worked ... and if the result is a list (neither matrix nor dataframe)
    #  Merge the results into a data.Frame
    if(      gotNet){
        #if needed, like when used in parallel, merge the results into a data.frame
        if(!is.data.frame(result) & !is.matrix(result)){result= data.frame(do.call(rbind,result))}
        #if LICORN actually did find some networks ... (meaning at least one GRN)
        if(ncol(result) >= 3 & nrow(result) >0){
            # Maybe LICORN did find somes nets, but not enough .. (for less then 5% of the genes)
            if(length(unique(result$Target)) < (0.05*length(GeneList))){
                gotNet=FALSE
            }
        }else{
            gotNet=FALSE
        }
        if(verbose){print(paste("got",nrow(result),"grn"))}
    }
    
    if(verbose){
        message(paste("adjusted thresh:", searchThresh))
    }
    
    
    #When done decrementing the threshold .... well if nothing was found maybe there is a probleme somewhere ...
    if(nrow(result) ==0 | ncol(result) <3){
        stop("Something went wrong. No GRN found.")
    }
    rownames(result)=NULL
    sigrns = coregnet(result)
    sigrns@inferenceParameters=list(minGeneSupport=minGeneSupport,maxCoreg=maxCoreg,minCoregSupport = minCoregSupport,searchThresh=searchThresh,nGRN=nGRN)
    return(sigrns)
    
    
}


oneGeneHLICORN = function(g,geneDiscExp,regDiscExp,coregs,transitemfreq,transRegBitData,searchThresh,regnexp,genexp,nresult){
    shift=ncol(geneDiscExp)
    #sample index with the target gene at 1 or -1
    pos =which( geneDiscExp[g,]==1)
    # negative samples are shifted because we are using a binary matrix with true false for ones in the first part
    # and true false for -1 in the second part
    neg=which( geneDiscExp[g,]== -1)+shift
    
    # select all the coregulators with a support of 50% minimum only in the samples with the target gene at ones or minus ones
    coact=coregs[which(support(transitemfreq, transRegBitData[c(pos,neg)])>= searchThresh )]
    pos =pos+shift
    neg=neg - shift
    corep=coregs[which(support(transitemfreq, transRegBitData[c(pos,neg)]) >= searchThresh )]
    
    # add empty coregulators to have the possibility to only have ativators or inhibitors
    corep = c(corep,list(""))
    coact = c(coact,list(""))
    
    # to have unique coregulators and a single vector of coreg (not a list)
    coactnames =unique(sapply(lapply(coact,sort),paste,collapse=" "))
    coact=strsplit(coactnames," ")
    coact[[which(coactnames=="")]]=""
    corepnames =unique(sapply(lapply(corep,sort),paste,collapse=" "))
    corep=strsplit(corepnames," ")
    corep[[which(corepnames=="")]]=""
    
    
    # merge merge expression of coregulator and corepressor
    coactexp = eand(coact,regDiscExp)
    corepexp = as.integer(eand(corep,regDiscExp))
    # active inhibitor has a stronger impact than naything else in licorn:
    corepexp[which(corepexp ==1)] = 2
    
    
    x= .C("combnLicorn",
    as.integer(t(coactexp)),as.integer(length(coact)),                     #expression and number of coactivators
    as.integer(t(corepexp)),as.integer(length(corep)),                     #expression and number of corepressors
    as.integer(geneDiscExp[g,]),      as.integer(ncol(geneDiscExp)),    #expression of gene, number of samples
    as.double(rep(-1,(length(coact))*(length(corep)))),             # vector to store MAE results
    as.integer(rep(-1,(length(coact))*(length(corep)))),            # vector to store index of coactivator
    as.integer(rep(-1,(length(coact))*(length(corep))))             # vector to store index of corepressor
    )
    
    # bad index will store all bad comparisons (could be done before computing .. right?)
    # then, no intersection between act and rep (includes both empty
    goodindex=which(apply(cbind(x[[8]],x[[9]]),1,function(y){
        return(length(intersect( coact[[y[1]]],corep[[y[2]]] )))
    })==0)
    
    
    selact = coactnames[x[[8]][goodindex]]
    selrep = corepnames[x[[9]][goodindex]]
    
    # all emty set of coregulators are set to NA
    selact[which(selact=="")]=NA
    selrep[which(selrep=="")]=NA
    
    mae = x[[7]][goodindex]
    if(!is.na(nresult)){
        # get 100 first ranks, if ties, might get more ...
        bestindex= which(rank(mae,ties.method="min")<=nresult)
        GRN = data.frame("Target"=rep(g, length(bestindex)),"coact"=selact[bestindex],   "corep"=selrep[bestindex] ,stringsAsFactors=FALSE)      
    }else{
        bestindex= which(rank(mae,ties.method="min")==1)
        GRN = data.frame("Target"=rep(g, length(bestindex)),"coact"=selact[bestindex],   "corep"=selrep[bestindex] ,stringsAsFactors=FALSE)      
        
    }
  
    # if no grn are found return NULL
    if(nrow(GRN)==0){
        return(NULL)
    }
    
    if(is.matrix(GRN) | is.data.frame(GRN)){
        linearmodels=.getEntry(apply(GRN,1,.fitGRN,genexp=t(genexp),regexp=t(regnexp),permut=FALSE),"numscores")
    }else{
        #    linearmodels=.linearCoregulationBootstrap(as.character(GRN),genexp=gexp,regnexp=regnexp,numBootstrap=numBootstrap)
        linearmodels=.fitGRN(as.character(GRN),genexp=t(genexp),regexp=t(regnexp),permut=FALSE)$numscores
    }
    
    numscores=data.frame(t(linearmodels),stringsAsFactors = FALSE)
    
    numscores[,3]=as.numeric(numscores[,3])
    numscores[,4]=as.numeric(numscores[,4])
    numscores[,5]=as.numeric(numscores[,5])
    numscores[,6]=as.numeric(numscores[,6])
    colnames(GRN)=c("Target","coact","corep")
    
    return(data.frame(GRN,numscores,stringsAsFactors = FALSE))
    
}



eand = function(coact,regDiscExp,multip=1){
    do.call(rbind,lapply(coact,function(co){
        if(co[1] == ""){
            return(rep(0,ncol(regDiscExp)))
        }else if(length(co) ==1){return(regDiscExp[co,])}
        n =length(co)
        x=apply(regDiscExp[co,],2,sum)
        y=x
        x[1:length(x)]=0
        x[which(y== - n )]=- multip
        x[which(y == n)] = multip
        return(x)
    }))
}











undirectedNetworkEnrichment = function(net1,net2,commonTF=NULL,commonGene=NULL,verbose=TRUE){
    
    if(class(net1)=="coregnet"){
        net1=coregulators(net1)
    }
    if(class(net2)=="coregnet"){
        net2=coregulators(net2)
    }
    
    if(is.data.frame(net1) | is.matrix(net1) )
    {
        x=as.character(net1[,1])
        y=as.character(net1[,2])        
        net1=unique(data.frame("A1"=c(x,y),"A2"=c(y,x),stringsAsFactors=FALSE))
            
    }else{
        stop("wrongformat")
    }
    if(is.data.frame(net2) | is.matrix(net2) )
    {
        x=as.character(net2[,1])
        y=as.character(net2[,2])
        net2=unique(data.frame("A1"=c(x,y),"A2"=c(y,x),stringsAsFactors=FALSE))

    }else{
        stop("wrong format ...")
    }
    
      
  
    
    common=intersect(c(net1$A1,net1$A2),c(net2$A1,net2$A2))
    
    
    net1 = net1[which(net1$A1 %in% common & !is.na(net1$A1) & net1$A2 %in% common & !is.na(net1$A2)& net1$A1 !=net1$A2),]
    net2 = net2[which(net2$A1 %in% common & !is.na(net2$A1) & net2$A2 %in% common & !is.na(net2$A2)& net2$A1 !=net2$A2),]
    int1 = apply(net1,1,paste,collapse="_")
    int2 = apply(net2,1,paste,collapse="_")
    comint=length(intersect(int1,int2)) /2
    nint1 = (length(int1)/2) - comint
    nint2 = (length(int2)/2) - comint
    totint = (length(common) *(length(common) -1))/2
    if(verbose){print(    matrix(c((comint),nint1,nint2,(totint)-nint1-nint2+comint),nrow=2))
    print(paste("Common",length(common)))}
    return(unlist(fisher.test(   matrix(c((comint),nint1,nint2,(totint)-nint1-nint2+comint),nrow=2) )[c("estimate","p.value")]))
    
}


directedNetworkEnrichment = function(net1,net2,commonTF=NULL,commonGene=NULL,verbose=FALSE){
  
  if(is.data.frame(net1) | is.matrix(net1) )
  {
    if(ncol(net1)==2){
      net1=data.frame(net1,stringsAsFactors=FALSE)
      colnames(net1)=c("Regulator","Target")
    }else{
      stop("net1 is a matrix or data.fram with more than 2 column. Need an edge or adjacency list")
    }
  }
  if(is.data.frame(net2) | is.matrix(net2) )
  {
    if(ncol(net2)==2){
      net2=data.frame(net2,stringsAsFactors=FALSE)
      colnames(net2)=c("Regulator","Target")
    }else{
      stop("net2 is a matrix or data.fram with more than 2 column. Need an edge or adjacency list")
    }
  }
  
  if(class(net1)=="coregnet"){
    net1= coregnetToDataframe(net1)
  }
  if(class(net2)=="coregnet"){
    net2=coregnetToDataframe(net2)
  }
  
  if( class(net1)=="list" ){
    net1=data.frame("Regulator"=unlist(lapply(names(net1),function(y){rep.int(y,length(net1[[y]]))})),
                    "Target"=unlist(net1,use.names=FALSE),stringsAsFactors=FALSE)
  }
  if(class(net2)=="list"){
    net2=data.frame("Regulator"=unlist(lapply(names(net2),function(y){rep.int(y,length(net2[[y]]))})),
                    "Target"=unlist(net2,use.names=FALSE),stringsAsFactors=FALSE)
    
  }
  
  if(is.null(commonTF)){
    commonTF=intersect((net1$Regulator),(net2$Regulator))
  }
  if(is.null(commonGene)){
    commonGene=intersect((net1$Target),(net2$Target))
  }
  net1 = net1[which(net1$Regulator %in% commonTF & !is.na(net1$Regulator) & net1$Target %in% commonGene & !is.na(net1$Target)),]
  net2 = net2[which(net2$Regulator %in% commonTF & !is.na(net2$Regulator) & net2$Target %in% commonGene & !is.na(net2$Target)),]
  int1 = apply(net1,1,paste,collapse="_")
  int2 = apply(net2,1,paste,collapse="_")
  comint=length(intersect(int1,int2))
  if(verbose){print(    matrix(c((comint),length(int1)-comint,length(int2)-comint,(length(commonTF)*length(commonGene))-length(int1)-length(int2)+comint),nrow=2))}
  return(unlist(fisher.test(    matrix(c((comint),length(int1)-comint,length(int2)-comint,
              (length(commonTF)*length(commonGene))-length(int1)-length(int2)+comint),nrow=2))[c("estimate","p.value")]))
  
}






.invlist = function(y){    w =rep.int(names(y),sapply(y,length))
                           z = unlist(y,use.names=FALSE)
                           return(tapply(w,z,c))}




.subsigrns = function(sigrns,genes)
{
  
  subtg  =intersect(names(sigrns$bygene),genes)
  subtf = intersect(names(sigrns$bytf),genes)
  subsigrns = list()
  
  subsigrns$bytf = lapply(sigrns$bytf[subtf],function(subr){
    subn = list()
    subn$act = intersect(   subr$act  , subtg   )
    subn$rep = intersect(subr$rep,subtg)
    return(subn)
  })
  
  subsigrns$bygene = lapply(sigrns$bygene[subtg],function(subr){
    subn = list()
    subn$act = intersect(   subr$act  , subtf   )
    subn$rep = intersect(subr$rep,subtf)
    return(subn)
  })
  return(subsigrns)
}
