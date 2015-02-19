# returns the size of the intersection of two network in the form of adjacency list with genes as entries (names of lists)
.netintersect = function(grn1,grn2){
    #list with target genes as entry
    genes = intersect(names(grn1),names(grn2))
    return(sum(sapply(genes,function(g){
        length(intersect(unlist(grn1[[g]]),unlist(grn2[[g]])))
    })))
    
}


# given a data frame containing interactions (regulator -> gene, one per row)
# score the GRN and add the data in the coregnet object
.addOneRegulatoryEvidence=function(geneRegulatoryNetwork,evidenceData,evname){
    ################################################################################################################
    #controlling input
    # in the end, evlist should be an adjacency list and evidenceData a 2 (or more) column data frame.
    if(class(geneRegulatoryNetwork) != "coregnet"){
        stop("Uses only coregnet network as an input") }
    if(is.data.frame(evidenceData) | is.matrix(evidenceData)){
        if(ncol(evidenceData)<2){
            stop("Evidence data should be either list or a two column data.frame (or matrix). Refer to the manual.")
        }else{
            # if evidence is in a matrix or a data frame, emove duplicate interactions, if possible by keeping the interaction with lowest "p.value"
            evidenceData[,1]=as.character(evidenceData[,1])
            evidenceData[,2]=as.character(evidenceData[,2])
            # removes duplicated evidence. If a pvalue column is present, keeps only the lowest pvalue. Otherwise only keeps the two
            evidenceData=.uniqueEvidenceDataFrameByPvalue(evidenceData)
            evlist = tapply(evidenceData[,1],evidenceData[,2],function(i){return(i)})  }
    }else{
        stop("Evidence data should be  a two column matrix or data.frame. Refer to the manual.")
    }
    
    ################################################################################################################
    # Intersection of gene and regulator names
    # keep only the evidences about genes and regs that are in the network
    originalGenes =length(unique(evidenceData[,2]))
    originalReg =length(unique(evidenceData[,1]))
    originalInt = nrow(evidenceData)
    genesandregs=c(names(geneRegulatoryNetwork@adjacencyList$bygene),names(geneRegulatoryNetwork@adjacencyList$bytf))
    commongenes = length(intersect(names(evlist),names(geneRegulatoryNetwork@adjacencyList$bygene)))
    #print(commongenes)
    commonreg = length(intersect(unique(unlist(evlist)),names(geneRegulatoryNetwork@adjacencyList$bytf)))
    if(commongenes==0 |  commonreg ==0){   return(NULL)}
    evlist = evlist[intersect(names(evlist),names(geneRegulatoryNetwork@adjacencyList$bygene))]
    geneNames = names(evlist)
    evlist = lapply(evlist,intersect,y=names(geneRegulatoryNetwork@adjacencyList$bytf))
    names(evlist)   =geneNames
    evlist =evlist[which(sapply(evlist,length)>0)]
    evidenceData = evidenceData[which(evidenceData[,1] %in% genesandregs & evidenceData[,2] %in% genesandregs),]
    
    ################################################################################################################
    # Scoring GRNs
    # count for each grn the number of known evidences and normalize.
    evcount = .parallelRegulationEvidence(geneRegulatoryNetwork@GRN,  evlist)
    if(sum(evcount) == 0){
        print(paste("No evidence from",evname,"were found in the inferred network."))
        return(NULL)
    }
    geneRegulatoryNetwork@GRN=data.frame(geneRegulatoryNetwork@GRN,evcount)
    colnames(geneRegulatoryNetwork@GRN)[ncol(geneRegulatoryNetwork@GRN)]=evname
    
    
    ################################################################################################################
    #
    # add evidence data to the object of type coregnet
    # evidenceData = cbind(evidenceData,matrix(NA,nrow=nrow(evidenceData),ncol=(4-ncol(evidenceData))))
    colnames(evidenceData)[1:2] = c("reg","target")
    geneRegulatoryNetwork@evidences[[evname]] = evidenceData
    
    desc=.descriptionUpdate(geneRegulatoryNetwork,evname,"regulatory")
    eviDesc = c("regulatory",originalGenes,originalReg,originalInt,desc)
    if(nrow(geneRegulatoryNetwork@evidenceDescription) == 0){
        geneRegulatoryNetwork@evidenceDescription = data.frame(t(eviDesc),stringsAsFactors=FALSE)
        colnames(geneRegulatoryNetwork@evidenceDescription)=c("evidenceType","originalGene","originalReg","originalEvidence",
        "commonGene","commonReg","evidences","commonEvidences", "enrichment","p.value")
        rownames(geneRegulatoryNetwork@evidenceDescription)[1] = evname
    }else{
        geneRegulatoryNetwork@evidenceDescription = rbind(geneRegulatoryNetwork@evidenceDescription, eviDesc)
    }
    rownames(geneRegulatoryNetwork@evidenceDescription)[nrow(geneRegulatoryNetwork@evidenceDescription)] = evname
    
    return(geneRegulatoryNetwork)
}


# given a data frame containing interactions (regulator -> gene, one per row)
# score the GRN and add the data in the coregnet object
.addOneCoRegulatoryEvidence=function(geneRegulatoryNetwork,evidenceData,evname){
    ################################################################################################################
    #controlling input
    # in the end, evlist should be an adjacency list and evidenceData a 2 (or more) column data frame.
    if(is.data.frame(evidenceData) | is.matrix(evidenceData)){
        if(ncol(evidenceData)<2){
            stop("Evidence data should be either list or a two column data.frame (or matrix). Refer to the manual.")
        }else{
            # if evidence is in a matrix or a data frame, emove duplicate interactions, if possible by keeping the interaction with lowest "p.value"
            evidenceData[,1]=as.character(evidenceData[,1])
            evidenceData[,2]=as.character(evidenceData[,2])
            # removes duplicated evidence. If a pvalue column is present, keeps only the lowest pvalue. Otherwise only keeps the two
            evidenceData=.uniqueEvidenceDataFrameByPvalue(evidenceData)
            # Need the set of interactions in the form of an adjancecy list.
            # Because the interaction are undirected and to ease the computations, all edges are listed twice (A to B and B to A)
            ppi = as.matrix(unique(evidenceData[,1:2]))
            ppi = unique(rbind(ppi[,1:2],ppi[,2:1]))
            ppi= unique(ppi[which(ppi[,1] != ppi[,2]),])
            evlist = tapply(ppi[,1],ppi[,2],function(i){return(i)})  }
    }else{
        stop("Evidence data should be  a two column matrix or data.frame. Refer to the manual.")
    }
    
    ################################################################################################################
    # Intersection of gene and regulator names
    # keep only the evidences about genes and regs that are in the network
    originalReg =length(unique(unlist(c(evidenceData[,1],evidenceData[,2]))))
    originalInt = nrow(ppi)  /2
    genesandregs=c(names(geneRegulatoryNetwork@adjacencyList$bygene),names(geneRegulatoryNetwork@adjacencyList$bytf))
    commonreg = length(intersect(unique(unlist(evlist)),names(geneRegulatoryNetwork@adjacencyList$bytf)))
    if(  commonreg ==0){   return(NULL)}
    evlist = evlist[intersect(names(evlist),names(geneRegulatoryNetwork@adjacencyList$bytf))]
    geneNames = names(evlist)
    evlist = lapply(evlist,intersect,y=names(geneRegulatoryNetwork@adjacencyList$bytf))
    names(evlist)   =geneNames
    evlist =evlist[which(sapply(evlist,length)>0)]
    evidenceData = evidenceData[which(evidenceData[,1] %in% genesandregs & evidenceData[,2] %in% genesandregs),]
    
    ################################################################################################################
    # Scoring GRNs
    # count for each grn the number of known evidences and normalize.
    evcount = .parallelRegulatorCooperativityEvidence(geneRegulatoryNetwork@GRN,  evlist)
    if(sum(evcount) == 0){
          print(paste("No evidence from",evname,"were found in the inferred network."))
          return(NULL)
    }
    
    geneRegulatoryNetwork@GRN=data.frame(geneRegulatoryNetwork@GRN,evcount)
    colnames(geneRegulatoryNetwork@GRN)[ncol(geneRegulatoryNetwork@GRN)]=evname
    
    
    ################################################################################################################
    #
    # add evidence data to the object of type coregnet
    # evidenceData = cbind(evidenceData,matrix(NA,nrow=nrow(evidenceData),ncol=(4-ncol(evidenceData))))
    colnames(evidenceData)[1:2] = c("reg1","reg2")
    geneRegulatoryNetwork@evidences[[evname]] = evidenceData
    
    desc=.descriptionUpdate(geneRegulatoryNetwork,evname,"coregulatory")
    eviDesc = c("coregulatory",NA,originalReg,originalInt,desc)
    if(nrow(geneRegulatoryNetwork@evidenceDescription) == 0){
        geneRegulatoryNetwork@evidenceDescription = data.frame(t(eviDesc),stringsAsFactors=FALSE)
        colnames(geneRegulatoryNetwork@evidenceDescription)=c("evidenceType","originalGene","originalReg","originalEvidence",
        "commonGene","commonReg","evidences","commonEvidences", "enrichment","p.value")
        rownames(geneRegulatoryNetwork@evidenceDescription)[1] = evname
    }else{
        geneRegulatoryNetwork@evidenceDescription = rbind(geneRegulatoryNetwork@evidenceDescription, eviDesc)
    }
    rownames(geneRegulatoryNetwork@evidenceDescription)[nrow(geneRegulatoryNetwork@evidenceDescription)] = evname
    
    return(geneRegulatoryNetwork)
}





.descriptionUpdate = function(regnet,evname,type="regulatory"){
    ev=regnet@evidences[[evname]]
    
    if(type=="regulatory"){  
      net=coregnetToDataframe(regnet)
      commongenes = intersect(unique(ev[,2]),unique(net[,2]))
      commonreg = intersect(unique(ev[,1]),unique(net[,1]))
      ev=ev[which(ev[,1] %in% commonreg & ev[,2] %in% commongenes),1:2]
      net=net[which(net[,1] %in% commonreg & net[,2] %in%  commongenes),1:2]
      
      evg=igraph::graph.data.frame(ev)
      g=igraph::graph.data.frame(net)
      n21=length(igraph::E(igraph::graph.difference(evg, g)))
      n12=length(igraph::E(igraph::graph.difference(g,evg )))
      n11=length(igraph::E(igraph::graph.intersection(g,evg)))
      n22=(length(commonreg)*length(commongenes))-n11-n12-n21      
      
      f=fisher.test(matrix(c(n11,n12,n21,n22),ncol=2),alternative="great")
         return(c(length(commongenes),length(commonreg),length(igraph::E(evg)),n11,f$estimate ,f$p.value))
    }else{      
      coreg=coregulators(regnet)
      alltfs = unique(c(as.character(unlist(coreg[,1])),as.character(unlist(coreg[,2]))))
      coregG =igraph::simplify(igraph::graph.edgelist(as.matrix(coreg[,1:2]),directed=FALSE))
      g=igraph::simplify(igraph::graph.edgelist(as.matrix(regnet@evidences[[evname]][,1:2]),directed=FALSE))
      
      TF = intersect(igraph::V(g)$name,alltfs)
      subg=igraph::induced.subgraph(g, which(V(g)$name %in%TF))
      subcoreg=igraph::induced.subgraph(coregG, which(igraph::V(coregG)$name %in%TF))
      
      n11=length(igraph::E(igraph::graph.intersection(subcoreg,subg)))    
      n21=length(igraph::E(igraph::graph.difference(subg, subcoreg)))   
      n12=length(igraph::E(igraph::graph.difference(subcoreg,subg )))      
      n22=((length(TF)*(length(TF)-1))/2)-n11-n12-n21
      #print(matrix(c(n11,n12,n21,n22),ncol=2))
      f=fisher.test(matrix(c(n11,n12,n21,n22),ncol=2),alternative="great")
        return(c(NA,length(TF),length(igraph::E(subg)),n11,f$estimate ,f$p.value))
    }
}


.uniqueEvidenceDataFrameByPvalue =function(evdataframe){
    # try possible col names that would correspond to pvalues
    
    possiblepvaluecolumnnames=c("adj.p.value","adj.pvalue","p.value")
    
    #uniforming column names, setting to lower case and changin - to .
    colnames(evdataframe) =tolower(colnames(evdataframe))
    colnames(evdataframe)=gsub("-",".",colnames(evdataframe))
    # if only two columns, should be two character column, first with Regulators second with target genes.
    if(ncol(evdataframe) ==2){
        return(unique(evdataframe))
        #if one of the additional columns
    }else if(ncol(evdataframe)>2 & sum(possiblepvaluecolumnnames %in% colnames(evdataframe) )>0 ){
        pvalcol=intersect(possiblepvaluecolumnnames,colnames(evdataframe))[1]
        interactionsonly = apply(evdataframe[,1:2],1,paste,collapse=" ")
        dupsint = interactionsonly[which(base::duplicated(interactionsonly))]
        if(length(dupsint)==0){
            return(evdataframe)
        }
        nondupinteractions = evdataframe[interactionsonly[which(!interactionsonly %in% dupsint)],]
        minpval=lapply(dupsint,function(oneInt){
            dupev=evdataframe[which(interactionsonly == oneInt),]
            return(dupev[which.min(dupev[,pvalcol])])
        })
        minpval=do.call(rbind,minpval)
        return(rbind(nondupinteractions,minpval))
    }else{
        # if no pvalue column, test if there are duplicates.
        if(nrow(evdataframe) > nrow(unique(evdataframe[,1:2]))){
            return(unique(evdataframe[,1:2]))
        }else{
            return(evdataframe)
        }
        
    }
}




.parallelRegulatorCooperativityEvidence  <- function(grns,ppi,targetGene = NULL )
{
    
    
    if(is.null(targetGene)) targetGene = unique(grns[,1]);
    
        evidence =suppressWarnings( pvec(targetGene,function(x,ppi,grns)    {
            
            ex =  .regulatorCooperativityEvidence(grns[which(grns[,1] %in% x),],ppi)
            return(unlist(ex))
        },   ppi=ppi,grns=grns))
        
    
    
    if(length(evidence) != nrow(grns)){stop("Evidence and GRN have different size. ERROR !")}
    return(unlist(evidence))
}

.regulatorCooperativityEvidence <- function(grns,ppi)
{
    res = apply(grns,1,.addRegulatorCooperativityEvidence,ppi=ppi)
    return(res)
}

.addRegulatorCooperativityEvidence <- function(x,ppi)
{
    
    activ=c()
    inhib=c()
    activators = as.character(x[2])
    activ = unlist(strsplit(activators, " "))
    inhibitors = as.character(x[3])
    inhib = unlist(strsplit(inhibitors, " "))
    
    if( length(activ )> 1) {
        activppiKnown =sum(unlist(ppi[activ],use.names=FALSE) %in% activ)
        possibleAPPI = ((length(activ)^2)-length(activ) )
    }else {
        activppiKnown=0
        activatorCorrectedEvidence=0
        possibleAPPI=0
    }
    
    if(length(inhib )> 1) {
        inhibppiKnown =sum(unlist(ppi[inhib],use.names=FALSE) %in% inhib)
        possibleIPPI=((length(inhib)^2)-length(inhib) )
    }else {
        inhibitorCorrectedEvidence=0
        inhibppiKnown=0
        possibleIPPI=0
    }
    
    if(possibleIPPI+possibleAPPI == 0){return(0)}else{
      return( (inhibppiKnown+activppiKnown)  / (possibleAPPI+possibleIPPI) )
    }

}




.parallelRegulationEvidence <- function(grns,tfgi,targetGene = NULL ,multi=TRUE)
{
    if(is.null(targetGene)) targetGene = unique(grns[,1]);
    
    if(multi){
        evidence =     suppressWarnings( pvec(targetGene,function(gene,tfgi,grns)
        {
            .regulationEvidence(grns[which(grns[,1] %in% gene),],tfgi[gene])
        }
        ,tfgi=tfgi,grns=grns))
    }else{
        evidence = lapply(targetGene,function(gene,tfgi,grns)
        {
            .regulationEvidence(grns[which(grns[,1] %in% gene),],tfgi[gene])
        }
        ,tfgi=tfgi,grns=grns)
    }
    if(length(evidence) != nrow(grns)){stop("Evidence and GRN have different size. ERROR !")}
    return(unlist(evidence))
}

.regulationEvidence <- function(grns,tfgi) apply(grns,1,.addRegulationEvidence,tfgi=tfgi);

.addRegulationEvidence <- function(x,tfgi)
{
    
    target = as.character(x[1])
    activ=c()
    inhib=c()
    activators = as.character(x[2])
    activ = unlist(strsplit(activators, " "))
    inhibitors = as.character(x[3])
    inhib = unlist(strsplit(inhibitors, " "))
    knownReg = unique(unlist(tfgi[[target]], use.names=FALSE))
    nreg=0
    
    if( ! "EMPTY" %in% activ | sum(is.na(activ))==0 ) 
    {
        activKnown = length(intersect(knownReg ,activ ))
        nreg=nreg+length(activ)
    }else 
    {
        nreg=nreg+0      
        activKnown=0  
    }
    
    if( ! "EMPTY" %in%inhib | sum(is.na(inhib))==0) 
    {
        inhibKnown = length(intersect(knownReg ,inhib ))
        nreg=nreg+length(inhib)
    }else {
        nreg=nreg+0  
        inhibKnown=0
    }
    
    return ((activKnown+inhibKnown)/nreg)
}



