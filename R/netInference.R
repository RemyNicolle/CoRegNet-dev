
# a function ...
.quicknonuniqgrnsTOSIGRNS = function(sigrns)
{
    act = strsplit(sigrns[,2]," ")
    rep = strsplit(sigrns[,3]," ")
    is = 1:length(rep)
    
    sigrns=unique(sigrns)
    goodgrnsIndex = mclapply(is,function(x,act,rep){
        if(length(intersect(act[[x]],rep[[x]])) == 0) {return (x);}
        else {return(0)}
    },act=act,rep=rep)
    goodgrnsIndex = unlist(goodgrnsIndex,use.names=FALSE)
    sigrns = sigrns[goodgrnsIndex,]
    
    SIGRNS = list()
    #SIGRNS$GRN = sigrns
    genes = unique(sigrns[,1]  )
    tf = unique(c(unique(unlist(strsplit(sigrns[,2]," "))),
    unique(unlist(strsplit(sigrns[,3]," ")))))
    tf = tf[which(tf!="EMPTY" & !is.na(tf) & tf != "NA")]
    bygene = mclapply(genes,function(x){
        y = sigrns[which(sigrns[,1] == x),]
        sub = list()
        sub$act = unique(unlist(strsplit(y[,2]," ")))
        sub$rep = unique(unlist(strsplit(y[,3]," ")))
        sub$act=sub$act[which(sub$act !="EMPTY" & !is.na(sub$act) & sub$act!="NA")]
        sub$rep=sub$rep[which(sub$rep !="EMPTY"& !is.na(sub$rep) & sub$rep!="NA")]
        return(sub)
    })
    names(bygene) = genes
    SIGRNS$bygene = bygene
    bytf = mclapply(tf,function(x){
        sub = list()
        sub$act=c()
        sub$rep = c()
        sub$act = unique(sigrns[grep(paste("\\b",x,"\\b",sep=""),sigrns[,2]),1])
        sub$rep = unique(sigrns[grep(paste("\\b",x,"\\b",sep=""),sigrns[,3]),1])
        return(sub)
    })
    names(bytf) = tf
    head(bytf)
    length(bytf)
    SIGRNS$bytf = bytf
    return(list("sigrns"=sigrns,"adjList"=SIGRNS))
}





.fitGRN=function(grn,genexp,regexp,permut=FALSE,CVlm=TRUE)
{
  
  #get all regulators (predictors) and gene (response) data and organise in a new dataframe
  #deal with several, one or no coregulators
  act = unique(unlist(strsplit(grn[2]," ")))
  act=act[which(act!="EMPTY" & !is.na(act))]
  if(permut){act = sample(colnames(regexp),length(act))}
  rep = unique(unlist(strsplit(grn[3]," ")))
  rep=rep[which(rep!="EMPTY" & !is.na(rep))]
  if(permut){rep = sample(colnames(regexp),length(rep))}
  X = regexp[,c(act,rep)]
  if(length(act)>1){
    coact = apply(regexp[,act],1,prod)
    X=data.frame(X,"coact"=as.numeric(coact))
  }
  if(length(rep)>1){
    corep = apply(regexp[,rep],1,prod)
    X=data.frame(X,"corep"=as.numeric(corep))
  }
  
  # da contains all the necessary data with the first column y as the gene, response variable
  Y=genexp[,grn[1]]
  da = data.frame("y"= Y,X)
  
  # get 10-cross validated coefficients and fitted values
  if(CVlm){
    Icol=sample(1:nrow(da))
    Ks=as.integer(cut(1:nrow(da),10))
    da=da[Icol,]
    cvlm=lapply(1:10,function(k){
      daL = da[which(Ks != k),]
      daP = da[which(Ks == k),]
      l=lm("y~.",daL)
      return(list("fitted"=predict(l,daP),"coefs"=coef(l)))
    })
    
    fitted=unlist( lapply(cvlm,function(x){return(x$fitted)}))[order(Icol)]
    coefs=apply(sapply(cvlm,function(x){return(x$coefs)}),1,mean)
    coefs= coefs[2:length(coefs)]
    
    R2=as.numeric(cor(Y,fitted)^2 )
    residues = Y-fitted
    RMSE = sqrt(sum(( (residues) ^2) )/length( fitted))
  }else{
    l=lm("y~.",da)
    fitted=l$fitted.values
    coefs=coef(l)
    coefs= coefs[2:length(coefs)]
    residues = Y-fitted
    RMSE = sqrt(sum(( (residues) ^2) )/length( fitted))
    R2 =  summary(l)$r.squared
  }
  if(permut){
    return(c(R2,RMSE))
  }else{
    numscores =c(rep.int(0,4),R2,RMSE)
    names(numscores)=c("Coef.Acts","Coef.Reps","Coef.coActs","Coef.coReps","R2","RMSE")
    
    if(length(act)==0){
      numscores["Coef.Acts"]=NA
      numscores["Coef.coActs"]=NA
    }else if(length(act)==1){
      numscores["Coef.Acts"]=coefs[1]
      numscores["Coef.coActs"]=NA
    }else{
      numscores["Coef.Acts"]=paste(coefs[1:length(act)],collapse= " ")
      numscores["Coef.coActs"]=coefs["coact"]
    }
    
    if(length(rep)==0){
      numscores["Coef.Reps"]=NA
      numscores["Coef.coReps"]=NA
    }else if(length(rep)==1){
      numscores["Coef.Reps"]=coefs[length(act)+1]
      numscores["Coef.coReps"]=NA
    }else{
      numscores["Coef.Reps"]=paste(coefs[(length(act)+1):length(rep)],collapse= " ")
      numscores["Coef.coReps"]=coefs["corep"]
    }
    return(list("fitted"=fitted,"residuals"=residues,"numscores"=numscores))
  }
}

#Very simple function, maybe already implemented in R. Simplifies getting data out of a list.
.getEntry = function(x,name){
    return(sapply(x,function(y){return(y[[name]])}))
}


.adaptiveCoregSupport = function(regexp,maxCoreg=floor(log(ncol(regexp),base=3)),thresholds=c(0.5,0.45,0.4,0.35,0.3,0.25,0.2,0.15,0.1),maxCoregSets=10^5){
    
    trans = list()
    for( i in 1:ncol(regexp)){
        trans[[paste(i,"+",sep="")]] = names(which(regexp[,i] ==1))
    }
    for( i in 1:ncol(regexp)){
        trans[[paste(i,"-",sep="")]] = names(which(regexp[,i] ==-1))
    }
    
    adapthresh=thresholds[1]
    for( th in thresholds){
        # again, the thresholds of the minimum number of samples (support) is divided by 2 because each of the samples are counted twice (positive TF and negative TF)
        if(length(suppressWarnings(apriori(trans,parameter=list(support =th/2,minlen=2,maxlen=maxCoreg,target="closed frequent itemsets")
        ,control=list(verbose=FALSE)))) > maxCoregSets){
            return(adapthresh)
        }else{
            adapthresh=th
        }
        
    }
    return(adapthresh)
    
}

