
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


# load("~/workspace/coRegNet/TCGA-bladder/BLCArnaseq.271.RData")

# load("~/workspace/coRegNet/CIT/CITtumsampandNormalsEXP.RData")
# expression = exp
# library(CoRegNet)
# library(parallel)
# data(HumanTF)
# TF =HumanTF
# nThreshQuantile = 100
# scaled=F
# verbose=T



automaticParameters = function(expression,TF,nThreshQuantile = 1000,scaled=FALSE,verbose=FALSE){
	#work only on TF
    subtf = expression[intersect(rownames(expression),TF),]
	# scale expression if not already done
    if(!scaled){
        subtf = t(scale(t(subtf),scale=F))
	}else{
			 subtf=as.matrix(subtf)
	}
	
    #generate random TF expression
		#     randomtf = mclapply(1:sampling,function(i){
		# return(t(apply(subtf,1,function(x){x[sample(1:ncol(subtf))]})))
		#     })
	#define several thresholds to test
	thresholds= quantile(abs(as.vector(subtf)),probs=(1:(nThreshQuantile-1))/nThreshQuantile)
    nco = ncol(subtf)
    nro = nrow(subtf)	
	#for each threshold, compare the random number of pairs in each sample versus the real number of pairs
	discTests=mclapply(thresholds,function(th){
		print(th)
	    realSup=meanSupport(subtf,th,nco,nro)
		prob= meanSupport(subtf,th,nco,nro,getProb=TRUE)
		estimeProb = ((prob)^2)*(((nro*(nro -1))/2))*(2*nco)
		return(c(pnorm(realSup,estimeProb,sqrt(estimeProb)), (realSup-estimeProb)/sqrt(estimeProb) ))
#		# all this is to compare to random data		
#	    return(ppois(realSup,estimeProb,lower.tail=F))
#		if(is.null(realSup)){return(NULL)}
#		randomSup=unlist(lapply(randomtf,meanSupport,th,nco,nro))
		#(realSup-mean(randomSup))/sd(randomSup)  Z-score, kind of...
#		return(c(realSup,estimeProb,max(randomSup),mean(randomSup),sd(randomSup),meanSupport(subtf,th,nco,nro,T)))
	})
	discTests=cbind(do.call(rbind,discTests),thresholds)
	
	# plot the distribution of the z-score of the number of pairs for each threshold
	if(verbose){
		plot(discTests[,3:2],type="l",ylab="z-score",xlab="threshold",main="Deviation from random",
		sub="blue : standard deviation, green 0.5 and 2 times SD.")
	    abline(v=sd(subtf),col="green")
		abline(v=0.5*sd(subtf),col="blue",lty="dotted")
		abline(v=2*sd(subtf),col="blue",lty="dotted")
	}
	
	discTh = discTests[which.max(discTests[,2]),3]
	# if the maximum z-score is too high or too low, just return the standard deviation
	if(discTh > 0.5* sd(subtf) & discTh < 2*sd(subtf)){
		return(discTh)
	}else{
		return(sd(subtf))
	}
	
}


#function to compute the mean Support (in the itemset sense, meaning the number of samples with a non zero value)
# of all pairs of TF given a discretiation threshold
meanSupport=function(centered,threshold,nco,nro,getMeanSup=F,getProb=F){
	regBitData = cbind(	matrix(  centered >= threshold   ,nrow=nro,ncol=nco),
	matrix(  centered <= (-threshold)   ,nrow=nro,ncol=nco))
	if(getProb){
		return(sum(regBitData)/(nro* (2*nco)))
	}
	if(getMeanSup){
		return(mean(apply(regBitData,1,sum)))
	}
	sum(suppressWarnings(apriori(	as(t(regBitData),"transactions"),parameter=list(support = 10^-100,maxlen=2,minlen=2,target="frequent itemsets")
    ,control=list(verbose=FALSE)))@quality[,1] * nco*2)#/ ((nro*(nro -1))/2)

}