
test_discretizeExpressionData <- function() {
  ex=matrix(seq(-2,1.8,by=0.2),nrow=2,dimnames=list(paste("gene",1:2,sep=""),paste("sample",1:10,sep="")))
  checkEquals(as.vector(discretizeExpressionData(ex, threshold=1)), c(rep(-1,6),rep(0,9),rep(1,5)))
}

test_regulatorInfluence<- function()
{
  acts=apply(matrix(rep(letters[1:4],7),nrow=2),2,paste,collapse=" ")[1:13]
  reps=apply(matrix(rep(letters[5:8],7),nrow=2),2,paste,collapse=" ")[1:13]
  grn=data.frame("Target"= LETTERS[1:26] ,"coact"=c(acts,reps),"corep"= c(reps,acts),"R2"=runif(26),stringsAsFactors=FALSE)
  co=coregnet(grn)
  samples= paste("S",1:2,sep="")
  expression=matrix(c(1+(0.01*1:34),1+(0.01*34:1)),ncol=2)
  dimnames(expression) = list(c(grn$Target,names(regulators(co))),samples)
  
  #Minimum number of targets is adjusted because of the small size of the network
  TFA = regulatorInfluence(co,expression,minTarg=4)
  
  
  checkEqualsNumeric(c( regulatorInfluence(co,expression,minTarg=4)), 
                     c(-5.629165 ,-5.629165, -6.017831 ,-6.017831  ,5.629165 ,
                        5.629165 , 6.017831 , 6.017831 , 5.629165 , 5.629165,
                         6.017831 , 6.017831 ,-5.629165, -5.629165, -6.017831 ,-6.017831),
                     tolerance = 0.1)
}

test_regulators <- function(){
  acts=apply(matrix(rep(letters[1:4],7),nrow=2),2,paste,collapse=" ")[1:13]
  reps=apply(matrix(rep(letters[5:8],7),nrow=2),2,paste,collapse=" ")[1:13]
  grn=data.frame("Target"= LETTERS[1:26] ,"coact"=c(acts,reps),"corep"= c(reps,acts),"R2"=runif(26),stringsAsFactors=FALSE)
  co=coregnet(grn)
  checkEquals(regulators(co),c("a"=14,"b"=14,"e"=14,"f"=14,"c"=12,"d"=12,"g"=12,"h"=12))
  checkEquals( regulators(co,target = "A"),c("a","b","e","f"))
  checkEquals( regulators(co,target = "a"),NA)
  
  checkEquals( activators(co,"A"),c("a","b"))
  checkEquals( activators(co,"A","coregulators"),"a b")
  
  checkEquals( repressors(co,"A"),c("e","f"))
  checkEquals( repressors(co,"A","coregulators"),"e f")
  
  checkEquals( targets(co),c("A","C","E","G","I","K","M","N","P","R","T","V","X","Z","B","D","F","H","J","L","O","Q","S","U","W","Y"))
  checkEquals( targets(co,"a"),c("A","C","E","G","I","K","M","N","P","R","T","V","X","Z"))
  checkEquals( targets(co,"a","reg"),c("A","C","E","G","I","K","M","N","P","R","T","V","X","Z"))
  checkEquals( targets(co,"a","act"),c("A","C","E","G","I","K","M"))
  checkEquals( targets(co,"a","rep"),c("N","P","R","T","V","X","Z"))
  checkEquals( targets(co,c("a","b"),"act"),c("A","C","E","G","I","K","M"))
  
}

test_coregulators <- function(){
  acts=apply(matrix(rep(letters[1:4],7),nrow=2),2,paste,collapse=" ")[1:13]
  reps=apply(matrix(rep(letters[5:8],7),nrow=2),2,paste,collapse=" ")[1:13]
  grn=data.frame("Target"= LETTERS[1:26] ,"coact"=c(acts,reps),"corep"= c(reps,acts),"R2"=runif(26),stringsAsFactors=FALSE)
  grn=coregnet(grn)
  co=coregulators(grn)
  checkEquals(co$Reg1,c("a","e","c","g"))
  checkEquals(co$Reg2,c("b","f","d","h"))
  checkEqualsNumeric(co$Support,c(0.2692308,0.2692308,0.2307692,0.2307692),tolerance = 0.001)
  checkEqualsNumeric(co$fisherTest,c(1.035443e-07,1.035443e-07,1.035443e-07,1.035443e-07),tolerance = 0.001)
  checkEqualsNumeric(co$adjustedPvalue,c(1.035443e-07,1.035443e-07,1.035443e-07,1.035443e-07),tolerance = 0.001)
}

test_addEvidences <- function(){
acts=apply(matrix(rep(letters[1:4],4),nrow=2),2,paste,collapse=" ")
reps=apply(matrix(rep(letters[5:8],4),nrow=2),2,paste,collapse=" ")
grn=data.frame("Target"= LETTERS[1:16] ,"coact"=c(acts,reps),"corep"= c(reps,acts),"R2"=runif(16),stringsAsFactors=FALSE)
co=coregnet(grn)
evidence1=unique(data.frame(tf=rep(letters[1:8],times=6),target=rep(LETTERS[1:16],each=3),stringsAsFactors = FALSE))

enrichco=addEvidences(co,evidence1)
checkEqualsNumeric(as.numeric(unlist(enrichco@evidenceDescription[,2:8])),c(16, 8, 48, 16, 8,48,24))


coregevidence1=unique(data.frame(tf1=rep(letters[1:8],times=4),tf2=rep(letters[1:8],each=4),stringsAsFactors = FALSE))
enrichcoop=addCooperativeEvidences(co,coregevidence1)

checkEqualsNumeric(as.numeric(unlist(enrichcoop@evidenceDescription[,2:8])),c(NA, 8, 22, NA, 8,22,4))

}

test_hLICORN <-function(){
    ex=matrix(seq(-3,3,length.out=140),ncol=20,dimnames=list(c(letters[1:2],LETTERS[1:5]),paste("sample",1:20,sep="")))
  ex[5:7,] =-0.6*ex[5:7,]
  dummyNet1=hLICORN(ex,TFlist =LETTERS[1:5] )@GRN
  dummyNet1
  checkEquals( dummyNet1$Target ,rep(c("a","b"),each=2))
  checkEquals( dummyNet1$coact ,c("B","A B","B","A B"))
  #  checkEquals( dummyNet1$corep ,rep( NA,6))
  checkEqualsNumeric(dummyNet1$R2,rep(1,4),tolerance=1.0e-4)
  checkEqualsNumeric(dummyNet1$RMSE,rep(0,4),tolerance=1.0e-4)
  dummyNet2=hLICORN(ex,TFlist =LETTERS[1:5]  ,parallel = "no" )@GRN
  checkEquals(dummyNet1[,1:3],dummyNet2[,1:3])
}

test_refine <-function(){
  acts=apply(matrix(rep(letters[1:4],8),nrow=2),2,paste,collapse=" ")
  reps=apply(matrix(rep(letters[5:8],8),nrow=2),2,paste,collapse=" ")
  grn=data.frame("Target"= rep(LETTERS[1:16],each=2) ,"coact"=c(acts,reps),"corep"= c(reps,acts),"R2"= seq(0.5,0.9,length.out = 32),stringsAsFactors=FALSE)
  co=coregnet(grn)
  evidence1=unique(data.frame(tf=rep(letters[1:8],times=6),target=rep(LETTERS[1:16],each=3),stringsAsFactors = FALSE))
  coregevidence1=unique(data.frame(tf1=rep(letters[1:8],times=4),tf2=rep(letters[1:8],each=4),stringsAsFactors = FALSE))
  
  enrichco=addEvidences(co,evidence1)
  enrichcoop=addCooperativeEvidences(enrichco,coregevidence1)
  
  unsupervisedNet=refine(enrichcoop)
  checkEquals(unsupervisedNet@GRN$Target,LETTERS[1:16])
  checkEquals(unsupervisedNet@GRN$Target,LETTERS[1:16])
  checkEquals(unsupervisedNet@GRN$coact,c("a b","a b","c d","c d","a b","a b","c d","c d","e f","e f","g h","g h","e f","e f","g h","g h"))
  checkEquals(unsupervisedNet@GRN$corep,c("e f","e f","g h","g h","e f","e f","g h","g h","a b","a b","c d","c d","a b","a b","c d","c d" ))
}



test_masterRegulator <- function(){
  acts=apply(matrix(rep(letters[1:8],26),nrow=2),2,paste,collapse=" ")
  reps=apply(matrix(rep(letters[9:16],26),nrow=2),2,paste,collapse=" ")
  grn=data.frame("Target"= rep(LETTERS[1:26],4) ,"coact"=c(acts,reps),"corep"= c(reps,acts),"R2"=runif(16),stringsAsFactors=FALSE)
  co=coregnet(grn)
  
  x=rep(0.6776398,16)
  names(x)=letters[1:16]
  checkEqualsNumeric(masterRegulator(co,LETTERS[1:6]),x, tolerance=1.0e-3)
  
  
  exampleWeight = seq(-3,3,length.out = 26)
  names(exampleWeight) = LETTERS[1:26]
  x=rep(0.7623108,16)
  names(x)=letters[1:16]
  checkEqualsNumeric(masterRegulator(co,exampleWeight,"list"),x, tolerance=1.0e-3)
  
  examplePvalue =c( 10^-(2:9),seq(0,1,length.out = 20) )
  names(examplePvalue) = LETTERS[1:26]
  x=c(rep(0,8),rep(1.390255e-15,8))
  names(x)=letters[1:16]
  checkEqualsNumeric(masterRegulator(co,examplePvalue,"merg"),x,tolerance=1.0e-3)
  
  
}
