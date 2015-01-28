

.DisplayEnv <- new.env()
assign("SELECTEDCOREG",  c(), envir = .DisplayEnv)



display <- function(coregnetwork,expressionData=NULL,TFA = NULL,alterationData=NULL ,clinicalData=NULL,TFnotes=NULL,
                    allTFplot =.heatplot ,oneTFplot=.tfPlot ){
  
  # testing input
  if(is.null(TFA) & is.null(expressionData)){
    stop("TFA or Expression is needed.")
  }
  if(is.null(TFA) & !is.null(expressionData)){
    TFA= regulatorInfluence(coregnetwork,expressionData)  
  }
  
  x=coregulators(coregnetwork,alpha = 10^-2)[1:4]
  if(is.null(nrow(x) )| nrow(x) ==0){
      stop("No co-regulators in the provided network. If it was inferred with the hLICORN function, try a lower minCoregSupport.")
  }
  initThreshCoreg=ceiling(0.01*max(x$nGRN))
  maxThreshCoreg = max(x$nGRN)
  x=x[which(x$nGRN >=initThreshCoreg),1:2]
  if(length( intersect(  unique(unlist( x  )),rownames(TFA)))<5){
      stop("Too few co-regulators have a measure of influence making the interactive visualization both useless and difficult")
  }
  
  assign("SELECTEDCOREG",   intersect(  unique(unlist( x  )),rownames(TFA)), envir = .DisplayEnv)
  TFA=TFA[get("SELECTEDCOREG",envir=.DisplayEnv),]
  
  ##############
  if(!is.null(alterationData)){
    alterationData = alterationData[intersect(rownames(alterationData),rownames(TFA))
                                    ,intersect(colnames(alterationData),colnames(TFA))]
    if(nrow(alterationData) <6 | ncol(alterationData) <10  ){
      warnings("Very few samples or regulators in alteration data which will not be used")
      alterationData=NULL
    }
  }
  ##############
  
  # will contain a factor vector of the sample (clinical) classification
  sampleClassif=c()
  ##############
  # This section processes the input clinical data
  # Clinical data is needed for two things. The color on the side of the heatmap and the color of nodes when selecting one of the info
  # Testing the clinical data. Can be a simple factor (or factorizable ) vector or a list of samples each entry corresponding to a class    
  if(!is.null(clinicalData)){  
    if(is.factor(clinicalData)& sum(names(clinicalData) %in% colnames(TFA)) > 0.5*ncol(TFA)){
      if(sum(names(clinicalData) %in% colnames(TFA))<length(clinicalData)){
        givenSample =names(clinicalData)
        missingSamples = setdiff(colnames(TFA),givenSample)
        sampleClassif = c(clinicalData,rep(NA,length(missingSamples)))
        names(sampleClassif) = c(givenSample,missingSamples)
        sampleClassif=sampleClassif[colnames(TFA)]
        sampleClassif=factor(sampleClassif)
      }else{
        sampleClassif=clinicalData[colnames(TFA)]
        sampleClassif=factor(sampleClassif)
      }
      clinicalData=split(names(sampleClassif),sampleClassif)                
      clinicalData = c(clinicalData,list("targ" =colnames(TFA) ))
      
    }else if(is.list(clinicalData)){
      factorized =rep(names(clinicalData),sapply(clinicalData,length))
      if(length(factorized) > 0 & length(factorized) <= ncol(TFA) ){
        sampleClassif = factorized
        names(sampleClassif)= unlist(clinicalData)
        sampleClassif=sampleClassif[colnames(TFA)]
        sampleClassif=factor(sampleClassif)
      }else{
        sampleClassif=NULL
        warning("Samples have more than one class.")
      }
      clinicalData = c(clinicalData,list("targ" =colnames(TFA) ))
    }else{
      clinicalData = list("targ" =colnames(TFA) )
      sampleClassif=NULL
    }
  }else{
    clinicalData = list("targ" =colnames(TFA) )
    sampleClassif=NULL
  }  
  clinicalInfo = names(clinicalData)
  names(clinicalInfo) = clinicalInfo
  names(clinicalInfo)[length(clinicalInfo)]="All"
  ##############
  # in the end, sampleClassif should be a factor, named with sample names, each level being a class of samples
  # and clinicalData should
  
  shinyApp= list(
    ui=navbarPage("CoRegNet",
                  tabPanel("Co-regulation",


                           sidebarLayout(position="left",
                                         sidebarPanel(
                                           ##### including javascripts and css
                                           

                                           includeScript(system.file('www/js/jquery-1.11.1.min.js', package = 'CoRegNet')),
                                           
                                           includeCSS(system.file('www/css/style.css', package = 'CoRegNet')),
                                           
                                           includeScript(system.file('www/js/cytoscape.min.js', package = 'CoRegNet')),
                                           includeScript(system.file('www/js/cy.js', package = 'CoRegNet')),
                                           

                                           includeCSS(system.file('www/js/jquery-slick/css/dcslick.css', package = 'CoRegNet')),
                 
                                           includeScript(system.file('www/js/jquery-slick/js/jquery.slick.2.1.js', package = 'CoRegNet')),
                                           
                                           selectInput("info", "Information:",clinicalInfo,selected="targ"),
                                           numericInput("mingrn", "Minimum number of shared GRN", initThreshCoreg, min = initThreshCoreg, max = maxThreshCoreg,step = 1),
                                           textInput("seltf", "Search", value = ""),textOutput("findInfo")
                                         ),
                                         mainPanel( 
                                           actionButton("layout", "Layout"),
                                           HTML( "<div id=\"cy\" class=\"cynetwork\"></div>"),
                                           HTML("<div id=\"legendNetwork_slick\" style=\"float: left; display: none;\">"), plotOutput("legendNetwork"), HTML("</div>"),
                                           HTML("<div id=\"legend_slick\" style=\"float: left; display: none;\">"), plotOutput("legend"), HTML("</div>")
                                         )
                           ),
                           plotOutput("plot")

                  ),tabPanel("Snapshot",
                             
                             HTML( '<img src="" id="snap" alt="Cytoscape needs two " width=100% height:100%>'),
                             textOutput(paste("A snapshot of the network.", 
                                              "Cytoscape always shows the previous network.", 
                                              "To get the latest snapshot, you need to go back to the coregnet tab and back here again." ,
                                              "Sorry for this inconvenience."))
                             ),id="tabset"),
    
    #SERVER
    server= function(input, output,session) { 

      
      .clinic <- clinicalData[["targ"]]
      
      # observe if clicked on the layout button to reperform layout
      observe({
        x=input$layout        
        output$cy <- redoLayout(function(){})
      })
      
      # observe going to snapshot
      observe({
        x=input$tabset             
        if(x== "Snapshot"){
          output$cy <-printNetwork(function(){})             
        }
      })
      
      # observe search field
      observe({
        x=input$seltf  
        
        if(x !=""){
        if(x %in% get("SELECTEDCOREG",envir=.DisplayEnv)){
          #print(x)
          output$cy <-selectNode(function(){},x)
          output$findInfo <- renderText("")
        }else{
          matched=grep(paste("^",x,sep=""),get("SELECTEDCOREG",envir=.DisplayEnv),ignore.case=TRUE)
          if(length(matched)==0){
            output$findInfo <- renderText("Nothing found.")
          }else if( length(matched)==1){
            x=get("SELECTEDCOREG",envir=.DisplayEnv)[matched]
            #print(x)
            output$cy <-selectNode(function(){},x)
            output$findInfo <- renderText(x)
          }else{
            Y=get("SELECTEDCOREG",envir=.DisplayEnv)[matched]            
            Y = sort(Y)
            output$findInfo <- renderText(paste(paste(Y[1:min(c(length(Y),3))],collapse=", "),"..."))
          }
        }
        }
      })
      
      
      # observe parameters of coregulators which changes the net
       observe({ 
         thresh=input$mingrn
         output$cy <- reactiveAdjacencyMatrix(function() {},coregnetwork,.clinic,TFA,alterationData,input$mingrn)  
        })

            
      output$legend <- renderPlot({legendPlot(sampleClassif,alterationData)   })
      
      # observe modification of the selected clinical data
      observe({  

        
        if(is.null(input$info)){
          .clinic<-clinicalData[["targ"]]
          output$legendNetwork <- renderPlot({networkLegendPlot(coregnetwork)})
        }else if(input$info == "targ"| input$info==""){
          .clinic<-clinicalData[["targ"]]
          
          output$legendNetwork <- renderPlot({networkLegendPlot(coregnetwork)})
          output$cy <- updateData(function() {return(NULL)},coregnetwork,get("SELECTEDCOREG",envir=.DisplayEnv), .clinic,TFA,alterationData)
        }else{
          .clinic<-clinicalData[[input$info]]  
      
          output$cy <- updateData(function() {return(NULL)},coregnetwork,get("SELECTEDCOREG",envir=.DisplayEnv), .clinic,TFA,alterationData)
          output$legendNetwork <- renderPlot({networkLegendPlot(coregnetwork,FCcol=c("red","blue"))})
        }                  
      },priority=1)
      
      #observe the cynetwork to see if we select or click some nodes    
      observe({
        clicked=input$cy
        output$plot <- renderPlot({
            
          if(is.null(clicked)[1] ){
                allTFplot(TFA[intersect(rownames(TFA),get("SELECTEDCOREG",envir=.DisplayEnv)),],clinical=sampleClassif,TFnotes=TFnotes)
          }else if(clicked[1] == "NULL" & length(clicked)<2){
              
            allTFplot(TFA[intersect(rownames(TFA),get("SELECTEDCOREG",envir=.DisplayEnv)),],clinical=sampleClassif,TFnotes=TFnotes)
          }else if(length(clicked)==1){
              
            oneTFplot(coreg=coregnetwork,TF=clicked,TFA=TFA,
                      expressionData=expressionData,alteration=alterationData,clinical=sampleClassif)
          }else{
                                                   allTFplot(TFA[intersect(intersect(rownames(TFA),clicked),get("SELECTEDCOREG",envir=.DisplayEnv)),],clinical=sampleClassif,TFnotes=TFnotes)
          }
        },height="auto")
      })
    }) 
  shiny::runApp(shinyApp)
}




.distfun=function(x) as.dist((1 - cor(t(x)))/2 )
.hclustfun = function(d) hclust(d,method="ward.D2")

.heatplot = function(subexpr,clinical=NULL,TFnotes=NULL){
  if(!is.null(clinical)){
    if(length(clinical) != ncol(subexpr)){
      stop("Wrong size of clinical data")
    }
    if(!is.factor(clinical)){clinical=factor(clinical) }       
    if(require(RColorBrewer)){
      
      mypalette<-RColorBrewer::brewer.pal(max(c(length(levels(clinical))),3),"Set3")
    }else{mypalette<-rainbow(length(levels(clinical)))}        
    clinical = mypalette[as.numeric(clinical)]
  }
  subexpr=t(subexpr)
  distf = function(x){as.dist((1-cor(t(x))))}
  hclustf = function(d){hclust(d,method="ward.D2")}
  
  bluetored=colorRampPalette(c("blue","white", "red"))
    if(require(gplots)){
    if(is.null(clinical) & is.null(TFnotes)){
      heatmap.2(subexpr, col=bluetored(75), scale="none",distfun=.distfun,hclustfun=.hclustfun,
                key=TRUE, symkey=TRUE, density.info="none", trace="none", cexRow=0.5,
                breaks= quantile(as.vector(subexpr),probs=seq(0,1,by=1/75)))
    }else if(is.null( clinical)){
      heatmap.2(subexpr, col=bluetored(75), scale="none",distfun=.distfun,hclustfun=.hclustfun,
                key=TRUE, symkey=TRUE, density.info="none", trace="none", cexRow=0.5,breaks= quantile(as.vector(subexpr),
                                                                                                   probs=seq(0,1,by=1/75)),
                ColSideColors=                    )
    }else if(is.null(TFnotes)){
      heatmap.2(subexpr, col=bluetored(75), scale="none",distfun=.distfun,hclustfun=.hclustfun,
                key=TRUE, symkey=TRUE, density.info="none", trace="none", cexRow=0.5,
                breaks= quantile(as.vector(subexpr),probs=seq(0,1,by=1/75))
                ,RowSideColors=clinical)
    }else{
      heatmap.2(subexpr, col=bluetored(75), scale="none",distfun=.distfun,hclustfun=.hclustfun,
                key=TRUE, symkey=TRUE, density.info="none", trace="none", cexRow=0.5,
                breaks= quantile(as.vector(subexpr),probs=seq(0,1,by=1/75))
                ,ColSideColors=clinical ,RowSideColors=clinical)
    }
    }   
}
  



.tfPlot = function(coreg,TF,TFA,expressionData,alteration=NULL, clinical=NULL){
  BR=colorRampPalette(c("blue","white", "red"))
  GR=colorRampPalette(c("green","black", "red"))  
  commonSamples=colnames(TFA)
  

  
  if(!is.null(clinical)){
    commonSamples=intersect(intersect(names(clinical),colnames(TFA)),colnames(expressionData))  
  }  
  if(!is.null(alteration)){
    if(TF %in% rownames(alteration)){    
      commonSamples=intersect(commonSamples,colnames(alteration))
    }
  }  
  if(!is.null(clinical)){
    clinical = clinical[commonSamples]  
    if(!is.factor(clinical)){clinical=factor(clinical) }       
    clinical = as.numeric(clinical)
    names(clinical)=commonSamples
    if(require(RColorBrewer)){    
        mypalette<-RColorBrewer::brewer.pal(max(c(max(clinical),3)),"Set3")
    }else{mypalette<-rainbow(length(levels(clinical)))}   
    
  }
   
  par(mfrow=c(4+(!is.null(alteration))+(!is.null(clinical)),1),mar=c(1,1,1,1),oma=c(0,0,4,0))
  
  tfA = TFA[TF,commonSamples]
  orderedSample=names(tfA)[order(tfA)]
  tfA=tfA[orderedSample]  
  tfexpr = as.numeric(expressionData[TF,orderedSample])  
  if(!is.null(clinical)){
  image(as.matrix(as.data.frame(as.numeric(clinical[orderedSample]))),axes=FALSE,main="Clinical",
        breaks = seq(0.5,max(clinical)+0.5,by=1),col=mypalette)
  mtext(TF, outer = TRUE, cex = 2)
  }

  if(!is.null(alteration)){
    if(TF %in% rownames(alteration)){ 
      tfCNA = as.numeric(alteration[TF,orderedSample])
      image(as.matrix(data.frame("a"=tfCNA)),main=paste(TF,"Copy Number"),breaks=seq(-2.5,2.5,by=1),
            col=c("#2b83ba", "#abdda4","#ffffbf","#fdae61","#d7191c"),axes=FALSE)    
    }
  }  
  image(as.matrix(data.frame("a"=tfexpr)),col=GR( 10  ),main=paste(TF,"Expression"),breaks=seq(min(tfexpr),max(tfexpr),by=(max(tfexpr)-min(tfexpr))/10),axes=FALSE)
  image(as.matrix(data.frame("a"=tfA)),col=BR(20),main=paste(TF,"Influence"),breaks=seq(min(tfA),max(tfA),by=(max(tfA)-min(tfA))/20),axes=FALSE)
  activated = scale(t(as.matrix(expressionData[targets(coreg,TF,"activating"),orderedSample])),scale=FALSE)
  repressed =scale(t(as.matrix(expressionData[targets(coreg,TF,"repressing"),orderedSample])),scale=FALSE)
  image((activated),col=GR(10),axes=FALSE,breaks=quantile(as.vector(as.numeric(unlist(activated))),probs=seq(0,1,by=1/10)),main="Activated genes expression")
  image((repressed),col=GR(10),axes=FALSE,breaks=quantile(as.vector(as.numeric(unlist(repressed))),probs=seq(0,1,by=1/10)),main="Repressed genes expression")
}

legendPlot = function(clinical=NULL,alteration=NULL)
{

  np =(sum(!is.null(clinical))>=1)  + (sum(!is.null(alteration))>=1) +2
  #par(mfrow=c(np,1))
  par(mfrow=c(4,1))
  BR=colorRampPalette(c("blue","white", "red"))
  n=75
  image(as.matrix(data.frame(1:n)),main="Transcriptional activity (influence)",breaks=seq(0.5,(n+0.5),by=1),
        col=BR(n),axes=FALSE)   
  axis(1, at = c(0,1), labels =c("loss","high"), tick=FALSE)
  
  
  GR=colorRampPalette(c("green","black", "red")) 
  image(as.matrix(data.frame(1:n)),main="Expression",breaks=seq(0.5,(n+0.5),by=1),
        col=GR(n),axes=FALSE)   
  axis(1, at = 0:1, labels =c("low","high"), tick=FALSE)
  
  
  
  if(!is.null(alteration)){
  image(as.matrix(data.frame(-2:2)),main="Copy Number",breaks=seq(-2.5,2.5,by=1),
        col=c("#2b83ba", "#abdda4","#ffffbf","#fdae61","#d7191c"),axes=FALSE)   
  axis(1, at = (0:4)/4, labels =c("deletion","loss","diploid","gain","amplification"), tick=FALSE)
  
  }
  
  if(!is.null(clinical)){
   n=length(levels(clinical))
    if(!is.factor(clinical)){clinical=factor(clinical) }       
    if(require(RColorBrewer)){
        mypalette<-RColorBrewer::brewer.pal(max(c(n,3)),"Set3")
    }else{mypalette<-rainbow(n)}           
    clinicalLevelsCol = mypalette[1:n]    
   
    image(as.matrix(data.frame(1:n)),main="Sample classification",breaks=seq(0.5,(n+0.5),by=1),
          col=clinicalLevelsCol,axes=FALSE)   
    axis(1, at = (0:(n-1))/(n-1), labels =levels(clinical), tick=FALSE)    
  }

}

networkLegendPlot=function(coreg,defCol="purple",sizeInfo="Number of targets",FCcol=NULL)
{
  narro=length(coreg@evidences)+1
  par(mar=c(0,0,0,0))
  plot(NULL,xlim=c(1,10),ylim=c(1,7+2*(narro)),bty="n",main="",xlab="",ylab="",yaxt="n",xaxt="n")
  
  symbols(1.5,5+2*narro,circles=0.2,fg="purple",bg="purple",add=TRUE,inches=FALSE)
  symbols(1.5,3+2*narro,circles=0.4,fg="purple",bg="purple",add=TRUE,inches=FALSE)
  if(is.null(FCcol)){
  
  }else{
  text(c(1,2.5,2.5), c((6+2*narro),5+2*narro,3+2*narro), c("Differential influence significance","low","high"),adj=0)  
    symbols(5.5,5+2*narro,circles=0.3,fg=FCcol[1],bg=FCcol[1],add=TRUE,inches=FALSE)
    symbols(5.5,3+2*narro,circles=0.3,fg=FCcol[2],bg=FCcol[2],add=TRUE,inches=FALSE)
    
    text(c(6,7.5,7.5), c((6+2*narro),5+2*narro,3+2*narro), c("Fold change","high","low"),adj=0)
    
  }
  
  if(length(coreg@evidences) > 0){
    if(require(RColorBrewer)){
        mypalette<-RColorBrewer::brewer.pal(length(coreg@evidences),"Set1")
    }else{ mypalette<-rainbow(length(coreg@evidences))}
    
    segments(x0=1,y0=2*narro+1 ,x1=4,y1=2*narro+1 ,col="#bdbdbd",lwd=4)
    text( 5,2*narro+1, "Co-regulation",adj=0)
    
    
    for( i in 1:length(coreg@evidences)){
      if(coreg@evidenceDescription$evidenceType[i] =="regulatory"){            
        arrows(x0=1,y0=(2*(narro)+1)-(2*i),x1=4,y1=(2*(narro)+1)-(2*i),col=mypalette[i],lwd=4)        
      }else{
        segments(x0=1,y0=(2*(narro)+1)-(2*i),x1=4,y1=(2*(narro)+1)-(2*i),col=mypalette[i],lwd=4)
      }
      text( 5,(2*(narro)+1)-(2*i),names(coreg@evidences)[i],adj=0)  
    }
  }
}

testCyto <- function(func){
  function(){x=func()
             
             return(list(names=list("PPARG"=list("ampli"=50,delet=20,diplo=10,gain=10,loss=10)),type="test"))
  }
}


selectNode <- function(func,tfname){
  function(){
    x=func()
    #print(tfname)
    return(list(tosel=tfname,type="select"))
  }
}


updateData <-function(func,coreg,regulators, clinicalDataSamples,DATA,CNV,absoluteNodeColor=TRUE){
  function(){
    x=func()
    regs=   regulators(coreg)[unique(unlist(regulators))]
    nodes=lapply(names(regs),function(r){
      x=list()
      x$color=0
      x$size=regs[r]/(max(regs))
      return(x)
    })
    names(nodes) = names(regs)

    if(length(clinicalDataSamples) < ncol(DATA)& length(clinicalDataSamples) > 3 & !absoluteNodeColor){
      ysamp=clinicalDataSamples
      xsamp=setdiff(colnames(DATA),ysamp)      
      forsize=-log10(apply(DATA,1,function(x){ return(t.test(x[ysamp],x[xsamp])$p.value)}) )
      forsize[which(forsize < 2)]=0
      forsize=forsize/max(forsize)
      forcolor = as.numeric(apply(DATA,1,function(x){ return(mean(x[ysamp])-mean(x[xsamp]))}))
      names(forcolor) = rownames(DATA)
      names(forsize) = rownames(DATA)
      nodes=lapply(names(nodes),function(nn){
        nx=nodes[[nn]]
        nx$size=unname(forsize[nn])
        nx$color=unname(forcolor[nn])
        return(nx)
      })
      names(nodes) = names(regs)
    }else if( length(clinicalDataSamples) >2){
      forsize=abs(as.numeric(apply(DATA[,clinicalDataSamples],1,mean)))
      forsize=forsize/max(forsize)
      forcolor = as.numeric(apply(DATA[,clinicalDataSamples],1,mean))
      names(forcolor) = rownames(DATA)
      names(forsize) = rownames(DATA)      
      nodes=lapply(names(nodes),function(nn){
        nx=nodes[[nn]]
        nx$size=unname(forsize[nn])
        nx$color=unname(forcolor[nn])
        return(nx)
      })
      names(nodes) = names(regs)
    }else if( length(clinicalDataSamples)==1 ){
      forsize=abs(as.numeric(DATA[,clinicalDataSamples]))
      forsize=forsize/max(forsize)
      forcolor = as.numeric(DATA[,clinicalDataSamples])
      names(forcolor) = rownames(DATA)
      names(forsize) = rownames(DATA)     
      nodes=lapply(names(nodes),function(nn){
        nx=nodes[[nn]]
        nx$size=unname(forsize[nn])
        nx$color=unname(forcolor[nn])
        return(nx)
      })
      names(nodes) = names(regs)
      
    }
    if(!is.null(CNV)){
      if((length(clinicalDataSamples) <= ncol(DATA)& length(clinicalDataSamples )> 3)){
        samplestouse = intersect(colnames(CNV),clinicalDataSamples)  
      }else{
        samplestouse = intersect(colnames(CNV),colnames(DATA))
      }      

  regnames=names(nodes)
      nodes=lapply(names(nodes),function(nn){
        if(nn %in% rownames(CNV)){
        nonNA = sum(!is.na(as.numeric(CNV[nn,samplestouse])))      
        nodecnv=list(
          delet=signif(((sum(as.numeric(CNV[nn,samplestouse]) == -2,na.rm=TRUE)/nonNA)*100),3)  ,
          loss=signif(((sum(as.numeric(CNV[nn,samplestouse]) == -1,na.rm=TRUE)/nonNA)*100),3)  ,
          diplo=signif(((sum(as.numeric(CNV[nn,samplestouse]) == 0,na.rm=TRUE)/nonNA)*100),3)  ,
          gain=signif(((sum(as.numeric(CNV[nn,samplestouse]) == 1,na.rm=TRUE)/nonNA)*100) ,3) ,
          ampli=signif(((sum(as.numeric(CNV[nn,samplestouse]) == 2,na.rm=TRUE)/nonNA)*100),3)  ,
          unknown=signif(((sum(is.na(as.numeric(CNV[nn,samplestouse])))/nonNA)*100)  ,3)
        )
        if(sum(unlist(nodecnv))>100){
          abov=sum(unlist(nodecnv))-100
          nodecnv[[which.max(unlist(nodecnv))]]=nodecnv[[which.max(unlist(nodecnv))]] - (abov)
        }
        return(c(nodes[[nn]],nodecnv))
        }else{
          return(c(nodes[[nn]],list(delet=NA,loss=NA,diplo=NA,gain=NA,ampli=NA,unknowl=100)))                      
        }
      })
  names(nodes)=regnames

    }

    return(list(names=nodes,type="update"))
  }
}

reactiveAdjacencyMatrix <- function(func,coreg,clinicalData,DATA,CNV,minGRN){
  function(){
      
      
      OKcoregulators=coregulators(coreg,alpha = 10^-2)
      OKcoregulators=OKcoregulators[which(OKcoregulators$nGRN >=minGRN),]
      
      assign("SELECTEDCOREG",  intersect(  unique(unlist( OKcoregulators[,1:2]  )),rownames(DATA)), envir = .DisplayEnv)
      
      regs=   regulators(coreg)[get("SELECTEDCOREG",envir=.DisplayEnv)]
      
      interactions = data.frame( OKcoregulators,"arrow"=rep("none",nrow(OKcoregulators)),color=rep("#bdbdbd",nrow(OKcoregulators)),stringsAsFactors=FALSE)
      interactions$nGRN =interactions$nGRN/max(interactions$nGRN)
      i=which(colnames(interactions) %in% c("fisherTest","adjustedPvalue"))
      interactions=interactions[,- i]
      
      if(length(coreg@evidences)>0){
          if(require(RColorBrewer)){
              mypalette<-RColorBrewer::brewer.pal(length(coreg@evidences),"Set1")
          }else{
              mypalette<-rainbow(length(coreg@evidences))
          }
    
          names(mypalette)=names(coreg@evidences)
          additionalInts=lapply(names(coreg@evidences),function(evname){
              ev=coreg@evidences[[evname]]
              #print(evname)
              addint = ev[which(ev[,1] %in% names(regs) & ev[,2] %in%names(regs) & ev[,1] != ev[,2]),1:2]
              addint= addint[which(sapply(1:nrow(addint),function(i){
                  return(length(intersect(which(OKcoregulators[,1] == addint[i,1] | OKcoregulators[,2] == addint[i,1] ),
                                          which(OKcoregulators[,1] == addint[i,2] | OKcoregulators[,2] == addint[i,2] )))>0)
              })),]
              addint=data.frame(addint,Support=rep(1,nrow(addint)),nGRN=rep(NA,nrow(addint)),
                                arrow=rep(ifelse(coreg@evidenceDescription[evname,"evidenceType"]=="regulatory","triangle","none"),nrow(addint)),
                                color=rep(mypalette[evname],nrow(addint)))
              colnames(addint)[1:2 ] = c("Reg1","Reg2")
              return(addint)
          })
          interactions=as.data.frame(rbind(interactions,do.call(rbind,additionalInts))      ,stringsAsFactors=FALSE)
      }

     assign("SELECTEDCOREG",  unique(c(interactions[,1],interactions[,2])), envir = .DisplayEnv)
     
    
    regs=   regulators(coreg)[unique(unlist(interactions[,1:2]))]
    nodes=lapply(names(regs),function(r){
      x=list()
      x$id=r
      x$name=r
      x$color=0
      x$size=unname(regs[r]/(max(regs)))
      return(x)
    })
    names(nodes) = names(regs)        
    interactions=data.frame("id"=paste("e",1:nrow(interactions),sep=""),interactions,stringsAsFactor=FALSE)      
    if(!is.null(CNV)){
      if((length(clinicalData) < ncol(DATA)& length(clinicalData )> 3)){
        samplestouse = intersect(colnames(CNV),clinicalData)  
      }else{
        samplestouse = intersect(colnames(CNV),colnames(DATA))
      }

      regnames=names(nodes)
      nodes=lapply(names(nodes),function(nn){
        if(nn %in% rownames(CNV)){
        nonNA = sum(!is.na(as.numeric(CNV[nn,samplestouse])))        
        return(c(nodes[[nn]],list(
          delet=signif(((sum(as.numeric(CNV[nn,samplestouse]) == -2,na.rm=TRUE)/nonNA)*100),3)  ,
          loss=signif(((sum(as.numeric(CNV[nn,samplestouse]) == -1,na.rm=TRUE)/nonNA)*100),3)  ,
          diplo=signif(((sum(as.numeric(CNV[nn,samplestouse]) == 0,na.rm=TRUE)/nonNA)*100),3)  ,
          gain=signif(((sum(as.numeric(CNV[nn,samplestouse]) == 1,na.rm=TRUE)/nonNA)*100) ,3) ,
          ampli=signif(((sum(as.numeric(CNV[nn,samplestouse]) == 2,na.rm=TRUE)/nonNA)*100),3)  ,
          unknown=signif(((sum(is.na(as.numeric(CNV[nn,samplestouse])))/nonNA)*100)  ,3)
        )))        
        }else{
          return(c(nodes[[nn]],list(delet=NA,loss=NA,diplo=NA,gain=NA,ampli=NA,unknown=100)))                      
        } 
      })
      edges=list(nodes=unname(nodes),links=interactions,type="pienetwork")
    }else{
      nodes=unname(nodes ) # necessary so that JS thinks it's an array and not an object with lots of attributes, each being a TF...
      edges=list(nodes=nodes,links=interactions,type="network")
      #print("return net")
    }
    return(edges)
  }
}

redoLayout <-function(func){function(){return(list(type="layout"))}}
printNetwork <-function(func){function(){  return(list(type="png"))  }}

