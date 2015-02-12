CoRegNet
========
CoRegNet is an R/Bioconductor package for gene regulatory network inference and analysis.

It uses a a Shiny/Cytoscape.js interface for the analysis of networks.


Install

From R :
```
source("http://bioconductor.org/biocLite.R")
biocLite("CoRegNet")
```

Quick user guide
----------------

1. reconsruct a large-scale regulatory network from a gene expression matrix EXP
  
  ```
	  GRN = hLICORN(EXP)
  ```
  
2. Infer transcription factor activity
  
  ```
  	influence = regulatorInfluence(GRN,EXP)
  ```
  
3. Retrieve inferred co-coregulators
  
  ```
  	coregs= coregulators(GRN,EXP)
  ```
  
4. Analyze the network of cooperative regulators using an interactive display
 ```
  display(GRN,EXP,influence)
```

