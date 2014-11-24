#include <R.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdio.h>	
#include <stdlib.h>



void combnLicorn(int * coactexp, int * ncoacts, int * corepexp , int * ncoreps, int * gexp, int * nsamples,  double * result , int * iact,int * irep){    
    int nsamp = (int)nsamples[0];
//    printf("%d \n",nsamp);
    int ncoact =(int)ncoacts[0];
    int ncorep = (int)ncoreps[0];
//        printf("%d \n",ncoact);
//        printf("%d \n",ncorep);
    int i =0;
    int j =0;
    int k=0;
    int  gexpred;
    int nonNormalizedMAE=0;
    for( i =0; i< ncoact;i++){
        for(  j =0; j< ncorep;j++){
            nonNormalizedMAE=0;
            for(  k =0; k <nsamp;k++){
                
                gexpred= coactexp[(i * nsamp)+k] - corepexp[(j * nsamp)+k];
                if (gexpred >1)  {
                    gexpred=1;
                } else if(gexpred < -1){
                    gexpred = -1;
                }
                nonNormalizedMAE += abs(gexp[k] -  gexpred);
//                printf("test %d; i %d j%d k%d: act %d , rep %d makes %d for %d\n",(i*ncoact)+j,i,j,k,coactexp[(i * nsamp)+k],corepexp[(j * nsamp)+k],gexpred,gexp[k]);
            }
            result[(i*ncorep)+j]= (nonNormalizedMAE *1.0)/(nsamp *1.0);
            iact[(i*ncorep)+j]=i+1;
            irep[(i*ncorep)+j]=j+1;
        }
    }
}



// internal function for loading in R
R_CMethodDef cMethods[] = {
    {"combnLicorn", (DL_FUNC) &combnLicorn, 9},
    {NULL, NULL, 0}
};

//R_CallMethodDef callMethods[] = {
//    {"combnLicorn", (DL_FUNC) &combnLicorn, 9},
//    {NULL, NULL, 0}
//};

// internal function for loading in R
void R_init_myLib(DllInfo *info)
{
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}

