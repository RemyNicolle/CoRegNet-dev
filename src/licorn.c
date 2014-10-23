#include <R.h>
#include <stdio.h>	
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define TAILLE_MOT  100// maximum size of the name of a coreg
#define TAILLE_OCTET 8
#define TAILLE_LIGNE 10000 //TAILLE_LIGNE=taille maximum des lignes dans le fichier des corégulateurs
#define TAILLE_SUPPORT_BIT 500 //TAILLE_SUPPORT > (nombre d'expérience*2)/8
#define TAILLE_SUPPORT_INT 4000 //TAILLE_SUPPORT_INT = nombre d'expérience*2

int nb_permut;
int bool_seuil;
float seuil;
int nb_gene;
int nb_coregs;
int taille_ouvert;
int taille_ouvert_pm;
int taille_fermer;
int taille_fermer_pm;
int taille_final;
int nb_experiences;
int taille_max_ouvert;
int taille_max_ouvert_pm;
int taille_max_fermer;
int taille_max_fermer_pm;
int taille_max_final;



#include "prototypes.h"

Gene gene_vide;
Noeud noeud_vide;

#include "structure.c"
#include "score.c"
#include "chargement.c"
#include "pvalue.c"


char ** LICORN(char ** coreg, char ** genexp, char ** result,double * threshold,int * ouvertferme, int * taille_max_ligne, int * taille_max_reg, int * nexp,int *nombre_coreg , int *nombre_gene,int * maxGRNpergene){
    
    
    nb_gene=*nombre_gene;
    nb_coregs=*nombre_coreg;
    nb_experiences=*nexp;
    
    
    nb_permut=0;
    bool_seuil=1;
    seuil=*threshold;
    
    taille_ouvert=0;
    taille_ouvert_pm=0;
    taille_fermer=0;
    taille_fermer_pm=0;
    taille_final=0;

    
    taille_max_ouvert=((int)*ouvertferme);
    taille_max_ouvert_pm=taille_max_ouvert;
    taille_max_fermer=((int)*ouvertferme);
    taille_max_fermer_pm=taille_max_fermer;
    taille_max_final=((int)*maxGRNpergene);
    
    
    gene_empty();
    noeud_empty();
    
    Liste l;
    init_liste(&l);
    
    Liste_gene lg;
    init_liste_gene(&lg);
    
    srand(time(NULL));
 
    chargement_noeudsSansfichier(coreg, &l, taille_max_ligne, taille_max_reg, nexp,nombre_coreg);
    chargement_genesSansfichier(genexp,&lg ,taille_max_ligne, taille_max_reg,nombre_gene);
    remplissage_noeuds(&l);
    a_star(l,lg,result);
    return(result);
}

