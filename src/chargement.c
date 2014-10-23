

/////////////////////////////////////////////////////////////////////////////
void parcours_ligne_coreg(char * ligne,char * nom, char * support_char,int *nb_mot);
void parcours_ligne_gene(char * ligne,char * nom, char * support_char);
void lire_support(char *tab_char,unsigned char *tab_int);
void lire_support_gene(char *tab_char,unsigned char *tab_int);
void support_int_vers_bit(unsigned char *tab,unsigned char *var);
int mot_appartient(char * mot1,char * mot2);
int est_fils(char *mot1, char *mot2);
void chargement_noeuds(FILE *fichier,Liste *l);
void chargement_genes(FILE *fichier,Liste_gene *lg);
void remplissage_noeuds(Liste *l);
/////////////////////////////////////////////////////////////////////////////
void parcours_ligne_coreg(char * ligne,char * nom, char * support_char,int *nb_mot){

	int i=0;
	int j=0;
	*nb_mot=1;


	while (	(i<TAILLE_LIGNE) && (*(ligne+i)!='\0') && (*(ligne+i)!=10) && (!((*(ligne+i)==' ')&&(*(ligne+i+1)-'0'>=0)&&(*(ligne+i+1)-'0'<=9)))){
		*(nom+i)=*(ligne+i);

		if(*(nom+i)==' ') *nb_mot=(*nb_mot)+1;
		i++;
	}
    
	*(nom+i)='\0';
	i++;
	while (i<TAILLE_LIGNE && *(ligne+i)!='\0'){
		*(support_char+j)=*(ligne+i);
		i++;
		j++;
	}
	*(support_char+j)='\0';
		
}
//////////////////////////////////////////////////////////////////////////////
void parcours_ligne_gene(char * ligne,char * nom, char * support_char){
	int i=0;
	int j=0;

	while (	(i<TAILLE_LIGNE) && (*(ligne+i)!='\0') && (*(ligne+i)!=10) && (!(*(ligne+i)==' '))){
		*(nom+i)=*(ligne+i);
		i++;
	}
	*(nom+i)='\0';
	i++;
	while (i<TAILLE_LIGNE && *(ligne+i)!='\0'){
		*(support_char+j)=*(ligne+i);
		i++;
		j++;
	}
	*(support_char+j)='\0';

}
/////////////////////////////////////////////////////////////////////////////
void lire_support(char *tab_char,unsigned char *tab_int){

	int i=0;
	int var=0;
	//initialisation à zéro de tout le support
	for(i=0;i<TAILLE_SUPPORT_INT;i++){
		*(tab_int+i)=0;
	}
	i=0;
	
	//verifier les limites!!!
	//tab_int[0] et tab_int[1] inutilisé , pas grave
	while (i<TAILLE_SUPPORT_INT && *(tab_char+i)!='\0'){
		//is chiffre
		if ((*(tab_char+i)-'0'>=0)&&(*(tab_char+i)-'0'<=9)){
			if ((*(tab_char+i+1)-'0'>=0)&&(*(tab_char+i+1)-'0'<=9)){
				var=(10*(*(tab_char+i)-'0'))+(*(tab_char+i+1)-'0');
				i++;
			}else{
			var=*(tab_char+i)-'0';
			}
			//var ok regarder + ou -
			if (*(tab_char+i+1)=='-'){
				*(tab_int+(2*var))=1;
			}else{
				*(tab_int+((2*var)+1))=1;
			}
	
		}
		i++;
	}

}
/////////////////////////////////////////////////////////////////////////////
void support_int_vers_bit(unsigned char *tab,unsigned char *var){
	int i=0;
	int j=0;
	
	for(i=0;i<TAILLE_SUPPORT_BIT;i++){
		*(var+i)=0;
	}
	
	i=0;
	while(i<TAILLE_SUPPORT_INT){
		*(var+j)=*(var+j) << 1;
		*(var+j)=*(var+j)+*(tab+i);


		i++;
		if(!(i%TAILLE_OCTET))	j++;
	}
}
/////////////////////////////////////////////////////////////////////////////
int mot_appartient(char * mot1,char * mot2){
	int i=0;
	int j=0;
	int saute=0;
	while(j<TAILLE_MOT && i<TAILLE_MOT && *(mot1+i)!='\0' && *(mot2+j)!='\0'){


		if(*(mot2+j)==' '){saute=0;i=0;
			}else{ if(saute==1){i=0;
				}else{ if(*(mot1+i)!=*(mot2+j)){saute=1;i=0;

					}else{if(*(mot1+i+1)=='\0' && !(*(mot2+j+1)==' ' || *(mot2+j+1)=='\0')){

						saute=1;i=0;
						}else{

						i++;
						
						}
					}
				}
			}
		j++;
	}

	return((*(mot1+i)=='\0'));
}

/////////////////////////////////////////////////////////////////////////////
int est_fils(char *mot1, char *mot2){
	char sous_mot[TAILLE_MOT];
	int i=0;
	int j=0;
	int appartient=0;
	while(i<TAILLE_MOT){

		if(*(mot1+i)==' ' || *(mot1+i)=='\0'){
			*(sous_mot+j)='\0';
			appartient=mot_appartient(sous_mot,mot2);
			if(!appartient) return(0);
			if(appartient && *(mot1+i)=='\0'){
				return(1);
			}else{
				j=-1;
			}

		}else{
			*(sous_mot+j)=*(mot1+i);
		}
		i++;
		j++;
	}
return(0);
}
///////////////////////////////////////////////////////////////////////
void chargement_noeudsSansfichier(char ** coreg, Liste *l, int * taille_max_ligne, int *   taille_max_reg, int*  nexp,int*  nombre_coreg){

    

	char nom[TAILLE_MOT];
	char support_char[TAILLE_LIGNE];
	unsigned char support_int[TAILLE_SUPPORT_INT];
	unsigned char support_bit[TAILLE_SUPPORT_BIT];

    int *nb_mot=malloc(sizeof(int));
    int n = 0; // nombre de lignes lues


    while(n < *nombre_coreg){
        parcours_ligne_coreg((*(coreg+n)),nom,support_char,nb_mot);
		lire_support(support_char,support_int);
		support_int_vers_bit(support_int,support_bit);
		//remplir un noeud
		nouveau_noeud(n+1,*nb_mot,nom,support_bit,l);
		n++;
    }


}


void chargement_noeuds(FILE *fichier,Liste *l){
	char ligne[TAILLE_LIGNE];
	char nom[TAILLE_MOT];
	char support_char[TAILLE_LIGNE];
	unsigned char support_int[TAILLE_SUPPORT_INT];
	unsigned char support_bit[TAILLE_SUPPORT_BIT];

	int *nb_mot=malloc(sizeof(int));
	int n = 0; // nombre de lignes lues 

	while (fgets( ligne, TAILLE_LIGNE, fichier)) {

		parcours_ligne_coreg(ligne,nom,support_char,nb_mot);
		lire_support(support_char,support_int);
		support_int_vers_bit(support_int,support_bit);
		//remplir un noeud
		nouveau_noeud(n+1,*nb_mot,nom,support_bit,l);
		n++;
	}

fclose(fichier);
}
//////////////////////////////////////////////////////////////////////

void chargement_genesSansfichier(char ** geneexp ,Liste_gene *lg,int * taille_max_ligne, int *taille_max_gene, int * nombre_gene){
  char nom[TAILLE_MOT];
	char support_char[TAILLE_LIGNE];
	unsigned char support_int[TAILLE_SUPPORT_INT];
	unsigned char support_bit[TAILLE_SUPPORT_BIT]; 
	int support_neg;
	int support_pos;
	int n=0;
	
	while ( n < * nombre_gene ){
		parcours_ligne_gene(    (*(geneexp+n)),nom,support_char);
		lire_support_gene(support_char,support_int);
		support_int_vers_bit(support_int,support_bit);
		support_neg=score_gene_neg(support_bit);
		support_pos=score_gene_pos(support_bit);
    			nouveau_gene(n,nom,support_bit,support_neg,support_pos,lg);
        

        
		n++;
	}

}




void chargement_genes(FILE *fichier,Liste_gene *lg){
	char ligne[TAILLE_LIGNE];
	char nom[TAILLE_MOT];
	char support_char[TAILLE_LIGNE];
	unsigned char support_int[TAILLE_SUPPORT_INT];
	unsigned char support_bit[TAILLE_SUPPORT_BIT];
	int support_neg;
	int support_pos;
	int n=0;
	//int i=0;



	while (fgets( ligne, TAILLE_LIGNE, fichier)){// && n<10) {

		parcours_ligne_gene(ligne,nom,support_char);
		lire_support_gene(support_char,support_int);
		support_int_vers_bit(support_int,support_bit);


			//remplir un noeud
			//**********************************//
		support_neg=score_gene_neg(support_bit);
		support_pos=score_gene_pos(support_bit);

		nouveau_gene(n,nom,support_bit,support_neg,support_pos,lg);


		n++;
	}

fclose(fichier);
}
//////////////////////////////////////////////////////////////////////
void remplissage_noeuds(Liste *l){
	Noeud *noeud_temp_pere=(l->deb);
	Noeud *noeud_temp_fils;

	while(noeud_temp_pere){
		noeud_temp_fils = noeud_temp_pere->suivant;
		while(noeud_temp_fils && ((noeud_temp_pere->nb_mot+1 == noeud_temp_fils->nb_mot)||(noeud_temp_pere->nb_mot == noeud_temp_fils->nb_mot) )){
			if(noeud_temp_pere->nb_mot+1 == noeud_temp_fils->nb_mot){
				if(est_fils(noeud_temp_pere->nom,noeud_temp_fils->nom)){
					ajouter_fils(noeud_temp_pere,noeud_temp_fils);


				}
			}
		noeud_temp_fils = noeud_temp_fils->suivant;
		}

		noeud_temp_pere = noeud_temp_pere->suivant;
		//i++;

	}


}
//////////////////////////////////////////////////////////////////////
void lire_support_gene(char *tab_char,unsigned char *tab_int){

	int i=0;
	int compteur=2;//le même décalage que pour les corégulateurs (2 bits)
	//initialisation à zéro de tout le support
	for(i=0;i<TAILLE_SUPPORT_INT;i++){
		*(tab_int+i)=0;
	}



	i=0;
	
	//verifier les limites!!!
	//tab_int[0] et tab_int[1] inutilisé , pas grave
	while (i<TAILLE_LIGNE && *(tab_char+i)!='\0'){
		//is chiffre


		if(*(tab_char+i)=='-'){
			*(tab_int+compteur)=1;
			compteur=compteur+2;
			i=i+3;
		}
		if(*(tab_char+i)=='1'){
			//*(tab_int+compteur)=0;
			*(tab_int+compteur+1)=1;
			compteur=compteur+2;
			i=i+2;
		}
		if(*(tab_char+i)=='0'){
			compteur=compteur+2;
			i=i+2;
		}
		i++;
	}

}
/////////////////////////////////////////////////////////////////
