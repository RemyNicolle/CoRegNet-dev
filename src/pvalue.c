
/////////////////////////////////////////////////////////////////////////////
int pif(int t){
	int r;
	r = (int)(rand() / (double)RAND_MAX * (t - 1));
return(r);
}
/////////////////////////////////////////////////////////////////////////////

int iemeRand(int * tab_bool,int r){

	int j=0;
	int nb=-1;
	int val=-1;

	while(nb<r){
		if(*(tab_bool+j)){val++;}
			else{	val++;
				nb++;
			}
		j++;
	}
return(val);
}

/////////////////////////////////////////////////////////////////////////////
int lire_ieme_bit(unsigned char * tab_gene,int i){
	int resultat;
	int zone;
	int rang;
	int decalage;
	int masque=3;
	zone = i / (TAILLE_OCTET/2);
	rang = i % (TAILLE_OCTET/2);
	decalage = 2*((TAILLE_OCTET/2 - 1) - rang);
	resultat = ( *(tab_gene+zone) >> decalage) & masque;
return(resultat);
}
/////////////////////////////////////////////////////////////////////////////
void ecrire_ieme_bit(unsigned char * tab_gene,int i,int val){
	int zone;
	int rang;
	int decalage;
	zone = i / (TAILLE_OCTET/2);
	rang = i % (TAILLE_OCTET/2);
	decalage = 2*((TAILLE_OCTET/2 - 1) - rang);
	*(tab_gene+zone) = (val << decalage) | *(tab_gene+zone);



}

/////////////////////////////////////////////////////////////////////////////
void melange_gene(unsigned char * tab_gene,unsigned char * tab_res, int nb_experiences){
	int i;
	int r;
	int j;
	int b;
	int dec=0;
	int tab_bool[nb_experiences];
	for(i=0;i<nb_experiences;i++){
		tab_bool[i]=0;
	}

	for(i=0;i<TAILLE_SUPPORT_BIT;i++){
		tab_res[i]=0;
	}

	for(i=nb_experiences+1;i>=0;i--){
		r=pif(i);
		j=iemeRand(tab_bool,r);
		b=lire_ieme_bit(tab_gene,j);
		tab_bool[j]=1;
		ecrire_ieme_bit(tab_res,dec,b);
		dec++;
	}
}
//////////////////////////////////////////////////////////////////////////////////////
