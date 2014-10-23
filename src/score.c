



//////////////////////////////////////////////////////////////////////////////////////////////////
float score(unsigned char * coreg,unsigned char *gene,int support_gene_neg,int support_gene_pos){

	unsigned char resultat[TAILLE_SUPPORT_BIT];
	unsigned char masque_neg=170;
	unsigned char masque_pos=85;
	
	int i=0;
	int compte=0;

	int j=TAILLE_OCTET;
	unsigned char bit;
	unsigned char temp;
	float score_pos;
	float score_neg;

	if(!support_gene_pos) score_pos=0;
		else{
			//masque de +
			for(i=0;i<TAILLE_SUPPORT_BIT;i++){
				*(resultat+i)=0;
				*(resultat+i)=*(coreg+i) & *(gene+i) & masque_pos;
			}
	
			//comptage numérateur des +
			for(i=0;i<TAILLE_SUPPORT_BIT;i++){
				temp=*(resultat+i);
				for(j=0;j<TAILLE_OCTET;j++){
					bit=temp & 1;
					temp=temp >> 1;
					compte=compte+(int)bit;
				}
			}
			score_pos=(float)compte/support_gene_pos;
		}
	if(!support_gene_neg) score_neg=0;
		else{
			//masque de -
			for(i=0;i<TAILLE_SUPPORT_BIT;i++){
				*(resultat+i)=*(coreg+i) & *(gene+i) & masque_neg;
			}

			//comptage numérateur des -
			compte=0;
			for(i=0;i<TAILLE_SUPPORT_BIT;i++){
				temp=*(resultat+i);
				for(j=0;j<TAILLE_OCTET;j++){
					bit=temp & 1;
					temp=temp >> 1;
					compte=compte+(int)bit;
				}
			}
			score_neg=(float)compte/support_gene_neg;
		}
	return(score_pos+score_neg);
}
///////////////////////////////////////////////////////////////////////////////////////////////
int score_gene_pos(unsigned char *gene){
	int i;
	int j;
	unsigned char bit;
	unsigned char temp;
	int compte;
	unsigned char resultat[TAILLE_SUPPORT_BIT];
	unsigned char masque_pos=85;

	for(i=0;i<TAILLE_SUPPORT_BIT;i++){
		*(resultat+i)=0;
		*(resultat+i)=*(gene+i) & masque_pos;
	}

	compte=0;
	for(i=0;i<TAILLE_SUPPORT_BIT;i++){
		temp=*(resultat+i);
		for(j=0;j<TAILLE_OCTET;j++){
			bit=temp & 1;
			temp=temp >> 1;
			compte=compte+(int)bit;
		}
	}
return(compte);
}
///////////////////////////////////////////////////////////////////////////////////////////////
int score_gene_neg(unsigned char *gene){
	int i;
	int j;
	unsigned char bit;
	unsigned char temp;
	int compte;
	unsigned char resultat[TAILLE_SUPPORT_BIT];
	unsigned char masque_neg=170;


	for(i=0;i<TAILLE_SUPPORT_BIT;i++){        
		*(resultat+i)=0;
		*(resultat+i)=*(gene+i) & masque_neg;        
	}

	compte=0;
	for(i=0;i<TAILLE_SUPPORT_BIT;i++){
		temp=*(resultat+i);
		for(j=0;j<TAILLE_OCTET;j++){
			bit=temp & 1;
			temp=temp >> 1;
			compte=compte+(int)bit;
		}
	}

return(compte);
}
///////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
float score_pm(unsigned char * coreg,unsigned char *gene,int support_gene_neg,int support_gene_pos){

	unsigned char resultat[TAILLE_SUPPORT_BIT];
	//unsigned char masque_neg=170;
	unsigned char masque_pos=85;
	
	int i=0;
	int compte=0;

	int j=TAILLE_OCTET;
	unsigned char bit;
	unsigned char temp;
	float score_pos;
	float score_neg;

	if(!support_gene_pos) score_pos=0;
		else{
			//masque de coreg- inter gene+
			for(i=0;i<TAILLE_SUPPORT_BIT;i++){
				*(resultat+i)=0;
				*(resultat+i)=((*(gene+i) & masque_pos) << 1) & *(coreg+i);
			}

			//comptage numérateur des +
			for(i=0;i<TAILLE_SUPPORT_BIT;i++){
				temp=*(resultat+i);
				for(j=0;j<TAILLE_OCTET;j++){
					bit=temp & 1;
					temp=temp >> 1;
					compte=compte+(int)bit;
				}
			}
			score_pos=(float)compte/support_gene_pos;
		}
	if(!support_gene_neg) score_neg=0;
		else{
			//masque de coreg+ inter gene-
			for(i=0;i<TAILLE_SUPPORT_BIT;i++){
				*(resultat+i)=0;
				*(resultat+i)=((*(coreg+i) & masque_pos) << 1) & *(gene+i);
			}

			//comptage numérateur des -
			compte=0;
			for(i=0;i<TAILLE_SUPPORT_BIT;i++){
				temp=*(resultat+i);
				for(j=0;j<TAILLE_OCTET;j++){
					bit=temp & 1;
					temp=temp >> 1;
					compte=compte+(int)bit;
				}
			}
			score_neg=(float)compte/support_gene_neg;
		}
	return(score_pos+score_neg);
}
//////////////////////////////////////////////////////////////////////////////////////////////////
float score2(unsigned char * coreg,unsigned char *gene,int support_gene_neg,int support_gene_pos){
	
	unsigned char resultat[TAILLE_SUPPORT_BIT];
	unsigned char masque_neg=170;
	unsigned char masque_pos=85;
	
	int i=0;
	int compte=0;
	int compte_pos=0;
	int compte_neg=0;
	int support_gene=0;
	

	int j=TAILLE_OCTET;
	unsigned char bit;
	unsigned char temp;
	


	//masque de +
	for(i=0;i<TAILLE_SUPPORT_BIT;i++){
		*(resultat+i)=0;
		*(resultat+i)=*(coreg+i) & *(gene+i) & masque_pos;
	}

	//comptage numérateur des +
	for(i=0;i<TAILLE_SUPPORT_BIT;i++){
		temp=*(resultat+i);
		for(j=0;j<TAILLE_OCTET;j++){
			bit=temp & 1;
			temp=temp >> 1;
			compte_pos=compte_pos+(int)bit;
		}
	}


	//masque de -
	for(i=0;i<TAILLE_SUPPORT_BIT;i++){
		*(resultat+i)=*(coreg+i) & *(gene+i) & masque_neg;
	}

	//comptage numérateur des -
	//compte=0;
	for(i=0;i<TAILLE_SUPPORT_BIT;i++){
		temp=*(resultat+i);
		for(j=0;j<TAILLE_OCTET;j++){
			bit=temp & 1;
			temp=temp >> 1;
			compte_neg=compte_neg+(int)bit;
		}
	}


	compte=compte_pos+compte_neg;
	support_gene=support_gene_pos+support_gene_neg;

	if(!support_gene) return(0);
		else
			return((float)compte/support_gene);
	
}
///////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
float score2_pm(unsigned char * coreg,unsigned char *gene,int support_gene_neg,int support_gene_pos){
	
	unsigned char resultat[TAILLE_SUPPORT_BIT];
	//unsigned char masque_neg=170;
	unsigned char masque_pos=85;
	
	int i=0;
	int compte=0;
	int compte_coreg_pos_gene_neg=0;
	int compte_coreg_neg_gene_pos=0;
	int support_gene=0;
	

	int j=TAILLE_OCTET;
	unsigned char bit;
	unsigned char temp;
	//float score_neg;
	


	//masque de coreg+ inter gene-
	for(i=0;i<TAILLE_SUPPORT_BIT;i++){
		*(resultat+i)=0;
		/////
		*(resultat+i)=((*(coreg+i) & masque_pos) << 1) & *(gene+i);
		/////

	}

	//comptage numérateur des coreg+ inter gene-
	for(i=0;i<TAILLE_SUPPORT_BIT;i++){
		temp=*(resultat+i);
		for(j=0;j<TAILLE_OCTET;j++){
			bit=temp & 1;
			temp=temp >> 1;
			compte_coreg_pos_gene_neg=compte_coreg_pos_gene_neg+(int)bit;
		}
	}


	//masque de coreg- inter gene+
	for(i=0;i<TAILLE_SUPPORT_BIT;i++){
		*(resultat+i)=0;
		*(resultat+i)=((*(gene+i) & masque_pos) << 1) & *(coreg+i);
	}

	//comptage numérateur des -
	//compte=0;
	for(i=0;i<TAILLE_SUPPORT_BIT;i++){
		temp=*(resultat+i);
		for(j=0;j<TAILLE_OCTET;j++){
			bit=temp & 1;
			temp=temp >> 1;
			compte_coreg_neg_gene_pos=compte_coreg_neg_gene_pos+(int)bit;
		}
	}
	

	compte=compte_coreg_pos_gene_neg+compte_coreg_neg_gene_pos;
	support_gene=support_gene_pos+support_gene_neg;

	if(!support_gene) return(0);
		else
			return((float)compte/support_gene);
	
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
int mae(unsigned char *gene,unsigned char *profile_gene){

	int i;
	int j;
	unsigned char bit;
	unsigned char temp;
	int compte;
	unsigned char resultat[TAILLE_SUPPORT_BIT];
	//unsigned char masque_pos=85;

	for(i=0;i<TAILLE_SUPPORT_BIT;i++){
		*(resultat+i)=0;
		*(resultat+i)=*(gene+i) ^ *(profile_gene+i);

	}

	compte=0;
	for(i=0;i<TAILLE_SUPPORT_BIT;i++){
		temp=*(resultat+i);
		for(j=0;j<TAILLE_OCTET;j++){
			bit=temp & 1;
			temp=temp >> 1;
			compte=compte+(int)bit;
		}
	}

return(compte);
}
///////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
void profile_predit(unsigned char *act,unsigned char *inib,unsigned char *resultat){

	char tempAct;
	char tempInib;
	char tempRes=0;
	int i;
	int j;
	
	for(i=0;i<TAILLE_SUPPORT_BIT;i++){
		*(resultat+i)=0;
		j=TAILLE_OCTET-2;
		while(j>=0){
			tempAct=*(act+i) >> j;
			tempInib=*(inib+i) >> j;
			tempAct=tempAct & 3;
			tempInib=tempInib & 3;

			if (!tempInib) tempRes=tempAct;
			if (tempInib==1) tempRes=2;
			if (tempInib==2) {
				if (tempAct==2) tempRes=0;
				else tempRes=1;
			}

			tempRes=tempRes << j;
			*(resultat+i)=*(resultat+i) | tempRes;
			j=j-2;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////




