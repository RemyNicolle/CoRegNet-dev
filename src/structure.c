


void fils_mauvais(Noeud *n);
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void Insert(Ouvert **sl,Noeud *n,float Val){
	static int taille_max=0;
        Ouvert *tmp = NULL;
	Ouvert *tmp2 = NULL;
        Ouvert *csl = *sl;
        Ouvert *elem = malloc(sizeof(Ouvert));
	int i=0;
	int nbnbmax=30;
        if(!elem) error("illegal elem value");
        elem->val = Val;
	elem->info=n;
	elem->suivant=NULL;
        while(csl && csl->val > Val && csl->info!=n){
		i++;
		tmp = csl;
        	csl = csl->suivant;
	}
	if(!csl || csl->info!=n){
       		elem->suivant = csl;
        	if(tmp)	{tmp->suivant = elem;
      		}else{
			*sl = elem;
		}
	}
	if(taille_max==nbnbmax){
		while(csl && csl->suivant){
		i++;
		tmp2 = csl;
		//if(!csl->suivant) 
		csl = csl->suivant;
		}
		if(!csl->suivant)
		free(csl);
		tmp2->suivant=NULL;
	}

}

///////////////////////////////////////////////////////////////////////////
//////////////////////////////////////

void PushFront(Liste_ouvert2 *l,Noeud *n,float val)
{

Ouvert2 *tmp =NULL;
Ouvert2 *csl = l->deb;

   Ouvert2 *nouv = malloc(sizeof(Ouvert2));
   if(!nouv) error("illegal nouv value");
   nouv->val = val;
   nouv->info=n;


while(csl && csl->val > val){
	tmp=csl;
	csl=csl->suivant;
}

nouv->suivant = csl;
if (csl) csl->precedent = nouv;
else l->fin=nouv;
if(tmp) {
	tmp->suivant = nouv;
	nouv->precedent = tmp;
	}else {
		if(l->deb) l->deb->precedent = nouv;
		else l->fin=nouv;
		nouv->precedent = NULL;
		l->deb = nouv;
	}

}
//////////////////////////////////////////////////////////////////

void Insert2(Liste_ouvert2 *l,Noeud *n,float val){
	//version avac ordre
	//version sans ordre
	if (taille_ouvert < taille_max_ouvert){

		PushFront(l,n,val);
		taille_ouvert++;
	}else if (val > l->fin->val && (!bool_seuil || val > seuil)){
		PushFront(l,n,val);
		PopBack_new(l);

	}
}
//////////////////////////////////////////////////////////////////

void Insert2_pm(Liste_ouvert2 *l,Noeud *n,float val){
	//version avec ordre
	//version sans ordre
	if (taille_ouvert_pm < taille_max_ouvert_pm){

		PushFront(l,n,val);
		taille_ouvert_pm++;
	}else if (val > l->fin->val && (!bool_seuil || val > seuil)){
		PushFront(l,n,val);
		PopBack_new(l);
		
	}
}
//////////////////////////////////////////////////////////////////

void Insert2_fermer_pm(Liste_ouvert2_fermer *l,Noeud *n,Gene *g,float val){
	if (taille_fermer_pm < taille_max_fermer_pm){
		PushFront_fermer(l,n,g,val);
		taille_fermer_pm++;
	}else if (val > l->fin->val){
		PushFront_fermer(l,n,g,val);
		PopBack_fermer_new(l);
	}
}
//////////////////////////////////////////////////////////////////
float PopBack_new(Liste_ouvert2 *l)
{
   float val;
   Ouvert2 *tmp = l->fin;
   if(!tmp) return -1;
   val = tmp->val;
   l->fin = tmp->precedent;
   if(l->fin) l->fin->suivant = NULL;
   else l->deb = NULL;
   free(tmp);
   return val;
}
//////////////////////////////////////////////////////////////////
float PopBack(Liste_ouvert2 *l,Noeud **n)
{
   float val;
   Ouvert2 *tmp = l->fin;
   if(!tmp) return -1;
	*n = malloc(sizeof(Noeud));
	*n = tmp->info;
   val = tmp->val;
   l->fin = tmp->precedent;
   if(l->fin) l->fin->suivant = NULL;
   else l->deb = NULL;
   free(tmp);
   return val;
}
//////////////////////////////////////////////////////////////////

float PopFront(Liste_ouvert2 *l,Noeud **n)
{
   float val;
   Ouvert2 *tmp = l->deb;
   if(!tmp) return -1;
   if(*n==NULL)	*n = malloc(sizeof(Noeud));
   *n = tmp->info;
   val = tmp->val;
   l->deb = tmp->suivant;
   if(l->deb)l->deb->precedent = NULL;
   else l->fin = NULL;
   free(tmp);
   taille_ouvert--;
   return val;
}
//////////////////////////////////////////////////////////////////

float PopFront_pm(Liste_ouvert2 *l,Noeud **n)
{
   float val;
   Ouvert2 *tmp = l->deb;
   if(!tmp) return -1;
   if(*n==NULL)	*n = malloc(sizeof(Noeud));
   *n = tmp->info;
   val = tmp->val;
   l->deb = tmp->suivant;
   if(l->deb)l->deb->precedent = NULL;
   else l->fin = NULL;
   free(tmp);
   taille_ouvert_pm--;
   return val;
}

//////////////////////////////////////////////////////////////////

void Insert2_fermer(Liste_ouvert2_fermer *l,Noeud *n,Gene *g,float val){
		if (taille_fermer < taille_max_fermer){
		PushFront_fermer(l,n,g,val);
		taille_fermer++;
	}else if (val > l->fin->val){
		PushFront_fermer(l,n,g,val);	//BUG i?i???!!!
		PopBack_fermer_new(l);
		///////////////////
	}
}

//////////////////////////////////////////////////////////////////

void Insert2_final(Liste_ouvert_final *l,Noeud *n1,Noeud *n2,float val){
	if (taille_final < taille_max_final){
		PushFront_final(l,n1,n2,val);
		taille_final++;
	}else if (val < l->fin->val){
		//final
		PushFront_final(l,n1,n2,val);
		PopBack_final(l);
	}
}

///////////////////////////////////////////////////////////////////////////
//////////////////////////////////////

void PushFront_final(Liste_ouvert_final *l,Noeud *n1,Noeud *n2,float val)
{

Ouvert_final *tmp =NULL;
Ouvert_final *csl = l->deb;

   Ouvert_final *nouv = malloc(sizeof(Ouvert_final));
   if(!nouv) error("illegal nouv value");
   nouv->val = val;
   nouv->info1=n1;
   nouv->info2=n2;

while(csl && csl->val < val){
	tmp=csl;
	csl=csl->suivant;
}
nouv->suivant = csl;
if (csl) csl->precedent = nouv;
else l->fin=nouv;
if(tmp) {
	tmp->suivant = nouv;
	nouv->precedent = tmp;
	}else {
		if(l->deb) l->deb->precedent = nouv;
		else l->fin=nouv;
		nouv->precedent = NULL;
		l->deb = nouv;
	}
}
///////////////////////////////////////////////////////////////////////////
//////////////////////////////////////

void PushFront_fermer(Liste_ouvert2_fermer *l,Noeud *n,Gene *g,float val)
{

Ouvert2_fermer *tmp =NULL;
Ouvert2_fermer *csl = l->deb;

   Ouvert2_fermer *nouv = malloc(sizeof(Ouvert2_fermer));
   if(!nouv) error("illegal nouv value");
   nouv->val = val;
   nouv->info=n;
   nouv->info_gene=g;

while(csl && csl->val > val){
	tmp=csl;
	csl=csl->suivant;
}

nouv->suivant = csl;
if (csl) csl->precedent = nouv;
else l->fin=nouv;
if(tmp) {
	tmp->suivant = nouv;
	nouv->precedent = tmp;
	}else {
		if(l->deb) l->deb->precedent = nouv;
		else l->fin=nouv;
		nouv->precedent = NULL;
		l->deb = nouv;
	}

}

//////////////////////////////////////////////////////////////////
float PopBack_final(Liste_ouvert_final *l)
{
   float val;
   Ouvert_final *tmp = l->fin;
   if(!tmp) return -1;
   val = tmp->val;
   l->fin = tmp->precedent;
   if(l->fin) l->fin->suivant = NULL;
   else l->deb = NULL;
   free(tmp);
   return val;
}


//////////////////////////////////////////////////////////////////
float PopBack_fermer(Liste_ouvert2_fermer *l,Noeud **n)
{
   float val;
   Ouvert2_fermer *tmp = l->fin;
   if(!tmp) return -1;
	*n = malloc(sizeof(Noeud));
	*n = tmp->info;
   val = tmp->val;
   l->fin = tmp->precedent;
   if(l->fin) l->fin->suivant = NULL;
   else l->deb = NULL;
   free(tmp);
   return val;
}
//////////////////////////////////////////////////////////////////
float PopBack_fermer_new(Liste_ouvert2_fermer *l)
{
   float val;
   Ouvert2_fermer *tmp = l->fin;
   if(!tmp) return -1;
   val = tmp->val;
   l->fin = tmp->precedent;
   if(l->fin) l->fin->suivant = NULL;
   else l->deb = NULL;
   free(tmp);
   return val;
}
//////////////////////////////////////////////////////////////////

float PopFront_fermer(Liste_ouvert2_fermer *l,Noeud **n)
{
   float val;
   Ouvert2_fermer *tmp = l->deb;
   if(!tmp) return -1;
	*n = malloc(sizeof(Noeud));
	*n = tmp->info;
   val = tmp->val;
   l->deb = tmp->suivant;
   if(l->deb)l->deb->precedent = NULL;
   else l->fin = NULL;
   free(tmp);
   taille_fermer--;
   return val;
}
//////////////////////////////////////////////////////////////////

float PopFront_fermer_pm(Liste_ouvert2_fermer *l,Noeud **n)
{
   float val;
   Ouvert2_fermer *tmp = l->deb;
   if(!tmp) return -1;
	*n = malloc(sizeof(Noeud));
	*n = tmp->info;
   val = tmp->val;
   l->deb = tmp->suivant;
   if(l->deb)l->deb->precedent = NULL;
   else l->fin = NULL;
   free(tmp);
   taille_fermer_pm--;
   return val;
}
///////////////////////////////////////////////////////////////////////////////
void Insert_fermer(Fermer **sl,Noeud *n,Gene *g,float Val){
        Fermer *tmp = NULL;
        Fermer *csl = *sl;
        Fermer *elem = malloc(sizeof(Fermer));
        if(!elem) error("illegal elem value");
        elem->val = Val;
	elem->info=n;
	elem->info_gene=g;
        while(csl && csl->val > Val){ // && csl->info_gene!=g){
		tmp = csl;
		csl = csl->suivant;
	}
        		elem->suivant = csl;
	        	if(tmp)	tmp->suivant = elem;
	        		else {*sl = elem;
				}

}

//////////////////////////////////////////////////////////////////////////////
float Pop(Ouvert **sl,Noeud **n){
	float Val;
	Ouvert *tmp;
	*n = malloc(sizeof(Noeud));
        if(!*sl) return -1;
	*n = (*sl)->info;
	tmp = (*sl)->suivant;
	Val = (*sl)->val;
	free(*sl);
	*sl = tmp;
	//compteur_ouvert--;
	return(Val);
}
//////////////////////////////////////////////////////////////////////////////
void fils_mauvais(Noeud *n){
	Fils *fils1 = n->premier_fils;
	while(fils1){
		fils1->info->marqueur=1;
		fils_mauvais(fils1->info);
		fils1=fils1->suivant;
	}
}
//////////////////////////////////////////////////////////////////////////////
void fils_mauvais_a_zero(Liste l){
	Noeud * n=l.deb;
	//Noeud *tmp=n;

	while(n){
		n->marqueur=0;
		n=n->suivant;
	}
}
//////////////////////////////////////////////////////////////////////////////
void passage_a_zero(Liste l){
	Noeud * n=l.deb;


	while(n){
		n->passage=0;
		n=n->suivant;
	}
}
///////////////////////////////////////////////////////////////////////////////
void init_a_star(Liste l,Liste_ouvert2 *lo,Gene *g){
	float resultat;
	while(l.deb && (l.deb)->nb_mot==1){
        
		resultat=score((l.deb)->support_bit,g->support_bit,g->score_neg,g->score_pos);
		if(resultat > seuil) Insert2(lo,l.deb,resultat);
		l.deb=(l.deb)->suivant;
	}
}
///////////////////////////////////////////////////////////////////////////////
void init_a_starpm(Liste l,Liste_ouvert2 *lo,Gene *g){
	float resultat;
	while(l.deb && (l.deb)->nb_mot==1){
		resultat=score_pm((l.deb)->support_bit,g->support_bit,g->score_neg,g->score_pos);
		if(resultat > seuil) Insert2_pm(lo,l.deb,resultat);
		l.deb=(l.deb)->suivant;
	}
}

void main_star(Liste l,Gene *lgdeb,Liste_ouvert_final *final_liste){

	Liste_ouvert2 lo2;
	init_liste_ouvert2(&lo2);

	Liste_ouvert2_fermer f;
	init_liste_ouvert2_fermer(&f);	

	float resultat;
	Noeud *n=NULL;
	Fils *fils1=NULL;
	//A*+-
	Liste_ouvert2 lo2_pm;
	init_liste_ouvert2(&lo2_pm);

	Liste_ouvert2_fermer f_pm;
	init_liste_ouvert2_fermer(&f_pm);

	float resultat_pm;
	Noeud *n_pm=NULL;
	Fils *fils1_pm=NULL;


		init_a_star(l,&lo2,lgdeb);
		init_a_starpm(l,&lo2_pm,lgdeb);

		while(lo2.deb){
			resultat=PopFront(&lo2,&n);
			Insert2_fermer(&f,n,lgdeb,resultat);
			
			fils1 = n->premier_fils;
			while(fils1){
                R_CheckUserInterrupt();
				if( fils1->info->marqueur==0 && fils1->info->passage==0 && (resultat=score(fils1->info->support_bit,(lgdeb)->support_bit,(lgdeb)->score_neg,(lgdeb)->score_pos)) > seuil ) {
					Insert2(&lo2,fils1->info,resultat);
					fils1->info->passage=1;
				}else{
					fils_mauvais(fils1->info);
				}
				fils1=fils1->suivant;
			}
		}

		fils_mauvais_a_zero(l);
		passage_a_zero(l);
		//
		while(lo2_pm.deb){
			resultat_pm=PopFront_pm(&lo2_pm,&n_pm);
			Insert2_fermer_pm(&f_pm,n_pm,lgdeb,resultat_pm);
			R_CheckUserInterrupt();
			fils1_pm = n_pm->premier_fils;
			while(fils1_pm){
				if( fils1_pm->info->marqueur==0 && fils1_pm->info->passage==0 && (resultat_pm=score_pm(fils1_pm->info->support_bit,(lgdeb)->support_bit,(lgdeb)->score_neg,(lgdeb)->score_pos)) > seuil ) {
					Insert2_pm(&lo2_pm,fils1_pm->info,resultat_pm);
					
					fils1_pm->info->passage=1;
				}else{
					fils_mauvais(fils1_pm->info);
				}
				fils1_pm=fils1_pm->suivant;
			}
		}

		fils_mauvais_a_zero(l);
		passage_a_zero(l);

		appariement(f,f_pm,lgdeb,final_liste);

		
		Clear_fermer2(&f);
		Clear_fermer2_pm(&f_pm);
		n=NULL;
		n_pm=NULL;
}

void a_star(Liste l,Liste_gene lg, char ** result){
	Liste_ouvert_final score_mae;
	Liste_ouvert_final score_mae_alea;
	int num_pvalue;
	int i;
	unsigned char gene_origine[TAILLE_SUPPORT_BIT];
    int grnNum=0;
	while(lg.deb){

		main_star(l,lg.deb,&score_mae);
		num_pvalue=0;

		for(i=0;i<TAILLE_SUPPORT_BIT;i++){
			gene_origine[i]=lg.deb->support_bit[i];
		}
		taille_final=0;

   
		for(i=0;i<nb_permut;i++){

			melange_gene(gene_origine,lg.deb->support_bit,nb_experiences);

			main_star(l,lg.deb,&score_mae_alea);///
			
			if (score_mae_alea.min < score_mae.min)	num_pvalue++;

			Clear_final(&score_mae_alea);
		}

		grnNum=merge_mae_result(lg.deb,&score_mae,num_pvalue, result,grnNum);

		Clear_final(&score_mae);

		lg.deb=(lg.deb)->suivant;
	}
}

//////////////////////////////////////////////////////////////////////////////

void gene_empty(void){
	int j;
	char nom[6]="EMPTY";

	gene_vide.id = -1;
	gene_vide.score_neg=0;
	gene_vide.score_pos=0;
	strcpy(gene_vide.nom,nom);
	for(j=0;j<TAILLE_SUPPORT_BIT;j++){
		gene_vide.support_bit[j]=0;
	}
	gene_vide.suivant=NULL;
}
void noeud_empty(void){
	int j=0;
	char nom[6]="EMPTY";

	noeud_vide.id=-1;
	noeud_vide.nb_mot=1;
	noeud_vide.marqueur=-1;
	noeud_vide.passage=-1;
	noeud_vide.score_local=-1;
	strcpy(noeud_vide.nom,nom);
	for(j=0;j<TAILLE_SUPPORT_BIT;j++){
		noeud_vide.support_bit[j]=0;
	}
	noeud_vide.suivant=NULL;
	noeud_vide.premier_fils=NULL;
}

///////////////////////////////////////////////////////////////////////////////
void appariement(Liste_ouvert2_fermer f,Liste_ouvert2_fermer f_pm,Gene *g,Liste_ouvert_final *final_liste){
	Ouvert2_fermer *element_fermer = f.deb;
	Ouvert2_fermer *element_fermer_pm = f_pm.deb;
	unsigned char gene_predit[TAILLE_SUPPORT_BIT];
	int final=0;
	init_liste_ouvert_final(final_liste);	

	
	while(element_fermer){
		element_fermer_pm=f_pm.deb;
		while(element_fermer_pm){
	
			if(element_fermer->info->id != element_fermer_pm->info->id){
				profile_predit(element_fermer->info->support_bit,element_fermer_pm->info->support_bit,gene_predit);

				final=mae(g->support_bit,gene_predit);
				
			
				Insert2_final(final_liste,element_fermer->info,element_fermer_pm->info,final);
			}
			element_fermer_pm = element_fermer_pm->suivant;
			
		}
		element_fermer = element_fermer->suivant;
	}


//	/**********************
		element_fermer_pm=f_pm.deb;
		while(element_fermer_pm){
			if(noeud_vide.id != element_fermer_pm->info->id){

				
				profile_predit(noeud_vide.support_bit,element_fermer_pm->info->support_bit,gene_predit);

				final=mae(g->support_bit,gene_predit);
				
			
				Insert2_final(final_liste,&noeud_vide,element_fermer_pm->info,final);
			}
			element_fermer_pm = element_fermer_pm->suivant;
			
		}


// ************************************
		element_fermer = f.deb;
		while(element_fermer){
			if(noeud_vide.id != element_fermer->info->id){

				

				profile_predit(element_fermer->info->support_bit,noeud_vide.support_bit,gene_predit);

				final=mae(g->support_bit,gene_predit);
				
			
				Insert2_final(final_liste,element_fermer->info,&noeud_vide,final);
			}
			element_fermer = element_fermer->suivant;
			
		}


// ************************

///////////////////////
	final_liste->nb_co_activateur=Nb_coregulateurs_fermer2(f)+1;
	final_liste->nb_co_inhibiteur=Nb_coregulateurs_fermer2(f_pm)+1;
	final_liste->min=mae_min(final_liste);

}
///////////////////////////////////////////////////////////////////////////////
int pas_vide(Ouvert *sl){
	if(sl){
		return 1;
	}else{
		return 0;
	}
}
///////////////////////////////////////////////////////////////////////////////
void View_fermer(Fermer *sl){
       while(sl)
          {

             sl = sl->suivant;
          }
}

///////////////////////////////////////////////////////////////////////////////
void View(Ouvert *sl){
       while(sl)
          {
             sl = sl->suivant;
          }
}
/////////////////////////////////////////////////////////////////////////////
int Nb_coregulateurs_fermer2(Liste_ouvert2_fermer f){
   Ouvert2_fermer *element_fermer = f.deb;
   int i=0;
   while(element_fermer){
	i++;
	element_fermer = element_fermer->suivant;
   }
   return(i);
}

/////////////////////////////////////////////////////////////////////////////
float mae_min(Liste_ouvert_final *f){

	float mae=-1;
	float mae_min=-1;

   Ouvert_final *element_final = f->deb;

	if(element_final){
		mae_min=(float)element_final->val/nb_experiences;
		element_final = element_final->suivant;
	}

   while(element_final){

	mae=(float)element_final->val/nb_experiences;
	if(mae<mae_min){mae_min=mae;}

	element_final = element_final->suivant;
   }
return(mae_min);
}
/////////////////////////////////////////////////////////////////////////////


int merge_mae_result(Gene *g,Liste_ouvert_final *f,int num_pvalue,char ** result ,int grnNum){


    Ouvert_final *element_final = f->deb;


    
    while(element_final){
        int len =strlen(g->nom) + strlen(element_final->info1->nom) + strlen(element_final->info1->nom)+46;

        char * ligne= malloc(  (len) *sizeof(char));
        char sep[1] =";";
        strcpy(ligne,g->nom);
        strncat(ligne,sep,1);
        strncat(ligne,element_final->info1->nom,strlen(element_final->info1->nom));
        strncat(ligne,sep,1);
        strncat(ligne,element_final->info2->nom,strlen(element_final->info2->nom));
        strncat(ligne,sep,1);
        if(nb_permut){
            char pc[20] = "pval";
             strncat(ligne,pc,strlen(pc));
        }else{
            char str[20];
            sprintf(str, "%f",  ((float)element_final->val/nb_experiences));
            strncat(ligne,str,strlen(str));
        }
       
        
        result[ grnNum++ ]=ligne;

        element_final = element_final->suivant;
    }

    return((int)grnNum);
}



/////////////////////////////////////////////////////////////////////////////

void View_final_affichage(Gene *g,Liste_ouvert2_fermer fe,Liste_ouvert2_fermer f_pm,Liste_ouvert_final *f){

   Ouvert_final *element_final = f->deb;
   while(element_final){
		
	element_final = element_final->suivant;
   }
}

/////////////////////////////////////////////////////////////////////////////

void View_final(Liste_ouvert_final f){
   Ouvert_final *element_final = f.deb;
   while(element_final){
	element_final = element_final->suivant;
   }
}
/////////////////////////////////////////////////////////////////////////////

void Clear_final(Liste_ouvert_final *f){
   Ouvert_final *tmp;
   Ouvert_final *pelem = f->deb;
	
   while(pelem){
	tmp = pelem;
	pelem = pelem->suivant;
	free(tmp);
	//taille_final--;
   }
   taille_final=0;
   f->fin = NULL;
   f->deb = NULL;
}

/////////////////////////////////////////////////////////////////////////////

void View_fermer2(Liste_ouvert2_fermer f){
   Ouvert2_fermer *element_fermer = f.deb;
   while(element_fermer){
	element_fermer = element_fermer->suivant;
   }
}
/////////////////////////////////////////////////////////////////////////////

void Clear_fermer2(Liste_ouvert2_fermer *f){
   Ouvert2_fermer *tmp;
   Ouvert2_fermer *pelem = f->deb;
	
   while(pelem){
	tmp = pelem;
	pelem = pelem->suivant;
	free(tmp);
	taille_fermer--;
   }
   f->fin = NULL;
   f->deb = NULL;
}
/////////////////////////////////////////////////////////////////////////////

void Clear_fermer2_pm(Liste_ouvert2_fermer *f){
   Ouvert2_fermer *tmp;
   Ouvert2_fermer *pelem = f->deb;
	
   while(pelem){
	tmp = pelem;
	pelem = pelem->suivant;
	free(tmp);
	taille_fermer_pm--;
   }
   f->fin = NULL;
   f->deb = NULL;
}

//////////////////////////////////////////////////////////////////////////////
void Clear_fermer(Fermer **p){
        Fermer *tmp;
        while(*p)
          {
             tmp = (*p)->suivant;
             free(*p);
             *p = tmp;
          }
}


/////////////////////////////////////////////////////////////////////////////

void nouveau_noeud2(int i,Noeud **graphe){
	Noeud *n = malloc(sizeof(Noeud));
	n->id = i;
	n->suivant=*graphe;
	n->premier_fils=NULL;
	*graphe = n;


}
/////////////////////////////////////////////////////////////////////////////
void nouveau_noeud(int i,int nb_mot,char *nom,unsigned char *support_bit,Liste *l){
	int j;
	Noeud *n = malloc(sizeof(Noeud));
	n->id = i;
	n->nb_mot = nb_mot;
	n->marqueur=0;
	n->passage=0;
	strcpy(n->nom,nom);
	for(j=0;j<TAILLE_SUPPORT_BIT;j++){
		n->support_bit[j]=support_bit[j];
	}
	n->suivant=NULL;
	if(!l->deb){	l->deb=n;
	}else{
			l->fin->suivant=n;
	}
	n->premier_fils=NULL;
	l->fin= n ;


}
/////////////////////////////////////////////////////////////////////////////
void nouveau_gene(int i,char *nom,unsigned char *support_bit,int support_neg,int support_pos,Liste_gene *lg){
	Gene *g = malloc(sizeof(Gene));
	g->id = i;
	g->score_neg=support_neg;
	g->score_pos=support_pos;

	strcpy(g->nom,nom);

	int j;

	for(j=0;j<TAILLE_SUPPORT_BIT;j++){
		g->support_bit[j]=support_bit[j];
	}

	g->suivant=NULL;
	if(!lg->deb){	lg->deb=g;
	}else{
			lg->fin->suivant=g;
	}
	lg->fin = g;


}
/////////////////////////////////////////////////////////////////////////////
void init_liste(Liste *l){
	l->deb=NULL;
	l->fin=NULL;
}
/////////////////////////////////////////////////////////////////////////////
void init_liste_gene(Liste_gene *lg){
	lg->deb=NULL;
	lg->fin=NULL;
}
/////////////////////////////////////////////////////////////////////////////
void init_liste_ouvert2(Liste_ouvert2 *l){
	l->deb=NULL;
	l->fin=NULL;
}
/////////////////////////////////////////////////////////////////////////////
void init_liste_ouvert2_fermer(Liste_ouvert2_fermer *l){
	l->deb=NULL;
	l->fin=NULL;
}
/////////////////////////////////////////////////////////////////////////////
void init_liste_ouvert_final(Liste_ouvert_final *l){
	l->deb=NULL;
	l->fin=NULL;
}
/////////////////////////////////////////////////////////////////////////////
void parcours_liste(Liste l){
	while(l.deb){
		parcours_liste_fils(l.deb);
		l.deb=(l.deb)->suivant;
	}
}
/////////////////////////////////////////////////////////////////////////////
void parcours_liste_gene(Liste_gene lg){
	while(lg.deb){
		lg.deb=(lg.deb)->suivant;
	}
}
/////////////////////////////////////////////////////////////////////////////
void parcours_liste_ouvert2(Liste_ouvert2 lo){
	while(lo.deb){
		lo.deb=(lo.deb)->suivant;
	}
}
/////////////////////////////////////////////////////////////////////////////
void parcours_liste_fils(Noeud *n){

	Fils *fils1;
	fils1 = n->premier_fils;
	while(fils1){
		fils1=fils1->suivant;
	}
}
/////////////////////////////////////////////////////////////////////////////
void ajouter_fils(Noeud *n, Noeud *fils){

	Fils *new_fils = malloc(sizeof(Fils));
	new_fils->info=fils;
	new_fils->suivant=n->premier_fils;
	n->premier_fils=new_fils;

}
/////////////////////////////////////////////////////////////////////////////

