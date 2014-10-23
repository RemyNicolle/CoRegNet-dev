/////////////////////////////////////////////////////////////////////////////
typedef struct tagNoeud {
	int id;
	int nb_mot;
	int marqueur;
	int passage;
	int score_local;
	char nom[TAILLE_MOT];
	unsigned char support_bit[TAILLE_SUPPORT_BIT];
	struct tagNoeud *suivant;
	struct tagFils *premier_fils;
}Noeud;

typedef struct tagGene {
	int id;
	char nom[TAILLE_MOT];
	unsigned char support_bit[TAILLE_SUPPORT_BIT];
	int score_neg;
	int score_pos;
	struct tagGene *suivant;
}Gene;


typedef struct tagFils {
	int cout;
	Noeud *info;
	struct tagFils *suivant;
}Fils;
typedef struct {
	Noeud * deb;
	Noeud * fin;
}Liste;

typedef struct tagOuvert {
	float val;
	Noeud *info;
	struct tagOuvert *suivant;
}Ouvert;


typedef struct tagOuvert2 {
	float val;
	Noeud *info;
	struct tagOuvert2 *precedent;
	struct tagOuvert2 *suivant;
}Ouvert2;

typedef struct tagOuvert_final {
	float val;
	Noeud *info1;
	Noeud *info2;
	//Gene *info_gene besoin d'un gene???? a priori non
	struct tagOuvert_final *precedent;
	struct tagOuvert_final *suivant;
}Ouvert_final;
//2->_final

typedef struct {
	Ouvert_final *deb;
	Ouvert_final *fin;
	int nb_co_activateur;
	int nb_co_inhibiteur;
	float min;
}Liste_ouvert_final;


typedef struct tagOuvert2_fermer {
	float val;
	Noeud *info;
	Gene *info_gene;
	struct tagOuvert2_fermer *precedent;
	struct tagOuvert2_fermer *suivant;
}Ouvert2_fermer;

typedef struct {
	Ouvert2 *deb;
	Ouvert2 *fin;
}Liste_ouvert2;

typedef struct {
	Ouvert2_fermer *deb;
	Ouvert2_fermer *fin;
}Liste_ouvert2_fermer;


typedef struct tagFermer {
	float val;
	Noeud *info;
	Gene *info_gene;
	struct tagFermer *suivant;
}Fermer;
typedef struct {
	Gene * deb;
	Gene * fin;
}Liste_gene;
/////////////////////////////////////////////////////////////////////////////

//faires les fonctions FREE MEMORY!!!!

/////////////////////////////////////////////////////////////////////////////
float PopBack_new(Liste_ouvert2 *l);
float PopBack_fermer_new(Liste_ouvert2_fermer *l);
void nouveau_noeud2(int i,Noeud **graphe);
void nouveau_noeud(int i,int nb_mot,char *nom,unsigned char *support_bit,Liste *l);
void nouveau_gene(int i,char *nom,unsigned char *support_bit,int support_neg,int support_pos,Liste_gene *lg);
void init_liste(Liste *l);
void init_liste_gene(Liste_gene *lg);
void init_liste_ouvert2(Liste_ouvert2 *l);
void parcours_liste(Liste l);
void parcours_liste_gene(Liste_gene lg);
void parcours_liste_ouvert2(Liste_ouvert2 lo);
void parcours_liste_fils(Noeud *n);
void ajouter_fils(Noeud *n, Noeud *fils);
void Insert(Ouvert **sl, Noeud *n, float Val);
void Insert2(Liste_ouvert2 *l,Noeud *n,float val);
float PopBack(Liste_ouvert2 *l,Noeud **n);
void Insert2_fermer(Liste_ouvert2_fermer *l,Noeud *n,Gene *g,float val);
void PushFront_fermer(Liste_ouvert2_fermer *l,Noeud *n,Gene *g,float val);
float PopBack_fermer(Liste_ouvert2_fermer *l,Noeud **n);
void init_liste_ouvert2_fermer(Liste_ouvert2_fermer *l);
void View_fermer2(Liste_ouvert2_fermer f);
void Clear_fermer2(Liste_ouvert2_fermer *f);
void Clear_fermer2_pm(Liste_ouvert2_fermer *f);
void View(Ouvert *sl);
int pas_vide(Ouvert *sl);
void View_fermer(Fermer *sl);
void appariement(Liste_ouvert2_fermer f,Liste_ouvert2_fermer f_pm,Gene *g,Liste_ouvert_final *final_liste);
void profile_predit(unsigned char *act,unsigned char *inib,unsigned char *resultat);
int mae(unsigned char *gene,unsigned char *profile_gene);
void init_liste_ouvert_final(Liste_ouvert_final *l);
void Insert2_final(Liste_ouvert_final *l,Noeud *n1,Noeud *n2,float val);
float PopBack_final(Liste_ouvert_final *l);
void PushFront_final(Liste_ouvert_final *l,Noeud *n1,Noeud *n2,float val);
void View_final(Liste_ouvert_final f);
void View_final_affichage(Gene *g,Liste_ouvert2_fermer fe,Liste_ouvert2_fermer f_pm,Liste_ouvert_final *f);
void Clear_final(Liste_ouvert_final *f);
int Nb_coregulateurs_fermer2(Liste_ouvert2_fermer f);
void View_mae_affichage(Gene *g,Liste_ouvert_final *f,int num_pvalue);
float mae_min(Liste_ouvert_final *f);
///////////////////////////////////////////////////////////////////////////////////////////////
float score(unsigned char * coreg,unsigned char *gene,int support_gene_neg,int support_gene_pos);
float score_pm(unsigned char * coreg,unsigned char *gene,int support_gene_neg,int support_gene_pos);
float score2(unsigned char * coreg,unsigned char *gene,int support_gene_neg,int support_gene_pos);
float score2_pm(unsigned char * coreg,unsigned char *gene,int support_gene_neg,int support_gene_pos);
int score_gene_pos(unsigned char *gene);
int score_gene_neg(unsigned char *gene);
///////////////////////////////////////////////////////////////////////////////////////////////


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

int pif(int t);
int iemeRand(int * tab_bool,int r);
int lire_ieme_bit(unsigned char * tab_gene, int i);
void ecrire_ieme_bit(unsigned char * tab_gene, int i,int val);
void melange_gene(unsigned char * tab_gene,unsigned char * tab_res, int nb_experiences);
///////////////////////////////////////////////////////////////////////////////
void gene_empty();
void noeud_empty();

int merge_mae_result(Gene *g,Liste_ouvert_final *f,int num_pvalue,char ** result ,int grnNum);
