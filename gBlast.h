//***************GBLAST.H************************

//

//Program Entry

//

//

//Author: Rong She

//Date: April 2007

//***********************************************





#ifndef GBLAST_H /* GBLAST_H */

#define GBLAST_H


#include <time.h>




#include "graph.h"





//parameters and their default values

float HSP_SK_ALPHA = 0.5f;

float GENE_MS_BETA = 0.5f;

float MAX_OVERLAP_PER = 0.5f; //no longer used

int	MAX_GAP_BTWN_HSP = 100000;

int TOP_RANK_NUM = INT_MAX;

float MIN_GENE_COVER_PER = 0.0f;

float MIN_GENE_SCORE;

bool USE_MIN_GENE_SCORE = false;



int MIN_EXON_SIZE = 2;

int MIN_INTRON_SIZE = 30; //no longer used?

int MAX_EXON_SIZE = INT_MAX;

int MAX_INTRON_SIZE = MAX_GAP_BTWN_HSP; //INT_MAX;

int MAX_EXON_HSP_DIFF = 4;

int MAX_INTRON_HSP_DIFF = 4;

//int MIN_INTRON_LEN = 15; //used by version 2.3 (=15/3, in terms of amino acids)
int START_CODON_CHECK_LEN = 15;
int MAX_NUM_SPLICE_SITES = 20; //used by version 2.3
int SPLICE_LEN = 5; //used by version 2.3, in terms of amino acids

//int MIN_INTERNAL_EXON_LEN = 20;


string PHASE1_VERSION = "1.1";

string PHASE1_OUTPUT = "c"; //'a': no overlap in base pair; 'b': allow everything; 'c': no overlap in hsp

string PHASE2_VERSION = "2.3"; //2.0: old; 2.1: new; 2.2: newest

int SPLICE_SEGMENT_VERSION = 1; //default: 1, optional: 2

float TREE_CLS_THRESHOLD = 0;

int TREE_DATA_SHIFT=0;

int TREE_DATA_DISTORT=0;//scheme 0: no distortion; scheme 1



//map<string, int> ALIGN_SCORE_MATRIX; //20 amino acid pairwise alignment scores, plus gap scores (opening and extended) (20*21/2+2)
int ALIGN_SCORE_MATRIX[28][28];

map<string, string> DNA_CODON_TBL;
map<string, char> DNA_CODON_TBL_SL;

string QUERY_SEQ_FILE;

int END_EXON_LEN = 1;

int MID_EXON_LEN = 6;

int EXPLORE_END_EXON_LEN = 1000;
float EXPLORE_END_EXON_PER = 0.1;


int REPAIR_HSP_MIN_INIT_SCORE = 0; //41; //Modified: now always try to find new HSP during repair?

int REPAIR_HSP_EXTEND_SCORE_DROP = 16;

bool EXTEND_HSP_BY_SCORE;

int REPAIR_HSP_AFTER_EXON = 1;
bool VERBOSE = false;
bool GENBLASTG_NEED_PID = false;

bool OUTPUT_Protein = false;
bool OUTPUT_GFF = false;
bool OUTPUT_cDNA = false;
string USER_ID = "";
bool phase1_only = false;


bool GENEWISE_COMMAND = false;
bool EXONERATE_COMMAND = false;

map<string, pair<string, set< pair<int, int> > > > WormBase_Exons;

int Total_WormBase_Exons;

int Total_WormBase_Exons_Effective;

int Total_GenBlast_Exons;

int	Total_Matching_Exons;

int Total_Matching_Genes;



int num_of_correct_exons;

int num_of_incorrect_exons;

int num_of_missed_exons;



int num_of_edges;



int num_of_groups_pruned_due_to_overlap;


float max_gblastg_final_align_pid;
float max_wormbase_final_align_pid;

float max_genewise_final_align_pid;

int chromosome_start_pos = 1;//keep track of the start position (on the chromosome sequence) of the current group region

clock_t gblast_start_time; //used in combination with duration_gblastg_without_comp_final_align
double duration_gblastg_without_comp_final_align;


#endif /* GBLAST_H */

