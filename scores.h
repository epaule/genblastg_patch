//***************SCORES.H************************

//

//Some utility functions used to compute scores, map functions, string functions

//

//

//Author: Rong She

//Date: April 2007

//***********************************************




#if defined(_MSC_VER)
	#pragma warning(disable: 4786)
	#pragma warning(disable: 4503)
#endif


#ifndef SCORES_H /* SCORES_H */

#define SCORES_H



#include "c45extern.h"

#include "struct.h"



#include <string>

#include <set>

#include <map>

#include <vector>

#include <fstream>

#include <iostream>
#include <algorithm>
using namespace std;



//not currently used parameters

const float HSP_ALPHA = 0.5f;

const float GENE_BETA = 0.5f;



const float GENE_OL_GAMA = 0;



//used parameters
extern float HSP_SK_ALPHA;

extern float GENE_MS_BETA;



extern float MAX_OVERLAP_PER; //no longer used

extern int MAX_GAP_BTWN_HSP;

extern int TOP_RANK_NUM;

extern float MIN_GENE_COVER_PER;

extern float MIN_GENE_SCORE;

extern bool USE_MIN_GENE_SCORE;



extern int MIN_EXON_SIZE;

extern int MIN_INTRON_SIZE;

extern int MAX_EXON_SIZE;

extern int MAX_INTRON_SIZE;

extern int MAX_EXON_HSP_DIFF;

extern int MAX_INTRON_HSP_DIFF;

//extern int MIN_INTRON_LEN;
extern int START_CODON_CHECK_LEN;
extern int MAX_NUM_SPLICE_SITES;
extern int SPLICE_LEN;

//extern int MIN_INTERNAL_EXON_LEN;

extern string PHASE1_VERSION;

extern string PHASE1_OUTPUT;

extern string PHASE2_VERSION;

extern int SPLICE_SEGMENT_VERSION;

extern float TREE_CLS_THRESHOLD;

extern int TREE_DATA_SHIFT;

extern int TREE_DATA_DISTORT;



extern int ALIGN_SCORE_MATRIX[28][28]; //map<string, int> ALIGN_SCORE_MATRIX;

extern string QUERY_SEQ_FILE;
extern int END_EXON_LEN;

extern int MID_EXON_LEN;

extern int EXPLORE_END_EXON_LEN;
extern float EXPLORE_END_EXON_PER;


extern int REPAIR_HSP_MIN_INIT_SCORE;

extern int REPAIR_HSP_EXTEND_SCORE_DROP;

extern bool EXTEND_HSP_BY_SCORE;

extern int REPAIR_HSP_AFTER_EXON;
extern bool VERBOSE;
extern bool GENBLASTG_NEED_PID;

extern bool OUTPUT_Protein, OUTPUT_GFF, OUTPUT_cDNA;
extern string USER_ID;
extern bool phase1_only;
extern bool GENEWISE_COMMAND;
extern bool EXONERATE_COMMAND;

extern 	map<string, multiset<long> > genewise_ug_exons;

extern 	map<string, multiset<long> > genewise_g_exons;

extern map< string, map< pair<long, long>, multiset<long> > > wormbase_exons;
extern ofstream of_perform;

extern bool GAPPED;

extern int	chromosome_start_pos; //keep track of the start position (on the chromosome sequence) of the current group region

//extern bool GBLASTG_FOUND;



//=================================================================================//

/* CONSTANT DEFINITIONS */

const int DELIMIT = '\t';

const string OUT_EXT = ".gblast";

const int MAX_CHAR_PER_LINE = 80;

const int MAX_CHAR_NAME = 32;

const int MAX_LINE = 1<<12; //4096;//each chromosome string has up to MAX_LINE chars (2^12)





const string START_CODON = "atg";

const string STOP_CODONS[3] = {"tag", "tga", "taa"};

const string START_CODON_REV = "gta";

const string STOP_CODONS_REV[3] = {"gat", "agt", "aat"};



const string DONOR_SITE[2] = {"gt", "gc"};

const string ACCEPTOR_SITE[1] = {"ag"};

const string DONOR_SITE_REV[2] = {"tg", "cg"};

const string ACCEPTOR_SITE_REV[1] = {"ga"};





//utility functions

extern float Score(int, int, float);

extern float Penalty_old(float HSP_skip_penalty, float gene_missing_penalty );

extern float Penalty_new(float HSP_skip_penalty, float gene_missing_penalty, float gene_overlap_penalty );



extern void MyMapInsert(const int pos, const int index, const bool isStart, map<int, vector<int> >& myMap);
extern void MyMapInsert( const string sort_key, const long value, map<string, multiset<long> >& myMap);

extern void MyMapInsert( const string sort_key, const long value, multimap<string, multiset<long> >& myMap);
extern void MyMapInsert( const int pos, const int pos_end, const float PID, const bool isStart, map<int, vector<pair<float, int> > >& myMap);
extern void MyMapInsert( const string sort_key, const string chr_seq, const int start, const int end, 
				 map<string, pair<string, set< pair<int, int> > > >& myMap);

extern int CountCharInStr(const string& str, char ch);

extern void GetSubstrFromVecStrs(vector<string>& chr_seq, bool isPosStrand, const int start_pos, const int len, string& resultStr);

extern int GetSubstrFromVecStrs_ForRepair(vector<string>& chr_seq, bool isPosStrand, const int start_pos, const int len, string& resultStr);

extern void GetSubstrFromVecStrs(vector<string>& chr_seq, bool isPosStrand, const int start_pos, const int len, string& resultStr, 
								 string& prev_2nt, string& next_2nt);
extern void GetSubstrFromVecStrs_ForRepair(vector<string>& chr_seq, bool isPosStrand, const int start_pos, const int len, string& resultStr, 
								 string& prev_2nt, string& next_2nt);

extern void GetSubstrFromVecStrs_NegRev(vector<string>& chr_seq, bool isPosStrand, const int start_pos, const int len, string& resultStr);
extern int GetSubstrFromVecStrs_NegRev_ForRepair(vector<string>& chr_seq, bool isPosStrand, const int start_pos, const int len, string& resultStr);

extern void StrToComplement(string& resultStr);

extern void StrToReverseComplement(string& resultStr);



extern void BackwardSearch(string& segStr, int start_pos, bool isPos, const string* matchStrs, int num_matchStrs, vector<int>& results, bool is_donor, ofstream& outFile);

extern void ForwardSearch(string& segStr, int start_pos, bool isPos, const string* matchStrs, int num_matchStrs, vector<int>& results, bool is_donor, ofstream& outFile);

//extern void BackwardSearch(string& segStr, int start_pos, bool isPos, const string* matchStrs, int num_matchStrs, vector<int>& results, bool is_donor);

//extern void ForwardSearch(string& segStr, int start_pos, bool isPos, const string* matchStrs, int num_matchStrs, vector<int>& results, bool is_donor);

extern void BackwardSearch_nPerBorder(string& segStr, int start_pos, bool isPos, const string* matchStrs, int num_matchStrs, vector<int>& results, bool is_donor, vector<string>& site_2nt, ofstream& outFile, 

									  string& prev_2nt, string& next_2nt);

extern void ForwardSearch_nPerBorder(string& segStr, int start_pos, bool isPos, const string* matchStrs, int num_matchStrs, vector<int>& results,

				   bool is_donor, vector<string>& site_2nt, ofstream& outFile, string& prev_2nt, string& next_2nt);

extern map<string, string> DNA_CODON_TBL;

extern map<string, char> DNA_CODON_TBL_SL;

extern void LoadCodonMap();

extern void LoadCodonMap_SingleLetterAA();

extern void TranslateDNACodons(string& dna, string& protein);

extern bool FindRealPos(string& alignStr, int align_start, int align_end, int search_start, int search_end,

						int& search_start_pos, int& search_end_pos, int& extra_front, int& extra_end);

extern void GetPenaltyFromMatchStr(string& alignStr, int align_start, int align_end, int search_start, int search_end,

						string& matchStr, int& exact_match, int& gap, int& pos_match,

						int query_start, string& queryStr, set< pair<int, int> >& query_match_pos, bool& end_is_match);

extern void GetPenaltyFromMatchStr(int search_start_pos, int search_end_pos, int extra_front, int extra_end,

							string& matchStr, int& exact_match, int& gap, int& pos_match);

extern void FindTgtStrPos(string& targetStr, int cur_pos, int search_start, int search_end, int& start_index, int& end_index);

//extern void FillFrontExtArea(string& matchStr, int exon_index, int& ext, vector<string>& chr_seq, int exon_pos, 

//				 vector<string>&  strChrs, string& queryStr, string& targetStr);

//extern void FillEndExtArea(string& matchStr, int exon_index, int& ext, vector<string>& chr_seq, int exon_pos,

//						   vector<string>&  strChrs, string& queryStr, string& targetStr);

//extern void FillExonArea(vector<string>&  strChrs, string& queryStr, string& matchStr, string& targetStr, vector<string>& chr_seq, 

//				  bool isPosStrand, int exon_pos, int front_gap, int frame, int start_index, int end_index);

//extern void FillFromTargetStr(vector<string>& chr_seq, bool isPosStrand, int start_pos, int length, bool revStr, vector<string>&  strChrs);

//extern void FillFrontFixedLenExt(vector<string>& chr_seq, int& ext, int exon_pos, vector<string>& strChrs);

//extern void FillEndFixedLenExt(vector<string>& chr_seq, int& ext, int exon_pos, vector<string>& strChrs);

extern void StrToLower(string& str);

extern int LenOfStrVec(vector<string>& str_vec);

extern void ProteinToDNAStyle(string& str);

extern int ProteinPosToDNAStylePos(int pos);
extern int DNAToProteinStylePos_StartSite(int pos);

extern int DNAToProteinStylePos_EndSite(int pos);


extern void GetStopCodonIndexes(string& cur_chr_seq, vector<int>& stop_indexes, const string* stop_codons);



extern int FindHSPNo(vector<int>& sites, vector<int>& site_hsp_no, int site);

extern int FindMax(vector<int>& hsp_no);

extern bool RegionOverlap(int cur_region_start, int cur_region_end, set< pair<int,int> >& regions);



extern int GetSingleStrFromVecStrs(vector<string>& chr_seq, bool isPosStrand, 

								   const int vecIndex, const int start_pos, string& resultStr, int& leftover);



extern bool IsSorted(vector<int>& vec);

extern bool IsUnique(vector<int>& vec);



extern void BackwardSearch(string& segStr, int start_pos, bool isPos, const string* matchStrs, int num_matchStrs, 

						   vector<int>& results, bool is_donor, ofstream& outFile, int seg_index, vector<int>& segs);

extern void ForwardSearch(string& segStr, int start_pos, bool isPos, const string* matchStrs, int num_matchStrs, 

						  vector<int>& results, bool is_donor, ofstream& outFile, int seg_index, vector<int>& segs);

//extern void CombineInFrameSpliceSites(int border, vector<int>& sites, vector<int>& final_sites, 

//									  vector<int>& site_regions, int region, bool is_acceptor, vector<int>& site_head_tail);
extern void CombineSpliceSites(int border, vector<int>& sites, vector<int>& final_sites, vector<int>& site_regions, 

							   int region, bool is_acceptor, vector<int>& site_head_tail, int frame_reference, 

							   vector<string>& site_2nt, vector<string>& final_site_2nt, ofstream& os);

//extern void CombineSpliceSites_S1(int border, vector<int>& sites, vector<int>& final_sites, vector<int>& site_regions, 

//							   int region, bool is_acceptor, vector<int>& site_head_tail, int frame_reference,

//							   int hsp_ID, vector<int>& acceptor_HSP_ID, vector<int>& donor_HSP_ID, 

//							   vector<string>& site_2nt, vector<string>& final_site_2nt);



extern bool FindRealPos_DNA(string& alignStr, int align_start, int align_end, int search_start, int search_end,

						int& search_start_pos, int& search_end_pos);



extern void GetGaps(const string& str, vector<int>& gap_starts, vector<int>& gap_ends);

extern int GetMatches(const string& str, vector<int>& match_starts, vector<int>& match_ends, vector<int>& id);

extern int ConvertToNoGappedPos(int pos, vector<int>& gap_starts, vector<int>& gap_ends, bool need_after_gap_pos);

extern int ConvertToGappedPos(int pos, vector<int>& gap_starts, vector<int>& gap_ends);

extern int ComputeBackId(int pos, vector<int>& match_starts, vector<int>& match_ends, vector<int>& match_id, int total_id);

extern int ComputeFrontId(int pos, vector<int>& match_starts, vector<int>& match_ends, vector<int>& match_id, int total_id);



extern void GetSegments(Tree T, vector<int>& gap_starts, vector<int>& gap_ends, int hsp_start, 

						//multimap<int, pair<int, int> >& donor_segments, 

						vector<SegmentsInThreeBounds>& donor_segments_pair, 

						vector< vector<int> >& donor_acceptor_HSP_ID, //vector<int>& donor_acceptor_HSP_ID, //vector<int>& acceptor_HSP_ID, 

						int hsp_ID, //int& prev_hsp_ID, 

						int& last_segment_start, int& last_segment_end, 

						bool& lastLeafIsExon); //, vector<int>& segment_hsp_start, bool& isFirstLeaf);



extern void ProcSegmentsByLens(multimap<int, pair<int, int> >& donor_segments, 

							   int& last_segment_start, vector<int>& segment_hsp_start, int& last_hsp_start);



extern bool GetAlignPairStr(string& curStr, char t, char q, bool open_gap);

extern void remove_trailing_cr_lf(string& str);
extern int remove_trailing_cr_lf(char* str, int len);

extern bool all_white_space(string& str);

extern void DNA2AA(string& target_seq, int start_index, string& result_seq);

extern bool TranslateTargetDNA(string& str);

extern bool HasInFrameStop(string& str, bool frame_start_at_str_beginning, int& stop_pos, bool ignore_last_codon);

extern void EraseAll(string& str, char ch);



extern bool HasStopCodon(int acceptor_site, int donor_site, vector<string>& seq, int& frame, //int& site_before_stop, 
						 int& site_after_stop, 
						 string& align, int prev_donor_site, string& left_chars, 
						 bool output_DNA, ofstream& DNA_os, bool temp_DNA, string& temp_DNA_string, 
						 bool is_last_exon, bool& found_stop);

extern void HasStopCodon_KeepStopCodon(int acceptor_site, int donor_site, vector<string>& seq, int& frame, 

				  string& align, string& left_chars);


extern void CombineConsecutiveSpaces(string& str);

extern int similarity_score(char a,char b, bool opengap, bool& opengap_cur);

#endif /* SCORES_H */


