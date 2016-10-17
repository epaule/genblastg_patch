//***************GRAPH.H************************

//

//Main functionality of the genBlast program is here

//Two graph classes are defined: one for first phase, one for second phase

//

//

//Author: Rong She

//Date: April 2007

//***********************************************



#ifndef GRAPH_H /* GRAPH_H */

#define GRAPH_H

#include <time.h>




#include "edge.h"

#include "data_manager.h"



extern map<string, pair<string, set< pair<int, int> > > > WormBase_Exons;

extern int Total_WormBase_Exons;

extern int Total_GenBlast_Exons;

extern int	Total_Matching_Exons;

extern int 	Total_WormBase_Exons_Effective;

extern int Total_Matching_Genes;



extern int num_of_correct_exons;

extern int num_of_incorrect_exons;

extern int num_of_missed_exons;



extern int num_of_edges;



extern int num_of_groups_pruned_due_to_overlap;


extern float max_gblastg_final_align_pid;

extern float max_wormbase_final_align_pid;

extern float max_genewise_final_align_pid;



extern clock_t gblast_start_time;

extern double duration_gblastg_without_comp_final_align;


//=================================================================================//

/*

//abstract base class

class BASE_GRAPH

{

public:

	int					num_nodes;

	vector<Edge>		edges;



	void Reset();

	void ShortestPath();

	virtual bool PrintPath(int* pred, bool* cut, float* penalty) = 0;



};

*/



class ACCP_DONR_graph //: public BASE_GRAPH

{

	int					num_nodes;

	vector<ACCP_DONR_Edge>		edges;

	vector<ACCP_DONR_Edge_w_Score> edges_w_score;

	vector<ACCP_DONR_Edge_Scores> edges_scores;

	DataManager&	dataMgr;



	//a set of following paras for each HSP group (candidate gene)

	int			start_site; //start position

	int			end_site; //end position

	//int			second_end_site; //this is only used in the new modified version where we test for the possibility of ending at the second-to-last exon during repair
	vector<ExonSiteInfo> second_end_temp_sites; //use this to record the bunch of temp_sites with second_end_site as gene end

	int			final_start_site; //for final gff output only, to hold the first exon start after repair and everything
	int			final_end_site; //for final gff output only, to hold the last exon end after repair etc.

	vector<int>	acceptors; //all acceptors between start and end, except those before the first donor

	vector<int> donors; //all donors between start and end, except those after the last acceptor



	vector<int> acceptor_region_index;
	vector<int> donor_region_index;

	vector<int> acceptor_region_index_copy; //keep the original region_index, ComputeExons() may change it, and we need the original to RemoveSegment etc.

	vector<int> donor_region_index_copy;

	int			hspID_added_for_repair;


	vector<ACCP_DONR_Edge> all_edges;



	//record the number of heading base pairs (number of base pairs before the HSP frame) for each acceptor

	//and the number of trailing base pairs for each donor

	vector<int>	acceptor_head;

	vector<int> donor_tail;



	//vector<int> gap_centers; //for version 2.3, no-intron case (donor/acceptor both in the center of gap)

	//multimap<int, pair<int, int> > donor_segments; //segments (areas) used for searching for donor/acceptors <segment_mid_point (where search starts), pair<front_end, back_end> >

	vector<SegmentsInThreeBounds> donor_segments_pair; 

	

	//vector<int> donor_acceptor_HSP_ID; //for S1, record one HSP for each donor_segment, this is the HSP to be used later in aligning donor and acceptor

	

	vector< vector<int> > donor_acceptor_HSP_ID;

	

	vector<int> donor_HSP_ID; //for S1, have to keep with each donor, in case the donor gets bumped up in its "region_index" in GenePredict4()

	vector<int> acceptor_HSP_ID; //for S1, for acceptors, same as donor_HSP_ID

	//however, for S2 algorithm, acceptor and donor that are in the same region may correspond to different HSPs! 

	//due to decision tree leaf class and lastleaf's class may be the same, thus we use donor_acceptor_HSP_ID only to record donor_hsp_ID, 

	//and we need acceptor_HSP_ID to record acceptor_hsp_ID separately

	//OR: for S2, we do NOT record hsp_IDs at all, we still get that by checking HSPs themselves!

	//WELL... that doesn't work very well, so we are back to recording everything

	//vector<int> acceptor_HSP_ID;



	map<int, int> hspID_to_hspIndex_map;
	set<int> hspID_due_to_stop; //stores the hsp_id that has in-frame stop between current hsp and prev hsp (was split due to "fixinframestop")
	

	vector<string > donor_prev_2nt;

	vector<string > acceptor_next_2nt;


	vector<int> next_acceptors; //store tmp acceptors that should be used in next HSP computation (for phase2 version 2 only)

//	vector<int> cur_acceptors;

	//int	last_exon_start; //2 special nodes that may need to be considered

	//int last_exon_end;

	//ACCP_DONR_Edge last_exon_edge; //revised: now is the first exon edge in prev_HSP_edges

	ACCP_DONR_Edge cur_exon_edge[3]; //last exon edge in curret solution

	//int cur_edge_penalty_index[3]; //is it necessary?



	//int last_edge_penalty_index[3]; //penalty index of the first exon edge in prev_HSP_edges

	vector<ACCP_DONR_Edge>	prev_HSP_edges; //store all previous edges in all possible paths in previous HSP

	vector<ACCP_DONR_Edge>	cur_HSP_edges;

	set<int> prev_HSP_edges_nodes;

	

	//int	prev_id[3];

	//int	prev_total[3];

	vector<int> prev_path[3]; //its index is the remainder that corresponds to each prev_path that results from the first HSP

	char prev_path_remainder[3]; 

	//AD_Edge_Penalty prev_dist[3];



	//int			special_site;	//the site that is always before all other sites (start_site-1)

								//(used only for constructing the first special "intron" edge for 2nd+ HSPs,

								//so that edges from this site are always at the front when edges are sorted)



	//int identity;

	//int totalpos;



	bool out_frame_merge; //signals if we need to do out_frame_merge

	bool in_frame_merge; //signals whether in_frame_merge has been finished



	//keep track of which hsp the acceptor/donor comes from

	vector<int> acceptor_hsp_no;

	vector<int> donor_hsp_no;

	

	//one group of HSPs for each frame (w.r.t. "start", start is 1)

	vector<HSP_Gene_Pair*>	HSP_frame[3];



	vector<HSP_Gene_Pair>	HSPs_dup;

	vector<HSP_Gene_Pair> HSPs_dup_copy;



	string query_seq;

	/*******************METHODS***********************/


	int CalcStartPos(int hsp_start, int gene_start, bool isPosStrand, 
		//pair<int, char*>& chr_seq); 
		vector<string>& chr_seq);

	int CalcEndPos(int hsp_end, int gene_end, bool isPosStrand, 
		//pair<int, char*>& chr_seq); 
		vector<string>& chr_seq);

	bool IsDonorAcceptorEmpty();

	bool GetSingleHSPDonorSegments_SplitOL_GapIsMinIntron(int cur_hsp_start, int& cur_hsp_end, int hsp_ID, 
		//multimap<int, pair<int, int> >& donor_segments, 
		int& last_segment_start, int& last_segment_hsp_ID, 
		int next_hsp_ID, int& next_hsp_start, int cur_query_start, int& cur_query_end, 
		int& next_query_start, int next_query_end, int& last_segment_end, bool& cur_hsp_exists);

	void Reset();

	void MergeOverlapHSPs(bool isPosStrand, vector<HSP_Gene_Pair*>& HSPs);

	bool TestMerge(vector<HSP_Gene_Pair*>::iterator ptrHSPIt, int& j, bool isPosStrand, 

		vector<HSP_Gene_Pair*>& HSPs_to_be_erased); // bool& remove_lastHSP);

	bool Merge(HSP_Gene_Pair* lastHSP, HSP_Gene_Pair* curHSP, int search_start, int cur_search_end, int last_search_end, 

		vector<HSP_Gene_Pair*>& HSPs_to_be_erased);

		//bool& remove_lastHSP);

	void CalcStartEndPos2(bool isPosStrand, vector<HSP_Gene_Pair*>& HSPs, 
		//pair<int, char*>& chr_seq);
		vector<string>& chr_seq);

	void PrintHeader(int rank, int count, int alt_count, multimap<Group_Info, vector<HSP_Gene_Pair*> >::iterator groupMapIt );

	bool GetSpliceSegments2(//multimap<int, pair<int, int> >& donor_segments, 

		vector<HSP_Gene_Pair*>& HSPs,

		int& first_segment_start, int& last_segment_start, int& last_segment_end);//, 

		//vector<int>& segment_hsp_start, int& last_hsp_start);


	void GetDonorsAcceptors_nPerBorder(int search_start, int search_end, 
		//pair<int, char*>& chr_seq,
		vector<string>& chr_seq, 
		vector<int>& accs, vector<int>& dons, vector<string>& acc_next_2nt, vector<string>& don_prev_2nt);


	int GetDonorHSP(vector<HSP_Gene_Pair*>& HSPs, int donor_site, int hsp_ID,  int segment_index, 

		Input_Alignment& donor_align, int& align_start, int& align_end, char& border_query_aa_front, char& border_query_aa_end,

		vector<string>& chr_seq);
		//pair<int, char*>& chr_seq);

	int GetAcceptorHSP(vector<HSP_Gene_Pair*>& HSPs, int acceptor_site, int hsp_ID,  int segment_index, 

		Input_Alignment& acceptor_align, int& align_start, int& align_end, char& border_query_aa_front, char& border_query_aa_end, 

		vector<string>& chr_seq);
		//pair<int, char*>& chr_seq);


	bool ComputeSpliceAlignment_SplitOL1(Input_Alignment& donor_align, Input_Alignment& acceptor_align, int donor_start, int donor_end, 

		int acceptor_start, int acceptor_end, string& splice_align, int& donor_splice_end, int& acceptor_splice_start, 

		int& new_donor_end, int& new_acceptor_start, //bool alignment_update, 

		//char missing_query_aa,

		string donor_2nt, string acceptor_2nt, int donor_tail_num, int acceptor_head_num);

	void ComputeSpliceAlignment_SplitOL(Input_Alignment& donor_align, Input_Alignment& acceptor_align, int donor_start, int donor_end, 

		int acceptor_start, int acceptor_end, string& splice_align, int& donor_splice_end, int& acceptor_splice_start, 

		int& new_donor_end, int& new_acceptor_start);



	void EraseSites(vector<int>& sites, vector<int>& site_regions, vector<int>& site_head_tail, vector<int>& site_HSP_ID, 

		vector<string>& site_2nt, int border, bool erase_before_border);

	void EraseSites(vector<int>& sites, vector<int>& site_regions, vector<int>& site_head_tail, 

		vector<string>& site_2nt, int border, bool erase_before_border);

	void EraseSpecificSites(vector<int>& sites, vector<int>& site_regions, 

								 vector<int>& site_head_tail, vector<int>& site_HSP_ID, vector<string>& site_2nt, 

								 vector<int>& sites_to_be_erased);



	int LoadHSPAlignScores(HSP_Gene_Pair* hsp_it, vector<int>& tgt_gap_starts, vector<int>& tgt_gap_ends);



	bool CutHeadOrTrailGapsInAlignment(Input_Alignment& align, int& gene_s, int& gene_e, int& hsp_s, int& hsp_e);



	int GetDonorHSP_S2(vector<HSP_Gene_Pair*>& HSPs, int donor_site, int segment_index, 

								  Input_Alignment& donor_align, int& align_start, int& align_end, 

								  char& border_query_aa_front, char& border_query_aa_end, 

								  //pair<int, char*>& chr_seq);
								  vector<string>& chr_seq);

	int GetAcceptorHSP_S2(vector<HSP_Gene_Pair*>& HSPs, int acceptor_site, int segment_index, 

									   Input_Alignment& acceptor_align, int& align_start, int& align_end, 

									   char& border_query_aa_front, char& border_query_aa_end, 

										//pair<int, char*>& chr_seq);
										vector<string>& chr_seq);


	void PrintHSPs(vector<HSP_Gene_Pair*>& HSPs);

	void PrintHSPs(vector<HSP_Gene_Pair>& HSPs);


	void FixHSPsWithInFrameStop(bool isPosStrand, 

		vector<HSP_Gene_Pair*>& newHSPs, vector<HSP_Gene_Pair*>& newHSP_ptrs, int& max_alignment_HSP_ID);



	bool RepairHeadTailExon(bool isHeadExon, int target_site, int query_site, int target_end_site, int query_end_site, 

			 //pair<int, char*>& chr_seq, 
			 vector<string>& chr_seq, 
			 int& max_alignment_HSP_ID, vector<HSP_Gene_Pair*>& newHSP_ptrs, 

			 vector<HSP_Gene_Pair*>& blastHSPs, int blast_hsp_index, 

			 vector<ExonSiteInfo>& temp_sites, 

			 vector<HSP_Gene_Pair>& temp_HSPs, //vector<HSP_Gene_Pair*>& temp_HSPs, 

			 int exon_start_frame, 

			 vector<ExonSiteInfo>& additional_temp_sites, 
			 bool is_last_exon);



	void PreProcHSPs(bool isPosStrand, vector<HSP_Gene_Pair*>& HSPs, 

								vector<HSP_Gene_Pair*>& newHSP_ptrs, int& max_alignment_HSP_ID);

	bool ComputeExons(vector<HSP_Gene_Pair*>& HSPs, 
		//pair<int, char*>& chr_seq, 
		vector<string>& chr_seq, 
		int start_frame, int end_frame, 
		vector< vector< ExonSiteInfo > >& all_alternative_acceptors, 
		vector< vector< ExonSiteInfo > >& all_alternative_donors, bool& has_possible_exon, 
		bool repair_only);



	bool RepairMidExon(ExonSiteInfo prev_exon_start, ExonSiteInfo prev_exon_end,

						ExonSiteInfo cur_exon_start, ExonSiteInfo cur_exon_end, 

						//pair<int, char*>& chr_seq, 
						vector<string>& chr_seq, 
						int& max_alignment_HSP_ID, vector<HSP_Gene_Pair*>& newHSP_ptrs, 

						vector<HSP_Gene_Pair*>& blastHSPs, 

						//int prev_start_blast_hsp_index, int prev_end_blast_hsp_index, 

						//int cur_start_blast_hsp_index, int cur_end_blast_hsp_index,

						vector<ExonSiteInfo>& temp_sites, 

						vector<HSP_Gene_Pair>& temp_HSPs, 

						vector<ExonSiteInfo>& additional_temp_sites, 
						bool is_last_exon);//vector<HSP_Gene_Pair*>& temp_HSPs 



	void GetSpliceSegments1_forGenePredict6(//multimap<int, pair<int, int> >& donor_segments, 

						vector<HSP_Gene_Pair*>& HSPs, 

						 int& first_segment_start, int& last_segment_start, int& last_segment_end);

	bool LoadData4_forGenePredict6(vector<HSP_Gene_Pair*>& HSPs, 
		//pair<int, char*>& chr_seq, 
		vector<string>& chr_seq, 
		bool start_site_fixed, bool end_site_fixed);


	int FindDonorHSPIndex_UseLargestHSPInSegment(int donor_site, vector<HSP_Gene_Pair*>& HSPs, int donor_segment_index, 

					   int& cur_hsp_start, int& cur_hsp_end, int& gene_start, int& gene_end, int& j);

	int FindAcceptorHSPIndex_UseLargestHSPInSegment(int acceptor_site, vector<HSP_Gene_Pair*>& HSPs, int donor_segment_index, 

						int& cur_hsp_start, int& cur_hsp_end, int& gene_start, int& gene_end, int& j);



	void CombineSpliceSites_S1(int border, vector<int>& sites, vector<int>& final_sites, vector<int>& site_regions, 

							   int region, bool is_acceptor, vector<int>& site_head_tail, 

							   //int frame_reference, int hsp_ID, 

							   vector<HSP_Gene_Pair*>& HSPs, //use this to compute hsp_ID

							   //int segment_index, //it's the same as "region"!

							   //pair<int, char*>& chr_seq, 
							   vector<string>& chr_seq, 

							   map<int, int>& hsp_donor_stop_pos, //after HSP's end, position of first stop

							   map<int, int>& hsp_acceptor_stop_pos, //before HSP's start, position of first stop codon

							   vector<int>& acceptor_HSP_ID, vector<int>& donor_HSP_ID, 

							   vector<string>& site_2nt, vector<string>& final_site_2nt);



	void PrintDonorAcceptorHSPID();

	void UseLargestHSPInSegment(int cur_seg, vector<HSP_Gene_Pair*>& HSPs, int search_start, int search_end);

	bool RemoveHSPBeforeStop(vector<HSP_Gene_Pair*>& HSPs, int stop_site, 
		//pair<int, char*>& chr_seq, 
		vector<string>& chr_seq, 
		bool& removed_hsp_is_new_hsp, bool remove_after_stop);



	bool RemoveSegments(int hsp_ID, int hsp_index, vector<HSP_Gene_Pair*>& HSPs, 
		//pair<int, char*>& chr_seq);
		vector<string>& chr_seq);

	void EraseSites_ByRegions(vector<int>& sites, vector<int>& site_regions, 

								 vector<int>& site_head_tail, vector<int>& site_HSP_ID, vector<string>& site_2nt, 

								 int region_erase_start, int region_erase_end);

	void IncrementSiteRegions(vector<int>& site_regions, int region_no, int inc_count);


	void GetNewSites(vector<SegmentsInThreeBounds>& newsites_donor_segments_pair, 

								  int start_region_index, 

									 vector<HSP_Gene_Pair*>& HSPs, 
									 //pair<int, char*>& chr_seq, 
									 vector<string>& chr_seq, 

									 vector<int>& newsites_acceptors, vector<int>& newsites_donors, 

									 vector<int>& newsites_acceptor_region_index, vector<int>& newsites_donor_region_index, 

									 vector<int>& newsites_acceptor_head, vector<int>& newsites_donor_tail, 

									 vector<int>& newsites_acceptor_HSP_ID, vector<int>& newsites_donor_HSP_ID, 

									 vector<string>& newsites_acceptor_next_2nt, vector<string>& newsites_donor_prev_2nt);

	void InsertAccSites_ByRegions(vector<int>& sites, vector<int>& site_regions, 

								 vector<int>& site_head_tail, vector<int>& site_HSP_ID, vector<string>& site_2nt, 

								 int region_start_insert);

	void InsertDonSites_ByRegions(vector<int>& sites, vector<int>& site_regions, 

								 vector<int>& site_head_tail, vector<int>& site_HSP_ID, vector<string>& site_2nt, 

								 int region_start_insert);

	void OutputDonorAcceptorWithAllInfo(bool isPosStrand);

	void UseOldExon(vector<ExonSiteInfo>& temp_sites, vector<ExonSiteInfo>& additional_temp_sites, 
								 vector<HSP_Gene_Pair>& temp_HSPs, 
								 int target_site, int query_site, int target_end_site, int query_end_site, 
								 int blast_hsp_index, vector<HSP_Gene_Pair*>& blastHSPs);
	void UseOldExon_Mid(vector<ExonSiteInfo>& temp_sites, vector<ExonSiteInfo>& additional_temp_sites, 
								 vector<HSP_Gene_Pair>& temp_HSPs, 
									 ExonSiteInfo& cur_exon_start, ExonSiteInfo& cur_exon_end, 
									 ExonSiteInfo& prev_exon_start, 
									 int cur_start_blast_hsp_index, int cur_end_blast_hsp_index, 
									 vector<HSP_Gene_Pair*>& blastHSPs);

	bool PostProcessSpliceSites(int first_segment_start, int last_segment_end, 
		//pair<int, char*>& chr_seq);
		vector<string>& chr_seq);

	void PrintGFFExons(bool isPosStrand, int exon_start, int exon_end, int exon_count, int alt_count);
	void CutHeadOrTrailGapsInHSPs(vector<HSP_Gene_Pair*>& HSPs); //this calls CutHeadOrTrailGapsInAlignment() to process one HSP at a time
	void SortHSPsToColinear(bool isPosStrand, vector<HSP_Gene_Pair*>& HSPs);
	bool CheckSglExon(vector<HSP_Gene_Pair*>& HSPs, 
		//pair<int, char*>& chr_seq, 
		vector<string>& chr_seq, 
		int start_frame, int end_frame, int cur_exon_start_site, 
		vector< vector< ExonSiteInfo > >& all_alternative_acceptors, 
		vector< vector< ExonSiteInfo > >& all_alternative_donors, 
		bool& has_possible_exon, bool repair_only);

	void GetPredGeneSeq(vector<ExonSiteInfo>& temp_sites, 
		//pair<int, char*>& chr_seq, 
		vector<string>& chr_seq, 
		string& final_alignment, //predicted protein
		string& temp_DNA_string, bool& found_stop); //predicted cDNA
	align_compare GetBetterPredGeneSeq(vector<ExonSiteInfo>& temp_sites, vector<ExonSiteInfo>& prev_temp_sites, 
						  //pair<int, char*>& chr_seq,
						  vector<string>& chr_seq, 
						string& final_alignment, string& prev_alignment, //predicted protein
						string& final_DNA_string, string& prev_DNA_string, //predicted cDNA
						Input_Alignment& newAlign, Input_Alignment& prev_align, 
						bool& found_stop); 

	bool RemoveHSPAndResetComputeExons(vector<HSP_Gene_Pair*>& HSPs, int stop_site, 
		//pair<int, char*>& chr_seq, 
		vector<string>& chr_seq, 
		bool remove_after_stop_preferred, 
		bool& has_possible_exon, 
		vector< vector< ExonSiteInfo > >& all_alternative_acceptors,
		vector< vector< ExonSiteInfo > >& all_alternative_donors);
	int ComputeSpliceAlignment_SplitOL1_Optimal(vector<HSP_Gene_Pair*>& HSPs, 
		//pair<int, char*>& chr_seq, 
		vector<string>& chr_seq, 
		IntronAndDonorTail& donor_acceptor_info, float& align_pid);
	bool GetHSPSpliceAlignment1(string& targetStr, string& queryStr, string& matchStr, 
		int cur_hsp_start, int cur_hsp_end, int gene_start, int gene_end, 
		int site_start, int site_end, //start and end of the splice segment to be obtained
		//the following FIVE are output from this function
		Input_Alignment& donor_align, int& align_start, int& align_end, char& border_query_aa_front, char& border_query_aa_end, 
		vector<string>& chr_seq); //border_query_aa records the char (protein amino acid, single letter) at the border (in case there's broken codon)
		//pair<int, char*>& chr_seq);
	bool ComputeSpliceAlignment_SplitOL1_ForGetHSPSpliceAlignment1(Input_Alignment& donor_align, Input_Alignment& acceptor_align,
		int donor_start, int donor_end, int acceptor_start, int acceptor_end, 
		string& splice_align, int& splice_align_score, 
		string& donor_2nt, string& acceptor_2nt, int donor_tail_num, int acceptor_head_num);
public:

	ACCP_DONR_graph(DataManager& dMgr): dataMgr(dMgr) 

	{

		num_nodes = 0;

		start_site = 0;

		end_site = 0;



		//identity = 0;

		//totalpos = 0;

		//cur_edge_penalty_index = -1; //signals no edge



	}



	void GenePredict6();


};



class HSP_graph //: public BASE_GRAPH

{

	int						num_nodes;

	vector<Edge>		edges; //list of edges, used for later SP algorithm



	list<int>		last_nodes;//used for creating ext/cut edges

	map<pair_dist, distance_info, pair_dist_compare>	all_pair_dist; //used for creating skip edges



	DataManager&			dataMgr;

	bool					isPosStrand;

	int						chr_index;



	map<int, float>			source_nodes; //all possible source nodes <node_id, penalty_before_node>

	multimap<float, vector<int> > local_short_paths; //store the local solutions for all possible source nodes



	float LocalShortestPath(int node_id, float dist_before_fst_node, vector<int>& path);

	void FindAllLocalShortestPaths();



	bool ExtCutPenalty(HSP_node* source, HSP_node* dest, float& gap_penalty, float& gap_penalty_nextgroup, float& overlap_penalty, vector<HSP_Gene_Pair>& curHSPs, 

		bool isSpecial, bool isSkipEdge);

	float GetCutPenalty_ToStart(int pos);

	float GetCutPenalty_ToEnd(int pos);

	float GetPenalty_Between(int startpos, int endpos, bool end_included);



	bool FindParent(int candidate, int child, vector<HSP_Gene_Pair>& curHSPs, set<int>& can_set);

	bool HasPath(int souce, int dest, vector<HSP_Gene_Pair>& curHSPs, map<pair<int,int>, bool>& map_HasPath);



	bool TooMuchOverlap(int source_start, int source_end, int dest_start, int dest_end, bool& pos_cut, bool& neg_cut);

	bool TooFar(int source_hsp_end, int dest_hsp_start);



	bool IsNextGeneSeg_old(HSP_node* source, HSP_node* dest); 

	bool IsNextGeneSeg_new(HSP_node* source, HSP_node* dest); 



	void AddExtCutEdge(HSP_node* source, HSP_node* dest, vector<HSP_Gene_Pair>& curHSPs, bool isSpecial);

	bool AddSkipEdge(HSP_node* source, HSP_node* dest, vector<HSP_Gene_Pair>& curHSPs);



	void ProcHSPwithNextSeg(int index, vector<HSP_Gene_Pair>& curHSPs);

	bool TrySkipEdge(vector<PID_HSPstart_Index>& curPIDs, int i, vector<HSP_Gene_Pair>& curHSPs, int& lastHSPstart);



	void CalcSkipDist_SP(int source_num, vector<HSP_Gene_Pair>& curHSPs);



	void CreateAllExtCutEdges(vector<HSP_Gene_Pair>& curHSPs);

	void CreateAllSkipEdges(vector<HSP_Gene_Pair>& curHSPs);

	void ShortestPath();







	bool PrintPath(int* pred, bool* cut, float* penalty, float* penalty_nextgroup);

	bool PrintPath_Backtrack(int* pred, bool* cut);



	void DestroyNodes(vector<HSP_Gene_Pair>& curHSPs)

	{

		for (int i=0; i<num_nodes; i++)

		{

			delete curHSPs[i].node;

			curHSPs[i].node = 0;

		}

	}



	//for debug

#ifdef	DEBUG

	ofstream				debugFile;

#endif

	int						cur_edge; //used for debug only

	void PrintEdges(ostream& os); //used for debug only



public:



	HSP_graph(DataManager& dMgr, bool isPos, int chr): dataMgr(dMgr), isPosStrand(isPos), chr_index(chr) 

	{

		num_nodes = 0;



		cur_edge = 0;





#ifdef	DEBUG


		string str = dMgr.inputFile;

		str += "_";

		str += dMgr.query_gene;

		str += "_";

		str += dMgr.HSP_chr[chr];

		if (isPos)

			str += "_posdebug";

		else

			str += "_negdebug";

		debugFile.open(str.c_str());

		if (!debugFile.is_open())

			cout << "debug file open error" << "\n";


#endif





	}



	~HSP_graph() { 		

		//delete all nodes

		if (isPosStrand)

			DestroyNodes(dataMgr.HSP_gene[chr_index]);

		else

			DestroyNodes(dataMgr.HSP_neg_gene[chr_index]);



#ifdef DEBUG

		if (debugFile.is_open())

			debugFile.close();

#endif

	}





	void GeneGrouping();



	int NodeCount() { return num_nodes - 2; } //for performance result output only, subtract two special nodes



	int EdgeCount() { return edges.size(); }



};



#endif /* GRAPH_H */

