//***************DATA_MANAGER.H************************
//
//Main class that controls the input/output and handles data storage and data conversions
//
//
//Author: Rong She
//Date: April 2007
//***********************************************

#ifndef DATA_MANAGER_H /* DATA_MANAGER_H */
#define DATA_MANAGER_H

//#define DEBUG_VERSION
//#define VERBOSE

#define COMMAND_WITH_BLAST //run it in command line mode (with blast integrated)
#define PERFORMANCE //record performance numbers in perform.txt

//#define TIMING

//GENEWISE_COMMAND no longer used in #define
//#define GENEWISE_COMMAND //run genewise from command line, requires target sequence file to be prepared from genBlastA result (one gene at a time!)
//#define GENEWISE_PERFORMANCE //record genewise command running time

#define GENBLASTG //do we need to run gblastg?
//#define GENBLASTG_NEED_PID //do we need to compute gblastg final alignment pid? (GENBLASTG must be defined)

//now replace GENBLASTG_NO_REPAIR by REPAIR_HSP_AFTER_EXON!
//#define GENBLASTG_NO_REPAIR //try this for EvsE set (when target genome is the same as query) (turn off REPAIR_HSP_AFTER_EXON)

//the following define is for evaluating against other predictions (not real truth), we compare their alignment PID
//#define GENEWISE //do we need to compute genewise final alignment pid?
//#define WORMBASE //do we need to compute wormbase final alignment pid?
//#define OTHERS_FINAL_ALIGN_KEEP_STOP //do we compute genewise/wormbase pid with internal stops kept?

//#define EVAL_ON_TRUE_EXONS //this is the case for EvsE wormbase exons (wormbase annotations are supposed to be REAL TRUTH)

//#define SKIP_VERSION

//#define DEBUG
//#define USE_CHAR_ARR_IN_STORE_SEQ

#define V23_EXTRA_IS_GAP
#define PID_BEFORE_SCORE
//#define USE_LAST_HSP_WHEN_HAS_STOP
//#define V23_ALL_HSP_GAP_SPLICE_SEGMENT
//#define COMPUT_EXON_FULL_STEP_BACK

#define MEMORY_LIMIT_SW 81920000 //80M, this roughly translates to 80M*4(sizeof(int))*3(num_of_arrays) ~ 960M memory, so about 1G memory is needed


#include "edge.h"
#include "scores.h"
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <list>
#include <numeric>
#include <sstream>
using namespace std;


#if defined(__GNUC__)
	#if defined(__APPLE__) && defined(__MACH__) /*mac os X*/
		#define IMPLICIT_LFS_64
		typedef fpos_t file_pos64; /*Mac does implicit large file support*/
	#else
		typedef fpos64_t file_pos64;
	#endif
#elif defined(_MSC_VER) || defined(__BORLANDC__)
	#define IMPLICIT_LFS_64
	typedef __int64 file_pos64;
#else
	#error This compiler is not supported (need GNU CC, MSC, or BORLANDC).
#endif


//=================================================================================//

/* the basic struct that stores HSPs */

struct HSP_Gene_Pair

{

	int					ID;

//	string				chr;



	int					gene_start;

	int					gene_end;

	int					HSP_start;

	int					HSP_end;

	float				pid;



	HSP_node*			node; //used for deleting the nodes

	bool				hasOverlap; //used for AddSkipEdge()



	HSP_Gene_Pair() { node = 0; }



	HSP_Gene_Pair(int i, int gstart, int gend, int hstart, int hend, float pid): ID(i), 

	  gene_start(gstart), gene_end(gend), HSP_start(hstart), HSP_end(hend), pid(pid), node(0) {

		  hasOverlap = false;

	  }



	//check whether two HSPs overlap

	bool Overlap(HSP_Gene_Pair* prev_HSP, int& gOL_start, int& gOL_end, 

		int& hOL_start, int& hOL_end, bool& g_align, bool& h_align);



	bool FitWithNeighborHSPs(HSP_Gene_Pair* neighborHSP, bool isPosStrand);



	bool HSPInRegion(int region_start, int region_end);



    bool operator<(const HSP_Gene_Pair& h) const

    { return HSP_start < h.HSP_start || (HSP_start == h.HSP_start && HSP_end < h.HSP_end); }



	friend ostream& operator<<(ostream& os, const HSP_Gene_Pair& h)

	{

		os << "HSP_ID["<< h.ID << "]:(" << h.HSP_start << "-" << h.HSP_end << ");" << "query:(" << h.gene_start << "-" << h.gene_end << "); pid: " << h.pid 

			<< "\n";



		return os;

	}

};





//=================================================================================//

/* stores computed HSP groups, for output purposes (used in DataManager) */

struct Group_Info

{

	float			score;

	int				HSP_start;

	int				HSP_end;

	int				gene_cover_len;



	int				chr_index;

	bool			isPosStrand;

	Group_Info(float s, int h_start, int h_end, int gene_cover, int chr, bool isPos): score(s), HSP_start(h_start), HSP_end(h_end), 
		gene_cover_len(gene_cover), chr_index(chr), isPosStrand(isPos) { }

	bool operator<(const Group_Info& g) const
	{	return score > g.score || (score == g.score && chr_index < g.chr_index ) 
		|| (score == g.score && chr_index == g.chr_index && HSP_start < g.HSP_start);
	}

	friend ostream& operator<<(ostream& os, const Group_Info& info)
	{
		os << "score: " << info.score << "; region: " << info.HSP_start << "-" << info.HSP_end 
			<< "; gene cover: " << info.gene_cover_len << "\n";
		return os;
	}

};

/*struct Group_Rank_Count
{
	//int	group_index;
	int count; //absolute position of the current group in all groups, start from 1
	int rank;

	Group_Rank_Count(int r, int c): rank(r), count(c) {}

};*/

struct Group_RegEnd_Count
{
	int region_end;
	vector<int> group_count;

	Group_RegEnd_Count(int e, vector<int>& v)
	{
		region_end = e;
		group_count = v;
	}
};

//struct Group_Rank_Count
struct Group_Count_RegStart_RegEnd
{
	//int	group_index;
	int count; //absolute position of the current group in all groups, start from 1
	//int rank;
	int region_start;
	int region_end;

	//Group_Rank_Count(int r, int c): rank(r), count(c) {}
	Group_Count_RegStart_RegEnd(int c, int s, int e): count(c), region_start(s), region_end(e) {}

	bool operator<(const Group_Count_RegStart_RegEnd& g) const
	{
		return (region_start < g.region_start) || (region_start == g.region_start && region_end < g.region_end) 
		|| (region_start == g.region_start && region_end == g.region_end && count < g.count);
	}
};

//=================================================================================//

/* stores input alignments from ".align" file  */

struct Input_Alignment

{

	string			query_align;

	string			target_align;

	string			match_align;



	Input_Alignment() {}



	Input_Alignment(string s1, string s2, string match): query_align(s1), target_align(s2), match_align(match) {}



	Input_Alignment(const Input_Alignment& align)

	{

		query_align = align.query_align;

		target_align = align.target_align;

		match_align = align.match_align;

	}

	int GetAlignmentScore(int start_index, int end_index)
	{
		int score = 0;
		bool opengap = false;
		bool opengap_cur = false;
		for (int i=start_index; i<end_index; i++)
		{
			score += similarity_score(query_align[i], target_align[i], opengap, opengap_cur);
			opengap = opengap_cur;
		}
		return score;
	}


	friend ostream& operator<<(ostream& os, const Input_Alignment& align)

	{

		os << "query:" << align.query_align << "\n";

		os << "match:" << align.match_align << "\n";

		os << "targt:" << align.target_align << "\n";

		return os;

	}

};


//******************************************************************//
struct align_compare {
	int align_score;
	float align_pid;

	align_compare(): align_score(0), align_pid(0) {}
	align_compare(int score, float pid): align_score(score), align_pid(pid) {}
	bool operator==(const align_compare& ac) const {return align_pid == ac.align_pid && align_score == ac.align_score; }

#ifdef PID_BEFORE_SCORE //DOESN'T WORK, this is defined in data_manager.h, which hasn't be defined yet! OUCH, darn circular declarations!
	bool operator>(const align_compare& ac) const {
		return align_pid > ac.align_pid || ((align_pid == ac.align_pid) && (align_score > ac.align_score));
	}
#else
	bool operator>(const align_compare& ac) const {
		return (align_score > ac.align_score) || ((align_score == ac.align_score) && (align_pid > ac.align_pid));
	}
#endif

	friend ostream& operator<<(ostream& os, const align_compare& ac)
	{
		os << "align_score:" << ac.align_score << "; align_pid:" << ac.align_pid << ";";
		return os;
	}
};


struct RankAlignment

{

	//float	pid;
	align_compare pid;

	int		len;



	int		match_minus_gap_mismatch;

//	float	pid_optimal;



	//Updated: this source_len ranking doesn't seem to be correct, so not used any more

//	int		source_len; //the total length of original alignments (donor_alignment and acceptor_alignment) (before splice align)



	//RankAlignment(float p, int l, int m)//, int sl)
	RankAlignment(align_compare p, int l, int m)
	{

		pid = p;

		len = l;

		match_minus_gap_mismatch = m;

//		source_len = sl;

	}

	bool operator==(const RankAlignment& ra) const

	{

		return pid == ra.pid && len == ra.len && match_minus_gap_mismatch == ra.match_minus_gap_mismatch;
			//&& source_len == ra.source_len;

	}



	bool operator<(const RankAlignment& ra) const

	{

//		return match_minus_gap_mismatch > ra.match_minus_gap_mismatch || 

//			(match_minus_gap_mismatch == ra.match_minus_gap_mismatch && pid > ra.pid) || 

//			(match_minus_gap_mismatch == ra.match_minus_gap_mismatch && pid == ra.pid && len < ra.len);

		return pid > ra.pid || (pid == ra.pid && match_minus_gap_mismatch > ra.match_minus_gap_mismatch) || 
			(pid == ra.pid && match_minus_gap_mismatch == ra.match_minus_gap_mismatch && len < ra.len); // || 
			//(pid == ra.pid && match_minus_gap_mismatch == ra.match_minus_gap_mismatch && len == ra.len && source_len < ra.source_len);

	}



	friend ostream& operator<<(ostream& os, const RankAlignment& ra)

	{

		os << "pid:" << ra.pid << "; align_len:" << ra.len << "; match-gap/mismatch:" << ra.match_minus_gap_mismatch;

			//<< "; source_len:" << ra.source_len;

		return os;

	}

};


//=================================================================================//

/* stores splice sites for each strand of each chromosome  */

/* positive strand use positive integers, negative strand use negative integers */

struct Splice_Sites

{

	vector<int>		acceptors; 

	vector<int>		donors;



	Splice_Sites() 

	{

		acceptors.clear();

		donors.clear();

	}



	void LoadSite(int site, bool isAcceptor)

	{

		if (isAcceptor)

			acceptors.push_back(site);

		else

			donors.push_back(site);

	}

};



//=================================================================================//



struct PairExonSiteInfo

{

	ExonSiteInfo prev_exon_start;

	ExonSiteInfo prev_exon_end;

	ExonSiteInfo cur_exon_start;

	ExonSiteInfo cur_exon_end;



	PairExonSiteInfo(ExonSiteInfo& p_exon_start, ExonSiteInfo& p_exon_end, ExonSiteInfo& c_exon_start, ExonSiteInfo& c_exon_end)

	{

		prev_exon_start = p_exon_start;

		prev_exon_end = p_exon_end;

		cur_exon_start = c_exon_start;

		cur_exon_end = c_exon_end;

	}



	bool operator==(const PairExonSiteInfo& pair_info) const

	{

		return prev_exon_start == pair_info.prev_exon_start && prev_exon_end == pair_info.prev_exon_end 

			&& cur_exon_start == pair_info.cur_exon_start && cur_exon_end == pair_info.cur_exon_end;

	}



	bool operator<(const PairExonSiteInfo& pair_info) const

	{

		return prev_exon_start < pair_info.prev_exon_start || 

			(prev_exon_start == pair_info.prev_exon_start && prev_exon_end < pair_info.prev_exon_end) || 

			(prev_exon_start == pair_info.prev_exon_start && prev_exon_end == pair_info.prev_exon_end && cur_exon_start < pair_info.cur_exon_start) || 

			(prev_exon_start == pair_info.prev_exon_start && prev_exon_end == pair_info.prev_exon_end && cur_exon_start == pair_info.cur_exon_start 

				&& cur_exon_end < pair_info.cur_exon_end);

	}

};



struct InfoAfterRepair

{

	bool		need_pop;

	vector<ExonSiteInfo>	temp_sites;

	vector<HSP_Gene_Pair>	temp_HSPs;



	InfoAfterRepair(bool pop, vector<ExonSiteInfo>& t_sites, vector<HSP_Gene_Pair>& t_HSPs)

	{

		need_pop = pop;



		int i;

		for (i=0; i<t_sites.size(); i++)

			temp_sites.push_back(t_sites[i]);

		for (i=0; i<t_HSPs.size(); i++)

			temp_HSPs.push_back(t_HSPs[i]);

	}

};



//=================================================================================//

/* THE MAIN CLASS THAT HANDLES ALL DATA INPUT/OUTPUT/CONVERSIONS(PREPARATIONS) */

class DataManager

{

private:

	/*******************Data Members*************************/
	//map<int, vector<int> >	gene_start_end_map; //used to compute query fragment scores

	map<int, vector<pair<float, int> > >	gene_start_end_map; //used to compute query fragment scores



	/*******************Methods*************************/
	//void ReadOneGeneSegment(int, int, int);

	void ReadOneGeneSegment(int, int, float);

	//void GetFragmentScores(vector<HSP_Gene_Pair>& curHSPs);

	void GetFragmentScores();

	void ProcHSPs(vector<HSP_Gene_Pair>& curHSPs, bool);

	void PrepareMaps(vector<HSP_Gene_Pair>& curHSPs, bool);
	
	void CompRegion(multimap<Group_Info, vector<HSP_Gene_Pair*> >::iterator groupIt, int& region_start, int& region_end);
	void UpdateGroupRegEnd(map<int, Group_RegEnd_Count>& cur_map, int region_start, int region_end, int group_count);
	void GetRegionSeq(vector<string>& seq, 
		//char* seq, int& seq_len
		char* line, int cur_line_start_pos, int cur_len);
	void StoreGroupSeq(vector<string>& seq, 
		//char* seq, int seq_len, 
		vector<int>& group, int reg_start_pos, int seq_count);
	bool GetRegion(int cur_pos, int reg_start_pos, int reg_end_pos,
							vector<string>& seq, 
							//char* seq, int& seq_len, 
							char* line, int line_len, 
							bool& reg_start_found, bool& reg_end_found);

public:


	/*******************Data Members*************************/

	string					inputFile; //*.report (blast report)

	string					outputFile;

	string					alignFile; //*.align (blast alignment)

	string					chrSeqFile; //chromosome (target) sequence file (FASTA format) (ONE big file for all chromosomes!)

	string					spliceFile; //result of genesplicer (result file header)

	string					queryFile; //query protein sequence



	ifstream				inputFile_is;

	int						headerIndex[9];

	bool					inputFile_Open;

	bool					cur_gene_start;

	bool					inputFile_Finish;

	string					next_query_gene;



	ifstream				alignFile_is;

	bool					alignFile_Open;

	bool					cur_align_start;

//	bool					alignFile_Finish;



	ofstream				outFile;

	ofstream				gff_os;
	ostringstream			gff_gene_str;
	ostringstream			gff_str;
	ostringstream			gff_geneinfo_str;

	ofstream				cDNA_os;
	ofstream				pred_protein_os;

	string					query_gene;

	int						query_len;



	int						chr_shortest_hsp_len; //use the shortest candidate gene length as reference for output

	int						chr_longest_cand_len;

	int						scale;



	multimap<Group_Info, vector<HSP_Gene_Pair*> > groups; //computed HSP groups
	
	//MODIFIED: finish all groups on the same chromosome, then move on to next chromosome
	//map<string, vector<Group_Rank_Count> > chr_groups;//keep track of the correpondence between chromosome name and groups, so only load one chromosome at a time
	//vector<string> chromosomes; //keep track of the processing order of chromosomes

	//MODIFIED: get all DNA regions for all groups at once
	map<int, pair<int, int> > group_start_seqno;
	vector< vector<string> > group_dna_regions;
	//vector<pair<int, char*> > group_dna_regions;
	map<string, file_pos64> chr_seq_map; //only load once for a target DNA fasta file

	string					cur_chr_name; //used only to store the chromosome name of current HSP group (current gene)
	char					cur_gene_id[MAX_CHAR_NAME]; //used only to store geneid (ID=gene*, * is rank), for GFF3 format output

	map<string, vector<string> > chr_name_seq; //map of <chromosome_name, chromosome_nucleotide_sequence>



	map<int, Input_Alignment>	input_alignments; // <HSP_realID, alignment> (sorted by HSP_ID)

	map<int, Input_Alignment>	input_alignments_HSPs_dup;



	map<string, Splice_Sites>	chr_splice_sites;

	

	

	vector<string>			HSP_chr; //chr info (header info)

	vector< vector<HSP_Gene_Pair> >	HSP_gene; //on positive strand

	vector< vector<HSP_Gene_Pair> >	HSP_neg_gene; //on negative strand





	map<int, float>			fragment_score_map; //<frag_start, frag_score>

	multimap<int, int>		gene_start_HSP_num_map; //<gene_start, HSP_num> (pos) <gene_end, HSP_num> (neg)



	map<string, string>		query_gene_seq;

	/*******************Methods*************************/


	void Reset()

	{

		gene_start_end_map.clear();

		query_len = 0;

		chr_shortest_hsp_len = INT_MAX;

		chr_longest_cand_len = 0;

		

		groups.clear();

		//chr_name_seq.clear();

		input_alignments.clear();

		chr_splice_sites.clear();

		HSP_chr.clear();

		HSP_gene.clear();

		HSP_neg_gene.clear();

		fragment_score_map.clear();

		gene_start_HSP_num_map.clear();

		gff_str.str(""); //clear gff_str content
		gff_gene_str.str("");
		gff_geneinfo_str.str("");
		
	}



	void ResetChrSeq()
	{

		chr_name_seq.clear();
		group_start_seqno.clear();
		for (int i = 0; i<group_dna_regions.size(); i++)
			//delete [] group_dna_regions[i].second;
			group_dna_regions[i].clear();
		group_dna_regions.clear();
	}



	DataManager(const string* fileNames, const char* outFilename, bool phase1_only) { //user can specify output filename

		inputFile = fileNames[0];

		chrSeqFile = fileNames[1];//chrFileName;

		alignFile = fileNames[2];//alignFileName;

		//spliceFile = spliceFileName;

		//queryFile = queryFileName;



		outputFile = outFilename;

		char buffer[256];

		if (SPLICE_SEGMENT_VERSION == 1)

		{

			sprintf(buffer, "_%s%s_%s_s%d_%d_%d_%d", PHASE1_VERSION.c_str(), PHASE1_OUTPUT.c_str(), 

				PHASE2_VERSION.c_str(), SPLICE_SEGMENT_VERSION, REPAIR_HSP_MIN_INIT_SCORE, REPAIR_HSP_EXTEND_SCORE_DROP, 

				REPAIR_HSP_AFTER_EXON);

		}

		else

		{

			if (EXTEND_HSP_BY_SCORE)

			sprintf(buffer, "_%s%s_%s_s%d_tdshift%d_tddis%d_tcls%.1f_m%d_score_i%d_d%d_%d", PHASE1_VERSION.c_str(), PHASE1_OUTPUT.c_str(), 

				PHASE2_VERSION.c_str(), SPLICE_SEGMENT_VERSION, TREE_DATA_SHIFT, TREE_DATA_DISTORT, TREE_CLS_THRESHOLD, MINOBJS, 

				REPAIR_HSP_MIN_INIT_SCORE, REPAIR_HSP_EXTEND_SCORE_DROP, REPAIR_HSP_AFTER_EXON);

			else

			sprintf(buffer, "_%s%s_%s_s%d_tdshift%d_tddis%d_tcls%.1f_m%d_pid_i%d_d%d_%d", PHASE1_VERSION.c_str(), PHASE1_OUTPUT.c_str(), 

				PHASE2_VERSION.c_str(), SPLICE_SEGMENT_VERSION, TREE_DATA_SHIFT, TREE_DATA_DISTORT, TREE_CLS_THRESHOLD, MINOBJS, 

				REPAIR_HSP_MIN_INIT_SCORE, REPAIR_HSP_EXTEND_SCORE_DROP, REPAIR_HSP_AFTER_EXON);



		}



		outputFile += buffer;



	

		outFile.open(outputFile.c_str()); 

		if (!outFile.is_open())

			cout << "output file open error!" << "\n";



		//if (hasGFF) //hasGFF is not used any more
		if (OUTPUT_GFF)

		{

		string gff_filename = outputFile;

		gff_filename += ".gff";

		gff_os.open(gff_filename.c_str());

		if (!gff_os.is_open())

			cout << gff_filename << " file open error!" << "\n";

		else

			gff_os << "##gff-version 3\n";

		if (!phase1_only)
		{
			if (OUTPUT_cDNA)
			{
			gff_filename = outputFile;
			gff_filename += ".DNA";
			cDNA_os.open(gff_filename.c_str());
			if (!cDNA_os.is_open())
				cout << gff_filename << " file open error!" << "\n";
			}

			if (OUTPUT_Protein)
			{
			gff_filename = outputFile;
			gff_filename += ".pro";
			pred_protein_os.open(gff_filename.c_str());
			if (!pred_protein_os.is_open())
				cout << gff_filename << " file open error!" << "\n";
			}
		}

		}




		query_len = 0; //query gene length, initialized to 0, used only to compute the percentage of gene_cover in output

		query_gene = "";

		next_query_gene = "";



		chr_shortest_hsp_len = INT_MAX;

		chr_longest_cand_len = 0;



		inputFile_Open = false;

		cur_gene_start = false;

		inputFile_Finish = false;



		alignFile_Open = false;

//		cur_align_start = false;

//		alignFile_Finish = false;



	}



	~DataManager() 

	{

		if (outFile.is_open())

			outFile.close(); 



		if (gff_os.is_open())

			gff_os.close();


		if (cDNA_os.is_open())
			cDNA_os.close();

		if (pred_protein_os.is_open())
			pred_protein_os.close();
	}



	bool ReadFile();

	bool ReadFile_Skip(); //for debug skip





	void PrepareData(bool isPosStrand, int chr_index);



	void PrintHSPs(int chr_index) //for  debug

	{

		cout << "pos HSPs: " << "\n";

		int i, j;

		j = HSP_gene[chr_index].size();

		for (i=0; i<j; i++)

			cout << HSP_gene[chr_index][i];



		cout << "neg HSPs: " << "\n";

		j = HSP_neg_gene[chr_index].size();

		for (i=0; i<j; i++)

			cout << HSP_neg_gene[chr_index][i];

	}



	void PrintGroups(bool printOverview);



	void PrintFragmentScores(ostream& os);





	int	GetHSPID(bool isPosStrand, int chr_index, int index);



	void PrepareOutputHSPs();

	void GetAlignments_Skip(); //for debug skip

	void GetChromosomes(map<string, map<int, Group_RegEnd_Count> >& chr_group_region_count);
	void GetChromosome(string& chr_name); //get one chromosome

	void GetAlignments(set<int>& HSP_IDs); //read ".align" file and store needed alignments

	void GetSpliceSites(map<string, pair<int, int> >& chr_coordinates);



	int ComputeScale();



	void PrintGroupsTxt(multimap<Group_Info, vector<HSP_Gene_Pair*> >::iterator mapIt, int rank, int count);

	void PrintGroupsOverview(multimap<Group_Info, vector<HSP_Gene_Pair*> >::iterator mapIt);



	void OutputStrsWithWrap(vector<string>& strs);

	void OutputStrsWithWrapHeader(vector<string>& strs, vector<string>& headers);





	

	void LoadQuerySeqData(const char* filename);

	void GetQuerySeq(string& query_str);

};



extern double GetLocalAlignment(string& seq_a, string& seq_b, int seq_a_start, int seq_b_start, int newId, 

					   HSP_Gene_Pair& newHSP, Input_Alignment& newAlign, bool close_to_start, ofstream& debugFile, 

					   int* align_pos);

extern int GetGlobalAlignment(string& seq_a, string& seq_b, Input_Alignment& newAlign, float& align_pid, ofstream& debugFile);
extern int GetGlobalAlignment_scorematrix(string& seq_a, string& seq_b, Input_Alignment& newAlign, float& align_pid, ofstream& debugFile);
extern int GetGlobalAlignment_PID_CompScore(string& seq_a, string& seq_b, Input_Alignment& newAlign, float& align_pid, ofstream& debugFile);

extern void ExtendAlignment(double max_score, string& query_seq_str, string& target_seq_left, string& target_seq_right, 

		//pair<int, char*>& chr_seq, 
		vector<string>& chr_seq, 
		bool isPosStrand,

		HSP_Gene_Pair* curHSP, Input_Alignment& curAlign, ofstream& debugFile,

		int target_left_limit, int target_right_limit);







#endif /* DATA_MANAGER_H */

