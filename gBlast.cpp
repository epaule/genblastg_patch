#define _CRTDBG_MAP_ALLOC

#if defined(linux) || defined(__linux)
	#include <sys/wait.h> //this only works on *nix systems, not on Windows
#endif 

#include "c45extern.h"
#include "gBlast.h"
#include <errno.h>

//#include <crtdbg.h>


//#include "FORTIFY.H"

//#include <direct.h>

#define SYNTAX_ERROR "genBlast command syntax error, please check genBlast usage by typing \"genblast\""

#define XD_FILE_ERROR "XDF file error"

#define WU_FILE_ERROR "wuBlast/Blast file error"

#define RPT_FILE_ERROR "intermediate file error, cannot parse wuBlast/Blast result"

#define FILE_READ_ERROR "file read error: "



map<string, multiset<long> > genewise_ug_exons, genewise_g_exons;

bool GAPPED=false;

map< string, map< pair<long, long>, multiset<long> > > wormbase_exons;

ofstream of_perform;

//bool GBLASTG_FOUND; //not used


int TestPara(bool use_score_or_pid, int HSP_init_score, int HSP_drop_score, int cls_sep_thr, int phase2_version);



void help()
{
	cout << "\ngenBlastG release v1.0.138\n";
	cout << "\nSYNOPSIS:\n";
	cout << "Given a list of query protein sequences and a target database that \n";

	cout << "consists of DNA sequences, this program runs wu-blast tblastn on the list \n";

	cout << "of sequences provided, then for each query, it predicts potential\n";

	cout << "genes that are homologous to the query. The output is ranked according\n";

	cout << "to their homology to the query.\n\n";

	cout << "Command line options:\n";
	cout << "\t-P\tSearch program used to produce HSPs,\n\t\tcan be either \"blast\" or \"wublast\", default is \"blast\",\n\t\toptional\n";
	cout << "\t-p\tspecifies the program option of genBlast: genblasta or genblastg\n";
	cout << "\t-q\tList of query sequences to blast, must be in fasta format,\n\t\trequired\n";

	cout << "\t-t\tThe target database of genomic sequences in fasta format,\n\t\trequired\n";

	//cout << "\t-p\tWhether query sequences are protein sequences (T/F)\n\t\t[default: T], optional\n";

	cout << "\t-e\tparameter for blast: The e-value, [default: 1e-2],\n\t\toptional\n";
	cout << "\t-g\tparameter for blast: Perform gapped alignment (T/F) \n\t\t[default: F], optional\n";
	cout << "\t-f\tparameter for blast: Perform filtering (T/F) [default: F],\n\t\toptional\n";
	cout << "\t-W\tparameter for blast: Set word size, 0 means using blast default [default: 0],\n\t\toptional\n";

	cout << "\t-a\tparameter for genBlast: weight of penalty for skipping HSPs,\n\t\tbetween 0 and 1 [default: 0.5], optional\n";

	//cout << "\t-b\tparameter for genBlast: weight of penalty for missing gene \n\t\tcoverage [default: 0.5], optional\n";

	cout << "\t-d\tparameter for genBlast: maximum allowed distance between HSPs \n\t\twithin the same gene, a non-negative integer [default: 100000],\n\t\toptional\n";

	cout << "\t-r\tparameter for genBlast: number of ranks in the output,\n\t\ta positive integer, optional\n";

	cout << "\t-c\tparameter for genBlast: minimum percentage of query gene \n\t\tcoverage in the output, between 0 and 1 ";

	cout << "(e.g. for 50%\n\t\tgene coverage, use \"0.5\"), optional\n";

	cout << "\t-s\tparameter for genBlast: minimum score of the HSP group in \n\t\tthe output, a real number, optional\n";
	cout << "\t-scodon\tThe number of base pairs to search for start codon within the region of HSP \
		group (inside the first HSP). If not specified, default is 15.\n";

	//for genblastg:
	//"genblastg (all genblasta options) [-i min-intron-length] [-x min-internal-exon-length] 
	//[-n max-num-splice-sites] 
	//[-v splice-region-version-s1ors2] [-h s2-tree-data-shift] [-j s2-tree-cls-threshold] 
	//[-norepair repair-or-not]
	//[-re end-exon-len] [-rm mid-exon-len] [-rl explore-end-exon-len]
	//[-rs repair-hsp-min-init-score] [-rd repair-hsp-extend-score-drop] 
	//[-gff] [-cdna] [-pro]
	cout << "\t-i\tparameter for genBlastG: minimum intron length, optional. If not specified, the default value is 15.\n";
	cout << "\t-x\tparameter for genBlastG: minimum internal exon length, optional. If not specified, default is 20.\n";
	cout << "\t-n\tparameter for genBlastG: maximum number of splice sites per region, optional. If not specified, default is 20.\n";
	cout << "\t-v\tgenBlastG splice region algorithm version: 1 or 2. Optional. If not specified, default is 1\n";
	cout << "\t-h\tparameter for genBlastG splice region algorithm2: data shift, optional. If not specified, default is 0.\n";
	cout << "\t-j\tparameter for genBlastG splice region algorithm2: class threshold, optional. If not specified, default is 0.\n";
	cout << "\t-norepair\tturn on the no-repair option of genBlastG\n";
	cout << "\t-re\tparameter for genBlastG repair process: minimum length of missing query region for repairing head or tail exon, ";
	cout << "optional. If not specified, default is 1.\n";
	cout << "\t-rm\tparameter for genBlastG repair process: minimum length of missing query region for repairing internal exon, ";
	cout << "optional. If not specified, default is 6.\n";
	cout << "\t-rl\tparameter for genBlastG repair process: length of DNA region before first exon or after last exon for searching ";
	cout << "additional alignments for repairing exons, optional. If not specified, default is 1000.\n";
	cout << "\t-rs\tparameter for genBlastG repair process: minimum initial score of alignment to be considered further (similar to BLAST), ";
	cout << "optional. If not specified, default is 41.\n";
	cout << "\t-rd\tparameter for genBlastG repair process: max allowed score reduction in extending initial alignment (similar to BLAST), ";
	cout << "optional. If not specified, default is 16.\n";
	cout << "\t-o\toutput filename, optional. If not specified, the output\n\t\twill be the same as the query filename with \".gblast\"\n\t\textension.\n";
	cout << "\t-gff\toutput options: turn on GFF output\n";
	cout << "\t-cdna\toutput options: turn on cDNA output\n";
	cout << "\t-pro\toutput options: output protein sequence of the predicted gene\n";
	cout << "\t-id\tThe GFF output user_id\n";
	cout << "\t-b\tTurn on the verbose on-screen output\n";
	cout << "\t-pid\tturn on final alignment PID computation (global alignment between predicted gene and query) in output.\n";

	cout << "\nExample:\n";

	cout << "genblast -p genblastg -q myquery -t mytarget -e 1e-2 -g T -f F -a 0.5 -d 100000 -r 10 -c 0.5 -s 0 ";
	cout << "-i 15 -x 20 -n 20 -v 2 -h 0 -j 3 -norepair -gff -cdna -pro ";
	cout << "-o myoutput\n";

	cout << "\n(Rong She\tLast updated: Nov. 2010)\n\n";
	//[-i min-intron-length] [-x min-internal-exon-length] 
	//[-n max-num-splice-sites] 
	//[-v splice-region-version-s1ors2] [-h s2-tree-data-shift] [-j s2-tree-cls-threshold] 
	//[-norepair repair-or-not]
	//[-re end-exon-len] [-rm mid-exon-len] [-rl explore-end-exon-len]
	//[-rs repair-hsp-min-init-score] [-rd repair-hsp-extend-score-drop] 
	//[-gff] [-cdna] [-pro]
}


int getOpt(int argc, char* argv[], int& a, string& str)
{
	a++;
	if (a >= argc)
	{
		cout << SYNTAX_ERROR << "\n";
		return -1;
	}
	str = argv[a];

	return 0;
}

int checkOpt_Int(int argc, char* argv[], int& a, int& val, int range_low, int range_high)
{
	int op, error;
	char* stopstring;
	char* option = argv[a]; //hold the "-x" option
	string tmpStr;
	op = getOpt(argc, argv, a, tmpStr); //a is incremented to point to the optionvalue
	if (op == 0)
	{
		error = errno;
		val = strtol(tmpStr.c_str(), &stopstring, 10);
		if (errno != error)
		{
			cout << "Cannot recognize option " << option << " " << argv[a] << "\n";
			op = -1;
		}
		else
			if (val < range_low || val > range_high)
			{
				cout << "invalid value for option " << option << " " << val << ", it must be between " 
					<< range_low << " and " << range_high << "\n";
				op = -1;
			}
	}
	return op;
}

int checkOpt_Float(int argc, char* argv[], int& a, float& val, float range_low, float range_high)
{
	int op, error;
	char* stopstring;
	char* option = argv[a]; //hold the "-x" option
	string tmpStr;
	op = getOpt(argc, argv, a, tmpStr); //a is incremented to point to the optionvalue
	if (op == 0)
	{
		error = errno;
		val = strtod(tmpStr.c_str(), &stopstring);
		if (errno != error)
		{
			cout << "Cannot recognize option " << option << " " << argv[a] << "\n";
			op = -1;
		}
		else
			if (val < range_low || val > range_high)
			{
				cout << "invalid value for option " << option << " " << val << ", it must be between " 
					<< range_low << " and " << range_high << "\n";
				op = -1;
			}
	}
	return op;
}

bool hasFile(const char* filename)

{

	bool flag = false;

	fstream fin;

	fin.open(filename, ios::in);

	if( fin.is_open() )

		flag=true;

	fin.close();



	return flag;

}



void LoadAlignScore_Error(string& line)

{

	cout << "unrecognized score in alignscore.txt: " << line << "\n";

	exit(-1);

}

//This whole function and probably the alignscore.txt should be modified for efficiency...
//but, it's a one-time load function, so not really big deal, let's use it for now :-)
void LoadAlignScore() //load scores from "alignscore.txt"
{
	ifstream scores("alignscore.txt");

	if (!scores.is_open())
	{
		cout << "alignscore.txt not found" << "\n";
		exit(-1);
	}

	string line;
	char *stopstring;
	int error;
	while (!scores.eof())
	{
		getline(scores, line);
		remove_trailing_cr_lf(line);

		if (line.empty() || line.find("#") == 0 || all_white_space(line)) //skip empty lines
			continue;

		int delimit_pos, count=0, cur_pos = 0;
		string curStr;
		int curScore;
		while ((delimit_pos = line.find(":", cur_pos)) != string::npos)
		{
			curStr = line.substr(cur_pos, delimit_pos);
			if (curStr.length() != 2) //matrix header must be a 2-char string
				LoadAlignScore_Error(line);

			error = errno;
			curScore = strtol(line.substr(delimit_pos+1).c_str(), &stopstring, 10 );//atoi(line.substr(delimit_pos+1).c_str() );
			if (errno != error) //error code set
				LoadAlignScore_Error(line);

			cur_pos = delimit_pos+1;
			count++;
		}

		if (count != 1) //make sure the line is in correct format
			LoadAlignScore_Error(line);

		//ALIGN_SCORE_MATRIX.insert(map<string, int>::value_type(curStr, curScore));
		int i;
		if (curStr == "-o")
			ALIGN_SCORE_MATRIX[26][0] = curScore;
		else
			if (curStr == "-e")
				ALIGN_SCORE_MATRIX[26][1] = curScore;
			else
				if (curStr == "Xx")
				{
					for (i=0; i<28; i++)
						ALIGN_SCORE_MATRIX['X'-'A'][i] = ALIGN_SCORE_MATRIX[i]['X'-'A'] = curScore;
				}
				else
					if (curStr == "*x")
					{
						for (i=0; i<27; i++)
						{
							ALIGN_SCORE_MATRIX[27][i] = ALIGN_SCORE_MATRIX[i][27] = curScore;
							//treat 'U' as '*' (stop codon)
							ALIGN_SCORE_MATRIX['U'-'A'][i] = ALIGN_SCORE_MATRIX[i]['U'-'A'] = curScore;
						}
					}
					else
						if (curStr == "**")
						{
							ALIGN_SCORE_MATRIX[27][27] = curScore;
							//treat 'U' as '*' (stop codon)
							ALIGN_SCORE_MATRIX['U'-'A']['U'-'A'] = ALIGN_SCORE_MATRIX['U'-'A'][27] = ALIGN_SCORE_MATRIX[27]['U'-'A'] = curScore;
						}
						else
							ALIGN_SCORE_MATRIX[curStr[0]-'A'][curStr[1]-'A'] = ALIGN_SCORE_MATRIX[curStr[1]-'A'][curStr[0]-'A'] = curScore;

		/*char ch;
		ch = curStr.at(0);
		curStr[0] = curStr[1];
		curStr[1] = ch;
		ALIGN_SCORE_MATRIX.insert(map<string, int>::value_type(curStr, curScore)); //add the reverse to map*/
	}
}





int main(int argc, char* argv[])
//int TestPara(bool use_score_or_pid, int HSP_init_score, int HSP_drop_score, int cls_sep_thr, int phase2_version)
{

	EXTEND_HSP_BY_SCORE = true;

//	REPAIR_HSP_MIN_INIT_SCORE = 41;

//	REPAIR_HSP_EXTEND_SCORE_DROP = 16;

//	TREE_CLS_THRESHOLD = -6;



/*	EXTEND_HSP_BY_SCORE = use_score_or_pid;

	REPAIR_HSP_MIN_INIT_SCORE = HSP_init_score;

	REPAIR_HSP_EXTEND_SCORE_DROP = HSP_drop_score;

	TREE_CLS_THRESHOLD = cls_sep_thr;

	SPLICE_SEGMENT_VERSION = phase2_version;

*/



//	GBLASTG_FOUND = false;



#ifdef PERFORMANCE
	Total_WormBase_Exons = 0;
	Total_WormBase_Exons_Effective = 0;
	Total_GenBlast_Exons = Total_Matching_Exons = 0;
	Total_Matching_Genes = 0;

	of_perform.open("perform.txt");

//	char perform_filename[64];
//	sprintf(perform_filename, "perform_S%d_%d_%d.txt", SPLICE_SEGMENT_VERSION, REPAIR_HSP_MIN_INIT_SCORE, REPAIR_HSP_EXTEND_SCORE_DROP);
//	of_perform.open(perform_filename);

	of_perform << "gene\tnum_of_HSPs\tnum_of_edges\ttime_of_1st_phase\tscore_of_1st_group\tnum_of_groups_pruned_due_to_overlap\t";
	of_perform << "num_of_edges_2nd_phase\ttime_of_2nd_phase\t"
		<< "time_2nd_phase_without_comp_pid\t#_of_correct_exons\t#_of_incorrect_exons\t#_of_missed_exons" << "\n";

	int num_of_HSPs_for_query;
	int num_of_edges_for_query;
#endif

	//syntax of genblast command: (genblasta)
	//"genblast -q <query protein, in fasta format> -t <target DNA database, fasta> [-e e-value(1e-2)] [-g T/F(T)] [-f T/F(F)] 
	//[-a alpha(0.5)] [-b beta(0.5)] [-d max-allowed-gap-distance-between-hsps-that-are-in-same-group(100000)] 
	//[-r num-of-top-ranks-in-output(INT_MAX)] [-c min-gene-cover-percentage-in-output-groups(0)] 
	//[-s min-score-of-output-groups(not used, 0 may be a choice)]"
	//now we must first read command line inputs, spit out options, check input errors (file format), 
	//create index file if it doesn't exist, pass to tblastn, then parse wublast result to create ".report", 
	//finally run genblast phase 1

	//for genblastg:
	//"genblastg (all genblasta options) [-i min-intron-length] [-x min-internal-exon-length] 
	//[-n max-num-splice-sites] 
	//[-v splice-region-version-s1ors2] [-h s2-tree-data-shift] [-j s2-tree-cls-threshold] 
	//[-norepair repair-or-not]
	//[-re end-exon-len] [-rm mid-exon-len] [-rl explore-end-exon-len]
	//[-rs repair-hsp-min-init-score] [-rd repair-hsp-extend-score-drop] 
	//[-gff] [-cdna] [-pro]

	string queryList="", targetDb="";
	string eValue="1e-2", gapSet="F", filterSet="F", proteinSeq = "T", tmpStr, outStr; //stuff we pass to tblastn
        int word_length=0;
	char *stopstring;
	int error;
	int op=-1;
	int a=1;
	bool SKIP = false;
	string skip_gene;
	//bool phase1_only = false;
//	int genewise_gene_index;
	string progStr;
	bool prog_is_wublast = false; //default is blast

	if (argc == 1) //"genBlast" only
	{
		help();
		return 0;
	}

	while (a<argc)
	{
		if (strcmp(argv[a], "-q")==0)
		{
			op = getOpt(argc, argv, a, queryList);
			outStr = queryList + OUT_EXT; //default output filename
		}
		else
		if (strcmp(argv[a], "-t")==0)
			op = getOpt(argc, argv, a, targetDb);
		else
		if (strcmp(argv[a], "-P")==0)
		{
			op = getOpt(argc, argv, a, progStr);
			if (op == 0)
			{
				std::transform(progStr.begin(), progStr.end(), progStr.begin(), ::tolower);
				if (progStr.compare("wublast") == 0)
				{
					prog_is_wublast = true;
				}
				else
				{
					if (progStr.compare("blast") == 0)
						prog_is_wublast = false;
					else
						op = -1; //invalid option?
				}
			}
		}
		else
/*		if (strcmp(argv[a], "-p")==0)
		{
			op = getOpt(argc, argv, a, proteinSeq);
			if (proteinSeq.compare("T") != 0 && proteinSeq.compare("F")!=0 )
			{
				cout << "Cannot recognize option -p " << proteinSeq << "\n";
				op = -1;
			}
		}
		else
*/		if (strcmp(argv[a], "-e")==0)
		{
			op = getOpt(argc, argv, a, eValue);
			//double e = strtod(eValue.c_str(), &stopstring); //wublast will verify this, so I'll skip testing
		}
		else
		if (strcmp(argv[a], "-g")==0)
			op = getOpt(argc, argv, a, gapSet);
		else
		if (strcmp(argv[a], "-f")==0)
			op = getOpt(argc, argv, a, filterSet);
		else
                if (strcmp(argv[a], "-W")==0) {
			op = checkOpt_Int(argc, argv, a, word_length, 0, INT_MAX);
		}
		else
		if (strcmp(argv[a], "-a")==0)
		{
			op = getOpt(argc, argv, a, tmpStr);
			if (op == 0)
			{
				error = errno;
				HSP_SK_ALPHA = strtod(tmpStr.c_str(), &stopstring); //this is my stuff, I need to verify the input
				if (errno != error) //error code set
				{
					cout << "Cannot recognize option -a " << argv[a] << "\n";
					op = -1;
				}
				else
					if (HSP_SK_ALPHA < 0 || HSP_SK_ALPHA > 1)
					{
						cout << "invalid value for option -a " << HSP_SK_ALPHA << ", it must be between 0 and 1" << "\n";
						op = -1;
					}
					else
						GENE_MS_BETA = 1 - HSP_SK_ALPHA;
			}
		}
		else
//		if (strcmp(argv[a], "-b")==0) //this option is not used anymore
//		{
//			op = getOpt(argc, argv, a, tmpStr);
//			if (op == 0)
//			{
//				error = errno;
//				GENE_MS_BETA = strtod(tmpStr.c_str(), &stopstring);
//				if (errno != error)
//				{
//					cout << "Cannot recognize option -b " << argv[a] << "\n";
//					op = -1;
//				}
//				else
//					if (GENE_MS_BETA < 0)
//					{
//						cout << "invalid value for option -b " << GENE_MS_BETA << ", it must be non-negative" << "\n";
//						op = -1;
//					}
//			}
//		}
//		else
		if (strcmp(argv[a], "-d")==0)
		{
			op = getOpt(argc, argv, a, tmpStr);
			if (op == 0)
			{
				error = errno;
				MAX_GAP_BTWN_HSP = strtol(tmpStr.c_str(), &stopstring, 10);
				if (errno != error)
				{
					cout << "Cannot recognize option -d " << argv[a] << "\n";
					op = -1;
				}
				else
					if (MAX_GAP_BTWN_HSP <= 0)
					{
						cout << "invalid value for option -d " << MAX_GAP_BTWN_HSP << ", it must be positive" << "\n";
						op = -1;
					}
			}
		}
		else
		if (strcmp(argv[a], "-r")==0)
		{
			op = getOpt(argc, argv, a, tmpStr);
			if (op == 0)
			{
				error = errno;
				TOP_RANK_NUM = strtol(tmpStr.c_str(), &stopstring, 10);
				if (errno != error)
				{
					cout << "Cannot recognize option -r " << argv[a] << "\n";
					op = -1;
				}
				else
					if (TOP_RANK_NUM <= 0) //make sure it's > 0
					{
						cout << "invalid value for option -r " << TOP_RANK_NUM << ", it must be positive" << "\n";
						op = -1;
					}
			}
		}
		else
		if (strcmp(argv[a], "-c")==0)
		{
			op = getOpt(argc, argv, a, tmpStr);
			if (op == 0)
			{
				error = errno;
				MIN_GENE_COVER_PER = strtod(tmpStr.c_str(), &stopstring);
				if (errno != error)
				{
					cout << "Cannot recognize option -c " << argv[a] << "\n";
					op = -1;
				}
				else
					if (MIN_GENE_COVER_PER < 0 || MIN_GENE_COVER_PER > 1)
					{
						cout << "invalid value for option -c " << MIN_GENE_COVER_PER << ", it must be between 0 and 1" << "\n";
						op = -1;
					}
			}
		}
		else
		if (strcmp(argv[a], "-s")==0)
		{
			op = getOpt(argc, argv, a, tmpStr);
			if (op == 0)
			{
				error = errno;
				MIN_GENE_SCORE = strtod(tmpStr.c_str(), &stopstring);
				if (errno != error)
				{
					cout << "Cannot recognize option -s " << argv[a] << "\n";
					op = -1;
				}
				else
					USE_MIN_GENE_SCORE = true;
			}
		}
		else
		if (strcmp(argv[a], "-scodon")==0)
			op = checkOpt_Int(argc, argv, a, START_CODON_CHECK_LEN, 0, 3000);
		else
		if (strcmp(argv[a], "-o")==0)
			op = getOpt(argc, argv, a, outStr); //output filename if specified
		else
		if (strcmp(argv[a], "-i")==0)
			op = checkOpt_Int(argc, argv, a, MIN_INTRON_LEN, 1, 300);
		else
		if (strcmp(argv[a], "-x")==0)
			op = checkOpt_Int(argc, argv, a, MIN_INTERNAL_EXON_LEN, 1, 300);
		else
		if (strcmp(argv[a], "-n")==0)
			op = checkOpt_Int(argc, argv, a, MAX_NUM_SPLICE_SITES, 1, 100);
		else
		if (strcmp(argv[a], "-v")==0)
			op = checkOpt_Int(argc, argv, a, SPLICE_SEGMENT_VERSION, 1, 2);
		else
		if (strcmp(argv[a], "-h")==0)
			op = checkOpt_Int(argc, argv, a, TREE_DATA_SHIFT, -20, 20);
		else
		if (strcmp(argv[a], "-j")==0)
			op = checkOpt_Float(argc, argv, a, TREE_CLS_THRESHOLD, -20, 20);
		else
		if (strcmp(argv[a], "-norepair")==0)
		{
			REPAIR_HSP_AFTER_EXON = 0;
			op = 0;
		}
		else
		if (strcmp(argv[a], "-b")==0)
		{
			VERBOSE = true;
			op = 0;
		}
		else
		if (strcmp(argv[a], "-pid")==0)
		{
			GENBLASTG_NEED_PID = true;
			op = 0;
		}
		else
		if (strcmp(argv[a], "-re")==0)
			op = checkOpt_Int(argc, argv, a, END_EXON_LEN, 1, 1000);
		else
		if (strcmp(argv[a], "-rm")==0)
			op = checkOpt_Int(argc, argv, a, MID_EXON_LEN, 1, 1000);
		else
		if (strcmp(argv[a], "-rl")==0)
			op = checkOpt_Int(argc, argv, a, EXPLORE_END_EXON_LEN, 100, 5000);
		else
		if (strcmp(argv[a], "-rp")==0) //percentage of padding length for repair (in terms of genBlastA group length)
			op = checkOpt_Float(argc, argv, a, EXPLORE_END_EXON_PER, 0.01, 1);
		else
		if (strcmp(argv[a], "-rs")==0)
			op = checkOpt_Int(argc, argv, a, REPAIR_HSP_MIN_INIT_SCORE, 1, 1000);
		else
		if (strcmp(argv[a], "-rd")==0)
			op = checkOpt_Int(argc, argv, a, REPAIR_HSP_EXTEND_SCORE_DROP, 1, 1000);
		else
		if (strcmp(argv[a], "-gff")==0)
		{
			OUTPUT_GFF = true;
			op = 0;
		}
		else
		if (strcmp(argv[a], "-cdna")==0)
		{
			OUTPUT_cDNA = true;
			op = 0;
		}
		else
		if (strcmp(argv[a], "-pro")==0)
		{
			OUTPUT_Protein = true;
			op = 0;
		}
		else
		if (strcmp(argv[a], "-id")==0)
		{
			op = getOpt(argc, argv, a, USER_ID);
		}
		else
		if (strcmp(argv[a], "-p")==0) //3 options here: genblasta, genblastg, genewise (hidden option)
		{
			string program_option;
			op = getOpt(argc, argv, a, program_option);
			if (program_option.compare("genblasta")==0)
				phase1_only = true;
			else
				if (program_option.compare("genewise")==0) //for experiment only, to compare with genewise (run genewise from within genblast codes, just for easier testing)
				{
					GENEWISE_COMMAND = true;
				}
				else
				if (program_option.compare("exonerate")==0) //for comparing with exonerate!
				{
					EXONERATE_COMMAND = true;
				}
				else
				{
					if (program_option.compare("genblastg") != 0) //phase1_only default is false
						op = -1;
				}
		}
/*		else
		if (strcmp(argv[a], "-gwgene") == 0)
		{
			string gw_gene_index_str;
			op = getOpt(argc, argv, a, gw_gene_index_str);
			genewise_gene_index = atoi(gw_gene_index_str.c_str());
		}
*/		else
		if (strcmp(argv[a], "-skip") == 0) //hidden option, for debug
		{
			SKIP = true;
			op = getOpt(argc, argv, a, skip_gene);
		}
		else
		{
			cout << "Cannot recognize option " << argv[a] << "\n";
			op = -1;
		}

		if (op == -1)
		{
			cout << SYNTAX_ERROR << "\n";
			return op;
		}

		a++;
	}

	//error checking: check whether we have stuff in query and target
	if (queryList.length() == 0 || targetDb.length() == 0)
	{
		cout << SYNTAX_ERROR << "\n";
		return -1;
	}

	char command[1024];

	//now ask wu-blastall to run
	string pathStr;
	char* gb_path = getenv("GBLAST_PATH");
	if (gb_path == NULL)
		pathStr = ".";
	else
		pathStr = gb_path;

	string xndfile, xnsfile, xntfile;
	if (prog_is_wublast)
	{
		xndfile = targetDb + ".xnd";
		xnsfile = targetDb + ".xns";
		xntfile = targetDb + ".xnt";
	}
	else
	{
		xndfile = targetDb + ".nhr";
		xnsfile = targetDb + ".nin";
		xntfile = targetDb + ".nsq";
	}

	if (!hasFile(xndfile.c_str()) || !hasFile (xnsfile.c_str()) || !hasFile (xntfile.c_str()))
	{		
		if (prog_is_wublast)
		{
			cout << "formatting target database for wublast..." << "\n";
			snprintf(command, 1024 ,"%s/xdformat -n %s", pathStr.c_str(), targetDb.c_str());
		}
		else
		{
			cout << "formatting target database for blast..." << "\n";
			snprintf(command, 1024 ,"%s/formatdb -i %s -p F", pathStr.c_str(), targetDb.c_str());
		}
		int i = system(command);
                if (i >=1)
                {
                  cerr << "error running" << command << "\n";
                  return -1;
                }
                     
	}
	if (!hasFile(xndfile.c_str()) || !hasFile (xnsfile.c_str()) || !hasFile (xntfile.c_str())) //test if xdformat is successful
	{
		cerr << XD_FILE_ERROR << "\n";
		return -1;
	}

	char wuFilename[256];
	char reportFilename[256];

	string targetDbFilename;
	int tgt_pos = targetDb.find_last_of("/\\");
	if (tgt_pos == string::npos)
		tgt_pos = -1;
	targetDbFilename = targetDb.substr(tgt_pos+1);
	if (prog_is_wublast)
		snprintf(wuFilename, 256 ,"%s_%s.wublast", queryList.c_str(), targetDbFilename.c_str());
	else
		snprintf(wuFilename, 256 ,"%s_%s.blast", queryList.c_str(), targetDbFilename.c_str());
	if (hasFile(wuFilename))
		cout << "blast already done, move forward" << "\n";
	else
	{
		if (prog_is_wublast)
		{
			cout << "running wublast..." <<"\n";
			if (proteinSeq.compare("T") == 0)
				snprintf(command, 1024 , "%s/wu-blastall -p tblastn -d %s -i %s -e %s -F %s -g %s -W %d -o %s warnings notes", pathStr.c_str(), 
					targetDb.c_str(), queryList.c_str(),
					eValue.c_str(), filterSet.c_str(), gapSet.c_str(), word_length, wuFilename);
			else
				snprintf(command, 1024 ,"%s/wu-blastall -p tblastx -d %s -i %s -e %s -F %s -g %s -W %d -o %s warnings notes", pathStr.c_str(), 
					targetDb.c_str(), queryList.c_str(),
					eValue.c_str(), filterSet.c_str(), gapSet.c_str(), word_length, wuFilename);
		}
		else
		{
			cout << "running blast..." << "\n";
			if (proteinSeq.compare("T") == 0)
				progStr = "tblastn";
			else
				progStr = "blastn";
			snprintf(command, 1024 ,"%s/blastall -p %s -d %s -i %s -e %s -F %s -g %s -W %d -o %s", pathStr.c_str(), progStr.c_str(), 
				targetDb.c_str(), queryList.c_str(), 
				eValue.c_str(), filterSet.c_str(), gapSet.c_str(), word_length, wuFilename);
		}
		int i = system(command);
                if (i >= 1)
                {
                   cerr << "error running " << command << "\n";
                   return -1;
                }
	}

	snprintf(reportFilename, 256 ,"%s.report", wuFilename);
	if (hasFile(reportFilename))
		cout << "blast.report already done, move forward" << "\n";
	else
	{
		//parse wublast output to get ".report" file
		ifstream wuFile;
		wuFile.open(wuFilename); //now this is my intermediate wublast file
		if (!wuFile.is_open())
		{
                        cerr << wuFilename << "\n";
			cerr << WU_FILE_ERROR << "\n";
			return -1;
		}
		if (prog_is_wublast)
			cout << "wublast done" << "\n";
		else
			cout << "blast done" << "\n";

		ofstream reportFile;
		reportFile.open(reportFilename);
		if (!reportFile.is_open())
		{
			cout << RPT_FILE_ERROR << "\n";
			return -1;
		}

	reportFile << "HSP_ID\tProtein_length\tHSP_Score\tE-value\tP-value\tP_Idn\tChr\tChr_start\tChr_end\tChr_length\tchr_strand\tquery_start\tquery_end\tquery_length\n";
	string wuline, chromosome;
	int cur_pos, cur_space_pos, hsp_count=1, query_len, chr_start, chr_end, query_start, query_end, exit_code;
	bool strandness, lastEnded = true, wublast_ok = true, cur_output = true;
	string cur_query_name;
	while (!wuFile.eof() && wublast_ok)
	{
		getline(wuFile, wuline);
		remove_trailing_cr_lf(wuline);

		if (wuline.find("Query=") == 0) //found query name
		{
			if (!cur_output && hsp_count > 1) //has not done previous output
			{
				reportFile << chr_start << "\t" << chr_end << "\t" << chr_end-chr_start+1 << "\t";
				if (strandness)
					reportFile << "1\t";
				else
					reportFile << "-1\t";
				reportFile << query_start << "\t" << query_end << "\t" << query_end-query_start+1 << "\n";
				cur_output = true;
			}
			//reset
			hsp_count = 1; 
			lastEnded = true;

			int fst_space_pos = wuline.find_first_not_of(" \t", 6);
			if (fst_space_pos != string::npos)
			{
				int snd_space_pos = wuline.find_first_of(" \t", fst_space_pos);
				if (snd_space_pos == string::npos)
					snd_space_pos = wuline.length();
				wuline.erase(snd_space_pos); //remove trailing part after the second white space
				cur_query_name = wuline.substr(fst_space_pos);
				reportFile << ">" << cur_query_name << "\n"; //get stuff after "Query=  " (8 chars including two trailing spaces)
			}
			else
			{
				cout << "check wuBlast/Blast file: Query= is not followed by query name" << "\n";
				exit(-1);
			}
			//if (cur_output || hsp_count == 1)
				cur_output = false; //reset to false at the beginning of every query
			//else
			//{
			//	cout << "previous query before " << cur_query_name << " not done properly, check wublast file" << "\n";
			//	exit(-1);
			//}

			//UPDATE:now read additional lines to get query protein length (assume it's always in the next line following "Query=")
			//read until hit "(* letters" in a line, * is some integer number
			bool query_len_found = false;
			while (!wuFile.eof() && wublast_ok && !query_len_found)
			{
				getline(wuFile, wuline);
				remove_trailing_cr_lf(wuline);
				if ((cur_pos = wuline.find("(")) != string::npos && (cur_space_pos = wuline.find(" letters")) != string::npos)
				{
					if (cur_space_pos > cur_pos)
					{
						error = errno;
						query_len = strtol(wuline.substr(cur_pos+1, cur_space_pos - cur_pos - 1).c_str(), &stopstring, 10);
						if (errno == error) //there is no error produced by strtol function
							query_len_found = true; //we are good to go!
					}
				}
			}
			if (!query_len_found)
			{
				cout << "wublast format has problem (after getting query name " << cur_query_name << "), check wublast file" << "\n";
				exit(-1);
			}
		}
		else
			/*if (wuline.find("record") != string::npos)//UPDATE:if there's only 1 query, wublast has no "record"!
			{
				cur_pos = wuline.find("(");
				cur_space_pos = wuline.find(" ", cur_pos);
				query_len = atoi(wuline.substr(cur_pos+1, cur_space_pos - cur_pos - 1).c_str());
			}
			else
			*/if (wuline.find(">") == 0) //format is: ">id description"
			{
				int desc_pos = wuline.find(" "); //find first space
				if (desc_pos == string::npos) //no description after id
					desc_pos = wuline.length();
				chromosome = wuline.substr(1, desc_pos-1);
			}
			else
			if (wuline.find("Score = ") != string::npos)
			{
					if (!lastEnded)
					{
						reportFile << chr_start << "\t" << chr_end << "\t" << chr_end-chr_start+1 << "\t";
						if (strandness)
							reportFile << "1\t";
						else
							reportFile << "-1\t";
						reportFile << query_start << "\t" << query_end << "\t" << query_end-query_start+1 << "\n";
					}

				query_start = INT_MAX; //reset
				query_end = 0;
				chr_start = INT_MAX;
				chr_end = 0;

				int pre_space_pos = wuline.find_first_not_of(" \t", 9); //hsp score
				cur_space_pos = wuline.find_first_of(" \t", pre_space_pos);
				if (cur_space_pos == string::npos)
					cur_space_pos = wuline.length();
				reportFile << hsp_count << "\t" << query_len << "\t" << wuline.substr(pre_space_pos, cur_space_pos - pre_space_pos) << "\t";
				
				//cout << "here1\n";
				hsp_count++; //increment count
				cur_pos = wuline.find("=", cur_space_pos); 
				if (cur_pos == string::npos)
				{
					cout << "wuBlast/Blast file format error" << "\n";
					exit(-1);
				}
				cur_pos = wuline.find_first_not_of(" \t", cur_pos+1);
				cur_space_pos = wuline.find_first_of(" \t,", cur_pos); //evalue
				if (cur_space_pos == string::npos)
					cur_space_pos = wuline.length();
				reportFile << wuline.substr(cur_pos, cur_space_pos - cur_pos) << "\t";

				//cout << "here2\n";
				cur_pos = wuline.find("=", cur_space_pos); //pvalue
				if (cur_pos != string::npos)
					reportFile << wuline.substr(cur_pos+2) << "\t";
				else
					reportFile << "unknown\t";

				//cout << "here3\n";
				lastEnded = false;
			}
			else
			/*if (wuline.find("EXIT CODE") == 0) //this is only for last HSP in each query //UPDATE:if there's only 1 query, there's no "EXIT CODE"!!!
			{
				error = errno;
				exit_code = strtol(wuline.substr(10).c_str(), &stopstring, 10);
				if (errno != error)
				{					
					wublast_ok = false;
				}
				else
				{
					if (exit_code != 0)
						wublast_ok = false;*/
					/*else //make sure exit code is 0
					{
						//if (hsp_count > 1) //if hsp_count > 1, means there's at least 1 HSP, otherwise nothing is outputted yet
						if (!cur_output && hsp_count > 1)
						{							
						reportFile << chr_start << "\t" << chr_end << "\t" << chr_end-chr_start+1 << "\t";
						if (strandness)
							reportFile << "1\t";
						else
							reportFile << "-1\t";
						reportFile << query_start << "\t" << query_end << "\t" << query_end-query_start+1 << "\n";
						cur_output = true;
						}
					}*/
				/*}
			}
			else
			if (wuline.find("Statistics:")==0) //UPDATE: for case when there's only 1 query
			{
				if (!cur_output && hsp_count > 1) //has not see the end of current output
				{
					reportFile << chr_start << "\t" << chr_end << "\t" << chr_end-chr_start+1 << "\t";
					if (strandness)
						reportFile << "1\t";
					else
						reportFile << "-1\t";
					reportFile << query_start << "\t" << query_end << "\t" << query_end-query_start+1 << "\n";
					cur_output = true;
				}
			}
			else*/
			if (wuline.find("Identities") != string::npos)
			{
				int identity, totalchar;
				cur_pos = wuline.find("/");
				identity = atoi(wuline.substr(14, cur_pos - 14).c_str());

				//cout << "here4\n";
				cur_space_pos = wuline.find_first_of(" \t,(", cur_pos);
				totalchar = atoi(wuline.substr(cur_pos+1, cur_space_pos - cur_pos - 1).c_str());

				//cout << "here5\n";
				double pid = (double)identity*100 / (double)totalchar; //pid
				reportFile << pid << "\t" << chromosome << "\t"; //chromosome
				if ((cur_pos = wuline.find("Frame = ")) != string::npos) //for wublast, Frame and Identities are on the same line
					if (wuline.substr(cur_pos+8, 1).compare("+") == 0) //strandness
						strandness = true;
					else
						strandness = false;
				//cout << "here6\n";
			}
			else
			if ((cur_pos = wuline.find("Frame = ")) != string::npos) //for blast
			{
				if (wuline.substr(cur_pos+8, 1).compare("+") == 0) //strandness
					strandness = true;
				else
					strandness = false;
			}
			else
			if (wuline.find("Query: ") != string::npos)
			{
				cur_pos = wuline.find_first_of("0123456789");
				cur_space_pos = wuline.find(" ", cur_pos);
				int s = atoi(wuline.substr(cur_pos, cur_space_pos - cur_pos).c_str());
				
				//cout << "here7\n";
				if (query_start > s)
					query_start = s;
				cur_pos = wuline.find_first_of("0123456789", cur_space_pos);
				int e = atoi(wuline.substr(cur_pos).c_str());

				//cout << "here8\n";
				if (query_end < e)
					query_end = e;
			}
			else
			if (wuline.find("Sbjct:") != string::npos)
			{
				cur_pos = wuline.find_first_of("0123456789");
				cur_space_pos = wuline.find(" ", cur_pos);
				int s = atoi(wuline.substr(cur_pos, cur_space_pos - cur_pos).c_str());

				//cout << "here9\n";
				cur_pos = wuline.find_first_of("0123456789", cur_space_pos);
				int e = atoi(wuline.substr(cur_pos).c_str());

				//cout << "here10\n";
				if (strandness)
				{
					if (chr_start > s)
						chr_start = s;
					if (chr_end < e)
						chr_end = e;
				}
				else
				{
					if (chr_end < s)
						chr_end = s;
					if (chr_start > e)
						chr_start = e;
				}
			}
	}
	//output last gene
	if (hsp_count > 1) //if hsp_count > 1, means there's at least 1 HSP, otherwise nothing is outputted yet
	{							
		reportFile << chr_start << "\t" << chr_end << "\t" << chr_end-chr_start+1 << "\t";
		if (strandness)
			reportFile << "1\t";
		else
			reportFile << "-1\t";
		reportFile << query_start << "\t" << query_end << "\t" << query_end-query_start+1 << "\n";
	}

	wuFile.close();
	reportFile.close();

	if (!wublast_ok)
	{
		cout << "unexpected blast report, please check " << wuFilename << "\n";
		return -1;
	}
	}

	//now do genBlast first phase
	string filenames[3];
	filenames[0] = reportFilename; //"seqList.txt.wublast.report"; //the report filename
	filenames[1] = targetDb; //target sequence file
	filenames[2] = wuFilename;//the wublast filename
	QUERY_SEQ_FILE = queryList; //the query sequence file
	MIN_INTRON_LEN_AA = MIN_INTRON_LEN / 3 + 1;

//	int gw_gene_count = 1; //damn! have to run it one by one due to segmentation fault after each gene?!

	//string outStr = "test.out";
	DataManager	dm(filenames, outStr.c_str(), phase1_only); //phase1_only: true for genBlastA, false is for genBlastG
	while (!dm.inputFile_Finish)
	{
		if (SKIP)
		{
			if (dm.next_query_gene.compare(skip_gene) != 0)
			{
				if (!dm.ReadFile_Skip())
					return -2;
				if (!phase1_only)
					dm.GetAlignments_Skip();
				dm.Reset();
				continue;
			}
			else
				SKIP = false;
		}

		if (!dm.ReadFile())
		{
			cout << FILE_READ_ERROR << filenames[0] << "\n";
			return -1;
		}

#ifdef PERFORMANCE
			of_perform << dm.query_gene << "\t";
			num_of_HSPs_for_query = num_of_edges_for_query = num_of_groups_pruned_due_to_overlap = 0;
#endif
		//now start phase 1
		clock_t start_time = clock();

		int hsp_chr_size = dm.HSP_chr.size();
		for (int i=0; i<hsp_chr_size; i++)
		{
			HSP_graph	hsp_graph(dm, true, i);
			hsp_graph.GeneGrouping();

			HSP_graph	hsp_neg_graph(dm, false, i);
			hsp_neg_graph.GeneGrouping();
		}

		dm.PrintGroups(false); //print text only
		//dm.PrintGroups(true); //print overview after text

		//record running time (phase 1 finished)
		double duration = (double)(clock() - start_time) / CLOCKS_PER_SEC;
		dm.outFile << "\n" << "genBlast phase 1 running time: " << duration << " seconds" << "\n";
#ifdef PERFORMANCE
			of_perform << num_of_HSPs_for_query << "\t" << num_of_edges_for_query << "\t" << duration;

			multimap<Group_Info, vector<HSP_Gene_Pair*> >::iterator mapIt = dm.groups.begin();

			if (mapIt != dm.groups.end())

				of_perform << "\t" << (*mapIt).first.score; //score of rank 1 group

			of_perform << "\t" << num_of_groups_pruned_due_to_overlap;
#endif


		//now start phase 2
		if (!phase1_only)
		{
		duration_gblastg_without_comp_final_align = 0;
		//dm.GetChromosomes();
		gblast_start_time = clock();
		dm.PrepareOutputHSPs(); //this function collects all useful HSP alignments

#ifdef TIMING
		dm.outFile << "after PrepareOutputHSPs: at " << (double)(clock() - gblast_start_time) / CLOCKS_PER_SEC << " seconds" << "\n";
#endif

		ACCP_DONR_graph ad_graph(dm);
		//execute phase2-version-2.3
		dm.LoadQuerySeqData(QUERY_SEQ_FILE.c_str());

#ifdef TIMING
		dm.outFile << "after LoadQuerySeqData: at " << (double)(clock() - gblast_start_time) / CLOCKS_PER_SEC << " seconds" << "\n";
#endif

		//if (SPLICE_SEGMENT_VERSION == 2)
		LoadAlignScore();
		LoadCodonMap_SingleLetterAA();
		//if (REPAIR_HSP_AFTER_EXON == 1)
			ad_graph.GenePredict6();
		//else
		//	ad_graph.GenePredict4();

		duration = (double)(clock() - gblast_start_time) / CLOCKS_PER_SEC;
		//dm.outFile << "\n" << "phase 2 running time: " << duration << " seconds" << "\n";
		dm.outFile << "\n" << "phase 2 running time: " << duration << " seconds\n";
			//phase 2 time without computing final pid: " << duration_gblastg_without_comp_final_align << " seconds" << "\n";

#ifdef PERFORMANCE
				of_perform << "\t" << num_of_edges << "\t" << duration << "\t" << duration_gblastg_without_comp_final_align;// << "\n";
				of_perform << "\t" << num_of_correct_exons << "\t" << num_of_incorrect_exons << "\t" << num_of_missed_exons;// << "\n";

		if (GENBLASTG_NEED_PID)
				of_perform << "\tGBG\t" << max_gblastg_final_align_pid;
#endif //#ifdef PERFORMANCE

		}

		cout << dm.query_gene << ": finished" << "\n";
#ifdef PERFORMANCE
			of_perform << "\n";
#endif

		dm.Reset();

	} //end while (for each query_gene in combined input report file)


#ifdef PERFORMANCE
	of_perform << "Total_WormBase_Exons: " << Total_WormBase_Exons << "\n";

	of_perform << "Total_WormBase_Exons_Effective: " << Total_WormBase_Exons_Effective << "\n";

	of_perform << "Total_GenBlast_Exons: " << Total_GenBlast_Exons << "\n";

	of_perform << "Total_Matching_Exons: " << Total_Matching_Exons << "\n";

	of_perform << "Total_Matching_Genes: " << Total_Matching_Genes << "\n";

	of_perform.close();
#endif

//	_CrtDumpMemoryLeaks();
	return 0;
}

