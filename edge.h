//***************EDGE.H************************

//

//Classes used by graphs (edge, node, penalty, distance etc.)

//

//

//Author: Rong She

//Date: April 2007

//***********************************************


#if defined(_MSC_VER)
	#pragma warning(disable: 4786)
#endif


#ifndef EDGE_H /* EDGE_H */

#define EDGE_H



#include "scores.h"



#include <string.h> //why?



#include <stack>

using namespace std;



//******************************************************************//

class Base_Edge

{

public:

	int						source_num;

	int						dest_num;

//	float					penalty;



	bool					isCut; //used for checking monotonity of paths when constructing skip edges



/*	Edge(int source, int dest, float p, bool isC): isCut(isC) 

	{

		source_num = source;

		dest_num = dest;

		penalty = p; 

	}

*/

    bool operator<(const Base_Edge& h) const

    { return (source_num < h.source_num) || ((source_num == h.source_num) && (dest_num < h.dest_num)); }



    //bool operator<(const Base_Edge& h) const //reverse order

    //{ return source_num > h.source_num || (source_num == h.source_num && dest_num > h.dest_num); }



	friend ostream& operator<<(ostream& os, const Base_Edge& h)

	{

		os << "edge:("<< h.source_num << "-" << h.dest_num  << "); "//penalty: " << h.penalty 

			<< "; isCut? "<< h.isCut << "\n";

		return os;

	}

};



class AD_Edge_Penalty

{

public:

	int						exact_match;

	int						gap;

	int						pos_match;



	set<pair<int, int> >	query_match_pos; //this is exact-matched regions (differ from AD_Edge_QueryCover)

	//vector<int>				query_match_pos;



	AD_Edge_Penalty()

	{

		exact_match = 0;

		gap = 0; //INT_MAX;

		pos_match = 0;

	}



	AD_Edge_Penalty(int e, int g, int p)

	{

		exact_match = e;

		gap = g;

		pos_match = p;

		query_match_pos.clear();

	}



	AD_Edge_Penalty(int e, int g, int p, //vector<int>& query_matches)

		set<pair<int, int> >& query_matches)

	{

		exact_match = e;

		gap = g;

		pos_match = p;

		//copy(query_matches.begin(), query_matches.end(), query_match_pos.begin());

		set< pair<int,int> >::iterator it = query_matches.begin();

		query_match_pos.clear();

		//vector<int>::iterator it = query_matches.begin();

		for (; it != query_matches.end(); it++)

			query_match_pos.insert( *it );

			//query_match_pos.push_back(*it);

	}



	AD_Edge_Penalty(const AD_Edge_Penalty& p) //copy constructor

	{

		exact_match = p.exact_match;

		gap = p.gap;

		pos_match = p.pos_match;

		set< pair<int,int> >::const_iterator it = p.query_match_pos.begin();

		query_match_pos.clear();

		//vector<int>::const_iterator it=p.query_match_pos.begin();

		for (; it != p.query_match_pos.end(); it++)

			query_match_pos.insert( *it );

			//query_match_pos.push_back(*it);

	}

	

	AD_Edge_Penalty operator=(const AD_Edge_Penalty& p)

	{

		exact_match = p.exact_match;

		gap = p.gap;

		pos_match = p.pos_match;

		set< pair<int,int> >::const_iterator it = p.query_match_pos.begin();

		query_match_pos.clear();

		for (; it != p.query_match_pos.end(); it++)

			query_match_pos.insert( *it );

		return *this;

	}



	bool IsSame(AD_Edge_Penalty& p)

	{

		return (exact_match == p.exact_match && gap == p.gap && pos_match == p.pos_match);

	}



	bool IsSame(int e, int g, int p)

	{

		return (exact_match == e && gap == g && pos_match == p);

	}



	void Reset() //only need to reset exact_match (used only in ComputePenalty() of edge, to skip creation of edges without matches)

	{

		exact_match = 0;

		//gap = 0; //INT_MAX;

		//pos_match = 0;

		query_match_pos.clear();

	}



	void ClearSpace()

	{

		query_match_pos.clear();

	}



	//for exon edges, if no exact_match, or (updated to reflect comparison criteria) if exact_match is not worth it

	bool NoExactMatch()

	{

		return exact_match == 0 || exact_match <= (gap+pos_match);

	}



	int TotalPos() const

	{

		return exact_match + gap + pos_match;

	}



	float PID() const

	{

		return (float)exact_match / TotalPos();

	}





	AD_Edge_Penalty operator+(const AD_Edge_Penalty& p);



	//used only for penalty comparison within single HSP area, since all penalties must not have overlapping query matches

	AD_Edge_Penalty SimpleAddition(const AD_Edge_Penalty& p)

	{

		AD_Edge_Penalty tmp(exact_match, gap, pos_match);//first, copy the first set into "tmp"

		tmp.exact_match += p.exact_match;

		tmp.gap += p.gap;

		tmp.pos_match += p.pos_match;

		return tmp;

	}



	int Score() const {return exact_match - gap - pos_match;}



	//comparison criteria

	bool operator<(const AD_Edge_Penalty& p) const

	{

//		return exact_match > p.exact_match || (exact_match == p.exact_match && gap < p.gap) 

//			|| (exact_match == p.exact_match && gap == p.gap && pos_match > p.pos_match);

//		return float(exact_match ) / float(exact_match+gap) > float(p.exact_match ) / float(p.exact_match+p.gap);



		//return exact_match+pos_match-gap > p.exact_match+p.pos_match-p.gap;

		//return (float)exact_match + (float)pos_match/2.0 - (float)gap > (float)p.exact_match + (float)p.pos_match/2.0 - (float)p.gap;

		return 

//			(exact_match - gap - pos_match > p.exact_match - p.gap - p.pos_match) ||

//			( (exact_match - gap - pos_match == p.exact_match - p.gap - p.pos_match) && (exact_match > p.exact_match) );



			//(exact_match - gap > p.exact_match - p.gap) || 

			//( (exact_match-gap == p.exact_match-p.gap) && (exact_match < p.exact_match) );

			(Score() > p.Score() ) || 

			( (Score() == p.Score()) && (exact_match < p.exact_match) );



//			( (exact_match-gap == p.exact_match-p.gap) && (pos_match < p.pos_match) ) || 

//			( (exact_match-gap == p.exact_match-p.gap) && (pos_match == p.pos_match) 

//			  && (exact_match > p.exact_match) );



			//( (exact_match-gap == p.exact_match-p.gap) && (exact_match > p.exact_match) ) || 

			//( (gap == p.gap) && (exact_match == p.exact_match) && (pos_match > p.pos_match) );



	}



	bool operator<=(const AD_Edge_Penalty& p) const

	{

		//return (exact_match+pos_match- gap >= p.exact_match + p.pos_match - p.gap); // || (exact_match == p.exact_match && pos_match == p.pos_match && gap == p.gap);

		//return (float)exact_match + (float)pos_match/2.0 - (float)gap >= (float)p.exact_match + (float)p.pos_match/2.0 - (float)p.gap;

		return (*this < p || (Score() == p.Score() && exact_match == p.exact_match)) ;

	}



	friend ostream& operator<<(ostream& os, const AD_Edge_Penalty& p)

	{

		os << "exact match:" << p.exact_match << "; gap:" << p.gap << "; pos match:" << p.pos_match;

		os << "; exact matches: ";

		set< pair<int, int> >::const_iterator  it = p.query_match_pos.begin();

		for (; it != p.query_match_pos.end(); it++)

			os << "[" << (*it).first << "-" << (*it).second << "); " ;

		os << "\n";

		return os;

	}

};



class QuerySeg_Start_End_PID

{

public:

	int first; //start position of segment

	int second; //ending position of segment

	float pid; //pid of this segment



	QuerySeg_Start_End_PID(int start, int end, float p)

	{

		first = start;

		second = end;

		pid = p;

	}



	bool operator<(const QuerySeg_Start_End_PID& q) const

	{

		return (first < q.first || (first == q.first && second < q.second));

	}

};





class AD_Edge_QueryCover

{

//1        10   8     14

//1234567890    8901234

//AAAAAAAAAA----AAAAAAA (query)

//|| || ||||    ||| |||

//BBBBB-BBBBBBBBBBBBBBB (target)

//12345 678901234567890

//          10        20

//(HSP[1]:pid1) (HSP[2]:pid2) (pid2 > pid1)

	//modified: use 3* query position (e.g. [1,10] is [1,30], query length 100 => 300 etc.)

public:

	int		query_cover; //number of query base pairs covered (14)

	//set<pair<int, int> > query_segments; //the actual segments of query covered (not necessarily matched regions)

	set<QuerySeg_Start_End_PID> query_segments; //[1,7,pid1], [8,14,pid2]

	float	qc_score; //7*pid1 + 7*pid2, computed according to query_segments

	int		target_len; //20 (total length, not counting gaps, i.e. dest-source+1)

	set<pair<int, int> > target_segments; //[1,20], used to compute target_len

	int		gap_len; //the length on target sequence that does not have any query correspondance

	//set<pair<int, int> > target_fragments; //[1,9], [14,20]

	bool	gap_at_end; //for operator+()



	AD_Edge_QueryCover()

	{

		query_cover = 0;

		qc_score = 0.0;

		target_len = 0;

		gap_len = 0;

		gap_at_end = true;

	}



	AD_Edge_QueryCover(int cover, set<QuerySeg_Start_End_PID>& query_matches, float score, 

		int t_len, set<pair<int, int> >& t_segments, int g_len=0, bool g_end=true) //set<pair<int, int> >& query_matches)

	{

		query_cover = cover;

		//set< pair<int,int> >::iterator it = query_matches.begin();

		set<QuerySeg_Start_End_PID>::iterator it = query_matches.begin();

		for (; it != query_matches.end(); it++)

			query_segments.insert( *it );

		qc_score = score;

		target_len = t_len;

		set< pair<int,int> >::iterator target_it = t_segments.begin();

		for (; target_it != t_segments.end(); target_it++)

			target_segments.insert( *target_it);

		gap_len = g_len;

		gap_at_end = g_end;

	}



	AD_Edge_QueryCover(const AD_Edge_QueryCover& qc) //copy constructor

	{

		query_cover = qc.query_cover;

		//set< pair<int,int> >::const_iterator it = qc.query_segments.begin();

		set<QuerySeg_Start_End_PID>::const_iterator it = qc.query_segments.begin();

		for (; it != qc.query_segments.end(); it++)

			query_segments.insert( *it );

		qc_score = qc.qc_score;

		target_len = qc.target_len;

		set< pair<int,int> >::const_iterator target_it = qc.target_segments.begin();

		for (; target_it != qc.target_segments.end(); target_it++)

			target_segments.insert( *target_it);

		gap_len = qc.gap_len;

		gap_at_end = qc.gap_at_end;

	}



/*	AD_Edge_QueryCover operator=(AD_Edge_QueryCover& qc) //copy constructor

	{

		query_cover = qc.query_cover;

		set< pair<int,int> >::iterator it = qc.query_segments.begin();

		for (; it != qc.query_segments.end(); it++)

			query_segments.insert( *it );



		return *this;

	}

*/

	AD_Edge_QueryCover operator+(const AD_Edge_QueryCover& qc);



	//comparison criteria

	bool operator<(const AD_Edge_QueryCover& qc) const

	{

		return query_cover < qc.query_cover;

	}



	bool operator<=(const AD_Edge_QueryCover& qc) const

	{

		return query_cover <= qc.query_cover;

	}



	friend ostream& operator<<(ostream& os, const AD_Edge_QueryCover& p)

	{

		os << "query cover: " << p.query_cover;

		//set<pair<int,int> >::const_iterator it=p.query_segments.begin();

		set<QuerySeg_Start_End_PID>::const_iterator it = p.query_segments.begin();

		for (; it != p.query_segments.end(); it++)

			os << "[" << (*it).first << "-" << (*it).second << "], pid=" << (*it).pid;

		os << "; score: " << p.qc_score;

		os << "\n";

		set<pair<int, int> >::const_iterator t_it = p.target_segments.begin();

		for (; t_it != p.target_segments.end(); t_it++)

			os << "[" << (*t_it).first << "-" << (*t_it).second << "], ";

		os << "t_len=" << p.target_len << "; gap_len=" << p.gap_len << "(gap_at_end:" << p.gap_at_end << ")" << "\n";

		return os;

	}

};



class AD_Edge_Score

{

public:

	AD_Edge_QueryCover	query_coverage;

	AD_Edge_Penalty		pid_penalty;

	//float				hsp_pid;



	AD_Edge_Score() 

	{ 

		//hsp_pid = 0.0;

	}



	AD_Edge_Score(AD_Edge_QueryCover& qc, AD_Edge_Penalty& p)

	{

		query_coverage = qc;

		pid_penalty = p;

		//hsp_pid = pid;

	}



	AD_Edge_Score(const AD_Edge_Score& s): query_coverage(s.query_coverage), pid_penalty(s.pid_penalty) { }



/*	void Set_Intron_Edge_Score(int source, int dest)

	{

		query_coverage.target_len = dest - source + 1;

		query_coverage.target_segments.insert(pair<int, int>(source, dest));

	}

*/	

	//"<" means better!!!

	bool operator<(const AD_Edge_Score& edge_score) const //scoring function (compare based on number of base pairs...)

	{

//		return pid_penalty.Score() + 0.5 * query_coverage.query_cover * 3 > 

//			edge_score.pid_penalty.Score() + 0.5 * edge_score.query_coverage.query_cover * 3;

//		return pid_penalty < 

//			edge_score.pid_penalty;

//		return (float) query_coverage.query_cover * pid_penalty.PID() > 

//			(float) edge_score.query_coverage.query_cover * edge_score.pid_penalty.PID();

//		return 3 * query_coverage.qc_score + pid_penalty.Score() > 

//			edge_score.query_coverage.qc_score + edge_score.pid_penalty.Score();

//		return query_coverage.query_cover*3 - (Base_Len() - query_coverage.query_cover*3) + pid_penalty.Score() > 

//			edge_score.query_coverage.query_cover*3 - (edge_score.Base_Len() - edge_score.query_coverage.query_cover*3) + edge_score.pid_penalty.Score();

//		return query_coverage.query_cover + pid_penalty.Score() >

//			edge_score.query_coverage.query_cover + edge_score.pid_penalty.Score();

		return query_coverage.qc_score + pid_penalty.Score() >

			edge_score.query_coverage.qc_score + edge_score.pid_penalty.Score();

//		if (pid_penalty < edge_score.pid_penalty)

//			return true;

//		else

//			//if (query_coverage.qc_score > edge_score.query_coverage.qc_score)

//			if (query_coverage.query_cover * pid_penalty.PID() > edge_score.query_coverage.query_cover * edge_score.pid_penalty.PID())

//				return true;

//			else

//				return false;

//		return query_coverage.query_cover*3 - (Base_Len() - query_coverage.query_cover*3) >

//			edge_score.query_coverage.query_cover*3 - (edge_score.Base_Len() - edge_score.query_coverage.query_cover*3);

//		return (float) query_coverage.query_cover / Base_Len() > 

//			(float) edge_score.query_coverage.query_cover / edge_score.Base_Len();

//		return query_coverage.query_cover*3 - (Base_Len() - query_coverage.query_cover*3) + pid_penalty.Score() >

//			edge_score.query_coverage.query_cover*3 - (edge_score.Base_Len() - edge_score.query_coverage.query_cover*3)

//			+ edge_score.pid_penalty.Score();

//		return query_coverage.query_cover >

//			edge_score.query_coverage.query_cover;

//		return query_coverage.qc_score - query_coverage.gap_len >

//			edge_score.query_coverage.qc_score - edge_score.query_coverage.gap_len;

//		return query_coverage.qc_score - query_coverage.gap_len * pid_penalty.PID() >

//			edge_score.query_coverage.qc_score - edge_score.query_coverage.gap_len * edge_score.pid_penalty.PID();

	}



	bool operator<=(const AD_Edge_Score& edge_score) const

	{

//		return pid_penalty.Score() + 0.5 * query_coverage.query_cover * 3 >= 

//			edge_score.pid_penalty.Score() + 0.5 * edge_score.query_coverage.query_cover * 3;

//		return pid_penalty <= 

//			edge_score.pid_penalty;

//		return (float) query_coverage.query_cover * pid_penalty.PID() >= 

//			(float) edge_score.query_coverage.query_cover * edge_score.pid_penalty.PID();

//		return query_coverage.qc_score >= edge_score.query_coverage.qc_score;

//		return 3 * query_coverage.qc_score + pid_penalty.Score() >= 

//			edge_score.query_coverage.qc_score + edge_score.pid_penalty.Score();

//		return query_coverage.query_cover*3 - (Base_Len() - query_coverage.query_cover*3) + pid_penalty.Score() >=

//			edge_score.query_coverage.query_cover*3 - (edge_score.Base_Len() - edge_score.query_coverage.query_cover*3)+ edge_score.pid_penalty.Score();

		return query_coverage.qc_score + pid_penalty.Score() >=

			edge_score.query_coverage.qc_score + edge_score.pid_penalty.Score();

/*		if (pid_penalty <= edge_score.pid_penalty)

			return true;

		else

			if (pid_penalty.Score() + query_coverage.query_cover * pid_penalty.PID() >= edge_score.query_coverage.query_cover * edge_score.pid_penalty.PID())

				return true;

			else

				return false;

*/

//		return query_coverage.query_cover*3 - (Base_Len() - query_coverage.query_cover*3) >=

//			edge_score.query_coverage.query_cover*3 - (edge_score.Base_Len() - edge_score.query_coverage.query_cover*3);

//		return (float) query_coverage.query_cover / Base_Len() >= 

//			(float) edge_score.query_coverage.query_cover / edge_score.Base_Len();

//		return query_coverage.query_cover*3 - (Base_Len() - query_coverage.query_cover*3) + pid_penalty.Score() >=

//			edge_score.query_coverage.query_cover*3 - (edge_score.Base_Len() - edge_score.query_coverage.query_cover*3)

//			+ edge_score.pid_penalty.Score();

//		return query_coverage.query_cover >=

//			edge_score.query_coverage.query_cover;

//		return query_coverage.qc_score - query_coverage.gap_len >=

//			edge_score.query_coverage.qc_score - edge_score.query_coverage.gap_len;

//		return query_coverage.qc_score - query_coverage.gap_len * pid_penalty.PID() >=

//			edge_score.query_coverage.qc_score - edge_score.query_coverage.gap_len * edge_score.pid_penalty.PID();

	}



	int Base_Len() const

	{

		//return pid_penalty.TotalPos() > query_coverage.query_cover*3 ? pid_penalty.TotalPos() : query_coverage.query_cover*3;

		return query_coverage.target_len > query_coverage.query_cover*3 ? query_coverage.target_len : query_coverage.query_cover*3;

	}



	AD_Edge_Score operator+(AD_Edge_Score& score)

	{

		AD_Edge_Score tmp;

		tmp.query_coverage = query_coverage + score.query_coverage;

		tmp.pid_penalty = pid_penalty + score.pid_penalty;

		//tmp.hsp_pid = (hsp_pid + score.hsp_pid) / 2;



		return tmp;

	}



	friend ostream& operator<<(ostream& os, const AD_Edge_Score& s)

	{

		os << s.pid_penalty;

		os << s.query_coverage;

		//os << s.hsp_pid;

		return os;

	}

};



class ACCP_DONR_Edge_w_Score : public Base_Edge //used by version 2.3

{

public:

	//score contains two parts: number of query base pairs that are non-gap, number of exact match-gap-mismatch

	AD_Edge_Score		score;



	//record frame position of the exon

	int		tail; //trailing part

	int		head; //leading part



	//for checking possible inter-exon stop codon

	char	start_bp[3];

	char	end_bp[3];



	ACCP_DONR_Edge_w_Score() {}



	ACCP_DONR_Edge_w_Score(int source, int dest, bool isC, AD_Edge_Score& s, int h, int t, const char s_bp[3], const char e_bp[3])

	{

		source_num = source;

		dest_num = dest;

		isCut = isC;

		score = s;

		head = h;

		tail = t;

		strcpy(start_bp, s_bp);

		strcpy(end_bp, e_bp);

	}



	ACCP_DONR_Edge_w_Score(const ACCP_DONR_Edge_w_Score& s): score(s.score) 

	{

		source_num = s.source_num;

		dest_num = s.dest_num;

		isCut = s.isCut;

		tail = s.tail;

		head = s.head;

		strcpy(start_bp, s.start_bp);

		strcpy(end_bp, s.end_bp);

	}



	friend ostream& operator<<(ostream& os, const ACCP_DONR_Edge_w_Score& e)

	{

		os << "edge[" << e.source_num << " .. " << e.dest_num << "]:isCut?" << e.isCut << "\nscore: " << e.score;

		os << "head: " << e.head << "(" << e.start_bp << "); tail: " << e.tail << "(" << e.end_bp << ")" << "\n";

		return os;

	}

};



class ACCP_DONR_Edge : public Base_Edge

{

public:

//	float					penalty[3]; //one penalty for each frame

	AD_Edge_Penalty			penalty[3];



	//used only for exons (record the number of nucleotides in a possibly 'broken' codon, must be [0,3)

	char						leftover[3]; //start positions corresponding to each penalty

	char						remainder[3]; //the rest of nucleotides remained from current exon

	

	//used only for introns

	char					start_bp[3]; //record the last two chars (base pair) in previous exon (last space for '\0')

	char					end_bp[3]; //record the first two chars in next exon



	//used for shortestpath3()

	bool					dest_is_acceptor; //because any edge between acceptor and donor is both "exon" and "intron" edge,

											//use this to distinguish which edge is valid in path (final exon must be from acc to donor...)



	//for exons (set default to true, so don't have to change method 2.0/2.1 edge functions)

	bool					dest_matched; //records whether the last base pair is an exact match or not, if yes, will have to deduct from shared base pair when computing path length



	ACCP_DONR_Edge() {}



	//ACCP_DONR_Edge(int source, int dest, bool isC, float* p, int* left, int* remain)

	ACCP_DONR_Edge(int source, int dest, bool isC, AD_Edge_Penalty* p, char* left, char* remain, 

		const char* sbp, const char* ebp)

	{

		source_num = source;

		dest_num = dest;

		isCut = isC;

		int i;

		for (i=0; i<3; i++)

		{

			penalty[i] = p[i];

			leftover[i] = left[i];

			remainder[i] = remain[i];

		}



		for (i=0; i<2; i++)

		{

			start_bp[i] = sbp[i];

			end_bp[i] = ebp[i];

		}

		start_bp[2] = end_bp[2] = '\0';

		dest_matched = false;

	}



	//used only by CreateEdges3()

	ACCP_DONR_Edge(int source, int dest, bool isC, AD_Edge_Penalty* p, char* left, char* remain, 

		const char* sbp, const char* ebp, bool dest_is_acc, bool dest_is_match)

	{

		source_num = source;

		dest_num = dest;

		isCut = isC;

		int i;

		for (i=0; i<3; i++)

		{

			penalty[i] = p[i];

			leftover[i] = left[i];

			remainder[i] = remain[i];

		}



		for (i=0; i<2; i++)

		{

			start_bp[i] = sbp[i];

			end_bp[i] = ebp[i];

		}

		start_bp[2] = end_bp[2] = '\0';



		dest_is_acceptor = dest_is_acc;

		dest_matched = dest_is_match;

	}



	//this constructor is only used to create special starting edges from "start_site-1", only need one penalty[0],remainder[0]

	ACCP_DONR_Edge(int source, int dest, bool isC, AD_Edge_Penalty& p, char remain)

	{

		source_num = source;

		dest_num = dest;

		isCut = isC;

			penalty[0] = p;

			remainder[0] = remain;

	}



	ACCP_DONR_Edge(const ACCP_DONR_Edge& e)

	{

		source_num = e.source_num;

		dest_num = e.dest_num;

		isCut = e.isCut;

		int i;

		for (i=0; i<3; i++)

		{

			penalty[i] = e.penalty[i];

			leftover[i] = e.leftover[i];

			remainder[i] = e.remainder[i];

		}



		for (i=0; i<2; i++)

		{

			start_bp[i] = e.start_bp[i];

			end_bp[i] = e.end_bp[i];

		}

		start_bp[2] = end_bp[2] = '\0';

		dest_is_acceptor = e.dest_is_acceptor;

		dest_matched = e.dest_matched;

	}



	bool operator==(ACCP_DONR_Edge& e)

	{

		return (source_num == e.source_num && dest_num == e.dest_num && isCut == e.isCut 

			&& dest_is_acceptor == e.dest_is_acceptor && dest_matched == e.dest_matched 

			&& penalty[0].IsSame(e.penalty[0]) && penalty[1].IsSame(e.penalty[1]) && penalty[2].IsSame(e.penalty[2]) 

			&& leftover[0] == e.leftover[0] && leftover[1] == e.leftover[1] && leftover[2] == e.leftover[2]

			&& remainder[0] == e.remainder[0] && remainder[1] == e.remainder[1] && remainder[2] == e.remainder[2] );

	}



	friend ostream& operator<<(ostream& os, const ACCP_DONR_Edge& e)

	{

		os << "edge[" << e.source_num << " .. " << e.dest_num << "]:isCut?" << e.isCut << "; dest_is_acceptor?" << e.dest_is_acceptor 

			<< "start:" << e.start_bp << "; end:" << e.end_bp << "; dest_matched:" << e.dest_matched << "\n";

		for (int i=0; i<3; i++)

			os << "\tstart_frame(" << (int)(e.leftover[i]) << "):" << e.penalty[i] << " (remain:" << (int)(e.remainder[i]) << "\n";

		return os;

	}

};



//3 frames

class ACCP_DONR_Edge_Scores : public Base_Edge

{

public:

//	float					penalty[3]; //one penalty for each frame

	AD_Edge_Score			score[3];



	//used only for exons (record the number of nucleotides in a possibly 'broken' codon, must be [0,3)

	char						leftover[3]; //start positions corresponding to each penalty

	char						remainder[3]; //the rest of nucleotides remained from current exon

	

	//used only for introns

	char					start_bp[3]; //record the last two chars (base pair) in previous exon (last space for '\0')

	char					end_bp[3]; //record the first two chars in next exon



	//used for shortestpath3()

	bool					dest_is_acceptor; //because any edge between acceptor and donor is both "exon" and "intron" edge,

											//use this to distinguish which edge is valid in path (final exon must be from acc to donor...)



	//for exons (set default to true, so don't have to change method 2.0/2.1 edge functions)

	bool					dest_matched; //records whether the last base pair is an exact match or not, if yes, will have to deduct from shared base pair when computing path length



	ACCP_DONR_Edge_Scores() {}



	//ACCP_DONR_Edge(int source, int dest, bool isC, float* p, int* left, int* remain)

	ACCP_DONR_Edge_Scores(int source, int dest, bool isC, AD_Edge_Score p[3], char* left, char* remain, 

		const char* sbp, const char* ebp)

	{

		source_num = source;

		dest_num = dest;

		isCut = isC;

		int i;

		for (i=0; i<3; i++)

		{

			score[i] = p[i];

			leftover[i] = left[i];

			remainder[i] = remain[i];

		}



		for (i=0; i<2; i++)

		{

			start_bp[i] = sbp[i];

			end_bp[i] = ebp[i];

		}

		start_bp[2] = end_bp[2] = '\0';

		dest_matched = false;

	}



	//used only by CreateEdges3()

	ACCP_DONR_Edge_Scores(int source, int dest, bool isC, AD_Edge_Score p[3], char* left, char* remain, 

		const char* sbp, const char* ebp, bool dest_is_acc, bool dest_is_match)

	{

		source_num = source;

		dest_num = dest;

		isCut = isC;

		int i;

		for (i=0; i<3; i++)

		{

			score[i] = p[i];

			leftover[i] = left[i];

			remainder[i] = remain[i];

		}



		for (i=0; i<2; i++)

		{

			start_bp[i] = sbp[i];

			end_bp[i] = ebp[i];

		}

		start_bp[2] = end_bp[2] = '\0';



		dest_is_acceptor = dest_is_acc;

		dest_matched = dest_is_match;

	}



	//this constructor is only used to create special starting edges from "start_site-1", only need one penalty[0],remainder[0]

	ACCP_DONR_Edge_Scores(int source, int dest, bool isC, AD_Edge_Score& p, char remain)

	{

		source_num = source;

		dest_num = dest;

		isCut = isC;

			score[0] = p;

			remainder[0] = remain;

	}



	ACCP_DONR_Edge_Scores(const ACCP_DONR_Edge_Scores& e)

	{

		source_num = e.source_num;

		dest_num = e.dest_num;

		isCut = e.isCut;

		int i;

		for (i=0; i<3; i++)

		{

			score[i] = e.score[i];

			leftover[i] = e.leftover[i];

			remainder[i] = e.remainder[i];

		}



		for (i=0; i<2; i++)

		{

			start_bp[i] = e.start_bp[i];

			end_bp[i] = e.end_bp[i];

		}

		start_bp[2] = end_bp[2] = '\0';

		dest_is_acceptor = e.dest_is_acceptor;

		dest_matched = e.dest_matched;

	}



	bool operator==(ACCP_DONR_Edge_Scores& e)

	{

		return (source_num == e.source_num && dest_num == e.dest_num && isCut == e.isCut 

			&& dest_is_acceptor == e.dest_is_acceptor && dest_matched == e.dest_matched 

			&& score[0].pid_penalty.IsSame(e.score[0].pid_penalty) 

			&& score[1].pid_penalty.IsSame(e.score[1].pid_penalty) 

			&& score[2].pid_penalty.IsSame(e.score[2].pid_penalty) 

			&& leftover[0] == e.leftover[0] && leftover[1] == e.leftover[1] && leftover[2] == e.leftover[2]

			&& remainder[0] == e.remainder[0] && remainder[1] == e.remainder[1] && remainder[2] == e.remainder[2] );

	}



	friend ostream& operator<<(ostream& os, const ACCP_DONR_Edge_Scores& e)

	{

		os << "edge[" << e.source_num << " .. " << e.dest_num << "]:isCut?" << e.isCut << "; dest_is_acceptor?" << e.dest_is_acceptor 

			<< "start:" << e.start_bp << "; end:" << e.end_bp << "; dest_matched:" << e.dest_matched << "\n";

		for (int i=0; i<3; i++)

			os << "\tstart_frame(" << (int)(e.leftover[i]) << "):" << e.score[i] << " (remain:" << (int)(e.remainder[i]) << "\n";

		return os;

	}

};



class Edge : public Base_Edge

{

public:

	float					penalty; //this is the total edge penalty; for cut edge, it's two parts: (-x|y-), previous group and next group

	float					penalty_nextgroup; //this records the 'y' part of the penalty from (-x|y-)



	Edge(int source, int dest, float p, float p_ng, bool isC)

	{

		source_num = source;

		dest_num = dest;

		penalty = p; 

		penalty_nextgroup = p_ng;

		isCut = isC;

	}



	friend ostream& operator<<(ostream& os, const Edge& h)

	{

		os << "edge:("<< h.source_num << "-" << h.dest_num  << "); penalty: " << h.penalty 

			<< "; nextpart penalty: " << h.penalty_nextgroup 

			<< "; isCut? "<< h.isCut << "\n";

		return os;

	}



};



//******************************************************************//

//each node is a HSP

class HSP_node

{

public:

	int						HSP_num; //the internal serial index of the HSP (in HSP_gene)

	float					HSP_score; //precompute score for each node



	vector<int>			parents; //used in Ext/Cut edge construction



	map<int, float>			gene_pid; //used only for skip distance computation

	int						last_skip_end;



	HSP_node(int num): HSP_num(num) {

		HSP_score = 0; //score initially 0

		last_skip_end = -1;

	}



	~HSP_node() {

	}



	void UpdateMap_For_SkipDist(int dest_start, int dest_end, float dest_pid);

	float CompSkipDist(int prev_start, int cur_start);



};



//******************************************************************//

//the following stuff used only for constructing skip edges 

//******************************************************************//



//***record skipping distance of HSP pairs***

//HSP pair header

struct pair_dist

{

	int						source_num;

	int						dest_num;



	pair_dist(int s, int d): source_num(s), dest_num(d) {}



//	bool operator<(const pair_dist& e)

//	{ return source_num < e.source_num || (source_num == e.source_num && dest_num < e.dest_num) ; }

};



//method to order the "pair_dist"

struct pair_dist_compare

{

	bool operator()(pair_dist const& d1, pair_dist const& d2) const

	{ return d1.source_num < d2.source_num || (d1.source_num == d2.source_num && d1.dest_num < d2.dest_num) ; }

};



//record some info about each HSP pair

struct distance_info

{

	float					penalty; //skipping distance

	bool					isMono; //true if at least one path is mono

	bool					allMono; //true if all paths are mono



	distance_info(float p, bool isM, bool allM) : penalty(p), isMono(isM), allMono(allM) {}



};

//*******************************************



//used in TrySkipEdge(), order the potential HSP nodes (the other end of skip edge) by pid and position

struct PID_HSPstart_Index

{

	float				pid;

	int					HSP_start;

	int					HSP_end;

	int					HSP_num;



	PID_HSPstart_Index(float p, int start, int end, int i): pid(p), HSP_start(start), HSP_end(end), HSP_num(i) {}



	bool operator<(const PID_HSPstart_Index& p) const//sort first by pid (big first), then by HSP_start (small first)

	//{ return pid > p.pid || (pid == p.pid && HSP_start < p.HSP_start); }

	{ return pid > p.pid || (pid == p.pid && HSP_end < p.HSP_end) || (pid == p.pid && HSP_end == p.HSP_end && HSP_start < p.HSP_start); } //non-strict version



};



//******************************************************************//



//used in ShortestPath3()

struct EdgeInfo_QuadLet

{

	int source_node_num;

	bool dest_is_acceptor;

	bool is_cut; //true: intron; false: exon

	int prev_remainder;

	char prev_bp[3]; //last position is for '\0'

	char next_bp[3]; 

	bool exonedge_end_is_match;



	EdgeInfo_QuadLet(int source, bool dest_is_acc, bool is_c, int prev_rem, const char* p_bp, const char* n_bp, 

		bool end_is_match)

	{

		source_node_num = source;

		dest_is_acceptor = dest_is_acc;

		is_cut = is_c;

		prev_remainder = prev_rem;

		int i;

		for (i=0; i<3; i++)

		{

			prev_bp[i] = p_bp[i];

			next_bp[i] = n_bp[i];



		}

		exonedge_end_is_match = end_is_match;

	}



	bool operator==(const EdgeInfo_QuadLet& edgeinfo) const

	{

		return source_node_num == edgeinfo.source_node_num 

			&& dest_is_acceptor == edgeinfo.dest_is_acceptor 

			&& is_cut == edgeinfo.is_cut && exonedge_end_is_match == edgeinfo.exonedge_end_is_match;

	}



	bool operator<(const EdgeInfo_QuadLet& edgeinfo) const

	{

		return (source_node_num < edgeinfo.source_node_num) 

			|| (source_node_num == edgeinfo.source_node_num && is_cut > edgeinfo.is_cut) //if same distance, intron is preferred

			|| (source_node_num == edgeinfo.source_node_num && is_cut == edgeinfo.is_cut 

			    && dest_is_acceptor < edgeinfo.dest_is_acceptor); //is this possible? same distance, same source, both intron/exon, 

																	//1 source is acceptor, 1 source is donor?

	}



	friend ostream& operator<<(ostream& os, const EdgeInfo_QuadLet& edgeinfo)

	{

		os << "source: " << edgeinfo.source_node_num <<"; isCut? " << edgeinfo.is_cut 

			<< "; dest_is_acceptor? " << edgeinfo.dest_is_acceptor << "; prev_remainder: " << edgeinfo.prev_remainder

			<< "; next_bp:" << edgeinfo.next_bp << "; prev_bp:" << edgeinfo.prev_bp 

			<< "; exonedge_end_is_match?" << edgeinfo.exonedge_end_is_match;



		return os;

	}

};




struct IntronAndDonorTail

{

	int first; //donor_site;

	int second; //acceptor_site;



	int donor_tail; //acceptor_head must be its 3-complement

	int donor_hsp_index; //the index(position) of the HSP (in "HSPs" vector) that was used to calc donor splice alignment

	int acceptor_hsp_index; //for acceptor



	//use the following for GenePredict6(), to keep track of query gaps resulted from using current donor and acceptor pair

	int donor_query_end;

	int acceptor_query_start;



	//use the following for ComputeExon(), to step back when there's in-frame stop and we must step back

	int donor_region;



/*	IntronAndDonorTail(int d, int a, int t, int d_hsp_index, int a_hsp_index)

	{

		//donor_site = d;

		first = d;

		//acceptor_site = a;

		second = a;

		donor_tail = t;

		donor_hsp_index = d_hsp_index;

		acceptor_hsp_index = a_hsp_index;

	}

*/

	IntronAndDonorTail(int d, int a, int t, int d_hsp_index, int a_hsp_index, int d_query_e = 0, int a_query_s = 0, 

		int d_region = -1)

	{

		//donor_site = d;

		first = d;

		//acceptor_site = a;

		second = a;

		donor_tail = t;

		donor_hsp_index = d_hsp_index;

		acceptor_hsp_index = a_hsp_index;

		donor_query_end = d_query_e;

		acceptor_query_start = a_query_s;

		donor_region = d_region;

	}



	friend ostream& operator<<(ostream& os, const IntronAndDonorTail& idt)

	{

		os << "donor:" << idt.first << "(tail:" << idt.donor_tail << ")(hsp_index:" << idt.donor_hsp_index 

			<< ") to acceptor:" << idt.second << "(hsp_index:" << idt.acceptor_hsp_index << ")" 

			<< "(donor query end:" << idt.donor_query_end << ";acc query start:" << idt.acceptor_query_start << ")" 

			<< "; donor_region:" << idt.donor_region

			<< "\n";

		return os;

	}

};



struct ExonSiteInfo //used only by GenePredict6, for exon repairs

{

	int first; //target_site;

	int second; //query_site;



	int hsp_index;



	int frame; //tail for donor, head for acceptor



	ExonSiteInfo() { first = second = 0; }



	ExonSiteInfo(int t, int q, int h_index, int f)

	{

		first = t;

		second = q;

		hsp_index = h_index;

		frame = f;

	}



	bool operator==(const ExonSiteInfo& exon_site) const

	{

		return first == exon_site.first && second == exon_site.second && frame == exon_site.frame;

	}



	bool operator<(const ExonSiteInfo& exon_site) const

	{

		return first < exon_site.first || (first == exon_site.first && second < exon_site.second);

	}



	friend ostream& operator<<(ostream& os, const ExonSiteInfo& exon_site)

	{

		os << "target:" << exon_site.first << "; query:" << exon_site.second << "; hsp_index:" << exon_site.hsp_index

			<<  "; frame:" << exon_site.frame ;

		return os;

	}

};





#endif /* EDGE_H */

