#include <assert.h>

#include "c45extern.h"
#include "graph.h"

#include <float.h>
#include <time.h>

#include <queue>
#include <iomanip>
using namespace std;





//****** compute cut/ext edge penalty ******

//case 1: the HSP is the last in group, compute the missing gene penalty

float HSP_graph::GetCutPenalty_ToEnd(int pos)

{

	float result = 0;

	

	map<int, float>::iterator it = dataMgr.fragment_score_map.find(pos);

	while ( it != dataMgr.fragment_score_map.end() )

	{

		it++; //advance it first, so score does not include "pos"

		if (it != dataMgr.fragment_score_map.end())

			result += (*it).second;

		//it++;

	}



	return result;



}



//case 2: the HSP is the first in group, compute the missing gene penalty

float HSP_graph::GetCutPenalty_ToStart(int pos)

{



	float result = 0;



	map<int, float>::iterator it = dataMgr.fragment_score_map.find(pos);

	while ( it != dataMgr.fragment_score_map.begin() )

	{

		it--;

		result += (*it).second;

	}



	return result;

}



//case 3: two HSPs have gap between, compute missing gene penalty corresponding to that gap

float HSP_graph::GetPenalty_Between(int startpos, int endpos, bool end_included)

{

	if (!end_included)

		if (endpos <= startpos + 1)

			return 0;



	float result = 0;

	map<int,float>::iterator sIt, tmpIt;

	if (startpos == -1 || startpos == 0)

		sIt = dataMgr.fragment_score_map.begin();

	else

		sIt = dataMgr.fragment_score_map.find(startpos);

	//sIt++; //advance sIt, so that the result does not include score of [start,start] (sIt cannot be the last iterator, since startpos is not last)

	tmpIt = sIt;

	map<int, float>::iterator endIt = dataMgr.fragment_score_map.find(endpos);



/*	if (isPosStrand && (endpos == INT_MAX || endpos == INT_MAX-1))

		endIt = dataMgr.fragment_score_map.end();

	if (!isPosStrand && (startpos == -1 || startpos == 0))

		sIt = dataMgr.fragment_score_map.begin();

*/	

	while ( sIt != endIt)

	{

		result += (*sIt).second;

		sIt++;

	}

	if (!end_included)

		result -= (*tmpIt).second; //subtract [start,start]

	else

		result += (*endIt).second; //add [end, end]



	return result;



}

//********************************************

bool HSP_graph::TooFar(int source_hsp_end, int dest_hsp_start)

{



	if (dest_hsp_start - source_hsp_end > MAX_GAP_BTWN_HSP)

		return true;

	else

		return false;

}







//return true if two HSPs overalp too much in query segments (defined by user input "MAX_OVERLAP_PER")

//Update: MAX_OVERLAP_PER is no longer used. This function now only provides pos_cut/neg_cut results

bool HSP_graph::TooMuchOverlap(int source_start, int source_end, int dest_start, int dest_end, bool& pos_cut, bool& neg_cut)

{

	bool overlap_exceed_threshold = false;

/*

	int source_len = source_end - source_start;

	int dest_len = dest_end - dest_start;

	int base_len = (source_len < dest_len ? source_len : dest_len);



	if ( dest_start < source_start )

	{

		pos_cut = true;

		if (dest_end > source_end)

			overlap_exceed_threshold = true;

		else

		{

			if (dest_end > source_start)

				overlap_exceed_threshold = ((float)(dest_end-source_start)/(float)base_len > MAX_OVERLAP_PER);

		}

			

	}

	else

	{

		neg_cut = true;

		if (dest_end < source_end)

			overlap_exceed_threshold = true;

		else

		{

			if (dest_start < source_end)

				overlap_exceed_threshold = ((float)(source_end-dest_start)/(float)base_len > MAX_OVERLAP_PER);

			

		}

	}

*/

	//if ( dest_start < source_start )

	//	pos_cut = true;

	//else

	//	neg_cut = true;

	if ( dest_start <= source_start || source_end >= dest_end )

		pos_cut = true;

	if ( dest_start >= source_start || source_end <= dest_end )

		neg_cut = true;



	return overlap_exceed_threshold;



}





//Note: cut edge is (for positive strand): when source.gene_end > dest.gene_end; 

//      or source.gene_end<dest.gene_end && source.gene_start > dest.gene_start!!!

//also: when source.hsp_end and dest.hsp_start are too far apart (more than MAX_GAP_BTWN_HSP base pairs)

bool HSP_graph::ExtCutPenalty(HSP_node* source, HSP_node* dest, float& gap_penalty, float& gap_penalty_nextgroup, 

							  float& overlap_penalty, vector<HSP_Gene_Pair>& curHSPs, bool isSpecial, bool isSkipEdge)

{

	int source_end = curHSPs[source->HSP_num].gene_end;

	int dest_start = curHSPs[dest->HSP_num].gene_start;



	int source_start = curHSPs[source->HSP_num].gene_start;

	int dest_end = curHSPs[dest->HSP_num].gene_end;



//	if (dataMgr.query_gene.compare("C53B7.5") == 0 && dataMgr.HSP_chr[chr_index].compare("X") == 0 && !isPosStrand)

//		int stophere2 = 1;



	bool pos_cut = false;

	bool neg_cut = false;

	bool overlap_exceed_threshold = TooMuchOverlap(source_start, source_end, dest_start, dest_end, pos_cut, neg_cut);

	bool hsp_too_far_apart = TooFar(curHSPs[source->HSP_num].HSP_end, curHSPs[dest->HSP_num].HSP_start);



	//******positive strand******

	if (isPosStrand) 

	{

		if (isSkipEdge)

		{

			if (hsp_too_far_apart || pos_cut) //if cut, then no edge constructed since skip edge must be extension edge

				return false; //this return value is used differently than ext/cut edge

			else

			{

				if (source_end < dest_start) //penalty only if there is gap between source_end and dest_start

					gap_penalty += GetPenalty_Between(source_end, dest_start, false);

				else //special case

					overlap_penalty += GetPenalty_Between(dest_start, source_end, true);

				return true; //this return value is used

			}

		}

		else

		{

		//if (source_end > dest_end || source_start > dest_start)//if (source_end > dest_start) //non-increasing, cut edge

		//if (!isSpecial && source_end > dest_start) //"stricter" cut edge

		if (!isSpecial && (hsp_too_far_apart || pos_cut )) //|| overlap_exceed_threshold ))

		{

			//compute 1st group missing fragments

			gap_penalty += GetCutPenalty_ToEnd(source_end);



			//2nd group

			gap_penalty_nextgroup = GetCutPenalty_ToStart(dest_start);

			gap_penalty += gap_penalty_nextgroup;

			//gap_penalty += GetCutPenalty_ToStart(dest_start);



//			gap_penalty = 0; //changed here!!!!!!!

			return true; //is cut edge

		}



		else //extension edge

		{

			if (source_end < dest_start) //penalty only if there is gap between source_end and dest_start

				gap_penalty += GetPenalty_Between(source_end, dest_start, false);

			else //special case

				overlap_penalty += GetPenalty_Between(dest_start, source_end, true);

		

			//gap_penalty += GetPenalty_Between(source_end, dest_start); //"stricter"



			return false; //is ext edge



		}

		}

	}



	//******negative strand******

	else

	{

		if (isSkipEdge)

		{

			if (hsp_too_far_apart || neg_cut)

				return false;

			else

			{

				if (source_start > dest_end) //penalty only if there is gap between source_end and dest_start

					gap_penalty += GetPenalty_Between(dest_end, source_start, false);

				else //special case

					overlap_penalty += GetPenalty_Between(source_start, dest_end, true);

				return true; //this return value is used

			}

		}

		else

		{

		//if (source_end < dest_end || source_start < dest_start)

		//if (!isSpecial && source_start < dest_end) //"stricter"

		if (!isSpecial && (hsp_too_far_apart || neg_cut )) // || overlap_exceed_threshold ))

		{

			//compute 1st group missing fragments

			gap_penalty += GetCutPenalty_ToStart(source_start);

		

			//2nd group

			gap_penalty_nextgroup = GetCutPenalty_ToEnd(dest_end);

			gap_penalty += gap_penalty_nextgroup;

			//gap_penalty += GetCutPenalty_ToEnd(dest_end);



//			gap_penalty = 0; //changed here!!!!!!!

			

			return true; //is cut edge

		}



		else //extension edge

		{

			if (source_start > dest_end) //penalty only if there is gap between source_end and dest_start

				gap_penalty += GetPenalty_Between(dest_end, source_start, false);

			else //special case

				overlap_penalty += GetPenalty_Between(source_start, dest_end, true);



			//gap_penalty += GetPenalty_Between(dest_end, source_start); //"stricter"



			return false; //is ext edge



		}

		}



	}

	

}



//physically add ext/cut edge

void HSP_graph::AddExtCutEdge(HSP_node* source, HSP_node* dest, vector<HSP_Gene_Pair>& curHSPs, bool isSpecial)

{

	dest->parents.push_back(source->HSP_num);



	float p = 0;

	float gap_p_nextgroup = 0;

	float gap_p = 0;

	float over_p = 0;

	bool isCut = true;



	//special case for 1st and last node

	if (source->HSP_num == 0) //1st node

		//edges.push_back(Edge(source->HSP_num, dest->HSP_num, 0.0f, false));

	{

		//2nd group

		if (isPosStrand)

			gap_p_nextgroup = p = GetCutPenalty_ToStart(curHSPs[dest->HSP_num].gene_start);

		else

			gap_p_nextgroup = p = GetCutPenalty_ToEnd(curHSPs[dest->HSP_num].gene_end);

			

		//edges.push_back(Edge(source->HSP_num, dest->HSP_num, p, p, true);



		source_nodes.insert(map<int, float>::value_type(dest->HSP_num, GENE_MS_BETA * gap_p_nextgroup));

	}

	else if (dest->HSP_num == num_nodes-1) //last node

		{

			//compute 1st group missing fragments

			if (isPosStrand)

				p = GetCutPenalty_ToEnd(curHSPs[source->HSP_num].gene_end);

			else

				p = GetCutPenalty_ToStart(curHSPs[source->HSP_num].gene_start);

			//gap_p_nextgroup = 0;

			//edges.push_back(Edge(source->HSP_num, dest->HSP_num, p, 0, true));

		}

		else

		{



			isCut = ExtCutPenalty(source, dest, gap_p, gap_p_nextgroup, over_p, curHSPs, isSpecial, false);

		

			p = Penalty_new(0, gap_p, over_p);



			//p -= source->HSP_score; //subtract the node score



			//edges.push_back(Edge(source->HSP_num, dest->HSP_num, p, isCut));

			//edges.push_back(Edge(source->HSP_num, dest->HSP_num, p, gap_p_nextgroup, isCut));

			

			if (isCut)

				source_nodes.insert(map<int, float>::value_type(dest->HSP_num, GENE_MS_BETA * gap_p_nextgroup));

		}



#ifdef DEBUG_VERSION

		dataMgr.outFile << "ext/cut edge from ID " << curHSPs[source->HSP_num].ID << "(" << curHSPs[source->HSP_num].HSP_start << "-" 

			<< curHSPs[source->HSP_num].HSP_end << ") to ID " << curHSPs[dest->HSP_num].ID << "(" 

			<< curHSPs[dest->HSP_num].HSP_start << "-" << curHSPs[dest->HSP_num].HSP_end << ")" 

			<< source->HSP_num << "-" << dest->HSP_num 

			<< "; Penalty_miss(" << gap_p 

			<< "); Penalty_overlap(" << over_p << "); P_total_before_subtract_node(" 

			<< p << "); node(" << source->HSP_score << ") isCut?" << isCut << "\n";

#endif

	p -= source->HSP_score; //subtract the node score



	//edges.push_back(Edge(source->HSP_num, dest->HSP_num, p, isCut));

	edges.push_back(Edge(source->HSP_num, dest->HSP_num, p, gap_p_nextgroup, isCut));



}



// OLD: (not used)

// if (gap_s <= gap_e), gap exists, maybe next segment; otherwise no gap means not next segment 

bool HSP_graph::IsNextGeneSeg_old(HSP_node* source, HSP_node* dest)

{

	int gap_s, gap_e; 



	if (isPosStrand)

	{

		gap_s = dataMgr.HSP_gene[chr_index][source->HSP_num].gene_end;

		gap_e = dataMgr.HSP_gene[chr_index][dest->HSP_num].gene_start;

	}

	else

	{

		gap_e = dataMgr.HSP_neg_gene[chr_index][source->HSP_num].gene_start;

		gap_s = dataMgr.HSP_neg_gene[chr_index][dest->HSP_num].gene_end;

	}





	if (gap_e < gap_s)

		return false;

	else

		if (gap_e == gap_s)

			return true;

		else //now gap_s < gap_e

		{



				map<int, float>::iterator it= dataMgr.fragment_score_map.find(gap_s);

				map<int, float>::iterator endIt = dataMgr.fragment_score_map.find(gap_e);

				for (;	it != endIt; it++)

				{

					if ( (*it).second > 0)

						return false;

				}



				return true;

		}



}



// NEW:

// for pos: if (dest->gene_start >= source->gene_start) && ( overlap || gap==0 || gap.gene_segment==0 ), is next segment

bool HSP_graph::IsNextGeneSeg_new(HSP_node* source, HSP_node* dest)

{

	int gap_s, gap_e, source_s, source_e, dest_s, dest_e; 



	if (isPosStrand)

	{

		gap_s = dataMgr.HSP_gene[chr_index][source->HSP_num].gene_end;

		gap_e = dataMgr.HSP_gene[chr_index][dest->HSP_num].gene_start;



		source_s = dataMgr.HSP_gene[chr_index][source->HSP_num].gene_start;

		dest_e = dataMgr.HSP_gene[chr_index][dest->HSP_num].gene_end;

	}

	else

	{

		gap_e = dataMgr.HSP_neg_gene[chr_index][source->HSP_num].gene_start;

		gap_s = dataMgr.HSP_neg_gene[chr_index][dest->HSP_num].gene_end;



		source_e = dataMgr.HSP_neg_gene[chr_index][source->HSP_num].gene_end;

		dest_s = dataMgr.HSP_neg_gene[chr_index][dest->HSP_num].gene_start;

	}





	if (gap_e < gap_s) //overlap

	{

		if (isPosStrand)

			if (source_s > gap_e || dest_e < gap_s)

				return false;

			else

				return true;

		else

			if (gap_s > source_e || gap_e < dest_s)

				return false;

			else

				return true;

	}

	else

		if (gap_e == gap_s)

			return true;

		else //now gap_s < gap_e (has gap)

		{



				map<int, float>::iterator it= dataMgr.fragment_score_map.find(gap_s);

				map<int, float>::iterator endIt = dataMgr.fragment_score_map.find(gap_e);

				for (;	it != endIt; it++)

				{

					if ( (*it).second > 0)

						return false;

				}



				return true;

		}



}



//constructing all ext/cut edges

void HSP_graph::CreateAllExtCutEdges(vector<HSP_Gene_Pair>& curHSPs)

{
	clock_t cur_start_time = clock();

//	cout << "CreateAllExtCutEdges" << "\n";



	map< pair<int, int>, bool>  map_HasPath;



	int i;

	for ( i=0; i<num_nodes; i++)

	{

		//construct node

		HSP_node* curNode = new HSP_node(i);

		int h_s = curHSPs[i].HSP_start;

		int h_e = curHSPs[i].HSP_end;

		curNode->HSP_score = Score(curHSPs[i].gene_start, curHSPs[i].gene_end, curHSPs[i].pid)*HSP_SK_ALPHA; //changed to factor in "GENE_MS_BETA"!?



		curHSPs[i].node = curNode;



#ifdef DEBUG_VERSION

		dataMgr.outFile << "node_" << i << ": " << curHSPs[i] << "; score:" << curNode->HSP_score << "\n";

		//PrintEdges(debugFile);

#endif



		//check curNode against "last_nodes"

		set<int>	parent_candidates; //list of candidates for edges (candidate->curNode)

		set<int>  counter_candidates; //record HSPs where an edge is added from, used to remove parent_candidates



		list<int>::iterator it = last_nodes.begin();

		while (it!=last_nodes.end())

		{

			//cout << "HSP " << (*it) << ": " << curHSPs[*it] << "\n";



			list<int>::iterator tmpIt = it; //hold it

			tmpIt++;



			int prev_h_e = curHSPs[(*it)].HSP_end;

			if (h_s >= prev_h_e) //do not overlap

			{

				AddExtCutEdge(curHSPs[(*it)].node, curNode, curHSPs, false);



				counter_candidates.insert(*it);



				last_nodes.erase(it); //does it affect the iterator "it"?



			}

			else //overlaps

			{

				//isSpecial: special case if curNode matches the immediate next gene segment of previous node

				if (h_e >= prev_h_e && IsNextGeneSeg_new(curHSPs[*it].node, curNode))

				{

					AddExtCutEdge(curHSPs[*it].node, curNode, curHSPs, true);



					counter_candidates.insert(*it);



					last_nodes.erase(it);



				}

				else

				{

					vector<int>::iterator parentIt = curHSPs[(*it)].node->parents.begin();

					for (; parentIt != curHSPs[(*it)].node->parents.end(); parentIt++)

					{

						//cout << "here1" << "\n";



						//find appropriate candidates first (can have edge to curNode)

						set<int> can_set;

						if (FindParent(*parentIt, curNode->HSP_num, curHSPs, can_set))

							can_set.insert(*parentIt);



						for (set<int>::iterator canSetIt = can_set.begin(); canSetIt != can_set.end(); canSetIt++)

						{

							//cout << "here2" << "\n";

							//check whether it's valid candidate (no conflict with anything in parent_candidates)

							bool overlap = true;

							set<int>::iterator setIt = parent_candidates.begin(); 

							set<int>::iterator tmpsetIt;

							while (overlap && setIt != parent_candidates.end())

							{

								//cout << "here3: from " <<  *canSetIt << " to " << *setIt << "\n";

								//debugFile << "here3: from " <<  *canSetIt << " to " << *setIt << "\n";



								//if (*canSetIt == 378 && *setIt == 247)

								//	int stophere = 1;

								tmpsetIt = setIt; //hold setIt

								tmpsetIt++;



								if (HasPath(*canSetIt, *setIt, curHSPs, map_HasPath))

									overlap = false;

								else

									if (HasPath(*setIt, *canSetIt, curHSPs, map_HasPath))

										parent_candidates.erase(setIt);



								setIt = tmpsetIt;

							}



							//if no conflict, insert it as a parent_candidate

							if (overlap)

								parent_candidates.insert(*canSetIt);

						}





					}

		

					curHSPs[(*it)].hasOverlap = true;

					curHSPs[i].hasOverlap = true;

				}



			}



			it = tmpIt;



		}

		

		//remove redundant candidates

		if (!counter_candidates.empty())

		{

		set<int>::iterator pIt = parent_candidates.begin(); 

		set<int>::iterator tmp_pIt;

		while (pIt != parent_candidates.end())

		{

			//cout << "source " << *pIt << "\n";



			tmp_pIt = pIt;

			tmp_pIt++;

			for (set<int>::iterator cIt = counter_candidates.begin(); cIt != counter_candidates.end(); cIt++)

			{

				if (HasPath(*pIt, *cIt, curHSPs, map_HasPath))

				{

					parent_candidates.erase(pIt);

					break;

				}

			}

			pIt = tmp_pIt;

		}

		}



		//now add edges physically

		for (set<int>::iterator pIt = parent_candidates.begin(); pIt != parent_candidates.end(); pIt++)

		{

			if (curHSPs[*pIt].HSP_end > curHSPs[i].HSP_start) //overlap, must be special case

				AddExtCutEdge(curHSPs[*pIt].node, curNode, curHSPs, true);

			else

				AddExtCutEdge(curHSPs[*pIt].node, curNode, curHSPs, false);

		}



		last_nodes.push_back(i);//add curNode to the list



/*		//for debug only

#ifdef VERBOSE

		debugFile << "last_nodes: ";

		for (it=last_nodes.begin(); it!=last_nodes.end(); it++)

			debugFile << (*it) << ",";

		debugFile << "\n";

		PrintEdges(debugFile);



#endif

*/

	}



	//sort edges so that we can easily find all edges from a source node

	sort(edges.begin(), edges.end());





	//preprocess for skip edge construction

	for (i=0; i<num_nodes-2; i++)

		CalcSkipDist_SP(i, curHSPs);


//	cout << "CreateAllExtCutEdges done" << "\n";
#ifdef DEBUG
	dataMgr.outFile << "time for CreateAllExtCutEdges:" << (double)(clock()-cur_start_time)/CLOCKS_PER_SEC << ";num_nodes:" << num_nodes << ";num_edges:" << edges.size() << "\n";
#endif
}



//check whether 'candidate' is an ancestor of 'child'

bool HSP_graph::FindParent(int candidate, int child, vector<HSP_Gene_Pair>& curHSPs, set<int>& can_set)

{

	if (curHSPs[candidate].HSP_end <= curHSPs[child].HSP_start)

		return true;



	if (curHSPs[child].HSP_end >= curHSPs[candidate].HSP_end && IsNextGeneSeg_new(curHSPs[candidate].node, curHSPs[child].node))

		return true;



	for (vector<int>::iterator it = curHSPs[candidate].node->parents.begin(); it != curHSPs[candidate].node->parents.end(); it++)

	{

		if (FindParent(*it, child, curHSPs, can_set))

			can_set.insert(*it);

	}



	return false;



}



//check whether 'source' has path to 'dest'

bool HSP_graph::HasPath(int source, int dest, vector<HSP_Gene_Pair>& curHSPs, map<pair<int,int>, bool>& map_HasPath)

{

	if (source == dest)

		return true;



	if (source > dest)

		return false;



	if (curHSPs[source].HSP_end <= curHSPs[dest].HSP_start)

		return true;



	map<pair<int,int>, bool>::iterator mapIt = map_HasPath.find(pair<int,int>(source, dest));

	if (mapIt != map_HasPath.end())

		return (*mapIt).second;



	bool has_path = false;

//	cout << "HasPath " << curHSPs[dest] << "\n";

	for (vector<int>::iterator it = curHSPs[dest].node->parents.begin(); it != curHSPs[dest].node->parents.end(); it++)

	{

//		cout << "HasPath? from " << source << " to " << *it << "\n";

		if (HasPath(source, *it, curHSPs, map_HasPath))

		{

			has_path = true;

			break;

		}		

	}



	map_HasPath.insert(map< pair<int, int>, bool>::value_type(pair<int, int>(source, dest), has_path));

	return has_path;

}





//if source->dest path (any one path) is monotically increasing, no edge is added; otherwise add edge

//if all paths are mono && dest node does not overlap with any other node, return false (no need to go further)

//otherwise return true (still need to go further)

bool HSP_graph::AddSkipEdge(HSP_node* source, HSP_node* dest, vector<HSP_Gene_Pair>& curHSPs)

{



	pair_dist curpair(source->HSP_num, dest->HSP_num);

	map<pair_dist, distance_info, pair_dist_compare>::iterator it= all_pair_dist.find(curpair);



	if (!(*it).second.isMono)

	{

		//construct edge physically, with penalty (it's never a cut edge)

		float s_penalty = (*it).second.penalty; //HSP skip penalty



		//add missing_gene_segment penalty

		float g_penalty = 0;

		float o_penalty = 0;

		float g_penalty_nextgroup = 0; //skip edge must be extension edge, so this must be zero/not used

		bool need_edge = ExtCutPenalty(source, dest, g_penalty, g_penalty_nextgroup, o_penalty, curHSPs, false, true); 



		if (need_edge)

		{

		float penalty = Penalty_new(s_penalty, g_penalty, o_penalty);



#ifdef DEBUG_VERSION

		dataMgr.outFile << "skip edge from ID " << curHSPs[source->HSP_num].ID << "("

			<< curHSPs[source->HSP_num].HSP_start << "-" << curHSPs[source->HSP_num].HSP_end 

			<< ") to ID " << curHSPs[dest->HSP_num].ID 

			<< "(" << curHSPs[dest->HSP_num].HSP_start << "-" << curHSPs[dest->HSP_num].HSP_end 

			<< ")" << source->HSP_num << "-" << dest->HSP_num 

			<< "; Penalty_skip(" << s_penalty << "); Penalty_miss(" << g_penalty 

			<< "); Penalty_overlap(" << o_penalty << "); P_total_before_subtract_node(" 

			<< penalty << "); node(" << source->HSP_score << ")" << "\n";

#endif



		penalty -= source->HSP_score; //subtract node score



		edges.push_back(Edge(source->HSP_num, dest->HSP_num, penalty, g_penalty_nextgroup, false));

		}

		return true;

	}

	else

		if ((*it).second.allMono && !curHSPs[dest->HSP_num].hasOverlap)

			return false;

		else

			return true;



}



//for skip edge construction

bool HSP_graph::TrySkipEdge(vector<PID_HSPstart_Index>&	curPIDs, int i, vector<HSP_Gene_Pair>& curHSPs, 

							//int& lastHSPstart)

							int& lastHSPend) //non-strict version

{

//		cout << "TrySkipEdge" << "\n";



		//process HSPs that with same gene_start (all > A, all after A)

		sort(curPIDs.begin(), curPIDs.end());



		vector<PID_HSPstart_Index>::iterator pidIt = curPIDs.begin();

		if (!AddSkipEdge(curHSPs[i].node, curHSPs[(*pidIt).HSP_num].node, curHSPs))  //always add 1st edge

			return false;



		float lastPID = curHSPs[(*pidIt).HSP_num].pid;

		//lastHSPstart = curHSPs[(*pidIt).HSP_num].HSP_start; //update lastHSPstart

		//int lastHSPend = curHSPs[(*pidIt).HSP_num].HSP_end;

		int lastHSPstart = curHSPs[(*pidIt).HSP_num].HSP_start; 

		lastHSPend = curHSPs[(*pidIt).HSP_num].HSP_end;//update lastHSPend



		pidIt++;



		while (pidIt != curPIDs.end() )

		{

			float curPID = (*pidIt).pid;

			int curHSPstart = (*pidIt).HSP_start;

			int curIndex = (*pidIt).HSP_num;

			int curHSPend = curHSPs[curIndex].HSP_end;



			if (curPID == lastPID && lastHSPend > curHSPstart ) //same pid and overlap with last

			{

				if (!AddSkipEdge(curHSPs[i].node, curHSPs[curIndex].node, curHSPs))

					return false;

				lastHSPend = curHSPend; //non-strict version

			}

			else

			//if (curPID<lastPID && curHSPend <= lastHSPstart)

			if (curPID<lastPID && curHSPend <= lastHSPend)

			{

				if (!AddSkipEdge(curHSPs[i].node, curHSPs[curIndex].node, curHSPs))

					return false;



				lastPID = curPID;

				lastHSPstart = curHSPstart; //update lastHSPstart

				lastHSPend = curHSPend;

			}

			pidIt++;

		}



//		cout << "TrySkipEdge done" << "\n";



		return true;



}



//construct skip edges:

//allow skip edge from one HSP(GeneSeg1,GeneSeg2) to another HSP(GeneSeg2,GeneSeg3)

void HSP_graph::ProcHSPwithNextSeg(int i, vector<HSP_Gene_Pair>& curHSPs)

{

//	cout << "Creating skip edges from node_" << i << "\n";





	int fstgene_end = curHSPs[i].gene_end;

	int fstgene_start = curHSPs[i].gene_start;



	multimap<int, int>::iterator it;

	if (isPosStrand)

		//B: smallest gene_start that is bigger than fstgene_start

		//smallest gene_start that is >= fstgene_end //"stricter"

		//it = dataMgr.gene_start_HSP_num_map.lower_bound(fstgene_end); 

		it = dataMgr.gene_start_HSP_num_map.upper_bound(fstgene_start); //non-strict version

	else

		//from this, travel backward (biggest gene_start that is smaller than fstgene_start)

		//biggest gene_end that's <= fstgene_start!!!

		//it = dataMgr.gene_start_HSP_num_map.upper_bound(fstgene_start); //now it points to smallest gene_end that's > fstgene_start

		it = dataMgr.gene_start_HSP_num_map.lower_bound(fstgene_end); //non-strict version

	

	int lastHSPstart = INT_MAX;

	bool cont = true;



//	bool pos_cut;//not used, only for calling TooMuchOverlap(...)

//	bool neg_cut;//not used, only for calling TooMuchOverlap(...)



	//******positive strand******

	if (isPosStrand)

	while (cont && it != dataMgr.gene_start_HSP_num_map.end())

	{

		pair< multimap<int, int>::iterator, multimap<int, int>::iterator >  rangeIt = dataMgr.gene_start_HSP_num_map.equal_range((*it).first);



		vector<PID_HSPstart_Index>	curPIDs;

		for ( ; rangeIt.first != rangeIt.second; rangeIt.first++)

		{

			int index = (*(rangeIt.first)).second;



			if ( ( curHSPs[index].HSP_start > curHSPs[i].HSP_end) //B is at A's succeeding path

				//&& (curHSPs[index].HSP_end <= lastHSPstart) //)  //B is before 1st node in previous round; UPDATE! "lastHSPstart" now is HSP_end of last HSP

				&& (curHSPs[index].HSP_start <= lastHSPstart) //non-strict version

				&& (curHSPs[index].gene_end > fstgene_end) //) //B.gene_end >= A.gene_end 

				//&& !TooMuchOverlap(fstgene_start, fstgene_end, curHSPs[index].gene_start, curHSPs[index].gene_end, pos_cut, neg_cut)

				&& !TooFar(curHSPs[i].HSP_end, curHSPs[index].HSP_start) )

				curPIDs.push_back(PID_HSPstart_Index(curHSPs[index].pid, curHSPs[index].HSP_start, curHSPs[index].HSP_end, index));



		}



		it = (rangeIt.second); //update iterator for next loop



		if (curPIDs.empty()) //no B at A's succeeding path, check next range

			continue;



		cont = TrySkipEdge(curPIDs, i, curHSPs, lastHSPstart);

		

	}

	//******negative strand******

	else



	while (cont && it != dataMgr.gene_start_HSP_num_map.begin())

	{

		it--;

		pair< multimap<int, int>::iterator, multimap<int, int>::iterator >  rangeIt = dataMgr.gene_start_HSP_num_map.equal_range((*it).first);

		it = (rangeIt.first); //update iterator for next loop



		vector<PID_HSPstart_Index>	curPIDs;

		for ( ; rangeIt.first != rangeIt.second; rangeIt.first++)

		{

			int index = (*(rangeIt.first)).second;

			

			if ( ( curHSPs[index].HSP_start > curHSPs[i].HSP_end) //B is at A's succeeding path

				//&& (curHSPs[index].HSP_end <= lastHSPstart) //)  //B is before 1st node in previous round

				&& (curHSPs[index].HSP_start <= lastHSPstart) //non-strict version

				&& (curHSPs[index].gene_start < fstgene_start) //) //B.gene_start <= A.gene_start

				//&& !TooMuchOverlap(fstgene_start, fstgene_end, curHSPs[index].gene_start, curHSPs[index].gene_end, pos_cut, neg_cut) 

				&& !TooFar(curHSPs[i].HSP_end, curHSPs[index].HSP_start) )

				curPIDs.push_back(PID_HSPstart_Index(curHSPs[index].pid, curHSPs[index].HSP_start, curHSPs[index].HSP_end, index));



		}



		if (curPIDs.empty()) //no B at A's succeeding path, check next range

			continue;



		cont = TrySkipEdge(curPIDs, i, curHSPs, lastHSPstart);

		

	}



//	cout << "Creating skip edges from node_" << i << " done" << "\n";



}



//construct all skip edges

void HSP_graph::CreateAllSkipEdges(vector<HSP_Gene_Pair>& curHSPs)

{
	clock_t cur_start_time = clock();

	for (int i=0; i<num_nodes-2; i++) //do not need to include 2nd last and last node (special node)!!!

	{

		ProcHSPwithNextSeg(i, curHSPs);



	}
#ifdef DEBUG
	dataMgr.outFile << "time for CreateAllSkipEdges:" << (double)(clock()-cur_start_time)/CLOCKS_PER_SEC << ";num_nodes:" << num_nodes << ";num_edges:" << edges.size() << "\n";
#endif



}



//compute skipping distance of all potential HSP pairs

//(linear time algorithm for DAG shortest path (our nodes are already topologically sorted :-)!)

//(when skiping consecutive nodes with same (overlapping) gene segments, penalty should be for only one copy with max pid)

void HSP_graph::CalcSkipDist_SP(int source_num, vector<HSP_Gene_Pair>& curHSPs)

{

//	cout << "ShortestPath from " << source_num << "\n";



	if (edges.empty())

		return;



	//initialize distance

	float* dist=new float[num_nodes-source_num];

	bool* allMono = new bool[num_nodes-source_num];

	bool* isMono = new bool[num_nodes-source_num];

	bool* hasPath = new bool[num_nodes-source_num];

	dist[0]=0;

	allMono[0]=true;

	isMono[0]=true;

	hasPath[0]=true;

	int i;

	for ( i=1; i<num_nodes-source_num; i++)

	{

		dist[i]=FLT_MAX;

		allMono[i] = true;

		isMono[i] = true;

		hasPath[i] = false;

	}



	//dynamic programming part

	//for each vertex v in sorted order, for each outgoing edge e(v,u), 

	//if dist(v) + weight(e) < dist(u), set dist(u)=dist(v) + weight(e) and the predecessor of u to v

	curHSPs[source_num].node->gene_pid.clear();

	curHSPs[source_num].node->gene_pid.insert(map<int, float>::value_type(curHSPs[source_num].gene_start, curHSPs[source_num].pid));

	curHSPs[source_num].node->gene_pid.insert(map<int, float>::value_type(curHSPs[source_num].gene_end, 0));



	if (isPosStrand)

		curHSPs[source_num].node->last_skip_end = curHSPs[source_num].gene_end;

	else

		curHSPs[source_num].node->last_skip_end = curHSPs[source_num].gene_start;

	float cur_score=0;

	

	for (i=source_num; i<num_nodes; i++)

	{

		if (hasPath[i-source_num])

		{

		vector<Edge>::iterator it=lower_bound(edges.begin(), edges.end(), Edge(i, 0, 0, 0, true));



		while (it != edges.end() && (*it).source_num==i)

		{

			hasPath[(*it).dest_num-source_num] = true;



			int cur_dest_start = curHSPs[(*it).dest_num].gene_start;

			int cur_dest_end = curHSPs[(*it).dest_num].gene_end;

			float cur_dest_pid = curHSPs[(*it).dest_num].pid;



			int last_skip_end = curHSPs[i].node->last_skip_end;

			

			if ((*it).isCut)

			{

				allMono[(*it).dest_num-source_num] = false;

				isMono[(*it).dest_num-source_num] = false;



				curHSPs[(*it).dest_num].node->gene_pid.clear();

				curHSPs[(*it).dest_num].node->gene_pid.insert(map<int, float>::value_type(cur_dest_start, cur_dest_pid));

				curHSPs[(*it).dest_num].node->gene_pid.insert(map<int, float>::value_type(cur_dest_end, 0));



				if (isPosStrand)

				{

					cur_score = curHSPs[i].node->CompSkipDist(last_skip_end, curHSPs[i].gene_end);

					curHSPs[(*it).dest_num].node->last_skip_end = cur_dest_start;

				}

				else

				{

					cur_score = curHSPs[i].node->CompSkipDist(curHSPs[i].gene_start, last_skip_end);

					curHSPs[(*it).dest_num].node->last_skip_end = cur_dest_end;

				}

				



			}

			else

			{

				allMono[(*it).dest_num-source_num] = allMono[i-source_num];

				isMono[(*it).dest_num-source_num] = isMono[i-source_num];



				curHSPs[(*it).dest_num].node->gene_pid = curHSPs[i].node->gene_pid; //copy previous map first

				curHSPs[(*it).dest_num].node->UpdateMap_For_SkipDist(cur_dest_start, cur_dest_end, cur_dest_pid);



				if (isPosStrand)

				{

					cur_score = curHSPs[(*it).dest_num].node->CompSkipDist(last_skip_end, cur_dest_start);

					curHSPs[(*it).dest_num].node->last_skip_end = (last_skip_end < cur_dest_start ? cur_dest_start : last_skip_end);

				}

				else

				{

					cur_score = curHSPs[(*it).dest_num].node->CompSkipDist(cur_dest_end, last_skip_end);

					curHSPs[(*it).dest_num].node->last_skip_end = (last_skip_end > cur_dest_end ? cur_dest_end : last_skip_end);

				}



			}



			float d = dist[i-source_num] + cur_score;

			if (d < dist[(*it).dest_num-source_num])

				dist[(*it).dest_num-source_num] = d;





			it++;

		}



		}

	}



	//fill up all_pair_dist 

	for (i=1; i<num_nodes-source_num; i++)

		all_pair_dist.insert(map<pair_dist, distance_info>::value_type(pair_dist(source_num, source_num+i), 

			distance_info(dist[i], isMono[i], allMono[i])));



	delete [] dist;

	delete [] isMono;

	delete [] allMono;

	delete [] hasPath;



//	cout << "ShortestPath from " << source_num << " done" << "\n";



}





//OLD output, not used

bool HSP_graph::PrintPath_Backtrack(int* pred, bool* cut)

{



	dataMgr.outFile << "//**************START***************//" << "\n";



	if (isPosStrand)

		dataMgr.outFile << "//   positive strand: "; 

	else

		dataMgr.outFile << "//   negative strand: ";

	dataMgr.outFile << "(backtrack)   //" << "\n";

	dataMgr.outFile << "//**********************************//" << "\n" ;

	int i=num_nodes-1; //last node, back track



	int count = 1;

	dataMgr.outFile << "\n" << "=== GROUP " << count << ":" << "\n";

	count++;

	while (pred[i]>0)

	{

		if (cut[i])

		{

			dataMgr.outFile << "\n" << "=== GROUP " << count << ":" << "\n";

			count++;

		}

		if (isPosStrand)

			dataMgr.outFile << "node_" << pred[i] << ": " << dataMgr.HSP_gene[chr_index][pred[i]];

		else

			dataMgr.outFile << "node_" << pred[i] << ": " << dataMgr.HSP_neg_gene[chr_index][pred[i]];

		i = pred[i];

	}

	dataMgr.outFile << "\n" << "//***************END****************//" << "\n" << "\n";



//	outFile.close();



	return true;



}



//THE OUTPUT FUNCTION: (store all groups)

//Groups should be ranked and output in ranked order;

//Each group should be summarized like this:

//Query_gene_id|Target_chr_id:Start..End|Strand|%Coverage_of_query_gene|Score|Ranking

bool HSP_graph::PrintPath(int* pred, bool* cut, float* penalty, float* penalty_nextgroup)

{

//	cout << "getting path info..." << "\n";



//	if (dataMgr.query_gene.compare("T26C11.2")==0)

//		int stop = 1;



	int i=num_nodes-1; //last node, back track

	float score = 0;

	vector<HSP_Gene_Pair*> indexes;

	//set<int>	HSP_IDs; //used to collect all HSP_IDs (the original ID in input file) that are included in the final output, sorted

	int HSP_start = 0;

	int HSP_end = 0;

	int gene_cover = 0;

	bool isFst = true;



	int gene_start = 0;

	int gene_end = 0;

	while (pred[i]>0)

	{

		if (cut[i] && (i < num_nodes-1) ) //skip the first round

		{

			score -= GENE_MS_BETA * penalty_nextgroup[i]; //last edge section of current group



			if (isPosStrand)

				HSP_start = dataMgr.HSP_gene[chr_index][i].HSP_start;

			else

				HSP_start = dataMgr.HSP_neg_gene[chr_index][i].HSP_start;



			//if (gene_cover >= MIN_GENE_COVER_PER * dataMgr.query_len) //only store the group with gene_cover exceed MIN_GENE_COVER_PER

			{

				Group_Info info(score, HSP_start, HSP_end, gene_cover, chr_index, isPosStrand);

				dataMgr.groups.insert(multimap<Group_Info, vector<HSP_Gene_Pair*> >::value_type(info, indexes));

			}



			indexes.clear();

			score = 0;

			isFst = true;

			gene_cover = 0;

		}



		//indexes.push_back(pred[i]);

		if (isPosStrand)

			indexes.push_back(&(dataMgr.HSP_gene[chr_index][pred[i]]));

		else

			indexes.push_back(&(dataMgr.HSP_neg_gene[chr_index][pred[i]]));



		//HSP_IDs.insert(dataMgr.GetHSPID(isPosStrand, chr_index, pred[i]));

		if (isPosStrand)

		{

			//score += dataMgr.HSP_gene[chr_index][pred[i]].node->HSP_score;

						

			if (isFst)

			{

				//score += dataMgr.HSP_gene[chr_index][pred[i]].node->HSP_score;

				score -= penalty[i] - GENE_MS_BETA * penalty_nextgroup[i]; //first edge section of current group

				gene_cover += dataMgr.HSP_gene[chr_index][pred[i]].gene_end - dataMgr.HSP_gene[chr_index][pred[i]].gene_start + 1;

				HSP_end = dataMgr.HSP_gene[chr_index][pred[i]].HSP_end;

				isFst = false;

			}

			else

			{

				score -= penalty[i]; //mid edges of current group

				if (dataMgr.HSP_gene[chr_index][pred[i]].gene_end >= gene_start)

					gene_cover += gene_start - dataMgr.HSP_gene[chr_index][pred[i]].gene_start;

				else

					gene_cover += dataMgr.HSP_gene[chr_index][pred[i]].gene_end - dataMgr.HSP_gene[chr_index][pred[i]].gene_start + 1;

			}

			gene_start = dataMgr.HSP_gene[chr_index][pred[i]].gene_start;

			

		}

		else

		{

			//score += dataMgr.HSP_neg_gene[pred[i]].node->HSP_score;

			

			if (isFst)

			{

				//score += dataMgr.HSP_neg_gene[chr_index][pred[i]].node->HSP_score;

				score -= penalty[i] - GENE_MS_BETA * penalty_nextgroup[i];

				gene_cover += dataMgr.HSP_neg_gene[chr_index][pred[i]].gene_end - dataMgr.HSP_neg_gene[chr_index][pred[i]].gene_start + 1;

				HSP_end = dataMgr.HSP_neg_gene[chr_index][pred[i]].HSP_end;

				isFst = false;

			}

			else

			{

				score -= penalty[i];

				if (dataMgr.HSP_neg_gene[chr_index][pred[i]].gene_start <= gene_end)

					gene_cover += dataMgr.HSP_neg_gene[chr_index][pred[i]].gene_end - gene_end;

				else

					gene_cover += dataMgr.HSP_neg_gene[chr_index][pred[i]].gene_end - dataMgr.HSP_neg_gene[chr_index][pred[i]].gene_start + 1;

			}

			gene_end = dataMgr.HSP_neg_gene[chr_index][pred[i]].gene_end;

		}



		



		i = pred[i];

	}

	//last group

	score -= GENE_MS_BETA * penalty_nextgroup[i];

	if (isPosStrand)

		HSP_start = dataMgr.HSP_gene[chr_index][i].HSP_start;

	else

		HSP_start = dataMgr.HSP_neg_gene[chr_index][i].HSP_start;

	//if (gene_cover >= MIN_GENE_COVER_PER * dataMgr.query_len) //only store the group with gene_cover exceed MIN_GENE_COVER_PER

	{

		Group_Info info(score, HSP_start, HSP_end, gene_cover, chr_index, isPosStrand);

		dataMgr.groups.insert(multimap<Group_Info, vector<HSP_Gene_Pair*> >::value_type(info, indexes));

	}



	//cut off everything below the TOP_RANK_NUM

	multimap<Group_Info, vector<HSP_Gene_Pair*> >::iterator groupMapIt, tmpGroupMapIt;

	int num_of_groups = dataMgr.groups.size();

	if (num_of_groups>TOP_RANK_NUM)

	{

		groupMapIt = dataMgr.groups.begin();

		for (i=0; i<TOP_RANK_NUM; i++)

			groupMapIt++;

		dataMgr.groups.erase(groupMapIt, dataMgr.groups.end());

	}



	//test for minimum gene coverage

	if (dataMgr.groups.size() > 1) //if there's only 1 group/gene, skip

	{

		groupMapIt=dataMgr.groups.begin();



		//keep the top-ranked one

		multimap<Group_Info, vector<HSP_Gene_Pair*> > tmpGroups;

		tmpGroups.insert(*groupMapIt);

		

		//groupMapIt++;

		//wrong: start from second group, so there's always at least 1 group in output

		while (groupMapIt != dataMgr.groups.end())

		{

			if ((*groupMapIt).first.gene_cover_len < MIN_GENE_COVER_PER * dataMgr.query_len)

			{

				//groupMapIt = dataMgr.groups.erase(groupMapIt); //groupMapIt points to the 1st element beyond the erased element (DARN: NON STANDARD!)

				tmpGroupMapIt = groupMapIt;

				tmpGroupMapIt++;

				dataMgr.groups.erase(groupMapIt);

				groupMapIt = tmpGroupMapIt;

			}

			else

				groupMapIt++;

		}



		//make sure we have at least 1 group left in output

		if (dataMgr.groups.empty())

		{

			groupMapIt = tmpGroups.begin();

			dataMgr.groups.insert(*groupMapIt);

		}

	}





//	cout << "stage 1: path stored" << "\n";





	return true;



}



void HSP_graph::ShortestPath()

{

//	cout << "num_nodes:" << num_nodes << "\n";



	if (edges.empty())

		return;



	//sort edges so that we can easily find all edges from a source node

	sort(edges.begin(), edges.end());



#ifdef DEBUG_VERSION

	PrintEdges(dataMgr.outFile);

#endif



	//initialize distance

	float* dist=new float[num_nodes];

	int* pred=new int[num_nodes];

	bool* cut = new bool[num_nodes];

	float* penalty = new float[num_nodes];

	float* penalty_nextgroup = new float[num_nodes];

	dist[0]=0;

	//pred[0]=-1;

	int i;

	for ( i=1; i<num_nodes; i++)

	{

		dist[i]=FLT_MAX;

		//pred[i]=-1;

	}



	//dynamic programming part

	//for each vertex v in sorted order, for each outgoing edge e(v,u), 

	//if dist(v) + weight(e) < dist(u), set dist(u)=dist(v) + weight(e) and the predecessor of u to v

	vector<Edge>::iterator it=edges.begin();

	for (i=0; i<num_nodes; i++)

	{

		//cout << "proc node " << i << "\n";

		while (it != edges.end() && (*it).source_num==i)

		{

			float d = dist[i]+(*it).penalty;

			if (d < dist[(*it).dest_num])

			{

				dist[(*it).dest_num] = d;

				pred[(*it).dest_num] = i;



				cut[(*it).dest_num] = (*it).isCut;

				penalty[(*it).dest_num] = (*it).penalty; //edge penalty = edge.cost - source_node.score

				penalty_nextgroup[(*it).dest_num] = (*it).penalty_nextgroup;

			}

			it++;

		}

	}



	//output path

	//PrintPath_Backtrack(pred, cut);

	PrintPath(pred, cut, penalty, penalty_nextgroup);



	delete [] dist;

	delete [] pred;

	delete [] cut;

	delete [] penalty;

	delete [] penalty_nextgroup;



//	cout << "result done" << "\n";



}



void HSP_graph::FindAllLocalShortestPaths()

{
	clock_t cur_start_time = clock();


	if (edges.empty())

		return;



	//sort edges so that we can easily find all edges from a source node

	sort(edges.begin(), edges.end());



#ifdef DEBUG_VERSION

	PrintEdges(dataMgr.outFile);

#endif



#ifdef DEBUG_VERSION

	dataMgr.outFile << "source nodes: " << "\n";

	map<int, float>::iterator setIt;

	for (setIt = source_nodes.begin(); setIt != source_nodes.end(); setIt++)

		dataMgr.outFile << (*setIt).first << "," << (*setIt).second << "; " ;

	dataMgr.outFile << "\n";

#endif



	map<int, float>::iterator it;

	float dist;

	//for (it = source_nodes.begin(); it!= source_nodes.end(); it++)

	while (!source_nodes.empty())

	{

		it = source_nodes.begin();



		//if ((*it).first == 90)

		//	dataMgr.outFile << "node 90 work start here" << "\n";



		vector<int> localpath;

		dist = LocalShortestPath((*it).first, (*it).second, localpath);

		//dist += (*it).second; //now add the weight before the source node (so path_dist is the total weight of the entire group)



		local_short_paths.insert(multimap<float, vector<int> >::value_type(dist, localpath));



#ifdef DEBUG_VERSION

		dataMgr.outFile << "source " << (*it).first << ": (path weight: " << dist << ")" << "\n";

		vector<int>::iterator pathIt;

		for (pathIt = localpath.begin(); pathIt != localpath.end(); pathIt++)

			dataMgr.outFile << "Node " << *pathIt << " - ";

		dataMgr.outFile << "\n";

#endif



		source_nodes.erase(it);

	}

//	dataMgr.outFile << "got all local paths" << "\n";



	//output path in weighted order

	int i, start_node, end_node, region_start, region_end, num_of_groups=0;

	set< pair<int, int> >  hsp_group_regions;

	vector<int>::iterator vecIt;

	//vector<int>::reverse_iterator vecRevIt, vecRevIt_tmp;

	multimap<float, vector<int> >::iterator mapIt = local_short_paths.begin();

	set<int> hsp_ids; //store the list of hsp_id that has already been used in the output

	float cur_score;
	set<float> rank_scores; //keep all encountered scores, for cut off by TOP_RANK_NUM
	while (mapIt != local_short_paths.end() )
	{
		cur_score = (*mapIt).first;
		if (rank_scores.size() >= TOP_RANK_NUM)
		{
			if (-cur_score < *(rank_scores.begin()) )
				break;
		}

		vector<int> path_tmp; //use this to collect real path (for option c)
		if (PHASE1_OUTPUT.compare("c") == 0)
		{

			//for (vecRevIt = (*mapIt).second.rbegin(); vecRevIt != (*mapIt).second.rend(); )

			for (vecIt = (*mapIt).second.begin(); vecIt != (*mapIt).second.end(); vecIt++)

			{

				if (isPosStrand)

					i = dataMgr.HSP_gene[chr_index][*vecIt].ID;

				else

					i = dataMgr.HSP_neg_gene[chr_index][*vecIt].ID;

				if (hsp_ids.find(i) != hsp_ids.end())

				{

					if (isPosStrand)

						cur_score +=(dataMgr.HSP_gene[chr_index][*vecIt].node)->HSP_score; //remove the HPS node and its score from the total score

					else

						cur_score +=(dataMgr.HSP_neg_gene[chr_index][*vecIt].node)->HSP_score;



					//(*mapIt).second.erase((++vecRevIt).base());



					//change to: erase all until the path end (they must all be overlapping)

/*					if (vecRevIt != (*mapIt).second.rbegin())

					{

						vecRevIt--;

						(*mapIt).second.erase((*mapIt).second.begin(), vecRevIt.base());

					}

					else

						(*mapIt).second.clear();

					break;

*/

				}

				else

				{

					hsp_ids.insert(i);

					path_tmp.push_back(*vecIt);

					//++vecIt;

				}

			}



			//if ((*mapIt).second.empty())

			if (path_tmp.empty())

			{

				mapIt++;

				continue;

			}

		}

		else //other options ('a' and 'b')

			 path_tmp = (*mapIt).second;

		if (rank_scores.size() >= TOP_RANK_NUM && -cur_score < *(rank_scores.begin()) )
			continue; //current group falls out of TOP_RANK_NUM, then skip it and keep checking next
		rank_scores.insert(-cur_score);


		//vecIt = (*mapIt).second.begin();

		vecIt = path_tmp.begin();

		end_node = *vecIt;

		//start_node = *((*mapIt).second.rbegin());

		start_node = *(path_tmp.rbegin());

		if (isPosStrand)

		{

			region_start = dataMgr.HSP_gene[chr_index][start_node].HSP_start;

			region_end = dataMgr.HSP_gene[chr_index][end_node].HSP_end;

		}

		else

		{

			region_start = dataMgr.HSP_neg_gene[chr_index][start_node].HSP_start;

			region_end = dataMgr.HSP_neg_gene[chr_index][end_node].HSP_end;

		}

#ifdef DEBUG_VERSION

		dataMgr.outFile << "current path weight: " << (*mapIt).first << " (source node " << start_node << ")" << "\n";

		dataMgr.outFile << "region start: " << region_start << "; region end: " << region_end << "\n";

#endif



		if (PHASE1_OUTPUT.compare("a")==0) 

		{

			if (RegionOverlap(region_start, region_end, hsp_group_regions))

			{

				num_of_groups_pruned_due_to_overlap++;

				mapIt++;

#ifdef DEBUG_VERSION

				dataMgr.outFile << "overlap" << "\n";

#endif

				continue;

			}



			hsp_group_regions.insert( pair<int,int>(region_start, region_end) );

		}



		int gene_start=INT_MAX, gene_end=0, gene_cover=0;

		vector<HSP_Gene_Pair*> indexes;

		//for (; vecIt != (*mapIt).second.end(); vecIt++)

		for (; vecIt != path_tmp.end(); vecIt++)

		{

			i = *vecIt;

#ifdef DEBUG_VERSION

			dataMgr.outFile << "node " << i << "-";

#endif



			if (isPosStrand)

			{

				indexes.push_back(&(dataMgr.HSP_gene[chr_index][i]));

				//if (dataMgr.HSP_gene[chr_index][i].gene_end >= gene_start)

					//gene_cover += gene_start - dataMgr.HSP_gene[chr_index][i].gene_start;

				//else

				//	gene_cover += dataMgr.HSP_gene[chr_index][i].gene_end - dataMgr.HSP_gene[chr_index][i].gene_start + 1;

				if (dataMgr.HSP_gene[chr_index][i].gene_end < gene_start)

					gene_cover += dataMgr.HSP_gene[chr_index][i].gene_end - dataMgr.HSP_gene[chr_index][i].gene_start + 1;

				else

					gene_cover += gene_start - dataMgr.HSP_gene[chr_index][i].gene_start;

				gene_start = dataMgr.HSP_gene[chr_index][i].gene_start;

			}

			else

			{

				indexes.push_back(&(dataMgr.HSP_neg_gene[chr_index][i]));

				if (dataMgr.HSP_neg_gene[chr_index][i].gene_start <= gene_end)

					gene_cover += dataMgr.HSP_neg_gene[chr_index][i].gene_end - gene_end;

				else

					gene_cover += dataMgr.HSP_neg_gene[chr_index][i].gene_end - dataMgr.HSP_neg_gene[chr_index][i].gene_start + 1;

				gene_end = dataMgr.HSP_neg_gene[chr_index][i].gene_end;

			}

		}

		//Group_Info info(-(*mapIt).first, region_start, region_end, gene_cover, chr_index, isPosStrand); //score is negative of path weight

		Group_Info info(-cur_score, region_start, region_end, gene_cover, chr_index, isPosStrand); //score is negative of path weight

		dataMgr.groups.insert(multimap<Group_Info, vector<HSP_Gene_Pair*> >::value_type(info, indexes));

#ifdef DEBUG_VERSION

		dataMgr.outFile << "\n";

		dataMgr.outFile << info;

#endif

		//num_of_groups++;



		mapIt++;

	}

//	dataMgr.outFile << "all paths stored" << "\n";



	//cut off everything below the TOP_RANK_NUM

	multimap<Group_Info, vector<HSP_Gene_Pair*> >::iterator groupMapIt; //, tmpGroupMapIt;

//	num_of_groups = dataMgr.groups.size();

//	if (num_of_groups>TOP_RANK_NUM)

//	{

#ifdef DEBUG

		groupMapIt = dataMgr.groups.begin();

		dataMgr.outFile << "all groups now: " << "\n";

		for (; groupMapIt != dataMgr.groups.end(); groupMapIt++)

			dataMgr.outFile << "cur score: " << (*groupMapIt).first.score << "\n";

#endif


		groupMapIt = dataMgr.groups.begin();

		dist = (*groupMapIt).first.score;



#ifdef DEBUG

		dataMgr.outFile << "top score is " << dist << "\n";

#endif

		num_of_groups = 1;

		groupMapIt++;

		while (groupMapIt != dataMgr.groups.end() && num_of_groups <= TOP_RANK_NUM)

		{

			//if ((*groupMapIt).first.score + Epsilon < dist)
			if ((*groupMapIt).first.score < dist)
			{

#ifdef DEBUG

				dataMgr.outFile << "lower score: " << (*groupMapIt).first.score << "\n";
#endif
				num_of_groups++;

				dist = (*groupMapIt).first.score;

				if (num_of_groups > TOP_RANK_NUM)

					break;

			}

#ifdef DEBUG

			else

				dataMgr.outFile << "same score: " << (*groupMapIt).first.score << " as " << dist << "\n";
#endif
			groupMapIt++;

		}

		if (groupMapIt != dataMgr.groups.end())

			dataMgr.groups.erase(groupMapIt, dataMgr.groups.end());

//	}



	//cut off groups that has low GENE_COVER than MIN_GENE_COVER_PER, and lower score than MIN_GENE_SCORE

	if (dataMgr.groups.size() > 1) //if there's only 1 group/gene, skip

	{

		groupMapIt=dataMgr.groups.begin();



		//keep the top-ranked one

		multimap<Group_Info, vector<HSP_Gene_Pair*> > tmpGroups;

		tmpGroups.insert(*groupMapIt);

		

		//groupMapIt++;

		//wrong: start from second group, so there's always at least 1 group in output

		multimap<Group_Info, vector<HSP_Gene_Pair*> >::iterator tmpGroupMapIt;

		while (groupMapIt != dataMgr.groups.end())

		{

			if ( (*groupMapIt).first.gene_cover_len < MIN_GENE_COVER_PER * dataMgr.query_len ||

				 (USE_MIN_GENE_SCORE && (*groupMapIt).first.score < MIN_GENE_SCORE) )

			{

				//groupMapIt = dataMgr.groups.erase(groupMapIt); //groupMapIt points to the 1st element beyond the erased element (DARN: NON STANDARD!)

				tmpGroupMapIt = groupMapIt;

				tmpGroupMapIt++;

				dataMgr.groups.erase(groupMapIt);

				groupMapIt = tmpGroupMapIt;

			}

			else

				groupMapIt++;

		}



		//make sure we have at least 1 group left in output

		if (dataMgr.groups.empty())

		{

			groupMapIt = tmpGroups.begin();

			dataMgr.groups.insert(*groupMapIt);

		}

	}
#ifdef DEBUG
	dataMgr.outFile << "time for FindAllLocalShortestPaths:" << (double)(clock()-cur_start_time)/CLOCKS_PER_SEC << ";num_nodes:" << num_nodes << ";num_edges:" << edges.size() << "\n";
#endif
}



float HSP_graph::LocalShortestPath(int node_id, float dist_before_fst_node, vector<int>& path)

{

	float path_dist = FLT_MAX;

	int path_end = INT_MAX; //if this doesn't change, then something's wrong?



	//initialize distance

	float* dist=new float[num_nodes - node_id];

	float* dist_ext = new float[num_nodes - node_id]; //keep track of node's distance when the node is used as an extension

	int* pred=new int[num_nodes - node_id];

	int* pred_ext = new int[num_nodes - node_id]; //used with dist_ext

//	bool* cut = new bool[num_nodes - node_id];

	dist[0]=0;

	dist_ext[0]=0;

	pred[0]=-1;

	pred_ext[0]=-1;

	int i;

	for ( i=node_id+1; i<num_nodes; i++)

	{

		dist[i-node_id]=FLT_MAX;

		dist_ext[i-node_id]=FLT_MAX;

	}



	set<int> nodes; //queue<int> nodes;

	nodes.insert(node_id); //nodes.push(node_id);

	vector<Edge>::iterator it;

	int cur_node;

	float d;

	while (!nodes.empty())

	{

#ifdef DEBUG_VERSION

		set<int>::iterator nodeIt = nodes.begin();

		dataMgr.outFile << "nodes in queue: ";

		for (; nodeIt != nodes.end(); nodeIt++)

			dataMgr.outFile << *nodeIt << "; ";

		dataMgr.outFile << "\n";

		for ( i=node_id; i<num_nodes; i++)

			if (dist[i-node_id]  < FLT_MAX)

				dataMgr.outFile << "node " << i << "dist: " << dist[i-node_id] << "; dist_ext: " << dist_ext[i-node_id] << "\n";

#endif



		//cur_node = nodes.front(); 

		cur_node = *(nodes.begin());

		//nodes.pop();

		nodes.erase(nodes.begin());

		

#ifdef DEBUG_VERSION

		dataMgr.outFile << "cur_node: " << cur_node << "\n";

#endif



		it = lower_bound(edges.begin(), edges.end(), Edge(cur_node, -1, 0, 0, true));



		while (it != edges.end() && (*it).source_num == cur_node)

		{

#ifdef DEBUG_VERSION

			dataMgr.outFile << "edge: " << *it << "\n";

#endif



			//float d = dist[cur_node-node_id]+(*it).penalty;

			//float d = dist[cur_node - node_id] + (*it).penalty - GENE_MS_BETA * (*it).penalty_nextgroup;

			d = dist_ext[cur_node - node_id] + (*it).penalty - GENE_MS_BETA * (*it).penalty_nextgroup;

			if (d < dist[(*it).dest_num-node_id])

			{

				dist[(*it).dest_num - node_id] = d;

				pred[(*it).dest_num - node_id] = cur_node;



//				cut[(*it).dest_num - node_id] = (*it).isCut;



			if (!(*it).isCut)

			{

				if (d < dist_ext[(*it).dest_num - node_id])

				{

					dist_ext[(*it).dest_num - node_id] = d;

					pred_ext[(*it).dest_num - node_id] = cur_node;

				}

			}

			else //the current node is a possible end of an extension path

			{

				//check if the current local path is a valid one (with only 1 cut between the path end and the source node)

				//if (dist[(*it).dest_num - node_id] < path_dist)

				float group_penalty = dist[(*it).dest_num - node_id]; //- GENE_MS_BETA * (*it).penalty_nextgroup;

				if ( group_penalty < path_dist)

				{

/*					int cut_count = 0;

					bool valid_path = true;

					i = (*it).dest_num;

					while (pred[i-node_id] >= node_id)

					{

						if (cut[i-node_id])

						{

							cut_count++;

							if (cut_count > 1)

							{

								valid_path = false;

								break;

							}

						}

						i = pred[i - node_id];

					}



					if (valid_path)

					{

*/						//path_dist = dist[(*it).dest_num - node_id];

						path_dist = group_penalty;

						path_end = (*it).dest_num;

#ifdef DEBUG_VERSION

						dataMgr.outFile << "current path end: " << path_end << "\n";

#endif

/*					}

					else

					{

#ifdef DEBUG_VERSION

						dataMgr.outFile << "path not valid (more than 1 cut)" << "\n";

#endif

					}

*/

				}

				else

				{

#ifdef DEBUG_VERSION

					dataMgr.outFile << "current path worse than previous" << "\n";

#endif

				}

			}

			}

			

			if (!(*it).isCut)

				nodes.insert((*it).dest_num);

				//nodes.push((*it).dest_num);



			it++;

		}

	}



	i = path_end;

	path.push_back(pred[i-node_id]); //last node (the end of extension path) must have cut edge before it, so should use pred[], not pred_ext[]

	i = pred[i-node_id];

	//while (pred[i-node_id] >= node_id)

	while (pred_ext[i-node_id] >= node_id)

	{

		//path.push_back(pred[i-node_id]);

		//i = pred[i-node_id];

		path.push_back(pred_ext[i-node_id]);

		i = pred_ext[i-node_id];

	}

	path_dist += dist_before_fst_node; //add the weight before the source node (so path_dist is the total weight of the entire group)



	//optimization: check the list of source_nodes, see if there's any other source node 

	//that is already included in the current shortest path. If that partial path is shorter,

	//include that partial path (both paths in the output, current shortest path will be cut down 

	//later, in post-processing output step); otherwise remove that source node (only current path 

	//in the output)

	map<int, float>::iterator mapIt;

	int k;

	for (i=0; i<path.size()-1; i++) //check all nodes on path except the current source node

	{

		mapIt = source_nodes.find(path[i]);

		if (mapIt != source_nodes.end())

		{

#ifdef DEBUG_VERSION

			dataMgr.outFile << "another source node contained in current shortest path: " << path[i] << "\n";

#endif

			d = dist[path_end-node_id] - dist_ext[path[i]-node_id] + (*mapIt).second;

			if (d < path_dist) //better than current path

			{

#ifdef DEBUG_VERSION

				dataMgr.outFile << "shorter path: ";

#endif

				vector<int> partial_path;

				for (k=0; k<=i; k++)

				{

#ifdef DEBUG_VERSION

					dataMgr.outFile << path[k] << " - ";

#endif

					partial_path.push_back(path[k]);

				}

#ifdef DEBUG_VERSION

				dataMgr.outFile << "\n";

#endif

				local_short_paths.insert(multimap<float, vector<int> >::value_type(d, partial_path));

			}



			source_nodes.erase(mapIt); //there is no need to compute for that source node again

		}

	}



	delete [] dist;

	delete [] pred;

//	delete [] cut;

	delete [] dist_ext;

	delete [] pred_ext;



	return path_dist;

}



//THE MAIN FUNCTION

void HSP_graph::GeneGrouping()

{



	dataMgr.PrepareData(isPosStrand, chr_index);


//#ifdef VERBOSE
	if (VERBOSE)
	cout << "prepare data done" << "\n";
//#endif


//	dataMgr.PrintHSPs(chr_index);



	if (isPosStrand)

	{

		//cout << "pos strand ";

		num_nodes = dataMgr.HSP_gene[chr_index].size();

		//cout << num_nodes << "\n";

		if (num_nodes == 2)

			return;

		

		CreateAllExtCutEdges(dataMgr.HSP_gene[chr_index]);

		//cout << "extcut edges done" << "\n";



#ifdef DEBUG
		PrintEdges(debugFile);
#endif



		CreateAllSkipEdges(dataMgr.HSP_gene[chr_index]);

		//cout << "skip edges done" << "\n";

#ifdef DEBUG
		PrintEdges(debugFile);

#endif



	}

	else

	{

		//cout << "neg strand ";

		num_nodes = dataMgr.HSP_neg_gene[chr_index].size();

		//cout << num_nodes << "\n";

		if (num_nodes == 2)

			return;



		CreateAllExtCutEdges(dataMgr.HSP_neg_gene[chr_index]);

		//cout << "extcut edges done" << "\n";

#ifdef DEBUG
		PrintEdges(debugFile);
#endif



		CreateAllSkipEdges(dataMgr.HSP_neg_gene[chr_index]);

		//cout << "skip edges done" << "\n";

#ifdef DEBUG
		PrintEdges(debugFile);
#endif



	}



	//cout << "starting ShortestPath()..." << "\n";



	//if (strcmp(PHASE1_VERSION, "1.1") == 0)

	if (PHASE1_VERSION.compare("1.1") == 0)

		FindAllLocalShortestPaths(); //v1.1

	else

		ShortestPath(); //v1.0



	//cout << "shortest path done" << "\n";

}



//output edges, for debug only

void HSP_graph::PrintEdges(ostream& os)

{

	os << "******EDGES******" << "\n";

	int i, j = edges.size();

	for ( i=cur_edge; i<j; i++)

	{

		if (isPosStrand)

			os << "source HSP_ID_" << dataMgr.HSP_gene[chr_index][edges[i].source_num].ID 

			<< " dest HSP_ID_" << dataMgr.HSP_gene[chr_index][edges[i].dest_num].ID;

		else

			os << "source HSP_ID_" << dataMgr.HSP_neg_gene[chr_index][edges[i].source_num].ID 

			<< " dest HSP_ID_" << dataMgr.HSP_neg_gene[chr_index][edges[i].dest_num].ID;



		os << " " << edges[i];

	}

	cur_edge = i;



}



//****************************************************************************************

//a special comparison function used only for sorting edges in "CreateEdges2" when creating edges within single HSP

bool less_edge_dest(const Base_Edge & e1, const Base_Edge & e2) {

        return e1.dest_num < e2.dest_num || (e1.dest_num == e2.dest_num && e1.source_num < e2.source_num);

}





//****************************************************************************************



bool ACCP_DONR_graph::IsDonorAcceptorEmpty()

{

	if (donors.empty())

	{

		acceptors.clear();

		donor_region_index.clear();

		acceptor_region_index.clear();

		return true;

	}

	else

	{

		if (acceptors.empty())

		{

			donors.clear();

			donor_region_index.clear();

			acceptor_region_index.clear();

			return true;

		}

		else

			return false;

	}



}






//modified: simply start searching from hsp_chr_end of last HSP, until hit "stop"

int ACCP_DONR_graph::CalcEndPos(int hsp_chr_end, int hsp_gene_end, bool isPosStrand, 
								//pair<int, char*>& chr_seq)
								vector<string>& chr_seq)

{

	int cur_pos, i, end;

	string cur_frame;

	bool stop=false;



	int max_pos = chromosome_start_pos-1 + 
		//chr_seq.first 
		LenOfStrVec(chr_seq)
		- 2; //make sure there's at least 2 chars following "cur_pos" in "chr_seq"

	if (isPosStrand)

	{

//		cur_pos = hsp_chr_end + 3*(dataMgr.query_len - hsp_gene_end) - 2;		

//		if (cur_pos > max_pos)

//			cur_pos = max_pos;

//		end = cur_pos + 2; //initial value

		//cur_pos = hsp_chr_end + 3*(dataMgr.query_len - hsp_gene_end) + 1;

		cur_pos = hsp_chr_end + 1;

		if (cur_pos > max_pos)
		{
			if (cur_pos > max_pos + 2)
				return cur_pos -1 - 3*((cur_pos - 1 - max_pos)/3);
			else
				return cur_pos-1; //this is the end
		}

		end = cur_pos - 1;

	}

	else

	{

		//cur_pos = hsp_chr_end - 3*(dataMgr.query_len - hsp_gene_end); //here "hsp_chr_end" is actually the hsp_chr_start

		//if (cur_pos < 1) //minimum 1

		//	cur_pos = 1;

		//end = cur_pos; //initial value

		//cur_pos = hsp_chr_end - 3*(dataMgr.query_len - hsp_gene_end)-3;

		cur_pos = hsp_chr_end - 3;

		//if (cur_pos < 1)
		if (cur_pos < chromosome_start_pos)
			//return cur_pos + 3*((3-cur_pos)/3);
			return chromosome_start_pos + (3-(chromosome_start_pos - cur_pos)%3)%3;

		end = cur_pos + 3;

	}



#ifdef DEBUG
	dataMgr.outFile << "searching end from " << end << "\n";
	dataMgr.outFile << "(hsp_chr_end at " << hsp_chr_end << "; hsp_gene_end at " << hsp_gene_end << ")" << "\n";
	dataMgr.outFile << "max_pos: " << max_pos << "\n";
#endif



	while (cur_pos > chromosome_start_pos-1 //0
		&& cur_pos <= max_pos && !stop) //go on from "end" position until hitting a stop codon or end of sequence

	{

		//cur_frame = chr_seq.substr(cur_pos-1, 3);

		GetSubstrFromVecStrs(chr_seq, isPosStrand, cur_pos-1, 3, cur_frame);

#ifdef DEBUG

		dataMgr.outFile << "cur_frame: " << cur_frame << "\n";
		dataMgr.outFile << "cur_pos: " << cur_pos << "\n";

#endif

		StrToLower(cur_frame);

		if (isPosStrand) //positive strand

		{

			for (i=0; i<3; i++)

			{

				if (cur_frame.compare(STOP_CODONS[i]) == 0)

				{

					//end = cur_pos - 1; //without stop codon! //cur_pos+2;
					end = cur_pos + 2;

					stop = true;

					break;

				}

			}



			cur_pos += 3;

		}

		else //negative strand

		{

			//dataMgr.outFile << cur_frame;

			for (i=0; i<3; i++)

			{

				if (cur_frame.compare(STOP_CODONS_REV[i])==0)

				{

					//end = cur_pos + 3; //cur_pos;
					end = cur_pos;

					stop = true;

					break;

				}

			}



			cur_pos -= 3;

		}



	}

//	dataMgr.outFile << "\n";



/*	int vecIndex=-1;

	int leftover, cur_index;

	int index;

	if (!isPosStrand)

		cur_pos += 2;

	while ( !stop && (vecIndex = GetSingleStrFromVecStrs(chr_seq, isPosStrand, vecIndex, cur_pos-1, cur_frame, leftover)) != -2)

	{

		if (isPosStrand) //positive strand

		{

			index = INT_MAX;

			bool cur_stop;

			for (i=0; i<3; i++)

			{

				cur_stop = false;

				cur_index = 0;

				while ((cur_index = cur_frame.find(STOP_CODONS[i], cur_index)) != string::npos && !cur_stop)

				{

					if (cur_index % 3 == 0)

					{

						if (cur_index<index)

							index = cur_index;



						cur_stop = true;

					}

					else

						cur_index++;

				}

			}





			if (cur_stop)

			{

				end = cur_pos + index - 1;

				stop = true;			

			}



			cur_pos += cur_frame.length();

		}

		else //negative strand

		{

			//dataMgr.outFile << cur_frame << "\n";



			index = -1;

			bool cur_stop;

			for (i=0; i<3; i++)

			{

				cur_stop = false;

				cur_index = string::npos;

				while ((cur_index = cur_frame.rfind(STOP_CODONS_REV[i], cur_index)) != string::npos && !cur_stop)

				{

					if (cur_index % 3 == 0)

					{

						if (cur_index>index)

							index = cur_index;

						cur_stop = true;

					}

					else

						cur_index--;

				}

			}



			cur_pos -= cur_frame.length();



			if (cur_stop)

			{

				end = cur_pos + 1 + index + 3;

				stop = true;			

			}



		}





	}

*/
#ifdef DEBUG
	dataMgr.outFile << "end_site: " << end << "\n";
#endif

	return end;

}



int ACCP_DONR_graph::CalcStartPos(int hsp_chr_start, int hsp_gene_start, bool isPosStrand, 
								  //pair<int, char*>& chr_seq) 
								  vector<string>& chr_seq)

{

	int cur_pos, start, i;

	string cur_frame;

#ifdef DEBUG
	dataMgr.outFile << "dump chr_seq:\n";
	for (i=0; i<chr_seq.size(); i++)
	{
		dataMgr.outFile << i << ":" << chr_seq[i].length() << "\n";
		dataMgr.outFile << chr_seq[i] << "\n";
	}
#endif

	int max_pos = chromosome_start_pos-1 + 
		//chr_seq.first 
		LenOfStrVec(chr_seq)
		- 2;

	int stop_pos=-1; //signals no stop

	if (isPosStrand)

	{

		cur_pos = hsp_chr_start-3*(hsp_gene_start-1); //the supposed mark "1" position

		//if (cur_pos < 1)
		if (cur_pos < chromosome_start_pos)
			cur_pos = chromosome_start_pos + ((3-(chromosome_start_pos - cur_pos)%3))%3; //1; //modified: make sure it's also in frame

#ifdef DEBUG

		dataMgr.outFile << "start: cur_pos:" << cur_pos << "; chr_seq size:" << chr_seq.size() 
			<< "; hsp_chr_start:" << hsp_chr_start << "\n";
#endif


		//check stop codon between "1" and hsp_chr_start

		GetSubstrFromVecStrs(chr_seq, isPosStrand, cur_pos-1, hsp_chr_start - cur_pos, cur_frame);

#ifdef DEBUG

		dataMgr.outFile << "cur_frame:" << cur_frame << "\n";

#endif
		StrToLower(cur_frame);

#ifdef DEBUG

		dataMgr.outFile << "cur_frame:" << cur_frame << "\n";

#endif
		int cur_stop_pos;

		for (i=0; i<3; i++)

		{

			cur_stop_pos = -1;

			while ((cur_stop_pos = cur_frame.find(STOP_CODONS[i], cur_stop_pos+1)) != string::npos)

				if ((cur_frame.length() - cur_stop_pos ) % 3 == 0)

				{

					if (stop_pos < (cur_pos + cur_stop_pos))

						stop_pos = cur_pos + cur_stop_pos;

				}

		}

#ifdef DEBUG

		dataMgr.outFile << "cur_stop_pos:" << cur_stop_pos << ";stop_pos:" << stop_pos << "\n";

#endif


		if (stop_pos != -1)

			cur_pos = stop_pos + 3;

	}

	else

	{

		cur_pos = hsp_chr_start+3*(hsp_gene_start-1)-2; //here "hsp_chr_start" is actually the hsp_chr_end

#ifdef DEBUG

		dataMgr.outFile << "start: cur_pos: " << cur_pos << "\n";

#endif

		if (cur_pos > max_pos)

			cur_pos = max_pos - (3-(cur_pos - max_pos)%3)%3; //modified: make sure it's also in frame



		GetSubstrFromVecStrs(chr_seq, isPosStrand, hsp_chr_start, cur_pos + 2 - hsp_chr_start, cur_frame);

		StrToLower(cur_frame);

		int cur_stop_pos;

		for (i=0; i<3; i++)

		{

			cur_stop_pos = -1;

			while ((cur_stop_pos = cur_frame.find(STOP_CODONS_REV[i], cur_stop_pos+1)) != string::npos)

			{

				if (cur_stop_pos % 3 == 0)

				{

					if (stop_pos == -1 || stop_pos > hsp_chr_start+1+cur_stop_pos)

					{

						stop_pos = hsp_chr_start+1+cur_stop_pos;

						break;

					}

				}

			}

		}

#ifdef DEBUG

		dataMgr.outFile << "stop_pos: " << stop_pos << "\n";

#endif


		if (stop_pos != -1)
			cur_pos = stop_pos - 3;


#ifdef DEBUG

		dataMgr.outFile << "cur_pos: " << cur_pos << "\n";

#endif

	}

	int mark_pos = cur_pos; //this is the mark "1" position


#ifdef DEBUG
		dataMgr.outFile << "mark_pos:" << mark_pos << "\n";
#endif



//	if (cur_pos > 0) //this must be the case for negative strand?

//	{

		//cur_frame=chr_seq.substr(cur_pos-1, 3);

		GetSubstrFromVecStrs(chr_seq, isPosStrand, cur_pos-1, 3, cur_frame);
#ifdef DEBUG
		dataMgr.outFile << "cur_frame:" << cur_frame << "\n";
#endif

		StrToLower(cur_frame);

		if (isPosStrand && (cur_frame.compare(START_CODON) == 0)) //found at the exact position

		{

			start = cur_pos;
#ifdef DEBUG
		dataMgr.outFile << "cur_frame:" << cur_frame << "\n";
#endif

		}

		else

		if (!isPosStrand && (cur_frame.compare(START_CODON_REV) == 0)) //found at the exact position

		{

			start = cur_pos + 2;

		}

		else //not found at exact position "1"

		{

			bool found=false;

			bool stop=false;

			if (isPosStrand)

				cur_pos -= 3;

			else

				cur_pos += 3;

			while (cur_pos > chromosome_start_pos-1 //0
				&& cur_pos <= max_pos && !stop && !found) //walk thru frames before position "1", until hit on either Start/Stop codons 

			{					

				//cur_frame=chr_seq.substr(cur_pos-1, 3);

				GetSubstrFromVecStrs(chr_seq, isPosStrand, cur_pos-1, 3, cur_frame);

				StrToLower(cur_frame);
#ifdef DEBUG
		dataMgr.outFile << "cur_frame:" << cur_frame << "\n";
#endif



				if (isPosStrand) //positive strand

				{

					if (cur_frame.compare(START_CODON) == 0) //hit start codon

					{

						found = true;

						start = cur_pos;
#ifdef DEBUG
		dataMgr.outFile << "cur_frame:" << cur_frame << "\n";
#endif

					}

					else

					{

						for (i=0; i<3; i++)

						{					

							if (cur_frame.compare(STOP_CODONS[i])==0) //hit stop codon

							{

								stop = true;

								break;

							}

						}

					}



					cur_pos -= 3;

				}

				else //negative strand

				{

					if (cur_frame.compare(START_CODON_REV) == 0)//hit start codon

					{

						found = true;

						start = cur_pos+2;

					}

					else

					{

						for (i=0; i<3; i++)

						{

							if (cur_frame.compare(STOP_CODONS_REV[i])==0) //hit stop codon

							{

								stop = true;

								break;

							}

						}

					}



					cur_pos += 3;

				}

					

			}



			if (!found) //didn't find start codon, now walk thru frames between "1" and hsp_chr_start

			{

				if (hsp_gene_start>1)

				{

					if (isPosStrand) //positive strand

					{

						//cur_pos = hsp_chr_start-3*(hsp_gene_start-2);

						cur_pos = mark_pos + 3;

						while (cur_pos < hsp_chr_start + START_CODON_CHECK_LEN)

						{

							//cur_frame = chr_seq.substr(cur_pos-1, 3);

							GetSubstrFromVecStrs(chr_seq, isPosStrand, cur_pos-1, 3, cur_frame);

							StrToLower(cur_frame);

							if (cur_frame.compare(START_CODON)==0)

								break;



							cur_pos += 3;

						}

						start = cur_pos;

					}

					else //negative strand

					{

						//cur_pos = hsp_chr_start+3*(hsp_gene_start-2)-2;

						cur_pos = mark_pos - 3;

						while (cur_pos > hsp_chr_start - START_CODON_CHECK_LEN)

						{

							//cur_frame = chr_seq.substr(cur_pos-1, 3);

							GetSubstrFromVecStrs(chr_seq, isPosStrand, cur_pos-1, 3, cur_frame);

							StrToLower(cur_frame);

							if (cur_frame.compare(START_CODON_REV)==0)

								break;



							cur_pos -= 3;

						}

						if (cur_pos < hsp_chr_start - START_CODON_CHECK_LEN)

							start = hsp_chr_start;

						else

							start = cur_pos+2;

					}

				}

				else

				{

					start = hsp_chr_start;

				}



			}



		}

//	}

/*	else //nothing before "1", (this can only happen for positive strand) now check only between "1" and hsp_chr_start

	{

		cur_pos = hsp_chr_start;

		start = cur_pos;

		while (cur_pos > 0 )

		{

			//cur_frame = chr_seq.substr(cur_pos-1, 3);

			GetSubstrFromVecStrs(chr_seq, isPosStrand, cur_pos-1, 3, cur_frame);

			StrToLower(cur_frame);

			if (cur_frame.compare(START_CODON)==0)

				start = cur_pos;



			cur_pos -= 3;

		}

	}

*/
#ifdef DEBUG
	dataMgr.outFile << "start(position only, no strandness):" << start << "\n";
#endif

	return start;



}





void ACCP_DONR_graph::Reset()

{

	num_nodes = 0;

	edges.clear();


	edges_w_score.clear();

}








//this also group HSPs based on frames
void ACCP_DONR_graph::MergeOverlapHSPs(bool isPosStrand, vector<HSP_Gene_Pair*>& HSPs)

{

	//initialize flags

	in_frame_merge = false;

	out_frame_merge = false;



	int i;



	vector<HSP_Gene_Pair*>::iterator ptrHSPIt, tmpIt;

	int j;



	vector<HSP_Gene_Pair*> HSPs_to_be_erased;


//	bool remove_lastHSP; //only used in TestMerge() -> Merge()


	if (isPosStrand)

	{

		do { //2 pass, first pass do in-frame merge, second pass do out-frame merge

		for (i=0; i<3; i++)

			HSP_frame[i].clear();



		ptrHSPIt = HSPs.end();

		ptrHSPIt--;



		while (ptrHSPIt != HSPs.begin())
		{

			//remove_lastHSP = false; //reset before every call to TestMerge()

			if (!TestMerge(ptrHSPIt, j, true, HSPs_to_be_erased)) //, remove_lastHSP))
			{

/*				if (remove_lastHSP) //if we need to remove lastHSP (lastHSP must be in frame with curHSP, so also in frame j)

				{

					tmpIt = ptrHSPIt;

					tmpIt++;

					HSPs.erase(tmpIt);

					//also remove lastHSP from HSP_frame!

					//(for in-frame merge, lastHSP must be HSP_frame[j].back())

					HSP_frame[j].pop_back();



#ifdef DEBUG

					dataMgr.outFile << "lastHSP removed" << "\n";

#endif

				}
*/
				HSP_frame[j].push_back(*ptrHSPIt);
				ptrHSPIt--;
			}
			else
			{
#ifdef DEBUG

					dataMgr.outFile << "curHSP removed" << "\n";
					dataMgr.outFile << "curHSP:" << *(*ptrHSPIt) << "\n";

#endif


				tmpIt = ptrHSPIt;

				tmpIt--;

#ifdef DEBUG
				dataMgr.outFile << "curHSP hold(" << *ptrHSPIt << ") to be erased:" << *(*ptrHSPIt) << "\n";
				dataMgr.outFile << "tmpIt hold:" << *tmpIt << "\n" << *(*tmpIt) << "\n";
#endif
				HSPs.erase(ptrHSPIt);

				ptrHSPIt = tmpIt;

			}

			//dataMgr.outFile << "all HSPs at this point: " << "\n";

			//vector<HSP_Gene_Pair*>::iterator all_hsp_it = HSPs.begin();

			//for (; all_hsp_it != HSPs.end(); all_hsp_it++)

			//	dataMgr.outFile << *(*all_hsp_it) << "\n";			

		}



		//now ptrHSPIt points to the first HSP

		//remove_lastHSP = false; //reset remove_lastHSP
		if (!TestMerge(ptrHSPIt, j, true, HSPs_to_be_erased)) //, remove_lastHSP))

		{

/*			if (remove_lastHSP)

			{

				tmpIt = ptrHSPIt;

				tmpIt++;

				HSPs.erase(tmpIt);

				HSP_frame[j].pop_back();

#ifdef DEBUG

					dataMgr.outFile << "lastHSP removed" << "\n";

#endif

			}

*/			

			HSP_frame[j].push_back(*ptrHSPIt);

		}
		else

		{
			HSPs.erase(ptrHSPIt);

#ifdef DEBUG

			dataMgr.outFile << "curHSP removed" << "\n";

#endif

		}

		for (ptrHSPIt = HSPs_to_be_erased.begin(); ptrHSPIt != HSPs_to_be_erased.end(); ptrHSPIt++)

			HSPs.erase(remove(HSPs.begin(), HSPs.end(), *ptrHSPIt), HSPs.end());



		if (!in_frame_merge)

			in_frame_merge = true; //in_frame_merge is done!

		else

			out_frame_merge = false; //in_frame_merge already true, must be second round already, set flag so we can exit

		} while (out_frame_merge); //out_frame_merge may have been updated by "TestMerge()"



/*		if (out_frame_merge)

		{

			map<int, vector<int> > hsp_end_hsp_ids;

			for ( ptrHSPIt = HSPs.begin(); ptrHSPIt != HSPs.end(); ptrHSPIt++)

			{

				int cur_hsp_start = (*ptrHSPIt)->HSP_start;

				int cur_hsp_end = (*ptrHSPIt)->HSP_end;

				MyMapInsert(cur_hsp_start, (*ptrHSPIt)->ID, true, hsp_end_hsp_ids);

				MyMapInsert(cur_hsp_end, (*ptrHSPIt)->ID, false, hsp_end_hsp_ids);

			}



			map<int, vector<int> >::iterator hsp_ids_it, hsp_ids_next_it;

			vector<int>::iterator ids_it;

			map<int, Input_Alignment>::iterator alignIt;



			hsp_ids_it = hsp_end_hsp_ids.begin();

			hsp_ids_next_it = hsp_end_hsp_ids.begin();

			hsp_ids_next_it++;

			while (hsp_ids_next_it != hsp_end_hsp_ids.end())

			{

				int cur_region_start = (*hsp_ids_it).first;

				int cur_region_end = (*hsp_ids_next_it).first;

				AD_Edge_Penalty cur_region_penalty(0, INT_MAX, 0);

				int cur_region_hsp_id=-1;

				for (ids_it = (*hsp_ids_it).second.begin(); ids_it != (*hsp_ids_it).second.end(); ids_it++)

				{

					alignIt = dataMgr.input_alignments.find(*ids_it);

					string& cur_targetStr = (*alignIt).second.target_align;

					string& cur_matchStr = (*alignIt).second.match_align;

					string& cur_queryStr = (*alignIt).second.query_align;



					int cur_search_start_pos, cur_search_end_pos, cur_extra_front, cur_extra_end;

					if (!FindRealPos(cur_targetStr, HSPs[*ids_it]->HSP_start, HSPs[*ids_it]->HSP_end,  

						cur_region_start, cur_region_end-1, 

						cur_search_start_pos, cur_search_end_pos, cur_extra_front, cur_extra_end))

						continue;



					int tmp_exact=0, tmp_gap=0, tmp_positive=0;

					GetPenaltyFromMatchStr(cur_search_start_pos, cur_search_end_pos, cur_extra_front, cur_extra_end, 

						cur_matchStr, tmp_exact, tmp_gap, tmp_positive);

					AD_Edge_Penalty tmp_cur_penalty(tmp_exact, tmp_gap, tmp_positive);



					if (tmp_cur_penalty < cur_region_penalty)

					{

						cur_region_penalty = tmp_cur_penalty;

						cur_region_hsp_id = *ids_it;

					}

				}

				hsp_ids_it++;

				hsp_ids_next_it++;

			}

		}

*/

	}

	else

	{

		do {

		for (i=0; i<3; i++)

			HSP_frame[i].clear();



		HSP_Gene_Pair* hsp_to_be_deleted = NULL;



		ptrHSPIt = HSPs.begin();

		while (ptrHSPIt != HSPs.end())

		{

			//remove_lastHSP = false; //reset remove_lastHSP

			if (!TestMerge(ptrHSPIt, j, false, HSPs_to_be_erased)) //, remove_lastHSP))

			{

/*				if (remove_lastHSP)

				{

					//tmpIt = ptrHSPIt;

					//tmpIt--;

					//HSPs.erase(tmpIt); //cannot erase here, because it will invalidate the ptrHSPIt iterator!

					HSPs_to_be_erased.push_back(HSP_frame[j].back());

					HSP_frame[j].pop_back();

#ifdef DEBUG

					dataMgr.outFile << "lastHSP will removed" << "\n";

#endif

				}

*/

				HSP_frame[j].push_back(*ptrHSPIt);

			}
			else
			{
				(*ptrHSPIt) = hsp_to_be_deleted; //mark this to be deleted
#ifdef DEBUG

				dataMgr.outFile << "curHSP will be removed" << "\n";

#endif

			}



			ptrHSPIt++;

		}



		HSPs.erase(remove(HSPs.begin(), HSPs.end(), hsp_to_be_deleted), HSPs.end()); //delete all "0"s

		

		for (ptrHSPIt = HSPs_to_be_erased.begin(); ptrHSPIt != HSPs_to_be_erased.end(); ptrHSPIt++)

			HSPs.erase(remove(HSPs.begin(), HSPs.end(), *ptrHSPIt), HSPs.end());



		if (!in_frame_merge)

			in_frame_merge = true; //in_frame_merge is done!

		else

			out_frame_merge = false; //in_frame_merge already true, must be second round already, reset flag so we can exit

		} while (out_frame_merge);

	}



	//in_frame_merge = false; //reset flag



/*	for (i=0; i<3; i++)

	{

		cout << "HSP_frame[" << i << "]:";

		for (ptrHSPIt = HSP_frame[i].begin(); ptrHSPIt != HSP_frame[i].end(); ptrHSPIt++)

			cout << *(*ptrHSPIt) << "\n";

	}

*/

}



bool ACCP_DONR_graph::Merge(HSP_Gene_Pair* lastHSP, HSP_Gene_Pair* curHSP, int search_start, int cur_search_end, int last_search_end, 

							vector<HSP_Gene_Pair*>& HSPs_to_be_erased)

							//bool& remove_lastHSP)
{

	//MODIFIED: if last_hsp is embeded in curHSP, remove lastHSP

	int last_hsp_start_coord;

	if (search_start > 0)

		last_hsp_start_coord = lastHSP->HSP_start;

	else

		last_hsp_start_coord = -lastHSP->HSP_end;

	if (last_hsp_start_coord > search_start)

	{

		dataMgr.outFile << "lastHSP should be removed" << "\n";

		HSPs_to_be_erased.push_back(lastHSP);

		return false;

	}



	int cur_id = curHSP->ID;



	int last_id = lastHSP->ID;



#ifdef DEBUG
	dataMgr.outFile << "resolving lastHSP [" << last_id << "]" << *lastHSP << " and curHSP [" 
		<< cur_id << "]" << *curHSP << "\n";
#endif

	map<int, Input_Alignment>::iterator alignIt = dataMgr.input_alignments.find(cur_id);



//	dataMgr.outFile << (*alignIt).second << "\n";

	string& cur_targetStr = (*alignIt).second.target_align;

	string& cur_matchStr = (*alignIt).second.match_align;

	string& cur_queryStr = (*alignIt).second.query_align;



	int cur_search_start_pos, cur_search_end_pos, cur_extra_front, cur_extra_end;

	//if (!FindRealPos(cur_targetStr, cur_start, cur_end,  cur_start, last_end, 

	//	cur_search_start_pos, cur_search_end_pos, cur_extra_front, cur_extra_end))

	if (search_start > 0)

	{

		if (!FindRealPos(cur_targetStr, curHSP->HSP_start, curHSP->HSP_end,  search_start, cur_search_end, 

			cur_search_start_pos, cur_search_end_pos, cur_extra_front, cur_extra_end))

		{

			cout << "something wrong1? " << "\n";

			exit(-1);
			//return false;

		}

	}

	else

	{

		if (!FindRealPos(cur_targetStr, -curHSP->HSP_end, -curHSP->HSP_start, search_start, cur_search_end, 

			cur_search_start_pos, cur_search_end_pos, cur_extra_front, cur_extra_end))

		{
			cout << "something wrong2?\n" << cur_targetStr << "\nhsp_start:" << -curHSP->HSP_end << ";hsp_end:" << -curHSP->HSP_start 

				<< ";search_start:" << search_start << ";search_end:" << cur_search_end << "\n";

			exit(-1);

			//return false;

		}

	}





	int tmp_exact=0, tmp_gap=0, tmp_positive=0;

	GetPenaltyFromMatchStr(cur_search_start_pos, cur_search_end_pos, cur_extra_front, cur_extra_end, 

		cur_matchStr, tmp_exact, tmp_gap, tmp_positive);

	AD_Edge_Penalty tmp_cur_penalty(tmp_exact, tmp_gap, tmp_positive);



	alignIt = dataMgr.input_alignments.find(last_id);

	string& last_targetStr = (*alignIt).second.target_align;

	string& last_matchStr = (*alignIt).second.match_align;

	string& last_queryStr = (*alignIt).second.query_align;



#ifdef DEBUG

	dataMgr.outFile << "before resolving:" << "\n";

	dataMgr.outFile << "last HSP: " << "\n";

	dataMgr.outFile << last_targetStr << "\n";

	dataMgr.outFile << last_matchStr << "\n";

	dataMgr.outFile << last_queryStr << "\n";

	dataMgr.outFile << "cur HSP: " << "\n";

	dataMgr.outFile << cur_targetStr << "\n";

	dataMgr.outFile << cur_matchStr << "\n";

	dataMgr.outFile << cur_queryStr << "\n";

#endif



	int last_search_start_pos, last_search_end_pos, last_extra_front, last_extra_end;

	//if (!FindRealPos(last_targetStr, last_start, last_end, cur_start, last_end, 

	//	last_search_start_pos, last_search_end_pos, last_extra_front, last_extra_end))

	if (search_start > 0)

	{

		if (!FindRealPos(last_targetStr, lastHSP->HSP_start, lastHSP->HSP_end, 

			search_start,

			last_search_end, 

			last_search_start_pos, last_search_end_pos, last_extra_front, last_extra_end))

		{
			cout << "something wrong3? " << "\n";

			exit(-1);

			//return false;

		}

	}

	else

	{

		if (!FindRealPos(last_targetStr,-lastHSP->HSP_end, -lastHSP->HSP_start,  

			search_start, 

			last_search_end, 

			last_search_start_pos, last_search_end_pos, last_extra_front, last_extra_end))

		{
			cout << "something wrong4? search_start:" << search_start << ", last_search_end:" << last_search_end 

				<< ", last_search_start_pos:" << last_search_start_pos << ", last_search_end_pos:" 

				<< last_search_end_pos << ", last_extra_front:" << last_extra_front 

				<< ", last_extra_end:" << last_extra_end << "\n";

			exit(-1);

			//return false;

		}

	}



	tmp_exact=tmp_gap=tmp_positive=0;

	GetPenaltyFromMatchStr(last_search_start_pos, last_search_end_pos, last_extra_front, last_extra_end, 

		last_matchStr, tmp_exact, tmp_gap, tmp_positive);

	AD_Edge_Penalty tmp_last_penalty(tmp_exact, tmp_gap, tmp_positive);

	if (in_frame_merge) //in_frame_merge has been done, doing out frame merge

	//because it's out-frame, two HSPs must have different starts/ends! So no need to worry about invalid situations after cut-off
	{

		string tmp_queryStr;

		int gap_add, len_advance;

		if (tmp_cur_penalty < tmp_last_penalty) //cut off last HSP

		{
#ifdef DEBUG
			dataMgr.outFile << "cut off portion of lastHSP...";
#endif


				if (search_start > 0) //positive strand

				{

					gap_add = (search_start - lastHSP->HSP_start) % 3 ; //make sure frame match!!!

					lastHSP->HSP_end = search_start-1- gap_add; 
//#ifdef VERBOSE
					if (VERBOSE)
					cout << "last start:" << lastHSP->HSP_start << "; end:" << lastHSP->HSP_end << "; pos:" << search_start 
						<< "; new_hsp_end:" << lastHSP->HSP_end << "\n";
//#endif
#ifdef DEBUG
					dataMgr.outFile << "last start:" << lastHSP->HSP_start << "; end:" << lastHSP->HSP_end << "; pos:" << search_start 
						<< "; new_hsp_end:" << lastHSP->HSP_end << "\n";
#endif


				}

				else //negative strand

				{

					//gap_add = (cur_search_end + lastHSP->HSP_end) % 3;//i.e. (cur_search_end - (-lastHSP->HSP_end)) % 3

					//lastHSP->HSP_start = -cur_search_end + 1 + gap_add;

					gap_add = (search_start + lastHSP->HSP_end) % 3;

					lastHSP->HSP_start = -search_start + 1 + gap_add;
//#ifdef VERBOSE
					if (VERBOSE)
					cout << "last start:" << -lastHSP->HSP_end << "; end:" << -lastHSP->HSP_start << "; pos_start:" << search_start 
						<< "; pos_end:" << cur_search_end << "; new_hsp_start:" << lastHSP->HSP_start << "\n";
//#endif
#ifdef DEBUG
					dataMgr.outFile << "last start:" << -lastHSP->HSP_end << "; end:" << -lastHSP->HSP_start << "; pos_start:" << search_start 
						<< "; pos_end:" << cur_search_end << "; new_hsp_start:" << lastHSP->HSP_start << "\n";
#endif


				}

				if (lastHSP->HSP_end < lastHSP->HSP_start) //invalid after cut, then delete

				{
#ifdef DEBUG
					dataMgr.outFile << "lastHSP should be removed" << "\n";
#endif
					HSPs_to_be_erased.push_back(lastHSP);

					return false;

				}



				tmp_queryStr = last_queryStr.substr(0, last_search_start_pos-1);

				lastHSP->gene_end = lastHSP->gene_start + last_search_start_pos-1 - CountCharInStr(tmp_queryStr, '-') - 1;

				//dataMgr.outFile << "outframe merge: lastHSP gene_end is updated to " << lastHSP->gene_end << "\n";



				const_cast< string& >(last_targetStr).erase(last_search_start_pos-1);

				const_cast< string& >(last_matchStr).erase(last_search_start_pos-1);

				const_cast< string& >(last_queryStr).erase(last_search_start_pos-1);



#ifdef DEBUG

				dataMgr.outFile << last_targetStr << "\n";

				dataMgr.outFile << last_matchStr << "\n";

				dataMgr.outFile << last_queryStr << "\n";

#endif

				return false;

		}

		else

		{
#ifdef DEBUG
			dataMgr.outFile << "cut off portion of curHSP...";
#endif


			if (cur_search_end == last_search_end)//cut off curHSP

			{

				if (search_start > 0)

				{

					gap_add = ( curHSP->HSP_end - last_search_end ) % 3;

					curHSP->HSP_start = last_search_end+1 + gap_add;
//#ifdef VERBOSE
					if (VERBOSE)
					cout << "cur start:" << curHSP->HSP_start << "; end:" << curHSP->HSP_end << "; pos:" << last_search_end
						<< "; new_hsp_start:" << curHSP->HSP_start << "\n";
//#endif
#ifdef DEBUG
					dataMgr.outFile << "cur start:" << curHSP->HSP_start << "; end:" << curHSP->HSP_end << "; pos:" << last_search_end
						<< "; new_hsp_start:" << curHSP->HSP_start << "\n";
#endif
				}

				else

				{

					gap_add = (-last_search_end - curHSP->HSP_start) % 3;

					curHSP->HSP_end = -last_search_end - 1 - gap_add;
//#ifdef VERBOSE
					if (VERBOSE)
					cout << "cur start:" << -curHSP->HSP_end << "; end:" << -curHSP->HSP_start << "; pos:" << last_search_end 
						<< "; new_hsp_end:" << curHSP->HSP_end << "\n";
//#endif
#ifdef DEBUG
					dataMgr.outFile << "cur start:" << -curHSP->HSP_end << "; end:" << -curHSP->HSP_start << "; pos:" << last_search_end 
						<< "; new_hsp_end:" << curHSP->HSP_end << "\n";
#endif
				}

				if (curHSP->HSP_end < curHSP->HSP_start) //invalid after cut, then delete

				{
#ifdef DEBUG
					dataMgr.outFile << "curHSP should be removed" << "\n";
#endif
					return true;

				}


				if (cur_extra_end > 0) //has extra_end, then need to advance 1 more (non-gapped positions)

				{

					len_advance = cur_search_end_pos + 2;

					while (cur_targetStr[len_advance-1] == '-')

						len_advance++;

				}
				else

				{
					len_advance = cur_search_end_pos + 1;

				}

				tmp_queryStr = cur_queryStr.substr(0, len_advance);

				curHSP->gene_start += len_advance - CountCharInStr(tmp_queryStr, '-'); //cur_search_end_pos - CountCharInStr(tmp_queryStr, '-');

				//dataMgr.outFile << "outframe merge: curHSP gene_start is updated to " << curHSP->gene_start << "\n";



				const_cast< string& >(cur_targetStr).erase(0, len_advance);//cur_search_end_pos);

				const_cast< string& >(cur_matchStr).erase(0, len_advance);//cur_search_end_pos);

				const_cast< string& >(cur_queryStr).erase(0, len_advance);//cur_search_end_pos);



#ifdef DEBUG

				dataMgr.outFile << cur_targetStr << "\n";

				dataMgr.outFile << cur_matchStr << "\n";

				dataMgr.outFile << cur_queryStr << "\n";

#endif

				return false;

			}

			else //just remove curHSP

			{

				return true;

			}

		}

	}

	else //in_frame_merge //updated: keep both HSPs! just make them into disjoint HSPs

	//in-frame merge may result in invalid start/end after cut-off, if HSPs' start/end were the same before
	{

	string tmp_queryStr;

	//update HSP_end and gene_end, no need to update pid since it's not used in computation later

	//(does it update the HSP_Gene_Pair it points to?)

/*	if (search_start > 0) //isPosStrand

	{

		lastHSP->HSP_end = curHSP->HSP_end; //cur_search_end;

		//dataMgr.outFile << "inframe merge: lastHSP hsp_end is updated to " << lastHSP->HSP_end << "\n";

	}

	else

	{

		lastHSP->HSP_start = curHSP->HSP_start; //-cur_search_end;

		//dataMgr.outFile << "inframe merge: lastHSP hsp_start is updated to " << lastHSP->HSP_start << "\n";

	}

	lastHSP->gene_end = curHSP->gene_end;

*/	//dataMgr.outFile << "inframe merge: lastHSP gene_end is updated to " << lastHSP->gene_end << "\n";



	if ( tmp_cur_penalty < tmp_last_penalty) //merge by replacing last portion of previous HSP by first portion of current HSP
	{

#ifdef DEBUG

		dataMgr.outFile << "merge by cutting off last portion of previous HSP, curHSP unchanged" << "\n";

#endif


		if (search_start > 0) //isPosStrand
			lastHSP->HSP_end = curHSP->HSP_start-1;
		else
			lastHSP->HSP_start = curHSP->HSP_end + 1;



		if (lastHSP->HSP_start > lastHSP->HSP_end) //invalid lastHSP start/end, must remove lastHSP altogether

		{

			//remove_lastHSP = true; //set up flag to remove lastHSP
#ifdef DEBUG
			dataMgr.outFile << "lastHSP needs to be removed" << "\n";
#endif
			HSPs_to_be_erased.push_back(lastHSP);

			return false; //do not remove or change curHSP

		}

		tmp_queryStr = last_queryStr.substr(0, last_search_start_pos);
		lastHSP->gene_end = lastHSP->gene_start + last_search_start_pos - CountCharInStr(tmp_queryStr, '-') - 1;

		//update alignment strings
		const_cast< string& >(last_targetStr).erase(last_search_start_pos);
		//const_cast< string& >(last_targetStr).append(cur_targetStr);
		const_cast< string& >(last_matchStr).erase(last_search_start_pos);
		//const_cast< string& >(last_matchStr).append(cur_matchStr);
		const_cast< string& >(last_queryStr).erase(last_search_start_pos);
		//const_cast< string& >(last_queryStr).append(cur_queryStr);

		//dataMgr.outFile << "replacing portion of last HSP by portion of cur HSP" << "\n";

	}
	else
	{
		//const_cast< string& >(last_targetStr).append(cur_targetStr.substr(cur_search_end_pos+1));
		//const_cast< string& >(last_matchStr).append(cur_matchStr.substr(cur_search_end_pos+1));
		//const_cast< string& >(last_queryStr).append(cur_queryStr.substr(cur_search_end_pos+1));
		//dataMgr.outFile << "appending portion of cur HSP to last HSP" << "\n";


#ifdef DEBUG

		dataMgr.outFile << "merge by cutting off first portion of current HSP, lastHSP unchanged" << "\n";

#endif


		if (search_start > 0)
			curHSP->HSP_start = lastHSP->HSP_end+1;
		else
			curHSP->HSP_end = lastHSP->HSP_start-1;



		//when both prev and cur HSP ends at same position, the new hsp_start/hsp_end will be invalid, must remove!

		if (curHSP->HSP_start > curHSP->HSP_end)

		{
#ifdef DEBUG
			dataMgr.outFile << "curHSP needs to be removed" << "\n";
#endif
			return true; //remove current HSP

		}

		tmp_queryStr = cur_queryStr.substr(0, cur_search_end_pos+1);

		curHSP->gene_start += cur_search_end_pos+1 - CountCharInStr(tmp_queryStr, '-');

		const_cast< string& >(cur_targetStr).erase(0, cur_search_end_pos+1);

		const_cast< string& >(cur_matchStr).erase(0, cur_search_end_pos+1);

		const_cast< string& >(cur_queryStr).erase(0, cur_search_end_pos+1);

	}



#ifdef DEBUG

	dataMgr.outFile << "last HSP: " << "\n";

	dataMgr.outFile << last_targetStr << "\n";

	dataMgr.outFile << last_matchStr << "\n";

	dataMgr.outFile << last_queryStr << "\n";

	dataMgr.outFile << "cur HSP: " << "\n";

	dataMgr.outFile << cur_targetStr << "\n";

	dataMgr.outFile << cur_matchStr << "\n";

	dataMgr.outFile << cur_queryStr << "\n";

#endif



	//return true;

	return false;

	}

}







bool ACCP_DONR_graph::TestMerge(vector<HSP_Gene_Pair*>::iterator ptrHSPIt, int& j, bool isPosStrand, 

								vector<HSP_Gene_Pair*>& HSPs_to_be_erased) //bool& remove_lastHSP)

{

/*	int k;

	for (k=0; k<3; k++)

	{

		dataMgr.outFile << "frame " << k << ":" << "\n";

		for (int kk=0; kk<HSP_frame[k].size(); kk++)

			dataMgr.outFile << *(HSP_frame[k][kk]) << "\n";

	}

*/

	HSP_Gene_Pair* lastHSP;


	int cur_start, cur_end, last_start, last_end;



	if (isPosStrand)

	{

		cur_start = (*ptrHSPIt)->HSP_start;

		cur_end = (*ptrHSPIt)->HSP_end;

	}

	else

	{

		cur_start = -(*ptrHSPIt)->HSP_end;

		cur_end = -(*ptrHSPIt)->HSP_start;

	}



	j= ((cur_start - start_site) % 3 + 3) % 3; //make sure j is non-negative (must be 0, 1, 2)!!!



	bool merged = false;

	int of[2];

	int i;


	vector<HSP_Gene_Pair*>::reverse_iterator hsp_it;


	if (!in_frame_merge) //if we are in the process of doing in frame merge

	{

	//dataMgr.outFile << "trying inframe_merge: " << cur_start << "-" << cur_end << "\n";



	if (!out_frame_merge) //not yet set, then check if out_frame_merge should be set

	{

	//check other frames, if overlap but in different frame, choose only one and delete the other one? keep partial HSPs?

	switch (j)

	{

	case 0:  {  of[0] = 1; of[1] = 2;  break;   }

	case 1:  {  of[0] = 0; of[1] = 2;  break;   }

	case 2:  {  of[0] = 0; of[1] = 1;  break;   }

	default: break;

	}

	for (i=0; i<2; i++)

	{

		if (!HSP_frame[of[i]].empty())

		{

			lastHSP = HSP_frame[of[i]].back();

			if (isPosStrand)

			{

				last_start = lastHSP->HSP_start;

				last_end = lastHSP->HSP_end;

			}

			else

			{

				last_start = -lastHSP->HSP_end;

				last_end = -lastHSP->HSP_start;

			}



			if (cur_start <= last_end) //overlap

			{

				//compare and select only one/remove other one?
#ifdef DEBUG
				dataMgr.outFile << "Overlapping HSP not in same frame..." << "\n";
#endif
				//cout << "Overlapping HSP not in same frame?! ..." << "\n";

				//exit(-1);

				out_frame_merge = true;

			}



		}





	}

	}



	//now doing its in_frame_merge

	//modified: must check more than just the very last one in HSP_frame[j], because there may be cases 

	//where curHSP overlaps with more than more previous HSPs!

	if (HSP_frame[j].empty())

		return false;



	//lastHSP = HSP_frame[j].back();



	hsp_it = HSP_frame[j].rbegin();

	while ( hsp_it != HSP_frame[j].rend() )

	{

		lastHSP = *hsp_it;

		if (isPosStrand)

		{

			last_start = lastHSP->HSP_start;

			last_end = lastHSP->HSP_end;

		}

		else

		{

			last_start = -lastHSP->HSP_end;

			last_end = -lastHSP->HSP_start;

		}



#ifdef DEBUG

		dataMgr.outFile << "cur_start: " << cur_start << "; last_end: " << last_end << "\n";

#endif



		if (cur_start > last_end) //no overlap, stop immediately

		{

#ifdef DEBUG

			dataMgr.outFile << "curHSP finished merging" << "\n";

#endif

			return false;

		}


	//if (cur_start <= last_end) //overlap, also in same frame, so merge
	//{

		//modified: cut off curHSP head portion if it's sticking out before last_start (because lastHSP has been cut off)

		//updated: now hold the original info of curHSP, in case it needs to be restored

		int old_cur_hsp_start = cur_start;

		int old_cur_gene_start = (*ptrHSPIt)->gene_start;

		string tmp_queryStr, tmp_matchStr, tmp_targetStr;

		map<int, Input_Alignment>::iterator alignIt;

		if (cur_start < last_start)

		{			

			cur_start = last_start;



			if (cur_end < cur_start) //invalid start/end, then remove curHSP (this should not happen?)

				return true;



			//change hsp_start (for positive strand) or hsp_end (for negative strand)

			if (isPosStrand)

				(*ptrHSPIt)->HSP_start = cur_start;

			else

				(*ptrHSPIt)->HSP_end = -cur_start;



			//change the alignment strings and gene_start

			int cur_id = (*ptrHSPIt)->ID;

			alignIt = dataMgr.input_alignments.find(cur_id);

			string& cur_targetStr = (*alignIt).second.target_align;

			string& cur_matchStr = (*alignIt).second.match_align;

			string& cur_queryStr = (*alignIt).second.query_align;



			int cur_search_start_pos, cur_search_end_pos, cur_extra_front, cur_extra_end; //because in frame, no extra_front or extra_end

			FindRealPos(cur_targetStr, old_cur_hsp_start, cur_end, cur_start, cur_start, cur_search_start_pos, cur_search_end_pos, 

				cur_extra_front, cur_extra_end);

			

			tmp_queryStr = cur_queryStr.substr(0, cur_search_start_pos);

			tmp_matchStr = cur_matchStr.substr(0, cur_search_start_pos);

			tmp_targetStr = cur_targetStr.substr(0, cur_search_start_pos);

			

			const_cast< string& >(cur_targetStr).erase(0, cur_search_start_pos);

			const_cast< string& >(cur_matchStr).erase(0, cur_search_start_pos);



			//change gene_start (for both positive and negative)
			//const string& fStr = cur_queryStr.substr(0, cur_search_start_pos);
			//(*ptrHSPIt)->gene_start += fStr.length() - CountCharInStr(fStr, ' ');

			(*ptrHSPIt)->gene_start += cur_search_start_pos - CountCharInStr(tmp_queryStr, '-');



			const_cast< string& >(cur_queryStr).erase(0, cur_search_start_pos);

		}



		if (last_end <= cur_end) //overlap region [cur_start, last_end]
		{
			//cout << "here1" << "\n";
			merged = Merge(lastHSP, *ptrHSPIt, cur_start, last_end, last_end, HSPs_to_be_erased); //remove_lastHSP);
		}
		else
		{

			//somthing wrong? (this shouldn't happen given how we get phase-1 HSPs)
			/*dataMgr.outFile << "something wrong???" << "\n";
			cout << "cur_start:" << cur_start << "<=last_end:" << last_end 
				<< "; && cur_end:" << cur_end << "<last_end:" << last_end << "! something wrong???" << "\n";
			exit(-1);*/
			return true; //remove curHSP
		}



		if (!merged)

		{

			if (cur_start != old_cur_hsp_start //curHSP has been cut off before Merge() and Merge() didn't change curHSP

			&& ( (isPosStrand && (*ptrHSPIt)->HSP_start == cur_start) || (!isPosStrand && (*ptrHSPIt)->HSP_end == -cur_start)))

			//restore curHSP back

			{

			if (isPosStrand)

				(*ptrHSPIt)->HSP_start = old_cur_hsp_start;

			else

				(*ptrHSPIt)->HSP_end = -old_cur_hsp_start;



			(*ptrHSPIt)->gene_start = old_cur_gene_start;



			string& cur_targetStr = (*alignIt).second.target_align;

			string& cur_matchStr = (*alignIt).second.match_align;

			string& cur_queryStr = (*alignIt).second.query_align;

			const_cast< string& >(cur_targetStr).insert(0, tmp_targetStr);

			const_cast< string& >(cur_matchStr).insert(0, tmp_matchStr);

			const_cast< string& >(cur_queryStr).insert(0, tmp_queryStr);



			cur_start = old_cur_hsp_start; //reset (in case curHSP was cut then restored, then need to check against next lastHSP)


#ifdef DEBUG
			dataMgr.outFile << "restore curHSP back:" << "\n";
			//dataMgr.outFile << *(*ptrHSPIt) << "\n";
#endif
			}

		}

		else //if curHSP needs to be removed, stop and return immediately!

			return true;

	//}

		hsp_it++;


		//if lastHSP is to be deleted, also remove it from HSP_frame[j], 

		//but first advance hsp_it so it won't affect that iterator

		if (!HSPs_to_be_erased.empty())

		{

			if (HSPs_to_be_erased.back() == lastHSP)

				HSP_frame[j].pop_back();

		}



#ifdef DEBUG

		for (vector<HSP_Gene_Pair*>::iterator debug_it = HSP_frame[j].begin(); debug_it != HSP_frame[j].end(); debug_it++)

			dataMgr.outFile << "all hsps in HSP_frame: " << *(*debug_it) << "\n";

		if (hsp_it != HSP_frame[j].rend())

			dataMgr.outFile << "lastHSP to be checked next: " << *(*hsp_it) << "\n";

		else

			dataMgr.outFile << "lastHSP ptr at end of list" << "\n";

#endif



	} //end while()



		return merged;

	} //end if()


	else //in frame merge has been finished, we must be doing out frame merge

	{

		//dataMgr.outFile << "doing outframe_merge: " << cur_start << "-" << cur_end << "\n";



		switch (j)

		{

		case 0:  {  of[0] = 1; of[1] = 2;  break;   }

		case 1:  {  of[0] = 0; of[1] = 2;  break;   }

		case 2:  {  of[0] = 0; of[1] = 1;  break;   }

		default: break;

		}

		for (i=0; i<2; i++)
		{

//NEED_CHECK:	if (!HSP_frame[of[i]].empty())

			hsp_it = HSP_frame[of[i]].rbegin();

			while (hsp_it != HSP_frame[of[i]].rend())
			{
				//lastHSP = HSP_frame[of[i]].back();

				lastHSP = *hsp_it;

				if (isPosStrand)
				{
					last_start = lastHSP->HSP_start;
					last_end = lastHSP->HSP_end;
				}
				else
				{
					last_start = -lastHSP->HSP_end;
					last_end = -lastHSP->HSP_start;
				}

				//get the cur_start and cur_end again, since curHSP (*ptrHSPIt) coordinates may have been updated by previous merging
				if (isPosStrand)
				{
					cur_start = (*ptrHSPIt)->HSP_start;
					cur_end = (*ptrHSPIt)->HSP_end;
				}
				else
				{
					cur_start = -(*ptrHSPIt)->HSP_end;
					cur_end = -(*ptrHSPIt)->HSP_start;
				}

				if (cur_start > last_end) //do not overlap
					break;

				if (cur_start < last_start && cur_end > last_end)
				{
					HSPs_to_be_erased.push_back(lastHSP); //remove lastHSP!
#ifdef DEBUG

					dataMgr.outFile << "lastHSP will be removed" << "\n";

#endif

					//vector<HSP_Gene_Pair*>::iterator to_be_erased = HSP_frame[of[i]].end();
					//to_be_erased--;
					//HSP_frame[of[i]].erase(to_be_erased);



					hsp_it++;

				

					HSP_frame[of[i]].pop_back();

					merged = false;

					//goto NEED_CHECK;

					continue;
				}
				else

				{
					if ( last_start < cur_start && last_end > cur_end)
					{
						//merged = true; //remove curHSP!

						return true; //return immediately

					}
					else

					{
						if ((cur_start <= last_end && cur_start >= last_start )
							|| (cur_end <= last_end && cur_end >= last_start) )//overlap
						{
							//compare and select only one/remove other one? (cut off appropriate portion if applicable)

							if (cur_end > last_end)
							{
								//cout << "here2" << "\n";
								//merged = Merge_Out_Frame((*ptrHSPIt), lastHSP, cur_start, last_end, last_end);
								merged = Merge(lastHSP, *ptrHSPIt, cur_start, last_end, last_end, HSPs_to_be_erased); //remove_lastHSP);
							}
							else
							{
								//cout << "here3" << "\n";
								//merged = Merge_Out_Frame((*ptrHSPIt), lastHSP, cur_start, cur_end, last_end);
								//merged = Merge(lastHSP, *ptrHSPIt, cur_start, cur_end, last_end);
								merged = Merge( *ptrHSPIt, lastHSP, last_start, cur_end, cur_end, HSPs_to_be_erased); //remove_lastHSP); //lastHSP actually has position after curHSP
							}

							if (merged)//if curHSP is already deemed redundant, no need to check further
								return true;//return immediately

							if (!HSPs_to_be_erased.empty())
							{
								if (HSPs_to_be_erased.back() == lastHSP)
								{
									hsp_it++;

									HSP_frame[of[i]].pop_back();

									continue;

								}

							}

						}
					}

				} 



				hsp_it++;
			} //end while()

			if (merged == true) //if curHSP is already deemed redundant, no need to check further
				break;

		} //end for()



		return merged;



	}



}





//the only difference between CalcStartEndPos and CalcStartEndPos2 is "LoadData" or not

void ACCP_DONR_graph::CalcStartEndPos2(bool isPosStrand, vector<HSP_Gene_Pair*>& HSPs, 
									   //pair<int, char*>& chr_seq)
									   vector<string>& chr_seq)

{

	//compute start/end positon

	int hsp_gene_start, hsp_gene_end, hsp_chr_start, hsp_chr_end;

	hsp_chr_start = HSPs.back()->HSP_start;

	hsp_chr_end = HSPs.front()->HSP_end;

	if (isPosStrand)

	{

		hsp_gene_start = HSPs.back()->gene_start;

		hsp_gene_end = HSPs.front()->gene_end;



		//start = CalcStartPos(hsp_chr_start, hsp_gene_start, true, chr_seq);

		//end = CalcEndPos(hsp_chr_end, hsp_gene_end, true, chr_seq);


//#ifdef VERBOSE
		if (VERBOSE)
		cout << "start pos" << "\n";
//#endif
		start_site = CalcStartPos(hsp_chr_start, hsp_gene_start, true, chr_seq); //(*chrMapIt).second);
//#ifdef VERBOSE
		if (VERBOSE)
		cout << "end pos" << "\n";
//#endif
		end_site = CalcEndPos(hsp_chr_end, hsp_gene_end, true, chr_seq); //(*chrMapIt).second);

	}

	else

	{

		hsp_gene_start = HSPs.front()->gene_start;

		hsp_gene_end = HSPs.back()->gene_end;



		//start = - CalcStartPos(hsp_chr_end, hsp_gene_start, false, chr_seq); //negative integer

		//end = - CalcEndPos(hsp_chr_start, hsp_gene_end, false, chr_seq);  //negative integer


//#ifdef VERBOSE
		if (VERBOSE)
		cout << "start pos" << "\n";
//#endif
		start_site = - CalcStartPos(hsp_chr_end, hsp_gene_start, false, chr_seq); //(*chrMapIt).second); //negative integer
//#ifdef VERBOSE
		if (VERBOSE)
		cout << "end pos" << "\n";
//#endif
		end_site = - CalcEndPos(hsp_chr_start, hsp_gene_end, false, chr_seq); //(*chrMapIt).second);  //negative integer

	}



}







void ACCP_DONR_graph::PrintHeader(int rank, int count, int alt_count, 
								  multimap<Group_Info, vector<HSP_Gene_Pair*> >::iterator groupMapIt)
{
		float per_cover = (float)(*groupMapIt).first.gene_cover_len/(float)dataMgr.query_len*100;

/*		dataMgr.gff_os << dataMgr.query_gene << "-Rank" << rank 
				<< ";Name=" << dataMgr.query_gene << ";Rank=" << rank << "_" << count 
				<< ";Cover=" << (*groupMapIt).first.gene_cover_len 
				<< "(" << per_cover << "%)\n";
*/
		sprintf(dataMgr.cur_gene_id, "%s-R%d-%d", dataMgr.query_gene.c_str(), rank, count); //for use in gff, pro, cDNA, raw output

		if (OUTPUT_GFF)
		{
			dataMgr.gff_gene_str << dataMgr.cur_chr_name  << "\tgenBlastG" << USER_ID << "\ttranscript\t";
			dataMgr.gff_geneinfo_str << "\t" << (*groupMapIt).first.score;
			if ((*groupMapIt).first.isPosStrand)
				dataMgr.gff_geneinfo_str << "\t+\t.\tID=";
			else
				dataMgr.gff_geneinfo_str << "\t-\t.\tID=";
			//dataMgr.gff_os << dataMgr.query_gene << "-R" << rank << "-" << count << "-A" << alt_count << "\n";			
			dataMgr.gff_geneinfo_str << dataMgr.cur_gene_id << "-A" << alt_count 
				<< ";Name=" << dataMgr.query_gene;
			if (GENBLASTG_NEED_PID)
				dataMgr.gff_geneinfo_str << ";PID=";
			else
				dataMgr.gff_geneinfo_str << "\n";
				//<< "\n";
		}

		dataMgr.outFile << (*groupMapIt).first.gene_cover_len << "(" << per_cover 
			<< "%)|score:" << (*groupMapIt).first.score << "|rank:" << rank << "\n";

		vector<HSP_Gene_Pair*>::iterator vecIt;
		for (vecIt = (*groupMapIt).second.begin(); vecIt != (*groupMapIt).second.end(); vecIt++)
			dataMgr.outFile << "HSP\t" << (*vecIt)->ID << "\t" << (*vecIt)->HSP_start << "-" << (*vecIt)->HSP_end 
			<< "\t[" << (*vecIt)->gene_start << "-" << (*vecIt)->gene_end << "]" << "\n";
}





















bool ACCP_DONR_graph::GetSingleHSPDonorSegments_SplitOL_GapIsMinIntron(int cur_hsp_start, int& cur_hsp_end, int hsp_ID, 



		//multimap<int, pair<int, int> >& donor_segments, 

		int& last_segment_start, int& last_segment_hsp_ID, 

		int next_hsp_ID, int& next_hsp_start, int cur_query_start, int& cur_query_end, 

		int& next_query_start, int next_query_end, int& last_segment_end, bool& cur_hsp_exists)

//[last_segment_start, last_segment_end] records the last "exon" segment from previous HSP

{

//	if (cur_hsp_start == -6746886 && cur_hsp_end == -6746869)

//		int stophere=1;



	int next_qStart = next_query_start;

	int cur_qEnd = cur_query_end;



	//find alignment

	map<int, Input_Alignment>::iterator alignIt = dataMgr.input_alignments.find(hsp_ID);

	string& matchStr = (*alignIt).second.match_align;

	string& targetStr = (*alignIt).second.target_align;

	string& queryStr = (*alignIt).second.query_align;



	int cur_len_ori = targetStr.length();

	

	vector<int> tgt_gap_starts, tgt_gap_ends;

	GetGaps(targetStr, tgt_gap_starts, tgt_gap_ends);



#ifdef DEBUG



	dataMgr.outFile << "queryStr: " << queryStr << "\n";

	dataMgr.outFile << "cur_hsp_start:" << cur_hsp_start << "; cur_hsp_end:" << cur_hsp_end 

		<< "cur_query_start:" << cur_query_start << "; cur_query_end:" << cur_query_end 

		<< "; last_segment_start:" << last_segment_start << "; last_segment_end:" << last_segment_end 

		<< "; next_hsp_start:" << next_hsp_start 

		<< "; next_query_start:" << next_query_start << "; next_query_end:" << next_query_end << "\n";



#endif



	int cur_tgt_pos, next_tgt_pos;



	//cur_hsp_exists = true;

	bool next_hsp_exists = true;

	map<int, Input_Alignment>::iterator next_alignIt;

	if (next_hsp_ID != -1) //has next HSP

	{		

		next_alignIt = dataMgr.input_alignments.find(next_hsp_ID); //next HSP

		string& nextQueryStr = (*next_alignIt).second.query_align; //nest HSP query string

		string& nextMatchStr = (*next_alignIt).second.match_align;

		string& nextTargetStr = (*next_alignIt).second.target_align;

	

		int next_len_ori = nextTargetStr.length();



		if (next_query_start <= cur_query_end) //check whether they overlap on query

		{

			vector<int> next_tgt_gap_starts, next_tgt_gap_ends;

			GetGaps(nextTargetStr, next_tgt_gap_starts, next_tgt_gap_ends);



			string hsps_align;

			int cur_splice_end, next_splice_start;

			ComputeSpliceAlignment_SplitOL((*alignIt).second, (*next_alignIt).second, cur_query_start, cur_query_end,

				next_query_start, next_query_end, hsps_align, cur_splice_end, next_splice_start, 

				cur_qEnd, next_qStart); //change the alignments themselves!!!



#ifdef DEBUG

			dataMgr.outFile << "splice align:\n" << hsps_align << "\n";

			dataMgr.outFile << "cur_splice_end: " << cur_splice_end << "; next_splice_start: " << next_splice_start << "\n";

			dataMgr.outFile << "possible new cur_qEnd:" << cur_qEnd << "; possible new next_qStart:" << next_qStart << "\n";

			dataMgr.outFile << "get segments: \nquery:" << queryStr << "\n" << "match:" << matchStr << "\n" << "targt:" << targetStr << "\n";

			dataMgr.outFile << "query:" << nextQueryStr << "\n" << "match:" << nextMatchStr << "\n" << "targt:" << nextTargetStr << "\n";

#endif

			if (cur_splice_end == -1) //curHSP is all cut off (nextHSP must be unchanged)

			{

				cur_hsp_exists = false;

				next_tgt_pos = ConvertToNoGappedPos(next_splice_start, next_tgt_gap_starts, next_tgt_gap_ends, true);



				if (last_segment_start == 0)

				{

					if (!donor_segments_pair.empty())

					{

						//multimap<int, pair<int, int> >::reverse_iterator seg_it = donor_segments.rbegin();

						//const_cast<int&>((*seg_it).second.second) = next_hsp_start + next_tgt_pos*3 - 1;

						donor_segments_pair[donor_segments_pair.size()-1].intron_seg_end = next_hsp_start + next_tgt_pos*3 - 1;

						//vector<int>::reverse_iterator gap_it = gap_centers.rbegin();

						//gap_centers.back() = ( (*seg_it).second.second + (*seg_it).first )/2;



#ifdef DEBUG

						dataMgr.outFile << "segment end updated: " << next_hsp_start + next_tgt_pos*3 - 1  << "\n";

#endif

					}

				}

				else

				{

					if (last_segment_end == 0)

					{

						//donor_segments.insert(multimap<int, pair<int, int> >::value_type(cur_hsp_start, 

						//	pair<int, int>(last_segment_start, next_hsp_start + next_tgt_pos*3 - 1)));

						donor_segments_pair.push_back( SegmentsInThreeBounds(last_segment_start, cur_hsp_start, next_hsp_start + next_tgt_pos*3 - 1) );

						//gap_centers.push_back((cur_hsp_start + next_hsp_start + next_tgt_pos*3-1)/2);



#ifdef DEBUG

						dataMgr.outFile << "segment inserted: " << last_segment_start << " -> " << cur_hsp_start 

							<< "->" << next_hsp_start + next_tgt_pos*3 - 1 << "\n";

#endif

					}

					else

					{

						//donor_segments.insert(multimap<int, pair<int, int> >::value_type(last_segment_end, 

						//	pair<int, int>(last_segment_start, next_hsp_start + next_tgt_pos*3 - 1)));

						donor_segments_pair.push_back( SegmentsInThreeBounds(last_segment_start, last_segment_end, next_hsp_start + next_tgt_pos*3 - 1) );

						//gap_centers.push_back( (last_segment_end + next_hsp_start + next_tgt_pos*3-1)/2);



#ifdef DEBUG

						dataMgr.outFile << "segment inserted: " << last_segment_start << "->" << last_segment_end 

							<< "->" << next_hsp_start + next_tgt_pos*3 - 1 << "\n";

#endif

						last_segment_end = 0;  //reset

					}

					last_segment_start = 0; //reset



					//donor_acceptor_HSP_ID.push_back(last_segment_hsp_ID);

				}



				return next_hsp_exists;

			}

			else //curHSP still exists

			{

				if (next_splice_start == nextQueryStr.length()) //nextHSP is cut off

				{

					next_hsp_exists = false;

					//last_segment_end = cur_hsp_start + tgt_start_pos*3 + 3;

				}				

				//else //otherwise both HSPs are left, then their alignments already updated by ComputeSpliceAlignment, so just do regular stuff

				cur_tgt_pos = ConvertToNoGappedPos(cur_splice_end, tgt_gap_starts, tgt_gap_ends, false);



			if (next_qStart != 0)

				next_query_start = next_qStart;

			if (cur_qEnd != 0)

				cur_query_end = cur_qEnd;



			if (cur_splice_end < cur_len_ori - 1)

			{

				cur_hsp_end = cur_hsp_start + (cur_tgt_pos+1)*3 - 1;

				tgt_gap_starts.clear();

				tgt_gap_ends.clear();

				GetGaps(targetStr, tgt_gap_starts, tgt_gap_ends); //after splice, targetStr is changed, so redo GetGaps

			}

			if (next_splice_start > 0 && next_splice_start < next_len_ori)

				next_hsp_start += ConvertToNoGappedPos(next_splice_start, next_tgt_gap_starts, next_tgt_gap_ends, true) * 3;



#ifdef DEBUG

			dataMgr.outFile << "now cur_query_end is: " << cur_query_end << "; next_query_start is: " << next_query_start << "\n";

#endif

			if (next_query_start > next_query_end)

				next_hsp_exists = false;

			if (cur_query_start > cur_query_end)

				cur_hsp_exists = false;



			if (!cur_hsp_exists)

				return next_hsp_exists;

			}

		}

		else //otherwise HSPs do not overlap on query, do regular stuff

		{

			//cur_tgt_pos = targetStr.length() - CountCharInStr(targetStr, '-') - 1; //no gapped index of last char in curHSP's targetStr

		}

	}

	//else //no next HSP, do regular stuff



	//int i;

	//string gapStr="";

	//for (i=0; i<MIN_INTRON_LEN; i++)

	//	gapStr += '-';

	//int min_intron_len_aa = MIN_INTRON_LEN / 3;//convert from DNA length to a.a. length

	string gapStr(MIN_INTRON_LEN_AA, '-'); 

	int start_ori;

	int start_pos = 0;

	if (last_segment_start == 0)

	{

		start_ori = cur_hsp_start;

		start_pos =queryStr.find_first_not_of("-"); //start from the first non-gap position

		start_ori += start_pos*3;

		//tgt_start_pos = ConvertToNoGappedPos(start_pos, tgt_gap_starts, tgt_gap_ends); //no need here, since first gaps on query must correspond to all non-gap on target

		//start_ori += tgt_start_pos*3;

		donor_acceptor_HSP_ID.push_back( vector<int>(1, hsp_ID) );

#ifdef DEBUG

		PrintDonorAcceptorHSPID();

#endif

	}

	else

	{

		if (last_segment_end == 0)

		{

			start_ori = last_segment_start; //last_segment_start is not zero, curHSP must start with non-gap, i.e. start_pos=0

			last_segment_start = 0; //reset



			if (!donor_acceptor_HSP_ID.empty())

			{

				if ( donor_acceptor_HSP_ID.back().back() != hsp_ID)

					donor_acceptor_HSP_ID.back().push_back(hsp_ID);

			}

			else

			{

				donor_acceptor_HSP_ID.push_back( vector<int>(1, hsp_ID) ); //is it possible?

				dataMgr.outFile << "damn: why is it possible? hsp_ID " << hsp_ID << "\n"; exit(-1);

			}

#ifdef DEBUG

		PrintDonorAcceptorHSPID();

#endif

		}

		else

		{

			start_ori = cur_hsp_start;

			start_pos = queryStr.find_first_not_of("-");

			start_ori += start_pos*3;

			//donor_segments.insert(multimap<int, pair<int, int> >::value_type(last_segment_end, 

			//	pair<int, int>(last_segment_start, start_ori-1)));

			donor_segments_pair.push_back( SegmentsInThreeBounds(last_segment_start,last_segment_end, start_ori-1) );

			//gap_centers.push_back( (last_segment_end + start_ori - 1)/2 );

			//donor_acceptor_HSP_ID.push_back(last_segment_hsp_ID);

			donor_acceptor_HSP_ID.push_back( vector<int>(1, hsp_ID) );



#ifdef DEBUG

			PrintDonorAcceptorHSPID();

			dataMgr.outFile << "start segment: " << last_segment_start << "->" << last_segment_end << "->" << start_ori-1 << "\n";

#endif

			last_segment_start = 0;

			last_segment_end = 0;

		}

	}

	

	//start_ori += start_pos*3;

#ifdef DEBUG

	dataMgr.outFile << "start_ori: " << start_ori << " (start_pos:" << start_pos << ")" << "\n";

#endif



	int tgt_start_pos, tgt_end_pos;



	int end_pos=0;

	while ((start_pos = queryStr.find(gapStr, start_pos)) != string::npos) //found a gap that is at least MIN_INTRON_LEN long

	{
#ifdef DEBUG
		dataMgr.outFile << "start_pos: " << start_pos;
#endif


		//end_pos cannot be string::npos, because HSP cannot end with gaps!?

		//NOT after merging/cutting off part of HSPs due to overlap merging!

		if (start_pos+MIN_INTRON_LEN_AA < queryStr.length())

		{

			end_pos = queryStr.find_first_not_of("-", start_pos+MIN_INTRON_LEN_AA); //look for non-gap position after the current gap

			
#ifdef DEBUG
			dataMgr.outFile << "; end_pos: " << end_pos << "\n";
#endif


			if (end_pos != string::npos) //normal segment, record it

			{

				donor_acceptor_HSP_ID.push_back( vector<int>(1, hsp_ID) );



				tgt_start_pos = ConvertToNoGappedPos(start_pos, tgt_gap_starts, tgt_gap_ends, false);

				tgt_end_pos = ConvertToNoGappedPos(end_pos, tgt_gap_starts, tgt_gap_ends, false);

				

				//donor_segments.insert(multimap<int, pair<int, int> >::value_type(cur_hsp_start+tgt_start_pos*3, //cur_hsp_start + start_pos*3, //cur_hsp_start+mid_pos, 

				//	pair<int, int>(start_ori, cur_hsp_start+tgt_end_pos*3-1)));//cur_hsp_start + end_pos*3 - 1)));

				donor_segments_pair.push_back( SegmentsInThreeBounds(start_ori, cur_hsp_start+tgt_start_pos*3, cur_hsp_start+tgt_end_pos*3-1) );

				//gap_centers.push_back((cur_hsp_start + tgt_end_pos*3 + cur_hsp_start + tgt_start_pos*3)/2);

				//donor_acceptor_HSP_ID.push_back(hsp_ID);



#ifdef DEBUG

		PrintDonorAcceptorHSPID();

				dataMgr.outFile << "segment: " << start_ori << "<-" << cur_hsp_start + tgt_start_pos*3 << "->" << cur_hsp_start + tgt_end_pos*3 - 1 << "\n";

#endif			

				//start_ori += end_pos*3;

				start_ori = cur_hsp_start + tgt_end_pos*3; //end_pos*3;

				start_pos = end_pos;

			}

			else //current HSP ends with gap

			{

				end_pos = queryStr.length();

				break;

			}

		}

		else //current HSP ends with gap (must due to merging and there must be another HSP immediately following it)

		{

			end_pos = queryStr.length();

			break;

		}

	}

	





	if (next_hsp_ID != -1) //has next HSP

	{

		last_segment_start = start_ori;

		last_segment_end = cur_hsp_end+1;//cur_hsp_start + cur_tgt_pos*3 + 3;



		if (next_hsp_exists)

		{

			string& nextQueryStr = (*next_alignIt).second.query_align; //next HSP query string

			int next_start_pos = nextQueryStr.find_first_not_of("-"); //this should be the same as next_tgt_pos?

#ifdef DEBUG

			dataMgr.outFile << "next_start_pos: " << next_start_pos << "\n";

#endif

/*

#ifdef V23_ALL_HSP_GAP_SPLICE_SEGMENT

		if (end_pos < queryStr.length()) //does not end with gap, record last_segment, for use with next HSP segment checking

		{

			if (next_hsp_start != cur_hsp_end + 1 || next_start_pos > 0)//next HSP is not immediate to cur HSP, or next HSP start with a gap

			//if (next_hsp_start + next_start_pos*3 > cur_hsp_end + MIN_INTRON_LEN) //the gap between cur_hsp_end and next HSP's first non-gap position is larger than MIN_INTRON_LEN

			{

				//donor_segments.insert(multimap<int, pair<int, int> >::value_type(last_segment_end, //cur_hsp_end+1,

				//	pair<int, int>(start_ori, next_hsp_start+next_start_pos*3-1)));

				donor_segments_pair.push_back( SegmentsInThreeBounds(start_ori, last_segment_end, next_hsp_start+next_start_pos*3-1) );

				//gap_centers.push_back((next_hsp_start+next_start_pos*3-1+last_segment_end)/2);

				donor_acceptor_HSP_ID.push_back(hsp_ID);



#ifdef DEBUG

			dataMgr.outFile << "segment1: " << start_ori << "<-" << last_segment_end << "->" << next_hsp_start+next_start_pos*3-1 << "\n";

#endif			

				last_segment_start = 0;

				last_segment_end = 0;

			}

			else //next HSP follows immediately and starts with non-gap, record the partial segment to extend into next HSP

			{

				//last_segment_start = start_ori;

				last_segment_end = 0;

			}

		}

		else //end with gap, then look for next_HSP non-gap start position

		{

			//if (next_hsp_start + next_start_pos*3 > cur_hsp_start + start_pos*3 -1 + MIN_INTRON_LEN)

			//{

			tgt_start_pos = ConvertToNoGappedPos(start_pos, tgt_gap_starts, tgt_gap_ends); //same as cur_tgt_pos?

#ifdef DEBUG

			dataMgr.outFile << "tgt_start_pos: " << tgt_start_pos << "\n";

#endif

			//donor_segments.insert(multimap<int, pair<int, int> >::value_type(cur_hsp_start + tgt_start_pos*3, //cur_hsp_start+mid_pos, 

			//	pair<int, int>(start_ori, next_hsp_start+next_start_pos*3-1)));

			donor_segments_pair.push_back( SegmentsInThreeBounds(start_ori, cur_hsp_start + tgt_start_pos*3, next_hsp_start+next_start_pos*3-1) );

			//gap_centers.push_back((next_hsp_start+next_start_pos*3+cur_hsp_start+tgt_start_pos*3)/2);

			donor_acceptor_HSP_ID.push_back(hsp_ID);



#ifdef DEBUG

			dataMgr.outFile << "segment2: " << start_ori << "<-" << cur_hsp_start + tgt_start_pos*3 << "->" << next_hsp_start+next_start_pos*3-1 << "\n";

#endif			

			last_segment_start = 0;

			last_segment_end = 0;

			

			//}

			//else

			//	last_segment_start = start_ori;

		}

#else

*/

		int cur_hsp_non_gap_end;

		if (end_pos < queryStr.length()) //does not end with gap, record last_segment, for use with next HSP segment checking

		{

			cur_hsp_non_gap_end = cur_hsp_end;

		}

		else//end with gap, then look for next_HSP non-gap start position

		{

			//cur_hsp_non_gap_end = cur_hsp_start + start_pos*3 - 1;

			tgt_start_pos = ConvertToNoGappedPos(start_pos, tgt_gap_starts, tgt_gap_ends, false);

			cur_hsp_non_gap_end = cur_hsp_start + tgt_start_pos*3 -1;

		}			

		if (next_hsp_start + next_start_pos*3 > cur_hsp_non_gap_end + MIN_INTRON_LEN //the gap between cur_hsp_end and next HSP's first non-gap position is larger than MIN_INTRON_LEN
			|| (next_hsp_start - cur_hsp_start) % 3 != 0 //or, cur hsp and next hsp not in same frame
			|| ((next_hsp_start - cur_hsp_start) % 3 == 0 && hspID_due_to_stop.find(next_hsp_ID)!= hspID_due_to_stop.end())) //or, the two HSPs are in frame, but there are in-frame stop between them

		{

			//donor_segments.insert(multimap<int, pair<int, int> >::value_type(cur_hsp_non_gap_end+1,

			//	pair<int, int>(start_ori, next_hsp_start+next_start_pos*3-1)));

			donor_segments_pair.push_back( SegmentsInThreeBounds(start_ori, cur_hsp_non_gap_end+1, next_hsp_start+next_start_pos*3-1) );

			//gap_centers.push_back((next_hsp_start+next_start_pos*3+cur_hsp_non_gap_end+1)/2);

			//donor_acceptor_HSP_ID.push_back(hsp_ID);



#ifdef DEBUG

			dataMgr.outFile << "segment3: " << start_ori << "<-" << cur_hsp_non_gap_end+1 << "->" << next_hsp_start+next_start_pos*3-1 << "\n";

#endif



			//donor_acceptor_HSP_ID.push_back( vector<int>(1, next_hsp_ID) );



			last_segment_start = 0;

			last_segment_end = 0;

		}

		else //next HSP follows immediately and starts with non-gap, record the partial segment to extend into next HSP

		{

			//last_segment_start = start_ori;

			last_segment_end = 0;



			//donor_acceptor_HSP_ID.back().push_back(next_hsp_ID);

		}

//#endif

		}



	}

	else //no next HSP

	{

		//must be ending with non-gap!

		last_segment_start = start_ori;



	}



	last_segment_hsp_ID = hsp_ID; //update last_segment_hsp_ID (curHSP must exist here)



	return next_hsp_exists;

}






void ACCP_DONR_graph::GetDonorsAcceptors_nPerBorder(int search_start, int search_end, 
													//pair<int, char*>& chr_seq, 
													vector<string>& chr_seq, 
		vector<int>& accs, vector<int>& dons, vector< string >& acc_next_2nt, vector< string >& don_prev_2nt)

{

#ifdef DEBUG



	dataMgr.outFile << "searching: " << search_start << " to " << search_end << "\n";



#endif



	int len, string_start; 



	string segStr, prev_2nt, next_2nt;



	if (search_start > 0) //positive strand



	{



		if (search_start > search_end) 



		{



			if (start_site >= search_start)



				return;







			//len = search_start - search_end + 1;



			//string_start = search_end-2 < start_site ? start_site : search_end-2;



			string_start = search_end-1 < start_site? start_site : search_end-1; //UPDATE: extend only 1 position at the front end







			len = search_start - string_start + 1;







			if (string_start < 1)



			{



				len -= 1-string_start;



				string_start = 1;



			}







			//GetSubstrFromVecStrs(chr_seq, true, string_start-1, len+4, segStr); //extends 2 positions (for double-char checking) at both ends



			GetSubstrFromVecStrs_ForRepair(chr_seq, true, string_start-1, len, segStr, prev_2nt, next_2nt); 



			StrToLower(segStr);



			



#ifdef DEBUG



			dataMgr.outFile << "seeking dnr 1(backward from "<< string_start << ";len:" << len << "): " << segStr << "\n";



#endif







			//dataMgr.outFile << "seeking dnr: " << segStr << "\n";



			BackwardSearch_nPerBorder(segStr, string_start, true, DONOR_SITE, 2, dons, true, don_prev_2nt, dataMgr.outFile, 

				prev_2nt, next_2nt); //backward search for donors



			//BackwardSearch(segStr, string_start, true, DONOR_SITE, 2, donors, true); //backward search for donors







			//dataMgr.outFile << "seeking acc: " << segStr << "\n";



			ForwardSearch_nPerBorder(segStr, string_start, true, ACCEPTOR_SITE, 1, accs, false, acc_next_2nt, dataMgr.outFile, 

				prev_2nt, next_2nt); //forward search for acceptors



			//ForwardSearch(segStr, string_start, true, ACCEPTOR_SITE, 1, acceptors, false); //forward search for acceptors



		}



		else 



		{



			if (start_site >= search_end)



				return;







			//len = search_end - search_start + 1;



			//string_start = search_start < start_site ? start_site : search_start;



			string_start = search_start-1 < start_site ? start_site : search_start-1;



			len = search_end - string_start + 1;







			if (string_start < 1)



			{



				len -= 1-string_start;



				string_start = 1;



			}







			//GetSubstrFromVecStrs(chr_seq, true, string_start-1, len, segStr);



			GetSubstrFromVecStrs_ForRepair(chr_seq, true, string_start-1, len, segStr, prev_2nt, next_2nt);



			StrToLower(segStr);







#ifdef DEBUG



			dataMgr.outFile << "seeking dnr 2(forward from " << string_start << ";len:"<< len << "): " << segStr << "\n";



#endif



			//dataMgr.outFile << "seeking dnr: " << segStr << "\n";



			ForwardSearch_nPerBorder(segStr, string_start, true, DONOR_SITE, 2, dons, true, don_prev_2nt, dataMgr.outFile, 

				prev_2nt, next_2nt);



			//ForwardSearch(segStr, string_start, true, DONOR_SITE, 2, donors, true);







			//dataMgr.outFile << "seeking acc: " << segStr << "\n";



			BackwardSearch_nPerBorder(segStr, string_start, true, ACCEPTOR_SITE, 1, accs, false, acc_next_2nt, dataMgr.outFile, 

				prev_2nt, next_2nt);



			//BackwardSearch(segStr, string_start, true, ACCEPTOR_SITE, 1, acceptors, false);



		}



	}



	else //negative strand



	{


#ifdef DEBUG
		dataMgr.outFile << "neg strand" << "\n";
#endif

		if (search_start > search_end) 



		{


			if (start_site >= search_start)



				return;







			//len = search_start - search_end + 1;



			//string_start = (-(search_start+2)) > (-start_site) ? (-start_site) : (-(search_start+2)) ;



			//string_start = (-(search_start+1)) > (-start_site) ? (-start_site) : (-(search_start+1)) ;



			string_start = -search_start > -start_site ? -start_site : -search_start;



			//len = -string_start - search_end + 1;



			len = -string_start - search_end + 2; //for negative strand, extend the other end, so length is 1 extra from string_start!


#ifdef DEBUG
			dataMgr.outFile << "string_start:" << string_start << ";len:" << len << "\n";
#endif





			if (string_start < 1)



			{



				len -= 1-string_start;



				string_start = 1;



			}




#ifdef DEBUG
			dataMgr.outFile << "again string_start:" << string_start << ";len:" << len << "\n";
#endif



			//GetSubstrFromVecStrs(chr_seq, false, string_start-1, len+4, segStr); //extends 2 positions (for double-char checking) at both ends



			GetSubstrFromVecStrs_ForRepair(chr_seq, false, string_start-1, len, segStr, prev_2nt, next_2nt);



			StrToLower(segStr);







#ifdef DEBUG



			dataMgr.outFile << "seeking dnr 3(forward from " << string_start << ";len:"<< len << "): " << segStr << "\n";



#endif



			ForwardSearch_nPerBorder(segStr, string_start, false, DONOR_SITE_REV, 2, dons, true, don_prev_2nt, dataMgr.outFile, 

				prev_2nt, next_2nt); //now it's forward search for donors



			//ForwardSearch(segStr, string_start, false, DONOR_SITE_REV, 2, donors, true); //now it's forward search for donors






#ifdef DEBUG
			dataMgr.outFile << "seeking acc: (backward)" << segStr << "\n";
#endif


			BackwardSearch_nPerBorder(segStr, string_start, false, ACCEPTOR_SITE_REV, 1, accs, false, acc_next_2nt, dataMgr.outFile, 

				prev_2nt, next_2nt); //backward search for acceptors



			//BackwardSearch(segStr, string_start, false, ACCEPTOR_SITE_REV, 1, acceptors, false); //backward search for acceptors



		}



		else 



		{



			if (start_site >= search_end)



				return;







			//len = search_end - search_start + 1;



			//string_start = -(search_end) > -start_site? -start_site : -(search_end);



			//string_start = -(search_end+1) > -start_site? -start_site : -(search_end+1);



			string_start = -search_end > -start_site? -start_site : -search_end;



			//len = -string_start - search_start + 1;



			len = -string_start - search_start + 2;







			if (string_start < 1)



			{



				len -= 1-string_start;



				string_start = 1;



			}







			//GetSubstrFromVecStrs(chr_seq, false, string_start-1, len, segStr);



			GetSubstrFromVecStrs_ForRepair(chr_seq, false, string_start-1, len, segStr, prev_2nt, next_2nt);



			StrToLower(segStr);



			//if (string_start == 463236)



			//int stop = 1;







#ifdef DEBUG



			dataMgr.outFile << "seeking dnr 4(backward from " << string_start << ";len:" << len << "): " << segStr << "\n";



#endif



			BackwardSearch_nPerBorder(segStr, string_start, false, DONOR_SITE_REV, 2, dons, true, don_prev_2nt, dataMgr.outFile, 

				prev_2nt, next_2nt);



			//BackwardSearch(segStr, string_start, false, DONOR_SITE_REV, 2, donors, true);







			//dataMgr.outFile << "seeking acc (forward): " << segStr << "\n";



			ForwardSearch_nPerBorder(segStr, string_start, false, ACCEPTOR_SITE_REV, 1, accs, false, acc_next_2nt, dataMgr.outFile, 

				prev_2nt, next_2nt);



			//ForwardSearch(segStr, string_start, false, ACCEPTOR_SITE_REV, 1, acceptors, false);



		}

	}



}








//"align_start" and "align_end" are start/end positions of the query segment in "alignment"

int ACCP_DONR_graph::GetDonorHSP(vector<HSP_Gene_Pair*>& HSPs, int donor_site, int hsp_ID, int segment_index, 

								  Input_Alignment& donor_align, int& align_start, int& align_end, 

								  char& border_query_aa_front, char& border_query_aa_end, 

								  vector<string>& chr_seq) //store results here
								  //pair<int, char*>& chr_seq)

{

	//int splice_start = donor_site - SPLICE_LEN*3 + 1; //the beginning position of exon-piece after splicing

	int splice_start;



	//find THE hsp to the left of donor

	int cur_hsp_start, cur_hsp_end, gene_start, gene_end, hsp_index=-1;

	//int tmp_start;

	int i;



	map<int, int>::iterator mapIt = hspID_to_hspIndex_map.find(hsp_ID);

	if (mapIt == hspID_to_hspIndex_map.end()) //not found, something wrong!

	{

		cout << "invalid donor site (hsp_ID)... why? " << donor_site << "(" << hsp_ID << ")" << "\n"; //because the last/first HSP was erased!

		dataMgr.outFile << "invalid donor site (hsp_ID)... why? " << donor_site << "(" << hsp_ID << ")" << "\n"; 

		exit(-1);

		//return -1;

	}



	i= (*mapIt).second;

	//i = hsp_ID; //MODIFIED: now hsp_ID is directly the index(relative pos) of HSP!

	if (donor_site > 0) //positive strand HSP

	{

				cur_hsp_start = HSPs[i]->HSP_start;

				cur_hsp_end = HSPs[i]->HSP_end; //(*vec_revIt)->HSP_end;

				gene_start = HSPs[i]->gene_start; //(*vec_revIt)->gene_start;

				gene_end = HSPs[i]->gene_end; //(*vec_revIt)->gene_end;

				hsp_index = hsp_ID; //(*vec_revIt)->ID;

	}

	else //negative HSP

	{

			cur_hsp_start = - HSPs[i]->HSP_end;

			cur_hsp_end = - HSPs[i]->HSP_start; //-(*vec_It)->HSP_start; //negate

			gene_start = HSPs[i]->gene_start; //(*vec_It)->gene_start;

			gene_end = HSPs[i]->gene_end; //(*vec_It)->gene_end;

			hsp_index = hsp_ID; //(*vec_It)->ID;

	}



#ifdef DEBUG

	dataMgr.outFile << "found hsp for donor " << donor_site << ": hsp_ID " << hsp_index << "\n";

#endif



		//make sure hsp_index != -1

		if (hsp_index == -1) //didn't find that hsp?

		//if (hsp_ID < 0 || hsp_ID >= HSPs.size())//MODIFIED: now hsp_ID is directly the index(relative pos) of HSP!

		{

			//cout << "invalid donor site... why? " << donor_site << "\n"; //because the last/first HSP was erased!

			//exit(-1);

			dataMgr.outFile << "invalid donor site ... why? " << donor_site << "(" << hsp_ID << ")" << "\n"; 

			exit(-1);

			//return -1;

		}



		if (cur_hsp_start > donor_site) //in case this hsp didn't pan out... back to our old method!

			return GetDonorHSP_S2(HSPs, donor_site, segment_index, 

									donor_align, align_start, align_end, 

									border_query_aa_front, border_query_aa_end, chr_seq);



		map<int, Input_Alignment>::iterator align_it = dataMgr.input_alignments.find(hsp_index);

		string& queryStr = (*align_it).second.query_align;

		string& matchStr = (*align_it).second.match_align;

		string& targetStr = (*align_it).second.target_align;



#ifdef DEBUG

		dataMgr.outFile << (*align_it).second << "\n";

#endif



		//ProteinToDNAStyle(queryStr);

		//ProteinToDNAStyle(matchStr);

		//ProteinToDNAStyle(targetStr);



		splice_start = cur_hsp_start;



		//if (donor_site == -1206537)

		//	int stophere = 1;



		if (!GetHSPSpliceAlignment1(targetStr, queryStr, matchStr, cur_hsp_start, cur_hsp_end, 

				gene_start, gene_end, splice_start, donor_site, donor_align, align_start, align_end, 

				border_query_aa_front, border_query_aa_end, chr_seq))

		{

#ifdef DEBUG

		dataMgr.outFile << "donorHSP has in-frame stop" << "\n";

#endif

				return -1;

		}



		return i;

}



int ACCP_DONR_graph::GetDonorHSP_S2(vector<HSP_Gene_Pair*>& HSPs, int donor_site, int segment_index, 

								  Input_Alignment& donor_align, int& align_start, int& align_end, 

								  char& border_query_aa_front, char& border_query_aa_end, 

								  vector<string>& chr_seq) //store results here
								  //pair<int, char*>& chr_seq)

{

	//int splice_start = donor_site - SPLICE_LEN*3 + 1; //the beginning position of exon-piece after splicing

	int splice_start;



	//find THE hsp to the left of donor

	int cur_hsp_start, cur_hsp_end, gene_start, gene_end;

	int j; //j is the index (position) of HSP in the vector

	int hsp_index = FindDonorHSPIndex_UseLargestHSPInSegment(donor_site, HSPs, segment_index, 

		cur_hsp_start, cur_hsp_end, gene_start, gene_end, j);



#ifdef DEBUG

	dataMgr.outFile << "found hsp for donor " << donor_site << ": hsp_ID " << hsp_index << "\n";

#endif



		//make sure hsp_index != -1

		if (hsp_index == -1) //didn't find that hsp?

		{

			//cout << "invalid donor site... why? " << donor_site << "\n"; //because the last/first HSP was erased!

			//exit(-1);

			dataMgr.outFile << "GetDonorHSP_S2: invalid donor site (still)... why? " << donor_site << "\n"; 

			exit(-1);

			//return -1;

		}



		map<int, Input_Alignment>::iterator align_it = dataMgr.input_alignments.find(hsp_index);

		string& queryStr = (*align_it).second.query_align;

		string& matchStr = (*align_it).second.match_align;

		string& targetStr = (*align_it).second.target_align;



#ifdef DEBUG

		dataMgr.outFile << (*align_it).second << "\n";

#endif



		//ProteinToDNAStyle(queryStr);

		//ProteinToDNAStyle(matchStr);

		//ProteinToDNAStyle(targetStr);



		splice_start = cur_hsp_start;



		//if (donor_site == -10062179)

		//	int stophere = 1;



		if (!GetHSPSpliceAlignment1(targetStr, queryStr, matchStr, cur_hsp_start, cur_hsp_end, 

				gene_start, gene_end, splice_start, donor_site, donor_align, align_start, align_end, 

				border_query_aa_front, border_query_aa_end, chr_seq))

		{

#ifdef DEBUG

		dataMgr.outFile << "donorHSP has in-frame stop" << "\n";

#endif

				return -1;

		}



		return j;

}






int ACCP_DONR_graph::GetAcceptorHSP(vector<HSP_Gene_Pair*>& HSPs, int acceptor_site, int hsp_ID, int segment_index, 

									 Input_Alignment& acceptor_align, int& align_start, int& align_end, 

									 char& border_query_aa_front, char& border_query_aa_end, 

									 vector<string>& chr_seq) //store results here
									 //pair<int, char*>& chr_seq)

{

	//modified: splice_end until the end of splicing segment

	//int splice_end = acceptor_site + SPLICE_LEN*3 - 1;

	int splice_end;



	//find THE hsp to the right of acceptor

	int cur_hsp_start, cur_hsp_end, gene_start, gene_end, hsp_index=-1;

	//int tmp_end;

	int i;



	map<int, int>::iterator mapIt = hspID_to_hspIndex_map.find(hsp_ID);

	if (mapIt == hspID_to_hspIndex_map.end()) //not found, something wrong!

	{

		cout << "invalid acceptor site (hsp_ID)... why? " << acceptor_site << "(" << hsp_ID << ")" << "\n"; //because the last/first HSP was erased!

		dataMgr.outFile << "invalid acceptor site (hsp_ID)... why? " << acceptor_site << "(" << hsp_ID << ")" << "\n"; 

		exit(-1);

		//return -1;

	}



	i= (*mapIt).second;

	if (acceptor_site > 0) //positive strand HSP

	{

				hsp_index = hsp_ID; //(*vec_It)->ID;

				cur_hsp_start = HSPs[i]->HSP_start; //(*vec_It)->HSP_start;

				cur_hsp_end = HSPs[i]->HSP_end;

				gene_start = HSPs[i]->gene_start; //(*vec_It)->gene_start;

				gene_end = HSPs[i]->gene_end; //(*vec_It)->gene_end;

	}

	else //negative HSP

	{

				hsp_index = hsp_ID; //(*vec_revIt)->ID;

				cur_hsp_start = - HSPs[i]->HSP_end; //-(*vec_revIt)->HSP_end; //negate

				cur_hsp_end = - HSPs[i]->HSP_start;

				gene_start = HSPs[i]->gene_start; //(*vec_revIt)->gene_start;

				gene_end = HSPs[i]->gene_end; //(*vec_revIt)->gene_end;		

	}

#ifdef DEBUG

	dataMgr.outFile << "found hsp for acc " << acceptor_site << ": hsp_ID " << hsp_index << "\n";

#endif



		//make sure hsp_index != -1

		if (hsp_index == -1) //didn't find that hsp?

		{

			//cout << "invalid acceptor site... why? " << acceptor_site << "\n";

			//exit(-1);

			dataMgr.outFile << "invalid acceptor site ... why? " << acceptor_site << "(" << hsp_ID << ")" << "\n"; 

			exit(-1);

			//return -1;

		}



		if (cur_hsp_end < acceptor_site) //in case this hsp does not work, use the old method to get a good hsp

			return GetAcceptorHSP_S2(HSPs, acceptor_site, segment_index, acceptor_align, align_start, align_end, 

							border_query_aa_front, border_query_aa_end, chr_seq);



		map<int, Input_Alignment>::iterator align_it = dataMgr.input_alignments.find(hsp_index);

		string& queryStr = (*align_it).second.query_align;

		string& matchStr = (*align_it).second.match_align;

		string& targetStr = (*align_it).second.target_align;



		//ProteinToDNAStyle(queryStr);

		//ProteinToDNAStyle(matchStr);

		//ProteinToDNAStyle(targetStr);



		splice_end = cur_hsp_end; //splice_end until the end of HSP?



		if (!GetHSPSpliceAlignment1(targetStr, queryStr, matchStr, cur_hsp_start, cur_hsp_end, 

				gene_start, gene_end, acceptor_site, splice_end, acceptor_align, align_start, align_end, 

				border_query_aa_front, border_query_aa_end, chr_seq))

		{

#ifdef DEBUG

		dataMgr.outFile << "acceptorHSP has in-frame stop" << "\n";

#endif

				return -1;

		}



		return i;

}



int ACCP_DONR_graph::GetAcceptorHSP_S2(vector<HSP_Gene_Pair*>& HSPs, int acceptor_site, int segment_index, 

									   Input_Alignment& acceptor_align, int& align_start, int& align_end, 

									   char& border_query_aa_front, char& border_query_aa_end, 
									   vector<string>& chr_seq) //store results here
									   //pair<int, char*>& chr_seq)

{

	//modified: splice_end until the end of splicing segment

	//int splice_end = acceptor_site + SPLICE_LEN*3 - 1;

	int splice_end;



	//find THE hsp to the right of acceptor

	int cur_hsp_start, cur_hsp_end, gene_start, gene_end;

	int j; //relative position (index) of HSP in the vector

	int hsp_index = FindAcceptorHSPIndex_UseLargestHSPInSegment(acceptor_site, HSPs, segment_index, cur_hsp_start, cur_hsp_end, gene_start, gene_end, j);



#ifdef DEBUG

	dataMgr.outFile << "found hsp for acc " << acceptor_site << ": hsp_ID " << hsp_index << "\n";

#endif



		//make sure hsp_index != -1

		if (hsp_index == -1) //didn't find that hsp?

		{

			//cout << "invalid acceptor site... why? " << acceptor_site << "\n";

			//exit(-1);

			dataMgr.outFile << "GetAcceptorHSP_S2: invalid acceptor site (still)... why? " << acceptor_site << "\n"; 

			exit(-1);

			//return -1;

		}



		map<int, Input_Alignment>::iterator align_it = dataMgr.input_alignments.find(hsp_index);

		string& queryStr = (*align_it).second.query_align;

		string& matchStr = (*align_it).second.match_align;

		string& targetStr = (*align_it).second.target_align;



		//ProteinToDNAStyle(queryStr);

		//ProteinToDNAStyle(matchStr);

		//ProteinToDNAStyle(targetStr);



		splice_end = cur_hsp_end; //splice_end until the end of HSP?



		if (!GetHSPSpliceAlignment1(targetStr, queryStr, matchStr, cur_hsp_start, cur_hsp_end, 

				gene_start, gene_end, acceptor_site, splice_end, acceptor_align, align_start, align_end, 

				border_query_aa_front, border_query_aa_end, chr_seq))

		{

#ifdef DEBUG

		dataMgr.outFile << "acceptorHSP has in-frame stop" << "\n";

#endif

				return -1;

		}



		return j;

}






void ACCP_DONR_graph::EraseSpecificSites(vector<int>& sites, vector<int>& site_regions, 

								 vector<int>& site_head_tail, vector<int>& site_HSP_ID, vector<string>& site_2nt, 

								 vector<int>& sites_to_be_erased)

{

	vector<int>::iterator setIt, regionIt, htIt, idIt;

	vector<string>::iterator site_2nt_It;

	setIt = sites.begin();

	regionIt = site_regions.begin();

	htIt = site_head_tail.begin();

	idIt = site_HSP_ID.begin();

	site_2nt_It = site_2nt.begin();



	int i;

	for (i=sites_to_be_erased.size()-1; i>=0; i--)

	{

		sites.erase(setIt+sites_to_be_erased[i]);

		site_regions.erase(regionIt+sites_to_be_erased[i]);

		site_head_tail.erase(htIt+sites_to_be_erased[i]);

		site_HSP_ID.erase(idIt+sites_to_be_erased[i]);

		site_2nt.erase(site_2nt_It+sites_to_be_erased[i]);

	}



}



void ACCP_DONR_graph::EraseSites(vector<int>& sites, vector<int>& site_regions, 

								 vector<int>& site_head_tail, vector<int>& site_HSP_ID, vector<string>& site_2nt, 

								 int border, bool erase_before_border)

{
#ifdef DEBUG
	dataMgr.outFile << "border: " << border << "(erase_before_border:" << erase_before_border << ") border itself is kept\n";
#endif

	if (erase_before_border)

	{

		vector<int>::iterator setIt, regionIt, htIt, idIt;

		vector<string>::iterator site_2nt_It;

		setIt = sites.begin();

		regionIt = site_regions.begin();

		htIt = site_head_tail.begin();

		idIt = site_HSP_ID.begin();

		site_2nt_It = site_2nt.begin();

		while (setIt != sites.end() && *setIt < border )

		{

			setIt++;

			regionIt++;

			htIt++;

			idIt++;

			site_2nt_It++;

		}

		if (setIt != sites.begin())

		{

			sites.erase(sites.begin(), setIt);

			site_regions.erase(site_regions.begin(), regionIt);

			site_head_tail.erase(site_head_tail.begin(), htIt);

			site_HSP_ID.erase(site_HSP_ID.begin(), idIt);

			site_2nt.erase(site_2nt.begin(), site_2nt_It);

		}

	}

	else

	{

		vector<int>::reverse_iterator setIt, regionIt, htIt, idIt;

		vector<int>::iterator eSetIt, eRegionIt, eHtIt, eIdIt;

		vector<string>::reverse_iterator site_2nt_It;

		vector<string>::iterator eSite_2nt_It;

		setIt = sites.rbegin();

		eSetIt = sites.end();

		regionIt = site_regions.rbegin();

		eRegionIt = site_regions.end();

		htIt = site_head_tail.rbegin();

		eHtIt = site_head_tail.end();

		idIt = site_HSP_ID.rbegin();

		eIdIt = site_HSP_ID.end();

		site_2nt_It = site_2nt.rbegin();

		eSite_2nt_It = site_2nt.end();

		while (setIt != sites.rend() && *setIt > border)

		{

			setIt++;

			regionIt++;

			htIt++;

			idIt++;

			site_2nt_It++;



			eSetIt--;

			eRegionIt--;

			eHtIt--;

			eIdIt--;

			eSite_2nt_It--;

		}

		if (setIt != sites.rbegin())

		{

			sites.erase(eSetIt, sites.end());

			site_regions.erase(eRegionIt, site_regions.end());

			site_head_tail.erase(eHtIt, site_head_tail.end());

			site_HSP_ID.erase(eIdIt, site_HSP_ID.end());

			site_2nt.erase(eSite_2nt_It, site_2nt.end());

		}

	}



}




//region_erase_start is the start region number that should be erased

//region_erase_end is the region number that is right after the erased region (itself should not be erased)

void ACCP_DONR_graph::EraseSites_ByRegions(vector<int>& sites, vector<int>& site_regions, 

								 vector<int>& site_head_tail, vector<int>& site_HSP_ID, vector<string>& site_2nt, 

								 int region_erase_start, int region_erase_end) 

{
#ifdef DEBUG
	dataMgr.outFile << "region_erase_start: " << region_erase_start << "; region_erase_end: " << region_erase_end << "\n";
#endif

	int i = 0;

	int index_start = -1, index_end = -1;

	

	while (i<site_regions.size())

	{

		if (site_regions[i] >= region_erase_start)
		{

			if (site_regions[i] < region_erase_end)
			{
				if (index_start == -1) //haven't changed (first time), then assign the value to it
					index_start = i;
			}
			else
			{
				index_end = i;
				break;
			}

		}

		i++;

	}

	if (index_start != -1)

	{

		if (index_end != -1)

		{

			sites.erase(sites.begin()+index_start, sites.begin()+index_end);

			site_regions.erase(site_regions.begin()+index_start, site_regions.begin()+index_end);

			site_head_tail.erase(site_head_tail.begin()+index_start, site_head_tail.begin()+index_end);

			site_HSP_ID.erase(site_HSP_ID.begin()+index_start, site_HSP_ID.begin()+index_end);

			site_2nt.erase(site_2nt.begin()+index_start, site_2nt.begin()+index_end);

		}

		else

		{

			sites.erase(sites.begin()+index_start, sites.end());

			site_regions.erase(site_regions.begin()+index_start, site_regions.end());

			site_head_tail.erase(site_head_tail.begin()+index_start, site_head_tail.end());

			site_HSP_ID.erase(site_HSP_ID.begin()+index_start, site_HSP_ID.end());

			site_2nt.erase(site_2nt.begin()+index_start, site_2nt.end());

		}

	}

	else

	{

		dataMgr.outFile << "nothing to erase(index_start is -1): region_erase_start:" << region_erase_start << "; region_erase_end:" << region_erase_end 

			<< "\n";

		//cout << "nothing to erase(index_start is -1)? region_erase_start:" << region_erase_start << "; region_erase_end:" << region_erase_end << "\n";

		//exit(-1);

	}

}



bool ACCP_DONR_graph::ComputeSpliceAlignment_SplitOL1(Input_Alignment& donor_align, Input_Alignment& acceptor_align,
	int donor_start, int donor_end, int acceptor_start, int acceptor_end, 
	string& splice_align, 
	//the following FOUR are no longer used
	int& donor_splice_end, int& acceptor_splice_start, 
	int& new_donor_end, int& new_acceptor_start, //bool alignment_update, //this function doesn't update alignment!
	//char missing_query_aa, 
	string donor_2nt, string acceptor_2nt, int donor_tail_num, int acceptor_head_num)
//donor_splice_end: index on the donor_align.matchStr of the position where the alignment is cut (index is the last position of useful string)
//acceptor_splice_start: index on acceptor_align.matchStr for the cut (index is the first position of useful string)
{

	string codonStr;

	if (donor_tail_num > 0) // acceptor_head_num must also > 0

	{

		codonStr = donor_2nt.substr(2-donor_tail_num) + acceptor_2nt.substr(0, acceptor_head_num);

		map<string, char>::iterator codon_tbl_it = DNA_CODON_TBL_SL.find(codonStr);

		if (codon_tbl_it != DNA_CODON_TBL_SL.end())
			if ((*codon_tbl_it).second == '*') //is stop codon
				return false;
	}



	new_acceptor_start = acceptor_start; //the possibly new acceptor_start position (nextHSP's gene_start)

	new_donor_end = donor_end; //possible new donor_end position (curHSP's gene_end)



	splice_align = "";



	if (donor_start == 0 || acceptor_start == 0) //if either segment is all gaps

	{

		splice_align += donor_align.match_align;

		splice_align += acceptor_align.match_align;

//		cout << "already done!" << "\n";

		donor_splice_end = donor_align.match_align.length()-1; //nothing is cut

		acceptor_splice_start = 0; //nothing is cut



#ifdef DEBUG

		dataMgr.outFile << "donor or acceptor segment is all gap (donor_start:" << donor_start 

			<< "; acceptor_start:" << acceptor_start << ")" << "\n";

#endif

		return true;

	}



	//int i;

	int search_start1, search_end1, search_start2, search_end2;

	if (donor_end < acceptor_start) //non overlap

	{

		//revised: now check the trailing / heading gaps in query_align before inserting more gaps!

		//revised again: if there is no overlap and the missing query part is exactly three base pairs,

		//then check donor_prev_nt with acceptor_next_nt against donor_aa_end/acceptor_aa_front, if it's a match, then revise splice_align

		//revised again: now do global alignment for all cases when there's gap on query (globally align query gap with target gap)

		//MODIFIED AGAIN: alignment now starts from the last non-matching position of donor-align to the first non-match of acceptor-align



		int total_gap_len = acceptor_start - donor_end - 1; //this is the total gap length that we should have



		int last_non_gap_pos_donor_align = donor_align.query_align.find_last_not_of("-");

		int first_non_gap_pos_acceptor_align = acceptor_align.query_align.find_first_not_of("-");

		if (last_non_gap_pos_donor_align == string::npos)

			last_non_gap_pos_donor_align = -1;

		if (first_non_gap_pos_acceptor_align == string::npos)

			first_non_gap_pos_acceptor_align = acceptor_align.query_align.length();



		int need_gap_len = total_gap_len - (donor_align.query_align.length() - 1 - last_non_gap_pos_donor_align + first_non_gap_pos_acceptor_align);

		

		int last_match_pos_donor_align = donor_align.match_align.find_last_not_of("+ ");

		int first_match_pos_acceptor_align = acceptor_align.match_align.find_first_not_of("+ ");

		if (last_match_pos_donor_align == string::npos)

			last_match_pos_donor_align = -1;

		if (first_match_pos_acceptor_align == string::npos)

			first_match_pos_acceptor_align = acceptor_align.query_align.length();



		//if (first_non_gap_pos_acceptor_align > 0 || last_non_gap_pos_donor_align < donor_align.query_align.length()-1)

		if (first_match_pos_acceptor_align > 0 || last_match_pos_donor_align < donor_align.query_align.length()-1)

		{

			//MODIFIED: do global alignment for the stuff in the gap region

			string target_seq_str = "";

			//target_seq_str must be a multiple of 3 (since donor/acceptor must have matching frame)

			int donor_tail_length = donor_align.target_align.length()-last_match_pos_donor_align-1-donor_tail_num;

			if (donor_tail_length > 0)

				target_seq_str = donor_align.target_align.substr(last_match_pos_donor_align+1, donor_tail_length);

			if (donor_tail_num > 0) // acceptor_head_num must also > 0

			{

				//string codonStr = donor_2nt.substr(2-donor_tail_num) + acceptor_2nt.substr(0, acceptor_head_num);

				target_seq_str += codonStr;

			}

			int acceptor_head_length = first_match_pos_acceptor_align - acceptor_head_num;

			if (acceptor_head_length > 0)

				target_seq_str += acceptor_align.target_align.substr(acceptor_head_num, acceptor_head_length);

			

			//special handling to lowercase letters (actual nucleotides)

			if (!TranslateTargetDNA(target_seq_str))

				return false; 



			//modified: ignore the very very long sequence, for now

			if (target_seq_str.length() > 4096)

			{

#ifdef DEBUG

				dataMgr.outFile << "length>4096, treat as gaps\n";

#endif

				goto TEMP_SOLUTION1;

			}



			//string query_seq_str = query_seq.substr(donor_end/3, (acceptor_start-donor_end-1)/3);

			string query_seg_missing_str = query_seq.substr(donor_end/3, (acceptor_start-donor_end-1)/3);

			ProteinToDNAStyle(query_seg_missing_str);



			string query_seq_str = donor_align.query_align.substr(last_match_pos_donor_align+1);

			EraseAll(query_seq_str, '-');

			query_seq_str += query_seg_missing_str;

			string query_acc_str = acceptor_align.query_align.substr(0, first_match_pos_acceptor_align);

			EraseAll(query_acc_str, '-');

			query_seq_str += query_acc_str;

			

			Input_Alignment gapAlign;

			float align_pid;

			GetGlobalAlignment(query_seq_str, target_seq_str, gapAlign, align_pid, dataMgr.outFile);

			splice_align += donor_align.match_align.substr(0, last_match_pos_donor_align+1);

			splice_align += gapAlign.match_align;

			splice_align += acceptor_align.match_align.substr(first_match_pos_acceptor_align);

		}

		else

		{

TEMP_SOLUTION1:		splice_align += donor_align.match_align;



/*		bool broken_codon_is_match = false;

#ifdef DEBUG

		dataMgr.outFile << "missing query char: " << missing_query_aa << "; donor_tail_num: " << donor_tail_num << "; acceptor_head_num: " << acceptor_head_num << "\n";

		dataMgr.outFile << "donor_tail string: " << donor_2nt << "; acceptor_2nt: " << acceptor_2nt << "\n";

#endif

		if (total_gap_len == 3 && donor_tail_num > 0 && acceptor_head_num > 0)

		{

			string codonStr = donor_2nt.substr(2-donor_tail_num) + acceptor_2nt.substr(0, acceptor_head_num);

#ifdef DEBUG

			dataMgr.outFile << "codonStr: " << codonStr << "\n";

#endif			

			map<string, char>::iterator codonIt = DNA_CODON_TBL_SL.find(codonStr);

			if (codonIt != DNA_CODON_TBL_SL.end() )

			{

#ifdef DEBUG				

				dataMgr.outFile << "amino acid: " << (*codonIt).second << "\n";

#endif

				if (((*codonIt).second) == missing_query_aa)

				{

					broken_codon_is_match = true;

					splice_align.replace(splice_align.size() - donor_tail_num, donor_tail_num, 

						string(3, missing_query_aa)); //replace last donor_tail_num of chars to 3 missing_query_aa!

				}

			}

		}

*/

		//for (i=0; i<acceptor_start-donor_end-1; i++)

		//	splice_align += ' ';

		if (need_gap_len > 0)

		{

			string tmpStr(need_gap_len, ' ');//(acceptor_start-donor_end-1, ' ');

			splice_align += tmpStr;

		}

//		if (broken_codon_is_match)

//			splice_align += acceptor_align.match_align.substr(acceptor_head_num);

//		else

			splice_align += acceptor_align.match_align;

		}



		donor_splice_end = donor_align.match_align.length()-1;

		acceptor_splice_start = 0;



#ifdef DEBUG

		dataMgr.outFile << "donor and acceptor segments no overlap, should have " << total_gap_len << " gaps; filled in " << need_gap_len << " gaps" << "\n";

			//tmpStr.length() << " gaps" << "\n";

#endif

	}

	else //overlap? ( acceptor_start ... donor_end)

	{

		int identity1, identity2;

		if (donor_start > acceptor_end) //cross align! so only one aligns, the other is turned to gaps 

			//(acceptor_start...acceptor_end...donor_start ... donor_end)

		{

			const string& qStr1 = donor_align.match_align;

			const string& qStr2 = acceptor_align.match_align;



			identity1 = qStr1.length() - CountCharInStr(qStr1, ' ') - CountCharInStr(qStr1, '+');

			identity2 = qStr2.length() - CountCharInStr(qStr2, ' ') - CountCharInStr(qStr2, '+');



			if (identity1 > identity2)

			{

				splice_align += qStr1;

				//for (i=0; i<qStr2.length(); i++)

				//	splice_align += ' ';

				string tmpStr(qStr2.length(), ' ');

				splice_align += tmpStr;



				donor_splice_end = qStr1.length()-1;

				acceptor_splice_start = qStr2.length(); //acceptor is all cut (to the end)



				dataMgr.outFile << "acceptor is turned to gap" << "\n";

			}

			else

			{

				//for (i=0; i<qStr1.length(); i++)

				//	splice_align += ' ';

				string tmpStr(qStr1.length(), ' ');

				splice_align += tmpStr;

				splice_align += qStr2;



				donor_splice_end = -1; //donor is all cut

				acceptor_splice_start = 0;



				dataMgr.outFile << "donor is turned to gap" << "\n";

			}



#ifdef DEBUG

			dataMgr.outFile << " (cross align, one is turned to gaps)" << "\n";

#endif

		}

		else //(...donor_start...acceptor_start...acceptor_end ... donor_end)

			//or (...donor_start...acceptor_start ... donor_end...acceptor_end)

			//or (...acceptor_start...donor_start...acceptor_end ... donor_end)

			//or (...acceptor_start...donor_start ... donor_end...acceptor_end)

		{

			if (donor_start < acceptor_start)//(...donor_start...acceptor_start...acceptor_end ... donor_end)

			//or (...donor_start...acceptor_start ... donor_end...acceptor_end)

			{

				

				FindRealPos_DNA(donor_align.query_align, donor_start, donor_end, acceptor_start, donor_end, 

					search_start1, search_end1);

				const string& qStr1 = donor_align.match_align.substr(search_start1, search_end1-search_start1+1);

				identity1 = qStr1.length() - CountCharInStr(qStr1, ' ') - CountCharInStr(qStr1, '+');

				int qStr1_left_len = donor_align.query_align.length() - search_end1 - 1; //length of stuff after search_end1 (gaps)



				if (acceptor_end >= donor_end) //donor_start...acceptor_start...donor_end...acceptor_end

				{

					//the overlapping part (query segment) is [acceptor_start, donor_end]

					FindRealPos_DNA(acceptor_align.query_align, acceptor_start, acceptor_end, acceptor_start, donor_end, 

						search_start2, search_end2);

					const string& qStr2 = acceptor_align.match_align.substr(search_start2, search_end2-search_start2+1);

					int qStr2_left_len = search_start2; //length of stuff before search_start2 (gaps)



					//this is only used here

					const string& donor_ol_queryStr = donor_align.query_align.substr(search_start1, search_end1-search_start1+1);

					const string& acceptor_ol_queryStr = acceptor_align.query_align.substr(search_start2, search_end2-search_start2+1);



					vector<int> donor_ol_queryStr_gap_starts, donor_ol_queryStr_gap_ends;

					vector<int> acceptor_ol_queryStr_gap_starts, acceptor_ol_queryStr_gap_ends;

					GetGaps(donor_ol_queryStr, donor_ol_queryStr_gap_starts, donor_ol_queryStr_gap_ends);

					GetGaps(acceptor_ol_queryStr, acceptor_ol_queryStr_gap_starts, acceptor_ol_queryStr_gap_ends);



#ifdef DEBUG

					dataMgr.outFile << "query:" << donor_ol_queryStr << "\n" << "match:" << qStr1 << "\n";

					dataMgr.outFile << "query:" << acceptor_ol_queryStr << "\n" << "match:" << qStr2 << "\n";

					vector<int>::iterator tmp_it1, tmp_it2=donor_ol_queryStr_gap_ends.begin();

					dataMgr.outFile << "donor gaps: ";

					for (tmp_it1 = donor_ol_queryStr_gap_starts.begin(); tmp_it1 != donor_ol_queryStr_gap_starts.end(); tmp_it1++, tmp_it2++)

						dataMgr.outFile << *tmp_it1 << "-" << *tmp_it2 << "; ";

					dataMgr.outFile << "\nacceptor gaps: ";

					tmp_it2 = acceptor_ol_queryStr_gap_ends.begin();

					for (tmp_it1 = acceptor_ol_queryStr_gap_starts.begin(); tmp_it1 != acceptor_ol_queryStr_gap_starts.end(); tmp_it1++, tmp_it2++)

						dataMgr.outFile << *tmp_it1 << "-" << *tmp_it2 << "; ";

					dataMgr.outFile << "\n";

#endif



					vector<int> donor_ol_matchStr_match_starts, donor_ol_matchStr_match_ends;

					vector<int> acceptor_ol_matchStr_match_starts, acceptor_ol_matchStr_match_ends;

					vector<int> donor_match_end_ids, acceptor_match_start_ids;



					int total_ol_id1 = GetMatches(qStr1, donor_ol_matchStr_match_starts, donor_ol_matchStr_match_ends, donor_match_end_ids);

					int total_ol_id2 = GetMatches(qStr2, acceptor_ol_matchStr_match_starts, acceptor_ol_matchStr_match_ends, acceptor_match_start_ids);



					//now try each match_end on donor_ol_matchStr (qStr1) to be used as split point

					vector<int>::iterator d_start_it=donor_ol_matchStr_match_starts.begin();

					vector<int>::iterator d_end_it = donor_ol_matchStr_match_ends.begin();

					int cur_match_count = 0; //donor side match count

					//int max_match_count=-1, max_match_donor_pos=-1, max_match_acceptor_pos=-1; //initialize

					int max_match_count=total_ol_id2, max_match_donor_pos=-1, max_match_acceptor_pos=-1; //initialize

					//first region_end = -1, region_start = 0 (to simulate all overlapping part uses acceptor)

					for (; d_start_it != donor_ol_matchStr_match_starts.end(); d_start_it++, d_end_it++)

					{

						int region_start = (*d_start_it);

						int region_end = (*d_end_it);

						cur_match_count += region_end - region_start + 1;



						int donor_queryPos = ConvertToNoGappedPos(region_end, donor_ol_queryStr_gap_starts, donor_ol_queryStr_gap_ends, false);

						int acceptor_matchPos = ConvertToGappedPos(donor_queryPos, acceptor_ol_queryStr_gap_starts, acceptor_ol_queryStr_gap_ends);



						int acceptor_match_count = ComputeBackId(acceptor_matchPos, acceptor_ol_matchStr_match_starts, 

							acceptor_ol_matchStr_match_ends, acceptor_match_start_ids, total_ol_id2);



						int total_match_count = cur_match_count + acceptor_match_count;

#ifdef DEBUG

						dataMgr.outFile << "donor matchPos (abs): " << region_end << "; converted to queryPos (all match): " << donor_queryPos;

						dataMgr.outFile << "; acceptor matchPos (abs): " << acceptor_matchPos << "\n";

						dataMgr.outFile << "cur max: " << max_match_count << "; this time: " << total_match_count << "\n";

#endif

						

						if (total_match_count > max_match_count)

						{

							max_match_count = total_match_count;

							max_match_donor_pos = region_end;

							max_match_acceptor_pos = acceptor_matchPos;

#ifdef DEBUG

							dataMgr.outFile << "cur donor pos: " << max_match_donor_pos  << "; cur acc pos: " << max_match_acceptor_pos << "\n";

#endif

						}

					}



					//now also compute from the other side (acceptor side), if both come to the same max_match_count but at different positions, keep matched stuff at both sides

					cur_match_count = 0; //acceptor side match count

					int max_match_count_acc=total_ol_id1, max_match_donor_pos_acc=donor_ol_queryStr.length(), max_match_acceptor_pos_acc=acceptor_ol_queryStr.length(); //initialize

					//first region_end = -1, region_start = 0 (to simulate all overlapping part uses acceptor)

					vector<int>::reverse_iterator a_start_it = acceptor_ol_matchStr_match_starts.rbegin(); 

					vector<int>::reverse_iterator a_end_it = acceptor_ol_matchStr_match_ends.rbegin();

					for (; a_start_it != acceptor_ol_matchStr_match_starts.rend(); a_start_it++, a_end_it++)

					{

						int region_start = (*a_start_it);

						int region_end = (*a_end_it);

						cur_match_count += region_end - region_start + 1;



						int acceptor_queryPos = ConvertToNoGappedPos(region_start, acceptor_ol_queryStr_gap_starts, acceptor_ol_queryStr_gap_ends, false);

						int donor_matchPos = ConvertToGappedPos(acceptor_queryPos, donor_ol_queryStr_gap_starts, donor_ol_queryStr_gap_ends);



						int donor_match_count = ComputeFrontId(donor_matchPos, donor_ol_matchStr_match_starts, 

							donor_ol_matchStr_match_ends, donor_match_end_ids, total_ol_id1);



						int total_match_count = cur_match_count + donor_match_count;

#ifdef DEBUG

						dataMgr.outFile << "acceptor matchPos (abs): " << region_start << "; converted to queryPos (all match): " << acceptor_queryPos;

						dataMgr.outFile << "; donor matchPos (abs): " << donor_matchPos << "\n";

						dataMgr.outFile << "cur max: " << max_match_count_acc << "; this time: " << total_match_count << "\n";

#endif

						

						if (total_match_count > max_match_count_acc)

						{

							max_match_count_acc = total_match_count;

							max_match_acceptor_pos_acc = region_start;

							max_match_donor_pos_acc = donor_matchPos;

#ifdef DEBUG

							dataMgr.outFile << "cur donor pos: " << max_match_donor_pos_acc  << "; cur acc pos: " << max_match_acceptor_pos_acc << "\n";

#endif

						}

					}



					int donor_gap_len = qStr1.length() - max_match_donor_pos - 1 + qStr1_left_len;

					int acceptor_gap_len = qStr2_left_len + max_match_acceptor_pos+1;

					int acceptor_gap_len_init = acceptor_gap_len; //hold this position, used later for splice_align

					//splice_align += donor_align.match_align.substr(0, search_start1+max_match_donor_pos+1);

#ifdef DEBUG

					dataMgr.outFile << "First: donor_gap_len: " << donor_gap_len << "; acceptor_gap_len: " << acceptor_gap_len << "\n";

#endif

					//string tmpStr(donor_gap_len + acceptor_gap_len, ' ');

					string tmpStr1;

					int strIt;

					if (donor_gap_len == 0)

						splice_align += donor_align.match_align.substr(0, search_start1+max_match_donor_pos+1);

					else if (donor_gap_len > 0)

					{

						tmpStr1 = donor_align.target_align.substr(search_start1+max_match_donor_pos+1);

						string tmpStr2 = donor_align.target_align.substr(0, search_start1+max_match_donor_pos+1);

						int don_valid_pos = tmpStr2.find_last_not_of('-');

						string tmpStr3;

						if (don_valid_pos == string::npos)

						{

							tmpStr3 = donor_align.query_align.substr(0, search_start1+max_match_donor_pos+1);

							don_valid_pos = tmpStr3.length();

						}

						else

						{

							if (don_valid_pos > 0)

								tmpStr3 = donor_align.query_align.substr(don_valid_pos+1, search_start1+max_match_donor_pos - don_valid_pos);

							splice_align += donor_align.match_align.substr(0, don_valid_pos+1);

						}



						if (donor_gap_len != tmpStr1.length()) //is this thing right? ... let's check

						{

							cout << "donor_gap_len is not tmpStr1.length()" << "\n";

							dataMgr.outFile << donor_align << "\n";

							dataMgr.outFile << "tmpStr1 from " << search_start1+max_match_donor_pos+1 << "\n";

							dataMgr.outFile << tmpStr1 << "\n";

							dataMgr.outFile << "donor_gap_len:" << donor_gap_len << " vs tmpStr.length:" << tmpStr1.length() << "\n";

							dataMgr.outFile << "qStr1:" << qStr1 << "\n";

							dataMgr.outFile << "search_start1:" << search_start1 << "; search_end1:" << search_end1 

								<< "; max_match_donor_pos:" << max_match_donor_pos << "\n";

							dataMgr.outFile << "qStr1_left_len:" << qStr1_left_len << "\n";

							

							exit(-1);

						}

						while ((strIt = tmpStr1.find('-')) != string::npos)

							tmpStr1.erase(strIt, 1);

						//donor_gap_len = tmpStr1.length();

						if (tmpStr3.length()>0 && tmpStr1.length()>0)//(don_valid_pos > 0 && 

						{

							Input_Alignment newAlign;

							float align_pid;

							GetGlobalAlignment(tmpStr1, tmpStr3, newAlign, align_pid, dataMgr.outFile);

							splice_align += newAlign.match_align;

						}

						else //tmpStr3 is empty

						{

							donor_gap_len = tmpStr1.length();

							string tmpStr(donor_gap_len, ' ');

							splice_align += tmpStr;

						}

					}



					if (acceptor_gap_len == 0)

						splice_align += acceptor_align.match_align;

					else if (acceptor_gap_len > 0)

					{

						tmpStr1 = acceptor_align.target_align.substr(0, acceptor_gap_len);

						string tmpStr2 = acceptor_align.target_align.substr(acceptor_gap_len);

						int acc_valid_pos = tmpStr2.find_first_not_of('-');

						string tmpStr3;

						if (acc_valid_pos == string::npos) //everything after acceptor_gap_len is "-"!

						{

							tmpStr3 = acceptor_align.query_align.substr(acceptor_gap_len);

							acc_valid_pos = tmpStr3.length();

						}

						else

						{

							if (acc_valid_pos > 0)

								tmpStr3 = acceptor_align.query_align.substr(acceptor_gap_len, acc_valid_pos);

						}



						if (acceptor_gap_len != tmpStr1.length()) //is this thing right? ... let's check

						{

							cout << "acceptor_gap_len is not tmpStr1.length()" << "\n";

							exit(-1);

						}

						while ((strIt = tmpStr1.find('-')) != string::npos)

							tmpStr1.erase(strIt);



						//acceptor_gap_len = tmpStr1.length();



						if (acc_valid_pos > 0)//align tmpStr1 with tmpStr3, to see if the gaps can be further reduced

						{

							Input_Alignment newAlign;

							float align_pid;

							GetGlobalAlignment(tmpStr1, tmpStr3, newAlign, align_pid, dataMgr.outFile);

							splice_align += newAlign.match_align;

						}

						else //acc_valid_pos must be 0

						{

							acceptor_gap_len = tmpStr1.length();

							string tmpStr(acceptor_gap_len, ' ');

							splice_align += tmpStr;

						}



						if (acceptor_gap_len_init+acc_valid_pos < acceptor_align.match_align.length())

							splice_align += acceptor_align.match_align.substr(acceptor_gap_len_init+acc_valid_pos);

					}

#ifdef DEBUG

					dataMgr.outFile << "Now: donor_gap_len: " << donor_gap_len << "; acceptor_gap_len: " << acceptor_gap_len << "\n";

#endif

				//	string tmpStr(donor_gap_len + acceptor_gap_len, ' ');

				//	splice_align += tmpStr;

				//	splice_align += acceptor_align.match_align.substr(acceptor_gap_len_init);



					//modified when max_match_count == max_match_count_acc (they must be equal?!), 

					//if the two splicing positions for each max_match_count are different:

					//donor_splice_end is at the later of the two positions (one for each max_match_count)

					//similarly for acceptor_splice_start



					int donor_pos_larger, acceptor_pos_smaller;

					max_match_donor_pos_acc--; //convert to the way of splicing using donor side computing?

					max_match_acceptor_pos_acc--;

					if (max_match_donor_pos == max_match_donor_pos_acc) //both methods splice at the same position

					{

					donor_pos_larger = max_match_donor_pos;

					acceptor_pos_smaller = max_match_acceptor_pos;

					//donor_splice_end = max_match_donor_pos;

					//acceptor_splice_start = max_match_acceptor_pos+1;

					donor_splice_end = search_start1 + max_match_donor_pos;

					acceptor_splice_start = acceptor_gap_len_init; //acceptor_gap_len;

					}

					else

					{

						donor_pos_larger = max_match_donor_pos < max_match_donor_pos_acc ? max_match_donor_pos_acc : max_match_donor_pos;

						acceptor_pos_smaller = max_match_acceptor_pos < max_match_acceptor_pos_acc ? max_match_acceptor_pos : max_match_acceptor_pos_acc;

						donor_splice_end = search_start1 + donor_pos_larger;

						acceptor_splice_start = qStr2_left_len + acceptor_pos_smaller + 1;

					}



					

/*					//also change the alignments themselves!

					if (alignment_update)

					{

					//const_cast<string&>(donor_align.query_align).replace(search_start1+max_match_donor_pos+1, 

					//	donor_gap_len, string(donor_gap_len, '-'));

					//const_cast<string&>(donor_align.match_align).replace(search_start1+max_match_donor_pos+1, 

					//	donor_gap_len, string(donor_gap_len, ' '));

					//const_cast<string&>(acceptor_align.query_align).replace(0, acceptor_gap_len, string(acceptor_gap_len, '-'));

					//const_cast<string&>(acceptor_align.match_align).replace(0, acceptor_gap_len, string(acceptor_gap_len, ' '));

						const_cast<string&>(donor_align.query_align).erase(donor_splice_end+1);//(search_start1+max_match_donor_pos+1);

						const_cast<string&>(donor_align.match_align).erase(donor_splice_end+1);//(search_start1+max_match_donor_pos+1);

						const_cast<string&>(donor_align.target_align).erase(donor_splice_end+1);//(search_start1+max_match_donor_pos+1);

						const_cast<string&>(acceptor_align.query_align).erase(0, acceptor_splice_start);//(0, acceptor_gap_len);

						const_cast<string&>(acceptor_align.match_align).erase(0, acceptor_splice_start);//(0, acceptor_gap_len);

						const_cast<string&>(acceptor_align.target_align).erase(0, acceptor_splice_start);//(0, acceptor_gap_len);



					//new_acceptor_start = acceptor_start + ConvertToNoGappedPos(max_match_acceptor_pos+1, acceptor_ol_queryStr_gap_starts, acceptor_ol_queryStr_gap_ends);

					new_acceptor_start = acceptor_start + ConvertToNoGappedPos(acceptor_pos_smaller+1, acceptor_ol_queryStr_gap_starts, acceptor_ol_queryStr_gap_ends, true);

					//new_donor_end = donor_end - (donor_ol_queryStr.length() - 1 - 

					//	ConvertToNoGappedPos(max_match_donor_pos, donor_ol_queryStr_gap_starts, donor_ol_queryStr_gap_ends) );

					new_donor_end = acceptor_start + 

						ConvertToNoGappedPos(donor_pos_larger, donor_ol_queryStr_gap_starts, donor_ol_queryStr_gap_ends, false) ;

					}



					identity2 = qStr2.length() - CountCharInStr(qStr2, ' ') - CountCharInStr(qStr2, '+');



					if (identity1 > identity2)

					{

						splice_align += donor_align.match_align;

						//for (i=0; i<qStr2.length(); i++)

						//	splice_align += ' ';

						string tmpStr(search_end2+1, ' '); //(qStr2.length(), ' '); //changed, because qStr2 may not be at the very end, if there's gaps at the end of query segment

						splice_align += tmpStr;

						splice_align += acceptor_align.match_align.substr(search_end2+1);

					}

					else

					{

						splice_align += donor_align.match_align.substr(0, search_start1);

						//for (i=0; i<qStr1.length(); i++)

						//	splice_align += ' ';

						string tmpStr(donor_align.match_align.length()-search_start1, ' ');//(qStr1.length(), ' ');

						splice_align += tmpStr;

						splice_align += acceptor_align.match_align;

					}

*/

				}

				else //donor_start...acceptor_start...acceptor_end...donor_end

				{

					const string& qStr2 = acceptor_align.match_align;

					identity2 = qStr2.length() - CountCharInStr(qStr2, ' ') - CountCharInStr(qStr2, '+');



					if (identity1 > identity2)

					{

						splice_align += donor_align.match_align;

						//for (i=0; i<qStr2.length(); i++)

						//	splice_align += ' ';

						//string tmpStr(acceptor_align.match_align.length(), ' ');//(qStr2.length(), ' ');

						int target_bp = acceptor_align.target_align.length() - CountCharInStr(acceptor_align.target_align, '-');

						string tmpStr(target_bp, ' ');

						splice_align += tmpStr;



						donor_splice_end = donor_align.match_align.length()-1;

						acceptor_splice_start = acceptor_align.match_align.length();

					}

					else

					{

						splice_align += donor_align.match_align.substr(0, search_start1);

						//for (i=0; i<qStr1.length(); i++)

						//	splice_align += ' ';

						string donor_tmpStr = donor_align.target_align.substr(search_start1);

						int donor_gap_len = donor_align.target_align.length() - search_start1 - CountCharInStr(donor_tmpStr, '-');

						//int donor_gap_len = donor_align.match_align.length()-search_start1;

						string tmpStr(donor_gap_len, ' ');//(qStr1.length(), ' ');

						splice_align += tmpStr;

						splice_align += qStr2;



						donor_splice_end = search_start1 - 1;

						acceptor_splice_start = 0;



						//also change the alignments themselves!

/*						if (alignment_update)

						{

						//const_cast<string&>(donor_align.query_align).replace(search_start1, 

						//	donor_gap_len, string(donor_gap_len, '-'));

						//const_cast<string&>(donor_align.match_align).replace(search_start1, 

						//	donor_gap_len, string(donor_gap_len, ' '));

							const_cast<string&>(donor_align.query_align).erase(search_start1);

							const_cast<string&>(donor_align.match_align).erase(search_start1);

							const_cast<string&>(donor_align.target_align).erase(search_start1);

						}

*/

						new_donor_end = acceptor_start - 1;

					}

				}

			}

			else // (donor_start >= acceptor_start), cross align?

				//(...acceptor_start...donor_start...acceptor_end ... donor_end)

			//or (...acceptor_start...donor_start ... donor_end...acceptor_end)

			{

				const string& qStr1 = donor_align.match_align;

				identity1 = qStr1.length() - CountCharInStr(qStr1, ' ') - CountCharInStr(qStr1, '+');



				int cur_end = acceptor_end < donor_end? acceptor_end : donor_end;

				FindRealPos_DNA(acceptor_align.query_align, acceptor_start, acceptor_end, acceptor_start, cur_end, 

					search_start2, search_end2);

				const string& qStr2 = acceptor_align.match_align.substr(0, search_end2+1);

				identity2 = qStr2.length() - CountCharInStr(qStr2, ' ') - CountCharInStr(qStr2, '+');



				if (identity1 > identity2)

				{

					splice_align += qStr1;

					//for (i=0; i<qStr2.length(); i++)

					//	splice_align += ' ';

					//string tmpStr(qStr2.length(), ' ');

					string acceptor_tmpStr = acceptor_align.target_align.substr(0, search_end2+1);

					string tmpStr(search_end2+1 - CountCharInStr(acceptor_tmpStr, '-'), ' ');

					splice_align += tmpStr;

					if (acceptor_end > cur_end)

					{

						splice_align += acceptor_align.match_align.substr(search_end2+1);

						acceptor_splice_start = search_end2+1;



/*						if (alignment_update)

						{

						//const_cast<string&>(acceptor_align.query_align).replace(0, search_end2+1, string(search_end2+1, '-'));

						//const_cast<string&>(acceptor_align.match_align).replace(0, search_end2+1, string(search_end2+1, ' '));

							const_cast<string&>(acceptor_align.query_align).erase(0, search_end2+1);

							const_cast<string&>(acceptor_align.match_align).erase(0, search_end2+1);

							const_cast<string&>(acceptor_align.target_align).erase(0, search_end2+1);

						}

*/

						new_acceptor_start = cur_end+1;

					}

					else

					{

						acceptor_splice_start = acceptor_align.match_align.length();

					}



					donor_splice_end = qStr1.length() - 1;

				}

				else

				{

					//for (i=0; i<qStr1.length(); i++)

					//	splice_align += ' ';

					//string tmpStr(qStr1.length(), ' ');

					string tmpStr(donor_align.target_align.length() - CountCharInStr(donor_align.target_align, '-'), ' ');

					splice_align += tmpStr;

					splice_align += acceptor_align.match_align;



					donor_splice_end = -1;

					acceptor_splice_start = 0;

				}

			}





		}



	}



	return true;

}



void ACCP_DONR_graph::ComputeSpliceAlignment_SplitOL(Input_Alignment& donor_align, Input_Alignment& acceptor_align, int donor_start, int donor_end, 

											 int acceptor_start, int acceptor_end, string& splice_align, 

											 int& donor_splice_end, int& acceptor_splice_start, 

											 int& new_donor_end, int& new_acceptor_start)

//donor_splice_end: index on the donor_align.matchStr of the position where the alignment is cut (index is the last position of useful string)

//acceptor_splice_start: index on acceptor_align.matchStr for the cut (index is the first position of useful string)

{

	new_acceptor_start = acceptor_start; //the possibly new acceptor_start position (nextHSP's gene_start)

	new_donor_end = donor_end; //possible new donor_end position (curHSP's gene_end)



	splice_align = "";



	if (donor_start == 0 || acceptor_start == 0) //if either segment is all gaps

	{

		splice_align += donor_align.match_align;

		splice_align += acceptor_align.match_align;

//		cout << "already done!" << "\n";

		donor_splice_end = donor_align.match_align.length()-1; //nothing is cut

		acceptor_splice_start = 0; //nothing is cut



#ifdef DEBUG

		dataMgr.outFile << "donor or acceptor segment is all gap (donor_start:" << donor_start 

			<< "; acceptor_start:" << acceptor_start << ")" << "\n";

#endif

		return;

	}



	//int i;

	int search_start1, search_end1, search_start2, search_end2;

	if (donor_end < acceptor_start) //non overlap

	{

		//revised: now check the trailing / heading gaps in query_align before inserting more gaps!

		//revised again: if there is no overlap and the missing query part is exactly three base pairs,

		//then check donor_prev_nt with acceptor_next_nt against donor_aa_end/acceptor_aa_front, if it's a match, then revise splice_align

		int total_gap_len = acceptor_start - donor_end - 1; //this is the total gap length that we should have

		int last_non_gap_pos_donor_align = donor_align.query_align.find_last_not_of("-");

		int first_non_gap_pos_acceptor_align = acceptor_align.query_align.find_first_not_of("-");

		if (last_non_gap_pos_donor_align == string::npos)

			last_non_gap_pos_donor_align = -1;

		if (first_non_gap_pos_acceptor_align == string::npos)

			first_non_gap_pos_acceptor_align = acceptor_align.query_align.length();



		int need_gap_len = total_gap_len - (donor_align.query_align.length() - 1 - last_non_gap_pos_donor_align + first_non_gap_pos_acceptor_align);



		splice_align += donor_align.match_align;



		//for (i=0; i<acceptor_start-donor_end-1; i++)

		//	splice_align += ' ';

		if (need_gap_len > 0)

		{

			string tmpStr(need_gap_len, ' ');//(acceptor_start-donor_end-1, ' ');

			splice_align += tmpStr;

		}

		splice_align += acceptor_align.match_align;



		donor_splice_end = donor_align.match_align.length()-1;

		acceptor_splice_start = 0;



#ifdef DEBUG

		dataMgr.outFile << "donor and acceptor segments no overlap, should have " << total_gap_len << " gaps; filled in " << need_gap_len << " gaps" << "\n";

			//tmpStr.length() << " gaps" << "\n";

#endif

	}

	else //overlap? ( acceptor_start ... donor_end)

	{

		int identity1, identity2;

		if (donor_start > acceptor_end) //cross align! so only one aligns, the other is turned to gaps 

			//(acceptor_start...acceptor_end...donor_start ... donor_end)

		{

			const string& qStr1 = donor_align.match_align;

			const string& qStr2 = acceptor_align.match_align;



			identity1 = qStr1.length() - CountCharInStr(qStr1, ' ') - CountCharInStr(qStr1, '+');

			identity2 = qStr2.length() - CountCharInStr(qStr2, ' ') - CountCharInStr(qStr2, '+');



			if (identity1 > identity2)

			{

				splice_align += qStr1;

				//for (i=0; i<qStr2.length(); i++)

				//	splice_align += ' ';

				string tmpStr(qStr2.length(), ' ');

				splice_align += tmpStr;



				donor_splice_end = qStr1.length()-1;

				acceptor_splice_start = qStr2.length(); //acceptor is all cut (to the end)



				dataMgr.outFile << "acceptor is turned to gap" << "\n";

			}

			else

			{

				//for (i=0; i<qStr1.length(); i++)

				//	splice_align += ' ';

				string tmpStr(qStr1.length(), ' ');

				splice_align += tmpStr;

				splice_align += qStr2;



				donor_splice_end = -1; //donor is all cut

				acceptor_splice_start = 0;



				dataMgr.outFile << "donor is turned to gap" << "\n";

			}



#ifdef DEBUG

			dataMgr.outFile << " (cross align, one is turned to gaps)" << "\n";

#endif

		}

		else //(...donor_start...acceptor_start...acceptor_end ... donor_end)

			//or (...donor_start...acceptor_start ... donor_end...acceptor_end)

			//or (...acceptor_start...donor_start...acceptor_end ... donor_end)

			//or (...acceptor_start...donor_start ... donor_end...acceptor_end)

		{

			if (donor_start < acceptor_start)//(...donor_start...acceptor_start...acceptor_end ... donor_end)

			//or (...donor_start...acceptor_start ... donor_end...acceptor_end)

			{

				

				FindRealPos_DNA(donor_align.query_align, donor_start, donor_end, acceptor_start, donor_end, 

					search_start1, search_end1);

				const string& qStr1 = donor_align.match_align.substr(search_start1, search_end1-search_start1+1);

				identity1 = qStr1.length() - CountCharInStr(qStr1, ' ') - CountCharInStr(qStr1, '+');

				int qStr1_left_len = donor_align.query_align.length() - search_end1 - 1; //length of stuff after search_end1 (gaps)



				if (acceptor_end >= donor_end) //donor_start...acceptor_start...donor_end...acceptor_end

				{

					//the overlapping part (query segment) is [acceptor_start, donor_end]

					FindRealPos_DNA(acceptor_align.query_align, acceptor_start, acceptor_end, acceptor_start, donor_end, 

						search_start2, search_end2);

					const string& qStr2 = acceptor_align.match_align.substr(search_start2, search_end2-search_start2+1);

					int qStr2_left_len = search_start2; //length of stuff before search_start2 (gaps)



					//this is only used here

					const string& donor_ol_queryStr = donor_align.query_align.substr(search_start1, search_end1-search_start1+1);

					const string& acceptor_ol_queryStr = acceptor_align.query_align.substr(search_start2, search_end2-search_start2+1);



					vector<int> donor_ol_queryStr_gap_starts, donor_ol_queryStr_gap_ends;

					vector<int> acceptor_ol_queryStr_gap_starts, acceptor_ol_queryStr_gap_ends;

					GetGaps(donor_ol_queryStr, donor_ol_queryStr_gap_starts, donor_ol_queryStr_gap_ends);

					GetGaps(acceptor_ol_queryStr, acceptor_ol_queryStr_gap_starts, acceptor_ol_queryStr_gap_ends);



#ifdef DEBUG

					dataMgr.outFile << "query:" << donor_ol_queryStr << "\n" << "match:" << qStr1 << "\n";

					dataMgr.outFile << "query:" << acceptor_ol_queryStr << "\n" << "match:" << qStr2 << "\n";

					vector<int>::iterator tmp_it1, tmp_it2=donor_ol_queryStr_gap_ends.begin();

					dataMgr.outFile << "donor gaps: ";

					for (tmp_it1 = donor_ol_queryStr_gap_starts.begin(); tmp_it1 != donor_ol_queryStr_gap_starts.end(); tmp_it1++, tmp_it2++)

						dataMgr.outFile << *tmp_it1 << "-" << *tmp_it2 << "; ";

					dataMgr.outFile << "\nacceptor gaps: ";

					tmp_it2 = acceptor_ol_queryStr_gap_ends.begin();

					for (tmp_it1 = acceptor_ol_queryStr_gap_starts.begin(); tmp_it1 != acceptor_ol_queryStr_gap_starts.end(); tmp_it1++, tmp_it2++)

						dataMgr.outFile << *tmp_it1 << "-" << *tmp_it2 << "; ";

					dataMgr.outFile << "\n";

#endif



					vector<int> donor_ol_matchStr_match_starts, donor_ol_matchStr_match_ends;

					vector<int> acceptor_ol_matchStr_match_starts, acceptor_ol_matchStr_match_ends;

					vector<int> donor_match_end_ids, acceptor_match_start_ids;



					int total_ol_id1 = GetMatches(qStr1, donor_ol_matchStr_match_starts, donor_ol_matchStr_match_ends, donor_match_end_ids);

					int total_ol_id2 = GetMatches(qStr2, acceptor_ol_matchStr_match_starts, acceptor_ol_matchStr_match_ends, acceptor_match_start_ids);



					//now try each match_end on donor_ol_matchStr (qStr1) to be used as split point

					vector<int>::iterator d_start_it=donor_ol_matchStr_match_starts.begin();

					vector<int>::iterator d_end_it = donor_ol_matchStr_match_ends.begin();

					int cur_match_count = 0; //donor side match count

					//int max_match_count=-1, max_match_donor_pos=-1, max_match_acceptor_pos=-1; //initialize

					int max_match_count=total_ol_id2, max_match_donor_pos=-1, max_match_acceptor_pos=-1; //initialize

					//first region_end = -1, region_start = 0 (to simulate all overlapping part uses acceptor)

					for (; d_start_it != donor_ol_matchStr_match_starts.end(); d_start_it++, d_end_it++)

					{

						int region_start = (*d_start_it);

						int region_end = (*d_end_it);

						cur_match_count += region_end - region_start + 1;



						int donor_queryPos = ConvertToNoGappedPos(region_end, donor_ol_queryStr_gap_starts, donor_ol_queryStr_gap_ends, false);

						int acceptor_matchPos = ConvertToGappedPos(donor_queryPos, acceptor_ol_queryStr_gap_starts, acceptor_ol_queryStr_gap_ends);



						int acceptor_match_count = ComputeBackId(acceptor_matchPos, acceptor_ol_matchStr_match_starts, 

							acceptor_ol_matchStr_match_ends, acceptor_match_start_ids, total_ol_id2);



						int total_match_count = cur_match_count + acceptor_match_count;

#ifdef DEBUG

						dataMgr.outFile << "donor matchPos (abs): " << region_end << "; converted to queryPos (all match): " << donor_queryPos;

						dataMgr.outFile << "; acceptor matchPos (abs): " << acceptor_matchPos << "\n";

						dataMgr.outFile << "cur max: " << max_match_count << "; this time: " << total_match_count << "\n";

#endif

						

						if (total_match_count > max_match_count)

						{

							max_match_count = total_match_count;

							max_match_donor_pos = region_end;

							max_match_acceptor_pos = acceptor_matchPos;

#ifdef DEBUG

							dataMgr.outFile << "cur donor pos: " << max_match_donor_pos  << "; cur acc pos: " << max_match_acceptor_pos << "\n";

#endif

						}

					}



					//now also compute from the other side (acceptor side), if both come to the same max_match_count but at different positions, keep matched stuff at both sides

					cur_match_count = 0; //acceptor side match count

					int max_match_count_acc=total_ol_id1, max_match_donor_pos_acc=donor_ol_queryStr.length(), max_match_acceptor_pos_acc=acceptor_ol_queryStr.length(); //initialize

					//first region_end = -1, region_start = 0 (to simulate all overlapping part uses acceptor)

					vector<int>::reverse_iterator a_start_it = acceptor_ol_matchStr_match_starts.rbegin(); 

					vector<int>::reverse_iterator a_end_it = acceptor_ol_matchStr_match_ends.rbegin();

					for (; a_start_it != acceptor_ol_matchStr_match_starts.rend(); a_start_it++, a_end_it++)

					{

						int region_start = (*a_start_it);

						int region_end = (*a_end_it);

						cur_match_count += region_end - region_start + 1;



						int acceptor_queryPos = ConvertToNoGappedPos(region_start, acceptor_ol_queryStr_gap_starts, acceptor_ol_queryStr_gap_ends, false);

						int donor_matchPos = ConvertToGappedPos(acceptor_queryPos, donor_ol_queryStr_gap_starts, donor_ol_queryStr_gap_ends);



						int donor_match_count = ComputeFrontId(donor_matchPos, donor_ol_matchStr_match_starts, 

							donor_ol_matchStr_match_ends, donor_match_end_ids, total_ol_id1);



						int total_match_count = cur_match_count + donor_match_count;

#ifdef DEBUG

						dataMgr.outFile << "acceptor matchPos (abs): " << region_start << "; converted to queryPos (all match): " << acceptor_queryPos;

						dataMgr.outFile << "; donor matchPos (abs): " << donor_matchPos << "\n";

						dataMgr.outFile << "cur max: " << max_match_count_acc << "; this time: " << total_match_count << "\n";

#endif

						

						if (total_match_count > max_match_count_acc)

						{

							max_match_count_acc = total_match_count;

							max_match_acceptor_pos_acc = region_start;

							max_match_donor_pos_acc = donor_matchPos;

#ifdef DEBUG

							dataMgr.outFile << "cur donor pos: " << max_match_donor_pos_acc  << "; cur acc pos: " << max_match_acceptor_pos_acc << "\n";

#endif

						}

					}



					int donor_gap_len = qStr1.length() - max_match_donor_pos - 1 + qStr1_left_len;

					int acceptor_gap_len = qStr2_left_len + max_match_acceptor_pos+1;

					splice_align += donor_align.match_align.substr(0, search_start1+max_match_donor_pos+1);

					string tmpStr(donor_gap_len + acceptor_gap_len, ' ');

					splice_align += tmpStr;

					splice_align += acceptor_align.match_align.substr(acceptor_gap_len);



					//modified when max_match_count == max_match_count_acc (they must be equal?!), 

					//if the two splicing positions for each max_match_count are different:

					//donor_splice_end is at the later of the two positions (one for each max_match_count)

					//similarly for acceptor_splice_start



					int donor_pos_larger, acceptor_pos_smaller;

					max_match_donor_pos_acc--; //convert to the way of splicing using donor side computing?

					max_match_acceptor_pos_acc--;

					if (max_match_donor_pos == max_match_donor_pos_acc) //both methods splice at the same position

					{

					donor_pos_larger = max_match_donor_pos;

					acceptor_pos_smaller = max_match_acceptor_pos;

					//donor_splice_end = max_match_donor_pos;

					//acceptor_splice_start = max_match_acceptor_pos+1;

					donor_splice_end = search_start1 + max_match_donor_pos;

					acceptor_splice_start = acceptor_gap_len;

					}

					else

					{

						donor_pos_larger = max_match_donor_pos < max_match_donor_pos_acc ? max_match_donor_pos_acc : max_match_donor_pos;

						acceptor_pos_smaller = max_match_acceptor_pos < max_match_acceptor_pos_acc ? max_match_acceptor_pos : max_match_acceptor_pos_acc;

						donor_splice_end = search_start1 + donor_pos_larger;

						acceptor_splice_start = qStr2_left_len + acceptor_pos_smaller + 1;

					}



					

					//also change the alignments themselves!

//					if (alignment_update)

//					{

					//const_cast<string&>(donor_align.query_align).replace(search_start1+max_match_donor_pos+1, 

					//	donor_gap_len, string(donor_gap_len, '-'));

					//const_cast<string&>(donor_align.match_align).replace(search_start1+max_match_donor_pos+1, 

					//	donor_gap_len, string(donor_gap_len, ' '));

					//const_cast<string&>(acceptor_align.query_align).replace(0, acceptor_gap_len, string(acceptor_gap_len, '-'));

					//const_cast<string&>(acceptor_align.match_align).replace(0, acceptor_gap_len, string(acceptor_gap_len, ' '));

						const_cast<string&>(donor_align.query_align).erase(donor_splice_end+1);//(search_start1+max_match_donor_pos+1);

						const_cast<string&>(donor_align.match_align).erase(donor_splice_end+1);//(search_start1+max_match_donor_pos+1);

						const_cast<string&>(donor_align.target_align).erase(donor_splice_end+1);//(search_start1+max_match_donor_pos+1);

						const_cast<string&>(acceptor_align.query_align).erase(0, acceptor_splice_start);//(0, acceptor_gap_len);

						const_cast<string&>(acceptor_align.match_align).erase(0, acceptor_splice_start);//(0, acceptor_gap_len);

						const_cast<string&>(acceptor_align.target_align).erase(0, acceptor_splice_start);//(0, acceptor_gap_len);



					//new_acceptor_start = acceptor_start + ConvertToNoGappedPos(max_match_acceptor_pos+1, acceptor_ol_queryStr_gap_starts, acceptor_ol_queryStr_gap_ends);

					new_acceptor_start = acceptor_start + ConvertToNoGappedPos(acceptor_pos_smaller+1, acceptor_ol_queryStr_gap_starts, acceptor_ol_queryStr_gap_ends, true);

					//new_donor_end = donor_end - (donor_ol_queryStr.length() - 1 - 

					//	ConvertToNoGappedPos(max_match_donor_pos, donor_ol_queryStr_gap_starts, donor_ol_queryStr_gap_ends) );

					new_donor_end = acceptor_start + 

						ConvertToNoGappedPos(donor_pos_larger, donor_ol_queryStr_gap_starts, donor_ol_queryStr_gap_ends, false) ;

//					}



/*					identity2 = qStr2.length() - CountCharInStr(qStr2, ' ') - CountCharInStr(qStr2, '+');



					if (identity1 > identity2)

					{

						splice_align += donor_align.match_align;

						//for (i=0; i<qStr2.length(); i++)

						//	splice_align += ' ';

						string tmpStr(search_end2+1, ' '); //(qStr2.length(), ' '); //changed, because qStr2 may not be at the very end, if there's gaps at the end of query segment

						splice_align += tmpStr;

						splice_align += acceptor_align.match_align.substr(search_end2+1);

					}

					else

					{

						splice_align += donor_align.match_align.substr(0, search_start1);

						//for (i=0; i<qStr1.length(); i++)

						//	splice_align += ' ';

						string tmpStr(donor_align.match_align.length()-search_start1, ' ');//(qStr1.length(), ' ');

						splice_align += tmpStr;

						splice_align += acceptor_align.match_align;

					}

*/

				}

				else //donor_start...acceptor_start...acceptor_end...donor_end

				{

					const string& qStr2 = acceptor_align.match_align;

					identity2 = qStr2.length() - CountCharInStr(qStr2, ' ') - CountCharInStr(qStr2, '+');



					if (identity1 > identity2)

					{

						splice_align += donor_align.match_align;

						//for (i=0; i<qStr2.length(); i++)

						//	splice_align += ' ';

						string tmpStr(acceptor_align.match_align.length(), ' ');//(qStr2.length(), ' ');

						splice_align += tmpStr;



						donor_splice_end = donor_align.match_align.length()-1;

						acceptor_splice_start = acceptor_align.match_align.length();

					}

					else

					{

						splice_align += donor_align.match_align.substr(0, search_start1);

						//for (i=0; i<qStr1.length(); i++)

						//	splice_align += ' ';

						int donor_gap_len = donor_align.match_align.length()-search_start1;

						string tmpStr(donor_gap_len, ' ');//(qStr1.length(), ' ');

						splice_align += tmpStr;

						splice_align += qStr2;



						donor_splice_end = search_start1 - 1;

						acceptor_splice_start = 0;



						//also change the alignments themselves!

//						if (alignment_update)

//						{

						//const_cast<string&>(donor_align.query_align).replace(search_start1, 

						//	donor_gap_len, string(donor_gap_len, '-'));

						//const_cast<string&>(donor_align.match_align).replace(search_start1, 

						//	donor_gap_len, string(donor_gap_len, ' '));

							const_cast<string&>(donor_align.query_align).erase(search_start1);

							const_cast<string&>(donor_align.match_align).erase(search_start1);

							const_cast<string&>(donor_align.target_align).erase(search_start1);

//						}



						new_donor_end = acceptor_start - 1;

					}

				}

			}

			else // (donor_start >= acceptor_start), cross align?

				//(...acceptor_start...donor_start...acceptor_end ... donor_end)

			//or (...acceptor_start...donor_start ... donor_end...acceptor_end)

			{

				const string& qStr1 = donor_align.match_align;

				identity1 = qStr1.length() - CountCharInStr(qStr1, ' ') - CountCharInStr(qStr1, '+');



				int cur_end = acceptor_end < donor_end? acceptor_end : donor_end;

				FindRealPos_DNA(acceptor_align.query_align, acceptor_start, acceptor_end, acceptor_start, cur_end, 

					search_start2, search_end2);

				const string& qStr2 = acceptor_align.match_align.substr(0, search_end2+1);

				identity2 = qStr2.length() - CountCharInStr(qStr2, ' ') - CountCharInStr(qStr2, '+');



				if (identity1 > identity2)

				{

					splice_align += qStr1;

					//for (i=0; i<qStr2.length(); i++)

					//	splice_align += ' ';

					string tmpStr(qStr2.length(), ' ');

					splice_align += tmpStr;

					if (acceptor_end > cur_end)

					{

						splice_align += acceptor_align.match_align.substr(search_end2+1);

						acceptor_splice_start = search_end2+1;



//						if (alignment_update)

//						{

						//const_cast<string&>(acceptor_align.query_align).replace(0, search_end2+1, string(search_end2+1, '-'));

						//const_cast<string&>(acceptor_align.match_align).replace(0, search_end2+1, string(search_end2+1, ' '));

							const_cast<string&>(acceptor_align.query_align).erase(0, search_end2+1);

							const_cast<string&>(acceptor_align.match_align).erase(0, search_end2+1);

							const_cast<string&>(acceptor_align.target_align).erase(0, search_end2+1);

//						}



						new_acceptor_start = cur_end+1;

					}

					else

					{

						acceptor_splice_start = acceptor_align.match_align.length();

					}



					donor_splice_end = qStr1.length() - 1;

				}

				else

				{

					//for (i=0; i<qStr1.length(); i++)

					//	splice_align += ' ';

					string tmpStr(qStr1.length(), ' ');

					splice_align += tmpStr;

					splice_align += acceptor_align.match_align;



					donor_splice_end = -1;

					acceptor_splice_start = 0;

				}

			}





		}



	}



}






bool ACCP_DONR_graph::CutHeadOrTrailGapsInAlignment(Input_Alignment& align, int& gene_s, int& gene_e, int& hsp_s, int& hsp_e)

{

		//cut off possible gapped heading/trailing portions of HSPs (trailing or heading gaps at either end of HSP)

		string& tStr = align.target_align;

		string& qStr = align.query_align;

		string& mStr = align.match_align;



		int t_pos_start, t_pos_end, q_pos_start, q_pos_end;

		t_pos_start = tStr.find_first_not_of("-*");
		if (t_pos_start == string::npos)
			return false;

		t_pos_end = tStr.find_last_not_of("-*");

		q_pos_start = qStr.find_first_not_of("-*");
		if (q_pos_start == string::npos)
			return false;

		q_pos_end = qStr.find_last_not_of("-*");



#ifdef DEBUG

		dataMgr.outFile << "tStr: " << tStr << "\n" << "t_pos_start: " << t_pos_start << "; t_pos_end: " << t_pos_end << "\n";

		dataMgr.outFile << "qStr: " << qStr << "\n" << "q_pos_start: " << q_pos_start << "; q_pos_end: " << q_pos_end << "\n";

#endif

		int start_len_erased = 0;

		if (q_pos_start < t_pos_start)
			start_len_erased = t_pos_start;
		else
			start_len_erased = q_pos_start;

		if (start_len_erased)
		{
			string tStr_cut = tStr.substr(0, start_len_erased);
			string qStr_cut = qStr.substr(0, start_len_erased);
			int count_t = start_len_erased - CountCharInStr(tStr_cut, '-'); //number of bp that needs to be advanced
			int count_q = start_len_erased - CountCharInStr(qStr_cut, '-'); //number of query a.a. that needs to be advanced
			if (start_site > 0)
				hsp_s += count_t*3;
			else
				hsp_e -= count_t*3;
			gene_s += count_q; //t_pos_start;
			const_cast<string&>(qStr).erase(0, start_len_erased);
			const_cast<string&>(tStr).erase(0, start_len_erased);
			const_cast<string&>(mStr).erase(0, start_len_erased);
		}

		int end_len_erased = 0;
		
		int remain_len = tStr.length();
		if (q_pos_end < t_pos_end)
			end_len_erased = remain_len - (q_pos_end - start_len_erased + 1);
		else
			end_len_erased = remain_len - (t_pos_end - start_len_erased + 1);

		if (end_len_erased)
		{
			string tStr_cut = tStr.substr(remain_len - end_len_erased, end_len_erased);
			string qStr_cut = qStr.substr(remain_len - end_len_erased, end_len_erased);
			int count_t = end_len_erased - CountCharInStr(tStr_cut, '-');
			int count_q = end_len_erased - CountCharInStr(qStr_cut, '-');
			if (start_site > 0)
				hsp_e -= count_t*3; //(tStr.length() - 1 - q_pos_end)*3;
			else
				hsp_s += count_t*3; //(tStr.length() - 1 - q_pos_end)*3;
			gene_e -= count_q; //q_pos_end - t_pos_end; //qStr.length() - 1- t_pos_end;
			const_cast<string&>(tStr).erase(remain_len - end_len_erased);
			const_cast<string&>(qStr).erase(remain_len - end_len_erased);
			const_cast<string&>(mStr).erase(remain_len - end_len_erased);
		}

	return true;


}





bool ACCP_DONR_graph::GetSpliceSegments2(//multimap<int, pair<int, int> >& donor_segments, 

										 vector<HSP_Gene_Pair*>& HSPs,

		int& first_segment_start, int& last_segment_start, int& last_segment_end)//, 

		//vector<int>& segment_hsp_start, int& last_hsp_start)

{

	bool lastLeafIsExon = false;



	int prev_hsp_ID = -1;



	hspID_to_hspIndex_map.clear();

	int cur_hsp_index;

	//int last_hsp_id;

	int last_hsp_end;



	if (start_site > 0) //positive strand HSP

	{

		cur_hsp_index = HSPs.size()-1;



		vector<HSP_Gene_Pair*>::reverse_iterator vec_revIt, last_revIt, next_revIt;

		vec_revIt = HSPs.rbegin();

		next_revIt = HSPs.rbegin(); next_revIt++;

		first_segment_start = (*vec_revIt)->HSP_start;

	

		last_revIt = HSPs.rend();

		last_revIt--;

		for (; vec_revIt != HSPs.rend(); vec_revIt++, next_revIt++) //for each HSP, in ascending order!

		//while (vec_revIt != HSPs.rend())

		{

			//hspID_to_hspIndex_map.insert(map<int, int>::value_type((*vec_revIt)->ID, cur_hsp_index--));

			//cout << "<" << (*vec_revIt)->ID << "," << cur_hsp_id-1 << "> inserted to hspID_to_hspIndex_map" << "\n";



			if (vec_revIt == HSPs.rbegin())

				FIRSTHSP = true;

			else

				FIRSTHSP = false;



			if (next_revIt == HSPs.rend())

			{

				LASTHSP = true;

				//last_hsp_id = (*vec_revIt)->ID;

				last_hsp_end = (*vec_revIt)->HSP_end;

			}

			else

				LASTHSP = false;			

				

			//last_hsp_start =  (*vec_revIt)->HSP_start;

			//if (vec_revIt != HSPs.rbegin() && vec_revIt != last_revIt)

				//MINOBJS = Min(MIN_INTRON_LEN, MIN_INTERNAL_EXON_LEN); //change MINOBJS for middle HSPs



			vector<int> tgt_gap_starts, tgt_gap_ends;



			MaxItem = LoadHSPAlignScores(*vec_revIt, tgt_gap_starts, tgt_gap_ends);

			InitialiseTreeData();

			InitialiseWeights();



			//Raw = new Tree[1];

			Raw = (Tree *)calloc(1, sizeof(Tree));

		    Raw[0] = FormTree(0, MaxItem);

#ifdef DEBUG

			PrintTree(Raw[0]);

			dataMgr.outFile << bufDecTree << "\n";

#endif

			

			Prune(Raw[0]);

#ifdef DEBUG

			PrintTree(Raw[0]);

			dataMgr.outFile << bufDecTree << "\n";

#endif



			bool isFirstLeaf = true;

			GetSegments(Raw[0], tgt_gap_starts, tgt_gap_ends, (*vec_revIt)->HSP_start, 

				//donor_segments, 

				donor_segments_pair, donor_acceptor_HSP_ID, //acceptor_HSP_ID, 

				(*vec_revIt)->ID, //prev_hsp_ID, 

				last_segment_start, last_segment_end, lastLeafIsExon);//, segment_hsp_start, isFirstLeaf);



#ifdef DEBUG

	//multimap<int, pair<int, int> >::iterator seg_map_It;

	vector<SegmentsInThreeBounds>::iterator seg_map_It;

	dataMgr.outFile << "FROM CUR HSP TREE: " << "\n";

	for (seg_map_It = donor_segments_pair.begin(); seg_map_It != donor_segments_pair.end(); seg_map_It++)

	{

		//dataMgr.outFile << "search: " << (*seg_map_It).first-1 << " to " << (*seg_map_It).second.first 

		//	<< "; " << (*seg_map_It).first << " to " << (*seg_map_It).second.second << "\n";

		dataMgr.outFile << "search: " << (*seg_map_It).exon_seg_end-1 << " to " << (*seg_map_It).exon_seg_start 

			<< "; " << (*seg_map_It).exon_seg_end << " to " << (*seg_map_It).intron_seg_end << "\n";

	}

	dataMgr.outFile << "last_segment_start: " << last_segment_start << " - last_segment_end: " << last_segment_end << "\n";

#endif



			if (next_revIt != HSPs.rend())

			{

				if (lastLeafIsExon)

				{

					if ((*next_revIt)->HSP_start > last_segment_end)

					{

						//donor_segments.insert(multimap<int, pair<int, int> >::value_type(last_segment_end, 

						//	pair<int, int>(last_segment_start, (*next_revIt)->HSP_start - 1)));

						donor_segments_pair.push_back( SegmentsInThreeBounds(last_segment_start,last_segment_end, (*next_revIt)->HSP_start - 1) );

						//donor_acceptor_HSP_ID.push_back((*vec_revIt)->ID);

						//donor_acceptor_HSP_ID.push_back(prev_hsp_ID);



						lastLeafIsExon = false;

					}

					//gap_centers.push_back(( last_segment_end + (*next_revIt)->HSP_start - 1 )/2 );



					//segment_hsp_start.push_back(last_hsp_start);



					//lastLeafIsExon = false;					

				}

				else

				{

					//if (!donor_segments.empty())

					if (!donor_segments_pair.empty())

					{

					//multimap<int, pair<int, int> >::reverse_iterator seg_it = donor_segments.rbegin();

					//const_cast<int&>((*seg_it).second.second) = (*next_revIt)->HSP_start - 1;

						donor_segments_pair[donor_segments_pair.size()-1].intron_seg_end = (*next_revIt)->HSP_start - 1;

					}

					//gap_centers.back() = ( (*seg_it).second.second + (*seg_it).first )/2;

				}				

			}

			else //if current HSP is the last one

			{

				if (!lastLeafIsExon) //if last leaf is intron, delete it

				{

					//if (!donor_segments.empty())

					if (!donor_segments_pair.empty())

					{

					//multimap<int, pair<int, int> >::iterator seg_it = donor_segments.end();

					//seg_it--;

					//donor_segments.erase(seg_it);

						donor_segments_pair.pop_back();

						//donor_acceptor_HSP_ID.pop_back();

					}

				}

			}



#ifdef DEBUG

	dataMgr.outFile << "AFTER CONNECTING WITH NEXT HSP: " << "\n";

	for (seg_map_It = donor_segments_pair.begin(); seg_map_It != donor_segments_pair.end(); seg_map_It++)

	{

		//dataMgr.outFile << "search: " << (*seg_map_It).first-1 << " to " << (*seg_map_It).second.first 

		//	<< "; " << (*seg_map_It).first << " to " << (*seg_map_It).second.second << "\n";

		dataMgr.outFile << "search: " << (*seg_map_It).exon_seg_end-1 << " to " << (*seg_map_It).exon_seg_start 

			<< "; " << (*seg_map_It).exon_seg_end << " to " << (*seg_map_It).intron_seg_end << "\n";

	}

	dataMgr.outFile << "last_segment_start: " << last_segment_start << " - last_segment_end: " << last_segment_end << "\n";

#endif



			//clean up

			for (int i=0; i<=MaxItem; i++)

				delete [] Item[i];

			delete [] Item;

			delete [] HSPAlignScores;



			//delete [] Raw;

			ReleaseTree(Raw[0]);



			prev_hsp_ID = (*vec_revIt)->ID;

		}



	}

	else //negative HSP

	{

		cur_hsp_index = 0;



		vector<HSP_Gene_Pair*>::iterator vec_It, next_It;

		vec_It = HSPs.begin();

		next_It = HSPs.begin(); next_It++;

		first_segment_start = -(*vec_It)->HSP_end;



		for (; vec_It != HSPs.end(); vec_It++, next_It++) //for each HSP, in ascending order (1st HSP is closest to start)!

		//while (vec_It != HSPs.end())

		{

			//hspID_to_hspIndex_map.insert(map<int, int>::value_type((*vec_It)->ID, cur_hsp_index++));

			//cout << "<" << (*vec_It)->ID << "," << cur_hsp_id-1 << "> inserted to hspID_to_hspIndex_map" << "\n";



			if (vec_It == HSPs.begin())

				FIRSTHSP = true;

			else

				FIRSTHSP = false;



			if (next_It == HSPs.end())

			{

				LASTHSP = true;

				//last_hsp_id = (*vec_It)->ID;

				last_hsp_end = -(*vec_It)->HSP_start;

			}

			else

				LASTHSP = false;



			//last_hsp_start = -(*vec_It)->HSP_end;



			vector<int> tgt_gap_starts, tgt_gap_ends;



			MaxItem = LoadHSPAlignScores(*vec_It, tgt_gap_starts, tgt_gap_ends);

			InitialiseTreeData();

			InitialiseWeights();



			Raw = (Tree *)calloc(1, sizeof(Tree));

		    Raw[0] = FormTree(0, MaxItem);

#ifdef DEBUG

			PrintTree(Raw[0]);

			dataMgr.outFile << bufDecTree << "\n";

#endif



			Prune(Raw[0]);

#ifdef DEBUG

			PrintTree(Raw[0]);

			dataMgr.outFile << bufDecTree << "\n";

#endif



			bool isFirstLeaf = true;

			GetSegments(Raw[0], tgt_gap_starts, tgt_gap_ends, -(*vec_It)->HSP_end, 

				//donor_segments, 

				donor_segments_pair, donor_acceptor_HSP_ID, //acceptor_HSP_ID, 

				(*vec_It)->ID, //prev_hsp_ID, 

				last_segment_start, last_segment_end, lastLeafIsExon);//, segment_hsp_start, isFirstLeaf);



#ifdef DEBUG

	//multimap<int, pair<int, int> >::iterator seg_map_It;

	vector<SegmentsInThreeBounds>::iterator seg_map_It;

	dataMgr.outFile << "FROM CUR HSP TREE: " << "\n";

	for (seg_map_It = donor_segments_pair.begin(); seg_map_It != donor_segments_pair.end(); seg_map_It++)

	{

		//dataMgr.outFile << "search: " << (*seg_map_It).first-1 << " to " << (*seg_map_It).second.first 

		//	<< "; " << (*seg_map_It).first << " to " << (*seg_map_It).second.second << "\n";

		dataMgr.outFile << "search: " << (*seg_map_It).exon_seg_end-1 << " to " << (*seg_map_It).exon_seg_start 

			<< "; " << (*seg_map_It).exon_seg_end << " to " << (*seg_map_It).intron_seg_end << "\n";

	}

	dataMgr.outFile << "last_segment_start: " << last_segment_start << " - last_segment_end: " << last_segment_end << "\n";

#endif



			if (next_It != HSPs.end())

			{

				if (lastLeafIsExon)

				{

					if ((*next_It)->HSP_start > last_segment_end)

					{

						//donor_segments.insert(multimap<int, pair<int, int> >::value_type(last_segment_end, 

						//	pair<int, int>(last_segment_start, -(*next_It)->HSP_end - 1)));

						donor_segments_pair.push_back(SegmentsInThreeBounds(last_segment_start, last_segment_end, -(*next_It)->HSP_end - 1));

						//donor_acceptor_HSP_ID.push_back((*vec_It)->ID);

						//donor_acceptor_HSP_ID.push_back(prev_hsp_ID);

					}

					//gap_centers.push_back(( last_segment_end - (*next_It)->HSP_end - 1 )/2 );

					//segment_hsp_start.push_back(last_hsp_start);



					lastLeafIsExon = false;

				}

				else

				{

					if (!donor_segments_pair.empty())

					{

					//multimap<int, pair<int, int> >::reverse_iterator seg_it = donor_segments.rbegin();

					//const_cast<int&>((*seg_it).second.second) = -(*next_It)->HSP_end - 1;

						donor_segments_pair[donor_segments_pair.size()-1].intron_seg_end = -(*next_It)->HSP_end - 1;

					}

					//gap_centers.back() = ( (*seg_it).second.second + (*seg_it).first )/2;

				}

			}

			else

			{

				if (!lastLeafIsExon)

				{

					if (!donor_segments_pair.empty())

					{

					//multimap<int, pair<int, int> >::iterator seg_it = donor_segments.end();

					//seg_it--;

					//donor_segments.erase(seg_it);

						donor_segments_pair.pop_back();

						//donor_acceptor_HSP_ID.pop_back();

					}

				}

			}



#ifdef DEBUG

	dataMgr.outFile << "AFTER CONNECTING WITH NEXT HSP: " << "\n";

	for (seg_map_It = donor_segments_pair.begin(); seg_map_It != donor_segments_pair.end(); seg_map_It++)

	{

		//dataMgr.outFile << "search: " << (*seg_map_It).first-1 << " to " << (*seg_map_It).second.first 

		//	<< "; " << (*seg_map_It).first << " to " << (*seg_map_It).second.second << "\n";

		dataMgr.outFile << "search: " << (*seg_map_It).exon_seg_end-1 << " to " << (*seg_map_It).exon_seg_start 

			<< "; " << (*seg_map_It).exon_seg_end << " to " << (*seg_map_It).intron_seg_end << "\n";

	}

	dataMgr.outFile << "last_segment_start: " << last_segment_start << " - last_segment_end: " << last_segment_end << "\n";

#endif



			//clean up

			for (int i=0; i<=MaxItem; i++)

				delete [] Item[i];

			delete [] Item;

			delete [] HSPAlignScores;



			//delete [] Raw;

			ReleaseTree(Raw[0]);



			prev_hsp_ID = (*vec_It)->ID;

		}



	}



	//donor_acceptor_HSP_ID.push_back(last_hsp_id);//need an extra one at the end, for last_segment

	//donor_acceptor_HSP_ID.push_back(prev_hsp_ID);

	//if (!donor_acceptor_HSP_ID.empty()) //no need, must have already been inserted!

	//	donor_acceptor_HSP_ID.back().insert(prev_hsp_ID);



	//MODIFIED: now we also need to check if the first HSP and last HSP still are valid (contain exon segments),

	//otherwise we need to re-look for gene_start and/or gene_end again

	int all_seg_start, fst_seg_hsp_ID, last_seg_hsp_ID;

	if (!donor_segments_pair.empty())

	{

		all_seg_start = donor_segments_pair.front().exon_seg_start;

		//fst_seg_hsp_ID = donor_acceptor_HSP_ID.front();

	}

	else

	{

		all_seg_start = last_segment_start;

		//fst_seg_hsp_ID = prev_hsp_ID; //the same as donor_acceptor_HSP_ID.back()

	}

	//last_seg_hsp_ID = prev_hsp_ID; //the same as donor_acceptor_HSP_ID.back()



	fst_seg_hsp_ID = donor_acceptor_HSP_ID.front().front();

	last_seg_hsp_ID = donor_acceptor_HSP_ID.back().back();



	vector<HSP_Gene_Pair*>::iterator vec_It, vec_next_It;

	if (all_seg_start != 0)

	{

		

		if (all_seg_start > first_segment_start) //need to fix HSPs

		{

			first_segment_start = donor_segments_pair.front().exon_seg_start;

			//fix HSPs

			if (start_site > 0)

			{

				vec_It = HSPs.end(); vec_It--; //points to the last HSP, must exist				

				while (vec_It != HSPs.begin() && (*vec_It)->ID != fst_seg_hsp_ID)

				{

					vec_next_It = vec_It; vec_next_It--; //get to the next position

					HSPs.erase(vec_It);

					vec_It = vec_next_It;

				}

				//now vec_It points to the first useful HSP, use this to calc the updated start_site

				//do this in LoadData4()

				//CalcStartPos((*vec_It)->HSP_start, (*vec_It)->gene_start, true, chr_seq);

			}

			else

			{

				vec_It = HSPs.begin();

				while (vec_It != HSPs.end() && (*vec_It)->ID != fst_seg_hsp_ID)

				{

					HSPs.erase(vec_It);

					vec_It = HSPs.begin(); //vec_It points to the next one after erasing

				}

			}



		}

	}

	else //no exon at all, we skip this gene

	{

		return false;

	}



	if (last_segment_end < last_hsp_end) //need to fix HSPs

	{

		//fix HSPs

		if (start_site > 0)

		{

			vec_It = HSPs.begin();

			while (vec_It != HSPs.end() && (*vec_It)->ID != last_seg_hsp_ID)

			{

				HSPs.erase(vec_It);

				vec_It = HSPs.begin();

			}

		}

		else

		{

			vec_It = HSPs.end(); vec_It--;

			while (vec_It != HSPs.begin() && (*vec_It)->ID != last_seg_hsp_ID)

			{

				vec_next_It = vec_It; vec_next_It--;

				HSPs.erase(vec_It);

				vec_It = vec_next_It;

			}

		}

	}



	//MODIFIED: now update hspID_to_hspIndex_map, since HSPs may have een erased as above!

	cur_hsp_index = 0;

	for (vec_It = HSPs.begin(); vec_It != HSPs.end(); vec_It++)

		hspID_to_hspIndex_map.insert(map<int, int>::value_type((*vec_It)->ID, cur_hsp_index++));



	return true;

}



int ACCP_DONR_graph::LoadHSPAlignScores(HSP_Gene_Pair* hsp_it, vector<int>& tgt_gap_starts, vector<int>& tgt_gap_ends)
{
	int hsp_ID = hsp_it->ID;

	//find alignment
	map<int, Input_Alignment>::iterator alignIt = dataMgr.input_alignments.find(hsp_ID);
	string& matchStr = (*alignIt).second.match_align;
	string& targetStr = (*alignIt).second.target_align;
	string& queryStr = (*alignIt).second.query_align;

	//record gapped regions, for later convert relative index (in Tree) to actual HSP position
	GetGaps(targetStr, tgt_gap_starts, tgt_gap_ends);

	int len = targetStr.length(); //this is MaxItem for C4.5! (may not be the number of a.a. due to possible gaps)
	Item = new Description[len];
	HSPAlignScores = new float[len];

	bool open_gap = false;
	int score;
	int i;
	for (i=0; i<len; i++)
	{
		Item[i] = new AttValue[1];		

		score = similarity_score(targetStr[i], queryStr[i], open_gap, open_gap);

		//first shirt score itself if necessary, then determine class using TREE_CLS_THRESHOLD
		score += TREE_DATA_SHIFT;

		if (score > TREE_CLS_THRESHOLD) //if (score > 0)
			Class(Item[i]) = 1; //positive, exon class
		else
			Class(Item[i]) = 0; //negative score, intron class

		HSPAlignScores[i] = (float)(abs(score)); //weight
	}

#ifdef DEBUG
	dataMgr.outFile << "weights of each position:" << "\n";
	for (i=0; i<len; i++)
	{
		if (Class(Item[i]) == 1)
			dataMgr.outFile << HSPAlignScores[i] << ";";
		else
			dataMgr.outFile << -HSPAlignScores[i] << ";";
	}
#endif

	return len - 1; //MaxItem is the index of last item
}









void ACCP_DONR_graph::PrintHSPs(vector<HSP_Gene_Pair*>& HSPs)

{

	vector<HSP_Gene_Pair*>::iterator hsp_it = HSPs.begin();

	for (; hsp_it != HSPs.end(); hsp_it++)

		dataMgr.outFile << "address:" << *hsp_it << ":" << *(*hsp_it) << "\n";

}

void ACCP_DONR_graph::PrintHSPs(vector<HSP_Gene_Pair>& HSPs)



{



	vector<HSP_Gene_Pair>::iterator hsp_it = HSPs.begin();



	for (; hsp_it != HSPs.end(); hsp_it++)



		dataMgr.outFile << (*hsp_it) << "\n";



}



//this is basically GenePredict4, except first getting potential exons with original HSPs, then try to repair HSPs and possibly 
//modify exons
void ACCP_DONR_graph::GenePredict6()
{
	duration_gblastg_without_comp_final_align += (double)(clock() - gblast_start_time) / CLOCKS_PER_SEC;

	dataMgr.outFile << "//**************START***************//" << "\n" ;

	if (GENBLASTG_NEED_PID)
		max_gblastg_final_align_pid	= 0;

#ifdef WORMBASE
	max_wormbase_final_align_pid = 0;
#endif

#ifdef GENEWISE
	max_genewise_final_align_pid = 0;
#endif

	num_of_correct_exons = 0;
	int max_pred_exons = 0;

#ifdef EVAL_ON_TRUE_EXONS
	bool need_test = true;
	map<string, pair<string, set< pair<int, int> > > >::iterator wb_exon_map_it = WormBase_Exons.find(dataMgr.query_gene);
	if (wb_exon_map_it == WormBase_Exons.end())
	{
		dataMgr.outFile << "NOT found in WormBase EXONS!" << "\n";
		//cout << "NOT found in WormBase genes!" << "\n";
		need_test = false;
	}
	else
	{
		Total_WormBase_Exons += (*wb_exon_map_it).second.second.size();
	}
#endif

	clock_t cur_start_time = clock();

	//get dataMgr.input_alignments' max HSP_ID, used to insert newHSPs (after repairing BLAST HSPs)
	int max_alignment_HSP_ID = 0;
	map<int, Input_Alignment>::iterator alignment_it = dataMgr.input_alignments.begin();
	for (; alignment_it != dataMgr.input_alignments.end(); alignment_it++)
		if (max_alignment_HSP_ID < (*alignment_it).first)
			max_alignment_HSP_ID = (*alignment_it).first;

	max_alignment_HSP_ID++;

	int last_dir_pos = dataMgr.chrSeqFile.find_last_of("\\/");
	if (last_dir_pos == string::npos)
		last_dir_pos = -1;

	string seqid="";

	multimap<Group_Info, vector<HSP_Gene_Pair*> >::iterator groupMapIt = dataMgr.groups.begin();
	
	int rank = 1;
	float dist = (*groupMapIt).first.score;
	int count=1; //for GFF file only, to distinguish genes of the same rank
	
	for (; groupMapIt != dataMgr.groups.end(); groupMapIt++)//, rank++) //for each gene (group)
	{
		if ((*groupMapIt).first.score < dist)
		{
			rank++;
			dist = (*groupMapIt).first.score;
		}

		//alignment pid, use to keep the max pid of current spliced alignment
#ifdef WORMBASE
		int wormbase_align_start, wormbase_align_end;
		wormbase_align_start = (*groupMapIt).first.HSP_start;
		wormbase_align_end = (*groupMapIt).first.HSP_end;
#endif
		//get corresponding chromosome sequence
		string cur_chr = dataMgr.HSP_chr[(*groupMapIt).first.chr_index];
#ifdef DEBUG
		dataMgr.outFile << "cur_chr:" << cur_chr <<"\n";
#endif

		map<int, pair<int, int> >::iterator group_chr_it = dataMgr.group_start_seqno.find(count);
		if (group_chr_it == dataMgr.group_start_seqno.end())
		{
			cout << "cannot find region for this group? " << count << "\n";
			exit(-1);
		}
		vector<string>& chr_seq = dataMgr.group_dna_regions[(*group_chr_it).second.second];
		chromosome_start_pos = (*group_chr_it).second.first;

		//output some info to gff file
		dataMgr.cur_chr_name = seqid + cur_chr;

		if (OUTPUT_GFF)
		{
			dataMgr.gff_str.str(""); //reset gff_str
			dataMgr.gff_gene_str.str("");
			dataMgr.gff_geneinfo_str.str("");
			dataMgr.gff_os << "##sequence-region\t" << dataMgr.cur_chr_name << "_group" << count << "\t1\t" 
				//<< chr_seq.first 
				<< LenOfStrVec(chr_seq) 
				//<< LenOfStrVec((*chrMapIt).second) 
				<< "\n";
		}
		final_start_site = INT_MAX;
		final_end_site = 0;
		//string query_seq;
		dataMgr.GetQuerySeq(query_seq);

		set< pair<int, int> > exon_ref; //this thing is also used in function OutputSites(), so need to be outside #ifdef
#ifdef EVAL_ON_TRUE_EXONS
		if (need_test)
		{
			string gb_chr = cur_chr;
			if ((*groupMapIt).first.isPosStrand)
				gb_chr += '+';
			else
				gb_chr += '-';
			if (gb_chr.compare((*wb_exon_map_it).second.first) == 0)
				exon_ref = (*wb_exon_map_it).second.second;
			else
				need_test = false;
		}
#endif

#ifdef GENBLASTG
		//try to repair missing HSPs (HSPs that Blast did not return):
		vector<HSP_Gene_Pair*> newHSP_ptrs; //collect the new HSP pointers, used for later deleting from memory

		PreProcHSPs((*groupMapIt).first.isPosStrand, (*groupMapIt).second, newHSP_ptrs, max_alignment_HSP_ID);

#ifdef TIMING
		dataMgr.outFile << "after PreProcHSPs: at " << (double)(clock() - gblast_start_time) / CLOCKS_PER_SEC << " seconds" << "\n";
#endif


SELECT_DONOR_ACCEPTOR6:

		//MODIFIED: clear dups only once, then after it's got HSPs, hold it, in case we have to come back 
		//and do SELECT_DONOR_ACCEPTOR again due to first exon having stop (also need to remove the first HSP in HSPs_dup then)
		if (SPLICE_SEGMENT_VERSION == 1)
		{
			HSPs_dup.clear();
			dataMgr.input_alignments_HSPs_dup.clear();
		}

		//LoadData to get all candidate acceptors and donors
		//CalStartEndPos2() are now inside LoadData4()
		if (!LoadData4_forGenePredict6((*groupMapIt).second, chr_seq, false, false))
		{
			for (vector<HSP_Gene_Pair*>::iterator newhsp_it = newHSP_ptrs.begin(); newhsp_it != newHSP_ptrs.end(); newhsp_it++)
				delete *newhsp_it;
			Reset();
			count++;
			continue; //something is wrong, just go on to do the next, for now!
		}

#ifdef TIMING
		dataMgr.outFile << "after LoadData4_forGenePredict6: at " << (double)(clock() - gblast_start_time) / CLOCKS_PER_SEC << " seconds" << "\n";
#endif


		if (VERBOSE)
			cout << "acceptor size: " << acceptors.size() << "; donor size: " << donors.size() << "\n";

		//print header
		dataMgr.outFile << "\n" << dataMgr.query_gene << "|" << dataMgr.HSP_chr[(*groupMapIt).first.chr_index] << ":"; 

		if ((*groupMapIt).first.isPosStrand)
			dataMgr.outFile << start_site << ".." << end_site << "|+|gene cover:";
		else
			dataMgr.outFile << -end_site << ".." << -start_site << "|-|gene cover:";
		PrintHeader(rank, count, 1, groupMapIt);


#ifdef DEBUG_VERSION
		dataMgr.outFile << "input alignments collected: " << dataMgr.input_alignments.size() << "\n";
		cout << "input alignments collected: " << dataMgr.input_alignments.size() << "\n";
#endif

		vector< vector< ExonSiteInfo > > all_alternative_acceptors, all_alternative_donors;

		bool has_possible_exon;
		hspID_added_for_repair = -1; //initialize to -1 (first time compute exon, before repair)

#ifdef COMPUT_EXON_FULL_STEP_BACK

		if (!ComputeExons((*groupMapIt).second, chr_seq, 0, 0, all_alternative_acceptors, all_alternative_donors, has_possible_exon, false))
			if (has_possible_exon)
				goto SELECT_DONOR_ACCEPTOR6;

#else

		while (!ComputeExons((*groupMapIt).second, chr_seq, 0, 0, all_alternative_acceptors, all_alternative_donors, has_possible_exon, false))
			if (!has_possible_exon)
				break;
#endif

#ifdef TIMING
		dataMgr.outFile << "after 1st-round ComputeExons: at " << (double)(clock() - gblast_start_time) / CLOCKS_PER_SEC << " seconds" << "\n";
#endif


		float align_max_pid=0; //no need?
		int align_start, align_end; //not used?
		int i, j;
		int exon_count = 1, alt_count = 1;
		vector<ExonSiteInfo> temp_sites;
		vector<HSP_Gene_Pair> temp_HSPs; //record the HSPs (may be new and old) that corresponds to the previous exon

		map<PairExonSiteInfo, InfoAfterRepair> ExonRepairInfo;

		for (i=0; i<all_alternative_acceptors.size(); i++)
		{
			if (OUTPUT_cDNA)
				dataMgr.cDNA_os << ">" << dataMgr.cur_gene_id << "-A" << alt_count << "\n";

			temp_sites.clear();

			if (REPAIR_HSP_AFTER_EXON == 0) //if -norepair option, simply copy sites to temp_sites for output later
			{
				for (j=0; j<all_alternative_acceptors[i].size(); j++)
				{
					temp_sites.push_back(all_alternative_acceptors[i][j]);
					temp_sites.push_back(all_alternative_donors[i][j]);
				}
			}
			else //for the repair option (default option)
			{
				bool need_last_exon = true;

				for (j=0; j<all_alternative_acceptors[i].size(); j++)
				{
					dataMgr.outFile << "EXON before HSP repair\t" << exon_count << "\t" << all_alternative_acceptors[i][j].first 
							<< "-" << all_alternative_donors[i][j].first << "\n";
					if (j==0 ) //first exon
					{
						ExonSiteInfo p_exon_start, p_exon_end;
						PairExonSiteInfo cur_info(p_exon_start, p_exon_end, all_alternative_acceptors[i][j], all_alternative_donors[i][j]);
						map<PairExonSiteInfo, InfoAfterRepair>::iterator exoninfo_it;

						if ((exoninfo_it = ExonRepairInfo.find(cur_info)) != ExonRepairInfo.end()) {
							temp_sites = (*exoninfo_it).second.temp_sites;
							temp_HSPs = (*exoninfo_it).second.temp_HSPs;
						} else {
							vector<ExonSiteInfo> additional_temp_sites;
							RepairHeadTailExon(true, all_alternative_acceptors[i][j].first, all_alternative_acceptors[i][j].second, 
								all_alternative_donors[i][j].first, all_alternative_donors[i][j].second, 
								chr_seq, max_alignment_HSP_ID, 
								newHSP_ptrs, (*groupMapIt).second, all_alternative_donors[i][j].hsp_index, //all_alternative_acceptors[i][j].hsp_index, 
								temp_sites, temp_HSPs, //collect sites (may be altered) in temp_sites
								all_alternative_acceptors[i][j].frame, additional_temp_sites, 
								j == all_alternative_acceptors[i].size()-1 ); 

							//...need to record the last exon (target, query) in temp_sites, for checking against the next exon (for RepairMidExon)
							//also need to record last bunch of HSPs (may contain newHSPs)!
							if (temp_sites.size() < 2 || temp_sites.size() % 2 > 0 ) //something wrong?!
							{
								cout << "error1: temp_sites has " << temp_sites.size() << "\n";
								exit(-1);
							}

							ExonRepairInfo.insert(map<PairExonSiteInfo, InfoAfterRepair>::value_type(cur_info, 
								InfoAfterRepair(false, temp_sites, temp_HSPs)));
						}
					}
					else //mid exons
					{
						PairExonSiteInfo cur_info(temp_sites[temp_sites.size()-2], temp_sites[temp_sites.size()-1], 
							all_alternative_acceptors[i][j], all_alternative_donors[i][j]);
						map<PairExonSiteInfo, InfoAfterRepair>::iterator exoninfo_it;
						if ((exoninfo_it = ExonRepairInfo.find(cur_info)) != ExonRepairInfo.end())
						{
							if ((*exoninfo_it).second.need_pop)
							{
								temp_sites.pop_back();
								temp_sites.pop_back();
							}

							temp_sites.insert(temp_sites.end(), (*exoninfo_it).second.temp_sites.begin(), (*exoninfo_it).second.temp_sites.end());
							temp_HSPs = (*exoninfo_it).second.temp_HSPs;
						} else {
							vector<ExonSiteInfo> additional_temp_sites;
							bool need_pop = RepairMidExon(temp_sites[temp_sites.size()-2], temp_sites[temp_sites.size()-1], 
								all_alternative_acceptors[i][j], all_alternative_donors[i][j], 
								chr_seq, max_alignment_HSP_ID, newHSP_ptrs, (*groupMapIt).second, 
								temp_sites, temp_HSPs, additional_temp_sites, 
								j == all_alternative_acceptors[i].size()-1 );

							if (temp_sites.size() < 2 || temp_sites.size() % 2 > 0 ) //something wrong?!
							{
								cout << "error2: temp_sites has " << temp_sites.size() << "\n";
								exit(-1);
							}
							ExonRepairInfo.insert(map<PairExonSiteInfo, InfoAfterRepair>::value_type(cur_info, 
							InfoAfterRepair(need_pop, additional_temp_sites, temp_HSPs)));
						}
					}

					exon_count++;
				}

				//last exon
				if (need_last_exon)
				{
					dataMgr.outFile << "last EXON before HSP repair\t" << exon_count << "\t" 
						<< temp_sites[temp_sites.size()-2] 
						<< "-" 
						<< temp_sites[temp_sites.size()-1] 
						<< "\n";

					ExonSiteInfo c_exon_start, c_exon_end;

					PairExonSiteInfo cur_info(temp_sites[temp_sites.size()-2], temp_sites[temp_sites.size()-1], 
						c_exon_start, c_exon_end);

					map<PairExonSiteInfo, InfoAfterRepair>::iterator exoninfo_it;

					if ((exoninfo_it = ExonRepairInfo.find(cur_info)) != ExonRepairInfo.end())
					{
						if ((*exoninfo_it).second.need_pop)
						{
							temp_sites.pop_back();
							temp_sites.pop_back();
						}
						temp_sites.insert(temp_sites.end(), (*exoninfo_it).second.temp_sites.begin(), (*exoninfo_it).second.temp_sites.end());
					}
					else
					{
						vector<ExonSiteInfo> additional_temp_sites;
						bool need_pop = RepairHeadTailExon(false, temp_sites[temp_sites.size()-2].first, temp_sites[temp_sites.size()-2].second, 
							temp_sites[temp_sites.size()-1].first, temp_sites[temp_sites.size()-1].second, 
							chr_seq, max_alignment_HSP_ID, 
							newHSP_ptrs, (*groupMapIt).second, temp_sites[temp_sites.size()-2].hsp_index, 
							temp_sites, temp_HSPs, //collect sites (may be altered) in temp_sites
							temp_sites[temp_sites.size()-2].frame, 
							additional_temp_sites, true); //this function also fills in second_end_temp_sites

						ExonRepairInfo.insert(map<PairExonSiteInfo, InfoAfterRepair>::value_type(cur_info, 
							InfoAfterRepair(need_pop, additional_temp_sites, temp_HSPs)));

					}

				}

#ifdef TIMING
		dataMgr.outFile << "after repairing exons: at " << (double)(clock() - gblast_start_time) / CLOCKS_PER_SEC << " seconds" << "\n";
#endif

}

				//MODIFIED: add checking stop codon
				if (!temp_sites.empty())
				{
					exon_count = 1;
					string final_alignment, final_DNA_string;
					bool found_stop;
					if (!second_end_temp_sites.empty())
					{
						Input_Alignment newAlign;
						string second_end_align, second_end_DNA_string;
						Input_Alignment prev_align;
						align_compare align_score_pid = GetBetterPredGeneSeq(temp_sites, second_end_temp_sites, chr_seq, 
							final_alignment, second_end_align, final_DNA_string, second_end_DNA_string, 
							newAlign, prev_align, found_stop);
						for (int k=0; k<temp_sites.size(); k+=2)
						{
							PrintGFFExons((*groupMapIt).first.isPosStrand, temp_sites[k].first, temp_sites[k+1].first, exon_count, alt_count);
							exon_count++;
						}
						if (OUTPUT_cDNA)
							dataMgr.cDNA_os << final_DNA_string << "\n";

						duration_gblastg_without_comp_final_align += (double)(clock() - cur_start_time) / CLOCKS_PER_SEC;

						dataMgr.outFile << newAlign << "\n";
						//dataMgr.outFile << "PID:" << align_score_pid.align_pid << "\n";
						dataMgr.outFile << "Gene:ID=" << dataMgr.cur_gene_id << "-A" << alt_count << "|"
							<< dataMgr.HSP_chr[(*groupMapIt).first.chr_index] << ":" 
							<< final_start_site << "-" << final_end_site;
						if (start_site > 0)
							dataMgr.outFile << "|+|";
						else
							dataMgr.outFile << "|-|";
						float output_pid = align_score_pid.align_pid * 100;
						dataMgr.outFile << "rank:" << rank << "|score:" << align_score_pid.align_score 
							<< "|PID:" << setprecision(4) << output_pid << "\n";
						float output_coverage = ((float)final_alignment.length() / dataMgr.query_len)*100;
						if (output_coverage > 100)
							output_coverage = 100;
						if (OUTPUT_GFF && GENBLASTG_NEED_PID)
							dataMgr.gff_geneinfo_str << setprecision(4) << output_pid << ";Coverage=" << output_coverage 
								<< ";Note=PID:" << output_pid << "-Cover:" << output_coverage << "\n";
						if (align_max_pid < align_score_pid.align_pid)
						{
							align_max_pid = align_score_pid.align_pid;
							align_start = start_site;
							align_end = temp_sites.back().first;
						}
						second_end_temp_sites.clear(); //reset this
					}
					else
					{
/*
#ifdef EVAL_ON_TRUE_EXONS
				int cur_matching_exons=0;
				int cur_pred_exons=temp_sites.size()/2;
				if (rank == 1)
					if (max_pred_exons < cur_pred_exons)
						max_pred_exons = cur_pred_exons;
#endif
*/						int frame=0;
						string left_chars;
						int prev_donor_site = 0;

						for (int k=0; k<temp_sites.size(); k+=2)
						{
							int site_before_stop;
							string tmp_str;

							if (!HasStopCodon(temp_sites[k].first, temp_sites[k+1].first, 
								chr_seq, frame, site_before_stop, final_alignment, prev_donor_site, left_chars, 
								OUTPUT_cDNA, dataMgr.cDNA_os, false, tmp_str, k==temp_sites.size()-2, found_stop ))
							{
								PrintGFFExons((*groupMapIt).first.isPosStrand, temp_sites[k].first, temp_sites[k+1].first, exon_count, alt_count);
/*
#ifdef EVAL_ON_TRUE_EXONS
	if (rank == 1 && need_test)
	{
						int last_coord_incl_stop;
						if (k < temp_sites.size()-2)
							last_coord_incl_stop = temp_sites[k+1].first;
						else
							last_coord_incl_stop = temp_sites[k+1].first+3;
						if (CheckAgainstTrueExons(exon_ref, temp_sites[k].first, last_coord_incl_stop))
							cur_matching_exons++;
	}
#endif
*/								prev_donor_site = temp_sites[k].first;
								exon_count++;
							}
							else
							{
								PrintGFFExons((*groupMapIt).first.isPosStrand, temp_sites[k].first, site_before_stop, exon_count, alt_count);
/*						
#ifdef EVAL_ON_TRUE_EXONS
	if (rank == 1 && need_test)
	{
						//wormbase exon coordinates include the final stop, so we need to add 3 to our coordinate
						if (CheckAgainstTrueExons(exon_ref, temp_sites[k].first, site_before_stop+3))
							cur_matching_exons++;
	}
#endif
*/								prev_donor_site = site_before_stop;
								break;
							}
						}

						duration_gblastg_without_comp_final_align += (double)(clock() - cur_start_time) / CLOCKS_PER_SEC;

					if (GENBLASTG_NEED_PID)
					{
						Input_Alignment newAlign;
						float align_pid;
						int align_score;
#ifdef PID_BEFORE_SCORE
						align_score = GetGlobalAlignment_PID_CompScore(query_seq, final_alignment, newAlign, align_pid, dataMgr.outFile);
#else
						align_score = GetGlobalAlignment_scorematrix(query_seq, final_alignment, newAlign, align_pid, dataMgr.outFile);
#endif
						dataMgr.outFile << newAlign << "\n";
						//dataMgr.outFile << "PID:" << align_pid << "; score:" << align_score << "\n";
						dataMgr.outFile << "Gene:ID=" << dataMgr.cur_gene_id << "-A" << alt_count << "|" 
							<< dataMgr.HSP_chr[(*groupMapIt).first.chr_index] << ":" 
							<< final_start_site << "-" << final_end_site;
						if (start_site > 0)
							dataMgr.outFile << "|+|";
						else
							dataMgr.outFile << "|-|";
						float output_pid = align_pid * 100;
						dataMgr.outFile << "rank:" << rank << "|score:" << align_score 
							<< "|PID:" << setprecision(4) << output_pid << "\n";
						float output_coverage = ((float)final_alignment.length() / dataMgr.query_len)*100;
						if (output_coverage > 100)
							output_coverage = 100;
						if (OUTPUT_GFF)
							dataMgr.gff_geneinfo_str << setprecision(4) << output_pid << ";Coverage=" << output_coverage 
								<< ";Note=PID:" << output_pid << "-Cover:" << output_coverage << "\n";
						if (align_max_pid < align_pid)
						{
							align_max_pid = align_pid;
							align_start = start_site;
							align_end = prev_donor_site;
						}
					}

					}		
					

					if (OUTPUT_Protein)
					{
						dataMgr.pred_protein_os << ">" << dataMgr.cur_gene_id << "-A" << alt_count << "\n";
						dataMgr.pred_protein_os << final_alignment;
						if (found_stop)
							dataMgr.pred_protein_os << '*';
						dataMgr.pred_protein_os << "\n";
					}

					if (OUTPUT_GFF) //final_start_site and final_end_site have been updated in PrintGFFExons()
						dataMgr.gff_os << dataMgr.gff_gene_str.str() << final_start_site << "\t" << final_end_site << dataMgr.gff_geneinfo_str.str() << dataMgr.gff_str.str();

					cur_start_time = clock();
				}

				exon_count = 1; //reset
				if (i+1 < all_alternative_acceptors.size())
				{
					//modified: at most certain number of alternatives (to reduce unnecessary alternatives that don't help much)?!
					if (alt_count >= 50)
						break;
					dataMgr.outFile << "ALTERNATIVE:" << "\n";
					alt_count++;
					if (OUTPUT_GFF)
					{
						dataMgr.gff_str.str(""); //reset
						dataMgr.gff_gene_str.str("");
						dataMgr.gff_geneinfo_str.str("");
						//if ((*groupMapIt).first.isPosStrand)
						//	dataMgr.gff_os << dataMgr.cur_chr_name  << "\tgenBlastG" << USER_ID << "\tgene\t";
						//<< start_site << "\t" << end_site 
						//		<< "\t" << (*groupMapIt).first.score << "\t+\t.\tID=";
						//else
						//	dataMgr.gff_os << dataMgr.cur_chr_name  << "\tgenBlastG" << USER_ID << "\tgene\t";
						//<< -end_site << "\t" << -start_site 
						//		<< "\t" << (*groupMapIt).first.score << "\t-\t.\tID=";						
					}
					final_start_site = INT_MAX;
					final_end_site = 0;
					
					dataMgr.outFile << "\n" << dataMgr.query_gene << "|" << dataMgr.HSP_chr[(*groupMapIt).first.chr_index] << ":"; 
					if ((*groupMapIt).first.isPosStrand)
						dataMgr.outFile << start_site << ".." << end_site << "|+|gene cover:";
					else
						dataMgr.outFile << -end_site << ".." << -start_site << "|-|gene cover:";
					PrintHeader(rank, count, alt_count, groupMapIt);
				}
			}
			
		//clean up memory
CLEAN_UP6:
		for (vector<HSP_Gene_Pair*>::iterator newhsp_it = newHSP_ptrs.begin(); newhsp_it != newHSP_ptrs.end(); newhsp_it++)
			delete *newhsp_it;
#endif //endif for #GENBLASTG

	duration_gblastg_without_comp_final_align += (double)(clock() - cur_start_time) / CLOCKS_PER_SEC;



		//compute alignment and pid from wormbase and genewise

	if (rank == 1)

	{

		if (GENBLASTG_NEED_PID)
		{
		//of_perform << "\tGBLASTG\t" << align_max_pid;
		if (align_max_pid > max_gblastg_final_align_pid)
			max_gblastg_final_align_pid = align_max_pid;
		//GBLASTG_FOUND = true;
		}

	}//end if (rank == 1)



		cur_start_time = clock();



		//reset graph for next construction

		Reset();

		dataMgr.outFile << "//---------------------------------//\n";

		cout << "current gene (rank " << rank << ") done\n";
		count++;
	//} //one chromosome done

	} //all chromosomes done

#ifdef EVAL_ON_TRUE_EXONS
	Total_Matching_Exons += num_of_correct_exons;
	Total_GenBlast_Exons += max_pred_exons;
#endif

	dataMgr.outFile << "\n" << "//***************END****************//" << "\n" << "\n";

}





//check for in-frame stop codon, if there is, then cut HSP into pieces that have no stops

void ACCP_DONR_graph::FixHSPsWithInFrameStop(bool isPosStrand, vector<HSP_Gene_Pair*>& newHSPs, 

											 vector<HSP_Gene_Pair*>& newHSP_ptrs, int& max_alignment_HSP_ID)

{
	hspID_due_to_stop.clear();

	int size_ori = newHSPs.size();



	map<int, vector<HSP_Gene_Pair*> > map_newHSPs;



	int i, t_pos_start, t_pos_end;

	for (i = 0; i < size_ori; i++)

	{

		Input_Alignment& align = (*(dataMgr.input_alignments.find(newHSPs[i]->ID))).second;

		string& qStr = align.query_align;

		string& mStr = align.match_align;

		string& tStr = align.target_align;



#ifdef USE_LAST_HSP_WHEN_HAS_STOP

		t_pos_start = tStr.find_last_of('*');

		if (t_pos_start != string::npos)

		{

			const_cast<string&>(tStr).erase(0, t_pos_start+1);

			const_cast<string&>(qStr).erase(0, t_pos_start+1);

			const_cast<string&>(mStr).erase(0, t_pos_start+1);

			newHSPs[i]->gene_start = newHSPs[i]->gene_end - (qStr.length() - CountCharInStr(qStr, '-')) + 1;

			if (isPosStrand)

				newHSPs[i]->HSP_start = newHSPs[i]->HSP_end - 3*(tStr.length() - CountCharInStr(tStr, '-')) + 1;

			else

				newHSPs[i]->HSP_end = newHSPs[i]->HSP_start + 3*(tStr.length() - CountCharInStr(tStr, '-')) - 1;

		}



#else

		//now check for '*' in targetStr

		string hold_qStr(qStr);

		string hold_tStr(tStr);

		string hold_mStr(mStr);

		int hold_q_s, hold_q_e, hold_h_s, hold_h_e;

		hold_q_e = newHSPs[i]->gene_end;

		hold_h_e = newHSPs[i]->HSP_end;

		bool first_portion_is_valid = true;

		if ((t_pos_start = hold_tStr.find('*')) != string::npos) //find first stop, modify existing HSP

		{
			//UPDATE: also remove the part at either end of HSP that has target or query aligned with gaps
			int t_pos_ori = t_pos_start; //hold this position
			while (t_pos_start - 1 >= 0)
				if (qStr[t_pos_start-1] == '-' || tStr[t_pos_start-1] == '-')
					t_pos_start--;
				else
					break;

			const_cast<string&>(qStr).erase(t_pos_start);

			const_cast<string&>(tStr).erase(t_pos_start);

			const_cast<string&>(mStr).erase(t_pos_start);

			newHSPs[i]->gene_end = newHSPs[i]->gene_start + t_pos_start - CountCharInStr(qStr, '-') - 1;

			if (isPosStrand)

				newHSPs[i]->HSP_end = newHSPs[i]->HSP_start + (t_pos_start - CountCharInStr(tStr, '-'))*3 - 1;

			else

				newHSPs[i]->HSP_start = newHSPs[i]->HSP_end - (t_pos_start - CountCharInStr(tStr, '-'))*3 + 1;



			if (mStr.length() == CountCharInStr(mStr, ' ')+CountCharInStr(mStr, '+')) //no perfect match, need to remove this HSP!

				first_portion_is_valid = false;



#ifdef DEBUG

				dataMgr.outFile << *(newHSPs[i]) << "(first portion is valid?" << first_portion_is_valid << ")" << "\n";

				dataMgr.outFile << align << "\n";

#endif

			t_pos_end = hold_tStr.find_first_not_of("*-", t_pos_ori); //forward to next non-stop position (also non-gap)

			if (t_pos_end != string::npos)
			{
				while (t_pos_end < hold_tStr.length()) //get to the next "good" base pair
					if (hold_tStr[t_pos_end] == '*' || hold_tStr[t_pos_end] == '-' || hold_qStr[t_pos_end] == '-')
						t_pos_end++;
					else
						break;

				string tmp_qStr = hold_qStr.substr(t_pos_start, t_pos_end - t_pos_start);

				string tmp_tStr = hold_tStr.substr(t_pos_start, t_pos_end - t_pos_start);

				hold_q_s = newHSPs[i]->gene_end + 1 + t_pos_end - t_pos_start - CountCharInStr(tmp_qStr, '-');

				if (isPosStrand)

					hold_h_s = newHSPs[i]->HSP_end + 1 + (t_pos_end - t_pos_start - CountCharInStr(tmp_tStr, '-'))*3; //forward, skip the position for stop codon

				else

					hold_h_e = newHSPs[i]->HSP_start - 1 - (t_pos_end - t_pos_start - CountCharInStr(tmp_tStr, '-'))*3;



#ifdef DEBUG

					dataMgr.outFile << "hold_q_s:" << hold_q_s << "\n";

					dataMgr.outFile << "hold_h_s:" << hold_h_s << "; hold_h_e:" << hold_h_e << "\n";

#endif



				t_pos_start = t_pos_end;



				vector<HSP_Gene_Pair*> tmp_newHSPs; //collect possible new HSPs here



				t_pos_end = hold_tStr.find('*', t_pos_start);
				if (t_pos_end == string::npos)
					t_pos_end = hold_tStr.length();

				//while ((t_pos_end = hold_tStr.find('*', t_pos_start)) != string::npos)
				while (true)
				{
					t_pos_ori = t_pos_end; //hold value
					while (t_pos_end-1 >= 0) //t_pos_end get to the earliest "bad" base pair
						if (hold_tStr[t_pos_end-1] == '*' || hold_tStr[t_pos_end-1] == '-' || hold_qStr[t_pos_end-1] == '-')
							t_pos_end--;
						else
							break;

					Input_Alignment newAlign;

					newAlign.match_align = hold_mStr.substr(t_pos_start, t_pos_end - t_pos_start);

					newAlign.target_align = hold_tStr.substr(t_pos_start, t_pos_end - t_pos_start);

					newAlign.query_align = hold_qStr.substr(t_pos_start, t_pos_end - t_pos_start);



					int cur_g_end, cur_h_start, cur_h_end; //use these to hold values

					cur_g_end = hold_q_s + t_pos_end - t_pos_start - CountCharInStr(newAlign.query_align, '-') - 1;



					HSP_Gene_Pair* newHSP = new HSP_Gene_Pair;

					newHSP->gene_start = hold_q_s;

					newHSP->gene_end = cur_g_end; //hold_q_s + t_pos_end - t_pos_start - CountCharInStr(newAlign.query_align, '-') - 1;

					if (isPosStrand)

					{

						newHSP->HSP_start = hold_h_s;

						newHSP->HSP_end = hold_h_s + (t_pos_end - t_pos_start - CountCharInStr(newAlign.target_align, '-'))*3 - 1;

					}

					else

					{

						newHSP->HSP_end = hold_h_e;

						newHSP->HSP_start = hold_h_e - (t_pos_end - t_pos_start - CountCharInStr(newAlign.target_align, '-'))* 3 + 1;

					}

					newHSP->pid = ((float)(t_pos_end - t_pos_start - CountCharInStr(newAlign.match_align, ' ') - CountCharInStr(newAlign.match_align, '+'))/

						(t_pos_end - t_pos_start))*100;



					cur_h_start = newHSP->HSP_start;

					cur_h_end = newHSP->HSP_end;



					if (newAlign.match_align.length() > (CountCharInStr(newAlign.match_align, ' ')+CountCharInStr(newAlign.match_align, '+')))

					{ //current portion is valid



					if (first_portion_is_valid)

					{

					newHSP->ID = max_alignment_HSP_ID;



					tmp_newHSPs.push_back(newHSP);

					newHSP_ptrs.push_back(newHSP); //only used to delete the memory later on

					dataMgr.input_alignments.insert(map<int, Input_Alignment>::value_type(max_alignment_HSP_ID, newAlign));

					max_alignment_HSP_ID++;



					//cur_g_end = newHSP->gene_end;

					}

					else //if not (first_portion_is_valid)

					{ //replace the first portion

						const_cast<string&>(tStr).assign(newAlign.target_align);

						const_cast<string&>(qStr).assign(newAlign.query_align);

						const_cast<string&>(mStr).assign(newAlign.match_align);



						int t_ID = newHSPs[i]->ID;

						*(newHSPs[i]) = *newHSP; //copy HSP

						newHSPs[i]->ID = t_ID;



						//cur_g_end = newHSPs[i]->gene_end;



						delete newHSP;

						first_portion_is_valid = true;

					}

					} //end if (current portion is valid)



					if (t_pos_ori == hold_tStr.length()) //this is already the last HSP, get out

						break;



					t_pos_start = hold_tStr.find_first_not_of("*-", t_pos_ori); //forward to next non-stop position
					if (t_pos_start == string::npos) //end with stop codon, get out
						break;
					while (t_pos_start < hold_tStr.length())
						if (hold_tStr[t_pos_start] == '*' || hold_tStr[t_pos_start] == '-' || hold_qStr[t_pos_start] == '-')
							t_pos_start++;
						else
							break;
					if (t_pos_start == hold_tStr.length()) //end, no more good HSPs, get out
						break;

					tmp_qStr = hold_qStr.substr(t_pos_end, t_pos_start - t_pos_end);

					tmp_tStr = hold_tStr.substr(t_pos_end, t_pos_start - t_pos_end);

					hold_q_s = cur_g_end + 1 + t_pos_start - t_pos_end - CountCharInStr(tmp_qStr, '-');



					if (isPosStrand)

						hold_h_s = cur_h_end + 1 + (t_pos_start - t_pos_end - CountCharInStr(tmp_tStr, '-'))*3; //forward, skip the position for stop codon

					else

						hold_h_e = cur_h_start - 1 - (t_pos_start - t_pos_end - CountCharInStr(tmp_tStr, '-'))*3;



#ifdef DEBUG

					dataMgr.outFile << "tmp_qStr:" << tmp_qStr << "\n";

					dataMgr.outFile << "tmp_tStr:" << tmp_tStr << "\n";

					dataMgr.outFile << "hold_q_s:" << hold_q_s << "\n";

					dataMgr.outFile << "hold_h_s:" << hold_h_s << "; hold_h_e:" << hold_h_e << "\n";

#endif



					t_pos_end = hold_tStr.find('*', t_pos_start); //go on to next stop

					if (t_pos_end == string::npos) //this signals the last HSP

						t_pos_end = hold_tStr.length();

				} //end while (true)



				if (!tmp_newHSPs.empty())

					map_newHSPs.insert(map<int, vector<HSP_Gene_Pair*> >::value_type(i, tmp_newHSPs));

			}



			

				

		}



#endif //#ifdef...else... USE_LAST_HSP_WHEN_HAS_STOP

	}



#ifndef USE_LAST_HSP_WHEN_HAS_STOP

	//put the newHSPs into the HSP vector, in correct order

	int j=0;

	map<int, vector<HSP_Gene_Pair*> >::iterator it = map_newHSPs.begin();

	for (; it != map_newHSPs.end(); it++)

	{

		i = (*it).first + j;	

		if (isPosStrand)

		{

			vector<HSP_Gene_Pair*>::reverse_iterator vec_it = (*it).second.rbegin();

			for (; vec_it != (*it).second.rend(); vec_it++)
			{
				newHSPs.insert(newHSPs.begin()+i, *vec_it);			
				hspID_due_to_stop.insert(newHSPs[i]->ID); //this is the id of the newly-inserted hsp
				i++;
			}			

		}

		else

		{

			vector<HSP_Gene_Pair*>::iterator vec_it = (*it).second.begin();

			for (; vec_it != (*it).second.end(); vec_it++)

			{

				newHSPs.insert(newHSPs.begin()+i+1, *vec_it);
				hspID_due_to_stop.insert(newHSPs[i+1]->ID); //this is the id of the newly-inserted hsp
				i++;

			}

		}

		j += (*it).second.size(); //keep count of how many HSPs have been inserted

	}

#endif


#ifdef DEBUG
	dataMgr.outFile << "after fixing internal stops in HSPs" << "\n";
	for (vector<HSP_Gene_Pair*>::iterator hsps_it = newHSPs.begin(); hsps_it != newHSPs.end(); hsps_it++)
	{
		dataMgr.outFile << "HSP[" << (*hsps_it)->ID << "]:" << (*hsps_it)->HSP_start << "-" << (*hsps_it)->HSP_end 
			<< "(" << (*hsps_it)->gene_start << "-" << (*hsps_it)->gene_end << ")" << "\n";

		dataMgr.outFile << "alignment:" << "\n";
		map<int, Input_Alignment>::iterator map_it = dataMgr.input_alignments.find((*hsps_it)->ID);
		if (map_it == dataMgr.input_alignments.end())
		{
			cout << "error3: cannot find the alignment for this HSP" << "\n";
			exit(-1);
		}
		dataMgr.outFile << (*map_it).second << "\n";
	}
#endif


}



//for either start or end exon, if the gene_start or gene_end is missing some query coverage, 

//try find an additional HSP and may modify exon:

//after finding the HSP, do the LoadData() thing, including finding the new start_site/end_site, finding splice segments from 

//the newHSP and previous start/end HSP, finding acceptors/donors, do genepredict() thing, test alignment with query,

//if the new alignment is better, use the new one, otherwise keep the old exon

//if there's next exon to be examined, 

//temp_sites record the "new" exons resulted from current repair (if it's indeed new, first 

//pop out the last exon in temp_sites, since it's replaced by the new exons);

//temp_HSPs record the "new" HSPs that correspond to the last exon resulted from current repair

bool ACCP_DONR_graph::RepairHeadTailExon(bool isHeadExon, int target_site, int query_site, int target_end_site, int query_end_site, 
										 //pair<int, char*>& chr_seq, 
										 vector<string>& chr_seq, 
										 int& max_alignment_HSP_ID, vector<HSP_Gene_Pair*>& newHSP_ptrs, 
										 vector<HSP_Gene_Pair*>& blastHSPs, int blast_hsp_index, //blast_hsp_index is the index (position) of the HSP which is used to get the first/last exon
										 vector<ExonSiteInfo>& temp_sites, //tmp_sites keep all unconfirmed sites?
										 vector<HSP_Gene_Pair>& temp_HSPs, //vector<HSP_Gene_Pair*>& temp_HSPs, 
										 int exon_start_frame, 
										 vector<ExonSiteInfo>& additional_temp_sites, 
										 bool is_last_exon)
{
#ifdef DEBUG
					dataMgr.outFile << "exon_start_frame:" << exon_start_frame;
					dataMgr.outFile << "; started with temp_sites: " << "\n";
					int j;
					for (j=0; j<temp_sites.size(); j++)
						dataMgr.outFile << temp_sites[j] << "\n";
					dataMgr.outFile << "started with temp_HSPs (for prev exon): " << "\n";
					for (j=0; j<temp_HSPs.size(); j++)
						dataMgr.outFile << (temp_HSPs[j]) << "\n";
#endif
	int exon_end_frame = (target_end_site - target_site + 1 - exon_start_frame) % 3;

	//string query_seq, target_seq, q_seq;
	//dataMgr.GetQuerySeq(query_seq);

	string target_seq, q_seq;



	vector<HSP_Gene_Pair> newHSPs; //used to store everything

	vector<HSP_Gene_Pair*> ptr_to_newHSPs;



	//align_pos[frame][]: all relative indexes!

	//1: start pos on query_seq; 2: end pos on query_seq; 3: start pos on target; 4: end pos on target

	int align_pos[3][4]; 

	int tgt_str_len;



		if (isHeadExon && query_site > END_EXON_LEN )

		{

#ifdef DEBUG

			dataMgr.outFile << "Repairing head exon " << target_site << "-" << target_end_site << " (query " << query_site 

				<<"-" << query_end_site << ")" << "\n";

#endif

				double max_score=0, cur_score; //MODIFIED: now use the alignment that has highest PID

				float max_pid = 0, cur_pid;


				if (target_site < 0) //neg strand

					tgt_str_len = GetSubstrFromVecStrs_NegRev_ForRepair(chr_seq, false, -target_site, EXPLORE_END_EXON_LEN, target_seq);

				else

					tgt_str_len = GetSubstrFromVecStrs_ForRepair(chr_seq, true, target_site-EXPLORE_END_EXON_LEN-1, EXPLORE_END_EXON_LEN, target_seq);

				StrToLower(target_seq);



				HSP_Gene_Pair newHSP[3];

				HSP_Gene_Pair* curHSP = new HSP_Gene_Pair;

				Input_Alignment newAlign[3];

				int j=0;

				string result_seq[3];

				q_seq = query_seq.substr(0, query_site - 1);

				for (int i=0; i<3; i++)

				{

					DNA2AA(target_seq, i, result_seq[i]);

#ifdef DEBUG

					dataMgr.outFile << "target_seq:" << "\n" << target_seq << "\n";

					dataMgr.outFile << result_seq[i] << "\n";

#endif					

					//if (target_site < 0) //same as positive?

					//	cur_score = GetLocalAlignment(q_seq, result_seq[i], 

					//		1, -(-target_site + EXPLORE_END_EXON_LEN - i), max_alignment_HSP_ID, newHSP[i], newAlign[i], 

					//		false, dataMgr.outFile, align_pos[i]);

					//else

						cur_score = GetLocalAlignment(q_seq, result_seq[i], 

							1, target_site-tgt_str_len+i, max_alignment_HSP_ID, newHSP[i], newAlign[i], 

							false, dataMgr.outFile, align_pos[i]);

					cur_pid = newHSP[i].pid;



					if (cur_score == 0)

						continue;

					

if (EXTEND_HSP_BY_SCORE)

{

					if (max_score < cur_score || 

						(target_site < 0 && max_score == cur_score && newHSP[i].HSP_start < curHSP->HSP_start) ||

						(target_site > 0 && max_score == cur_score && newHSP[i].HSP_end > curHSP->HSP_end) )

					{

						j = i;

						max_score = cur_score;

						max_pid = cur_pid;

						*curHSP = newHSP[i];

					}

}

else

{

					if (max_pid < cur_pid || 

						(target_site < 0 && max_pid == cur_pid && newHSP[i].HSP_start < curHSP->HSP_start) || 

						(target_site > 0 && max_pid == cur_pid && newHSP[i].HSP_end > curHSP->HSP_end) )

					{

						j = i;

						max_score = cur_score;

						max_pid = cur_pid;

						*curHSP = newHSP[i];

					}

}

				}



				//dataMgr.HSP_neg_gene[chr_index].push_back(newHSP[j]);

				//newHSPs.push_back(&(dataMgr.HSP_neg_gene[chr_index].back()));

				if (max_score == 0 || max_score < REPAIR_HSP_MIN_INIT_SCORE) //<= 0)

				{
						UseOldExon(temp_sites, additional_temp_sites, temp_HSPs, 
								 target_site, query_site, target_end_site, query_end_site, 
								 blast_hsp_index, blastHSPs);

#ifdef DEBUG

					dataMgr.outFile << "no extra alignment, use old exon: temp_HSPs for this exon: " << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

						return false; //no good alignment that can be used to adjust exons, just use the old exon

				}

				else

				{

					//if the score of init alignment is big enough, try extend this alignment to both directions, just like BLAST!

					//extend curHSP to both directions if possible, stop extending when the score of alignment is 

					//lower than ("max_score" - REPAIR_HSP_EXTEND_SCORE_DROP), "max_score" is the running max score

					//if (max_score >= REPAIR_HSP_MIN_INIT_SCORE)
					//{

						string tmpStr = "";

						if (align_pos[j][3]+1 < result_seq[j].length())

							tmpStr = result_seq[j].substr(align_pos[j][3]+1);


						string result_begin_str = result_seq[j].substr(0, align_pos[j][2]);

						if (target_site < 0)
						ExtendAlignment(max_score, query_seq, 
							//result_seq[j].substr(0, align_pos[j][2]), //this may also be empty string
							result_begin_str,

							tmpStr, 

							chr_seq, false, curHSP, newAlign[j], dataMgr.outFile, 

							//INT_MAX, -target_end_site);

							INT_MAX, blastHSPs.front()->HSP_start);
						else

						ExtendAlignment(max_score, query_seq, 

							//result_seq[j].substr(0, align_pos[j][2]), //this may also be empty string
							result_begin_str,

							tmpStr, 

							chr_seq, true, curHSP, newAlign[j], dataMgr.outFile, 

							//0, target_end_site);

							0, blastHSPs.back()->HSP_end);

#ifdef DEBUG

						dataMgr.outFile << "extended alignment: " << "\n";

						dataMgr.outFile << *curHSP << "\n";

						dataMgr.outFile << newAlign[j] << "\n";

#endif

					//}



					newHSPs.push_back(*curHSP);

					newHSP_ptrs.push_back(curHSP);

					dataMgr.input_alignments.insert(map<int, Input_Alignment>::value_type(max_alignment_HSP_ID, newAlign[j]));
					
					hspID_added_for_repair = max_alignment_HSP_ID; //set hspID_added_for_repair to the new HSP
#ifdef DEBUG
					dataMgr.outFile << "hspID_added_for_repair is: HSP[" << hspID_added_for_repair << "]\n";
#endif
					max_alignment_HSP_ID++;

					



					//now redo PreProcHSPs and LoadData4 etc.

					if (target_site < 0)

					{

						for (j=0; j<=blast_hsp_index; j++)

						{

							newHSPs.push_back(*(blastHSPs[j]));

						newHSPs.back().ID = max_alignment_HSP_ID; //also replicate the alignment so it won't be changed on the original blast alignment

						dataMgr.input_alignments.insert(map<int, Input_Alignment>::value_type(max_alignment_HSP_ID, 

							(*(dataMgr.input_alignments.find(blastHSPs[j]->ID))).second));

						max_alignment_HSP_ID++;

						}

					}

					else

					{

						for (j=blastHSPs.size()-1; j>=blast_hsp_index; j--)

						{

							newHSPs.insert(newHSPs.begin(), *(blastHSPs[j]));

						newHSPs.front().ID = max_alignment_HSP_ID; //also replicate the alignment so it won't be changed on the original blast alignment

						dataMgr.input_alignments.insert(map<int, Input_Alignment>::value_type(max_alignment_HSP_ID, 

							(*(dataMgr.input_alignments.find(blastHSPs[j]->ID))).second));

						max_alignment_HSP_ID++;

						}

					}

#ifdef DEBUG

					dataMgr.outFile << "newHSPs for repairing head exon: " << "\n";

					for (j=0; j<newHSPs.size(); j++)

						dataMgr.outFile << (newHSPs[j]) << "\n";

#endif

					for (j=0; j<newHSPs.size(); j++)

						ptr_to_newHSPs.push_back(&(newHSPs[j]));

					PreProcHSPs((target_site>0), ptr_to_newHSPs, newHSP_ptrs, max_alignment_HSP_ID);

					//need to update newHSPs? since PreProcHSPs may have inserted new HSPs...nah, just use ptr_to_newHSPs!

/*					int k=0;

					for (j=0; j<ptr_to_newHSPs.size(); j++) 

					{

						if (ptr_to_newHSPs[j]->ID != newHSPs[k].ID)

							newHSPs.insert(newHSPs.begin()+k, *(ptr_to_newHSPs[j]));

						k++;

					}

*/

#ifdef DEBUG

					dataMgr.outFile << "newHSPs after preproc: " << "\n";

					for (j=0; j<ptr_to_newHSPs.size(); j++)

						dataMgr.outFile << *(ptr_to_newHSPs[j]) << "\n";

#endif



					end_site = target_end_site; //fix end_site first

SELECT_DONOR_ACCEPTOR6_HEADEXON:

					if (SPLICE_SEGMENT_VERSION == 1)

					{

						HSPs_dup.clear();

						dataMgr.input_alignments_HSPs_dup.clear();

					}

					if (!LoadData4_forGenePredict6(ptr_to_newHSPs, chr_seq, false, true))
					{
						UseOldExon(temp_sites, additional_temp_sites, temp_HSPs, 
								 target_site, query_site, target_end_site, query_end_site, 
								 blast_hsp_index, blastHSPs);

#ifdef DEBUG

					dataMgr.outFile << "use old exon: temp_HSPs for this exon: " << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

						return false; //something is wrong, just use the old exon

					}



					vector< vector<ExonSiteInfo> > all_alternative_acceptors, all_alternative_donors;

					bool has_possible_exon;
#ifdef COMPUT_EXON_FULL_STEP_BACK

					if (!ComputeExons(ptr_to_newHSPs, chr_seq, 0, exon_end_frame, all_alternative_acceptors, 
						all_alternative_donors, has_possible_exon, true))
						if (has_possible_exon)
							goto SELECT_DONOR_ACCEPTOR6_HEADEXON;

#else

					while (!ComputeExons(ptr_to_newHSPs, chr_seq, 0, exon_end_frame, all_alternative_acceptors, 
						all_alternative_donors, has_possible_exon, true))
						if (!has_possible_exon)
							break;

#endif

#ifdef DEBUG
					dataMgr.outFile << "computed new exon candidates (" << has_possible_exon << ")" << "\n";
#endif
					if (!has_possible_exon)
					{
						UseOldExon(temp_sites, additional_temp_sites, temp_HSPs, 
								 target_site, query_site, target_end_site, query_end_site, 
								 blast_hsp_index, blastHSPs);
						return false;
					}

					//multiply all_alternative_acceptors/donors into the old all_alternative_acceptors/donors?
					int frame = 0;
					string final_alignment, left_chars;
					int prev_donor_site = 0;
					//temp_sites.clear();
					int site_before_stop;

					for (j=0; j<all_alternative_acceptors[0].size(); j++) //use the first one, for now
					{
						string tmp_str;
						bool found_stop;
						if (HasStopCodon(all_alternative_acceptors[0][j].first, all_alternative_donors[0][j].first, 
							chr_seq, frame, site_before_stop, final_alignment, prev_donor_site, left_chars, 
							false, dataMgr.cDNA_os, false, tmp_str, (is_last_exon && j==all_alternative_acceptors[0].size()-1), found_stop ))
						{ //if there's stop in new exons, discard them
							UseOldExon(temp_sites, additional_temp_sites, temp_HSPs, 
								 target_site, query_site, target_end_site, query_end_site, 
								 blast_hsp_index, blastHSPs);

#ifdef DEBUG

					dataMgr.outFile << "new exons have stop, use old exon: temp_HSPs for this exon: " << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

							return false;

						}

						else

						{

							prev_donor_site = all_alternative_donors[0][j].first;

						}

					}

					Input_Alignment newAlign;

					float align_pid;
					string query_begin_str = query_seq.substr(0, query_end_site);
					int align_score;
#ifdef PID_BEFORE_SCORE
					align_score = GetGlobalAlignment_PID_CompScore(query_begin_str, final_alignment, newAlign, align_pid, dataMgr.outFile);
#else
					align_score = GetGlobalAlignment_scorematrix(//query_seq.substr(0, query_end_site), //G++ compile error
						query_begin_str, final_alignment, newAlign, align_pid, dataMgr.outFile);
#endif
					align_compare align_new(align_score, align_pid);

					//dataMgr.outFile << newAlign << "\n";

					//dataMgr.outFile << "PID:" << align_pid << "\n";

					

					frame = 0;

					final_alignment = "";

					left_chars= "";

					prev_donor_site = 0;

					int num_of_new_exons = all_alternative_acceptors[0].size();

#ifdef DEBUG

					dataMgr.outFile << "num_of_new_exons:" << num_of_new_exons << "\n";

#endif

					string tmp_str;
					bool found_stop;
					if (HasStopCodon(target_site, target_end_site, chr_seq, frame, site_before_stop, final_alignment, prev_donor_site, left_chars, 
						false, dataMgr.cDNA_os, false, tmp_str, is_last_exon, found_stop))

					{

						for (j=0; j<num_of_new_exons; j++) //use the first one, for now

						{

							temp_sites.push_back(all_alternative_acceptors[0][j]);

							temp_sites.push_back(all_alternative_donors[0][j]);



							additional_temp_sites.push_back(all_alternative_acceptors[0][j]);

							additional_temp_sites.push_back(all_alternative_donors[0][j]);

							//output some additional info, for inspection purpose only
							if (all_alternative_donors[0][j].first < target_site)
								dataMgr.outFile << "NEW EXON after repair:" << all_alternative_acceptors[0][j].first << "-" << all_alternative_donors[0][j].first << "\n";
						}

						//temp_HSPs record the HSPs used for the last exon among [all_alternative_acceptors, all_alternative_donors]

						temp_HSPs.clear();

						if (target_site < 0)

						{

							for (j=all_alternative_acceptors[0][num_of_new_exons-1].hsp_index; 

								j<=all_alternative_donors[0][num_of_new_exons-1].hsp_index; j++)

									temp_HSPs.push_back(*(ptr_to_newHSPs[j]));

						}

						else

						{

							for (j=all_alternative_donors[0][num_of_new_exons-1].hsp_index;

								j<=all_alternative_acceptors[0][num_of_new_exons-1].hsp_index; j++)

									temp_HSPs.push_back(*(ptr_to_newHSPs[j]));



						}

#ifdef DEBUG

					dataMgr.outFile << "old exon has stop, use new exons, HSPs for this last exon:" << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

					return true;

					}

					else

					{

						float ori_align_pid;
						string query_begin_str = query_seq.substr(0, query_end_site);
						int ori_align_score;
#ifdef PID_BEFORE_SCORE
						ori_align_score = GetGlobalAlignment_PID_CompScore(query_begin_str, final_alignment, newAlign, ori_align_pid, dataMgr.outFile);
#else
						ori_align_score = GetGlobalAlignment_scorematrix(//query_seq.substr(0, query_end_site), 
							query_begin_str, final_alignment, newAlign, ori_align_pid, dataMgr.outFile);
#endif
						align_compare align_old(ori_align_score, ori_align_pid);

#ifdef DEBUG
							dataMgr.outFile << "align_new:" << align_new << ";align_old:" << align_old << "\n";
#endif

						//if (ori_align_score < align_score || (ori_align_score == align_score && ori_align_pid < align_pid))
						if (align_new > align_old)
						{

							for (j=0; j<num_of_new_exons; j++) //use the first one, for now

							{

								temp_sites.push_back(all_alternative_acceptors[0][j]);

								temp_sites.push_back(all_alternative_donors[0][j]);



								additional_temp_sites.push_back(all_alternative_acceptors[0][j]);

								additional_temp_sites.push_back(all_alternative_donors[0][j]);

							//output some additional info, for inspection purpose only
							if (all_alternative_donors[0][j].first < target_site)
								dataMgr.outFile << "NEW EXON after repair:" << all_alternative_acceptors[0][j].first << "-" << all_alternative_donors[0][j].first << "\n";
							}

							temp_HSPs.clear();

							if (target_site < 0)

							{

								for (j=all_alternative_acceptors[0][num_of_new_exons-1].hsp_index; 

									j<=all_alternative_donors[0][num_of_new_exons-1].hsp_index; j++)

										temp_HSPs.push_back(*(ptr_to_newHSPs[j]));

							}

							else

							{

								for (j=all_alternative_donors[0][num_of_new_exons-1].hsp_index;

									j<=all_alternative_acceptors[0][num_of_new_exons-1].hsp_index; j++)

										temp_HSPs.push_back(*(ptr_to_newHSPs[j]));



							}

#ifdef DEBUG

					dataMgr.outFile << "both have no stop, use new exons, HSPs for this last exon:" << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

					return true;

						}

						else
						{
							UseOldExon(temp_sites, additional_temp_sites, temp_HSPs, 
								 target_site, query_site, target_end_site, query_end_site, 
								 blast_hsp_index, blastHSPs);

#ifdef DEBUG

					dataMgr.outFile << "both have no stop, use old exon, HSPs for this last exon:" << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

						return false;

						}

					}

					

				}

		}

		else //if (isHeadExon && query_site > END_EXON_LEN )

		if (!isHeadExon && query_end_site < dataMgr.query_len - END_EXON_LEN) //last exon

		{

#ifdef DEBUG

			dataMgr.outFile << "Repairing tail exon " << target_site << "-" << target_end_site << " (query " << query_site 

				<<"-" << query_end_site << ")" << "\n";

#endif
			//Modified: now also check for the possibility of ending at the second-to-last exon
			//this is the end of second-to-last exon, extend it to the next stop codon and store it 
			//for later comparing with the other options
			second_end_temp_sites.clear();
			if (temp_sites.size()>3)
			{
			int second_end_site;
			int second_end_index = temp_sites.size()-3;
			if (target_site > 0) //pos strand
				second_end_site = CalcEndPos(temp_sites[second_end_index].first-temp_sites[second_end_index].frame,
					0, true, chr_seq); //gene_end is not used any more in CalcEndPos, so just leave 0 here
			else
				second_end_site = -CalcEndPos(-temp_sites[second_end_index].first+temp_sites[second_end_index].frame, 
					0, false, chr_seq);
			for (int s=0; s<temp_sites.size()-3; s++)
				second_end_temp_sites.push_back(temp_sites[s]);
			second_end_temp_sites.push_back(ExonSiteInfo(second_end_site,0,0,0)); //store the bunch of temp_sites with second_end_site as its gene end
			}

			//do the usual repair stuff (try new extension, compare with existing exon, keep the better one, etc.)
					double max_score=0, cur_score;

					float max_pid = 0, cur_pid;



					if (target_site > 0)

						tgt_str_len = GetSubstrFromVecStrs_ForRepair(chr_seq, true, target_end_site, EXPLORE_END_EXON_LEN, target_seq);

					else

						tgt_str_len = GetSubstrFromVecStrs_NegRev_ForRepair(chr_seq, false, -target_end_site-EXPLORE_END_EXON_LEN-1, EXPLORE_END_EXON_LEN, target_seq);

					StrToLower(target_seq);



					HSP_Gene_Pair newHSP[3];

					HSP_Gene_Pair* curHSP = new HSP_Gene_Pair;

					Input_Alignment newAlign[3];

					int j=0;

					string result_seq[3];

					q_seq = query_seq.substr(query_end_site);

					for (int i=0; i<3; i++)

					{

						DNA2AA(target_seq, i, result_seq[i]);

#ifdef DEBUG

					dataMgr.outFile << "target_seq:" << "\n" << target_seq << "\n";

					dataMgr.outFile << result_seq[i] << "\n";

#endif



						

						//if (target_site > 0)

							cur_score = GetLocalAlignment(q_seq, result_seq[i], 

								query_end_site+1, target_end_site + 1 + i, max_alignment_HSP_ID, newHSP[i], newAlign[i], 

								true, dataMgr.outFile, align_pos[i]);

						//else //same?

						//	cur_score = GetLocalAlignment(q_seq, result_seq[i], 

						//		query_end_site+1, -(-target_end_site - 1 - i), max_alignment_HSP_ID, newHSP[i], newAlign[i],

						//		true, dataMgr.outFile, align_pos[i]);



						cur_pid = newHSP[i].pid;



						if (cur_score == 0)

							continue;

					

if (EXTEND_HSP_BY_SCORE)

{

						if (max_score < cur_score || 

							(target_site > 0 && max_score == cur_score && newHSP[i].HSP_start < curHSP->HSP_start) || 

							(target_site < 0 && max_score == cur_score && newHSP[i].HSP_end > curHSP->HSP_end))

						{

							j = i;

							max_score = cur_score;

							max_pid = cur_pid;

							*curHSP = newHSP[i];

						}

}

else

{

						if (max_pid < cur_pid || 

							(target_site > 0 && max_pid == cur_pid && newHSP[i].HSP_start < curHSP->HSP_start) || 

							(target_site < 0 && max_pid == cur_pid && newHSP[i].HSP_end > curHSP->HSP_end) )

						{

							j = i;

							max_score = cur_score;

							max_pid = cur_pid;

							*curHSP = newHSP[i];

						}

}

					}



					//dataMgr.HSP_gene[chr_index].push_back(newHSP[j]);

					//newHSPs.push_back(&(dataMgr.HSP_gene[chr_index].back()));

					if (max_score == 0 || max_score < REPAIR_HSP_MIN_INIT_SCORE) //<= 0)

					{

						temp_HSPs.clear();

#ifdef DEBUG

					dataMgr.outFile << "use old exon: temp_HSPs for this exon: " << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

					return false; //no valid extra alignment, just use the old exon

					}

					else

					{

					//if the score of init alignment is big enough, try extend this alignment to both directions, just like BLAST!

					//extend curHSP to both directions if possible, stop extending when the score of alignment is 

					//lower than ("max_score" - REPAIR_HSP_EXTEND_SCORE_DROP), "max_score" is the running max score

					//if (max_score >= REPAIR_HSP_MIN_INIT_SCORE)
					//{

						string tmpStr = "";

						if (align_pos[j][3]+1 < result_seq[j].length())

							tmpStr = result_seq[j].substr(align_pos[j][3]+1);


						string result_begin_str = result_seq[j].substr(0, align_pos[j][2]);

						if (target_site > 0)

						ExtendAlignment(max_score, query_seq, 

							//result_seq[j].substr(0, align_pos[j][2]), //this may also be empty string
							result_begin_str, 

							tmpStr, 

							chr_seq, true, curHSP, newAlign[j], dataMgr.outFile, 

							//target_end_site, INT_MAX);

							temp_HSPs.front().HSP_start, INT_MAX);

						else

						ExtendAlignment(max_score, query_seq, 

							//result_seq[j].substr(0, align_pos[j][2]), //this may also be empty string
							result_begin_str, 

							tmpStr, 

							chr_seq, false, curHSP, newAlign[j], dataMgr.outFile, 

							//-target_end_site, 0);

							temp_HSPs.back().HSP_end, 0);

#ifdef DEBUG

						dataMgr.outFile << "extended alignment: " << "\n";

						dataMgr.outFile << *curHSP << "\n";

						dataMgr.outFile << newAlign[j] << "\n";

#endif

					//}



					newHSPs.push_back(*curHSP);

					newHSP_ptrs.push_back(curHSP);

					dataMgr.input_alignments.insert(map<int, Input_Alignment>::value_type(max_alignment_HSP_ID, newAlign[j]));
					hspID_added_for_repair = max_alignment_HSP_ID; //set hspID_added_for_repair to the new HSP
#ifdef DEBUG
					dataMgr.outFile << "hspID_added_for_repair is: HSP[" << hspID_added_for_repair << "]\n";
#endif
					max_alignment_HSP_ID++;



					//now redo PreProcHSPs and LoadData4 etc.

					if (target_site > 0)

					{

						//for (j=0; j<=blast_hsp_index; j++)

						//	newHSPs.push_back(blastHSPs[j]);

						for (j=0; j<temp_HSPs.size(); j++)

						{

							newHSPs.push_back(temp_HSPs[j]);

						newHSPs.back().ID = max_alignment_HSP_ID; //also replicate the alignment so it won't be changed on the original blast alignment

						dataMgr.input_alignments.insert(map<int, Input_Alignment>::value_type(max_alignment_HSP_ID, 

							(*(dataMgr.input_alignments.find(temp_HSPs[j].ID))).second));

						max_alignment_HSP_ID++;

						}

					}

					else

					{

						//for (j=blastHSPs.size()-1; j>=blast_hsp_index; j--)

						//	newHSPs.insert(newHSPs.begin(), blastHSPs[j]);

						for (j=temp_HSPs.size()-1; j>=0; j--)

						{

							newHSPs.insert(newHSPs.begin(), temp_HSPs[j]);

						newHSPs.front().ID = max_alignment_HSP_ID; //also replicate the alignment so it won't be changed on the original blast alignment

						dataMgr.input_alignments.insert(map<int, Input_Alignment>::value_type(max_alignment_HSP_ID, 

							(*(dataMgr.input_alignments.find(temp_HSPs[j].ID))).second));

						max_alignment_HSP_ID++;

						}

					}

#ifdef DEBUG

					dataMgr.outFile << "newHSPs for repairing tail exon: " << "\n";

					for (j=0; j<newHSPs.size(); j++)

						dataMgr.outFile << (newHSPs[j]) << "\n";

#endif

					for (j=0; j<newHSPs.size(); j++)

						ptr_to_newHSPs.push_back(&(newHSPs[j]));

					PreProcHSPs((target_site>0), ptr_to_newHSPs, newHSP_ptrs, max_alignment_HSP_ID);

					//need to update newHSPs? since PreProcHSPs may have inserted new HSPs...nah, just use ptr_to_newHSPs!

/*					int k=0;

					for (j=0; j<ptr_to_newHSPs.size(); j++) 

					{

						if (ptr_to_newHSPs[j]->ID != newHSPs[k].ID)

							newHSPs.insert(newHSPs.begin()+k, *(ptr_to_newHSPs[j]));

						k++;

					}

*/



#ifdef DEBUG

					dataMgr.outFile << "newHSPs after preproc: " << "\n";

					for (j=0; j<ptr_to_newHSPs.size(); j++)

						dataMgr.outFile << *(ptr_to_newHSPs[j]) << "\n";

#endif



					start_site = target_site; //fix start_site first

SELECT_DONOR_ACCEPTOR6_TAILEXON:

					if (SPLICE_SEGMENT_VERSION == 1)

					{

						HSPs_dup.clear();

						dataMgr.input_alignments_HSPs_dup.clear();

					}

					if (!LoadData4_forGenePredict6(ptr_to_newHSPs, chr_seq, true, false))

					{ //do not need to record now, since it's the last exon (has been recorded previously)!

						//temp_sites.push_back(target_site);

						//temp_sites.push_back(target_end_site);

						temp_HSPs.clear();

#ifdef DEBUG

					dataMgr.outFile << "use old exon: temp_HSPs for this exon: " << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif
					return false; //something is wrong, just use the old exon

					}



					vector< vector<ExonSiteInfo> > all_alternative_acceptors, all_alternative_donors;

					bool has_possible_exon;
#ifdef COMPUT_EXON_FULL_STEP_BACK

					if (!ComputeExons(ptr_to_newHSPs, chr_seq, exon_start_frame, 0, //exon_end_frame, //tail exon must end at frame 0
						all_alternative_acceptors, all_alternative_donors, has_possible_exon, true))
						if (has_possible_exon)
							goto SELECT_DONOR_ACCEPTOR6_TAILEXON;

#else

					while (!ComputeExons(ptr_to_newHSPs, chr_seq, exon_start_frame, 0, 
						all_alternative_acceptors, all_alternative_donors, has_possible_exon, true))
						if (!has_possible_exon)
							break;

#endif



#ifdef DEBUG

					dataMgr.outFile << "computed new exon candidates" << "\n";

#endif
					if (!has_possible_exon)
					{
						temp_HSPs.clear();
						return false;
					}

					//multiply all_alternative_acceptors/donors into the old all_alternative_acceptors/donors?

					int frame = exon_start_frame; //for last exon, its start frame may not be 0!

					string final_alignment;

					string left_chars((3-frame)%3, 'x'); //initialize to some 'x's, so the broken frame won't be stop codon

					int prev_donor_site = 0;

					//temp_sites.clear();

					int site_before_stop;

					for (j=0; j<all_alternative_acceptors[0].size(); j++) //use the first one, for now

					{
						string tmp_str;
						bool found_stop;
						if (HasStopCodon(all_alternative_acceptors[0][j].first, all_alternative_donors[0][j].first, 

							chr_seq, frame, site_before_stop, final_alignment, prev_donor_site, left_chars, 
							false, dataMgr.cDNA_os, false, tmp_str, j==all_alternative_acceptors[0].size()-1, found_stop))

						{ //if there's stop in new exons, discard them

							//temp_sites.push_back(target_site);

							//temp_sites.push_back(target_end_site);

							temp_HSPs.clear();

#ifdef DEBUG

					dataMgr.outFile << "new exons have stop, use old exon: temp_HSPs for this exon: " << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

							return false;

						}

						else

						{

							prev_donor_site = all_alternative_donors[0][j].first;

						}

					}

					Input_Alignment newAlign;

					float align_pid;

					string query_end_str = query_seq.substr(query_site-1);
					int align_score;
#ifdef PID_BEFORE_SCORE
					align_score = GetGlobalAlignment_PID_CompScore(query_end_str, final_alignment, newAlign, align_pid, dataMgr.outFile);
#else
					align_score = GetGlobalAlignment_scorematrix(//query_seq.substr(query_site-1), 
						query_end_str, final_alignment, newAlign, align_pid, dataMgr.outFile);
#endif
					align_compare align_new(align_score, align_pid);

					dataMgr.outFile << "New End spliced alignment:\n" << newAlign << "\n";
					dataMgr.outFile << "New End PID:" << align_pid << "\n";

					

					//reset stuff

					frame = exon_start_frame;

					final_alignment = "";

					left_chars.assign((3-frame)%3, 'x');

					prev_donor_site = 0;

					int num_of_new_exons = all_alternative_acceptors[0].size();

#ifdef DEBUG

					dataMgr.outFile << "num_of_new_exons:" << num_of_new_exons << "\n";

#endif

					string tmp_str;
					bool found_stop;
					if (HasStopCodon(target_site, target_end_site, chr_seq, frame, site_before_stop, final_alignment, prev_donor_site, left_chars, 
						false , dataMgr.cDNA_os, false, tmp_str, true, found_stop))

					{

						temp_sites.pop_back();

						temp_sites.pop_back();

						for (j=0; j<num_of_new_exons; j++) //use the first one, for now

						{

							temp_sites.push_back(all_alternative_acceptors[0][j]);

							temp_sites.push_back(all_alternative_donors[0][j]);

							//no need to record temp_HSPs, since there's no next exon to check

							additional_temp_sites.push_back(all_alternative_acceptors[0][j]);

							additional_temp_sites.push_back(all_alternative_donors[0][j]);

							//output some additional info, for inspection purpose only
							if (all_alternative_acceptors[0][j].first > target_end_site)
								dataMgr.outFile << "NEW EXON after repair:" << all_alternative_acceptors[0][j].first << "-" << all_alternative_donors[0][j].first << "\n";
						}

						temp_HSPs.clear();

#ifdef DEBUG

					dataMgr.outFile << "old exon has stop, use new exons, HSPs for this last exon:" << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

					return true;

					}

					else

					{

						float ori_align_pid;

						string query_end_str = query_seq.substr(query_site-1);
						int ori_align_score;
#ifdef PID_BEFORE_SCORE
						ori_align_score = GetGlobalAlignment_PID_CompScore(query_end_str, final_alignment, newAlign, ori_align_pid, dataMgr.outFile);
#else
						ori_align_score = GetGlobalAlignment_scorematrix(//query_seq.substr(query_site-1), 
							query_end_str, final_alignment, newAlign, ori_align_pid, dataMgr.outFile);
#endif
						align_compare align_old(ori_align_score, ori_align_pid);

						dataMgr.outFile << "Initial End spliced alignment:\n" << newAlign << "\n";
						dataMgr.outFile << "Initial End PID:" << ori_align_pid << "\n";
#ifdef DEBUG
							dataMgr.outFile << "align_new:" << align_new << ";align_old:" << align_old << "\n";
							dataMgr.outFile << "align_new is bigger?" << (align_new > align_old) 
								<< "(align_new:pid is bigger?" << (align_new.align_pid > align_old.align_pid)
								<< ")" << "\n";
#endif

						//if (ori_align_score < align_score || (ori_align_score == align_score && ori_align_pid < align_pid))
						if (align_new > align_old)
						{
							temp_sites.pop_back();

							temp_sites.pop_back();

							for (j=0; j<num_of_new_exons; j++) //use the first one, for now

							{

								temp_sites.push_back(all_alternative_acceptors[0][j]);

								temp_sites.push_back(all_alternative_donors[0][j]);

								additional_temp_sites.push_back(all_alternative_acceptors[0][j]);

								additional_temp_sites.push_back(all_alternative_donors[0][j]);
							//output some additional info, for inspection purpose only
							if (all_alternative_acceptors[0][j].first > target_end_site)
								dataMgr.outFile << "NEW EXON after repair:" << all_alternative_acceptors[0][j].first << "-" << all_alternative_donors[0][j].first << "\n";

							}

							temp_HSPs.clear();

#ifdef DEBUG

					dataMgr.outFile << "both have no stop, use new exons, HSPs for this last exon:" << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

					return true;

						}

						else

						{

							//temp_sites.push_back(target_site);

							//temp_sites.push_back(target_end_site);

							temp_HSPs.clear();

#ifdef DEBUG

					dataMgr.outFile << "both have no stop, use old exon, HSPs for this last exon:" << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

					return false;

						}

					}

					

					

					}



			}

			else //no need to repair, old exon is good already!

			{

				int j;

					if (isHeadExon)

					{

						temp_sites.push_back( ExonSiteInfo(target_site, query_site, -1, 0)); //hsp_index doesn't matter?

						temp_sites.push_back( ExonSiteInfo(target_end_site, query_end_site, -1, 

							((target_end_site - target_site + 1)% 3 )));



						//additional_temp_sites.push_back( ExonSiteInfo(target_site, query_site, -1, 0)); //hsp_index doesn't matter?

						//additional_temp_sites.push_back( ExonSiteInfo(target_end_site, query_end_site, -1, 

						//	((target_end_site - target_site + 1)% 3 )));

						

						temp_HSPs.clear(); //now update temp_HSPs... no need! since it's the last exon already!

						if (target_site < 0)

						{

							for ( j=0; j<=blast_hsp_index; j++)

								temp_HSPs.push_back(*(blastHSPs[j]));

						}

						else

						{

							//for (j=blastHSPs.size()-1; j>=blast_hsp_index; j--)

								//temp_HSPs.insert(temp_HSPs.begin(), blastHSPs[j]);

							for ( j=blast_hsp_index; j<blastHSPs.size(); j++)

								temp_HSPs.push_back(*(blastHSPs[j]));

						}

					}

					else

					{

						temp_HSPs.clear();

					}



#ifdef DEBUG

					dataMgr.outFile << "no need to repair, use old exon: temp_HSPs for this exon: " << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

					return false;

			}





}



//pre-processing steps for GenePredict, with HSPs, clean them

void ACCP_DONR_graph::PreProcHSPs(bool isPosStrand, vector<HSP_Gene_Pair*>& HSPs, 

								vector<HSP_Gene_Pair*>& newHSP_ptrs, int& max_alignment_HSP_ID)

{

		//compute start and end positions (start of 1st exon, end of last exon)

//		CalcStartEndPos2((*groupMapIt).first.isPosStrand, (*groupMapIt).second, (*chrMapIt).second);

		//MODIFIED: now just use a temporary start_site and end_site...

		if (isPosStrand)

		{

			start_site = HSPs.back()->HSP_start;

			end_site = HSPs.front()->HSP_end;

		}

		else

		{

			start_site = -HSPs.front()->HSP_end;

			end_site = -HSPs.back()->HSP_start;

		}



		//merge overlapping HSPs (if they're in same frame) first

		//"vector<HSP_Gene_Pair*> (*groupMapIt).second" may be changed by this

		//this also groups HSPs on frames
//#ifdef VERBOSE
		if (VERBOSE)
		cout << "merging..." << "\n";
//#endif
		MergeOverlapHSPs(isPosStrand, HSPs);
//#ifdef VERBOSE
		if (VERBOSE)
		cout << "merge done" << "\n";
//#endif

#ifdef DEBUG
		dataMgr.outFile << "after merging HSPs" << "\n";
		for (vector<HSP_Gene_Pair*>::iterator hsps_it = HSPs.begin(); hsps_it != HSPs.end(); hsps_it++)
			dataMgr.outFile << "HSP[" << (*hsps_it)->ID << "]:" << (*hsps_it)->HSP_start << "-" << (*hsps_it)->HSP_end 
				<< "(" << (*hsps_it)->gene_start << "-" << (*hsps_it)->gene_end << ")" << "\n";
#endif

	//ADDED: re-order the HSPs, since they may be out of order because of merge/cutoffs

	map<int, HSP_Gene_Pair*> HSP_sort_map; //sort HSPs by their HSP_start (after merge, every HSP will be non-overlapping)
	vector<HSP_Gene_Pair*>::iterator ptrHSPIt;

	for (ptrHSPIt = HSPs.begin(); ptrHSPIt != HSPs.end(); ptrHSPIt++)

		HSP_sort_map.insert(map<int, HSP_Gene_Pair*>::value_type((*ptrHSPIt)->HSP_start, *ptrHSPIt));

	HSPs.clear();

	for (map<int, HSP_Gene_Pair*>::reverse_iterator mapHSPIt = HSP_sort_map.rbegin(); mapHSPIt != HSP_sort_map.rend(); mapHSPIt++)

		HSPs.push_back((*mapHSPIt).second);	

	CutHeadOrTrailGapsInHSPs(HSPs);

	FixHSPsWithInFrameStop(isPosStrand, HSPs, newHSP_ptrs, max_alignment_HSP_ID); //this now also makes sure there's no gapped head or trail

	SortHSPsToColinear(isPosStrand, HSPs);//clean up HSPs to eliminate non-colinear ones (due to overlap merging and stop fixing)	
	CutHeadOrTrailGapsInHSPs(HSPs); //do this again, since the sort may introduce gap again
}

void ACCP_DONR_graph::SortHSPsToColinear(bool isPosStrand, vector<HSP_Gene_Pair*>& HSPs)
{
	bool sorted = false;
	while (!sorted && HSPs.size()>1)
	{
		sorted = true;
		int i, j;
		if (isPosStrand)
		{
			i = HSPs.size()-1;
			j=HSPs.size()-2;
		}
		else
		{
			j = HSPs.size()-1;
			i = HSPs.size()-2;
		}
		while (j>=0 && i>=0)
		{
			int k;
			if (isPosStrand)
				k = i;
			else
				k = j;

			if (HSPs[i]->gene_start >= HSPs[j]->gene_start || HSPs[i]->gene_end >= HSPs[j]->gene_end) //not colinear
			{
				sorted = false;
				//compare HSPs[i] and HSPs[j], and delete/cut one with less match
				int i_gene_start, i_gene_end, j_gene_start, j_gene_end;
				i_gene_end = HSPs[i]->gene_end;
				if (HSPs[i]->gene_end < HSPs[j]->gene_end)
					j_gene_end = HSPs[i]->gene_end;
				else
					j_gene_end = HSPs[j]->gene_end;
				j_gene_start = HSPs[j]->gene_start;
				if (HSPs[i]->gene_start < HSPs[j]->gene_start)
					i_gene_start = HSPs[j]->gene_start;
				else
					i_gene_start = HSPs[i]->gene_start;
#ifdef DEBUG
				dataMgr.outFile << "HSPs[i]:" << *(HSPs[i]) << "\n";
				dataMgr.outFile << "HSPs[j]:" << *(HSPs[j]) << "\n";
				dataMgr.outFile << "i_gene_start:" << i_gene_start << ";i_gene_end:" << i_gene_end << "\n";
				dataMgr.outFile << "j_gene_start:" << j_gene_start << ";j_gene_end:" << j_gene_end << "\n";
#endif
				map<int, Input_Alignment>::iterator it;
				it = dataMgr.input_alignments.find(HSPs[i]->ID);
				Input_Alignment& hsp_i_align = (*it).second;
				it = dataMgr.input_alignments.find(HSPs[j]->ID);
				Input_Alignment& hsp_j_align = (*it).second;
#ifdef DEBUG
				dataMgr.outFile << "HSPs[i]:\n" << hsp_i_align << "\n";
				dataMgr.outFile << "HSPs[j]:\n" << hsp_j_align << "\n";
#endif

				int i_gene_start_index, i_gene_end_index, j_gene_start_index, j_gene_end_index;
				if (!FindRealPos_DNA(hsp_i_align.query_align, HSPs[i]->gene_start, HSPs[i]->gene_end, 
					i_gene_start, i_gene_end, i_gene_start_index, i_gene_end_index))
				{
					cout << "cannot find hsp_i query overlap region properly" << "\n";
					exit(-1);
				}
				if (!FindRealPos_DNA(hsp_j_align.query_align, HSPs[j]->gene_start, HSPs[j]->gene_end, 
					j_gene_start, j_gene_end, j_gene_start_index, j_gene_end_index))
				{
					cout << "cannot find hsp_j query overlap region properly" << "\n";
					exit(-1);
				}
#ifdef DEBUG
				dataMgr.outFile << "i_gene_start_index:" << i_gene_start_index << ";i_gene_end_index:" << i_gene_end_index << "\n";
				dataMgr.outFile << "j_gene_start_index:" << j_gene_start_index << ";j_gene_end_index:" << j_gene_end_index << "\n";
#endif

				string str_i = hsp_i_align.match_align.substr(i_gene_start_index, i_gene_end_index - i_gene_start_index + 1);
				string str_j = hsp_j_align.match_align.substr(j_gene_start_index, j_gene_end_index - j_gene_start_index + 1);
				int i_identity = str_i.length() - CountCharInStr(str_i, ' ') - CountCharInStr(str_i, '+');
				int j_identity = str_j.length() - CountCharInStr(str_j, ' ') - CountCharInStr(str_j, '+');
				if ( i_identity < j_identity || (i_identity == j_identity && str_i.length() > str_j.length()))
				{
					if (i_gene_start > HSPs[i]->gene_start) //keep the part before i_gene_start, cut off rest of i
					{
						HSPs[i]->gene_end = HSPs[j]->gene_start-1;
						const_cast<string&>(hsp_i_align.target_align).erase(i_gene_start_index);
						const_cast<string&>(hsp_i_align.query_align).erase(i_gene_start_index);
						const_cast<string&>(hsp_i_align.match_align).erase(i_gene_start_index);
						if (isPosStrand)
							HSPs[i]->HSP_end = HSPs[i]->HSP_start + 3*(i_gene_start_index - CountCharInStr(hsp_i_align.target_align, '-')) - 1;
						else
							HSPs[i]->HSP_start = HSPs[i]->HSP_end - 3*(i_gene_start_index - CountCharInStr(hsp_i_align.target_align, '-')) + 1;
						int i_id = i_gene_start_index - CountCharInStr(hsp_i_align.match_align, ' ') - CountCharInStr(hsp_i_align.match_align, '+');
						HSPs[i]->pid = (float)i_id / i_gene_start_index;
#ifdef DEBUG
						dataMgr.outFile << "HSPs[i]:" << *(HSPs[i]) << "\n";
#endif
					}
					else //delete i altogether
					{
						HSPs.erase(HSPs.begin()+i);
#ifdef DEBUG
						dataMgr.outFile << "erased HSPs[i]" << "\n";
#endif
					}
				}
				else
				//if (j_identity < i_identity)
				{
					if (j_gene_end < HSPs[j]->gene_end)
					{
						HSPs[j]->gene_start = HSPs[i]->gene_end+1;
						const_cast<string&>(hsp_j_align.target_align).erase(0, j_gene_end_index+1);
						const_cast<string&>(hsp_j_align.query_align).erase(0, j_gene_end_index+1);
						const_cast<string&>(hsp_j_align.match_align).erase(0, j_gene_end_index+1);
						if (isPosStrand)
							HSPs[j]->HSP_start = HSPs[j]->HSP_end - 3*(hsp_j_align.target_align.length()-CountCharInStr(hsp_j_align.target_align, '-'))+1;
						else
							HSPs[j]->HSP_end = HSPs[j]->HSP_start + 3*(hsp_j_align.target_align.length()-CountCharInStr(hsp_j_align.target_align, '-'))-1;
						int j_id = hsp_j_align.match_align.length() - CountCharInStr(hsp_j_align.match_align, ' ') - CountCharInStr(hsp_j_align.match_align, '+');
						HSPs[j]->pid = (float)j_id / hsp_j_align.match_align.length();
#ifdef DEBUG
						dataMgr.outFile << "HSPs[j]:" << *(HSPs[j]) << "\n";
#endif
					}
					else
					{
						HSPs.erase(HSPs.begin()+j);
#ifdef DEBUG
						dataMgr.outFile << "erased HSPs[j]" << "\n";
#endif
					}
				}
			}

			if (isPosStrand)
			{
				i = k-1;
				j = i-1;
			}
			else
			{
				j = k-1;
				i = j-1;
			}
		}
	}
}

void ACCP_DONR_graph::CutHeadOrTrailGapsInHSPs(vector<HSP_Gene_Pair*>& HSPs)
{
	vector<HSP_Gene_Pair*>::iterator ptrHSPIt;
	vector< vector<HSP_Gene_Pair*>::iterator > ptrHSP_to_be_deleted;	
	for ( ptrHSPIt = HSPs.begin(); ptrHSPIt != HSPs.end(); ptrHSPIt++ )
	{
		//cut off possible gapped heading/trailing portions of HSPs (trailing or heading gaps at either end of HSP)
		Input_Alignment& align = (*(dataMgr.input_alignments.find((*ptrHSPIt)->ID))).second;
#ifdef DEBUG
		dataMgr.outFile << "left hsp before processing: " << *(*ptrHSPIt) << "\n";
		dataMgr.outFile << align;
#endif
		if (!CutHeadOrTrailGapsInAlignment(align, (*ptrHSPIt)->gene_start, (*ptrHSPIt)->gene_end, 
			(*ptrHSPIt)->HSP_start, (*ptrHSPIt)->HSP_end))
			ptrHSP_to_be_deleted.push_back(ptrHSPIt);

#ifdef DEBUG
		dataMgr.outFile << "left hsp after processing: " << *(*ptrHSPIt) << "\n";
		dataMgr.outFile << (*(dataMgr.input_alignments.find((*ptrHSPIt)->ID))).second;
#endif
	}
	for (int i = ptrHSP_to_be_deleted.size()-1; i>=0; i--)
		HSPs.erase(ptrHSP_to_be_deleted[i]);
}


bool ACCP_DONR_graph::ComputeExons(vector<HSP_Gene_Pair*>& HSPs, 
								   //pair<int, char*>& chr_seq, 
								   vector<string>& chr_seq, 
								   int start_frame, int end_frame, 
								   vector< vector< ExonSiteInfo > >& all_alternative_acceptors, 
								   vector< vector< ExonSiteInfo > >& all_alternative_donors, bool& has_possible_exon, 
								   bool repair_only)
{
	has_possible_exon = true;

#ifdef DEBUG
	OutputDonorAcceptorWithAllInfo(start_site>0);
#endif

	int tmp;

	acceptor_region_index_copy.clear();

	for (tmp=0; tmp<acceptor_region_index.size(); tmp++)

		acceptor_region_index_copy.push_back(acceptor_region_index[tmp]);

	donor_region_index_copy.clear();

	for (tmp=0; tmp<donor_region_index.size(); tmp++)

		donor_region_index_copy.push_back(donor_region_index[tmp]);



		//int prev_intron_end = start_site;

		vector< ExonSiteInfo > prev_intron_end;



		vector<int> prev_intron_start_donor_index;

		//vector< pair<int, int> > neg_exons; //temporary store negative strand splice sites

		vector< ExonSiteInfo > all_prev_donors; //used to try to possibly skip previously-computed introns!

		vector<int> all_prev_donors_tail;

		vector< ExonSiteInfo > all_prev_acceptors;

		int last_used_donor_index;

		vector<int> count_alternative_splice_pairs; //used to keep count of how many alternative ways of splicing happened for each splicing segment region



		int start_site_hsp_rel_index=0, end_site_hsp_rel_index=0; //the RELATIVE index of HSP that corresponds to start_site/end_site



		int i=0;		

		//int prev_intron_end = start_site;

		prev_intron_end.clear();

		HSP_Gene_Pair *startHSP, *endHSP;

		if (start_site > 0)

		{

			startHSP = HSPs.back();

			endHSP = HSPs.front();

			prev_intron_end.push_back( ExonSiteInfo(start_site, startHSP->gene_start, HSPs.size()-1, start_frame) ); //initialize to one start_site

		}

		else

		{

			startHSP = HSPs.front();

			endHSP = HSPs.back();

			prev_intron_end.push_back( ExonSiteInfo(start_site, startHSP->gene_start, 0, start_frame) ); //initialize to one start_site

		}

		all_prev_donors_tail.clear();

		all_prev_donors_tail.push_back((3-start_frame)%3); //for RepairExon functions, "start_site" may not start at frame 0!



		//bunch of declarations, put here only because otherwise the 1st "goto CLEAN_UP" will result in compile error

		int ai, di;

		//vector< vector< ExonSiteInfo > > all_alternative_acceptors, all_alternative_donors;

		/*UPDATE: the single exon case will be checked later (when checking last exon)
		if (donors.empty() ) //there is only start_site and end_site
			if (!CheckSglExon(HSPs, chr_seq, start_frame, end_frame, start_site, 
				all_alternative_acceptors, all_alternative_donors, has_possible_exon, repair_only))
				return false; //if CheckSglExon returned true, then go on to next statement
		*/

		bool break_loop = false;
		while (i<donors.size())

		{

			last_used_donor_index = i;

#ifdef DEBUG

			dataMgr.outFile << "donor " << i << " (region: " << donor_region_index[i] << ")" << "\n";

#endif

			//if (donors[i] <= prev_intron_end || 

			//	(prev_intron_end > start_site && donors[i] - prev_intron_end + 1 < MIN_INTERNAL_EXON_LEN) ) //once previous intron is fixed, the following must follow it

			if (donors[i] <= prev_intron_end.back().first || 

				(prev_intron_end.front().first > start_site && 

				donors[i] - prev_intron_end.front().first + 1 < MIN_INTERNAL_EXON_LEN) )//once previous intron is fixed, the following must follow it

			{

#ifdef DEBUG

			dataMgr.outFile << "donor " << i << " distance from previous acceptor is not valid" << "\n";

#endif

				i++;

				continue;

			}



			int match_count=0;



			int donor_region = donor_region_index[i];

			int acc_region = donor_region+1;

			vector<int>::iterator acc_pos = lower_bound(acceptor_region_index.begin(), acceptor_region_index.end(), acc_region);

			if (acc_pos == acceptor_region_index.end()) //no acceptor after this point, we stop

				break;

			int acc_start_index = acc_pos - acceptor_region_index.begin();

			//MODIFIED: there is no acceptor in this <donor,acceptor>pair region, immediately 
			//increase donor region and go to next round
			if (acceptor_region_index[acc_start_index] > acc_region)
			{
				//change the current bunch of donors to next region
				while (i < donors.size() && donor_region_index[i] == donor_region)
					i++;
				for (int di = last_used_donor_index; di < i; di++)
					donor_region_index[di]++;

				//reset the iterator i back to last_used_donor_index
				i = last_used_donor_index;
				continue;
			}

			Input_Alignment donor_align, acceptor_align;

			int donor_align_start, donor_align_end, acceptor_align_start, acceptor_align_end;

			string splice_align;
			int splice_align_score;

			//multimap<RankAlignment, pair<int,int> > result_intron;

			multimap<RankAlignment, IntronAndDonorTail> result_intron;

			//result_intron.clear();



//			char missing_query_aa = '?';



			int j;
			while (i < donors.size() && donor_region_index[i] == donor_region)

			{

#ifdef DEBUG

				dataMgr.outFile << "donor " << i << ":" << donors[i] << " (region: " << donor_region_index[i] << ")" << "\n";

#endif

				//UPDATED: check donor against previous intron end for in-frame stop codon

				//UPDATED: also check if current donor matches frame with prev_intron_end (in case they are from different HSPs... possible?)

				int prev_intron_end_head = (3 - all_prev_donors_tail.back())%3;

				if ((donors[i] - prev_intron_end.back().first + 1 - donor_tail[i] - prev_intron_end_head)%3 > 0) //if frame doesn't match, skip

				{

#ifdef DEBUG

			dataMgr.outFile << "donor " << i << " frame doesn't match previous acceptor" << "\n";

#endif

					i++;

					continue;

				}



				string testStr, codonStr;

				if (donors[i] > 0)

					GetSubstrFromVecStrs(chr_seq, true, prev_intron_end.back().first-1, 

						donors[i]-prev_intron_end.back().first+1-donor_tail[i], testStr);

				else

					GetSubstrFromVecStrs_NegRev(chr_seq, false, -donors[i]-1+donor_tail[i], 

						donors[i]-prev_intron_end.back().first+1-donor_tail[i], testStr);

				StrToLower(testStr);

/*				int tmp_start_pos = testStr.length();

				bool has_stop = false;

				while (!has_stop && tmp_start_pos >= 3)

				{

					codonStr = testStr.substr(tmp_start_pos-3, 3);

					map<string, char>::iterator codon_it = DNA_CODON_TBL_SL.find(codonStr);

					if (codon_it != DNA_CODON_TBL_SL.end())

					{

						if ((*codon_it).second == '*')

							has_stop = true;

					}

					tmp_start_pos -= 3;

				}

				if (has_stop) //if there's in-frame stop, the rest of donors will also disqualify, so we'll need to skip them

*/

				int stop_pos;

				if (HasInFrameStop(testStr, false, stop_pos, false))

				{

#ifdef DEBUG

			dataMgr.outFile << "donor " << i << " has in-frame stop from previous acceptor" << "\n";

#endif

					//if there's supposed to be more donor_regions after this, we try fix it by removing a HSP and go back do everything again

					//also: only do this when result_intron is empty! (there's no current valid result)

					//if (donor_region == 0 || donor_region == 1 || 

					//	(result_intron.empty() && donor_region < donor_region_index.back()-1))

					if (result_intron.empty())						
					{
						//if (donor_region < donor_region_index.back()-1 || donor_region == 0)
						//ALWAYS remove 1 hsp before stop?
						//{
						return RemoveHSPAndResetComputeExons(HSPs, prev_intron_end.back().first+stop_pos, chr_seq, 
							stop_pos > (end_site-prev_intron_end.back().first+1)/2, has_possible_exon, 
							all_alternative_acceptors, all_alternative_donors);
						/*}
						else
						{
#ifdef DEBUG
							dataMgr.outFile << "break_loop entirely, has possible exon\n";
#endif
							break_loop = true; //end the whole loop
							break;
						}*/
					}
					else //already has some results, then just skip some donors
					{
#ifdef DEBUG
						dataMgr.outFile << "skip rest of donors, has possible exon\n";
#endif
						while (i < donors.size() && donor_region_index[i] == donor_region)
							i++; //forward the pointer to skip rest of donors in this region
						break;
					}
				}

				j = acc_start_index;



				//if (!GetDonorHSP((*groupMapIt).second, donors[i], donor_align, donor_align_start, donor_align_end)) //donors[i]

				int cur_donor_hsp_index;

				char donor_border_query_aa_front=' ', donor_border_query_aa_end=' ';

				//if (SPLICE_SEGMENT_VERSION == 1)

					cur_donor_hsp_index = GetDonorHSP(HSPs, donors[i], donor_HSP_ID[i], donor_region, //donor_acceptor_HSP_ID[donor_region], 

						donor_align, donor_align_start, donor_align_end, donor_border_query_aa_front, donor_border_query_aa_end, 

						chr_seq);

				//else

				//	cur_donor_hsp_index = GetDonorHSP_S2((*groupMapIt).second, donors[i], 

				//		donor_align, donor_align_start, donor_align_end, donor_border_query_aa_front, donor_border_query_aa_end);

				if (cur_donor_hsp_index == -1)

				{

#ifdef DEBUG

			dataMgr.outFile << "donor " << i << " something wrong with hsp_index" << "\n";

#endif

					i++;

					continue;

				}



				



				while (j < acceptors.size() && acceptor_region_index[j] == acc_region)

				{

					//for each viable donor/acceptor pair

					if ( acceptors[j] > donors[i] && acceptor_head[j] == (3 - donor_tail[i]) % 3 

						&& (acceptors[j] - donors[i]) > MIN_INTRON_LEN) //if frame matches && intron size is at least MIN_INTRON_LEN

					{

#ifdef DEBUG

						dataMgr.outFile << "matching acceptor " << j << ":" << acceptors[j] << " (region: " << acceptor_region_index[j] << ")" << "\n";

#endif
/*						//MODIFIED: now if acceptor is in the last region, also check acceptor with the end_site to see if they match!
						if (acc_region == acceptor_region_index.back())
						{
							if ((end_site - acceptors[j] + 1 - acceptor_head[j] - end_frame) % 3 > 0)
							{
								j++;
								continue;
							}
						}
*/
						int cur_acceptor_hsp_index;

						char acceptor_border_query_aa_front = ' ', acceptor_border_query_aa_end = ' ';

						//if (SPLICE_SEGMENT_VERSION == 1)

							cur_acceptor_hsp_index = GetAcceptorHSP(HSPs, acceptors[j], acceptor_HSP_ID[j], acc_region, //donor_acceptor_HSP_ID[acc_region], 

								acceptor_align, acceptor_align_start, acceptor_align_end, acceptor_border_query_aa_front, acceptor_border_query_aa_end, 

								chr_seq);

						//else

						//	cur_acceptor_hsp_index = GetAcceptorHSP_S2((*groupMapIt).second, acceptors[j], 

						//		acceptor_align, acceptor_align_start, acceptor_align_end, acceptor_border_query_aa_front, acceptor_border_query_aa_end);



						if (cur_acceptor_hsp_index == -1) //acceptors[j]

						{

#ifdef DEBUG

			dataMgr.outFile << "something wrong with acceptor " << j << "'s hsp_index" << "\n";

#endif

							j++;

							continue;

						}



#ifdef DEBUG

						dataMgr.outFile << "donor " << donors[i] << " (query:" << donor_align_start << "-" << donor_align_end 

							<< ") (" << donor_prev_2nt[i] << "=" << acceptor_next_2nt[j] << ") align:\n" << donor_align << "\n";

						dataMgr.outFile << "acceptor " << acceptors[j] << " (query:" << acceptor_align_start << "-" << acceptor_align_end

							<< ") align:\n" << acceptor_align << "\n";

#endif

						//ComputeSpliceAlignment(donor_align, acceptor_align, donor_align_start, donor_align_end, 

						//	acceptor_align_start, acceptor_align_end, splice_align);

						//int non_used;

/*						if (donor_border_query_aa_end != ' ')

							missing_query_aa = donor_border_query_aa_end;

						else

							if (acceptor_border_query_aa_front != ' ')

								missing_query_aa = acceptor_border_query_aa_front;

*/

						//if (donors[i] == 749659 && acceptors[j] == 749701)

						//	int stophere=1;



						/*if (!ComputeSpliceAlignment_SplitOL1(donor_align, acceptor_align, donor_align_start, donor_align_end, 

							acceptor_align_start, acceptor_align_end, splice_align, non_used, non_used, non_used, non_used, 

							//missing_query_aa, 

							donor_prev_2nt[i], acceptor_next_2nt[j], donor_tail[i], acceptor_head[j]))*/
						if (!ComputeSpliceAlignment_SplitOL1_ForGetHSPSpliceAlignment1(donor_align, acceptor_align, 
							donor_align_start, donor_align_end, acceptor_align_start, acceptor_align_end, 
							splice_align, splice_align_score, 
							donor_prev_2nt[i], acceptor_next_2nt[j], donor_tail[i], acceptor_head[j]))

						{

#ifdef DEBUG

			dataMgr.outFile << "something wrong with splice alignment" << "\n";

#endif

							j++;

							continue;

						}

#ifdef DEBUG

						dataMgr.outFile << "splice align:\n" << splice_align << "\n";

#endif						

						int identity = splice_align.length()-CountCharInStr(splice_align, ' ') - CountCharInStr(splice_align, '+');

						int match_minus_gap_mismatch = identity - (splice_align.length() - identity);

						float cur_pid = (float)identity / splice_align.length();

						//result_intron.insert(multimap<RankAlignment, pair<int,int> >::value_type(RankAlignment(cur_pid, splice_align.length(), match_minus_gap_mismatch), 

						//	pair<int, int>(donors[i], acceptors[j])));

						result_intron.insert(multimap<RankAlignment, IntronAndDonorTail>::value_type(

							RankAlignment(align_compare(splice_align_score, cur_pid), splice_align.length(), match_minus_gap_mismatch), //, donor_align.match_align.length()+acceptor_align.match_align.length()), 

							IntronAndDonorTail(donors[i], acceptors[j], donor_tail[i], cur_donor_hsp_index, cur_acceptor_hsp_index, 
								//DNAToProteinStylePos_EndSite(donor_align_end), DNAToProteinStylePos_StartSite(acceptor_align_start),
								donor_align_end, acceptor_align_start, 
								donor_region)));



						match_count++;

					}

					j++;

				}

				i++;

			}

			if (break_loop)
				break;

			/*if (acc_region == acceptor_region_index.back() && result_intron.empty()) //if it's the last region already
			{
				has_possible_exon = false;
				return false;
			}*/

			int count_result = 0;

			if (match_count)

			{

				multimap<RankAlignment, IntronAndDonorTail>::iterator result_it, next_result_it;

#ifdef DEBUG

				dataMgr.outFile << "all possible introns before gap_center: " << "\n";

				for (result_it = result_intron.begin(); result_it != result_intron.end(); result_it++)

					dataMgr.outFile << (*result_it).first 

						<< "; " << (*result_it).second << "\n";

#endif

				

				//if find the best donor/acceptor pair 

				//(compare it with the case of both donor/acceptor is at gap_centers[donor_region])

				//if (donor_region < gap_centers.size())		
				//MODIFIED: donor_segments_pair was modified to include the last segment, so the last segment should be excluded from gap_center consideration
				if (donor_region < donor_segments_pair.size() - 1)
				{				

				//dataMgr.outFile << "gap_centers[" << donor_region << "]: " << gap_centers[donor_region] << "\n";





				result_it = result_intron.begin(); //compute gap_center alignment according to the HSP used for the so-far-best result

				//int cur_gap_center = ((*result_it).second.first + (*result_it).second.second ) / 2;

				//?revised: compute gap_center according to current donor_segments (the center of current intron region)

				int cur_gap_center = (donor_segments_pair[donor_region].exon_seg_end + donor_segments_pair[donor_region].intron_seg_end) / 2;



				//if (GetDonorHSP((*groupMapIt).second, gap_centers[donor_region], donor_align, donor_align_start, donor_align_end)) //gap_centers[donor_region]

				char donor_border_query_aa_front, donor_border_query_aa_end, acceptor_border_query_aa_front, acceptor_border_query_aa_end;

				//GetDonorSpecificHSP((*groupMapIt).second, (*result_it).second.donor_hsp_index, cur_gap_center, //gap_centers[donor_region], 

				//	donor_align, donor_align_start, donor_align_end, donor_border_query_aa_front, donor_border_query_aa_end);

				//int non_used;

				//int gap_center_hspID = FindDonorHSPIndex(cur_gap_center, HSPs, donor_region, 
				//		non_used, non_used, non_used, non_used, non_used);
				int gap_center_hspID = HSPs[(*result_it).second.donor_hsp_index]->ID;

				int gap_center_hspIndex_don, gap_center_hspIndex_acc;
				if ((gap_center_hspIndex_don = GetDonorHSP(HSPs, cur_gap_center, //donor_acceptor_HSP_ID[donor_region], 

					gap_center_hspID, donor_region, 

					donor_align, donor_align_start, donor_align_end, donor_border_query_aa_front, donor_border_query_aa_end, 

					chr_seq)) == -1)

					goto AFTER_GAP_CENTER6;
#ifdef DEBUG
				dataMgr.outFile << "donor_align:\n" << donor_align << "\n";
#endif

//				{

					//if (GetAcceptorHSP((*groupMapIt).second, gap_centers[donor_region]+1, acceptor_align, acceptor_align_start, acceptor_align_end)) //gap_centers[donor_region]

					//GetAcceptorSpecificHSP((*groupMapIt).second, (*result_it).second.acceptor_hsp_index, cur_gap_center+1, //gap_centers[donor_region]+1, 

					//	acceptor_align, acceptor_align_start, acceptor_align_end, acceptor_border_query_aa_front, acceptor_border_query_aa_end);

					//gap_center_hspID = FindAcceptorHSPIndex(cur_gap_center, HSPs, acc_region, 
					//	non_used, non_used, non_used, non_used, non_used);
					gap_center_hspID = HSPs[(*result_it).second.acceptor_hsp_index]->ID;

					int cur_gap_center_acc = cur_gap_center+1;
					if ((gap_center_hspIndex_acc = GetAcceptorHSP(HSPs, 
						cur_gap_center_acc, //cur_gap_center, //donor_acceptor_HSP_ID[acc_region], 
						gap_center_hspID, acc_region, 
						acceptor_align, acceptor_align_start, acceptor_align_end, acceptor_border_query_aa_front, acceptor_border_query_aa_end, 
						chr_seq)) == -1)

						goto AFTER_GAP_CENTER6;
#ifdef DEBUG
				dataMgr.outFile << "acceptor_align:\n" << acceptor_align << "\n";
#endif

					//if these two HSPs (for gap center computation) are not in same frame, discard gap center!
					if ((HSPs[gap_center_hspIndex_acc]->HSP_start - HSPs[gap_center_hspIndex_don]->HSP_start) % 3 != 0)
						goto AFTER_GAP_CENTER6;

					//if these two HSPs have stop in between, then there must be intron between them, discard gap center!?
					//MODIFIED: still compute PID for gap center, but if this PID turns out to be the highest, do removeHSP stuff and ComputeExons again
					/*if (hspID_due_to_stop.find(HSPs[gap_center_hspIndex_acc]->ID) != hspID_due_to_stop.end())
						goto AFTER_GAP_CENTER6;*/

//					{

					//ComputeSpliceAlignment(donor_align, acceptor_align, donor_align_start, donor_align_end, 

					//	acceptor_align_start, acceptor_align_end, splice_align);

					int donor_tail_start_pos = donor_align.target_align.length()-2;

					if (donor_tail_start_pos < 0)

						donor_tail_start_pos = 0;



					/*if (!ComputeSpliceAlignment_SplitOL1(donor_align, acceptor_align, donor_align_start, donor_align_end, 

							acceptor_align_start, acceptor_align_end, splice_align, non_used, non_used, non_used, non_used, 

							donor_align.target_align.substr(donor_tail_start_pos), acceptor_align.target_align.substr(0, 2), 

							donor_align.target_align.length() % 3, acceptor_align.target_align.length() % 3))*/
					int gap_center_donor_tail_num, gap_center_acceptor_head_num;
					string donor_tail_str, acceptor_head_str;
					if (start_site > 0)
					{
						gap_center_donor_tail_num = (cur_gap_center - HSPs[gap_center_hspIndex_don]->HSP_start + 1)%3;
						gap_center_acceptor_head_num = (HSPs[gap_center_hspIndex_acc]->HSP_end - cur_gap_center_acc + 1)%3;
						GetSubstrFromVecStrs(chr_seq, true, cur_gap_center-2, 2, donor_tail_str);
						GetSubstrFromVecStrs(chr_seq, true, cur_gap_center_acc-1, 2, acceptor_head_str);
					}
					else
					{
						gap_center_donor_tail_num = (HSPs[gap_center_hspIndex_don]->HSP_end + cur_gap_center + 1)%3;
						gap_center_acceptor_head_num = (-cur_gap_center_acc - HSPs[gap_center_hspIndex_acc]->HSP_start + 1)%3;
						GetSubstrFromVecStrs_NegRev(chr_seq, false, -cur_gap_center-1, 2, donor_tail_str);
						GetSubstrFromVecStrs_NegRev(chr_seq, false, -cur_gap_center_acc-2, 2, acceptor_head_str);
					}
					if (!ComputeSpliceAlignment_SplitOL1_ForGetHSPSpliceAlignment1(donor_align, acceptor_align, 
						donor_align_start, donor_align_end, acceptor_align_start, acceptor_align_end, 
						splice_align, splice_align_score, 
						//donor_align.target_align.substr(donor_tail_start_pos), acceptor_align.target_align.substr(0, 2), 
						donor_tail_str, acceptor_head_str, 
						gap_center_donor_tail_num, gap_center_acceptor_head_num))
							goto AFTER_GAP_CENTER6;
#ifdef DEBUG
					dataMgr.outFile << "splice_align:\n" << splice_align << "\n" 
						//<< donor_align.target_align.substr(donor_tail_start_pos) << "=" << acceptor_align.target_align.substr(0, 2)
						<< donor_tail_str << "=" << acceptor_head_str 
						<< "\n";
#endif

					int identity = splice_align.length()-CountCharInStr(splice_align, ' ') - CountCharInStr(splice_align, '+');

					int match_minus_gap_mismatch = identity - (splice_align.length() - identity);

					float cur_pid = (float)identity / splice_align.length();

					//MODIFIED: check if the highest PID is the "no-intron" case and the two HSPs have stop between them. then do removeHSP...
					if (gap_center_hspIndex_acc != gap_center_hspIndex_don && hspID_due_to_stop.find(HSPs[gap_center_hspIndex_acc]->ID) != hspID_due_to_stop.end())
					{
						const align_compare& top_align_compare = ((*(result_intron.begin())).first).pid;
						if (splice_align_score > top_align_compare.align_score || 
							(splice_align_score == top_align_compare.align_score && cur_pid >= top_align_compare.align_pid))
						{
							string testStr;
							if (start_site > 0)
								GetSubstrFromVecStrs(chr_seq, true, prev_intron_end.back().first-1, 
									donor_segments_pair[donor_region].intron_seg_end-prev_intron_end.back().first+1, testStr);
							else
								GetSubstrFromVecStrs_NegRev(chr_seq, false, -donor_segments_pair[donor_region].intron_seg_end-1, 
									donor_segments_pair[donor_region].intron_seg_end-prev_intron_end.back().first+1, testStr);
							StrToLower(testStr);
							int stop_pos;
							if (HasInFrameStop(testStr, false, stop_pos, false))
							{
								return RemoveHSPAndResetComputeExons(HSPs, prev_intron_end.back().first+stop_pos, chr_seq, 
									stop_pos > (end_site-prev_intron_end.back().first+1)/2, has_possible_exon, 
									all_alternative_acceptors, all_alternative_donors);
							}
							else
							{
								cout << "stop codon not found between two HSPs that were originally separated by stop?" << "\n";
								exit(-1);
							}
						}
						else
							goto AFTER_GAP_CENTER6;
					}

					//result_intron.insert(multimap<RankAlignment, pair<int,int> >::value_type(RankAlignment(cur_pid, splice_align.length(), match_minus_gap_mismatch), 

					//			pair<int, int>(0, 0))); //"0" signals "no intron"

					result_intron.insert(multimap<RankAlignment, IntronAndDonorTail>::value_type(

						RankAlignment(align_compare(splice_align_score, cur_pid), splice_align.length(), match_minus_gap_mismatch), // donor_align.match_align.length()+acceptor_align.match_align.length()), 
							IntronAndDonorTail(0, 0, 0, gap_center_hspIndex_don, gap_center_hspIndex_acc, //-1, -1, 
							//DNAToProteinStylePos_EndSite(donor_align_end), DNAToProteinStylePos_StartSite(acceptor_align_start), 
							donor_align_end, acceptor_align_start, 
							donor_region))); //"0" signals "no intron"
//					}

//				}



				}



				//If there's intron, output it immediately

				//for debug checking

				//multimap<RankAlignment, pair<int, int> >::iterator result_it;

AFTER_GAP_CENTER6:



#ifdef DEBUG

				dataMgr.outFile << "all possible introns: " << "\n";

				for (result_it = result_intron.begin(); result_it != result_intron.end(); result_it++)

					dataMgr.outFile << (*result_it).first 

						<< "; " << (*result_it).second.first << "-" << (*result_it).second.second << "\n";

#endif

				//use the top ranked one (what if there's more than one top-ranked?)

				//keep all top-ranked results! all_prev_donors_tail is no longer needed...



				//modified!! now don't just immediately output (may be skipped later!!!)
				result_it = result_intron.begin();
				RankAlignment curRank = (*result_it).first;	

				//MODIFIED: if there're multiple pairs of the same highest rank, then compute "optimal PID" to try to separate them
				multimap<RankAlignment, IntronAndDonorTail>::iterator results_to_be_erased_it = (result_intron.equal_range(curRank)).second;
				if (results_to_be_erased_it != result_intron.end())
					result_intron.erase(results_to_be_erased_it, result_intron.end());
#ifdef DEBUG
				dataMgr.outFile << "all possible introns:(result_intron size:" << result_intron.size() << ")" << "\n";
				for (result_it = result_intron.begin(); result_it != result_intron.end(); result_it++)
					dataMgr.outFile << (*result_it).first 
						<< "; " << (*result_it).second.first << "-" << (*result_it).second.second << "\n";
#endif
				align_compare op_align_max(INT_MIN,0);
				align_compare* op_align_cur = new align_compare[result_intron.size()];
				int c=0; 
				op_align_cur[0].align_score = INT_MIN; //initialize the first one
				op_align_cur[0].align_pid = 0;
				if (result_intron.size() > 1)
				{
					for (result_it = result_intron.begin(); result_it != result_intron.end(); result_it++, c++)
					{
						op_align_cur[c].align_score = ComputeSpliceAlignment_SplitOL1_Optimal(HSPs, chr_seq, (*result_it).second, 
							op_align_cur[c].align_pid);
						if (op_align_cur[c] > op_align_max)
							op_align_max = op_align_cur[c];
#ifdef DEBUG
						dataMgr.outFile << "current op_align_max:" << op_align_max << "; this align:" << op_align_cur[c] << "\n";
#endif
					}
				}

				//reset result_it
				result_it = result_intron.begin(); c=0;

				//now do the usual output (with possible alternatives)
				next_result_it = result_intron.begin();
				next_result_it++;
				bool fstResult = true;
				vector< ExonSiteInfo > cur_prev_intron_end(prev_intron_end); //hold the value
				//if ((*result_it).second.first != 0)
				while (result_it != result_intron.end() && 
					//(*result_it).first == curRank && 
					(*result_it).second.first != 0 )
				{
#ifdef DEBUG
					dataMgr.outFile << "current op_align_cur[" << c << "]:" << op_align_cur[c] << "\n";
#endif
					//if (op_align_pid[c] < op_align_pid_max) //skip those with lower optimal PIDs
					if (op_align_max > op_align_cur[c])
					{
#ifdef DEBUG
						dataMgr.outFile << "skip current op_align_cur[" << c << "], op_align_max:" << op_align_max << "\n";
#endif
						c++;
						result_it++;
						next_result_it++;
						continue;
					}						
					
						bool good_exon = true; //, good_intron = true;

						//check if it gives a good exon

						if (!all_prev_acceptors.empty() && //not the first exon

							//donor_region < donor_region_index.back()  && //not the last exon -> nah, no need!

							((*result_it).second.first - cur_prev_intron_end.front().first + 1) < MIN_INTERNAL_EXON_LEN) 

							good_exon = false;



						//MODIFIED: for first or last exon, check for in-frame stop codon, if there is, then we must

						//revise start_site (or end_site?)!

/*						if (all_prev_acceptors.empty()) //if it's the first exon

						{

							bool possible_valid_exons;

							if (!FirstExonIsValid((*chrMapIt).second, (*result_it).second.first, (*groupMapIt).second, 

								start_site_hsp_rel_index, possible_valid_exons))

							{

								if (possible_valid_exons)

									goto SELECT_DONOR_ACCEPTOR6;

								else

									goto CLEAN_UP6;

							}

						}

*/

						

						//if good_exon, then check if it gives a good intron

//						if (((*result_it).second.second - (*result_it).second.first - 1) < MIN_INTRON_LEN)

//							good_intron = false;



						if (good_exon) 

						{

							if (count_result == 0)

							{

								//once there's a good result, we push all prev_intron_ends to store them, and just this once

								for (int intron_end_count = 0; intron_end_count < cur_prev_intron_end.size(); intron_end_count++)

									all_prev_acceptors.push_back(cur_prev_intron_end[intron_end_count]);

								prev_intron_end.clear();

								all_prev_donors_tail.clear();



								//all_prev_donors_tail.push_back((*result_it).second.donor_tail); //may not need it

								//prev_intron_end = (*result_it).second.second;								

								prev_intron_start_donor_index.push_back( last_used_donor_index );

							}



							//all_prev_acceptors.push_back(cur_prev_intron_end);

							all_prev_donors.push_back(ExonSiteInfo((*result_it).second.first, 

								(*result_it).second.donor_query_end, (*result_it).second.donor_hsp_index, 

								(*result_it).second.donor_tail));

						

							//push all good results' intron_ends to store them, for every good result

							prev_intron_end.push_back(ExonSiteInfo((*result_it).second.second,

								(*result_it).second.acceptor_query_start, (*result_it).second.acceptor_hsp_index, 

								(3 - (*result_it).second.donor_tail) % 3));

							all_prev_donors_tail.push_back((*result_it).second.donor_tail);



							count_result++;

						}

						else//if !good_exon //now may not need this, since we've already checked for MIN_LENGTHs when computing alignments in the first place

						{

							if (count_result == 0) //if count_result>0, then skip this bad result, go on checking next

							{



							//last candidate of top-ranked intron in this region, still no good result, then maybe backtrack

							if (next_result_it == result_intron.end()) 

							{

							dataMgr.outFile << "bad result: exon below length limit, backtrack and try again" << "\n";

							cout << "bad result? exon below length limit" << "\n";

							exit(-1);



							}

							}

						}



					result_it++;

					next_result_it++;
					c++;

				}

				delete [] op_align_cur;
			}

			

			//else //there's no matching donor-acceptor pair, then reconsider these donors in the next round!

			if (count_result == 0)

			{

				//change the last bunch of donors to next region

				for (int di = last_used_donor_index; di < i; di++)

					donor_region_index[di]++;



				//reset the iterator i back to last_used_donor_index

				i = last_used_donor_index;

			}

			else

			{

				count_alternative_splice_pairs.push_back(count_result);

			}



		}



		//output
		//first, push in start_site
		if (start_site > 0)

			all_alternative_acceptors.push_back( vector< ExonSiteInfo >(1, ExonSiteInfo(start_site, 

				startHSP->gene_start, HSPs.size()-1, start_frame)) );

		else

			all_alternative_acceptors.push_back( vector< ExonSiteInfo >(1, ExonSiteInfo(start_site, 

				startHSP->gene_start, 0, start_frame)) );

		all_alternative_donors.push_back( vector< ExonSiteInfo >() );

		int j, k, cur_size;

		ai = 0; di=0;

		//UPDATE: always check the last exon for stop codon etc.
		//UPDATE: EXCEPT the last stop codon! Now there will always be a stop codon at the end since stop is included!
		if (!CheckSglExon(HSPs, chr_seq, prev_intron_end.back().frame, end_frame, prev_intron_end.back().first, 
			all_alternative_acceptors, all_alternative_donors, has_possible_exon, repair_only))
			return false; //if CheckSglExon returned true, then go on to next statement
		
		if (all_prev_acceptors.empty()) //has predicted only one exon (start_site to end_site)
		{
			//if in repair mode, make sure the frame of start and end matches
			if (repair_only && (end_site - start_site + 1) % 3 != 0) {
				has_possible_exon = false;
				all_alternative_acceptors.clear();
				all_alternative_donors.clear();
				return false;
			}
		}
		else //has predicted more than one exons
		{
			all_prev_acceptors.erase(all_prev_acceptors.begin()); //the first one must be start_site, and it's the only one that's not a pair

			//add prev_intron_ends to all_prev_acceptors, so they're all paired up

			int count_last_set_of_accs = prev_intron_end.size();
			if (repair_only)
			{
				int count_all_set_of_dons = all_prev_donors.size();
				int count_erased = 0;
				for (i=0; i<count_last_set_of_accs; i++)
				{
					if ( (end_site - prev_intron_end[i].first + 1 - prev_intron_end[i].frame - end_frame) % 3 == 0)
						all_prev_acceptors.push_back(prev_intron_end[i]);
					else //does not match end_site
					{
						all_prev_donors.erase(all_prev_donors.begin()+count_all_set_of_dons-count_last_set_of_accs
							+i-count_erased);
						count_alternative_splice_pairs.back()--;
						count_erased++;
					}
				}

				if (count_alternative_splice_pairs.back() == 0)
				{
					has_possible_exon = false;
					all_alternative_acceptors.clear();
					all_alternative_donors.clear();
					return false;
				}
			}
			else
			{
				for (i=0; i<count_last_set_of_accs; i++)
					all_prev_acceptors.push_back(prev_intron_end[i]);
			}


		for (i=0; i<count_alternative_splice_pairs.size(); i++)

		{

			int count_result = count_alternative_splice_pairs[i];

			cur_size = all_alternative_acceptors.size();

			for (j=1; j<count_result; j++)

			{

				for (k=0; k<cur_size; k++)

				{

					all_alternative_acceptors.push_back( all_alternative_acceptors[k] );

					all_alternative_donors.push_back(all_alternative_donors[k]);

				}

			}

			int inc = 0;				

			for (j=0; j<count_result; j++)

			{

				for (k=0; k<cur_size; k++)

				{

					all_alternative_donors[inc].push_back(all_prev_donors[ai]);

					all_alternative_acceptors[inc].push_back(all_prev_acceptors[ai]);

					inc++;

				}

				ai++;

			}

		}

		}



		//add end_site to all_alternative_donors

		cur_size = all_alternative_donors.size();

		if (start_site > 0)

		{

			for (j=0; j<cur_size; j++)

				all_alternative_donors[j].push_back( ExonSiteInfo(end_site, endHSP->gene_end, 0, 0) ); //end_site's frame doesn't matter?

		}

		else

		{

			for (j=0; j<cur_size; j++)

				all_alternative_donors[j].push_back( ExonSiteInfo(end_site, endHSP->gene_end, HSPs.size()-1, 0) );

		}

		

			

#ifdef DEBUG

		for (i=0; i<all_prev_acceptors.size(); i++)

			dataMgr.outFile << "all_prev_acceptors: " << all_prev_acceptors[i].first << "(" << all_prev_acceptors[i].second << ")" << "\n";

		for (i=0; i<all_prev_donors.size(); i++)

			dataMgr.outFile << "all_prev_donors: " << all_prev_donors[i].first << "(" << all_prev_donors[i].second << ")" << "\n";

		for (i=0; i<prev_intron_end.size(); i++)

			dataMgr.outFile << "prev_intron_end: " << prev_intron_end[i].first << "(" << prev_intron_end[i].second << ")" << "\n";



		for (i=0; i<all_alternative_acceptors.size(); i++)

			for (j=0; j<all_alternative_acceptors[i].size(); j++)

				dataMgr.outFile << i << " acceptor: " << all_alternative_acceptors[i][j].first << "\n";

		for (i=0; i<all_alternative_donors.size(); i++)

			for (j=0; j<all_alternative_donors[i].size(); j++)

				dataMgr.outFile << i << " donor: " << all_alternative_donors[i][j].first << "\n";

#endif



	return true;

}




bool ACCP_DONR_graph::RepairMidExon(//int prev_exon_start, int prev_exon_end, int prev_query_start, int prev_query_end, 
									//int cur_exon_start, int cur_exon_end, int cur_query_start, int cur_query_end, 
									ExonSiteInfo prev_exon_start, ExonSiteInfo prev_exon_end,
									ExonSiteInfo cur_exon_start, ExonSiteInfo cur_exon_end, 
									 //pair<int, char*>& chr_seq, 
									 vector<string>& chr_seq, 
									 int& max_alignment_HSP_ID, vector<HSP_Gene_Pair*>& newHSP_ptrs, 
									 vector<HSP_Gene_Pair*>& blastHSPs, 
									 //int prev_start_blast_hsp_index, int prev_end_blast_hsp_index, 
									 //int cur_start_blast_hsp_index, int cur_end_blast_hsp_index,
									 vector<ExonSiteInfo>& temp_sites,  //tmp_sites keep all unconfirmed sites?
									 vector<HSP_Gene_Pair>& temp_HSPs, 
									 vector<ExonSiteInfo>& additional_temp_sites, 
									 bool is_last_exon) //vector<HSP_Gene_Pair*>& temp_HSPs)
{
	int j;

#ifdef DEBUG
	dataMgr.outFile << "all old HSPs:" << "\n";
	PrintHSPs(blastHSPs);
	dataMgr.outFile << "RepairMidExon: start with: " << "\n";
	for (j=0; j<temp_HSPs.size(); j++)
		dataMgr.outFile << (temp_HSPs[j]) << "\n";
	for (j=0; j<temp_sites.size(); j++)
		dataMgr.outFile << temp_sites[j] << "\n";
#endif

	string target_seq, q_seq;



	//vector<HSP_Gene_Pair*> newHSPs; //used to store everything

	vector<HSP_Gene_Pair> newHSPs; //used to store everything, do a copy of HSPs, do NOT change the old BlastHSPs!

	vector<HSP_Gene_Pair*> ptr_to_newHSPs;



	//align_pos[frame][]: all relative indexes!

	//1: start pos on query_seq; 2: end pos on query_seq; 3: start pos on target; 4: end pos on target

	int align_pos[3][4]; 

	int prev_exon_start_frame = prev_exon_start.frame;
	int cur_exon_end_frame = cur_exon_end.frame;
	//int prev_start_blast_hsp_index = prev_exon_start.hsp_index;
	//int prev_end_blast_hsp_index = prev_exon_end.hsp_index;
	int cur_start_blast_hsp_index = cur_exon_start.hsp_index;
	int cur_end_blast_hsp_index = cur_exon_end.hsp_index;

/*	//UPDATE: now check the second exon, if it has in-frame stop, then remove the HSP before that stop and do it again
	string exonStr;
	if (start_site > 0)
		//GetSubstrFromVecStrs(chr_seq, true, start_site-1, end_site-start_site+1-end_frame, testStr);
		GetSubstrFromVecStrs(chr_seq, true, cur_exon_start.first-1+cur_exon_start.frame, cur_exon_end.first-cur_exon_start.first+1-cur_exon_end_frame-cur_exon_start.frame, exonStr);
	else
		//GetSubstrFromVecStrs_NegRev(chr_seq, false, -end_site-1+end_frame, end_site-start_site+1-end_frame, testStr);
		GetSubstrFromVecStrs_NegRev(chr_seq, false, -cur_exon_end.first-1+cur_exon_end_frame, cur_exon_end.first-cur_exon_start.first+1-cur_exon_end_frame-cur_exon_start.frame, exonStr);
	StrToLower(testStr);
#ifdef DEBUG
	dataMgr.outFile << "cur_exon_start:" << cur_exon_start.first << ";cur_exon_end:" << cur_exon_end.first << "\n";
	dataMgr.outFile << "exonStr:" << exonStr << "\n";
#endif

	int stop_pos;
	if (HasInFrameStop(exonStr, true, stop_pos))
	{
#ifdef DEBUG
		dataMgr.outFile << "end_site has in-frame stop from exon_start_site: stop_pos:" << stop_pos << "\n";
#endif
		bool removed_hsp_is_new_hsp;
		if (RemoveHSPBeforeStop(HSPs, cur_exon_start.first+cur_exon_start.frame+stop_pos, chr_seq, removed_hsp_is_new_hsp))
		{
			if (HSPs.size() > 0) //there's still HSPs left after removal
			{
#ifdef COMPUT_EXON_FULL_STEP_BACK
	#ifdef DEBUG
				dataMgr.outFile << "removed a HSP, still HSP left, now going back one full step" << "\n";
	#endif
#else
	#ifdef DEBUG								
				dataMgr.outFile << "encounted a stop, adjusted segments and go back (not full step back)" << "\n";
	#endif
#endif
			}
			else
			{
#ifdef DEBUG
				dataMgr.outFile << "no more HSPs, no possible exon\n";
#endif
			}
		}
		else
		{
#ifdef DEBUG
			dataMgr.outFile << "did not find the HSP to erase, no possible exon\n";
#endif
		}
	}
*/

	int hsp_gap = cur_exon_start.first - prev_exon_end.first - 1; //distance between the two exons on target



			if (!(hsp_gap > 0 && cur_exon_start.second - prev_exon_end.second > MID_EXON_LEN)) //no need for repair

			{

				temp_sites.push_back(cur_exon_start);

				temp_sites.push_back(cur_exon_end);



				additional_temp_sites.push_back(cur_exon_start);

				additional_temp_sites.push_back(cur_exon_end);



					temp_HSPs.clear(); //now update temp_HSPs

					if (prev_exon_start.first > 0)

					{

						for (j=cur_end_blast_hsp_index; j<=cur_start_blast_hsp_index; j++)

							temp_HSPs.push_back(*(blastHSPs[j]));

					}

					else

					{

						for (j=cur_start_blast_hsp_index; j<=cur_end_blast_hsp_index; j++)

							temp_HSPs.push_back(*(blastHSPs[j]));

					}

#ifdef DEBUG

	dataMgr.outFile << "all old HSPs:" << "\n";

	PrintHSPs(blastHSPs);

					dataMgr.outFile << "no need to repair, use old exon: temp_HSPs for this exon: " << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

					return false; //no need to pop

			}

			else //need repair

			{

				double max_score=0, cur_score;

				float max_pid = 0, cur_pid;



				if (prev_exon_start.first < 0)

					GetSubstrFromVecStrs_NegRev(chr_seq, false, -cur_exon_start.first, hsp_gap, target_seq);

				else

					GetSubstrFromVecStrs(chr_seq, true, prev_exon_end.first, hsp_gap, target_seq);

				StrToLower(target_seq);



				HSP_Gene_Pair newHSP[3];

				HSP_Gene_Pair* curHSP = new HSP_Gene_Pair;

				Input_Alignment newAlign[3];

				int j=0;

				string result_seq[3];

				int query_gap = cur_exon_start.second - prev_exon_end.second - 1;

				q_seq = query_seq.substr(prev_exon_end.second, query_gap);

				for (int i=0; i<3; i++)

				{

					DNA2AA(target_seq, i, result_seq[i]);

#ifdef DEBUG

					dataMgr.outFile << "target_seq:" << "\n" << target_seq << "\n";

					dataMgr.outFile << result_seq[i] << "\n";

#endif



/*					if (prev_exon_start < 0) //same as positive?

						cur_score = GetLocalAlignment(q_seq, result_seq[i], 

							prev_query_end+1, -(-prev_exon_end-1-i), max_alignment_HSP_ID, newHSP[i], newAlign[i],

							true, dataMgr.outFile, align_pos[i]); //close_to_start or end, doesn't matter?

					else

*/						cur_score = GetLocalAlignment(q_seq, result_seq[i], 

							prev_exon_end.second+1, prev_exon_end.first+1+i, max_alignment_HSP_ID, newHSP[i], newAlign[i], 

							true, dataMgr.outFile, align_pos[i]); //close_to_start or end, doesn't matter?

					cur_pid = newHSP[i].pid;



					if (cur_score == 0)

						continue;

					

if (EXTEND_HSP_BY_SCORE)

{

					if (max_score < cur_score )

					{

						j = i;

						max_score = cur_score;

						max_pid = cur_pid;

						*curHSP = newHSP[i];

					}

}

else

{

					if (max_pid < cur_pid || (max_pid == cur_pid && max_score < cur_score))

					{

						j = i;

						max_score = cur_score;

						max_pid = cur_pid;

						*curHSP = newHSP[i];

					}

}

				}



				//dataMgr.HSP_neg_gene[chr_index].push_back(newHSP[j]);

				//newHSPs.push_back(&(dataMgr.HSP_neg_gene[chr_index].back()));

				if (max_score == 0 || max_score < REPAIR_HSP_MIN_INIT_SCORE) //<= 0)

				{

						temp_sites.push_back(cur_exon_start);

						temp_sites.push_back(cur_exon_end);



						additional_temp_sites.push_back(cur_exon_start);

						additional_temp_sites.push_back(cur_exon_end);



						temp_HSPs.clear();

						if (prev_exon_start.first > 0)

						{

							for (j=cur_end_blast_hsp_index; j<=cur_start_blast_hsp_index; j++)

								temp_HSPs.push_back(*(blastHSPs[j]));

						}

						else

						{

							for (j=cur_start_blast_hsp_index; j<=cur_end_blast_hsp_index; j++)

								temp_HSPs.push_back(*(blastHSPs[j]));

						}

#ifdef DEBUG

	dataMgr.outFile << "all old HSPs:" << "\n";

	PrintHSPs(blastHSPs);

					dataMgr.outFile << "loaddata return false, use old exon: temp_HSPs for this exon: " << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

						return false; //no valid extra alignment, just use the old exon



				}

				else

				{

					//if the score of init alignment is big enough, try extend this alignment to both directions, just like BLAST!

					//extend curHSP to both directions if possible, stop extending when the score of alignment is 

					//lower than ("max_score" - REPAIR_HSP_EXTEND_SCORE_DROP), "max_score" is the running max score

					//if (max_score >= REPAIR_HSP_MIN_INIT_SCORE)
					//{

						string tmpStr = "";

						if (align_pos[j][3]+1 < result_seq[j].length())

							tmpStr = result_seq[j].substr(align_pos[j][3]+1);


						string result_begin_str = result_seq[j].substr(0, align_pos[j][2]);

						if (prev_exon_start.first < 0)

						ExtendAlignment(max_score, query_seq, 

							//result_seq[j].substr(0, align_pos[j][2]), //this may also be empty string
							result_begin_str, 

							tmpStr, 

							chr_seq, false, curHSP, newAlign[j], dataMgr.outFile, 

							//(*hsp_it)->HSP_end, (*next_hsp_it)->HSP_start);

							//-prev_exon_start, -cur_exon_end);

							temp_HSPs.back().HSP_end, blastHSPs[cur_start_blast_hsp_index]->HSP_start);

						else

						ExtendAlignment(max_score, query_seq, 

							//result_seq[j].substr(0, align_pos[j][2]), //this may also be empty string
							result_begin_str, 

							tmpStr, 

							chr_seq, true, curHSP, newAlign[j], dataMgr.outFile, 

							//(*next_hsp_it)->HSP_start, (*hsp_it)->HSP_end);

							//prev_exon_start, cur_exon_end);

							temp_HSPs.front().HSP_start, blastHSPs[cur_start_blast_hsp_index]->HSP_end);



#ifdef DEBUG

						dataMgr.outFile << "extended alignment: " << "\n";

						dataMgr.outFile << *curHSP << "\n";

						dataMgr.outFile << newAlign[j] << "\n";

#endif

					//}



				//before we accept curHSP, still need to check it against one previous BlastHSP and one next BlastHSP, 

				//to make sure it is valid (doesn't entirely include or be included by other HSP)

				HSP_Gene_Pair *prev_hsp, *next_hsp;

				if (prev_exon_start.first < 0)

				{

					prev_hsp = &(temp_HSPs.back());

					//next_hsp = blastHSPs[cur_start_blast_hsp_index];

				}

				else

				{

					prev_hsp = &(temp_HSPs.front());

					//next_hsp = blastHSPs[cur_start_blast_hsp_index];

				}

				next_hsp = blastHSPs[cur_start_blast_hsp_index]; //same for both positive and negative

#ifdef DEBUG

				dataMgr.outFile << "prev_hsp:" << *prev_hsp << "\n";

				dataMgr.outFile << "next_hsp:" << *next_hsp << "\n";

#endif

				//MODIFIED: check if prev_hsp is the same as next_hsp

				bool fit_with_neighbors = true;

				if (!(curHSP->FitWithNeighborHSPs(prev_hsp, start_site>0)))

					fit_with_neighbors = false;

				else

				{

					if ((next_hsp->ID != prev_hsp->ID) && !(curHSP->FitWithNeighborHSPs(next_hsp, start_site>0)))

						fit_with_neighbors = false;

				}


				if (!fit_with_neighbors)
				{

					delete curHSP;

					temp_sites.push_back(cur_exon_start);

					temp_sites.push_back(cur_exon_end);



					additional_temp_sites.push_back(cur_exon_start);

					additional_temp_sites.push_back(cur_exon_end);



					temp_HSPs.clear();

					if (prev_exon_start.first > 0)

					{

						for (j=cur_end_blast_hsp_index; j<=cur_start_blast_hsp_index; j++)

							temp_HSPs.push_back(*(blastHSPs[j]));

					}

					else

					{

						for (j=cur_start_blast_hsp_index; j<=cur_end_blast_hsp_index; j++)

							temp_HSPs.push_back(*(blastHSPs[j]));

					}

#ifdef DEBUG

	dataMgr.outFile << "all old HSPs:" << "\n";

	PrintHSPs(blastHSPs);

					dataMgr.outFile << "new HSP doesn't fit, use old exon: temp_HSPs for this exon: " << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

					return false;

				}



				newHSPs.push_back(*curHSP);

				newHSP_ptrs.push_back(curHSP);

				dataMgr.input_alignments.insert(map<int, Input_Alignment>::value_type(max_alignment_HSP_ID, newAlign[j]));
				hspID_added_for_repair = max_alignment_HSP_ID; //set hspID_added_for_repair to the new HSP
#ifdef DEBUG
				dataMgr.outFile << "hspID_added_for_repair is: HSP[" << hspID_added_for_repair << "]\n";
#endif
				max_alignment_HSP_ID++;



				//now redo PreProcHSPs and LoadData4 etc.

				set<int> collected_hsp_IDs; //used to collect the old HSPs (IDs) that have been added to newHSPs, so we don't add the same HSP multiple times!
				if (prev_exon_start.first < 0) //negative

				{

					//for (j=prev_end_blast_hsp_index; j>=prev_start_blast_hsp_index; j--)

					//	newHSPs.insert(newHSPs.begin(), blastHSPs[j]);

					for (j=temp_HSPs.size()-1; j>=0; j--)

					{

#ifdef DEBUG

						dataMgr.outFile << "before insertion: " << "\n";

						PrintHSPs(newHSPs);

#endif

						if (collected_hsp_IDs.find(temp_HSPs[j].ID) != collected_hsp_IDs.end())

							continue;

						else

							collected_hsp_IDs.insert(temp_HSPs[j].ID);

						newHSPs.insert(newHSPs.begin(), temp_HSPs[j]);

						newHSPs.front().ID = max_alignment_HSP_ID; //also replicate the alignment so it won't be changed on the original blast alignment

						dataMgr.input_alignments.insert(map<int, Input_Alignment>::value_type(max_alignment_HSP_ID, 

							(*(dataMgr.input_alignments.find(temp_HSPs[j].ID))).second));

						max_alignment_HSP_ID++;
#ifdef DEBUG

						dataMgr.outFile << "after insertion: " << "\n";

						PrintHSPs(newHSPs);

#endif


					}

					for (j=cur_start_blast_hsp_index; j<=cur_end_blast_hsp_index; j++)

					{
#ifdef DEBUG

						dataMgr.outFile << "before extension: " << "\n";

						PrintHSPs(newHSPs);

#endif

						if (collected_hsp_IDs.find(blastHSPs[j]->ID) != collected_hsp_IDs.end())

							continue;

						else

							collected_hsp_IDs.insert(blastHSPs[j]->ID);


						newHSPs.push_back(*(blastHSPs[j]));

						newHSPs.back().ID = max_alignment_HSP_ID; //also replicate the alignment so it won't be changed on the original blast alignment

						dataMgr.input_alignments.insert(map<int, Input_Alignment>::value_type(max_alignment_HSP_ID, 

							(*(dataMgr.input_alignments.find(blastHSPs[j]->ID))).second));

						max_alignment_HSP_ID++;

#ifdef DEBUG

						dataMgr.outFile << "after extension: " << "\n";

						PrintHSPs(newHSPs);

#endif


					}

				}

				else //positive

				{

					for (j=cur_start_blast_hsp_index; j>=cur_end_blast_hsp_index; j--)

					{
#ifdef DEBUG

						dataMgr.outFile << "before insertion: " << "\n";

						PrintHSPs(newHSPs);

#endif

						if (collected_hsp_IDs.find(blastHSPs[j]->ID) != collected_hsp_IDs.end())

							continue;

						else

							collected_hsp_IDs.insert(blastHSPs[j]->ID);


						newHSPs.insert(newHSPs.begin(), *(blastHSPs[j]));

						newHSPs.front().ID = max_alignment_HSP_ID; //also replicate the alignment so it won't be changed on the original blast alignment

						dataMgr.input_alignments.insert(map<int, Input_Alignment>::value_type(max_alignment_HSP_ID, 

							(*(dataMgr.input_alignments.find(blastHSPs[j]->ID))).second));

						max_alignment_HSP_ID++;

#ifdef DEBUG

						dataMgr.outFile << "after insertion: " << "\n";

						PrintHSPs(newHSPs);

#endif

					}

					//for (j=prev_end_blast_hsp_index; j<=prev_start_blast_hsp_index; j++)

					//	newHSPs.push_back(blastHSPs[j]);

					for (j=0; j<temp_HSPs.size(); j++)

					{
#ifdef DEBUG

						dataMgr.outFile << "before extension: " << "\n";

						PrintHSPs(newHSPs);

#endif

						if (collected_hsp_IDs.find(temp_HSPs[j].ID) != collected_hsp_IDs.end())

							continue;

						else

							collected_hsp_IDs.insert(temp_HSPs[j].ID);



						newHSPs.push_back(temp_HSPs[j]);

						newHSPs.back().ID = max_alignment_HSP_ID; //also replicate the alignment so it won't be changed on the original blast alignment

						dataMgr.input_alignments.insert(map<int, Input_Alignment>::value_type(max_alignment_HSP_ID, 

							(*(dataMgr.input_alignments.find(temp_HSPs[j].ID))).second));

						max_alignment_HSP_ID++;

#ifdef DEBUG

						dataMgr.outFile << "after extension: " << "\n";

						PrintHSPs(newHSPs);

#endif

					}

				}

#ifdef DEBUG

	dataMgr.outFile << "all old HSPs:" << "\n";

	PrintHSPs(blastHSPs);

					dataMgr.outFile << "newHSPs for repairing: " << "\n";

					for (j=0; j<newHSPs.size(); j++)

						dataMgr.outFile << (newHSPs[j]) << "\n";

#endif



				for (j=0; j<newHSPs.size(); j++)

					ptr_to_newHSPs.push_back(&(newHSPs[j])); //record the address for these new HSPs

				PreProcHSPs((prev_exon_start.first>0), ptr_to_newHSPs, newHSP_ptrs, max_alignment_HSP_ID);

					//need to update newHSPs? since PreProcHSPs may have inserted new HSPs...nah, just use ptr_to_newHSPs!

/*					int k=0;

					for (j=0; j<ptr_to_newHSPs.size(); j++) 

					{

						if (ptr_to_newHSPs[j]->ID != newHSPs[k].ID)

							newHSPs.insert(newHSPs.begin()+k, *(ptr_to_newHSPs[j]));

						k++;

					}

*/



#ifdef DEBUG

	dataMgr.outFile << "all old HSPs:" << "\n";

	PrintHSPs(blastHSPs);

					dataMgr.outFile << "newHSPs after preproc: " << "\n";

					PrintHSPs(ptr_to_newHSPs);

					//for (j=0; j<ptr_to_newHSPs.size(); j++)

					//	dataMgr.outFile << *(ptr_to_newHSPs[j]) << "\n";

#endif



					end_site = cur_exon_end.first; //fix both start_site and end_site first

					start_site = prev_exon_start.first;

SELECT_DONOR_ACCEPTOR6_MIDEXON:

					if (SPLICE_SEGMENT_VERSION == 1)

					{

						HSPs_dup.clear();

						dataMgr.input_alignments_HSPs_dup.clear();

					}

					if (!LoadData4_forGenePredict6(ptr_to_newHSPs, chr_seq, true, true))

					{

						temp_sites.push_back(cur_exon_start);

						temp_sites.push_back(cur_exon_end);



						additional_temp_sites.push_back(cur_exon_start);

						additional_temp_sites.push_back(cur_exon_end);



						temp_HSPs.clear();

						if (prev_exon_start.first > 0)

						{

							for (j=cur_end_blast_hsp_index; j<=cur_start_blast_hsp_index; j++)

								temp_HSPs.push_back(*(blastHSPs[j]));

						}

						else

						{

							for (j=cur_start_blast_hsp_index; j<=cur_end_blast_hsp_index; j++)

								temp_HSPs.push_back(*(blastHSPs[j]));

						}

#ifdef DEBUG

	dataMgr.outFile << "all old HSPs:" << "\n";

	PrintHSPs(blastHSPs);

					dataMgr.outFile << "loaddata return false, use old exon: temp_HSPs for this exon: " << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

						return false; //something is wrong, just use the old exon

					}

#ifdef DEBUG

	dataMgr.outFile << "after LoadData, all old HSPs:" << "\n";

	PrintHSPs(blastHSPs);

	dataMgr.outFile << "newHSPs:" << "\n";

	PrintHSPs(ptr_to_newHSPs);

#endif

					vector< vector<ExonSiteInfo> > all_alternative_acceptors, all_alternative_donors;

					bool has_possible_exon;
#ifdef COMPUT_EXON_FULL_STEP_BACK

					if (!ComputeExons(ptr_to_newHSPs, chr_seq, prev_exon_start_frame, cur_exon_end_frame, 
						all_alternative_acceptors, all_alternative_donors, has_possible_exon, true))
						if (has_possible_exon)
							goto SELECT_DONOR_ACCEPTOR6_MIDEXON;

#else

					while (!ComputeExons(ptr_to_newHSPs, chr_seq, prev_exon_start_frame, cur_exon_end_frame, 
						all_alternative_acceptors, all_alternative_donors, has_possible_exon, true))
						if (!has_possible_exon)
							break;

#endif



#ifdef DEBUG

	dataMgr.outFile << "all old HSPs:" << "\n";

	PrintHSPs(blastHSPs);

					dataMgr.outFile << "new exons computed(" << has_possible_exon << ")" << "\n";

#endif
					if (!has_possible_exon)
					{
						UseOldExon_Mid(temp_sites, additional_temp_sites, temp_HSPs, 
								 cur_exon_start, cur_exon_end, prev_exon_start, 
								 cur_start_blast_hsp_index, cur_end_blast_hsp_index, 
								 blastHSPs);
						return false;
					}



					//multiply all_alternative_acceptors/donors into the old all_alternative_acceptors/donors?

					int frame = prev_exon_start_frame;

					string final_alignment, left_chars;

					left_chars.assign((3-frame)%3, 'x');

					int prev_donor_site = 0;

					//temp_sites.clear();

					int site_before_stop;

					for (j=0; j<all_alternative_acceptors[0].size(); j++) //use the first one, for now

					{

						string tmp_str;
						bool found_stop;
						if (HasStopCodon(all_alternative_acceptors[0][j].first, all_alternative_donors[0][j].first, 

							chr_seq, frame, site_before_stop, final_alignment, prev_donor_site, left_chars, 
							false, dataMgr.cDNA_os, false, tmp_str, is_last_exon && j==all_alternative_acceptors[0].size()-1, found_stop))

						{ //if there's stop in new exons, discard them
							UseOldExon_Mid(temp_sites, additional_temp_sites, temp_HSPs, 
									 cur_exon_start, cur_exon_end, prev_exon_start, 
									 cur_start_blast_hsp_index, cur_end_blast_hsp_index, 
									 blastHSPs);

#ifdef DEBUG

	dataMgr.outFile << "all old HSPs:" << "\n";

	PrintHSPs(blastHSPs);

					dataMgr.outFile << "new exon has stop, use old exon: temp_HSPs for this exon: " << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << temp_HSPs[j] << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

							return false;

						}

						else

						{

							prev_donor_site = all_alternative_donors[0][j].first;

						}

					}

					Input_Alignment newAlign;

					float align_pid;

					string query_mid_str = query_seq.substr(prev_exon_start.second-1, cur_exon_end.second - prev_exon_start.second + 1);
					int align_score;
#ifdef PID_BEFORE_SCORE
					align_score = GetGlobalAlignment_PID_CompScore(query_mid_str, final_alignment, newAlign, align_pid, dataMgr.outFile);
#else
					align_score = GetGlobalAlignment_scorematrix(//query_seq.substr(prev_exon_start.second-1, cur_exon_end.second - prev_exon_start.second + 1), 
						query_mid_str, final_alignment, newAlign, align_pid, dataMgr.outFile);
#endif
					align_compare align_new(align_score, align_pid);

					//dataMgr.outFile << newAlign << "\n";

					//dataMgr.outFile << "PID:" << align_pid << "\n";

					

					frame = prev_exon_start_frame;

					final_alignment = "";

					left_chars.assign((3-frame)%3, 'x');

					prev_donor_site = 0;

					int num_of_new_exons = all_alternative_acceptors[0].size();

#ifdef DEBUG

					dataMgr.outFile << "num_of_new_exons:" << num_of_new_exons << "\n";

#endif

					vector<ExonSiteInfo> old_exon_sites;

					old_exon_sites.push_back(prev_exon_start);

					old_exon_sites.push_back(prev_exon_end);

					old_exon_sites.push_back(cur_exon_start);

					old_exon_sites.push_back(cur_exon_end);

					for (j=0; j<4; j+=2)

					{

						string tmp_str;
						bool found_stop;
						if (HasStopCodon(old_exon_sites[j].first, old_exon_sites[j+1].first, 

							chr_seq, frame, site_before_stop, final_alignment, prev_donor_site, left_chars, 
							false, dataMgr.cDNA_os, false, tmp_str, is_last_exon, found_stop))

						{

						temp_sites.pop_back(); //pop out one previous exon

						temp_sites.pop_back();

						for (j=0; j<num_of_new_exons; j++) //use the first one, for now

						{

							temp_sites.push_back(all_alternative_acceptors[0][j]);

							temp_sites.push_back(all_alternative_donors[0][j]);



							additional_temp_sites.push_back(all_alternative_acceptors[0][j]);

							additional_temp_sites.push_back(all_alternative_donors[0][j]);

							//output some additional info, for inspection purpose only
							if (all_alternative_acceptors[0][j].first > prev_exon_end.first && all_alternative_donors[0][j].first < cur_exon_start.first )
								dataMgr.outFile << "NEW EXON after repair:" << all_alternative_acceptors[0][j].first << "-" << all_alternative_donors[0][j].first << "\n";

						}

						temp_HSPs.clear();

						if (prev_exon_start.first > 0)

						{

							for (j=all_alternative_donors[0][num_of_new_exons-1].hsp_index; 

								j<=all_alternative_acceptors[0][num_of_new_exons-1].hsp_index; j++)

									if (ptr_to_newHSPs[j]->HSPInRegion(all_alternative_acceptors[0][num_of_new_exons-1].first, 

										all_alternative_donors[0][num_of_new_exons-1].first))

										temp_HSPs.push_back(*(ptr_to_newHSPs[j]));

						}

						else

						{

							for (j=all_alternative_acceptors[0][num_of_new_exons-1].hsp_index; 

								j<=all_alternative_donors[0][num_of_new_exons-1].hsp_index; j++)

									if (ptr_to_newHSPs[j]->HSPInRegion(-all_alternative_donors[0][num_of_new_exons-1].first, 

										-all_alternative_acceptors[0][num_of_new_exons-1].first))

										temp_HSPs.push_back(*(ptr_to_newHSPs[j]));

						}

#ifdef DEBUG

	dataMgr.outFile << "all old HSPs:" << "\n";

	PrintHSPs(blastHSPs);

					dataMgr.outFile << "old exon has stop, use new exon: temp_HSPs for this exon: " << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

						return true;

						}

						else //has no stop codon

						{

							prev_donor_site = old_exon_sites[j+1].first;

						}

					} //end for (old_exon_sites checking and got target sequence portion)



						float ori_align_pid;
						int ori_align_score;
#ifdef PID_BEFORE_SCORE
						ori_align_score = GetGlobalAlignment_PID_CompScore(query_mid_str, final_alignment, newAlign, ori_align_pid, dataMgr.outFile);
#else
						ori_align_score = GetGlobalAlignment_scorematrix(//query_seq.substr(prev_exon_start.second-1, cur_exon_end.second - prev_exon_start.second + 1), 
							query_mid_str, final_alignment, newAlign, ori_align_pid, dataMgr.outFile);
#endif
						align_compare align_old(ori_align_score, ori_align_pid);

						//if (ori_align_score < align_score || (ori_align_score == align_score && ori_align_pid < align_pid))
						if (align_new > align_old)
						{

							temp_sites.pop_back();

							temp_sites.pop_back();

							for (j=0; j<num_of_new_exons; j++) //use the first one, for now

							{

								temp_sites.push_back(all_alternative_acceptors[0][j]);

								temp_sites.push_back(all_alternative_donors[0][j]);



								additional_temp_sites.push_back(all_alternative_acceptors[0][j]);

								additional_temp_sites.push_back(all_alternative_donors[0][j]);

							//output some additional info, for inspection purpose only
							if (all_alternative_acceptors[0][j].first > prev_exon_end.first && all_alternative_donors[0][j].first < cur_exon_start.first )
								dataMgr.outFile << "NEW EXON after repair:" << all_alternative_acceptors[0][j].first << "-" << all_alternative_donors[0][j].first << "\n";
							}

			

							temp_HSPs.clear();

							if (prev_exon_start.first > 0)

							{

								for (j=all_alternative_donors[0][num_of_new_exons-1].hsp_index; 

									j<=all_alternative_acceptors[0][num_of_new_exons-1].hsp_index; j++)

									if (ptr_to_newHSPs[j]->HSPInRegion(all_alternative_acceptors[0][num_of_new_exons-1].first, 

										all_alternative_donors[0][num_of_new_exons-1].first))

										temp_HSPs.push_back(*(ptr_to_newHSPs[j]));

							}

							else

							{

								for (j=all_alternative_acceptors[0][num_of_new_exons-1].hsp_index; 

									j<=all_alternative_donors[0][num_of_new_exons-1].hsp_index; j++)

									if (ptr_to_newHSPs[j]->HSPInRegion(-all_alternative_donors[0][num_of_new_exons-1].first, 

										-all_alternative_acceptors[0][num_of_new_exons-1].first))

										temp_HSPs.push_back(*(ptr_to_newHSPs[j]));

							}

#ifdef DEBUG

	dataMgr.outFile << "all old HSPs:" << "\n";

	PrintHSPs(blastHSPs);

					dataMgr.outFile << "both have no stop, use new exon: temp_HSPs for this exon: " << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

					return true;

						}

						else

						{

						UseOldExon_Mid(temp_sites, additional_temp_sites, temp_HSPs, 
								 cur_exon_start, cur_exon_end, prev_exon_start, 
								 cur_start_blast_hsp_index, cur_end_blast_hsp_index, 
								 blastHSPs);
#ifdef DEBUG

	dataMgr.outFile << "all old HSPs:" << "\n";

	PrintHSPs(blastHSPs);

					dataMgr.outFile << "both have no stop, use old exon: temp_HSPs for this exon: " << "\n";

					for (j=0; j<temp_HSPs.size(); j++)

						dataMgr.outFile << (temp_HSPs[j]) << "\n";

					for (j=0; j<temp_sites.size(); j++)

						dataMgr.outFile << temp_sites[j] << "\n";

#endif

					return false;

						}

				}

			}



}



void ACCP_DONR_graph::GetSpliceSegments1_forGenePredict6(//multimap<int, pair<int, int> >& donor_segments, 

										 vector<HSP_Gene_Pair*>& HSPs, 

										 int& first_segment_start, int& last_segment_start, int& last_segment_end)

{

	//duplicate HSPs and their alignment first (to keep an original copy of these)

	vector<HSP_Gene_Pair*>::iterator HSPs_ori_it;

	if (HSPs_dup.empty()) //record HSPs in HSPs_dup in the first run only

	//if (REPAIR_HSP_AFTER_EXON == 1 || HSPs_dup.empty())

	{

		HSPs_dup.clear();

	for (HSPs_ori_it = HSPs.begin(); HSPs_ori_it != HSPs.end(); HSPs_ori_it++)

	{

		HSPs_dup.push_back(*(*HSPs_ori_it));

		dataMgr.input_alignments_HSPs_dup.insert(map<int, Input_Alignment>::value_type((*HSPs_ori_it)->ID, 

			(*(dataMgr.input_alignments.find((*HSPs_ori_it)->ID))).second));



#ifdef DEBUG

		dataMgr.outFile << "HSPs_dup recorded:" << "\n";

		dataMgr.outFile << *(*HSPs_ori_it) << "\n";

		dataMgr.outFile << (*(dataMgr.input_alignments_HSPs_dup.find((*HSPs_ori_it)->ID))).second << "\n";

#endif

	}

	}



	int cur_hsp_start, cur_hsp_end;

	bool next_exists, cur_hsp_exists;

	vector<HSP_Gene_Pair*> HSPs_to_be_erased;

	list<int> HSPs_dup_to_be_erased;

	int dummy_query_start = INT_MAX, dummy_query_end = 0, dummy_hsp_start = 0;

	int cur_query_start, cur_query_end, next_query_start, next_query_end;



	int last_segment_hsp_ID = -1; //initialize to -1 (in fact, this is the previous hsp_ID that corresponds to the current [last_segment_start, last_segment_end])



	if (start_site > 0) //positive strand HSP

	//if (isPosStrand)

	{



		vector<HSP_Gene_Pair*>::reverse_iterator vec_revIt, vec_revIt_next;

		vec_revIt_next = HSPs.rbegin(); vec_revIt_next++; //vec_revIt_next points to next HSP

		//vector<HSP_Gene_Pair>::reverse_iterator vec_revIt_next_dup = HSPs_dup.rbegin(); vec_revIt_next_dup++;

		int dup_i=HSPs_dup.size()-1, dup_next_i=HSPs_dup.size()-2;



		vec_revIt = HSPs.rbegin();

		first_segment_start = (*vec_revIt)->HSP_start;

		//vector<HSP_Gene_Pair>::reverse_iterator vec_revIt_dup = HSPs_dup.rbegin();		

		

		//for (; vec_revIt != HSPs.rend(); vec_revIt++) //for each HSP, in ascending order!

		while (vec_revIt != HSPs.rend())

		{

			//determine first segment (special treatment required)

			cur_hsp_start = (*vec_revIt)->HSP_start;

			cur_hsp_end = (*vec_revIt)->HSP_end;



			cur_query_start = (*vec_revIt)->gene_start;

			cur_query_end = (*vec_revIt)->gene_end;



			cur_hsp_exists = true;



			if (vec_revIt_next != HSPs.rend())

			{

				next_query_start = (*vec_revIt_next)->gene_start;

				next_query_end = (*vec_revIt_next)->gene_end;



				next_exists = GetSingleHSPDonorSegments_SplitOL_GapIsMinIntron(cur_hsp_start, cur_hsp_end, (*vec_revIt)->ID, 

					//donor_segments, 

					last_segment_start, last_segment_hsp_ID, 

					(*vec_revIt_next)->ID, (*vec_revIt_next)->HSP_start, 

					cur_query_start, cur_query_end, next_query_start, next_query_end, last_segment_end, cur_hsp_exists);

				

				if (next_query_start != (*vec_revIt_next)->gene_start)

					(*vec_revIt_next)->gene_start = next_query_start;

				if (cur_query_end != (*vec_revIt)->gene_end)

					(*vec_revIt)->gene_end = cur_query_end;

				if (cur_hsp_end != (*vec_revIt)->HSP_end)

					(*vec_revIt)->HSP_end = cur_hsp_end;



				//GetDonorSegments(cur_hsp_start, cur_hsp_end, prev_hsp_end, cur_start, (*vec_revIt)->ID, 

				//	donor_segments, prev_segment_end); //used for matchStr gap region

				if (!cur_hsp_exists)

				{

					HSPs_to_be_erased.push_back(*vec_revIt);

					HSPs_dup_to_be_erased.push_back(dup_i);//(*vec_revIt_dup);

				}



				if (!next_exists) //next time, use next HSP since it's still there (otherwise advance extra vec_revIt,vec_revIt_next)

				{

					HSPs_to_be_erased.push_back(*vec_revIt_next);

					HSPs_dup_to_be_erased.push_back(dup_next_i);//(*vec_revIt_next_dup);



					vec_revIt++;

					vec_revIt++;

					//vec_revIt_dup++;

					//vec_revIt_dup++;

					dup_i-=2;



					vec_revIt_next++;

					//vec_revIt_next_dup++;

					dup_next_i--;

					if (vec_revIt_next != HSPs.rend())

					{

						vec_revIt_next++;

						///vec_revIt_next_dup++;

						dup_next_i--;

					}

				}

				else

				{

					vec_revIt++;

					//vec_revIt_dup++;

					dup_i--;

					vec_revIt_next++;

					//vec_revIt_next_dup++;

					dup_next_i--;

				}

			}

			else //there's no next HSP

			{

				GetSingleHSPDonorSegments_SplitOL_GapIsMinIntron(cur_hsp_start, cur_hsp_end, (*vec_revIt)->ID, 

					//donor_segments, 

					last_segment_start, last_segment_hsp_ID, 

					-1, dummy_hsp_start, 0, dummy_query_end, dummy_query_start, INT_MAX, last_segment_end, cur_hsp_exists);

				//GetDonorSegments(cur_hsp_start, cur_hsp_end, prev_hsp_end, cur_start, (*vec_revIt)->ID, 

				//	donor_segments, prev_segment_end); //used for matchStr gap region

				

				if (!cur_hsp_exists)

				{

					HSPs_to_be_erased.push_back(*vec_revIt);

					HSPs_dup_to_be_erased.push_back(dup_i);//(*vec_revIt_dup);

				}



				vec_revIt++; //since there's no next HSP, always advance vec_revIt to stop while loop

				//vec_revIt_dup++;

				dup_i--;

				last_segment_end = cur_hsp_end;

			}

		}



	}

	else //negative HSP

	{



		vector<HSP_Gene_Pair*>::iterator vec_It, vec_It_next;

		vec_It_next = HSPs.begin(); vec_It_next++;

		//vector<HSP_Gene_Pair>::iterator vec_It_dup, vec_It_next_dup;

		//vec_It_next_dup = HSPs_dup.begin(); vec_It_next_dup++;

		int dup_i=0, dup_next_i=1;



		vec_It = HSPs.begin();

		first_segment_start = -(*vec_It)->HSP_end;

		//vec_It_dup = HSPs_dup.begin();



		//for (; vec_It != HSPs.end(); vec_It++) //for each HSP, in ascending order (1st HSP is closest to start)!

		while (vec_It != HSPs.end())

		{

			cur_hsp_start = -(*vec_It)->HSP_end; //negate

			cur_hsp_end = -(*vec_It)->HSP_start; //negate



			cur_query_start = (*vec_It)->gene_start;

			cur_query_end = (*vec_It)->gene_end;



			cur_hsp_exists = true;



			if (vec_It_next != HSPs.end())

			{

				next_query_start = (*vec_It_next)->gene_start;

				next_query_end = (*vec_It_next)->gene_end;



				int next_hsp_start = -(*vec_It_next)->HSP_end;



				next_exists = GetSingleHSPDonorSegments_SplitOL_GapIsMinIntron(cur_hsp_start, cur_hsp_end, (*vec_It)->ID, 

					//donor_segments, 

					last_segment_start, last_segment_hsp_ID, 

					(*vec_It_next)->ID, next_hsp_start, 

					cur_query_start, cur_query_end, next_query_start, next_query_end, last_segment_end, cur_hsp_exists);



				if (next_query_start != (*vec_It_next)->gene_start)

					(*vec_It_next)->gene_start = next_query_start;

				if (cur_query_end != (*vec_It)->gene_end)

					(*vec_It)->gene_end = cur_query_end;

				if (cur_hsp_end != -(*vec_It)->HSP_start)

					(*vec_It)->HSP_start = -cur_hsp_end;

				if (next_hsp_start != -(*vec_It_next)->HSP_end)

					(*vec_It_next)->HSP_end = -next_hsp_start;



				//GetDonorSegments(cur_hsp_start, cur_hsp_end, prev_hsp_end, cur_start, (*vec_It)->ID, 

				//	donor_segments, prev_segment_end); //used for matchStr gap region

				if (!cur_hsp_exists)

				{

					HSPs_to_be_erased.push_back(*vec_It);

					HSPs_dup_to_be_erased.push_front(dup_i);//(*vec_It_dup);

				}



				if (!next_exists)

				{

					HSPs_to_be_erased.push_back(*vec_It_next);

					HSPs_dup_to_be_erased.push_front(dup_next_i);//(*vec_It_next_dup);



					vec_It++;

					vec_It++;

					//vec_It_dup++;

					//vec_It_dup++;

					dup_i+=2;



					vec_It_next++;

					//vec_It_next_dup++;

					dup_next_i++;

					if (vec_It_next != HSPs.end())

					{

						vec_It_next++;

						//vec_It_next_dup++;

						dup_next_i++;

					}

				}

				else

				{

					vec_It++;

					vec_It_next++;

					//vec_It_dup++;

					//vec_It_next_dup++;

					dup_i++;

					dup_next_i++;

				}

			}

			else //there's no next HSP

			{

				GetSingleHSPDonorSegments_SplitOL_GapIsMinIntron(cur_hsp_start, cur_hsp_end, (*vec_It)->ID, 

					//donor_segments, 

					last_segment_start, last_segment_hsp_ID, 

					-1, dummy_hsp_start, 0, dummy_query_end, dummy_query_start, INT_MAX, last_segment_end, cur_hsp_exists);

				//GetDonorSegments(cur_hsp_start, cur_hsp_end, prev_hsp_end, cur_start, (*vec_It)->ID, 

				//	donor_segments, prev_segment_end); //used for matchStr gap region

				if (!cur_hsp_exists)

				{

					HSPs_to_be_erased.push_back(*vec_It);

					HSPs_dup_to_be_erased.push_front(dup_i);//(*vec_It_dup);

				}



				vec_It++;

				//vec_It_dup++;

				dup_i++;



				last_segment_end = cur_hsp_end;

			}

		}



	}



	//donor_acceptor_HSP_ID.push_back(last_segment_hsp_ID); //need an extra one at the end, for last_segment



	//vector<HSP_Gene_Pair*>::iterator ptrHSPIt;

	list<int>::iterator ptrHSPIt;

/*	for (ptrHSPIt = HSPs_to_be_erased.begin(); ptrHSPIt != HSPs_to_be_erased.end(); ptrHSPIt++)

	{

#ifdef DEBUG

		dataMgr.outFile << "erasing HSPs: " << "\n";

		dataMgr.outFile << *(*ptrHSPIt) << "\n";

#endif

		HSPs.erase(remove(HSPs.begin(), HSPs.end(), *ptrHSPIt), HSPs.end());

	}

*/



/*	for (ptrHSPIt = HSPs.begin(); ptrHSPIt != HSPs.end(); ptrHSPIt++)

	{

#ifdef DEBUG

		dataMgr.outFile << "left hsp before processing: " << *(*ptrHSPIt) << "\n";

		dataMgr.outFile << (*(dataMgr.input_alignments.find((*ptrHSPIt)->ID))).second;

#endif

		//cut off possible gapped heading/trailing portions of HSPs (trailing or heading gaps at either end of HSP)

		Input_Alignment& align = (*(dataMgr.input_alignments.find((*ptrHSPIt)->ID))).second;

		CutHeadOrTrailGapsInAlignment(align, (*ptrHSPIt)->gene_start, (*ptrHSPIt)->gene_end, 

			(*ptrHSPIt)->HSP_start, (*ptrHSPIt)->HSP_end);



#ifdef DEBUG

		dataMgr.outFile << "left hsp after processing: " << *(*ptrHSPIt) << "\n";

		dataMgr.outFile << (*(dataMgr.input_alignments.find((*ptrHSPIt)->ID))).second;

#endif

	}

*/

	for (ptrHSPIt = HSPs_dup_to_be_erased.begin(); ptrHSPIt != HSPs_dup_to_be_erased.end(); ptrHSPIt++)

	{

#ifdef DEBUG

		dataMgr.outFile << "erasing dup HSPs: " << "\n";

		dataMgr.outFile << HSPs_dup[*ptrHSPIt] << "\n";

#endif

		HSPs_dup.erase(HSPs_dup.begin()+(*ptrHSPIt));

	}



	//recover from HSPs_dup!

	hspID_to_hspIndex_map.clear(); //use hspID_to_hspIndex_map, so later we can easily get from HSP_ID to the actual HSP in HSPs vector

	int i=0;

#ifdef DEBUG

	dataMgr.outFile << "recover from HSPs_dup:" << "\n";

#endif

	vector<HSP_Gene_Pair>::iterator HSPs_dup_it;

	for (HSPs_dup_it = HSPs_dup.begin(); HSPs_dup_it != HSPs_dup.end(); HSPs_dup_it++)

	{

		*(HSPs[i]) = *HSPs_dup_it; //assign the content

		map<int, Input_Alignment>::iterator align_it = dataMgr.input_alignments.find((*HSPs_dup_it).ID);

		const_cast<Input_Alignment&>((*align_it).second) = (*(dataMgr.input_alignments_HSPs_dup.find((*HSPs_dup_it).ID))).second;

		hspID_to_hspIndex_map.insert(map<int, int>::value_type((*HSPs_dup_it).ID, i++));

	}

	if (HSPs.size() > i)

		HSPs.erase(HSPs.begin()+i, HSPs.end());

#ifdef DEBUG

	PrintHSPs(HSPs);

#endif



/*

	HSPs.clear();

	//dataMgr.input_alignments.clear();

	HSPs_dup_copy.clear();

	HSPs_dup_copy.assign(HSPs_dup.begin(), HSPs_dup.end()); //need to make a copy, so we work on this copy later without changing HSPs_dup

	//for (HSPs_dup_it = HSPs_dup.begin(); HSPs_dup_it != HSPs_dup.end(); HSPs_dup_it++)

	for (HSPs_dup_it = HSPs_dup_copy.begin(); HSPs_dup_it != HSPs_dup_copy.end(); HSPs_dup_it++)

	{

		HSPs.push_back(&(*HSPs_dup_it)); //this is the address of the HSP, so later may get the HSP changed

		map<int, Input_Alignment>::iterator align_it = dataMgr.input_alignments.find(HSPs_dup_it->ID);

		const_cast<Input_Alignment&>((*align_it).second) = (*(dataMgr.input_alignments_HSPs_dup.find(HSPs_dup_it->ID))).second;

		hspID_to_hspIndex_map.insert(map<int, int>::value_type((*HSPs_dup_it).ID, i++));



#ifdef DEBUG

		dataMgr.outFile << (*HSPs_dup_it) << "\n";

		dataMgr.outFile << (*(dataMgr.input_alignments_HSPs_dup.find((HSPs_dup_it)->ID))).second << "\n";

#endif

	}

*/

}



bool ACCP_DONR_graph::LoadData4_forGenePredict6(vector<HSP_Gene_Pair*>& HSPs, 
												//pair<int, char*>& chr_seq, 
												vector<string>& chr_seq, 
												bool start_site_fixed, bool end_site_fixed)

{



	donors.clear();

	acceptors.clear();



	acceptor_region_index.clear();

	donor_region_index.clear();



	acceptor_head.clear();

	donor_tail.clear();



	//gap_centers.clear();



	donor_segments_pair.clear();

	donor_acceptor_HSP_ID.clear();

	donor_HSP_ID.clear();

	acceptor_HSP_ID.clear();



	donor_prev_2nt.clear();

	acceptor_next_2nt.clear();



//	int cur_start = start_site;

//	int prev_hsp_end = start_site-1;



	int last_segment_start=0; //used to record the possible last partial segment from last HSP

	int last_segment_end=0, first_segment_start;



//	vector<pair<int, int> > prev_segment_end; //used for matchStr gap region



	//vector<int> segment_hsp_start; //record the hsp_start for each donor_segment, used to compute acceptor_head/donor_tail frames

	//int last_hsp_start; //corresponding to [last_segment_start, last_segment_end]



	if (SPLICE_SEGMENT_VERSION == 1)

	{

		//HSPs_dup.clear();

		//dataMgr.input_alignments_HSPs_dup.clear();

		//GetSpliceSegments1(donor_segments, HSPs, first_segment_start, last_segment_start, last_segment_end);

		GetSpliceSegments1_forGenePredict6(HSPs, first_segment_start, last_segment_start, last_segment_end);

	}

	else

	{

		//GetSpliceSegments2(donor_segments, HSPs, first_segment_start, last_segment_start, last_segment_end, 

		if (!GetSpliceSegments2(HSPs, first_segment_start, last_segment_start, last_segment_end))//segment_hsp_start, last_hsp_start);

			return false;



//		ProcSegmentsByLens(donor_segments, last_segment_start, segment_hsp_start, last_hsp_start);

	}

	//MODIFIED: now put CalcStartEndPos2() here

	if (start_site > 0)

	{

		if (!start_site_fixed)

		{

			if (!end_site_fixed)

				CalcStartEndPos2( true, HSPs, chr_seq); //do the regular stuff

			else

				start_site = CalcStartPos(HSPs.back()->HSP_start, HSPs.back()->gene_start, true, chr_seq);

		}

		else //start_site already fixed, must be needing end_site

		{

			if (!end_site_fixed)

				end_site = CalcEndPos(HSPs.front()->HSP_end, HSPs.front()->gene_end, true, chr_seq);

		}

	}

	else

	{

		if (!start_site_fixed)

		{

			if (!end_site_fixed)

				CalcStartEndPos2( false, HSPs, chr_seq); //do the regular stuff

			else

				start_site = - CalcStartPos(HSPs.front()->HSP_end, HSPs.front()->gene_start, false, chr_seq);

		}

		else //start_site already fixed, must be needing end_site

		{

			if (!end_site_fixed)

				end_site = -CalcEndPos(HSPs.back()->HSP_start, HSPs.back()->gene_end, false, chr_seq);

		}

	}



	//remember: last segment after all HSPs(special treatment)!



	//now start searching for donors/acceptors

	//multimap<int, pair<int, int> >::iterator seg_map_It;

	vector< SegmentsInThreeBounds >::iterator seg_map_It;


	int t=0;

#ifdef DEBUG
	//for debug only, print out all segments
	//for (seg_map_It = donor_segments.begin(); seg_map_It != donor_segments.end(); seg_map_It++)
	for (seg_map_It = donor_segments_pair.begin(); seg_map_It != donor_segments_pair.end(); seg_map_It++, t++)
	{
		//dataMgr.outFile << "search: " << (*seg_map_It).first-1 << " to " << (*seg_map_It).second.first 
		//	<< "; " << (*seg_map_It).first << " to " << (*seg_map_It).second.second << "\n";
		dataMgr.outFile << "search: " << (*seg_map_It).exon_seg_end - 1 << " to " << (*seg_map_It).exon_seg_start 
			<< "; " << (*seg_map_It).exon_seg_end << " to " << (*seg_map_It).intron_seg_end;
		dataMgr.outFile << "; hspIDs:";
		for (int tt=0; tt<donor_acceptor_HSP_ID[t].size(); tt++)
			dataMgr.outFile << donor_acceptor_HSP_ID[t][tt] << ",";
//		if (SPLICE_SEGMENT_VERSION == 2)
//			dataMgr.outFile << "; segment_hsp_start (might be wrong!): " << segment_hsp_start[t];
		dataMgr.outFile << "\n";
	}
	dataMgr.outFile << "last_segment_start: " << last_segment_start << " - last_segment_end: " << last_segment_end;
	dataMgr.outFile << "; hspIDs:";
	for (int tt=0; tt<donor_acceptor_HSP_ID[t].size(); tt++)
		dataMgr.outFile << donor_acceptor_HSP_ID[t][tt] << ",";
	dataMgr.outFile << "\n";
#endif

	//MODIFIED: adjust donor_segments_pair (start_site/end_site may be in the middle, due to repair exons)
	t = 0;
	while ( t < donor_segments_pair.size() )
	{
		if (donor_segments_pair[t].exon_seg_end <= start_site)
		{
			donor_segments_pair.erase(donor_segments_pair.begin()+t);
			donor_acceptor_HSP_ID.erase(donor_acceptor_HSP_ID.begin()+t);
			continue;
		}

		if (donor_segments_pair[t].exon_seg_start >= end_site)
		{
			donor_segments_pair.erase(donor_segments_pair.begin()+t, donor_segments_pair.end());
			donor_acceptor_HSP_ID.erase(donor_acceptor_HSP_ID.begin()+t, donor_acceptor_HSP_ID.end());
			if (donor_segments_pair.size() > 0)
			{
				last_segment_start = donor_segments_pair.back().exon_seg_start;
				last_segment_end = donor_segments_pair.back().exon_seg_end-1;
				donor_segments_pair.pop_back();
			}
			else
			{
				return false;
			}
			break;
		}

		t++;
	}

#ifdef DEBUG
	OutputDonorAcceptorWithAllInfo(start_site > 0);
#endif

	int search_start, search_end;

	//int seg_map_size = donor_segments.size();

	int seg_map_size = donor_segments_pair.size();

	int cur_seg = 0;

	vector<int> accs, dons;

	vector< string > acc_next_2nt, don_prev_2nt;



	map<int, int> hsp_donor_stop_pos, hsp_acceptor_stop_pos;

	

	//for (seg_map_It = donor_segments.begin(); seg_map_It != donor_segments.end(); seg_map_It++)

	for (seg_map_It = donor_segments_pair.begin(); cur_seg < seg_map_size; seg_map_It++, cur_seg++) //do all except last one

	{

		search_start = (*seg_map_It).exon_seg_end - 1; //(*seg_map_It).first-1;



		search_end = (*seg_map_It).exon_seg_start; //(*seg_map_It).second.first;



		//MODIFIED: use the largest HSP in each splice segment as the only HSP there (used for donor_acceptor_HSP_ID)

		UseLargestHSPInSegment(cur_seg, HSPs, search_start, search_end);



		if (search_start > search_end)

		{

			//if (search_start == 749824 && search_end == 749699)

			//	int stophere=1;



			//GetDonorsAcceptors_1PerFrame(search_start, search_end, chr_seq, accs, dons);
#ifdef DEBUG
			dataMgr.outFile << "GetDonorsAcceptors_nPerBorder search_start:" << search_start << ";search_end:" << search_end << "\n";
#endif
			GetDonorsAcceptors_nPerBorder(search_start, search_end, chr_seq, accs, dons, acc_next_2nt, don_prev_2nt);
#ifdef DEBUG
			dataMgr.outFile << "GetDonorsAcceptors_nPerBorder ended" << "\n";
#endif



			//Added: remove any acceptor beyond search_start, or any donor beyond search_end?

			vector<int>::reverse_iterator check_acc_it;

			while (!accs.empty()) 

			{

				check_acc_it = accs.rbegin();

				if (*check_acc_it > search_start)

					accs.pop_back();

				else

					break;

			}

/*			vector<int>::iterator check_don_it = dons.begin();

			while (*check_don_it < search_end)

			{

				dons.erase(check_don_it);

				check_don_it = dons.begin();

			}

*/		}



#ifdef DEBUG		



		dataMgr.outFile << "now: searched " << search_start << "-" << search_end << "\n";

		dataMgr.outFile << "accs:" << "\n";

		vector<int>::iterator t_it;

		int t=0;

		for (t_it = accs.begin(); t_it != accs.end(); t_it++, t++)

			dataMgr.outFile << *t_it << ", " << acc_next_2nt[t] << ";";

		dataMgr.outFile << "\n";



		dataMgr.outFile << "dons:" << "\n";

		t = 0;

		for (t_it = dons.begin(); t_it != dons.end(); t_it++, t++)

			dataMgr.outFile << *t_it << ", " << don_prev_2nt[t] << ";";

		dataMgr.outFile << "\n";



#endif



		if (seg_map_It != donor_segments_pair.begin()) //first segment doesn't need acceptor

		{

			//combine same frame acceptors, so that there's at most 1 acceptor for each frame around this segment border

			//CombineInFrameSpliceSites(search_end, accs, acceptors, acceptor_region_index, cur_seg, true, acceptor_head);

//			if (SPLICE_SEGMENT_VERSION == 1)

//			{

				//int hsp_index = (*(hspID_to_hspIndex_map.find(donor_acceptor_HSP_ID[cur_seg]))).second;

				//if (start_site > 0)

					CombineSpliceSites_S1(search_end, accs, acceptors, acceptor_region_index, cur_seg, true, acceptor_head, 

						//HSPs[hsp_index]->HSP_start, //search_end, 

						//donor_acceptor_HSP_ID[cur_seg], 

						HSPs, chr_seq, hsp_donor_stop_pos, hsp_acceptor_stop_pos, 

						acceptor_HSP_ID, donor_HSP_ID, acc_next_2nt, acceptor_next_2nt);

				//else

				//	CombineSpliceSites_S1(search_end, accs, acceptors, acceptor_region_index, cur_seg, true, acceptor_head, 

				//		-HSPs[hsp_index]->HSP_end, //search_end, 

				//		donor_acceptor_HSP_ID[cur_seg], acceptor_HSP_ID, donor_HSP_ID, acc_next_2nt, acceptor_next_2nt);

					

/*

				map<int, int>::iterator mapIt = hspID_to_hspIndex_map.find(donor_acceptor_HSP_ID[cur_seg]);

				if (mapIt == hspID_to_hspIndex_map.end()) //not found, something wrong!

				{

					cout << "invalid hsp_ID... why? " << donor_acceptor_HSP_ID[cur_seg] << "\n"; //because the last/first HSP was erased!

					exit(-1);

				}

				if (start_site > 0)

					CombineSpliceSites_S1(search_end, accs, acceptors, acceptor_region_index, cur_seg, true, acceptor_head, 

						HSPs[(*mapIt).second]->HSP_start, //search_end, 

						donor_acceptor_HSP_ID[cur_seg], acceptor_HSP_ID, donor_HSP_ID);

				else

					CombineSpliceSites_S1(search_end, accs, acceptors, acceptor_region_index, cur_seg, true, acceptor_head, 

						-(HSPs[(*mapIt).second]->HSP_end), //search_end, 

						donor_acceptor_HSP_ID[cur_seg], acceptor_HSP_ID, donor_HSP_ID);

*/

//			}

//			else

//				CombineSpliceSites(search_end, accs, acceptors, acceptor_region_index, cur_seg, true, acceptor_head, 

//					search_end, //segment_hsp_start[cur_seg], 

//					acc_next_2nt, acceptor_next_2nt, dataMgr.outFile);

		}

		else

		{

			accs.erase(accs.begin(), accs.end());

			acc_next_2nt.erase(acc_next_2nt.begin(), acc_next_2nt.end());

		}



#ifdef DEBUG		



		dataMgr.outFile << "after combine: searched " << search_start << "-" << search_end << "\n";



		//PrintDonorAcceptor();
		OutputDonorAcceptorWithAllInfo(start_site > 0);



#endif







		search_start++;



		search_end = (*seg_map_It).intron_seg_end; //(*seg_map_It).second.second;



		if (search_start < search_end) //now this is the gap region

		{

			//GetDonorsAcceptors_1PerFrame(search_start, search_end, chr_seq, accs, dons);

			GetDonorsAcceptors_nPerBorder(search_start, search_end, chr_seq, accs, dons, acc_next_2nt, don_prev_2nt);



			//gap_centers.push_back((search_start + search_end)/2); //used for matchStr gap region



			//Added: remove any acceptor beyond search_start, or any donor beyond search_end?

			vector<int>::reverse_iterator check_don_it;

			while (!dons.empty()) 

			{

				check_don_it = dons.rbegin();

				if (*check_don_it > search_end)

					dons.pop_back();

				else

					break;

			}

/*			vector<int>::iterator check_acc_it = accs.begin();

			while (*check_acc_it < search_start)

			{

				accs.erase(check_acc_it);

				check_acc_it = accs.begin();

			}

*/		}



#ifdef DEBUG		



		dataMgr.outFile << "now: searched " << search_start << "-" << search_end << "\n";

		dataMgr.outFile << "accs:" << "\n";



		t=0;

		for (t_it = accs.begin(); t_it != accs.end(); t_it++,t++)

			dataMgr.outFile << *t_it << ", " << acc_next_2nt[t] << ";";

		dataMgr.outFile << "\n";



		dataMgr.outFile << "dons:" << "\n";

		t=0;

		for (t_it = dons.begin(); t_it != dons.end(); t_it++,t++)

			dataMgr.outFile << *t_it << ", " << don_prev_2nt[t] << ";";

		dataMgr.outFile << "\n";



#endif



		//if (search_start == -400392 && search_end == -400412)
		//	int stophere = 1;



		//combine same frame donors, so that there's at most 1 donor for each frame

		//CombineInFrameSpliceSites(search_start-1, dons, donors, donor_region_index, cur_seg, false, donor_tail);

//		if (SPLICE_SEGMENT_VERSION == 1)

//		{

			//int hsp_index = (*(hspID_to_hspIndex_map.find(donor_acceptor_HSP_ID[cur_seg]))).second;

			//if (start_site > 0)

				CombineSpliceSites_S1(search_start-1, dons, donors, donor_region_index, cur_seg, false, donor_tail, 

					//HSPs[hsp_index]->HSP_end, //search_start-1, 

					//donor_acceptor_HSP_ID[cur_seg], 

					HSPs, chr_seq, hsp_donor_stop_pos, hsp_acceptor_stop_pos, 

					acceptor_HSP_ID, donor_HSP_ID, don_prev_2nt, donor_prev_2nt);

			//else

			//	CombineSpliceSites_S1(search_start-1, dons, donors, donor_region_index, cur_seg, false, donor_tail, 

			//		-HSPs[hsp_index]->HSP_start, //search_start-1, 

			//		donor_acceptor_HSP_ID[cur_seg], acceptor_HSP_ID, donor_HSP_ID, don_prev_2nt, donor_prev_2nt);



/*			map<int, int>::iterator mapIt = hspID_to_hspIndex_map.find(donor_acceptor_HSP_ID[cur_seg]);

			if (mapIt == hspID_to_hspIndex_map.end()) //not found, something wrong!

			{

				cout << "invalid hsp_ID... why? " << donor_acceptor_HSP_ID[cur_seg] << "\n"; //because the last/first HSP was erased!

				exit(-1);

			}

			if (start_site > 0)

				CombineSpliceSites_S1(search_start-1, dons, donors, donor_region_index, cur_seg, false, donor_tail, 

					HSPs[(*mapIt).second]->HSP_end, //search_start-1, 

					donor_acceptor_HSP_ID[cur_seg], acceptor_HSP_ID, donor_HSP_ID);

			else

				CombineSpliceSites_S1(search_start-1, dons, donors, donor_region_index, cur_seg, false, donor_tail, 

					-(HSPs[(*mapIt).second]->HSP_start), //search_start-1, 

					donor_acceptor_HSP_ID[cur_seg], acceptor_HSP_ID, donor_HSP_ID);

*/

//		}

//		else

//			CombineSpliceSites(search_start-1, dons, donors, donor_region_index, cur_seg, false, donor_tail, 

//				search_start-1, //segment_hsp_start[cur_seg]+2, 

//				don_prev_2nt, donor_prev_2nt, dataMgr.outFile);





#ifdef DEBUG



		dataMgr.outFile << "after combine: searched " << search_start << "-" << search_end << "\n";



		//PrintDonorAcceptor();
		OutputDonorAcceptorWithAllInfo(start_site > 0);



#endif



	}



	//now process last partial segment

	//search_start = end_site;

	search_start = last_segment_end;

	search_end = last_segment_start;



	UseLargestHSPInSegment(cur_seg, HSPs, search_start, search_end);



	//used for matchStr gap region

/*	if (prev_segment_end.size() == 1) //because HSPs now do no overlap (after merging), 

		//there can be only one segment in "last_segment_end", and there must be one!

	{

		search_start = prev_segment_end.back().second; 

		search_end = prev_segment_end.back().first;

	}

	else

	{

		cout << "wrong last segment end size " << prev_segment_end.size() << "\n";

		exit(-1);

	}

*/

	//GetDonorsAcceptors_1PerFrame(search_start, search_end, chr_seq, accs, dons);

	GetDonorsAcceptors_nPerBorder(search_start, search_end, chr_seq, accs, dons, acc_next_2nt, don_prev_2nt);



			//Added: remove any acceptor beyond search_start, or any donor beyond search_end?

			vector<int>::reverse_iterator check_acc_it;

			while (!accs.empty()) 

			{

				check_acc_it = accs.rbegin();

				if (*check_acc_it > search_start)

					accs.pop_back();

				else

					break;

			}

/*			vector<int>::iterator check_don_it = dons.begin();

			while (*check_don_it < search_end)

			{

				dons.erase(check_don_it);

				check_don_it = dons.begin();

			}

*/

#ifdef DEBUG		



		dataMgr.outFile << "now: searched " << search_start << "-" << search_end << "\n";

		dataMgr.outFile << "accs:" << "\n";

		vector<int>::iterator t_it;

		for (t_it = accs.begin(); t_it != accs.end(); t_it++)

			dataMgr.outFile << *t_it << ", " ;

		dataMgr.outFile << "\n";



		dataMgr.outFile << "dons:" << "\n";

		for (t_it = dons.begin(); t_it != dons.end(); t_it++)

			dataMgr.outFile << *t_it << ", " ;

		dataMgr.outFile << "\n";



#endif



	//CombineInFrameSpliceSites(search_end, accs, acceptors, acceptor_region_index, cur_seg, true, acceptor_head);

//	if (SPLICE_SEGMENT_VERSION == 1)

		//int hsp_index = (*(hspID_to_hspIndex_map.find(donor_acceptor_HSP_ID[cur_seg]))).second;

		//if (start_site > 0)

			CombineSpliceSites_S1(search_end, accs, acceptors, acceptor_region_index, cur_seg, true, acceptor_head, 

				//HSPs[hsp_index]->HSP_start, //search_end, 

				//donor_acceptor_HSP_ID[cur_seg], 

				HSPs, chr_seq, hsp_donor_stop_pos, hsp_acceptor_stop_pos, 

				acceptor_HSP_ID, donor_HSP_ID, acc_next_2nt, acceptor_next_2nt);

		//else

		//	CombineSpliceSites_S1(search_end, accs, acceptors, acceptor_region_index, cur_seg, true, acceptor_head, 

		//		-HSPs[hsp_index]->HSP_end, //search_end, 

		//		donor_acceptor_HSP_ID[cur_seg], acceptor_HSP_ID, donor_HSP_ID, acc_next_2nt, acceptor_next_2nt);



//	else

//		CombineSpliceSites(search_end, accs, acceptors, acceptor_region_index, cur_seg, true, acceptor_head, 

//			search_end, //last_hsp_start, 

//			acc_next_2nt, acceptor_next_2nt, dataMgr.outFile);



#ifdef DEBUG



		dataMgr.outFile << "after combine: searched " << search_start << "-" << search_end << "\n";

		//PrintDonorAcceptor();
		OutputDonorAcceptorWithAllInfo(start_site > 0);



#endif


	//MODIFIED: now combine last_segment_start ~ last_segment_end to donor_segments_pair, so that we have everything stored 
	//(in case we need to adjust splice segments later, due to finding stop codon when computing exons)
	donor_segments_pair.push_back(SegmentsInThreeBounds(last_segment_start, last_segment_end+1, last_segment_end+1));
//#ifdef VERBOSE
	if (VERBOSE)
	cout << "load data done" << "\n";
//#endif
	return PostProcessSpliceSites(first_segment_start, last_segment_end, chr_seq);
}


int ACCP_DONR_graph::FindDonorHSPIndex_UseLargestHSPInSegment(int donor_site, vector<HSP_Gene_Pair*>& HSPs, int donor_segment_index, 

					   int& cur_hsp_start, int& cur_hsp_end, int& gene_start, int& gene_end, 

					   int& j) //j is the relative index, hsp_index is the HSP_ID

{

	int hsp_index=-1;

	int tmp_start;

	int i;

	vector<int>::reverse_iterator segment_HSP_ID_revIt;

	segment_HSP_ID_revIt = donor_acceptor_HSP_ID[donor_segment_index].rbegin(); 

	if (donor_site > 0) //positive strand HSP

	{

			i = (*(hspID_to_hspIndex_map.find(*segment_HSP_ID_revIt))).second;



			tmp_start = HSPs[i]->HSP_start; //(*vec_revIt)->HSP_start;



			if (tmp_start <= donor_site)

			{

				cur_hsp_start = tmp_start;

				cur_hsp_end = HSPs[i]->HSP_end; //(*vec_revIt)->HSP_end;

				gene_start = HSPs[i]->gene_start; //(*vec_revIt)->gene_start;

				gene_end = HSPs[i]->gene_end; //(*vec_revIt)->gene_end;

				hsp_index = HSPs[i]->ID; //(*vec_revIt)->ID;

				j=i;

			}

	}

	else //negative HSP

	{



			i = (*(hspID_to_hspIndex_map.find(*segment_HSP_ID_revIt))).second;



			tmp_start = - HSPs[i]->HSP_end; //-(*vec_It)->HSP_end; //negate



			if (tmp_start <= donor_site)

			{

			cur_hsp_start = tmp_start;

			cur_hsp_end = - HSPs[i]->HSP_start; //-(*vec_It)->HSP_start; //negate

			gene_start = HSPs[i]->gene_start; //(*vec_It)->gene_start;

			gene_end = HSPs[i]->gene_end; //(*vec_It)->gene_end;

			hsp_index = HSPs[i]->ID; //(*vec_It)->ID;

			j=i;

			}

	}



	return hsp_index;

}




//find the HSP according to the last HSP_ID in donor_acceptor_HSP_ID (the largest HSP), assume UseLargestHSPInSegment() is used!

int ACCP_DONR_graph::FindAcceptorHSPIndex_UseLargestHSPInSegment(int acceptor_site, vector<HSP_Gene_Pair*>& HSPs, int donor_segment_index, 

						  int& cur_hsp_start, int& cur_hsp_end, int& gene_start, int& gene_end, 

						  int& j)

{

	int hsp_index=-1;

	int tmp_end;

	int i;

	//vector<int>::iterator segment_HSP_ID_it;

	vector<int>::reverse_iterator segment_HSP_ID_it; //chose the last HSP!

	segment_HSP_ID_it = donor_acceptor_HSP_ID[donor_segment_index].rbegin(); //the last one in vector!

	if (acceptor_site > 0) //positive strand HSP

	{	

			i = (*(hspID_to_hspIndex_map.find(*segment_HSP_ID_it))).second;



			tmp_end = HSPs[i]->HSP_end; //(*vec_It)->HSP_end;



			//if (cur_hsp_start >= acceptor_site || (cur_hsp_start < acceptor_site && cur_hsp_end > acceptor_site))

			//if (cur_hsp_end >= acceptor_site)

			if (tmp_end >= acceptor_site)

			{

				hsp_index = HSPs[i]->ID; //(*vec_It)->ID;

				cur_hsp_start = HSPs[i]->HSP_start; //(*vec_It)->HSP_start;

				cur_hsp_end = tmp_end;

				gene_start = HSPs[i]->gene_start; //(*vec_It)->gene_start;

				gene_end = HSPs[i]->gene_end; //(*vec_It)->gene_end;

				j = i;

			}



	}

	else //negative HSP

	{

		i = (*(hspID_to_hspIndex_map.find(*segment_HSP_ID_it))).second;



			tmp_end = - HSPs[i]->HSP_start; //-(*vec_revIt)->HSP_start; //negate



			if (tmp_end >= acceptor_site)

			{

				hsp_index = HSPs[i]->ID; //(*vec_revIt)->ID;

				cur_hsp_start = - HSPs[i]->HSP_end; //-(*vec_revIt)->HSP_end; //negate

				cur_hsp_end = tmp_end;

				gene_start = HSPs[i]->gene_start; //(*vec_revIt)->gene_start;

				gene_end = HSPs[i]->gene_end; //(*vec_revIt)->gene_end;		

				j = i;

			}

	}



	return hsp_index;



}





//MODIFIED: do not use old frame_reference any more, compute site_head_tail based on the actual HSP start/end (determined by HSP_ID)

//MODIFIED AGAIN: do not use the "hsp_ID", compute the HSP that needs to be used based on the closest distance!

void ACCP_DONR_graph::CombineSpliceSites_S1(int border, vector<int>& sites, vector<int>& final_sites, 

											vector<int>& site_regions, 

							   int region, bool is_acceptor, vector<int>& site_head_tail, 

							   //int frame_reference, int hsp_ID, 

							   vector<HSP_Gene_Pair*>& HSPs, //use this to compute hsp_ID

							   //int segment_index, 

							   vector<string>& chr_seq, 
							   //pair<int, char*>& chr_seq, 

							   map<int, int>& hsp_donor_stop_pos, //after HSP's end, position of first stop

							   map<int, int>& hsp_acceptor_stop_pos, //before HSP's start, position of first stop codon

							   vector<int>& acceptor_HSP_ID, vector<int>& donor_HSP_ID, 

							   vector<string>& site_2nt, vector<string>& final_site_2nt)

{

	multimap<int,int> tmpSitesMap;



	//vector<int>::iterator site_it;

	int site_i;

	//for (site_it = sites.begin(); site_it != sites.end(); site_it++)

	for (site_i = 0; site_i < sites.size(); site_i++)

		//tmpSitesMap.insert(multimap<int, int>::value_type(abs(*site_it - border), *site_it)); //sort sites by their distance from border

		tmpSitesMap.insert(multimap<int, int>::value_type(abs(sites[site_i] - border), site_i));

	multimap<int, int>::iterator map_it = tmpSitesMap.begin();
/*
	//vector<int> tmpSites;
	map<int, string> tmpSites;
	//int count = 0;
	//while (count < MAX_NUM_SPLICE_SITES && map_it != tmpSitesMap.end())
	while (map_it != tmpSitesMap.end()) //sort the sites!?
	{
		//tmpSites.push_back((*map_it).second); //get MAX_NUM_SPLICE_SITES sites that are closest around border
		tmpSites.insert(map<int,string>::value_type(sites[(*map_it).second], site_2nt[(*map_it).second]));

		map_it++;
		//count++;
	}
	//sort(tmpSites.begin(), tmpSites.end()); //sort so we can store to final_sites
	//vector<int>::iterator tmp_site_it;
	map<int, string>::iterator tmp_site_it;
*/
	
	map<int, Triplet > site_info; //hold the current valid sites and its info (sort them using map)

	int cur_site;
	int frame_reference, hsp_ID, not_used, j;
	int count = 0;
	bool found;
//	tmp_site_it = tmpSites.begin();
//	while ( tmp_site_it != tmpSites.end() && count < MAX_NUM_SPLICE_SITES)
	while (map_it != tmpSitesMap.end() && count < MAX_NUM_SPLICE_SITES)
	{
		//cur_site = *tmp_site_it;
		//cur_site = (*tmp_site_it).first;

		cur_site = sites[(*map_it).second];

		if (is_acceptor) //for acceptor, record head
		{

			//frame_reference is cur_hsp_start!
			hsp_ID = FindAcceptorHSPIndex_UseLargestHSPInSegment(cur_site, HSPs, region, frame_reference, not_used, not_used, not_used, j);
			if (hsp_ID == -1)
			{
				//tmp_site_it++;
				map_it++;
				continue;
			}

/*			vector<HSP_Gene_Pair*> cur_hsps;

			found = FindAcceptorHSPs(cur_site, HSPs, segment_index, cur_hsps);

			if (!found)

			{

				tmp_site_it++;

				continue;

			}



			found = false; //use this for other purpose now

			for (j=0; j<cur_hsps.size(); j++)

			{

				if (start_site > 0)

					frame_reference = cur_hsps[j]->HSP_start;

				else

					frame_reference = -cur_hsps[j]->HSP_end;



				hsp_ID = cur_hsps[j]->ID;



*/
				//check for in-frame stop here! MODIFIED: ONLY check acceptor, not donor?! ComputeExons will check donors
				if (cur_site < frame_reference)
				{
					map<int, int>::iterator pos_it;
					if ((pos_it=hsp_acceptor_stop_pos.find(hsp_ID)) == hsp_acceptor_stop_pos.end())
					{
						string tmpStr3;
						if (cur_site > 0)
							GetSubstrFromVecStrs(chr_seq, true, cur_site-1, frame_reference-cur_site, tmpStr3);
						else
							GetSubstrFromVecStrs_NegRev(chr_seq, false, -frame_reference, frame_reference-cur_site, tmpStr3);
						StrToLower(tmpStr3);
						int stop_pos;
						if (HasInFrameStop(tmpStr3, false, stop_pos, false))
						{
							hsp_acceptor_stop_pos.insert(map<int,int>::value_type(hsp_ID, cur_site+stop_pos));
							//tmp_site_it++;
							map_it++;
							continue;
						}
					}
					else
					{
						if (cur_site <= (*pos_it).second)
						{
							//tmp_site_it++;
							map_it++;
							continue;
						}
						else //since we recorded the actual first stop codon position, cur_site now is past that, so we should be safe? NOT unless sites are sorted
						{
							string tmpStr3;
							if (cur_site > 0)
								GetSubstrFromVecStrs(chr_seq, true, cur_site-1, frame_reference-cur_site, tmpStr3);
							else
								GetSubstrFromVecStrs_NegRev(chr_seq, false, -frame_reference, frame_reference-cur_site, tmpStr3);
							StrToLower(tmpStr3);
							int stop_pos;
							if (HasInFrameStop(tmpStr3, false, stop_pos, false))
							{
								const_cast<int&>((*pos_it).second) = cur_site+stop_pos; //record the largest stop position
								//tmp_site_it++;
								map_it++;
								continue;
							}
						}
					}
				}
		
/*				//if (cur_site <= border)
				if (cur_site <= frame_reference)
					//site_head_tail.push_back((border - cur_site)%3);
					site_head_tail.push_back((frame_reference - cur_site)%3);
				else
					//site_head_tail.push_back((3 - (cur_site - border) % 3) % 3);
					site_head_tail.push_back((3 - (cur_site - frame_reference) % 3) % 3);

				acceptor_HSP_ID.push_back(hsp_ID);
				final_sites.push_back(cur_site);
				//final_site_2nt.push_back((*tmp_site_it).second);
				final_site_2nt.push_back(site_2nt[(*map_it).second]);
				site_regions.push_back(region);
*/
				site_info.insert(map<int, Triplet >::value_type(cur_site, 
					Triplet(hsp_ID, frame_reference, (*map_it).second)));
				found = true;
//			}
			

		}

		else //for donor, record tail
		{
			//if (cur_site == -1206537)
			//	int stophere = 1;

			hsp_ID = FindDonorHSPIndex_UseLargestHSPInSegment(cur_site, HSPs, region, not_used, frame_reference, not_used, not_used, j);
			if (hsp_ID == -1)
			{
				//tmp_site_it++;
				map_it++;
				continue;
			}

/*			vector<HSP_Gene_Pair*> cur_hsps;
			found = FindDonorHSPs(cur_site, HSPs, segment_index, cur_hsps);
			if (!found)
			{
				tmp_site_it++;
				continue;
			}



			found = false;

			for (j=0; j<cur_hsps.size(); j++)

			{

				hsp_ID = cur_hsps[j]->ID;

				if (start_site > 0)

					frame_reference = cur_hsps[j]->HSP_end;

				else

					frame_reference = -cur_hsps[j]->HSP_start;

*/

			/*if (cur_site > frame_reference)
			{
			map<int, int>::iterator pos_it;
			if ((pos_it=hsp_donor_stop_pos.find(hsp_ID)) == hsp_donor_stop_pos.end())
			{
				string tmpStr3;
				if (cur_site > 0)
					GetSubstrFromVecStrs(chr_seq, true, frame_reference, cur_site-frame_reference, tmpStr3);
				else
					GetSubstrFromVecStrs_NegRev(chr_seq, false, -cur_site-1, cur_site-frame_reference, tmpStr3);
				StrToLower(tmpStr3);
				int stop_pos;
				if (HasInFrameStop(tmpStr3, true, stop_pos))
				{
					hsp_donor_stop_pos.insert(map<int,int>::value_type(hsp_ID, frame_reference+stop_pos+1));
					tmp_site_it++;
					continue;
				}
			}
			else
			{
				if (cur_site >= (*pos_it).second)
				{
					tmp_site_it++;
					continue;
				}
				else//sites are sorted, is it possible to have cur_site smaller?
				{
					cout << "eh?" << "\n";
					exit(-1);
				}
			}
			}*/

/*			donor_HSP_ID.push_back(hsp_ID);
			//if (cur_site <= border)
			if (cur_site <= frame_reference)
				//site_head_tail.push_back((3 - (border - cur_site)%3) % 3);
				site_head_tail.push_back((3 - (frame_reference - cur_site)%3) % 3);
			else
				//site_head_tail.push_back( (cur_site - border) % 3);
				site_head_tail.push_back( (cur_site - frame_reference) % 3);
			final_sites.push_back(cur_site);
			//final_site_2nt.push_back((*tmp_site_it).second);
			final_site_2nt.push_back(site_2nt[(*map_it).second]);
			site_regions.push_back(region);
*/
			site_info.insert(map<int, Triplet >::value_type(cur_site, 
				Triplet(hsp_ID, frame_reference, (*map_it).second)));

			found = true;
//			}
		}

		if (found)
			count++;		

		//tmp_site_it++;
		map_it++;
	}

	map<int, Triplet>::iterator site_info_it;
	for (site_info_it = site_info.begin(); site_info_it != site_info.end(); site_info_it++)
	{
		cur_site = (*site_info_it).first;
		hsp_ID = (*site_info_it).second.hsp_ID;
		frame_reference = (*site_info_it).second.frame_reference;
		if (is_acceptor)
		{
			acceptor_HSP_ID.push_back(hsp_ID);
			if (cur_site <= frame_reference)
				site_head_tail.push_back((frame_reference - cur_site)%3);
			else
				site_head_tail.push_back((3 - (cur_site - frame_reference) % 3) % 3);
		}
		else
		{
			donor_HSP_ID.push_back(hsp_ID);
			if (cur_site <= frame_reference)
				site_head_tail.push_back((3 - (frame_reference - cur_site)%3) % 3);
			else
				site_head_tail.push_back( (cur_site - frame_reference) % 3);
		}
		final_sites.push_back(cur_site);
		final_site_2nt.push_back(site_2nt[(*site_info_it).second.site_index]);
		site_regions.push_back(region);
	}

	sites.clear(); //reset sites
	site_2nt.clear();
}



void ACCP_DONR_graph::PrintDonorAcceptorHSPID()

{

	dataMgr.outFile << "donor_acceptor_HSP_ID: " << "\n";

	int i, j;

	for (i = 0; i<donor_acceptor_HSP_ID.size(); i++)

	{

		dataMgr.outFile << "region " << i << ":";

		for (j=0; j<donor_acceptor_HSP_ID[i].size(); j++)

				dataMgr.outFile << donor_acceptor_HSP_ID[i][j] << ",";

		dataMgr.outFile << "\n";

	}

}






//MODIFIED: put the largest HSP at the end of vector, leave the rest HSPs unchanged (the largest HSP will occur twice in vector)!

void ACCP_DONR_graph::UseLargestHSPInSegment(int cur_seg, vector<HSP_Gene_Pair*>& HSPs, int search_start, int search_end)

{

#ifdef DEBUG

	dataMgr.outFile << "UseLargestHSPInSegment: donor_acceptor_HSP_ID[" << cur_seg << "]:";
	int k;
	for (k = 0; k < donor_acceptor_HSP_ID[cur_seg].size(); k++)

		dataMgr.outFile << donor_acceptor_HSP_ID[cur_seg][k] << ",";

	dataMgr.outFile << "\n";

#endif



		int seg_hsp_size = donor_acceptor_HSP_ID[cur_seg].size();

		if (seg_hsp_size > 1)

		{

			int max_overlap_len = 0, overlap_len, max_overlap_hspID=-1;

			for (int j=0; j<seg_hsp_size; j++)

			{

				int hsp_index = (*(hspID_to_hspIndex_map.find(donor_acceptor_HSP_ID[cur_seg][j]))).second;

				int cur_start_pos, cur_end_pos;

				if (start_site > 0)

				{

					cur_start_pos = HSPs[hsp_index]->HSP_start;

					cur_end_pos = HSPs[hsp_index]->HSP_end;

				}

				else

				{

					cur_start_pos = -HSPs[hsp_index]->HSP_end;

					cur_end_pos = -HSPs[hsp_index]->HSP_start;

				}

				if (cur_start_pos < search_end)

				{

					if (cur_end_pos >= search_end)

					{

						if (cur_end_pos < search_start)

							overlap_len = cur_end_pos - search_end;

						else

							overlap_len = search_start - search_end;

					}

				}

				else

				{

					if (cur_start_pos <= search_start)

					{

						if (cur_end_pos < search_start)

							overlap_len = cur_end_pos - cur_start_pos;

						else

							overlap_len = search_start - cur_start_pos;

					}

				}

				if (max_overlap_len < overlap_len)

				{

					max_overlap_len = overlap_len;

					max_overlap_hspID = donor_acceptor_HSP_ID[cur_seg][j];

				}

			}

			if (max_overlap_hspID == -1)

			{

				cout << "region: " << cur_seg << " max_overlap_hspID is -1?" << "\n";

				exit(-1);

			}

			//donor_acceptor_HSP_ID[cur_seg].clear();

			donor_acceptor_HSP_ID[cur_seg].push_back(max_overlap_hspID);



		}

		else //there's only 1 HSP_ID

		{

			if (donor_acceptor_HSP_ID[cur_seg].empty())

			{

				cout << "region " << cur_seg << " has no HSP_ID?" << "\n";

				exit(-1);

			}

			donor_acceptor_HSP_ID[cur_seg].push_back(donor_acceptor_HSP_ID[cur_seg].front());

		}



#ifdef DEBUG

	dataMgr.outFile << "after UseLargestHSPInSegment: donor_acceptor_HSP_ID[" << cur_seg << "]:";

	for (k = 0; k < donor_acceptor_HSP_ID[cur_seg].size(); k++)

		dataMgr.outFile << donor_acceptor_HSP_ID[cur_seg][k] << ",";

	dataMgr.outFile << "\n";

#endif



}





//stop_site is the stop codon start coordinate (is negative for negative strand)
//remove the HSP that's JUST BEFORE stop and is in-frame with that stop
//UPDATE: now if the stop_site is close to the end, remove the HSP that's after stop!
//UPDATE AGAIN: now compare the hsp to the left and the hsp to the right of stop_site, remove lower-quality one
bool ACCP_DONR_graph::RemoveHSPBeforeStop(vector<HSP_Gene_Pair*>& HSPs, int stop_site, 
										  //pair<int, char*>& chr_seq, 
										  vector<string>& chr_seq, 
										  bool& removed_hsp_is_new_hsp, bool remove_after_stop_preferred)

{
	int i;

#ifdef DEBUG
	dataMgr.outFile << "BEFORE REMOVING HSP --- HSPs:(stop_site:"<< stop_site << ")" << "\n";
	for (i=0; i<HSPs.size(); i++)
		dataMgr.outFile << *(HSPs[i]);
	dataMgr.outFile << "\n" << "splice segments:" << "\n";
	for (i=0; i<donor_segments_pair.size(); i++)
	{
		for (int t=0; t<donor_acceptor_HSP_ID[i].size(); t++)
			dataMgr.outFile << donor_acceptor_HSP_ID[i][t] << ",";
		dataMgr.outFile << donor_segments_pair[i];
	}
	dataMgr.outFile << "\n";
	OutputDonorAcceptorWithAllInfo(start_site>0);
#endif

	i=0;
	int hspIndex_after_stop, hspIndex_before_stop;
	bool hsp_after_stop_found = false, hsp_before_stop_found = false;
	while (i<HSPs.size())
	{
		if (stop_site > 0) //pos strand
		{
			if (!hsp_before_stop_found && HSPs[i]->HSP_end < stop_site && (stop_site - HSPs[i]->HSP_start) % 3 == 0)
			{
				hspIndex_before_stop = i;
				hsp_before_stop_found = true;
				break; //no need to go on
			}
			if (HSPs[i]->HSP_start > stop_site && (HSPs[i]->HSP_start - stop_site) % 3 == 0)
			{
				hspIndex_after_stop = i;
				hsp_after_stop_found = true; //go on
			}
		}
		else
		{
			if (!hsp_after_stop_found && HSPs[i]->HSP_end < -stop_site && (-stop_site - HSPs[i]->HSP_end) % 3 == 0 )
			{
				hspIndex_after_stop = i;
				hsp_after_stop_found = true;
				break;
			}
			if (HSPs[i]->HSP_start > -stop_site && (HSPs[i]->HSP_end + stop_site) % 3 == 0)
			{
				hspIndex_before_stop = i;
				hsp_before_stop_found = true;
			}
		}
		i++;
	}

	bool remove_after_stop = remove_after_stop_preferred;
	if (hsp_before_stop_found && hsp_after_stop_found)
	{
		Input_Alignment& align_before_stop = (*(dataMgr.input_alignments.find(HSPs[hspIndex_before_stop]->ID))).second;
		Input_Alignment& align_after_stop = (*(dataMgr.input_alignments.find(HSPs[hspIndex_after_stop]->ID))).second;
		int b = align_before_stop.match_align.length() - CountCharInStr(align_before_stop.match_align, ' ') - CountCharInStr(align_before_stop.match_align, '+');
		int a = align_after_stop.match_align.length() - CountCharInStr(align_after_stop.match_align, ' ') - CountCharInStr(align_after_stop.match_align, '+');
		if (b > a)
			remove_after_stop = true;
		else
			if (b < a)
				remove_after_stop = false;
	}
	
	removed_hsp_is_new_hsp = false;

	bool found = false;
	if (remove_after_stop && hsp_after_stop_found)
	{
		found = true;
		i = hspIndex_after_stop;
	}
	else
		if (!remove_after_stop && hsp_before_stop_found)
		{
			found = true;
			i = hspIndex_before_stop;
		}
		else
			return false;		

	if (HSPs[i]->ID == hspID_added_for_repair)
	{
		removed_hsp_is_new_hsp = true;
		return true;
	}

	bool first_hsp_erased = false;
	bool last_hsp_erased = false;

#ifdef DEBUG
	dataMgr.outFile << "removing HSP[" << HSPs[i]->ID << "]\n";
#endif

#ifndef COMPUT_EXON_FULL_STEP_BACK
	first_hsp_erased = RemoveSegments(HSPs[i]->ID, i, HSPs, chr_seq);
	if ((start_site>0 && i==0) || (start_site<0 && i==HSPs.size()-1))
		last_hsp_erased = true;
#endif

	HSPs.erase(HSPs.begin()+i);

#ifndef COMPUT_EXON_FULL_STEP_BACK
	hspID_to_hspIndex_map.clear();
	for (int k=0; k<HSPs.size(); k++)
		hspID_to_hspIndex_map.insert(map<int,int>::value_type(HSPs[k]->ID, k));
#endif


/*	if (start_site > 0) //positive strand
	{
		i = 0;
		while (i<HSPs.size())
		{
			if (HSPs[i]->HSP_end < stop_site && ((stop_site - HSPs[i]->HSP_start) % 3 == 0) ) 
			{
#ifdef DEBUG
				dataMgr.outFile << "removing HSP[" << HSPs[i]->ID << "]\n";
#endif
				if (HSPs[i]->ID == hspID_added_for_repair)
				{
					removed_hsp_is_new_hsp = true;
					return true;
				}
				found = true;

#ifndef COMPUT_EXON_FULL_STEP_BACK
				first_hsp_erased = RemoveSegments(HSPs[i]->ID, i, HSPs, chr_seq);
#endif
				HSPs.erase(HSPs.begin()+i);

#ifndef COMPUT_EXON_FULL_STEP_BACK
				hspID_to_hspIndex_map.clear();
				for (int k=0; k<HSPs.size(); k++)
					hspID_to_hspIndex_map.insert(map<int,int>::value_type(HSPs[k]->ID, k));
#endif
				break;
			}
			i++;
		}
	}
	else //negative strand
	{
		i = HSPs.size()-1;
		while (i>=0)
		{
			if (-HSPs[i]->HSP_start < stop_site && ((stop_site + HSPs[i]->HSP_end) % 3 == 0) )
			{
#ifdef DEBUG
				dataMgr.outFile << "removing HSP[" << HSPs[i]->ID << "]\n";
#endif
				if (HSPs[i]->ID == hspID_added_for_repair)
				{
					removed_hsp_is_new_hsp = true;
					return true;
				}
				found = true;
#ifndef COMPUT_EXON_FULL_STEP_BACK
				first_hsp_erased = RemoveSegments(HSPs[i]->ID, i, HSPs, chr_seq);
#endif
				HSPs.erase(HSPs.begin()+i);

#ifndef COMPUT_EXON_FULL_STEP_BACK
				hspID_to_hspIndex_map.clear();
				for (int k=0; k<HSPs.size(); k++)
					hspID_to_hspIndex_map.insert(map<int,int>::value_type(HSPs[k]->ID, k));
#endif
				break;
			}
			i--;
		}
	}

#ifndef COMPUT_EXON_FULL_STEP_BACK
	//check if we need a new start_site
	if (first_hsp_erased && HSPs.size() > 0)
	{
		if (start_site > 0)
			start_site = CalcStartPos(HSPs.back()->HSP_start, HSPs.back()->gene_start, true, chr_seq);
		else
			start_site = - CalcStartPos(HSPs.front()->HSP_end, HSPs.front()->gene_start, false, chr_seq);
	}
#endif
*/

/*	if ((start_site > 0) != remove_after_stop) //positive strand and before stop (or, neg strand and after stop)
	{
		i = 0;
		while (i<HSPs.size())
		{
			bool found_first_eligible_hsp;
			if (remove_after_stop) //must be neg strand
				found_first_eligible_hsp = (HSPs[i]->HSP_end < -stop_site) && ((-stop_site - HSPs[i]->HSP_end) % 3 == 0);
			else //must be pos strand
				found_first_eligible_hsp = (HSPs[i]->HSP_end < stop_site) && ((stop_site - HSPs[i]->HSP_start) % 3 == 0);
				
			if  (found_first_eligible_hsp)
			{
#ifdef DEBUG
				dataMgr.outFile << "removing HSP[" << HSPs[i]->ID << "]\n";
#endif
				if (HSPs[i]->ID == hspID_added_for_repair)
				{
					removed_hsp_is_new_hsp = true;
					return true;
				}
				found = true;

#ifndef COMPUT_EXON_FULL_STEP_BACK
				first_hsp_erased = RemoveSegments(HSPs[i]->ID, i, HSPs, chr_seq);
				if ((start_site>0 && i==0) || (start_site<0 && i==HSPs.size()-1))
					last_hsp_erased = true;
#endif
				HSPs.erase(HSPs.begin()+i);

#ifndef COMPUT_EXON_FULL_STEP_BACK
				hspID_to_hspIndex_map.clear();
				for (int k=0; k<HSPs.size(); k++)
					hspID_to_hspIndex_map.insert(map<int,int>::value_type(HSPs[k]->ID, k));
#endif
				break;
			}
			i++;
		}
	}
	else //negative strand and before stop (pos strand && after stop)
	{
		i = HSPs.size()-1;
		while (i>=0)
		{
			bool found_first_eligible_hsp;
			if (remove_after_stop) //must be pos strand
				found_first_eligible_hsp = (HSPs[i]->HSP_start > stop_site) && ((HSPs[i]->HSP_start - stop_site) % 3 == 0);
			else //must be neg strand
				found_first_eligible_hsp = (-HSPs[i]->HSP_start < stop_site) && ((stop_site + HSPs[i]->HSP_end) % 3 == 0);
			if (found_first_eligible_hsp )
			{
#ifdef DEBUG
				dataMgr.outFile << "removing HSP[" << HSPs[i]->ID << "]\n";
#endif
				if (HSPs[i]->ID == hspID_added_for_repair)
				{
					removed_hsp_is_new_hsp = true;
					return true;
				}
				found = true;
#ifndef COMPUT_EXON_FULL_STEP_BACK
				first_hsp_erased = RemoveSegments(HSPs[i]->ID, i, HSPs, chr_seq);
				if ((start_site > 0 && i == 0) || (start_site < 0 && i == HSPs.size()-1))
					last_hsp_erased = true;				
#endif
				HSPs.erase(HSPs.begin()+i);

#ifndef COMPUT_EXON_FULL_STEP_BACK
				hspID_to_hspIndex_map.clear();
				for (int k=0; k<HSPs.size(); k++)
					hspID_to_hspIndex_map.insert(map<int,int>::value_type(HSPs[k]->ID, k));
#endif
				break;
			}
			i--;
		}
	}
*/

#ifndef COMPUT_EXON_FULL_STEP_BACK
	//check if we need a new start_site
	if (HSPs.size() > 0)
	{
		if (first_hsp_erased)
			if (start_site > 0)
				start_site = CalcStartPos(HSPs.back()->HSP_start, HSPs.back()->gene_start, true, chr_seq);
			else
				start_site = - CalcStartPos(HSPs.front()->HSP_end, HSPs.front()->gene_start, false, chr_seq);
		
		if (last_hsp_erased)
			if (start_site > 0)
				end_site = CalcEndPos(HSPs.front()->HSP_end, HSPs.front()->gene_end, true, chr_seq);
			else
				end_site = -CalcEndPos(HSPs.back()->HSP_start, HSPs.back()->gene_end, false, chr_seq);

	}
#endif

	//make sure start_site and end_site are still legit
	if (start_site >= end_site)
		return false;

#ifdef DEBUG
	dataMgr.outFile << "After removing HSP before/after stop:"<< "\n";
	for (i=0; i<HSPs.size(); i++)
		dataMgr.outFile << *(HSPs[i]);
	dataMgr.outFile << "\n" << "splice segments:" << "\n";
	for (i=0; i<donor_segments_pair.size(); i++)
	{
		for (int t=0; t<donor_acceptor_HSP_ID[i].size(); t++)
			dataMgr.outFile << donor_acceptor_HSP_ID[i][t] << ",";
		dataMgr.outFile << donor_segments_pair[i];
	}
	dataMgr.outFile << "\n";
	OutputDonorAcceptorWithAllInfo(start_site>0);
#endif

	return found;

}



bool ACCP_DONR_graph::RemoveSegments(int hsp_ID, int hsp_index, vector<HSP_Gene_Pair*>& HSPs, 
									 //pair<int, char*>& chr_seq)
									 vector<string>& chr_seq)

{

/*	int hsp_start, hsp_end;

	if (start_site > 0)

	{

		hsp_start = HSPs[hsp_index]->HSP_start;

		hsp_end = HSPs[hsp_index]->HSP_end;

	}

	else

	{

		hsp_start = -HSPs[hsp_index]->HSP_end;

		hsp_end = -HSPs[hsp_index]->HSP_start;

	}

*/

	int tmp;

	acceptor_region_index.clear();

	for (tmp=0; tmp<acceptor_region_index_copy.size(); tmp++)

		acceptor_region_index.push_back(acceptor_region_index_copy[tmp]);

	donor_region_index.clear();

	for (tmp=0; tmp<donor_region_index_copy.size(); tmp++)

		donor_region_index.push_back(donor_region_index_copy[tmp]);



	int i=0;

	int j;

	bool found = false;

	bool hsp_in_mid_segment = false, seg_check_done = false;

	bool first_hsp_erased = false;



	int donor_erase_region_start = -1, donor_erase_region_end = -1;

	int acceptor_erase_region_start = -1, acceptor_erase_region_end = -1;



	vector<SegmentsInThreeBounds> newsites_donor_segments_pair; //for getting new sites!
	//bool newsites_donor_segments_pair_last = false; //signals if the new segments are all done

	int start_region_index;

	

	vector<int> donor_acceptor_HSP_ID_to_be_erased;

	int cur_first_segment = 0;
	int num_segment_erased = 0;
	while (i<donor_acceptor_HSP_ID.size())
	{	
		bool cur_found = false;
		int cur_size = donor_acceptor_HSP_ID[i].size();
		for (j = 0; j < cur_size-1; j++) //skip the last one, since it's a duplicate!
		{
			if (donor_acceptor_HSP_ID[i][j] == hsp_ID) //found this hsp
			{
				found = true;
				cur_found = true;
				if (j == 0)
				{
					if (j+1 == cur_size-1) //j is the one and only HSP in this segment
					{
						if (i > cur_first_segment)
						{
							donor_segments_pair[i-num_segment_erased-1].intron_seg_end = donor_segments_pair[i-num_segment_erased].intron_seg_end;
						}
						else
						{
							first_hsp_erased = true;
							cur_first_segment++;
						}
						//update newsites_donor_segments_pair
						if (!newsites_donor_segments_pair.empty())
							newsites_donor_segments_pair.back().intron_seg_end = donor_segments_pair[i-num_segment_erased].intron_seg_end;
						else
						{
							if (i > cur_first_segment)
								newsites_donor_segments_pair.push_back(donor_segments_pair[i-num_segment_erased-1]);
						}

						donor_segments_pair.erase(donor_segments_pair.begin()+i-num_segment_erased);
						num_segment_erased++;
						donor_acceptor_HSP_ID_to_be_erased.push_back(i);

						if (donor_erase_region_start == -1) //first time
						{
							donor_erase_region_start = i-1;
							acceptor_erase_region_start = i;
						}
						donor_erase_region_end = i;
						acceptor_erase_region_end = i+1;
					}
					else //j is the first one, but not the only one
					{
						int next_hsp_index = (*(hspID_to_hspIndex_map.find(donor_acceptor_HSP_ID[i][j+1]))).second;
						int next_hsp_start;
						if (start_site > 0)
							next_hsp_start = HSPs[next_hsp_index]->HSP_start;
						else
							next_hsp_start = -HSPs[next_hsp_index]->HSP_end;

						if (i>cur_first_segment)
						{
							donor_segments_pair[i-num_segment_erased-1].intron_seg_end = next_hsp_start-1;
						}
						else
						{
							first_hsp_erased = true;
							//cur_first_segment++; //no need since we are breaking out anyway
						}
						donor_segments_pair[i-num_segment_erased].exon_seg_start = next_hsp_start;

						if (donor_acceptor_HSP_ID[i].back() == donor_acceptor_HSP_ID[i].front())
						{
							donor_acceptor_HSP_ID[i].erase(donor_acceptor_HSP_ID[i].begin());
							donor_acceptor_HSP_ID[i].pop_back();
							UseLargestHSPInSegment(i, HSPs, donor_segments_pair[i-num_segment_erased].exon_seg_end-1, 
								donor_segments_pair[i-num_segment_erased].exon_seg_start);
						}
						else
						{
							donor_acceptor_HSP_ID[i].erase(donor_acceptor_HSP_ID[i].begin());
						}

						//update newsites_donor_segments_pair
						if (!newsites_donor_segments_pair.empty())
							newsites_donor_segments_pair.back().intron_seg_end = next_hsp_start-1;
						else
						{
							if (i>cur_first_segment)
								newsites_donor_segments_pair.push_back(donor_segments_pair[i-num_segment_erased-1]);
						}
						newsites_donor_segments_pair.push_back(donor_segments_pair[i-num_segment_erased]);
						//newsites_donor_segments_pair_last = true;

						if (donor_erase_region_start == -1)
						{
							donor_erase_region_start = i-1;
							acceptor_erase_region_start = i;
						}
						donor_erase_region_end = i;
						acceptor_erase_region_end = i; //same as start

						seg_check_done = true;
						break;
					}
				}
				else
				{
					if (j+1==cur_size-1) //j is the last one, and not the only one
					{
						int prev_hsp_index = (*(hspID_to_hspIndex_map.find(donor_acceptor_HSP_ID[i][j-1]))).second;
						int prev_hsp_end;

						if (start_site > 0)
							prev_hsp_end = HSPs[prev_hsp_index]->HSP_end;
						else
							prev_hsp_end = -HSPs[prev_hsp_index]->HSP_start;

						if (i<donor_acceptor_HSP_ID.size()-1)
							donor_segments_pair[i].exon_seg_end = prev_hsp_end+1;

						if (donor_acceptor_HSP_ID[i].back() == donor_acceptor_HSP_ID[i][j])
						{
							donor_acceptor_HSP_ID[i].pop_back();
							donor_acceptor_HSP_ID[i].pop_back();
							UseLargestHSPInSegment(i, HSPs, donor_segments_pair[i].exon_seg_end-1, 
								donor_segments_pair[i].exon_seg_start);
						}
						else
						{
							donor_acceptor_HSP_ID[i].erase(donor_acceptor_HSP_ID[i].begin()+j);
						}

						if (i>0)
							newsites_donor_segments_pair.push_back(SegmentsInThreeBounds(donor_segments_pair[i-1].exon_seg_end, 
								donor_segments_pair[i-1].exon_seg_end, donor_segments_pair[i-1].intron_seg_end));
						newsites_donor_segments_pair.push_back(donor_segments_pair[i]);

						donor_erase_region_start = i;
						acceptor_erase_region_start = i;
						donor_erase_region_end = i;//same as start
						acceptor_erase_region_end = i+1;

					}
					else //j is not the first one, not the last one, must be the only segment that involved this HSP
					{
						int next_hsp_index = (*(hspID_to_hspIndex_map.find(donor_acceptor_HSP_ID[i][j+1]))).second;
						int next_hsp_start;
						int prev_hsp_index = (*(hspID_to_hspIndex_map.find(donor_acceptor_HSP_ID[i][j-1]))).second;
						int prev_hsp_end;
						if (start_site > 0)
						{
							next_hsp_start = HSPs[next_hsp_index]->HSP_start;
							prev_hsp_end = HSPs[prev_hsp_index]->HSP_end;
						}
						else
						{
							next_hsp_start = -HSPs[next_hsp_index]->HSP_end;
							prev_hsp_end = -HSPs[prev_hsp_index]->HSP_start;
						}

						//insert an extra segment
						donor_segments_pair.insert(donor_segments_pair.begin()+i+1, //donor_segments_pair.begin()+i+1-num_segment_erased, //no need, since num_segment_erased must be 0 here
							SegmentsInThreeBounds(next_hsp_start, donor_segments_pair[i].exon_seg_end, 
							donor_segments_pair[i].intron_seg_end));
						
						vector<int> newseg_hspIDs;
						j++;						
						while (j<cur_size-1)
						{
							donor_acceptor_HSP_ID[i].pop_back();
							newseg_hspIDs.push_back(donor_acceptor_HSP_ID[i][j]);
							j++;
						}
						donor_acceptor_HSP_ID[i].pop_back();
						donor_acceptor_HSP_ID[i].pop_back();

						donor_acceptor_HSP_ID.insert(donor_acceptor_HSP_ID.begin()+i+1, newseg_hspIDs);
						UseLargestHSPInSegment(i+1, HSPs, donor_segments_pair[i+1].exon_seg_end-1, 
							donor_segments_pair[i+1].exon_seg_start);

						donor_segments_pair[i].exon_seg_end = prev_hsp_end+1;
						donor_segments_pair[i].intron_seg_end = next_hsp_start-1;
						UseLargestHSPInSegment(i, HSPs, donor_segments_pair[i].exon_seg_end-1, 
							donor_segments_pair[i].exon_seg_start);

						//get newsites_donor_segments_pair here (3 segments if applicable)
						if (i>0)
							newsites_donor_segments_pair.push_back(SegmentsInThreeBounds(donor_segments_pair[i-1].exon_seg_end, 
								donor_segments_pair[i-1].exon_seg_end, donor_segments_pair[i-1].intron_seg_end));
						newsites_donor_segments_pair.push_back(donor_segments_pair[i]);
						newsites_donor_segments_pair.push_back(donor_segments_pair[i+1]);
						//newsites_donor_segments_pair_last = true;

						donor_erase_region_start = donor_erase_region_end = i;
						acceptor_erase_region_start = acceptor_erase_region_end = i;

						hsp_in_mid_segment = true;
						break;
					}
				}
			}
		}

		if (seg_check_done || hsp_in_mid_segment)
			break;

		if (found && !cur_found)
		{
			newsites_donor_segments_pair.push_back(SegmentsInThreeBounds(donor_segments_pair[i-num_segment_erased].exon_seg_start, 
				donor_segments_pair[i-num_segment_erased].exon_seg_end, donor_segments_pair[i-num_segment_erased].exon_seg_end));
			break;
		}

		i++;
	}



	vector<int>::iterator id_it = donor_acceptor_HSP_ID_to_be_erased.begin();

	j = 0;

	for (; id_it != donor_acceptor_HSP_ID_to_be_erased.end(); id_it++)

	{

		int k = (*id_it) - j;

		donor_acceptor_HSP_ID.erase(donor_acceptor_HSP_ID.begin()+k);

		j++; //j records the number of segments that's already erased

	}

#ifdef DEBUG
	dataMgr.outFile << "after adjusted segments:" << "\n";
	OutputDonorAcceptorWithAllInfo(start_site > 0);
#endif

	donor_erase_region_end++;

	acceptor_erase_region_end++;

	EraseSites_ByRegions(donors, donor_region_index, donor_tail, donor_HSP_ID, donor_prev_2nt, donor_erase_region_start, donor_erase_region_end);

	EraseSites_ByRegions(acceptors, acceptor_region_index, acceptor_head, acceptor_HSP_ID, acceptor_next_2nt, acceptor_erase_region_start, acceptor_erase_region_end);

#ifdef DEBUG
	dataMgr.outFile << "after erase affected sites:" << "\n";
	OutputDonorAcceptorWithAllInfo(start_site > 0);
#endif


	if (j > 0) //there're segments got erased, need to decrement the region_index for sites after the erased segments

	{

		IncrementSiteRegions(donor_region_index, donor_erase_region_end, -j);

		IncrementSiteRegions(acceptor_region_index, acceptor_erase_region_end, -j);

	}

	else

	{

		if (hsp_in_mid_segment) //had one extra segment, need to increment the region_index!

		{

			IncrementSiteRegions(donor_region_index, donor_erase_region_end, 1);

			IncrementSiteRegions(acceptor_region_index, acceptor_erase_region_end, 1);

		}

	}

#ifdef DEBUG
	dataMgr.outFile << "after adjusted site region indexes:" << "\n";
	OutputDonorAcceptorWithAllInfo(start_site > 0);
#endif


	//finally, insert new sites!	
	if (donor_erase_region_start < acceptor_erase_region_start)
		start_region_index = donor_erase_region_start < 0 ? 0 : donor_erase_region_start;
	else
		start_region_index = (donor_erase_region_start - 1) < 0 ? 0 : (donor_erase_region_start - 1);
/*	if (donor_erase_region_start < acceptor_erase_region_start)
	{
		i = donor_erase_region_start;
		if (i>=0)
		{
			start_region_index = i;
			newsites_donor_segments_pair.push_back(SegmentsInThreeBounds(donor_segments_pair[i].exon_seg_start, 
				donor_segments_pair[i].exon_seg_end, donor_segments_pair[i].intron_seg_end));
			i++;
		}
		else
		{
			i = 0;
			start_region_index = i;
		}
		if (donor_erase_region_end < acceptor_erase_region_end)
			newsites_donor_segments_pair.push_back(SegmentsInThreeBounds(donor_segments_pair[i].exon_seg_start, 
				donor_segments_pair[i].exon_seg_end, donor_segments_pair[i].exon_seg_end));
		else
			newsites_donor_segments_pair.push_back(SegmentsInThreeBounds(donor_segments_pair[i].exon_seg_start, 
				donor_segments_pair[i].exon_seg_end, donor_segments_pair[i].intron_seg_end));
	}
	else
	{
		i = donor_erase_region_start - 1;
		if (i >= 0)
		{
			start_region_index = i;
			newsites_donor_segments_pair.push_back(SegmentsInThreeBounds(donor_segments_pair[i].exon_seg_end-1,
				donor_segments_pair[i].exon_seg_end, donor_segments_pair[i].intron_seg_end));
			i++;
		}
		else
		{
			i = 0;
			start_region_index = i;
		}
		newsites_donor_segments_pair.push_back(SegmentsInThreeBounds(donor_segments_pair[i].exon_seg_start,
			donor_segments_pair[i].exon_seg_end, donor_segments_pair[i].intron_seg_end));
		i++;
		if (donor_erase_region_end < acceptor_erase_region_end)
			newsites_donor_segments_pair.push_back(SegmentsInThreeBounds(donor_segments_pair[i].exon_seg_start, 
				donor_segments_pair[i].exon_seg_end, donor_segments_pair[i].exon_seg_end));
		else
			newsites_donor_segments_pair.push_back(SegmentsInThreeBounds(donor_segments_pair[i].exon_seg_start, 
				donor_segments_pair[i].exon_seg_end, donor_segments_pair[i].intron_seg_end));
	}
*/
	vector<int> newsites_acceptors, newsites_donors, newsites_acceptor_region_index, newsites_donor_region_index;
	vector<int> newsites_acceptor_head, newsites_donor_tail, newsites_acceptor_HSP_ID, newsites_donor_HSP_ID;
	vector<string> newsites_acceptor_next_2nt, newsites_donor_prev_2nt;
	GetNewSites(newsites_donor_segments_pair, start_region_index, HSPs, chr_seq, 
		newsites_acceptors, newsites_donors, newsites_acceptor_region_index, newsites_donor_region_index, 
		newsites_acceptor_head, newsites_donor_tail, newsites_acceptor_HSP_ID, newsites_donor_HSP_ID, 
		newsites_acceptor_next_2nt, newsites_donor_prev_2nt);

	if (!newsites_acceptors.empty())
		InsertAccSites_ByRegions(newsites_acceptors, newsites_acceptor_region_index, newsites_acceptor_head, 
			newsites_acceptor_HSP_ID, newsites_acceptor_next_2nt, newsites_acceptor_region_index.front());

	if (!newsites_donors.empty())
		InsertDonSites_ByRegions(newsites_donors, newsites_donor_region_index, newsites_donor_tail, 
			newsites_donor_HSP_ID, newsites_donor_prev_2nt, newsites_donor_region_index.front());

#ifdef DEBUG
	dataMgr.outFile << "after inserted new sites:" << "\n";
	OutputDonorAcceptorWithAllInfo(start_site > 0);
#endif

	//finally do the post processing stuff on splice sites! (make sure there's no odd sites etc.)
/*	int first_start, last_end;
	if (start_site > 0)
	{
		first_start = HSPs[HSPs.size()-1]->HSP_start;
		last_end = HSPs[0]->HSP_end;
	}
	else
	{
		first_start = -HSPs[0]->HSP_end;
		last_end = -HSPs[HSPs.size()-1]->HSP_start;
	}
	PostProcessSpliceSites(first_start, last_end, chr_seq);
*/
	PostProcessSpliceSites(start_site, end_site, chr_seq);

	return first_hsp_erased;
}



void ACCP_DONR_graph::IncrementSiteRegions(vector<int>& site_regions, int region_no, int inc_count)

{

	int i = 0;

	while (i<site_regions.size())

	{

		if (site_regions[i] >= region_no)

			break;

		i++;

	}



	if (i<site_regions.size())

	{

		for (int j = i; j<site_regions.size(); j++)

			site_regions[j] += inc_count;

	}

	else

	{

#ifdef DEBUG

		dataMgr.outFile << "did not find anything with region_no >= site_region_no " << region_no << "\n";

#endif

	}

}







void ACCP_DONR_graph::GetNewSites(vector<SegmentsInThreeBounds>& newsites_donor_segments_pair, 

								  int start_region_index, 

									 vector<HSP_Gene_Pair*>& HSPs, 
									 //pair<int, char*>& chr_seq, 
									 vector<string>& chr_seq, 
									 vector<int>& newsites_acceptors, vector<int>& newsites_donors, 

									 vector<int>& newsites_acceptor_region_index, vector<int>& newsites_donor_region_index, 

									 vector<int>& newsites_acceptor_head, vector<int>& newsites_donor_tail, 

									 vector<int>& newsites_acceptor_HSP_ID, vector<int>& newsites_donor_HSP_ID, 

									 vector<string>& newsites_acceptor_next_2nt, vector<string>& newsites_donor_prev_2nt)

{

#ifdef DEBUG

	dataMgr.outFile << "GetNewSites...newsites_segments:" << "\n";

	for (int k=0; k<newsites_donor_segments_pair.size(); k++)

		dataMgr.outFile << newsites_donor_segments_pair[k] << "\n";

	dataMgr.outFile << "end of segments" << "\n";

	int t;

	vector<int>::iterator t_it;

#endif



	int search_start, search_end;

	

	int cur_seg = start_region_index;



	vector<int> accs, dons;

	vector< string > acc_next_2nt, don_prev_2nt;

	map<int, int> hsp_donor_stop_pos, hsp_acceptor_stop_pos;



	bool comp_start = false;



	vector<SegmentsInThreeBounds>::iterator seg_map_It;

	for (seg_map_It = newsites_donor_segments_pair.begin(); seg_map_It != newsites_donor_segments_pair.end(); 

		seg_map_It++, cur_seg++)

	{

		search_start = (*seg_map_It).exon_seg_end - 1; //(*seg_map_It).first-1;

		search_end = (*seg_map_It).exon_seg_start; //(*seg_map_It).second.first;



		//use search_start == search_end to skip unwanted segment regions!



		if (search_start > search_end)

		{

			GetDonorsAcceptors_nPerBorder(search_start, search_end, chr_seq, accs, dons, acc_next_2nt, don_prev_2nt);



			//Added: remove any acceptor beyond search_start, or any donor beyond search_end?

			vector<int>::reverse_iterator check_acc_it;

			while (!accs.empty()) 

			{

				check_acc_it = accs.rbegin();

				if (*check_acc_it > search_start)

					accs.pop_back();

				else

					break;

			}



#ifdef DEBUG		



		dataMgr.outFile << "now: searched " << search_start << "-" << search_end << "\n";

		dataMgr.outFile << "accs:" << "\n";

		

		t=0;

		for (t_it = accs.begin(); t_it != accs.end(); t_it++, t++)

			dataMgr.outFile << *t_it << ", " << acc_next_2nt[t] << ";";

		dataMgr.outFile << "\n";



		dataMgr.outFile << "dons:" << "\n";

		t = 0;

		for (t_it = dons.begin(); t_it != dons.end(); t_it++, t++)

			dataMgr.outFile << *t_it << ", " << don_prev_2nt[t] << ";";

		dataMgr.outFile << "\n";



#endif



			if (comp_start)

			{

				CombineSpliceSites_S1(search_end, accs, newsites_acceptors, newsites_acceptor_region_index, 

					cur_seg, true, newsites_acceptor_head, 

					HSPs, chr_seq, hsp_donor_stop_pos, hsp_acceptor_stop_pos, 

					newsites_acceptor_HSP_ID, newsites_donor_HSP_ID, acc_next_2nt, newsites_acceptor_next_2nt);



#ifdef DEBUG		



		dataMgr.outFile << "after combine: searched " << search_start << "-" << search_end << "\n";



		//PrintDonorAcceptor();



#endif

			}

			else //first time

			{

				accs.clear();

				acc_next_2nt.clear();

				comp_start = true;

			}



		}



		search_start++;

		search_end = (*seg_map_It).intron_seg_end;



		if (search_start < search_end) //now this is the gap region

		{

			GetDonorsAcceptors_nPerBorder(search_start, search_end, chr_seq, accs, dons, acc_next_2nt, don_prev_2nt);



			//Added: remove any acceptor beyond search_start, or any donor beyond search_end?

			vector<int>::reverse_iterator check_don_it;

			while (!dons.empty()) 

			{

				check_don_it = dons.rbegin();

				if (*check_don_it > search_end)

					dons.pop_back();

				else

					break;

			}



#ifdef DEBUG		



		dataMgr.outFile << "now: searched " << search_start << "-" << search_end << "\n";

		dataMgr.outFile << "accs:" << "\n";



		t=0;

		for (t_it = accs.begin(); t_it != accs.end(); t_it++,t++)

			dataMgr.outFile << *t_it << ", " << acc_next_2nt[t] << ";";

		dataMgr.outFile << "\n";



		dataMgr.outFile << "dons:" << "\n";

		t=0;

		for (t_it = dons.begin(); t_it != dons.end(); t_it++,t++)

			dataMgr.outFile << *t_it << ", " << don_prev_2nt[t] << ";";

		dataMgr.outFile << "\n";



#endif

			if (comp_start)

			{

				CombineSpliceSites_S1(search_start-1, dons, newsites_donors, newsites_donor_region_index, 

					cur_seg, false, newsites_donor_tail, 

					HSPs, chr_seq, hsp_donor_stop_pos, hsp_acceptor_stop_pos, 

					newsites_acceptor_HSP_ID, newsites_donor_HSP_ID, don_prev_2nt, newsites_donor_prev_2nt);



#ifdef DEBUG



		dataMgr.outFile << "after combine: searched " << search_start << "-" << search_end << "\n";



		//PrintDonorAcceptor();



#endif

			}

			else

			{

				dons.clear();

				don_prev_2nt.clear();

				comp_start = true;

			}

		}

	}



}



void ACCP_DONR_graph::InsertAccSites_ByRegions(vector<int>& sites, vector<int>& site_regions, 

								 vector<int>& site_head_tail, vector<int>& site_HSP_ID, vector<string>& site_2nt, 

								 int region_start_insert)

{

#ifdef DEBUG

	dataMgr.outFile << "InsertAccSites_ByRegions..." << "\n";

#endif



	//if (sites.empty())

	//	return;



	int i = 0;

	while (i<acceptor_region_index.size())

	{

		if (acceptor_region_index[i] < region_start_insert)

			i++;

		else

			break;

	}

	

	acceptors.insert(acceptors.begin()+i, sites.begin(), sites.end());

	acceptor_region_index.insert(acceptor_region_index.begin()+i, site_regions.begin(), site_regions.end());

	acceptor_head.insert(acceptor_head.begin()+i, site_head_tail.begin(), site_head_tail.end());

	acceptor_HSP_ID.insert(acceptor_HSP_ID.begin()+i, site_HSP_ID.begin(), site_HSP_ID.end());

	acceptor_next_2nt.insert(acceptor_next_2nt.begin()+i, site_2nt.begin(), site_2nt.end());



}



void ACCP_DONR_graph::InsertDonSites_ByRegions(vector<int>& sites, vector<int>& site_regions, 

								 vector<int>& site_head_tail, vector<int>& site_HSP_ID, vector<string>& site_2nt, 

								 int region_start_insert)

{

#ifdef DEBUG

	dataMgr.outFile << "InsertDonSites_ByRegions..." << "\n";

#endif



	//if (sites.empty())

	//	return;



	int i = 0;

	while (i<donor_region_index.size())

	{

		if (donor_region_index[i] < region_start_insert)

			i++;

		else

			break;

	}

	

	donors.insert(donors.begin()+i, sites.begin(), sites.end());

	donor_region_index.insert(donor_region_index.begin()+i, site_regions.begin(), site_regions.end());

	donor_tail.insert(donor_tail.begin()+i, site_head_tail.begin(), site_head_tail.end());

	donor_HSP_ID.insert(donor_HSP_ID.begin()+i, site_HSP_ID.begin(), site_HSP_ID.end());

	donor_prev_2nt.insert(donor_prev_2nt.begin()+i, site_2nt.begin(), site_2nt.end());



}

void ACCP_DONR_graph::OutputDonorAcceptorWithAllInfo(bool isPosStrand)
{
	dataMgr.outFile << "===============================\nOutputDonorAcceptorWithAllInfo\n";

	int i;

		if (isPosStrand)



		{



			dataMgr.outFile << "acceptor\t0\t" << start_site << "-" << start_site+100 << "\n";



			dataMgr.outFile << "donor\t0\t" << end_site << "-" << end_site+100 << "\n";



			//output acceptor/donors between start and end



			for (i=0; i<acceptors.size(); i++)

			{

				dataMgr.outFile << "acceptor\t" << i+1 << "\t" << acceptors[i] << "-" << acceptors[i]+100 //<< "\n";

					<< "; region[" << acceptor_region_index[i] << "]" << " head:" << acceptor_head[i];

				//if (SPLICE_SEGMENT_VERSION == 1)

					dataMgr.outFile << "; hsp_ID:" << acceptor_HSP_ID[i];

				dataMgr.outFile << "\n";

			}

			for (i=0; i<donors.size(); i++)

			{

				dataMgr.outFile << "donor\t" << i+1 << "\t" << donors[i] << "-" << donors[i]+100 //<< "\n";

					<< "; region[" << donor_region_index[i] << "]" << " tail:" << donor_tail[i];

				//if (SPLICE_SEGMENT_VERSION == 1)

					dataMgr.outFile << "; hsp_ID:" << donor_HSP_ID[i];

				dataMgr.outFile << "\n";

			}



		}



		else



		{



			dataMgr.outFile << "acceptor\t0\t" << -start_site << "-" << -start_site+100 << "\n";



			dataMgr.outFile << "donor\t0\t" << -end_site << "-" << -end_site+100 << "\n";



			for (i=0; i<acceptors.size(); i++)

			{

				dataMgr.outFile << "acceptor\t" << i+1 << "\t" << -acceptors[i] << "-" << -acceptors[i]+100 //<< "\n";

					<< "; region[" << acceptor_region_index[i] << "]" << " head:" << acceptor_head[i];

				//if (SPLICE_SEGMENT_VERSION == 1)

					dataMgr.outFile << "; hsp_ID:" << acceptor_HSP_ID[i];

				dataMgr.outFile << "\n";

			}



			for (i=0; i<donors.size(); i++)

			{

				dataMgr.outFile << "donor\t" << i+1 << "\t" << -donors[i] << "-" << -donors[i]+100 //<< "\n";

					<< "; region[" << donor_region_index[i] << "]" << " tail:" << donor_tail[i];

				//if (SPLICE_SEGMENT_VERSION == 1)

					dataMgr.outFile << "; hsp_ID:" << donor_HSP_ID[i];

				dataMgr.outFile << "\n";

			}



		}

	dataMgr.outFile << "acceptor_region_index:" << "\n";
	for (i=0; i<acceptor_region_index.size(); i++)
		dataMgr.outFile << acceptor_region_index[i] << ",";
	dataMgr.outFile << "acceptor_region_index DONE" << "\n";

	dataMgr.outFile << "donor_region_indexes:" << "\n";
	for (i=0; i<donor_region_index.size(); i++)
		dataMgr.outFile << donor_region_index[i] << ",";
	dataMgr.outFile << "donor_region_indexes DONE" << "\n";

	dataMgr.outFile << "OutputDonorAcceptorWithAllInfo DONE\n===============================\n";
}

void ACCP_DONR_graph::UseOldExon(vector<ExonSiteInfo>& temp_sites, vector<ExonSiteInfo>& additional_temp_sites, 
								 vector<HSP_Gene_Pair>& temp_HSPs, 
								 int target_site, int query_site, int target_end_site, int query_end_site, 
								 int blast_hsp_index, vector<HSP_Gene_Pair*>& blastHSPs)
{
	int j;

	temp_sites.push_back( ExonSiteInfo(target_site, query_site, -1, 0)); //hsp_index doesn't matter?
	temp_sites.push_back( ExonSiteInfo(target_end_site, query_end_site, -1, 
		((target_end_site - target_site + 1)% 3 )));

	additional_temp_sites.push_back( ExonSiteInfo(target_site, query_site, -1, 0)); //hsp_index doesn't matter?
	additional_temp_sites.push_back( ExonSiteInfo(target_end_site, query_end_site, -1, 
		((target_end_site - target_site + 1)% 3 )));

	temp_HSPs.clear(); //now update temp_HSPs

	if (target_site < 0)
	{
		for (j=0; j<=blast_hsp_index; j++)
			temp_HSPs.push_back(*(blastHSPs[j]));
	}
	else
	{
		//for (j=blastHSPs.size()-1; j>=blast_hsp_index; j--)
		//temp_HSPs.insert(temp_HSPs.begin(), blastHSPs[j]);
		for (j=blast_hsp_index; j<blastHSPs.size(); j++)
			temp_HSPs.push_back(*(blastHSPs[j]));
	}		
	
}

void ACCP_DONR_graph::UseOldExon_Mid(vector<ExonSiteInfo>& temp_sites, vector<ExonSiteInfo>& additional_temp_sites, 
								 vector<HSP_Gene_Pair>& temp_HSPs, 
									 ExonSiteInfo& cur_exon_start, ExonSiteInfo& cur_exon_end, 
									 ExonSiteInfo& prev_exon_start, 
									 int cur_start_blast_hsp_index, int cur_end_blast_hsp_index, 
									 vector<HSP_Gene_Pair*>& blastHSPs)
{
#ifdef DEBUG
	dataMgr.outFile << "UseOldExon_Mid (prev_exon_start:" << prev_exon_start << ") cur_start_blast_hsp_index:"
		<< cur_start_blast_hsp_index << "; cur_end_blast_hsp_index:" << cur_end_blast_hsp_index 
		<< "; blastHSPs.size():" << blastHSPs.size() 
		<< "\n";
#endif

	int j;
							temp_sites.push_back(cur_exon_start);
							temp_sites.push_back(cur_exon_end);

							additional_temp_sites.push_back(cur_exon_start);
							additional_temp_sites.push_back(cur_exon_end);

							temp_HSPs.clear();
							if (prev_exon_start.first > 0)
							{
								for (j=cur_end_blast_hsp_index; j<=cur_start_blast_hsp_index; j++)
									temp_HSPs.push_back(*(blastHSPs[j]));
							}
							else
							{
								for (j=cur_start_blast_hsp_index; j<=cur_end_blast_hsp_index; j++)
									temp_HSPs.push_back(*(blastHSPs[j]));
							}
#ifdef DEBUG
	dataMgr.outFile << "UseOldExon_Mid DONE" << "\n";
#endif

}

bool ACCP_DONR_graph::PostProcessSpliceSites(int first_segment_start, int last_segment_end, 
											 //pair<int, char*>& chr_seq)
											 vector<string>& chr_seq)
{
	if (IsDonorAcceptorEmpty())
		return true;

	if (!IsSorted(acceptors))
	{
		cout << "acceptors not sorted" << "\n";
		PrintDonorAcceptorHSPID();
		exit(-1);
		return false;
	}

	if (!IsSorted(donors))
	{
		cout << "donors not sorted" << "\n";
		PrintDonorAcceptorHSPID();
		exit(-1);
		return false;
	}

/*
	if (!IsUnique(acceptors))
	{
		cout << "acceptors not unique" << "\n";
		PrintDonorAcceptor();
		exit(-1);
		return false;
	}
	if (!IsUnique(donors))
	{
		cout << "donors not unique" << "\n";
		PrintDonorAcceptor();
		exit(-1);
		return false;
	}
*/

#ifdef DEBUG
	int i;
	for (i=0; i<acceptors.size(); i++)
		dataMgr.outFile << "acceptor\t" << i+1 << "\t" << acceptors[i] << "-" << acceptors[i]+100 << "\n";
	for (i=0; i<donors.size(); i++)
		dataMgr.outFile << "donor\t" << i+1 << "\t" << donors[i] << "-" << donors[i]+100 << "\n";
#endif


	//filter out acceptors before 1st donor, and donors after last acceptor;
	//AND also:
	//remove acceptor/donor before first_segment_start or after last_segment_end!
	vector<int>::iterator setIt;
	setIt = donors.begin(); 
	while (setIt != donors.end() && (*setIt <= start_site ) )
		setIt++;
	int fst_donor=0;
	if (setIt != donors.end())
		fst_donor = *setIt; //first donor after start_site

	if ( acceptors.front() >= end_site) //if all acceptors >= end, all must be erased
	{
		acceptors.clear();
		donors.clear();
		return true;
	}

	setIt = acceptors.end();
	setIt--;

	while (*setIt >= end_site ) //here, there must be at least 1 acceptor < end, so setIt won't go beyond acceptors.begin()
		setIt--;
	int last_acceptor = *setIt; //last acceptor before end_site

#ifdef DEBUG
	dataMgr.outFile << "last_acceptor to be kept: " << last_acceptor << "\n";
#endif

	int first_start = first_segment_start < fst_donor ? fst_donor : first_segment_start; //get whichever is later/bigger
	int last_end = last_segment_end < last_acceptor ? last_segment_end : last_acceptor; //get whichever is earlier/smaller

//	if (SPLICE_SEGMENT_VERSION == 1)
//	{
		EraseSites(donors, donor_region_index, donor_tail, donor_HSP_ID, donor_prev_2nt, first_start, true);
		EraseSites(acceptors, acceptor_region_index, acceptor_head, acceptor_HSP_ID, acceptor_next_2nt, first_start, true);
		EraseSites(donors, donor_region_index, donor_tail, donor_HSP_ID, donor_prev_2nt, last_end, false);
		EraseSites(acceptors, acceptor_region_index, acceptor_head, acceptor_HSP_ID, acceptor_next_2nt, last_end, false);

//	}
//	else
//	{
//		EraseSites(donors, donor_region_index, donor_tail, donor_prev_2nt, first_segment_start, true);
//		EraseSites(acceptors, acceptor_region_index, acceptor_head, acceptor_next_2nt, first_segment_start, true);
//		EraseSites(donors, donor_region_index, donor_tail, donor_prev_2nt, last_segment_end, false);
//		EraseSites(acceptors, acceptor_region_index, acceptor_head, acceptor_next_2nt, last_segment_end, false);
//	}

	if (IsDonorAcceptorEmpty())
		return true;

	//UPDATED: check for in-frame stop codon within the MIN_EXON_LENGTH area of acceptor/donor
	int last_acc_region_index = acceptor_region_index.back();
	string testStr, codonStr;

	int cur_site, site_index;

	vector<int> acceptors_to_be_erased, donors_to_be_erased;
	for (site_index = 0; site_index < acceptors.size(); site_index++)
	{
		//do not check last bunch of acceptors, since there's no min_exon_len for last exon
		if (acceptor_region_index[site_index] == last_acc_region_index) 
			break;

		cur_site = acceptors[site_index];
		if (cur_site > 0 && LenOfStrVec(chr_seq)+chromosome_start_pos-1 > cur_site + MIN_INTERNAL_EXON_LEN - 2)
			GetSubstrFromVecStrs(chr_seq, true, cur_site-1, MIN_INTERNAL_EXON_LEN, testStr);
		else
			if (cur_site < 0 && -cur_site-MIN_INTERNAL_EXON_LEN+1 >= chromosome_start_pos)
				GetSubstrFromVecStrs_NegRev(chr_seq, false, -cur_site-MIN_INTERNAL_EXON_LEN, MIN_INTERNAL_EXON_LEN, testStr);
			else //this acceptor does not have enough DNA length to make an internal exon
			{
				acceptors_to_be_erased.push_back(site_index);
				continue;
			}

		StrToLower(testStr);
		int tmp_start_pos = acceptor_head[site_index];
		bool has_stop = false;
		while (!has_stop && tmp_start_pos+3 <= testStr.length())
		{
			codonStr = testStr.substr(tmp_start_pos, 3);
			map<string, char>::iterator it = DNA_CODON_TBL_SL.find(codonStr);
			if (it != DNA_CODON_TBL_SL.end())
			{
				if ((*it).second == '*')
					has_stop = true;
			}

			tmp_start_pos += 3;
		}

		if (has_stop)
			acceptors_to_be_erased.push_back(site_index);
	}
	for (site_index = 0; site_index < donors.size(); site_index++)
	{
		//do not check first bunch of donors, since there's no min_exon_len for first exon
		if (donor_region_index[site_index] == 0) 
			continue;

		cur_site = donors[site_index];
		if (cur_site > 0 && cur_site-MIN_INTERNAL_EXON_LEN+1 >= chromosome_start_pos)
			GetSubstrFromVecStrs(chr_seq, true, cur_site-MIN_INTERNAL_EXON_LEN, 
				MIN_INTERNAL_EXON_LEN-donor_tail[site_index], testStr);
		else
			if (cur_site < 0 && LenOfStrVec(chr_seq)+chromosome_start_pos-1 > -cur_site-2+MIN_INTERNAL_EXON_LEN)
				GetSubstrFromVecStrs_NegRev(chr_seq, false, -cur_site-1+donor_tail[site_index], 
					MIN_INTERNAL_EXON_LEN-donor_tail[site_index], testStr);
			else
			{
				donors_to_be_erased.push_back(site_index);
				continue;
			}

		StrToLower(testStr);
		
		int tmp_start_pos = testStr.length();
		bool has_stop = false;
		while (!has_stop && tmp_start_pos >= 3)
		{
			codonStr = testStr.substr(tmp_start_pos-3, 3);
			map<string, char>::iterator it = DNA_CODON_TBL_SL.find(codonStr);
			if (it != DNA_CODON_TBL_SL.end())
			{
				if ((*it).second == '*')
					has_stop = true;
			}
			tmp_start_pos -= 3;
		}
		if (has_stop)
			donors_to_be_erased.push_back(site_index);
	}
	EraseSpecificSites(donors, donor_region_index, donor_tail, donor_HSP_ID, donor_prev_2nt, donors_to_be_erased);
	EraseSpecificSites(acceptors, acceptor_region_index, acceptor_head, acceptor_HSP_ID, acceptor_next_2nt, acceptors_to_be_erased);
	if (IsDonorAcceptorEmpty())
		return true;

	//remove acceptor/donor before first_segment_start or after last_segment_end!
/*	if (SPLICE_SEGMENT_VERSION == 1)
	{
		EraseSites(donors, donor_region_index, donor_tail, donor_HSP_ID, first_segment_start, true);
		EraseSites(acceptors, acceptor_region_index, acceptor_head, acceptor_HSP_ID, first_segment_start, true);
		EraseSites(donors, donor_region_index, donor_tail, donor_HSP_ID, last_segment_end, false);
		EraseSites(acceptors, acceptor_region_index, acceptor_head, acceptor_HSP_ID, last_segment_end, false);
	}
	else
	{
		EraseSites(donors, donor_region_index, donor_tail, first_segment_start, true);
		EraseSites(acceptors, acceptor_region_index, acceptor_head, first_segment_start, true);
		EraseSites(donors, donor_region_index, donor_tail, last_segment_end, false);
		EraseSites(acceptors, acceptor_region_index, acceptor_head, last_segment_end, false);
	}
*/

/*
	//filter out acceptors before 1st donor, and donors after last acceptor;
	vector<int>::iterator setIt, regionIt;
	setIt = donors.begin(); 
	regionIt = donor_region_index.begin();
	while (setIt != donors.end() && (*setIt <= start_site ) )
	{
		setIt++;
		regionIt++;
	}
	int fst_donor=0;
	if (setIt != donors.end())
		fst_donor = *setIt; //first donor after start_site
	if (setIt != donors.begin())
	{
		donors.erase(donors.begin(), setIt);
		donor_region_index.erase(donor_region_index.begin(), regionIt);
	}

	if (IsDonorAcceptorEmpty())
		return true;
*/


/*	setIt=acceptors.begin();
	regionIt = acceptor_region_index.begin();
	while (setIt != acceptors.end())
	{
		if (*setIt > fst_donor)
			break;
		setIt++;
		regionIt++;
	}
	if (setIt != acceptors.begin())
	{
		acceptors.erase(acceptors.begin(), setIt);
		acceptor_region_index.erase(acceptor_region_index.begin(), regionIt);
	}

	if (IsDonorAcceptorEmpty())
		return true;
*/

/*	if ( acceptors.front() >= end_site) //if all acceptors >= end, all must be erased
	{
		acceptors.clear();
		donors.clear();
		return true;
	}
	setIt = acceptors.end();
	setIt--;
	regionIt = acceptor_region_index.end();
	regionIt--;
	while (*setIt >= end_site ) //here, there must be at least 1 acceptor < end, so setIt won't go beyond acceptors.begin()
	{
		setIt--;
		regionIt--;
	}
	int last_acceptor = *setIt; //last acceptor before end_site
	setIt++;
	regionIt++;
	if (setIt != acceptors.end())
	{
		acceptors.erase(setIt, acceptors.end());
		acceptor_region_index.erase(regionIt, acceptor_region_index.end());
	}
	if (IsDonorAcceptorEmpty())
		return true;
*/	

/*	setIt = donors.end();
	regionIt = donor_region_index.end();
	while (setIt != donors.begin())
	{
		setIt--;
		regionIt--;
		if (*setIt < last_acceptor)
			break;
	}
	if (*setIt < last_acceptor)
	{
		setIt++;
		regionIt++;
	}
	if (setIt != donors.end())
	{
		donors.erase(setIt, donors.end());
		donor_region_index.erase(regionIt, donor_region_index.end());
	}

	if (IsDonorAcceptorEmpty())
		return true;
*/

	return true;
}


void ACCP_DONR_graph::PrintGFFExons(bool isPosStrand, int exon_start, int exon_end, int exon_count, int alt_count)
{
	//make sure exon_start < exon_end
	if (exon_start > exon_end)
	{
		dataMgr.outFile << "stop EXON (" << exon_count << ")?" << exon_start << "to" << exon_end << "\n";
		//return;
		int e = exon_start;
		exon_start = exon_end;
		exon_end = e;
	}

	if (OUTPUT_GFF)
	//dataMgr.gff_os << dataMgr.cur_chr_name << "\tgenBlastG" << USER_ID << "\texon\t";
		dataMgr.gff_str << dataMgr.cur_chr_name << "\tgenBlastG" << USER_ID << "\tcoding_exon\t";

	if (isPosStrand)
	{
		dataMgr.outFile << "EXON\t" << exon_count << "\t" << exon_start  
			<< "-" << exon_end << "\n";
		if (OUTPUT_GFF)
		{
		//dataMgr.gff_os << exon_start << "\t" << exon_end << "\t.\t+\t.\t";
			dataMgr.gff_str << exon_start << "\t" << exon_end << "\t.\t+\t.\t";
		}
		if (final_start_site > exon_start )
			final_start_site = exon_start;
		if (final_end_site < exon_end)
			final_end_site = exon_end;
	}
	else
	{
		dataMgr.outFile << "EXON\t" << exon_count << "\t" << -exon_end 
			<< "-" << -exon_start << "\n";
		if (OUTPUT_GFF)
		{
		//dataMgr.gff_os << -exon_end << "\t" << -exon_start << "\t.\t-\t.\t";
			dataMgr.gff_str << -exon_end << "\t" << -exon_start << "\t.\t-\t.\t";
		}
		if (final_start_site > -exon_end )
			final_start_site = -exon_end;
		if (final_end_site < -exon_start)
			final_end_site = -exon_start;
	}
	if (OUTPUT_GFF)
	//dataMgr.gff_os << "ID=" << dataMgr.cur_gene_id << "-A" << alt_count << "-E" << exon_count 
	//	<< ";Parent=" << dataMgr.cur_gene_id << "-A" << alt_count << "\n";
		dataMgr.gff_str << "ID=" << dataMgr.cur_gene_id << "-A" << alt_count << "-E" << exon_count 
			<< ";Parent=" << dataMgr.cur_gene_id << "-A" << alt_count << "\n";
}

bool ACCP_DONR_graph::CheckSglExon(vector<HSP_Gene_Pair*>& HSPs, 
								   //pair<int, char*>& chr_seq, 
								   vector<string>& chr_seq, 
								   int start_frame, int end_frame, int cur_exon_start_site, 
								   vector< vector< ExonSiteInfo > >& all_alternative_acceptors, 
								   vector< vector< ExonSiteInfo > >& all_alternative_donors, 
								   bool& has_possible_exon, 
								   bool repair_only)
{
#ifdef DEBUG
	dataMgr.outFile << "start_frame:" << start_frame << ";end_frame:" << end_frame << "; repair_only:" << repair_only 
		<< "; end_site:" << end_site << "; cur_exon_start_site:" << cur_exon_start_site << "\n";
#endif
	
	//check frame (only for "repair_only" step?)
	if (repair_only && ((end_site - cur_exon_start_site + 1 - start_frame - end_frame) % 3 != 0))
	{
		has_possible_exon = false;
		all_alternative_acceptors.clear();
		all_alternative_donors.clear();
		return false;
	}

	//check stop
	string testStr, codonStr;
	if (start_site > 0)
		//GetSubstrFromVecStrs(chr_seq, true, start_site-1, end_site-start_site+1-end_frame, testStr);
		GetSubstrFromVecStrs(chr_seq, true, cur_exon_start_site+start_frame-1, end_site-cur_exon_start_site+1-end_frame-start_frame, testStr);
	else
		//GetSubstrFromVecStrs_NegRev(chr_seq, false, -end_site-1+end_frame, end_site-start_site+1-end_frame, testStr);
		GetSubstrFromVecStrs_NegRev(chr_seq, false, -end_site-1+end_frame, end_site-cur_exon_start_site+1-end_frame-start_frame, testStr);
	StrToLower(testStr);
#ifdef DEBUG
	dataMgr.outFile << "testStr:" << testStr << "\n";
#endif

	int stop_pos;
	//if (HasInFrameStop(testStr, false, stop_pos))
	if (HasInFrameStop(testStr, true, stop_pos, true))
	{
#ifdef DEBUG
		dataMgr.outFile << "end_site has in-frame stop from cur_exon_start_site: stop_pos:" << stop_pos << "\n";
#endif
		if (repair_only)
		{
			has_possible_exon = false;
			all_alternative_acceptors.clear();
			all_alternative_donors.clear();
			return false;
		}

		bool removed_hsp_is_new_hsp;
		if (RemoveHSPBeforeStop(HSPs, cur_exon_start_site+start_frame+stop_pos, chr_seq, removed_hsp_is_new_hsp, stop_pos > testStr.length()/2))
		{
			if (removed_hsp_is_new_hsp) //only happens during repair!
			{
				has_possible_exon = false;
			all_alternative_acceptors.clear();
			all_alternative_donors.clear();
				return false;
			}
			if (HSPs.size() > 0) //there's still HSPs left after removal
			{
#ifdef COMPUT_EXON_FULL_STEP_BACK
	#ifdef DEBUG
				dataMgr.outFile << "removed a HSP, still HSP left, now going back one full step" << "\n";
	#endif
#else
				//for now, just restart ComputeExons()
	#ifdef DEBUG								
				dataMgr.outFile << "encounted a stop, adjusted segments and go back (not full step back)" << "\n";
	#endif
#endif
				//all_alternative_acceptors.clear();
				//all_alternative_donors.clear();
			}
			else
			{
				has_possible_exon = false;
#ifdef DEBUG
				dataMgr.outFile << "no more HSPs, no possible exon\n";
#endif
			}
		}
		else
		{
#ifdef DEBUG
			dataMgr.outFile << "did not find the HSP to erase, no possible exon\n";
#endif
			has_possible_exon = false;
		}

		all_alternative_acceptors.clear();
		all_alternative_donors.clear();
								
		return false;	

	}

	return true;
}

void ACCP_DONR_graph::GetPredGeneSeq(vector<ExonSiteInfo>& temp_sites, 
									 //pair<int, char*>& chr_seq, 
									 vector<string>& chr_seq, 
						string& final_alignment, //predicted protein
						string& temp_DNA_string, 
						bool& found_stop) //predicted cDNA
{
	int frame=0;
	string left_chars;
	int prev_donor_site = 0;
	int exon_count = 1;

	int j;
	for (j=0; j<temp_sites.size(); j+=2)
	{
		int site_before_stop;
		if (!HasStopCodon(temp_sites[j].first, temp_sites[j+1].first, 
				chr_seq, frame, site_before_stop, final_alignment, prev_donor_site, left_chars, 
				false, dataMgr.cDNA_os, true, temp_DNA_string, j==temp_sites.size()-2, found_stop ))
		{
#ifdef DEBUG
			dataMgr.outFile << "Checking EXON\t" << exon_count << "\t" << temp_sites[j].first
				<< "-" << temp_sites[j+1].first << "\n";
#endif
			prev_donor_site = temp_sites[j+1].first;
			exon_count++;
		}
		else
		{
#ifdef DEBUG
			dataMgr.outFile << "Stop at EXON\t" << exon_count << "\t" << temp_sites[j].first
				<< "-" << site_before_stop << "\n";
#endif
			prev_donor_site = site_before_stop;
			//adjust the temp_sites here!
			for (int i = j+1; i<temp_sites.size(); i++)
				temp_sites.pop_back();
			temp_sites.push_back(ExonSiteInfo(prev_donor_site,0,0,0));
			break;
		}
	}
}

align_compare ACCP_DONR_graph::GetBetterPredGeneSeq(vector<ExonSiteInfo>& temp_sites, vector<ExonSiteInfo>& prev_temp_sites, 
						   //pair<int, char*>& chr_seq, 
						   vector<string>& chr_seq, 
						string& final_alignment, string& prev_alignment, //predicted protein
						string& final_DNA_string, string& prev_DNA_string, //predicted cDNA
						Input_Alignment& newAlign, Input_Alignment& prev_align, 
						bool& found_stop)
{
	bool found_stop_prev;
	GetPredGeneSeq(prev_temp_sites, chr_seq, prev_alignment, prev_DNA_string, found_stop_prev);
	float prev_align_pid;
	int prev_align_score;
#ifdef PID_BEFORE_SCORE
	prev_align_score = GetGlobalAlignment_PID_CompScore(query_seq, prev_alignment, prev_align, prev_align_pid, dataMgr.outFile);
#else
	prev_align_score = GetGlobalAlignment_scorematrix(query_seq, prev_alignment, prev_align, prev_align_pid, dataMgr.outFile);
#endif

	align_compare align_p(prev_align_score, prev_align_pid);
//#ifdef DEBUG
	dataMgr.outFile << "Second-to-end alignment(to be compared):\n" << prev_align << "\n";
	dataMgr.outFile << "Second-to-end PID:" << prev_align_pid << "\n";
//#endif

	GetPredGeneSeq(temp_sites, chr_seq, final_alignment, final_DNA_string, found_stop);
	float align_pid;
	int align_score;
#ifdef PID_BEFORE_SCORE
	align_score = GetGlobalAlignment_PID_CompScore(query_seq, final_alignment, newAlign, align_pid, dataMgr.outFile);
#else
	align_score = GetGlobalAlignment_scorematrix(query_seq, final_alignment, newAlign, align_pid, dataMgr.outFile);
#endif
	align_compare align_c(align_score, align_pid);
//#ifdef DEBUG
	dataMgr.outFile << "End alignment(to be compared):\n" << newAlign << "\n";
	dataMgr.outFile << "End PID:" << align_pid << "\n";
//#endif	

	//if (prev_align_score > align_score || (prev_align_score == align_score && prev_align_pid > align_pid))
	if (align_p > align_c)
	{
		final_alignment = prev_alignment;
		final_DNA_string = prev_DNA_string;
		//align_pid = prev_align_pid;
		//align_score = prev_align_score;
		align_c = align_p;
		newAlign = prev_align;
		temp_sites = prev_temp_sites;
		found_stop = found_stop_prev;
#ifdef DEBUG
		dataMgr.outFile << "use prev_temp_sites" << "\n";
		for (int i=0; i<temp_sites.size(); i+=2)
			dataMgr.outFile << temp_sites[i].first << " to " << temp_sites[i+1].first << "\n";
#endif
	}
	return align_c; //align_compare(align_score, align_pid); //the final better option is stored in final_alignment, final_DNA_string, temp_sites
}

bool ACCP_DONR_graph::RemoveHSPAndResetComputeExons(vector<HSP_Gene_Pair*>& HSPs, int stop_site, 
													//pair<int, char*>& chr_seq, 
													vector<string>& chr_seq, 
													bool remove_after_stop_preferred, 
													bool& has_possible_exon, 
													vector< vector< ExonSiteInfo > >& all_alternative_acceptors,
													vector< vector< ExonSiteInfo > >& all_alternative_donors)
{
	bool removed_hsp_is_new_hsp;
	if (RemoveHSPBeforeStop(HSPs, stop_site, chr_seq, removed_hsp_is_new_hsp, remove_after_stop_preferred))
	{
		if (removed_hsp_is_new_hsp) //only happens during repair!
		{
			has_possible_exon = false;
			return false;
		}

		if (HSPs.size() > 0) //there's still HSPs left after removal
		{
#ifdef COMPUT_EXON_FULL_STEP_BACK
	#ifdef DEBUG
			dataMgr.outFile << "removed a HSP, still HSP left, now going back one full step" << "\n";
	#endif
#else
			//for now, just restart ComputeExons()
	#ifdef DEBUG								
			dataMgr.outFile << "encounted a stop, adjusted segments and go back (not full step back)" << "\n";
	#endif
#endif
			all_alternative_acceptors.clear();
			all_alternative_donors.clear();
		}
		else
		{
			has_possible_exon = false;
#ifdef DEBUG
			dataMgr.outFile << "no more HSPs, no possible exon\n";
#endif
		}
	}
	else
	{
#ifdef DEBUG
		dataMgr.outFile << "did not find the HSP to erase, no possible exon\n";
#endif
		has_possible_exon = false;
	}
								
	return false;	
}



//only used when trying to separate donor-acceptor pairs with same highest local PIDs (based on existing HSPs)
//this function computes the actual optimal alignment between the two corresponding subsequences from target and query
//it returns the alignment score, and also updates the align_pid
int ACCP_DONR_graph::ComputeSpliceAlignment_SplitOL1_Optimal(vector<HSP_Gene_Pair*>& HSPs, 
															 //pair<int, char*>& chr_seq, 
															 vector<string>& chr_seq, 
															  IntronAndDonorTail& donor_acceptor_info, float& align_pid)
{
	//get target subsequence
	string donorStr, acceptorStr;
	if (start_site > 0)
	{
		if (donor_acceptor_info.first == 0) //signals gap center
		{
			GetSubstrFromVecStrs(chr_seq, true, HSPs[donor_acceptor_info.donor_hsp_index]->HSP_start-1, 
				HSPs[donor_acceptor_info.acceptor_hsp_index]->HSP_end - HSPs[donor_acceptor_info.donor_hsp_index]->HSP_start + 1, 
				donorStr);
		}
		else
		{
			GetSubstrFromVecStrs(chr_seq, true, HSPs[donor_acceptor_info.donor_hsp_index]->HSP_start-1, 
				donor_acceptor_info.first - HSPs[donor_acceptor_info.donor_hsp_index]->HSP_start + 1, donorStr);
			GetSubstrFromVecStrs(chr_seq, true, donor_acceptor_info.second-1, 
				HSPs[donor_acceptor_info.acceptor_hsp_index]->HSP_end - donor_acceptor_info.second + 1, acceptorStr);
			donorStr += acceptorStr;
		}
	}
	else
	{
		if (donor_acceptor_info.first == 0) //signals gap center
		{
			GetSubstrFromVecStrs_NegRev(chr_seq, false, HSPs[donor_acceptor_info.acceptor_hsp_index]->HSP_start-1, 
				HSPs[donor_acceptor_info.donor_hsp_index]->HSP_end - HSPs[donor_acceptor_info.acceptor_hsp_index]->HSP_start + 1, 
				donorStr);
		}
		else
		{
			GetSubstrFromVecStrs_NegRev(chr_seq, false, -donor_acceptor_info.first-1, 
				HSPs[donor_acceptor_info.donor_hsp_index]->HSP_end + donor_acceptor_info.first + 1, donorStr);
			GetSubstrFromVecStrs_NegRev(chr_seq, false, HSPs[donor_acceptor_info.acceptor_hsp_index]->HSP_start-1, 
				-donor_acceptor_info.second-HSPs[donor_acceptor_info.acceptor_hsp_index]->HSP_start + 1, acceptorStr);
			donorStr += acceptorStr;
		}
	}
	//now donorStr contains the entire target subsequence to be aligned
	StrToLower(donorStr);
	
	string tStr;
	DNA2AA(donorStr, 0, tStr);

	//get query subsequence
	string qStr = query_seq.substr(HSPs[donor_acceptor_info.donor_hsp_index]->gene_start-1, 
		HSPs[donor_acceptor_info.acceptor_hsp_index]->gene_end - HSPs[donor_acceptor_info.donor_hsp_index]->gene_start + 1);

	Input_Alignment align;
	//float align_pid;
#ifdef PID_BEFORE_SCORE
	return GetGlobalAlignment_PID_CompScore(qStr, tStr, align, align_pid, dataMgr.outFile);
#else
	return GetGlobalAlignment_scorematrix(qStr, tStr, align, align_pid, dataMgr.outFile);
#endif
}

bool ACCP_DONR_graph::GetHSPSpliceAlignment1(string& targetStr, string& queryStr, string& matchStr, 
	int cur_hsp_start, int cur_hsp_end, int gene_start, int gene_end, 
	int site_start, int site_end, //start and end of the splice segment to be obtained
	//the following FIVE are output from this function
	Input_Alignment& donor_align, int& align_start, int& align_end, char& border_query_aa_front, char& border_query_aa_end, 
	vector<string>& chr_seq) //border_query_aa records the char (protein amino acid, single letter) at the border (in case there's broken codon)
	//pair<int, char*>& chr_seq)
{
	//we only need query_align and match_align (two strings)
	//UPDATE: now we need all three alignment strings (target_align used to check for gaps of target, if there're gaps to be inserted to query, in ComputeSpliceAlignment)
	donor_align.query_align = "";
	donor_align.match_align = "";
	donor_align.target_align = "";

	int search_start, search_end, extra_front, extra_end;
	if (site_start < cur_hsp_start)
	{
		if (site_end >= cur_hsp_start)
		{
			//first, fill in the gap before cur_hsp_start
			//string tmpStr3(cur_hsp_start-site_start, 'N'); //N represent some base pair, but unknown to us, doesn't matter, since it's aligning with '-'
			//MODIFIED: need to obtain the actual chars, since they might be needed for later alignment
			//MODIFIED AGAIN: need to check for in-frame stop codons, if there is, then return false (signalling unqualified donor/acceptor)
			string tmpStr3;

			if (site_start > 0)
				GetSubstrFromVecStrs(chr_seq, true, site_start-1, cur_hsp_start-site_start, tmpStr3);
			else
				GetSubstrFromVecStrs_NegRev(chr_seq, false, -cur_hsp_start, cur_hsp_start-site_start, tmpStr3);

			StrToLower(tmpStr3);

			int notused;

			if (HasInFrameStop(tmpStr3, false, notused, false))
				return false;

			int target_front_extra = (cur_hsp_start - site_start) % 3;
			//if (target_front_extra > 0)
			//	border_target_front = tmpStr3.substr(0, target_front_extra); //record the broken codon

			DNA2AA(tmpStr3, target_front_extra, donor_align.target_align); //convert the rest to AA style

			string tmpStr1((cur_hsp_start-site_start-target_front_extra)/3, '-');
			string tmpStr2((cur_hsp_start-site_start-target_front_extra)/3, ' ');
			donor_align.query_align += tmpStr1;
			donor_align.match_align += tmpStr2;

			//then fill in the rest
			if (site_end <= cur_hsp_end) //...site_start...cur_hsp_start...site_end...cur_hsp_end...
			{
				FindRealPos(targetStr, cur_hsp_start, cur_hsp_end, cur_hsp_start, site_end, 
					search_start, search_end, extra_front, extra_end);

				border_query_aa_end = queryStr[search_end+1];
				
				string qStr = queryStr.substr(0, search_end+1);

				donor_align.query_align += qStr;
				donor_align.match_align += matchStr.substr(0, search_end+1);
				donor_align.target_align += targetStr.substr(0, search_end+1);

				/*if (extra_end > 0) //get the actual couple of extra nucleotides from chromosome sequence (info not available from HSP)
				{
					if (site_start > 0)
						GetSubstrFromVecStrs(chr_seq, true, site_end-extra_end, extra_end, border_target_end);
					else
						GetSubstrFromVecStrs_NegRev(chr_seq, false, -site_end-1, extra_end, border_target_end);
					StrToLower(border_target_end);
				}*/

				align_start = gene_start;
				align_end = align_start + qStr.length() - CountCharInStr(qStr, '-') - 1;
			}
			else //...splice_start...cur_hsp_start...cur_hsp_end...donor...
			{
				donor_align.query_align += queryStr;
				donor_align.match_align += matchStr;
				donor_align.target_align += targetStr;

				//string tmpStr3(site_end-cur_hsp_end, 'N');
				if (site_end > 0)
					GetSubstrFromVecStrs(chr_seq, true, cur_hsp_end, site_end-cur_hsp_end, tmpStr3);
				else
					GetSubstrFromVecStrs_NegRev(chr_seq, false, -site_end-1, site_end-cur_hsp_end, tmpStr3);
				StrToLower(tmpStr3);

				int notused;
				if (HasInFrameStop(tmpStr3, true, notused, false))
					return false;

				int target_end_extra = (site_end - cur_hsp_end) % 3;
				//if (target_end_extra > 0)
				//	border_target_end = tmpStr3.substr(site_end - cur_hsp_end - target_end_extra);

				string tmpStr3_aa;
				DNA2AA(tmpStr3, 0, tmpStr3_aa); //convert the rest to AA style

				donor_align.target_align += tmpStr3_aa;

				string tmpStr1((site_end-cur_hsp_end-target_end_extra)/3, '-');
				string tmpStr2((site_end-cur_hsp_end-target_end_extra)/3, ' ');
				donor_align.query_align += tmpStr1;
				donor_align.match_align += tmpStr2;

				align_start = gene_start;
				align_end = gene_end;

			}

		}
		else //...splice_start...donor...cur_hsp_start...cur_hsp_end... (will this happen? NO)
		{
			cout << "wrong site: splice_start...donor...cur_hsp_start...cur_hsp_end:" 
				<< site_start << "," << site_end << "," << cur_hsp_start << "," << cur_hsp_end << "\n";
			exit(-1);
		}
	}
	else
	{
			if (site_start <= cur_hsp_end)
			{
				if (site_end <= cur_hsp_end) //...cur_hsp_start...splice_start...donor...cur_hsp_end...
				{
					FindRealPos(targetStr, cur_hsp_start, cur_hsp_end, site_start, site_end,
						search_start, search_end, extra_front, extra_end);

					border_query_aa_front = queryStr[search_start-1];
					border_query_aa_end = queryStr[search_end+1];
					/*if (extra_front > 0)
					{
						if (site_start > 0)
							GetSubstrFromVecStrs(chr_seq, true, site_start-1, extra_front, border_target_front);
						else
							GetSubstrFromVecStrs_NegRev(chr_seq, false, -site_start-extra_front, extra_front, border_target_front);
						StrToLower(border_target_front);
					}
					if (extra_end > 0)
					{
						if (site_start > 0)
							GetSubstrFromVecStrs(chr_seq, true, site_end-extra_end, extra_end, border_target_end);
						else
							GetSubstrFromVecStrs_NegRev(chr_seq, false, -site_end-1, extra_end, border_target_end);
						StrToLower(border_target_end);
					}*/

					string qStr;
					if (search_start < targetStr.length())
					{
						qStr = queryStr.substr(search_start, search_end - search_start + 1);
						donor_align.match_align = matchStr.substr(search_start, search_end - search_start + 1);
						donor_align.target_align = targetStr.substr(search_start, search_end - search_start + 1);
					}
					donor_align.query_align = qStr;

					string fStr = queryStr.substr(0, search_start);
					
					align_start = gene_start + fStr.length() - CountCharInStr(fStr, '-');
					align_end = align_start + qStr.length() - CountCharInStr(qStr, '-') - 1;
				}
				else //...cur_hsp_start...splice_start...cur_hsp_end...donor...
				{
					FindRealPos(targetStr, cur_hsp_start, cur_hsp_end, site_start, cur_hsp_end, 
						search_start, search_end, extra_front, extra_end);

					border_query_aa_front = queryStr[search_start-1];
					/*if (extra_front > 0)
					{
						if (site_start > 0)
							GetSubstrFromVecStrs(chr_seq, true, site_start-1, extra_front, border_target_front);
						else
							GetSubstrFromVecStrs_NegRev(chr_seq, false, -site_start-extra_front, extra_front, border_target_front);
						StrToLower(border_target_front);
					}*/

					string qStr;
					if (search_start < targetStr.length())
					{
						qStr = queryStr.substr(search_start, search_end - search_start + 1);
						donor_align.match_align = matchStr.substr(search_start, search_end - search_start + 1);
						donor_align.target_align = targetStr.substr(search_start, search_end - search_start + 1);
					}
					donor_align.query_align = qStr;

					string fStr = queryStr.substr(0, search_start);

					string tmpStr3;
					if (site_end > 0)
						GetSubstrFromVecStrs(chr_seq, true, cur_hsp_end, site_end-cur_hsp_end, tmpStr3);
					else
						GetSubstrFromVecStrs_NegRev(chr_seq, false, -site_end-1, site_end-cur_hsp_end, tmpStr3);
					StrToLower(tmpStr3);
					int notused;
					if (HasInFrameStop(tmpStr3, true, notused, false))
						return false;

					int target_end_extra = (site_end - cur_hsp_end) % 3;
					//if (target_end_extra > 0)
					//	border_target_end = tmpStr3.substr(site_end - cur_hsp_end - target_end_extra);

					string tmpStr3_aa;
					DNA2AA(tmpStr3, 0, tmpStr3_aa);

					donor_align.target_align += tmpStr3_aa;

					string tmpStr1((site_end - cur_hsp_end - target_end_extra)/3, '-');
					string tmpStr2((site_end - cur_hsp_end - target_end_extra)/3, ' ');
					donor_align.query_align += tmpStr1;
					donor_align.match_align += tmpStr2;

					align_start = gene_start + fStr.length() - CountCharInStr(fStr, '-');
					align_end = gene_end;
				}
			}
			else //...cur_hsp_start...cur_hsp_end...splice_start...donor...(impossible?)
			{
				cout << "wrong site: cur_hsp_start...cur_hsp_end...splice_start...donor:" 
					<< cur_hsp_start << "," << cur_hsp_end << "," << site_start << "," << site_end << "\n";
				exit(-1);
			}
	}

//	cout << "getHSPsplicealignment done" << "\n";
	return true;
}

bool ACCP_DONR_graph::ComputeSpliceAlignment_SplitOL1_ForGetHSPSpliceAlignment1(Input_Alignment& donor_align, Input_Alignment& acceptor_align,
	int donor_start, int donor_end, int acceptor_start, int acceptor_end, 
	string& splice_align, int& splice_align_score, 
	string& donor_2nt, string& acceptor_2nt, int donor_tail_num, int acceptor_head_num)
{
#ifdef DEBUG
	dataMgr.outFile << "donor_tail_num:" << donor_tail_num << ";acceptor_head_num:" << acceptor_head_num << "\n";
#endif

	string codonStr;
	if (donor_tail_num > 0) // acceptor_head_num must also > 0
	{
		codonStr = donor_2nt.substr(2-donor_tail_num) + acceptor_2nt.substr(0, acceptor_head_num);
		map<string, char>::iterator codon_tbl_it = DNA_CODON_TBL_SL.find(codonStr);
		if (codon_tbl_it != DNA_CODON_TBL_SL.end())
		{
			if ((*codon_tbl_it).second == '*') //is stop codon
				return false;
			else
				codonStr = (*codon_tbl_it).second; //convert to AA
		}
		else
			codonStr = "X";
	}

	splice_align = "";
	splice_align_score = 0;

	if (donor_start == 0 || acceptor_start == 0) //if either segment is all gaps (impossible?)
	{
		cout <<"impossible?" << donor_start << ";" << acceptor_start << "\n";
		exit(-1);
	}

	int search_start1, search_end1, search_start2, search_end2;
	if (donor_end < acceptor_start) //no overlap
	{
		//revised: now check the trailing / heading gaps in query_align before inserting more gaps!
		//revised again: if there is no overlap and the missing query part is exactly three base pairs,
		//then check donor_prev_nt with acceptor_next_nt against donor_aa_end/acceptor_aa_front, if it's a match, then revise splice_align
		//revised again: now do global alignment for all cases when there's gap on query (globally align query gap with target gap)
		//MODIFIED AGAIN: alignment now starts from the last non-matching position of donor-align to the first non-match of acceptor-align

		int total_gap_len = acceptor_start - donor_end - 1; //this is the total gap length that we should have

		int last_non_gap_pos_donor_align = donor_align.query_align.find_last_not_of("-");

		int first_non_gap_pos_acceptor_align = acceptor_align.query_align.find_first_not_of("-");

		if (last_non_gap_pos_donor_align == string::npos)
			last_non_gap_pos_donor_align = -1;

		if (first_non_gap_pos_acceptor_align == string::npos)
			first_non_gap_pos_acceptor_align = acceptor_align.query_align.length();

		int need_gap_len = total_gap_len - (donor_align.query_align.length() - 1 - last_non_gap_pos_donor_align + first_non_gap_pos_acceptor_align);

		int last_match_pos_donor_align = donor_align.match_align.find_last_not_of("+ ");

		int first_match_pos_acceptor_align = acceptor_align.match_align.find_first_not_of("+ ");

		if (last_match_pos_donor_align == string::npos)
			last_match_pos_donor_align = -1;

		if (first_match_pos_acceptor_align == string::npos)
			first_match_pos_acceptor_align = acceptor_align.query_align.length();

		//if (first_non_gap_pos_acceptor_align > 0 || last_non_gap_pos_donor_align < donor_align.query_align.length()-1)
		//if (first_match_pos_acceptor_align > 0 || last_match_pos_donor_align < donor_align.query_align.length()-1)
		//{
		
		//MODIFIED: do global alignment for the stuff in the gap region
		string target_seq_str = "";

		//target_seq_str must be a multiple of 3 (since donor/acceptor must have matching frame)
		int donor_tail_length = donor_align.target_align.length()-last_match_pos_donor_align-1;
		if (donor_tail_length > 0)
			target_seq_str = donor_align.target_align.substr(last_match_pos_donor_align+1, donor_tail_length);

		if (donor_tail_num > 0) // acceptor_head_num must also > 0
		{
			//string codonStr = donor_2nt.substr(2-donor_tail_num) + acceptor_2nt.substr(0, acceptor_head_num);
			target_seq_str += codonStr;
		}

		int acceptor_head_length = first_match_pos_acceptor_align;
		if (acceptor_head_length > 0)
			target_seq_str += acceptor_align.target_align.substr(0, acceptor_head_length);

		//modified: ignore the very long sequence, just treat it as all gap alignment
		if (target_seq_str.length() > 4096 || target_seq_str.length() == 0)
		{
#ifdef DEBUG
			dataMgr.outFile << "length>4096?" << (target_seq_str.length() > 4096) << ", treat as gaps\n";
#endif
			splice_align += donor_align.match_align;

			if (need_gap_len > 0)
			{
				string tmpStr(need_gap_len, ' ');//(acceptor_start-donor_end-1, ' ');
				splice_align += tmpStr;
			}
			splice_align += acceptor_align.match_align;
			splice_align_score = donor_align.GetAlignmentScore(0, donor_align.query_align.length());
			bool opengap = false, opengap_cur;
			int gap_score = 0;
			for (int i=0; i<need_gap_len; i++) //get the score of all_gap_alignment of length "need_gap_len" 
			{
				gap_score += similarity_score('-', 'x', opengap, opengap_cur);
				opengap = true;
			}
			splice_align_score += gap_score;
			splice_align_score += acceptor_align.GetAlignmentScore(0, acceptor_align.query_align.length());
		}
		else
		{
			//now get the corresponding query segment to be aligned with target_seq_str
			string query_seg_missing_str = query_seq.substr(donor_end, acceptor_start-donor_end-1);

			string query_seq_str = donor_align.query_align.substr(last_match_pos_donor_align+1);
			EraseAll(query_seq_str, '-');

			query_seq_str += query_seg_missing_str;

			string query_acc_str = acceptor_align.query_align.substr(0, first_match_pos_acceptor_align);
			EraseAll(query_acc_str, '-');

			query_seq_str += query_acc_str;
		
			Input_Alignment gapAlign;
			float align_pid;
			int gap_score;
#ifdef PID_BEFORE_SCORE
			gap_score = GetGlobalAlignment_PID_CompScore(query_seq_str, target_seq_str, gapAlign, align_pid, dataMgr.outFile);
#else
			gap_score = GetGlobalAlignment_scorematrix(query_seq_str, target_seq_str, gapAlign, align_pid, dataMgr.outFile);
#endif
			splice_align_score = donor_align.GetAlignmentScore(0, last_match_pos_donor_align+1);
			splice_align_score += gap_score;
			splice_align_score += acceptor_align.GetAlignmentScore(first_match_pos_acceptor_align, acceptor_align.query_align.length());
			splice_align += donor_align.match_align.substr(0, last_match_pos_donor_align+1);
			splice_align += gapAlign.match_align;
			splice_align += acceptor_align.match_align.substr(first_match_pos_acceptor_align);
		}

#ifdef DEBUG
		dataMgr.outFile << "donor and acceptor segments no overlap, should have " << total_gap_len << " gaps; filled in " << need_gap_len << " gaps" << "\n";
			//tmpStr.length() << " gaps" << "\n";
#endif
	}
	else //overlap? ( acceptor_start ... donor_end)
	{
		int identity1, identity2;
		if (donor_start > acceptor_end) //cross align! so only one aligns, the other is turned to gaps 
			//(acceptor_start...acceptor_end...donor_start ... donor_end) (impossible?)
		{
			const string& qStr1 = donor_align.match_align;
			const string& qStr2 = acceptor_align.match_align;

			identity1 = qStr1.length() - CountCharInStr(qStr1, ' ') - CountCharInStr(qStr1, '+');
			identity2 = qStr2.length() - CountCharInStr(qStr2, ' ') - CountCharInStr(qStr2, '+');

			if (identity1 >= identity2)
			{
				splice_align += qStr1;
				string tmpStr(qStr2.length()+(donor_tail_num+acceptor_head_num)/3, ' ');
				splice_align += tmpStr;
				splice_align_score = donor_align.GetAlignmentScore(0,donor_align.query_align.length());
				bool opengap = false, opengap_cur;
				for (int i=0; i<qStr2.length()+(donor_tail_num+acceptor_head_num)/3; i++) //get the score of all_gap_alignment
				{
					splice_align_score += similarity_score('-', 'x', opengap, opengap_cur);
					opengap = true;
				}

				dataMgr.outFile << "acceptor is turned to gap" << "\n";
			}
			else
			{
				string tmpStr(qStr1.length()+(donor_tail_num+acceptor_head_num)/3, ' ');
				splice_align += tmpStr;
				splice_align += qStr2;
				bool opengap = false, opengap_cur;
				int gap_score = 0;
				for (int i=0; i<qStr1.length()+(donor_tail_num+acceptor_head_num)/3; i++) //get the score of all_gap_alignment
				{
					gap_score += similarity_score('-', 'x', opengap, opengap_cur);
					opengap = true;
				}
				splice_align_score = gap_score + acceptor_align.GetAlignmentScore(0, acceptor_align.query_align.length());

				dataMgr.outFile << "donor is turned to gap" << "\n";
			}

#ifdef DEBUG
			dataMgr.outFile << " (cross align, one is turned to gaps)" << "\n";
#endif
		}
		else //(...donor_start...acceptor_start...acceptor_end ... donor_end)
			//or (...donor_start...acceptor_start ... donor_end...acceptor_end)
			//or (...acceptor_start...donor_start...acceptor_end ... donor_end)
			//or (...acceptor_start...donor_start ... donor_end...acceptor_end)
		{
			if (donor_start < acceptor_start)//(...donor_start...acceptor_start...acceptor_end ... donor_end)
			//or (...donor_start...acceptor_start ... donor_end...acceptor_end)
			{
				FindRealPos_DNA(donor_align.query_align, donor_start, donor_end, acceptor_start, donor_end, 
					search_start1, search_end1);

				const string& qStr1 = donor_align.match_align.substr(search_start1, search_end1-search_start1+1);
				identity1 = qStr1.length() - CountCharInStr(qStr1, ' ') - CountCharInStr(qStr1, '+');
				int qStr1_left_len = donor_align.query_align.length() - search_end1 - 1; //length of stuff after search_end1 (gaps)

				if (acceptor_end >= donor_end) //donor_start...acceptor_start...donor_end...acceptor_end
				{
					//the overlapping part (query segment) is [acceptor_start, donor_end]
					FindRealPos_DNA(acceptor_align.query_align, acceptor_start, acceptor_end, acceptor_start, donor_end, 
						search_start2, search_end2);

					const string& qStr2 = acceptor_align.match_align.substr(search_start2, search_end2-search_start2+1);
					int qStr2_left_len = search_start2; //length of stuff before search_start2 (gaps)

					//this is only used here
					const string& donor_ol_queryStr = donor_align.query_align.substr(search_start1, search_end1-search_start1+1);
					const string& acceptor_ol_queryStr = acceptor_align.query_align.substr(search_start2, search_end2-search_start2+1);

					vector<int> donor_ol_queryStr_gap_starts, donor_ol_queryStr_gap_ends;
					vector<int> acceptor_ol_queryStr_gap_starts, acceptor_ol_queryStr_gap_ends;

					GetGaps(donor_ol_queryStr, donor_ol_queryStr_gap_starts, donor_ol_queryStr_gap_ends);
					GetGaps(acceptor_ol_queryStr, acceptor_ol_queryStr_gap_starts, acceptor_ol_queryStr_gap_ends);

#ifdef DEBUG
					dataMgr.outFile << "query:" << donor_ol_queryStr << "\n" << "match:" << qStr1 << "\n";
					dataMgr.outFile << "query:" << acceptor_ol_queryStr << "\n" << "match:" << qStr2 << "\n";

					vector<int>::iterator tmp_it1, tmp_it2=donor_ol_queryStr_gap_ends.begin();
					dataMgr.outFile << "donor gaps: ";
					for (tmp_it1 = donor_ol_queryStr_gap_starts.begin(); tmp_it1 != donor_ol_queryStr_gap_starts.end(); tmp_it1++, tmp_it2++)
						dataMgr.outFile << *tmp_it1 << "-" << *tmp_it2 << "; ";

					dataMgr.outFile << "\nacceptor gaps: ";
					tmp_it2 = acceptor_ol_queryStr_gap_ends.begin();
					for (tmp_it1 = acceptor_ol_queryStr_gap_starts.begin(); tmp_it1 != acceptor_ol_queryStr_gap_starts.end(); tmp_it1++, tmp_it2++)
						dataMgr.outFile << *tmp_it1 << "-" << *tmp_it2 << "; ";
					dataMgr.outFile << "\n";
#endif

					vector<int> donor_ol_matchStr_match_starts, donor_ol_matchStr_match_ends;
					vector<int> acceptor_ol_matchStr_match_starts, acceptor_ol_matchStr_match_ends;
					vector<int> donor_match_end_ids, acceptor_match_start_ids;

					int total_ol_id1 = GetMatches(qStr1, donor_ol_matchStr_match_starts, donor_ol_matchStr_match_ends, donor_match_end_ids);
					int total_ol_id2 = GetMatches(qStr2, acceptor_ol_matchStr_match_starts, acceptor_ol_matchStr_match_ends, acceptor_match_start_ids);

					//now try each match_end on donor_ol_matchStr (qStr1) to be used as split point
					vector<int>::iterator d_start_it=donor_ol_matchStr_match_starts.begin();
					vector<int>::iterator d_end_it = donor_ol_matchStr_match_ends.begin();
					int cur_match_count = 0; //donor side match count
					//int max_match_count=-1, max_match_donor_pos=-1, max_match_acceptor_pos=-1; //initialize
					int max_match_count=total_ol_id2, max_match_donor_pos=-1, max_match_acceptor_pos=-1; //initialize
					//first region_end = -1, region_start = 0 (to simulate all overlapping part uses acceptor)
					for (; d_start_it != donor_ol_matchStr_match_starts.end(); d_start_it++, d_end_it++)
					{
						int region_start = (*d_start_it);
						int region_end = (*d_end_it);
						cur_match_count += region_end - region_start + 1;

						int donor_queryPos = ConvertToNoGappedPos(region_end, donor_ol_queryStr_gap_starts, donor_ol_queryStr_gap_ends, false);
						int acceptor_matchPos = ConvertToGappedPos(donor_queryPos, acceptor_ol_queryStr_gap_starts, acceptor_ol_queryStr_gap_ends);
						int acceptor_match_count = ComputeBackId(acceptor_matchPos, acceptor_ol_matchStr_match_starts, 
							acceptor_ol_matchStr_match_ends, acceptor_match_start_ids, total_ol_id2);

						int total_match_count = cur_match_count + acceptor_match_count;
#ifdef DEBUG
						dataMgr.outFile << "donor matchPos (abs): " << region_end << "; converted to queryPos (all match): " << donor_queryPos;
						dataMgr.outFile << "; acceptor matchPos (abs): " << acceptor_matchPos << "\n";
						dataMgr.outFile << "cur max: " << max_match_count << "; this time: " << total_match_count << "\n";
#endif
						if (total_match_count > max_match_count)
						{
							max_match_count = total_match_count;
							max_match_donor_pos = region_end;
							max_match_acceptor_pos = acceptor_matchPos;
#ifdef DEBUG
							dataMgr.outFile << "cur donor pos: " << max_match_donor_pos  << "; cur acc pos: " << max_match_acceptor_pos << "\n";
#endif
						}
					}

					//now also compute from the other side (acceptor side), if both come to the same max_match_count but at different positions, keep matched stuff at both sides
					cur_match_count = 0; //acceptor side match count
					int max_match_count_acc=total_ol_id1, max_match_donor_pos_acc=donor_ol_queryStr.length(), max_match_acceptor_pos_acc=acceptor_ol_queryStr.length(); //initialize
					//first region_end = -1, region_start = 0 (to simulate all overlapping part uses acceptor)
					vector<int>::reverse_iterator a_start_it = acceptor_ol_matchStr_match_starts.rbegin(); 
					vector<int>::reverse_iterator a_end_it = acceptor_ol_matchStr_match_ends.rbegin();
					for (; a_start_it != acceptor_ol_matchStr_match_starts.rend(); a_start_it++, a_end_it++)
					{
						int region_start = (*a_start_it);
						int region_end = (*a_end_it);
						cur_match_count += region_end - region_start + 1;

						int acceptor_queryPos = ConvertToNoGappedPos(region_start, acceptor_ol_queryStr_gap_starts, acceptor_ol_queryStr_gap_ends, false);
						int donor_matchPos = ConvertToGappedPos(acceptor_queryPos, donor_ol_queryStr_gap_starts, donor_ol_queryStr_gap_ends);

						int donor_match_count = ComputeFrontId(donor_matchPos, donor_ol_matchStr_match_starts, 
							donor_ol_matchStr_match_ends, donor_match_end_ids, total_ol_id1);

						int total_match_count = cur_match_count + donor_match_count;
#ifdef DEBUG
						dataMgr.outFile << "acceptor matchPos (abs): " << region_start << "; converted to queryPos (all match): " << acceptor_queryPos;
						dataMgr.outFile << "; donor matchPos (abs): " << donor_matchPos << "\n";
						dataMgr.outFile << "cur max: " << max_match_count_acc << "; this time: " << total_match_count << "\n";
#endif

						if (total_match_count > max_match_count_acc)
						{
							max_match_count_acc = total_match_count;
							max_match_acceptor_pos_acc = region_start;
							max_match_donor_pos_acc = donor_matchPos;
#ifdef DEBUG
							dataMgr.outFile << "cur donor pos: " << max_match_donor_pos_acc  << "; cur acc pos: " << max_match_acceptor_pos_acc << "\n";
#endif
						}
					}

					int donor_gap_len = qStr1.length() - max_match_donor_pos - 1 + qStr1_left_len + (donor_tail_num+acceptor_head_num)/3;
					int acceptor_gap_len = qStr2_left_len + max_match_acceptor_pos+1;
					int acceptor_gap_len_init = acceptor_gap_len; //hold this position, used later for splice_align

#ifdef DEBUG
					dataMgr.outFile << "First: donor_gap_len: " << donor_gap_len << "; acceptor_gap_len: " << acceptor_gap_len << "\n";
#endif

					string tmpStr1;
					int strIt;
					bool opengap_cur = false;
					if (donor_gap_len == 0)
					{
						splice_align += donor_align.match_align.substr(0, search_start1+max_match_donor_pos+1);
						splice_align_score += donor_align.GetAlignmentScore(0, search_start1+max_match_donor_pos+1);
						/*if ((donor_tail_num+acceptor_head_num)/3 > 0)
						{
							splice_align += ' ';							
							splice_align_score += similarity_score('-', 'x', false, opengap_cur);
						}*/
					}
					else if (donor_gap_len > 0)
					{
						tmpStr1 = donor_align.target_align.substr(search_start1+max_match_donor_pos+1);
						tmpStr1 += codonStr;

						string tmpStr2 = donor_align.target_align.substr(0, search_start1+max_match_donor_pos+1);
						int don_valid_pos = tmpStr2.find_last_not_of('-');

						string tmpStr3;
						if (don_valid_pos == string::npos)
						{
							tmpStr3 = donor_align.query_align.substr(0, search_start1+max_match_donor_pos+1);
							don_valid_pos = tmpStr3.length();
						}
						else
						{
							if (don_valid_pos > 0)
								tmpStr3 = donor_align.query_align.substr(don_valid_pos+1, search_start1+max_match_donor_pos - don_valid_pos);
							splice_align += donor_align.match_align.substr(0, don_valid_pos+1);
							splice_align_score += donor_align.GetAlignmentScore(0, don_valid_pos+1);
						}

						if (donor_gap_len != tmpStr1.length()) //is this thing right? ... let's check
						{
							cout << "donor_gap_len is not tmpStr1.length()" << "\n";
							dataMgr.outFile << donor_align << "\n";
							dataMgr.outFile << "tmpStr1 from " << search_start1+max_match_donor_pos+1 << "\n";
							dataMgr.outFile << tmpStr1 << "\n";
							dataMgr.outFile << "donor_gap_len:" << donor_gap_len << " vs tmpStr.length:" << tmpStr1.length() << "\n";
							dataMgr.outFile << "qStr1:" << qStr1 << "\n";
							dataMgr.outFile << "search_start1:" << search_start1 << "; search_end1:" << search_end1 
								<< "; max_match_donor_pos:" << max_match_donor_pos << "\n";
							dataMgr.outFile << "qStr1_left_len:" << qStr1_left_len << "\n";
							exit(-1);
						}

						while ((strIt = tmpStr1.find('-')) != string::npos)
							tmpStr1.erase(strIt, 1);

						//donor_gap_len = tmpStr1.length();
						if (tmpStr3.length()>0 && tmpStr1.length()>0)//(don_valid_pos > 0 && 
						{
							Input_Alignment newAlign;
							float align_pid;
#ifdef PID_BEFORE_SCORE
							splice_align_score += GetGlobalAlignment_PID_CompScore(tmpStr1, tmpStr3, newAlign, align_pid, dataMgr.outFile);
#else
							splice_align_score += GetGlobalAlignment_scorematrix(tmpStr1, tmpStr3, newAlign, align_pid, dataMgr.outFile);
#endif
							splice_align += newAlign.match_align;
						}
						else //tmpStr3 is empty
						{
							donor_gap_len = tmpStr1.length();
							string tmpStr(donor_gap_len, ' ');
							splice_align += tmpStr;
							for (int i=0; i<donor_gap_len; i++)
								splice_align_score += similarity_score('-', 'x', opengap_cur, opengap_cur);
						}
					}

					if (acceptor_gap_len == 0)
					{
						splice_align += acceptor_align.match_align;
						splice_align_score += acceptor_align.GetAlignmentScore(0, acceptor_align.query_align.length());
					}
					else if (acceptor_gap_len > 0)
					{
						tmpStr1 = acceptor_align.target_align.substr(0, acceptor_gap_len);

						string tmpStr2 = acceptor_align.target_align.substr(acceptor_gap_len);
						int acc_valid_pos = tmpStr2.find_first_not_of('-');

						string tmpStr3;
						if (acc_valid_pos == string::npos) //everything after acceptor_gap_len is "-"!
						{
							tmpStr3 = acceptor_align.query_align.substr(acceptor_gap_len);
							acc_valid_pos = tmpStr3.length();
						}
						else
						{
							if (acc_valid_pos > 0)
								tmpStr3 = acceptor_align.query_align.substr(acceptor_gap_len, acc_valid_pos);
						}

						if (acceptor_gap_len != tmpStr1.length()) //is this thing right? ... let's check
						{
							cout << "acceptor_gap_len is not tmpStr1.length()" << "\n";
							exit(-1);
						}
						while ((strIt = tmpStr1.find('-')) != string::npos)
							tmpStr1.erase(strIt);

						//acceptor_gap_len = tmpStr1.length();
						if (acc_valid_pos > 0)//align tmpStr1 with tmpStr3, to see if the gaps can be further reduced
						{
							Input_Alignment newAlign;
							float align_pid;
#ifdef PID_BEFORE_SCORE
							splice_align_score += GetGlobalAlignment_PID_CompScore(tmpStr1, tmpStr3, newAlign, align_pid, dataMgr.outFile);
#else
							splice_align_score += GetGlobalAlignment_scorematrix(tmpStr1, tmpStr3, newAlign, align_pid, dataMgr.outFile);
#endif
							splice_align += newAlign.match_align;
						}
						else //acc_valid_pos must be 0
						{
							acceptor_gap_len = tmpStr1.length();
							string tmpStr(acceptor_gap_len, ' ');
							splice_align += tmpStr;
							for (int i=0; i<acceptor_gap_len; i++)
								splice_align_score += similarity_score('-', 'x', opengap_cur, opengap_cur);
						}

						if (acceptor_gap_len_init+acc_valid_pos < acceptor_align.match_align.length())
						{
							splice_align += acceptor_align.match_align.substr(acceptor_gap_len_init+acc_valid_pos);
							splice_align_score += acceptor_align.GetAlignmentScore(acceptor_gap_len_init+acc_valid_pos, acceptor_align.query_align.length());
						}
					}
#ifdef DEBUG
					dataMgr.outFile << "Now: donor_gap_len: " << donor_gap_len << "; acceptor_gap_len: " << acceptor_gap_len << "\n";
#endif
				}
				else //donor_start...acceptor_start...acceptor_end...donor_end
				{
					const string& qStr2 = acceptor_align.match_align;
					identity2 = qStr2.length() - CountCharInStr(qStr2, ' ') - CountCharInStr(qStr2, '+');
					
					if (identity1 > identity2)
					{
						splice_align += donor_align.match_align;
						splice_align_score += donor_align.GetAlignmentScore(0, donor_align.query_align.length());
						int target_bp = acceptor_align.target_align.length() - CountCharInStr(acceptor_align.target_align, '-')
							+(donor_tail_num+acceptor_head_num)/3;
						string tmpStr(target_bp, ' ');
						splice_align += tmpStr;
						bool opengap = false;
						for (int i=0; i<target_bp; i++)
							splice_align_score += similarity_score('-', 'x', opengap, opengap);
					}
					else
					{
						splice_align += donor_align.match_align.substr(0, search_start1);
						splice_align_score += donor_align.GetAlignmentScore(0, search_start1);

						string donor_tmpStr = donor_align.target_align.substr(search_start1);
						int donor_gap_len = donor_align.target_align.length() - search_start1 - CountCharInStr(donor_tmpStr, '-')
							+(donor_tail_num+acceptor_head_num)/3;
						string tmpStr(donor_gap_len, ' ');//(qStr1.length(), ' ');
						splice_align += tmpStr;
						bool opengap = false;
						for (int i=0; i<donor_gap_len; i++)
							splice_align_score += similarity_score('-', 'x', opengap, opengap);

						splice_align += qStr2;
						splice_align_score += acceptor_align.GetAlignmentScore(0, acceptor_align.query_align.length());
					}

				}

			}
			else // (donor_start >= acceptor_start), cross align?
				//(...acceptor_start...donor_start...acceptor_end ... donor_end)
			//or (...acceptor_start...donor_start ... donor_end...acceptor_end)
			{
				const string& qStr1 = donor_align.match_align;
				identity1 = qStr1.length() - CountCharInStr(qStr1, ' ') - CountCharInStr(qStr1, '+');

				int cur_end = acceptor_end < donor_end? acceptor_end : donor_end;
				FindRealPos_DNA(acceptor_align.query_align, acceptor_start, acceptor_end, acceptor_start, cur_end, 
					search_start2, search_end2);

				const string& qStr2 = acceptor_align.match_align.substr(0, search_end2+1);
				identity2 = qStr2.length() - CountCharInStr(qStr2, ' ') - CountCharInStr(qStr2, '+');

				if (identity1 > identity2)
				{
					splice_align += qStr1;
					splice_align_score += donor_align.GetAlignmentScore(0, donor_align.query_align.length());

					string acceptor_tmpStr = acceptor_align.target_align.substr(0, search_end2+1);
					int gap_len = search_end2+1 - CountCharInStr(acceptor_tmpStr, '-')+(donor_tail_num+acceptor_head_num)/3;
					string tmpStr(gap_len, ' ');
					splice_align += tmpStr;
					bool opengap = false;
					for (int i=0; i<gap_len; i++)
						splice_align_score += similarity_score('-', 'x', opengap, opengap);

					if (acceptor_end > cur_end)
					{
						splice_align += acceptor_align.match_align.substr(search_end2+1);
						splice_align_score += acceptor_align.GetAlignmentScore(search_end2+1, acceptor_align.query_align.length());
					}
				}
				else
				{
					int gap_len = donor_align.target_align.length() - CountCharInStr(donor_align.target_align, '-')+(donor_tail_num+acceptor_head_num)/3;
					string tmpStr(gap_len, ' ');
					splice_align += tmpStr;
					bool opengap = false;
					for (int i=0; i<gap_len; i++)
						splice_align_score += similarity_score('-', 'x', opengap, opengap);

					splice_align += acceptor_align.match_align;
					splice_align_score += acceptor_align.GetAlignmentScore(0, acceptor_align.query_align.length());
				}
			}
		}
	}

	return true;

}
