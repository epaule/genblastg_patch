#include "edge.h"

//update current node's gene_pid (map structure)
void HSP_node::UpdateMap_For_SkipDist(int dest_start, int dest_end, float dest_pid)
{
	map<int, float>::iterator startIt = gene_pid.lower_bound(dest_start);

	if (startIt == gene_pid.end()) //everything is smaller than dest_start
	{
		gene_pid.insert(map<int,float>::value_type(dest_start, dest_pid));
		gene_pid.insert(map<int, float>::value_type(dest_end, 0));
	}

	else
	{
		if ((*startIt).first > dest_start)
		{
			float prev_pid = 0;
			if (startIt != gene_pid.begin())
			{
				startIt--;
				prev_pid = (*startIt).second;
			}

			gene_pid.insert(map<int,float>::value_type(dest_start, prev_pid));
		}

		//now process dest_end
		map<int, float>::iterator endIt = gene_pid.lower_bound(dest_end);
		if (endIt == gene_pid.end())
			gene_pid.insert(map<int, float>::value_type(dest_end, 0));
		else
		{
			if ((*endIt).first > dest_end) //now endIt can not point to the first element of the map (we have processed dest_start already)
			{
				endIt--;
				gene_pid.insert(map<int, float>::value_type(dest_end, (*endIt).second));
			}
		}

		//now update every pid between dest_start and dest_end
		map<int, float>::iterator mapIt = gene_pid.find(dest_start);
		for (; mapIt != gene_pid.find(dest_end); mapIt++)
		{
			if ((*mapIt).second < dest_pid)
				(*mapIt).second= dest_pid;
		}

	}


}


float HSP_node::CompSkipDist(int prev_start, int cur_start)
{
	if (prev_start >= cur_start)
		return 0;

	float result = 0;
	
	map<int, float>::iterator it = gene_pid.find(prev_start);
	int prev_bound = prev_start;
	float prev_pid = (*it).second;

	it++;
	for (; it != gene_pid.find(cur_start); it++)
	{
		result += Score(prev_bound, (*it).first, prev_pid);
		prev_bound = (*it).first;
		prev_pid = (*it).second;
	}
	result += Score(prev_bound, (*it).first, prev_pid);

	return result;
}


AD_Edge_Penalty AD_Edge_Penalty::operator+(const AD_Edge_Penalty& p)
{
		AD_Edge_Penalty tmp(exact_match, gap, pos_match, query_match_pos); //first, copy the first set into "tmp"
		//tmp.exact_match = exact_match + p.exact_match;
		//tmp.gap = gap + p.gap;
		//tmp.pos_match = pos_match + p.pos_match;

		set< pair<int,int> >::const_iterator it;
		set< pair<int,int> >::iterator tmp_it;// hold_it;
		set< pair<int,int> >::reverse_iterator rev_it;
		//vector<int>::const_iterator it;
		//vector<int>::iterator tmp_it;
		//vector<int>::reverse_iterator rev_it, tmp_rev_it;

		//also copy p.exact_match, gap
		int exact_match_p = p.exact_match;
		int gap_p = p.gap;

		//now use "tmp" for further operation
		int count, count_transfer;
		it = p.query_match_pos.begin();
		while (it!= p.query_match_pos.end() && !tmp.query_match_pos.empty())
		{
			count = 0;
			//hold_it = it;
			//hold_it++;

			int curseg_start = (*it).first;
			int curseg_end = (*it).second;
			//int curseg_start = *it;
			//it++;
			//int curseg_end = *it;
			int curseg_len = curseg_end - curseg_start;
			rev_it = tmp.query_match_pos.rbegin(); //last element
			tmp_it = tmp.query_match_pos.end();
			tmp_it--; //used only to match up rev_it position, for erase()!

			if (curseg_start > (*rev_it).second)
			//if (curseg_start > *rev_it)
				break;
			else
			{
				if (curseg_start == (*rev_it).second)
				//if ( curseg_start == *rev_it)
				{
					const_cast< pair<int,int>& >(*rev_it).second = curseg_end; //merge
					//*rev_it = curseg_end;
					//p.query_match_pos.erase( it );
					it++;
					break;
				}
				else
				{					
					count += curseg_len;
					bool pos_found = false;
					while (rev_it != tmp.query_match_pos.rend() && count>0)
					{
						//tmp_rev_it = rev_it;
						//rev_it++;
						//tmp_it--;

						if ((*rev_it).first <= curseg_start)
						//if (*rev_it <= curseg_start)
						{
							pos_found = true;
							if ((*rev_it).second > curseg_start)
								count -= (*rev_it).second - curseg_start;
							//if (*tmp_rev_it > curseg_start)
							//	count -= *tmp_rev_it - curseg_start;								
							break;
						}
						else
							count -= (*rev_it).second - (*rev_it).first;
							//count -= *tmp_rev_it - *rev_it;
						rev_it++; 
						tmp_it--; 
					}
					if (count > 0)
					{
						count_transfer = (curseg_len - count)*3;
						tmp.exact_match -= count_transfer;
						tmp.gap += count_transfer;
						if (pos_found)
						{
							if ((*rev_it).second >= curseg_start) //merge
							//if (*tmp_rev_it >= curseg_start)
							{
								//exact_match += curseg_len;
								const_cast< pair<int,int>& >(*rev_it).second = curseg_end;
								//*tmp_rev_it = curseg_end;
								//p.query_match_pos.erase( it );
								//p.exact_match -= curseg_len;
							}
							tmp_it++;
						}
						else //must be rev_it(tmp_it) is beyond the beginning of container (rend())
							tmp_it = tmp.query_match_pos.begin();
						//tmp_it++;
						if (tmp_it != tmp.query_match_pos.end())
							tmp.query_match_pos.erase(tmp_it, tmp.query_match_pos.end());
					}
					else
					{
						count_transfer = curseg_len*3;
						exact_match_p -= count_transfer;
						gap_p += count_transfer;
						//p.query_match_pos.erase( it );
					}
				}
			}

			//it = hold_it;
			it++;
		}
		//copy(query_match_pos.begin(), query_match_pos.end(), tmp.query_match_pos.begin());
		//copy(p.query_match_pos.begin(), p.query_match_pos.end(), tmp.query_match_pos.end());
		//for (tmp_it = query_match_pos.begin(); tmp_it != query_match_pos.end(); tmp_it++)
		//	tmp.query_match_pos.insert( *tmp_it );
		//for (tmp_it = p.query_match_pos.begin(); tmp_it != p.query_match_pos.end(); tmp_it++)
		//	tmp.query_match_pos.insert( *tmp_it );

		while (it != p.query_match_pos.end())
		{
			tmp.query_match_pos.insert( *it );
			it++;
		}
		tmp.exact_match += exact_match_p;
		tmp.gap += gap_p;
		tmp.pos_match += p.pos_match; //unchanged
		//tmp.exact_match = exact_match + p.exact_match;
		//tmp.gap = gap + p.gap;
		//tmp.pos_match = pos_match + p.pos_match;
		return tmp;
}


AD_Edge_QueryCover AD_Edge_QueryCover::operator+(const AD_Edge_QueryCover& qc)
{
	AD_Edge_QueryCover tmp=*this;

/*		//set< pair<int,int> >::const_iterator it;
		//set< pair<int,int> >::iterator tmp_it, tmp_end_it;// hold_it;
		//set< pair<int,int> >::reverse_iterator rev_it;
*/
	set<QuerySeg_Start_End_PID>::const_iterator it;
	set<QuerySeg_Start_End_PID>::iterator tmp_it, tmp_end_it;
	set<QuerySeg_Start_End_PID>::reverse_iterator rev_it;

		//examine each query segment covered by "qc" (the later operand) against each segment in "tmp" (backward)
		//if they overlap, resolve by merging them properly
		it = qc.query_segments.begin();
		while (it!= qc.query_segments.end() && !tmp.query_segments.empty())
		{
			int curseg_start = (*it).first;
			int curseg_end = (*it).second;
			rev_it = tmp.query_segments.rbegin(); //last element
			tmp_it = tmp.query_segments.end();
			tmp_it--; //used only to match up rev_it position, for erase()! (start iterator for erase)

			if (curseg_start > (*rev_it).second)
				break;
			else
			{
/*				if (curseg_start == (*rev_it).second)
				{
					const_cast< pair<int,int>& >(*rev_it).second = curseg_end; //merge into "tmp"
					tmp.query_cover += curseg_end - curseg_start;
					it++; //skip this segment in "qc" since it's already contained in "tmp"
					break;
				}
				else
				{					
*/					bool pos_found = false;
					while (rev_it != tmp.query_segments.rend() )//&& !pos_found)
					{
						tmp.query_cover -= (*rev_it).second - (*rev_it).first + 1;
						if ((*rev_it).first <= curseg_start) //find the "tmp.segment"  that is just before "qc.segment.start"
						{
							pos_found = true;
							break;
						}						
						rev_it++;
						tmp_it--;
					}
					if (!pos_found) //must be rev_it(tmp_it) is beyond the beginning of container (rend())
						tmp_it = tmp.query_segments.begin();
					//map<int, vector<int> > segment_map;
					map<int, vector<pair<float, int> > > segment_map;
					for (tmp_end_it = tmp_it; tmp_end_it != tmp.query_segments.end(); tmp_end_it++)
					{
						//MyMapInsert((*tmp_end_it).first, (*tmp_end_it).second, true, segment_map);
						//MyMapInsert((*tmp_end_it).second, 0, false, segment_map);
						MyMapInsert((*tmp_end_it).first, (*tmp_end_it).second, (*tmp_end_it).pid, true, segment_map);
						MyMapInsert((*tmp_end_it).second, 0, (*tmp_end_it).pid, false, segment_map);
					}
					//MyMapInsert((*it).first, (*it).second, true, segment_map);
					//MyMapInsert((*it).second, 0, false, segment_map);
					MyMapInsert((*it).first, (*it).second, (*it).pid, true, segment_map);
					MyMapInsert((*it).second, 0, (*it).pid, false, segment_map);

					if (tmp_it != tmp.query_segments.end()) //now we can erase first
						tmp.query_segments.erase(tmp_it, tmp.query_segments.end());

					//map<int, vector<int> >::iterator map_it, map_next_it;
					map<int, vector<pair<float, int> > >::iterator map_it, map_next_it;
					int start=-1, end=-1;
					for (map_it=segment_map.begin(); map_it != segment_map.end(); map_it++)
					{
						map_next_it=map_it;
						map_next_it++;
						if (!(*map_it).second.empty()) //non-empty range, means it's covered, there must be something next
						{
							//vector<int>::iterator vec_it;
							vector<pair<float, int> >::iterator vec_it;
							float cur_pid = 0.0;
							for (vec_it = (*map_it).second.begin(); vec_it != (*map_it).second.end(); vec_it++)
							{
								if (map_next_it != segment_map.end())
								//if ((*vec_it) > (*map_next_it).first)
								if ((*vec_it).second > (*map_next_it).first)
									(*map_next_it).second.push_back(*vec_it);
								
								if ((*vec_it).first > cur_pid)
									cur_pid = (*vec_it).first; //cur_pid is the largest pid for current segment
							}

							//if (start == -1)
							//	start = (*map_it).first;
							if (map_next_it != segment_map.end())
							if ((*map_next_it).second.empty()) //after we try extend the (*it) into (*next_it), if it's still empty
							{
								//end = (*map_next_it).first;
								//tmp.query_segments.insert(pair<int, int>(start, end));
								//tmp.query_cover += end-start+1;
								//start = -1; //reset start
								tmp.query_segments.insert(QuerySeg_Start_End_PID((*map_it).first, (*map_next_it).first, cur_pid));
								tmp.query_cover += (*map_next_it).first - (*map_it).first + 1;
							}
							else
							{
								tmp.query_segments.insert(QuerySeg_Start_End_PID((*map_it).first, (*map_next_it).first-1, cur_pid));
								tmp.query_cover += (*map_next_it).first - (*map_it).first ;
							}
							else //no next and this is already the last section, this section must be "start == end"
							{
								tmp.query_segments.insert(QuerySeg_Start_End_PID((*map_it).first, (*map_it).first, cur_pid));
								tmp.query_cover++;
							}
						}
						
					}

				//}
			}

			it++;
		}

		while (it != qc.query_segments.end())
		{
			tmp.query_segments.insert( *it );
			tmp.query_cover += (*it).second - (*it).first + 1;
			it++;
		}

		//merge adjacent segments (modified: that also has same pid)
		tmp.qc_score = 0;
		tmp_it = tmp.query_segments.begin();
		while (tmp_it != tmp.query_segments.end())
		{
			tmp.qc_score += (*tmp_it).pid/100 * ((*tmp_it).second - (*tmp_it).first + 1); //compute score
			tmp_end_it = tmp_it;
			tmp_end_it++;
			if (tmp_end_it != tmp.query_segments.end())
				if ((*tmp_it).second == (*tmp_end_it).first - 1 && (*tmp_it).pid == (*tmp_end_it).pid)
				{
					const_cast< QuerySeg_Start_End_PID& >(*tmp_it).second = (*tmp_end_it).second;
					tmp.query_segments.erase(tmp_end_it);
					continue; //do not forward tmp_it, since *tmp_it has just been modified
				}
			tmp_it++;
		}


		//========================================
		//process target_segments
		set< pair<int,int> >::const_iterator target_it;
		set< pair<int,int> >::iterator target_tmp_it, target_tmp_end_it;// hold_it;
		set< pair<int,int> >::reverse_iterator target_rev_it;

		//compute gap_len
		if (tmp.target_segments.empty())
		{
			tmp.gap_len += qc.gap_len;
			tmp.gap_at_end = qc.gap_at_end;
		}
		else
		{
			if (!qc.target_segments.empty())
			{
				target_tmp_it = tmp.target_segments.end();
				target_tmp_it--;
				target_it = qc.target_segments.begin();
				if ((*target_tmp_it).second == (*target_it).first && tmp.gap_at_end)
					tmp.gap_len += qc.gap_len - 1;
				else
					tmp.gap_len += qc.gap_len;

				tmp.gap_at_end = qc.gap_at_end;
			}
		}

		//examine each query segment covered by "qc" (the later operand) against each segment in "tmp" (backward)
		//if they overlap, resolve by merging them properly
		target_it = qc.target_segments.begin();
		while (target_it!= qc.target_segments.end() && !tmp.target_segments.empty())
		{
			int curseg_start = (*target_it).first;
			int curseg_end = (*target_it).second;
			target_rev_it = tmp.target_segments.rbegin(); //last element
			target_tmp_it = tmp.target_segments.end();
			target_tmp_it--; //used only to match up rev_it position, for erase()! (start iterator for erase)

			if (curseg_start > (*target_rev_it).second)
				break;
			else
			{
/*				if (curseg_start == (*rev_it).second)
				{
					const_cast< pair<int,int>& >(*rev_it).second = curseg_end; //merge into "tmp"
					tmp.query_cover += curseg_end - curseg_start;
					it++; //skip this segment in "qc" since it's already contained in "tmp"
					break;
				}
				else
				{					
*/					bool pos_found = false;
					while (target_rev_it != tmp.target_segments.rend() )
					{
						tmp.target_len -= (*target_rev_it).second - (*target_rev_it).first + 1;
						if ((*target_rev_it).first <= curseg_start) //find the "tmp.segment"  that is just before "qc.segment.start"
						{
							pos_found = true;
							break; //do not forward iterators
						}
						target_rev_it++;
						target_tmp_it--;
					}
					if (!pos_found) //must be rev_it(tmp_it) is beyond the beginning of container (rend())
						target_tmp_it = tmp.target_segments.begin();
					map<int, vector<int> > segment_map;
					for (target_tmp_end_it = target_tmp_it; target_tmp_end_it != tmp.target_segments.end(); target_tmp_end_it++)
					{
						MyMapInsert((*target_tmp_end_it).first, (*target_tmp_end_it).second, true, segment_map);
						MyMapInsert((*target_tmp_end_it).second, 0, false, segment_map);
					}
					MyMapInsert((*target_it).first, (*target_it).second, true, segment_map);
					MyMapInsert((*target_it).second, 0, false, segment_map);

					if (target_tmp_it != tmp.target_segments.end()) //now we can erase first
						tmp.target_segments.erase(target_tmp_it, tmp.target_segments.end());

					map<int, vector<int> >::iterator map_it, map_next_it;
					int start=-1, end=-1;
					for (map_it=segment_map.begin(); map_it != segment_map.end(); map_it++)
					{
						map_next_it=map_it;
						map_next_it++;
						if (!(*map_it).second.empty()) //non-empty range, means it's covered, there must be something next
						{
							vector<int>::iterator vec_it;
							for (vec_it = (*map_it).second.begin(); vec_it != (*map_it).second.end(); vec_it++)
							{
								if (map_next_it != segment_map.end())
								if ((*vec_it) > (*map_next_it).first)
									(*map_next_it).second.push_back(*vec_it);
							}

							if (start == -1)
								start = (*map_it).first;
							if (map_next_it != segment_map.end())
							{
							if ((*map_next_it).second.empty()) //after we try extend the (*it) into (*next_it), if it's still empty
							{
								end = (*map_next_it).first;
								tmp.target_segments.insert(pair<int, int>(start, end));
								tmp.target_len += end-start+1;
								start = -1; //reset start
							}
							}
							else //no next and this is already the last section, this section must be "start == end"
							{
								tmp.target_segments.insert(pair<int, int>((*map_it).first, (*map_it).first));
								tmp.target_len++;
							}
						}
						
					}

				//}
			}

			target_it++;
		}

		while (target_it != qc.target_segments.end())
		{
			tmp.target_segments.insert( *target_it );
			tmp.target_len += (*target_it).second - (*target_it).first + 1;
			target_it++;
		}

		//merge adjacent segments 
		target_tmp_it = tmp.target_segments.begin();
		while (target_tmp_it != tmp.target_segments.end())
		{
			target_tmp_end_it = target_tmp_it;
			target_tmp_end_it++;
			if (target_tmp_end_it != tmp.target_segments.end())
				if ((*target_tmp_it).second == (*target_tmp_end_it).first - 1 )
				{
					const_cast< pair<int, int>& >(*target_tmp_it).second = (*target_tmp_end_it).second;
					tmp.target_segments.erase(target_tmp_end_it);
					continue; //do not forward tmp_it, since *tmp_it has just been modified
				}
			target_tmp_it++;
		}


		return tmp;

}
