#include "data_manager.h"



bool HSP_Gene_Pair::FitWithNeighborHSPs(HSP_Gene_Pair* neighborHSP, bool isPosStrand) //neighborHSP is the reference

{

	int cur_hStart = HSP_start;

	int cur_hEnd = HSP_end;



	int nb_hStart = neighborHSP->HSP_start;

	int nb_hEnd = neighborHSP->HSP_end;

	//check compatibility of genomic positions

	if (!((cur_hStart < nb_hStart && cur_hEnd < nb_hEnd) || (cur_hStart > nb_hStart && cur_hEnd > nb_hEnd)))

		return false;



	//check compatibility of query correspondence

	if (isPosStrand)

		return ((cur_hStart < nb_hStart && gene_start < neighborHSP->gene_start && gene_end < neighborHSP->gene_end) || 

			(cur_hStart > nb_hStart && gene_start > neighborHSP->gene_start && gene_end > neighborHSP->gene_end));

	else

		return ((cur_hStart < nb_hStart && gene_start > neighborHSP->gene_start && gene_end > neighborHSP->gene_end) || 

			(cur_hStart > nb_hStart && gene_start < neighborHSP->gene_start && gene_end < neighborHSP->gene_end));



}



bool HSP_Gene_Pair::HSPInRegion(int region_start, int region_end)

{

	if (HSP_start <= region_start)

	{

		if (HSP_end >= region_end)

		{

			return true;

		}

		else

		{

			if (HSP_end >= region_start)

				return true;

			else

				return false;

		}

	}

	else

	{

		if (HSP_start > region_end)

			return false;

		else

			return true;

	}

}



bool HSP_Gene_Pair::Overlap(HSP_Gene_Pair* prev_HSP,  int& gOL_start, int& gOL_end,

							 int& hOL_start, int& hOL_end, bool& g_align, bool& h_align)

{

	bool OL_gene = false;

	bool OL_hsp = false;



	int cur_gStart = gene_start;

	int cur_gEnd = gene_end;

	int cur_hStart = HSP_start;

	int cur_hEnd = HSP_end;



	int prev_gEnd = prev_HSP->gene_end;

	int prev_hStart = prev_HSP->HSP_start;

	int prev_hEnd = prev_HSP->HSP_end;

		

	if (cur_gStart <= prev_gEnd)

	{

		gOL_start = cur_gStart;

		if (prev_gEnd < cur_gEnd)

			gOL_end = prev_gEnd;

		else

			gOL_end = cur_gEnd;

		OL_gene = true;

	}

	else

		if (cur_gStart == prev_gEnd + 1)

			g_align = true;



	if (cur_hStart >= prev_hStart) //positive strand HSP

	{

		if (cur_hStart <= prev_hEnd)

		{

			hOL_start = cur_hStart;

			if (prev_hEnd < cur_hEnd)

				hOL_end = prev_hEnd;

			else

				hOL_end = cur_hEnd;

			OL_hsp = true;

		}

		else

			if (cur_hStart == prev_hEnd + 1)

				h_align = true;

	}

	else //negatve strand HSP

	{

		if (prev_hStart <= cur_hEnd)

		{

			hOL_start = prev_hStart;

			if (cur_hEnd < prev_hEnd)

				hOL_end = cur_hEnd;

			else

				hOL_end = prev_hEnd;

			OL_hsp = true;

		}

		else

			if (prev_hStart == cur_hEnd + 1)

				h_align = true;

	}



	return OL_gene || OL_hsp;

}





//read input file, store HSPs

bool DataManager::ReadFile()

{

	int count=0;

	int neg_count = 0;

	int id, start, end, HSP_start, HSP_end;

	float pid;



	//ifstream inFile(inputFile.c_str());

	if (!inputFile_Open)

	{

		inputFile_is.open(inputFile.c_str());

		if (!inputFile_is.is_open())

		{

			cout << "input file open error" << "\n";

			return false;

		}

		inputFile_Open = true;

	}

	else

		query_gene = next_query_gene;

	

	//bool cur_gene_finish = false;

	string line;

	vector<string> items;

	vector<string> headers;

	//int headerIndex[9];

	while (!inputFile_is.eof())	//read in file, parse, 

	{

		items.clear();

		headers.clear();

		//get one set of "ID, start, end, pid, HSP_start, HSP_end"

		getline(inputFile_is, line);

		if (line.find("HSP_ID") != string::npos) //header line

		{

			int pos = 0;

			int start_pos = 0;

			while (( pos = line.find(DELIMIT, start_pos) )!= string::npos) //the last one is ignored???

			{

				while (pos == start_pos) //extra delimiter???!!!

				{

					start_pos = pos+1;

					pos = line.find(DELIMIT, start_pos);

				}

				headers.push_back(line.substr(start_pos, pos-start_pos));

				start_pos = pos+1;

			}



			for (int i=0; i<headers.size(); i++)

			{

				//cout << "header " << i << ":" << headers[i] << ";"<< "\n";

				if (headers[i].compare("HSP_ID") == 0)

					headerIndex[0] = i;

				else if (headers[i].compare("query_start")==0)

					headerIndex[1] = i;

				else if (headers[i].compare("query_end")==0)

					headerIndex[2] = i;

				else if (headers[i].compare("Chr_start")==0)

					headerIndex[3] = i;

				else if (headers[i].compare("Chr_end")==0)

					headerIndex[4] = i;

				else if (headers[i].compare("P_Idn")==0)

					headerIndex[5] = i;

				else if (headers[i].compare("chr_strand")==0)

					headerIndex[6] = i;

				else if (headers[i].compare("Chr")==0)

					headerIndex[7] = i;

				else if (headers[i].compare("Protein_length")==0)

					headerIndex[8] = i;

			}

			

		}

		else

		if (line[0] == '>') //gene name line

		{

			if (!cur_gene_start)

			{

				cur_gene_start = true;

				query_gene = line.substr(1, line.length()-1);

				//cout << "line: " << line << "\nquery_gene: " << query_gene << "\n";

			}

			else

			{

				//cur_gene_finish = true;

				next_query_gene = line.substr(1, line.length()-1);

				break;

			}

		}			

		else //now it must be a data line

		if (line.find(DELIMIT) != string::npos) //make sure there is '\t'

		{

			int pos = 0;

			int start_pos = 0;

			while (( pos = line.find(DELIMIT, start_pos) )!= string::npos) //the last one is ignored???

			{

				while (pos == start_pos) //extra delimiter???!!!

				{

					start_pos = pos+1;

					pos = line.find(DELIMIT, start_pos);

				}

				items.push_back(line.substr(start_pos, pos-start_pos));

				start_pos = pos+1;

			}



			id = atoi(items[headerIndex[0]].c_str());

			start = atoi(items[headerIndex[1]].c_str());

			end = atoi(items[headerIndex[2]].c_str());

			HSP_start = atoi(items[headerIndex[3]].c_str());

			HSP_end = atoi(items[headerIndex[4]].c_str());

			pid = atof(items[headerIndex[5]].c_str());



			if (!query_len) //if still zero

				query_len = atoi(items[headerIndex[8]].c_str()); //3 times protein length?



			int chr_total = HSP_chr.size();

			int c;

			for (c = 0; c < chr_total; c++)

			{

				if (HSP_chr[c].compare(items[headerIndex[7]]) == 0)

					break;

			}

			if (c == chr_total) //not found

			{

				HSP_chr.push_back(items[headerIndex[7]]);

				HSP_gene.push_back(vector<HSP_Gene_Pair>());

				HSP_neg_gene.push_back(vector<HSP_Gene_Pair>());

			}



			int strand = atoi(items[headerIndex[6]].c_str());



			switch (strand)

			{

			case 1:

				{

					//store everything first

					HSP_gene[c].push_back(HSP_Gene_Pair(id,  start,end,HSP_start,HSP_end,pid));



					count++;

//					cout << "read " << count << "HSP" << "\n";

				}

				break;

			case -1:

				{

					HSP_neg_gene[c].push_back(HSP_Gene_Pair(id, start,end,HSP_start,HSP_end,pid));



					neg_count++;

//					cout << "read " << neg_count << "neg HSP" << "\n";

				}

				break;

			default:

				cout << "something is wrong" << "\n";

				break;

			}

			

			

		}

	}



	gene_start_end_map.clear();

	fragment_score_map.clear();

	int i,j;

	for (j=0; j<HSP_gene.size(); j++)

		for (i=0; i<HSP_gene[j].size(); i++)

			ReadOneGeneSegment(HSP_gene[j][i].gene_start, HSP_gene[j][i].gene_end, HSP_gene[j][i].pid);

	for (j=0; j<HSP_neg_gene.size(); j++)

		for (i=0; i<HSP_neg_gene[j].size(); i++)

			ReadOneGeneSegment(HSP_neg_gene[j][i].gene_start, HSP_neg_gene[j][i].gene_end, HSP_neg_gene[j][i].pid);



	GetFragmentScores();



	if (inputFile_is.eof())

	{

		inputFile_is.close();

		inputFile_Open = false; //reset for next file

		cur_gene_start = false;

		inputFile_Finish = true;

	}



	return true;

}



//read input file, store HSPs

bool DataManager::ReadFile_Skip()

{

	int count=0;

	int neg_count = 0;

//	int id, start, end, HSP_start, HSP_end;

//	float pid;



	//ifstream inFile(inputFile.c_str());

	if (!inputFile_Open)

	{

		inputFile_is.open(inputFile.c_str());

		if (!inputFile_is.is_open())

		{

			cout << "input file open error" << "\n";

			return false;

		}

		inputFile_Open = true;

	}

	else

		query_gene = next_query_gene;

	

	//bool cur_gene_finish = false;

	string line;

	vector<string> items;

	vector<string> headers;

	//int headerIndex[9];

	while (!inputFile_is.eof())	//read in file, parse, 

	{

		items.clear();

		headers.clear();

		//get one set of "ID, start, end, pid, HSP_start, HSP_end"

		getline(inputFile_is, line);

		if (line.find("HSP_ID") != string::npos) //header line

		{

			int pos = 0;

			int start_pos = 0;

			while (( pos = line.find(DELIMIT, start_pos) )!= string::npos) //the last one is ignored???

			{

				while (pos == start_pos) //extra delimiter???!!!

				{

					start_pos = pos+1;

					pos = line.find(DELIMIT, start_pos);

				}

				headers.push_back(line.substr(start_pos, pos-start_pos));

				start_pos = pos+1;

			}



			for (int i=0; i<headers.size(); i++)

			{

				if (headers[i].compare("HSP_ID") == 0)

					headerIndex[0] = i;

				else if (headers[i].compare("query_start")==0)

					headerIndex[1] = i;

				else if (headers[i].compare("query_end")==0)

					headerIndex[2] = i;

				else if (headers[i].compare("Chr_start")==0)

					headerIndex[3] = i;

				else if (headers[i].compare("Chr_end")==0)

					headerIndex[4] = i;

				else if (headers[i].compare("P_Idn")==0)

					headerIndex[5] = i;

				else if (headers[i].compare("chr_strand")==0)

					headerIndex[6] = i;

				else if (headers[i].compare("Chr")==0)

					headerIndex[7] = i;

				else if (headers[i].compare("Protein_length")==0)

					headerIndex[8] = i;

			}

			

		}

		else

		if (line[0] == '>') //gene name line

		{

			if (!cur_gene_start)

			{

				cur_gene_start = true;

				query_gene = line.substr(1, line.length()-1);

			}

			else

			{

				//cur_gene_finish = true;

				next_query_gene = line.substr(1, line.length()-1);

				break;

			}

		}			

/*		else //now it must be a data line

		if (line.find(DELIMIT) != string::npos) //make sure there is '\t'

		{

			int pos = 0;

			int start_pos = 0;

			while (( pos = line.find(DELIMIT, start_pos) )!= string::npos) //the last one is ignored???

			{

				items.push_back(line.substr(start_pos, pos-start_pos));

				start_pos = pos+1;

			}



			id = atoi(items[headerIndex[0]].c_str());

			start = atoi(items[headerIndex[1]].c_str());

			end = atoi(items[headerIndex[2]].c_str());

			HSP_start = atoi(items[headerIndex[3]].c_str());

			HSP_end = atoi(items[headerIndex[4]].c_str());

			pid = atof(items[headerIndex[5]].c_str());



			if (!query_len) //if still zero

				query_len = atoi(items[headerIndex[8]].c_str()); //3 times protein length?



			int chr_total = HSP_chr.size();

			int c;

			for (c = 0; c < chr_total; c++)

			{

				if (HSP_chr[c].compare(items[headerIndex[7]]) == 0)

					break;

			}

			if (c == chr_total) //not found

			{

				HSP_chr.push_back(items[headerIndex[7]]);

				HSP_gene.push_back(vector<HSP_Gene_Pair>());

				HSP_neg_gene.push_back(vector<HSP_Gene_Pair>());

			}



			int strand = atoi(items[headerIndex[6]].c_str());



			switch (strand)

			{

			case 1:

				{

					//store everything first

					HSP_gene[c].push_back(HSP_Gene_Pair(id,  start,end,HSP_start,HSP_end,pid));



					count++;

//					cout << "read " << count << "HSP" << "\n";

				}

				break;

			case -1:

				{

					HSP_neg_gene[c].push_back(HSP_Gene_Pair(id, start,end,HSP_start,HSP_end,pid));



					neg_count++;

//					cout << "read " << neg_count << "neg HSP" << "\n";

				}

				break;

			default:

				cout << "something is wrong" << "\n";

				break;

			}

			

			

		}

*/

	}



	if (inputFile_is.eof())

	{

		inputFile_is.close();

		inputFile_Open = false; //reset for next file

		cur_gene_start = false;

		inputFile_Finish = true;

	}



	return true;

}



/*

void DataManager::PrepareData(bool isPosStrand, int chr_index)

{

	

	if (isPosStrand)

		PrepareMaps(HSP_gene[chr_index], isPosStrand);

	else

		PrepareMaps(HSP_neg_gene[chr_index], isPosStrand);



}



void DataManager::PrepareMaps(vector<HSP_Gene_Pair>& curHSPs, bool isPosStrand)

{

	gene_start_end_map.clear();

	for (int i=0; i<curHSPs.size(); i++)

		ReadOneGeneSegment(curHSPs[i].gene_start, curHSPs[i].gene_end, i);



	fragment_score_map.clear();

	GetFragmentScores(curHSPs);



	gene_start_HSP_num_map.clear();

	ProcHSPs(curHSPs, isPosStrand);



}

*/





void DataManager::PrepareData(bool isPosStrand, int chr_index)

{

	

	gene_start_HSP_num_map.clear();

	if (isPosStrand)

		ProcHSPs(HSP_gene[chr_index], isPosStrand);

	else

		ProcHSPs(HSP_neg_gene[chr_index], isPosStrand);



}



/*

//fill up gene_start_end_map

void DataManager::ReadOneGeneSegment(int start, int end, int ind)

{



	MyMapInsert( start, ind, true, gene_start_end_map);

	MyMapInsert( end, ind, false, gene_start_end_map);



}

*/



void DataManager::ReadOneGeneSegment(int start, int end, float PID)

{

	MyMapInsert( start, end, PID, true, gene_start_end_map);

	MyMapInsert( end, end, PID, false, gene_start_end_map);



}



/*

void DataManager::MyMapInsert( const int pos, const int count, const bool isStart)

{

	map<int, vector<int> >::iterator sit = gene_start_end_map.lower_bound(pos);

	if (sit != gene_start_end_map.end() && !gene_start_end_map.key_comp()(pos, (*sit).first)) //already has this key

	{

		if (isStart)

			(*sit).second.push_back(count); //record which lines start from this position

	}

	else

	{

		if (isStart)

			gene_start_end_map.insert(sit, map<int, vector<int> >::value_type(pos, vector<int>(1, count)));//construct the vector with 1 integer: count

		else

			gene_start_end_map.insert(sit, map<int, vector<int> >::value_type(pos, vector<int>()));

	}

}



*/





/*

//using gene_start_end_map, compute fragment scores, fill up fragment_score_map

void DataManager::GetFragmentScores(vector<HSP_Gene_Pair>& curHSPs)

{

	if (gene_start_end_map.empty())

		return;



	int start, end;

	float pid;



	

	map<int, vector<int> >::iterator it=gene_start_end_map.begin();



	//add special segment so it starts from 1 (no use: it's 0!)

//	if ((*it).first > 1)

//		fragment_score_map.insert(map<int, float>::value_type(1, Score(1, (*it).first, 0)));



	map<int, vector<int> >::iterator end_it=gene_start_end_map.begin();

	end_it++;





	while (end_it!=gene_start_end_map.end())

	{

		start = (*it).first; //fragment start

		end = (*end_it).first; //fragment end



		pid =0;

		int count = (*it).second.size();

		vector<int>::iterator vecIt = (*it).second.begin(); //computing average pid

		while (vecIt != (*it).second.end() )

		{

			pid += curHSPs[(*vecIt)].pid;

			if (curHSPs[(*vecIt)].gene_end > end)

				(*end_it).second.push_back(*vecIt);



			vecIt++;

		}

		if (count)

			pid = pid / (float)count;



		//store results

		fragment_score_map.insert(map<int, float>::value_type(start, Score(start, end, pid)));

		

		it++;

		end_it++;





	}





}

*/



//using gene_start_end_map, compute fragment scores, fill up fragment_score_map

void DataManager::GetFragmentScores()

{

	if (gene_start_end_map.empty())

		return;



	int start, end;

	float pid;



	

	map<int, vector<pair<float, int> > >::iterator it=gene_start_end_map.begin();



	map<int, vector<pair<float, int> > >::iterator end_it=gene_start_end_map.begin();

	end_it++;





	while (end_it!=gene_start_end_map.end())

	{

		start = (*it).first; //fragment start

		end = (*end_it).first; //fragment end



		pid =0;

		int count = (*it).second.size();

		vector<pair<float, int> >::iterator vecIt = (*it).second.begin(); //computing average pid

		while (vecIt != (*it).second.end() )

		{

			pid += (*vecIt).first;

			if ((*vecIt).second > end)

				(*end_it).second.push_back(*vecIt);



			vecIt++;

		}

		if (count)

			pid = pid / (float)count;



		//store results

		//fragment_score_map.insert(map<int, float>::value_type(start, Score(start, end, pid))); 

		fragment_score_map.insert(map<int,float>::value_type(start, Score(start, start, pid))); //[start,start]

		if (start+1 < end)

			fragment_score_map.insert(map<int, float>::value_type(start+1, Score(start+1, end-1, pid))); //(start, end)

		fragment_score_map.insert(map<int, float>::value_type(end, Score(end, end, pid))); //[end,end]

		

		it++;

		end_it++;





	}



#ifdef DEBUG_VERSION

	//cout << "fragment_score_map size: " << fragment_score_map.size() << "\n";

	PrintFragmentScores(outFile);

#endif



}



void DataManager::PrintFragmentScores(ostream& os)

{

	map<int, float>::iterator it = fragment_score_map.begin();

	while (it != fragment_score_map.end())

	{

		os << "query " << (*it).first << "->: score: " << (*it).second << "\n";

		it++;

	}

}





//some preprocessing stuff, for graph construction

void DataManager::ProcHSPs(vector<HSP_Gene_Pair>& curHSPs, bool isPosStrand)

{



	//fill up gene_start_HSP_num_map, to be used later for constructing skip edges

	//(skip edge: the node with smallest gene_start > curNode->gene_end && HSP_start > curNode->HSP_end)

	if (isPosStrand)

	{

		//adding 2 special HSPs

		curHSPs.push_back(HSP_Gene_Pair(-1, -1,0,-1,0,0));

		curHSPs.push_back(HSP_Gene_Pair(-2, INT_MAX-1,INT_MAX,INT_MAX-1,INT_MAX,0));





		//fragment scores have been computed, now we can sort by HSP_start

		sort(curHSPs.begin(), curHSPs.end());

		//***HSP_num from now on is fixed***

		

		int total = curHSPs.size();



		for (int i=1; i<total; i++)

			gene_start_HSP_num_map.insert(multimap<int,int>::value_type(curHSPs[i].gene_start, i));

	}

	else

	{

		//adding 2 special HSPs

		curHSPs.push_back(HSP_Gene_Pair(-1, INT_MAX-1,INT_MAX,-1,0,0));

		curHSPs.push_back(HSP_Gene_Pair(-2, -1,0,INT_MAX-1,INT_MAX,0));



		//fragment scores have been computed, now we can sort by HSP_start

		sort(curHSPs.begin(), curHSPs.end());

		//***HSP_num from now on is fixed***



		int total = curHSPs.size();



		for (int i=1; i<total; i++)

			gene_start_HSP_num_map.insert(multimap<int,int>::value_type(curHSPs[i].gene_end, i)); //record gene_end for negative strand

	}



}



//the output function

void DataManager::PrintGroups(bool printOverview)

{

	if (printOverview)

		scale = ComputeScale();



	outFile << "//*****************START******************//" << "\n";

	outFile << "//   candidate genes: (in ranked order)   //" << "\n"; //44 chars on one line

	outFile << "//for query: " << query_gene;

	for (int i=0; i<44-13-2-query_gene.length(); i++)

		outFile << " ";

	outFile << "//" << "\n";

	outFile << "//****************************************//" << "\n" ;


	if (groups.size() == 0)

	{

		outFile << "NONE" << "\n";

	}

	else

	{


	int rank = 1;
	int count = 1;
	multimap<Group_Info, vector<HSP_Gene_Pair*> >::iterator mapIt = groups.begin();

	float dist = (*mapIt).first.score;

	for (; mapIt != groups.end(); mapIt++)

	{
		if ((*mapIt).first.score < dist)

		{

			rank++;

			dist = (*mapIt).first.score;

		}



		//print header

		outFile << "\n";


		cur_chr_name = HSP_chr[(*mapIt).first.chr_index];

		outFile << "\n" << query_gene << "|" << cur_chr_name << ":" << (*mapIt).first.HSP_start << ".." << (*mapIt).first.HSP_end;
		
		if (OUTPUT_GFF && phase1_only)
		{
		gff_os << "##sequence-region\t" << cur_chr_name << "\t" 
			<< (*mapIt).first.HSP_start << "\t" << (*mapIt).first.HSP_end << "\n";
			//<< LenOfStrVec((*(chr_name_seq.find(cur_chr_name))).second) << "\n";
		gff_os << cur_chr_name << "\tgenBlastA" << USER_ID << "\tregion\t" << (*mapIt).first.HSP_start << "\t" << (*mapIt).first.HSP_end
			<< "\t" << (*mapIt).first.score << "\t" ;
		}
		
		if ((*mapIt).first.isPosStrand)
		{
			outFile << "|+|gene cover:";
			if (OUTPUT_GFF && phase1_only)
				gff_os << "+\t.\tID=";
		}
		else
		{
			outFile << "|-|gene cover:";
			if (OUTPUT_GFF && phase1_only)
				gff_os << "-\t.\tID=";
		}
		if (OUTPUT_GFF && phase1_only)
			gff_os << query_gene << "-R" << rank << "-" << count << "\n";

		outFile << (*mapIt).first.gene_cover_len << "(" << (float)(*mapIt).first.gene_cover_len/(float)query_len*100 

			<< "%)|score:" << (*mapIt).first.score << "|rank:" << rank << "\n";


		if (printOverview)

			PrintGroupsOverview(mapIt);

		else

			PrintGroupsTxt(mapIt, rank, count);



		count++;



	}



	}



	outFile << "\n" << "//******************END*******************//" << "\n" << "\n";





}



void DataManager::PrintGroupsTxt(multimap<Group_Info, vector<HSP_Gene_Pair*> >::iterator mapIt, int rank, int count)

{

		int cur_hsp_len;

		vector<HSP_Gene_Pair*>::iterator vecIt;

		for (vecIt = (*mapIt).second.begin(); vecIt != (*mapIt).second.end(); vecIt++)

		{

			//print text

			outFile << *(*vecIt);

			if (OUTPUT_GFF && phase1_only)
			{
			gff_os << cur_chr_name << "\tgenBlastA" << USER_ID << "\tmatch\t" << (*vecIt)->HSP_start << "\t" << (*vecIt)->HSP_end 
				<< "\t" << (*vecIt)->pid << "\t" ;
			if ((*mapIt).first.isPosStrand)
				gff_os << "+\t.\tName=HSP";
			else
				gff_os << "-\t.\tName=HSP";
			gff_os << (*vecIt)->ID << ";Parent="<< query_gene << "-R" << rank << "-" << count << "\n";
			}

			//update chr_shortest_hsp_len (to be used for PrintGroupsOverview)

			cur_hsp_len = (*vecIt)->HSP_end - (*vecIt)->HSP_start + 1;

			if (chr_shortest_hsp_len > cur_hsp_len)

				chr_shortest_hsp_len = cur_hsp_len;



		}



		//update chr_longest_cand_len (to be used for PrintGroupsOverview)

		int cur_start = (*mapIt).first.HSP_start;

		int cur_end = (*mapIt).first.HSP_end;

		int cur_total_len = cur_end - cur_start + 1;

		if (chr_longest_cand_len < cur_total_len)

			chr_longest_cand_len = cur_total_len;



}



void DataManager::PrintGroupsOverview(multimap<Group_Info, vector<HSP_Gene_Pair*> >::iterator mapIt)

{	

	vector<string> strChrs(3, ""); //three strings: 1-entire ">>>"; 2-hsp_id; 3-hsp_marker

	

	//print entire ">>>" first

	int start = (*mapIt).first.HSP_start;

	int end = (*mapIt).first.HSP_end;

	int total_len = end - start + 1;



	int i;

	for (i=0; i<(total_len + scale -1)/scale; i++)

		//outFile << ">";

		strChrs[0] += '>';

	

	//collect HSPs in reverse order (so always starts with HSPs with smallest coordinate)

	int cur_pos = start; //keep track of multiples of scales

	int next_hsp_start = start; //actual HSP_start position

	vector<HSP_Gene_Pair*>::reverse_iterator vecIt;

	char tmp[32];

	string tmpStr;

	int cur_id_pos = start;



	for (vecIt = (*mapIt).second.rbegin(); vecIt != (*mapIt).second.rend(); vecIt++)

	{

		if ((*vecIt)->HSP_start < next_hsp_start) //overlap

		{

			//outFile << '|';

			//strChrs[2] += '|';

			strChrs[1] += '|';

			cur_pos += scale;

		}

		else

		{

			//outFile << '>';

			if ((*vecIt)->HSP_start > cur_pos)

			{

				int num_gap1 = ((*vecIt)->HSP_start-cur_pos + scale - 1)/scale;

				for (i=0; i<num_gap1; i++)

					//outFile << ' ';

					//strChrs[2] += ' ';

					strChrs[1] += ' ';

				cur_pos += num_gap1*scale;

			}



			if ((*vecIt)->HSP_start > cur_id_pos)

			{

				int num_gap2 = ((*vecIt)->HSP_start - cur_id_pos + scale -1)/scale;

				for (i=0; i<num_gap2; i++)

					//strChrs[1] += ' ';

					strChrs[2] += ' ';

				cur_id_pos += num_gap2*scale;

			}

		}

		

		int len = (*vecIt)->HSP_end - cur_pos + 1;

		

		int num_frag = (len + scale - 1)/scale; //must be at least 2 fragments?



		//print HSP_ID

		//_itoa((*vecIt)->ID, tmp, 10);

		sprintf(tmp, "%d", (*vecIt)->ID);		

		//strChrs[1] += '[';

		//strChrs[1] += tmp;

		//strChrs[1] += ']';

		strChrs[2] += '[';

		strChrs[2] += tmp;

		strChrs[2] += ']';

		tmpStr = tmp;

		int size = tmpStr.size()+2;

		cur_id_pos += size*scale;



		if (!(*mapIt).first.isPosStrand)

		{

			//outFile << '<';

			//strChrs[2] += '<';

			strChrs[1] += '<';

		}



		if (num_frag > 0)

		{

			for (i=0; i<num_frag-1; i++)

				//outFile << '=';

				//strChrs[2] += '=';

				strChrs[1] += '=';

			cur_pos += (num_frag-1)*scale;

		}



		if ((*mapIt).first.isPosStrand)

			//outFile << '>';

			//strChrs[2] += '>';

			strChrs[1] += '>';

		cur_pos += scale;

		

		if (cur_id_pos < cur_pos)

		{

			for (i=0; i<(cur_pos - cur_id_pos)/scale; i++)

				//strChrs[1] += ' ';

				strChrs[2] += ' ';

			cur_id_pos = cur_pos;

		}



		next_hsp_start = (*vecIt)->HSP_end+1;

	}





	OutputStrsWithWrap(strChrs);

}



int DataManager::ComputeScale()

{

	int frag_scale = (chr_shortest_hsp_len+1) / 2;



	if ((chr_longest_cand_len+frag_scale-1) / frag_scale < MAX_CHAR_PER_LINE)

		frag_scale = (chr_longest_cand_len+MAX_CHAR_PER_LINE-1)/MAX_CHAR_PER_LINE;



	return frag_scale;



}



void DataManager::OutputStrsWithWrap(vector<string>& strs)

{

	int i;

	int size_max = 0;

	for (i=0; i<strs.size(); i++)

		if (strs[i].size() > size_max)

			size_max = strs[i].size();



	int j=0;

	while (j<size_max)

	{

		for (i=0; i<strs.size(); i++)

		{

			if (j<strs[i].size())

				outFile << strs[i].substr(j, MAX_CHAR_PER_LINE);

			outFile << "\n";

		}

		outFile << "\n" << "\n";

		j += MAX_CHAR_PER_LINE;

	}



}





void DataManager::OutputStrsWithWrapHeader(vector<string>& strs, vector<string>& headers)

{

	int i;

	int size_max = 0;

	for (i=0; i<strs.size(); i++)

		if (strs[i].size() > size_max)

			size_max = strs[i].size();



	int j=0;

	while (j<size_max)

	{

		for (i=0; i<strs.size(); i++)

		{

			if (j<strs[i].size())

			{

				outFile << headers[i] << "\t";

				outFile << strs[i].substr(j, MAX_CHAR_PER_LINE);

			}

			outFile << "\n";

		}

		outFile << "\n" << "\n";

		j += MAX_CHAR_PER_LINE;

	}



}

void DataManager::CompRegion(multimap<Group_Info, vector<HSP_Gene_Pair*> >::iterator groupIt, 
							 int& region_start, int& region_end)
{
	int len = 2 * EXPLORE_END_EXON_LEN;

	int hsp_chr_start, hsp_gene_start, hsp_chr_end;
	int cur_pos;
	if ((*groupIt).first.isPosStrand)
	{
		hsp_chr_start = (*groupIt).first.HSP_start;
		hsp_gene_start = (*groupIt).second.back()->gene_start;
		hsp_chr_end = (*groupIt).first.HSP_end;

		cur_pos = hsp_chr_start-3*(hsp_gene_start-1); //the supposed mark "1" position
		if (cur_pos < 1)
			cur_pos = 1;
		int len1 = (int)((hsp_chr_end - cur_pos + 1) * EXPLORE_END_EXON_PER);
		if (len < len1)
			len = len1;
		cur_pos -= len;
		if (cur_pos < 1)
			cur_pos = 1;
		region_start = cur_pos;
		region_end = hsp_chr_end + len; //region_end may exceed the chromosome length, will adjust when reading chromosomes
	}
	else
	{
		hsp_gene_start = (*groupIt).second.front()->gene_start;
		hsp_chr_start = (*groupIt).second.back()->HSP_start;
		hsp_chr_end = (*groupIt).second.front()->HSP_end;

		cur_pos = hsp_chr_end+3*(hsp_gene_start-1)-2;
		int len1 = (int)((cur_pos - hsp_chr_start + 1) * EXPLORE_END_EXON_PER);
		if (len < len1)
			len = len1;
		region_end = cur_pos + len;
		hsp_chr_start -= len;
		if (hsp_chr_start < 1)
			hsp_chr_start = 1;
		region_start = hsp_chr_start;
	}

}


void DataManager::PrepareOutputHSPs()
{
	//reset for each query gene
	//chr_groups.clear();
	//chromosomes.clear();
	group_start_seqno.clear();
	group_dna_regions.clear();

	cout << "current gene: " << query_gene << "\n";



	//collect all HSP_IDs (the original ID in input file) that are included in the final output, sorted

	set<int>	HSP_IDs; 

	//collect chromosome names that are in the final output
	//set<string> chr_names;

	//collect <chr,(start,end)> so later we know which genesplicer output to read

//	map<string, pair<int, int> > chr_coordinates;

	multimap<Group_Info, vector<HSP_Gene_Pair*> >::iterator groupMapIt = groups.begin();

	//collect rank and count here, to be used in GenePredict6(), where we simply extract groups with certain rank and count
	//int rank = 1;
	//float dist = (*groupMapIt).first.score;
	int count=1; //for GFF file only, to distinguish genes of the same rank

	map<string, vector<Group_Count_RegStart_RegEnd> > chr_group_region;
	int region_start, region_end;
	for (; groupMapIt != groups.end(); groupMapIt++)
	{
		//outFile << "group:" << "\n";
		/*if ((*groupMapIt).first.score < dist)
		{
			rank++;
			dist = (*groupMapIt).first.score;
		}*/

		string cur_chr = HSP_chr[(*groupMapIt).first.chr_index];
		CompRegion(groupMapIt, region_start, region_end);
		/*map<string, vector<Group_Rank_Count> >::iterator chr_it = chr_groups.find(cur_chr);
		if (chr_it == chr_groups.end())
		{
			vector<Group_Rank_Count> cur_vec;
			cur_vec.push_back(Group_Rank_Count(rank, count));
			chr_groups.insert(map<string, vector<Group_Rank_Count> >::value_type(cur_chr, cur_vec) );
			chromosomes.push_back(cur_chr);
		}
		else
		{
			(*chr_it).second.push_back(Group_Rank_Count(rank, count));
		}*/
		map<string, vector<Group_Count_RegStart_RegEnd> >::iterator chr_it = chr_group_region.find(cur_chr);
		if (chr_it == chr_group_region.end())
		{
			vector<Group_Count_RegStart_RegEnd> cur_vec;
			cur_vec.push_back(Group_Count_RegStart_RegEnd(count, region_start, region_end));
			chr_group_region.insert(map<string, vector<Group_Count_RegStart_RegEnd> >::value_type(cur_chr, cur_vec));
		}
		else
		{
			(*chr_it).second.push_back(Group_Count_RegStart_RegEnd(count, region_start, region_end));
		}

		vector<HSP_Gene_Pair*>::iterator hspVecIt = (*groupMapIt).second.begin();
		for (; hspVecIt != (*groupMapIt).second.end(); hspVecIt++)
		{
			//outFile << "hsp_" << (*hspVecIt)->ID << "\n";
			HSP_IDs.insert((*hspVecIt)->ID);
		}

		count++;
	}



	//outFile << "HSP_IDs size:" << HSP_IDs.size() << "\n";

	GetAlignments(HSP_IDs);


	//collect all acceptors and donors from genesplicer output
//	GetSpliceSites(chr_coordinates);

	//UPDATED: now examine the collected regions, merge overlapping ones, organize them and store into chr_group_region_count
	map<string, map<int, Group_RegEnd_Count> > chr_group_region_count;
	map<string, vector<Group_Count_RegStart_RegEnd> >::iterator group_region_it;
	for (group_region_it = chr_group_region.begin(); group_region_it != chr_group_region.end(); group_region_it++)
	{
		string cur_chr = (*group_region_it).first;
		map<int, Group_RegEnd_Count> cur_map; //empty one
	
		vector<Group_Count_RegStart_RegEnd>& regions = (*group_region_it).second;
		sort(regions.begin(), regions.end());
		int i=0, j=1;
		while (j<=regions.size())
		{
			if (j < regions.size() && regions[i].region_end >= regions[j].region_start) //region i and j overlap
			{
				if (regions[j].region_end > regions[i].region_end)
					regions[i].region_end = regions[j].region_end;
				UpdateGroupRegEnd(cur_map, regions[i].region_start, regions[i].region_end, regions[j].count);//record j
				regions.erase(regions.begin()+j);
			}
			else
			{
				UpdateGroupRegEnd(cur_map, regions[i].region_start, regions[i].region_end, regions[i].count);//record i
				i++;
				j=i+1;
			}
		}
		chr_group_region_count.insert(map<string, map<int, Group_RegEnd_Count> >::value_type(cur_chr, cur_map));
	}

	//now use chr_group_region_count to collect necessary DNA regions from chromosome sequence file
	GetChromosomes(chr_group_region_count);
}

void DataManager::UpdateGroupRegEnd(map<int, Group_RegEnd_Count>& cur_map, int region_start, int region_end, int group_count)
{
	map<int, Group_RegEnd_Count>::iterator cur_map_it = cur_map.find(region_start);
	if (cur_map_it == cur_map.end())
	{
		vector<int> cur_vec;
		cur_vec.push_back(group_count); //record j
		Group_RegEnd_Count reg(region_end, cur_vec);
		cur_map.insert(map<int, Group_RegEnd_Count>::value_type(region_start, reg));
	}
	else
	{
		int prev_end = (*cur_map_it).second.region_end;
		if (prev_end < region_end)
			const_cast<int&>((*cur_map_it).second.region_end) = region_end;
		const_cast<vector<int>&>((*cur_map_it).second.group_count).push_back(group_count);
	}
}



//not used anymore (genesplicer: out!)

void DataManager::GetSpliceSites(map<string, pair<int, int> >& chr_coordinates)

{

	cout << "get splice sites" << "\n";



	map<string, pair<int, int> >::iterator mapIt;



	for (mapIt = chr_coordinates.begin(); mapIt != chr_coordinates.end(); mapIt++)

	{

		int start = (*mapIt).second.first;

		int end = (*mapIt).second.second;



		//genesplicer output filename format: "spliceFile"_"chr_name", where "spliceFile" is 

		string filename = spliceFile;

		filename += "_";

		filename += (*mapIt).first;



		ifstream inSpliceFile(filename.c_str());

		if (!inSpliceFile.is_open())

		{

			cout << "splice site file open error" << "\n";

			return;

		}



		//first set up a new entry (empty)

		Splice_Sites curSites;



		string line;

		while (!inSpliceFile.eof())

		{

			getline(inSpliceFile, line);

			if (line.empty())

				continue;



			int pos1= line.find(' ');

			if (pos1 == string::npos)

				continue;



			int first=atoi(line.substr(0, pos1).c_str());

			if (first < start) //skip those not between start and end

				continue;



			int pos2 = line.find(' ', pos1+1);

			int second=atoi(line.substr(pos1+1, pos2-pos1-1).c_str());

			if (second > end) //skip those not between start and end

				continue;



			int pos3 = line.find(' ', pos2+1);

			int pos4 = line.find(' ', pos3+1);



			string name = line.substr(pos4+1); //acceptor or donor



			//update entry by adding new sites to it

			if (name.compare("acceptor") == 0)

				if (first < second) //positive acceptor

					curSites.LoadSite(second, true);

				else //negative acceptor

					curSites.LoadSite(-second, true); //use negative coordinates for negative strand

			else

				if (first < second) //positive donor

					curSites.LoadSite(first, false);

				else //negative donor

					curSites.LoadSite(-first, false); //use negative coordinates for negative strand

		}



		//(pos and neg all in the same list, sort them together)

		//positive strand order: from small to big (splice result already ordered)

		//negative strand order: absolute value from big to small (because of negative integers)

		sort(curSites.acceptors.begin(), curSites.acceptors.end()); 

		sort(curSites.donors.begin(), curSites.donors.end());





		//add this entry to our map

		chr_splice_sites.insert(map<string, Splice_Sites>::value_type((*mapIt).first, curSites));



		inSpliceFile.close();



	}





}



//UPDATED: get DNA sequence regions for all groups (for the same query gene) from the sequence file
void DataManager::GetChromosomes(map<string, map<int, Group_RegEnd_Count> >& chr_group_region_count)
{
	if (VERBOSE)
		cout << "get chromosome sequences" << "\n";

	bool build_seq_index;
	if (chr_seq_map.empty())
		build_seq_index = true;
	else
		build_seq_index = false;

#ifdef DEBUG
	outFile << "build_seq_index is " << build_seq_index << "\n";
#endif

	int num_of_chr = chr_group_region_count.size();
	int chr_count = 0;

	string filename = chrSeqFile;
	//ifstream inSeqFile(filename.c_str());
	FILE* inSeqFile;
	inSeqFile = fopen(filename.c_str(), "r");

	//if (!inSeqFile.is_open())
	if (inSeqFile == NULL)
	{
		cout << "target sequence file open error: " << filename << "\n";
		return;
	}

	string curChrName; //line;
	char* line = new char[4096];
	vector<string> seq;
	//seq.push_back("");
	//char* seq;

	bool found = false;

	//char* whitespaces = " \t\n";
	unsigned int ws_pos;
	int seq_count = 0;
	int reg_start_pos, reg_end_pos, cur_pos;//, cur_len;
	map<string, map<int, Group_RegEnd_Count> >::iterator region_it;
	map<int, Group_RegEnd_Count>::iterator cur_reg_it;
	bool reg_start_found, reg_end_found;
	
	int line_len;
	//int seq_len;
	file_pos64 inSeqFile_pos;
	if (build_seq_index)
	{
#ifdef DEBUG
	outFile << "building seq_index and process first query gene" << "\n";
#endif
	//while (!inSeqFile.eof())
	while (!feof(inSeqFile))
	{
		//getline(inSeqFile, line);
		if (!fgets(line, 4096, inSeqFile)) //check if EOF and no char is read, then break
			break;
		//remove_trailing_cr_lf(line);
		line_len = remove_trailing_cr_lf(line, strlen(line)); //for dos/windows text files used in unix/linux
		//if (line.empty())
		if (line_len == 0)
			continue;

		//if (line.find('>') == 0)
		if (line[0] == '>')
		{
			if (found)
			{
				//chr_name_seq.insert(map<string, vector<string> >::value_type(curChrName, seq));
				StoreGroupSeq(seq, //seq_len, 
					(*cur_reg_it).second.group_count, reg_start_pos, seq_count);
				seq_count++;
				//seq.clear();
				//seq.push_back("");

#ifdef DEBUG
				outFile << "cur chr sequence stored:" << curChrName << ";seq_count:" << seq_count 
					<< ";chr_count:" << chr_count << "\n";
#endif
				chr_count++;
				if (chr_count == num_of_chr)
				{
					found = false;
					//break;
				}
			}

			//build_seq_index
#if !defined(IMPLICIT_LFS_64)
			fgetpos64(inSeqFile, &inSeqFile_pos);	
#else
			fgetpos(inSeqFile, &inSeqFile_pos);
#endif

			//ignore anything after whitespace chars
			/*ws_pos = line.find_first_of(whitespaces);
			if (ws_pos == string::npos)
				curChrName = line.substr(1, line.length() -1);
			else
				curChrName = line.substr(1, ws_pos-1);
			*/
			ws_pos = 0;
			while (ws_pos < line_len && line[ws_pos] != ' ' && line[ws_pos] != '\t')
				ws_pos++;
			line[ws_pos]='\0';//ignore whatever is after space
			string curChrName(line+1); //copy from the 2nd position to skip '>'
			chr_seq_map.insert(map<string,file_pos64>::value_type(curChrName, inSeqFile_pos));
#ifdef DEBUG
			char pos_str[64];
			sprintf(pos_str, "%lld", inSeqFile_pos);
			outFile << "curChrName is " << curChrName << "; position: " << pos_str << "\n";
#endif

			region_it = chr_group_region_count.find(curChrName);
			if (region_it != chr_group_region_count.end()) //need this chromosome
			{
				found = true;
				cur_reg_it = (*region_it).second.begin();
				reg_start_pos = (*cur_reg_it).first; //the first region start on this chromosome
				reg_end_pos = (*cur_reg_it).second.region_end;
				//seq = new char[reg_end_pos - reg_start_pos + 2]; //C string, with null terminator
				//seq[0]='\0'; //initialize to an empty C string with proper size
				//seq_len = 0;
				seq.clear();
				seq.push_back("");
				reg_start_found = false;
				reg_end_found = false;
				cur_pos = 0;
				outFile << "chr:" << curChrName << ";cur_region_start:" << reg_start_pos << ";cur_region_end:" << reg_end_pos << "\n";
			}
			else
				found = false;
		}
		else
			if (found)
			{
				//cur_pos += line.length();
				cur_pos += line_len;
				if (!GetRegion(cur_pos, reg_start_pos, reg_end_pos, seq, //seq_len,
					line, line_len, reg_start_found, reg_end_found))
					continue;

				if (reg_end_found) //region ended
				{
					//first, update group_dna_regions and group_start_seqno
					StoreGroupSeq(seq, //seq_len, 
						(*cur_reg_it).second.group_count, reg_start_pos, seq_count);
					seq_count++;
					//seq.clear();
					//seq.push_back("");

#ifdef DEBUG
				outFile << "cur chr sequence stored:" << curChrName << ";seq_count:" << seq_count 
					<< ";chr_count:" << chr_count << "\n";
#endif
					//reset for next region
					cur_reg_it++;
					if (cur_reg_it != (*region_it).second.end()) //next region in same chromosome
					{
						reg_start_pos = (*cur_reg_it).first;
						reg_end_pos = (*cur_reg_it).second.region_end;
						//seq = new char[reg_end_pos - reg_start_pos + 2]; //C string, with null terminator
						//seq[0] = '\0';
						//seq_len = 0;
						seq.clear();
						seq.push_back("");
						reg_start_found = false;
						reg_end_found = false;
						outFile << "chr:" << curChrName << ";cur_region_start:" << reg_start_pos << ";cur_region_end:" << reg_end_pos << "\n";
					}
					else
					{
						found = false;
						chr_count++;
						//if (chr_count == num_of_chr)
							//break;
					}
				}
			}
	}
	//if (VERBOSE)
	//	cout << "number of seq strings : " << seq.size() << "\n";
	//chr_name_seq.insert(map<string, vector<string> >::value_type(curChrName, seq));
	if (found)
	{
		StoreGroupSeq(seq, //seq_len, 
		(*cur_reg_it).second.group_count, reg_start_pos, seq_count);
#ifdef DEBUG
		outFile << "cur chr sequence stored:" << curChrName << ";seq_count:" << seq_count 
			<< ";chr_count:" << chr_count << "\n";
#endif
	}
	}
	else //if build_seq_index is already done, use it to set file position directly
	{
#ifdef DEBUG
	outFile << "seq_index already built, process subsequent query genes" << "\n";
#endif
		map<string,file_pos64>::iterator chr_seq_map_it;
		for (region_it = chr_group_region_count.begin(); region_it != chr_group_region_count.end(); region_it++)
		{
			chr_seq_map_it = chr_seq_map.find((*region_it).first);
			if (chr_seq_map_it == chr_seq_map.end())
			{
				cout << "sequence " << (*region_it).first << " not found in fasta file" << "\n";
				exit(-1);
			}
			inSeqFile_pos = (*chr_seq_map_it).second;
#if !defined(IMPLICIT_LFS_64)
			fsetpos64(inSeqFile, &inSeqFile_pos);	
#else
			fsetpos(inSeqFile, &inSeqFile_pos); //seek to the proper position of that specific chromosome
#endif
#ifdef DEBUG
			char pos_str[64];
			sprintf(pos_str, "%lld", inSeqFile_pos);
			outFile << "curChrName is " << (*region_it).first << "; position: " << pos_str << "\n";
#endif
			cur_pos = 0;
			for (cur_reg_it = (*region_it).second.begin(); cur_reg_it != (*region_it).second.end(); cur_reg_it++)
			{
#ifdef DEBUG
				outFile << "get region..." << "\n";
#endif
				reg_start_pos = (*cur_reg_it).first;
				reg_end_pos = (*cur_reg_it).second.region_end;
				//seq = new char[reg_end_pos - reg_start_pos + 2]; //C string, with null terminator
				//seq[0] = '\0';
				//seq_len = 0;
				seq.clear();
				seq.push_back("");
				reg_start_found = false;
				reg_end_found = false;

				while (!feof(inSeqFile))
				{
					if (!fgets(line, 4096, inSeqFile))
						break;
					line_len = remove_trailing_cr_lf(line, strlen(line)); //for dos/windows text files used in unix/linux
					if (line_len == 0)
						continue;
					if (line[0] == '>')
						break;
					cur_pos += line_len;
					if (!GetRegion(cur_pos, reg_start_pos, reg_end_pos, seq, //seq_len, 
						line, line_len, reg_start_found, reg_end_found))
						continue;
					if (reg_end_found) //region ended
						break;
				}
				StoreGroupSeq(seq, //seq_len, 
					(*cur_reg_it).second.group_count, reg_start_pos, seq_count);
				seq_count++;
#ifdef DEBUG
				outFile << "cur chr sequence stored:" << (*region_it).first << ";seq_count:" << seq_count << "\n";
#endif
			}
		}
	}

	delete [] line;

	//inSeqFile.close();
	fclose(inSeqFile);

/*	//read the chromosome file, get required chromosomes (fill in chr_name_seq map)
	ifstream inSeqFile(chrSeqFile.c_str());
	if (!inSeqFile.is_open())
	{
		cout << "target sequence file open error" << "\n";
		return;
	}

	string line;
	string name;
	string seq="";
	bool found=false;
	while (!inSeqFile.eof())
	{
		getline(inSeqFile, line);

//		cout << "line: " << line << "\n";

		if (line.empty())
			continue;

		if (line.find('>') == 0)
		{
			if (found)
			{
				cout << "chr_name: " << name << "\n";
				chr_name_seq.insert(map<string, string>::value_type(name, seq));
			}
			found = false;
			seq = "";
			name = line.substr(1);
			if (chr_names.find(name) != chr_names.end())
				found = true;
		}
		else
		{
			if (found)
			{
				seq += line;
//				cout << "string size: " << seq.size() << "\n";
			}
		}
	}

	inSeqFile.close();
*/
}

void DataManager::StoreGroupSeq(vector<string>& seq, 
								//char* seq, int seq_len, 
								vector<int>& group, int reg_start_pos, int seq_count)
{
	group_dna_regions.push_back(seq);//(pair<int, char*>(seq_len, seq));
	//vector<int>& group = (*cur_reg_it).second.group_count;
	for (int i=0; i<group.size(); i++)
	{
		group_start_seqno.insert(map<int, pair<int, int> >::value_type(group[i], pair<int,int>(reg_start_pos, seq_count)));
		outFile << "group_start_seqno:" << group[i] << ",<reg_start:" << reg_start_pos << ",seq#:" << seq_count << ">\n";
		outFile << "number of seq strings : " << seq.size() << "(" << LenOfStrVec(seq) << ")" << "\n";
		//outFile << "seq size:" << seq_len << "\n";
	}
}

bool DataManager::GetRegion(int cur_pos, int reg_start_pos, int reg_end_pos,
							vector<string>& seq, 
							//char* seq, int& seq_len, 
							char* line, int line_len, 
							bool& reg_start_found, bool& reg_end_found)
{
	int cur_len;
	if (!reg_start_found) //looking for region start
	{
		if (cur_pos < reg_start_pos)
			return false;
		else
		{
			reg_start_found = true;
			//int cur_line_start_pos = line.length()-1-(cur_pos - reg_start_pos);
			int cur_line_start_pos = line_len-1-(cur_pos - reg_start_pos);
			//cur_len = line.length() - cur_line_start_pos;
			cur_len = line_len - cur_line_start_pos;
			if (cur_pos >= reg_end_pos)
			{
				cur_len -= cur_pos - reg_end_pos;
				reg_end_found = true;
			}
			GetRegionSeq(seq, 
				//seq_len, 
				line, cur_line_start_pos, cur_len);
		}
	}
	else //region start already found
	{
		if (!reg_end_found) //looking for region end
		{
			int cur_line_start_pos = 0;
			if (cur_pos < reg_end_pos)
				//cur_len = line.length();
				cur_len = line_len;
			else
			{
				reg_end_found = true;
				//cur_len = line.length() - (cur_pos - reg_end_pos);
				cur_len = line_len - (cur_pos - reg_end_pos);
			}
			GetRegionSeq(seq, 
				//seq_len, 
				line, cur_line_start_pos, cur_len);
		}
	}
	return true;
}

void DataManager::GetRegionSeq(vector<string>& seq, 
							   //char* seq, int& seq_len, 
							   char* line, int cur_line_start_pos, int cur_len)
{
/*	char tmp_str[4096];
	int i=0, j=cur_line_start_pos;
	while (i<cur_len)
		tmp_str[i++] = line[j++];
	tmp_str[i] = '\0';
	strcat(seq, tmp_str);
	seq_len += cur_len;
*/
	string& last_string = seq.back();
	int leftover = last_string.length() + cur_len - MAX_LINE;
	if (leftover > 0)
	{
		//seq.back() += line.substr(cur_line_start_pos, cur_len-leftover);
		seq.back().append(line+cur_line_start_pos, cur_len-leftover);
		while (leftover > MAX_LINE)
		{
			//seq.push_back(line.substr(cur_line_start_pos+cur_len - leftover, MAX_LINE));
			string tmp_string1;
			tmp_string1.assign(line+cur_line_start_pos+cur_len-leftover, MAX_LINE);
			seq.push_back(tmp_string1);
			leftover -= MAX_LINE;
		}
		//seq.push_back(line.substr(cur_line_start_pos+cur_len-leftover));
		string tmp_string2;
		tmp_string2.assign(line+cur_line_start_pos+cur_len-leftover);
		seq.push_back(tmp_string2);
	}
	else
	{
		//seq.back() += line.substr(cur_line_start_pos, cur_len);
		seq.back().append(line+cur_line_start_pos, cur_len);
	}
}


void DataManager::GetChromosome(string& chr_name)
//modified: get 1 set of target/chromosome sequence for all input report files that correspond to that target
{
	if (chr_name_seq.find(chr_name) != chr_name_seq.end()) //if the previously-loaded chr is the one we need, return
		return;

	chr_name_seq.clear(); //reset chr_name_seq

	if (VERBOSE)
		cout << "get current chromosome sequence: " << chr_name << "\n";

	string filename = chrSeqFile;
	ifstream inSeqFile(filename.c_str());

	if (!inSeqFile.is_open())
	{
		cout << "target sequence file open error: " << filename << "\n";
		return;
	}

	string line, curChrName;
	vector<string> seq;
	seq.push_back("");
	bool found = false;

	const char* whitespaces = " \t\n";
	int ws_pos;

#ifdef DEBUG
	unsigned long line_no = 0;
#endif
	while (!inSeqFile.eof())
	{
		getline(inSeqFile, line);
#ifdef DEBUG
		line_no++;
#endif
		remove_trailing_cr_lf(line); //for dos/windows text files used in unix/linux

		if (line.empty())
			continue;

		if (line.find('>') == 0)
		{
			if (!found)
			{
				//ignore anything after whitespace chars
				ws_pos = line.find_first_of(whitespaces);
				if (ws_pos == string::npos)
					curChrName = line.substr(1, line.length() -1);
				else
					curChrName = line.substr(1, ws_pos-1);

				if (curChrName.compare(chr_name) == 0)
				{
					found = true;
#ifdef DEBUG
					outFile << "found sequence:" << curChrName << "\n";
#endif
				}
				else
				{
#ifdef DEBUG
					outFile << "current sequence:" << curChrName << " at line " << line_no << "\n";
#endif
				}

				continue;
			}
			else
			{
#ifdef DEBUG
				outFile << "sequence finished at line " << line_no << "\n";
#endif
				//break out of the loop immediately, as we are good after finding the single chromosome we want
				break;
			}
		}

		if (found)
		{
			string& last_string = seq.back();

			int leftover = last_string.length() + line.length() - MAX_LINE;
			if (leftover > 0)
			{
				seq.back() += line.substr(0, line.length()-leftover);

				//Modified (bug fix)! Consider the possibility that single input line exceeds MAX_LINE!
				while (leftover > MAX_LINE)
				{
					seq.push_back(line.substr(line.length()-leftover, MAX_LINE));
					leftover -= MAX_LINE;
				}
				seq.push_back(line.substr(line.length()-leftover));
			}
			else
				seq.back() += line;
		}
	}

	if (VERBOSE)
		cout << "number of seq strings : " << seq.size() << "\n";

#ifdef DEBUG
	outFile << "number of seq strings : " << seq.size() << "\n";
#endif

	if (seq[0].length() == 0)
	{
		cout << "get chromosome sequence (" << chr_name << ") failed, check your sequence file" << "\n";
		exit(-1);
	}

	chr_name_seq.insert(map<string, vector<string> >::value_type(chr_name, seq)); //only one chromosome

	inSeqFile.close();
}


//now read ".align" file and store alignment info for HSPs in HSP_IDs

void DataManager::GetAlignments(set<int>& HSP_IDs)

{

//	cout << "get blast alignments" << "\n";



	if (!alignFile_Open)

	{

		alignFile_is.open(alignFile.c_str());

		if (!alignFile_is.is_open())

		{

			cout << "input alignment file open error" << "\n";

			return;

		}

		alignFile_Open = true;

	}



	string line;

	int id=1; //HSP_ID start with "1", has one-to-one correspondence with HSP in report file

	int pos1, pos2;

	getline(alignFile_is, line); //skip 1 line (first time, it skips the 1st TBLASTN line; later times, skips a useless line)

	//cur_align_start = true;

	while (!alignFile_is.eof())	//read in file, parse, 

	{

		getline(alignFile_is, line);
		remove_trailing_cr_lf(line);

		//if (query_gene.compare("F39C12.2a")==0)

		//	outFile << line << "\n";



		//if (cur_align_start && line.find("TBLASTN") == 0) //signal of next batch of HSPs

		if (line.find("TBLASTN") == 0) //signal of next batch of HSPs

		{

			//cur_align_start = false;			

			break;

		}



		int pos = line.find("Query:");

		if (pos == 0) //header line, starts with '>'

		{

			//if (query_gene.compare("F39C12.2a")==0)

			//	outFile << "alignment: id" << id << "\n";



			//int tab_pos = line.find('\t');

			//int id = atoi(line.substr(pos+1, tab_pos-pos-1).c_str());

			

			if (HSP_IDs.find(id) != HSP_IDs.end()) //found this ID, we need to store its alignment then

			{

				//if (query_gene.compare("F39C12.2a")==0)

				//	outFile << "needed: id" << id << "\n";



				string qStr="", matchStr="", tStr="";

				while (pos == 0)

				{

				//getline(alignFile_is, line);

				pos1 = line.find_first_of("0123456789"); //digit

				pos1 = line.find_first_not_of(" 0123456789", pos1); //skip digits, spaces

				pos2 = line.find(' ', pos1);

				//string tStr = line; //first line is target string

				qStr += line.substr(pos1, pos2-pos1);

				getline(alignFile_is, line);
				remove_trailing_cr_lf(line);

				//if (query_gene.compare("F39C12.2a")==0)

				//	outFile << "matchStr: " << line << "\n";



				//string matchStr = line; //similarity between two strings

				matchStr += line.substr(pos1, pos2-pos1);

				getline(alignFile_is, line);
				remove_trailing_cr_lf(line);


				//if (query_gene.compare("F39C12.2a")==0)

				//	outFile << "targetStr: " << line << "\n";



				//string qStr = line; //second line is query string

				tStr += line.substr(pos1, pos2-pos1);

				

				getline(alignFile_is, line); //empty line

				getline(alignFile_is, line);
				remove_trailing_cr_lf(line);

				pos = line.find("Query:");

				}

				input_alignments.insert(map<int, Input_Alignment>::value_type(id, Input_Alignment(qStr, tStr, matchStr)));



				//outFile << "alignment collected for hsp_" << id << "\n";

			}

			else //this HSP is not used, skip it

			{

				while (pos == 0)

				{

					getline(alignFile_is, line); //match line

					getline(alignFile_is, line); //target line

					getline(alignFile_is, line); //empty line

					getline(alignFile_is, line); //next possible query line
					remove_trailing_cr_lf(line);

					pos = line.find("Query:");

				}



/*				getline(alignFile_is, line);

				while ( line.find("Score") == string::npos)

					getline(alignFile_is, line);

*/			}

			id++;

		}

	}



	if (alignFile_is.eof())

	{

		alignFile_is.close();

		alignFile_Open = false;

		//cur_align_start = false;

		//alignFile_Finish = true;

	}



	cout << "obtained current blast alignments" << "\n";



}



//now read ".align" file and store alignment info for HSPs in HSP_IDs

void DataManager::GetAlignments_Skip()

{

	cout << "current gene: " << query_gene << "\n";



//	cout << "get blast alignments" << "\n";



	if (!alignFile_Open)

	{

		alignFile_is.open(alignFile.c_str());

		if (!alignFile_is.is_open())

		{

			cout << "input alignment file open error" << "\n";

			return;

		}

		alignFile_Open = true;

	}



	string line;

	getline(alignFile_is, line); //skip 1 line (first time, it skips the 1st TBLASTN line; later times, skips a useless line)

	//cur_align_start = true;

	while (!alignFile_is.eof())	//read in file, parse, 

	{

		getline(alignFile_is, line);
		remove_trailing_cr_lf(line);

//		if (query_gene.compare("F08D12.9")==0)

//			outFile << line << "\n";



		//if (cur_align_start && line.find("TBLASTN") == 0) //signal of next batch of HSPs

		if (line.find("TBLASTN") == 0) //signal of next batch of HSPs

		{

			//cur_align_start = false;			

			break;

		}



	}



	if (alignFile_is.eof())

	{

		alignFile_is.close();

		alignFile_Open = false;

		//cur_align_start = false;

		//alignFile_Finish = true;

	}



	cout << "skipped current blast alignments" << "\n";



}



int	DataManager::GetHSPID(bool isPosStrand, int chr_index, int index)

{

	if (isPosStrand)

		return HSP_gene[chr_index][index].ID;

	else

		return HSP_neg_gene[chr_index][index].ID;

}



void DataManager::LoadQuerySeqData(const char* filename) //filename must contain query sequences, in FASTA format

{

	query_gene_seq.clear();



	ifstream inFile(filename);

	if (!inFile.is_open())

	{

		cout << "cannot find query sequence file: " << filename << "\n";

		exit(-1);

	}



	string line, gene, seq;

	gene = "";

	seq = "";

	int para_pos;

	while (!inFile.eof())

	{

		getline(inFile, line);
		remove_trailing_cr_lf(line);


		if (line.empty() || all_white_space(line)) //skip empty lines

			continue;

		//cout << "1:" << line << ":" << "\n";

		if ((para_pos = line.find_first_of(" \t\n\r")) != string::npos)
			line.erase(para_pos); //remove trailing white spaces
		//CombineConsecutiveSpaces(line);

		//cout << "2:" << line << ":" << "\n";

		if (line.find('>') == 0)

		{

			if (seq.length() > 0)
			{
				query_gene_seq.insert(map<string, string>::value_type(gene, seq));
				//cout << "obtained query_gene:" << gene << ":\n" << seq << "\n";
			}

			gene = line.substr(1);

			seq = ""; //reset

		}

		else

			seq += line;



	}

	if (seq.length()>0)

		query_gene_seq.insert(map<string, string>::value_type(gene, seq)); //last one



	inFile.close();

}



void DataManager::GetQuerySeq(string& query_str)

{
	map<string, string>::iterator it;
	
	//for (it = query_gene_seq.begin(); it != query_gene_seq.end(); it++)
	//	cout << "gene:" << (*it).first << ":seq:" << (*it).second << "\n";

	it = query_gene_seq.find(query_gene);


	if (it == query_gene_seq.end())

	{

		cout << "cannot find this query gene in query sequences: " << query_gene << " in " << QUERY_SEQ_FILE << "\n";

		exit(-1);

	}



	query_str = (*it).second;



}



/*

void DataManager::Print_Alignments(ostream& os, Alignment* am, vector<HSP_Gene_Pair*>& HSPs, bool isPosStrand)

{

		//print query_markers

		am->Print_Markers(os, am->query_markers, '<', '>');



		//print query_string

		os << "Query\t" ;

		char tmp[32];

		_itoa(HSPs.front()->gene_start, tmp, 10);

		os << tmp;



		string tmp_str = tmp;

		int size = tmp_str.size();

		int i;

		for (i=size; i<10; i++)

			os << " ";



		am->Print_Align_Str(os, am->query_layout, am->query_overlap_align_gap, am->target_overlap_align_gap);

		os << "\t" << HSPs.back()->gene_end << "\n";





		//print target_string

		os << "Target\t";

		int pos1, pos2;

		if (isPosStrand)

		{

			pos1 = HSPs.front()->HSP_start;

			pos2 = HSPs.back()->HSP_end;

		}

		else

		{

			pos1 = HSPs.front()->HSP_end;

			pos2 = HSPs.back()->HSP_start;

		}

		_itoa(pos1, tmp, 10);

		os << tmp;

		tmp_str = tmp;

		size = tmp_str.size();

		for (i=size; i<10; i++)

			os << " ";



		am->Print_Align_Str(os, am->target_layout, am->target_overlap_align_gap, am->query_overlap_align_gap);

		os << "\t" << pos2 << "\n";



		//print target_markers

		am->Print_Markers(os, am->target_markers, '[', ']');

}



//TODO:

//each marker should have enough space (adjust query_layout and target_layout if necessary)

//each marker should record the relative start and end position w.r.t. its current alignment string

//each marker should record its original coordinates

//use ".." for gaps between query_layout / target_layout (adjust marker position accordingly)

//use "//" before each copy of previous segments (adjust marker position accordingly)

//mark position for every start and end (above all markers)

Alignment* DataManager::ConstGroupAlignment(vector<HSP_Gene_Pair*>& HSPs)

{

	Alignment* am = new Alignment();

	int size = HSPs.size();

	int i;

	int query_marker_index = 1;

	int target_marker_index = 1;

	for (i=0; i < size; i++) //for each HSP_Gene_Pair in HSPs

	{

		int id = HSPs[i]->ID;



		if (id == 89)

			int stop =1;



		map<int, Input_Alignment>::iterator mapIt  = input_alignments.find(id);

		if (mapIt == input_alignments.end())

		{

			cout << "cannot find HSP[" << id << "] in input alignments" << "\n";

			continue;

		}



		Input_Alignment inAlign = mapIt->second;

		int curhsp_size = inAlign.query_align.size();



		am->query_layout.push_back(inAlign.query_align);

		am->target_layout.push_back(inAlign.target_align);



		int gOL_start = 0;

		int gOL_end = 0;

		int hOL_start = 0;

		int hOL_end = 0;

		bool g_align = false;

		bool h_align = false;



		int prev_index = i-1;

		if (prev_index >=0 )

		{

			//set up flags for "..._overlap_align_gap"

			bool prev_ol = HSPs[i]->Overlap(HSPs[prev_index], gOL_start, gOL_end, hOL_start, hOL_end, g_align, h_align);



			//set up next_start_pos

			if (g_align && h_align)

				am->next_start_pos.push_back(am->next_start_pos[prev_index]+curhsp_size);

			else

			{

				am->next_start_pos[prev_index] += 2;//adding "//" or ".."

				am->next_start_pos.push_back(am->next_start_pos[prev_index]+curhsp_size); 

			}



			if (gOL_end != 0) //query overlap, need markers

			{

				am->query_overlap_align_gap.push_back(-1);



				AssignMarker(am, i, prev_index, true, HSPs[prev_index]->gene_end, gOL_start, gOL_end, query_marker_index);

				query_marker_index++;

			}

			else //no need for markers

				if (g_align)

					am->query_overlap_align_gap.push_back(0);

				else

					am->query_overlap_align_gap.push_back(1);



			if (hOL_end != 0) //target overlap, need markers

			{

				am->target_overlap_align_gap.push_back(-1);



				AssignMarker(am, i, prev_index, false, HSPs[prev_index]->HSP_end, hOL_start, hOL_end, target_marker_index);

				target_marker_index++;

			}

			else //no need for markers

				if (h_align)

					am->target_overlap_align_gap.push_back(0);

				else

					am->target_overlap_align_gap.push_back(1);



			//check previous ones if there is still overlap

			prev_index--;

			while (prev_index >=0 && prev_ol)

			{

				gOL_start = 0;

				gOL_end = 0;

				hOL_start = 0;

				hOL_end = 0;

				//g_align = false;

				//h_align = false;



				prev_ol = HSPs[i]->Overlap(HSPs[prev_index],  gOL_start, gOL_end, 

					hOL_start, hOL_end, g_align, h_align);



				//only need to check for overlap (cannot be align or gap)

				if (gOL_end != 0)

				{

					AssignMarker(am, i, prev_index, true, HSPs[prev_index]->gene_end, gOL_start, gOL_end, query_marker_index);

					query_marker_index++;

				}



				if (hOL_end != 0)

				{

					AssignMarker(am, i, prev_index, false, HSPs[prev_index]->HSP_end, hOL_start,hOL_end, target_marker_index);

					target_marker_index++;

				}

				

				prev_index--;



			}



		}

		else //first hsp

			am->next_start_pos.push_back(curhsp_size);



	}



	//arrange line_index in markers

	am->Arrange_Markers();



	return am;

}





void DataManager::AssignMarker(Alignment* am, int cur_index, int prev_index, bool isQuery, int prev_HSP_end,

							   int ol_start, int ol_end, int cur_marker_index)

{

	//fix prev_align marker

	int prev_align_len;

	if (isQuery)

		prev_align_len = am->query_layout[prev_index].size();

	else

		prev_align_len = am->target_layout[prev_index].size();



	int e; // = prev_align_len-1 ; //the end of overlap must be always the last char (e is the relative distance from the start of this alignment)

	int i=0;

	int step = isQuery ? 3 : 1;

	int j=prev_align_len - 1;

	while (i < ol_end - prev_HSP_end)

	{

		char cur_char;

		if (isQuery)

			cur_char = am->query_layout[prev_index].at(j);

		else

			cur_char = am->target_layout[prev_index].at(j);

		if ((cur_char <= 'z' && cur_char >= 'a') || (cur_char <='Z' && cur_char >= 'A')) //skip all non alphbet chars

			i++;



		j -= step;

	}

	e = j;



	i=0;

	if (isQuery)

		j -= 2;

	while ( i<ol_end-ol_start )

	{

		j -= step;

		char cur_char;

		if (isQuery)

			cur_char = am->query_layout[prev_index].at(j);

		else

			cur_char = am->target_layout[prev_index].at(j);

		if ((cur_char <= 'z' && cur_char >= 'a') || (cur_char <='Z' && cur_char >= 'A')) //skip all non alphbet chars

			i++;



	}

	//int s = e - ol_end + ol_start-(prev_align_len - (prev_end-prev_start+1));//actual end of the prev alignment string

	int s = j; //actual start of the overlap



				

	char marker_index_chars[32];

	string marker_str;

	_itoa(cur_marker_index, marker_index_chars, 10);

	marker_str = marker_index_chars;

	int total_len = marker_str.size()+2; //this is the string length of the marker_index string



	int dif1 = total_len - 1 - (e-s);

	if (dif1 > 0)

	{

		e += dif1;

	

		int k;

		for (k=prev_index; k<=cur_index; k++)

			am->next_start_pos[k] += dif1;



	}

	else

		dif1 = 0; 



	//now fix the cur_align marker

	int cur_s, cur_e;



	cur_s = 0; //must always be the start of current alignment

	i=0;

	j=-1;

	while ( i<ol_end-ol_start+1 )

	{

		j+=step;

		char cur_char;

		if (isQuery)

			cur_char = am->query_layout[cur_index].at(j);

		else

			cur_char = am->target_layout[cur_index].at(j);

		if ((cur_char <= 'z' && cur_char >= 'a') || (cur_char <='Z' && cur_char >= 'A')) //skip all non alphbet chars

			i++;

	}

	//cur_e = cur_s + ol_end - ol_start + dif;

	cur_e = j;



	int dif2 = total_len - 1 - (cur_e-cur_s);

	if (dif2 > 0)

	{

		cur_e += dif2;

	

		am->next_start_pos[cur_index] += dif2;



	}

	else

		dif2 = 0; 





	//now stick spaces at the back of prev_align and at the front of cur_align if necessary

	//also fix the alignment_markers finally

	for (j=0; j<dif1; j++)

	{

		am->query_layout[prev_index] += " ";

		am->target_layout[prev_index] += " ";

	}

	for (j=0; j<dif2; j++)

	{

		am->query_layout[cur_index].insert(0, " ");

		am->target_layout[cur_index].insert(0, " ");

	}



	if (isQuery)

	{

		am->query_markers.push_back(Align_Marker(prev_index, s, e, cur_marker_index, ol_start, ol_end));

		am->query_markers.push_back(Align_Marker(cur_index, cur_s, cur_e, cur_marker_index, ol_start, ol_end));

	}

	else

	{

		am->target_markers.push_back(Align_Marker(prev_index, s, e, cur_marker_index, ol_start, ol_end));

		am->target_markers.push_back(Align_Marker(cur_index, cur_s, cur_e, cur_marker_index, ol_start, ol_end));

	}

}





void Alignment::Arrange_Markers()

{

	//sort before doing anything further (now every marker has line_index 0)

	sort(query_markers.begin(), query_markers.end());

	sort(target_markers.begin(), target_markers.end());



	//now fix line_index in Align_Marker

	vector<int> line_indexes;

	vector<Align_Marker>::iterator it;

	for (it=query_markers.begin(); it != query_markers.end(); it++)

		//it->line_index = Last_line_Of_Overlap_Marker( query_markers, it) + 1;

		line_indexes.push_back (Last_line_Of_Overlap_Marker( query_markers, it) + 1); //use this, just to make sure markers aren't sorted again automatically

	int i=0;

	for (it=query_markers.begin(); it != query_markers.end(); it++, i++)

		it->line_index = line_indexes[i];



	line_indexes.clear();

	for (it=target_markers.begin(); it != target_markers.end(); it++)

		//it->line_index = Last_line_Of_Overlap_Marker( target_markers, it) + 1;

		line_indexes.push_back(Last_line_Of_Overlap_Marker( target_markers, it) + 1);

	i=0;

	for (it=target_markers.begin(); it != target_markers.end(); it++, i++)

		it->line_index = line_indexes[i];



	//sort again, this time with line_index filled

	sort(query_markers.begin(), query_markers.end());

	sort(target_markers.begin(), target_markers.end());



}



int Alignment::Last_line_Of_Overlap_Marker(vector<Align_Marker>& markers, vector<Align_Marker>::iterator cur_it)

{

	int line_index = -1;



	int cur_align_index = cur_it->align_index;

	int cur_end_pos = cur_it->end_pos;

	int cur_start_pos = cur_it->start_pos;



	vector<Align_Marker>::iterator lower_it = lower_bound(markers.begin(), markers.end(), 

		Align_Marker(cur_align_index,0,0,0,0,0 ));

//	vector<Align_Marker>::iterator upper_it = lower_bound(markers.begin(), markers.end(), 

//		Align_Marker(cur_align_index,cur_end_pos+1, 0,0,0,0)); 



	while (lower_it != cur_it)

	{

		//everything in this range must have start_pos<=cur_start_pos

		if ( lower_it->end_pos >= cur_start_pos  ) //for each overlapping marker

			line_index = lower_it->line_index > line_index ? lower_it->line_index : line_index; //find the largest line_index



		lower_it++;

	}



	return line_index;



}

*/

