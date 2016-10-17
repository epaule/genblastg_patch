#include "scores.h"

#include "c45extern.h"


#include <stack>

using namespace std;



//gene score / HSP score (node->score)

float Score(int start, int end, float pid)

{

//	return pid*pid/10000.0f*(float)(end-start+1);

	return pid/100.0f*(float)(end-start+1); //*HSP_SK_ALPHA; //GENE_MS_BETA

}



//edge penalty

float Penalty_old(float HSP_skip_penalty, float gene_missing_penalty)

{

	return (HSP_ALPHA*HSP_skip_penalty+GENE_BETA*gene_missing_penalty );



}



//edge penalty

//Sep. 8, 2007 Note: on one edge, it can be only either gene missing or overlapping, not both, 

//so right now we treat both with same GENE_MS_BETA parameter

float Penalty_new(float HSP_skip_penalty, float gene_missing_penalty, float gene_overlap_penalty )

{

//	return (HSP_SK_ALPHA*HSP_skip_penalty+GENE_MS_BETA*gene_missing_penalty + GENE_OL_GAMA*gene_overlap_penalty);

	//return (HSP_SK_ALPHA*HSP_skip_penalty+GENE_MS_BETA*gene_missing_penalty + HSP_SK_ALPHA*gene_overlap_penalty);



	return (HSP_SK_ALPHA*HSP_skip_penalty+GENE_MS_BETA*gene_missing_penalty + GENE_MS_BETA*gene_overlap_penalty);

	//return (HSP_skip_penalty+GENE_MS_BETA*gene_missing_penalty + GENE_MS_BETA*gene_overlap_penalty); //skip_penalty already computed based on Score(), which takes account HSP_SK_ALPHA

}





/* utility function used by DataManger.gene_start_end_map and ACCP_DONR_graph.cur_start_HSP_map */

void MyMapInsert( const int pos, const int index, const bool isStart, map<int, vector<int> >& myMap)

{

	map<int, vector<int> >::iterator sit = myMap.lower_bound(pos);

	if (sit != myMap.end() && !myMap.key_comp()(pos, (*sit).first)) //already has this key

	{

		if (isStart)

			(*sit).second.push_back(index); //record which lines start from this position

	}

	else

	{

		if (isStart)

			myMap.insert(sit, map<int, vector<int> >::value_type(pos, vector<int>(1, index)));//construct the vector with 1 integer: count

		else

			myMap.insert(sit, map<int, vector<int> >::value_type(pos, vector<int>()));

	}

}



void MyMapInsert( const string sort_key, const long value, map<string, multiset<long> >& myMap)
{
	map<string, multiset<long> >::iterator sit = myMap.lower_bound(sort_key);
	if (sit != myMap.end() && !myMap.key_comp()(sort_key, (*sit).first)) //already has this key
	{
		(*sit).second.insert(value); 
	}
	else
	{
		multiset<long> temp_set;
		temp_set.insert(value);
		myMap.insert(sit, map<string, multiset<long> >::value_type(sort_key, temp_set));//construct the vector with 1 integer: count
	}
}

void MyMapInsert( const string sort_key, const long value, multimap<string, multiset<long> >& myMap)
{
	multimap<string, multiset<long> >::iterator sit = myMap.lower_bound(sort_key);
	if (sit != myMap.end() && !myMap.key_comp()(sort_key, (*sit).first)) //already has this key
	{
		(*sit).second.insert(value); 
	}
	else
	{
		multiset<long> temp_set;
		temp_set.insert(value);
		myMap.insert(sit, multimap<string, multiset<long> >::value_type(sort_key, temp_set));//construct the vector with 1 integer: count
	}
}

void MyMapInsert( const int pos, const int pos_end, const float PID, const bool isStart, map<int, vector<pair<float, int> > >& myMap)
{

	map<int, vector<pair<float, int> > >::iterator sit = myMap.lower_bound(pos);

	if (sit != myMap.end() && !myMap.key_comp()(pos, (*sit).first)) //already has this key

	{

		if (isStart)

			(*sit).second.push_back(pair<float, int>(PID, pos_end) ); //record which lines start from this position

	}

	else

	{

		if (isStart)

			myMap.insert(sit, map<int, vector<pair<float, int> > >::value_type(pos, vector<pair<float, int> >(1, pair<float, int>(PID, pos_end) )));//construct the vector with 1 pair: count

		else

			myMap.insert(sit, map<int, vector<pair<float, int> > >::value_type(pos, vector<pair<float, int> >())); //empty vector

	}

}

//for WormBase_Exons map
void MyMapInsert( const string sort_key, const string chr_seq, const int start, const int end, 
				 map<string, pair<string, set< pair<int, int> > > >& myMap)
{
	map<string, pair<string, set< pair<int, int> > > >::iterator sit = myMap.lower_bound(sort_key);
	if (sit != myMap.end() && !myMap.key_comp()(sort_key, (*sit).first)) //already has this key
	{
		(*sit).second.second.insert( pair<int, int>(start, end) ); 
	}
	else
	{
		set< pair<int, int> > temp_set;
		temp_set.insert( pair<int,int>(start, end) );
		myMap.insert(sit, map<string, pair<string, set< pair<int, int> > > >::value_type(sort_key, 
			pair<string, set< pair<int, int> > >(chr_seq, temp_set)));
	}
}


int CountCharInStr(const string& str, char ch)

{

	int count = 0;

	int pos = str.find(ch);

	while (pos != string::npos)

	{

		count++;

		pos = str.find(ch, pos+1);

	}

	return count;

}





void GetSubstrFromVecStrs(vector<string>& chr_seq, bool isPosStrand, const int start_pos, const int len, string& resultStr) 
//start_pos is index (starting from 0)
{

//	cout << "getting from " << start_pos << " for " << len << "(" << isPosStrand << ")" << "\n";

	resultStr = "";

	if (len <= 0)

		return;

	int cur_start_pos = start_pos;
	cur_start_pos -= chromosome_start_pos - 1; //UPDATE cur_start_pos by adjusting using chromosome_start_pos

	int len_need = len;
	if (cur_start_pos <0)
	{
		cout << "error4: start_pos is negative (" << start_pos << ") " << cur_start_pos << "\n";
		return;
		//len_need += start_pos; //this should only be done during repair?
		//cur_start_pos = 0;
	}



	//vector<string>::iterator it = chr_seq.begin();

	//int i=0;

	int count = 0;


	int size = chr_seq[0].length();

	int i = cur_start_pos / size;

	int j = cur_start_pos % size;

	if (i >= chr_seq.size()) //error checking

	{

		cout << "error5: start_pos (" << start_pos << ") is bigger than sequence length " << cur_start_pos << "\n";

		return;

	}

	int left = chr_seq[i].length() - j;

	bool error = false;

	while (left < len_need)

	{

		resultStr += chr_seq[i].substr(j, left);

		len_need -= left;

		i++;

		if (i >= chr_seq.size()) //error checking again

		{

			cout << "GetSubstrFromVecStrs error: not enough (" << len << ") chars after start_pos (" << start_pos << ")" 
				<< cur_start_pos << ";" 
				<< i << ">=" << chr_seq.size() << "\n";

			exit(-1);

			error = true;

			break;

		}

		j=0;

		left = chr_seq[i].length();

	}

	if (!error)

		resultStr += chr_seq[i].substr(j, len_need);



	//while (it != chr_seq.end())

/*	while (i<chr_seq.size())

	{

		//size = (*it).size();

		size = chr_seq[i].size();

		count += size;

		if (count <= cur_start_pos)

		{			

			//it++;

			i++;

			continue;

		}



		int left = count - cur_start_pos ;

		if (left >= len_need)

		{

			//resultStr += (*it).substr(size-left, len_need);

			resultStr += chr_seq[i].substr(size-left, len_need);

			break;

		}

		else

		{

			//resultStr += (*it).substr(size-left, left);

			resultStr += chr_seq[i].substr(size-left, left);

			cur_start_pos += left;

			len_need -= left;

			//it++;

			i++;

		}

	}

	if (count <= start_pos)

	{

		cout << "error: start_pos (" << start_pos << ") is bigger than sequence length (" << count << ")" << "\n";

		return;

	}

*/

	//for negative genes, complement the resultStr

	if (!isPosStrand)

		StrToComplement(resultStr);

		//StrToReverseComplement(resultStr); //Update!!!

}

//allows to get string that is not up to "len" long
int GetSubstrFromVecStrs_ForRepair(vector<string>& chr_seq, bool isPosStrand, const int start_pos, const int len, string& resultStr) 
//start_pos is index (starting from 0)
{
//	cout << "getting from " << start_pos << " for " << len << "(" << isPosStrand << ")" << "\n";
	resultStr = "";
	if (len <= 0)
		return 0; //0 indicates there's no valid string returned

	int cur_start_pos = start_pos;
	cur_start_pos -= chromosome_start_pos - 1;

	int len_need = len;
	if (cur_start_pos <0)
	{
		//cout << "error2: start_pos is negative (" << start_pos << ")" << "\n";
		//return 0;
		len_need += cur_start_pos;
		cur_start_pos = 0;		
	}

	int count = 0;
	int size = chr_seq[0].length();
	int i = cur_start_pos / size;
	int j = cur_start_pos % size;
	if (i >= chr_seq.size()) //error checking
	{
		//cout << "error6: start_pos (" << start_pos << ") is bigger than sequence length" << "\n";
		return 0;
	}
	int left = chr_seq[i].length() - j;
	bool error = false;
	while (left < len_need)
	{
		resultStr += chr_seq[i].substr(j, left);
		len_need -= left;
		i++;
		if (i >= chr_seq.size()) //error checking again
		{
			//cout << "error: not enough (" << len << ") chars after start_pos (" << start_pos << ")" 
			//	<< i << ">=" << chr_seq.size() << "\n";
			//exit(-1);
			error = true;
			break;
		}
		j=0;
		left = chr_seq[i].length();
	}
	if (!error)
		resultStr += chr_seq[i].substr(j, len_need);

	//for negative genes, complement the resultStr
	if (!isPosStrand)
		StrToComplement(resultStr);
		//StrToReverseComplement(resultStr); //Update!!!

	return resultStr.length();
}


void GetSubstrFromVecStrs(vector<string>& chr_seq, bool isPosStrand, const int start_pos, const int len, string& resultStr, 
						  string& prev_2nt, string& next_2nt) //start_pos is index (starting from 0)
{

//	cout << "getting from " << start_pos << " for " << len << "(" << isPosStrand << ") 2nts" << "\n";

	resultStr = "";
	if (len <= 0)
		return;

	int count = 0;
	int cur_start_pos = start_pos;
	cur_start_pos -= chromosome_start_pos - 1;

	if (cur_start_pos <0)
	{
		cout << "error7: start_pos is negative (" << start_pos << ")" << cur_start_pos << "\n";
		return;
	}


	int size;

	int len_need = len;



	size = chr_seq[0].length();

	int i = cur_start_pos / size;

	int j = cur_start_pos % size;

	if (i >= chr_seq.size()) //error checking

	{

		cout << "error8: start_pos (" << start_pos << ") is bigger than sequence length)" << cur_start_pos << "\n";

		return;

	}

	int left = chr_seq[i].length() - j;

	bool error = false;

	if (j >= 2)

		prev_2nt = chr_seq[i].substr(j-2, 2);

	else

		if (i > 0) //has chr_seq[i-1]

		{

			if (j == 1)

			{

				prev_2nt = chr_seq[i-1][chr_seq[i-1].length()-1];

				prev_2nt += chr_seq[i][0];

			}

			else

			{

				prev_2nt = chr_seq[i-1].substr(chr_seq[i-1].length()-2, 2);

			}

		}

		else

		{

			if (j == 1)

			{

				prev_2nt = ' ';

				prev_2nt += chr_seq[i][0];

			}

			else

			{

				prev_2nt = "  ";

			}

		}



	while (left < len_need)

	{

		resultStr += chr_seq[i].substr(j, left);

		len_need -= left;

		i++;

		if (i >= chr_seq.size()) //error checking again

		{

			cout << "GetSubstrFromVecStrs (2 extras) error: not enough (" << len << ") chars after start_pos (" << start_pos << ")" 
				<< cur_start_pos << ";" 
				<< i << ">=" << chr_seq.size() << "\n";

			error = true;

			exit(-1);

			break;

		}

		j=0;

		left = chr_seq[i].length();

	}

	if (!error)

	{

		resultStr += chr_seq[i].substr(j, len_need);

		if (j+len_need+2 <= chr_seq[i].length())

			next_2nt = chr_seq[i].substr(j+len_need, 2);

		else

			if (j+len_need == chr_seq[i].length() - 1)

			{

				next_2nt = chr_seq[i][chr_seq[i].length()-1];

				if (i < chr_seq.size() - 1)

					next_2nt += chr_seq[i+1][0];

				else

					next_2nt += ' ';

			}

			else

			{

				if (i < chr_seq.size() - 1)

					if (chr_seq[i+1].length() >= 2)

					{

						next_2nt = chr_seq[i+1].substr(0, 2);

					}

					else

					{

						next_2nt = chr_seq[i+1][0];

						next_2nt += ' ';

					}

				else

					next_2nt = "  ";

			}



	}



	//while (it != chr_seq.end())

/*	while (i<chr_seq.size())

	{

		//size = (*it).size();

		size = chr_seq[i].size();

		count += size;

		if (count <= cur_start_pos)

		{			

			//it++;

			i++;

			continue;

		}



		int left = count - cur_start_pos ;

		if (left >= len_need)

		{

			//resultStr += (*it).substr(size-left, len_need);

			resultStr += chr_seq[i].substr(size-left, len_need);

			break;

		}

		else

		{

			//resultStr += (*it).substr(size-left, left);

			resultStr += chr_seq[i].substr(size-left, left);

			cur_start_pos += left;

			len_need -= left;

			//it++;

			i++;

		}

	}

	if (count <= start_pos)

	{

		cout << "error: start_pos (" << start_pos << ") is bigger than sequence length (" << count << ")" << "\n";

		return;

	}

*/

	//for negative genes, complement the resultStr

	if (!isPosStrand)

	{

		StrToComplement(resultStr);

		StrToComplement(prev_2nt);

		StrToComplement(next_2nt);

		//StrToReverseComplement(resultStr); //Update!!!

	}

}


void GetSubstrFromVecStrs_ForRepair(vector<string>& chr_seq, bool isPosStrand, const int start_pos, const int len, string& resultStr, 
						  string& prev_2nt, string& next_2nt) //start_pos is index (starting from 0)
{
//	cout << "getting from " << start_pos << " for " << len << "(" << isPosStrand << ") 2nts" << "\n";
	resultStr = "";
	if (len <= 0)
		return;
	int cur_start_pos = start_pos;
	cur_start_pos -= chromosome_start_pos - 1;

	int len_need = len;
	if (cur_start_pos <0)
	{
		//cout << "error4: start_pos is negative (" << start_pos << ")" << "\n";
		//return;
		len_need += cur_start_pos;
		cur_start_pos = 0;
	}

	//vector<string>::iterator it = chr_seq.begin();
	//int i=0;
	int count = 0;
	int size = chr_seq[0].length();
	int i = cur_start_pos / size;
	int j = cur_start_pos % size;
	if (i >= chr_seq.size()) //error checking
	{
		//cout << "error9: start_pos (" << start_pos << ") is bigger than sequence length" << "\n";
		return;
	}
	int left = chr_seq[i].length() - j;

	bool error = false;
	if (j >= 2)
		prev_2nt = chr_seq[i].substr(j-2, 2);
	else
		if (i > 0) //has chr_seq[i-1]
		{
			if (j == 1)
			{
				prev_2nt = chr_seq[i-1][chr_seq[i-1].length()-1];
				prev_2nt += chr_seq[i][0];
			}
			else
			{
				prev_2nt = chr_seq[i-1].substr(chr_seq[i-1].length()-2, 2);
			}
		}
		else
		{
			if (j == 1)
			{
				prev_2nt = ' ';
				prev_2nt += chr_seq[i][0];
			}
			else
			{
				prev_2nt = "  ";
			}
		}

	while (left < len_need)
	{
		resultStr += chr_seq[i].substr(j, left);
		len_need -= left;
		i++;
		if (i >= chr_seq.size()) //error checking again
		{
			//cout << "GetSubstrFromVecStrs (2 extras) error: not enough (" << len << ") chars after start_pos (" << start_pos << ")" 
			//	<< i << ">=" << chr_seq.size() << "\n";
			error = true;
			//exit(-1);
			break;
		}
		j=0;
		left = chr_seq[i].length();
	}
	if (!error)
	{
		resultStr += chr_seq[i].substr(j, len_need);
		if (j+len_need+2 <= chr_seq[i].length())
			next_2nt = chr_seq[i].substr(j+len_need, 2);
		else
			if (j+len_need == chr_seq[i].length() - 1)
			{
				next_2nt = chr_seq[i][chr_seq[i].length()-1];
				if (i < chr_seq.size() - 1)
					next_2nt += chr_seq[i+1][0];
				else
					next_2nt += ' ';
			}
			else
			{
				if (i < chr_seq.size() - 1)
					if (chr_seq[i+1].length() >= 2)
					{
						next_2nt = chr_seq[i+1].substr(0, 2);
					}
					else
					{
						next_2nt = chr_seq[i+1][0];
						next_2nt += ' ';
					}
				else
					next_2nt = "  ";
			}
	}
	else //resultStr is at the very end of chromosome
	{
		next_2nt = "  ";
	}

	//for negative genes, complement the resultStr
	if (!isPosStrand)
	{
		StrToComplement(resultStr);
		StrToComplement(prev_2nt);
		StrToComplement(next_2nt);
		//StrToReverseComplement(resultStr); //Update!!!
	}
}


//no longer used
//if vecIndex==-1, it signals that it's the first time doing this, so getting the rest of string that contains "start_pos"
//otherwise getting the string at vecIndex
//-2 signals stop...
int GetSingleStrFromVecStrs(vector<string>& chr_seq, bool isPosStrand, const int vecIndex, const int start_pos, string& resultStr, int& leftover) //start_pos is index (starting from 0)

{

	int i, len;

	if (vecIndex == -1)

	{

	resultStr = "";

	if (start_pos <0 || chr_seq.empty())

	{

		cout << "error10: start_pos is negative (" << start_pos << ") or chr_seq is empty" << "\n";

		return -2;

	}



	int size = chr_seq[0].size(); //all chr_seq[i] is of the same size, except the last string

	i = start_pos / size;

	int j = start_pos % size;

	if (i < chr_seq.size() && j < chr_seq[i].length())

	{

		if (isPosStrand)

		{

			len = chr_seq[i].length() - j;

			leftover = len % 3;

			resultStr = chr_seq[i].substr(j, len-leftover);

		}

		else

		{

			leftover = (j+1) % 3;

			resultStr = chr_seq[i].substr(leftover, j+1-leftover);



			//for negative genes, complement the resultStr

			StrToComplement(resultStr);

		}



		return i;

	}

	else

	{

		cout << "error11: start_pos (" << start_pos << ") is bigger than sequence length" << "\n";

		return -2;

	}



	}

	else

	{

		if (isPosStrand)

		{

			i = vecIndex+1;

			if (i<chr_seq.size())

			{

				resultStr = chr_seq[i-1].substr(chr_seq[i-1].length()-leftover);

				len = chr_seq[i].length()+leftover;

				leftover = len % 3; //new leftover

				resultStr += chr_seq[i].substr(0, chr_seq[i].length() - leftover);

	

				return i;

			}

			else

			{

				return -2;

			}

		}

		else //negative strand

		{

			i = vecIndex-1;

			if (i>=0)

			{

				int last_leftover = leftover;

				len = leftover + chr_seq[i].length();

				leftover = len % 3; //new leftover

				resultStr = chr_seq[i].substr(leftover);

				resultStr += chr_seq[i+1].substr(0, last_leftover);



				//for negative genes, complement the resultStr

				StrToComplement(resultStr);



				return i;

			}

			else

			{

				return -2;

			}





		}

	}



}



void GetSubstrFromVecStrs_NegRev(vector<string>& chr_seq, bool isPosStrand, const int start_pos, const int len, string& resultStr) //start_pos is index (starting from 0)

{

//	cout << "NegRev: getting from " << start_pos << " for " << len << "(" << isPosStrand << ")" << "\n";



	resultStr = "";

	if (len <= 0)

		return;



	//vector<string>::iterator it = chr_seq.begin();

	int count = 0;

	int cur_start_pos = start_pos;
	cur_start_pos -= chromosome_start_pos - 1;

	if (cur_start_pos <0)

	{

		cout << "error12: start_pos is negative (" << start_pos << ")" << cur_start_pos << "\n";

		return;

	}

	int size;

	int len_need = len;



	size = chr_seq[0].length();

	int i = cur_start_pos / size;

	int j = cur_start_pos % size;

	if (i >= chr_seq.size()) //error checking

	{

		cout << "error13: start_pos (" << start_pos << ") is bigger than sequence length)" << cur_start_pos << "\n";

		return;

	}

	int left = chr_seq[i].length() - j;

	bool error = false;

	while (left < len_need)

	{

		resultStr += chr_seq[i].substr(j, left);

		len_need -= left;

		i++;

		if (i >= chr_seq.size()) //error checking again

		{

			cout << "GetSubstrFromVecStrs_NegRev error: not enough (" << len << ") chars after start_pos (" << start_pos << ")" 
				<< cur_start_pos << ";" 
				<< i << ">=" << chr_seq.size() << "\n";

			exit(-1);

			error = true;

			break;

		}

		j=0;

		left = chr_seq[i].length();

	}

	if (!error)

		resultStr += chr_seq[i].substr(j, len_need);



/*	while (it != chr_seq.end())

	{

		size = (*it).size();

		count += size;

		if (count <= cur_start_pos)

		{			

			it++;

			continue;

		}



		int left = count - cur_start_pos ;

		if (left >= len_need)

		{

			resultStr += (*it).substr(size-left, len_need);

			break;

		}

		else

		{

			resultStr += (*it).substr(size-left, left);

			cur_start_pos += left;

			len_need -= left;

			it++;

		}

	}

	if (count <= start_pos)

	{

		cout << "error: start_pos (" << start_pos << ") is bigger than sequence length (" << count << ")" << "\n";

		return;

	}

*/

	//for negative genes, complement the resultStr

	if (!isPosStrand)

		StrToReverseComplement(resultStr); //Update!!!

}

//this function allows to get string that is not up to "len" long
int GetSubstrFromVecStrs_NegRev_ForRepair(vector<string>& chr_seq, bool isPosStrand, const int start_pos, const int len, string& resultStr) //start_pos is index (starting from 0)
{
//	cout << "NegRev: getting from " << start_pos << " for " << len << "(" << isPosStrand << ")" << "\n";

	resultStr = "";
	if (len <= 0)
		return 0; //0 means there's no valid string returned
	int cur_start_pos = start_pos;
	cur_start_pos -= chromosome_start_pos - 1;

	int len_need = len;
	if (cur_start_pos <0)
	{
		//cout << "error7: start_pos is negative (" << start_pos << ")" << "\n";
		//return 0;
		len_need += cur_start_pos;
		cur_start_pos = 0;
	}

	//vector<string>::iterator it = chr_seq.begin();
	int count = 0;
	int size = chr_seq[0].length();
	int i = cur_start_pos / size;
	int j = cur_start_pos % size;
	if (i >= chr_seq.size()) //error checking
	{
		//cout << "error14: start_pos (" << start_pos << ") is bigger than sequence length" << "\n";
		return 0;
	}
	int left = chr_seq[i].length() - j;
	bool error = false;
	while (left < len_need)
	{
		resultStr += chr_seq[i].substr(j, left);
		len_need -= left;
		i++;
		if (i >= chr_seq.size()) //error checking again
		{
			//cout << "error: not enough (" << len << ") chars after start_pos (" << start_pos << ")" 
			//	<< i << ">=" << chr_seq.size() << "\n";
			//exit(-1);
			error = true;
			break;
		}
		j=0;
		left = chr_seq[i].length();
	}
	if (!error)
		resultStr += chr_seq[i].substr(j, len_need);

	//for negative genes, complement the resultStr
	if (!isPosStrand)
		StrToReverseComplement(resultStr); //Update!!!

	return resultStr.length();
}



void StrToComplement(string& resultStr)

{

	set<char> non_recog_chars;

	set<char>::iterator setIt;



	non_recog_chars.clear();



	int i;//, j;

	for (i=0; i<resultStr.length(); i++)

	{

		//j = resultStr[i];



		switch (tolower(resultStr[i]))

		{

		case 'a':

			resultStr[i] = 't';

			break;

		case 't':

			resultStr[i] = 'a';

			break;

		case 'g':

			resultStr[i] = 'c';

			break;

		case 'c':

			resultStr[i] = 'g';

			break;

		default:

			setIt = non_recog_chars.find(resultStr[i]);

			if (setIt == non_recog_chars.end())

			{

				//cout << "warning: sequence contains non-recognizable char " << resultStr[i] << " - code:" << j << "\n";

				non_recog_chars.insert(resultStr[i]);

			}

			break;

		}

	}



}



void StrToReverseComplement(string& resultStr)

{

	int i, k, m;//, j;



	k= resultStr.length();

	m = k*2;

	resultStr.resize(m); //allocate empty space



	set<char> non_recog_chars;

	set<char>::iterator setIt;



	non_recog_chars.clear();



	for (i=0; i<k; i++)

	{

		//j = resultStr[i];



		switch (tolower(resultStr[i]))

		{

		case 'a':

			resultStr[m-1-i] = 't';

			break;

		case 't':

			resultStr[m-1-i] = 'a';

			break;

		case 'g':

			resultStr[m-1-i] = 'c';

			break;

		case 'c':

			resultStr[m-1-i] = 'g';

			break;

		default:

			setIt = non_recog_chars.find(resultStr[i]);

			if (setIt == non_recog_chars.end())

			{

				//cout << "warning: sequence contains non-recognizable char " << resultStr[i] << " - code:" << j << "\n";

				non_recog_chars.insert(resultStr[i]);

			}

			resultStr[m-1-i] = resultStr[i];

			break;

		}

	}



	resultStr.erase(0, k);



}



void BackwardSearch(string& segStr, int start_pos, bool isPos, const string* matchStrs, int num_matchStrs, vector<int>& results, bool is_donor, ofstream& outFile)

//void BackwardSearch(string& segStr, int start_pos, bool isPos, const string* matchStrs, int num_matchStrs, vector<int>& results, bool is_donor)

{

	int i;

	int pos=-1; //string::npos;

	int match_pos=-1;//0;

	int count=0;

	bool frame[3] = {false, false, false};



	int tmpResults[3];

	while (count < 3)

	{

		int cur_match_pos;

		for (i=0; i<num_matchStrs; i++)

		{

			if ((cur_match_pos = segStr.rfind(matchStrs[i], pos)) != string::npos)

				if (match_pos < cur_match_pos)

					match_pos = cur_match_pos; //match_pos keeps the last matching position

		}



		if (match_pos == -1)//0)

			break;



		//compute frame

		int cur_frame = match_pos % 3;

		if (!frame[cur_frame])

		{

			frame[cur_frame] = true;

			if (is_donor)  //donors

				if (isPos) //positive strand

				{

					//outFile << "pos_donor: " << start_pos+match_pos-1 << "\n";

					//results.push_back(start_pos+match_pos-1);

					tmpResults[count]=start_pos+match_pos-1;

				}

				else //negative strand

				{

					//outFile << "neg_donor: " << -(start_pos+match_pos+2) << "\n";

					results.push_back(-(start_pos+match_pos+2));

					//tmpResults[count] = -(start_pos+match_pos+2);

				}

			else//acceptors

				if (isPos)

				{

					//outFile << "pos_acceptor: " << start_pos+match_pos+2 << "\n";

					//results.push_back(start_pos+match_pos+2); 

					tmpResults[count] = start_pos+match_pos+2;

				}

				else

				{

					//outFile << "neg_acceptor: " << -(start_pos+match_pos-1) << "\n";

					results.push_back(-(start_pos+match_pos-1));

					//tmpResults[count] = -(start_pos+match_pos-1);

				}

		

			count++;

		}



		if (match_pos == 0)

			break;

		pos = match_pos - 1;

		match_pos = -1;//0;

	}



	if (count) //now store tmpResults to results, in correct order!

	{

		if (isPos)

		{

			while (count>0)

			{

				results.push_back(tmpResults[count-1]);

				count--;

			}

		}

	}

}



//n is "MAX_NUM_SPLICE_SITES"

void BackwardSearch_nPerBorder(string& segStr, int start_pos, bool isPos, const string* matchStrs, int num_matchStrs, 

							   vector<int>& results, bool is_donor, vector<string>& site_2nt, ofstream& outFile, 

							   string& prev_2nt, string& next_2nt)

{

	int i, strLen = segStr.size();

	int pos;

	//int match_pos=-1;//0;

	int count=0;

	//bool frame[3] = {false, false, false};



	//vector<int> tmpResults;

	vector< pair<int, string> > tmpResults;



	//set<int> tmpSites;

	map<int, string> tmpSites;

	int cur_match_pos;

	for (i=0; i<num_matchStrs; i++)

	{

		pos = segStr.length()-1; //string::npos;

		int local_count = 0;

		while (local_count < MAX_NUM_SPLICE_SITES && pos >= 0) //find MAX_NUM_SPLICE_SITES for each matchStr

		{

			if ((cur_match_pos = segStr.rfind(matchStrs[i], pos)) != string::npos)

			{

				//tmpSites.insert(cur_match_pos);

				string tmpStr = "??";

				if (is_donor)

				{

					if (isPos)

					{

						//tmpStr = segStr[cur_match_pos-2];

						//tmpStr += segStr[cur_match_pos-1];						

						if (cur_match_pos > 1)

							tmpStr = segStr.substr(cur_match_pos-2, 2);

						else

							if (cur_match_pos == 1)

							{

								tmpStr = prev_2nt[1];

								tmpStr += segStr[0];

							}

							else

							{

								tmpStr = prev_2nt;

							}

					}

					else

					{

						if (cur_match_pos < strLen - 3)

						{

							tmpStr = segStr[cur_match_pos+3];

							tmpStr += segStr[cur_match_pos+2];

						}

						else

						{

							if (cur_match_pos == strLen - 3)

							{

								tmpStr = next_2nt[0];

								tmpStr += segStr[strLen-1];

							}

							else

							{

								tmpStr = next_2nt[1];

								tmpStr += next_2nt[0];

							}

						}

					}

				}

				else

				{

					if (isPos)

					{

						if (cur_match_pos < strLen - 3)

							tmpStr = segStr.substr(cur_match_pos+2, 2);

						else

							if (cur_match_pos == strLen - 3)

							{

								tmpStr = segStr[strLen-1];

								tmpStr += next_2nt[0];

							}

							else

							{

								tmpStr = next_2nt;

							}



					}

					else

					{

						if (cur_match_pos > 1)

						{

							tmpStr = segStr[cur_match_pos-1];

							tmpStr += segStr[cur_match_pos-2];

						}

						else

						{

							if (cur_match_pos == 1)

							{

								tmpStr = segStr[0];

								tmpStr += prev_2nt[1];

							}

							else

							{

								tmpStr = prev_2nt[1];

								tmpStr += prev_2nt[0];

							}

						}

					}

				}

				tmpSites.insert(map<int, string>::value_type(cur_match_pos, tmpStr));



				local_count++;

				pos = cur_match_pos - 1;

			}

			else

				break;

		}

	}



	if (!tmpSites.empty())

	{

		//set<int>::reverse_iterator site_pos = tmpSites.rbegin();

		map<int, string>::reverse_iterator site_pos = tmpSites.rbegin();

		while (count < MAX_NUM_SPLICE_SITES && site_pos != tmpSites.rend())

		{

			if (is_donor)

				if (isPos)

				{

					tmpResults.push_back(pair<int, string>(start_pos + (*site_pos).first - 1, (*site_pos).second));

				}

				else

				{

					results.push_back(-(start_pos + (*site_pos).first + 2));

					site_2nt.push_back((*site_pos).second);

				}

			else

				if (isPos)

				{

					tmpResults.push_back(pair<int, string>(start_pos+(*site_pos).first+2, (*site_pos).second));

				}

				else

				{

					results.push_back(-(start_pos+(*site_pos).first-1));

					site_2nt.push_back((*site_pos).second);

				}



			site_pos++;

			count++;

		}



		if (count) //now store tmpResults to results, in correct order!

		{

			if (isPos)

			{

				while (count>0)

				{

					results.push_back(tmpResults[count-1].first);

					site_2nt.push_back(tmpResults[count-1].second);

					count--;

				}

			}

		}

	}

}



void ForwardSearch(string& segStr, int start_pos, bool isPos, const string* matchStrs, int num_matchStrs, vector<int>& results,

				   bool is_donor, ofstream& outFile)

//void ForwardSearch(string& segStr, int start_pos, bool isPos, const string* matchStrs, int num_matchStrs, vector<int>& results, bool is_donor)

{

	int i, strLen = segStr.size();

	int pos=0;

	int match_pos=strLen; //string::npos; (<== this turns out to be -1!)

	int count=0;

	bool frame[3] = {false, false, false};



	int tmpResults[3];

	while (count < 3)

	{

		int cur_match_pos;

		for (i=0; i<num_matchStrs; i++)

		{

			if ((cur_match_pos = segStr.find(matchStrs[i], pos)) != string::npos)

				if (match_pos > cur_match_pos)

					match_pos = cur_match_pos; //match_pos keeps the first matching position

		}



		if (match_pos == strLen) //string::npos)

			break;



		//compute frame

		int cur_frame = match_pos % 3;

		if (!frame[cur_frame])

		{

			frame[cur_frame] = true;

			if (is_donor)  //donors

				if (isPos) //positive strand

				{

					//outFile << "pos_donor: " << start_pos+match_pos-1 << "\n";

					results.push_back(start_pos+match_pos-1);

				}

				else //negative strand

				{

					//outFile << "neg_donor: " << -(start_pos+match_pos+2) << "\n";

					//results.push_back(-(start_pos+match_pos+2));

					tmpResults[count] = -(start_pos+match_pos+2);

				}

			else//acceptors

				if (isPos)

				{

					//outFile << "pos_acceptor: " << start_pos+match_pos+2 << "\n";

					results.push_back(start_pos+match_pos+2); 

				}

				else

				{

					//outFile << "neg_acceptor: " << -(start_pos+match_pos-1) << "\n";

					//results.push_back(-(start_pos+match_pos-1));

					tmpResults[count] = -(start_pos+match_pos-1);

				}



			count++;

		}



		pos = match_pos + 1;

		match_pos = strLen; //string::npos;

	}



	if (count)

		if (!isPos)

			while (count>0)

			{

				results.push_back(tmpResults[count-1]);

				count--;

			}



}



void ForwardSearch_nPerBorder(string& segStr, int start_pos, bool isPos, const string* matchStrs, int num_matchStrs, vector<int>& results,

				   bool is_donor, vector<string>& site_2nt, ofstream& outFile, string& prev_2nt, string& next_2nt)

{

	int i, strLen = segStr.size();

	int pos;

	//int match_pos=strLen; //string::npos; (<== this turns out to be -1!)

	int count=0;

	//bool frame[3] = {false, false, false};



	//vector<int> tmpResults;

	vector< pair<int, string> > tmpResults;



	//set<int> tmpSites;

	map<int, string> tmpSites;

	int cur_match_pos;

	for (i=0; i<num_matchStrs; i++)

	{

		pos = 0;

		int local_count = 0;

		while (local_count < MAX_NUM_SPLICE_SITES && pos < strLen) //find MAX_NUM_SPLICE_SITES for each matchStr

		{

			if ((cur_match_pos = segStr.find(matchStrs[i], pos)) != string::npos)

			{

				//tmpSites.insert(cur_match_pos);

				string tmpStr="??";

				if (is_donor)

				{

					if (isPos)

					{

						//tmpStr = segStr[cur_match_pos-2];

						//tmpStr += segStr[cur_match_pos-1];						

						if (cur_match_pos > 1)

							tmpStr = segStr.substr(cur_match_pos-2, 2);

						else

							if (cur_match_pos == 1)

							{

								tmpStr = prev_2nt[1];

								tmpStr += segStr[0];

							}

							else

							{

								tmpStr = prev_2nt;

							}

					}

					else

					{

						if (cur_match_pos < strLen - 3)

						{

							tmpStr = segStr[cur_match_pos+3];

							tmpStr += segStr[cur_match_pos+2];

						}

						else

						{

							if (cur_match_pos == strLen - 3)

							{

								tmpStr = next_2nt[0];

								tmpStr += segStr[strLen-1];

							}

							else

							{

								tmpStr = next_2nt[1];

								tmpStr += next_2nt[0];

							}

						}

					}

				}

				else //acceptor

				{

					if (isPos)

					{

						if (cur_match_pos < strLen - 3)

							tmpStr = segStr.substr(cur_match_pos+2, 2);

						else

							if (cur_match_pos == strLen - 3)

							{

								tmpStr = segStr[strLen-1];

								tmpStr += next_2nt[0];

							}

							else

							{

								tmpStr = next_2nt;

							}



					}

					else

					{

						if (cur_match_pos > 1)

						{

							tmpStr = segStr[cur_match_pos-1];

							tmpStr += segStr[cur_match_pos-2];

						}

						else

						{

							if (cur_match_pos == 1)

							{

								tmpStr = segStr[0];

								tmpStr += prev_2nt[1];

							}

							else

							{

								tmpStr = prev_2nt[1];

								tmpStr += prev_2nt[0];

							}

						}

					}

				}

				tmpSites.insert(map<int, string>::value_type(cur_match_pos, tmpStr));



				local_count++;

				pos = cur_match_pos + 1;

			}

			else

				break;

		}

	}



	if (!tmpSites.empty())

	{

		//set<int>::iterator site_pos = tmpSites.begin();

		map<int, string>::iterator site_pos = tmpSites.begin();

		while (count < MAX_NUM_SPLICE_SITES && site_pos != tmpSites.end())

		{

			if (is_donor)

				if (isPos)

				{

					//results.push_back(start_pos+(*site_pos)-1);

					results.push_back(start_pos+(*site_pos).first-1);

					site_2nt.push_back((*site_pos).second);

				}

				else

				{

					tmpResults.push_back(pair<int, string>(-(start_pos + (*site_pos).first + 2), (*site_pos).second));

				}

			else

				if (isPos)

				{

					results.push_back(start_pos+(*site_pos).first+2);

					site_2nt.push_back((*site_pos).second);

				}

				else

				{

					tmpResults.push_back(pair<int, string>(-(start_pos+(*site_pos).first-1), (*site_pos).second));

				}

			site_pos++;

			count++;

		}



		if (count) //now store tmpResults to results, in correct order!

		{

			if (!isPos)

				while (count>0)

				{

					results.push_back(tmpResults[count-1].first);

					site_2nt.push_back(tmpResults[count-1].second);

					count--;

				}

		}

	}

}





void LoadCodonMap()

{



	DNA_CODON_TBL.insert(map<string, string>::value_type("att", "Ile"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("atc", "Ile"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("ata", "Ile"));



	DNA_CODON_TBL.insert(map<string, string>::value_type("ctt", "Leu"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("ctc", "Leu"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("cta", "Leu"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("ctg", "Leu"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("tta", "Leu"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("ttg", "Leu"));



	DNA_CODON_TBL.insert(map<string, string>::value_type("gtt", "Val"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("gtc", "Val"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("gta", "Val"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("gtg", "Val"));



	DNA_CODON_TBL.insert(map<string, string>::value_type("ttt", "Phe"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("ttc", "Phe"));



	DNA_CODON_TBL.insert(map<string, string>::value_type("atg", "Met"));



	DNA_CODON_TBL.insert(map<string, string>::value_type("tgt", "Cys"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("tgc", "Cys"));



	DNA_CODON_TBL.insert(map<string, string>::value_type("gct", "Ala"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("gcc", "Ala"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("gca", "Ala"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("gcg", "Ala"));





	DNA_CODON_TBL.insert(map<string, string>::value_type("ggt", "Gly"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("ggc", "Gly"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("gga", "Gly"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("ggg", "Gly"));



	DNA_CODON_TBL.insert(map<string, string>::value_type("cct", "Pro"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("ccc", "Pro"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("cca", "Pro"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("ccg", "Pro"));





	DNA_CODON_TBL.insert(map<string, string>::value_type("act", "Thr"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("acc", "Thr"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("aca", "Thr"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("acg", "Thr"));



	DNA_CODON_TBL.insert(map<string, string>::value_type("tct", "Ser"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("tcc", "Ser"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("tca", "Ser"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("tcg", "Ser"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("agt", "Ser"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("agc", "Ser"));



	DNA_CODON_TBL.insert(map<string, string>::value_type("tat", "Tyr"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("tac", "Tyr"));



	DNA_CODON_TBL.insert(map<string, string>::value_type("tgg", "Trp"));



	DNA_CODON_TBL.insert(map<string, string>::value_type("caa", "Gln"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("cag", "Gln"));



	DNA_CODON_TBL.insert(map<string, string>::value_type("aat", "Asn"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("aac", "Asn"));



	DNA_CODON_TBL.insert(map<string, string>::value_type("cat", "His"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("cac", "His"));



	DNA_CODON_TBL.insert(map<string, string>::value_type("gaa", "Glu"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("gag", "Glu"));





	DNA_CODON_TBL.insert(map<string, string>::value_type("gat", "Asp"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("gac", "Asp"));



	DNA_CODON_TBL.insert(map<string, string>::value_type("aaa", "Lys"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("aag", "Lys"));



	DNA_CODON_TBL.insert(map<string, string>::value_type("cgt", "Arg"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("cgc", "Arg"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("cga", "Arg"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("cgg", "Arg"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("aga", "Arg"));

	DNA_CODON_TBL.insert(map<string, string>::value_type("agg", "Arg"));







}



void LoadCodonMap_SingleLetterAA()

{

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("tag", '*')); //stop codon

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("tga", '*')); //stop codon

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("taa", '*')); //stop codon



	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("att", 'I'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("atc", 'I'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("ata", 'I'));



	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("ctt", 'L'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("ctc", 'L'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("cta", 'L'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("ctg", 'L'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("tta", 'L'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("ttg", 'L'));



	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("gtt", 'V'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("gtc", 'V'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("gta", 'V'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("gtg", 'V'));



	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("ttt", 'F'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("ttc", 'F'));



	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("atg", 'M'));



	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("tgt", 'C'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("tgc", 'C'));



	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("gct", 'A'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("gcc", 'A'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("gca", 'A'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("gcg", 'A'));





	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("ggt", 'G'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("ggc", 'G'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("gga", 'G'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("ggg", 'G'));



	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("cct", 'P'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("ccc", 'P'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("cca", 'P'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("ccg", 'P'));





	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("act", 'T'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("acc", 'T'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("aca", 'T'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("acg", 'T'));



	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("tct", 'S'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("tcc", 'S'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("tca", 'S'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("tcg", 'S'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("agt", 'S'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("agc", 'S'));



	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("tat", 'Y'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("tac", 'Y'));



	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("tgg", 'W'));



	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("caa", 'Q'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("cag", 'Q'));



	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("aat", 'N'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("aac", 'N'));



	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("cat", 'H'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("cac", 'H'));



	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("gaa", 'E'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("gag", 'E'));





	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("gat", 'D'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("gac", 'D'));



	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("aaa", 'K'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("aag", 'K'));



	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("cgt", 'R'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("cgc", 'R'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("cga", 'R'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("cgg", 'R'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("aga", 'R'));

	DNA_CODON_TBL_SL.insert(map<string, char>::value_type("agg", 'R'));







}



void TranslateDNACodons(string& dna, string& protein)

{

	protein = "";

	int i, j;

	map<string, string>::iterator it;

	

	int dna_len = dna.length();

	for (i=0; i<dna_len; i+=3)

	{

		string codon = dna.substr(i, 3);

		StrToLower(codon);

		it = DNA_CODON_TBL.find(codon);



		int codon_len = codon.length();

		if (it == DNA_CODON_TBL.end())

		{

			for (j=0; j<codon_len; j++)

				protein += ' '; //up to 3 spaces

		}

		else

		{

			protein += (*it).second;

		}

	}

}





//find the real position (index) of [search_start, search_end] in the area [align_start, align_end] on the alignStr (targetStr)

//search_start/end, align_start/end are ABSOLUTE positions, and this returns the RELATIVE index of search_start/end on alignStr

bool FindRealPos(string& alignStr, int align_start, int align_end, int search_start, int search_end,

						int& search_start_pos, int& search_end_pos, int& extra_front, int& extra_end)

//modified: each char in "alignStr" now represent 3 chars on chromosome sequence (FindRealPos is only used by ComputePenalty)

{

/*	if (search_start < align_start || search_end > align_end)

	{

		search_start_pos = 0;

		search_end_pos = -1; //indicate no search area

		return;

	}







	int count = align_start;

	int cur_pos = 0;

		

	while (count <= search_start)

	//while (count < search_start)

	{

		if (alignStr[cur_pos] != '-')

			//count++;

			count += 3;

		cur_pos++;

	}

	search_start_pos = cur_pos-1;

	//search_start_pos = cur_pos;



	while (count <= search_end)

	//while (count < search_end)

	{

		if (alignStr[cur_pos] != '-')

			//count++;

			count += 3;

		cur_pos++;

	}

	search_end_pos = cur_pos-1;

	//search_end_pos = cur_pos;

*/

	if (alignStr.empty() || search_start < align_start || search_end > align_end) //no str

	{

		//indicate no search area

		return false;

	}



	int count = align_start;

	int cur_pos = 0;

		

	while (count < search_start)

	{

		if (alignStr[cur_pos] != '-')

			//count++;

			count += 3;

		cur_pos++;

	}

	search_start_pos = cur_pos;

	extra_front = count - search_start;



	while (count < search_end)

	{

		if (alignStr[cur_pos] != '-')

			//count++;

			count += 3;

		cur_pos++;

	}

	search_end_pos = cur_pos-1;

	switch (count - search_end)

	{

	case 0:

		extra_end = 1;

		break;

	case 1:

		extra_end = 0;

		break;

	case 2:

		search_end_pos--;

		extra_end = 2;

		break;

	default:

		break;

	}



	return true;



}



//exact_match, gap, pos_match must have already been initialized to 0 before calling this function

void GetPenaltyFromMatchStr(string& alignStr, int align_start, int align_end, int search_start, int search_end,

						string& matchStr, int& exact_match, int& gap, int& pos_match, 

						int query_start, string& queryStr, set< pair<int,int> >& query_match_pos, bool& end_is_match)

{

	int search_start_pos, search_end_pos;

	int extra_front, extra_end;



/*	if (search_start < align_start || search_end > align_end) //no str

	{

		exact_match = 0;

		gap = 0;

		pos_match = 0;

		return;

	}



	int count = align_start;

	int cur_pos = 0;

		

	while (count < search_start)

	{

		if (alignStr[cur_pos] != '-')

			//count++;

			count += 3;

		cur_pos++;

	}

	search_start_pos = cur_pos;

	extra_front = count - search_start;



	while (count < search_end)

	{

		if (alignStr[cur_pos] != '-')

			//count++;

			count += 3;

		cur_pos++;

	}

	search_end_pos = cur_pos-1;

	switch (count - search_end)

	{

	case 0:

		extra_end = 1;

		break;

	case 1:

		extra_end = 0;

		break;

	case 2:

		search_end_pos--;

		extra_end = 2;

		break;

	default:

		break;

	}

*/

	if (!FindRealPos(alignStr, align_start, align_end, search_start, search_end, 

		search_start_pos, search_end_pos, extra_front, extra_end))

	{

/*		exact_match = 0;

		gap = 0;

		pos_match = 0;

*/		return;

	}



	//added, record query positions

	int query_count=query_start;

	int i;

	for (i=0; i<search_start_pos; i++)

		if (queryStr[i] <= 'Z' && queryStr[i] >= 'A')

			query_count ++;



	//string read_str = matchStr.substr(search_start_pos, search_end_pos - search_start_pos + 1);

	//gap = CountCharInStr(read_str, ' ')*3;

	//pos_match = CountCharInStr(read_str, '+')*3;

	//exact_match = read_str.length()*3 - gap - pos_match;

	bool matched = false;

	int match_start;

	for (i=search_start_pos; i<=search_end_pos; i++)

	{

		switch (matchStr[i])

		{

		case ' ':

			gap += 3;

			if (matched)

			{

				matched = false;

				query_match_pos.insert( pair<int,int>(match_start, query_count));

			}

			if (queryStr[i] <= 'Z' && queryStr[i] >= 'A')

				query_count++;

			break;

		case '+':

			pos_match += 3;

			if (matched)

			{

				matched = false;

				query_match_pos.insert( pair<int,int>(match_start, query_count));

			}

			if (queryStr[i] <= 'Z' && queryStr[i] >= 'A')

				query_count++;

			break;

		default: //queryStr[i] here must be a letter (>='A' and <='Z')

			exact_match += 3;			

			if (!matched)

			{

				matched = true;				

				match_start = query_count;

			}

			query_count++;

			break;

		}

	}

	if (matched)

		query_match_pos.insert( pair<int,int>(match_start, query_count) );

	if (extra_end == 0)

	{

		if (i == search_end_pos + 1)

		{

			if (matchStr[i-1] == ' ' || matchStr[i-1] == '+')

				end_is_match = false;

			else

				end_is_match = true;

		}

		else

			end_is_match = false; //impossible case?

	}





	char ch;

	if (extra_front > 0)

	{

		ch = matchStr[search_start_pos-1];

		if (ch == ' ')

			gap += extra_front;

		else

			if (ch == '+')

				pos_match += extra_front;

			else

				exact_match += extra_front;

	}

	if (extra_end > 0)

	{

		ch = matchStr[search_end_pos+1];

		if (ch == ' ')

		{

			gap += extra_end;

			end_is_match = false;

		}

		else

			if (ch == '+')

			{

				pos_match += extra_end;

				end_is_match = false;

			}

			else

			{

				exact_match += extra_end;

				end_is_match = true;

			}

	}



}



void GetPenaltyFromMatchStr(int search_start_pos, int search_end_pos, int extra_front, int extra_end,

							string& matchStr, int& exact_match, int& gap, int& pos_match)

{

		

	int i;



	//string read_str = matchStr.substr(search_start_pos, search_end_pos - search_start_pos + 1);

	//gap = CountCharInStr(read_str, ' ')*3;

	//pos_match = CountCharInStr(read_str, '+')*3;

	//exact_match = read_str.length()*3 - gap - pos_match;

	for (i=search_start_pos; i<=search_end_pos; i++)

	{

		switch (matchStr[i])

		{

		case ' ':

			gap += 3;

			break;

		case '+':

			pos_match += 3;

			break;

		default: //queryStr[i] here must be a letter (>='A' and <='Z')

			exact_match += 3;			

			break;

		}

	}



	char ch;

	if (extra_front > 0)

	{

		ch = matchStr[search_start_pos-1];

		if (ch == ' ')

			gap += extra_front;

		else

			if (ch == '+')

				pos_match += extra_front;

			else

				exact_match += extra_front;

	}

	if (extra_end > 0)

	{

		ch = matchStr[search_end_pos+1];

		if (ch == ' ')

			gap += extra_end;

		else

			if (ch == '+')

				pos_match += extra_end;

			else

				exact_match += extra_end;

	}



}





//roll back target string (used only for roll back overlapping HSP )

//search_start and search_end both < cur_pos

void FindTgtStrPos(string& targetStr, int cur_pos, int search_start, int search_end, int& start_index, int& end_index)

{

	int count = cur_pos;

	int cur_index = targetStr.length() - 1;



	while (count > search_end)

	{

		if (targetStr[cur_index] != '-')

			count--;

		cur_index--;

	}

	end_index = cur_index;



	while (count > search_start)

	{

		if (targetStr[cur_index] != '-')

			count--;

		cur_index--;

	}

	start_index = cur_index;

}



/*

//update ext, update strChrs

void FillFrontExtArea(string& matchStr, int exon_index, int& ext, vector<string>& chr_seq, int exon_pos, 

				 vector<string>&  strChrs, string& queryStr, string& targetStr)

{

	if (exon_index < 0)

	{

		FillFrontFixedLenExt(chr_seq, ext, exon_pos, strChrs);

	}

	else

	{



		if (matchStr[exon_index] != ' ')

		{

			int m = matchStr.rfind(' ', exon_index);

			if (m != string::npos)

				m++;

			else

				m=0;

			

			if (m < exon_index)

			{

				ext = exon_index-m;



				int i;

				for (i=m; i<exon_index; i++)

				{

					strChrs[0] += queryStr[i];

					strChrs[1] += matchStr[i];

					strChrs[2] += ' ';

					strChrs[3] += toupper(targetStr[i]);

				}

			}

			else //if HSP happens to start from exon_index

			{

				FillFrontFixedLenExt(chr_seq, ext, exon_pos, strChrs);

			}



		}

		else //otherwise display 3 base pairs intron

		{

			FillFrontFixedLenExt(chr_seq, ext, exon_pos, strChrs);

		}

	}



}







void FillEndExtArea(string& matchStr, int exon_index, int& ext, vector<string>& chr_seq, int exon_pos,

						   vector<string>&  strChrs, string& queryStr, string& targetStr)

{

	int i, m;



	if (matchStr[exon_index] != ' ')

	{

		m = matchStr.find(' ', exon_index);

		if (m == string::npos)

			m=matchStr.length();

		

		if (m > exon_index+1)

		{

			ext = m-exon_index-1;



			for (i=exon_index+1; i<m; i++)

			{

				strChrs[0] += queryStr[i];

				strChrs[1] += matchStr[i];

				strChrs[2] += ' ';

				strChrs[3] += toupper(targetStr[i]);

			}

		}

		else

			FillEndFixedLenExt(chr_seq, ext, exon_pos, strChrs);

	}

	else //otherwise display 3 base pairs intron

	{

		FillEndFixedLenExt(chr_seq, ext, exon_pos, strChrs);

	}







}



void FillFrontFixedLenExt(vector<string>& chr_seq, int& ext, int exon_pos, vector<string>& strChrs)

{

	ext = 3;

	if (exon_pos > 0) //pos strand

		FillFromTargetStr(chr_seq, true, exon_pos-4, ext, false, strChrs);

	else

		FillFromTargetStr(chr_seq, false, -exon_pos, ext, true, strChrs);

}



void FillEndFixedLenExt(vector<string>& chr_seq, int& ext, int exon_pos, vector<string>& strChrs)

{

	ext = 3;



	if (exon_pos > 0) //pos strand

		FillFromTargetStr(chr_seq, true, exon_pos, ext, false, strChrs);

	else

		FillFromTargetStr(chr_seq, false, -exon_pos-4, ext, true, strChrs);

}





void FillFromTargetStr(vector<string>& chr_seq, bool isPosStrand, int start_pos, int length, bool revStr, vector<string>& strChrs)

{

	string tmpStr;

	GetSubstrFromVecStrs(chr_seq, isPosStrand, start_pos, length, tmpStr);

	StrToLower(tmpStr);



	if (revStr)

	{

		//reverse tmpStr

		char tmpChar = tmpStr[0];

		tmpStr[0] = tmpStr[2];

		tmpStr[2] = tmpChar;

	}



	int i;

	for (i=0; i<3; i++)

	{

		strChrs[0] += ' ';

		strChrs[1] += ' ';

		strChrs[2] += ' ';

		strChrs[3] += toupper(tmpStr[i]);

	}

}



void FillExonArea(vector<string>& strChrs, string& queryStr, string& matchStr, string& targetStr, vector<string>& chr_seq, 

				  bool isPosStrand, int exon_pos, int front_gap, int frame, int start_index, int end_index)

{

	int i;

	string tmpStr, revStr="", proStr, passStr;

	if (front_gap > 0)

	{

		

		if (exon_pos > 0)

		{

			GetSubstrFromVecStrs(chr_seq, isPosStrand, exon_pos-1, front_gap, tmpStr);

			revStr = tmpStr;

			for (i=0; i<front_gap; i++)

			{

				strChrs[0] += ' ';

				strChrs[1] += ' ';

			}

		}

		else

		{

			GetSubstrFromVecStrs(chr_seq, isPosStrand, -exon_pos-front_gap, front_gap, tmpStr);

			//reverse tmpStr

			for (i=front_gap-1; i>=0; i--)

			{

				strChrs[0] += ' ';

				strChrs[1] += ' ';

				revStr += tmpStr[i];

			}

		}

			

		passStr = revStr.substr(frame);

		TranslateDNACodons(passStr, proStr);

		for (i=0; i<frame; i++)

			strChrs[2] += ' ';

		strChrs[2] += proStr;

		strChrs[3] += revStr;

	}



	strChrs[0] += queryStr.substr(start_index, end_index-start_index+1);

	strChrs[1] += matchStr.substr(start_index, end_index-start_index+1);

	strChrs[3] += targetStr.substr(start_index, end_index-start_index+1);



	if (front_gap > 0) //start_index must be 0

	{

		proStr = "";

		if (end_index >= start_index)

		{

			passStr = targetStr.substr(start_index, end_index-start_index+1);

			TranslateDNACodons(passStr, proStr);

		}

	}

	else

	{

		passStr = targetStr.substr(start_index+frame, end_index-start_index-frame+1);

		TranslateDNACodons(passStr, proStr);

		for (i=0; i<frame; i++)

			strChrs[2] += ' ';

	}

	strChrs[2] += proStr;



}

*/



void StrToLower(string& str)

{
#if defined(linux) || defined(__linux)
	transform(str.begin(), str.end(), str.begin(), (int(*)(int))std::tolower); //GNU G++ need (int(*)(int))std::tolower to compile correctly
#elif defined(__APPLE__) && defined(__MACH__) /*mac os X*/
	std::transform(str.begin(), str.end(), str.begin(), (int(*)(int))std::tolower);
#else
	transform(str.begin(), str.end(), str.begin(), tolower); //for visual c++
#endif

/*	unsigned int i;

	for (i=0; i<str.length(); i++)

		str[i] = tolower(str[i]);

*/

}


//this function is only used to check length of specific type of vector<string>, chr_seq, 
//which store chromosome sequences (see GetChromosomes() function in data_manager), 
//where each string is the same length (MAX_LINE) except the very last string in the vector!
int LenOfStrVec(vector<string>& str_vec)
{

/*	int len = 0;
	vector<string>::iterator vecIt = str_vec.begin();
	for (; vecIt != str_vec.end(); vecIt++)
		len += (*vecIt).length();
	return len;
*/
	if (str_vec.empty())
		return 0;
	else
		return (str_vec.size()-1)*MAX_LINE + str_vec.back().size();
}



//VERY COSTLY function, do not use

void ProteinToDNAStyle(string& str)
{
	string tmpStr="";
	unsigned int i, j;
	for (i=0; i<str.length(); i++)
	{
		for (j=0; j<3; j++)
			tmpStr += str[i];
	}

	str = tmpStr;
}

int ProteinPosToDNAStylePos(int pos) //for query_start only; for query_end, must add additional 2
{
	return (pos - 1) * 3 + 1;
}

int DNAToProteinStylePos_StartSite(int pos)
{
	return (pos - 1) / 3 + 1;
}

int DNAToProteinStylePos_EndSite(int pos)
{
	return pos / 3;
}


//No longer used
void GetStopCodonIndexes(string& cur_chr_seq, vector<int>& stop_indexes, const string* stop_codons)

{

	int pos = 0;

	int i;



	for (i=0; i<3; i++) //3 stop_codons

	{

			while ((pos = cur_chr_seq.find(stop_codons[i], pos)) != string::npos)

			{

				stop_indexes.push_back(pos);

				pos++;

			}

	}



}



int FindHSPNo(vector<int>& sites, vector<int>& site_hsp_no, int site)

{

	int i;

	vector<int>::iterator it = site_hsp_no.begin();

	for (i=0; i<sites.size(); i++, it++)

		if (sites[i] == site)

			return *it;

	return -1; //failure? start_site/end_site?

}



int FindMax(vector<int>& hsp_no)

{

	if (hsp_no.empty())

		return -1;



	int max = hsp_no[0];

	int i;

	for (i=1; i< hsp_no.size(); i++)

		if (max < hsp_no[i])

			max = hsp_no[i];



	return max;

}



bool RegionOverlap(int start, int end, set< pair<int,int> >& regions)

{

	set< pair<int,int> >::iterator it = regions.begin();

	for (; it != regions.end(); it++)

	{

		if ((*it).second >= start )

		{

			if ((*it).first <= end)

				return true;

		}

	}



	return false;

}



//stupid stl... my simple implementation of "is_sorted"

bool IsSorted(vector<int>& vec)

{

	if (vec.size() < 2) //size is 1 or 0

		return true;



	int i;

	for (i=0; i<vec.size()-1; i++)

		if (vec[i] > vec[i+1])

			return false;



	return true;

}



//based on a sorted vector!!!

bool IsUnique(vector<int>& vec)

{

	if (vec.size() < 2) //size is 1 or 0

		return true;



	int i;

	for (i=0; i<vec.size()-1; i++)

		if (vec[i] == vec[i+1])

			return false;



	return true;

}



void BackwardSearch(string& segStr, int start_pos, bool isPos, const string* matchStrs, int num_matchStrs, vector<int>& results, 

					bool is_donor, ofstream& outFile, int seg_index, vector<int>& segs)

//void BackwardSearch(string& segStr, int start_pos, bool isPos, const string* matchStrs, int num_matchStrs, vector<int>& results, bool is_donor)

{

	int i;

	int pos=-1; //string::npos;

	int match_pos=-1;//0;

	int count=0;

	bool frame[3] = {false, false, false};



	int tmpResults[3];

	while (count < 3)

	{

		int cur_match_pos;

		for (i=0; i<num_matchStrs; i++)

		{

			if ((cur_match_pos = segStr.rfind(matchStrs[i], pos)) != string::npos)

				if (match_pos < cur_match_pos)

					match_pos = cur_match_pos; //match_pos keeps the last matching position

		}



		if (match_pos == -1)//0)

			break;



		//compute frame

		int cur_frame = match_pos % 3;

		if (!frame[cur_frame])

		{

			frame[cur_frame] = true;

			if (is_donor)  //donors

				if (isPos) //positive strand

				{

					//outFile << "pos_donor: " << start_pos+match_pos-1 << "\n";

					//results.push_back(start_pos+match_pos-1);

					tmpResults[count]=start_pos+match_pos-1;

				}

				else //negative strand

				{

					//outFile << "neg_donor: " << -(start_pos+match_pos+2) << "\n";

					results.push_back(-(start_pos+match_pos+2));

					segs.push_back(seg_index);

					//tmpResults[count] = -(start_pos+match_pos+2);

				}

			else//acceptors

				if (isPos)

				{

					//outFile << "pos_acceptor: " << start_pos+match_pos+2 << "\n";

					//results.push_back(start_pos+match_pos+2); 

					tmpResults[count] = start_pos+match_pos+2;

				}

				else

				{

					//outFile << "neg_acceptor: " << -(start_pos+match_pos-1) << "\n";

					results.push_back(-(start_pos+match_pos-1));

					segs.push_back(seg_index);

					//tmpResults[count] = -(start_pos+match_pos-1);

				}

		

			count++;

		}



		if (match_pos == 0)

			break;

		pos = match_pos - 1;

		match_pos = -1;//0;

	}



	if (count) //now store tmpResults to results, in correct order!

	{

		if (isPos)

		{

			while (count>0)

			{

				results.push_back(tmpResults[count-1]);

				segs.push_back(seg_index);

				count--;

			}

		}

	}

}



void ForwardSearch(string& segStr, int start_pos, bool isPos, const string* matchStrs, int num_matchStrs, vector<int>& results, 

				   bool is_donor, ofstream& outFile, int seg_index, vector<int>& segs)

//void ForwardSearch(string& segStr, int start_pos, bool isPos, const string* matchStrs, int num_matchStrs, vector<int>& results, bool is_donor)

{

	int i, strLen = segStr.size();

	int pos=0;

	int match_pos=strLen; //string::npos; (<== this turns out to be -1!)

	int count=0;

	bool frame[3] = {false, false, false};



	int tmpResults[3];

	while (count < 3)

	{

		int cur_match_pos;

		for (i=0; i<num_matchStrs; i++)

		{

			if ((cur_match_pos = segStr.find(matchStrs[i], pos)) != string::npos)

				if (match_pos > cur_match_pos)

					match_pos = cur_match_pos; //match_pos keeps the first matching position

		}



		if (match_pos == strLen) //string::npos)

			break;



		//compute frame

		int cur_frame = match_pos % 3;

		if (!frame[cur_frame])

		{

			frame[cur_frame] = true;

			if (is_donor)  //donors

				if (isPos) //positive strand

				{

					//outFile << "pos_donor: " << start_pos+match_pos-1 << "\n";

					results.push_back(start_pos+match_pos-1);

					segs.push_back(seg_index);

				}

				else //negative strand

				{

					//outFile << "neg_donor: " << -(start_pos+match_pos+2) << "\n";

					//results.push_back(-(start_pos+match_pos+2));

					tmpResults[count] = -(start_pos+match_pos+2);

				}

			else//acceptors

				if (isPos)

				{

					//outFile << "pos_acceptor: " << start_pos+match_pos+2 << "\n";

					results.push_back(start_pos+match_pos+2); 

					segs.push_back(seg_index);

				}

				else

				{

					//outFile << "neg_acceptor: " << -(start_pos+match_pos-1) << "\n";

					//results.push_back(-(start_pos+match_pos-1));

					tmpResults[count] = -(start_pos+match_pos-1);

				}



			count++;

		}



		pos = match_pos + 1;

		match_pos = strLen; //string::npos;

	}



	if (count)

		if (!isPos)

			while (count>0)

			{

				results.push_back(tmpResults[count-1]);

				segs.push_back(seg_index);

				count--;

			}



}



/*

void CombineInFrameSpliceSites(int border, vector<int>& sites, vector<int>& final_sites, vector<int>& site_regions, 

							   int region, bool is_acceptor, vector<int>& site_head_tail)

{

	//combine same frame acceptors, so that there's at most 1 acceptor for each frame around this segment border

	vector<int>::iterator site_i1, site_i2;



	sort(sites.begin(), sites.end());

	vector<int> real_sites;	

			

			while ( !sites.empty() )

			{

				site_i1=sites.begin();

				bool inframe = false;

				site_i2 = site_i1;

				site_i2++;

				for (; site_i2 != sites.end(); site_i2++)

				{



					if ( (*site_i1 - *site_i2) % 3 == 0) //same frame

					{

						int d1=abs(*site_i1-border);

						int d2=abs(*site_i2-border);

						if ( d1<d2 )

						{

							real_sites.push_back(*site_i1);

						}

						else

						{

							if (d1>d2)

								real_sites.push_back(*site_i2);

							else

							{

								real_sites.push_back(*site_i1);

								real_sites.push_back(*site_i2);

							}

						}

						sites.erase(site_i2);

						sites.erase(site_i1);

						inframe = true;

						break;

					}

				}

				if (!inframe)

				{

					real_sites.push_back(*site_i1);

					sites.erase(site_i1);

				}

			}

			sort(real_sites.begin(), real_sites.end());

			for (site_i1 = real_sites.begin(); site_i1 != real_sites.end(); site_i1++)

			{

				final_sites.push_back(*site_i1);

				site_regions.push_back(region);

				if (is_acceptor) //for acceptor, record head

				{

					if (*site_i1 <= border)

						site_head_tail.push_back((border - *site_i1)%3);

					else

						site_head_tail.push_back((3 - (*site_i1 - border) % 3) % 3);

				}

				else //for donor, record tail

				{

					if (*site_i1 <= border)

						site_head_tail.push_back((3 - (border - *site_i1)%3) % 3);

					else

						site_head_tail.push_back( (*site_i1 - border) % 3);

				}

			}





}

*/



void CombineSpliceSites(int border, vector<int>& sites, vector<int>& final_sites, vector<int>& site_regions, 

							   int region, bool is_acceptor, vector<int>& site_head_tail, int frame_reference, 

							   vector<string>& site_2nt, vector<string>& final_site_2nt, ofstream& os)

{

	multimap<int,int> tmpSitesMap;



	//vector<int>::iterator site_it;

	int site_i;

	//for (site_it = sites.begin(); site_it != sites.end(); site_it++)

	for (site_i = 0; site_i < sites.size(); site_i++)

		//tmpSitesMap.insert(multimap<int, int>::value_type(abs(*site_it - border), *site_it)); //sort sites by their distance from border

		tmpSitesMap.insert(multimap<int, int>::value_type(abs(sites[site_i] - border), site_i));



	//vector<int> tmpSites;

	map<int, string> tmpSites;

	int count = 0;

	multimap<int, int>::iterator map_it = tmpSitesMap.begin();

	while (count < MAX_NUM_SPLICE_SITES && map_it != tmpSitesMap.end())

	{

		//tmpSites.push_back((*map_it).second); //get MAX_NUM_SPLICE_SITES sites that are closest around border

		tmpSites.insert(map<int, string>::value_type(sites[(*map_it).second], site_2nt[(*map_it).second]));



		map_it++;

		count++;

	}



	//sort(tmpSites.begin(), tmpSites.end()); //sort so we can store to final_sites



	sites.clear(); //reset sites

	site_2nt.clear();



	//vector<int>::iterator tmp_site_it;

	map<int, string>::iterator tmp_site_it;

	int cur_site;

	for (tmp_site_it = tmpSites.begin(); tmp_site_it != tmpSites.end(); tmp_site_it++)

	{

		//cur_site = *tmp_site_it;

		cur_site = (*tmp_site_it).first;

		final_sites.push_back(cur_site);

		final_site_2nt.push_back((*tmp_site_it).second);

		site_regions.push_back(region);



#ifdef DEBUG

		os << "cur_site:" << cur_site << "; site_2nt:" << (*tmp_site_it).second << "; site_head_tail:";

#endif

		if (is_acceptor) //for acceptor, record head

		{

			//if (cur_site <= border)

			if (cur_site <= frame_reference)

				//site_head_tail.push_back((border - cur_site)%3);

				site_head_tail.push_back((frame_reference - cur_site)%3);

			else

				//site_head_tail.push_back((3 - (cur_site - border) % 3) % 3);

				site_head_tail.push_back((3 - (cur_site - frame_reference) % 3) % 3);

#ifdef DEBUG

			os << "(acceptor)" << site_head_tail.back() << "(frame_reference:" << frame_reference << ")" << "\n";

#endif



		}

		else //for donor, record tail

		{

			//if (cur_site <= border)

			if (cur_site <= frame_reference)

				//site_head_tail.push_back((3 - (border - cur_site)%3) % 3);

				site_head_tail.push_back((3 - (frame_reference - cur_site)%3) % 3);

			else

				//site_head_tail.push_back( (cur_site - border) % 3);

				site_head_tail.push_back( (cur_site - frame_reference) % 3);

#ifdef DEBUG

			os << "(donor)" << site_head_tail.back() << "(frame_reference:" << frame_reference << ")" << "\n";

#endif

		}



	}

}







//find the real position (index) of [search_start, search_end] in the area [align_start, align_end] on the alignStr (targetStr)

//search_start/end, align_start/end are ABSOLUTE positions, and this returns the RELATIVE index of search_start/end on alignStr

//alignStr is now a DNA string

bool FindRealPos_DNA(string& alignStr, int align_start, int align_end, int search_start, int search_end,

						int& search_start_pos, int& search_end_pos)

{

//	if (search_start == 12739337)

//		int stop=1;



	if (alignStr.empty() || search_start < align_start || search_end > align_end) //no str

	{

		//indicate no search area

		return false;

	}



/*

#ifdef DEBUG

	cout << alignStr << "\n" << "string start: " << align_start << "; string end: " << align_end 

		<< "; search_start:" << search_start << "; search_end:" << search_end << "\n";

#endif

*/

	int count = align_start;

	int cur_pos = alignStr.find_first_not_of("+- ");//0;

		

	while (count < search_start)

	{

		if (alignStr[cur_pos] != '-')

			count++;

		cur_pos++;

	}

	search_start_pos = cur_pos;



	while (count < search_end)

	{

		if (alignStr[cur_pos] != '-')

			count++;

		cur_pos++;

	}

	search_end_pos = cur_pos;



	return true;



}



void GetGaps(const string& str, vector<int>& gap_starts, vector<int>& gap_ends)

{

	int pos = 0;

	

	while ((pos = str.find('-', pos)) != string::npos)

	{

		gap_starts.push_back(pos);

		

		pos = str.find_first_not_of('-', pos); //find the end of gap

		if (pos != string::npos)

		{

			gap_ends.push_back(pos-1);

		}

		else //the end of str is gap

		{

			gap_ends.push_back(str.length()-1);

			break;

		}

	}

}



int GetMatches(const string& str, vector<int>& match_starts, vector<int>& match_ends, vector<int>& id)

{

	int pos = 0;

	int count = 0;



	int tmp;



	while ((pos = str.find_first_not_of("+ ", pos)) != string::npos)

	{

		tmp = pos;



		match_starts.push_back(pos);



		pos = str.find_first_of("+ ", pos);

		if (pos != string::npos)

		{

			match_ends.push_back(pos-1);

			count += pos-tmp; //length of the match

			id.push_back(count); //the total identity at the end of current match (from the beginning of str)

		}

		else

		{

			int len = str.length();

			match_ends.push_back(len-1);

			count += len - tmp;

			id.push_back(count);

			break;

		}

	}



	return count; //total number of identities

}



// need_after_gap_pos default is false

int ConvertToNoGappedPos(int pos, vector<int>& gap_starts, vector<int>& gap_ends, bool need_after_gap_pos)

{

	int queryPos = pos;



	vector<int>::iterator end_it = gap_ends.begin();

	vector<int>::iterator start_it = gap_starts.begin();



	int gap_count = 0;



	bool pos_in_gap = false;

	while ( start_it != gap_starts.end() )

	{

		if (*start_it <= pos)

		{

			if (*end_it < pos)

			{

				gap_count += *end_it - *start_it + 1;

			}

			else

			{

				gap_count += pos - *start_it + 1;

				pos_in_gap = true;

			}

		}

		else

		{

			break;

		}



		start_it++;

		end_it++;

	}



	queryPos -= gap_count;



	if ( need_after_gap_pos && pos_in_gap)

		queryPos++;



	return queryPos;

}



int ConvertToGappedPos(int pos, vector<int>& gap_starts, vector<int>& gap_ends)

{

	int matchPos = pos;



	vector<int>::iterator start_it = gap_starts.begin();

	vector<int>::iterator end_it = gap_ends.begin();



	while (start_it != gap_starts.end())

	{

		if (*start_it <= matchPos)

		{

			if (*end_it <= matchPos)

			{

				matchPos += *end_it - *start_it + 1;

			}

			else

			{

				matchPos += matchPos - *start_it + 1;

				break;

			}

		}

		else

			break;



		start_it++;

		end_it++;

	}



	return matchPos;

}



//pos here is the "matchPos", we need to compute the identities from [pos+1] to end of "matchStr"

int ComputeBackId(int pos, vector<int>& match_starts, vector<int>& match_ends, vector<int>& match_id, int total_id)

{

	if (match_ends.empty())

		return 0;



	vector<int>::iterator end_it = lower_bound(match_ends.begin(), match_ends.end(), pos);

	

	if (end_it == match_ends.end()) //if pos is bigger than all "match_ends"

		return 0;



	int index = end_it-match_ends.begin();



	if (*end_it == pos)

		return total_id - match_id[index];

	else //*end_it must be > pos

	{

		if (match_starts[index] < pos)

			return total_id - (match_id[index] - (*end_it - pos));

		else

			return total_id - (match_id[index] - (*end_it - match_starts[index] + 1));

	}

}



//pos here is the "matchPos", we need to compute the identities from [pos-1] to the front of "matchStr"

int ComputeFrontId(int pos, vector<int>& match_starts, vector<int>& match_ends, vector<int>& match_id, int total_id)

{

	if (match_ends.empty()) //no matched region at all

		return 0;



	vector<int>::iterator end_it = lower_bound(match_ends.begin(), match_ends.end(), pos);

	if (end_it == match_ends.end()) //if pos is bigger than all "match_ends"

		return total_id;



	int index = end_it - match_ends.begin();



	if (*end_it == pos)

		return match_id[index] - 1;

	else

	{

		if (match_starts[index] < pos)

			return match_id[index] - (*end_it - pos + 1);

		else

			return (index > 0 ? match_id[index-1] : 0);

	}



}





/*************************************************************************/

/*                                                              	 */

/*  Classify a case description using the given subtree by adjusting	 */

/*  the value ClassSum for each class					 */

/*                                                              	 */

/*************************************************************************/





//MODIFIED: now donor_acceptor_HSP_ID needs to collect hsp_ID (may contain multiple for each donor_segment, 

//and now hsp_ID is not the actual hsp->ID, but the relative index of the HSP!

void GetSegments(Tree T, vector<int>& gap_starts, vector<int>& gap_ends, int hsp_start, 

				 //multimap<int, pair<int, int> >& donor_segments, 

				 vector<SegmentsInThreeBounds>& donor_segments_pair, 

				 vector< vector<int> >& donor_acceptor_HSP_ID, //vector<int>& donor_acceptor_HSP_ID, //vector<int>& acceptor_HSP_ID, 

				 int hsp_ID, //int& prev_hsp_ID, //now this prev_hsp_ID is the HSP_ID for the last_segment (not necessarily the previous HSP)

				 int& last_segment_start, int& last_segment_end, 

				 bool& lastLeafIsExon)//, vector<int>& segment_hsp_start, bool& isFirstLeaf)

{

    DiscrValue v;

    ClassNo c, BestClass;



	float	*ClassSum = new float[MaxClass+1];		/* ClassSum[c] = total weight of class c */

    for(c=0; c<=MaxClass; c++)

		ClassSum[c] = 0;



	bool isExon;

	int start_pos, end_pos;

    switch ( T->NodeType )

    {

        case 0:  /* leaf */



	    if ( T->Items > 0 )

	    {

		/*  Update from ALL classes  */



		ForEach(c, 0, MaxClass)

		{

		    if ( T->ClassDist[c] )

		    {

			ClassSum[c] +=  T->ClassDist[c] / T->Items; //Weight * T->ClassDist[c] / T->Items;

		    }

		}

	    BestClass = 0;

		ForEach(c, 0, MaxClass)

			if ( ClassSum[c] > ClassSum[BestClass] ) BestClass = c;



	    }

	    else

	    {

			BestClass = T->Leaf;

	    }



		//convert to insert entry into donor_segments

		start_pos = ConvertToNoGappedPos(T->Start, gap_starts, gap_ends, true);

		end_pos = ConvertToNoGappedPos(T->End, gap_starts, gap_ends, false);



		isExon = BestClass; //class can only be 0 or 1



		if (isExon)

		{

			if (!lastLeafIsExon) //is this possible?

			{

				last_segment_start = hsp_start + start_pos*3; //this is the start of a new donor_segment

				donor_acceptor_HSP_ID.push_back( vector<int>(1, hsp_ID) );

			}

			//else //last leaf is also exon, so they'll be combined into one exon segment, and prev_hsp_ID should be recorded

				//donor_acceptor_HSP_ID.push_back( vector<int>(1, prev_hsp_ID) );

				//donor_acceptor_HSP_ID.back().push_back(prev_hsp_ID);



			last_segment_end = hsp_start + end_pos*3 + 3; //2;



			//prev_hsp_ID = hsp_ID; //update prev_hsp_ID so it always corresponds to the current last_segment,

			if (donor_acceptor_HSP_ID.back().back() != hsp_ID)

				donor_acceptor_HSP_ID.back().push_back(hsp_ID); //only insert when there's no such HSP_ID already

			//so there's no need for isFirstLeaf any more



			//if (isFirstLeaf)

			//	isFirstLeaf = false;

		}

		else

		{

			if (lastLeafIsExon)

			{

				//donor_segments.insert(multimap<int, pair<int, int> >::value_type(last_segment_end, 

				//	pair<int, int>(last_segment_start, hsp_start + end_pos*3 + 2)));

				donor_segments_pair.push_back(SegmentsInThreeBounds(last_segment_start, last_segment_end, hsp_start + end_pos*3 + 2 ) );

				//if (isFirstLeaf)

				//{

					//donor_acceptor_HSP_ID.push_back(prev_hsp_ID);

					//donor_acceptor_HSP_ID.back().insert(prev_hsp_ID);

				//	isFirstLeaf = false;

				//}

				//else

				//{

				//	donor_acceptor_HSP_ID.push_back(hsp_ID);

				//}

				//gap_centers.push_back((last_segment_end + hsp_start + end_pos*3 + 2)/2);



				//segment_hsp_start.push_back(hsp_start); segment_hsp_start not used any more

			}

			else

			{

				if (!donor_segments_pair.empty())

				{

					//multimap<int, pair<int, int> >::reverse_iterator seg_it = donor_segments.rbegin();

					//const_cast<int&>((*seg_it).second.second) = hsp_start + end_pos*3 + 2;

					donor_segments_pair[donor_segments_pair.size()-1].intron_seg_end = hsp_start + end_pos*3 + 2;

					//gap_centers.back() = ( (*seg_it).second.second + (*seg_it).first )/2;

				}

				//if (isFirstLeaf)

				//	isFirstLeaf = false;

			}

		}

		lastLeafIsExon = isExon;



	    break;



	case ThreshContin:  /* test of continuous attribute */



		//v = ( Cv <= T->Cut ? 1 : 2 );

		for (v=1; v<=2; v++)

			GetSegments(T->Branch[v], gap_starts, gap_ends, hsp_start, 

				//donor_segments, 

				donor_segments_pair, donor_acceptor_HSP_ID, //acceptor_HSP_ID, 

				hsp_ID, //prev_hsp_ID, 

				last_segment_start, last_segment_end, lastLeafIsExon);//, segment_hsp_start, isFirstLeaf);



	    break;



    } 



	delete [] ClassSum;

}



//not used

void ProcSegmentsByLens(multimap<int, pair<int, int> >& donor_segments, int& last_segment_start, 

						vector<int>& segment_hsp_start, int& last_hsp_start) //MIN_INTERNAL_EXON_LEN is currently not used!

{

	if (donor_segments.empty())

		return;



	int s, m, e;

	//int last_s, last_m, last_e;



	vector< multimap<int, pair<int, int> >::iterator > seg_to_be_erased;

	vector< vector<int>::iterator > seg_hsp_start_to_be_erased;

	

	bool init_erase = true;

	while (init_erase || seg_to_be_erased.size() > 0)

	{

		init_erase = false;

		seg_to_be_erased.clear();

		seg_hsp_start_to_be_erased.clear();



	multimap<int, pair<int, int> >::iterator seg_it, prev_seg_it = donor_segments.begin();

	

	bool good_exon, good_intron, last_good_intron = true, last_good_exon = true;

	//vector< pair<int, int> > exon_segments, intron_segments;

	vector<int>::iterator seg_hs_it = segment_hsp_start.begin(), prev_seg_hs_it = segment_hsp_start.begin();

	int i = 0;

	for (seg_it = donor_segments.begin(); seg_it != donor_segments.end(); seg_it++, seg_hs_it++, i++ )

	{

		cout << "checking segment " << i << ": " << "\n";



		s = (*seg_it).second.first;

		m = (*seg_it).first;

		e = (*seg_it).second.second;



		if (seg_it != donor_segments.begin() && m - s < MIN_INTERNAL_EXON_LEN) //internal exon test failed

			good_exon = false;

		else

			good_exon = true;



		if (e - m + 1 < MIN_INTRON_LEN) //intron test failed

			good_intron = false;

		else

			good_intron = true;



		if (good_exon)

		{

			if (!last_good_intron && last_good_exon)

			{

				const_cast<int&>((*prev_seg_it).first) = m;

				const_cast<int&>((*prev_seg_it).second.second) = e;

				seg_to_be_erased.push_back(seg_it);

				seg_hsp_start_to_be_erased.push_back(seg_hs_it);

				

				cout << "will erase last_intron " << m << " to " << e << "\n";

			}

			else

			{

				prev_seg_it = seg_it;

				prev_seg_hs_it = seg_hs_it;

			}



			last_good_exon = true;

			last_good_intron = good_intron;

		}

		else //!good_exon

		{

			//if !good_exon && both introns around the bad exon are good, then remove exon, otherwise keep everything

			if (last_good_intron && good_intron)

			{

				const_cast<int&>((*prev_seg_it).second.second) = e;

				seg_to_be_erased.push_back(seg_it);

				seg_hsp_start_to_be_erased.push_back(seg_hs_it);



				cout << "will erase cur_exon " << s << " to " << m << "\n";



				//last_good_intron = true;

				last_good_exon = (prev_seg_it == donor_segments.begin() || 

					((*prev_seg_it).first - (*prev_seg_it).second.first >= MIN_INTERNAL_EXON_LEN));

			}

			else

			{

				last_good_intron = good_intron;

				last_good_exon = good_exon;



				prev_seg_it = seg_it;

				prev_seg_hs_it = seg_hs_it;

			}

		}

	}

	if (last_good_exon && !last_good_intron)// && last_segment_end - last_segment_start + 1 >= MIN_INTERNAL_EXON_LEN)

	{

		last_segment_start = s;



		seg_it--; //seg_it now point to last segment

		if (prev_seg_it == seg_it) //last segment has not been deleted previously

		{

			seg_to_be_erased.push_back(prev_seg_it);

			seg_hsp_start_to_be_erased.push_back(prev_seg_hs_it);

		}

		else //otherwise, prev_seg_it must be before seg_it, need to insert to proper place!

		{

			vector< multimap<int, pair<int, int> >::iterator >::iterator tmp_seg_it = seg_to_be_erased.end();

			tmp_seg_it--;

			seg_to_be_erased.insert(tmp_seg_it, prev_seg_it);

			

			vector< vector<int>::iterator >::iterator tmp_seg_start_it = seg_hsp_start_to_be_erased.end();

			tmp_seg_start_it--;

			seg_hsp_start_to_be_erased.insert(tmp_seg_start_it, prev_seg_hs_it);

		}



		cout << "will erase last_intron_before_last_segment " << (*prev_seg_it).first << " to " << (*prev_seg_it).second.second << "\n";

	}



	for (i=seg_to_be_erased.size()-1; i>=0; i--)

	{

		donor_segments.erase(seg_to_be_erased[i]);

		cout << "erased segment [" << i << "]" << "\n";



		segment_hsp_start.erase(seg_hsp_start_to_be_erased[i]);

		cout << "erased segment_hsp_start " << *(seg_hsp_start_to_be_erased[i]) << "\n";

	}



	}



}


//no longer used
bool GetAlignPairStr(string& curStr, char t, char q, bool open_gap)
{
	bool open_gap_cur;
	//MODIFIED: treat 'U' as '*' (stop codon)
	t = (t=='U')? '*' : t;
	q = (q=='U')? '*' : q;

		if (t == '-' || q == '-' )
		{
			if (!open_gap)
				curStr = "-o";
			else
				curStr = "-e";
			open_gap_cur = true;
		}
		else
		{
			open_gap_cur = false;

			if (t == 'X' || t == 'x' || 
				q == 'X' || q == 'x')
				curStr = "Xx";
			else
			{
				//MODIFIED: also treat 'U' as '*'
				if (t == '*' || q == '*' )//|| t == 'U' || q == 'U') //stop codon
				{
					if (q != t)
						curStr = "*x";
					else
						curStr = "**";
				}
				else
				{
					curStr = t;
					curStr += q;			
				}
			}
		}

	return open_gap_cur;
}



//string trim function
void remove_trailing_cr_lf(string& str)
{
	int len = str.length();
	if (len == 0)
		return;
	const char* whitespaces = "\r\n"; //CR/LF
	int last_non_space = str.find_last_not_of(whitespaces);
	if (last_non_space == string::npos) //string::npos is -1 in all platforms?
		last_non_space=-1;
	if (last_non_space < len-1)
		str.erase(last_non_space+1);
}

int remove_trailing_cr_lf(char* str, int len)
{
	int i = len - 1;
	while (i>=0 && (str[i] == '\r' || str[i] == '\n'))
		i--;
	str[i+1] = '\0'; //no need? since we're returning the exact length of string
	return i+1;
}

bool all_white_space(string& str)
{
	const char* whitespaces = " \t\r\n"; //space and tab and newline
	int fst_non_space = str.find_first_not_of(whitespaces);
	if (fst_non_space == string::npos)
		return true;
	else
		return false;
}





void DNA2AA(string& target_seq, int start_index, string& result_seq)

{

	result_seq = "";



	map<string, char>::iterator it;

	string tmp;

	

	for (int i=start_index; i+3<=target_seq.length(); i+=3)

	{

		tmp = target_seq.substr(i, 3);		

		it = DNA_CODON_TBL_SL.find(tmp);

		if (it == DNA_CODON_TBL_SL.end())

			result_seq += 'X';

		else

			result_seq += (*it).second;	

		

	}

}



bool TranslateTargetDNA(string& str)

{

	string codonStr;

	map<string, char>::iterator codonIt;

	for (int i=0; i<str.length(); i+=3)

	{

		codonStr = str.substr(i, 3);

		if (codonStr[0] >= 'a') //lowercase, we convert to protein (upper case must already be protein, don't need to change)

		{

			codonIt = DNA_CODON_TBL_SL.find(codonStr);

			if (codonIt != DNA_CODON_TBL_SL.end() )

			{

				if ((*codonIt).second == '*') //stop codon

					return false;

				else

					str.replace(i, 3, 3, (*codonIt).second); //replace codonStr part with the amino acid letter (repeated 3 times)

			}
			else //if not in codon table?
				str.replace(i, 3, 3, 'X');

		}

	}

	return true;

}


//stop_pos points to the first stop codon position (relative index)
bool HasInFrameStop(string& str, bool frame_start_at_str_beginning, int& stop_pos, bool ignore_last_codon)

{

	int i;

	string tempStr;

	if (frame_start_at_str_beginning)
		i = 0;
	else //frame start from the end of string
		i = str.length() % 3;

	int len = str.length();
	if (ignore_last_codon)
		len -= 3;
		
	while (i+3 <= len )
	{
		tempStr = str.substr(i, 3);
		map<string, char>::iterator codon_it = DNA_CODON_TBL_SL.find(tempStr);
		if (codon_it != DNA_CODON_TBL_SL.end())
		{
			if ((*codon_it).second == '*')
			{
				stop_pos = i;
				return true;
			}
		}
		i+=3;
	}

	return false;
}


//erase all occurrences of "ch" in "str"
void EraseAll(string& str, char ch)
{

	int str_pos_index = str.length()-1;

	for (; str_pos_index >= 0; str_pos_index--) 

		if (str[str_pos_index] == ch)

			str.erase(str_pos_index, 1);

}



bool HasStopCodon(int acceptor_site, int donor_site, vector<string>& seq, int& frame, //int& site_before_stop, //frame here is the acceptor_head
				  int& site_after_stop, 
				  string& align, int prev_donor_site, string& left_chars, 
				  bool output_DNA, ofstream& DNA_os, bool temp_DNA, string& temp_DNA_string, //if needs to output DNA sequence, print the sequence to DNA_os file
				  bool is_last_exon, bool& found_stop) 
{
	found_stop = false;

	string exon_str;
	if (acceptor_site > 0)
		GetSubstrFromVecStrs_NegRev(seq, true, acceptor_site-1, donor_site-acceptor_site+1, exon_str);
	else
		GetSubstrFromVecStrs_NegRev(seq, false, -donor_site-1, donor_site-acceptor_site+1, exon_str);
	StrToLower(exon_str);

/*	if (output_DNA)
		DNA_os << exon_str;
	if (temp_DNA)
		temp_DNA_string += exon_str;
*/
	string cur_frame_str;
	if (frame > 0)
	{
		cur_frame_str = left_chars;
		cur_frame_str += exon_str.substr(0, frame);
		map<string, char>::iterator it = DNA_CODON_TBL_SL.find(cur_frame_str);
		if (it != DNA_CODON_TBL_SL.end())
		{
			
			if ((*it).second == '*') //stop codon
			{
				//site_before_stop = prev_donor_site - (3-frame);
				site_after_stop = acceptor_site + frame - 1;
				found_stop = true;

				if (output_DNA)
					DNA_os << cur_frame_str;
				if (temp_DNA)
					temp_DNA_string += cur_frame_str;

				return true;
			}
			else
			{
				align += (*it).second;
			}
		}
		else
		{
			//cout << "cannot find protein codon for " << cur_frame_str << "!" << "\n";
			//exit(-1);
			align += 'X';
		}

		if (output_DNA)
			DNA_os << cur_frame_str;
		if (temp_DNA)
			temp_DNA_string += cur_frame_str;
	}

	int t = frame;
	while (t+3 <= exon_str.length() )
	{
		cur_frame_str = exon_str.substr(t, 3);
		map<string, char>::iterator it = DNA_CODON_TBL_SL.find(cur_frame_str);
		if (it != DNA_CODON_TBL_SL.end())
		{
			if ((*it).second == '*') //stop codon
			{
				if (output_DNA)
					DNA_os << cur_frame_str;
				if (temp_DNA)
					temp_DNA_string += cur_frame_str;

				site_after_stop = acceptor_site + t + 2; //this is the site after stop!
				//site_before_stop = acceptor_site + t - 1; //this is the site just before stop
				found_stop = true;

				if (is_last_exon && t == exon_str.length() - 3) //if it's the last 3bp of last exon
					return false; //signals that this is a valid situation (gene end at stop codon)
				else
					return true;
			}
			else
			{
				align += (*it).second;
			}
		}
		else
		{
			//cout << "cannot find protein codon for " << cur_frame_str << "!" << "\n";
			//exit(-1);
			align += 'X';
		}

		if (output_DNA)
			DNA_os << cur_frame_str;
		if (temp_DNA)
			temp_DNA_string += cur_frame_str;

		t += 3;
	}
	frame = (3 - (exon_str.length()-t)) % 3;
	if (frame == 1)
		left_chars = exon_str.substr(exon_str.length()-2);
	else
		if (frame == 2)
			left_chars = exon_str.substr(exon_str.length()-1);

	if (output_DNA)
		DNA_os << "\n";
	if (temp_DNA)
		temp_DNA_string += "\n";
	return false;
}

void HasStopCodon_KeepStopCodon(int acceptor_site, int donor_site, vector<string>& seq, int& frame, 

				  string& align, string& left_chars)

{

	string exon_str;

	if (acceptor_site > 0)

		GetSubstrFromVecStrs_NegRev(seq, true, acceptor_site-1, donor_site-acceptor_site+1, exon_str);

	else

		GetSubstrFromVecStrs_NegRev(seq, false, -donor_site-1, donor_site-acceptor_site+1, exon_str);

	StrToLower(exon_str);

							

	string cur_frame_str;

	if (frame > 0)

	{

		cur_frame_str = left_chars;

		cur_frame_str += exon_str.substr(0, frame);

		map<string, char>::iterator it = DNA_CODON_TBL_SL.find(cur_frame_str);

		if (it != DNA_CODON_TBL_SL.end())

		{

			align += (*it).second;

		}

		else

		{

			//cout << "cannot find protein codon for " << cur_frame_str << "!" << "\n";

			//exit(-1);

			align += 'X';

		}

	}



	int t = frame;

	while (t+3 <= exon_str.length() )

	{

		cur_frame_str = exon_str.substr(t, 3);

		map<string, char>::iterator it = DNA_CODON_TBL_SL.find(cur_frame_str);

		if (it != DNA_CODON_TBL_SL.end())

		{

			align += (*it).second;

		}

		else

		{

			//cout << "cannot find protein codon for " << cur_frame_str << "!" << "\n";

			//exit(-1);

			align += 'X';

		}

		t += 3;

	}

	frame = (3 - (exon_str.length()-t)) % 3;

	if (frame == 1)

		left_chars = exon_str.substr(exon_str.length()-2);

	else

		if (frame == 2)

			left_chars = exon_str.substr(exon_str.length()-1);



}

void CombineConsecutiveSpaces(string& str)
{
	string new_str ="";
	bool prev_is_space = false;
	int i;
	for (i=0; i<str.length(); i++)
		if (str[i] == ' ' || str[i] == '\t' || str[i] == '\n')
		{
			if (prev_is_space)
				continue;
			else
			{
				new_str += ' ';
				prev_is_space = true;
			}
		}
		else
			new_str += str[i];

	str = new_str;
}




