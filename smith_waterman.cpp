#if defined(_MSC_VER)
	#pragma warning(disable: 4786)
	#pragma warning(disable: 4503)
#endif


#include <iostream>

#include <fstream>

#include <cstdlib>

#include <string>

#include <cmath>

#include <time.h> //<sys/time.h>

#include <map>

using namespace std;



#include "data_manager.h"

#include "Hirschberg.h"



int similarity_score(char a,char b, bool opengap, bool& opengap_cur);
int sim_score(char a,char b);

template <class T>
T find_array_max(T array[],int length);

//void insert_at(char arr[], int n, int idx, char val);

//void checkfile(int open, char filename[]);

//string read_sequence(ifstream& f);

//void ReverseString(char *ptr );



extern bool GetAlignPairStr(string& curStr, char t, char q, bool open_gap);

extern void GetSubstrFromVecStrs(//pair<int, char*>& chr_seq, 
								 vector<string>& chr_seq, 
								 bool isPosStrand, const int start_pos, const int len, string& resultStr);

extern void GetSubstrFromVecStrs_NegRev(//pair<int, char*>& chr_seq, 
										vector<string>& chr_seq, 
										bool isPosStrand, const int start_pos, const int len, string& resultStr);



int ind;

double mu,delta;

extern int ALIGN_SCORE_MATRIX[28][28]; //map<string, int> ALIGN_SCORE_MATRIX;


//alignment score is computed based on: match is 1, mismatch/deletion/insertion is 0
void AssignScores_std(int* temp, int H_im1_jm1, int H_im1_j, int H_i_jm1, char seq_a_im1, char seq_b_jm1)
{
  temp[0] = H_im1_jm1+sim_score(seq_a_im1,seq_b_jm1);
  temp[1] = H_im1_j;
  temp[2] = H_i_jm1;
}

//use score matrix to compute score
template <class T>
void AssignScores_scorematrix(T* temp, T H_im1_jm1, T H_im1_j, T H_i_jm1, char seq_a_im1, char seq_b_jm1, 
	bool& opengap_im1_jm1, bool& opengap_im1_j, bool& opengap_i_jm1, bool* temp_opengap)
{
	temp[0] = H_im1_jm1+similarity_score(seq_a_im1,seq_b_jm1, opengap_im1_jm1, temp_opengap[0]);
	temp[1] = H_im1_j+similarity_score(seq_a_im1, '-', opengap_im1_j, temp_opengap[1]);
	temp[2] = H_i_jm1+similarity_score(seq_b_jm1, '-', opengap_i_jm1, temp_opengap[2]);
	//temp[3] = 0.;
}


//This is the Smith-Waterman algirthm (a variation of Needleman-Wunsch algorithm), finds one optimal local alignment

//close_to_start indicates which alignment to take, in case there're multiple local alignments with same score

//we may take either the one closest to seq start, or closest to seq end

double GetLocalAlignment(string& seq_a, string& seq_b, int seq_a_start, int seq_b_start, int newId, 

					   HSP_Gene_Pair& newHSP, Input_Alignment& newAlign, bool close_to_start, ofstream& debugFile, 

					   int* align_pos)

{

//	cout << "seq_a:(" << seq_a_start << "..)" << seq_a << "\n";

//	cout << "seq_b:(" << seq_b_start << "..)" << seq_b << "\n";



  // string s_a=seq_a,s_b=seq_b;

  int N_a = seq_a.length();                     // get the actual lengths of the sequences

  int N_b = seq_b.length();

 

  ////////////////////////////////////////////////



  // initialize H

  //double H[N_a+1][N_b+1];     

  int i, j;

  double** H;

  int **I_i, **I_j;

  H = new double*[N_a+1];

  I_i = new int*[N_a+1];

  I_j = new int*[N_a+1];



  bool **opengap; //use this to keep track whether a gap is opengap

  opengap = new bool*[N_a+1];



  for(i=0;i<=N_a;i++){

	  H[i] = new double[N_b+1];

	  I_i[i] = new int[N_b+1];

	  I_j[i] = new int[N_b+1];

	  opengap[i] = new bool[N_b+1];



	  H[i][0]=0; //only need to initialize first row and first column

	  opengap[i][0] = false;

//    for(j=0;j<=N_b;j++){

//      H[i][j]=0.;

//    }

  } 

  for (j=0;j<=N_b;j++){

	  H[0][j]=0; //initialize first column

	  opengap[0][j] = false;

  }



  double temp[4];

  bool temp_opengap[3];

  //int I_i[N_a+1][N_b+1],I_j[N_a+1][N_b+1];     // Index matrices to remember the 'path' for backtracking

  

  // here comes the actual algorithm



  for(i=1;i<=N_a;i++){

    for(j=1;j<=N_b;j++){

      //temp[0] = H[i-1][j-1]+similarity_score(seq_a[i-1],seq_b[j-1]); 

      //temp[1] = H[i-1][j]-delta;                  

      //temp[2] = H[i][j-1]-delta;

	  //cout << "i:" << i << ";j:" << j << "\n";

	  AssignScores_scorematrix<double>(temp, H[i-1][j-1], H[i-1][j], H[i][j-1], seq_a[i-1], seq_b[j-1], 
		  opengap[i-1][j-1], opengap[i-1][j], opengap[i][j-1], temp_opengap);
	  /*temp[0] = H[i-1][j-1]+similarity_score(seq_a[i-1],seq_b[j-1], opengap[i-1][j-1], temp_opengap[0]);
      temp[1] = H[i-1][j]+similarity_score(seq_a[i-1], '-', opengap[i-1][j], temp_opengap[1]);
      temp[2] = H[i][j-1]+similarity_score(seq_b[j-1], '-', opengap[i][j-1], temp_opengap[2]);*/

      temp[3] = 0.;

      H[i][j] = find_array_max<double>(temp,4);

      switch(ind){

      case 0:                                  // score in (i,j) stems from a match/mismatch

   		I_i[i][j] = i-1;

		I_j[i][j] = j-1;

		opengap[i][j] = temp_opengap[0];

		break;

      case 1:                                  // score in (i,j) stems from a deletion in sequence A

     	I_i[i][j] = i-1;

		I_j[i][j] = j;

		opengap[i][j] = temp_opengap[1];

		break;

      case 2:                                  // score in (i,j) stems from a deletion in sequence B

      	I_i[i][j] = i;

		I_j[i][j] = j-1;

		opengap[i][j] = temp_opengap[2];

		break;

      case 3:                                  // (i,j) is the beginning of a subsequence

      	I_i[i][j] = i;

		I_j[i][j] = j;

		opengap[i][j] = false;

		break;

      }

    }

  }



/*  // Print the matrix H to the console

  cout<<"**********************************************"<<"\n";

  cout<<"The scoring matrix is given by  "<<"\n"<<"\n";

  for(int i=1;i<=N_a;i++){

    for(int j=1;j<=N_b;j++){

      cout<<H[i][j]<<" ";

    }

    cout<<"\n";

    }

*/

/*#ifdef DEBUG

  debugFile <<"**********************************************"<<"\n";

  debugFile<<"The scoring matrix is given by  "<<"\n"<<"\n";

  for( i=1;i<=N_a;i++){

    for( j=1;j<=N_b;j++){

      debugFile<<H[i][j]<<" ";

    }

    debugFile<<"\n";

    }

  debugFile <<"**********************************************"<<"\n";

#endif

*/

  // search H for the maximal score

  double H_max = 0.;

  int i_max=0,j_max=0;

  for( i=1;i<=N_a;i++){

    for( j=1;j<=N_b;j++){

		if (close_to_start)

		{

			if(H[i][j]>H_max){

				H_max = H[i][j];

				i_max = i;

				j_max = j;

			}

		}

		else

		{

			if(H[i][j]>=H_max){ //use >= instead of >, so we update to the later position even if both are same

				H_max = H[i][j];

				i_max = i;

				j_max = j;

			}

		}

    }

  }



  //what if H_max is 0?

  if (H_max == 0)

  {

#ifdef DEBUG

	  debugFile << seq_a << "\n" << seq_b << "\n" << "H_max is 0, no alignment!" << "\n";

#endif

	  return 0.;

  }

  

  //cout<<H_max<<"\n";

#ifdef DEBUG

  debugFile << "H_max:" << H_max << "(" << i_max << "," << j_max << ")" << "\n";

#endif



  align_pos[1] = i_max-1;

  align_pos[3] = j_max-1;



     // Backtracking from H_max

  int current_i=i_max,current_j=j_max;

  int next_i=I_i[current_i][current_j];

  int next_j=I_j[current_i][current_j];

//  int tick=0;

  //char consensus_a[N_a+N_b+2],consensus_b[N_a+N_b+2];

  string consensus_a(""), consensus_b(""), match("");



  int seq_a_align_start, seq_b_align_start, seq_a_align_end, seq_b_align_end;

  seq_a_align_end = seq_a_start + i_max - 1;

  seq_a_align_start = seq_a_align_end;



  seq_b_align_start = seq_b_start + 3*(j_max - 1);

  seq_b_align_end = seq_b_align_start + 2; //3 base pairs is one letter

  newAlign.match_align = newAlign.target_align = newAlign.query_align = "";

  int num_of_match = 0;



  while(((current_i!=next_i) || (current_j!=next_j)) && (next_j!=0) && (next_i!=0)){



    if(next_i==current_i)

	{

		consensus_a += '-'; //consensus_a[tick] = '-';                  // deletion in A

		match += ' ';

		consensus_b += seq_b[current_j-1];	//b must be some actual char, cannot be '-' aigns with '-'!

		seq_b_align_start-=3;

	}

    else

	{

		consensus_a += seq_a[current_i-1]; //consensus_a[tick] = seq_a[current_i-1];   // match/mismatch in A

		seq_a_align_start--;

	    if(next_j==current_j)

		{

			consensus_b += '-'; //consensus_b[tick] = '-';                  // deletion in B

			match += ' ';

		}

		else

		{

			consensus_b += seq_b[current_j-1]; //consensus_b[tick] = seq_b[current_j-1];   // match/mismatch in B

			seq_b_align_start-=3;

			if (seq_a[current_i-1] == seq_b[current_j-1])

			{

				match += seq_a[current_i-1];

				num_of_match++;

			}

			else

			{

				match += ' ';

			}

		}

	}



    current_i = next_i;

    current_j = next_j;

    next_i = I_i[current_i][current_j];

    next_j = I_j[current_i][current_j];

//    tick++;

    }



    //record the last one

//    if ((next_i==current_i) && (next_j == current_j))

//	{

		consensus_a += seq_a[current_i-1];

		consensus_b += seq_b[current_j-1];

		if (seq_a[current_i-1] == seq_b[current_j-1])

		{

			match += seq_a[current_i-1];

			num_of_match++;

		}

		else

		{

			match += ' ';

		}



		align_pos[0] = current_i-1;

		align_pos[2] = current_j-1;

//	}

//	else

//	{

//		seq_a_align_start++; //forwarded one extra position, so now get it back

//		seq_b_align_start+=3;

//	}



  newHSP.gene_start = seq_a_align_start;

  newHSP.gene_end = seq_a_align_end;

  if (seq_b_start > 0)

  {

	newHSP.HSP_start = seq_b_align_start;

	newHSP.HSP_end = seq_b_align_end;

  }

  else //neg strand

  {

	  newHSP.HSP_start = -seq_b_align_end;

	  newHSP.HSP_end = -seq_b_align_start;

  }

  newHSP.pid = ((float)(num_of_match) / match.length())*100;



  newHSP.ID = newId;



/*  for (i=match.length()-1; i>=0; i--) //reverse strings

  {

	  newAlign.query_align += consensus_a[i];

	  newAlign.match_align += match[i];

	  newAlign.target_align += consensus_b[i];

  }

*/

  reverse(consensus_a.begin(), consensus_a.end());

  reverse(consensus_b.begin(), consensus_b.end());

  reverse(match.begin(), match.end());

  newAlign.query_align += consensus_a;

  newAlign.match_align += match;

  newAlign.target_align += consensus_b;



 // Output of the consensus motif to the console

/*  cout<<"\n"<<"***********************************************"<<"\n";

  cout<<"The alignment of the sequences"<<"\n"<<"\n";

  for( i=0;i<N_a;i++){cout<<seq_a[i];}; cout<<"  and"<<"\n";

  for( i=0;i<N_b;i++){cout<<seq_b[i];}; cout<<"\n"<<"\n";

  cout<<"is for the parameters  mu = "<<mu<<" and delta = "<<delta<<" given by"<<"\n"<<"\n";  

  for( i=tick-1;i>=0;i--) cout<<consensus_a[i]; 

  cout<<"\n";

  for( j=tick-1;j>=0;j--) cout<<consensus_b[j];

  cout<<"\n";

*/

#ifdef DEBUG

  debugFile <<"The alignment of the sequences"<<"\n";

	debugFile << "(" << seq_a_start << "..)" << "\n" << seq_a << "\n";

	debugFile << "(" << seq_b_start << "..)" << "\n" << seq_b << "\n";

	debugFile << newAlign << "\n";

	debugFile << newHSP << "\n";

#endif



  //clean up memory

  for(i=0;i<=N_a;i++){

	  delete [] H[i];

	  delete [] I_i[i];

	  delete [] I_j[i];

	  delete [] opengap[i];

  } 

  delete [] H;

  delete [] I_i;

  delete [] I_j;

  delete [] opengap;



  return H_max;

} // END of main







/////////////////////////////////////////////////////////////////////////////

// auxiliary functions used by main:

/////////////////////////////////////////////////////////////////////////////





void checkfile(int open, char filename[]){



  if (open){cout << "Error: Can't open the file "<<filename<<"\n";exit(1);}

  else cout<<"Opened file "<<filename<<"\n";

}



/////////////////////////////////////////////////////////////////////////////
int score_index(char a)
{
  if (a <= 'Z' && a >= 'A')
	  return a - 'A';
  else
	  if (a <= 'z' && a >= 'a')
		return a - 'a';
	  else
		if (a == '-')
			return 26;
		else
		{
			if (a == '*')
				return 27;
			else
			{
				cout << "cannot find score when running SW local alignment: " << "char a:" << a << "\n";
				exit(-1);
			}
		}
}


int similarity_score(char a,char b, bool opengap, bool& opengap_cur){
  /*string tmp;
  opengap_cur = GetAlignPairStr(tmp, a, b, opengap);
  map<string, int>::iterator score_it = ALIGN_SCORE_MATRIX.find(tmp);
  if (score_it == ALIGN_SCORE_MATRIX.end())
  {
	  cout << "cannot find score when running SW local alignment: " << tmp << "(char a:" << a << ",char b:" << b << ")" << "\n";
	  exit(-1);
  }
  return (*score_it).second;*/

  //ALIGN_SCORE_MATRIX[28][28] stores the align scores, 26 Uppercase letters, '-', '*'
  int ai = score_index(a);
  int bi = score_index(b);
  if (ai != 26 && bi != 26)
  {
	  opengap_cur = false;
	  return ALIGN_SCORE_MATRIX[ai][bi];
  }
  else
  {
	  opengap_cur = true;
	  if (opengap) //previous is gap
		  return ALIGN_SCORE_MATRIX[26][1]; //gap extension
	  else
		  return ALIGN_SCORE_MATRIX[26][0]; //gap opening
  }
		
}



/////////////////////////////////////////////////////////////////////////////


template <class T>
T find_array_max(T array[],int length){

  T max = array[0];            // start with max = first element
  ind = 0;

  for(int i = 1; i<length; i++){
      if(array[i] > max){
		max = array[i];
		ind = i; 
      }
  }
  return max;                    // return highest value in array
}

/////////////////////////////////////////////////////////////////////////////


//not used
string read_sequence(ifstream& f)

{

  // overflows.

  string seq;

  char line[5000];

  while( f.good() )

    {

      f.getline(line,5000);

      // 		cout << "Line:" << line << "\n";

      if( line[0] == 0 || line[0]=='#' )

	continue;

      for(int i = 0; line[i] != 0; ++i)

	{

	  int c = toupper(line[i]);

	  if( c != 'A' && c != 'G' && c != 'C' && c != 'T' )

	    continue;

	  //cout << char(c);

	  //seq.push_back(char(c));

	  seq += char(c);

	}

    }

  return seq;

}



/////////////////////////////////////////////////////////////////////////////



int sim_score(char a,char b){



	if (a == b)

		return 1;

	else

		return 0;



}



//////////////////////////////////////////////////////////////////////////////
//This is the Needleman-Wunsch algorithm, which returns the optimal global alignment
int GetGlobalAlignment(string& seq_a, string& seq_b, Input_Alignment& newAlign, float& align_pid, ofstream& debugFile)
{

//	cout << "seq_a:(" << seq_a_start << "..)" << seq_a << "\n";

//	cout << "seq_b:(" << seq_b_start << "..)" << seq_b << "\n";



  // string s_a=seq_a,s_b=seq_b;

  int N_a = seq_a.length();                     // get the actual lengths of the sequences

  int N_b = seq_b.length();



//For long sequences when we don't have enough space to hold the H arrays, we'll need to use 

//Hirschberg's algorithm (a modification to Needleman-Wunsch). 

//It also returns optimal global alignment, but with only O(m+n) space and 

//at the cost of roughly twice as much time as Needleman-Wunsch?

  if (N_a * N_b > MEMORY_LIMIT_SW) //80M, this roughly translates to 80M*4(sizeof(int))*3(num_of_arrays) ~ 960M memory)

	  return (ComputeEditScript(const_cast<char*>(seq_a.c_str()), N_a, const_cast<char*>(seq_b.c_str()), N_b, 

		newAlign, align_pid));



  ////////////////////////////////////////////////



  // initialize H

  //double H[N_a+1][N_b+1];     

  int i, j;

  int** H;

  int **I_i, **I_j;

  H = new int*[N_a+1];

  I_i = new int*[N_a+1];

  I_j = new int*[N_a+1];



  for(i=0;i<=N_a;i++){

	  H[i] = new int[N_b+1];

	  I_i[i] = new int[N_b+1];

	  I_j[i] = new int[N_b+1];



	  H[i][0]=0; //only need to initialize first row and first column

//    for(j=0;j<=N_b;j++){

//      H[i][j]=0.;

//    }

  } 

  for (j=0;j<=N_b;j++){

	  H[0][j]=0; //initialize first column

  }



  int temp[3];

  //int I_i[N_a+1][N_b+1],I_j[N_a+1][N_b+1];     // Index matrices to remember the 'path' for backtracking

  

  // here comes the actual algorithm



  for(i=1;i<=N_a;i++){

    for(j=1;j<=N_b;j++){

      //temp[0] = H[i-1][j-1]+similarity_score(seq_a[i-1],seq_b[j-1]); 

      //temp[1] = H[i-1][j]-delta;                  

      //temp[2] = H[i][j-1]-delta;

	  //cout << "i:" << i << ";j:" << j << "\n";

	  temp[0] = H[i-1][j-1]+sim_score(seq_a[i-1],seq_b[j-1]);

      temp[1] = H[i-1][j];

      temp[2] = H[i][j-1];

      H[i][j] = find_array_max<int>(temp,3);

      switch(ind){

      case 0:                                  // score in (i,j) stems from a match/mismatch

   		I_i[i][j] = i-1;

		I_j[i][j] = j-1;

		break;

      case 1:                                  // score in (i,j) stems from a deletion in sequence A

     	I_i[i][j] = i-1;

		I_j[i][j] = j;

		break;

      case 2:                                  // score in (i,j) stems from a deletion in sequence B

      	I_i[i][j] = i;

		I_j[i][j] = j-1;

		break;

      }

    }

  }



/*  // Print the matrix H to the console

  cout<<"**********************************************"<<"\n";

  cout<<"The scoring matrix is given by  "<<"\n"<<"\n";

  for(int i=1;i<=N_a;i++){

    for(int j=1;j<=N_b;j++){

      cout<<H[i][j]<<" ";

    }

    cout<<"\n";

    }

*/

/*

#ifdef DEBUG

  debugFile <<"**********************************************"<<"\n";

  debugFile<<"The scoring matrix is given by  "<<"\n"<<"\n";

  for( i=1;i<=N_a;i++){

    for( j=1;j<=N_b;j++){

      debugFile<<H[i][j]<<" ";

    }

    debugFile<<"\n";

    }

  debugFile <<"**********************************************"<<"\n";

#endif

*/



  // search H for the maximal score

  int H_max = H[N_a][N_b];

  int i_max=N_a,j_max=N_b;



  //what if H_max is 0?

//  if (H_max == 0)

//  {

//	  debugFile << "global H_max is 0, no alignment!" << "\n";

//	  return 0.;

//  }

  

  //cout<<H_max<<"\n";

#ifdef DEBUG

  debugFile << "H_max:" << H_max << "(" << i_max << "," << j_max << ")" << "\n";

#endif



     // Backtracking from H_max

  int current_i=i_max,current_j=j_max;

  int next_i=I_i[current_i][current_j];

  int next_j=I_j[current_i][current_j];

//  int tick=0;

  //char consensus_a[N_a+N_b+2],consensus_b[N_a+N_b+2];

  string consensus_a(""), consensus_b(""), match("");



  newAlign.match_align = newAlign.target_align = newAlign.query_align = "";

  int num_of_match = 0;



  //while(next_j > 0 && next_i > 0){

  while(current_j > 0 && current_i > 0){

//	  debugFile << "current_i: " << current_i << "; current_j: " << current_j << "\n";



    if(next_i==current_i)

	{

		consensus_a += '-'; //consensus_a[tick] = '-';                  // deletion in A

		match += ' ';

		consensus_b += seq_b[current_j-1];	//b must be some actual char, cannot be '-' aigns with '-'!

	}

    else

	{

		consensus_a += seq_a[current_i-1]; //consensus_a[tick] = seq_a[current_i-1];   // match/mismatch in A

	    if(next_j==current_j)

		{

			consensus_b += '-'; //consensus_b[tick] = '-';                  // deletion in B

			match += ' ';

		}

		else

		{

			consensus_b += seq_b[current_j-1]; //consensus_b[tick] = seq_b[current_j-1];   // match/mismatch in B

			if (seq_a[current_i-1] == seq_b[current_j-1])

			{

				match += seq_a[current_i-1];

				num_of_match++;

			}

			else

			{

				match += ' ';

			}

		}

	}



    current_i = next_i;

    current_j = next_j;

    next_i = I_i[current_i][current_j];

    next_j = I_j[current_i][current_j];

//    tick++;

    }



	while (current_i > 0)

	{

//		debugFile << "current_i: " << current_i  << "\n";



		consensus_a += seq_a[current_i-1];

		consensus_b += '-';

		match += ' ';



		current_i--;

	}



	while (current_j > 0)

	{

//		debugFile << "current_j: " << current_j << "\n";



		consensus_a += '-';

		consensus_b += seq_b[current_j-1];

		match += ' ';



		current_j--;

	}

			

/*  for (i=match.length()-1; i>=0; i--) //reverse strings

  {

	  newAlign.query_align += consensus_a[i];

	  newAlign.match_align += match[i];

	  newAlign.target_align += consensus_b[i];

  }

*/

  reverse(consensus_a.begin(), consensus_a.end());

  reverse(consensus_b.begin(), consensus_b.end());

  reverse(match.begin(), match.end());

  newAlign.query_align += consensus_a;

  newAlign.match_align += match;

  newAlign.target_align += consensus_b;



  align_pid = (float)num_of_match / newAlign.match_align.length();



 // Output of the consensus motif to the console

/*  cout<<"\n"<<"***********************************************"<<"\n";

  cout<<"The alignment of the sequences"<<"\n"<<"\n";

  for( i=0;i<N_a;i++){cout<<seq_a[i];}; cout<<"  and"<<"\n";

  for( i=0;i<N_b;i++){cout<<seq_b[i];}; cout<<"\n"<<"\n";

  cout<<"is for the parameters  mu = "<<mu<<" and delta = "<<delta<<" given by"<<"\n"<<"\n";  

  for( i=tick-1;i>=0;i--) cout<<consensus_a[i]; 

  cout<<"\n";

  for( j=tick-1;j>=0;j--) cout<<consensus_b[j];

  cout<<"\n";

*/

#ifdef DEBUG

  debugFile <<"The global alignment of the sequences"<<"\n";

	debugFile << seq_a << "\n";

	debugFile << seq_b << "\n";

	debugFile << newAlign << "\n";

#endif



  //clean up memory

  for(i=0;i<=N_a;i++){

	  delete [] H[i];

	  delete [] I_i[i];

	  delete [] I_j[i];

  } 

  delete [] H;

  delete [] I_i;

  delete [] I_j;



  return H_max;

} // END of main


int GetGlobalAlignment_scorematrix(string& seq_a, string& seq_b, Input_Alignment& newAlign, float& align_pid, ofstream& debugFile)
{
  int N_a = seq_a.length();                     // get the actual lengths of the sequences
  int N_b = seq_b.length();

//For long sequences when we don't have enough space to hold the H arrays, we'll need to use 
//Hirschberg's algorithm (a modification to Needleman-Wunsch). 
//It also returns optimal global alignment, but with only O(m+n) space and 
//at the cost of roughly twice as much time as Needleman-Wunsch?
  if (N_a * N_b > MEMORY_LIMIT_SW) //80M, this roughly translates to 80M*4(sizeof(int))*3(num_of_arrays) ~ 960M memory)
	  return (ComputeEditScript(const_cast<char*>(seq_a.c_str()), N_a, const_cast<char*>(seq_b.c_str()), N_b, 
		newAlign, align_pid));
  ////////////////////////////////////////////////
  // initialize H
  int i, j;
  int** H;
  int **I_i, **I_j;
  H = new int*[N_a+1];
  I_i = new int*[N_a+1];
  I_j = new int*[N_a+1];
  bool **opengap; //use this to keep track whether a gap is opengap
  opengap = new bool*[N_a+1];

  //initialize first row and first column (score matrix based gap penalty)
  bool init_gap = false;
  H[0][0]=0; 
  for (i=1;i<=N_a;i++)
	  H[i][0] = similarity_score('-', 'x', init_gap, init_gap);
  init_gap = false;
  for (j=1;j<=N_b;j++)
	  H[0][j] = similarity_score('-', 'x', init_gap, init_gap);

  for(i=0;i<=N_a;i++){
	  H[i] = new int[N_b+1];
	  I_i[i] = new int[N_b+1];
	  I_j[i] = new int[N_b+1];
	  //H[i][0]=0; //only need to initialize first row and first column
	  opengap[i] = new bool[N_b+1];
	  opengap[i][0] = false;
  } 
  for (j=0;j<=N_b;j++){
	  //H[0][j]=0; //initialize first column
	  opengap[0][j] = false;
  }

  int temp[3];
  bool temp_opengap[3];

  for(i=1;i<=N_a;i++){
    for(j=1;j<=N_b;j++){
	  /*temp[0] = H[i-1][j-1]+sim_score(seq_a[i-1],seq_b[j-1]);
      temp[1] = H[i-1][j];
      temp[2] = H[i][j-1];*/
	  AssignScores_scorematrix<int>(temp, H[i-1][j-1], H[i-1][j], H[i][j-1], seq_a[i-1], seq_b[j-1], 
		  opengap[i-1][j-1], opengap[i-1][j], opengap[i][j-1], temp_opengap);
      H[i][j] = find_array_max<int>(temp,3);
      switch(ind){
      case 0:                                  // score in (i,j) stems from a match/mismatch
   		I_i[i][j] = i-1;
		I_j[i][j] = j-1;
		opengap[i][j] = temp_opengap[0];
		break;
      case 1:                                  // score in (i,j) stems from a deletion in sequence A
     	I_i[i][j] = i-1;
		I_j[i][j] = j;
		opengap[i][j] = temp_opengap[1];
		break;
      case 2:                                  // score in (i,j) stems from a deletion in sequence B
      	I_i[i][j] = i;
		I_j[i][j] = j-1;
		opengap[i][j] = temp_opengap[2];
		break;
      }
    }
  }

  // search H for the maximal score
  int H_max = H[N_a][N_b];
  int i_max=N_a,j_max=N_b;

#ifdef DEBUG
  debugFile << "H_max:" << H_max << "(" << i_max << "," << j_max << ")" << "\n";
#endif

  // Backtracking from H_max
  int current_i=i_max,current_j=j_max;
  int next_i=I_i[current_i][current_j];
  int next_j=I_j[current_i][current_j];

  //char consensus_a[N_a+N_b+2],consensus_b[N_a+N_b+2];
  string consensus_a(""), consensus_b(""), match("");

  newAlign.match_align = newAlign.target_align = newAlign.query_align = "";
  int num_of_match = 0;

  //while(next_j > 0 && next_i > 0){
  while(current_j > 0 && current_i > 0){
//	  debugFile << "current_i: " << current_i << "; current_j: " << current_j << "\n";
    if(next_i==current_i)
	{
		consensus_a += '-'; //consensus_a[tick] = '-';                  // deletion in A
		match += ' ';
		consensus_b += seq_b[current_j-1];	//b must be some actual char, cannot be '-' aigns with '-'!
	}
    else
	{
		consensus_a += seq_a[current_i-1]; //consensus_a[tick] = seq_a[current_i-1];   // match/mismatch in A
	    if(next_j==current_j)
		{
			consensus_b += '-'; //consensus_b[tick] = '-';                  // deletion in B
			match += ' ';
		}
		else
		{
			consensus_b += seq_b[current_j-1]; //consensus_b[tick] = seq_b[current_j-1];   // match/mismatch in B
			if (seq_a[current_i-1] == seq_b[current_j-1])
			{
				match += seq_a[current_i-1];
				num_of_match++;
			}
			else
			{
				match += ' ';
			}
		}
	}

    current_i = next_i;
    current_j = next_j;
    next_i = I_i[current_i][current_j];
    next_j = I_j[current_i][current_j];
//    tick++;
  }

  while (current_i > 0)
  {
//		debugFile << "current_i: " << current_i  << "\n";
		consensus_a += seq_a[current_i-1];
		consensus_b += '-';
		match += ' ';

		current_i--;
  }

  while (current_j > 0)
  {
//		debugFile << "current_j: " << current_j << "\n";
		consensus_a += '-';
		consensus_b += seq_b[current_j-1];
		match += ' ';

		current_j--;
  }
			
/*  for (i=match.length()-1; i>=0; i--) //reverse strings
  {
	  newAlign.query_align += consensus_a[i];
	  newAlign.match_align += match[i];
	  newAlign.target_align += consensus_b[i];
  }
*/

  reverse(consensus_a.begin(), consensus_a.end());
  reverse(consensus_b.begin(), consensus_b.end());
  reverse(match.begin(), match.end());
  newAlign.query_align += consensus_a;
  newAlign.match_align += match;
  newAlign.target_align += consensus_b;

  align_pid = (float)num_of_match / newAlign.match_align.length();

  // Output of the consensus motif to the console
/*  cout<<"\n"<<"***********************************************"<<"\n";
  cout<<"The alignment of the sequences"<<"\n"<<"\n";
  for( i=0;i<N_a;i++){cout<<seq_a[i];}; cout<<"  and"<<"\n";
  for( i=0;i<N_b;i++){cout<<seq_b[i];}; cout<<"\n"<<"\n";
  cout<<"is for the parameters  mu = "<<mu<<" and delta = "<<delta<<" given by"<<"\n"<<"\n";  
  for( i=tick-1;i>=0;i--) cout<<consensus_a[i]; 
  cout<<"\n";
  for( j=tick-1;j>=0;j--) cout<<consensus_b[j];
  cout<<"\n";
*/

#ifdef DEBUG
  debugFile <<"The global alignment of the sequences"<<"\n";
  debugFile << seq_a << "\n";
  debugFile << seq_b << "\n";
  debugFile << newAlign << "\n";
#endif

  //clean up memory
  for(i=0;i<=N_a;i++){
	  delete [] H[i];
	  delete [] I_i[i];
	  delete [] I_j[i];
	  delete [] opengap[i];
  } 
 
  delete [] H;
  delete [] I_i;
  delete [] I_j;
  delete [] opengap;

  return H_max;
} // END of main

int GetGlobalAlignment_PID_CompScore(string& seq_a, string& seq_b, Input_Alignment& newAlign, float& align_pid, ofstream& debugFile)
{
//	cout << "QUICK DIRT CHECK" << "\n";

  int N_a = seq_a.length();                     // get the actual lengths of the sequences
  int N_b = seq_b.length();

//For long sequences when we don't have enough space to hold the H arrays, we'll need to use 
//Hirschberg's algorithm (a modification to Needleman-Wunsch). 
//It also returns optimal global alignment, but with only O(m+n) space and 
//at the cost of roughly twice as much time as Needleman-Wunsch?
  if (N_a * N_b > MEMORY_LIMIT_SW) //80M, this roughly translates to 80M*4(sizeof(int))*3(num_of_arrays) ~ 960M memory)
	  return (ComputeEditScript(const_cast<char*>(seq_a.c_str()), N_a, const_cast<char*>(seq_b.c_str()), N_b, 
		newAlign, align_pid));

  ////////////////////////////////////////////////
  // initialize H
  //double H[N_a+1][N_b+1];     
  int i, j;
  int** H;
  int **I_i, **I_j;
  H = new int*[N_a+1];
  I_i = new int*[N_a+1];
  I_j = new int*[N_a+1];

  for(i=0;i<=N_a;i++){
	  H[i] = new int[N_b+1];
	  I_i[i] = new int[N_b+1];
	  I_j[i] = new int[N_b+1];

	  H[i][0]=0; //only need to initialize first row and first column
//    for(j=0;j<=N_b;j++){
//      H[i][j]=0.;
//    }
  } 
  for (j=0;j<=N_b;j++){
	  H[0][j]=0; //initialize first column
  }

  int temp[3];
  //int I_i[N_a+1][N_b+1],I_j[N_a+1][N_b+1];     // Index matrices to remember the 'path' for backtracking
  
  // here comes the actual algorithm

  for(i=1;i<=N_a;i++){
    for(j=1;j<=N_b;j++){
      //temp[0] = H[i-1][j-1]+similarity_score(seq_a[i-1],seq_b[j-1]); 
      //temp[1] = H[i-1][j]-delta;                  
      //temp[2] = H[i][j-1]-delta;
	  //cout << "i:" << i << ";j:" << j << "\n";
	  temp[0] = H[i-1][j-1]+sim_score(seq_a[i-1],seq_b[j-1]);
      temp[1] = H[i-1][j];
      temp[2] = H[i][j-1];
      H[i][j] = find_array_max<int>(temp,3);
      switch(ind){
      case 0:                                  // score in (i,j) stems from a match/mismatch
   		I_i[i][j] = i-1;
		I_j[i][j] = j-1;
		break;
      case 1:                                  // score in (i,j) stems from a deletion in sequence A
     	I_i[i][j] = i-1;
		I_j[i][j] = j;
		break;
      case 2:                                  // score in (i,j) stems from a deletion in sequence B
      	I_i[i][j] = i;
		I_j[i][j] = j-1;
		break;
      }
    }
  }
  
  // search H for the maximal score
  int H_max = H[N_a][N_b];
  int i_max=N_a,j_max=N_b;

  //what if H_max is 0?
//  if (H_max == 0)
//  {
//	  debugFile << "global H_max is 0, no alignment!" << "\n";
//	  return 0.;
//  }
  
#ifdef DEBUG
  debugFile << "H_max:" << H_max << "(" << i_max << "," << j_max << ")" << "\n";
#endif

   // Backtracking from H_max
  int current_i=i_max,current_j=j_max;
  int next_i=I_i[current_i][current_j];
  int next_j=I_j[current_i][current_j];

  //char consensus_a[N_a+N_b+2],consensus_b[N_a+N_b+2];
  string consensus_a(""), consensus_b(""), match("");
  newAlign.match_align = newAlign.target_align = newAlign.query_align = "";
  int num_of_match = 0;

  bool opengap = false;
  int score = 0;
  //while(next_j > 0 && next_i > 0){
  while(current_j > 0 && current_i > 0){
    if(next_i==current_i)
	{
		consensus_a += '-'; //consensus_a[tick] = '-';                  // deletion in A
		match += ' ';
		consensus_b += seq_b[current_j-1];	//b must be some actual char, cannot be '-' aigns with '-'!
		score += similarity_score('-', seq_b[current_j-1], opengap, opengap);
	}
    else
	{
		consensus_a += seq_a[current_i-1]; //consensus_a[tick] = seq_a[current_i-1];   // match/mismatch in A
	    if(next_j==current_j)
		{
			consensus_b += '-'; //consensus_b[tick] = '-';                  // deletion in B
			match += ' ';
			score += similarity_score(seq_a[current_i-1], '-', opengap, opengap);
		}
		else
		{
			consensus_b += seq_b[current_j-1]; //consensus_b[tick] = seq_b[current_j-1];   // match/mismatch in B
			score += similarity_score(seq_a[current_i-1], seq_b[current_j-1], opengap, opengap);
			if (seq_a[current_i-1] == seq_b[current_j-1])
			{
				match += seq_a[current_i-1];
				num_of_match++;
			}
			else
			{
				match += ' ';
			}
		}
	}	

    current_i = next_i;
    current_j = next_j;
    next_i = I_i[current_i][current_j];
    next_j = I_j[current_i][current_j];
    }

	while (current_i > 0)
	{
//		debugFile << "current_i: " << current_i  << "\n";
		consensus_a += seq_a[current_i-1];
		consensus_b += '-';
		match += ' ';
		score += similarity_score(seq_a[current_i-1], '-', opengap, opengap);

		current_i--;
	}

	while (current_j > 0)
	{
//		debugFile << "current_j: " << current_j << "\n";
		consensus_a += '-';
		consensus_b += seq_b[current_j-1];
		match += ' ';
		score += similarity_score('-', seq_b[current_j-1], opengap, opengap);

		current_j--;
	}
			
  reverse(consensus_a.begin(), consensus_a.end());
  reverse(consensus_b.begin(), consensus_b.end());
  reverse(match.begin(), match.end());
  newAlign.query_align += consensus_a;
  newAlign.match_align += match;
  newAlign.target_align += consensus_b;

  align_pid = (float)num_of_match / newAlign.match_align.length();

#ifdef DEBUG
  debugFile <<"The global alignment of the sequences"<<"\n";
	debugFile << seq_a << "\n";
	debugFile << seq_b << "\n";
	debugFile << newAlign << "\n";
	debugFile << "align_pid:" << align_pid << "(" << num_of_match << "/" << newAlign.match_align.length() 
		<< ");align_score:" << score << "\n";
#endif

  //clean up memory
	for(i=0;i<=N_a;i++){
	  delete [] H[i];
	  delete [] I_i[i];
	  delete [] I_j[i];
  } 
  delete [] H;
  delete [] I_i;
  delete [] I_j;

  //return H_max;
  return score;
} // END of main

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

bool TryExtAlign(char cur_query_char, char cur_target_char, 

				 double& cur_score, double& cur_max_score, bool& opengap)

{

	//now try to align cur_query_char with cur_target_char

	double temp[3];

	bool temp_opengap[3];



	temp[0] = similarity_score(cur_query_char, cur_target_char, opengap, temp_opengap[0]);

	temp[1] = similarity_score('-', cur_target_char, opengap, temp_opengap[1]);

	temp[2] = similarity_score(cur_query_char, '-', opengap, temp_opengap[2]);

	cur_score += find_array_max<double>(temp, 3);

	if (cur_score > cur_max_score)

		cur_max_score = cur_score;

	else

		if (cur_score < cur_max_score - REPAIR_HSP_EXTEND_SCORE_DROP)

			return false;



	opengap = temp_opengap[ind];

	return true;

}



//extend curHSP to both directions if possible, stop extending when the score of alignment is 

//lower than ("max_score" - REPAIR_HSP_EXTEND_SCORE_DROP), "max_score" is the running max score

void ExtendAlignment(double max_score, string& query_seq, string& target_seq_left, string& target_seq_right, 
				//pair<int, char*>& chr_seq, 
				vector<string>& chr_seq,
				bool isPosStrand, 

				HSP_Gene_Pair* curHSP, Input_Alignment& curAlign, ofstream& debugFile, 

				int target_left_limit, int target_right_limit)

{

	//cout << "query_seq:" << query_seq << ";\n";

	//cout << "target_seq_left:" << target_seq_left << ";\ntarget_seq_right:" << target_seq_right << ";\n";

#ifdef DEBUG

	debugFile << "Extend alignment: " << "\n";

	debugFile << query_seq << "\n";

	debugFile << "target_seq_left:" << target_seq_left << ";\ntarget_seq_right:" << target_seq_right << ";\n";

#endif





	int retrieve_length = 600;//retrieve the next 600 bp, must be multiple of 3



	int left_length = target_seq_left.length();



	int num_of_match = (int)(curAlign.match_align.length() * (curHSP->pid/100));



	double cur_max_score = max_score; //keep the running max score

	double cur_score = max_score; //keep the current score



	int target_total_len = chromosome_start_pos-1 + 
		//chr_seq.first; 
		LenOfStrVec(chr_seq);

#ifdef DEBUG
	debugFile << "target_left_limit:" << target_left_limit << "; target_right_limit:" << target_right_limit 
		<< "; target_total_len:" << target_total_len << "\n";
#endif

	if (isPosStrand)
	{
		if (target_left_limit < chromosome_start_pos)
			target_left_limit = chromosome_start_pos;
		if (target_right_limit > target_total_len)
			target_right_limit = target_total_len;
	}
	else
	{
		if (target_right_limit < chromosome_start_pos)
			target_right_limit = chromosome_start_pos;
		if (target_left_limit > target_total_len)
			target_left_limit = target_total_len;
	}

#ifdef DEBUG
	debugFile << "target_left_limit:" << target_left_limit << "; target_right_limit:" << target_right_limit 
		<< "; target_total_len:" << target_total_len << "\n";
#endif


	//extend to left

	int i, j;

	char cur_query_char, cur_target_char;

	bool opengap = false; //initial alignment must have no gap at the end

	Input_Alignment align_left; //for the left side of alignment extension, store in the reverse order first

	string tempStr;

	if (isPosStrand)

	{

		i = curHSP->gene_start - 1;

		j = curHSP->HSP_start - 1;
		
		int target_left_limit_in_frame = target_left_limit + ((curHSP->HSP_start - target_left_limit) % 3);

		//while (i > 0 && j > 0)
		//while (i > 0 && j > target_left_limit) //MODIFIED: make sure the left/right_limit is in frame so the extended alignment won't go over the original limit
		while (i > 0 && j > target_left_limit_in_frame)
		{			

			//string resultStr;

			//GetSubstrFromVecStrs(chr_seq, isPosStrand, j-1, 1, resultStr);

			if (left_length == 0)

			{

				if (j < chromosome_start_pos - 1 + retrieve_length)

					left_length = (j-chromosome_start_pos+1) - (j-chromosome_start_pos+1)%3; 

				else

					left_length = retrieve_length;
				if (left_length < 3)
					break;

				GetSubstrFromVecStrs(chr_seq, isPosStrand, j-left_length, left_length, tempStr);

				StrToLower(tempStr);

				DNA2AA(tempStr, 0, target_seq_left);

#ifdef DEBUG
				debugFile << "j:" << j << ";left_length:" << left_length << "\n";

					debugFile << "target_seq(tempStr):" << "\n" << tempStr << "\n";

					debugFile << target_seq_left << "\n";

#endif

				left_length = left_length/3;

				//cout << "target_seq_left:" << target_seq_left << ";\n";

			}

			cur_query_char = query_seq[i-1];

			cur_target_char = target_seq_left[left_length-1];

			//now try to align cur_query_char with cur_target_char

			if (!TryExtAlign(cur_query_char, cur_target_char, cur_score, cur_max_score, opengap))

				break;



			switch (ind) {

			case 0:

				curHSP->gene_start--;

				curHSP->HSP_start-=3;

				if (cur_query_char == cur_target_char)

				{

					num_of_match++;

					align_left.match_align += cur_query_char;

				}

				else

				{

					align_left.match_align += ' ';

				}

				align_left.query_align += cur_query_char;

				align_left.target_align += cur_target_char;				

				i--;

				j-=3;

				left_length--;

				break;

			case 1:

				curHSP->HSP_start-=3;

				align_left.query_align += '-';

				align_left.match_align += ' ';

				align_left.target_align += cur_target_char;

				j-=3;

				left_length--;

				break;

			case 2:

				curHSP->gene_start--;

				align_left.query_align += cur_query_char;

				align_left.match_align += ' ';

				align_left.target_align += '-';

				i--;

				break;

			}

		}

	}

	else //negative strand

	{

		i = curHSP->gene_start - 1;

		j = curHSP->HSP_end + 1;
		int target_left_limit_in_frame = target_left_limit - ((target_left_limit - curHSP->HSP_end) % 3);

		//while (i > 0 && j <= target_total_len)

		while (i > 0 && j < target_left_limit_in_frame)

		{

			//string resultStr;

			//GetSubstrFromVecStrs(chr_seq, isPosStrand, j-1, 1, resultStr);

			if (left_length == 0)

			{
				int remain_len = target_total_len - j + 1;

				if (remain_len < retrieve_length)

					left_length = remain_len - remain_len%3; 

				else

					left_length = retrieve_length;

				if (left_length < 3) //check left_length to make sure it's non zero
					break;

				GetSubstrFromVecStrs_NegRev(chr_seq, isPosStrand, j-1, left_length, tempStr);

				StrToLower(tempStr);

				DNA2AA(tempStr, 0, target_seq_left);

#ifdef DEBUG

					debugFile << "target_seq(tempStr):" << "\n" << tempStr << "\n";

					debugFile << target_seq_left << "\n";

#endif

				left_length = left_length/3;

				//cout << "target_seq_left:" << target_seq_left << ";\n";

			}

			cur_query_char = query_seq[i-1];

			cur_target_char = target_seq_left[left_length-1];

			if (!TryExtAlign(cur_query_char, cur_target_char, cur_score, cur_max_score, opengap))

				break;



			switch (ind) {

			case 0:

				curHSP->gene_start--;

				curHSP->HSP_end+=3;

				if (cur_query_char == cur_target_char)

				{

					num_of_match++;

					align_left.match_align += cur_query_char;

				}

				else

				{

					align_left.match_align += ' ';

				}

				align_left.query_align += cur_query_char;

				align_left.target_align += cur_target_char;				

				i--;

				j+=3;

				left_length--;

				break;

			case 1:

				curHSP->HSP_end+=3;

				align_left.query_align += '-';

				align_left.match_align += ' ';

				align_left.target_align += cur_target_char;

				j+=3;

				left_length--;

				break;

			case 2:

				curHSP->gene_start--;

				align_left.query_align += cur_query_char;

				align_left.match_align += ' ';

				align_left.target_align += '-';

				i--;

				break;

			}

		}



	}

	//now reverse align_left and add to curAlign

	reverse(align_left.query_align.begin(), align_left.query_align.end());

	reverse(align_left.match_align.begin(), align_left.match_align.end());

	reverse(align_left.target_align.begin(), align_left.target_align.end());

	curAlign.query_align = align_left.query_align + curAlign.query_align;

	curAlign.match_align = align_left.match_align + curAlign.match_align;

	curAlign.target_align = align_left.target_align + curAlign.target_align;



	//now do the right side

	//reset stuff

	cur_max_score = max_score; //keep the running max score

	cur_score = max_score; //keep the current score

	opengap = false;

	left_length = target_seq_right.length();



	int query_total_len = query_seq.length();

	if (isPosStrand)

	{

		i = curHSP->gene_end + 1;

		j = curHSP->HSP_end + 1;

		int target_right_limit_in_frame = target_right_limit - ((target_right_limit - curHSP->HSP_end)%3);

		int k = 0;

		//while (i <= query_total_len && j <= target_total_len)

		while (i <= query_total_len && j < target_right_limit_in_frame)

		{

			if (left_length == 0)

			{
				int remain_len = target_total_len - j + 1;

				if (remain_len < retrieve_length)

					left_length = remain_len - remain_len%3; 

				else

					left_length = retrieve_length;
				if (left_length < 3)
					break;

				GetSubstrFromVecStrs(chr_seq, isPosStrand, j-1, left_length, tempStr);

				StrToLower(tempStr);

				DNA2AA(tempStr, 0, target_seq_right);

#ifdef DEBUG

					debugFile << "target_seq(tempStr):" << "\n" << tempStr << "\n";

					debugFile << target_seq_right << "\n";

#endif

				k = 0;

				left_length = left_length/3;

				//cout << "target_seq_right:" << target_seq_right << ";\n";

			}

			cur_query_char = query_seq[i-1];

			cur_target_char = target_seq_right[k];

			//now try to align cur_query_char with cur_target_char

			if (!TryExtAlign(cur_query_char, cur_target_char, cur_score, cur_max_score, opengap))

				break;



			switch (ind) {

			case 0:

				curHSP->gene_end++;

				curHSP->HSP_end+=3;

				if (cur_query_char == cur_target_char)

				{

					num_of_match++;

					curAlign.match_align += cur_query_char;

				}

				else

				{

					curAlign.match_align += ' ';

				}

				curAlign.query_align += cur_query_char;

				curAlign.target_align += cur_target_char;				

				i++;

				j+=3;

				left_length--;

				k++;

				break;

			case 1:

				curHSP->HSP_end+=3;

				curAlign.query_align += '-';

				curAlign.match_align += ' ';

				curAlign.target_align += cur_target_char;

				j+=3;

				left_length--;

				k++;

				break;

			case 2:

				curHSP->gene_end++;

				curAlign.query_align += cur_query_char;

				curAlign.match_align += ' ';

				curAlign.target_align += '-';

				i++;

				break;

			}

		}



	}

	else //negative strand

	{

		i = curHSP->gene_end + 1;

		j = curHSP->HSP_start - 1;
		int target_right_limit_in_frame = target_right_limit + ((curHSP->HSP_start-target_right_limit)%3);

		int k = 0;

		//while (i <= query_total_len && j > 0)

		while (i <= query_total_len && j > target_right_limit_in_frame)

		{

			if (left_length == 0)

			{

				if ( j < chromosome_start_pos - 1 + retrieve_length)

					left_length = (j-chromosome_start_pos+1) - (j-chromosome_start_pos+1)%3; 

				else

					left_length = retrieve_length;
				if (left_length < 3)
					break;

				GetSubstrFromVecStrs_NegRev(chr_seq, isPosStrand, j-left_length, left_length, tempStr);

				StrToLower(tempStr);

				DNA2AA(tempStr, 0, target_seq_right);

#ifdef DEBUG

					debugFile << "target_seq(tempStr):" << "\n" << tempStr << "\n";

					debugFile << target_seq_right << "\n";

#endif

				k = 0;

				left_length = left_length/3;

				//cout << "target_seq_right:" << target_seq_right << ";\n";

			}

			cur_query_char = query_seq[i-1];

			cur_target_char = target_seq_right[k];

			//now try to align cur_query_char with cur_target_char

			if (!TryExtAlign(cur_query_char, cur_target_char, cur_score, cur_max_score, opengap))

				break;



			switch (ind) {

			case 0:

				curHSP->gene_end++;

				curHSP->HSP_start-=3;

				if (cur_query_char == cur_target_char)

				{

					num_of_match++;

					curAlign.match_align += cur_query_char;

				}

				else

				{

					curAlign.match_align += ' ';

				}

				curAlign.query_align += cur_query_char;

				curAlign.target_align += cur_target_char;				

				i++;

				j-=3;

				left_length--;

				k++;

				break;

			case 1:

				curHSP->HSP_start-=3;

				curAlign.query_align += '-';

				curAlign.match_align += ' ';

				curAlign.target_align += cur_target_char;

				j-=3;

				left_length--;

				k++;

				break;

			case 2:

				curHSP->gene_end++;

				curAlign.query_align += cur_query_char;

				curAlign.match_align += ' ';

				curAlign.target_align += '-';

				i++;

				break;

			}

		}



	}



	curHSP->pid = ((float) num_of_match / curAlign.match_align.length()) * 100;

}



//not used

/* This string reversal might be faster and more efficient -- uses minimal 

memory great for embedded systems BY: STEVE KING Tucson Arizona */ 

//Da... only works for char*, what about string?

void ReverseString(char *ptr ) 

{ 

	char saved_char; // need to save off a character 

	char *saved_ptr; // need to save off original ptr address 

	saved_ptr = ptr; // save the actual ptr address 

	// I only need to go through half the string to swap positions 

	// so when the ratio of savedADDRESS to incrementedADDRESS is 2, then we quit 

	while(strlen(saved_ptr)/strlen(ptr) != 2) 

	{ 

		saved_char = *(saved_ptr+strlen(ptr)-1); // save off a char 

		*(saved_ptr+strlen(ptr)-1) = *ptr; // swap chars 

		*ptr = saved_char; // swap chars 

		ptr++; // increment to next address 

	}

	ptr = saved_ptr; // need to return the ptr address to original

} 





