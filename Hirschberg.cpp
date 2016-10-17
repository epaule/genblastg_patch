/* Hirschberg's algorithm a non-recursive version.
 *
 * vamsik@engr	11/29/2008
 **/
//#include<stdio.h>
//#include<string.h>
//#include<stdlib.h>
//#include<assert.h>
#include "Hirschberg.h"

#include <queue>
using namespace std;

#define INS 1
#define DEL 1 
#define CHG 1

/*Space Used by the Entire Algorithm*/
unsigned long *D[2] = {NULL,NULL};
unsigned long *Dr[2] = {NULL,NULL};
/*We find the edit script relative to S2. So */
char *EditScriptS2 = NULL; /*Store Deletes and Changes in this*/
char *EditScriptS1 = NULL; /*Stores Inserts*/

void CreateWorkSpace(unsigned long s1len,unsigned long s2len){
	unsigned long i;
	for(i=0;i<2;i++){
		D[i] = new unsigned long[s2len+1];
		Dr[i] = new unsigned long[s2len+1];
	}
	/*Apply to the Second String*/
	EditScriptS2 = new char[s2len];
	EditScriptS1 = new char[s1len];
}

void FreeWorkSpace(){
	unsigned long i;
	for(i=0;i<2;i++){
		delete [] D[i];
		delete [] Dr[i];
	}
	delete [] EditScriptS1;
	delete [] EditScriptS2;
}

unsigned long min3(unsigned long a,unsigned long b,unsigned long c){
	return (a>b)?((b>c)?c:b):(a>c)?c:a;
}

inline unsigned long CELING(unsigned long a,unsigned long b){
	return (a%b)?(a/b)+1:(a/b);
}

inline char GetOp(unsigned int i,unsigned int j,char c1,char c2){
	char op = 'I';
	if(i){ /*Figure out D or C*/
		op = (D[i-1][j-1]+((c1 == c2)?0:CHG)== D[i][j])?'C':
			(D[i-1][j]+DEL == D[i][j])?'D':'I';
	}
	return op;
}

int ComputeEditScript(char *s1, unsigned long n1,char *s2,unsigned long n2, Input_Alignment& align, float& pid){
	
	CreateWorkSpace(n1, n2);

	unsigned long q;
	char *sub_s1,*sub_s2;
	char op;
	unsigned long sub_n1,sub_n2,i,j;
	queue<ESSubProblem> bfs_list; 
	ESSubProblem es_sub;

	/*Put the toplevel subproblem into the CQueue*/
	es_sub.sub_s1 = s1; es_sub.sub_n1 = n1;
	es_sub.sub_s2 = s2; es_sub.sub_n2 = n2;
	bfs_list.push(es_sub);
	
	while( !bfs_list.empty() ){
		es_sub = bfs_list.front();
		bfs_list.pop();

		sub_s1 = es_sub.sub_s1; sub_s2 = es_sub.sub_s2;
		sub_n1 = es_sub.sub_n1; sub_n2 = es_sub.sub_n2;

		if(sub_n1 <=1){ 
			/*If sub_n1 is 1 or 0 we cannot split it any more*/
			EditDistanceLS(sub_s1,sub_n1,sub_s2,sub_n2);
			i = sub_n1; j = sub_n2;
			/*Figure out the edit operations*/
			while(i || j){
				op = GetOp(i,j,(i?sub_s1[i-1]:'\0'),
					(j?sub_s2[j-1]:'\0'));
				switch(op){
					case 'I':
							EditScriptS2[(&sub_s2[j-1]) - s2] = 'I';
							j--;
							break;
					case 'D':
							EditScriptS1[(&sub_s1[i-1]) - s1] = 'D';
							i--;
							break;
					case 'C':
							EditScriptS1[(&sub_s1[i-1]) - s1] =
								(sub_s1[i-1] == sub_s2[j-1])?':':'C';
								i--; j--;
				}
			}
		}else {
			/*Add the two subproblems*/
			q = FindMinSplit(sub_s1,sub_n1,sub_s2,sub_n2);
			/*SUB PROBLEM 1*/
			es_sub.sub_n1 = CELING(sub_n1,2);
			es_sub.sub_n2 = q;
			bfs_list.push(es_sub);

			/*SUB PROBLEM 2*/
			es_sub.sub_s1 = &(sub_s1[CELING(sub_n1,2)]);
			es_sub.sub_s2 = &(sub_s2[q]);
			es_sub.sub_n1 = sub_n1/2;
			es_sub.sub_n2 = sub_n2-q;
			bfs_list.push(es_sub);
		}
	}

	int score = PrintAlignment(s1, n1, s2, n2, EditScriptS2,EditScriptS1, align, pid);
	
	FreeWorkSpace();

	return score;
}

/*Given 2 strings S1[0......n1-1] and S2[0.......n2-1], Find
 *a 1<=q<=n2 such that D[floor(n1/2)][q]+Dr[ceil(n1/2)][n2-q]
 *( S2[0......q-1] )            (S2[q.....n2-1])
 *( S1[0......floor(n1/2)-1] ) (S1[floor(n1/2)....n1-1])
 */
unsigned long FindMinSplit(char *s1,unsigned long n1,char *s2, unsigned long n2){
	unsigned long i,j,chg,chgr;
	unsigned long emin,qmin,n1_split,q;
	/*If S1=(null) then INS, If S2=(null) then DEL. 
	* If S1!=S2 (both single chars) then CHG*/
	for(j=0;j<=n2;j++){
		D[0][j] = j*(INS); 
		Dr[0][j] = j*(INS);
	} 

	n1_split = CELING(n1,2);
	for(i=1;i<=n1_split;i++){
		D[i%2][0] = i*DEL; Dr[i%2][0] = i*DEL;
		for(j=1;j<=n2;j++){
			chg = (s1[i-1]==s2[j-1])?0:CHG;
			chgr = (s1[n1-i] == s2[n2-j])?0:CHG;
			/*Forward DP*/	
			D[i%2][j] = min3(D[(i+1)%2][j-1]+chg,
				D[i%2][j-1]+INS,D[(i+1)%2][j]+DEL);
			/*Reverse DP*/
			Dr[i%2][j] = min3(Dr[(i+1)%2][j-1]+chgr,
				Dr[i%2][j-1]+INS,Dr[(i+1)%2][j]+DEL);
		}
	}
	/*emin = D[(n1_split+(n1)%2)%2][0]+Dr[(n1_split)%2][n2]; */
	emin = D[(n1_split)%2][0] + Dr[(n1_split+(n1)%2)%2][n2];
	qmin = 0;
	for(q=1;q<=n2;q++){
		if((Dr[(n1_split+(n1)%2)%2][n2-q]+D[(n1_split)%2][q]) < emin){
			emin = Dr[(n1_split+(n1)%2)%2][n2-q]+D[(n1_split)%2][q];
			qmin = q;
		}
	}
//#ifdef VERBOSE
	if (VERBOSE)
	printf("The minimum edit distance is %lu \n",emin);
//#endif
	return qmin;
}

/*EditDistanceLS: A Linear Space EditDistance computation.*/
unsigned long EditDistanceLS(char *s1,unsigned long n1,char *s2,unsigned long n2){
	unsigned long i,j,chg;//,chgr;
	/*[Begin] Initialization*/
	for(j=0;j<=n2;j++){
		D[0][j] = j*(INS); 
	}
	/*[End] Initialization*/

	/*Linear Space DP*/
	for(i=1;i<=n1;i++){
		D[i%2][0] = i*DEL; 
		for(j=1;j<=n2;j++){
			chg = (s1[i-1]==s2[j-1])?0:CHG;
			D[i%2][j] = min3(D[(i+1)%2][j-1]+chg,
				D[i%2][j-1]+INS,D[(i+1)%2][j]+DEL);
		}
	}
//#ifdef VERBOSE
	if (VERBOSE)
	printf("The edit distance between S1=(%s) and S2=(%s) is %lu D[n1][n2]\n",
		s1,s2,D[(i+1)%2][n2]);
//#endif
	return D[(i+1)%2][n2];
}

/*The Edit Script will be applied to S2 to make it S1*/
int PrintAlignment(char *s1,unsigned long n1,char *s2,unsigned long n2,char *es2,char *es1, 
					Input_Alignment& align, float& pid){
	unsigned long i=0,j=0;
	int k=0;
//	char buf2[n1+n2+1];
//	char buf1[n1+n2+1];
//	char buf3[n1+n2+1];

	while(i < n1 || j < n2){
		if(es2[j] == 'I'){
//			buf2[k] = s2[j++];
//			buf1[k] = '_';
//			buf3[k++] = ' ';
			align.query_align += '-';
			align.target_align += s2[j++];
			align.match_align += ' ';
		}else if(es1[i] == ':'){
//			buf2[k] = s2[j++];
//			buf1[k] = s1[i++];
//			buf3[k++] = '|';
			align.query_align += s1[i++];
			align.target_align += s2[j++];
			align.match_align += '|';
			k++;
		}else if(es1[i] == 'C'){
//			buf2[k] = s2[j++];
//			buf1[k] = s1[i++];
//			buf3[k++] = ':';
			align.query_align += s1[i++];
			align.target_align += s2[j++];
			align.match_align += ' ';
		}else if(es1[i] == 'D'){
//			buf2[k] = '_';
//			buf1[k] = s1[i++];
//			buf3[k++] = ' ';
			align.query_align += s1[i++];
			align.target_align += '-';
			align.match_align += ' ';
		}else{
			cout << "FATAL: Invalid Hirschberg output...\n";
			exit(-1);
		}
	}

	pid = (float) k / align.match_align.length();
	return k; //the score assumption is based on sim_score()! same as in GetGlobalAlignment() with Needleman algorithm

//	buf1[k] = '\0'; buf2[k] = '\0'; buf3[k] = '\0';
//	printf("%s\n",buf1);
//	printf("%s\n",buf3);
//	printf("%s\n",buf2);
}

/*
int main(int argc,char **argv){
	unsigned long i;
	if(!((s1 = ReadString(&s1len)) && 
		(s2 = ReadString(&s2len)))){
		fprintf(stderr,"WARNING: Unable to Reading Strings Properly\n");
	}
	CreateWorkSpace(s1len,s2len);
#ifdef VERBOSE
	printf("The min split is %lu \n",FindMinSplit(s1,s1len,s2,s2len));
#endif
	printf("The edit distance is %lu\n",EditDistanceLS(s1,s1len,s2,s2len));
	ComputeEditScript(s1,s1len,s2,s2len);
	printf("Printing the Edit Script\n");
	PrintAlignment(s1,s1len,s2,s2len,EditScriptS2,EditScriptS1);
	FreeWorkSpace();
	printf("\n");

}
*/
