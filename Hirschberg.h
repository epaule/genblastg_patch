#ifndef _HIRSCHBERG_H
#define _HIRSCHBERG_H

#include "data_manager.h"

typedef struct _EditScriptSubProblem_{
	char *sub_s1;
	char *sub_s2;
	unsigned long sub_n1;
	unsigned long sub_n2;
}ESSubProblem;

void CreateSpace(unsigned long s1len,unsigned long s2len);
unsigned long min3(unsigned long a,unsigned long b,unsigned long c);
int ComputeEditScript(char *s1,unsigned long n1,char *s2,unsigned long n2, Input_Alignment& align, float& pid);
unsigned long FindMinSplit(char *s1,unsigned long n1,char *s2, unsigned long n2);
unsigned long EditDistanceLS(char *s1,unsigned long n1,char *s2,unsigned long n2);
int PrintAlignment(char *s1,unsigned long n1,char *s2,unsigned long n2,char *es2,char *es1, 
					Input_Alignment& align, float& pid);
#endif
