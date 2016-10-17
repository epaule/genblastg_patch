#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef C45_DEFNS_H
#define C45_DEFNS_H
/*************************************************************************/
/*									 */
/*		Definitions used in C4.5				 */
/*              ------------------------				 */
/*									 */
/*************************************************************************/


#include <stdio.h>
#include <stdlib.h>

#define	 Eof			EOF             /*char read on end of file*/
#define	 Nil			0               /*null pointer*/
#define	 false			0 
#define	 true			1 
#define	 None			-1
#define	 Epsilon                1E-3

#define	 Max(a,b)               ((a)>(b) ? a : b) 
#define	 Min(a,b)               ((a)<(b) ? a : b) 
#define	 Round(x)		((int) (x+0.5))

#define	 Log2			0.69314718055994530942
#define	 Log(x)			((x) <= 0 ? 0.0 : log((float)x) / Log2)

#define	 Bit(b)			(1 << (b))
#define	 In(b,s)		((s[(b) >> 3]) & Bit((b) & 07))
#define	 ClearBits(n,s)		memset(s,0,n)
#define	 CopyBits(n,f,t)	memcpy(t,f,n)
#define	 SetBit(b,s)		(s[(b) >> 3] |= Bit((b) & 07))

#define	 ForEach(v,f,l)		for(v=f ; v<=l ; ++v) 

#define	 Verbosity(d)		if(VERBOSITY >= d)

#define	 Check(v,l,h)\
	     if ( v<l||v>h ) {printf("\t** illegal value **\n"); exit(1);}

/*************************************************************************/
/*									 */
/*		Type definitions for C4.5				 */
/*              -------------------------				 */
/*									 */
/*************************************************************************/


typedef  char	Boolean, *String, *Set;

//typedef  int	ItemNo;		/* data item number */
//typedef  float	ItemCount;	/* count of (partial) items */

typedef  short	ClassNo,	/* class number, 0..MaxClass */
		DiscrValue;	/* discrete attribute value (0 = ?) */
typedef  short	Attribute;	/* attribute number, 0..MaxAtt */

typedef  union  _attribute_value
	 {
	    DiscrValue	_discr_val;
	    float	_cont_val;
	 }
	 	AttValue, *Description;

//CHANGED: since there's no need for any attribute value in Description, 
//just store class in Case[0] (first attribute value of current case)

//#define  CVal(Case,Attribute)   Case[Attribute]._cont_val
//#define  DVal(Case,Attribute)   Case[Attribute]._discr_val
#define  Class(Case)		Case[MaxAtt]._discr_val //Case[MaxAtt+1]._discr_val 


//#define  Unknown  -999		/* unknown value for continuous attribute */


#define  BrDiscr	1	/* node types:	branch */
#define  ThreshContin	2	/*		threshold cut */
#define  BrSubset	3	/*		subset test */

typedef  struct _tree_record *Tree;
typedef  struct _tree_record
	 {
	    short	NodeType;	/* 0=leaf 1=branch 2=cut 3=subset */
	    ClassNo	Leaf;		/* most frequent class at this node */
	    float	Items,		/* no of items at this node */
			*ClassDist,	/* class distribution of items */
	    		Errors;		/* no of errors at this node */
//	    Attribute	Tested; 	/* attribute referenced in test */
	    short	Forks;		/* number of branches at this node */
	    float	Cut,		/* threshold for continuous attribute */
		  	Lower,		/* lower limit of soft threshold */
		  	Upper;		/* upper limit ditto */
	    //Set         *Subset;	/* subsets of discrete values  */
	    Tree	*Branch;	/* Branch[x] = (sub)tree for outcome x */

		//Tree	Parent; //ADDED: the root node's parent is "Nil"
		int		Start; //start position of items covered by this node (index) (inclusive)
		int		End; //end position (inclusive)
	 }
		TreeRec;

//#define  IGNORE		1	/* special attribute status: do not use */
//#define  DISCRETE	2	/* ditto: collect values as data read */


#endif

#ifdef __cplusplus
}
#endif
