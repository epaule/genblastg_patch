#include "c45defns.h"

#ifdef __cplusplus
extern "C" {
#endif

/*************************************************************************/
/*									 */
/*		Global data for C4.5					 */
/*		--------------------					 */
/*									 */
/*************************************************************************/


extern short MaxAtt; /* max att number */
extern short MaxClass;	/* max class number */
extern short MaxDiscrVal;	/* max discrete values for any att */

extern int		MaxItem;	/* max data item number */
extern Description	*Item;		/* data items */
extern float* HSPAlignScores;

//extern DiscrValue	*MaxAttVal;	/* number of values for each att */

//extern String		*ClassName,	/* class names */
//		  	*AttName,	/* att names */
//		  	**AttValName,	/* att value names */
//			FileName;	/* family name of files */

//extern Boolean		AllKnown;	/* true if there have been no splits
//					   on atts with missing values above
//					   the current position in the tree */


/*************************************************************************/
/*									 */
/*		Global parameters for C4.5				 */
/*		--------------------------				 */
/*									 */
/*************************************************************************/


extern short		VERBOSITY;	/* verbosity level (0 = none) */

extern int		MINOBJS;	/* minimum items each side of a cut */

extern Boolean		GAINRATIO;	/* true=gain ratio, false=gain */

//extern float		CF;		/* confidence limit for tree pruning */

extern Boolean		FIRSTHSP, LASTHSP;

extern int MIN_INTRON_LEN, MIN_INTERNAL_EXON_LEN, MIN_INTRON_LEN_AA, MIN_INTERNAL_EXON_LEN_AA; //amino acid length
extern int MINOBJS;

/*************************************************************************/
/*									 */
/*	  Global data for C4.5 used for building decision trees		 */
/*	  -----------------------------------------------------		 */
/*									 */
/*************************************************************************/


extern float
	*Weight,	/* Weight[i]  = current fraction of item i */
	**Freq,		/* Freq[x][c] = no. items of class c with outcome x */
	*ValFreq,	/* ValFreq[x]   = no. items with outcome x */
	*ClassFreq;	/* ClassFreq[c] = no. items of class c */

extern float
	*Gain,		/* Gain[a] = info gain by split on att a */
	*Info,		/* Info[a] = potential info of split on att a */
	*Bar;		/* Bar[a]  = best threshold for contin att a */
//	*UnknownRate;	/* UnknownRate[a] = current unknown rate for att a */

//char
//	*Tested;	/* Tested[a] = true if att a already tested */

extern Tree	*Raw;
extern char bufDecTree[65536];

/*************************************************************************/
/*									 */
/*	  extern functions in C4.5		 */
/*	  -----------------------------------------------------		 */
/*									 */
/*************************************************************************/

extern void InitialiseTreeData();
extern void InitialiseWeights();
extern Tree FormTree(int Fp, int Lp);
extern Boolean Prune(Tree T);
extern void ReleaseTree(Tree Node);
extern void PrintTree(Tree T);


#ifdef __cplusplus
}
#endif
