/*************************************************************************/
/*								 	 */
/*    Central tree-forming algorithm incorporating all criteria  	 */
/*    ---------------------------------------------------------	 	 */
/*								 	 */
/*************************************************************************/

//CHANGES TO C4.5:
//only 1 attribute, treated as "continuous" attribute (no sorting!)
//every time 1 cut is found, two branches are descende, then test the same attribute again
//(no need: discrete attributes, the "Tested[]" flag)


//#include "c45defns.h"
#include "c45extern.h"

/*  global data */

float
	*Weight,	/* Weight[i]  = current fraction of item i */ 
	**Freq,		/* Freq[x][c] = no. items of class c with outcome x */
	*ValFreq,	/* ValFreq[x]   = no. items with outcome x */
	*ClassFreq;	/* ClassFreq[c] = no. items of class c */

float
	*Gain,		/* Gain[a] = info gain by split on att a */
	*Info,		/* Info[a] = potential info of split on att a */
	*Bar;		/* Bar[a]  = best threshold for contin att a */

short	MaxAtt = 0; //only 1 attribute
short	MaxClass = 1; //intron-exon (2 classes)
short	MaxDiscrVal = 2; //do we need this?

//initialize MINOBJS to 2 as regular c4.5 stuff
int	MINOBJS = 2; //Min(MIN_INTRON_LEN, MIN_INTERNAL_EXON_LEN)?

int MIN_INTRON_LEN = 11;//15;
int MIN_INTERNAL_EXON_LEN = 20;
int MIN_INTRON_LEN_AA; // = MIN_INTRON_LEN / 3 + 1;//5; //used by gblastg version 2.3 (=15/3, in terms of amino acids)
int MIN_INTERNAL_EXON_LEN_AA = 6; //20/3


Boolean GAINRATIO = true;

short VERBOSITY = 0;

int		MaxItem; //read from HSP (length of HSP)!
Description	*Item; //read from HSP (each item: weight, class)
float* HSPAlignScores; //read from HSP (according to alignscore!)

Tree	*Raw;

Boolean FIRSTHSP = false, LASTHSP = false;

/*  External variables initialised here  */

extern float
	*SplitGain,	/* SplitGain[i] = gain with att value of item i as threshold */
	*SplitInfo;	/* SplitInfo[i] = potential info ditto */



extern void EvalContinuousAtt(Attribute Att, int Fp, int Lp);
extern void ContinTest(Tree Node, Attribute Att);
extern float Worth(float ThisInfo, float ThisGain, float MinGain);


/*************************************************************************/
/*								 	 */
/*		Allocate space for tree tables			 	 */
/*								 	 */
/*************************************************************************/


void InitialiseTreeData()
/*  ------------------  */
{ 
    DiscrValue v;
//    Attribute a;


    Gain	= (float *) calloc(MaxAtt+1, sizeof(float));
    Info	= (float *) calloc(MaxAtt+1, sizeof(float));
    Bar		= (float *) calloc(MaxAtt+1, sizeof(float));

    SplitGain = (float *) calloc(MaxItem+1, sizeof(float));
    SplitInfo = (float *) calloc(MaxItem+1, sizeof(float));

    Weight = (float *) calloc(MaxItem+1, sizeof(float));

    Freq  = (float **) calloc(MaxDiscrVal+1, sizeof(float *));
    ForEach(v, 0, MaxDiscrVal)
    {
	Freq[v]  = (float *) calloc(MaxClass+1, sizeof(float));
    }

    ValFreq = (float *) calloc(MaxDiscrVal+1, sizeof(float));
    ClassFreq = (float *) calloc(MaxClass+1, sizeof(float));

}



/*************************************************************************/
/*								 	 */
/*		Initialise the weight of each item		 	 */
/*								 	 */
/*************************************************************************/


void InitialiseWeights()
/*  -----------------  */
{
    int i;

    ForEach(i, 0, MaxItem)
    {
        Weight[i] = HSPAlignScores[i];//1.0;
    }
}

/*************************************************************************/
/*								 	 */
/*  Build a decision tree for the cases Fp through Lp:		 	 */
/*								 	 */
/*  - if all cases are of the same class, the tree is a leaf and so	 */
/*      the leaf is returned labelled with this class		 	 */
/*								 	 */
/*  - for each attribute, calculate the potential information provided 	 */
/*	by a test on the attribute (based on the probabilities of each	 */
/*	case having a particular value for the attribute), and the gain	 */
/*	in information that would result from a test on the attribute	 */
/*	(based on the probabilities of each case with a particular	 */
/*	value for the attribute being of a particular class)		 */
/*								 	 */
/*  - on the basis of these figures, and depending on the current	 */
/*	selection criterion, find the best attribute to branch on. 	 */
/*	Note:  this version will not allow a split on an attribute	 */
/*	unless two or more subsets have at least MINOBJS items. 	 */
/*								 	 */
/*  - try branching and test whether better than forming a leaf	 	 */
/*								 	 */
/*************************************************************************/


Tree FormTree(Fp, Lp)
/*   ---------  */
    int Fp, Lp; 
{ 
    int i, Kp, Ep; //Group();
    float Cases, NoBestClass, KnownCases, CountItems();
    float BestVal, Val, AvGain=0; //Factor, Worth();
    Attribute Att, BestAtt, Possible=0;
    ClassNo c, BestClass;
    Tree Node, Leaf();
    DiscrValue v;
//    Boolean PrevAllKnown;

    Cases = CountItems(Fp, Lp);

    /*  Generate the class frequency distribution  */

    ForEach(c, 0, MaxClass)
    {
	ClassFreq[c] = 0;
    }
    ForEach(i, Fp, Lp)
    { 
	ClassFreq[ Class(Item[i]) ] += Weight[i];
    } 

    /*  Find the most frequent class  */

    BestClass = 0;
    ForEach(c, 0, MaxClass)
    {
	if ( ClassFreq[c] > ClassFreq[BestClass] )
	{
	    BestClass = c;
		//printf("best class: %d", BestClass);
	}
    }
    NoBestClass = ClassFreq[BestClass];

    Node = Leaf(ClassFreq, BestClass, Cases, Cases - NoBestClass);
	//Node->Parent = Nil; //root node's parent is Nil
	Node->Start = Fp;
	Node->End = Lp;

    /*  If all cases are of the same class or there are not enough
	cases to divide, the tree is a leaf  */

    if ( NoBestClass == Cases  || ((Lp-Fp+1) < 2*MINOBJS) ) //Cases < 2 * MINOBJS )
    { 
	return Node;
    } 

    Verbosity(1)
    	printf("\n%d items, total weight %.1f\n", Lp - Fp + 1, Cases);

    /*  For each available attribute, find the information and gain  */

    ForEach(Att, 0, MaxAtt) 
    { 
	Gain[Att] = -Epsilon;

	    /*  continuous attribute  */

	    EvalContinuousAtt(Att, Fp, Lp);

		/*  Update average gain, excluding attributes with very many values  */

//		if ( Gain[Att] > -Epsilon &&
//			 ( MultiVal || MaxAttVal[Att] < 0.3 * (MaxItem + 1) ) )
//		{
			Possible++;
			AvGain += Gain[Att];
//		}
    } 

    /*  Find the best attribute according to the given criterion  */

    BestVal = -Epsilon;
    BestAtt = None;
    AvGain  = ( Possible ? AvGain / Possible : 1E6 );

    Verbosity(2)
    {
	if ( AvGain < 1E6 ) printf("\taverage gain %.3f\n", AvGain);
    }

    ForEach(Att, 0, MaxAtt) 
    { 
	if ( Gain[Att] > -Epsilon )
	{ 
	    Val = Worth(Info[Att], Gain[Att], AvGain);
	    if ( Val > BestVal ) 
	    { 
	        BestAtt  = Att; 
	        BestVal = Val;
	    } 
	} 
    } 

    /*  Decide whether to branch or not  */ 

    if ( BestAtt != None )
    { 
	Verbosity(1)
	{
	    //printf("\tbest attribute %s", AttName[BestAtt]);
	    //if ( ! MaxAttVal[BestAtt] )
	    //{
		printf(" cut %.3f", Bar[BestAtt]);
	    //}
	    printf(" inf %.3f gain %.3f val %.3f\n",
		   Info[BestAtt], Gain[BestAtt], BestVal);
	}	

	/*  Build a node of the selected test  */

	    /*  Continuous attribute  */

	    ContinTest(Node, BestAtt);

	/*  Remove unknown attribute values  */

//	PrevAllKnown = AllKnown;

	Kp = Fp; //Group(0, Fp, Lp, Node) + 1;
//	if ( Kp != Fp ) AllKnown = false;
	KnownCases = Cases;// - CountItems(Fp, Kp-1);
//	UnknownRate[BestAtt] = (Cases - KnownCases) / (Cases + 0.001);

	/*  Recursive divide and conquer  */

//	++Tested[BestAtt];

	Ep = Kp - 1;
	Node->Errors = 0;

	ForEach(v, 1, Node->Forks)
	{
		//Ep = Group(v, Kp, Lp, Node);
		if (v == 1) //Branch[1] is from Kp to Node->Cut
		{
			Ep = (int)Node->Cut;
		}
		else //Branch[2] is from Node->Cut+1 to Lp
		{
			Kp = (int)Node->Cut + 1;
			Ep = Lp;
		}

	    if ( Kp <= Ep )
	    {
/*		Factor = CountItems(Kp, Ep) / KnownCases;

		ForEach(i, Fp, Kp-1)
		{
		    Weight[i] *= Factor;
		}
*/
		//Node->Branch[v] = FormTree(Fp, Ep); //recursion
		Node->Branch[v] = FormTree(Kp, Ep);
		Node->Errors += Node->Branch[v]->Errors;

/*		Group(0, Fp, Ep, Node);
		ForEach(i, Fp, Kp-1)
		{
		    Weight[i] /= Factor;
		}
*/	    }
	    else
	    {
		Node->Branch[v] = Leaf(Node->ClassDist, BestClass, 0.0, 0.0);
	    }
	}

//	--Tested[BestAtt];
//	AllKnown = PrevAllKnown;

	/*  See whether we would have been no worse off with a leaf  */

	if ( Node->Errors >= Cases - NoBestClass - Epsilon )
	{ 
	    Node->NodeType = 0;
	} 
    }

    return Node; 
} 



/*************************************************************************/
/*								 	 */
/*  Group together the items corresponding to branch V of a test 	 */
/*  and return the index of the last such			 	 */
/*								 	 */
/*  Note: if V equals zero, group the unknown values		 	 */
/*								 	 */
/*************************************************************************/

/*
int Group(V, Fp, Lp, TestNode)
    DiscrValue V;
    int Fp, Lp;
    Tree TestNode;
{
    int i, last_index;
//    Attribute Att;
    float Thresh;
//    Set SS;
//    void Swap();

//    Att = TestNode->Tested;

	last_index = Fp - 1;
    if ( V )
    {
//  Group items on the value of attribute Att, and depending on the type of branch

	switch ( TestNode->NodeType )
	{
//	    case BrDiscr:

//		ForEach(i, Fp, Lp)
//		{
//		    if ( DVal(Item[i], Att) == V ) Swap(Fp++, i);
//		}
//		break;

	    case ThreshContin:

		Thresh = TestNode->Cut;
		ForEach(i, Fp, Lp)
		{
		    //if ( (CVal(Item[i], Att) <= Thresh) == (V == 1) ) Swap(Fp++, i);
			if ( ((float)i <= Thresh) == (V == 1) ) last_index = i; //Swap(Fp++, i);
		}
		break;

//	    case BrSubset:

//		SS = TestNode->Subset[V];
//		ForEach(i, Fp, Lp)
//		{
//		    if ( In(DVal(Item[i], Att), SS) ) Swap(Fp++, i);
//		}
//		break;
	}
    }
    else
    {
	//  Group together unknown values 

//	switch ( TestNode->NodeType )
//	{
//	    case BrDiscr:
//	    case BrSubset:

//		ForEach(i, Fp, Lp)
//		{
//		    if ( ! DVal(Item[i], Att) ) Swap(Fp++, i);
//		}
//		break;

//	    case ThreshContin:

//		ForEach(i, Fp, Lp)
//		{
//		    if ( CVal(Item[i], Att) == Unknown ) Swap(Fp++, i);
//		}
//		break;
//	}
  }

    //return Fp - 1;
	return last_index;
}

*/

/*************************************************************************/
/*								 	 */
/*	Return the total weight of items from Fp to Lp		 	 */
/*								 	 */
/*************************************************************************/


float CountItems(Fp, Lp)
/*        ----------  */
    int Fp, Lp;
{

    register float Sum=0.0, *Wt, *LWt;

//    if ( AllKnown ) return Lp - Fp + 1;

    for ( Wt = Weight + Fp, LWt = Weight + Lp ; Wt <= LWt ; )
    {
	Sum += *Wt++;
    }

    return Sum;

}



/*************************************************************************/
/*                                                               	 */
/*		Exchange items at a and b			 	 */
/*									 */
/*************************************************************************/


void Swap(a,b)
/*   ----  */
    int a, b;
{
    register Description Hold;
    register float HoldW;

    Hold = Item[a];
    Item[a] = Item[b];
    Item[b] = Hold;

    HoldW = Weight[a];
    Weight[a] = Weight[b];
    Weight[b] = HoldW;
}


