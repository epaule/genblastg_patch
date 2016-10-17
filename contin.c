/*************************************************************************/
/*                                                                	 */
/*	Evaluation of a test on a continuous valued attribute	  	 */
/*	-----------------------------------------------------	  	 */
/*								  	 */
/*************************************************************************/


#include "c45extern.h"


float
	*SplitGain,	/* SplitGain[i] = gain with att value of item i as threshold */
	*SplitInfo;	/* SplitInfo[i] = potential info ditto */

extern void ResetFreq(DiscrValue MaxVal);
extern void Sprout(Tree Node, DiscrValue Branches);
extern float ComputeGain(float BaseInfo, float UnknFrac, DiscrValue MaxVal, float TotalItems, int i, int Sp, int Ep);
extern float TotalInfo(float V[], DiscrValue MinVal, DiscrValue MaxVal);
extern float Worth(float ThisInfo, float ThisGain, float MinGain);



/*************************************************************************/
/*								  	 */
/*  Continuous attributes are treated as if they have possible values	 */
/*	0 (unknown), 1 (less than cut), 2(greater than cut)	  	 */
/*  This routine finds the best cut for items Fp through Lp and sets	 */
/*  Info[], Gain[] and Bar[]						 */
/*								  	 */
/*************************************************************************/


void EvalContinuousAtt(Att, Fp, Lp)
/*  -----------------  */ 
    Attribute Att;
    int Fp, Lp; 
{ 
    int i, BestI, Xp, Tries=0;
    float Items, KnownItems, LowItems, MinSplit, CountItems(), TotalLowWeights; 
    ClassNo c;
    float AvGain=0, Val, BestVal, BaseInfo, ThreshCost; //,
	//ComputeGain(), TotalInfo(), Worth();
//    void Swap();

//    Verbosity(2) printf("\tAtt %s", AttName[Att]);
    Verbosity(3) printf("\n");

    ResetFreq(2);

    /*  Omit and count unknown values */

    Items = CountItems(Fp, Lp);
    Xp = Fp;

    ValFreq[0] = 0;
    KnownItems = Items;// - ValFreq[0];

//    Quicksort(Xp, Lp, Att, Swap); //important! skip this and we can use the rest codes as for a "continuous" attribute!

    /*  Count base values and determine base information  */
    ForEach(i, Xp, Lp)
    {
	Freq[ 2 ][ Class(Item[i]) ] += Weight[i];
	SplitGain[i] = -Epsilon;
	SplitInfo[i] = 0;
    }

    BaseInfo = TotalInfo(Freq[2], 0, MaxClass) / KnownItems;

    /*  Try possible cuts between items i and i+1, and determine the
	information and gain of the split in each case.  We have to be wary
	of splitting a small number of items off one end, as we can always
	split off a single item, but this has little predictive power.  */

    //MinSplit = 0.10 * KnownItems / (MaxClass + 1);
	MinSplit = 0.10 * (Lp-Fp+1) / (MaxClass + 1);
    if ( MinSplit <= MINOBJS ) MinSplit = MINOBJS;
    else
    if ( MinSplit > 25 ) MinSplit = 25; //heuristics, I guess?

    LowItems = 0; //keep track of number of items (count)
	TotalLowWeights = 0; //keep track of weights of all previous items
    ForEach(i, Xp, Lp - 1)
    {
	c = Class(Item[i]);
	LowItems   += 1;//Weight[i];
	TotalLowWeights += Weight[i];
	Freq[1][c] += Weight[i];
	Freq[2][c] -= Weight[i];

	if ( LowItems < MinSplit ) continue;
	else
	//if ( LowItems > KnownItems - MinSplit ) break;
	if ( LowItems > (Lp-Fp+1) - MinSplit ) break;

	//if ( CVal(Item[i],Att) < CVal(Item[i+1],Att) - 1E-5 ) //what's this? => if: there's a significant gap between this value and its next value
	//{
	    ValFreq[1] = TotalLowWeights; //LowItems;
	    ValFreq[2] = KnownItems - TotalLowWeights; //LowItems;
	    SplitGain[i] = ComputeGain(BaseInfo, 0, 2, KnownItems, i, Xp, Lp-1); //UnknownRate[Att], 2, KnownItems);
	    SplitInfo[i] = TotalInfo(ValFreq, 0, 2) / Items;
	    AvGain += SplitGain[i];
	    Tries++;

	    Verbosity(3)
	    {	printf("\t\tCut at %.1f  (gain %.3f, val %.3f):",
	               ((float)i + i+1 )/2, //( CVal(Item[i],Att) + CVal(Item[i+1],Att) ) / 2,
	    	       SplitGain[i],
	    	       Worth(SplitInfo[i], SplitGain[i], Epsilon));
	    	       //PrintDistribution(Att, 2, true);
	    }
	//}
    }

    /*  Find the best attribute according to the given criterion  */

    ThreshCost = Log(Tries) / Items;

    BestVal = 0;
    BestI   = None;
    ForEach(i, Xp, Lp - 1)
    {
	if ( (Val = SplitGain[i] - ThreshCost) > BestVal )
	{
	    BestI   = i;
	    BestVal = Val;
	}
    }

    /*  If a test on the attribute is able to make a gain,
	set the best break point, gain and information  */ 

    if ( BestI == None )
    {
	Gain[Att] = -Epsilon;
	Info[Att] = 0.0;

	Verbosity(2) printf("\tno gain\n");
    }
    else
    {
	//instead of cutting by a value, we are cutting by a position!
	//note the "position" is relative to the beginning of HSP
	Bar[Att]  = (float)BestI; //(CVal(Item[BestI],Att) + CVal(Item[BestI+1],Att)) / 2;
	Gain[Att] = BestVal;
	Info[Att] = SplitInfo[BestI];

	Verbosity(2)
	    printf("\tcut=%.3f, inf %.3f, gain %.3f\n",
		   Bar[Att], Info[Att], Gain[Att]);
    }
} 



/*************************************************************************/
/*                                                                	 */
/*  Change a leaf into a test on a continuous attribute           	 */
/*                                                                	 */
/*************************************************************************/


void ContinTest(Node, Att)
/*  ----------  */
    Tree Node;
    Attribute Att;
{
    float Thresh, GreatestValueBelow();
    float CountItems();

    Sprout(Node, 2);

    Thresh = Bar[Att];//GreatestValueBelow(Att, Bar[Att]);

    Node->NodeType	= ThreshContin;
//    Node->Tested	= Att;
    Node->Cut		=
    Node->Lower		=
    Node->Upper		= Thresh;
    Node->Errors        = 0;

/*	//3 branches, we only care about Branch[1] and Branch[2]
	Node->Branch[1]->Start = Node->Start; //add the actual split information, for later recover donor_segments from it
	Node->Branch[1]->End = (int)Thresh;
	Node->Branch[2]->Start = Node->Branch[1]->End + 1;
	Node->Branch[2]->End = Node->End;
*/
}



/*************************************************************************/
/*                                                                	 */
/*  Return the greatest value of attribute Att below threshold t  	 */
/*                                                                	 */
/*************************************************************************/

/*
float GreatestValueBelow(Att, t)
    Attribute Att;
    float t;
{
    int i;
    float v, Best;
    Boolean NotYet=true;

    ForEach(i, 0, MaxItem)
    {
	v = CVal(Item[i], Att);
	if ( v != Unknown && v <= t && ( NotYet || v > Best ) )
	{
	    Best = v;
	    NotYet = false;
	}
    }

    return Best;
}
*/
