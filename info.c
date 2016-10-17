/*************************************************************************/
/*									 */
/*	Calculate information, information gain, and print dists	 */
/*	--------------------------------------------------------	 */
/*									 */
/*************************************************************************/


#include "c45extern.h"

float TotalInfo(float V[], DiscrValue MinVal, DiscrValue MaxVal);

ClassNo	GetBestClass(int Fp, int Lp)
{
	ClassNo c, BestClass;
    int i;
	float *ClsFreq;

    /*  Generate the class frequency distribution  */
	ClsFreq = (float *) calloc(MaxClass+1, sizeof(float));
//    ForEach(c, 0, MaxClass)
//    {
//	ClassFreq[c] = 0;
//    }

    ForEach(i, Fp, Lp)
    { 
	ClsFreq[ Class(Item[i]) ] += Weight[i];
    } 

    /*  Find the most frequent class  */

    BestClass = 0;
    ForEach(c, 0, MaxClass)
    {
	if ( ClsFreq[c] > ClsFreq[BestClass] )
	{
	    BestClass = c;
	}
    }

	free(ClsFreq);

	return BestClass;

}


/*************************************************************************/
/*									 */
/*  Determine the worth of a particular split according to the		 */
/*  operative criterion							 */
/*									 */
/*	    Parameters:							 */
/*		SplitInfo:	potential info of the split		 */
/*		SplitGain:	gain in info of the split		 */
/*		MinGain:	gain above which the Gain Ratio		 */
/*				may be used				 */
/*									 */
/*  If the Gain criterion is being used, the information gain of	 */
/*  the split is returned, but if the Gain Ratio criterion is		 */
/*  being used, the ratio of the information gain of the split to	 */
/*  its potential information is returned.				 */
/*									 */
/*************************************************************************/


float Worth(float ThisInfo, float ThisGain, float MinGain)
/*    -----  
    float ThisInfo, ThisGain, MinGain;
*/
{
    if ( GAINRATIO )
    {
	if ( ThisGain >= MinGain - Epsilon && ThisInfo > Epsilon )
	{
	    return ThisGain / ThisInfo;
	}
	else
	{
	    return -Epsilon;
	}
    }
    else
    {
	return ( ThisInfo > 0 && ThisGain > -Epsilon ? ThisGain : -Epsilon );
    }
}



/*************************************************************************/
/*									 */
/*  Zero the frequency tables Freq[][] and ValFreq[] up to MaxVal	 */
/*									 */
/*************************************************************************/


void ResetFreq(MaxVal)
/*  ---------  */
    DiscrValue MaxVal;
{
    DiscrValue v;
    ClassNo c;

    ForEach(v, 0, MaxVal)
    { 
	ForEach(c, 0, MaxClass)
	{
	    Freq[v][c] = 0;
	}
	ValFreq[v] = 0;
    } 
}



/*************************************************************************/
/*									 */
/*  Given tables Freq[][] and ValFreq[], compute the information gain.	 */
/*									 */
/*	    Parameters:							 */
/*		BaseInfo:	average information for all items with	 */
/*				known values of the test attribute	 */
/*		UnknownRate:	fraction of items with unknown ditto	 */
/*		MaxVal:		number of forks				 */
/*		TotalItems:	number of items with known values of	 */
/*				test att				 */
/*									 */
/*  where Freq[x][y] contains the no. of cases with value x for a	 */
/*  particular attribute that are members of class y,			 */
/*  and ValFreq[x] contains the no. of cases with value x for a		 */
/*  particular attribute						 */
/*									 */
/*************************************************************************/


float ComputeGain(float BaseInfo, float UnknFrac, DiscrValue MaxVal, float TotalItems, int i, int Sp, int Ep)
/*    -----------  
    float BaseInfo, UnknFrac;
    DiscrValue MaxVal;
    float TotalItems;
	int i, Sp, Ep;
*/
{
    DiscrValue v;
    float ThisInfo=0.0, ThisGain;// TotalInfo();
    short ReasonableSubsets=0;

    /*  Check whether all values are unknown or the same  */

    if ( ! TotalItems ) return -Epsilon;

    /*  There must be at least two subsets with MINOBJS items  */

//    ForEach(v, 1, MaxVal)
//    {
	//if ( ValFreq[v] >= MINOBJS ) ReasonableSubsets++;
//    }
	
/*	if (i - Sp + 1 >= MINOBJS) ReasonableSubsets++;
	if (Ep - i >= MINOBJS) ReasonableSubsets++;
*/
	if (FIRSTHSP && Sp == 0) //fist hsp's first segment
	{
		if (i - Sp + 1 >= MINOBJS) ReasonableSubsets++;

		if (GetBestClass(i+1, Ep) == 0) //second branch is intron
		{
			if (Ep - i >= MIN_INTRON_LEN_AA) ReasonableSubsets++;
		}
		else
		{
			if (Ep - i >= MINOBJS) ReasonableSubsets++;
		}
	}
	else
	{
		if (LASTHSP && Ep == MaxItem) //last hsp's last segment
		{
			if (GetBestClass(Sp, i) == 0) //first branch is intron
			{
				if (i - Sp + 1 >= MIN_INTRON_LEN_AA) ReasonableSubsets++;
			}
			else
			{
				if (i - Sp + 1 >= MINOBJS) ReasonableSubsets++;
			}

			if (Ep - i >= MINOBJS) ReasonableSubsets++;
		}
		else
		{
			if (GetBestClass(Sp, i) == 0) //first branch is intron
			{
				if (i - Sp + 1 >= MIN_INTRON_LEN_AA) ReasonableSubsets++;
				if (Ep - i >= MIN_INTERNAL_EXON_LEN_AA) ReasonableSubsets++;
			}
			else //first branch is exon
			{
				if (i - Sp + 1 >= MIN_INTERNAL_EXON_LEN_AA) ReasonableSubsets++;
				if (Ep - i >= MIN_INTRON_LEN_AA) ReasonableSubsets++;
			}
		}
	}	

    if ( ReasonableSubsets < 2 ) return -Epsilon;

    /*  Compute total info after split, by summing the
	info of each of the subsets formed by the test  */

    ForEach(v, 1, MaxVal)
    {
	ThisInfo += TotalInfo(Freq[v], 0, MaxClass);
    }

    /*  Set the gain in information for all items, adjusted for unknowns  */

    ThisGain = (1 - UnknFrac) * (BaseInfo - ThisInfo / TotalItems);

    Verbosity(5)
        printf("ComputeThisGain: items %.1f info %.3f base %.3f unkn %.3f result %.3f\n",
    		TotalItems + ValFreq[0], ThisInfo, BaseInfo, UnknFrac, ThisGain);

    return ThisGain;
}



/*************************************************************************/
/*									 */
/*  Compute the total information in V[ MinVal..MaxVal ]		 */
/*									 */
/*************************************************************************/


float TotalInfo(float V[], DiscrValue MinVal, DiscrValue MaxVal)
/*    ---------  
    float V[];
    DiscrValue MinVal, MaxVal;
*/
{
    DiscrValue v;
    float Sum=0.0;
    float N, TotalItems=0;

    ForEach(v, MinVal, MaxVal)
    {
	N = V[v];

	Sum += N * Log(N);
	TotalItems += N;
    }

    return TotalItems * Log(TotalItems) - Sum;
}



/*************************************************************************/
/*									 */
/*	Print distribution table for given attribute			 */
/*									 */
/*************************************************************************/

/*
    PrintDistribution(Att, MaxVal, ShowNames)

	Attribute Att;
    DiscrValue MaxVal;
    Boolean ShowNames;
{
    DiscrValue v;
    ClassNo c;
    String Val;

    printf("\n\t\t\t ");
    ForEach(c, 0, MaxClass)
    {
	printf("%7.6s", ClassName[c]);
    }
    printf("\n");

    ForEach(v, 0, MaxVal)
    {
	if ( ShowNames )
	{
	    Val = ( !v ? "unknown" :
		    MaxAttVal[Att] ? AttValName[Att][v] :
		    v == 1 ? "below" : "above" );
	    printf("\t\t[%-7.7s:", Val);
	}
	else
	{
	    printf("\t\t[%-7d:", v);
	}

	ForEach(c, 0, MaxClass)
	{
	    printf(" %6.1f", Freq[v][c]);
	}

	printf("]\n");
    }
}
*/
