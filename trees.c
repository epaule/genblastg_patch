/*************************************************************************/
/*									 */
/*	Routines for displaying, building, saving and restoring trees	 */
/*	-------------------------------------------------------------	 */
/*									 */
/*************************************************************************/


//#include "c45defns.h"
#include "c45extern.h"
#include <memory.h>
#include <string.h>

#define	Tab		"|   "
#define	TabSize		4
#define	Width		80	/* approx max width of printed trees */

	/*  If lines look like getting too long while a tree is being
	    printed, subtrees are broken off and printed separately after
	    the main tree is finished	 */

short	Subtree;		/* highest subtree to be printed */
Tree	Subdef[100];		/* pointers to subtrees */

//FILE	*TRf = 0, *fopen();	/* file pointer for tree i/o */
//char	Fn[500];		/* file name */

char	bufDecTree[65536]; //64k buffer, used to store strings of decision tree output

void Show(Tree T, short Sh);

void AppendBuffer(const char* newStr)
{
	int size;

	size = strlen(bufDecTree);
	strcat(bufDecTree, newStr);
}

/*************************************************************************/
/*									 */
/*	Find the maximum single line size for non-leaf subtree St.	 */
/*	The line format is						 */
/*			<attribute> <> X.xx:[ <class (<Items>)], or	 */
/*			<attribute> = <DVal>:[ <class> (<Items>)]	 */
/*									 */
/*************************************************************************/


short MaxLine(St)

    Tree St;
{
//    Attribute a;
    DiscrValue v, MaxV, Next;
    short Ll, MaxLl=0;

//    a = St->Tested;

    MaxV = St->Forks;
    ForEach(v, 1, MaxV)
    {
	//Ll = ( St->NodeType == 2 ? 4 : strlen(AttValName[a][v]) ) + 1;
		Ll = 4;

	//  Find the appropriate branch  

        Next = v;

	if ( ! St->Branch[Next]->NodeType )
	{
//	    Ll += strlen(ClassName[St->Branch[Next]->Leaf]) + 6;
		Ll += 6;
	}
	MaxLl = Max(MaxLl, Ll);
    }

    //return strlen(AttName[a]) + 4 + MaxLl;
	return 4 + MaxLl;
}



/*************************************************************************/
/*								   	 */
/*	Indent Sh columns					  	 */
/*								  	 */
/*************************************************************************/

void  Indent(Sh, Mark)

    short Sh;
    char *Mark;
{
    //printf("\n");
	char buf[1024];
	sprintf(buf, "\n");
	AppendBuffer(buf);
    while ( Sh-- ) sprintf(buf, "%s", Mark);//printf("%s", Mark);
	AppendBuffer(buf);
}



/*************************************************************************/
/*									 */
/*	Print a node T with offset Sh, branch value v, and continue	 */
/*									 */
/*************************************************************************/


void ShowBranch(Sh, T, v)

    short Sh;
    Tree T;
    DiscrValue v;
{
    DiscrValue Pv, Last;
    Attribute Att;
    Boolean FirstValue;
    short TextWidth, Skip, Values=0, i;
	char buf[1024];
    
//    Att = T->Tested;

    switch ( T->NodeType )
    {
//	case BrDiscr:

//	    Indent(Sh, Tab);

//	    printf("%s = %s:", AttName[Att], AttValName[Att][v]);
//	    break;

	case ThreshContin:

	    Indent(Sh, Tab);

	    //printf("%s %s %g ", AttName[Att],
		//printf(" %s %g ",
		sprintf(buf, " %s %g ",
			( v == 1 ? "<=" : ">" ), T->Cut);
		AppendBuffer(buf);

	    if ( T->Lower != T->Upper )
	    {
		//printf("[%g,%g]", T->Lower, T->Upper);
			sprintf(buf, "[%g,%g]", T->Lower, T->Upper);
			AppendBuffer(buf);
	    }

	    //printf(":");
		sprintf(buf, ":");
		AppendBuffer(buf);
	    break;

/*	case BrSubset:

	    //  Count values at this branch  

	    ForEach(Pv, 1, MaxAttVal[Att])
	    {
		if ( In(Pv, T->Subset[v]) )
		{
		    Last = Pv;
		    Values++;
		}
	    }
	    if ( ! Values ) return;

	    Indent(Sh, Tab);

	    if ( Values == 1 )
	    {
		printf("%s = %s:", AttName[Att], AttValName[Att][Last]);
		break;
	    }

	    printf("%s in {", AttName[Att]);
	    FirstValue = true;
	    Skip = TextWidth = strlen(AttName[Att]) + 5;

	    ForEach(Pv, 1, MaxAttVal[Att])
	    {
		if ( In(Pv, T->Subset[v]) )
		{
		    if ( ! FirstValue &&
			 TextWidth + strlen(AttValName[Att][Pv]) + 11 > Width )
		    {
		  	Indent(Sh, Tab);
			ForEach(i, 1, Skip) putchar(' ');

			TextWidth = Skip;
			FirstValue = true;
		    }

		    printf("%s%c", AttValName[Att][Pv], Pv == Last ? '}' : ',');
		    TextWidth += strlen(AttValName[Att][Pv]) + 1;
		    FirstValue = false;
		}
	    }
	    putchar(':');
*/    }

    Show(T->Branch[v], Sh+1);
}

/*************************************************************************/
/*									 */
/*	Display the tree T with offset Sh				 */
/*									 */
/*************************************************************************/


void Show(T, Sh)
    Tree T;
    short Sh;
{
    DiscrValue v, MaxV;
	char buf[1024];
    short MaxLine();

    if ( T->NodeType )
    {
	//  See whether separate subtree needed  

	if ( T != Nil && Sh && Sh * TabSize + MaxLine(T) > Width )
	{
	    if ( Subtree < 99 )
	    {
		Subdef[++Subtree] = T;
		//printf("[S%d]", Subtree);
		sprintf(buf, "[S%d]", Subtree);
		AppendBuffer(buf);
	    }
	    else
	    {
		//printf("[S??]");
			sprintf(buf, "[S??]");
			AppendBuffer(buf);
	    }
	}
	else
	{
	    MaxV = T->Forks;

	    //  Print simple cases first 

	    ForEach(v, 1, MaxV)
	    {
		if ( ! T->Branch[v]->NodeType )
		{
		    ShowBranch(Sh, T, v);
		}
	    }

	    //  Print subtrees 

	    ForEach(v, 1, MaxV)
	    {
		if ( T->Branch[v]->NodeType )
		{
		    ShowBranch(Sh, T, v);
		}
	    }
	}
    }
    else
    {
	//printf(" %s (%.1f", ClassName[T->Leaf], T->Items);
	//printf(" %d (%.1f", T->Leaf, T->Items);
		sprintf(buf, " %d (%.1f", T->Leaf, T->Items);
		AppendBuffer(buf);
	if ( T->Errors > 0 ) sprintf(buf, "/%.1f", T->Errors); //printf("/%.1f", T->Errors);
	AppendBuffer(buf);
	//printf(")");
	sprintf(buf, ")");
	AppendBuffer(buf);
    }
}




/*************************************************************************/
/*									 */
/*	Display entire decision tree T					 */
/*									 */
/*************************************************************************/


void PrintTree(T)
    Tree T;
{
    short s;
	char buf[1024];

	bufDecTree[0] = '\0'; //reset bufDecTree

    Subtree=0;
    //printf("Decision Tree:\n");
	sprintf(bufDecTree, "Decision Tree:\n");
    Show(T, 0);
    //printf("\n");
	sprintf(buf, "\n");
	AppendBuffer(buf);

    ForEach(s, 1, Subtree)
    {
	//printf("\n\nSubtree [S%d]\n", s);
		sprintf(buf, "\n\nSubtree [S%d]\n", s);
		AppendBuffer(buf);
	Show(Subdef[s], 0);
	sprintf(buf, "\n");//printf("\n");
	AppendBuffer(buf);
    }
    //printf("\n");
	sprintf(buf, "\n");
	AppendBuffer(buf);
}




/*************************************************************************/
/*									 */
/*	Save entire decision tree T in file with extension Extension	 */
/*									 */
/*************************************************************************/

/*
    SaveTree(T, Extension)
    Tree T;
    String Extension;
{
    static char *LastExt="";

    if ( strcmp(LastExt, Extension) )
    {
	LastExt = Extension;

	if ( TRf ) fclose(TRf);

	strcpy(Fn, FileName);
	strcat(Fn, Extension);
	if ( ! ( TRf = fopen(Fn, "w") ) )
	    Error(0, Fn, " for writing");
    }

    putc('\n', TRf);
    OutTree(T);

    SaveDiscreteNames();
}
*/


/*************************************************************************/
/*									 */
/*	Save tree T as characters					 */
/*									 */
/*************************************************************************/

/*
    OutTree(T)
    Tree T;
{
    DiscrValue v;
    int Bytes;

    StreamOut((char *) &T->NodeType, sizeof(short));
    StreamOut((char *) &T->Leaf, sizeof(ClassNo));
    StreamOut((char *) &T->Items, sizeof(float));
    StreamOut((char *) &T->Errors, sizeof(float));
    StreamOut((char *) T->ClassDist, (MaxClass + 1) * sizeof(float));

    if ( T->NodeType )
    {
	StreamOut((char *) &T->Tested, sizeof(Attribute));
	StreamOut((char *) &T->Forks, sizeof(short));

	switch ( T->NodeType )
	{
	    case BrDiscr:
		break;

	    case ThreshContin:
		StreamOut((char *) &T->Cut, sizeof(float));
		StreamOut((char *) &T->Lower, sizeof(float));
		StreamOut((char *) &T->Upper, sizeof(float));
		break;

	    case BrSubset:
		Bytes = (MaxAttVal[T->Tested]>>3) + 1;
		ForEach(v, 1, T->Forks)
		{
		    StreamOut((char *) T->Subset[v], Bytes);
		}
		break;
	}

	ForEach(v, 1, T->Forks)
	{
	    OutTree(T->Branch[v]);
	}
    }
}
*/


/*************************************************************************/
/*									 */
/*	Retrieve entire decision tree with extension Extension		 */
/*									 */
/*************************************************************************/

/*
Tree GetTree(Extension)
    String Extension;
{
    Tree Hold, InTree();
    static char *LastExt="";

    if ( strcmp(LastExt, Extension) )
    {
	LastExt = Extension;

	if ( TRf ) fclose(TRf);

	strcpy(Fn, FileName);
	strcat(Fn, Extension);
	if ( ! ( TRf = fopen(Fn, "r") ) ) Error(0, Fn, "");
    }

    if ( ! TRf || getc(TRf) == EOF ) return Nil;

    Hold = InTree();

    RecoverDiscreteNames();

    return Hold;
}
*/


/*************************************************************************/
/*									 */
/*	Retrieve tree from saved characters				 */
/*									 */
/*************************************************************************/

/*
Tree InTree()
{
    Tree T;
    DiscrValue v;
    int Bytes;

    T = (Tree) malloc(sizeof(TreeRec));

    StreamIn((char *) &T->NodeType, sizeof(short));
    StreamIn((char *) &T->Leaf, sizeof(ClassNo));
    StreamIn((char *) &T->Items, sizeof(float));
    StreamIn((char *) &T->Errors, sizeof(float));

    T->ClassDist = (float *) calloc(MaxClass+1, sizeof(float));
    StreamIn((char *) T->ClassDist, (MaxClass + 1) * sizeof(float));

    if ( T->NodeType )
    {
	StreamIn((char *) &T->Tested, sizeof(Attribute));
	StreamIn((char *) &T->Forks, sizeof(short));

	switch ( T->NodeType )
	{
	    case BrDiscr:
		break;

	    case ThreshContin:
		StreamIn((char *) &T->Cut, sizeof(float));
		StreamIn((char *) &T->Lower, sizeof(float));
		StreamIn((char *) &T->Upper, sizeof(float));
		break;

	    case BrSubset:
		T->Subset = (Set *) calloc(T->Forks + 1, sizeof(Set));

		Bytes = (MaxAttVal[T->Tested]>>3) + 1;
		ForEach(v, 1, T->Forks)
		{
		    T->Subset[v] = (Set) malloc(Bytes);
		    StreamIn((char *) T->Subset[v], Bytes);
		}
	}

	T->Branch = (Tree *) calloc(T->Forks + 1, sizeof(Tree));
	ForEach(v, 1, T->Forks)
	{
	    T->Branch[v] = InTree();
	}
    }

    return T;
}
*/


/*************************************************************************/
/*									 */
/*	Stream characters to/from file TRf from/to an address		 */
/*									 */
/*************************************************************************/

/*
    StreamOut(s, n)

    String s;
    int n;
{
    while ( n-- ) putc(*s++, TRf);
}



    StreamIn(s, n)

	String s;
    int n;
{
    while ( n-- ) *s++ = getc(TRf);
}

*/

/*************************************************************************/
/*									 */
/*	Free up space taken up by tree Node				 */
/*									 */
/*************************************************************************/


void ReleaseTree(Node)
/*  -------  */
    Tree Node;
{
    DiscrValue v;

    if ( Node->NodeType )
    {
	ForEach(v, 1, Node->Forks)
	{
	    ReleaseTree(Node->Branch[v]);
	}

	//cfree(Node->Branch);
	free(Node->Branch);

    }

    //cfree(Node->ClassDist);
	free(Node->ClassDist);
    //cfree(Node);
	free(Node);
}



/*************************************************************************/
/*									 */
/*	Construct a leaf in a given node				 */
/*									 */
/*************************************************************************/


Tree Leaf(ClassFreq, NodeClass, Cases, Errors)
/*   ----  */
    float *ClassFreq;
    ClassNo NodeClass;
    float Cases, Errors;
{
    Tree Node;

    Node = (Tree) calloc(1, sizeof(TreeRec));

    Node->ClassDist = (float *) calloc(MaxClass+1, sizeof(float));
    memcpy(Node->ClassDist, ClassFreq, (MaxClass+1) * sizeof(float));
    
    Node->NodeType	= 0; 
    Node->Leaf		= NodeClass;
    Node->Items		= Cases;
    Node->Errors	= Errors;

    return Node; 
}



/*************************************************************************/
/*									 */
/*	Insert branches in a node 	                 		 */
/*									 */
/*************************************************************************/


void Sprout(Node, Branches)
/*  ------  */
    Tree Node;
    DiscrValue Branches;
{
    Node->Forks = Branches;
    
    Node->Branch = (Tree *) calloc(Branches+1, sizeof(Tree));

}



/*************************************************************************/
/*									 */
/*	Count the nodes in a tree					 */
/*									 */
/*************************************************************************/

	
int TreeSize(Node)
/*  --------  */
    Tree Node;
{
    int Sum=0;
    DiscrValue v;

    if ( Node->NodeType )
    {
	ForEach(v, 1, Node->Forks)
	{
	    Sum += TreeSize(Node->Branch[v]);
	}
    }

    return Sum + 1;
}



/*************************************************************************/
/*									 */
/*	Return a copy of tree T						 */
/*									 */
/*************************************************************************/


Tree CopyTree(T)
/*   ---------  */
    Tree T;
{
    DiscrValue v;
    Tree New;

    New = (Tree) malloc(sizeof(TreeRec));
    memcpy(New, T, sizeof(TreeRec));

    New->ClassDist = (float *) calloc(MaxClass+1, sizeof(float));
    memcpy(New->ClassDist, T->ClassDist, (MaxClass + 1) * sizeof(float));

    if ( T->NodeType )
    {
	New->Branch = (Tree *) calloc(T->Forks + 1, sizeof(Tree));
	ForEach(v, 1, T->Forks)
	{
	    New->Branch[v] = CopyTree(T->Branch[v]);
	}
    }

    return New;
}



/*************************************************************************/
/*									 */
/*	Save attribute values read with "discrete N"			 */
/*									 */
/*************************************************************************/

/*
    SaveDiscreteNames()

{
    Attribute Att;
    DiscrValue v;
    int Length;

    ForEach(Att, 0, MaxAtt)
    {
	//if ( SpecialStatus[Att] != DISCRETE ) continue;

	StreamOut((char *) &MaxAttVal[Att], sizeof(int));

	ForEach(v, 1, MaxAttVal[Att])
	{
	    Length = strlen(AttValName[Att][v]) + 1;

	    StreamOut((char *) &Length, sizeof(int));
	    StreamOut((char *) AttValName[Att][v], Length);
	}
    }
}
*/


/*************************************************************************/
/*									 */
/*	Recover attribute values read with "discrete N"			 */
/*									 */
/*************************************************************************/

/*
    RecoverDiscreteNames()

{
    Attribute Att;
    DiscrValue v;
    int Length;

    ForEach(Att, 0, MaxAtt)
    {
	//if ( SpecialStatus[Att] != DISCRETE ) continue;

	StreamIn(&MaxAttVal[Att], sizeof(int));

	ForEach(v, 1, MaxAttVal[Att])
	{
	    StreamIn(&Length, sizeof(int));

	    AttValName[Att][v] = (char *) malloc(Length);
	    StreamIn(AttValName[Att][v], Length);
	}
    }
}
*/
