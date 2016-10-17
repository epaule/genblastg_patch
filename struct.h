//***************STRUCT.H************************

//

//misc data structures

//

//

//Author: Rong She

//Date: Nov 2008

//***********************************************



#if defined(_MSC_VER)
	#pragma warning(disable: 4786)
	#pragma warning(disable: 4503)
#endif


#ifndef STRUCT_H /* STRUCT_H */



#define STRUCT_H



#include <iostream>
#include <string>
using namespace std;



struct SegmentsInThreeBounds

{

	int exon_seg_start;

	int exon_seg_end;

	int intron_seg_end;



	SegmentsInThreeBounds(int start, int center, int end)

	{

		exon_seg_start = start;

		exon_seg_end = center;

		intron_seg_end = end;

	}



	friend ostream& operator<<(ostream& os, const SegmentsInThreeBounds& seg)

	{

		os << seg.exon_seg_start << "->" << seg.exon_seg_end << "->" << seg.intron_seg_end << "\n";

		return os;

	}

};

struct Triplet //3 integers
{
	int	hsp_ID;
	int frame_reference;
	int	site_index;

	Triplet(int hspID, int fref, int sitei)
	{
		hsp_ID = hspID;
		frame_reference = fref;
		site_index = sitei;
	}
};

struct DNAString_ProteinStyle
{
	string	center_PString;
	string	front_DString;
	string	end_DString;

	friend ostream& operator<<(ostream& os, const DNAString_ProteinStyle& dstr)
	{
		os << dstr.front_DString << "<" << dstr.center_PString << ">" << dstr.end_DString << "\n";
		return os;
	}
};

#endif /* STRUCT_H */

