#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include "api/api_global.h"
#include "api/BamReader.h"

using namespace std;
using namespace BamTools;

int main(int argc, char *argv[])
{
    //unsigned short count;
    int32_t chrLength, thresHold, inCounter = 0, outCounter = 0, inIdx = -1, outIdx = -1;
    int32_t currMatePos, currPos;
	  std::string filename = argv[1], chr = argv[2];
    BamTools::BamReader reader;
    const BamTools::RefVector *references;
    std::stringstream(argv[3]) >> thresHold;

    //open BAM and its index
    if (!reader.Open(filename)){
        cerr << "Could not open BAM file!"  << endl;
        return -1;
        }
    if ( !reader.LocateIndex() ){
        if ( !reader.CreateIndex() )
            cerr << "Could not open or create index file!" << endl;
        return -1;
        }
    int RefID = reader.GetReferenceID(chr);
    references = &reader.GetReferenceData();
    chrLength = (*references)[RefID].RefLength;
    BamTools::BamRegion region(RefID, 1, RefID, chrLength);

    reader.SetRegion(region);
    references = &reader.GetReferenceData();
    BamTools::BamAlignment al;
    //do the job
    while (reader.GetNextAlignment(al)){
      //cout << al.Length << endl;
      if (inIdx != -1 && outIdx != -1) break;
      currPos = al.Position;

      //go right
      if (al.AlignmentFlag == 97 && al.MateRefID != RefID &&al.MapQuality > 20 && outIdx== -1){
        currMatePos = al.MatePosition;
        reader.Jump(al.MateRefID, al.MatePosition);
        while (reader.GetNextAlignment(al)){
          if (al.Position == currMatePos && al.AlignmentFlag == 145){
            outCounter++;
            if (outCounter > thresHold) {
              outIdx = al.RefID;
            }
            reader.SetRegion(RefID, currPos, RefID, chrLength);
            break;
          }
        }
      }
      //go left
      if (al.AlignmentFlag == 145 && al.MateRefID != RefID &&al.MapQuality > 20 && inIdx == -1){
        currMatePos = al.MatePosition;
        reader.Jump(al.MateRefID, al.MatePosition-al.Length-50);
        while (reader.GetNextAlignment(al)){
          if (al.Position == currMatePos && al.AlignmentFlag == 97){
            inCounter++;
            if (inCounter> thresHold){
              inIdx = al.RefID;
            }
            reader.SetRegion(RefID, currPos, RefID, chrLength);
            break;
          }
        }
      }
    }
	if (inIdx != -1)
		std::cout <<  (*references)[inIdx].RefName << "\t" << chr << "\t" << inCounter << endl;
	if (outIdx != -1)
    std::cout << chr  << "\t" << (*references)[outIdx].RefName << "\t" << outCounter << endl;
}
