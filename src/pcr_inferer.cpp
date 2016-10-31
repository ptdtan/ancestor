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
    int32_t uRstart, uRend, uLen, uWall;
    int count=0;
    unsigned short sumInsertsize=0;
    bool flag = false;
    std::string filename = argv[1], chr = argv[2];
    BamTools::BamReader reader;

    std::stringstream(argv[3]) >> uRstart;
    std::stringstream(argv[4]) >> uWall;
    std::stringstream(argv[5]) >> uLen;
    std::stringstream(argv[6]) >> uMean;
    uRend = uWall + uLen;

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
    BamTools::BamRegion region(RefID, uRstart, RefID, uWall);

    reader.SetRegion(region);

    BamTools::BamAlignment al;
    //do the job
    while ( reader.GetNextAlignment(al) ){
        if (al.MatePosition > uWall && al.MatePosition < uRend){
            sumInsertsize+=abs(al.InsertSize)-uGap;
            count++;
        }
    }
    std::cout << std::fixed << sumInsertsize/count << endl;
}
