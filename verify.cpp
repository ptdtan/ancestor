#include <iostream>
#include <sstream>
#include "api/api_global.h"
#include "api/BamReader.h"

using namespace std;
using namespace BamTools;

int main(int argc, char *argv[])
{
    unsigned short count;
    //std::string filename, idxname;
    int32_t uGstart, uGend;
    std::string filename = argv[1], chr = argv[2];
    BamTools::BamReader reader;

    std::stringstream str1(argv[3]);
    str1 >> uGstart;
    std::stringstream str2(argv[4]);
    str2 >> uGend;

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
    BamTools::BamRegion region(RefID, uGstart, RefID, uGend);

    reader.SetRegion(region);

    BamTools::BamAlignment al;
    //do the job
    while ( reader.GetNextAlignment(al) ){
		if (al.isProperPair())
        std::cout << abs(al.InsertSize) << "\n";
