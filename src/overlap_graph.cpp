#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include "api/api_global.h"
#include "api/BamReader.h"

using namespace std;
using namespace BamTools;

struct node{
  char *in;
  char *out;
  int weight;
}

int main(int argc, char *argv[])
{
    //unsigned short count;
    int32_t uStart, uEnd
    int count=0;
    std::string filename = argv[1]
    BamTools::BamReader reader;

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
    //Get all reference sequence name
    int RefID = reader.GetReferenceID(chr);
    BamTools::BamRegion region(RefID, uRstart, RefID, uWall);

    reader.SetRegion(region);

    BamTools::BamAlignment al;
    //do the job
    while ( reader.GetNextAlignment(al) ){
		if (al.MatePosition > uWall && al.MatePosition < uRend){
            //cout  << (unsigned long)al.Position << "\t" <<(unsigned long)al.MatePosition << endl;
			//sumInsertsize+=abs(al.InsertSize)-uGap;
            //count++;
			cout << "True" << endl;
			flag = true;
			break;
        }
	}
	if (!flag)
		cout << "False" << endl;
		//std::cout << std::fixed << sumInsertsize/count;
}
