#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include "api/api_global.h"
#include "api/BamReader.h"

using namespace std;
using namespace BamTools;

float stdev(int32_t *arrInserts){
    int32_t len = sizeof(*arrInserts)/sizeof(int32_t), sum;
    int i;
    float mean, sumstd;
    for(i=0; i<len; i++)
        sum+=arrInserts[i];
    mean = sum/len;
    for(i=0; i<len; i++)
        sumstd+=(arrInserts[i]-mean)*(arrInserts[i]-mean);
    return sqrt(sumstd/len)

int main(int argc, char *argv[])
{
    //unsigned short count;
    int32_t uRstart, uLen, uRend, uGap, uMean;
    int count=0;
    unsigned short sumInsertsize=0;
    std::string filename = argv[1], chr = argv[2];
    BamTools::BamReader reader;

    std::stringstream(argv[3]) >> uRstart;
    std::stringstream(argv[4]) >> uLen;
    uRend = uRstart + uLen;

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
    BamTools::BamRegion region(RefID, uRstart, RefID, uRend);

    reader.SetRegion(region);

    BamTools::BamAlignment al;
    //do the job
    while ( reader.GetNextAlignment(al) ){
		if (al.MatePosition > uRend){

			sumInsertsize+=abs(al.InsertSize)-uGap;
            count++;
        }
	}
	std::cout << std::fixed << sumInsertsize/count;
}

