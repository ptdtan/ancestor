#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <map>
#include <string>
#include <utility>
#include <cmath>
#include "api/api_global.h"
#include "api/BamReader.h"

using namespace std;
using namespace BamTools;

float *stdev(int32_t *arrInserts, int len){
    //int32_t len = sizeof(*arrInserts)/sizeof(int32_t), sum;
	int32_t sum=0;
    //cout << len;
	int i;
    float *ret, mean, sumstd=0.0;
    ret = (float *)calloc(2,sizeof(float));
	for(i=0; i<len; i++){
		//cout << arrInserts[i] << endl;
        sum+=arrInserts[i];
	}
    mean = sum/len;
	ret[0] =mean;
    for(i=0; i<len; i++)
        sumstd+=(float)(arrInserts[i]-mean)*(float)(arrInserts[i]-mean);
	ret[1] = sqrt(sumstd/len);
    return ret;
}
int main(int argc, char *argv[])
{
    //unsigned short count;
    int32_t uRstart, uLen, uRend, uGap, uMean, *arrInserts;
    int count=0;
    unsigned short sumInsertsize=0;
	float *info;
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
    arrInserts = (int32_t *)calloc(1, sizeof(int32_t));
    while ( reader.GetNextAlignment(al) ){
		if (al.IsProperPair()){
			arrInserts[count]=(int32_t)abs(al.InsertSize);
            count++;
            arrInserts = (int32_t *)realloc(arrInserts, (count+1)*sizeof(int32_t));
        }
	}
	info = stdev(arrInserts, count);
	cout << "Mean: " << info[0] << "\t" << "Std: " << info[1] << "\t" << count << endl;
}
