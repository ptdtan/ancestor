#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include "api/api_global.h"
#include "api/BamReader.h"

using namespace std;
using namespace BamTools;

struct node{
	int inIdx;
	int outIdx;
  //double win, wout;
};

void initialize(std::vector<node> *nodes){
	for(std::vector<node>::iterator it=nodes->begin(); it!=nodes->end(); it++){
		it->inIdx = -1;
		it->outIdx = -1;
	}
}

int main(int argc, char *argv[])
{
    //unsigned short count;
    int count=0, RefID;
	int32_t currPos, currMatePos;
	int thresHold;
    //std::vector<node> nodes;
	const BamTools::RefVector *references;
    string filename = argv[1];
	std::stringstream(argv[2]) >> thresHold;
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
    references = &reader.GetReferenceData();
	std::vector<int> counter(references->size(), 0);
	std::vector<node> nodes(references->size());
	//initialize
	initialize(&nodes);
	//for(int i=0; i<counter.size(); i++)
	//	cout << counter[i] << endl;
	//for(int i=0; i<nodes.size(); i++)
	//	cout << nodes[i].inIdx << endl;
    //cout << (*references).size() << endl;
    //Get all reference sequence name
	for(std::vector<BamTools::RefData>::const_iterator it=references->begin(); it!=references->end(); it++){
		count = 0;
		cout << "Infering: " << it->RefName << endl;
		RefID = reader.GetReferenceID(it->RefName);
		//cout << RefID << endl;
		BamTools::BamRegion region(RefID,0,RefID,it->RefLength);
		reader.SetRegion(region);
		BamTools::BamAlignment al;

		//go right
		while (reader.GetNextAlignment(al)){
			cout << al.RefID << endl;
			currPos = al.Position;
			if (count > thresHold){
				cout << (*references)[al.RefID].RefName << endl;
				cout << count << endl;
				break;
			}
			if (al.AlignmentFlag == 97 && al.MateRefID != RefID){
				currMatePos = al.MatePosition;
				//cout << al.MateRefID << endl;
				reader.Jump(al.MateRefID, al.MatePosition);
				while (reader.GetNextAlignment(al)){
					if (al.Position == currMatePos && al.AlignmentFlag == 145){
					//cout <<al.Position << endl;
					cout << (*references)[al.RefID].RefName << endl;
					//reader.Jump(RefID, currPos);
					reader.Jump(RefID, currPos);
					count++;
					break;
					}
				}
			}
		}
	}
    /*int RefID = reader.GetReferenceID(chr);
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
		//std::cout << std::fixed << sumInsertsize/count;*/
}
