#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include "api/api_global.h"
#include "api/BamReader.h"

using namespace std;
using namespace BamTools;

struct node{
	int inIdx, outIdx;
};

void initialize_nodes(std::vector<node> *nodes){
	for(std::vector<node>::iterator it=nodes->begin(); it!=nodes->end(); it++){
		it->inIdx = -1;
		it->outIdx = -1;
	}
}

int main(int argc, char *argv[])
{
    //unsigned short count;
    int RefSize, RefID;
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
	RefSize = references->size();
	std::vector<node> nodes(RefSize);
	initialize_nodes(&nodes);
	std::vector<uint32_t> inCounter(RefSize,0), outCounter(RefSize, 0);
	//for(int i=0; i<counter.size(); i++)
	//	cout << counter[i] << endl;
	//for(int i=0; i<nodes.size(); i++)
	//	cout << nodes[i].inIdx << endl;
    //cout << (*references).size() << endl;
    //Get all reference sequence name
	for(std::vector<BamTools::RefData>::const_iterator it=references->begin(); it!=references->end(); it++){
		std::fill(inCounter.begin(), inCounter.end(), 0);
		std::fill(outCounter.begin(), outCounter.end(), 0);
		cout << "Infering: " << it->RefName << endl;
		RefID = reader.GetReferenceID(it->RefName);
		BamTools::BamRegion region(RefID,0,RefID,it->RefLength);
		reader.SetRegion(region);
		BamTools::BamAlignment al;

		while (reader.GetNextAlignment(al)){
			//cout << al.Length << endl;
			if (nodes[RefID].inIdx != -1 && nodes[RefID].outIdx != -1) break;
			currPos = al.Position;
			//cout << currPos << endl;
			//go right
			if (al.AlignmentFlag == 97 && al.MateRefID != RefID &&al.MapQuality > 50 && nodes[RefID].outIdx == -1){
				currMatePos = al.MatePosition;
				//cout << al.MateRefID << endl;
				reader.Jump(al.MateRefID, al.MatePosition);
				while (reader.GetNextAlignment(al)){
					if (al.Position == currMatePos && al.AlignmentFlag == 145){
					//cout <<al.Position << endl;
					//cout << (*references)[al.RefID].RefName << endl;
						outCounter[al.RefID]++;
						if (outCounter[al.RefID] > thresHold) {
							nodes[RefID].outIdx = al.RefID;
							//cout << RefID <<;
							cout << (*references)[al.RefID].RefName << endl;
						}
					reader.SetRegion(RefID, currPos, RefID, it->RefLength);
					break;
					}
				}
			}
			//go left
			if (al.AlignmentFlag == 145 && al.MateRefID != RefID &&al.MapQuality > 50 && nodes[RefID].inIdx == -1){
				currMatePos = al.MatePosition;
				reader.Jump(al.MateRefID, al.MatePosition-al.Length-50);
				while (reader.GetNextAlignment(al)){
					if (al.Position == currMatePos && al.AlignmentFlag == 97){
						inCounter[al.RefID]++;
						if (inCounter[al.RefID] > thresHold){
							nodes[RefID].inIdx = al.RefID;
							cout << (*references)[al.RefID].RefName << endl;
						}
						reader.SetRegion(RefID, currPos, RefID, it->RefLength);
						break;
					}
				}
			}

		}
	}
	cout << "\nAdjacency" << endl;
	for(int i = 0; i < nodes.size(); i++){
		if ( nodes[i].inIdx != -1)
			cout << (*references)[nodes[i].inIdx].RefName << "\t" <<  (*references)[i].RefName << endl;
		if ( nodes[i].outIdx != -1)
			cout << (*references)[i].RefName << "\t" << (*references)[nodes[i].outIdx].RefName << endl;
			
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

