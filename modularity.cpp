#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <cstring>
#include <fstream>
#include <cstdio>

#define omega 0.5
#define gamma 0.5

using namespace std;

vector<int> clusterMapping(vector< vector<int> > clusters, int n) {
	int i,j;
	vector<int> clusterMap(n);
	
	for(i=0;i<clusters.size();i++)
	{
		for(j=0;j<clusters[i].size();j++)
		{
			clusterMap[clusters[i][j]] = i;
		}
	}
	
	return clusterMap;
}

int numSlices(vector< pair<int, int> > edgeList, int n) {
	int s = 0;
	for(int i=0; i<edgeList.size(); i++) {
		if(s < edgeList[i].first)
			s = edgeList[i].first;
		if(s < edgeList[i].second)
			s = edgeList[i].second;
	}
	return s/n;
}
double quality(vector< pair<int, int> > edgeList, int num, vector<int> clusterMap, int numSlice) {
	
	vector< vector<int> > K(num, vector<int> (numSlice, 0));
	vector< vector<int> > C(num, vector<int> (numSlice, 0));
	vector<int> m(numSlice, 0);
	double mu2=0;
	for( int i=0; i < edgeList.size(); i++){

		if(edgeList[i].first%num == edgeList[i].second%num){

			C[edgeList[i].first%num][edgeList[i].second/num]+=omega;
			C[edgeList[i].second%num][edgeList[i].first/num]+=omega;
			mu2 += 2*omega;
		}
		else{
			K[edgeList[i].first%num][edgeList[i].second/num]++;
			K[edgeList[i].second%num][edgeList[i].first/num]++;
			m[edgeList[i].first/num]+=2;
			mu2 += 2;
		}
	}

	
	double Q=0, x;

	for( int i=0; i < edgeList.size(); i++) {	
		if(clusterMap[edgeList[i].first%num] == clusterMap[edgeList[i].second%num]) {
			if(edgeList[i].first/num == edgeList[i].second/num) {
				x = (1-(gamma*K[edgeList[i].first%num][edgeList[i].first/num]*K[edgeList[i].second%num][edgeList[i].first/num])/(1.0*m[edgeList[i].first/num]));
				Q = Q + 2*x;
			}
			else{
				Q += omega;
			}
		}
	}

	//cout << "Q = " << Q << " twomu = " << mu2 << "\n";
	Q /= mu2;
	
	return Q;

}

vector< pair<int, int> > readEdges(int numSlice, char* slices[], int numNodes) {

	ifstream fin;
	int x, y;
	vector< pair<int, int> > edgeList;
	for (int i = 0; i < numSlice; ++i)
	{	
		int M[150][150];
		for (int i1 = 0; i1 < numNodes; ++i1) {
			for (int j1 = 0; j1 < numNodes; ++j1) {
				M[i1][j1]=0;
			}
		}
		int c=0;
		printf("%s\n",slices[i]);
		fin.open(slices[i],ios::in);
		c=0;
		while(fin) {
			fin>> x >> y;
			if(M[x][y]!=1){
				edgeList.push_back(make_pair(x+numNodes*i,y+numNodes*i));				
				c++;
			}
			M[x][y]=1;
			M[y][x]=1;
		}
		fin.close();
	}
	//cout<<"before: "<<2*edgeList.size() + numSlice*(numSlice-1)*numNodes*0.5<<"\n";
	for (int i = 0; i < numSlice; ++i)
	{
		for (int j = i+1; j < numSlice; ++j)
		{
			for (int k = 0; k < numNodes; ++k)
			{
				edgeList.push_back(make_pair(k+numNodes*i,k+j*numNodes));
			}
			
		}
		
	}
	//cout<<2*edgeList.size()<<"\n";
	return edgeList;
}

vector<int> readClusters(char* clusterFile) {
	ifstream fin;
	fin.open(clusterFile,ios::in);

	vector<int> clusterMap;
	int clustnum;

	while(fin) {
		fin>>clustnum;
		clusterMap.push_back(clustnum);
	}
	fin.close();
	// clusters.erase(clusters.begin() + clusters.size() - 1);
	return clusterMap;
}

int main(int argc, char* argv[])
{
	vector< pair<int,int> > edgeList;// = readEdges(argv[1]);
	vector<int> clusters;//readClusters(argv[2]);

	if(argc < 6){
		printf("How to run:\nmodularity.exe <numNodes> <numClusters> <clusterFilePath> <numSlice> <slice1EdgeList> [<slice2EdgeList> ... ] ");
		exit(-1);
	}

	int numNodes = atoi(argv[1]);
	int numCluster = atoi(argv[2]);
	clusters = readClusters(argv[3]);
	int numSlice = atoi(argv[4]);
	//printf("numSlices: %d\n",numSlice);
	edgeList = readEdges(numSlice,argv+5,numNodes);
	
	cout<<"\n\nmodularity = "<<quality(edgeList,numNodes,clusters,numSlice);



}
