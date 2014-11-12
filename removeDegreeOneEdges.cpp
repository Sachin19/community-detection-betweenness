#include <cstdio>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
using namespace std;
int main(int argc, char* argv[]){
	ifstream fin;
	ofstream fout;
	vector< vector<int> > G;
	vector< vector<int> > G2;

	G.resize(390680);
	G2.resize(390680);
	int x,y;
	for (int i = 0; i < argc; i+=2)
	{
		fin.open(argv[i+1],ios::in);
		while(fin){
			fin>>x>>y;
			G[x].push_back(y);
			//G[y].push_back(x);
		}
		fin.close();
		fin.open(argv[i+2],ios::in);
		while(fin){
			fin>>x>>y;
			G2[x].push_back(y);
			//G[y].push_back(x);
		}
		fin.close();
		vector<int> marked(G.size(), 0);
		for (int i = 0; i < G.size(); ++i)
		{
			if(G[i].size() + G2[i].size() <= 1) {
				marked[i]=1;
			}
		}
		fin.open(argv[i+1],ios::in);
		while(fin){
			fin>>x>>y;
			if(marked[x] == 1 || marked[y] == 1) {
				marked[x] = 1;
				marked[y] = 1;
			}
			//G[y].push_back(x);
		}
		fin.close();

		fin.open(argv[i+2],ios::in);
		while(fin){
			fin>>x>>y;
			if(marked[x] == 1 || marked[y] == 1) {
				marked[x] = 1;
				marked[y] = 1;
			}
			//G[y].push_back(x);
		}
		fin.close();

		fin.open(argv[i+1],ios::in);
		char newname[100];
		sprintf(newname, "%s_small.txt",argv[i+1]);
		fout.open(newname,ios::out);
		while(fin){
			fin>>x>>y;
			if(marked[x] == 0 && marked[y] == 0) {
				fout << x << " " << y << endl;
			}
			//G[y].push_back(x);
		}
		fout.close();
		fin.close();

		fin.open(argv[i+2],ios::in);
		sprintf(newname, "%s_small.txt",argv[i+2]);
		fout.open(newname,ios::out);
		while(fin){
			fin>>x>>y;
			if(marked[x] == 0 && marked[y] == 0) {
				fout << x << " " << y << endl;
			}
			//G[y].push_back(x);
		}
		fout.close();
		fin.close();
		// char newname[100];
		// sprintf(newname, "%s_adj.txt",argv[i+1]);
		// fout.open(newname,ios::out);
		// for (int i = 0; i < 390680; ++i)
		// {
		// 	fout<<i<<" ";
		// 	for (int j = 0; j < G[i].size(); ++j)
		// 	{
		// 		fout<<G[i][j]<<" ";
		// 	}
		// 	fout<<"\n";
		// }
		// fout.close();
	}
}