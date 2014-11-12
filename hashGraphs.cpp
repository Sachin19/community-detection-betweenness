#include <iostream>
#include <fstream>
#include <map>
using namespace std;
int main(int argc, char* argv[]){
	map<int,int> m;
	ifstream fin;
	ofstream fout;
	int i=0,a,b;
	int ind=0;
	while(i+1<argc){
		printf("%d %s\n",i+1, argv[i+1]);
		char write[100];
		sprintf(write,"%s_ind.txt",argv[i+1]);
		fin.open(argv[i+1],ios::in);
		fout.open(write,ios::out);
		while(fin){
			fin>>a>>b;
			if(m.count(a) == 0){
				m[a] = ind;
				ind++;
			}
			if(m.count(b) == 0){
				m[b] = ind;
				ind++;
			}
			fout<<m[a]<<" "<<m[b]<<"\n";
		}
		i++;
		fin.close();
		fout.close();
	}
	cout << ind << "\n";
}