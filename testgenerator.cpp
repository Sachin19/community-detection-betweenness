#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <map>
#include <set>
#include <iostream>
#include <cstring>
using namespace std;

int G[1000][1000];

int main(int argc, char* argv[]){
	int t, n, len, l, x,y, start, end,j, maxe, edgenum;

	char ch;
	char s[60];
	string ee;
	map<int, string> mp;
	set<string> st;
	set<string> nodes;
	set<pair<int, int> > edges;
	t = 10;
	n = atoi(argv[1]);
	//int L = atoi(argv[2]);
	srand((unsigned int)time(NULL));
	//printf("%d\n",t);
	
	st.clear();
	edges.clear();
					
	maxe = (n*(n-1))/2;
	//printf("maxe %d\n", maxe);

	edgenum = rand()%maxe+1;
	//cout << edgenum << endl;
	//cout << mp[start] <<" "<<mp[end] << endl;
	//printf("edgenum %d\n",edgenum);
	int ll = 0;
	while(ll < edgenum)
	{
		do {
			x = rand()%n;
			y = rand()%n;
		}while(edges.find(make_pair(x,y)) != edges.end() || edges.find(make_pair(y,x)) != edges.end() || x==y);
		
		edges.insert(make_pair(x,y));
		edges.insert(make_pair(y,x));
		printf("%d %d\n",x,y);
		printf("%d %d\n",y,x);
		ll++;
	}
}
