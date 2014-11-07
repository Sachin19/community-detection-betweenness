#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <cstring>
#include <stack>
#include <set>
#include <fstream>
#include <queue>

#define gamma 0.5
#define threshold 0.65

using namespace std;

struct component{
	double dleft;
	double dright;
	vector<int> leftComp;
	vector<int> rightComp;
};


class Graph
{
	private:
		int V;
		int L;
		int N;
		vector< vector<int> > adj;
		vector< bool > deleted;
	public:
		Graph(int v, int l) {
			V = v*l;
			L = l;
			N = v;
			adj.resize(V);
			deleted.resize(V);
			for(int i=0; i<V; i++) {
				deleted[i] = false;
			}
			//cout<<deleted.size()<<"\n";
		}
		Graph(const Graph &g) {
			this->V = g.V;
			this->adj = g.adj;
		}
		void addEdge(pair<int,int> edge) {
			adj[edge.first].push_back(edge.second);
			adj[edge.second].push_back(edge.first);
		}
		void printGraph() {
			for(int i=0;i<V;i++) {
				if(isDeleted(i)){
					continue;
				}
				cout << i <<": ";
				for (int j = 0; j < adj[i].size(); j++) {
					cout<<adj[i][j]<<" ";
				}
				cout<<"\n";
			}
		}
		vector<int> getNeighbours(int i) {
			return adj[i];
		}
		bool isEmpty(){
			//cout<<V;
			for (int i = 0; i < V; ++i)
			{
				//cout<<deleted[i];
				if(!deleted[i]){

					return false;
				}
			}
			return true;
		}
		bool isDeleted(int i){
			return deleted[i];
		}
		double getDegree(int i) {
			return (double) adj[i].size();
		}
		void removeEdge(pair<int,int> edge) {
			adj[edge.first].erase(remove(adj[edge.first].begin(), adj[edge.first].end(), edge.second), adj[edge.first].end());
			adj[edge.second].erase(remove(adj[edge.second].begin(), adj[edge.second].end(), edge.first), adj[edge.second].end());
		}
		bool isEdgePresent(pair<int,int> edge) {
			//cout << "looking for edge between : "<< edge.first<<","<<edge.second<<"\n";
			std::vector<int>::iterator it,it2;
			it = find(adj[edge.first].begin(),adj[edge.first].end(),edge.second);
			//it2 = find(adj[edge.second].begin(),adj[edge.second].end(),edge.first);
			if(it!=adj[edge.first].end()  ) {
				return true;
			}
			else {
				return false;
			}
		}
		int numNodes() {
			int ans = 0;

			for (int i = 0; i < N*L; ++i) {
				if(!deleted[i])
				ans++;
			}
			
			return ans;
		}
		int numEdges() {
			int edges = 0;
			
			for(int i=0;i<adj.size();i++)
				if(!deleted[i])
					edges+=adj[i].size();
			
			return edges/2;
		}

		void removeSubgraphFromAllLevels(vector<int> nodes) {
			//cout << "Removing : " << nodes.size() << "\n";
			for(int l = 0; l < L; l++) {				
				for(int i = 0; i < adj.size(); i++) {
					for(int j = 0; j < nodes.size(); j++) {
						adj[i].erase(remove(adj[i].begin(), adj[i].end(), nodes[i]+l*N), adj[i].end());
					}
				}
				for(int i = 0; i < nodes.size(); i++) {
					//cout << "deleting " << nodes[i]<<" "<<nodes[i]%N + l*N<<"\n";
					deleted[nodes[i]+l*N] = true;
				}
			}			
		}
};

vector<double> betweennessCentralityNodes(Graph &g, int n, int l) {


	vector<double>Cb(n,0.0);

	
	
	double *sigma = new double[n*l];
	double *sigmaM = new double[n];
	int *D = new int[n*l];
	int *Dm = new int[n*l];
	vector< vector<int> > Vsame(n);
	


	for (int i = 0; i < n; i++)
	{
		if(g.isDeleted(i)){
			continue;
		}
		if(i%100 == 0)
			cout << i <<"\n";
		memset(sigma, 0, n*l*sizeof(double));
		memset(sigmaM, 0, n*sizeof(double));
		memset(D, -1, n*l*sizeof(int));
		memset(Dm, -1, n*l*sizeof(int));
		Vsame.clear();
		queue<int> Q;
		stack<int>S;
		vector< stack<int> > P(n*l);


		// fill(P.begin(), P.end(), stack<int>());

		for (int j = i; j < n*l; j += n)
		{
			sigma[j] = 1.0;
		}

		for (int j = i; j < n*l; j += n)
		{
			D[j] = 0;
		}

		for (int j = i; j < n*l; j += n)
		{
			Dm[j] = 0;
		}

		//fill(Vsame.begin(),Vsame.end(),vector<int>());
		Q.push(i);

		while(!Q.empty()) {
			int v = Q.front();
			Q.pop();
			if(g.isDeleted(v))
				continue;	

			S.push(v);

			if(v!=i) {
				vector<int> vNeighbours = g.getNeighbours(v);
				for(int j=0;j<vNeighbours.size();j++) {
					//W.insert(vNeighbours[j]);

					int w = vNeighbours[j];
					if(D[w] < 0) {
						Q.push(w);
						D[w] = D[v] + 1;

						if(Dm[w%n] < 0 || Dm[w%n] == D[w]) {
							Dm[w%n] = D[w];
							Vsame[w%n].push_back(w);
						}
					}
					if(D[w] == D[v] + 1) {
						sigma[w] += sigma[v];
						P[w].push(v);
					}

				}
			}
			else {
				for (int j = i; j < n*l; j += n)
				{
					vector<int> jNeighbours = g.getNeighbours(j);
					for(int k=0;k<jNeighbours.size();k++) {
						//W.insert(jNeighbours[k]);

						int w = jNeighbours[k];
						if(D[w] < 0) {
							Q.push(w);
							D[w] = D[v] + 1;

							if(Dm[w%n] < 0 || Dm[w%n] == D[w]) {
								Dm[w%n] = D[w];
								Vsame[w%n].push_back(w);
							}
						}
						if(D[w] == D[v] + 1) {
							sigma[w] += sigma[v];
							P[w].push(v);
						}
					}
				}
			}
		}

		// cout << "donebfs\n";
		// for(int w =0; w<n; w++) {
		// 	sigmaM[w] = 0.0;
		// 	for(int j=0; j<Vsame[w].size(); j++) {
		// 		int v = Vsame[w][j];
		// 		sigmaM[w] += sigma[v];
		// 	}
		// }

		// vector<double> delta(n*l,0.0);		

		// while(!S.empty()) {
		// 	int w = S.top();
		// 	S.pop();

		// 	while(!P[w].empty()) {
		// 		int v = P[w].top();
		// 		P[w].pop();

		// 		vector<int>::iterator p;
		// 		p = find(Vsame[w%n].begin(),Vsame[w%n].end(),w);

		// 		if(p != Vsame[w%n].end()) {
		// 			if(sigma[w] > 0 && sigmaM[w%n]>0)
		// 				delta[v] += (1-gamma)*(sigma[v]/sigma[w])*((sigma[w]/sigmaM[w%n]) + delta[w]);
		// 		}
		// 		else {
		// 			if(sigma[w] > 0 && sigmaM[w%n]>0)
		// 				delta[v] += gamma*(sigma[v]/sigma[w])*delta[w];
		// 		}
		// 		if(w!=i) {
		// 			Cb[w%n] += delta[w];
		// 		}
		// 	}
		// }
	}

	return Cb;
}

int main(int argc, char const *argv[]) {
	int numNodes, numSlice;
    numNodes = atoi(argv[1]);
    numSlice = atoi(argv[2]);

    char slices[100][100];
    for (int i = 0; i < numSlice; i++)
    {
        strcpy(slices[i], argv[i + 3]);
    }

    ifstream fin;
    int x, y;
    vector< pair<int, int> > edgeList;
    Graph g(numNodes, numSlice);
    Graph singleSliceGraph(numNodes, 1);

    for (int i = 0; i < numSlice; ++i)
    {
        fin.open(slices[i], ios::in);

        while (fin)
        {
            fin >> x >> y;

            if (!g.isEdgePresent(make_pair(x + numNodes * i, y + numNodes * i)) && !g.isEdgePresent(make_pair(y + numNodes * i, x + numNodes * i)))
                g.addEdge(make_pair(x + numNodes * i, y + numNodes * i));

            if (!singleSliceGraph.isEdgePresent(make_pair(x, y)) && !singleSliceGraph.isEdgePresent(make_pair(y, x)))
            {
                singleSliceGraph.addEdge(make_pair(x, y));
            }
        }

        fin.close();
    }
	
    vector<double> Cb = betweennessCentralityNodes(singleSliceGraph,numNodes,1);

    for(int i=0;i<Cb.size();i++)
    	cout<<Cb[i]<<"\n";
}