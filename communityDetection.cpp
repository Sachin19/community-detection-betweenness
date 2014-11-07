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
#define threshold 0.2

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
						adj[i].erase(remove(adj[i].begin(), adj[i].end(), nodes[j]+l*N), adj[i].end());
					}
				}
				for(int i = 0; i < nodes.size(); i++) {
					//cout << "deleting " << nodes[i]<<" "<<nodes[i]%N + l*N<<"\n";
					deleted[nodes[i]+l*N] = true;
				}
			}			
		}
};

// vector<double> betweennessCentralityNodes(Graph &g, int n, int l) {


// 	vector<double>Cb(n,0.0);

	
	
// 	double *sigma = new double[n*l];
// 	double *sigmaM = new double[n];
// 	int *D = new int[n*l];
// 	int *Dm = new int[n*l];
// 	vector< vector<int> > Vsame(n);
	


// 	for (int i = 0; i < n; i++)
// 	{
// 		if(g.isDeleted(i)){
// 			continue;
// 		}
// 		if(i%100 == 0)
// 			cout << i <<"\n";
// 		memset(sigma, 0, n*l*sizeof(double));
// 		memset(sigmaM, 0, n*sizeof(double));
// 		memset(D, -1, n*l*sizeof(int));
// 		memset(Dm, -1, n*l*sizeof(int));
// 		Vsame.clear();
// 		queue<int> Q;
// 		stack<int>S;
// 		vector< stack<int> > P(n*l);


// 		// fill(P.begin(), P.end(), stack<int>());

// 		for (int j = i; j < n*l; j += n)
// 		{
// 			sigma[j] = 1.0;
// 		}

// 		for (int j = i; j < n*l; j += n)
// 		{
// 			D[j] = 0;
// 		}

// 		for (int j = i; j < n*l; j += n)
// 		{
// 			Dm[j] = 0;
// 		}

// 		//fill(Vsame.begin(),Vsame.end(),vector<int>());
// 		Q.push(i);

// 		while(!Q.empty()) {
// 			int v = Q.front();
// 			Q.pop();
// 			if(g.isDeleted(v))
// 				continue;	

// 			S.push(v);

// 			if(v!=i) {
// 				vector<int> vNeighbours = g.getNeighbours(v);
// 				for(int j=0;j<vNeighbours.size();j++) {
// 					//W.insert(vNeighbours[j]);

// 					int w = vNeighbours[j];
// 					if(D[w] < 0) {
// 						Q.push(w);
// 						D[w] = D[v] + 1;

// 						if(Dm[w%n] < 0 || Dm[w%n] == D[w]) {
// 							Dm[w%n] = D[w];
// 							Vsame[w%n].push_back(w);
// 						}
// 					}
// 					if(D[w] == D[v] + 1) {
// 						sigma[w] += sigma[v];
// 						P[w].push(v);
// 					}

// 				}
// 			}
// 			else {
// 				for (int j = i; j < n*l; j += n)
// 				{
// 					vector<int> jNeighbours = g.getNeighbours(j);
// 					for(int k=0;k<jNeighbours.size();k++) {
// 						//W.insert(jNeighbours[k]);

// 						int w = jNeighbours[k];
// 						if(D[w] < 0) {
// 							Q.push(w);
// 							D[w] = D[v] + 1;

// 							if(Dm[w%n] < 0 || Dm[w%n] == D[w]) {
// 								Dm[w%n] = D[w];
// 								Vsame[w%n].push_back(w);
// 							}
// 						}
// 						if(D[w] == D[v] + 1) {
// 							sigma[w] += sigma[v];
// 							P[w].push(v);
// 						}
// 					}
// 				}
// 			}
// 		}

// 		for(int w =0; w<n; w++) {
// 			sigmaM[w] = 0.0;
// 			for(int j=0; j<Vsame[w].size(); j++) {
// 				int v = Vsame[w][j];
// 				sigmaM[w] += sigma[v];
// 			}
// 		}

// 		vector<double> delta(n*l,0.0);		

// 		while(!S.empty()) {
// 			int w = S.top();
// 			S.pop();

// 			while(!P[w].empty()) {
// 				int v = P[w].top();
// 				P[w].pop();

// 				vector<int>::iterator p;
// 				p = find(Vsame[w%n].begin(),Vsame[w%n].end(),w);

// 				if(p != Vsame[w%n].end()) {
// 					if(sigma[w] > 0 && sigmaM[w%n]>0)
// 						delta[v] += (1-gamma)*(sigma[v]/sigma[w])*((sigma[w]/sigmaM[w%n]) + delta[w]);
// 				}
// 				else {
// 					if(sigma[w] > 0 && sigmaM[w%n]>0)
// 						delta[v] += gamma*(sigma[v]/sigma[w])*delta[w];
// 				}
// 				if(w!=i) {
// 					Cb[w%n] += delta[w];
// 				}
// 			}
// 		}
// 	}

// 	return Cb;
// }

//function improved. working atleast 100 times faster. \m/
vector<double> betweennessCentralityNodes(Graph &g, int n, int l) {
    vector<double> Cb(n,0.0);

    double *sigma = new double[n*l];
    double *sigmaM = new double[n];
    int *D = new int[n*l];
    int *Dm = new int[n*l];
    vector< set<int> > Vsame(n);
    
    for (int i = 0; i < n; i++)
    {
        if(g.isDeleted(i)){
            continue;
        }

        memset(sigma, 0, n*l*sizeof(double));
        memset(sigmaM, 0, n*sizeof(double));
        memset(D, -1, n*l*sizeof(int));
        memset(Dm, -1, n*l*sizeof(int));
        queue<int> Q;
        stack<int>S;
        Vsame.resize(n, set<int>());
        map<int, vector<int> > P;

        Q.push(i);

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
                    if(P.count(w) == 0) {
                        P[w] = vector<int>();
                    }
                    if(D[w] < 0) {
                        Q.push(w);
                        D[w] = D[v] + 1;

                        if(Dm[w%n] < 0 || Dm[w%n] == D[w]) {
                            Dm[w%n] = D[w];
                            Vsame[w%n].insert(w);
                        }
                    }
                    if(D[w] == D[v] + 1) {
                        sigma[w] += sigma[v];
                        sigmaM[w%n] += sigma[v];
                        P[w].push_back(v);
                    }

                }
            }
            else {
                for (int j = i; j < n*l; j += n)
                {
                    vector<int> jNeighbours = g.getNeighbours(j);
                    for(int k=0;k<jNeighbours.size();k++) {

                        int w = jNeighbours[k];
                        if(P.count(w) == 0) {
                            P[w] = vector<int>();
                        }
                        if(D[w] < 0) {
                            Q.push(w);
                            D[w] = D[v] + 1;

                            if(Dm[w%n] < 0 || Dm[w%n] == D[w]) {
                                Dm[w%n] = D[w];
                                Vsame[w%n].insert(w);
                            }
                        }
                        if(D[w] == D[v] + 1) {
                            sigma[w] += sigma[v];
                            sigmaM[w%n] += sigma[v];
                            P[w].push_back(v);
                        }
                    }
                }
            }
        }

        //following piece of code is not required anymore.
        // cout << "donebfs\n";
        // for(int w =0; w<n; w++) {
        //  sigmaM[w] = 0.0;
        //  for(int j=0; j<Vsame[w].size(); j++) {
        //      int v = Vsame[w][j];
        //      sigmaM[w] += sigma[v];
        //  }
        // }
        //useless code ends here

        vector<double> delta(n*l,0.0);       

        while(!S.empty()) {
            int w = S.top();
            S.pop();

            for( int ll = 0; ll < P[w].size(); ll++) {
                int v = P[w][ll];

                if(Vsame[w%n].find(w) != Vsame[w%n].end()) {
                    if(sigma[w] > 0 && sigmaM[w%n]>0)
                        delta[v] += (1-gamma)*(sigma[v]/sigma[w])*((sigma[w]/sigmaM[w%n]) + delta[w]);
                }
                else {
                    if(sigma[w] > 0 && sigmaM[w%n]>0)
                        delta[v] += gamma*(sigma[v]/sigma[w])*delta[w];
                }
                if(w!=i) {
                    Cb[w%n] += delta[w];
                }
            }
        }
    }

    return Cb;
}

int commonElements(vector<int> a, vector<int> b) {
	set<int> inter;
	for (int i = 0; i < a.size(); ++i)
	{
		inter.insert(a[i]);
	}
	int c=0;
	for (int i = 0; i < b.size(); ++i)
	{
		if(inter.find(b[i]) != inter.end()) {
			c++;
		}
	}
	return c;
}

pair<int, int> topBetweennessCentralityInEdges(Graph &g, vector<double> Cb, int v, int l)
{
    int n = v * l;

    double maxCentrality = 0.0;
    double centrality;
    int e1 = -1, e2 = -1;
    int cnbrs;
    for (int i = 0; i < n; i++)
    {
        //cout<<"ta\n";
        if (g.isDeleted(i))
            continue;
        //printf("hello\n");
        vector<int> iNeighbours = g.getNeighbours(i);
        //cout << "size: "<<iNeighbours.size()<<"\n";
        for (int j = 0; j < iNeighbours.size(); j++)
        {
            if (g.isDeleted(iNeighbours[j]))
                continue;
            vector<int> jNeighbours = g.getNeighbours(iNeighbours[j]);

            //vector<int> commonNeighbours(n);
            // vector<int> kNeighbours = iNeighbours;
            // sort(kNeighbours.begin(), kNeighbours.end());
            // sort(jNeighbours.begin(), jNeighbours.end());
            
            // std::vector<int>::iterator it;
            // it = std::set_intersection (kNeighbours.begin(), kNeighbours.end(), jNeighbours.begin(),
            //                             jNeighbours.end(), commonNeighbours.begin());

            //commonNeighbours.resize(it - commonNeighbours.begin());
            cnbrs = commonElements(iNeighbours, jNeighbours);

            double di = g.getDegree(i);
            double dj = g.getDegree(iNeighbours[j]);

            centrality = (di * Cb[i%v] + dj * Cb[iNeighbours[j]%v] ) / ((di + dj) * (cnbrs + 1.0));

            if (centrality > maxCentrality)
            {
                maxCentrality = centrality;
                e1 = i;
                e2 = iNeighbours[j];
            }
        }
    }

    cout << e1 << ", " << e2 << ": " << maxCentrality << "\n";
    return make_pair(e1, e2);
}

component getComponents(int left, int right, Graph &G, int n)
{
    queue<int> Q;
    int f;
    Q.push(left);
    int edgeleft = 0, edgeright = 0;
    vector<int> leftComp;
    vector<int> rightComp;
    vector<int> visited(n, 0);
    visited[left] = 1;

    component Comp;

    //G.printGraph();
    int flag = 0;

    //cout<<"h\n"<<left<<" "<<right<<"\nBFS: ";
    while (!Q.empty())
    {
        f = Q.front();
        //cout << f <<" ";
        if (f == right)
        {
            leftComp.resize(0);
            rightComp.resize(0);
            visited.resize(0);
            //cout <<"YES\n";
            Comp.dleft = -1;
            Comp.dright = -1;
            //cout<<"NO\n";
            flag = 1;
            break;
        }
        leftComp.push_back(f);
        Q.pop();
        vector<int> nbrs = G.getNeighbours(f);
        for (int i = 0; i < nbrs.size(); ++i)
        {
            edgeleft++;
            if (visited[nbrs[i]] == 0)
            {
                Q.push(nbrs[i]);
                visited[nbrs[i]] = 1;
            }
        }
    }
    if (flag == 0)
    {
        //cout<<"\n";
        //cout<<"i\n";
        edgeleft /= 2;

        Q.push(right);
        fill(visited.begin(), visited.end(), 0);
        visited[right] = 1;
        while (!Q.empty())
        {
            f = Q.front();
            //cout << f <<"\n";
            if (f == left)
            {
                leftComp.resize(0);
                rightComp.resize(0);
                visited.resize(0);
                Comp.dleft = -1;
                Comp.dright = -1;
                flag = 1;
                break;
            }
            rightComp.push_back(f);
            Q.pop();
            //visited[f]=1;
            vector<int> nbrs = G.getNeighbours(f);
            for (int i = 0; i < nbrs.size(); ++i)
            {
                edgeright++;
                if (visited[nbrs[i]] == 0)
                {
                    Q.push(nbrs[i]);
                    visited[nbrs[i]] = 1;
                }
            }
        }

        if (flag == 0)
        {
            edgeright /= 2;
            //cout<<"j\n";

            int nleft = leftComp.size();
            int nright = rightComp.size();
            double dleft = (nleft > 1) ? edgeleft * 2.0 / (nleft * (nleft - 1)) : 1;
            double dright = (nright > 1) ? edgeright * 2.0 / (nright * (nright - 1)) : 1;

            Comp.dleft = dleft;
            Comp.dright = dright;
            Comp.leftComp = leftComp;
            Comp.rightComp = rightComp;

            leftComp.resize(0);
            rightComp.resize(0);
            visited.resize(0);
        }
    }

    return Comp;
}

vector< vector<int> >getCommunities(Graph &singleSliceGraph, Graph &g, int n, int l)
{

    vector< vector<int> >communities;
    int s = 0;
    int iter = 0;
    while (!g.isEmpty())
    {
        //cout<<"a\n";
        iter++;
        vector<double> Cb = betweennessCentralityNodes(g, n, l);
        cout << "iteration " << iter << "\n"; 
        pair<int, int> topEdge = topBetweennessCentralityInEdges(g, Cb, n, l);

        if (topEdge.first == -1 && topEdge.second == -1)
        {
            //singleSliceGraph.printGraph();
            cout << "g: " << g.numEdges() << "  " << g.numNodes() << "\n";
            break;
        }
        //cout << "NODES: " << topEdge.first << " " << topEdge.second <<"\n";
        int node1 = topEdge.first % n;
        int node2 = topEdge.second % n;



        for (int i = 0; i < l; i++)
        {
            g.removeEdge(make_pair(node1 + i * n, node2 + i * n));
        }

        //cout<<"Before : " << singleSliceGraph.numEdges() << "\n";
        singleSliceGraph.removeEdge(make_pair(node1, node2));
        //cout << "After: " << node1 << " " << singleSliceGraph.getDegree(node1) << " " << node2 << " " << singleSliceGraph.getDegree(node2) << "\n";

        component components = getComponents(node1, node2, singleSliceGraph, n);

        //cout << components.dleft << " " << components.dright << "\n";

        if (components.dleft > 0)
        {
            //cout << components.leftComp.size() << " " << components.rightComp.size() << "\n";

            if (components.dleft > threshold )
            {
                //cout << "making communties left " << components.leftComp.size() << "\n";

                communities.push_back(components.leftComp);
                singleSliceGraph.removeSubgraphFromAllLevels(components.leftComp);
                g.removeSubgraphFromAllLevels(components.leftComp);
            }

            if (components.dright > threshold )
            {
                //cout << "making communties right " << components.rightComp.size() << "\n";
                communities.push_back(components.rightComp);
                singleSliceGraph.removeSubgraphFromAllLevels(components.rightComp);
                g.removeSubgraphFromAllLevels(components.rightComp);
            }
        }
        //cout <<"out\n";
    }

    return communities;
}

int main(int argc, char const *argv[])
{
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
    // singleSliceGraph.printGraph();
    // cout << "\n\n";

    vector< vector<int> > communities = getCommunities(singleSliceGraph, g, numNodes, numSlice);
    cout << "Number of communities" << communities.size() << "\n";
    for (int i = 0; i < communities.size(); i++)
    {
        for (int j = 0; j < communities[i].size(); j++)
        {
            cout << communities[i][j] << " ";
        }
        cout << "\n\n";
    }


    return 0;
}
