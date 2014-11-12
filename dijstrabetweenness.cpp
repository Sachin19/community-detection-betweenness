#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <queue>
#include <map>
#include <cmath>
#include <cstring>
#include <functional>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <cstring>
#include <stack>
#include <set>
#include <fstream>
#include <queue>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/graph_traits.hpp>

typedef boost::property<boost::edge_weight_t, int> EdgeWeightProperty;
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, boost::no_property, EdgeWeightProperty> mygraph;
typedef boost::property_map< mygraph, boost::vertex_index_t>::type VertexIndexMap;

using namespace std;

#define gamma 0.5
#define threshold 0.04

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

int main(int argc, char const *argv[])
{
  	int numNodes,numSlice;
  	mygraph bg;
  	boost::shared_array_property_map<double, boost::property_map<mygraph, boost::vertex_index_t>::const_type> centrality_map(num_vertices(bg), get(boost::vertex_index, bg));
	numNodes = atoi(argv[1]);
	numSlice = atoi(argv[2]);

  char slices[100][100];
  for(int i=0; i<numSlice; i++) {
    strcpy(slices[i],argv[i+3]);
  }

    ifstream fin;
  int x,y;
  vector< pair<int,int> > edgeList;
  Graph g(numNodes, numSlice);
  Graph singleSliceGraph(numNodes,1);
  
  for (int i = 0; i < numSlice; ++i)
  {
    fin.open(slices[i],ios::in);
    
    while(fin) {
      fin>> x >> y;

      if(!g.isEdgePresent(make_pair(x+numNodes*i,y+numNodes*i)) && !g.isEdgePresent(make_pair(y+numNodes*i, x+numNodes*i)))
        g.addEdge(make_pair(x+numNodes*i,y+numNodes*i));
      
      if(!singleSliceGraph.isEdgePresent(make_pair(x,y)) && !singleSliceGraph.isEdgePresent(make_pair(y,x))) {
        singleSliceGraph.addEdge(make_pair(x,y));
        boost::add_edge(x,y,1,bg);
        //gg.add_edge(y, x);
      }
    }

    fin.close();
  }

  //brandes_betweenness_centrality(bg, centrality_map);
  	vector<double> centrality2(numNodes);
    brandes_betweenness_centrality(bg, make_iterator_property_map(centrality2.begin(), get(boost::vertex_index, bg),double()));

    for (int i = 0; i < centrality2.size(); ++i)
    {
    	cout << centrality2[i] << "\n";
    }

  return 0;
}
