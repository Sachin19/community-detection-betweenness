//
// Betweenness centrality of undirected unweighted graph (Brandes)
//
// Description:
// 
//   Compute betweenness centrality, defined by
//     f(u) := \sum_{u,t \eq v} |s-t shortest paths that contains v|/|s-t shortest paths|
//
// Algorithm:
//
//   Brandes's algorithm, O(nm) time, O(m) space.
//
// References:
//
//   U. Brandes (2001): A faster algorithm for betweenness centrality.
//   Journal of Mathematical Sociology, vol.25, pp.163–177.

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

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct edge {
  size_t src, dst;
};
struct graph {
  vector<edge> edges;
  void add_edge(size_t src, size_t dst) {
    edges.push_back({src, dst});
  }
  size_t n;
  vector<vector<edge>> adj;
  void make_graph(int n_ = 0) {
    n = n_;
    for (auto e: edges) 
      n = max(n, max(e.src, e.dst)+1);
    adj.resize(n);
    for (auto e: edges) {
      adj[e.src].push_back(e);
      swap(e.src, e.dst);
      adj[e.src].push_back(e);
    }
  }

  vector<double> betweeness_centrality() {
    vector<double> centrality(n);

    for (size_t s = 0; s < n; ++s) {
      vector<size_t> S;
      vector<double> sigma(n); sigma[s] = 1;
      vector<int> dist(n, -1); dist[s]  = 0;
      queue<size_t> que;       que.push(s);
      while (!que.empty()) {
        size_t u = que.front();
        S.push_back(u);
        que.pop();
        for (auto e: adj[u]) {
          if (dist[e.dst] < 0) {
            dist[e.dst] = dist[e.src] + 1;
            que.push(e.dst);
          }
          if (dist[e.dst] == dist[e.src] + 1) {
            sigma[e.dst] += sigma[e.src];
          }
        }
      }
      vector<double> delta(n);
      while (!S.empty()) {
        size_t u = S.back();
        S.pop_back();
        for (auto e: adj[u]) {
          if (dist[e.dst] == dist[e.src] + 1) {
            delta[e.src] += sigma[e.src] / sigma[e.dst] * (1 + delta[e.dst]);
          }
        }
        if (u != s) centrality[u] += delta[u];
      }
    }
    return centrality;
  }
};



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

//Graph singleSliceGraph;

Graph getGraph(int numSlice, char slices[][100], int numNodes) {
  // ifstream fin;
  // int x,y;
  // vector< pair<int,int> > edgeList;
  // Graph g(numNodes, numSlice);
  // singleSliceGraph = Graph(numNodes,1);
  
  // for (int i = 0; i < numSlice; ++i)
  // {
  //  fin.open(slices[i],ios::in);
    
  //  while(fin) {
  //    fin>> x >> y;

  //    g.addEdge(make_pair(x+numNodes*i,y+numNodes*i));
      
  //    if(!singleSliceGraph.isEdgePresent(make_pair(x,y))) {
  //      singleSliceGraph.addEdge(make_pair(x,y));
  //    }
  //  }

  //  fin.close();
  // }

  // for (int i = 0; i < numSlice; ++i)
  // {
  //  for (int j = i+1; j < numSlice; ++j)
  //  {
  //    for (int k = 0; k < numNodes; ++k)
  //    {
  //      g.addEdge(make_pair(k+numNodes*i,k+j*numNodes));
  //    }
  //  }
  // }


  //return g;
}

vector<double> betweennessCentralityNodes(Graph &g, int n, int l) {


  vector<double>Cb(n,0.0);

  for (int i = 0; i < n; i++)
  {
    if(g.isDeleted(i)){
      continue;
    }
    //cout<<"dsf\n";
    stack<int>S;
    vector< stack<int> > P(n*l);
    vector<double> sigma(n*l,0.0);
    vector<double> sigmaM(n,0.0);
    vector<int> D(n*l,-1);
    vector<int> Dm(n*l,-1);
    vector< vector<int> > Vsame(n);
    queue<int> Q;
    set<int> W;

    fill(P.begin(), P.end(), stack<int>());

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

    fill(Vsame.begin(),Vsame.end(),vector<int>());
    Q.push(i);

    while(!Q.empty()) {
      int v = Q.front();
      Q.pop();
      if(g.isDeleted(v))
        continue;
      

      S.push(v);

      W.clear();

      if(v!=i) {
        vector<int> vNeighbours = g.getNeighbours(v);
        for(int j=0;j<vNeighbours.size();j++) {
          W.insert(vNeighbours[j]);
        }
      }
      else {
        for (int j = i; j < n*l; j += n)
        {
          vector<int> jNeighbours = g.getNeighbours(j);
          for(int k=0;k<jNeighbours.size();k++) {
            W.insert(jNeighbours[k]);
          }
        }
      }

      set<int>::iterator it;

      for(it = W.begin(); it!=W.end(); it++) {
        int w = *it;
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

    for(int w =0; w<n; w++) {
      sigmaM[w] = 0.0;
      for(int j=0; j<Vsame[w].size(); j++) {
        int v = Vsame[w][j];
        sigmaM[w] += sigma[v];
      }
    }

    vector<double> delta(n*l,0.0);    

    while(!S.empty()) {
      int w = S.top();
      S.pop();

      while(!P[w].empty()) {
        int v = P[w].top();
        P[w].pop();

        vector<int>::iterator p;
        p = find(Vsame[w%n].begin(),Vsame[w%n].end(),w);

        if(p != Vsame[w%n].end()) {
          //if(sigma[w]==0 || sigmaM[w]==0) {
          //  exit(1);
          //}
          if(sigma[w] > 0 && sigmaM[w%n]>0)
            delta[v] += (1-gamma)*(sigma[v]/sigma[w])*((sigma[w]/sigmaM[w%n]) + delta[w]);
          else
            cout <<"load";
          //cout<<delta[v]<<" ";
        }
        else {
          //if(sigma[w]==0) {
          //  exit(1);
          //}
          if(sigma[w] > 0 && sigmaM[w%n]>0)
            delta[v] += gamma*(sigma[v]/sigma[w])*delta[w];
          else
            cout<<"sdfd";
          //cout<<delta[v]<<" ";
        }
        if(w!=i) {
          Cb[w%n] += delta[w];
        }
      }
    }
  }

  return Cb;
}

pair<int,int> topBetweennessCentralityInEdges(Graph &g, vector<double> Cb, int v, int l) {
  int n = v*l;

  double maxCentrality = 0.0;
  double centrality;
  int e1=-1,e2=-1;

  for(int i=0; i<n; i++) {
    //cout<<"ta\n";
    if(g.isDeleted(i))
      continue;
    vector<int> iNeighbours = g.getNeighbours(i);
    //cout << "size: "<<iNeighbours.size()<<"\n";
    for(int j=0; j<iNeighbours.size(); j++) {
      if(g.isDeleted(iNeighbours[j]))
        continue;
      vector<int> jNeighbours = g.getNeighbours(iNeighbours[j]);

      vector<int> commonNeighbours(n);
      sort(iNeighbours.begin(),iNeighbours.end());
      sort(jNeighbours.begin(),jNeighbours.end());
      //cout<<"tb\n";
      std::vector<int>::iterator it;
      it=std::set_intersection (iNeighbours.begin(), iNeighbours.end(), jNeighbours.begin(), 
        jNeighbours.end(), commonNeighbours.begin());

      commonNeighbours.resize(it-commonNeighbours.begin());
      //cout<<"tc\n";
      double di = g.getDegree(i);
      // for(int k = i; k<n; k+=v)
      //  di += g.getDegree(k);

      double dj = g.getDegree(iNeighbours[j]);
      // for(int k = j; k<n; k+=v)
      //  dj += g.getDegree(k);
      // cout<<"sdf"<<(di+dj)*(commonNeighbours.size()+1)<<"\n";
      // cout<<"sdfd"<<i%n<<" "<<iNeighbours[j]<<"\n";
      centrality = (di*Cb[i] + dj*Cb[iNeighbours[j]] )/((di+dj)*(commonNeighbours.size()+1.0));
      //cout<<"c: "<<centrality<<"\n";
      if(centrality > maxCentrality) {
        maxCentrality = centrality;
        e1 = i;
        e2 = iNeighbours[j];
      }
    }
  }

  cout << e1 << ", " << e2 << ": " << maxCentrality<<"\n";
  return make_pair(e1,e2);
}

component getComponents(int left, int right, Graph &G, int n){
  queue<int> Q;
  int f;
  Q.push(left);
  int edgeleft=0,edgeright=0;
  vector<int> leftComp;
  vector<int> rightComp;
  vector<int> visited(n,0);
  visited[left]=1;

  component Comp;

  //G.printGraph();
  int flag=0; 

  //cout<<"h\n"<<left<<" "<<right<<"\nBFS: ";
  while(!Q.empty()) {
    f = Q.front();
    //cout << f <<" ";
    if(f==right){
      leftComp.resize(0);
      rightComp.resize(0);
      visited.resize(0);
      //cout <<"YES\n";
      Comp.dleft = -1;
      Comp.dright = -1;
      //cout<<"NO\n";
      flag=1;
      break;
    }
    leftComp.push_back(f);
    Q.pop();
    vector<int> nbrs = G.getNeighbours(f);
    for (int i = 0; i < nbrs.size(); ++i)
    {
      edgeleft++;
      if(visited[nbrs[i]]==0){
        Q.push(nbrs[i]);
        visited[nbrs[i]]=1;
      }
    }
  }
  if(flag==0){
    //cout<<"\n";
    //cout<<"i\n";
    edgeleft/=2;

    Q.push(right);
    fill(visited.begin(),visited.end(),0);
    visited[right]=1;
    while(!Q.empty()){
      f = Q.front();
      //cout << f <<"\n";
      if(f==left){
        leftComp.resize(0);
        rightComp.resize(0);
        visited.resize(0);
        Comp.dleft = -1;
        Comp.dright = -1;
        flag=1;
        break;
      }
      rightComp.push_back(f);
      Q.pop();
      //visited[f]=1;
      vector<int> nbrs = G.getNeighbours(f);
      for (int i = 0; i < nbrs.size(); ++i)
      {
        edgeright++;
        if(visited[nbrs[i]]==0){
          Q.push(nbrs[i]);
          visited[nbrs[i]]=1;
        }
      }
    }

    if(flag==0){
      edgeright/=2;
      //cout<<"j\n";

      int nleft = leftComp.size();
      int nright = rightComp.size();
      double dleft = (nleft > 1) ? edgeleft*2.0/(nleft*(nleft-1)) : 1;
      double dright = (nright > 1) ? edgeright*2.0/(nright*(nright-1)) : 1;

      Comp.dleft = dleft;
      Comp.dright = dright;
      Comp.leftComp = leftComp;
      Comp.rightComp = rightComp;

      //cout<<"k\n";
      leftComp.resize(0);
      rightComp.resize(0);
      visited.resize(0);
    }
  }

  return Comp;
  
  //return result;
}

vector< vector<int> >getCommunities(Graph &singleSliceGraph, Graph &g, int n, int l) {

  vector< vector<int> >communities;
  int s=0;
  
  while(!g.isEmpty()) {
    //cout<<"a\n";
    vector<double> Cb = betweennessCentralityNodes(g,n,l);
    pair<int,int> topEdge = topBetweennessCentralityInEdges(g,Cb,n,l);

    if(topEdge.first == -1 && topEdge.second == -1){
      //singleSliceGraph.printGraph();
      cout <<"g: "<<g.numEdges()<<"  "<<g.numNodes()<<"\n";
      break;
    }
    //cout << "NODES: " << topEdge.first << " " << topEdge.second <<"\n";
    int node1 = topEdge.first % n;
    int node2 = topEdge.second % n;

    

    for(int i=0;i<l;i++) {
      g.removeEdge(make_pair(node1+i*n,node2+i*n));
    }

    //cout<<"Before : " << singleSliceGraph.numEdges() << "\n";
    singleSliceGraph.removeEdge(make_pair(node1,node2));
    cout<<"After: " << node1 << " "<< singleSliceGraph.getDegree(node1)<<" "<< node2 << " "<< singleSliceGraph.getDegree(node2) << "\n";

    component components = getComponents(node1,node2,singleSliceGraph, n);
    
    //cout << components.dleft << " " << components.dright << "\n";
    
    if(components.dleft > 0) {
      //cout << components.leftComp.size() << " " << components.rightComp.size() << "\n";
      
      if(components.dleft > threshold ) {
        cout << "making communties left " << components.leftComp.size() << "\n";

        communities.push_back(components.leftComp);
        singleSliceGraph.removeSubgraphFromAllLevels(components.leftComp);
        g.removeSubgraphFromAllLevels(components.leftComp);
      }

      if(components.dright > threshold ) {
        cout << "making communties right " << components.rightComp.size() << "\n";
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
  int numNodes,numSlice;
  graph gg;
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
        gg.add_edge(x, y);
        //gg.add_edge(y, x);
      }
    }

    fin.close();
  } 
  gg.make_graph();
  vector<double> Cb = gg.betweeness_centrality();

    for(int i=0;i<Cb.size();i++)
      if(Cb[i] > 0)
        cout<<i << ": " << Cb[i]<<"\n";


  return 0;
}
