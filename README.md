community-detection-betweenness
===============================

1. Community Detection in a multiplex network using betweenness centrality

The input graph should contain edge lists for each layer in seperate files and should be 0 indexed.

Compile the file communityDetection.cpp as "g++ communityDetection.cpp" and run as "./communityDetection num_nodes num_layers layer1_path [layer2_path ... ]

The output of this program contains number of communities  C on the first line, followed by C lines containing space separated integers with ith line denoting the nodes in the ith community

2. Other programs for performing analysis,

a) removeDegreeOneNodes.cpp

	This program removes degree one nodes from the citation network.

	./removeDegreeOneNodes firstLayerPath secondLayerPath

	the output files are name firstLayerPath_small.txt and secondLayerPath_small.txt

b) hashGraphs.cpp
	
	This program takes as input all the layers of the graph and outputs 0 indexed graphs, along with the number of nodes in the graph

	./hashGraphs firstLayerEdgeListPath [secondLayerEdgeListPath ... ]

	the output files are named firstLayerEdgeListPath_ind.txt ...

c) betweennessCentrality.cpp

	Outputs the betweennessCentrality using our implementation (variation of brandes) of each node in the graph. Input same as communityDetection.cpp

d) dijkstrabetweenness.cpp
	
	Outputs the betweennessCentrality using betweennessCentrality of each node in the graph. Input same as communityDetection.cpp. Note that it uses boostlibrary to compile

e) modularity.cpp

	It calculates modularity of communities in a multilayered graph. 

	./modularity <numNodes> <numClusters> <clusterFilePath> <numSlice> <slice1EdgeList> [<slice2EdgeList> ... ]

	clusterFile should contain numNodes lines where ith line contains a single integer indicating the cluster to which ith node belongs 

