/******************************************************************************
 * graph_io.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <math.h>
#include <sstream>
#include "graph_io.h"

#define MIN(A,B) ((A)<(B))?(A):(B)


graph_io::graph_io() {

}

graph_io::~graph_io() {

}

int graph_io::writeGraph_HMetisFormat(graph_access & G, const std::string & filename) {
        std::ofstream f(filename.c_str());
//        f << G.number_of_edges()/2 << " " << G.number_of_nodes() << " 11" <<  std::endl;
        f << G.number_of_edges()/2 << " " << G.number_of_nodes() <<  std::endl;

        forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
			if (node > G.getEdgeTarget(e)) {
//				f << G.getEdgeWeight(e)  <<  " " << (node+1) << " " << (G.getEdgeTarget(e)+1) <<  "\n";
				f << (node+1) << " " << (G.getEdgeTarget(e)+1) <<  "\n";
			}
                } endfor
        } endfor

//        forall_nodes(G, node) {
//                f <<  G.getNodeWeight(node)  <<  "\n";
//        } endfor

        f.close();
        return 0;
}

int graph_io::writeGraphWeighted(graph_access & G, const std::string & filename) {
        std::ofstream f(filename.c_str());
        f << G.number_of_nodes() <<  " " <<  G.number_of_edges()/2 <<  " 11" <<  std::endl;

        forall_nodes(G, node) {
                f <<  G.getNodeWeight(node) ;
                forall_out_edges(G, e, node) {
                        f << " " <<   (G.getEdgeTarget(e)+1) <<  " " <<  G.getEdgeWeight(e) ;
                } endfor
                f <<  "\n";
        } endfor

        f.close();
        return 0;
}

int graph_io::writeGraph(graph_access & G, const std::string & filename) {
        std::ofstream f(filename.c_str());
        f << G.number_of_nodes() <<  " " <<  G.number_of_edges()/2 << std::endl;

        forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        f <<   (G.getEdgeTarget(e)+1) << " " ;
                } endfor
                f <<  "\n";
        } endfor

        f.close();
        return 0;
}

int graph_io::readEdgeStreamToGraph(graph_access & G, const std::string & filename, bool relabel_nodes) {
        std::string line;
        long nmbNodes;
        long nmbEdges;
	int i, x, y;
	NodeID target;
	int num_edges = 0;
	int num_nodes = 0;
	int real_x, real_y;

        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in) {
                std::cerr << "Error opening " << filename << std::endl;
                return 1;
        }

        std::getline(in,line);
        //skip comments
        while( line[0] == '%' || line[0] == '#' ) {
                std::getline(in, line);
        }

        std::stringstream ss(line);
        ss >> nmbNodes;
        ss >> nmbEdges;

	std::unordered_map<int,int> map_virtual_real;

	std::vector<std::set<int>> neighbors;
	neighbors.resize(nmbNodes);
        while(  std::getline(in, line)) {
                if (line[0] == '%' || line[0] == '#' ) {
                        continue;
                }
		std::stringstream ss(line);
		ss >> x;
		ss >> y;

		if (x == y) {
			continue;
		}

		if (relabel_nodes) {
			if (map_virtual_real.find(x) == map_virtual_real.end()) {
				map_virtual_real[x] = num_nodes++;
			}
			if (map_virtual_real.find(y) == map_virtual_real.end()) {
				map_virtual_real[y] = num_nodes++;
			}
			real_x = map_virtual_real[x];
			real_y = map_virtual_real[y];
		} else {
			real_x = x;
			real_y = y;
			if (real_x >= nmbNodes || real_y >= nmbNodes) {
				std::cerr <<  "Nodes need to be relabeled for this input."  << std::endl;
				exit(0);
			}
		}	

//std::cout << x << " " << y << "\n";
		if (neighbors[real_x].find(real_y) == neighbors[real_x].end()) {
			neighbors[real_x].insert(real_y);
			num_edges++;
		}
		if (neighbors[real_y].find(real_x) == neighbors[real_y].end()) {
			neighbors[real_y].insert(real_x);
			num_edges++;
		}
	} 

	nmbEdges = num_edges;
        if( nmbEdges > std::numeric_limits<int>::max() || nmbNodes > std::numeric_limits<int>::max() ) {
                std::cerr <<  "The graph is too large. Currently only 32bit supported!"  << std::endl;
                exit(0);
        }
	
	G.start_construction(nmbNodes, nmbEdges);
	for( i = 0 ; i < nmbNodes ; i++) {
                NodeID node = G.new_node(); 
                G.setPartitionIndex(node, 0);
                NodeWeight weight = 1;
                G.setNodeWeight(node, weight);
		for (std::set<int>::iterator it=neighbors[node].begin(); it!=neighbors[node].end(); ++it) {
			target = *it;
                        EdgeWeight edge_weight = 1;
                        EdgeID e = G.new_edge(node, target);
                        G.setEdgeWeight(e, edge_weight);
                }
        }
	return 0;
}

int graph_io::readMatrixToGraph(PartitionConfig & config, graph_access & G, const std::string & filename) {
        std::string line;
        long nmbNodes;
        long nmbEdges;

        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in) {
                std::cerr << "Error opening " << filename << std::endl;
                return 1;
        }

        std::getline(in,line);
        //skip comments
        while( line[0] == '%' ) {
                std::getline(in, line);
        }

	int i, x, y;
	int *rowindices, *colindices;
	NodeID target;

        std::stringstream ss(line);
	// #rows #columns #nonzeros
        ss >> config.matrix_m;
        ss >> config.matrix_n;
        ss >> config.matrix_nnz;

        if( config.matrix_m+config.matrix_n+config.matrix_nnz > std::numeric_limits<int>::max() ) {
                std::cerr <<  "The graph is too large. Currently only 32bit supported!"  << std::endl;
                exit(0);
        }

	std::vector<std::vector<int>> neighbors;
	neighbors.resize(config.matrix_nnz + config.matrix_m + config.matrix_n);
	rowindices = (int *) malloc(sizeof(int) * (config.matrix_nnz+1));
	colindices = (int *) malloc(sizeof(int) * (config.matrix_nnz+1));

	for( i = 0 ; i < config.matrix_nnz ; i++) {
		std::getline(in, line);
                if (line[0] == '%') { // a comment in the file
                        continue;
                }

		std::stringstream ss(line);
		// rowIndex columnIndex value
		ss >> rowindices[i];
		ss >> colindices[i];
		// we just ignore the value of the matrix
		rowindices[i] --;
		colindices[i] --;

		x = config.matrix_nnz + colindices[i];
		neighbors[i].push_back(x);
		neighbors[x].push_back(i);
		
		y = config.matrix_nnz + config.matrix_n + rowindices[i];
		neighbors[i].push_back(y);
		neighbors[y].push_back(i);
	} 

	nmbNodes = config.matrix_nnz + config.matrix_m + config.matrix_n;
	nmbEdges = 4*config.matrix_nnz;
	G.start_construction(nmbNodes, nmbEdges);
	for( i = 0 ; i < nmbNodes ; i++) {
                NodeID node = G.new_node(); 
                G.setPartitionIndex(node, 0);
                NodeWeight weight = 1;
                G.setNodeWeight(node, weight);
		for(int j : neighbors[node]) {
			target = j;
                        EdgeWeight edge_weight = 1;
                        EdgeID e = G.new_edge(node, target);
                        G.setEdgeWeight(e, edge_weight);
                }
        }
	
	free(rowindices);
	free(colindices);
	return;
}


int graph_io::readPartition(graph_access & G, const std::string & filename) {
        std::string line;

        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in) {
                std::cerr << "Error opening file" << filename << std::endl;
                return 1;
        }

        PartitionID max = 0;
        forall_nodes(G, node) {
                // fetch current line
                std::getline(in, line);
                if (line[0] == '%') { //Comment
                        node--;
                        continue;
                }

                // in this line we find the block of Node node
                G.setPartitionIndex(node, (PartitionID) atol(line.c_str()));

                if(G.getPartitionIndex(node) > max)
                        max = G.getPartitionIndex(node);
        } endfor

        G.set_partition_count(max+1);
        in.close();

        return 0;
}


int graph_io::readGraphWeighted(graph_access & G, const std::string & filename) {
        std::string line;

        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in) {
                std::cerr << "Error opening " << filename << std::endl;
                return 1;
        }

        long nmbNodes;
        long nmbEdges;

        std::getline(in,line);
        //skip comments
        while( line[0] == '%' ) {
                std::getline(in, line);
        }

        int ew = 0;
        std::stringstream ss(line);
        ss >> nmbNodes;
        ss >> nmbEdges;
        ss >> ew;

        if( 2*nmbEdges > std::numeric_limits<int>::max() || nmbNodes > std::numeric_limits<int>::max()) {
                std::cerr <<  "The graph is too large. Currently only 32bit supported!"  << std::endl;
                exit(0);
        }

        bool read_ew = false;
        bool read_nw = false;

        if(ew == 1) {
                read_ew = true;
        } else if (ew == 11) {
                read_ew = true;
                read_nw = true;
        } else if (ew == 10) {
                read_nw = true;
        }
        nmbEdges *= 2; //since we have forward and backward edges

        NodeID node_counter   = 0;
        EdgeID edge_counter   = 0;
        long long total_nodeweight = 0;

        G.start_construction(nmbNodes, nmbEdges);

        while(  std::getline(in, line)) {

                if (line[0] == '%') { // a comment in the file
                        continue;
                }

                NodeID node = G.new_node(); node_counter++;
                G.setPartitionIndex(node, 0);

                std::stringstream ss(line);

                NodeWeight weight = 1;
                if( read_nw ) {
                        ss >> weight;
                        total_nodeweight += weight;
                        if( total_nodeweight > (long long) std::numeric_limits<NodeWeight>::max()) {
                                std::cerr <<  "The sum of the node weights is too large (it exceeds the node weight type)."  << std::endl;
                                std::cerr <<  "Currently not supported. Please scale your node weights."  << std::endl;
                                exit(0);
                        }
                }
                G.setNodeWeight(node, weight);

                NodeID target;
                while( ss >> target ) {
                        //check for self-loops
                        if(target-1 == node) {
                                std::cerr <<  "The graph file contains self-loops. This is not supported. Please remove them from the file."  << std::endl;
                        }

                        EdgeWeight edge_weight = 1;
                        if( read_ew ) {
                                ss >> edge_weight;
                        }
                        edge_counter++;
                        EdgeID e = G.new_edge(node, target-1);

                        G.setEdgeWeight(e, edge_weight);
                }

                if(in.eof()) {
                        break;
                }
        }

        if( edge_counter != (EdgeID) nmbEdges ) {
                std::cerr <<  "number of specified edges mismatch"  << std::endl;
                std::cerr <<  edge_counter <<  " " <<  nmbEdges  << std::endl;
                exit(0);
        }

        if( node_counter != (NodeID) nmbNodes) {
                std::cerr <<  "number of specified nodes mismatch"  << std::endl;
                std::cerr <<  node_counter <<  " " <<  nmbNodes  << std::endl;
                exit(0);
        }


        G.finish_construction();
        return 0;
}


void graph_io::writePartition(graph_access & G, const std::string & filename) {
        std::ofstream f(filename.c_str());
        std::cout << "writing partition to " << filename << " ... " << std::endl;

        forall_nodes(G, node) {
                f << G.getPartitionIndex(node) <<  "\n";
        } endfor

        f.close();
}


void graph_io::writeSpMxVPartition(PartitionConfig & config, graph_access & G, const std::string & filename) {
	int node=0;

        std::ofstream f1(filename+"NZ");
        std::cout << "writing partition to " << filename << "NZ" << " ... " << std::endl;
	while (node < config.matrix_nnz) {
                f1 << G.getPartitionIndex(node++) <<  "\n";
        } 
        f1.close();

        std::ofstream f2(filename+"IN");
        std::cout << "writing partition to " << filename << "IN" << " ... " << std::endl;
	while (node < config.matrix_nnz + config.matrix_n) {
                f2 << G.getPartitionIndex(node++) <<  "\n";
        }
        f2.close();

        std::ofstream f3(filename+"OUT");
        std::cout << "writing partition to " << filename << "OUT" << " ... " << std::endl;
	while (node < config.matrix_nnz + config.matrix_n + config.matrix_m) {
                f3 << G.getPartitionIndex(node++) <<  "\n";
        }
        f3.close();

}


