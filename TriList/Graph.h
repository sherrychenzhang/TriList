#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"
#include "Timer.h"
#include "HashSet.h"

using namespace std;

class Graph {
private:
	string dir; //input graph directory
	ui n, m; //number of nodes and edges of the graph

	ui *pstart; //offset of neighbors of nodes
	ui *edges; //adjacent ids of edges

#ifdef _HASHSET_
	HashSet **hs;
	//unordered_set<int> **hs;
#else
	dense_hash_set<int> **hs;
#endif

public:
	Graph(const char *_dir) ;
	~Graph() ;

	void read_graph() ;

	void ni_oriented() ;
	void ei_full_e() ;
	void ei_full_v() ;
	void ei_full_vo() ;
	void ei_full_vo2() ;
	void ei_full_c() ;
	void ei_full_co() ;
	void ei_oriented_v() ;
	void ei_roriented_v() ;
	void ei_oriented_c() ;

private:
	inline ui binary_search(const ui *array, ui e, ui val) ;
	inline void cross_link(ui *reverse) ;
	void build_degree_oriented_graph() ;
	void build_degree_oriented_graph_reverse() ;
	void build_least_oriented_graph() ;

	void build_hashset() ;
	int find_hashset(ui a, ui b) ;
};

#endif
