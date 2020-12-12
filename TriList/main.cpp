/*****************************************
 *
 * Triangle Listing Algorithms
 * Author: Lijun Chang
 * Date: 05/10/2015
 * Email: ljchang@outlook.com
 *
 ******************************************
 *
 * Three categories of algorithms
 * 1) no-extra space: EdgeIterator (computing common neighbors by sort-merge join)
 * 		i) degree_oriented (optimal complexity)
 * 		ii) full_graph (baseline algorithm)
 * 		iii) least_degree_first (optimal complexity)
 *
 * 2) hash structure (O(m) space): NodeIterator (Implement min {d(u), d(v)})
 * 		i) degree_oriented (optimal complexity)
 * 		ii) full_graph (optimal complexity)
 * 		iii) least_degree_first (optimal complexity)
 *
 * 3) arrayset (O(n) space): EdgeIterator
 * 		i) degree_oriented (optimal complexity)
 * 		ii) full_graph (optimal complexity)
 *		iii) id_oriented (ad hoc algorithm)
 * 		iv) least_degree_first (optimal complexity)
 *
 * 4) compressed input (byte or nibble): EdgeIterator
 * 		i) degree_oriented (optimal complexity)
 * 		ii) full_graph (optimal complexity)
 * 		iii) id_oriented (ad hoc algorithm)
 *
 *****************************************/

#include "Timer.h"
#include "Graph.h"

void print_usage() {
	printf("Usage: [1]exe [2]algorithm [3]graph-dir\n");
	printf("Algorithms:\n");
	printf("\t 1. EdIt-OC\n");
	printf("\t 2. EdIt-OV\n");
	printf("\t 3. EdIt-rOV\n");
	printf("\t 4. NoIt-OE\n");
	printf("\t 5. EdIt-FC\n");
	printf("\t 6. EdIt-FCo\n");
	printf("\t 7. EdIt-FV\n");
	printf("\t 8. EdIt-FVo\n");
	printf("\t 9. EdIt-FE\n");
}

int main(int argc, char *argv[]) {
	if(argc < 3) {
		print_usage();
		return 0;
	}

#ifndef NDEBUG
	printf("**** Triangle Listing (Debug): %s %s *** ", argv[1], argv[2]);
#else
	printf("**** Triangle Listing (Release): %s %s *** ", argv[1], argv[2]);
#endif
	printf("\n");

	Timer t;

	Graph *graph = new Graph(argv[2]);
	graph->read_graph();
	//printf("\t*** Finished reading graph\n");

	long long mtime1 = t.elapsed();
	t.restart();

	//if(strcmp(argv[1], "NoIt-OE") == 0||strcmp(argv[1], "ni-least") == 0) graph->ni_oriented();
	if(strcmp(argv[1], "EdIt-FE") == 0) graph->ei_full_e();
	else if(strcmp(argv[1], "EdIt-FV") == 0) graph->ei_full_v();
	else if(strcmp(argv[1], "EdIt-FVo") == 0) graph->ei_full_vo();
	else if(strcmp(argv[1], "EdIt-FVo2") == 0) graph->ei_full_vo2();
	else if(strcmp(argv[1], "EdIt-FC") == 0) graph->ei_full_c();
	//else if(strcmp(argv[1], "EdIt-FCo") == 0) graph->ei_full_co();
	//else if(strcmp(argv[1], "EdIt-rOV") == 0) graph->ei_roriented_v();
	//else if(strcmp(argv[1], "EdIt-OV") == 0||strcmp(argv[1], "ei-least-v") == 0||strcmp(argv[1], "ei-most-v") == 0) graph->ei_oriented_v();
	//else if(strcmp(argv[1], "EdIt-OC") == 0||strcmp(argv[1], "ei-least-c") == 0) graph->ei_oriented_c();
	else print_usage();

	// printf("\t*** Finished triangulation\n");

	return 0;
}
