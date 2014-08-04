#include <iostream>
#include <fstream>
#include <algorithm>
#include <atomic>
#include "graph.h"
using namespace std;

int optimal8node[8][8] =

{ { 0, 1, 0, 0, 1, 0, 0, 1 },
{ 1, 0, 1, 0, 0, 0, 1, 0 },
{ 0, 1, 0, 1, 0, 1, 0, 0 },
{ 0, 0, 1, 0, 1, 0, 0, 1 },
{ 1, 0, 0, 1, 0, 1, 0, 0 },
{ 0, 0, 1, 0, 1, 0, 1, 0 },
{ 0, 1, 0, 0, 0, 1, 0, 1 },
{ 1, 0, 0, 1, 0, 0, 1, 0 } };

int optimal16node[16][16] = { { 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1 },
							  { 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 },
							  { 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0 },
							  { 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 },
							  { 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0 },
							  { 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0 },
							  { 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0 },
							  { 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1 },
							  { 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0 },
							  { 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0 },
							  { 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0 },
							  { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1 },
							  { 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0 },
							  { 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0 },
							  { 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1 },
							  { 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0 } };



int * load_file(char * filename, int nodes){
	ifstream file(filename);
	int * x = (int*)malloc(nodes * nodes * sizeof(int));
	int k = 0;
	char c;
	while (!file.eof() && k < nodes * nodes){
		file >> noskipws >> c;
		if (c == '0') x[k++] = 0;
		else if (c == '1') x[k++] = 1;
	}
	return x;
}

bool * hypercube(int size){
	int nodes = 1 << size;
	bool * x = (bool*)malloc(nodes * nodes);
	memset(x, false, nodes * nodes * sizeof(bool));
	if (size == 0) {
		x[0] = false;
		return x;
	}
	bool * prev = hypercube(size - 1);
	for (int i = 0; i < nodes / 2; i++){
		for (int j = 0; j < nodes / 2; j++){
			x[i * nodes + j] = prev[i * (nodes / 2) + j];
			int a = i + nodes / 2, b = j + nodes / 2;
			x[a * nodes + b] = prev[i * (nodes / 2) + j];
		}
	}
	for (int i = 0; i < nodes / 2; i++){
		int j = i + nodes / 2;
		x[i * nodes + j] = true;
		x[j * nodes + i] = true;
	}
	return x;
}

float * matrix_mult(float ** matrix, float * v, int size){
	float * result = new float[size];
	for (int i = 0; i < size; i++){
		float sum = 0;
		for (int j = 0; j < size; j++){
			sum += matrix[i][j] * v[j];
		}
		result[i] = sum;
	}
	return result;
}

void dist_analysis(bool * x, int nodes){
	graph g(x, nodes);
	cout << "Node\tAverage\tDiameter\n";
	for (int i = 0; i < nodes; i++){
		cout << i << "\t" << g.node_average(i) << "\t" << g.node_diameter(i) << "\n";
	}
	cout << "Diameter: " << g.diameter() << "\n";
	cout << "Average:  " << g.average() << "\n";
}

void connecting_optimal_8_nodes(){
	bool x[16][16];
	for (int i = 0; i < 8; i++){
		for (int j = 0; j < 8; j++){
			if (optimal8node[i][j]){
				x[i][j] = true;
				x[i + 8][j + 8] = true;
			}
			else {
				x[i][j] = false;
				x[i + 8][j + 8] = false;

			}
		}
	}
	int pos[8];
	for (int i = 0; i < 8; i++){
		pos[i] = i;
	}
	do {
		for (int i = 0; i < 8; i++){
			for (int j = 0; j < 8; j++){
				x[8 + i][j] = false;
				x[i][8 + j] = false;
			}
		}
		for (int i = 0; i < 8; i++){
			x[8 + i][pos[i]] = true;
			x[i][pos[i] + 8] = true;
		}
		graph g(&x[0][0], 16);
		if (g.diameter() != 3) cout << "Incorrect Diameter!\n";
		else if (g.average() <= 1.8){
			cout << "Found graph with average " << g.average() << "\n";
			for (int i = 0; i < 8; i++) cout << pos[i];
		}
	} while (std::next_permutation(pos, pos + 8));
	cout << "Done!\n";
	getchar();
}

int max_degree(int * g, int nodes){
	int x = 0;
	for (int i = 0; i < nodes; i++){
		int sum = 0;
		for (int j = 0; j < nodes; j++){
			if (g[nodes * i + j]) sum++;
		}
		if (sum > x) x = sum;
	}
	return x;
}

bool is_regular(int * g, int nodes, int k){
	for (int i = 0; i < nodes; i++){
		int sum = 0;
		for (int j = 0; j < nodes; j++){
			if (g[nodes * i + j]) sum++;
		}
		if (sum != k) return false;
	}
	return true;
}

bool check_degree_regular(int * deg, int n, int k){
	while (--n > 0 && deg[n] == k);
	return n == 0;
}

std::mutex changing_diameter;

void connect_graphs(int * source, int source_size, int nodes, int k, char * file_output){
	int * x = combineGraphs(source, source_size, nodes);
	//For the fist 16x4 64 graphs
	//x[19 * nodes +  4] = 1; x[ 4 * nodes + 19] = 1;
	//x[40 * nodes + 18] = 1; x[18 * nodes + 40] = 1;
	//x[60 * nodes + 40] = 1; x[40 * nodes + 60] = 1;

	srand(0);
	for (int i = 0; i < source_size - 1; i++){
		int j = i + 1;
		int u = rand() % source_size;
		int v = rand() % source_size;
		x[i * nodes * source_size + j * source_size + u * nodes + v] = 1;
		x[j * nodes * source_size + i * source_size + v * nodes + u] = 1;
	}
	int v1, v2;
	int best_diameter = 100000000;
	float best_average = 100000000.0f;
	int count = 0;
	int * deg = new int[nodes];
	for (int i = 0; i < nodes; i++){
		int asdf = 0;
		for (int j = 0; j < nodes; j++){
			asdf += x[i * nodes + j];
		}
		deg[i] = asdf;
	}
	vector<int> searching_nodes, searchable_nodes;
	for (int i = 0; i < nodes; i++) searchable_nodes.push_back(i);
	while (!check_degree_regular(deg, nodes, k)){
		v1 = -1; v2 = -1;
		searching_nodes.clear();
		for (int i = 0; i < EDGE_SEARCH_SIZE && i < searchable_nodes.size();){
			int j = searchable_nodes[rand() % searchable_nodes.size()];
			bool asdf = true;
			for (auto r = searching_nodes.begin(); r != searching_nodes.end(); r++)
				if (*r == j){
					asdf = false;
					break;
				}
			if (asdf){
				searching_nodes.push_back(j);
				i++;
			}
		}
		int k = 0;
		vector<std::thread> threads;
		for (vector<int>::iterator a = searching_nodes.begin(); a != searching_nodes.end(); a++){
			if (k < MAX_THREADS){
				threads.push_back(std::thread(compare_graph, x, nodes, source_size,
									  		  *a, searching_nodes,
											  &best_diameter, &best_average,
											  &v1, &v2));
				k++;
			}
			else {
				for (vector<std::thread>::iterator t = threads.begin(); t != threads.end(); t++){
					t->join();
				}
				threads.clear();
				k = 0;
			}
		}
		if (v1 < 0){
			//find the two nodes of irregular degree and connect them
			for (int i = 0; i < nodes; i++){
				if (deg[i] < k){
					v1 = i;
					break;
				}
			}
			for (int i = v1 + 1; i < nodes; i++){
				if (deg[i] < k){
					v2 = i;
					break;
				}
			}
		}
		x[v1 * nodes + v2] = 1;
		x[v2 * nodes + v1] = 1;
		deg[v1]++;
		deg[v2]++;
		if (deg[v1] == k){
			for (vector<int>::iterator i = searchable_nodes.begin(); i != searchable_nodes.end(); i++){
				if (*i == v1){
					searchable_nodes.erase(i);
					break;
				}
			}
		}
		if (deg[v2] == k){
			for (vector<int>::iterator i = searchable_nodes.begin(); i != searchable_nodes.end(); i++){
				if (*i == v2){
					searchable_nodes.erase(i);
					break;
				}
			}
		}
		count++;
		std::cout << "Added an edge\n";
	}
	/*ofstream f(file_output);
	for (int i = 0; i < nodes; i++){
		for (int j = 0; j < nodes; j++){
			f << x[nodes * i + j];
		}
		f << "\n";
	}*/
	graph gr(x, nodes);
	//f << "Diameter: " << gr.diameter() << "\nAverage: " << gr.average() << "\n";
//	f.close();
	std::cout << source_size << "x" << source_size << " graph coonected to form a " << nodes << "x" << nodes << "graph\n";
	std::cout << "Diameter: " << gr.diameter() << "\nAverage: " << gr.average() << "\n";
	delete[] deg;
	delete[] x;
}

void compare_graph(int * gr, int nodes, int source_size,
	int i, vector<int> search_nodes,
	int * best_diameter, float * best_average,
	int * v1, int * v2){
	graph * g;
	int * x = new int[nodes * nodes];
	memcpy(x, gr, nodes * nodes * sizeof(int));
	for (vector<int>::iterator a = search_nodes.begin(); a != search_nodes.end(); a++){
		int j = *a;
		if (!x[nodes * i + j] &&
			i - (i % source_size) != j - (j % source_size)){
			x[nodes * i + j] = 1;
			x[nodes * j + i] = 1;
			g = new graph(x, nodes);
			g->calcDist();
			changing_diameter.lock();
			if (g->diameter() < *best_diameter){
				*best_diameter = g->diameter();
				*best_average = g->average();
				*v1 = i;
				*v2 = j;
			}
			else if (g->diameter() == *best_diameter && g->average() < *best_average){
				*best_average = g->average();
				*v1 = i;
				*v2 = j;
			}
			changing_diameter.unlock();
			delete g;
			x[nodes * i + j] = 0;
			x[nodes * j + i] = 0;
		}
	}
	delete[] x;
}

int main(){
	//connect_graphs(&optimal8node[0][0], 8, 64, 6, "Graph8x8Redone.txt");
	connect_graphs(&optimal16node[0][0], 16, 256, 8, "Graph16x16.txt");
	int * x = load_file("NewGraph.txt", 64);
	connect_graphs(x, 64, 256, 8, "Graph64x4.txt");
	getchar();
}

void random_walk(bool * x, int nodes){
	float ** rw = new float*[nodes];
	for (int i = 0; i < nodes; i++) rw[i] = new float[nodes];
	for (int i = 0; i < nodes; i++){
		for (int j = 0; j < nodes; j++){
			if (x[nodes * i + j]) rw[i][j] = 1.0 / 1.22;
			else rw[i][j] = 0;
		}
	}

	float * f = new float[nodes];
	for (int i = 0; i < nodes; i++) f[i] = 1.0 / nodes;
	int * iterations = new int[nodes];
	for (int i = 0; i < nodes; i++) iterations[i] = 0;
	for (int i = 0; i < nodes; i++){
		float * v = new float[nodes];
		for (int n = 0; n < nodes; n++) v[n] = 0.0f;
		v[i] = 1.0f;
		float error = 1000000.0f;
		int k = 0;
		while (error > .0001){
			float * temp = matrix_mult(rw, v, nodes);
			delete[] v;
			v = temp;
			error = 0;
			for (int j = 0; j < nodes; j++){
				error += (v[j] - f[j]) * (v[j] - f[j]);
			}
			iterations[i]++;
		}
	}
	for (int i = 0; i < nodes; i++){
		cout << "Iterations for node " << i << " " << iterations[i] << "\n";
	}
	delete[] iterations;
}

/*void create_graph(int nodes)
{
	bool x[nodes][nodes];
	for (int i = 0; i < nodes; i++){
		for (int j = 0; j < nodes; j++){
			x[i][j] = true;
		}
		x[i][i] = false;
	}
	int deg[nodes];
	for (int i = 0; i < nodes; i++){
		deg[i] = nodes - 1;
	}
	int max_degree = nodes - 1;
	while (true){
		int best_edge[2] = { -1, -1 };
		int best_diameter = 1000000000;
		float best_average = 1000000000;
		for (int i = 0; i < nodes; i++){
			for (int j = 0; j < nodes; j++){
				if (x[i][j] && deg[i] > K && deg[j] > K &&
					deg[i] >= max_degree && deg[j] >= max_degree){
					x[i][j] = false;
					x[j][i] = false;
					graph temp(&x[0][0]);
					int diameter = temp.diameter();
					if (diameter < best_diameter){
						best_diameter = diameter;
						best_edge[0] = i;
						best_edge[1] = j;
					}
					else if(diameter == best_diameter) {
						float average = temp.average();
						if (average < best_average){
							best_average = average;
							best_edge[0] = i;
							best_edge[1] = j;
						}
					}
					x[i][j] = true;
					x[j][i] = true;
				}
			}
		}
		if (best_edge[0] == -1){
			int a1, a2, b1, b2;
			for (int i = 0; i < nodes; i++){
				if (deg[i] == max_degree){
					a1 = i;
					break;
				}
			}
			for (int i = 0; i < nodes; i++){
				if (deg[i] == max_degree && i != a1){
					b2 = i;
					break;
				}
			}
			//select the other nodes as two random nodes connected to the proper nodes
			do {
				b1 = rand() % nodes;
			} while (b1 == a1 || b1 == b2 || !x[a1][b1]);
			do{
				a2 = rand() % nodes;
			} while (a2 == a1 || a2 == b1 || a2 == b2 || !x[a2][b2]);
			//perform the double adjancency switch
			x[a1][b1] = false; x[b1][a1] = false;
			x[a2][b2] = false; x[b2][a2] = false;
			x[a1][b2] = true ; x[b2][a1] = true ;
			x[b1][a2] = true ; x[a2][b1] = true ;
			best_edge[0] = a1;
			best_edge[1] = b2;
		}
		x[ best_edge[0] ][ best_edge[1] ] = false;
		x[ best_edge[1] ][ best_edge[0] ] = false;
		deg[best_edge[0]]--;
		deg[best_edge[1]]--;
		bool done = true;
		max_degree = 0;
		for (int i = 0; i < nodes; i++){
			if (deg[i] > max_degree) max_degree = deg[i];
		}
		if (max_degree <= K) break;
	}
	graph g(&x[0][0]);
	std::cout << "Diameter: " << g.diameter() << "\n";
	std::cout << "Average: " << g.average() << "\n";
	getchar();
}*/