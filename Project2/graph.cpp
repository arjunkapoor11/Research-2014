#include <iostream>
#include <stdlib.h>
#include <cstring>
#include "graph.h"

int * combineGraphs(int * graph, int size, int nodes) {
	int * g = (int*)malloc(nodes * nodes * sizeof(int));
	memset(g, 0, nodes * nodes * sizeof(int));
	for (int x = 0; x < nodes / size; x++){
		for (int i = 0; i < size; i++){
			for (int j = 0; j < size; j++){
				g[nodes * size * x + size * x + nodes * i + j] = graph[size * i + j];
			}
		}
	}
	return g;
}

void das(bool * g, int a1, int a2, int b1, int b2, int nodes){
	g[a1 * nodes + b1] = false;
	g[b1 * nodes + a1] = false;
	g[a2 * nodes + b2] = false;
	g[b2 * nodes + a2] = false;
	g[a1 * nodes + b2] = true;
	g[b2 * nodes + a1] = true;
	g[a2 * nodes + b1] = true;
	g[b1 * nodes + a2] = true;
}

void graph::calcDist(){
	if (distsCalculated) return;
	for (int i = 0; i < nodes; i++){
		for (int j = 0; j < nodes; j++){
			dist[i][j] = 10000000;
			dist[j][i] = 10000000;
		}
		dist[i][i] = 0;
	}
	//Q contains the yet to be searched vertecies
	vector<int> neighbors;
	vector<int> next_dist_nodes, current_dist_nodes;
	for (int r = 0; r < nodes; r++){
		int current_dist = 1;
		current_dist_nodes.push_back(r);
		while (true){
			neighbors.clear();
			for (vector<int>::iterator i = current_dist_nodes.begin(); i != current_dist_nodes.end(); i++){
				for (int j = 0; j < nodes; j++){
					if (adj[*i][j]) neighbors.push_back(j);
				}
			}
			for (vector<int>::iterator i = neighbors.begin(); i != neighbors.end(); i++){
				if (dist[r][*i] > current_dist){
					dist[r][*i] = current_dist;
					next_dist_nodes.push_back(*i);
				}
			}
			current_dist_nodes.clear();
			current_dist_nodes = next_dist_nodes;
			next_dist_nodes.clear();
			current_dist++;
			if (current_dist_nodes.empty()) break;
		}
	}
	for (int i = 0; i < nodes; i++){
		for (int j = i; j < nodes; j++){
			dist[i][j] = dist[j][i];
		}
	}
	distsCalculated = true;
}

int graph::diameter(){
	calcDist();
	int d = 0;
	for (int i = 0; i < nodes; i++){
		for (int j = i; j < nodes; j++){
			if (dist[i][j] > d) d = dist[i][j];
		}
	}
	return d;
}

float graph::average(){
	calcDist();
	int sum = 0;
	for (int i = 0; i < nodes; i++){
		for (int j = 0; j < nodes; j++){
			sum += dist[i][j];
		}
	}
	return ((float)sum) / ((float)nodes) / ((float)nodes - 1);
}

float graph::node_average(int node){
	calcDist();
	float sum = 0;
	for (int i = 0; i < nodes; i++){
		sum += dist[node][i];
	}
	return sum / (float)(nodes - 1);
}

int graph::node_diameter(int node){
	calcDist();
	int d = 0;
	for (int i = 0; i < nodes; i++){
		if (dist[node][i] > d) d = dist[node][i];
	}
	return d;
}

void graph::correct(){
	calcDist();
	int v1, v2, d;
	d = diameter();
	bool done = false;
	for (int i = 0; i < nodes; i++){
		for (int j = 0; j < nodes; j++){
			if (dist[i][j] = d){
				v1 = i;
				v2 = j;
				done = true;
				break;
			}
		}
		if (done) break;
	}
	vector<int> v1Neighbors, v2Neighbors;
	for (int i = 0; i < nodes; i++){
		if (adj[v1][i]) v1Neighbors.push_back(i);
		if (adj[v2][i]) v2Neighbors.push_back(i);
	}
	bool * x = new bool[nodes * nodes];
	bool * best = new bool[nodes * nodes];
	int best_d = 10000000;
	float best_a = 1000000.0f;
	for (vector<int>::iterator i = v1Neighbors.begin(); i != v1Neighbors.end(); i++){
		for (vector<int>::iterator j = v2Neighbors.begin(); j != v2Neighbors.end(); j++){
			for (int r = 0; r < nodes; r++) for (int c = 0; c < nodes; c++) x[r * nodes + c] = adj[r][c];
			das(x, v1, *j, *i, v2, nodes);
			graph g(x, nodes);
			if (g.diameter() <= best_d && g.average() < best_a){
				memcpy(best, x, nodes * nodes * sizeof(bool));
				best_d = g.diameter();
				best_a = g.average();
			}
		}
	}
	delete[] x;
	memcpy(&adj[0][0], best, nodes * nodes * sizeof(bool));
	delete[] best;
}

void graph::print(){
	for (int i = 0; i < nodes; i++){
		for (int j = 0; j < nodes; j++){
			if (adj[i][j]) std::cout << "1";
			else std::cout << "0";
		}
		std::cout << "\n";
	}
}