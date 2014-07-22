#include <iostream>
#include <stdlib.h>
#include <cstring>
#include "graph.h"

int optimal8node[8][8] = { { 0, 1, 1, 1, 0, 0, 0, 0 },
{ 1, 0, 0, 0, 1, 0, 1, 0 },
{ 1, 0, 0, 0, 1, 1, 0, 0 },
{ 1, 0, 0, 0, 0, 1, 0, 1 },
{ 0, 1, 1, 0, 0, 0, 0, 1 },
{ 0, 0, 1, 1, 0, 0, 1, 0 },
{ 0, 1, 0, 0, 0, 1, 0, 1 },
{ 0, 0, 0, 1, 1, 0, 1, 0 } };

int * combineGraphs(int * graph) {
	int * g = (int*)malloc(64 * 64 * sizeof(int));
	memset(g, 0, 64 * 64 * sizeof(int));
	for (int x = 0; x < 8; x++){
		for (int i = 0; i < 8; i++){
			for (int j = 0; j < 8; j++){
				g[64 * 8 * x + 8 * x + 64 * i + j] = graph[8 * i + j];
			}
		}
	}
	return g;
}

void graph::calcDist(){
	if (distsCalculated) return;
	for (int i = 0; i < 8; i++){
		for (int j = 0; j < 8; j++){
			dist[i][j] = 10000000;
		}
	}
	//Q contains the yet to be searched vertecies
	vector<int> Q, neighbors;
	int u, v;
	int altDistance;
	for (int r = 0; r < 8; r++){
		dist[r][r] = 0;
		Q.clear();
		for (int k = 0; k < r; k++) Q.push_back(k);
		while (Q.size() > 0){
			u = Q[0];
			for (std::vector<int>::iterator k = Q.begin(); k != Q.end(); k++){
				if (dist[r][*k] < dist[r][u]){
					u = *k;
				}
			}
			for (std::vector<int>::iterator i = Q.begin(); i != Q.end(); i++)
			if (*i == u) {
				Q.erase(i);
				break;
			}
			neighbors.clear();
			for (int k = 0; k < 8; k++) if (adj[u][k]) neighbors.push_back(k);
			while (neighbors.size() > 0){
				v = neighbors.back();
				neighbors.pop_back();
				altDistance = dist[r][u] + 1;
				if (altDistance < dist[r][v]) dist[r][v] = altDistance;
			}
		}
	}
	distsCalculated = true;
}