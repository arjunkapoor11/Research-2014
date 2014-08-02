#include <vector>
using std::vector;

int * combineGraphs(int * graph, int, int);
void das(bool *, int, int, int, int, int);

class graph
{
private:
	int nodes;
	bool ** adj;
	int ** dist;
	bool distsCalculated;
public:
	graph(bool * a, int numNodes){
		distsCalculated = false;
		nodes = numNodes;
		adj = new bool *[nodes];
		dist = new int *[nodes];
		for (int i = 0; i < nodes; i++){
			adj[i] = new bool[nodes];
			dist[i] = new int[nodes];
		}
		for (int i = 0; i < nodes; i++){
			for (int j = 0; j < nodes; j++){
				adj[i][j] = a[i * nodes + j];
			}
		}
	}
	graph(int * a, int numNodes){
		distsCalculated = false;
		nodes = numNodes;
		adj = new bool *[nodes];
		dist = new int *[nodes];
		for (int i = 0; i < nodes; i++){
			adj[i] = new bool[nodes];
			dist[i] = new int[nodes];
		}
		for (int i = 0; i < nodes; i++)
		for (int j = 0; j < nodes; j++){
			if (a[nodes * i + j]) adj[i][j] = true;
			else adj[i][j] = false;
		}
	}

	void calcDist();
	int diameter();
	float average();
	float node_average(int);
	int node_diameter(int);

	bool operator== (graph a);
	void print();
	void correct();
};