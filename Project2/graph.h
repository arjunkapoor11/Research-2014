#include <vector>
#include <atomic>
#include <thread>
#include <mutex>
using std::vector;
using namespace std;

#define EDGE_SEARCH_SIZE 100
#define MAX_THREADS 32

int * combineGraphs(int * graph, int, int);
void das(bool *, int, int, int, int, int);
void compare_graph(int * gr, int nodes, int source_size,
	int i, vector<int> search_nodes,
	int * best_diameter, float * best_average,
	int * v1, int * v2);

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

	~graph(){
		for (int i = 0; i < nodes; i++){
			delete[] adj[i];
			delete[] dist[i];
		}
		delete[] adj;
		delete[] dist;
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