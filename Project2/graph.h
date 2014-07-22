#include <vector>

int * combineGraphs(int * graph);

class graph
{
private:
	bool adj[64][64];
	int dist[64][64];
	bool distsCalculated;
public:
	graph(bool * a){
		distsCalculated = false;
		for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			adj[i][j] = a[8 * i + j];
	}
	void calcDist();
	int diameter();
	float average();

	bool operator== (graph a);
	void print();
};