#include <vector>
#include <fstream>
#include <cassert>

using namespace std;

class mGraph {
    public:
        int n, m;
        string basename;

        vector< vector<int> > g;
        vector< vector<int> > gT;
        vector< vector<double> > prob;

        void readNM() {
            string f = basename + "_att.txt";

            ifstream infile(f);

            infile >> n >> m;
        }

        template<bool fixed>
        void readGraph(double p_val){
            int u,v;
            double p;

            string f = basename + "_graph.txt";

            ifstream infile(f);

            g.clear();
            g.resize(n);
            gT.clear();
            gT.resize(n);
            prob.clear();
            prob.resize(n);

            int edges_found = 0;
            while (infile >> u >> v >> p) {
                if (fixed)
                    prob[u].push_back(p_val);
                else
                    prob[u].push_back(p);
                g[u].push_back(v);
                gT[v].push_back(u);
                edges_found++;
            }
            assert(m == edges_found);
        }

        mGraph(string the_basename, string prob_dist, double p_val): basename(the_basename) {
            readNM();
            if (prob_dist == string("fixed")) 
                readGraph<true>(p_val);
            else
                readGraph<false>(0.0);
        }

};
