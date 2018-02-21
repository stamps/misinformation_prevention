#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <queue>
#include <ctime>
#include <cmath>
#include <map>
#include <algorithm>

#include "sfmt/SFMT.h"
#include "mGraph.h"

using namespace std;

const double eps = .1;

class InfluenceGraph {

private:
	mGraph G;
	int& n;
	int& m;

	vector<int> c_visited;
    vector<int> visited;
    vector<int> depth;
    vector<int> has_beta;
    vector<int> beta;
	vector<int> node_sketch_count;

	deque<int> q;

    map<long,int> sketch_index;

	int k;
    int v_c;
    int curr_sketch_index;

	// hypergraph
	vector< vector<int> > h;
	vector< vector<int> > hT;

public:
	sfmt_t sfmtSeed;

	InfluenceGraph(string basename, int k_val, int v_val, string prob_dist, double p_val): G(basename, prob_dist, p_val), n(G.n), m(G.m), c_visited(n), visited(n), depth(n), has_beta(n), beta(n), node_sketch_count(n), k(k_val), v_c(v_val) {
		cout << "------------------" << endl;
		cout << "Initializing variables." << endl;

		cout << "n = " << n << "\tm = " << m << endl;

        sfmt_init_gen_rand(&sfmtSeed, int(time(NULL)));

        cout << "v_c = " << v_c << endl;
	}

	void clear_marked() {
        fill(c_visited.begin(), c_visited.end(), 0);
        fill(visited.begin(), visited.end(), 0);
        fill(has_beta.begin(), has_beta.end(), 0);
	}

	//
    //  Standard BFS traversal
    //
	int BFS(int v, int sketch) {
		int bfs_size = 0;

		q.clear();
		q.push_back(v);
		c_visited[v] = sketch + 1;
        has_beta[v] = sketch + 1;
        beta[v] = 0;
        bfs_size++;

		while (!q.empty()) {
			int u = q.front(); q.pop_front();

			for(int j = 0; j < G.g[u].size(); j++) {

				int w = G.g[u][j];
				double flip = double(sfmt_genrand_uint32(&sfmtSeed))/double(RAND_MAX)/2;

                if (flip > G.prob[u][j]) continue;
                if (c_visited[w] == sketch+1) continue;

                q.push_back(w);
                c_visited[w] = sketch + 1;
                has_beta[w] = sketch + 1;
                beta[w] = beta[u] + 1;
                bfs_size++;
			}
		}

		return bfs_size;
	}

    //
    //  Algorithm 2 - modified BFS for constructing an RRC set
    //
    template<bool addHypEdge>
    int mBFS(int u, int sketch, int distance) {
        int bfs_size = 0;

        q.clear();
        q.push_back(u);
        visited[u] = sketch + 1;
        depth[u] = 0;

        if(addHypEdge) {
            hT.push_back(vector<int>());
            sketch_index[sketch] = curr_sketch_index;
            hT[curr_sketch_index].push_back(u);
            node_sketch_count[u]++;
        }

        while (!q.empty()) {
            int w = q.front(); q.pop_front();

            if (depth[w] == distance) continue;
            if (has_beta[w] == sketch+1 && beta[w] == 0) continue;

            for(int z : G.gT[w]) {
                bfs_size++;

                if (visited[z] == sketch+1) continue;

                if (has_beta[w] == sketch+1 && has_beta[z] == sketch+1) {
                    if (beta[w] - 1 < beta[z]) {
                        beta[z] = beta[w] - 1;
                    }
                } else if (has_beta[w] == sketch+1 && has_beta[z] < sketch+1) {
                    has_beta[z] = sketch+1;
                    beta[z] = beta[w] - 1;
                }

                q.push_back(z);
                visited[z] = sketch + 1;
                depth[z] = depth[w] + 1;

                if(addHypEdge) {
                    hT[curr_sketch_index].push_back(z);
                    node_sketch_count[z]++;
                }
            }
        }

        if(addHypEdge)
            curr_sketch_index++;

        return bfs_size;
    }

	//
    //  Main Sampling Loop
    //
	void generate_hypergraph(long theta) {

        long successes = 0;
        long avg_bfs_size = 0;
        long avg_success_bfs_size = 0;
        long avg_distance = 0;
        long num_empty = 0;

		// initialize hypergraph
	    h.clear();
		h.resize(n);

		hT.clear();
		//hT.resize(theta);

        sketch_index.clear();
        curr_sketch_index = 0;

		// reset sketch counter
        fill(node_sketch_count.begin(), node_sketch_count.end(), 0);

		clock_t start = clock();

		clear_marked();

		// generate theta RRC sets
		for(int i = 0; i < theta; i++) {
            long b = BFS(v_c, i);
            avg_bfs_size += b;
			
            // generate a random node u to start mBFS from
			int u = sfmt_genrand_uint32(&sfmtSeed) % n;
            while(u == v_c) u = sfmt_genrand_uint32(&sfmtSeed) % n;

            // if u not in I(v_c) then it is not eligible to be saved -- continue
            if (c_visited[u] < i+1) {
                num_empty++;
                continue;
            }

            successes++;
            avg_success_bfs_size += b;

            // determine distance from v_c to u
            int distance = beta[u];
            avg_distance += distance;

			mBFS<true>(u, i, distance);
		}

		cout << "Time to generate hypergraph: " << (clock() - start) / (double)(CLOCKS_PER_SEC) << " s" << endl;
		start = clock();

        long avg_mbfs_size = 0;
        for ( const auto &p : sketch_index ) {
            int curr_sketch_loc = p.second;
            for(int t : hT[curr_sketch_loc]) {
                h[t].push_back(p.first);
                avg_mbfs_size++;
            }
        }

		cout << "Time to invert hypergraph: " << (clock() - start) / (double)(CLOCKS_PER_SEC) << " s" << endl;

		//cout << "hT size = " << hT.size() << endl;
		cout << "Succeses = " << successes << endl;

		avg_bfs_size /= (double)(theta);
        avg_success_bfs_size /= (double)(successes);
        avg_mbfs_size /= (double)(successes);
        avg_distance /= (double)(successes);
        double pct_empty = (double)(num_empty) / theta;
		cout << "Average BFS Size = " << avg_bfs_size << endl;
        cout << "Average success BFS Size = " << avg_success_bfs_size << endl;
        cout << "Average mBFS Size = " << avg_mbfs_size << endl;
        cout << "Average distance = " << avg_distance << endl;
        //cout << "Number of empty RRC sets = " << num_empty << endl;
        cout << "Percentage of empty RRC sets = " << pct_empty << endl;
	}

	//
    //  Algorithm 1 - NodeSelection
    //
	double node_selection(long num_sketches) {

		clock_t start = clock();
        
        //vector<int> sketches_removed(num_sketches, 0);
        map<long,int> sketches_removed;
        for ( const auto &p : sketch_index )
            sketches_removed[p.first] = 0;

        vector<int> max_nodes;
        vector<long> prevention;
        int max_node;
        long max_node_count;

        max_nodes.clear();
        prevention.clear();

        node_sketch_count[v_c] = 0;

        for (int i = 0; i < k; i++) {
            // 1. identify max node
            max_node = -1;
            max_node_count = -1;
            for (int i = 0; i < n; i++) {
                if (node_sketch_count[i] > max_node_count) {
                    max_node = i; 
                    max_node_count = node_sketch_count[i];
                }
            }
            //cout << "Max node: " << max_node << "\tCount: " << max_node_count << endl;
            max_nodes.push_back(max_node);
            prevention.push_back(max_node_count);

            // 2. identify sketches that max node participates in via inverted index
            for (int curr_sketch : h[max_node]) {
                if (!sketches_removed[curr_sketch]) {
                    int curr_sketch_loc = sketch_index[curr_sketch];
                    // 3. remove corresponding sketches: decrement node_sketch_count of each node in the sketch
                    for (int curr_node : hT[curr_sketch_loc]) {
                        node_sketch_count[curr_node]--;
                    }
                    // mark sketch as removed
                    sketches_removed[curr_sketch] = 1;
                }
            }
        } 

        cout << "Time to get seed set: " << (clock() - start) / (double)(CLOCKS_PER_SEC) << " s" << endl;

        // compute expected prevention
        double expected_prevention = 0.0;
        for (long inf : prevention) {
            expected_prevention += inf;
        }
        // E[pi(S)]= n * F_R
        expected_prevention = (expected_prevention / num_sketches) * n;
        cout << "Expected prevention: " << expected_prevention << endl;

        cout << "Seed set: ";
        for(auto item : max_nodes)
            cout << item << " ";
        cout << endl;

        return expected_prevention;
	}

	//
    //  Algorithm 3 - KptEstimation
    //
    double kpt_estimation() {

        int outer_loop, c_i;
        double sum, subwidth, lhs, rhs;

        int bfs_num = 0;

        outer_loop = log(n) / log(2);
        for (int i = 1; i < outer_loop; i++) {
            sum = 0.0;
            c_i = (6 * log(n) + 6 * log(log(n)/ log(2)) ) * pow(2,i);
            for (int j = 0; j < c_i; j++) {
                BFS(v_c, bfs_num);
                // generate a random node u to start mBFS from
                int u = sfmt_genrand_uint32(&sfmtSeed) % n;
                while(u == v_c) u = sfmt_genrand_uint32(&sfmtSeed) % n;
                // if u not in I(v_c) then it is not eligible to be saved -- continue
                if (c_visited[u] < bfs_num+1) {
                    bfs_num++;
                    continue;
                }
                // determine distance from v_c to u
                int distance = beta[u];
                subwidth = (double)mBFS<false>(u, bfs_num++, distance);
                sum += 1.0 - pow( (1.0 - subwidth/m), k );
            }

            lhs = sum / c_i;
            rhs = 1.0 / pow(2,i);
            if (lhs > rhs) {
                return n * (sum / (2.0 * c_i));
            }
        }

        return 1.0;
    }

    //
    //  Algorithm 3 - KptEstimation
    //
    double refine_kpt(double kpt_star) {

        double eps_prime = 5.0 * pow( (eps * eps)/(k+1), 1.0/3.0 );

        cout << "eps_prime = " << eps_prime << endl;

        long theta_prime = (2.0 + eps_prime) * ( n * log(n) ) / ( eps_prime * eps_prime * kpt_star);

        cout << "theta_prime = " << theta_prime << endl;
        
        generate_hypergraph(theta_prime);

        double kpt_prime = node_selection(theta_prime);
        
        kpt_prime /= 1.0 + eps_prime;

        cout << "kpt_prime = " << kpt_prime << endl;
        
        return max(kpt_prime, kpt_star);
    }

    double logcnk(int n, int k){
        double val = 0;
        for(int i = n-k+1; i <= n; i++) {
            val += log(i);
        }
        for(int i = 1; i <= k; i++) {
            val -= log(i);
        }
        return val;
    }

    long compute_theta(double kpt_plus) {
    	return (8.0 + 2.0 * eps) * n * ( log(n) + log(2) + logcnk(n, k) ) / ( eps * eps * kpt_plus);
    }
};

int select_adversary(string basename, int topk) {
    string f = basename + "_top" + to_string(topk) + ".txt";
    ifstream infile(f);

    int v;
    vector<int> candidates;
    while (infile >> v) {
        candidates.push_back(v);
    }

    sfmt_t sfmtSeed;
    sfmt_init_gen_rand(&sfmtSeed, int(time(NULL)));

    int c = sfmt_genrand_uint32(&sfmtSeed) % candidates.size();
    return candidates[c];
}

int main(int argc, char ** argv) {
	if (argc < 5){
		cout << "Usage: " << argv[0] << " <basename> <k> <topk> <prob_dist> <p>" << endl;
		return 1;
	}

	string basename = string(argv[1]);
	int k = atoi(argv[2]);
    int topk = atoi(argv[3]);
    string prob_dist = string(argv[4]);
    double p_val = 0.0;
    if(argv[4] == string("fixed")) p_val = atof(argv[5]);

	clock_t start = clock();

    int v_c = select_adversary(basename, topk);

	InfluenceGraph IG(basename, k, v_c, prob_dist, p_val);

	// KPT Estimation
	cout << "------------------" << endl;
	cout << "KPT Estimation" << endl;
    clock_t kpt_start = clock();
	double kpt_star = IG.kpt_estimation();
    cout << "Time to do kpt estimation: " << (clock() - kpt_start) / (double)(CLOCKS_PER_SEC) << " s" << endl;
	cout << "kpt_star = " << kpt_star << endl;

	// Refine KPT
	cout << "------------------" << endl;
	cout << "Refine KPT" << endl;
	double kpt_plus = IG.refine_kpt(kpt_star);
    cout << "------------------" << endl;
	cout << "kpt_plus = " << kpt_plus << endl;

	// Compute Theta
    cout << "------------------" << endl;
	long theta = IG.compute_theta(kpt_plus);
	cout << "theta = " << theta << endl;

	// Generate Hypergraph
	cout << "------------------" << endl;
	cout << "Generating Hypergraph" << endl;
	IG.generate_hypergraph(theta);

	// Compute Seed Set
	cout << "------------------" << endl;
	cout << "Computing Seed Set with k=" << k << endl;
	double prevention = IG.node_selection(theta);

	cout << "Total Time: " << (clock() - start) / (double)(CLOCKS_PER_SEC) << " s" << endl;

}
