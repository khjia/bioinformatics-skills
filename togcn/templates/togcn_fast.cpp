// Algorithm-equivalent reimplementation of TO-GCN_single (adjacency-list BFS).
// Semantics copied from TO-GCN-master/Single_Time-series_data/source_code/TO-GCN_single.cpp:
//   - Pearson r between every pair of nodes (i != j)
//   - directed edge query_i -> target_j when r >= cutoff
//   - BFS: seeds = level 1, neighbours level 2, ...  (output = level - 1)
//   - edges are consumed once (both directions), i.e. BFS runs on the undirected graph
// Output files identical to the original: Node_level.csv, Node_relation.csv
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

int main(int argc, char **argv) {
    if (argc != 5) {
        printf("\nUsage: togcn_fast No_of_features expression_table_file seed_table_file Cutoff\n\n");
        return 1;
    }
    int K = atoi(argv[1]);
    const char *expr_file = argv[2];
    const char *seed_file = argv[3];
    double cutoff = atof(argv[4]);

    // ---- read matrix ----
    FILE *fp = fopen(expr_file, "r");
    if (!fp) { printf("cannot open %s\n", expr_file); return 1; }
    std::vector<std::string> ids;
    std::vector<double> raw;
    char buf[64];
    std::vector<double> tmp(K);
    while (fscanf(fp, "%63s", buf) == 1) {
        ids.push_back(std::string(buf));
        for (int i = 0; i < K; i++) {
            if (fscanf(fp, "%lf", &tmp[i]) != 1) { printf("bad matrix\n"); return 1; }
        }
        raw.insert(raw.end(), tmp.begin(), tmp.end());
    }
    fclose(fp);
    int N = (int)ids.size();
    printf("Number of nodes: %d\n", N);
    printf("Number of features: %d\n", K);

    // ---- standardise rows (genes) -> z ----
    std::vector<double> Z((size_t)N * K, 0.0);
    for (int i = 0; i < N; i++) {
        const double *x = &raw[(size_t)i * K];
        double s = 0.0;
        bool allzero = true;
        for (int k = 0; k < K; k++) { s += x[k]; if (x[k] != 0.0) allzero = false; }
        if (allzero) continue;                       // reference: SUM_X == 0 -> r = 0
        double m = s / K;
        double v = 0.0;
        for (int k = 0; k < K; k++) { double d = x[k] - m; v += d * d; }
        double sd = sqrt(v / (K - 1));
        if (!(sd > 0.0)) continue;                   // constant row -> r = 0
        for (int k = 0; k < K; k++) Z[(size_t)i * K + k] = (x[k] - m) / sd;
    }
    raw.clear(); raw.shrink_to_fit();

    // ---- pairwise Pearson, keep r >= cutoff ----
    std::vector<int> eu, ev;      // undirected edge list (u < v)
    std::vector<float> er;
    eu.reserve(1 << 20); ev.reserve(1 << 20); er.reserve(1 << 20);
    long long npos = 0, nneg = 0;
    for (int i = 0; i < N; i++) {
        const double *zi = &Z[(size_t)i * K];
        for (int j = i + 1; j < N; j++) {
            const double *zj = &Z[(size_t)j * K];
            double d = 0.0;
            for (int k = 0; k < K; k++) d += zi[k] * zj[k];
            double r = d / (K - 1);
            if (r >= cutoff) { eu.push_back(i); ev.push_back(j); er.push_back((float)r); npos++; }
            else if (r <= -cutoff) nneg++;
        }
    }
    printf("Edges (r >= %.2f): %lld ; negative pairs (r <= -%.2f): %lld\n", cutoff, npos, cutoff, nneg);

    // ---- CSR adjacency (undirected) ----
    std::vector<int> deg(N + 1, 0);
    for (size_t e = 0; e < eu.size(); e++) { deg[eu[e] + 1]++; deg[ev[e] + 1]++; }
    for (int i = 0; i < N; i++) deg[i + 1] += deg[i];
    std::vector<int> adj(deg[N]);
    std::vector<char> visited(N, 0);
    {
        std::vector<int> pos(deg.begin(), deg.end() - 1);
        for (size_t e = 0; e < eu.size(); e++) { adj[pos[eu[e]]++] = ev[e]; adj[pos[ev[e]]++] = eu[e]; }
    }

    // ---- seeds ----
    std::vector<int> seeds;
    fp = fopen(seed_file, "r");
    if (!fp) { printf("cannot open %s\n", seed_file); return 1; }
    while (fscanf(fp, "%63s", buf) == 1) {
        int found = -1;
        for (int i = 0; i < N; i++) if (ids[i] == buf) { found = i; break; }
        if (found >= 0) seeds.push_back(found);
        else printf("seed %s not found\n", buf);
    }
    fclose(fp);
    printf("Number of seeds: %d\n", (int)seeds.size());
    printf("Cutoffs for correlation: %.3f\n\n", cutoff);

    // ---- BFS ----
    std::vector<int> lvl(N, -1);
    std::vector<int> frontier, next;
    for (size_t s = 0; s < seeds.size(); s++) if (lvl[seeds[s]] < 0) { lvl[seeds[s]] = 1; frontier.push_back(seeds[s]); }
    int level = 1;
    while (!frontier.empty()) {
        level++;
        next.clear();
        for (size_t f = 0; f < frontier.size(); f++) {
            int u = frontier[f];
            for (int p = deg[u]; p < deg[u + 1]; p++) {
                int v = adj[p];
                if (lvl[v] < 0) { lvl[v] = level; next.push_back(v); }
            }
        }
        frontier.swap(next);
    }

    // ---- output ----
    fp = fopen("Node_level.csv", "w");
    fprintf(fp, "TF gene ID,level in GCN\n");
    int assigned = 0, maxl = 0;
    for (int i = 0; i < N; i++) {
        if (lvl[i] > 0) { fprintf(fp, "%s,%d\n", ids[i].c_str(), lvl[i] - 1); assigned++; if (lvl[i] - 1 > maxl) maxl = lvl[i] - 1; }
    }
    fclose(fp);

    fp = fopen("Node_relation.csv", "w");
    fprintf(fp, "node_1 ID, node_2 ID, PCC_value\n");
    for (size_t e = 0; e < eu.size(); e++) {
        if (lvl[eu[e]] > 0 && lvl[ev[e]] > 0)
            fprintf(fp, "%s,%s,%lf\n", ids[eu[e]].c_str(), ids[ev[e]].c_str(), (double)er[e]);
    }
    fclose(fp);
    printf("Assigned nodes: %d / %d, levels 0..%d\nDone\n", assigned, N, maxl);
    return 0;
}
