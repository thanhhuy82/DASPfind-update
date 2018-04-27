#include<bits/stdc++.h>
#define min_drug_similarity 0.5
#define min_target_similarity 0.5
#define alpha 2.26

using namespace std;

typedef pair<int, double> pid;
typedef vector<vector<double> > vvd;

vvd drugs, targets, interactions, result;
vector<pid> g[100005];
vector<int> visited;
int nDrugs, nTargets, nVertex;

int numberOfDrugs() {
    ifstream fi("NR_Drugs_IDs.txt");
    string s;
    getline(fi, s);
    int ret=0;
    while (!fi.eof()) {
        getline(fi, s);
        ret++;
    }
    return ret;
}

int numberOfTargets() {
    ifstream fi("NR_Targets_IDs.txt");
    string s;
    getline(fi, s);
    int ret=0;
    while (!fi.eof()) {
        getline(fi, s);
        ret++;
    }
    return ret;
}

void init() {
    nDrugs = numberOfDrugs();
    nTargets = numberOfTargets();
    nVertex = nDrugs + nTargets;
    drugs.resize(nDrugs, vector<double>(nDrugs, 0));
    targets.resize(nTargets, vector<double>(nTargets, 0));
    interactions.resize(nDrugs, vector<double>(nTargets, 0));
    result.resize(nDrugs, vector<double>(nTargets, 0));
    visited.resize(nVertex);
}

void readDrugsSimilarity() {
    
    ifstream fi("NR_Drugs_Similarity.txt");
    for (size_t i = 0; i < nDrugs; i++) {
        for (size_t j = 0; j < nDrugs; j++) {
            fi >> drugs[i][j];
            
        }
        
    }
}

void readTargetsSimilarity() {
    
    ifstream fi("NR_Targets_Similarity.txt");
    for (size_t i = 0; i < nTargets; i++) {
        for (size_t j = 0; j < nTargets; j++) {
            fi >> targets[i][j];
        }
    }
}

void readDrugTargetInteractions() {
    
    ifstream fi("NR_DrugTarget_interactions.txt");
    for (size_t i = 0; i < nDrugs; i++) {
        for (size_t j = 0; j < nTargets; j++) {
            fi >> interactions[i][j];
        }
    }
}

void write_result(vvd _vvd) {
    ofstream fo("output.txt");
    for (size_t i = 0; i < _vvd.size(); i++) {
        for (size_t j = 0; j < _vvd[i].size() ; j++) {
            fo << setprecision(18) << fixed<< _vvd[i][j] << ' ';
        }
        fo << endl;
    }
}

void buildGraph() {
    for (size_t i = 0; i < nDrugs; i++) {
        for (size_t j = 0; j < nDrugs; j++) {
            if (i != j && drugs[i][j] > min_drug_similarity) {
                g[i].push_back(make_pair(j, drugs[i][j]));
            }
        }
    }
    
    for (size_t i = nDrugs; i < nDrugs+nTargets; i++) {
        for (size_t j = nDrugs; j < nDrugs+nTargets; j++) {
            if (i != j && targets[i-nDrugs][j-nDrugs] > min_target_similarity) {
                g[i].push_back(make_pair(j, targets[i-nDrugs][j-nDrugs]));
            }
        }
    }
    
    for (size_t i = 0; i < nDrugs; i++) {
        for (size_t j = 0; j < nTargets; j++) {
            if (interactions[i][j] > 0) {
                g[i].push_back(make_pair(j+nDrugs, 1.0));
                g[j+nDrugs].push_back(make_pair(i, 1.0));
            }
        }
    }
}

double score=0.0;
void traverse(int u, int des, int cnt, double sc) {
    if (cnt > 3){
        return;
    }
    if (u==des) {
        score+=pow(sc, alpha*cnt);
        return;
    }
    for (auto e : g[u]) {
        int v = e.first;
        double cost = e.second;
        if (!visited[v]) {
            visited[v] = 1;
            traverse(v, des, cnt+1, sc*cost);
            visited[v] = 0;
        }
    }
}

void mainTask() {
    for (int i = 0; i < nDrugs; i++) {
        for (int j = nDrugs; j < nDrugs + nTargets; j++) {
            for (int k = 0; k < nVertex; k++) {
                visited[k] = 0;
            }
            visited[i] = 1;
            traverse(i, j, 0, 1.0);
            result[i][j-nDrugs] = score;
            score = 0.0;
        }
    }
}

int main() {
    init();
    readDrugsSimilarity();
    readTargetsSimilarity();
    readDrugTargetInteractions();
    buildGraph();
    mainTask();
    write_result(result);
    return 0;
}
