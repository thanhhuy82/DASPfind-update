#include<bits/stdc++.h>
#define min_drug_similarity 0.5
#define min_target_similarity 0.5
#define alpha 2.26
#define fdrugID "NR_Drugs_IDs.txt"
#define ftargetID "NR_Targets_IDs.txt"
#define fdtInterraction "NR_DrugTarget_interactions.txt"
#define fdrugSimilarity "NR_Drugs_Similarity.txt"
#define ftargetSimilarity "NR_Targets_Similarity.txt"
using namespace std;

typedef pair<int, double> pid;
typedef vector<vector<double> > vvd;
typedef pair<pair<int,int>,double> ppd;
vvd drugs, targets, interactions, result;
vector<pid> g[100005];
vector<unsigned int> visited,outcome;
vector<ppd> lstpredicton;
vector<string> targetsname,drugnames;
int nDrugs, nTargets, nVertex;
//get numberOfDrugs
int numberOfDrugs() {
    ifstream fi(fdrugID);
    string s;
    getline(fi, s);
    int ret=0;
    while (!fi.eof()) {
        getline(fi, s);
        ret++;
        drugnames.push_back(s);
    }
    return ret;
}
//get numberOfTargets
int numberOfTargets() {
    ifstream fi(ftargetID);
    string s;
    getline(fi, s);
    int ret=0;
    while (!fi.eof()) {
        getline(fi, s);
        ret++;
        targetsname.push_back(s);
    }
    return ret;
}
// init all state and vector
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
//read DrugsSimilarity
void readDrugsSimilarity() {

    ifstream fi(fdrugSimilarity);
    for (size_t i = 0; i < nDrugs; i++) {
        for (size_t j = 0; j < nDrugs; j++) {
            fi >> drugs[i][j];

        }

    }
}
//read TargetsSimilarity
void readTargetsSimilarity() {

    ifstream fi(ftargetSimilarity);
    for (size_t i = 0; i < nTargets; i++) {
        for (size_t j = 0; j < nTargets; j++) {
            fi >> targets[i][j];
        }
    }
}
//read DrugTargetInteractions
void readDrugTargetInteractions() {

    ifstream fi(fdtInterraction);
    for (size_t i = 0; i < nDrugs; i++) {
        for (size_t j = 0; j < nTargets; j++) {
            fi >> interactions[i][j];
        }
    }
}
void findnewInteraction(vvd _result,vvd _interaction)
{

    lstpredicton.resize(nDrugs+nTargets);
    outcome.resize(nDrugs+nTargets);
    for(size_t i=0;i<nDrugs+nTargets; i++)
        outcome[i]=0;
    for(size_t i=0; i<_interaction.size(); i++)
    {
        for(size_t j=0; j<_interaction[i].size(); j++)
        {

                lstpredicton.push_back(make_pair(make_pair(i,j),_result[i][j]));
        }
    }

}
//write to matrix output
void write_result(vvd _vvd) {
    ofstream fo("output.txt");
    for (size_t i = 0; i < _vvd.size(); i++) {
        for (size_t j = 0; j < _vvd[i].size() ; j++) {
            fo << setprecision(18) << fixed<< _vvd[i][j] << ' ';
        }
        fo << endl;
    }
}
// make linkList graph
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
// traverse graph DFS
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
//write to CVS format
void write_CVS(vvd _vvd)
{
    size_t i;
    ofstream fo("result.csv");
    fo << "\t,";
    for(i=0; i<nTargets-1; i++)
        fo << targetsname[i]<<",";
    fo << targetsname[i];
    fo << endl;
    //write line
     for (size_t i = 0; i < _vvd.size(); i++) {
        fo << drugnames[i]<<",";
        for (size_t j = 0; j < _vvd[i].size() ; j++) {
            fo << setprecision(18) << fixed<< _vvd[i][j] << ',';
        }
        fo << endl;
    }
    fo.close();
}

void mainTask() {
    for (size_t i = 0; i < nDrugs; i++) {
        for (size_t j = nDrugs; j < nDrugs + nTargets; j++) {
            for (size_t k = 0; k < nVertex; k++) {
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
    write_CVS(result);
    findnewInteraction(result,interactions);
    return 0;
}
