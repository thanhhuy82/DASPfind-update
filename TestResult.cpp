#include<bits/stdc++.h>
#define EPSILON 0.000001
using namespace std;

bool AreSame(double a, double b)
{
    return fabs(a - b) < EPSILON;
}

int main() {
    double x, y;
    ifstream fi("1001_results.txt");
    ifstream fi1("output.txt");
    bool result=true;
    while (!fi.eof()) {
        fi >> x;
        fi1 >> y;
        if (!AreSame(x, y)) {
            result = false;
            cout << x << ' ' << y << endl;
        }
    }
    if (result) {
        cout << "True";
    } else {
        cout << "False";
    }
    return 0;
}
