#include<bits/stdc++.h>
using namespace std;
bool isEigenVector(vector<vector<long double>> A, vector<long double> v, long double lambda) {
    int n = A.size();
    vector<long double> Av(n);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            Av[i] += A[i][j] * v[j];
        }
        if(abs(Av[i] - lambda * v[i]) > 1e-4) {
            return false;
        }
    }
    return true;
}

int main()
{
    vector<vector<long double>> A = {{0.68112222, -0.03900667, 1.26519111, 0.51345778},
                                     {-0.03900667, 0.18675067, -0.319568, -0.11719467},
                                     {1.26519111, -0.319568, 3.09242489, 1.28774489},
                                     {0.51345778, -0.11719467, 1.28774489, 0.57853156}};
    
    vector<long double> v = {-0.47971899,  0.07252408 , 0.1757674,0.85657211};
    long double lambda = 4.196675163197978;
    if(isEigenVector(A, v, lambda))
        cout << "v is an eigenvector of A with eigenvalue lambda." << endl;
    else
        cout << "v is not an eigenvector of A with eigenvalue lambda." << endl;
}