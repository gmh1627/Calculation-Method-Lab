#include<bits/stdc++.h>
using namespace std;
bool isEigenVector(vector<vector<long double>> A, vector<long double> v, long double lambda) {
    int n = A.size();
    vector<long double> Av(n);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            Av[i] += A[i][j] * v[j];
        }
        if(abs(Av[i] - lambda * v[i]) > 1e-7) {
            return false;
        }
    }
    return true;
}

int main()
{
    vector<vector<long double>> A = {{1.71851, 1.66005, 1.14971, 0.792285}, {1.66005, 1.65417, 1.09037, 0.759718}, {1.14971, 1.09037, 0.906686, 0.636949}, {0.792285, 0.759718, 0.636949, 0.450517}};
    vector<long double> v = {2.05987, 2.00456, 1.43646, 1};
    long double lambda = 4.52037;

    if(isEigenVector(A, v, lambda))
        cout << "v is an eigenvector of A with eigenvalue lambda." << endl;
    else
        cout << "v is not an eigenvector of A with eigenvalue lambda." << endl;
}