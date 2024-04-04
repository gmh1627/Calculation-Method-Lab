#include<iostream>
#include<cmath>
#include <iomanip>
using namespace std;
#define n 3
long double x[n], A[n][n], b[n];

long double* Multiplie(long double (*Inv)[n], long double* x)
{
    long double* res = new long double[n]();
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
            res[i] += Inv[i][j] * x[j];
    }
    return res;
}

long double Norm(long double* x1, long double* x2)
{
    long double norm = 0;
    for(int i = 0; i < n; i++)
        if(abs(x1[i] - x2[i]) > norm)
            norm = abs(x1[i] - x2[i]);
    return norm;
}

long double* Gauss_Seidel(long double (*S)[n], long double (*Inv)[n])
{
    long double* x1 = new long double[n]();
    long double tempx[n];
    for(int i = 0; i < n; i++)
        tempx[i] = x[i];
    long double *temp1 = Multiplie(S, tempx);
    long double *temp2 = Multiplie(Inv, b);
    for(int i = 0; i < n; i++)    
        x1[i] = temp1[i] + temp2[i];
    while(Norm(tempx, x1) > 1e-6)
    {
        for(int i = 0; i < n; i++)
            cout << tempx[i] << " ";
        for(int i = 0; i < n; i++) 
            tempx[i] = x1[i];
        long double *temp1 = Multiplie(S, tempx);
        long double *temp2 = Multiplie(Inv, b);
        for(int i = 0; i < n; i++)    
            x1[i] = temp1[i] + temp2[i];
        cout << endl;
    }
    return x1;
}

int main()
{
    long double S[n][n], Inv[n][n];
    fill(&S[0][0], &S[0][0] + sizeof(S) / sizeof(S[0][0]), 0);
    fill(&Inv[0][0], &Inv[0][0] + sizeof(Inv) / sizeof(Inv[0][0]), 0);
    x[0] = 0, x[1] = 0, x[2] = 0;
    b[0] = 0, b[1] = -21, b[2] = -20;
    S[0][1] = 0.2, S[0][2] = 0.1;
    S[1][1] = 0.04, S[1][2] = 0.12;
    S[2][1] = 0.056, S[2][2] = 0.0428;
    Inv[0][0] = 0.1;
    Inv[1][0] = 0.02, Inv[1][1] = 0.1;
    Inv[2][0] = 0.028, Inv[2][1] = 0.04, Inv[2][2] = 0.2;
    long double* Seidel = Gauss_Seidel(S, Inv);
    for(int i = 0; i < n; i++)
        cout << Seidel[i] << " ";
}