#include<bits/stdc++.h>
using namespace std;

long double a = 1.0 / 2, n = 100, h = 1.0 / n;
vector<long double> x(99), b(99);
vector<vector<long double>> A(99, vector<long double>(99));

//Inv矩阵与x向量相乘
vector<long double> Multiplie(vector<vector<long double>>& Inv, vector<long double>& x)
{
    vector<long double> res(99);
    for(int i = 0; i < n - 1; i++)
        for(int j = 0; j < n - 1; j++)
            res[i] += Inv[i][j] * x[j];
    return res;
}

// 计算两个向量的差向量的无穷范数
long double Norm(vector<long double>& x1, vector<long double>& x2)
{
    long double norm = 0;
    for(int i = 0; i < n - 1; i++)
        if(abs(x1[i] - x2[i]) > norm)
            norm = abs(x1[i] - x2[i]);
    return norm;
}

// Gauss-Seidel迭代法求解线性方程组
vector<long double> Gauss_Seidel(vector<vector<long double>>& S, vector<vector<long double>>& Inv)
{
    vector<long double> x1(99), tempx(99);
    for(int i = 0; i < n - 1; i++)
        tempx[i] = x[i];
    vector<long double> temp1 = Multiplie(S, tempx);
    vector<long double> temp2 = Multiplie(Inv, b);
    for(int i = 0; i < n - 1; i++)    
        x1[i] = temp1[i] + temp2[i];
    while(Norm(tempx, x1) > 1e-6)
    {
        for(int i = 0; i < n - 1; i++) 
            tempx[i] = x1[i];
        temp1 = Multiplie(S, tempx);
        temp2 = Multiplie(Inv, b);
        for(int i = 0; i < n - 1; i++)    
            x1[i] = temp1[i] + temp2[i];
    }
    return x1;
}

// 计算精确解
vector<long double> PreciseSol(vector<long double>& x , long double epsilon)
{
    vector<long double> y(99);
    for(int i = 0; i < n - 1; i++)
        y[i] = (1 - a) * (1 - exp(-x[i] / epsilon)) / (1 - exp(-1 / epsilon)) + a * x[i];
    return y;
}

// 计算误差并输出
void Calculate(long double epsilon)
{
    for(int i = 0; i < n - 1; i++)
        b[i] = a * h * h;
    b[98] -= epsilon + h;
    vector<long double> Presol = PreciseSol(x , epsilon);
    for(int i = 0; i < n - 2; i++)
    {
        A[i][i + 1] = epsilon + h;
        A[i + 1][i] = epsilon;
    }
    for(int i = 0; i < n - 1; i++)
        A[i][i] = -(2 * epsilon + h);
    vector<vector<long double>> S(99, vector<long double>(99));
    vector<vector<long double>> Inv(99, vector<long double>(99));
    fill(&S[0][0], &S[0][0] + sizeof(S) / sizeof(S[0][0]), 0);
    fill(&Inv[0][0], &Inv[0][0] + sizeof(Inv) / sizeof(Inv[0][0]), 0);
    long double init1 = (epsilon + h) / (2 * epsilon + h);
    long double temp = epsilon / (2 * epsilon + h);
    for(int i = 0; i < n - 1; i++)//初始化S=-(D+L)^(-1)U
    {
        for(int j = i , k = 1; j < n -1 && k < n - 1; j++, k++)
        {
            S[j][k] = init1;
        }
        init1 *= temp; 
    }
    long double init2 = -1 / (2 * epsilon + h);
    for(int i = 0; i < n - 1; i++)//初始化(D+L)^(-1)
    {
        for(int j = i, k = 0; j < n - 1 && k < n - 1; j++, k++)
        {
            Inv[j][k] = init2;
        }
        init2 *= temp; 
    }
    vector<long double> Seidel = Gauss_Seidel(S, Inv);
    long double errorSeidel = 0;
    for(int i = 0; i < n - 1; i++)
    	errorSeidel += abs(Seidel[i] - Presol[i]) / 99;
    cout << "epsilon=" << epsilon << "时，Gauss_Seidel迭代法的解为：" ;
    for(int i = 0; i < n - 1; i++)
        cout << Seidel[i] << " ";
    cout << endl;
    cout << "epsilon=" << epsilon << "时，Gauss_Seidel迭代法的相对误差为"  << errorSeidel << endl;
}

int main()
{
    for(int i = 0; i < n - 1; i++)
        x[i] = (i + 1) * h;
    Calculate(1);
    Calculate(0.1);
    Calculate(0.01);
    Calculate(0.0001);
}