#include<bits/stdc++.h>
using namespace std;

// 列主元消去法求解线性方程组
vector<long double> Column_Elimination(vector<vector<long double>> A, vector<long double> b);

int main()
{
    ifstream file("point.txt");
    if (!file) {
        cerr << "Unable to open file point.txt";
        return 1;
    }

    vector<long double> x, y, lambda, d, h, miu;
    long double a, b;
    while (file >> a >> b)
    {
        x.push_back(a);
        y.push_back(b);
        lambda.push_back(0);
        d.push_back(0);
        h.push_back(0);
        miu.push_back(0);
    }
    int len = x.size();
    int n = len - 1;
    for(int i = 0; i < n; i++)
        h[i] = x[i + 1] - x[i];
    for(int i = 1; i < n; i++)
    {
        lambda[i] = h[i] / (h[i] + h[i - 1]);
        miu[i] = 1 - lambda[i];
        d[i] = 6 / (h[i] + h[i - 1]) * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
    }
    vector<vector<long double>> A(n - 1, vector<long double>(n - 1, 0));
    for(int i = 0; i < n - 1; i++)
    {
        A[i][i] = 2;
        if(i != 0)
            A[i][i - 1] = miu[i + 1];
        if(i != n - 2)
            A[i][i + 1] = lambda[i + 1];
    }
    vector<long double> B(n - 1, 0);
    for(int i = 0; i < n - 1; i++)
        B[i] = d[i + 1];
    vector<long double> M = Column_Elimination(A, B);
    for(int i = 0; i < len; i++)
        cout << M[i] << " ";
    file.close();
    return 0;
}

// 列主元消去法求解线性方程组
vector<long double> Column_Elimination(vector<vector<long double>> A, vector<long double> b)
{
    int n = A.size();
    vector<long double> x(n + 2, 0);
    vector<vector<long double>> matrix(n, vector<long double>(n + 1, 0));

    for(int i = 0; i < n; i++)
    {
        matrix[i][n] = b[i];
        for(int j = 0; j < n; j++)
            matrix[i][j] = A[i][j];
    }
        
    for(int col = 0; col < n; col++)//找到列主元
    {
        long double maxnum = abs(matrix[col][col]);
        int maxrow = col;
        for(int row = col + 1; row < n - 1; row++)
        {
            if(abs(matrix[row][col]) > maxnum)
            {
                maxnum = abs(matrix[row][col]);
                maxrow = row;
            }
        }
        swap(matrix[col], matrix[maxrow]);
        for(int row = col + 1; row < n; row++)
        {
            long double res = matrix[row][col] / matrix[col][col];
            for(int loc = col; loc <= n - 1; loc++)
                matrix[row][loc] -= matrix[col][loc] * res; 
        }
    }
    for(int row = n - 1; row >= 0; row--)//回代
    {
        for(int col = row + 1; col < n; col++)
        {
            matrix[row][n] -= matrix[col][n] * matrix[row][col] / matrix[col][col];
            matrix[row][col] = 0;
        }
        matrix[row][n] /= matrix[row][row];
        matrix[row][row] = 1;
    }
    for(int i = 0; i < n; i++)
        x[i+1] = matrix[i][n];
    return x;
}