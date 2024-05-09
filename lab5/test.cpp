#include<bits/stdc++.h>
using namespace std;

// 列主元消去法求解线性方程组
vector<long double> Column_Elimination(vector<vector<long double>> A, vector<long double> b);

int main()
{
    vector<vector<long double>> A(3, vector<long double>(3, 0));
    A[0][0] = 2, A[0][1] = 4, A[0][2] = -3;
    A[1][0] = 4, A[1][1] = -3, A[1][2] = 1;
    A[2][0] = 1, A[2][1] = 3, A[2][2] = -2;
    vector<long double> b(3, 1);
    vector<long double> X = Column_Elimination(A, b);
    for(int i = 0; i < 3; i++)
        cout << X[i] << " ";
    return 0;
}

// 列主元消去法求解线性方程组
vector<long double> Column_Elimination(vector<vector<long double>> A, vector<long double> b)
{
    int n = A.size();
    vector<long double> x(n, 0);
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
        for(int row = col + 1; row < n; row++)
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
            for(int loc = col; loc <= n; loc++)
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
        x[i] = matrix[i][n];
    return x;
}