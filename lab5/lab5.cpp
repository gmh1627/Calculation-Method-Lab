#include<bits/stdc++.h>
using namespace std;

// 列主元消去法求解线性方程组
vector<long double> Column_Elimination(vector<vector<long double>>& A)
{
    vector<long double> x(99);
    vector<vector<long double>> matrix(99, vector<long double>(100));
    fill(matrix.begin(), matrix.end(), vector<long double>(100, 0));
    for(int i = 0; i < n - 1; i++)
        for(int j = 0; j < n - 1; j++)
            matrix[i][j] = A[i][j];
    for(int i = 0; i < n - 1; i++)
        matrix[i][99] = b[i];
    for(int col = 0; col < n - 1; col++)//找到列主元
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
        for(int row = col + 1; row < n - 1; row++)
        {
            long double res = matrix[row][col] / matrix[col][col];
            for(int loc = col; loc <= n - 1; loc++)
                matrix[row][loc] -= matrix[col][loc] * res; 
        }
    }
    for(int row = n - 2; row >= 0; row--)//回代
    {
        for(int col = row + 1; col < n - 1; col++)
        {
            matrix[row][99] -= matrix[col][99] * matrix[row][col] / matrix[col][col];
            matrix[row][col] = 0;
        }
        matrix[row][99] /= matrix[row][row];
        matrix[row][row] = 1;
    }
    for(int i = 0; i < n - 1; i++)
        x[i] = matrix[i][99];
    return x;
}

int main()
{
    
}