#include<bits/stdc++.h>
using namespace std;

long double a = 1.0 / 2, n = 100, h = 1.0 / n;
vector<long double> x(99), b(99);
vector<vector<long double>> A(99, vector<long double>(99));

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
    vector<long double> Column = Column_Elimination(A);
    long double errorColumn = 0;
    for(int i = 0; i < n - 1; i++)
        errorColumn += abs(Column[i] - Presol[i]) / 99;
    cout << "epsilon=" << epsilon << "时，列主元消元法的解为：" ;
    for(int i = 0; i < n - 1; i++)
        cout << Column[i] << " ";
    cout << endl;
    cout << "epsilon=" << epsilon << "时，列主元消元法的相对误差为"  << errorColumn << endl;
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