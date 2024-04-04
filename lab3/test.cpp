#include<iostream>
#include<cmath>
#include <vector>
#include <iomanip>
using namespace std;

//计算向量的无穷范数，即找到向量中绝对值最大的元素
long double Norm(long double* x1, int n)
{
    long double norm = 0;
    for(int i = 0; i < n - 1; i++)
        if(fabsl(x1[i]) > norm)
            norm = fabsl(x1[i]);
    return norm;
}


long double* Column_Elimination(vector<std::vector<long double>>& A, long double* b, int n)
{
    long double* x = new long double[n]();
    vector<vector<long double>> matrix(n, vector<long double>(n+1, 0)); 
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            matrix[i][j] = A[i][j];
    for(int i = 0; i < n; i++)
        matrix[i][n] = b[i];
    for(int col = 0; col < n; col++)
    {
        long double maxnum = fabsl(matrix[col][col]);
        int maxrow = col;
        for(int row = col + 1; row < n; row++)
        {
            if(fabsl(matrix[row][col]) > maxnum)
            {
                maxnum = fabsl(matrix[row][col]);
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
    for(int row = n - 1; row >= 0; row--)
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

void Iteration(vector<vector<long double>>& A, long double* x, int n)
{
    long double* xtemp = new long double[n]();
    cout << "X(0):";
    for (int i = 0; i < n; i++)
    {
        cout << x[i] << " ";
        xtemp[i] = x[i];
    }
    cout << endl;
    long double norm1 = Norm(xtemp, n);
    cout << norm1 << endl;
    long double* x1 = nullptr;
    x1 = Column_Elimination(A, xtemp, n);
    cout << "X(1):";
    for (int i = 0; i < n; i++)
        cout << x1[i] << " ";
    cout << endl;
    long double norm2 = Norm(x1, n);
    cout << norm2 << endl;
    int k = 1;
    long double* y = new long double[n]();
    
    
    while (fabsl(norm1 - norm2) > 1e-5)
    {
        //long double norm = Norm(x1, n);
        for(int i = 0; i < n; i++)
            y[i] = x1[i] / norm2;
        cout << "Y(" << k << "):";
        for (int i = 0; i < n; i++)
            cout << y[i] << " ";
        cout << endl;
        if (xtemp != x1) 
            delete[] xtemp;
        //xtemp = x1;
        xtemp = Column_Elimination(A, y, n);
        cout << "X(" << k + 1 << "):";
        for (int i = 0; i < n; i++)
            cout << xtemp[i] << " ";
        cout << endl;
        norm1 = norm2;
        norm2 = Norm(xtemp, n);
        cout << norm2 << endl;
        k++;
        x1 = xtemp;
    }
    cout << "lambda=" << 1 / Norm(x1, n) << endl;
    cout << "特征向量为：";
    for(int i = 0; i < n; i++)
        cout << y[i] << " ";
    cout << endl << "迭代次数为：" << k << endl << endl;
    delete[] y;
    delete[] x1;
}

int main()
{
    int n = 5; 
    vector<std::vector<long double>> A1(n, std::vector<long double>(n));
    for(int i = 0; i < 5; i++)
        for(int j = 0; j < 5; j++)
            A1[i][j] = 1.0 / (9 - i - j);
    long double b1[5] = {1, 1, 1, 1, 1};
    Iteration(A1, b1, 5);
    n = 4; 
    vector<std::vector<long double>> A2(n, std::vector<long double>(n));
    A2[0][0] = 4; A2[0][1] = -1; A2[0][2] = 1; A2[0][3] = 3;
    A2[1][0] = 16; A2[1][1] = -2; A2[1][2] = -2; A2[1][3] = 5;
    A2[2][0] = 16; A2[2][1] = -3; A2[2][2] = -1; A2[2][3] = 7;
    A2[3][0] = 6; A2[3][1] = -4; A2[3][2] = 2; A2[3][3] = 9;
    long double b2[4] = {1, 1, 1, 1};
    Iteration(A2, b2, 4);
}
