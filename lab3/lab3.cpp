#include<bits/stdc++.h>
using namespace std;

//计算向量的无穷范数，即找到向量中绝对值最大的元素
long double Norm(const vector<long double>& x1)
{
    long double norm = 0;
    int n = x1.size();
    for(int i = 0; i < n; i++)
        if(fabsl(x1[i]) > norm)
            norm = fabsl(x1[i]);
    return norm;
}

//Doolittle分解解线性方程组
vector<long double> Doolittle(vector<vector<long double>>& A, const vector<long double>& b)
{
    int n = A.size();
    vector<long double> x(n, 0);
    vector<long double> temp(n, 0);
    vector<vector<long double>> matrix1(n, vector<long double>(n+1, 0)); 
    vector<vector<long double>> matrix2(n, vector<long double>(n+1, 0)); 
    vector<vector<long double>> Lmatrix(n, vector<long double>(n, 0));
    vector<vector<long double>> Umatrix(n, vector<long double>(n, 0));
    for(int i = 0; i < n; i++)
    {
        Lmatrix[i][i] = 1;
        Umatrix[0][i] = A[0][i];
        Lmatrix[i][0] = A[i][0] / Umatrix[0][0];
    }
    for(int i = 1; i < n; i++)
    {
        for(int j = i; j < n; j++)
        {
            Umatrix[i][j] = A[i][j];
            for(int k = 0; k < i; k++)
                Umatrix[i][j] -= Lmatrix[i][k] * Umatrix[k][j];
        }
        for(int j = i + 1; j < n; j++)
        {
            Lmatrix[j][i] = A[j][i];
            for(int k = 0; k < i; k++)
                Lmatrix[j][i] -= Lmatrix[j][k] * Umatrix[k][i];
            Lmatrix[j][i] /= Umatrix[i][i];
        }
    }

    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            matrix1[i][j] = Lmatrix[i][j];
    for(int i = 0; i < n; i++)
        matrix1[i][n] = b[i];
    for(int row = 0; row < n; row++)
    {
        for(int col = 0; col < row; col++)
        {
            matrix1[row][n] -= matrix1[col][n] * matrix1[row][col] / matrix1[col][col];
            matrix1[row][col] = 0;
        }
        matrix1[row][n] /= matrix1[row][row];
        matrix1[row][row] = 1;
    }
    for(int i = 0; i < n; i++)
        temp[i] = matrix1[i][n];
    
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            matrix2[i][j] = Umatrix[i][j];
    for(int i = 0; i < n; i++)
        matrix2[i][n] = temp[i];
    for(int row = n - 1; row >= 0; row--)
    {
        for(int col = row + 1; col < n; col++)
        {
            matrix2[row][n] -= matrix2[col][n] * matrix2[row][col] / matrix2[col][col];
            matrix2[row][col] = 0;
        }
        matrix2[row][n] /= matrix2[row][row];
        matrix2[row][row] = 1;
    }
    for(int i = 0; i < n; i++)
        x[i] = matrix2[i][n];
    return x;
}

//迭代求解特征值和特征向量
void Iteration(vector<vector<long double>>& A, vector<long double> x)
{
    int n = A.size();
    vector<long double> xtemp = x;
    cout << "X(0):";
    for (int i = 0; i < n; i++)
        cout << x[i] << " ";
    cout << endl;
    long double norm1 = Norm(xtemp);
    cout << norm1 << endl;
    vector<long double> x1 = Doolittle(A, xtemp);
    cout << "X(1):";
    for (int i = 0; i < n; i++)
        cout << x1[i] << " ";
    cout << endl;
    long double norm2 = Norm(x1);
    cout << norm2 << endl;
    int k = 1;
    vector<long double> y(n, 0);
    
    while (fabsl(norm1 - norm2) > 1e-5)
    {
        for(int i = 0; i < n; i++)
            y[i] = x1[i] / norm2;
        cout << "Y(" << k << "):";
        for (int i = 0; i < n; i++)
            cout << y[i] << " ";
        cout << endl;
        xtemp = Doolittle(A, y);
        cout << "X(" << k + 1 << "):";
        for (int i = 0; i < n; i++)
            cout << xtemp[i] << " ";
        cout << endl;
        norm1 = norm2;
        norm2 = Norm(xtemp);
        cout << norm2 << endl;
        k++;
        x1 = xtemp;
    }
    cout << "lambda=" << 1 / Norm(x1) << endl;
    cout << "特征向量为：";
    for(int i = 0; i < n; i++)
        cout << y[i] << " ";
    cout << endl << "迭代次数为：" << k << endl << endl;
}

int main()
{
    int n = 5;
    vector<vector<long double>> A1(n, vector<long double>(n));
    for(int i = 0; i < 5; i++)
        for(int j = 0; j < 5; j++)
            A1[i][j] = 1.0 / (9 - i - j);
    vector<long double> b1 = {1, 1, 1, 1, 1};
    Iteration(A1, b1);

    n = 4;
    vector<vector<long double>> A2(n, vector<long double>(n));
    A2[0][0] = 4; A2[0][1] = -1; A2[0][2] = 1; A2[0][3] = 3;
    A2[1][0] = 16; A2[1][1] = -2; A2[1][2] = -2; A2[1][3] = 5;
    A2[2][0] = 16; A2[2][1] = -3; A2[2][2] = -1; A2[2][3] = 7;
    A2[3][0] = 6; A2[3][1] = -4; A2[3][2] = 2; A2[3][3] = 9;
    vector<long double> b2 = {1, 1, 1, 1};
    Iteration(A2, b2);
}