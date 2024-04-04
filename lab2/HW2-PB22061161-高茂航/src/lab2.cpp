#include<iostream>
#include<cmath>
#include <iomanip>
using namespace std;

long double a = 1.0 / 2, n = 100, h = 1.0 / n;
long double x[99], A[99][99], b[99];

// 交换两个数组的元素
void swap(long double* x, long double* y)
{
    for(int i = 0; i <= n - 1; i++)
    {
        long double temp = x[i];
        x[i] = y[i];
        y[i] = temp;
    }
    return ;
}

// 矩阵和向量的乘法，用于求Gauss-Seidel迭代的f向量
long double* Multiplie(long double (*Inv)[99], long double* x)
{
    long double* res = new long double[99]();
    for(int i = 0; i < n - 1; i++)
        for(int j = 0; j < n - 1; j++)
            res[i] += Inv[i][j] * x[j];
    return res;
}

// 列主元消去法求解线性方程组
long double* Column_Elimination(long double (*A)[99])
{
    long double* x = new long double[99]();
    long double (*matrix)[100] = new long double[99][100];//增广矩阵
    fill_n(&matrix[0][0], 99*100, 0);
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

// 计算两个向量的差向量的无穷范数
long double Norm(long double* x1, long double* x2)
{
    long double norm = 0;
    for(int i = 0; i < n - 1; i++)
        if(abs(x1[i] - x2[i]) > norm)
            norm = abs(x1[i] - x2[i]);
    return norm;
}

// Gauss-Seidel迭代法求解线性方程组
long double* Gauss_Seidel(long double (*S)[99], long double (*Inv)[99])
{
    long double* x1 = new long double[99]();
    long double tempx[99];
    for(int i = 0; i < n - 1; i++)
        tempx[i] = x[i];
    long double *temp1 = Multiplie(S, tempx);
    long double *temp2 = Multiplie(Inv, b);
    for(int i = 0; i < n - 1; i++)    
        x1[i] = temp1[i] + temp2[i];
    while(Norm(tempx, x1) > 1e-6)
    {
        for(int i = 0; i < n - 1; i++) 
            tempx[i] = x1[i];
        long double *temp1 = Multiplie(S, tempx);
        long double *temp2 = Multiplie(Inv, b);
        for(int i = 0; i < n - 1; i++)    
            x1[i] = temp1[i] + temp2[i];
    }
    return x1;
}

// 计算精确解
long double* PreciseSol(long double* x , long double epsilon)
{
    long double* y = new long double[99];
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
    long double* Presol = PreciseSol(x , epsilon);
    for(int i = 0; i < n - 2; i++)
    {
        A[i][i + 1] = epsilon + h;
        A[i + 1][i] = epsilon;
    }
    for(int i = 0; i < n - 1; i++)
        A[i][i] = -(2 * epsilon + h);
    long double S[99][99], Inv[99][99];
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
    long double* Column = Column_Elimination(A);
    long double* Seidel = Gauss_Seidel(S, Inv);
    long double errorColumn = 0, errorSeidel = 0;
    for(int i = 0; i < n - 1; i++){
        errorColumn += abs(Column[i] - Presol[i]) / 99;
        errorSeidel += abs(Seidel[i] - Presol[i]) / 99;
    }
    cout << "epsilon=" << epsilon << "时，列主元消元法的解为：" ;
    for(int i = 0; i < n - 1; i++)
        cout << Column[i] << " ";
    cout << endl;
    cout << "epsilon=" << epsilon << "时，Gauss_Seidel迭代法的解为：" ;
    for(int i = 0; i < n - 1; i++)
        cout << Seidel[i] << " ";
    cout << endl;
    cout << "epsilon=" << epsilon << "时，列主元消元法的相对误差为"  << errorColumn <<"，Gauss_Seidel迭代法的相对误差为"  << errorSeidel << endl;
    delete[] Presol;
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