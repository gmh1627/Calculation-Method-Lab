#include<iostream>
#include<cmath>
#include <iomanip>
#define n 4
using namespace std;

long double a = 1.0 / 2, h = 1.0 / n;
long double x[n], A[n][n], b[n];

void swap(long double* x, long double* y)
{
    for(int i = 0; i <= n; i++)
    {
        long double temp = x[i];
        x[i] = y[i];
        y[i] = temp;
    }
    return ;
}

long double* Column_Elimination(long double (*A)[n])
{
    long double* x = new long double[n]();
    long double (*matrix)[n + 1] = new long double[n][n + 1];//增广矩阵
    fill_n(&matrix[0][0], n * (n + 1), 0);
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            matrix[i][j] = A[i][j];
    for(int i = 0; i < n; i++)
        matrix[i][n] = b[i];
    for(int col = 0; col < n; col++)
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
        //cout << matrix[1][0] <<  " "  << matrix[1][1] << " " << matrix[1][2] << endl;
        for(int row = col + 1; row < n; row++)
        {
            long double res = matrix[row][col] / matrix[col][col];
            for(int loc = col; loc <= n; loc++)
                matrix[row][loc] -= matrix[col][loc] * res; 
        }
        //cout << matrix[1][0] <<  " "  << matrix[1][1] << " " << matrix[1][2] << endl;
    }
    for(int row = n - 1; row >= 0; row--)
    {
        //matrix[row][3] /= matrix[row][row];
        //matrix[row][row] = 1;
        for(int col = row + 1; col < n; col++)
        {
            matrix[row][n] -= matrix[col][n] * matrix[row][col] / matrix[col][col];
            matrix[row][col] = 0;
        }
        matrix[row][n] /= matrix[row][row];
        matrix[row][row] = 1;
        //cout << matrix[1][0] <<  " "  << matrix[1][1] << " " << matrix[1][2] << endl;
    }
    for(int i = 0; i < n; i++)
        x[i] = matrix[i][n];
    return x;
}

int main()
{
    /*b[0] = 37.27;
    b[1] = -87.15;
    A[0][0] = 4.453, A[0][1] = -7.26, A[1][0] = -0.002, A[1][1] = -87.13; */
    A[0][0] = 7.2, A[0][1] = 2.3, A[0][2] = -4.4, A[0][3] = 0.5;
    A[1][0] = 1.3, A[1][1] = 6.3, A[1][2] = -3.5, A[1][3] = 2.8;
    A[2][0] = 5.6, A[2][1] = 0.9, A[2][2] = 8.1, A[2][3] = -1.3;
    A[3][0] = 1.5, A[3][1] = 0.4, A[3][2] = 3.7, A[3][3] = 5.9;
    b[0] = 15.1, b[1] = 1.8, b[2] = 16.6, b[3] = 36.9;
    long double* x = Column_Elimination(A);
    for(int i = 0; i < n; i++)
        cout << x[i] << endl;
}