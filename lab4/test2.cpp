#include<bits/stdc++.h>
using namespace std;

vector<vector<long double>> Column_Elimination(vector<vector<long double>> A)
{
    int n = A.size();
    vector<vector<long double>> Temp(n, vector<long double>(n));
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            Temp[i][j] = A[i][j];
    for(int col = 0; col < n; col++)
    {
        long double maxnum = abs(Temp[col][col]);
        int maxrow = col;
        for(int row = col + 1; row < n; row++)
        {
            if(abs(Temp[row][col]) > maxnum)
            {
                maxnum = abs(Temp[row][col]);
                maxrow = row;
            }
        }
        swap(Temp[col], Temp[maxrow]);
        for(int row = col + 1; row < n; row++)
        {
            long double res = Temp[row][col] / Temp[col][col];
            for(int loc = col; loc < n; loc++)
                Temp[row][loc] -= Temp[col][loc] * res; 
        }
    }
    return Temp;
}
/*vector<vector<long double>> A = vector<vector<long double>>(3, vector<long double>(4));
    A[0][0] = 3, A[0][1] = 1, A[0][2] = 2, A[0][3] = 1;
    A[1][0] = 1, A[1][1] = 3, A[1][2] = 4, A[1][3] = 1;
    A[2][0] = 2, A[2][1] = 4, A[2][2] = 6, A[2][3] = 1;
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 4; j++)
            cout << A[i][j] << " ";
        cout << endl;
    }*/
//
vector<long double> SolveUpperTriangle(vector<vector<long double>> A, vector<long double> b)
{
    int n = A.size();
    vector<long double> x(n);
    vector<vector<long double>> Temp(n, vector<long double>(n+1));
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            Temp[i][j] = A[i][j];
    for(int i = 0; i < n; i++)
        Temp[i][n] = b[i];
    for(int row = n-1; row >= 0; row--)
    {
        for(int col = row + 1; col < n; col++)
        {
            Temp[row][n] -= Temp[col][n] * Temp[row][col] / Temp[col][col];
            Temp[row][col] = 0;
        }
        Temp[row][n] /= Temp[row][row];
        Temp[row][row] = 1;
    }
    for(int i = 0; i < n; i++)
        x[i] = Temp[i][n];
    return x;
}

//
vector<vector<long double>> solve(vector<vector<long double>> A, int cnt)
{
    int n = A.size();
    vector<vector<long double>> x(cnt, vector<long double>(n));
    vector<vector<long double>> Temp(n-cnt, vector<long double>(n-cnt));
    vector<long double> Tempb(n-cnt);
    for(int i = 0; i < cnt; i++)
    {
        for(int j = n - 1; j >= n - cnt; j--)
        {
            if(j >= n - i)
                x[i][j] = 0;
            else
                x[i][j] = 1;
        }
    }

    /*for(int i = 0; i < cnt; i++)
        for(int j = 0; j < n; j++)
            cout << x[i][j] << " ";*/
    for(int i = 0; i < n - cnt; i++)
        for(int j = 0; j < n - cnt; j++)
            Temp[i][j] = A[i][j];
    for(int i = 0; i < cnt; i++)
    {
        for(int j = n - cnt - 1; j >=  0; j--)
        {
            Tempb[j] = 0;
            for(int k = 0; k < cnt; k++)
                Tempb[j] -= A[j][n- cnt + k] * x[i][n- cnt + k];
        }
        vector<long double> res = SolveUpperTriangle(Temp, Tempb);
        for(int j = 0; j < n - cnt; j++)
            x[i][j] = res[j];
    }
    return x;
}

int main()
{
    vector<vector<long double>> A = {{1, 2, 3}, {1, 2, 3}, {7, 8, 9}};
    int n = A.size();
    vector<vector<long double>> B = Column_Elimination(A);
    int cnt = 0;//记录消元后全为0的行数
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            if(fabsl(B[i][j]) > 1e-7)
                break;
            else if(j == n - 1)
                cnt++;
        }
    }
    cout << cnt << endl;
    vector<vector<long double>> x = solve(B, cnt);
    for(int i = 0; i < cnt; i++)
    {
        for(int j = 0; j < n; j++)
            cout << x[i][j] << " ";
        cout << endl;
    }
    cout << endl;
    return 0;
}