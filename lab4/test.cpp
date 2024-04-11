#include<bits/stdc++.h>
using namespace std;

vector<vector<long double>> generateRandomMatrix(int rows, int cols) {
    random_device rd;//随机数种子
    mt19937 gen(rd());//使用 rd 生成的随机数来初始化一个mt19937类型的随机数生成器 gen
    uniform_real_distribution<> dis(0.0, 1.0);//创建一个均匀实数分布 dis

    vector<std::vector<long double>> matrix(rows, std::vector<long double>(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            matrix[i][j] = dis(gen);

    return matrix;
}

pair<int, int> chooseMax(vector<vector<long double>>& A)
{
    long double max = 0;
    pair<int, int> maxPos;
    for(int i = 0; i < A.size(); i++)
        for(int j = 0; j < A.size(); j++)
            if(i != j && fabsl(A[i][j]) > max)
            {
                max = fabsl(A[i][j]);
                maxPos = make_pair(i, j);
            }
    return maxPos;
}

vector<vector<long double>> calAT(vector<vector<long double>>& A)
{
    vector<vector<long double>> AT(A.size(), vector<long double>(A.size()));
    for(int i = 0; i < A.size(); i++)
        for(int j = 0; j < A.size(); j++)
            AT[i][j] = A[j][i];
    return AT;
}

long double calculateElement(const vector<vector<long double>>& A, int i, int j, long double p, long double q, long double t, long double c, long double d) {
    if (i == p && j == p)
        return A[p][p] - t * A[p][q];
    else if (i == q && j == q)
        return A[q][q] + t * A[p][q];
    else if ((i == p && j == q) || (i == q && j == p))
        return 0;
    else if (i != q && i != p && (j == p || j == q))
        return (j == p ? c : d) * A[p][i] - (j == p ? d : (-c)) * A[q][i];
    else if ((i == p || i == q) && j != q && j != p)
        return (i == p ? c : d) * A[p][j] - (i == p ? d : (-c)) * A[q][j];
    else
        return A[i][j];
}

vector<vector<long double>> calQTAQ(vector<vector<long double>>& A)
{
    pair<int, int> maxPos = chooseMax(A);
    int row = maxPos.first;
    int col = maxPos.second;
    
    long double s = (A[col][col] - A[row][row]) / (2 * A[row][col]);
    long double t = 0;
    if(s == 0)
        t = 1;
    else if(abs(-s + sqrt(1 + s * s)) <= abs(-s - sqrt(1 + s * s)))
        t = -s + sqrt(1 + s * s);
    else
        t = -s - sqrt(1 + s * s);

    long double c = 1 / sqrt(1 + t * t);
    long double d = t * c;

    vector<vector<long double>> QTAQ(A.size(), vector<long double>(A.size()));
    for(int i = 0; i < A.size(); i++)
        for(int j = 0; j < A.size(); j++)
            QTAQ[i][j] = calculateElement(A, i, j, row, col, t, c, d);

    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
            cout << QTAQ[i][j] << " ";
        cout << endl;
    }
    cout << endl;
    return QTAQ;
}

int judgeEnd(vector<vector<long double>>& A)
{
    int i, j;
    for(i = 0; i < A.size(); i++)
        for(j = 0; j < A.size(); j++)
            if(i != j && fabsl(A[i][j]) >= 1e-7)
                return 0;
    if(i == A.size() && j == A.size()) 
        return 1;
}

vector<long double> calEigenValue(vector<vector<long double>> A)
{
    vector<long double> eigenValue(A.size());
    vector<vector<long double>> QTAQ= calQTAQ(A);
    int i, j;
    while(!judgeEnd(QTAQ))
        QTAQ = calQTAQ(QTAQ);
    for(i = 0; i < A.size(); i++)
        eigenValue[i] =QTAQ[i][i];
    return eigenValue;
}   

int main()
{
    vector<vector<long double>> A = vector<vector<long double>>(4, vector<long double>(4));
    A[0][0] = 1.31336, A[0][1] = 2.02598, A[0][2] = 1.08122, A[0][3] = 0;
    A[1][0] = 1.22657, A[1][1] = 1.94428, A[1][2] = 1.00043, A[1][3] = 0;
    A[2][0] = 0.829463, A[2][1] = 1.31122, A[2][2] = 0.841318, A[2][3] = 0;
    A[3][0] = 0.556764, A[3][1] = 0.895583, A[3][2] = 0.584515, A[3][3] = 0;
    
    /*vector<vector<long double>> A = {
        {1, 2, 1},
        {1, 2, 2},
        {3, 8, 3},
    };*/
    vector<long double> x =calEigenValue(A);
    cout << "EigenValue: ";
    for(int i = 0; i < A.size(); i++)
        cout << x[i] << " ";
}