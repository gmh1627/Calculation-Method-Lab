#include<bits/stdc++.h>
using namespace std;

//计算向量的无穷范数，即找到向量中绝对值最大的元素
/*long double Norm(long double* x1, int n)
{
    long double norm = 0;
    for(int i = 0; i < n - 1; i++)
        if(fabsl(x1[i]) > norm)
            norm = fabsl(x1[i]);
    return norm;
}*/
#include<bits/stdc++.h>
using namespace std;

// Function prototypes
vector<vector<long double>> generateRandomMatrix(int rows, int cols);
pair<int, int> chooseMax(vector<vector<long double>>& A);
vector<vector<long double>> calAT(vector<vector<long double>> A);
vector<vector<long double>> calAAT(vector<vector<long double>>& A);
long double calculateElement(const vector<vector<long double>>& A, int i, int j, long double p, long double q, long double t, long double c, long double d);
vector<vector<long double>> calQTAQ(vector<vector<long double>> A);
int judgeEnd(vector<vector<long double>> A);
vector<long double> calEigenValue(vector<vector<long double>>& A);
vector<vector<long double>> Column_Elimination(vector<vector<long double>> A);
vector<long double> SolveUpperTriangle(vector<vector<long double>> A, vector<long double> b);
vector<vector<long double>> solve(vector<vector<long double>> A, int cnt);
vector<vector<long double>> calEigenVector(vector<vector<long double>> A, vector<long double> eigenValue);
vector<vector<long double>> calSigma(vector<long double> x);

int main()
{
    vector<vector<long double>> A = generateRandomMatrix(4, 3);
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
    vector<vector<long double>> AT = calAT(A);
    vector<vector<long double>> AAT = calAAT(A);
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
            cout << AAT[i][j] << " ";
        cout << endl;
    }
    vector<long double> x =calEigenValue(AAT);
    cout << "EigenValue: ";
    for(int i = 0; i < x.size(); i++)
        cout << x[i] << " ";
    vector<vector<long double>> Sigma1 = calSigma(x);
    vector<long double> x1;
    x1.reserve(x.size());
    unique_copy(x.begin(), x.end(), back_inserter(x1));
    vector<vector<long double>> eigenVector1 = calEigenVector(AAT, x1);
    cout << endl << "EigenVector: " << endl;
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
            cout << eigenVector1[i][j] << " ";
        cout << endl;
    }
    cout << endl;
/*
    vector<vector<long double>> ATA = calAAT(AT);
    vector<long double> y =calEigenValue(ATA);
    cout << "EigenValue: ";
    for(int i = 0; i < 3; i++)
        cout << y[i] << " ";
    vector<long double> y1;
    x1.reserve(y.size());
    unique_copy(y.begin(), y.end(), back_inserter(y1));
    vector<vector<long double>> eigenVector2 = calEigenVector(ATA, y1);
    cout << endl << "EigenVector: " << endl;
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
            cout << eigenVector2[i][j] << " ";
        cout << endl;
    }*/
    return 0;
}

//Generate a random matrix with given rows and columns
vector<vector<long double>> generateRandomMatrix(int rows, int cols) {
    random_device rd;//随机数种子
    mt19937 gen(rd());//使用 rd 生成的随机数来初始化一个mt19937类型的随机数生成器 gen
    uniform_real_distribution<> dis(0.0, 1.0);//创建一个均匀实数分布 dis

    vector<vector<long double>> matrix(rows, vector<long double>(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            matrix[i][j] = dis(gen);

    return matrix;
}

// Find the position of the maximum absolute value on the non diagonal of the matrix
pair<int, int> chooseMax(vector<vector<long double>>& A)
{
    long double max = 0;
    pair<int, int> maxPos;
    int n = A.size();
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            if(i != j && fabsl(A[i][j]) > max)
            {
                max = fabsl(A[i][j]);
                maxPos = make_pair(i, j);
            }
    return maxPos;
}

//Calculate the transpose of the matrix
vector<vector<long double>> calAT(vector<vector<long double>> A)
{
    int n = A.size();
    vector<vector<long double>> AT(n, vector<long double>(n));
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            AT[i][j] = A[j][i];
    return AT;
}

//Calculate the product of A and A^T
vector<vector<long double>> calAAT(vector<vector<long double>>& A)
{
    int n = A.size();
    vector<vector<long double>> AAT(n, vector<long double>(n));
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            for(int k = 0; k < n; k++)
                AAT[i][j] += A[i][k] * A[j][k];
    return AAT;
}

//Calculate the element of the matrix Q^T * A * Q
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

//Calculate the product of Q^T * A * Q
vector<vector<long double>> calQTAQ(vector<vector<long double>> A)
{
    int n = A.size();
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

    vector<vector<long double>> QTAQ(n, vector<long double>(n));
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            QTAQ[i][j] = calculateElement(A, i, j, row, col, t, c, d);

    /*for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
            cout << QTAQ[i][j] << " ";
        cout << endl;
    }
    cout << endl;*/
    return QTAQ;
}

//Judge whether the iteration is end
int judgeEnd(vector<vector<long double>> A)
{
    int i, j;
    int n = A.size();
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            if(i != j && fabsl(A[i][j]) >= 1e-7)
                return 0;
    if(i == n && j == n) 
        return 1;
}

//Calculate the eigenvalue of the matrix A with Jacobi method
vector<long double> calEigenValue(vector<vector<long double>>& A)
{
    int n = A.size();
    vector<long double> eigenValue(n);
    vector<vector<long double>> QTAQ= calQTAQ(A);
    int i, j;
    while(!judgeEnd(QTAQ))
        QTAQ = calQTAQ(QTAQ);
    for(i = 0; i < n; i++)
        eigenValue[i] =QTAQ[i][i];
    return eigenValue;
}

//Column principal element elimination method
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

//Calculate the eigenvector of the matrix A
vector<vector<long double>> calEigenVector(vector<vector<long double>> A, vector<long double> eigenValue)
{
    int n = A.size();
    int num = 0;
    vector<vector<long double>> x(n, vector<long double>(n));
    vector<vector<long double>> tempMartix(n, vector<long double>(n));
    vector<vector<long double>> eigenVector(n, vector<long double>(n));
    for(int k = 0; k < n; k++)
    {
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                i == j ? tempMartix[i][j] = A[i][j] - eigenValue[k] : tempMartix[i][j] = A[i][j];
        /*cout << "tempMartix: " << endl;
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
                cout << tempMartix[i][j] << " ";
            cout << endl;
        }*/
        
        vector<vector<long double>> B = Column_Elimination(tempMartix);
        /*cout << "B: " << endl;
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
                cout << B[i][j] << " ";
            cout << endl;
        }*/
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
        vector<vector<long double>> result = solve(B, cnt);
        for(int i = 0; i < cnt; i++)
            copy(result[i].begin(), result[i].end(), x[num + i].begin());
        num += cnt;
    }
    return x;
}

//Calculate the matrix Sigma with the eigenvalue
vector<vector<long double>> calSigma(vector<long double> x)
{
    int n = x.size();
    vector<vector<long double>> Sigma(n, vector<long double>(n));
    for(int i = 0; i < n; i++)
        Sigma[i][i] = sqrt(x[i]);
    return Sigma;
}