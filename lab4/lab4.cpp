#include<bits/stdc++.h>
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

//交换两个向量
void swap(vector<long double>& x, vector<long double>& y)
{
    vector<long double> temp = x;
    x = y;
    y = temp;
}

//Generate a random matrix with given rows and columns
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

// Find the position of the maximum absolute value on the non diagonal of the matrix
pair<int, int> chooseMax(vector<vector<long double>>& A)
{
    long double max = 0;
    pair<int, int> maxPos;
    for(int i = 0; i < A.size(); i++)
        for(int j = 0; j < A.size(); j++)
            if(i != j && fabsl(A[i][j]) > max)
            {
                max = fabsl(A[i][j]);
                maxPos = std::make_pair(i, j);
            }
    return maxPos;
}

//Calculate the transpose of the matrix
vector<vector<long double>> calAT(vector<vector<long double>>& A)
{
    vector<vector<long double>> AT(A.size(), vector<long double>(A.size()));
    for(int i = 0; i < A.size(); i++)
        for(int j = 0; j < A.size(); j++)
            AT[i][j] = A[j][i];
    return AT;
}

//Calculate the product of A and A^T
vector<vector<long double>> calAAT(vector<vector<long double>>& A)
{
    vector<vector<long double>> AAT(A.size(), vector<long double>(A.size()));
    for(int i = 0; i < A.size(); i++)
        for(int j = 0; j < A.size(); j++)
            for(int k = 0; k < A.size(); k++)
                AAT[i][j] += A[i][k] * A[k][j];
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

//Judge whether the iteration is end
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

//Calculate the eigenvalue of the matrix A with Jacobi method
vector<long double> calEigenValue(vector<vector<long double>>& A)
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

//Column principal element elimination method
vector<long double> Column_Elimination(vector<vector<long double>>& A)
{
    vector<long double> x(99);
    int n = A.size();
    for(int col = 0; col < n; col++)
    {
        long double maxnum = abs(A[col][col]);
        int maxrow = col;
        for(int row = col + 1; row < n; row++)
        {
            if(abs(A[row][col]) > maxnum)
            {
                maxnum = abs(A[row][col]);
                maxrow = row;
            }
        }
        swap(A[col], A[maxrow]);
        for(int row = col + 1; row < n; row++)
        {
            long double res = A[row][col] / A[col][col];
            for(int loc = col; loc <= n; loc++)
                A[row][loc] -= A[col][loc] * res; 
        }
    }
    for(int row = n - 1; row >= 0; row--)
    {
        for(int col = row + 1; col < n; col++)
        {
            A[row][99] -= A[col][99] * A[row][col] / A[col][col];
            A[row][col] = 0;
        }
        A[row][99] /= A[row][row];
        A[row][row] = 1;
    }
    
    return x;
}

//Calculate the eigenvector of the matrix A
vector<vector<long double>> calEigenVector(vector<vector<long double>>& A, vector<long double>& eigenValue)
{
    vector<vector<long double>> tempMartix(A.size(), vector<long double>(A.size()));
    vector<vector<long double>> eigenVector(A.size(), vector<long double>(A.size()));
    for(int i = 0; i < A.size(); i++)
        for(int j = 0; j < A.size(); j++)
            i == j ? tempMartix[i][j] = A[i][j] - eigenValue[i] : tempMartix[i][j] = A[i][j];
}

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
    vector<long double> x =calEigenValue(AAT);
    cout << "EigenValue: ";
    for(int i = 0; i < 3; i++)
        cout << x[i] << " ";
    vector<vector<long double>> eigenVector1 = calEigenVector(AAT, x);
    cout << endl << "EigenVector: " << endl;
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
            cout << eigenVector1[i][j] << " ";
        cout << endl;
    }
    vector<vector<long double>> ATA = calAAT(AT);
    vector<long double> y =calEigenValue(ATA);
    cout << "EigenValue: ";
    for(int i = 0; i < 3; i++)
        cout << y[i] << " ";    
    vector<vector<long double>> eigenVector2 = calEigenVector(ATA, y);
    cout << endl << "EigenVector: " << endl;
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
            cout << eigenVector2[i][j] << " ";
        cout << endl;
    }
    return 0;
}