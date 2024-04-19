#include<bits/stdc++.h>
using namespace std;

// 生成一个具有指定行数和列数的随机矩阵
vector<vector<long double>> generateRandomMatrix(int rows, int cols);

// 从矩阵 A 非对角元中选择最大的元素，并返回其位置
pair<int, int> chooseMax(vector<vector<long double>> A);

// 计算矩阵 A 的转置
vector<vector<long double>> calAT(vector<vector<long double>> A);

// 计算两个矩阵的乘积
vector<vector<long double>> multiplyMatrices(const vector<vector<long double>> A, const vector<vector<long double>> B);

// 计算矩阵Q^T * A * Q的每个元素，使用给定的参数 p, q, t, c, d
long double calculateElement(const vector<vector<long double>> A, int i, int j, long double p, long double q, long double t, long double c, long double d);

// 计算矩阵 Q^T * A * Q
vector<vector<long double>> calQTAQ(vector<vector<long double>> A);

// 判断Jacobi方法是否满足结束条件
int judgeEnd(vector<vector<long double>> A);

// Jacobi方法计算矩阵 A 的特征值
vector<long double> calEigenValue(vector<vector<long double>> A);

// 对矩阵 A 进行列主元化成上三角
vector<vector<long double>> Column_Elimination(vector<vector<long double>> A);

// 求解系数矩阵为上三角矩阵A的线性方程组
vector<long double> SolveUpperTriangle(vector<vector<long double>> A, vector<long double> b);

// 解系数矩阵为上三角矩阵 A 的线性方程组，且A全为0的行数为 cnt
vector<vector<long double>> solve(vector<vector<long double>> A, int cnt);

// 计算矩阵 A 的特征向量，使用给定的特征值
vector<vector<long double>> calEigenVector(vector<vector<long double>> A, vector<long double> eigenValue);

// 计算 Sigma 矩阵，使用给定的特征值 x 和矩阵的行数 n1 和列数 n2
vector<vector<long double>> calSigma(vector<long double> x,int n1, int n2);

// 计算向量 x 的欧几里得范数
long double EuclideanNorm(vector<long double> x);

// 对矩阵 A 进行归一化
vector<vector<long double>> Normalization(vector<vector<long double>> A);

// 计算 VT 矩阵
vector<vector<long double>> calculateVT(const vector<vector<long double>>& U, const vector<vector<long double>>& A, const vector<vector<long double>>& Sigma1, int n2);

int main()
{
    vector<vector<long double>> A = generateRandomMatrix(4, 3);
    cout << "A: " << endl;
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 3; j++)
            cout << A[i][j] << " ";
        cout << endl;
    }
    cout << endl;
    vector<vector<long double>> AT = calAT(A);
    vector<vector<long double>> AAT = multiplyMatrices(A, AT);
    int n1 = AAT.size();
    int n2 = A[0].size();
    vector<long double> x =calEigenValue(AAT);
    sort(x.begin(), x.end());
    reverse(x.begin(), x.end());

    vector<vector<long double>> U = calAT(Normalization(calEigenVector(AAT, x)));
    cout << "U: " << endl;
    for(int i = 0; i < n1; i++)
    {
        for(int j = 0; j < n1; j++)
            cout << U[i][j] << " ";
        cout << endl;
    }
    cout << endl;
    
    vector<vector<long double>> Sigma1 = calSigma(x, n1, n2);
    cout << "Sigma: " << endl;
    for(int i = 0; i < n1; i++)
    {
        for(int j = 0; j < n2; j++)
            cout << Sigma1[i][j] << " ";
        cout << endl;
    }
    cout << endl;

    vector<vector<long double>> VT = calculateVT(U, A, Sigma1, n2);
    cout << "VT: " << endl;
    for(int i = 0; i < n2; i++)
    {
        for(int j = 0; j < n2; j++)
            cout << VT[i][j] << " ";
        cout << endl;
    }
    return 0;
}

// 生成一个具有指定行数和列数的随机矩阵
vector<vector<long double>> generateRandomMatrix(int rows, int cols) {
    random_device rd;//随机数种子
    mt19937 gen(rd());//使用 rd 生成的随机数来初始化一个mt19937类型的随机数生成器 gen
    uniform_real_distribution<> dis(0.0, 1.0);//创建一个均匀实数分布 dis

    vector<vector<long double>> matrix(rows, vector<long double>(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            matrix[i][j] = dis(gen);//生成随机数
    return matrix;
}

// 找到矩阵 A 中非对角元中绝对值最大的元素，并返回其位置
pair<int, int> chooseMax(vector<vector<long double>> A)
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

// 计算矩阵 A 的转置
vector<vector<long double>> calAT(vector<vector<long double>> A)
{
    int n1 = A.size();
    int n2 = A[0].size();
    vector<vector<long double>> AT(n2, vector<long double>(n1));
    for(int i = 0; i < n1; i++)
        for(int j = 0; j < n2; j++)
            AT[j][i] = A[i][j];
    return AT;
}

// 计算两个矩阵的乘积
vector<vector<long double>> multiplyMatrices(const vector<vector<long double>> A, const vector<vector<long double>> B) {
    int n1 = A.size();
    int n2 = B[0].size();
    int n3 = A[0].size();
    vector<vector<long double>> result(n1, vector<long double>(n2, 0.0));

    for(int i = 0; i < n1; i++) {
        for(int j = 0; j < n2; j++) {
            for(int k = 0; k < n3; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

// 计算矩阵Q^T * A * Q的每个元素，使用给定的参数 p, q, t, c, d
long double calculateElement(const vector<vector<long double>> A, int i, int j, long double p, long double q, long double t, long double c, long double d) {
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

// 计算矩阵 Q^T * A * Q
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
    return QTAQ;
}

// 判断Jacobi方法是否满足结束条件
int judgeEnd(vector<vector<long double>> A)
{
    int i, j;
    int n = A.size();
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            if(i != j && fabsl(A[i][j]) >= 1e-6)
                return 0;
    if(i == n && j == n) 
        return 1;
}

// Jacobi方法计算矩阵 A 的特征值
vector<long double> calEigenValue(vector<vector<long double>> A)
{
    int n = A.size();
    vector<long double> eigenValue(n);
    vector<vector<long double>> QTAQ= calQTAQ(A);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++)
            cout << QTAQ[i][j] << " ";
            cout << endl;
        }
        cout << endl;
    int i, j;
    while(!judgeEnd(QTAQ)){
        QTAQ = calQTAQ(QTAQ);
        for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            cout << QTAQ[i][j] << " ";   
        }
        cout << endl;
        }
        cout << endl;
    }
    for(i = 0; i < n; i++)
        eigenValue[i] =QTAQ[i][i];
    return eigenValue;
}

// 对矩阵A进行列主元化成上三角
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

// 求解系数矩阵为上三角矩阵A，且A的行列式不为0的线性方程组
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

// 解系数矩阵为上三角矩阵 A 的线性方程组，且A全为0的行数为 cnt
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

// 使用给定的特征值计算矩阵 A 的特征向量
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

        vector<vector<long double>> B = Column_Elimination(tempMartix);
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

// 使用给定的特征值 x 和矩阵的行数 n1 和列数 n2，计算 Sigma 矩阵
vector<vector<long double>> calSigma(vector<long double> x, int n1, int n2)
{
    vector<vector<long double>> Sigma(n1, vector<long double>(n2));
    for(int i = 0; i < min(n1, n2); i++)
        Sigma[i][i] = sqrt(x[i]);
    return Sigma;
}

// 计算向量 x 的欧几里得范数
long double EuclideanNorm(vector<long double> x)
{
    long double norm = 0;
    for(int i = 0; i < x.size(); i++)
        norm += x[i] * x[i];
    return sqrt(norm);
}

// 对矩阵 A 进行归一化
vector<vector<long double>> Normalization(vector<vector<long double>> A)
{
    int rows = A.size();
    for(int i = 0; i < rows; i++)
    {
        long double norm = EuclideanNorm(A[i]);
        int cols = A[i].size();
        for(int j = 0; j < cols; j++)
            A[i][j] /= norm;
    }
    return A;
}

// 计算 VT 矩阵
vector<vector<long double>> calculateVT(const vector<vector<long double>>& U, const vector<vector<long double>>& A, const vector<vector<long double>>& Sigma1, int n2) {
    vector<vector<long double>> UTA = multiplyMatrices(calAT(U), A);
    vector<vector<long double>> VT(n2, vector<long double>(n2));
    for(int i = 0; i < n2; i++)
        for(int j = 0; j < n2; j++)
            VT[i][j] = UTA[i][j]/Sigma1[i][i];
    return VT;
}