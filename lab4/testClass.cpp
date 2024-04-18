#include <iostream>
#include <vector>
#include <random>
#include <cmath>

using namespace std;

template<typename T>
class Matrix {
private:
    vector<vector<T>> A;

    // Function to find the maximum element in the matrix
    pair<int, int> chooseMax(vector<vector<T>>& A) {
        T max = 0;
        pair<int, int> maxPos;
        for(int i = 0; i < A.size(); i++)
            for(int j = 0; j < A.size(); j++)
                if(i != j && fabs(A[i][j]) > max)
                {
                    max = fabs(A[i][j]);
                    maxPos = make_pair(i, j);
                }
        return maxPos;
    }

    // Function to calculate the transpose of the matrix
    vector<vector<T>> calAT(vector<vector<T>>& A) {
        vector<vector<T>> AT(A.size(), vector<T>(A.size()));
        for(int i = 0; i < A.size(); i++)
            for(int j = 0; j < A.size(); j++)
                AT[i][j] = A[j][i];
        return AT;
    }

    // Function to calculate the element of the matrix
    T calculateElement(const vector<vector<T>>& A, int i, int j, T p, T q, T t, T c, T d) {
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

    // Function to calculate the QTAQ of the matrix
    vector<vector<T>> calQTAQ(vector<vector<T>>& A) {
        pair<int, int> maxPos = chooseMax(A);
        int row = maxPos.first;
        int col = maxPos.second;
        
        T s = (A[col][col] - A[row][row]) / (2 * A[row][col]);
        T t = 0;
        if(s == 0)
            t = 1;
        else if(abs(-s + sqrt(1 + s * s)) <= abs(-s - sqrt(1 + s * s)))
            t = -s + sqrt(1 + s * s);
        else
            t = -s - sqrt(1 + s * s);

        T c = 1 / sqrt(1 + t * t);
        T d = t * c;

        vector<vector<T>> QTAQ(A.size(), vector<T>(A.size()));
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

    // Function to check if the calculation should end
    int judgeEnd(vector<vector<T>>& A) {
        int i, j;
        for(i = 0; i < A.size(); i++)
            for(j = 0; j < A.size(); j++)
                if(i != j && fabs(A[i][j]) >= 1e-7)
                    return 0;
        if(i == A.size() && j == A.size()) 
            return 1;
    }

public:
    // Constructor
    Matrix(vector<vector<T>> matrix) : A(matrix) {}

    // Function to calculate the eigenvalue of the matrix
    vector<long double> calEigenValue() {
        vector<long double> eigenValue(A.size());
        vector<vector<long double>> QTAQ= calQTAQ(A);
        int i, j;
        while(!judgeEnd(QTAQ))
            QTAQ = calQTAQ(QTAQ);
        for(i = 0; i < A.size(); i++)
            eigenValue[i] = QTAQ[i][i];
        return eigenValue;
    }

    // Function to print the eigenvalues of the matrix
    void printEigenValues() {
        vector<long double> x = calEigenValue();
        cout << "EigenValue: ";
        for(int i = 0; i < A.size(); i++)
            cout << x[i] << " ";
    }
};

int main() {
    vector<vector<long double>> A = {
        {1.31336, 2.02598, 1.08122, 0},
        {1.22657, 1.94428, 1.00043, 0},
        {0.829463, 1.31122, 0.841318, 0},
        {0.556764, 0.895583, 0.584515, 0}
    };

    // Create a matrix and print its eigenvalues
    Matrix<long double> matrix(A);
    matrix.printEigenValues();
    return 0;
}