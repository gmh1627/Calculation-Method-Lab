#include<bits/stdc++.h>
using namespace std;

int satisfiedCount;
int M;
long double eps = 1e-6, proportion;
vector<long double> Vx(100), Vy(100);

void Calculate(int M);
long double ax(long double t);
long double ay(long double t);
long double vx(long double t);
long double vy(long double t);
long double romberg(function<long double(long double)> f, long double a, long double b, long double eps, int M, bool isX);// Perform the Romberg integration

int main() 
{
    Calculate(4);
    Calculate(8);
    Calculate(12);
    Calculate(16);
    Calculate(20);

    return 0;
}

void Calculate(int M) {
    satisfiedCount = 0; 
    cout << "M = " << M << endl;
    ofstream outFile("trajectory.txt", ios::app);
    for (long double t = 0.1; t < 10.1; t += 0.1) 
    { 
        Vx[10*t-1] = romberg(ax, 0, t, eps, M, 0);
        Vy[10*t-1] = romberg(ay, 0, t, eps, M, 0);
        long double x = romberg(vx, 0, t, eps, M, 1);
        long double y = romberg(vy, 0, t, eps, M, 0);
        cout << fixed << setprecision(1) << "At t = " << t << ", vx = " << setprecision(6) << Vx[10*t-1] << ", vy = " << setprecision(6) << Vy[10*t-1] << ", " << "(x, y) = (" << setprecision(6) << x << ", " << setprecision(6) << y << ")" << endl;

        if (M == 8)
            outFile << fixed << setprecision(6) << x << " " << y << "\n";
    }
    long double proportion = (long double)satisfiedCount / 100;
    cout << "At M = " << M << ", proportion of times the error requirement of x was satisfied: " << proportion << endl;
}

long double ax(long double t) 
{
    return sin(t) / (1 + sqrt(t));
}

long double ay(long double t) 
{
    return log(t + 1) / (t + 1);
}

long double vx(long double t) 
{
    return t == 0 ? 0 : Vx[ceil(10*t)-1];
}

long double vy(long double t) 
{
    return t == 0 ? 0 : Vy[ceil(10*t)-1];
}

// Perform the Romberg integration
long double romberg(function<long double(long double)> f, long double a, long double b, long double eps, int M, bool isX) {
    long double h[M+1], r[M+1][M+1];
    h[1] = b - a;
    r[1][1] = 0.5 * h[1] * (f(a) + f(b));

    for (int k = 2; k <= M; k++) 
    {
        h[k] = 0.5 * h[k-1];
        long double sum = 0;
        for (int i = 1; i <= pow(2, k-2); i++)
            sum += f(a + (2*i-1) * h[k]);
        r[k][1] = 0.5 * (r[k-1][1] + h[k-1] * sum);

        for (int j = 2; j <= k; j++)
            r[k][j] = r[k][j-1] + (r[k][j-1] - r[k-1][j-1]) / (pow(4, j-1) - 1);

        if (k > 2 && fabs(r[k][k] - r[k-1][k-1]) < eps)
        {
            if(isX)
                satisfiedCount++;
            return r[k][k];
        }
    }
    return r[M][M];
}