#include<bits/stdc++.h>
using namespace std;

// Define the acceleration functions ax and ay
long double ax(long double t) {
    return sin(t) / (1 + sqrt(t));
}

long double ay(long double t) {
    return log(t + 1) / (t + 1);
}

// Function to perform the Romberg integration
long double romberg(long double (*f)(long double), long double a, long double b, long double eps, int maxIter) {
    long double h[maxIter], r[maxIter][maxIter];
    h[0] = b - a;
    r[0][0] = 0.5 * h[0] * (f(a) + f(b));

    for (int i = 1; i < maxIter; i++) {
        h[i] = 0.5 * h[i-1];
        long double sum = 0;
        for (int k = 0; k < pow(2, i-1); k++)
            sum += f(a + (2*k+1) * h[i]);
        r[i][0] = 0.5 * r[i-1][0] + h[i] * sum;

        for (int j = 1; j <= i; j++)
            r[i][j] = r[i][j-1] + (r[i][j-1] - r[i-1][j-1]) / (pow(4, j) - 1);

        if (i > 1 && fabs(r[i][i] - r[i-1][i-1]) < eps)
            return r[i][i];
    }

    return r[maxIter-1][maxIter-1];
}

// Main function
int main() {
    long double eps = 1e-6;
    int maxIter = 8;

    for (long double t = 0.1; t <= 10; t += 0.1) {
        long double x = romberg(ax, 0, t, eps, maxIter);
        long double y = romberg(ay, 0, t, eps, maxIter);
        cout << fixed << setprecision(6) << "At t = " << t << ", (x, y) = (" << x << ", " << y << ")" << endl;
    }

    return 0;
}