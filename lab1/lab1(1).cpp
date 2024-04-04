#include<bits/stdc++.h>
using namespace std;

const long double PI = 3.14159265358979323;

long double calculateEquation(long double Px, long double Qx, long double Qy, long double theta) {
    long double eq1 = Px * sin(theta);
    long double eq2 = sqrt(Qx * Qx + Qy * Qy +1 - 2 * Qx * cos(theta) - 2 * Qy * sin(theta));
    long double eq3 = (Qx * sin(theta) - Qy * cos(theta));
    long double eq4 = sqrt(Px * Px -2 * Px * cos(theta) + 1);
    return eq1 * eq2 + eq3 * eq4;
}

int main(int argc, char *argv[]) {
    if(argc != 4) {
        cout << "Usage: " << argv[0] << " <Px> <Qx> <Qy>\n";
        return 1;
    }

    long double Px = stold(argv[1]);
    long double Qx = stold(argv[2]);
    long double Qy = stold(argv[3]);

    if(Px >= -1 || Qx >= 0 || Qy <= 0 || Qx * Qx + Qy * Qy <= 1) {
        cout << "输入错误，Px应该小于-1，Qx应该小于0，Qy应该大于0，Qx^2+Qy^2应该大于1\n";
        return 1;
    }

    long double Tx = 0;
    long double Ty = 0;
    long double theta = 0;
    long double low = PI, high = PI / 2;
    long double k = 0;

    while(1) {
        theta = (low + high) / 2;
        long double res = calculateEquation(Px, Qx, Qy, theta);
        long double absres = abs(res);
        if(absres <= 1e-7) {
            Tx = cos(theta);
            Ty = sin(theta);
            k = 1 / tan(PI - theta);
            break;
        }
        else 
            res > 0 ? low = theta : high = theta;
    }

    long double eq5 = 2 * Qy - Qx * tan(theta) + k * (2 * Tx - Qx) - 2 * Ty;
    long double eq6 = k - tan(theta); 
    long double Rx = eq5 / eq6;
    long double Ry = Qy - tan(theta) * (Qx - Rx);
    cout << "T = (" << fixed << setprecision(6) << Tx << " , " << fixed << setprecision(6) << Ty << ") , R = (" << fixed << setprecision(6) << Rx << " , " << fixed << setprecision(6) << Ry << ")";
    return 0;
}