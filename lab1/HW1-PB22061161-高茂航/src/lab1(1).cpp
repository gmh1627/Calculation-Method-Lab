#include<iostream>
#include<cmath>
#include <iomanip>
using namespace std;

int main(int argc, char *argv[])
{
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
    long double low = 3.14159265358979323 , high = 1.570796326794896;
    long double k = 0;
    while(1)
    {
        theta = (low + high) / 2;
        long double eq1 = Px * sin(theta);
        long double eq2 = sqrt(Qx * Qx + Qy * Qy +1 - 2 * Qx * cos(theta) - 2 * Qy * sin(theta));
        long double eq3 = (Qx * sin(theta) - Qy * cos(theta));
        long double eq4 = sqrt(Px * Px -2 * Px * cos(theta) + 1);
        long double res = eq1 * eq2 + eq3 * eq4;
        long double absres = abs(res);
        if(absres <= 1e-7)
        {
            Tx = cos(theta);
            Ty = sin(theta);
            k = 1 / tan(3.14159265358979323 - theta);
            break;
        }
        else if(res > 0)
        {
            low = theta;
        }
        else
        {
            high = theta;
        }
    }
    long double eq5 = 2 * Qy - Qx * tan(theta) + k * (2 * Tx - Qx) - 2 * Ty;
    long double eq6 = k - tan(theta); 
    long double Rx = eq5 / eq6;
    long double Ry = Qy - tan(theta) * (Qx - Rx);
    cout << "T = (" << fixed << setprecision(6) << Tx << " , " << fixed << setprecision(6) << Ty << ") , R = (" << fixed << setprecision(6) << Rx << " , " << fixed << setprecision(6) << Ry << ")";
    return 0;
}