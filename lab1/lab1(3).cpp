#include<iostream>
#include<cmath>
#include <iomanip>
using namespace std;

int main()
	{
	long double Px,Qx,Qy;
	cin>>Px>>Qx>>Qy;
/*int main(int argc, char *argv[])
{
    if(argc != 4) {
        cout << "Usage: " << argv[0] << " <Px> <Qx> <Qy>\n";
        return 1;
    }

    long double Px = stod(argv[1]);
    long double Qx = stod(argv[2]);
    long double Qy = stod(argv[3]);*/

    if(Px >= -1 || Qx >= 0 || Qy <= 0 || Qx * Qx + Qy * Qy <= 1) {
        cout << "输入错误，Px应该小于-1，Qx应该小于0，Qy应该大于0，Qx^2+Qy^2应该大于1\n";
        return 1;
    }
    
    long double k = Qy / (Qx - Px);
    long double Tx = 0;
    long double Ty = 0;

    for(long double x = -1 ; -x >= 1e-7 ; x += 1e-7)
    {
        long double ySquare = 1 - x * x , y = sqrt(1 - x * x);
        long double eq1 = Qy * Qy * (k * x - y) * (k * x - y)  + k * k * Px * Px * ySquare - 2 * k * Qy * Px * y * (k * x - y);
        long double eq2 = k * k * Px * Px * ySquare;
        long double eq3 = Qy * Qy + Qx * Qx + 1 - 2 * y * Qy - 2 * x *Qx;
        long double eq4 = 1 + Px * Px - 2 * x * Px;
        long double res = abs(eq1 * eq4 - eq2 * eq3);
        
        if(res <= 1e-7)
        {
            Tx = x;
            Ty = y;
            break;
        }
    }

    long double eq5 = Qy - Qx * Ty / Tx - Tx / Ty * (2 * Tx - Qx) - 2 * Ty + Qy;
    long double eq6 = -Tx / Ty - Ty / Tx; 
    long double Rx = eq5 / eq6;
    long double Ry = Qy - Ty / Tx * (Qx - Rx);
    cout << "T = (" << fixed << setprecision(6) << Tx << " , " << fixed << setprecision(6) << Ty << ") , R = (" << fixed << setprecision(6) << Rx << " , " << fixed << setprecision(6) << Ry << ")";
    return 0;
}