#include<bits/stdc++.h>

using namespace std;

const double PI = acos(-1);

complex<double> f1(double t);
complex<double> f2(double t);
void fft(vector<complex<double>> &a);
void ifft(vector<complex<double>> &a);

int main() {
    ofstream fft_file;
    ofstream discrete;
    ofstream ifft_file;
    cout << "Function 1:\n";
    for (int n : {16, 128}) {
        discrete.open(n == 16 ? "discrete_f1(1).txt" : "discrete_f1(2).txt");
        vector<complex<double>> g1(n);
        for (int i = 0; i < n; i++) {
            g1[i] = f1(i / (double)n);
            discrete << i / (double)n << ' ' << g1[i].real() << endl;
        }
        discrete.close();

        fft(g1);
        
        fft_file.open(n == 16 ? "fft_f1(1).txt" : "fft_f1(2).txt");
        for (int i = 0; i < n; i++)
            fft_file << i << ' ' << abs(g1[i]) << endl;
        fft_file.close();
        cout << "FFT of function for n = " << n << endl;
        cout << "Real part:\n";
        for (const auto& x : g1)
            cout << x.real() << ' ';
        cout << endl;
        cout << "Imaginary part:\n";
        for (const auto& x : g1)
            cout << x.imag() << ' ';
        cout << endl << endl;

        vector<complex<double>> inv1 = g1;
        ifft(inv1);
        ifft_file.open(n == 16 ? "ifft_f1(1).txt" : "ifft_f1(2).txt");
        cout << "IFFT of function for n = " << n << ":\n";
        for (int i = 0; i < n; i++) {
            cout << inv1[i].real() << ' ';
            ifft_file << i / (double)n << ' ' << inv1[i].real() << endl;
        }
        ifft_file.close();
        cout << endl << endl;
    }
    fft_file.close();


    cout << "Function 2:\n";
    vector<complex<double>> g2(128);
    discrete.open("discrete_f2.txt");
    for (int i = 0; i < 128; i++) {
        g2[i] = f2((double)i / 128);
        discrete << i / 128.0 << ' ' << g2[i].real() << endl;
    }
    discrete.close();

    fft(g2);

    fft_file.open("fft_f2.txt");
    for (int i = 0; i < 128; i++)
        fft_file << i << ' ' << abs(g2[i]) << endl;
    cout << "FFT of function for n = " << 128 << endl;
    cout << "Real part:\n";
    for (const auto& x : g2)
        cout << x.real() << ' ';
    cout << endl;
    cout << "Imaginary part:\n";
    for (const auto& x : g2)
        cout << x.imag() << ' ';
    cout << endl;
    
    vector<complex<double>> inv2 = g2;
    ifft(inv2);
    ifft_file.open("ifft_f2(1).txt");
    cout << "IFFT of function for n = " << 128 << ":\n";
    for (int i = 0; i < 128; i++) {
        cout << inv2[i].real() << ' ';
        ifft_file << (double)i / 128 << ' ' << inv2[i].real() << endl;
    }
    cout << endl << "imaginary part:\n";
    for (int i = 0; i < 128; i++)
        cout << inv2[i].imag() << ' ';
    ifft_file.close();

    vector<complex<double>> inv3(g2.begin(), g2.begin() + 32);
    inv3.resize(128, complex<double>(0, 0));
    ifft(inv3);
    ifft_file.open("ifft_f2(2).txt");
    cout << endl << "IFFT of the first 25% of frequencies" << endl;
    for (int i = 0; i < 128; i++) {
        cout << inv3[i].real() << ' ';
        ifft_file << (double)i / 128 << ' ' << inv3[i].real() << endl;
    }
    cout << endl << "imaginary part:\n";
    for (int i = 0; i < 128; i++)
        cout << inv3[i].imag() << ' ';
    ifft_file.close();
}

complex<double> f1(double t) {
    return 0.7 * sin(2 * PI * 2 * t) + sin(2 * PI * 5 * t);
}

complex<double> f2(double t) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0, 1);
    return 0.7 * sin(2 * PI * 2 * t) + sin(2 * PI * 5 * t) + 0.3 * dis(gen);
}

void fft(vector<complex<double>> &a) {
    int n = a.size();
    if (n == 1) 
        return;

    vector<complex<double>> a0(n / 2), a1(n / 2);
    for (int i = 0; 2 * i < n; i++) {
        a0[i] = a[2 * i];
        a1[i] = a[2 * i + 1];
    }

    fft(a0);
    fft(a1);
    
    double ang = 2 * PI / n;
    complex<double> w(1), wn(cos(ang), sin(ang));
    for (int i = 0; 2 * i < n; i++) {
        a[i] = a0[i] + w * a1[i];
        a[i + n / 2] = a0[i] - w * a1[i];
        w *= wn;
    }
}

/*void ifft(vector<complex<double>> &a) {
    int n = a.size();
    reverse(a.begin() + 1, a.end());

    fft(a);

    for (auto& x : a)
        x /= n;
}*/

void ifft(vector<complex<double>> &a) {
    int n = a.size();
    if (n == 1)
        return;

    vector<complex<double>> a0(n / 2), a1(n / 2);
    for (int i = 0; 2 * i < n; i++) {
        a0[i] = a[2 * i];
        a1[i] = a[2 * i + 1];
    }

    ifft(a0);
    ifft(a1);

    double ang = -2 * PI / n; // negative sign
    complex<double> w(1), wn(cos(ang), sin(ang));
    for (int i = 0; 2 * i < n; i++) {
        a[i] = a0[i] + w * a1[i];
        a[i + n / 2] = a0[i] - w * a1[i];
        w *= wn;
    }

    for (auto& x : a)
        x /= 2;
}