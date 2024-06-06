#include<iostream>
#include<complex>
#include<vector>
#include<cmath>
#include<bits/stdc++.h>
using namespace std;
#define PI 3.14159265358979323846
typedef double Double;

vector<complex<Double>> Generate_f1(int n){
    vector<complex<Double>> f(n);
    int i;
    double t;
    for(i=0;i<n;i++){
        t=static_cast<Double>(i)/n;
        f[i]=complex<Double>(0.7*sin(2*PI*2*t)+sin(2*PI*5*t),0.0);
    }
    return f;
}

vector<complex<Double>> Generate_f2(int n){
    vector<complex<Double>> f(n);
    int i;
    double t;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0,1);
    for(i=0;i<n;i++){
        t=static_cast<Double>(i)/n;
        f[i]=complex<Double>(0.7*sin(2*PI*2*t)+sin(2*PI*5*t),0.0)+0.3*dis(gen);
    }
    return f;
}

vector<complex<Double>> FFT(const vector<complex<Double>>& f){
    int n=f.size();
    if(n==1)
        return f;
    complex<Double> wn=exp(complex<Double>(0.0,-2*PI/n));
    complex<Double> w(1,0);
    int i,k;
    vector<complex<Double>> f0(n/2),f1(n/2);
    for(i=0;i<n/2;i++){
        f0[i]=f[2*i];
        f1[i]=f[2*i+1];
    }
    vector<complex<Double>> g0=FFT(f0);
    vector<complex<Double>> g1=FFT(f1);
    vector<complex<Double>> g(n);
    for(k=0;k<n/2;k++){
        g[k]=g0[k]+w*g1[k];
        g[k+n/2]=g0[k]-w*g1[k];
        w*=wn;
    }
    return g;
}

vector<complex<Double>> IFFT(const vector<complex<Double>>& f){
    int n=f.size();
    if(n==1)
        return f;
    complex<Double> wn=exp(complex<Double>(0.0,2*PI/n));
    complex<Double> w(1,0);
    int i,k;
    vector<complex<Double>> f0(n/2),f1(n/2);
    for(i=0;i<n/2;i++){
        f0[i]=f[2*i];
        f1[i]=f[2*i+1];
    }
    vector<complex<Double>> g0=IFFT(f0);
    vector<complex<Double>> g1=IFFT(f1);
    vector<complex<Double>> g(n);
    for(k=0;k<n/2;k++){
        g[k]=g0[k]+w*g1[k];
        g[k+n/2]=g0[k]-w*g1[k];
        w*=wn;
    }
    for (i= 0;i<n;i++) {
        g[i]/=2; // IFFT需要除以2
    }
    return g;
}

vector<complex<Double>> Low_Frequency(vector<complex<Double>> &g, int n){
    int i;
    vector<complex<Double>> g1(n);
    for(i=0;i<n;i++){
        if(i<n*0.25)
            g1.push_back(g[i]);
        else
            g1.push_back(complex<Double>(0.0,0.0));
    }
    return g1;
}

int main(){
    int i;
    //f1，当n=16时
    vector<complex<Double>> f1_1=Generate_f1(16);
    vector<complex<Double>> g1_1=FFT(f1_1);
    cout<<"Fuction 1:"<<endl;
    cout<<"FFT of function for n=16:"<<endl;
    cout<<"Real part:"<<endl;
    for(i=0;i<16;i++)
        cout<<g1_1[i].real()<<"\t";
    cout<<endl;
    cout<<"Imaginary part:"<<endl;
    for(i=0;i<16;i++)
        cout<<g1_1[i].imag()<<"\t";
    cout<<endl;
    vector<complex<Double>> h1_1=IFFT(g1_1);
    cout<<"Fuction 1:"<<endl;
    cout<<"IFFT of function for n=16:"<<endl;
    for(i=0;i<16;i++)
        cout<<h1_1[i]<<"\t";
    cout<<endl;
    //f1，当n=128时
    vector<complex<Double>> f1_2=Generate_f1(128);
    vector<complex<Double>> g1_2=FFT(f1_2);
    cout<<"FFT of function for n=128:"<<endl;
    cout<<"Real part:"<<endl;
    for(i=0;i<128;i++)
        cout<<g1_2[i].real()<<"\t";
    cout<<endl;
    cout<<"Imaginary part:"<<endl;
    for(i=0;i<128;i++)
        cout<<g1_2[i].imag()<<"\t";
    cout<<endl;
    vector<complex<Double>> h1_2=IFFT(g1_2);
    cout<<"Fuction 1:"<<endl;
    cout<<"IFFT of function for n=128:"<<endl;
    for(i=0;i<128;i++)
        cout<<h1_2[i]<<"\t";
    cout<<endl;
    //f2
    vector<complex<Double>> f2=Generate_f1(128);
    vector<complex<Double>> g2=FFT(f2);
    cout<<"Fuction 2:"<<endl;
    cout<<"FFT of fusnction for n=128:"<<endl;
    cout<<"Real part:"<<endl;
    for(i=0;i<128;i++)
        cout<<g2[i].real()<<"\t";
    cout<<endl;
    cout<<"Imaginary part:"<<endl;
    for(i=0;i<128;i++)
        cout<<g2[i].imag()<<"\t";
    cout<<endl;
    vector<complex<Double>> h2=IFFT(Low_Frequency(g2,128));
    cout<<"Fuction 2:"<<endl;
    cout<<"IFFT of function for n=128:"<<endl;
    for(i=0;i<128;i++) 
        cout<<h2[i]<<"\t";
    cout<<endl;
}