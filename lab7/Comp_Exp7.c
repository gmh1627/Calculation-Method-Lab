#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define printf __mingw_printf

long double Romberg(long double a, long double b, long double e, int M, int m);

long double Compute(int m, long double x, long double e, int M){
    if(m==1){
        return sin(x)/(sqrt(x)+1);
    }
    else if(m==2){
        return log(x+1)/(x+1);
    }
    else if(m==3){
        return Romberg(0, x, e, M ,1);
    }
    else if(m==4){
        return Romberg(0, x, e, M ,2);
    }
}

long double Romberg(long double a, long double b, long double e, int M, int m){
    int n=1;
    long double h=b-a, sum;
    int i,j,k;
    long double R_1[M], R_2[M], hk[M];//用R_1表示R_{k-1,j},R_2表示R_{k,j}
    for(k=1;k<=M;k++){
        hk[k-1]=h/(pow(2, k-1));
    }
    for(k=2;k<=M;k++){
        if(k==2){
            R_1[0]=(Compute(m, a, e, M)+Compute(m, b, e, M))*h/2;
        }
        sum=0;
        for(i=1;i<=pow(2, k-2);i++){
            //printf("%Lf,%Lf,%Lf\t",a+(2*i-1)*hk[k-1],Compute(m, a+(2*i-1)*hk[k-1], e, M),Compute(2, a+(2*i-1)*hk[k-1], e, M));
            sum+=Compute(m, a+(2*i-1)*hk[k-1], e, M);
        }
        //printf("%Lf,%d",sum,m);
        R_2[0]=(R_1[0]+hk[k-2]*sum)/2;
        //printf("%Lf",R_2[0]);
        for(j=2;j<=k;j++){
            R_2[j-1]=R_2[j-2]+(R_2[j-2]-R_1[j-2])/(pow(4, j-1)-1);
            //printf("%Lf",R_2[j-1]);
            //if(k==j) printf("\n"); 
        }
        if(k>2&&fabs(R_2[k-1]-R_1[k-2])<e&&(m==3||m==4)){
            //printf("%Lf\t%Lf",R_2[k-1],R_1[k-2]);
            break;
        }
        for(i=1;i<=k;i++){
            R_1[i-1]=R_2[i-1];
        }
    }
    //printf("%Lf",R_2[M-1]);
    return R_2[k-1];
}

int main(){
    long double i;
    for(i=0.1;i<=10;i+=0.1){
        printf("(%Lf,%Lf)\n",Romberg(0, i, 1e-6, 12, 3),Romberg(0, i, 1e-6, 12, 4));
    }
}