#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double A1[5][5];
double A2[4][4];
double LU1[5][5];
double LU2[4][4];

void InitialA1(void){
    int i,j,k=9;
    for(i=0;i<5;i++){
        for(j=0;j<5;j++){
            A1[i][j]=1.0/(k-j);
        }
        k--;
    }
}

void InitialA2(void){
    int i,j;
    A2[0][0]=4;
    A2[1][0]=A2[2][0]=16;
    A2[3][0]=6;
    A2[0][1]=A2[2][2]=-1;
    A2[0][2]=1;
    A2[0][3]=3;
    A2[1][1]=A2[1][2]=-2;
    A2[1][3]=5;
    A2[2][1]=-3;
    A2[2][3]=7;
    A2[3][1]=-4;
    A2[3][2]=2;
    A2[3][3]=9;
}

void Doolittle(void){
    int i,j,k,r;
    double temp;
//求矩阵A1的LU分解
    for(k=0;k<5;k++){
        LU1[k][k]=1;
    }
    for(k=0;k<5;k++){
        for(j=k;j<5;j++){
            for(temp=0.0,r=0;r<=k-1;r++){
                temp+=LU1[k][r]*LU1[r][j];
            }
            LU1[k][j]=A1[k][j]-temp;
        }
        for(i=k+1;i<5;i++){
            for(temp=0.0,r=0;r<=k-1;r++){
                temp+=LU1[i][r]*LU1[r][k];
            }
            LU1[i][k]=(A1[i][k]-temp)/LU1[k][k];
        }
    }
//求矩阵A2的LU分解
    for(k=0;k<4;k++){
        LU2[k][k]=1;
    }
    for(k=0;k<4;k++){
        for(j=k;j<4;j++){
            for(temp=0.0,r=0;r<=k-1;r++){
                temp+=LU2[k][r]*LU2[r][j];
            }
            LU2[k][j]=A2[k][j]-temp;
        }
        for(i=k+1;i<4;i++){
            for(temp=0.0,r=0;r<=k-1;r++){
                temp+=LU2[i][r]*LU2[r][k];
            }
            LU2[i][k]=(A2[i][k]-temp)/LU2[k][k];
        }
    }
}

double Norm(double x[],int n){
    int i;
    double norm=fabs(x[0]);
    for(i=1;i<n;i++){
        if(fabs(x[i])>norm){
            norm=fabs(x[i]);
        }
    }
    return norm;
}

void Eigenvalue1(void){
    double y[5]={1,1,1,1,1};
    double z[5];
    double x[5];
    double p_x[5];
    double temp;
    int i,j,k=0;
    //解方程
    for(i=0;i<5;i++){
        for(temp=0.0,j=0;j<=i-1;j++){
            temp+=LU1[i][j]*z[j];
        }
        z[i]=y[i]-temp;
    }
    for(i=4;i>=0;i--){
        for(temp=0.0,j=i+1;j<5;j++){
            temp+=LU1[i][j]*p_x[j];
        }
        p_x[i]=(z[i]-temp)/LU1[i][i];
    }
    printf("y(%d)=(%Lf,%Lf,%Lf,%Lf,%Lf)\tx(%d)=(%Lf,%Lf,%Lf,%Lf,%Lf)\tlamda=%Lf\n",k,y[0],y[1],y[2],y[3],y[4],k+1,p_x[0],p_x[1],p_x[2],p_x[3],p_x[4],Norm(p_x,5));
    //---
    k++;
    for(i=0;i<5;i++){
        y[i]=p_x[i]/Norm(p_x,5);
    }
    //解方程
    for(i=0;i<5;i++){
        for(temp=0.0,j=0;j<=i-1;j++){
            temp+=LU1[i][j]*z[j];
        }
        z[i]=y[i]-temp;
    }
    for(i=4;i>=0;i--){
        for(temp=0.0,j=i+1;j<5;j++){
            temp+=LU1[i][j]*x[j];
        }
        x[i]=(z[i]-temp)/LU1[i][i];
    }
    printf("y(%d)=(%Lf,%Lf,%Lf,%Lf,%Lf)\tx(%d)=(%Lf,%Lf,%Lf,%Lf,%Lf)\tlamda=%Lf\n",k,y[0],y[1],y[2],y[3],y[4],k+1,x[0],x[1],x[2],x[3],x[4],Norm(x,5));
    while(fabs(Norm(p_x,5)-Norm(x,5))>1e-5){
        k++;
        for(i=0;i<5;i++){
            p_x[i]=x[i];
        }
        for(i=0;i<5;i++){
            y[i]=p_x[i]/Norm(p_x,5);
        }
        //解方程
        for(i=0;i<5;i++){
            for(temp=0.0,j=0;j<=i-1;j++){
                temp+=LU1[i][j]*z[j];
            }
            z[i]=y[i]-temp;
        }
        for(i=4;i>=0;i--){
            for(temp=0.0,j=i+1;j<5;j++){
                temp+=LU1[i][j]*x[j];
            }
            x[i]=(z[i]-temp)/LU1[i][i];
        }
        printf("y(%d)=(%Lf,%Lf,%Lf,%Lf,%Lf)\tx(%d)=(%Lf,%Lf,%Lf,%Lf,%Lf)\tlamda=%Lf\n",k,y[0],y[1],y[2],y[3],y[4],k+1,x[0],x[1],x[2],x[3],x[4],Norm(x,5));
    }
    printf("特征向量：v=(%Lf,%Lf,%Lf,%Lf,%Lf)\n迭代次数：%d\n",y[0],y[1],y[2],y[3],y[4],k);
}

void Eigenvalue2(void){
    double y[4]={1,1,1,1};
    double z[4];
    double x[4];
    double p_x[4];
    double temp;
    int i,j,k=0;
    //解方程
    for(i=0;i<4;i++){
        for(temp=0.0,j=0;j<=i-1;j++){
            temp+=LU2[i][j]*z[j];
        }
        z[i]=y[i]-temp;
    }
    for(i=3;i>=0;i--){
        for(temp=0.0,j=i+1;j<4;j++){
            temp+=LU2[i][j]*p_x[j];
        }
        p_x[i]=(z[i]-temp)/LU2[i][i];
    }
    printf("y(%d)=(%Lf,%Lf,%Lf,%Lf)\tx(%d)=(%Lf,%Lf,%Lf,%Lf)\tlamda=%Lf\n",k,y[0],y[1],y[2],y[3],k+1,p_x[0],p_x[1],p_x[2],p_x[3],Norm(p_x,4));
    //---
    k++;
    for(i=0;i<4;i++){
        y[i]=p_x[i]/Norm(p_x,4);
    }
    //解方程
    for(i=0;i<4;i++){
        for(temp=0.0,j=0;j<=i-1;j++){
            temp+=LU2[i][j]*z[j];
        }
        z[i]=y[i]-temp;
    }
    for(i=3;i>=0;i--){
        for(temp=0.0,j=i+1;j<4;j++){
            temp+=LU2[i][j]*x[j];
        }
        x[i]=(z[i]-temp)/LU2[i][i];
    }
    printf("y(%d)=(%Lf,%Lf,%Lf,%Lf)\tx(%d)=(%Lf,%Lf,%Lf,%Lf)\tlamda=%Lf\n",k,y[0],y[1],y[2],y[3],k+1,x[0],x[1],x[2],x[3],Norm(x,4));
    while(fabs(Norm(p_x,4)-Norm(x,4))>1e-5){
        k++;
        for(i=0;i<4;i++){
            p_x[i]=x[i];
        }
        for(i=0;i<4;i++){
            y[i]=p_x[i]/Norm(p_x,4);
        }
        //解方程
        for(i=0;i<4;i++){
            for(temp=0.0,j=0;j<=i-1;j++){
                temp+=LU2[i][j]*z[j];
            }
            z[i]=y[i]-temp;
        }
        for(i=3;i>=0;i--){
            for(temp=0.0,j=i+1;j<4;j++){
                temp+=LU2[i][j]*x[j];
            }
            x[i]=(z[i]-temp)/LU2[i][i];
        }
        printf("y(%d)=(%Lf,%Lf,%Lf,%Lf)\tx(%d)=(%Lf,%Lf,%Lf,%Lf)\tlamda=%Lf\n",k,y[0],y[1],y[2],y[3],k+1,x[0],x[1],x[2],x[3],Norm(x,4));
    }
    printf("特征向量：v=(%Lf,%Lf,%Lf,%Lf)\n迭代次数：%d\n",y[0],y[1],y[2],y[3],k);
}

int main(){
    InitialA1();
    InitialA2();
    Doolittle();
    Eigenvalue1();
    Eigenvalue2();
}