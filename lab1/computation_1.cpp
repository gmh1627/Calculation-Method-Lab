#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double function1(double xp, double xq, double yq, double xt, double yt){
    double r1=xt*xt*xt*yq;
    double r2=xt*xt*xp*yq;
    double r3=xt*yt*xp*xq;
    double r4=yt*yt*yt*(xq+xp);
    double r5=xt*xt*yt*(xp+xq);
    double r6=yt*yt*xp*yq;
    double r7=xt*yt*yt*yq;
    double r=r1-r2+2*r3-r4-r5+r6+r7;
    return r;
}

double ReflectingPoint(double xp, double xq, double yq){
    double xt1,yt1,xt2;
    double theta=atan(yq/(-xq));
    xt1=-cos(theta);
    xt2=-1.0;
    while((xt1-xt2)>1e-6){
        yt1=sqrt(1-xt1*xt1);
        double mid_xt=(xt1+xt2)/2;
        double mid_yt=sqrt(1-mid_xt*mid_xt);
        double mid_re=function1(xp,xq,yq,xt1,yt1)*function1(xp,xq,yq,mid_xt,mid_yt);
        printf("%lf",mid_re);
        if(mid_re==0){
            return mid_xt;
        }
        else if(mid_re<0){
            xt2=mid_xt;
        }
        else{
            xt1=mid_xt;
        }
    }
    return (xt1+xt2)/2;
}

int main(){
    double xp, xq, yq;
    scanf("%Lf,%Lf,%Lf",&xp,&xq,&yq);
    printf("%Lf",ReflectingPoint(xp,xq,yq));
}