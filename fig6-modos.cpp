#include <stdio.h>
#include <math.h>

double q12,q22,gama,w,wL,w2,Omega,PassoDoppler;
double dd,h,t,passo,delta,VarreduraDoppler;
double a10,a20,x,y,cosseno,coswl,cosndeltat,LarguraDoppler,z;
int pontos,i,j,k,m,n,modos,quociente;
double a[4],b[4],c[4],k1[4],k2[4],k3[4],k4[4];
double Pi=3.141592653589793;

double f(double a1,double a2,double a12,double b12,double t,int i)
{
    if (i==1) return            +2*Omega*b12*cosseno     +q22*a2;
    if (i==2) return            -2*Omega*b12*cosseno     -q22*a2;
    if (i==3) return  -(w2-dd)*b12                       -q12*a12;
    if (i==4) return   (w2-dd)*a12+Omega*(a2-a1)*cosseno -q12*b12;    
}      

main() 
{
    FILE *arquivo;
    arquivo=fopen("dados.dat","w");
    
    //23/10/2012
    //modificações:
    //sistema de dois níveis
    //muitos campos cw
    //sem onda girante
    //dados
    
    gama=0.0025*1e9*0;
    q22=(2*Pi)*5e6;
    q12=0.5*q22;
    modos=50;
    Omega=10.1*q22;
    passo=10e-12*0.1;
    pontos=40000*10;             
    a10=1;                  
    a20=0;                  
    w2= (2*Pi)*0;
    wL= (2*Pi)*0;
    delta = (2*Pi)*100e6;
    dd=0;
    
    t=0;
    a[1]=a10; a[2]=a20; a[3]=0;  a[4]=0;
    b[1]=0;   b[2]=0;   b[3]=0;  b[4]=0;
    
        for (k=1;k<=pontos;k++)
        {
             t = t + passo;
             h = passo;
             z = delta*t - 1;//(1)
             
             cosndeltat = sin((modos+0.5)*z)/sin(z*0.5);//(1)
             coswl = cos(wL*t);
             cosseno = cosndeltat*coswl;      
                                        
             for (j=1;j<=4;j++)
                 k1[j]=f(a[1],a[2],a[3],a[4],t,j);
                                                        
             for (j=1;j<=4;j++)
                 k2[j]=f(a[1]+k1[1]*h/2,a[2]+k1[2]*h/2,a[3]+k1[3]*h/2,a[4]+k1[4]*h/2,t+h/2,j);
                  
             for (j=1;j<=4;j++)
                 k3[j]=f(a[1]+k2[1]*h/2,a[2]+k2[2]*h/2,a[3]+k2[3]*h/2,a[4]+k2[4]*h/2,t+h/2,j);
                  
             for (j=1;j<=4;j++)           
                 k4[j]=f(a[1]+k3[1]*h,a[2]+k3[2]*h,a[3]+k3[3]*h,a[4]+k3[4]*h,t+h,j);

             for (j=1;j<=4;j++)
                 b[j]=a[j]+h*(k1[j]/6+k2[j]/3+k3[j]/3+k4[j]/6);
                  
             for (m=1;m<=4;m++)
                 a[m]=b[m];     
             
             if (k % 10 == 0)
             {
                   w=b[1]+b[2];  
                   printf("\n");
        
                   printf("%d %12.10f %12.10f %12.10f %12.10f", 
                        k,b[1],b[2],sqrt(b[3]*b[3]+b[4]*b[4]),w);

                   fprintf(arquivo,"%16.14f %16.14f %16.14f %16.14f %16.14f %16.14f",
                        t*1e9,b[1],b[2],sqrt(b[3]*b[3]+b[4]*b[4]),b[3],b[4]);   
                   fprintf(arquivo,"\n");
              }
        
     }      

     //scanf ("%d",&i);
     fclose(arquivo);
     printf("\a");
}
