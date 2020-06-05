#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>
float sinx(float x){
    return sin(x);
}

float MonteCarlos(float a, float b, float (*f)(float), int n){
    srand((unsigned int)time(NULL));
    float h = b - a;
    float soma = 0;
    for (int i = 0; i < n; i++){
        soma += f(((float)rand()/(float)(RAND_MAX)) * h);
    }
    return soma/n;
}

float Lagrange(float x, float xi[], float yi[], int n){
    float p1, p2, soma = 0;
    for (int i = 0; i < n; i++){
        p1 = 1;
        p2 = 1;
        for(int j=0; j<n; j++)
        {
            if(j!=i)
            {
                p1 = p1 *(x - xi[j]);
                // printf("p1 %.3f\n", p1);
                p2 = p2 *(xi[i] - xi[j]);
                // printf("p2 %.3f\n", p2);
            }
        }
        soma += ((p1/p2)*yi[i]);
    }
   // printf("soma %.3f\n", soma);
    return soma;
}
//Integração no intervalo a, b com n trapezios calculando separadamente a interpolacao de F(x) e Cos(x)
float TrapezioComposto(float *xi, float *fi, float *theta, int n, int n_dados){
    float b = xi[n_dados - 1];
    float a = xi[0];
    float h = (b - a)/n;
    printf("a = %.3f b = %.3f h = %.3f\n", a, b, h);
    float total = 0;
    for(int i = 0; i < n; i++){
        //calculo de um trapezio de lados f(x0 + i*h)*cos(x0 + i*h) e f(x0 + (i+1)*h)*cos(x0 + (i+1)*h)
        total += (Lagrange(a + i*h, xi, fi, n_dados)*cos(Lagrange(a + i*h, xi, theta, n_dados))
                + Lagrange(a + (i+1)*h, xi, fi, n_dados)*cos(Lagrange(a + (i+1)*h, xi, theta, n_dados)))*h/2;
    }
    printf("Trapezio %.3f\n", total);
    return total;
}

float SimpsonComposto(float *xi, float *fi, float *theta, int n, int n_dados){
    if (n%2 != 0){
        printf("n precisa ser par\n");
        return 0;
    }
    float b = xi[n_dados - 1];
    float a = xi[0];
    float h = (b - a)/n;
    float total = (h/3)*Lagrange(a, xi, fi, n_dados)*cos(Lagrange(a, xi, theta, n_dados));
    for(int i = 1; i < n; i++){
        if (i%2 == 0)
            total += 2*(h/3)*Lagrange(a + i*h, xi, fi, n_dados)*cos(Lagrange(a + i*h, xi, theta, n_dados));
        else
            total += 4*(h/3)*Lagrange(a + i*h, xi, fi, n_dados)*cos(Lagrange(a + i*h, xi, theta, n_dados));
    }
    total += (h/3)*Lagrange(b, xi, fi, n_dados)*cos(Lagrange(b, xi, theta, n_dados));
    printf("Simpson %.3f\n", total);
    return total;
}


int main(){
    float xi[7] = {0, 5, 10, 15, 20, 25, 30};
    float f[7] = {0, 9, 13, 14, 10.5, 12, 5};
    float theta[7] = {0.5, 1.4, 0.75, 0.9, 1.3, 1.48, 1.5};
    float fcos[7] = {0, 1.5297, 9.5120, 8.7025, 2.8087, 1.0881, 0.3537};
    //foi optado interpolar o f(x) e cos(x) separadamente, pois pelos
    //teste foi observado que calculando-os separadamente se obtem uma intepolacao melhor de f(x)cos(x)
    TrapezioComposto(xi,f, theta, 24, 7);
    SimpsonComposto(xi,f, theta, 24, 7);
    printf("MonteCarlos sin(x) entre 0 e 1 %.5f", MonteCarlos(0, 1, sinx, 10000));
    return 0;
}