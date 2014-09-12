#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "RngStream.h"
 
#define SEARCH_AGENTS 30
#define dim 30
#define lb -100
#define ub 100
#define sInf 3425666547544
#define MI 500
#define PI 3.1415926535897932384
#define NS 1000000

typedef struct {
    double pos[dim];
    double fitness;
} Wolf;

Wolf wolf_pack[SEARCH_AGENTS], alpha, beta, delta;

RngStream g1;

double cH[4]={1, 1.2, 3, 3.2};

double aH1[4][6]={{10, 3, 17, 3.5, 1.7, 8},
                 {0.05, 10, 17, 0.1, 8, 14},
                 {3, 3.5, 1.7, 10, 17, 8},
                 {17, 8, 0.05, 10, 0.1, 14}};
double pH1[4][6]={{0.131, 0.169, 0.556, 0.012, 0.828, 0.588},
               {0.232, 0.413, 0.830, 0.373, 0.100, 0.999},
               {0.234, 0.141, 0.352, 0.288, 0.304, 0.665},
               {0.404, 0.882, 0.873, 0.574, 0.109, 0.038}};

double aH[4][3]={{3, 10, 30},
               {0.1, 10, 35},
               {3, 10, 30},
               {0.1, 10, 30}};
double pH[4][3]={{0.3689, 0.117, 0.2673},
              {0.4699, 0.4387, 0.747},
              {0.1091, 0.8732, 0.5547},
              {0.03815, 0.5743, 0.8828}};

void verLU(int index){
    for (int i = 0; i < dim; ++i){
        if (wolf_pack[index].pos[i] > ub || wolf_pack[index].pos[i] < lb){
            wolf_pack[index].pos[i]= (alpha.pos[i]+beta.pos[i]+delta.pos[i])/3;
         } 
    }
}

double sin2(double x){
    return (1-cos(2*x))/2;
}

void showPack(){
        for (int j = 0; j < dim; ++j){
            //printf(" %f ", alpha.pos[j]);        
        }
        //printf("\n");

    printf("Fitness de alpha:\t %.60f\n", alpha.fitness);
    //printf("Fitness de beta:\t %.60f\n", beta.fitness);
    //printf("Fitness de delta:\t %.60f\n", delta.fitness);
    //printf("\n\n");
}

void initializeAgents(){
    alpha.fitness=sInf;
    beta.fitness=sInf;
    delta.fitness=sInf;

    for(int i=0; i<SEARCH_AGENTS;i++){
        for (int j = 0; j < dim; j++){
            wolf_pack[i].pos[j]=RngStream_RandU01(g1)*(ub-lb)+lb;
        }
    }
}

double function(int i, int j){
    return wolf_pack[i].pos[j] * wolf_pack[i].pos[j];
}

double function2(){
    double fitness=0;

    for (int i = 0; i < SEARCH_AGENTS; ++i){
        for(int k=0; k<dim; k++){
            fitness = fitness + fabs(wolf_pack[i].pos[k]);              
        }
    }

    double mult=1;

    for (int i = 0; i < SEARCH_AGENTS; ++i){
        for (int m = 0; m < dim;m++){
            mult = mult*fabs(wolf_pack[i].pos[m]);
        }
    }

    return fitness + mult;
}

double function3(int j, int i){
    double sum=0;
    for (int z = 0; z < i; ++z){
        sum = sum + wolf_pack[j].pos[z];
    }

    return (sum * sum);
}

double function6(int j, int i){
    return pow(wolf_pack[j].pos[i] + 0.5,2);
}

double function7(int j, int i){
    return i * pow(wolf_pack[j].pos[i],4);
}

double function8(int j, int i){
    return wolf_pack[i].pos[j]*sin(sqrt(fabs(wolf_pack[i].pos[j])));
}

double function16(int i){
    double sum1=0, sum2=0, sum3=0;
    for(int j=0; j<dim; j++){
        sum1=sum1+pow(sin(wolf_pack[i].pos[j]),2);
        sum2=sum2+pow(wolf_pack[i].pos[j],2);
        sum3=sum3+sin2(sqrt(fabs(wolf_pack[i].pos[j])));
    }

    return (sum1-exp(-sum2)) * exp(-sum3);
}

double function19(int i){
    double sum=0;
    for (int j = 0; j < dim; ++j){
        sum = sum + aH[i][j]*pow(wolf_pack[i].pos[j]-pH[i][j],2);
    }
    return exp(-sum);
}

double function20(int i){
    double sum=0;
    for (int j = 0; j < dim; ++j){
        sum = sum + aH1[i][j]*pow(wolf_pack[i].pos[j]-pH1[i][j],2);
    }
    return exp(-sum);
}

double GWO(){
    initializeAgents();

    int l=0;
    while(l<MI){

        for(int i=0; i<SEARCH_AGENTS;i++){
            double sumF = 0;

            verLU(i);

            for (int j = 0; j < dim; j++){
                sumF = (sumF + function6(i,j));   
            }

            //sumF=function16(i);

            wolf_pack[i].fitness = sumF;

            if (wolf_pack[i].fitness < alpha.fitness) {
                memcpy(alpha.pos, wolf_pack[i].pos, sizeof(wolf_pack[i].pos));
                alpha.fitness=wolf_pack[i].fitness;
            }

            if (wolf_pack[i].fitness < beta.fitness && wolf_pack[i].fitness > alpha.fitness){
                memcpy(beta.pos, wolf_pack[i].pos, sizeof(wolf_pack[i].pos)); 
                beta.fitness=wolf_pack[i].fitness;   
            }

            if (wolf_pack[i].fitness > alpha.fitness && wolf_pack[i].fitness > beta.fitness && wolf_pack[i].fitness < delta.fitness){
                memcpy(delta.pos, wolf_pack[i].pos, sizeof(wolf_pack[i].pos));        
                delta.fitness=wolf_pack[i].fitness;
            }
        }
        
        double a=2-l*((2)/MI);

        for(int i=0; i<SEARCH_AGENTS;i++){

            double A1[dim],A2[dim],A3[dim],D_alpha[dim],D_beta[dim],D_delta[dim],C1[dim],C2[dim],C3[dim],X1[dim],X2[dim],X3[dim];

            for (int z = 0; z < dim; ++z){

                A1[z]=2*a*RngStream_RandU01(g1)-a;
                C1[z]=2*RngStream_RandU01(g1);

                A2[z]=2*a*RngStream_RandU01(g1)-a;
                C2[z]=2*RngStream_RandU01(g1);

                A3[z]=2*a*RngStream_RandU01(g1)-a;
                C3[z]=2*RngStream_RandU01(g1);
                
                D_alpha[z]=fabs(C1[z]*alpha.pos[z]-wolf_pack[i].pos[z]);
                D_beta[z]=fabs(C2[z]*beta.pos[z]-wolf_pack[i].pos[z]);
                D_delta[z]=fabs(C3[z]*delta.pos[z]-wolf_pack[i].pos[z]);
                
                X1[z]=alpha.pos[i]-A1[z]*D_alpha[z];
                X2[z]=beta.pos[i]-A2[z]*D_beta[z];
                X3[z]=delta.pos[i]-A3[z]*D_delta[z];

                wolf_pack[i].pos[z]=(X1[z]+X2[z]+X3[z])/3;
            }
        }
        l++;
    }
}
 
int main (int argc, char const *argv[])
{
    srand(time(0));
    for (int i = 0; i < 30; ++i){
        unsigned long germe[6] = { rand(), rand(), rand(), rand(), rand(), rand() };
        RngStream_SetPackageSeed (germe);
        g1 = RngStream_CreateStream ("Laplace");
        GWO();
        showPack();
        RngStream_DeleteStream (&g1);   
    }
    return 0;
}