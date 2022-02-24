#include <cmath>

#include "mosquito.h"

double mosquito::SporozoiteRate(double a, double x, double P, double n){
//calculate the proportion of anopheline mosquitoes with sporozites in 
//their salivary glands which are actually infective.
//the formular was proposed by MacDonald(1957):
//    S = (a x P^n)/(a x - ln(p))
//    S - the proportion
//    P - the probability of mosquito survival through 1 day
//    n - the duration, in days, of the extrinsic cycle of 
//        the parasite in the mosquito
//    a - average number of blood meals on man per day
//    x - the proportion of bites infective to man
    double S; //the proportion (return value)
    
    S = ( a * x * pow(P,n) ) / ( a * x - log(P) );
    
    return S;      
}


double mosquito::InoculationRate(double m, double a, double b, double S){
//Calculate the inoculation rate, or the mean daily number of bites(h)
//received by sporozoite-bearing mosquitoes.
//the formula for h is
//    h = m a b S
//    h - the inoculation rate
//    m - anopheline density in relation to man
//    a - average number of blood meals on man per day
//    b - proportion of bites that are infectious
//    S - sporozoite rate

      double h ; //inoculation rate (return value)
      
      h = m * a * b * S;
      
      return h;
}


double mosquito::ReproductiveRate(double m,double a,double b,double x,double P,double n,double z){
//Calculate the reproductive rate of the infection (r) or the number of 
//secondary cases resulting from a primary case.
//the formula is 
//    r = {(m a^2 b P^n)/(-z ln(P) )}{1 - ((a x) /(a x - ln(P)))}
//    r - the reproductive rate 
//    m - anopheline density in relation to man
//    a - average number of blood meals on man per day
//    b - proportion of bites that are infectious
//    x - the proportion of bites infective to man
//    P - the probability of mosquito survival through 1 day
//    n - the duration, in days, of the extrinsic cycle of 
//        the parasite in the mosquito
//    z - the recovery rate, or the reciprocal of the duration of human infectivity

      double r; // reproductive rate (return value)
      
      r = ((m * a*a * b * pow(P,n))/(-z * log(P) ))*(1 - ((a * x) /(a * x - log(P))));
      
      return r;
}


double mosquito::ReproductiveRate0(double m,double a,double b,double P,double n,double z){
//calculate the reproductive rate when the transmission is low(i.e., x ->0)
//the formula is
//    r0 = (m a^2 b P^n)/(-z ln(P))
//    r0 - the reproductive rate at low transmission
//    m - anopheline density in relation to man
//    a - average number of blood meals on man per day
//    b - proportion of bites that are infectious
//    P - the probability of mosquito survival through 1 day
//    n - the duration, in days, of the extrinsic cycle of 
//        the parasite in the mosquito
//    z - the recovery rate, or the reciprocal of the duration of human infectivity

      double r0; //reproductive rate at low transmission (return value)
      
      r0 = (m * a*a * b * pow(P,n))/(-z * log(P));
      
      return r0;     
}

