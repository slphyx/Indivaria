#ifndef _MOSQUITO_
#define _MOSQUITO_

class mosquito
{
public:
    double SporozoiteRate(double a,double x,double P,double n);
    double InoculationRate(double m,double a,double b,double S);  
    double ReproductiveRate(double m,double a,double b,double x,double P,double n,double z); 
    double ReproductiveRate0(double m,double a,double b,double P,double n,double z);
    
};

#endif
