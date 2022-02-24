//
// parasites class 
//
#ifndef _PARASITE_
#define _PARASITE_

class Parasite
{
      public:
             
//////////////Kwiatkowski Nowak model///////             
    int KwiatkowskiNowak; //0=disable;1=2stages,2=4stages;
    void Numbers_KwiatkowskiNowak2(double r,double d1,double d2,double p,double s,double x0,double y0,int day);
    void Numbers_KwiatkowskiNowak4(double r,double h,double d1,double d2,double d3,double d4,
         double s1,double s2,double s3,double s4,double x10,double x20,double x30,double x40,int day);         
    double KwiatkowskiNowakr;
    double KwiatkowskiNowakd1;
    double KwiatkowskiNowakd2;
    double KwiatkowskiNowakd3;
    double KwiatkowskiNowakd4;
    double KwiatkowskiNowakp;
    double KwiatkowskiNowaks;
    double KwiatkowskiNowakh;
    double KwiatkowskiNowakx0;
    double KwiatkowskiNowaky0;
    double KwiatkowskiNowaks1;       
    double KwiatkowskiNowaks2;       
    double KwiatkowskiNowaks3;       
    double KwiatkowskiNowaks4;       
    double KwiatkowskiNowakx10;
    double KwiatkowskiNowakx20;
    double KwiatkowskiNowakx30;
    double KwiatkowskiNowakx40;

///////////////////////////////////////////////////

//////////////////AndersonMayGupta////////////////
    int AndersonMayGupta; //0=disable;1=enable
    void Numbers_AndersonMayGupta(double x0,double y0,double s0,int day);
    void AndersonMayGuptaEquations(double x,double *y, double *yp);
    double AndersonMayGuptaLambda;
    double AndersonMayGuptamu;
    double AndersonMayGuptabeta;
    double AndersonMayGuptaalpha;
    double AndersonMayGuptadelta;
    double AndersonMayGuptax0;
    double AndersonMayGuptay0;
    double AndersonMayGuptas0;
    double AndersonMayGuptar;
    void AndersonMayGuptaRK4(int n, double x, double *y, double h, double *y1);
    
//////////////////////////////////////////////////

///////////////Gravenor et al.(JTB 2002) 217, 137-148/////////
//    int Gravenor;  //0=disable; 1=enable
//    void GravenorEquations(int n,double *y, double *yp);
//    void GravenorRK4n(int n, double x, double *y, double h, double *y1);
//    void Numbers_Gravenor();    
//    int GravenorCompartments;  //number of compartments
//    int GravenorR;             
//    double GravenorLambda;
//    double *Gravenormu;
//    double *Gravenor0;  //initial conditions  
    
//////////////////////////////////////////////////////////////

//////////////////S. Saralamba////////////////////////
    int initSPZ;  //initial number of sporozoites


/////////////////////////////////////////////////////    
                      
};



#endif
