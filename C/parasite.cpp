#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "parasite.h"
#include "parameter.h"


void Parasite::Numbers_KwiatkowskiNowak2(double r,double d1,double d2,double p,double s,double x0,double y0,int day){
// this function calculates the number of parasites
//based on the methods from the work of Dominic Kwiatkowski and Martin Nowak
// Kwiatkowski, D. & Nowak, M.(1991). Periodic and chaotic host parasite interactions in human malaria.
// Proceedings of the National Academy of Sciences, USA 88, 5111-5113.
//
// x_t - the number of young parasites(0-24 hr postinvasion) on day t 
// y_t - the number of old parasites (24-48 hr postinvasion)
//     y_t+1 = d1 x_t f(x_t)
//     x_t+1 = r d2 y_t g(x_t)
//           r - represents the number of newly infected erythocytes that arises from a single
//               rupturing schizont.
//          d1,d2 - the probabilities that parasite survives the first and second half of its life cycle
//                  in the absence of fever.
//          f(x_t),g(x_t) - the fraction of young and old parasites that survive the damaging effect of the fever.
//
// from their paper, they chose f(x)=exp(-p*x) and g(x)=exp(-s*x)
// and the the range of human body temperature F(x) = 36.5+4(1-exp(-0.0003*x))
//
       
    double x_n, y_n, total, ret;
    double Fx; // the temperature
       
    y_n = d1 * x0 * exp(-1 * p * x0);
    x_n = r * d2 * y0 * exp(-1 * s * x0);                      
    
    KwiatkowskiNowakx0 = x_n; //initial value for the next step
    KwiatkowskiNowaky0 = y_n; //initial value for the next step
    
    total = y_n + x_n ;
    
    Fx = 36.5 + 4 * (1 - exp(-1 * 0.0003 * x_n));
    
    FILE *out;   
    if((out=fopen("paraKN.txt","a"))==NULL){
             perror("\nOH! I CAN'T WRITE THE DATA TO paraKN.txt\n");
             exit(0);
    }         
    
    ret = log10(total);
    
    fprintf(out, "%d\t %lf\t %lf\t \n",day,ret,Fx);
    printf("\n\tLog10 PARASITEMIA: %lf\n",ret);
    
    fclose(out);      
}
/////////////////////////////////////////////////////////

void Parasite::Numbers_KwiatkowskiNowak4(double r,double h,double d1,double d2,double d3,double d4,
         double s1,double s2,double s3,double s4,double x10,double x20,double x30,double x40,int day){
// this function calculates the number of parasites
//based on the methods from the work of Dominic Kwiatkowski and Martin Nowak
// Kwiatkowski, D. & Nowak, M.(1991). Periodic and chaotic host parasite interactions in human malaria.
// Proceedings of the National Academy of Sciences, USA 88, 5111-5113.
//
//the parasite population is divided into four stages of development:
//    x1 (0-12 hr postinvastion), x2 (12-24 hr postinvasion)
//    x3 (24-36 hr postinvasion), x4 (36-48 hr postinvasion)
//

  double x1t_05,x2t_05,x3t_05,x4t_05;
  double x1t,x2t,x3t,x4t;
  double total,Fx,ret;
  
  x1t_05 = (1-h)*r*d4*x40*exp(-1*s1*x10)+h*d1*x10*exp(-1*s1*x10);
  x2t_05 = (1-h)*d1*x10*exp(-1*s1*x10)+h*d2*x20*exp(-1*s2*x10);
  x3t_05 = (1-h)*d2*x20*exp(-1*s2*x10)+h*d3*x30*exp(-1*s3*x10);
  x4t_05 = (1-h)*d3*x30*exp(-1*s3*x10)+h*d4*x40*exp(-1*s4*x10);       

  x1t = (1-h)*r*d4*x4t_05*exp(-1*s1*x1t_05)+h*d1*x1t_05*exp(-1*s1*x1t_05);
  x2t = (1-h)*d1*x1t_05*exp(-1*s1*x1t_05)+h*d2*x2t_05*exp(-1*s2*x1t_05);
  x3t = (1-h)*d2*x2t_05*exp(-1*s2*x1t_05)+h*d3*x3t_05*exp(-1*s3*x1t_05);
  x4t = (1-h)*d3*x3t_05*exp(-1*s3*x1t_05)+h*d4*x4t_05*exp(-1*s4*x1t_05);       

    total = x1t + x2t + x3t + x4t;
    ret = log10(total);
//printf("\t%lf \t%lf \t%lf \t%lf \n",x10,x20,x30,x40);
//printf("\t%ef \t%lf \n",total,ret);

 
    KwiatkowskiNowakx10 = x1t; //initial value for the next step
    KwiatkowskiNowakx20 = x2t; //initial value for the next step
    KwiatkowskiNowakx30 = x3t; //initial value for the next step
    KwiatkowskiNowakx40 = x4t; //initial value for the next step
        
    Fx = 36.5 + 4 * (1 - exp(-1 * 0.0003 * x1t));
    
    FILE *out;
    
    if((out=fopen("paraKN.txt","a"))==NULL){
             perror("\nOH! I CAN'T WRITE THE DATA TO paraKN.txt\n");
             exit(0);
    }         
       
    fprintf(out, "%d\t %lf\t %lf\t \n",day,ret,Fx);
    printf("\n\tLog10 PARASITEMIA: %lf\n",ret);
    
    fclose(out);      
}
///////////////////////////////////////////////////////

void Parasite::Numbers_AndersonMayGupta(double x0,double y0,double s0,int day){
// this function does the integration of Anderson, May and Gupta's ODEs.
// Anderson,R. M.,May, R.M. & Gupta,S.(1989). Non-linear phenomena in host-parasite
// interactions. Parasitology 99 (Suppl.), S59-S79.
// This model attemps to address the blood stage asexual cycle of P.falciparum, by following
// the invasion of erythroctyes(RBC) by merozoites(M), and the resulting infected red blood cells(IRBC).
// These IRBCs rupture after 48 h, releasing merozoites into the bloodstream to cause further invasion.
// the equations are as below.
//        x' = Lambda - mu*x - beta*x*s
//        y' = beta*x*s - alpha*y;
//        s' = alpha*r*y - delta*s - beta*x*y;
//        
//   x - uninfected RBC concentration(per ul)
//   y - IRBC
//   s - merozoites        
//   Lambda - Bone-marrow erythropoiesis rate
//   mu - normal decay of RBCs rate (first-order process)
//   alpha - rate of merozoites released from rupturing IRBCs 
//   r - multiplication rate of merozoites released from rupturing IRBCs
//   beta - loss of infectivity
//   delta - culminating with death rate 

    double X,Y[3],Y1[3];
    double xend,h;
    
    int    finesse,i,k,kl,l;

    FILE *out;   

    if((out=fopen("paraAMG.txt","a"))==NULL){
             perror("\nOH! I CAN'T WRITE THE DATA TO paraAMG.txt\n");
             exit(0);
    }         
    
//  n=3;           // size of DE system 
//  xend=1.0;        // ending x
//  kl=100;         // number of calculated points
//  finesse=20;
//  h=(xl-x0)/kl/finesse;  // elementary integration step
//  h=xend/kl/finesse;
  h=1.0/2000;
  X=0;
  
  /*initial conditions*/
    Y[0] = x0;
    Y[1] = y0;
    Y[2] = s0;
    
  AndersonMayGuptaRK4(3,X,Y,h,Y1);
// fprintf(out,"%lf\t %lf\t %lf\t %lf\t\n", X,Y[0],Y[1],Y[2]);
  for (k=1; k<=2000; k++){
    X += h;
    for (i=0; i<3; i++) Y[i]=Y1[i];
    
    AndersonMayGuptaRK4(3,X,Y,h,Y1);
    
    if (k%100==0){   //Tabulate point
      fprintf(out,"%d\t %E\t %E\t %E\t\n",day,log10(Y[0]),log10(Y[1]),log10(Y[2]));
      printf("\n\t Log10 Uninfected RBCs: %E\t\n",log10(Y[0]));
    }
  }
  /*return the initial values*/
  AndersonMayGuptax0=Y[0]; AndersonMayGuptay0=Y[1]; AndersonMayGuptas0=Y[2];

  fclose(out);
                
}

//////////////////////////////////////////////////////
void Parasite::AndersonMayGuptaEquations(double x,double *y, double *yp){
//y[0]=x, y[1]=y, y[2]=s
  yp[0] = AndersonMayGuptaLambda-AndersonMayGuptamu*y[0]-AndersonMayGuptabeta*y[0]*y[2];
  yp[1] = AndersonMayGuptabeta*y[0]*y[2]-AndersonMayGuptaalpha*y[1];
  yp[2] = AndersonMayGuptaalpha*AndersonMayGuptar*y[1]-AndersonMayGuptadelta*y[2]-AndersonMayGuptabeta*y[0]*y[2];
     
}
///////////////////////////////////////////////////////
void Parasite::AndersonMayGuptaRK4(int n,double x,double *y,double h,double *y1){
    double k1[3],k2[3],k3[3],k4[3],yy[3],h2; 
    int i;
    
    AndersonMayGuptaEquations(x,y,k1);

    h2=h/2.0;
    
    for (i=0; i<n; i++) yy[i]=y[i]+h2*k1[i];
    
    AndersonMayGuptaEquations(x+h2,yy,k2);
    
    for (i=0; i<n; i++) yy[i]=y[i]+h2*k2[i];
    
    AndersonMayGuptaEquations(x+h2,yy,k3);
    
    for (i=0; i<n; i++) yy[i]=y[i]+h*k3[i];
    
    AndersonMayGuptaEquations(x+h,yy,k4);
    
    for (i=0; i<n; i++)
        y1[i]=y[i]+h*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])/6.0;
}
////////////////////////////////////////////////////////

/*void Parasite::GravenorRK4n(int n, double x, double *y, double h, double *y1) {
    double c1[SIZE],c2[SIZE],c3[SIZE],c4[SIZE],yy[SIZE],h2; 
    int i;
    
    GravenorEquations(x,y,c1);
     
    h2=h/2.0;
    
    for (i=0; i<n; i++) yy[i]=y[i]+h2*c1[i];
    GravenorEquations(x+h2,yy,c2);
    
    for (i=0; i<n; i++) yy[i]=y[i]+h2*c2[i];
    GravenorEquations(x+h2,yy,c3);
    
    for (i=0; i<n; i++) yy[i]=y[i]+h*c3[i];
    GravenorEquations(x+h,yy,c4);
    
    for (i=0; i<n; i++)
        y1[i]=y[i]+h*(c1[i]+2.0*c2[i]+2.0*c3[i]+c4[i])/6.0;

}
///////////////////////////////////////////////////////////

void Parasite::GravenorEquations(int n,double *y, double *yp){
    int i;

    yp[0] = GravenorR*GravenorLamda*y[n-1] - (GravenorLambda + Gravenormu[0])*y[0];
    
    for(i=1;i<n;i++)
        yp[i] = GravenorLambda*y[i-1] - (GravenorLambda + Gravenormu[i])*y[i];

}
////////////////////////////////////////////////////////

void Parasite::Numbers_Gravenor(double *){
          
     
}
 */
