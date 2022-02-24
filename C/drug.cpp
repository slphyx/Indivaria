#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "drug.h"


//Drug.DrugName = {"Artemisinin","Chloroquine","Mefloquine",
//              "Quinine","Primaquine","Sulfadoxine","Pyrimethamine","Proguanil",
//              "Hydroxychloroquine"};

double Drug::ConcentrationFirstOrder(int t)
//calculate the drug concentration at the given time
//the formula for the 1st order process: C(t)=C(0)exp(-k t)
//    k     - the first order elimination rate
//    t     - time (days)
//    C(t)  - drug concentration at time t
//    C(0)  - initial drug concentration
{
      int DayToHours=24;  
      double Conc;        

      FILE *out;
    
      if((out=fopen("DrugConc.txt","a"))==NULL){
             perror("\nOH! I CAN'T WRITE THE DATA TO DrugConc.txt\n");
             exit(0);
      }    

      Conc = initConc*exp(-1*FirstOrderRate*t*DayToHours);    

      fprintf(out,"%d\t%lf\n",t,Conc);
      fclose(out);

      return Conc;
}

double Drug::ConcentrationZeroOrder(int t)
//the zero order process of drug elimination
//    rate = constant
{
       return(ZeroOrderRate);            
}

