#ifndef _DRUG_
#define _DRUG_

#define Artemisinin     1
#define Quinine         2
#define Sulfadoxine     3
#define Pyrimethamine   4              

//char   *DrugNames[20]={"Artemisinin","Chloroquine","Mefloquine","Quinine",
//      "Primaquine","Sulfadoxine","Pyrimethamine","Proguanil",
//      "Hydroxychloroquine"};

class Drug
{
public:
      
      int EliminationOrder; // 0=zoro order,1=first order,2=second order
       
      double initConc; //initial concentration C(0)
      double ZeroOrderRate;
      double FirstOrderRate; //the first order elimination rate
      
      double PTP; //Post Treatment Prophylaxis duration 
      double MIC; //Minimum Inhibitory Concentration
      double MPC; //Minimum Parasiticidal Concentration
      
            
      double ConcentrationZeroOrder(int t);

      double ConcentrationFirstOrder(int t);

};

#endif
