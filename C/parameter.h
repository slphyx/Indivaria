#ifndef _PARAMETERS_
#define _PARAMETERS_


typedef struct parm {

//variables for general information        
    int Sex; // M=Male=1, F=Femal=2, M2F=3, F2M=4
    int Age; // >=1 years
    int StateOfInfection; //0=susceptible, 1=infected
    double Weight;        
      
    int Year;
    int Date;             //from 1 to 365
    int Day;              //number of days for running the model
 
    int DateOfInfection;  //date of infection
    int DateOfDrugUse;    //date of drug use
    int DateOfBitten;     //date that the host is bitten
    int NumberOfBitten;   //bitten number in a day

 //variables for drug part
    int DrugName;          //the name of drug used in treatment
    int EliminationOrder;
    
    double initConc;
    double ZeroOrderRate;
    double FirstOrderRate;    
    

//variables for parasite part

/////KwiatkowskiNowak/////
    int KwiatkowskiNowak; //0 = disable; 1=2stages; 2=4stages
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
/////////////////////////

////Anderson May Gupta model/////
    int AndersonMayGupta;
    double AndersonMayGuptaLambda;
    double AndersonMayGuptamu;
    double AndersonMayGuptabeta;
    double AndersonMayGuptaalpha;
    double AndersonMayGuptadelta;
    double AndersonMayGuptax0;
    double AndersonMayGuptay0;
    double AndersonMayGuptas0;
    double AndersonMayGuptar;  
//////////////////////////////

    void ReadParm(char *filename);
    void ShowReadParm(void);      //use this function after Readparm!!                  

    void add_double(char *buf,int index); // store a value to a double type variable
    void add_int(char *buf,int index); // store a value to an integer type variable
      
};



#endif
