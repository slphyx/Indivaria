/////////////////////////////////////
// Indivalia
// Author : S. Saralamba
// Email: sompob@tropmedres.ac
/////////////////////////////////////

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cstring>

#include "drug.h"
#include "parameter.h"
#include "xxxlib.h"
#include "mosquito.h"
#include "parasite.h"



using namespace std;

int main(int argc, char *argv[])
{
    int i,tmp;
    Drug drug1;
    char *filename;
    parm SL;
    Parasite para1;
    
    FILE *parmfile;
    
    if(argc!=2){
        if((parmfile=fopen("indivaria.par","r"))!=NULL){
            filename="indivaria.par";
            fclose(parmfile); 
        }else{
              cout <<"\n====================Indivaria==================\n"
                   <<"version: 0.01"
                   <<"\nAuthor: S. Saralamba <sompob@tropmedres.ac>\n"
                   <<"=Mahidol-Oxford Tropical Medicine Research Unit=\n"
                   <<"Usage: "<<argv[0]<<" parameterfile\n"
                   <<"================================================\n";
              exit(0);
       }   
    }  
    if(strcasecmp(filename,"indivaria.par")!=0)
        filename=argv[1];     //name of parameter file         
        
    
    //read parameters from the input file
    SL.ReadParm(filename);
    SL.ShowReadParm();
    
    del_outputfile();  //delete the existed output files
    write_header();   //preparing the output files
    
    SL.StateOfInfection=0;

//read drug parameters
    drug1.EliminationOrder = SL.EliminationOrder;
    drug1.FirstOrderRate = SL.FirstOrderRate;
    drug1.ZeroOrderRate = SL.ZeroOrderRate;
    drug1.initConc = SL.initConc;


//read parasite parameters

///Kwiatkowski, Nowak model///
    para1.KwiatkowskiNowak = SL.KwiatkowskiNowak;
    if(para1.KwiatkowskiNowak==1){
        para1.KwiatkowskiNowakd1 = SL.KwiatkowskiNowakd1;
        para1.KwiatkowskiNowakd2 = SL.KwiatkowskiNowakd2;
        para1.KwiatkowskiNowakr = SL.KwiatkowskiNowakr;
        para1.KwiatkowskiNowakp = SL.KwiatkowskiNowakp;
        para1.KwiatkowskiNowaks = SL.KwiatkowskiNowaks;
        para1.KwiatkowskiNowakx0 = SL.KwiatkowskiNowakx0;
        para1.KwiatkowskiNowaky0 = SL.KwiatkowskiNowaky0;      
    }
    if(para1.KwiatkowskiNowak==2){
        para1.KwiatkowskiNowakd1 = SL.KwiatkowskiNowakd1;
        para1.KwiatkowskiNowakd2 = SL.KwiatkowskiNowakd2;
        para1.KwiatkowskiNowakd3 = SL.KwiatkowskiNowakd3;
        para1.KwiatkowskiNowakd4 = SL.KwiatkowskiNowakd4;
        para1.KwiatkowskiNowakh = SL.KwiatkowskiNowakh;
        para1.KwiatkowskiNowakr = SL.KwiatkowskiNowakr;
        para1.KwiatkowskiNowaks1 = SL.KwiatkowskiNowaks1;
        para1.KwiatkowskiNowaks2 = SL.KwiatkowskiNowaks2;
        para1.KwiatkowskiNowaks3 = SL.KwiatkowskiNowaks3;
        para1.KwiatkowskiNowaks4 = SL.KwiatkowskiNowaks4;
        para1.KwiatkowskiNowakx10 = SL.KwiatkowskiNowakx10;
        para1.KwiatkowskiNowakx20 = SL.KwiatkowskiNowakx20;          
        para1.KwiatkowskiNowakx30 = SL.KwiatkowskiNowakx30;
        para1.KwiatkowskiNowakx40 = SL.KwiatkowskiNowakx40;
    }

////AndersonMayGupta model//////
    para1.AndersonMayGupta = SL.AndersonMayGupta;
    if(para1.AndersonMayGupta==1){
        para1.AndersonMayGuptaLambda = SL.AndersonMayGuptaLambda;
        para1.AndersonMayGuptamu = SL.AndersonMayGuptamu;
        para1.AndersonMayGuptabeta = SL.AndersonMayGuptabeta;
        para1.AndersonMayGuptaalpha = SL.AndersonMayGuptaalpha;
        para1.AndersonMayGuptadelta = SL.AndersonMayGuptadelta;
        para1.AndersonMayGuptax0 = SL.AndersonMayGuptax0;
        para1.AndersonMayGuptay0 = SL.AndersonMayGuptay0;
        para1.AndersonMayGuptas0 = SL.AndersonMayGuptas0;
        para1.AndersonMayGuptar = SL.AndersonMayGuptar;                                                                
    }
    
    srand(time(NULL));  
    tmp=500;  
    
    double KNtmpX,KNtmpY;
    double KNtmpX1,KNtmpX2,KNtmpX3,KNtmpX4;
    KNtmpX = para1.KwiatkowskiNowakx0; KNtmpY = para1.KwiatkowskiNowaky0;  
    KNtmpX1 = para1.KwiatkowskiNowakx10; KNtmpX2 = para1.KwiatkowskiNowakx20;
    KNtmpX3 = para1.KwiatkowskiNowakx30; KNtmpX4 = para1.KwiatkowskiNowakx40;
    
    double AMGtmpX,AMGtmpY,AMGtmpS;
    AMGtmpX = para1.AndersonMayGuptax0; AMGtmpY = para1.AndersonMayGuptay0; 
    AMGtmpS = para1.AndersonMayGuptas0;
    
    for (i=1;i<=SL.Day;i++){
        cout<<"\nDAY: "<<i<<"\tBITTEN: "<<rand()%5+1<<"\tSTATUS: "<<SL.StateOfInfection<<"\n";

        if(i==SL.DateOfDrugUse||SL.StateOfInfection==1){
            SL.StateOfInfection=1;
            tmp=tmp-1;
 
            if(drug1.EliminationOrder==0)
                cout<<"\tDRUG CONCENTRATION: "<<drug1.ConcentrationZeroOrder(i);
 
            if(drug1.EliminationOrder==1)
                cout<<"\tDRUG CONCENTRATION: "<<drug1.ConcentrationFirstOrder(i);
            
            if(para1.KwiatkowskiNowak==1){
                para1.Numbers_KwiatkowskiNowak2(para1.KwiatkowskiNowakr,para1.KwiatkowskiNowakd1,para1.KwiatkowskiNowakd2,
                para1.KwiatkowskiNowakp,para1.KwiatkowskiNowaks,KNtmpX,KNtmpY,i);
                KNtmpX = para1.KwiatkowskiNowakx0; KNtmpY = para1.KwiatkowskiNowaky0;  
            }else if (para1.KwiatkowskiNowak==2){
                  para1.Numbers_KwiatkowskiNowak4(para1.KwiatkowskiNowakr,para1.KwiatkowskiNowakh,para1.KwiatkowskiNowakd1,
                  para1.KwiatkowskiNowakd2,para1.KwiatkowskiNowakd3,para1.KwiatkowskiNowakd4,para1.KwiatkowskiNowaks1,
                  para1.KwiatkowskiNowaks2,para1.KwiatkowskiNowaks3,para1.KwiatkowskiNowaks4,KNtmpX1,KNtmpX2,KNtmpX3,KNtmpX4,i);  
                  KNtmpX1 = para1.KwiatkowskiNowakx10; KNtmpX2 = para1.KwiatkowskiNowakx20;
                  KNtmpX3 = para1.KwiatkowskiNowakx30; KNtmpX4 = para1.KwiatkowskiNowakx40;
            }

            if(para1.AndersonMayGupta==1){
                para1.Numbers_AndersonMayGupta(AMGtmpX,AMGtmpY,AMGtmpS,i);
                AMGtmpX = para1.AndersonMayGuptax0; AMGtmpY = para1.AndersonMayGuptay0; AMGtmpS = para1.AndersonMayGuptas0;
                           
            }
            
            if(tmp==0){
                tmp=500;
                SL.StateOfInfection=0;
            }
        }   
        wait(0.1); 
    }

/*
    mosquito mos1;
    cout << "\nsporozoite rate:"<<mos1.SporozoiteRate(0.25,5.36,0.055,1)<<"\n inoculation rate:"
         << mos1.InoculationRate(2,3.5,0.5,mos1.SporozoiteRate(0.25,5.36,0.055,1))<<"\n";
    cout << "reproductive rate:"<<mos1.ReproductiveRate(1,2,3,4,5,6,7)<<"\n";     
//    system("PAUSE");
 */

    return EXIT_SUCCESS;
}
