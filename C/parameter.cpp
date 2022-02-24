#include <cstdio>
#include <iostream>
#include <cstring>

#include "parameter.h"
#include "xxxlib.h"

using namespace std;

//read parameters file
void parm::ReadParm(char *filename){
    FILE *pFile;
    char buffer[512];  //  Buffer to store each line of the file
    char first_word[512];  //  First word of the current line
 
    double TMP;

//list of integer variables
    enum var_list_int{age,sex,stateofinfection,year,date,
     day,dateofinfection,dateofdruguse,dateofbitten,
     numberofbitten,drugname,eliminationorder,kwiatkowskinowak,
     andersonmaygupta};

//list of double variables
    enum var_list_double{weight,zeroorderrate,firstorderrate,initconc,
     kwiatkowskinowakr,kwiatkowskinowakd1,kwiatkowskinowakd2,
     kwiatkowskinowakp,kwiatkowskinowaks,kwiatkowskinowakx0,kwiatkowskinowaky0,
     kwiatkowskinowakd3,kwiatkowskinowakd4,kwiatkowskinowakh,kwiatkowskinowaks1,
     kwiatkowskinowaks2,kwiatkowskinowaks3,kwiatkowskinowaks4,kwiatkowskinowakx10,
     kwiatkowskinowakx20,kwiatkowskinowakx30,kwiatkowskinowakx40,andersonmayguptaLambda,
    andersonmayguptamu,andersonmayguptabeta,andersonmayguptaalpha,andersonmayguptadelta,
    andersonmayguptax0,andersonmayguptay0,andersonmayguptas0,andersonmayguptar};
   
    if((pFile = fopen(filename,"r"))==NULL){
        perror("ARE YOU KIDDING ME!? \nI CAN'T READ YOUR PARAMETER FILE.\n");
        exit(0);
    }

    //Keep reading in lines until we hit the EOF 
    cout<<"\nI'M READING YOUR FILE..."<<filename<<". PLEASE WAIT...\n";
//    wait(1);
    while(read_line(pFile,buffer) != -1){
        //Get the first word of the line                                                            
        find_first_word(buffer,first_word);

        if(strcasecmp(first_word,"age") == 0)
            add_int(buffer,age);

        if(strcasecmp(first_word,"sex") == 0)
            add_int(buffer,sex);

        if(strcasecmp(first_word,"stateofinfection") == 0)
            add_int(buffer,stateofinfection);

        if(strcasecmp(first_word,"year") == 0)
            add_int(buffer,year);

        if(strcasecmp(first_word,"date") == 0)
            add_int(buffer,date);

        if(strcasecmp(first_word,"day") == 0)
            add_int(buffer,day);

        if(strcasecmp(first_word,"weight") == 0)
            add_double(buffer,weight);

        if(strcasecmp(first_word,"dateofinfection") == 0)
            add_int(buffer,dateofinfection);

        if(strcasecmp(first_word,"dateofbitten") == 0)
            add_int(buffer,dateofbitten);

        if(strcasecmp(first_word,"numberofbitten") == 0)
            add_int(buffer,numberofbitten);


//read parameters of Drug part
        if(strcasecmp(first_word,"dateofdruguse") == 0)
            add_int(buffer,dateofdruguse);
        if(strcasecmp(first_word,"drugname") == 0)
            add_int(buffer,drugname);
        if(strcasecmp(first_word,"zeroorderrate") == 0)
            add_double(buffer,zeroorderrate);
        if(strcasecmp(first_word,"firstorderrate") == 0)
            add_double(buffer,firstorderrate);
        if(strcasecmp(first_word,"initconc") == 0)
            add_double(buffer,initconc);          
        if(strcasecmp(first_word,"eliminationorder") ==0)
            add_int(buffer,eliminationorder);

//read parasite parameters

//////KwiatkowskiNowak////////////
        if(strcasecmp(first_word,"kwiatkowskinowak")==0)
            add_int(buffer,kwiatkowskinowak);
        if(strcasecmp(first_word,"kwiatkowskinowakr")==0)
            add_double(buffer,kwiatkowskinowakr);
        if(strcasecmp(first_word,"kwiatkowskinowakd1")==0)
            add_double(buffer,kwiatkowskinowakd1);
        if(strcasecmp(first_word,"kwiatkowskinowakd2")==0)
            add_double(buffer,kwiatkowskinowakd2);          
        if(strcasecmp(first_word,"kwiatkowskinowakd3")==0)
            add_double(buffer,kwiatkowskinowakd3);          
        if(strcasecmp(first_word,"kwiatkowskinowakd4")==0)
            add_double(buffer,kwiatkowskinowakd4);          
        if(strcasecmp(first_word,"kwiatkowskinowakp")==0)
            add_double(buffer,kwiatkowskinowakp);           
        if(strcasecmp(first_word,"kwiatkowskinowaks")==0)
            add_double(buffer,kwiatkowskinowaks);                          
        if(strcasecmp(first_word,"kwiatkowskinowaks1")==0)
            add_double(buffer,kwiatkowskinowaks1);                          
        if(strcasecmp(first_word,"kwiatkowskinowaks2")==0)
            add_double(buffer,kwiatkowskinowaks2);                          
        if(strcasecmp(first_word,"kwiatkowskinowaks3")==0)
            add_double(buffer,kwiatkowskinowaks3);                          
        if(strcasecmp(first_word,"kwiatkowskinowaks4")==0)
            add_double(buffer,kwiatkowskinowaks4);                          
        if(strcasecmp(first_word,"kwiatkowskinowakh")==0)
            add_double(buffer,kwiatkowskinowakh);                          
        if(strcasecmp(first_word,"kwiatkowskinowakx0")==0)
            add_double(buffer,kwiatkowskinowakx0);                          
        if(strcasecmp(first_word,"kwiatkowskinowaky0")==0)
            add_double(buffer,kwiatkowskinowaky0);                          
        if(strcasecmp(first_word,"kwiatkowskinowakx10")==0)
            add_double(buffer,kwiatkowskinowakx10);                          
        if(strcasecmp(first_word,"kwiatkowskinowakx20")==0)
            add_double(buffer,kwiatkowskinowakx20);                          
        if(strcasecmp(first_word,"kwiatkowskinowakx30")==0)
            add_double(buffer,kwiatkowskinowakx30);                          
        if(strcasecmp(first_word,"kwiatkowskinowakx40")==0)
            add_double(buffer,kwiatkowskinowakx40);                          

////////AndersonMayGupta model/////////////////////////
    if(strcasecmp(first_word,"andersonmaygupta")==0)
        add_int(buffer,andersonmaygupta);                                                           
    if(strcasecmp(first_word,"andersonmayguptaLambda")==0)
        add_double(buffer,andersonmayguptaLambda);
    if(strcasecmp(first_word,"andersonmayguptamu")==0)
        add_double(buffer,andersonmayguptamu);                                                           
    if(strcasecmp(first_word,"andersonmayguptabeta")==0)
        add_double(buffer,andersonmayguptabeta);                                                           
    if(strcasecmp(first_word,"andersonmayguptaalpha")==0)
        add_double(buffer,andersonmayguptaalpha);                                                           
    if(strcasecmp(first_word,"andersonmayguptadelta")==0)
        add_double(buffer,andersonmayguptadelta);                                                           
    if(strcasecmp(first_word,"andersonmayguptax0")==0)
        add_double(buffer,andersonmayguptax0);                                                           
    if(strcasecmp(first_word,"andersonmayguptay0")==0)
        add_double(buffer,andersonmayguptay0);                                                           
    if(strcasecmp(first_word,"andersonmayguptas0")==0)
        add_double(buffer,andersonmayguptas0);                                                           
    if(strcasecmp(first_word,"andersonmayguptar")==0)
        add_double(buffer,andersonmayguptar);                                                           

////////////////////////////////

   }
    fclose(pFile);
}

void parm::ShowReadParm(void){
    cout<<"AGE: "<<Age<<"\n"
        <<"SEX: "<<Sex<<"\n"
        <<"WEIGHT: "<<Weight<<"\n" 
        <<"DATE OF DRUG USE: "<<DateOfDrugUse<<"\n"    
        <<"ELIMINATION ORDER: "<<EliminationOrder<<"\n"
        <<"FIRST ORDER RATE: "<<FirstOrderRate<<"\n"
        <<"INITIAL CONCENTRATION: "<<initConc<<"\n";
        
        
     
}

void parm::add_int(char *buf,int index){
    int read_count; //count from sscanf
    char *tmp;      //temp for the keywords
    enum var_list_int{age,sex,stateofinfection,year,date,
     day,dateofinfection,dateofdruguse,dateofbitten,
     numberofbitten,drugname,eliminationorder,kwiatkowskinowak,
     andersonmaygupta};
     
    if(index==age) 
        read_count = sscanf(buf,"%s %d",tmp,&Age);
    if(index==sex) 
        read_count = sscanf(buf,"%s %d",tmp,&Sex);
    if(index==stateofinfection) 
        read_count = sscanf(buf,"%s %d",tmp,&StateOfInfection);
    if(index==year) 
        read_count = sscanf(buf,"%s %d",tmp,&Year);
    if(index==date) 
        read_count = sscanf(buf,"%s %d",tmp,&Date);
    if(index==day) 
        read_count = sscanf(buf,"%s %d",tmp,&Day);
    if(index==dateofinfection) 
        read_count = sscanf(buf,"%s %d",tmp,&DateOfInfection);
    if(index==dateofdruguse) 
        read_count = sscanf(buf,"%s %d",tmp,&DateOfDrugUse);
    if(index==dateofbitten) 
        read_count = sscanf(buf,"%s %d",tmp,&DateOfBitten);
    if(index==numberofbitten) 
        read_count = sscanf(buf,"%s %d",tmp,&NumberOfBitten);
    if(index==drugname) 
        read_count = sscanf(buf,"%s %d",tmp,&DrugName);
    if(index==eliminationorder)
        read_count = sscanf(buf,"%s %d",tmp,&EliminationOrder);
    if(index==kwiatkowskinowak)
        read_count = sscanf(buf,"%s %d",tmp,&KwiatkowskiNowak);    
    if(index==andersonmaygupta)
        read_count = sscanf(buf,"%s %d",tmp,&AndersonMayGupta);
            
}

void parm::add_double(char *buf,int index){
    int read_count; //count from sscanf
    char *tmp;      //temp for the keywords

    enum var_list_double{weight,zeroorderrate,firstorderrate,initconc,
     kwiatkowskinowakr,kwiatkowskinowakd1,kwiatkowskinowakd2,
     kwiatkowskinowakp,kwiatkowskinowaks,kwiatkowskinowakx0,kwiatkowskinowaky0,
     kwiatkowskinowakd3,kwiatkowskinowakd4,kwiatkowskinowakh,kwiatkowskinowaks1,
     kwiatkowskinowaks2,kwiatkowskinowaks3,kwiatkowskinowaks4,kwiatkowskinowakx10,
     kwiatkowskinowakx20,kwiatkowskinowakx30,kwiatkowskinowakx40,andersonmayguptaLambda,
    andersonmayguptamu,andersonmayguptabeta,andersonmayguptaalpha,andersonmayguptadelta,
    andersonmayguptax0,andersonmayguptay0,andersonmayguptas0,andersonmayguptar};
    
    if(index==weight) 
        read_count = sscanf(buf,"%s %lf",tmp,&Weight);
    if(index==zeroorderrate)
        read_count = sscanf(buf,"%s %lf",tmp,&ZeroOrderRate);
    if(index==firstorderrate)
        read_count = sscanf(buf,"%s %lf",tmp,&FirstOrderRate);
    if(index==initconc)
        read_count = sscanf(buf,"%s %lf",tmp,&initConc);
    if(index==kwiatkowskinowakr)
        read_count = sscanf(buf,"%s %lf",tmp,&KwiatkowskiNowakr);      
    if(index==kwiatkowskinowakd1)
        read_count = sscanf(buf,"%s %lf",tmp,&KwiatkowskiNowakd1);
    if(index==kwiatkowskinowakd2)
        read_count = sscanf(buf,"%s %lf",tmp,&KwiatkowskiNowakd2);
    if(index==kwiatkowskinowakd3)
        read_count = sscanf(buf,"%s %lf",tmp,&KwiatkowskiNowakd3);
    if(index==kwiatkowskinowakd4)
        read_count = sscanf(buf,"%s %lf",tmp,&KwiatkowskiNowakd4);
    if(index==kwiatkowskinowakp)
        read_count = sscanf(buf,"%s %lf",tmp,&KwiatkowskiNowakp);
    if(index==kwiatkowskinowaks)
        read_count = sscanf(buf,"%s %lf",tmp,&KwiatkowskiNowaks);      
    if(index==kwiatkowskinowaks1)
        read_count = sscanf(buf,"%s %lf",tmp,&KwiatkowskiNowaks1);      
    if(index==kwiatkowskinowaks2)
        read_count = sscanf(buf,"%s %lf",tmp,&KwiatkowskiNowaks2);      
    if(index==kwiatkowskinowaks3)
        read_count = sscanf(buf,"%s %lf",tmp,&KwiatkowskiNowaks3);      
    if(index==kwiatkowskinowaks4)
        read_count = sscanf(buf,"%s %lf",tmp,&KwiatkowskiNowaks4);      
    if(index==kwiatkowskinowakh)
        read_count = sscanf(buf,"%s %lf",tmp,&KwiatkowskiNowakh);
    if(index==kwiatkowskinowakx0)
        read_count = sscanf(buf,"%s %lf",tmp,&KwiatkowskiNowakx0);
    if(index==kwiatkowskinowaky0)
        read_count = sscanf(buf,"%s %lf",tmp,&KwiatkowskiNowaky0);
    if(index==kwiatkowskinowakx10)
        read_count = sscanf(buf,"%s %lf",tmp,&KwiatkowskiNowakx10);
    if(index==kwiatkowskinowakx20)
        read_count = sscanf(buf,"%s %lf",tmp,&KwiatkowskiNowakx20);
    if(index==kwiatkowskinowakx30)
        read_count = sscanf(buf,"%s %lf",tmp,&KwiatkowskiNowakx30);
    if(index==kwiatkowskinowakx40)
        read_count = sscanf(buf,"%s %lf",tmp,&KwiatkowskiNowakx40);
    if(index==andersonmayguptaLambda)
        read_count = sscanf(buf,"%s %lf",tmp,&AndersonMayGuptaLambda);
    if(index==andersonmayguptamu)
        read_count = sscanf(buf,"%s %lf",tmp,&AndersonMayGuptamu);
    if(index==andersonmayguptabeta)
        read_count = sscanf(buf,"%s %lf",tmp,&AndersonMayGuptabeta);
    if(index==andersonmayguptaalpha)
        read_count = sscanf(buf,"%s %lf",tmp,&AndersonMayGuptaalpha);
    if(index==andersonmayguptadelta)
        read_count = sscanf(buf,"%s %lf",tmp,&AndersonMayGuptadelta);
    if(index==andersonmayguptax0)
        read_count = sscanf(buf,"%s %lf",tmp,&AndersonMayGuptax0);
    if(index==andersonmayguptay0)
        read_count = sscanf(buf,"%s %lf",tmp,&AndersonMayGuptay0);
    if(index==andersonmayguptas0)
        read_count = sscanf(buf,"%s %lf",tmp,&AndersonMayGuptas0);
    if(index==andersonmayguptar)
        read_count = sscanf(buf,"%s %lf",tmp,&AndersonMayGuptar);


}
