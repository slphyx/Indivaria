#include <fstream>

class Data
{
public:
       void Read(char *FileName);
       
};

void Data::Read(char *FileName)
{
     ifstream InData;
    
     InData.open("indivaria.par"); // opens the file
     if(!InData) { // if file couldn't be opened
          cerr << "Error: indivaria.par could not be opened" << endl;
          exit(1);
     }

}
