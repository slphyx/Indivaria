#include <cstdio>
#include <cctype>
#include <ctime>
#include <cstdlib>

//find the first word from a line
void find_first_word(char *source, char *word)
{
	int i=0;	//  Position within source
	int j=0;	//  Position within word

	/*  Skip leading white space					*/
	while((source[i] != '\0') && isspace(source[i]))
		i++;

	/*  Copy the word						*/
	while((source[i] != '\0') && !isspace(source[i]))
	{
		word[j]=source[i];
		i++;
		j++;
	}

	word[j]='\0';

	return;
}

//read in a line from a file
int read_line(FILE *fd, char *buf)
{
	int i=0;	//  Current position in buf
	int c;		//  Character read in from file

	/*  Loop and read characters until we get either an EOF or a    */
	/*  newline							*/
	while ( ((c=fgetc(fd)) != EOF) && (c != '\n') )
	{
		/*  If we encounter a bracketed comment, skip it.  This */
		/*  basically means read EVERYTHING until the next } and*/
		/*  throw it into the big bit bucket			*/
		if (c == '{')
		{
			while ( ((c=fgetc(fd)) != EOF) && (c!='}') )
			{
			}
			if (c==EOF)
			{
				char err_msg[512];

				sprintf(err_msg, "ABNORMAL EOF FOUND - buffer=*%s*\n", buf);

			}
			continue;
		}

		/*  Also, use this little piece of logic to remove      */
		/*  any leading spaces from the line			*/
		if ((i>0) || !isspace(c))
		{
			buf[i] = c;
	
			i++;
		}
	}

	/*  NULL terminate the string					*/
	buf[i]='\0';

	/*  Check for an EOF in the middle of a line			*/
	if ((c==EOF) && (i!=0))
	{
		char err_msg[512];

		sprintf(err_msg, "ABNORMAL EOF FOUND - buffer=*%s*\n", buf);
		
	}

	/*  Return the appropriate value				*/
	if (c==EOF)
		return(-1);
	else
		return(0);
}

//wait for x seconds
void wait(float seconds)
{
  clock_t endwait;
  endwait = clock () + seconds * CLOCKS_PER_SEC ;
  while(clock() < endwait){} //loop !!
}

//writhe the header messages in the output files
void write_header(void)
{       
    FILE *out;   
    if((out=fopen("DrugConc.txt","a"))==NULL){
        perror("\nOH! I CAN'T WRITE THE HEADER FILE TO DrugConc.txt\n");
        exit(0);
    }  
    fprintf(out,"::::INDIVARIA::::\n");
    fprintf(out,"DRUG CONCENTRATION\n");
    fprintf(out,"DATE\tCONCENTRATION\n");
    
    fclose(out);

/* ******************************* */
    if((out=fopen("paraKN.txt","a"))==NULL){
        perror("\nOH! I CAN'T WRITE THE HEADER FILE TO paraKN.txt\n");
        exit(0);
    }
    fprintf(out,"::::INDIVARIA::::\n");
    fprintf(out,"TOTAL PARASITE\n");
    fprintf(out,"DATE\t PARASITE\t TEMPERATURE\n");
     
     
    fclose(out); 

/* ******************************* */
    if((out=fopen("paraAMG.txt","a"))==NULL){
        perror("\nOH! I CAN'T WRITE THE HEADER FILE TO paraAMG.txt\n");
        exit(0);
    }
    fprintf(out,"::::INDIVARIA::::\n");
    fprintf(out,"RED BLOOD CELLS\n");
    fprintf(out,"DATE\t UninfectedRBCs\t InfectedRBCs\t Merozoites\n");
     
    fclose(out); 
     
}



void del_outputfile(void){

     int number_of_files = 3,
         i = 0;

     if(remove("DrugConc.txt")!=0)
         perror("Error deleting file DrugConc.txt.");
     else{
         puts("DrugConc.txt is ready to be written."); i=i+1;
         }
     
     if(remove("paraKN.txt")!=0)
         perror("Error deleting file paraKN.txt.");
     else{
         puts("paraKN.txt is ready to be written."); i=i+1;
         }
      
     if(remove("paraAMG.txt")!=0)
         perror("Error deleting file paraAMG.txt.");
     else{
         puts("paraAMG.txt is ready to be written."); i=i+1;
         }
     
     if(i==number_of_files)   
         puts("Output files are ready to be written!\n");
          
}
