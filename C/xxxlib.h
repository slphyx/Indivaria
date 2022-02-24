#ifndef _XXX_
#define _XXX_

int read_line(FILE *fd, char *buf); // read in a line from a file 
void find_first_word(char *source, char *word); //find the first word of each line
void wait(float seconds); //wait for x seconds
void del_outputfile(void); //delete output files
void write_header(void); //write the header messages in the output files

#endif
