#include <stdio.h>
#include <stdlib.h>


int main () {
   float x,y,z;
   FILE * fp;

   fp = fopen ("horizons_results_venus","r");
   printf("File was opened\n");
  
   fscanf(fp, "%E %E  %E", &x, &y, &z);
   printf("file was scanned\n");
   
   printf("Read String1 |%e|\n", x );
   printf("Read String2 |%e|\n", y );
   printf("Read String3 |%e|\n", z );

   fclose(fp);
   
   return(0);
}
