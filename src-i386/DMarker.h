#include "stdio.h"
#include "time.h"

int init_dmarker( char **rfile, char **rBG, char **rBP, char **rUG, char **rUP, int *sBG, int *sBP, int *sUG, int *sUP, int *k1, int *k2 )
{
  FILE *P;
  int i = 0;
  time_t ti;
  time(&ti);
  P = fopen(rfile[0], "w");
  fprintf(P, "%s\nThe results of prediction is writing to the file : %s .\nDMarker will predict the", ctime(&ti), rfile[0]);
  if(*k2 == 1)  fprintf(P, " GENE ");
  if(*k2 == 2)  fprintf(P, " PROTEIN ");
  fprintf(P, "whether can be detected");
  if(*k1 == 1)  fprintf(P, " both in blood and urine .\n");
  if(*k1 == 2)  fprintf(P, " in blood .\n");
  if(*k1 == 3)  fprintf(P, " in urine .\n");
  fprintf(P, "\n");
  if( (*k1 == 1) || (*k1 == 2) )
  {
    if(*sBG  == 0)  fprintf(P, "No GENE or PROTEIN can be detected in BLOOD.\n");
    if(*sBG > 0)
    {
       fprintf(P, "%d GENEs and %d PROTEINs can be detected in BLOOD.\n[%d GENE]:\t", *sBG, *sBP, *sBG);
       for(i = 0; i < *sBG; i++)  fprintf(P, "%s\t", rBG[i]);
       fprintf(P, "\n[%d PROTEIN]:\t", *sBP);
       for(i = 0; i < *sBP; i++)  fprintf(P, "%s\t", rBP[i]);
       fprintf(P, "\n");
    }
  }
  if( (*k1 == 1) || (*k1 == 3) )
  {
    if(*sUG  == 0)  fprintf(P, "No GENE or PROTEIN can be detected in URINE.\n");
    if(*sUG > 0)
    {
       fprintf(P, "%d GENEs and %d PROTEINs can be detected in URINE.\n[%d GENE]:\t", *sUG, *sUP, *sUG);
       for(i = 0; i < *sUG; i++)  fprintf(P, "%s\t", rUG[i]);
       fprintf(P, "\n[%d PROTEIN]:\t", *sUP);
       for(i = 0; i < *sUP; i++)  fprintf(P, "%s\t", rUP[i]);
       fprintf(P, "\n");
    }
  }
  fclose(P);
  return 1;
}
