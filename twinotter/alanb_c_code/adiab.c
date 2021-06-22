/* adiab - programs to write to stdout and plot adiabatic L and T */

#include <stdio.h>   /* libraries to include for subroutines */
#include <math.h>
#include <strings.h>
#include "cdfhdr.h"
#include "pgraf.h"

#define MAX 1000
#define T0    288.0
#define P0    1013.0
#define G     9.81
#define RD    287.0
#define ALR   0.0065

float trev();    /* subroutine to do the adiabatic work */
float ptz();

FILE *fdate;  /* files */

main(argc,argv)
int argc;
char *argv[];
{
  float press[MAX],temp[MAX],la[MAX];
  float al,at,z;
  float tbase,pbase;
  float tmax,tmin;
  float pmin;
  char lab1[120];
  char timestr[120];
  char data[LINE],head[LINE];
  char realdate[LINE];
  int i,kk;
  int ttck,ptck,pm;

/* check input */
  if (argc != 1) {
    fprintf(stderr,"usage: adiab < sounding > print_file\n");
    exit(1);
  }
/* open the graphics file */
  gopen();
  badset(BADLIM);

/* read in the cloud base from the sounding - header first. */
  fgets(head,LINE,stdin);
  fgets(data,LINE,stdin);
  sscanf(data,"%f %f",&tbase,&pbase);

/* change units */
  tbase += 273.16;
  pbase *= 100.;
/* loop on pressure */
  for (i=0;i<150;i++) {
    press[i] = (float) (pbase - 100. * (i * 5));
    trev(pbase,tbase,press[i],&at,&al);
    z = ptz(press[i]/100.);
    fprintf(stdout,"alt = %5.1f m, pres = %5.1f mb, La = %6.2f g/m3, Ta = %6.1f\n",
                   z,press[i]/100.,al,at);
    temp[i] = at;
    la[i] = al;
    press[i] /= 100.;
  } 

/* plot here */
  gclear();
  window(1,3.0,12.0,9.0,10.0);
  window(2,0.0,7.5,0.0,9.0);
  window(3,7.5,15.0,0.0,9.0);
  pm = (int) pbase / 10000 * 100;
  pmin = (float) pm + 100.;
  ptck = (int) (pmin - 300) / 100 + 1; 
  scale(100,10.,temp,&tmax,&tmin);
  ttck = (int) (tmax - tmin) / 10. + 1;
  axes(2,0.0,5.0,6,"La (g/m3)",pmin,300.,ptck,"pressure (mb)");
  axes(3,tmin,tmax,ttck,"Ta (C)",pmin,300.,ptck,"pressure (mb)");
/* labels */
  sprintf(lab1,"Cloud base:  %4.0f mb %4.1f C",pbase/100.0,tbase-273.15);
  label(1,.1,.7,lab1);

/*  write current date and time */
  system("date '+%h %d 19%y' > tempdate");
  if ((fdate = fopen("tempdate","r")) == NULL) {
    fprintf(stderr,"can't open tempdate\n");
    exit(0);
  }
  fgets(realdate,LINE,fdate);
  fclose(fdate);
  realdate[strlen(realdate)-1] = 0;
  system("rm tempdate");
  sprintf(timestr,"adiab output as of %s",realdate); 
  label(1,.1,.3,timestr);

/* plot adiabatic lwc and temperature */
  line(2,1,100,la,press);
  line(3,1,100,temp,press);
  gpause();
  gclose();
}

/* function to calculate altitude for a given pressure */
/* calculation uses a constant lapse rate of 6.5 C/km (see Hess pp 82-83) */
/* constants are defined at the top of this program */
float ptz(p)
float p;
{
  float ex,ptzr;

  ex = (RD*ALR)/G;
  ptzr = T0*(1.0-pow((p/P0),ex))/ALR;
  return(ptzr);
}
