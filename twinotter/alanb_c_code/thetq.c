/* thetq */
/* usage: thetq amult < sounding_data                              */
/* compile:  cc -g -o thetq thetq.c -lcdf -lpgc -lm                */
/* march 11 83. modified may 85 to incorporate nstr & 2D criteria. */
/* purpose: to plot in-cloud and environmental thetaq.             */
/* inputs: sounding file, data file, time file are disc files.     */
/* data is high rate. Changed to C May 88.                         */

#include <stdio.h>   /* libraries to include for subroutines */
#include <math.h>
#include <strings.h>
#include "cdfhdr.h"
#include "pgraf.h"

#define MAX 1500

float sat();     /* goff-gratch formula for water saturation vapour pressure */
float thesq();   /* calculate saturated thetaq */
float theuq();   /* calculate undersaturated thetaq */
float clbas();   /* finds temp and pres for an undersaturated parcel */
float convert(); /* change HHMMSS to secs after midnight */

char buffer[HBMAX][LINE]; /* buffer for the cdf header */
float *sbuff;             /* buffer for the slice data */


FILE *dfile,*tfile,*fdate;  /* files */

main(argc,argv)
int argc;
char *argv[];
{
  float pmb[MAX],trf[MAX],sum;
  float ajw[MAX],qtot[MAX],amult;
  float thetaq[MAX],twodc[MAX];
  float qt[50],thetuq[50];
  float xclbt[2],yclbt[2];
  float xclbq[2],yclbq[2];
  float tt,td,pres,t,apres;
  float pc,tc,es,pd,ws,ql;
  float tbase,pbase,esbase,tqbase;
  float pdbase,wbase,qbase,q2base;
  float time,tstart,tstop;
  float time1,time2;
  float amax,amin;
  float tqdif,tqmin,tqmax;
  float tmax,tmin,qtmax,qtmin;
  float *conc2c,*atrf,*plwcf;
  float *strobe,*tptime,*psfdc;
  char lab1[120],lab2[120];
  char lab4[12],data[LINE];
  char timestr[120];
  char head[80],date[3][3];
  char dname[LINE],lwname[LINE];
  char realdate[LINE];
  long nel,stat;
  int press[50],nstr;
  int itqck,iqtck,dsize;
  int no,i,loop,loop1,indx;
  int flag1,flag2,level;
  struct field *fp;

/* set flags */
  flag1=flag2=0;
/* check input */
  if (argc != 3) {
    fprintf(stderr,"usage: thetq lwc_field amult < sounding\n");
    fprintf(stderr,"       lwc_field is plwcf, plwccz, lwccz, ... etc\n");
    exit(1);
  }
/* get the value of amult from the command line */
  amult = atof(argv[2]);
  strcpy(lwname,argv[1]);
/* open the graphics file */
  gopen();
  badset(BADLIM);
/* open the time file */
  if ((tfile=fopen("katim.d","r"))==NULL) {
    fprintf(stderr,"thetq: time file katim.d could not be opened\n");
    exit(1);
  }
/* read in the name to the data file */
  fgets(dname,LINE,tfile);
  level=strlen(dname);
  dname[level-1]='\0';
/* open the data file */
  if ((dfile=fopen(dname,"r"))==NULL) {
    fprintf(stderr,"thetq: datafile %s could not be opened\n",dname);
    exit(1);
  }
/* get the header from the file format */
  gethdr(dfile,buffer,HBMAX);
/* check the format, must be float */
  if (getfmt(buffer,HBMAX)!='f') {
    fprintf(stderr,"thetq: data file is not in float format!\n");
    exit(1);
  }
/* get date from header */
  sscanf(buffer[3],"%s %s %s",date[2],date[1],date[0]);

/*  strcpy(date,buffer[3]);
 * level=strlen(date);
 * date[level-1]='\0';
 */
/* read in static slice */
  nel=elemcnt(buffer,HBMAX,'s');
  sbuff=getbuff(nel);
  getslice(dfile,nel,sbuff);
/* set up memory for variable slice */
  free(sbuff);
  nel=elemcnt(buffer,HBMAX,'v');
  sbuff=getbuff(nel);
/* get pointers to the required data */
  tptime=getptr(buffer,HBMAX,sbuff,'v','d',"tptime");
  conc2c=getptr(buffer,HBMAX,sbuff,'v','d',"conc2c");
  atrf=getptr(buffer,HBMAX,sbuff,'v','d',"atrf");
/* note: inorder to change the liquid water content being read change the */
/* name in quotes. e.g. for plwccf change "plwcf" to "plwccf" */
  psfdc=getptr(buffer,HBMAX,sbuff,'v','d',"psfdc");
  plwcf=getptr2(buffer,HBMAX,sbuff,'v','d',lwname,&fp);
  dsize=fp->dsize1;
  strobe=getptr(buffer,HBMAX,sbuff,'v','d',"strobe");
/* read in the cloud base and sounding - header first. */
  fgets(head,LINE,stdin);
  fgets(data,LINE,stdin);
  sscanf(data,"%f %f",&tbase,&pbase);
/* continue */
  no=0;
  while (fgets(data,LINE,stdin)!=NULL) {
    sscanf(data,"%f %f %f",&pres,&tt,&td);
    press[no]=pres;
    pres=pres*100.0;
    tt=tt+273.15;
    td=td+273.15;
/* get tc & pc */
    clbas(pres,tt,td,&pc,&tc);
    es=sat(td);
    pd=pres-es;
    ws=.622*es/pd;
    qt[no]=ws*1000.0;
    thetuq[no++]=theuq(tt,tc,pd,pc,ws);
  }
/* put cloud base T & p in correct units */
  tbase=tbase+273.15;
  pbase=pbase*100.0;
/* calculate thetaq, qtot for these points */
  esbase=sat(tbase);
  pdbase=pbase-esbase;
  wbase=0.622*esbase/pdbase;
  qbase=wbase;
  q2base=qbase*1000.0;
  tqbase=thesq(tbase,pdbase,wbase,wbase);
  fprintf(stderr,"cb tq = %f\n",tqbase);
/* loop for the in-cloud points */
/* read in penetration times */
  while (fgets(data,LINE,tfile)!=NULL) {
    indx=0;
    sscanf(data,"%d %f %f",&level,&time1,&time2);
    if (time1 > time2) {
      fprintf(stderr,"thetq: tstart > tstop, closing files\n");
      gclose();
      exit(1);
    }
/* convert the times (given as HHMMSS) to secs after midnight */
    tstart=convert(time1-60000);
    tstop=convert(time2-60000);
/* find number of records to skip to beginning time segment */
    for (;;) {
      stat=getslice(dfile,nel,sbuff);
      if (stat == EOF) {
        fprintf(stderr,"thetq: unexpected End Of File. Closing pgraf.out\n");
        gclose();
        exit(1);
      }
      if (tptime[0] >= tstart-1.0)
        break;
    }
/* loop thru the records between the times records */
    loop1=0;
    while (loop1<1) {
      for (;;) {
        stat = getslice(dfile,nel,sbuff);
        if (stat == EOF) {
          fprintf(stderr,"thetq: unexpected End Of File. Closing pgraf.out\n");
          gclose();
          exit(1);
        }
        time=tptime[0];
        if (time >= tstart)
          break;
      }
/* record type 5 for j2 liquid water content */
      sum=0;
      for (loop=0;loop<dsize;loop++)
        sum = sum+plwcf[loop];
      ajw[indx]=sum*amult/dsize;
/* if cloud not hit yet then go back */
      if ((ajw[indx]<=0)&&(flag1==0))
        continue;
/* if cloud not hit yet and the time is past tstop+2.0 then get next  */
/* times from tfile. Tstop+2.0 is to ensure that indeed the cloud has */
/* been passed and to prevent the program from going any further into */
/* data with invalid times */
      if (time > tstop+2.0) {
/* set loop1 to a value greater than one */
        loop1=10;
        break;
      }
/* record type 0 for 2D C concentration. */
      sum=0;
      for (loop=0;loop<20;loop++)
        sum = sum+conc2c[loop];
      twodc[indx]=sum/20;
/* record type 3 for pressure. */
      sum=0;
      for (loop=0;loop<20;loop++)
        sum = sum+psfdc[loop];
/* convert pressure in mb to n/m2 */
      pmb[indx]=(sum/20)*100;
/* record type 9 for reverse flow temperature. */
      sum=0;
      for (loop=0;loop<20;loop++)
        sum = sum+atrf[loop];
      trf[indx]=sum/20.0+273.15;
      for (loop=0;loop<50;loop++) {
        nstr = (int) strobe[loop];
        if (nstr==0) {
          ajw[indx]=0.0;
          break;
        }
      }
/* set flag1=1 when into cloud for first time. */
      ++indx;
      if (indx==1)
        flag1=1;
/* back for more data else continue to go thru cloud */
      if (time >= tstop)
        break;
    }
/* check value of loop1, if greater than 1 then last times were invalid */
/* (i.e. the plane never encounterd an ajw value > 0) */
    if (loop1 > 1)
      continue;
/******** end of data loop  *********/
/* calculation for in-cloud points. */
/* first some units                 */
/* pressure in n/m2 = 100*p in mb   */
/* temp in k = 273.15+t in c        */
/* water vapour mix ratio in g/kg   */
/* l in g/kg = g/m3 *rt/p           */
/* q in g/kg = l+wat vap mix rat.   */
/* sat vapour pressure in n/m2      */
/* constants:                       */
/*          rd=287      j/kg/k      */
/*          cp=1.0046e3 j/kg/k      */
/*          cw=4.218e3  j/kg/k      */
/*          l =2.5e6    j/kg        */
/************************************/
    for (i=3;i<indx;i++) {
      if ((ajw[i-3]>0.0)&&(ajw[i-2]>0.0)&&(ajw[i-1]>0.0)&&
           (ajw[i]>0.0)&&(twodc[i]<=2.0)) {
        ql=ajw[i]*287.0*trf[i]/pmb[i];
/* water saturation vapour pressure */
        t=trf[i];
        es=sat(t);
/* dry air pressure */
        pd=pmb[i]-es;
/* saturation mixing ratio in kg/kg */
        ws=0.622*es/pd;
/* total water content mixing ratio in kg/kg */
        qtot[i]=ql*(1.e-3)+ws;
        thetaq[i]=thesq(t,pd,ws,qtot[i]);
        ws=ws*1.e3;
        qtot[i]=qtot[i]*1.e3;
      }
      else {
        thetaq[i]=BAD;
        qtot[i]=BAD;
      }
    }
/* average pressure in mb */
    sum=0;
    for (loop=0;loop<indx;loop++)
      sum=sum+pmb[loop];
    apres=sum/(100*indx);
/* plot here */
    gclear();
/* labels */
    window(1,0.0,15.0,9.0,10.0);
    sprintf(lab1," KA data collected: %s%s%s;  CB:",date[0],date[1],date[2]);
    sprintf(data,"%4.0f mb %4.1f C",pbase/100.0,tbase-273.15);
    strcat(lab1,data);
    sprintf(lab2," Time: %6.0f - %6.0f;  <p> = %4.0f mb",time1,
            time2,apres);
    label(1,.01,.7,lab1);
    label(1,.01,.3,lab2);

/* write lwc variable */
    sprintf(lab1,"lwc var: %s",lwname);
    label(1,.75,.7,lab1);
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
    sprintf(timestr,"thetq output as of %s",realdate); 
    label(1,.6,.3,timestr);

/* this section only on first penetration */
    ++flag2;
    if (flag2==1)  {
/* find limits for thetaq and Q from sounding values. */
      scale(no,5.0,thetuq,&amax,&amin);
      tqdif = amax - amin;
      tqmin=amin;
      if (tqdif < 15)
        tqmax = tqmin + 20.0;
      else if (tqdif < 20)
        tqmax = tqmin + 30.0;
      else
        tqmax = tqmin + 40.0;
      scale(no,5.,qt,&tmax,&tmin);
      qtmin = 0.0;
      if (tmax < 8)
        qtmax = 10.0;
      else if (tmax < 10)
        qtmax = 12.0;
      else
        qtmax = 16.0;
/* tick marks */
      itqck = (tqmax - tqmin) / 5 + 1;
      iqtck = qtmax / 2 + 1;
    }
/* axes */
    window(2,0.0,15.0,0.0,9.0);
    axes(2,tqmin,tqmax,itqck,"thetaq (k)",0.,qtmax,iqtck,"q (g/kg)");
/* sounding */
    line(2,1,no,thetuq,qt);
/* place pressure values next to points on sounding */
    for (i=0;i<no;i++) {
      sprintf(lab4,"%d",press[i]);
      label(2,thetuq[i]-2.0,qt[i],lab4);
    }
/* mark cloud base point. */
    mark(2,1,tqbase,q2base,1.2);
/* plot in-cloud points. */
    for (i=0;i<indx;i++)
      mark(2,1,thetaq[i],qtot[i],0.3);
    gpause();
/* zero thetaq and qtot. */
    for (i=0;i<indx;i++) {
      thetaq[i]=0.0;
      qtot[i]=0.0;
    }
/* return for another penetration. */
    flag1=0;
  }
/* end of time loop. close files */
  fclose(dfile);
  fclose(tfile);
  gclose();
}

/* goff-gratch formula for water saturation vapour pressure */
float sat(tfp)
float tfp;
{
  float t,e,es,satr;

  t=tfp;
  e = -7.90298*(373.16/t-1.0)+5.02808*log10(373.16/t)
         -1.3816e-7*(pow(10.0,11.344*(1.0-t/373.16))-1.0)
         +8.1328e-3*(pow(10.0,3.49149*(1.0-373.16/t))-1.0);
  es=1013.246*pow(10.0,e);
/* change to n/m2 */
  satr=es*1.e2;
  return(satr);
}

/* calculates saturated thetaq */
float thesq(tt,pp,ww,qq)
float tt,pp,ww,qq;
{
  float thesqr;

  thesqr = tt*pow((1.0e5/pp),(287.0/(1004.6+4218.0*qq)))
           *exp(ww*2.5e6/(tt*(1004.6+4218.0*qq)));
  return(thesqr);
}


/* calculates undersaturated thetaq */
float theuq(tt,tc,pp,pc,ww)
float tt,tc,pp,pc,ww;
{
  float theuqr;

  theuqr=tt*pow((1.0e5/pp),(.286*(1.0-.26*ww)))
      *exp(ww*2.5e6/(tc*(1004.6+4218.0*ww)))*pow((1.0e5/pc),(-1.12*ww));
  return(theuqr);
}

/* the temperature and pressure at the lifting condensation */
/* level is found for an undersaturated parcel with given   */
/* pressure (mb), dry temperature (k) and mixing ratio.     */
float clbas(p,t,td,pc,tc)
float p,t,td,*pc,*tc;
{
  float es,wc,ws;
  int i;

/* prepare iteration */
  *tc=t;
/* find wc (sat vap at td) */
  es=sat(td);
  wc=.622*es/(p-es);
/* iterate to find lcl */
/* that is, match the iterated vapour mixing ratio ws to wc */
  for (i=0;i<200;i++) {
    es=sat(( *tc));
    *pc=p*pow(((*tc)/t),(1.0046e3/287.0));
    ws=.622*es/((*pc)-es);
    *tc = *tc-(ws-wc)/(2.5e6*ws/(461.5*(*tc)*(*tc))-ws/(*tc));
    if ((ws > 0.9999*wc)&&(ws < 1.0001*wc))
      break;
  }
  return(OK);
}

float convert(t1)
float t1;
{
  int hh,mm,ss;  /* hours mins and secs */
  float sam;     /* sec after midnight  */
  

  hh = (int) t1/10000;
  mm = (int) (t1-hh*10000)/100;
  ss = (int) t1%100;
  sam = (float) hh*3600+mm*60+ss;
  return(sam);
}
