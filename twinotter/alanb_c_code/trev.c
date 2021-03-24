#include <stdio.h>
#include <math.h>

#define A 17.269
#define B 35.86
#define C 610.78
#define D 273.16
#define K 17.67
#define KK 243.5
#define A0 6.107799961
#define A1 4.436518521e-1
#define A2 1.428945805e-2
#define A3 2.650648471e-4
#define A4 3.031240396e-6
#define A5 2.034080948e-8
#define A6 6.136820929e-11

/* calculates temperature t and liquid water content alwc from */
/* cloud base pressure po and temperature to, for adiabatic    */
/* ascent to the pressure p */
/* -> input: cloud base pressure po and temperature to *
/*            and pressure at observation level */
/* -> output: "adiabatic" temperature t and liquid water content alwc */
/* pres in N/m2, T in K */
/* return T in C and L in g/m3 */

float lf(),bolton(),teten(),vapor();

float trev(po,to,p,t,alwc)
float po,to,p;
float *t,*alwc;
{
  float eps,cpd,cw,rd,alhv;
  float tk,r,e,cpt,t1;
  float alref,tref,cpv;
  float rv,thetaq,tw;
  int i;

  eps=0.622;
  cpd=1.0057e3;
  cw=4.18e3;
  cpv=1.875e3;
  rd=287.04;
  alref=2.43e06;
  tref=303.15;

  tk=to;
  e=vapor(tk);
  r=eps*e/(po-e);
  cpt=cpd+r*cw;
/*  alhv=alref+(cpv-cw)*(tk-tref);*/
  alhv = (2.501 - 0.00237 * (tk - 273.15)) * 1.0e6;
  thetaq=tk*pow((1.e5/(po-e)),(rd/cpt))*exp(alhv*r/(cpt*tk));
/*  1st approx  */
  t1=tk;
  e=vapor(t1);
  rv=eps*e/(p-e);
  t1=thetaq/(pow((1.e5/(p-e)),(rd/cpt))*exp(alhv*rv/(cpt*t1)));
/*  successive approximations  */
  for (i=0;i<100;i++) {
    e=vapor(t1);
    rv=eps*e/(p-e);
/*    alhv=alref+(cpv-cw)*(t1-tref);*/
    alhv = (2.501 - 0.00237 * (t1 - 273.15)) * 1.0e6;
    t1=(thetaq/(pow((1.e5/(p-e)),(rd/cpt))*exp(alhv*rv/(cpt*t1)))+t1)/2.0;
  }
  *t = t1 - 273.15;
/*  get lwc  */
  e=vapor(t1);
  rv=eps*e/(p-e);
  tw=r-rv;
  *alwc=tw*p*28.9644/(8.314e7*t1)*1.e7;
  return(0);
}

/* goff-gratch formula for water saturation vapour pressure */
float vapor(tfp)
float tfp;
{
  float t,e,es,satr;

  t=tfp;
  e = -7.90298*(373.16/t-1.0)+5.02808*log10(373.16/t)
         -1.3816e-7*(pow(10.0,11.344*(1.0-t/373.16))-1.0)
         +8.1328e-3*(pow(10.0,3.49149*(1.0-373.16/t))-1.0);
  es=1013.246*pow(10.0,e);
  satr=100.*es;
  return(satr);
}

/* teten's formula for water saturation vapor pressure */

float teten(t)
float t;
{
  return(C*exp(A*(t - D)/(t - B)));
}

/* Bolton's formula for water saturation (Bolton, Mon. Weather Rev, 1980) */

float bolton(t)
float t;
{
  t -= 273.15;
  return(611.2*exp(K*t/(t+KK)));
}

/* Lowe and Ficke's formula for saturation (see Pruppacher & Klett, p625)*/

float lf(t)
float t;
{
  t -= 273.15;
  return(100.*(A0+t*(A1+t*(A2+t*(A3+t*(A4+t*(A5+A6*t))))))); 
}
