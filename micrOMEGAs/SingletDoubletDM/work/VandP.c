#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../CalcHEP_src/include/extern.h"
#include "../../CalcHEP_src/include/VandP.h"
#include "autoprot.h"
extern int  FError;
/*  Special model functions  */
extern int access(const char *pathname, int mode);

int nModelParticles=21;
static ModelPrtclsStr ModelPrtcls_[21]=
{
  {"a","a",1, 22, "0","0",2,1,2,0}
, {"Z","Z",1, 23, "MZ","WZ",2,1,3,0}
, {"W+","W-",0, 24, "MW","WW",2,1,3,3}
, {"g","g",1, 21, "0","0",2,8,16,0}
, {"ve","ve~",0, 12, "0","0",1,1,2,0}
, {"vm","vm~",0, 14, "0","0",1,1,2,0}
, {"vt","vt~",0, 16, "0","0",1,1,2,0}
, {"e-","e+",0, 11, "Me","0",1,1,2,-3}
, {"mu-","mu+",0, 13, "MMU","0",1,1,2,-3}
, {"ta-","ta+",0, 15, "MTA","0",1,1,2,-3}
, {"u","u~",0, 2, "MU","0",1,3,6,2}
, {"c","c~",0, 4, "MC","0",1,3,6,2}
, {"t","t~",0, 6, "MT","WT",1,3,6,2}
, {"d","d~",0, 1, "MD","0",1,3,6,-1}
, {"s","s~",0, 3, "MS","0",1,3,6,-1}
, {"b","b~",0, 5, "MB","0",1,3,6,-1}
, {"H","H",1, 25, "MH","WH",0,1,1,0}
, {"~N1","~N1",1, 9000005, "MN1","0",1,1,2,0}
, {"~N2","~N2",1, 9000006, "MN2","0",1,1,2,0}
, {"~N3","~N3",1, 9000007, "MN3","0",1,1,2,0}
, {"~Chi","~Chi~",0, 9000008, "MChi","0",1,1,2,3}
};
ModelPrtclsStr *ModelPrtcls=ModelPrtcls_; 
int nModelVars=42;
int nModelFunc=32;
static int nCurrentVars=41;
int*currentVarPtr=&nCurrentVars;
static char*varNames_[74]={
 "aEWM1","Gf","aS","ymdo","ymup","yms","ymc","ymb","ymt","yme"
,"ymm","ymtau","cabi","mS","mD","YChi","ThetaY","MN1","MN2","MN3"
,"U11","U21","U31","U12","U22","U32","U13","U23","U33","MZ"
,"Me","MMU","MTA","MU","MC","MT","MD","MS","MB","MH"
,"E","Pi","aEW","MW","sw2","EE","cw","sw","gw","g1"
,"vev","lam","muH","ye","ym","ytau","yup","yc","yt","ydo"
,"ys","yb","CKM1x1","CKM1x2","CKM1x3","CKM2x1","CKM2x2","CKM2x3","CKM3x1","CKM3x2"
,"CKM3x3","Y1","Y2","MChi"};
char**varNames=varNames_;
static REAL varValues_[74]={
   1.279000E+02,  1.166370E-05,  1.184000E-01,  5.040000E-03,  2.550000E-03,  1.010000E-01,  1.270000E+00,  4.700000E+00,  1.720000E+02,  5.110000E-04
,  1.056600E-01,  1.777000E+00,  2.277360E-01,  1.500000E+03,  1.000000E+03,  1.000000E+00,  7.800000E-01,  1.000000E+03, -1.012067E+03,  1.512067E+03
,  0.000000E+00, -7.071068E-01,  7.071068E-01, -6.914131E-02,  7.054146E-01,  7.054146E-01,  9.976069E-01,  4.889029E-02,  4.889029E-02,  9.118760E+01
,  5.110000E-04,  1.056600E-01,  1.777000E+00,  2.550000E-03,  1.270000E+00,  1.720000E+02,  5.040000E-03,  1.010000E-01,  4.700000E+00,  1.250000E+02
,  2.718282E+00,  3.141593E+00};
REAL*varValues=varValues_;
int calcMainFunc(void)
{
   int i;
   static REAL * VV=NULL;
   static int iQ=-1;
   static int cErr=1;
   REAL *V=varValues;
   FError=0;
   if(VV && cErr==0)
   { for(i=0;i<nModelVars;i++) if(i!=iQ && VV[i]!=V[i]) break;
     if(i==nModelVars)     return 0;
   }
  cErr=1;
   nCurrentVars=42;
   V[42]=Pow(V[0],-1);
   if(!isfinite(V[42]) || FError) return 42;
   nCurrentVars=43;
   V[43]=Pow(Pow(V[29],2)/(2.)+Pow(Pow(V[29],4)/(4.)-V[42]*V[41]*Pow(2,-0.5)*Pow(V[1],-1)*Pow(V[29],2),0.5),0.5);
   if(!isfinite(V[43]) || FError) return 43;
   nCurrentVars=44;
   V[44]=1-Pow(V[43],2)*Pow(V[29],-2);
   if(!isfinite(V[44]) || FError) return 44;
   nCurrentVars=45;
   V[45]=2*Pow(V[42],0.5)*Pow(V[41],0.5);
   if(!isfinite(V[45]) || FError) return 45;
   nCurrentVars=46;
   V[46]=Pow(1-V[44],0.5);
   if(!isfinite(V[46]) || FError) return 46;
   nCurrentVars=47;
   V[47]=Pow(V[44],0.5);
   if(!isfinite(V[47]) || FError) return 47;
   nCurrentVars=48;
   V[48]=V[45]*Pow(V[47],-1);
   if(!isfinite(V[48]) || FError) return 48;
   nCurrentVars=49;
   V[49]=V[45]*Pow(V[46],-1);
   if(!isfinite(V[49]) || FError) return 49;
   nCurrentVars=50;
   V[50]=2*V[43]*V[47]*Pow(V[45],-1);
   if(!isfinite(V[50]) || FError) return 50;
   nCurrentVars=51;
   V[51]=Pow(V[39],2)*Pow(V[50],-2)/(2.);
   if(!isfinite(V[51]) || FError) return 51;
   nCurrentVars=52;
   V[52]=Pow(V[51]*Pow(V[50],2),0.5);
   if(!isfinite(V[52]) || FError) return 52;
   nCurrentVars=53;
   V[53]=V[9]*Pow(2,0.5)*Pow(V[50],-1);
   if(!isfinite(V[53]) || FError) return 53;
   nCurrentVars=54;
   V[54]=V[10]*Pow(2,0.5)*Pow(V[50],-1);
   if(!isfinite(V[54]) || FError) return 54;
   nCurrentVars=55;
   V[55]=V[11]*Pow(2,0.5)*Pow(V[50],-1);
   if(!isfinite(V[55]) || FError) return 55;
   nCurrentVars=56;
   V[56]=V[4]*Pow(2,0.5)*Pow(V[50],-1);
   if(!isfinite(V[56]) || FError) return 56;
   nCurrentVars=57;
   V[57]=V[6]*Pow(2,0.5)*Pow(V[50],-1);
   if(!isfinite(V[57]) || FError) return 57;
   nCurrentVars=58;
   V[58]=V[8]*Pow(2,0.5)*Pow(V[50],-1);
   if(!isfinite(V[58]) || FError) return 58;
   nCurrentVars=59;
   V[59]=V[3]*Pow(2,0.5)*Pow(V[50],-1);
   if(!isfinite(V[59]) || FError) return 59;
   nCurrentVars=60;
   V[60]=V[5]*Pow(2,0.5)*Pow(V[50],-1);
   if(!isfinite(V[60]) || FError) return 60;
   nCurrentVars=61;
   V[61]=V[7]*Pow(2,0.5)*Pow(V[50],-1);
   if(!isfinite(V[61]) || FError) return 61;
   nCurrentVars=62;
   V[62]=Cos(V[12]);
   if(!isfinite(V[62]) || FError) return 62;
   nCurrentVars=63;
   V[63]=Sin(V[12]);
   if(!isfinite(V[63]) || FError) return 63;
   nCurrentVars=64;
   V[64]=0;

   nCurrentVars=65;
   V[65]=-Sin(V[12]);
   if(!isfinite(V[65]) || FError) return 65;
   nCurrentVars=66;
   V[66]=Cos(V[12]);
   if(!isfinite(V[66]) || FError) return 66;
   nCurrentVars=67;
   V[67]=0;

   nCurrentVars=68;
   V[68]=0;

   nCurrentVars=69;
   V[69]=0;

   nCurrentVars=70;
   V[70]=1;

   nCurrentVars=71;
   V[71]=V[15]*Cos(V[16]);
   if(!isfinite(V[71]) || FError) return 71;
   nCurrentVars=72;
   V[72]=V[15]*Sin(V[16]);
   if(!isfinite(V[72]) || FError) return 72;
   nCurrentVars=73;
   V[73]=V[14];

   if(VV==NULL) 
   {  VV=malloc(sizeof(REAL)*nModelVars);
      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) iQ=i;
   }
   for(i=0;i<nModelVars;i++) VV[i]=V[i];
   cErr=0;
   nCurrentVars++;
   return 0;
}
