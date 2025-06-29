#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include "../../include/VandP.h"
#include "../../include/VandP_size.h"
#include"SLHAplus.h"
#include"../ntools/include/vegas.h"
#include"../ntools/include/1d_integration.h"
#include"../num/include/alphas2.h"

#include"dynamic_cs.h"
#include "vp.h"

#define VVmassGap  70


extern int displayPlot(char * title, char*xName, double xMin, double xMax ,  int lScale, int N, ...);

extern  char * trim(char *);

int ForceUG=0;
int useSLHAwidth=1;
decayTableStr* decayTable=NULL;

/*=============   decayPcm   and decayPcmW ================*/      

REAL  decayPcm(REAL am0,  REAL  am1,  REAL am2)
{
  REAL  summ, diffm;
  summ = am1 + am2;
  diffm = am1 - am2;
  if(am0<summ) return 0;
  return Sqrt((am0-summ)*(am0+ summ)*(am0-diffm)*(am0+diffm))/(am0*2);
}

static double ME(double m0,double m1,double m2)
{ 
  return  m2*m2*m1*m1 +(m0*m0 -m2*m2 -m1*m1)*(m0*m0 -m2*m2 -m1*m1)/8;
}


typedef struct{  double w1,w2,m0,m1,m2,m1p,m2p, y11,y12,y21,y22;  int nGauss;} sVarW;


#define w1_      sVar->w1 
#define w2_      sVar->w2
#define m0_      sVar->m0
#define m1_      sVar->m1
#define m2_      sVar->m2
#define m1_p     sVar->m1p
#define m2_p     sVar->m2p
#define y_11      sVar->y11
#define y_12      sVar->y12
#define y_21      sVar->y21
#define y_22      sVar->y22
#define n_Gauss  sVar->nGauss
#define nn       3

//static double w1_,w2_,m0_,m1_,m2_,m1_p,m2_p;
//static int n_Gauss=0;
//static double y_11,y_12,y_21,y_22;

static double intDecay2(double y2,sVarW*sVar)
{ double m2;
  m2= m2_p*m2_p + w2_*m2_p*tan(y2);
  if(m2<=0) return 0;
  m2=sqrt(m2);
  return  ME(m0_,m1_,m2)*decayPcm(m0_,m1_,m2);
}   

static double intDecay2_(double x,void*sVarVoid)
{ if(x<=0 || x>=1) return 0;
  sVarW*sVar=(sVarW*)sVarVoid;
  return intDecay2(y_21+pow(x,nn)*(nn+1-nn*x)*(y_22-y_21),sVar)*nn*(nn+1)*pow(x,nn-1)*(1-x)*(y_22-y_21);
}

static double  intDecay1(double y1,sVarW*sVar)
{ 
   m1_=m1_p*m1_p + w1_*m1_p*tan(y1);
   if(m1_<0) return 0;
   m1_=sqrt(m1_);
   
   y_21=atan(-m2_p/w2_);
   y_22=atan( ((m0_-m1_)*(m0_-m1_)-m2_p*m2_p)/(m2_p*w2_));
   if(n_Gauss) return  gauss_arg(intDecay2_,sVar,0,1,n_Gauss);
   else       return  simpson_arg(intDecay2_,sVar,0,1,1.E-5,NULL);
//   return  simpson(intDecay2, atan(-m2_p/w2_), atan( ((m0_-m1_)*(m0_-m1_)-m2_p*m2_p)/(m2_p*w2_)), 1.E-3);
}

static double intDecay1_(double x,void*sVarVoid)
{ //if(x<=0 || x>=1) return 0;
//   return intDecay1(y_11+pow(x,nn)*(1+nn*(1-x))*(y_12-y_11),sVar)*nn*(nn+1)*pow(x,nn-1)*(1-x)*(y_12-y_11);
   sVarW*sVar=(sVarW*)sVarVoid;
   double res= intDecay1(y_11+pow(x,nn)*(y_12-y_11),sVar)*nn*pow(x,nn-1)*(y_12-y_11);
//   printf("x =%e res=%e\n", x,res);

   return res;
}
  

double decayPcmW(double m0,double m1,double m2,double w1,double w2, int N)
{  

  sVarW sVarStr;
  sVarW*sVar=&sVarStr;

  n_Gauss=N;
  m0_=m0;
  if(w1==0 && w2==0) return decayPcm(m0,m1,m2);
  else if(w1==0)
  { if(m1>m0) return 0;
    m1_=m1;
    m2_p=m2;
    w2_=w2;
    y_21=atan(-m2/w2);
    y_22=atan( ((m0-m1)*(m0-m1)  -m2*m2)/(m2*w2));
    if(n_Gauss) return gauss_arg( intDecay2_,sVar,0,1,n_Gauss)/M_PI/ME(m0,m1,m2);
     else      return simpson_arg( intDecay2_,sVar,0,1,1.E-3,NULL)/M_PI/ME(m0,m1,m2); 
  }
  else if(w2==0)
  { if(m2>m0) return 0;
    m1_=m2;
    m2_p=m1;
    w2_=w1;
    y_21=atan(-m1/w1);
    y_22=atan( ((m0-m2)*(m0-m2)-m1*m1)/(m1*w1));
    if(n_Gauss) return gauss_arg( intDecay2_,sVar,0,1,n_Gauss)/M_PI/ME(m0,m1,m2);  
    else       return simpson_arg(intDecay2_,sVar,0,1,1.E-3,NULL)/M_PI/ME(m0,m1,m2);
  }
  else 
  {
    if(w1<=w2)  { w1_=w1; w2_=w2;  m1_p=m1;  m2_p=m2; }
     else       { w1_=w2; w2_=w1;  m1_p=m2;  m2_p=m1; }
   y_11=atan(-m1/w1); y_12=atan( (m0*m0-m1*m1)/(m1*w1));

//((m0-m2)^2-m1_p*m1_p)/(w1_*m1_p) =tan(y1);

// m1_=m1_p*m1_p + w1_*m1_p*tan(y1);
 
 
//double y_1m=atan( ((m0-m2+3*w2 )*(m0-m2+3*w2)-m1*m1)/(m1*w1));

//printf("y_11 y_1m y_12 = %E %E %E\n", y_11, y_1m, y_12);
//displayPlot("intDecay1_1","y",y_11,y_1m,0,1,"",0,intDecay1,sVar);
//displayPlot("intDecay1_2","y",y_1m,y_12,0,1,"",0,intDecay1,sVar);

//y_11=atan( ((m0-m2+3*w2 )*(m0-m2)-m1*m1)/(m1*w1));

//y_11+pow(x,nn)*(y_12-y_11)=y_1m 

//double x0=pow( (y_1m-y_11)/(y_12-y_11), 1./nn); 

//printf("x0=%e\n",x0);
//displayPlot("intDecay1_","x",0,1,0,1,"",0,intDecay1_,sVar);
//double sum=0;
//for(int i=0;i<10;i++) sum+=intDecay1_((i+0.5)/10.,sVar)/10;
//printf("=====================\n");
double res=simpson_arg(intDecay1_,sVar,0,1,1E-5,NULL);
//printf("sum=%E res=%E\n",sum,res);
return res/(M_PI*M_PI)/ME(m0,m1,m2);

  }
}

#undef w1_    
#undef w2_    
#undef m0_    
#undef m1_    
#undef m2_    
#undef m1_p   
#undef m2_p   
#undef y11    
#undef y12    
#undef y21    
#undef y22    
#undef n_Gauss


numout* xVtoxll(int Nin,int Nout,char**name,int *pdg, int lV, double *wV,  double *br)
{
  int i;   
  char* e_=NULL,*E_=NULL,*ne_=NULL,*Ne_=NULL,*m_=NULL,*M_=NULL,*nm_=NULL,*Nm_=NULL,*W_=NULL,*Z_=NULL;
  char processX3[50],plib13[50],exclude[20];
  int lV_=2*Nin+1-lV;
  char *c;
  numout*ccx3;
  
  double ww,wz,wBrE,wBrM,zBrEn,zBrMn;
   
  if(pdg[lV]!=23 &&  abs(pdg[lV])!=24) return NULL;
  
//if(Nin==2)printf("%s %s -> %s %s\n", name[0], name[1],name[2],name[3]);  

  for(i=0;i<nModelParticles &&!(e_&&ne_&&m_&&nm_&&e_&&Ne_&&M_&&Nm_&&W_&&Z_ )   ;i++) 
  switch(ModelPrtcls[i].NPDG)
  {  
       case  11: e_ =ModelPrtcls[i].name;  E_=ModelPrtcls[i].aname; break;
       case -11: e_ =ModelPrtcls[i].aname; E_=ModelPrtcls[i].name;  break;
       case  12: ne_=ModelPrtcls[i].name; Ne_=ModelPrtcls[i].aname; break;
       case -12: ne_=ModelPrtcls[i].aname;Ne_=ModelPrtcls[i].name;  break;
       case  13: m_ =ModelPrtcls[i].name;  M_=ModelPrtcls[i].aname; break;
       case -13: m_ =ModelPrtcls[i].aname; M_=ModelPrtcls[i].name;  break;
       case  14: nm_=ModelPrtcls[i].name; Nm_=ModelPrtcls[i].aname; break;
       case -14: nm_=ModelPrtcls[i].aname;Nm_=ModelPrtcls[i].name;  break;
       case  24: W_=ModelPrtcls[i].name;                            break;
       case -24:                           W_=ModelPrtcls[i].aname; break;
       case  23: Z_=ModelPrtcls[i].name;                            break;       
  }
    
  if(!(e_&&ne_&&m_&&nm_&&e_&&Ne_&&M_&&Nm_&&W_&&Z_)) return NULL;
   
  {  txtList wDlist,zDlist;
       char txt[20]; 
       ww=pWidth(W_,&wDlist); 
       wz=pWidth(Z_,&zDlist);       
       sprintf(txt,"%s,%s",E_,ne_); wBrE=findBr(wDlist,txt);
       sprintf(txt,"%s,%s",M_,nm_); wBrM=findBr(wDlist,txt);
       sprintf(txt,"%s,%s",ne_,Ne_);  zBrEn=findBr(zDlist,txt);         
       sprintf(txt,"%s,%s",nm_,Nm_);  zBrMn=findBr(zDlist,txt);
  } 
       
  sprintf(processX3,"%s",name[0]); 
  if(Nin>1) sprintf(processX3+strlen(processX3),",%s",name[1]);
  sprintf(processX3+strlen(processX3),"->%s,",name[lV_]);
  c=processX3+strlen(processX3);
             
  if(abs(pdg[lV_])==11 || abs(pdg[lV_])==12) switch(pdg[lV])
  {   case -24: sprintf(c,"%s,%s",m_,Nm_); *wV=ww; *br=wBrM; break;
      case  24: sprintf(c,"%s,%s",M_,nm_); *wV=ww; *br=wBrM; break;
      case  23: sprintf(c,"%s,%s",nm_,Nm_);  *wV=wz; *br=zBrMn; break;
  } else 
  {  
      switch(pdg[lV])
      { case -24: sprintf(c,"%s,%s",e_,Ne_); *wV=ww; *br=wBrE; break;
        case  24: sprintf(c,"%s,%s",E_,ne_); *wV=ww; *br=wBrE; break;
        case  23: sprintf(c,"%s,%s",ne_,Ne_);  *wV=wz; *br=zBrEn; break; 
       }
  }
  process2Lib(processX3,plib13); 
  if(Nin==2) sprintf(exclude,"%s","%Z+W<1"); else
  { if(pdg[lV]==23) sprintf(exclude,"%s<1",Z_); else sprintf(exclude,"%s<1",W_);}
  strcat(plib13,"V");    
  ccx3=getMEcode(0,ForceUG,processX3,exclude,"",plib13);
  if(ccx3)passParameters(ccx3);
  return ccx3;
}




/*================  kinematic 1->3 and  1->4  =======================*/

static double kinematic_1_3(REAL *pmass, int i3, double m12, double xcos, REAL * P)
{ 
  double factor;
  REAL pout,chY,shY,xsin, E1,P12,P13,E2,P22,P23, m0,m1,m2,m3;
  int i,i1,i2;
  
  for(i=1;i<4;i++)if(i3!=i) {i1=i; break;}
  for(i++;i<4;i++)if(i3!=i) {i2=i; break;}
  
  m0=pmass[0];
  m1=pmass[i1];
  m2=pmass[i2];
  m3=pmass[i3];

  if(m12<=m1+m2) return 0;
  for(i=0;i<16;i++) P[i]=0;

  P[0]=m0; 
  factor=1/(64*M_PI*M_PI*M_PI*m0*m0);
  
  pout=decayPcm(m0,m12,m3);
  if(!pout) return 0;  
  P[i3*4]=Sqrt(pout*pout+m3*m3); P[i3*4+3]=-pout; 

  factor*=pout;  
  
  shY=pout/m12;
  chY=Sqrt(1+shY*shY);  
  pout=decayPcm(m12,m1,m2);
  if(!pout) return 0;
  factor*=pout;
  xsin=Sqrt(1-xcos*xcos);
  E1=Sqrt(m1*m1+pout*pout);    E2=Sqrt(m2*m2+pout*pout);
  P13=xcos*pout;               P23=-P13;
  P12=xsin*pout;               P22=-P12;
  
  P[4*i1]  =chY*E1 + shY*P13;  P[4*i2]  =chY*E2 + shY*P23;
  P[4*i1+3]=shY*E1 + chY*P13;  P[4*i2+3]=shY*E2 + chY*P23;
  P[4*i1+2]=P12;               P[4*i2+2]=P22;
  
  return factor;
}

static double kinematic_1_4(REAL *pmass, double xm1, double xm2, double xcos1, double xcos2,double fi2,  REAL * P)
{ 
  double factor,M1,M2,Pcm,p1cm,p2cm,chY,shY;
  int i,j;

  factor= 1./(pow(2*M_PI,8)*2*pmass[0]*16);
  
  M1= pmass[1]+pmass[2]+ xm1*(pmass[0]-pmass[1]-pmass[2]-pmass[3]-pmass[4]);
  M2= pmass[3]+pmass[4]+ xm2*(pmass[0]-M1-pmass[3]-pmass[4]); 

  factor*=(pmass[0]-pmass[1]-pmass[2]-pmass[3]-pmass[4])*(pmass[0]-M1-pmass[3]-pmass[4]);
  Pcm=decayPcm(pmass[0],M1,M2);            factor*=4*M_PI*Pcm/pmass[0];
  p1cm=decayPcm(M1,pmass[1],pmass[2]);     
  p2cm=decayPcm(M2,pmass[3],pmass[4]);     factor*=2*M_PI*p1cm*p2cm;
  
  
  P[0]=pmass[0]; P[1]=P[2]=P[3]=0;
  
  P[4+0]=Sqrt(pmass[1]*pmass[1]+p1cm*p1cm);       P[8+0]=Sqrt(pmass[2]*pmass[2]+p1cm*p1cm);
  P[4+1]=0;   
  P[4+2]=p1cm*Sqrt(1-xcos1*xcos1);
  P[4+3]=p1cm*xcos1;
  
          
  P[12+0]=Sqrt(pmass[3]*pmass[3]+p2cm*p2cm);       P[16+0]=Sqrt(pmass[4]*pmass[4]+p2cm*p2cm);
  P[12+1]=Sin(fi2)*p2cm*Sqrt(1-xcos2*xcos2);
  P[12+2]=Cos(fi2)*p2cm*Sqrt(1-xcos2*xcos2); 
  P[12+3]=p2cm*xcos2;  
  
  for(i=1;i<=2;i++) for(j=1;j<=3;j++) P[i*8+j]=-P[i*8-4+j];
  shY=Pcm/M1;
  chY=Sqrt(1+shY*shY);
  for(i=1;i<3;i++)
  { double  p0=P[4*i];
    double  p3=P[4*i+3];
    P[4*i]=chY*p0+shY*p3;
    P[4*i+3]=shY*p0+chY*p3;
  }

  shY=-Pcm/M2;            
  chY=Sqrt(1+shY*shY);  
  for(i=3;i<5;i++)
  { double  p0=P[4*i];
    double  p3=P[4*i+3];
    P[4*i]=chY*p0+shY*p3;
    P[4*i+3]=shY*p0+chY*p3;
  }

//printf("Energy conservation\n");
/*
for(i=0;i<4;i++)
{ double sum=P[0+i]-P[4+i]-P[8+i]-P[12+i]-P[16+i];
  if(fabs(sum/P[0]) > 1.E-4)
  { printf("No Energy conservation %E i=%d  \n",sum/P[0],i);
    exit(22);
  }  
}    
*/
  for(i=0;i<5;i++)
  { double m;
    m=Sqrt(Fabs(P[4*i]*P[4*i]-P[4*i+1]*P[4*i+1]-P[4*i+2]*P[4*i+2]-P[4*i+3]*P[4*i+3]));
    if(Fabs(m-pmass[i])>pmass[0]*1.E-5) { printf("wrong mass %d (%E != %E) \n",i,m,(double)pmass[i]); exit(33);}
  }
  return factor;
}


/* ===========  Intergration ================ */

double (*sqme)(int nsub,double GG, REAL *pvect, REAL*cb_coeff, int * err_code)=NULL;
static int  nsub_stat;
static REAL*Q=NULL;
static REAL Pmass[5];
static int i3_;
static double M_;
static double GG=1.2;

static double dWidthdCos(double xcos)
{
  double factor;
  REAL P[16];
  int err_code=0;


  factor=kinematic_1_3(Pmass,i3_,M_,xcos, P);
  
  if(factor==0) return 0;
  
//printf("xcos=%e factor=%E sqme=%e\n", xcos,factor,(*sqme)(nsub_stat,GG,P,&err_code));  
  return  factor*(*sqme)(nsub_stat,GG,P,NULL,&err_code);

}

static int simpsonNeed=1;

static double dWidthdM(double M)
{ 
  M_=M; 
  if(simpsonNeed) return simpson(dWidthdCos,-1.,1.,0.5E-2,NULL);
  else          return peterson21(dWidthdCos,-1,1,NULL);
}

static double width13(numout * cc, int nsub, int * err) 
{
  int i;

  if(passParameters(cc)){ *err=4; return 0;}
  int pdg[4];
  for(i=0;i<4;i++) cc->interface->pinf(nsub,1+i,Pmass+i,pdg+i);
  if(cc->SC ) GG=*(cc->SC); else  GG=sqrt(4*M_PI*alphaQCD(Pmass[0]));
  *err=0;  

  double mm=Pmass[1];  
  i3_=1;  
  for(int i=2;i<=3;i++) if(Pmass[i]>mm) { i3_=i; mm=Pmass[i];}
  
// printf(" mass_1=%E\n", mm);

  double pp[3]={0,0,0};
  
  for(int n=1;;n++)  
  {  int m,w,pnum;
     char*s=cc->interface->den_info(nsub,n,&m,&w,&pnum);
     if(!s) break; 
     double mass=0;
     if(m) mass=fabs(cc->interface->va[m]);
     if(mass>0)
     {  double  mMax=Pmass[0]-Pmass[s[1]-1];
        double  mMin=Pmass[1]+Pmass[2]+Pmass[3]-Pmass[s[1]-1]; 
        if( mMin*0.9<mass && mass<mMax*1.1) pp[s[1]-2]=1;
     }   
  }
   
//  printf("=%E %E %E\n",pp[0],pp[1],pp[2]);            

  simpsonNeed=0;
  for(int i=1;i<=3;i++) if(i!=i3_ && pp[i-1]>0.7 ) simpsonNeed=1; 
     
  double ms=0;for(int i=1;i<=3;i++) if(i!=i3_) ms+=Pmass[i]; 
  sqme=cc->interface->sqme;
  nsub_stat=nsub; 

//  displayPlot("width", "M23",ms, Pmass[0]-mm, 0,1,"",0,dWidthdM,NULL);

  for(int i=1;i<=3;i++)  if(i!=i3_&& (abs(pdg[i])==1 || abs(pdg[i])==2) && Pmass[i] < 0.25 ) ms+=0.25-Pmass[i]; 
  if(ms>=Pmass[0]-mm) return 0;
  
    
  if(pp[i3_-1]<0.7 ) return peterson21(dWidthdM,ms, Pmass[0]-mm,NULL); 
  else               return  simpson(dWidthdM,ms, Pmass[0]-mm,1.E-2,NULL);
  
}


static double wInt14(double *x, double w)
{
   REAL pvect[20];
   int err_code=0;
   double res;
   res=kinematic_1_4(Pmass,x[0], x[1], 2*(x[2]-0.5),2*(x[3]-0.5),2*M_PI*x[4],pvect); 
   if(res==0) return 0;  
   res*= (*sqme)(1,GG, pvect,NULL, &err_code);
   if(err_code) return 0;
   return res*8*M_PI;
}

static double width14(numout * cc, int * err) 
{
  int i;
  vegasGrid * vegPtr;
  long ncall0=10000;                     /* number of integrand calls */
  double alph=1.5;                       /* rate of grid improvement  */
  double ti;                             /* integral estimation */
  double tsi;                            /* standard deviation */
                          
  *err=0;

  if(passParameters(cc)){ *err=4; return 0;}
    
  for(i=0;i<5;i++) cc->interface->pinf(1,1+i,Pmass+i,NULL);  
  if(cc->SC ) GG=*(cc->SC); else  GG=sqrt(4*M_PI*alphaQCD(Pmass[0]));
  *err=0;

  sqme=cc->interface->sqme;
  
  vegPtr=vegas_init(5,wInt14,50);
  
  vegas_int(vegPtr,ncall0,alph,nPROCSS,&ti,&tsi); 
  vegas_int(vegPtr,ncall0,alph,nPROCSS,&ti,&tsi);
  vegas_finish(vegPtr);
   
  return ti;
}


                       
int pname2lib(char*pname, char * libname)
{
  int n,p;
  char buff[30];
  strcpy(buff,pname);
  p=strlen(buff)-1;
  if(buff[p]=='%') {buff[p]=0;p=1;} else p=0;
  n=pTabPos(buff);
  if(!n) {printf("Wrong particle name '''%s'''\n",pname); libname[0]=0; return 1;}
  if(p) { if(n>0) sprintf(libname,"pp%d",n); else sprintf(libname,"ap%d",-n);}
  else  { if(n>0) sprintf(libname,"p%d",n); else sprintf(libname,"a%d",-n);}   
  return 0;
}


static int decodeProcess(char *txt,int*inList,int*outList)
{ char name[20];
   char *ch_,*ch;
   int i,p;   
   ch_=strstr(txt,"->");
   if(!ch_) { inList[0]=0; ch=txt;}
   else
   { 
     for(p=0,ch=txt;; )
     { sscanf(ch," %[^,]",name);
       ch_=strstr(name,"->");
       if(ch_) *ch_=0;
       for(i=strlen(name)-1; i>=0 && name[i]==' '; i--) name[i]=0;
       inList[p]=pTabPos(name);
       if(!inList[p]) return -(p+1);
       p++;
       if(ch_) break;     
       ch=strchr(ch,',');
       if(!ch) break; else ch++;
     }  
     inList[p]=0; 
     ch=strstr(txt,"->")+2; 
   }  

   for(p=0;ch; )
   { sscanf(ch," %[^,]",name);
     for(i=strlen(name)-1;i>=0 && name[i]==' '; i--) name[i]=0;
     outList[p]=pTabPos(name);
     if(!outList[p]) return p+1;
     p++;    
     ch=strchr(ch,',');
     if(!ch) break;
     ch++;
   }
   outList[p]=0;
   return 0;
}


 
void massFilter(double M, txtList * List)
{
  txtList lold=*List, lnew=NULL,lnext;

  while(lold)
  {  double Msum=0,dM;
     char *ch;
     ch=strstr(lold->txt,"->")+2;
     for( ; ch; ch=strchr(ch,','))
     { char buff[10];
       int n;
       char *nm;
        
       ch++;
       sscanf(ch,"%[^,]",buff); 
       
       n=pTabPos(buff);
       nm=ModelPrtcls[abs(n)-1].mass;
       if(nm[0]=='0')   
       {  dM=0;
          switch(abs(ModelPrtcls[abs(n)-1].NPDG))
          {
            case 1: case 2: dM=0.07; break;
            case 3: dM=0.3; break;
            case 4: dM=1.5; break;
            case 5: dM=5. ; break;
          }
       }else if(abs(ModelPrtcls[abs(n)-1].NPDG)==6) dM=MtEff(0);
       else  dM = fabs(*(varAddress(nm))); 
       switch(abs(ModelPrtcls[abs(n)-1].NPDG))
       {
          case 1: case 2: if(dM<0.07)  dM=0.07; break;
          case 3:         if(dM<0.3)   dM=0.3;  break;
          case 4:         if(dM<1.5)   dM=1.5;  break;
          case 5:         if(dM<5)     dM=5. ; break;
      }
      
       Msum+=dM;
     } 
     lnext=lold->next;
     if(M>Msum) {lold->next=lnew; lnew=lold;}
     else {free(lold->txt); free(lold);}
     lold=lnext;
 } 
 *List=lnew;
}

void gammaGluFilter(txtList * List)
{
  txtList lold=*List, lnew=NULL,lnext;

  while(lold)
  {  int del=0,code;
     char *ch;
     ch=strstr(lold->txt,"->")+2;
     for( ; !del && ch; ch=strchr(ch,','))
     { char buff[10];
       ch++;
       sscanf(ch,"%[^,]",buff); 
       code=pNum(buff);
       if(code==22 || code ==21) { del=1;}
     } 
     lnext=lold->next;
     if(del) {free(lold->txt); free(lold);}
     else    {lold->next=lnew; lnew=lold;}
  
     lold=lnext;
 } 
 *List=lnew;
}

int process2Lib(char * process,char * lib)
{ 
  char * ch, *ch1,*c;
  char bufflib[20];
  char *process_;
  int err=0,nX=0;
  process_=malloc(strlen(process)+1);
  strcpy(process_,process);
  int nIn=1;
  char * chAr= strstr(process_,"->");
  chAr[0]=' '; chAr[1]=' ';

  for(;;)
  { c=strchr(process_,','); 
    if(!c) break; else
    { c[0]=' ';
      if(c< chAr) nIn++;
    }
  } 

  char p[7][30];
  int n=sscanf(process_,"%s %s %s %s %s %s %s", p[0],p[1],p[2],p[3],p[4],p[5],p[6]);
  lib[0]=0;
  for(int i=0;i<n;i++) 
  {  if(strstr(p[i],"*x") || strstr(p[i],"*X")) { sscanf(p[i],"%d",&nX); continue;}
     err=pname2lib(p[i],bufflib);
     if(err) return i+1;
     strcat(lib,bufflib);
     if(i==nIn-1) strcat(lib,"_");
  }   
  if(nX) sprintf(lib+strlen(lib),"x%d",nX);
  free(process_);
  return 0;
}            
    


void process2Mass(char * process,double * mass)
{ 
  char * ch;
  char *process_;

  process_=malloc(strlen(process)+1);
  strcpy(process_,process);
    
  ch= strstr(process_,"->");
  ch[0]=' '; ch[1]=' ';
  for(ch=process_; ch[0]; ch++) if(ch[0]==',') ch[0]=' ';
  
  char p[7][30];
  int n=sscanf(process_,"%s %s %s %s %s %s %s", p[0],p[1],p[2],p[3],p[4],p[5],p[6]);
  for(int i=0;i<n;i++) mass[i]=pMass(p[i]);
     
  free(process_);
}


/*======================  1->2 decay ==================*/

double pWidth2(numout * cc, int nsub)
{
  REAL pvect[12];
  double width=0.;
  REAL m1,m2,m3; 
  int pdg1,pdg2,pdg3;
  int ntot,nin,nout;
  double GG;
  procInfo1(cc,&ntot,&nin,&nout);
  if(nsub<1 ||  nsub>ntot|| nin!=1||nout !=2)  return 0;
       
  if(passParameters(cc)) return -1;
  
  cc->interface->pinf(nsub,1,&m1,&pdg1);
  cc->interface->pinf(nsub,2,&m2,&pdg2); 
  cc->interface->pinf(nsub,3,&m3,&pdg3);
  if(cc->SC) GG=*(cc->SC); else GG=sqrt(4*M_PI*alphaQCD(m1));
  
  if(abs(pdg2)<=2 && m2<0.14) m2=0.14;
  if(abs(pdg3)<=2 && m3<0.14) m3=0.14;
  if(abs(pdg2)==3 && m2<0.5 ) m2=0.5;
  if(abs(pdg3)==3 && m3<0.5 ) m3=0.5;
   
  
  if(m1 >m2 + m3)
  {   int i,err_code=0; 
      double md=m2-m3;
      double ms=m2+m3;
      double pRestOut=Sqrt((m1*m1 - ms*ms)*(m1*m1-md*md))/(2*m1);
      double totcoef= pRestOut/(8. * M_PI * m1*m1);
           
      for(i=1;i<12;i++) pvect[i]=0;
      pvect[0]=m1;
      pvect[7]=pRestOut;
      pvect[4]=Sqrt(pRestOut*pRestOut+m2*m2);
      pvect[11]=-pRestOut;
      pvect[8]=Sqrt(pRestOut*pRestOut+m3*m3);
      width = totcoef * (cc->interface->sqme)(nsub,GG,pvect,NULL,&err_code);
  }
  return width;
}

 
double decay2Info(char * pname, FILE* f)
{ int i,j,ntot;
  numout * cc;
  double wtot;
  char pname2[30],process[30],plib[30];
  char * dname[8];

  for(i=0,j=0;pname[i];i++)
  if(pname[i]!=' ') pname2[j++]=pname[i];
  pname2[j]=0;
  strcpy(plib,"2width_");
  pname2lib(pname2,plib+7);
  sprintf(process,"%s->2*x",pname2);
  cc=getMEcode(0,ForceUG,process,NULL,"",plib);
  if(!cc) return -1; 
  procInfo1(cc,&ntot,NULL,NULL); 
  if(f) fprintf(f,"\n Partial width for %s->2x decays in GeV\n",pname2); 
  for(wtot=0,i=1;i<=ntot;i++)
  { double w;
    procInfo2(cc,i,dname,NULL);
    w=pWidth2(cc,i);
    if(w!=0)
    { wtot+=w;
      if(f) fprintf(f,"%3.3s %3.3s  %.2E\n",dname[1],dname[2],w); 
    }
  }
  if(f) fprintf(f," Total width %.2E GeV\n",wtot);
  return  wtot;
}

static int chOpen(numout*cc, int k)
{  double s; 
   int pdg[3],j;
   char*name[3]; 
   for(j=0;j<3;j++) name[j]=cc->interface->pinf(k,j+1,NULL,pdg+j);
   s=pMass(name[0]);
   for(j=1;j<3;j++) s-=pMass(name[j]);
   if( (pdg[0]!=23) &&  (abs(pdg[0])!=24))
   for(j=1;j<3;j++) if(pdg[j]==23 && VZdecay) s-=6; else if(abs(pdg[j])==24 && VWdecay) s-=5;
   if(s>0) return 1; else return 0;
}   
//=========================================================



static double decay22List(const char * pname, txtList *LL)
{ int i,j,ntot,no22;
  numout * cc;
  double wtot,w;
  char pname2[20],process[20],plib[120];
  char * dname[8];
  txtList L=NULL,L_;
  char buff[100];
  int pN=pNum(pname);

  if(pN==0) return 0;
  for(i=0,j=0;pname[i];i++)
  if(pname[i]!=' ') pname2[j++]=pname[i];
  pname2[j]=0;
  strcpy(plib,"2width_");
  pname2lib(pname2,plib+7);
  sprintf(process,"%s->2*x",pname2);
  cc=getMEcode(0,ForceUG,process,NULL,"",plib);
  if(!cc) { if(LL) *LL=NULL; return 0;}  // !!! 
  passParameters(cc);
  procInfo1(cc,&ntot,NULL,NULL);
  for(wtot=0,i=1;i<=ntot;i++)  if(chOpen(cc,i))
  {     
    w=pWidth2(cc,i);
    if(w!=0)
    {  
      procInfo2(cc,i,dname,NULL);    
      if(LL)
       { L_=malloc(sizeof(txtListStr));
         L_->next=L;
         L=L_;
         sprintf(buff,"%E  %s -> %s,%s",w,pname2,dname[1],dname[2]);
         L_->txt=malloc(20+strlen(buff));
         strcpy(L_->txt,buff);
       } 
       wtot+=w;
    }
  }

  no22= L?0:1;
  
  if(pN!=23 && abs(pN)!=24)
  { int k,l;
    REAL m[5];  
    int pdg[5]; 
    char*name[5];
    for(k=1;k<=ntot;k++) if(!chOpen(cc,k))
    { 
      int vd[3]={0,0,0};
      for(i=0;i<3;i++) name[i]=cc->interface->pinf(k,i+1,m+i,pdg+i);
      if(no22 && m[0]<= m[1]+m[2]) continue;
      
      if(pdg[0]==23 ||  abs(pdg[0])==24) continue;

      for(i=1;i<3;i++) vd[i]= (abs(pdg[i])==24 && VWdecay) || (pdg[i]==23 && VZdecay);

      if(vd[1]||vd[2])
      {
         if(vd[1] && vd[2]) { if(m[1]>m[2]) l=1; else l=2;} else {if(vd[1]) l=1; else l=2;} 
         int l_=3-l;
         if(m[0]>m[l_] + 5)    
         {  numout * cc13;
            double wV,brV; 
          
            if( vd[l_] && m[l]<m[l_]) {l=l_; l_=3-l;}
            if( (pdg[l_]==23 || abs(pdg[l_])==24) && m[l]<m[l_]) {l=l_; l_=3-l;}
            cc13=xVtoxll(1,2,name,pdg, l, &wV, &brV);
            if(cc13)         
            {               
              passParameters(cc13);
              *(cc13->interface->BWrange)=20;
              for(i3_=1;i3_<4;i3_++) if(strcmp(cc13->interface->pinf(1,i3_+1,NULL,NULL),name[l_])==0)break;
              for(i=0;i<4;i++) cc13->interface->pinf(1,1+i,Pmass+i,NULL);
              sqme=cc13->interface->sqme;
              
              nsub_stat=1;         
              {  double Mmin,Mmax;
                 Mmax=Pmass[0]-Pmass[i3_];
                 for(Mmin=0,j=1;j<4;j++) if(j!=i3_) Mmin+=Pmass[j]; 
                 w=simpson(dWidthdM, Mmin*1.00001 , Mmax*0.9999,1.E-3,NULL);
//printf("w=%E  %s -> %s %s %s\n", w, cc13->interface->pinf(1,1,NULL,NULL), cc13->interface->pinf(1,2,NULL,NULL),cc13->interface->pinf(1,3,NULL,NULL), cc13->interface->pinf(1,4,NULL,NULL));                 
//             if(findVal(ModelPrtcls[abs(pTabPos(name[l]))-1].width,&wVt)) K=1;else K=wVt/wV; 
                 w/=brV;
              }
              if( vd[l_] )   
              {  double w2= pWidth(name[l_], NULL);
                 w*=decayPcmW(m[0],m[l],m[l_],wV,w2,0)/decayPcmW(m[0],m[l],m[l_],wV,0,0);
                 if(pdg[l_]==pdg[l])  w/=2;
              }
                                                                
              if(w!=0)
              {  if(LL)
                 { 
/*                 
                   L_=malloc(sizeof(txtListStr));
                   L_->next=L;
                   L=L_;
                   sprintf(buff,"%E  %s -> %s,%s",w,name[0],name[1],name[2]);
                   L_->txt=malloc(20+strlen(buff));
                   strcpy(L_->txt,buff);
*/           
                   txtList Lv; 
                   double br;                   
                   char n0[30],n1[30],n2[30];
                   if(pdg[1]==pdg[2]) 
                   {  pWidth(name[1],&Lv);
                      for(;Lv;Lv=Lv->next)
                      { 
                         L_=malloc(sizeof(txtListStr));
                         L_->next=L;
                         L=L_;
                         sscanf(Lv->txt,"%lf %s -> %[^,],%[^,]",&br,n0,n1,n2);
                         sprintf(buff,"%E  %s -> %s,%s,%s",w*br,name[0],name[2],n1,n2);
                         L_->txt=malloc(20+strlen(buff));
                         strcpy(L_->txt,buff);                           
                      }
                   } else if(pdg[1]==-pdg[2])
                   {  pWidth(name[1],&Lv);
                      for(;Lv;Lv=Lv->next)
                      { 
                        L_=malloc(sizeof(txtListStr)); 
                        L_->next=L;
                        L=L_;
                        sscanf(Lv->txt,"%lf %s -> %[^,],%[^,]",&br,n0,n1,n2);
                        sprintf(buff,"%E  %s -> %s,%s,%s",w*br/2,name[0],name[2],n1,n2);  
                        L_->txt=malloc(20+strlen(buff));
                        strcpy(L_->txt,buff);   
                      }
                      pWidth(name[2],&Lv);
                      for(;Lv;Lv=Lv->next)
                      {  L_=malloc(sizeof(txtListStr));
                         L_->next=L;
                         L=L_;
                         sscanf(Lv->txt,"%lf %s -> %[^,],%[^,]",&br,n0,n1,n2);
                         sprintf(buff,"%E  %s -> %s,%s,%s",w*br/2,name[0],name[1],n1,n2);
                         L_->txt=malloc(20+strlen(buff));
                         strcpy(L_->txt,buff);  
                      }  
                   }
                   else                   
                   {
                      pWidth(name[l],&Lv);
                      for(;Lv;Lv=Lv->next)
                      {  L_=malloc(sizeof(txtListStr));
                         L_->next=L;
                         L=L_;
                         sscanf(Lv->txt,"%lf %s -> %[^,],%[^,]",&br,n0,n1,n2);
                         sprintf(buff,"%E  %s -> %s,%s,%s",w*br,name[0],name[l_],n1,n2);
                         L_->txt=malloc(20+strlen(buff));
                         strcpy(L_->txt,buff);   
                      }  
                   }
                 }           
                 wtot+=w; 
              }
            }    
         }
         
      }
    }       
  }
  if(LL) *LL=L;
/*  
  {  for(L_=L;L_;L_=L_->next)
     { 
       sscanf(L_->txt,"%lf %[^\n]",&w,buff);
       sprintf(L_->txt,"%E %s",w/wtot,buff);
     }   
    *LL=L;
  }
*/  
  return  wtot;
}

static txtList conBrList(txtList BrList)
{ 

txtList out=NULL;
  char buff[100];
  double br;
  int inCode[10], outCode[10],i;
  for(;BrList;BrList=BrList->next)
  { txtList new=malloc(sizeof(txtListStr));
    new->next=out;out=new;
    sscanf(BrList->txt,"%lf %[^\n]",&br,buff); 
    decodeProcess(buff,inCode,outCode);
    if(inCode[0]>0) sprintf(buff,"%E  %s -> ",br,ModelPrtcls[inCode[0]-1].aname);   
    else            sprintf(buff,"%E  %s -> ",br,ModelPrtcls[-inCode[0]-1].name);
    for(i=0;outCode[i];i++)
    { if(i) strcat(buff,",");
      if(outCode[i]>0) strcat(buff,ModelPrtcls[outCode[i]-1].aname);   
      else             strcat(buff,ModelPrtcls[-outCode[i]-1].name);  
    }
    new->txt=malloc(strlen(buff)+1);
    strcpy(new->txt,buff);
  }
  return out;  
}


void setQforParticle(REAL *Q,const char*pname)
{
  char *nm;
  REAL*ma;  
  int n,i,cdim;
  int pdg;
  if(!Q) return;
 
  n=pTabPos(pname);
  if(!n){printf("Wrong particle name '%s'\n",pname); return ;}
  nm=ModelPrtcls[abs(n)-1].mass;
  if(nm[0]=='0') return ; else ma=varAddress(nm);

  cdim=abs(ModelPrtcls[abs(n)-1].cdim);
  pdg=abs(ModelPrtcls[abs(n)-1].NPDG);
  
  if(cdim==1){ calcMainFunc(); *Q=fabs(*ma);  calcMainFunc(); return;}
  switch(pdg)
  { case 1:case 2:case 3: *Q=1; return;
    case 4: *Q=1.5; break;
    case 5: *Q=5;   break;
    case 6: *Q=175; break;
  }
  calcMainFunc();
  for(i=0;i<10;i++) 
  { 
    if( fabs(*Q-fabs(*ma)) < 1E-2*(*Q)) break;
    *Q=fabs(*ma);
    calcMainFunc();
  }
} 


static void sortWidthTxtList(txtList L)
{
  if(!L) return;
  for(;;)
  {
    txtList l=L;
    for(l=L; l && l->next; l=l->next)
    { double n1,n2;
      sscanf(l->txt,"%lf",&n1);
      sscanf(l->next->txt,"%lf",&n2);
      if(n1<n2) { char*m=l->txt; l->txt=l->next->txt; l->next->txt=m; break;} 
    }
    if(l->next==NULL) break;   
  } 
}


static double pWidthSTD(const char * name, txtList * LL)
{  
  double width=0;
  txtList Lout=NULL,L;
  
  for(int nout=2; nout<5 && width==0; nout++)
  {  
    if(nout==2)  width=decay22List(name,&Lout); else 
    {  char libName[100];
       L=makeDecayList(name,nout);
       massFilter(pMass(name),&L);
       gammaGluFilter(&L);
       if(L==NULL) continue;
     
       for(txtList l=L;l;l=l->next)                                                    
       { numout* cc;
         double dwidth;                                                                 
         int err=0;                                                                  
         txtList newr; 
         process2Lib(l->txt ,libName);                                               
         cc=getMEcode(0,ForceUG,l->txt,NULL,"",libName);                             
         if(!cc) continue;                                                           
         dwidth=pWidthCC(cc,&err);       
         if(dwidth >0)                                                                
         {                                                                           
           width+=dwidth;                                                               
           newr=malloc(sizeof(txtListStr));                                          
           newr->next=Lout;                                                          
           Lout=newr;                                                                
           newr->txt=malloc(strlen(l->txt)+20);                                      
           sprintf(newr->txt,"%E  %s",dwidth,l->txt);                                 
         }
       }
       
       if(nout==3 && L) {  cleanTxtList(L); break;} 
       cleanTxtList(L); 
    }
  }                                                                    
    
  if(Lout)
  for(L=Lout;L;L=L->next)
  { char buff[100];
    double dwidth;
    sscanf(L->txt,"%lf %[^\n]",&dwidth,buff);
    sprintf(L->txt,"%E %s",dwidth/width,buff);  
  }   
 
  *LL=Lout;
  return width;
}


static double pWidthPlus(const char *name, txtList * LL)
{
  txtList L,l,Lout;
  char libName[100];
  double width;
  int i0,j0,nout;
  REAL Qstat;
  REAL*Q=NULL;

  i0=pTabPos(name);
  if(!i0)  { printf("%s out of model particles\n",name); if(LL) *LL=NULL; return 0; }  
  if(i0>0) { i0--; j0=0;} else { i0=-i0-1; j0=1;} 
  
  if(pMass(name)==0) { if(LL) *LL=NULL; return 0;}
  
  for(int i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0){ Q= varValues+i; break;}
  if(Q) { Qstat=*Q; setQforParticle(Q,name);}  //ok 
    
  width=0;
  Lout=NULL;
  for(nout=2; nout<5 && width==0; nout++)
  {  
    if(nout==2)  width=decay22List(name,&Lout); else 
    { 
       L=makeDecayList(name,nout);
       massFilter(pMass(name),&L);
       gammaGluFilter(&L);
       if(L==NULL) continue;
     
       for(l=L;l;l=l->next)                                                    
       { numout* cc;
         double dwidth;                                                                 
         int err=0;                                                                  
         txtList newr;                                                               
         process2Lib(l->txt ,libName);                                               
         cc=getMEcode(0,ForceUG,l->txt,NULL,"",libName);                             
         if(!cc) continue;                                                           
         dwidth=pWidthCC(cc,&err);       
         if(dwidth >0)                                                                
         {                                                                           
           width+=dwidth;                                                               
           newr=malloc(sizeof(txtListStr));                                          
           newr->next=Lout;                                                          
           Lout=newr;                                                                
           newr->txt=malloc(strlen(l->txt)+20);                                      
           sprintf(newr->txt,"%E  %s",dwidth,l->txt);                                 
         }
       } 
       cleanTxtList(L); 
    }
  }

  if(nout<5) 
  { 
       L= makeDecayList(name,nout);
       massFilter(pMass(name),&L);
       
       gammaGluFilter(&L);
       for(l=L;l;l=l->next)                                                    
       { numout* cc;
         double dwidth;                                                                 
         int err=0;                                                                  
         txtList newr;                                                               
         process2Lib(l->txt ,libName);                                               
         cc=getMEcode(0,ForceUG,l->txt,NULL,"",libName);
           
         if(!cc) continue;  
         passParameters(cc);
         { int n,m,w;
           char*s;
           REAL pmass[6];
           int pdg[6];
//           for(int i=0;i<nout+1;i++) printf("%s ",  cc->interface->pinf(1,i+1,pmass+i,pdg+i)); printf("\n");
           int onSh=0,pnum;
           for(n=1;(s=cc->interface->den_info(1,n,&m,&w,&pnum));n++)  if(m!=0)
           { double m0=pmass[0], mu=fabs(cc->interface->va[m]),mx1=0,mx2=0;   // mx1 -sum  mass of decay product, mx2 - mass of rest 
             if(strchr(s,1)) { for(int i=1; i<1+nout; i++) if(strchr(s,i+1)) mx2+=pmass[i]; else mx1+=pmass[i];}
             else            { for(int i=1; i<1+nout; i++) if(strchr(s,i+1)) mx1+=pmass[i]; else mx2+=pmass[i];} 
             if(m0>mu+mx2 && mu>mx1) {onSh=1;break;}
             if(nout==3)
             { 
//               printf(" mass=%E pnum=%d %d\n", (double)cc->interface->va[m], pnum, ModelPrtcls[pnum].NPDG );
               if(VWdecay && abs(ModelPrtcls[pnum].NPDG)==24) { onSh=1;break;}
               if(VZdecay &&  ModelPrtcls[pnum].NPDG==23 )    { onSh=1;break;}
             }  
           } // else  printf("s[0]=%d s[1]=%d s[2]=%d m=%E\n", s[0],s[1],s[2],(double)cc->interface->va[m] );
          // printf("onSh=%d\n", onSh);
           if(onSh) continue;   
         }                                                        
         dwidth=pWidthCC(cc,&err);       
         if(dwidth >0)                                                                
         {                                                                           
           width+=dwidth;                                                               
           newr=malloc(sizeof(txtListStr));                                          
           newr->next=Lout;                                                          
           Lout=newr;                                                                
           newr->txt=malloc(strlen(l->txt)+20);                                      
           sprintf(newr->txt,"%E  %s",dwidth,l->txt);                                 
         }   
       }
       cleanTxtList(L);
  }
         
  if(Lout)
  {
     for(L=Lout;L;L=L->next)
     { char buff[100];
       double pwidth; 
       sscanf(L->txt,"%lf %[^\n]",&pwidth,buff);
       sprintf(L->txt,"%E %s",pwidth/width,buff);  
     }   
     sortWidthTxtList(Lout);
     if(LL) *LL=Lout;
  } else *LL=NULL;
  
  if(Q) { *Q=Qstat; calcMainFunc();}
  return width;
}


static double  pWidthSLHA(int i0, int j0,txtList *LL, int *err)                                                                                    
{  int pdg,pdg0,Len,decay[10];                                                                                                                     
   int i,j;                                                                                                                                        
   double br;                                                                                                                                      
   double width;                                                                                                                                   
                                                                                                                                                   
   pdg0=ModelPrtcls[i0].NPDG;  
   if(j0) pdg0*=-1;
   for(i=1; allDecays(i,0,&pdg,&Len,decay,&width,&br) ;i++)                                                                       
   {                                                                                                                                               
      if(abs(pdg)==abs(pdg0))                                                                                                                      
      {  txtListStr*l,*L=NULL;                                                                                                                     
         l=malloc(sizeof(txtListStr));                                                                                                             
         decayTable[i0].width=width;                                                                                                               
         decayTable[i0].status=1;                                                                                                                  
         for(j=1; allDecays(i,j,&pdg,&Len,decay,&width,&br) ;j++) if(br>0)                                                                         
         { int k;                                                                                                                                  
           char*ch;                                                                                                                                
           l=malloc(sizeof(txtListStr));                                                                                                           
           l->txt=malloc(100);                                                                                                                     
           l->next=L;                                                                                                                              
           ch=pdg2name(pdg);   if(ch) sprintf(l->txt,"%E  %s -> ",br,ch); else sprintf(l->txt,"%E  #%d -> ",br,pdg);                               
           ch=pdg2name(decay[0]); if(ch)  sprintf(l->txt+strlen(l->txt),"%s",ch); else  sprintf(l->txt+strlen(l->txt),"#%d",decay[0]);             
           for(k=1;k<Len;k++)                                                                                                                      
           { ch=pdg2name(decay[k]);                                                                                                                
             if(ch)sprintf(l->txt+strlen(l->txt),", %s",ch); else sprintf(l->txt+strlen(l->txt),", #%d",decay[k]);                                 
           }                                                                                                                                       
           L=l;                                                                                                                                    
         }                     
         sortWidthTxtList(L);
         if(strcmp(ModelPrtcls[i0].name,ModelPrtcls[i0].aname)==0)  decayTable[i0].pdList[0]=L;
         else 
         {
            if(pdg0==pdg)                                                                                                                             
            {  decayTable[i0].pdList[j0]=L;                                                                                                                                                                                            
               decayTable[i0].pdList[1-j0]=conBrList(L); 
            } else                                                                                                                                    
            {  decayTable[i0].pdList[1-j0]=L;                                                                                                                                                                                            
               decayTable[i0].pdList[j0]=conBrList(L);                                                                                   
            }
         }                                                                                                                                            
         if(LL) *LL=decayTable[i0].pdList[j0];                                                                                                     
         *err=0;                                                                                                                                   
         return width;                                                                                                                             
      }                                                                                                                                            
   }                                                                                                                                               
   *err=1;                                                                                                                                         
   return 0;                                                                                                                                       
}                                                                                                                                                  


double pWidth(const char *name, txtList * LL)
{
  txtList L,l,Lout;
  char libName[100];
  double width;
  int i0,j0,nout;
  REAL Qstat;
  REAL*Q=NULL;
//  Check decay Table 
  i0=pTabPos(name);
  if(i0==0)  { printf("%s out of model particles\n",name); if(LL) *LL=NULL; return 0; }  
  if(i0>0)   { i0--; j0=0;} else { i0=-i0-1; j0=1;}
  
  if(pMass(name)==0) { if(LL) *LL=NULL; return 0;}

  if(decayTable[i0].status==1)
  {     if(LL) *LL=decayTable[i0].pdList[j0];
        return decayTable[i0].width;
  } else if(decayTable[i0].status==-1)
  {     if(LL) *LL=NULL;
        return 0;
  } 

  int pref=decayTable[i0].pref;

  int err=1;
  if(pref==2 || pref==3 || (pref==4 && useSLHAwidth)) width=pWidthSLHA(i0, j0,&L,&err);
  if(err==0)
  { if(LL) *LL=L;
    return width;
  }  

  if(dynamic_cs_mutex) 
  {   pthread_mutex_lock(dynamic_cs_mutex); 
      if(decayTable[i0].status==1)   // calculated while was locked 
      {        
        if(LL) *LL=decayTable[i0].pdList[j0];
        pthread_mutex_unlock(dynamic_cs_mutex);
        return decayTable[i0].width;
      }   
  }
  decayTable[i0].status=-1;
  for(int i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0){ Q= varValues+i; break;}
  if(Q) { Qstat=*Q; setQforParticle(Q,name);}  //ok 
  
  if(pref==0 || pref==2 || pref==4) width=pWidthSTD(name,&L); else width=pWidthPlus(name,&L);
  decayTable[i0].width=width;
  decayTable[i0].status=1;
  decayTable[i0].pdList[j0]=L;

  if(strcmp(ModelPrtcls[i0].name,ModelPrtcls[i0].aname)) decayTable[i0].pdList[1-j0]=conBrList(L);
  
  
  if(dynamic_cs_mutex) 
  {   pthread_mutex_lock(dynamic_cs_mutex); 
      if(decayTable[i0].status==1)   // calculated while was locked 
      {        
        if(LL) *LL=decayTable[i0].pdList[j0];
        pthread_mutex_unlock(dynamic_cs_mutex);
        return decayTable[i0].width;
      }   
  }

  if(Q) { *Q=Qstat; calcMainFunc();} 
  
  decayTable[i0].status=1;
  
  if(LL) *LL=L;
   if(dynamic_cs_mutex){ pthread_mutex_unlock(dynamic_cs_mutex);  /* printf("pWidth(%s) unlock =%e\n",name,sum);*/ }

  return width;
}




int  pWidthPref(const char *name, int pref)
{  if(pref<0 || pref>4) { printf("Unlegal second argument of pWidthPref. It has to be a number betweem 0 and 4\n"); return 2;}
   int i0=pTabPos(name);
   if(i0==0)  { printf("%s out of model particles\n",name); return 1;} 
   i0=abs(i0)-1;
   if(!decayTable) cleanDecayTable();
   if(decayTable[i0].pref==pref) return 0;
   for(int j=0;j<2;j++) if(decayTable[i0].pdList[j])  cleanTxtList(decayTable[i0].pdList[j]); 
   decayTable[i0].width=0;
   decayTable[i0].status=0;
   decayTable[i0].pref=pref;  
   return 0;     
}




double aWidth(char *name) {  return pWidth(name,NULL); }

static int pListEq(char * txt1, char * txt2)  
{  char buff[100];
   char rd1[10][10];
   char rd2[10][10];
   int n1,n2,i1,i2;
   char *ch;
    
   strcpy(buff,txt1); while((ch=strchr(buff,','))) ch[0]=' ';
   
   n1=sscanf(buff,"%s %s %s %s %s %s %s %s %s %s",
   rd1[0],rd1[1],rd1[2],rd1[3],rd1[4],rd1[5],rd1[6],rd1[7],rd1[8],rd1[9]); 
   int star=0;
   for(int i=0;i<n1;i++) 
   { if(strcmp(rd1[i],"*")==0) star++; 
     if(star && i+star<n1) strcpy(rd1[i],rd1[i+star]);
   }  
   n1-=star;
      
   strcpy(buff,txt2); while((ch=strchr(buff,','))) ch[0]=' ';
   
   n2=sscanf(buff,"%s %s %s %s %s %s %s %s %s %s",
   rd2[0],rd2[1],rd2[2],rd2[3],rd2[4],rd2[5],rd2[6],rd2[7],rd2[8],rd2[9]); 
   
   if(n1!=n2 && star==0) return 0;
   for(i1=0;i1<n1;i1++)
   { for(i2=0;i2<n2;i2++) if(strcmp(rd1[i1],rd2[i2])==0){rd2[i2][0]=0; break;}
     if(i2==n2) return 0;
   } 
   return 1;
}      

double findBr(txtList L, char * pattern)
{ char buff[100];
  char *ch;
  double width,widthTot=0;
  
  for(;L;L=L->next)
  { 
     sscanf(L->txt,"%lf %[^\n]",&width,buff);
     ch=strstr(buff,"->");
     ch+=2;
     if( pListEq(pattern,ch)) widthTot+=width;
  }
  return widthTot;   
}

/* =============  ProcInfo ================ */

int  procInfo1(numout*cc, int *nsub, int * nin, int *nout)
{
  if(nin) *nin=cc->interface->nin;
  if(nout)*nout=cc->interface->nout;
  if(nsub)*nsub=cc->interface->nprc;
  return 0;
}

int procInfo2(numout*cc,int nsub,char**name,REAL*mass)
{
  int i;
  int ntot=cc->interface->nin+cc->interface->nout;
    
  if(nsub<1 || nsub> cc->interface->nprc) return 2;

  if(name)for(i=0;i<ntot ;i++) 
  name[i]=(cc->interface->pinf)(nsub,i+1,NULL,NULL);

  if(mass)
  {  
    if(passParameters(cc)){ return 4;}
    for(i=0;i<ntot ;i++) cc->interface->pinf(nsub,i+1,mass+i,NULL);     
  }
  return 0;
}

/* ======================decayTable ===================*/
 static int nPrtcls_old=0; 
 void cleanDecayTable(void)
 { int i,j;
 
   if(decayTable) for(i=0;i<nPrtcls_old;i++) for(j=0;j<2;j++) if(decayTable[i].pdList[j]) 
      cleanTxtList(decayTable[i].pdList[j]);    
   decayTable=realloc(decayTable, nModelParticles*sizeof(decayTableStr));
   for(i=0;i<nModelParticles;i++)
   { for(j=0;j<2;j++) decayTable[i].pdList[j]=NULL;
     decayTable[i].width=0;
     decayTable[i].status=0;
     if( nPrtcls_old!=nModelParticles)  decayTable[i].pref=4;
   }
   
   nPrtcls_old=nModelParticles;
   cleanHiggs_AA_GG();
 }

/*============= Export of parameters ============*/
int passParameters(numout*cc)
{
// printf("passParameters(%p  %s %s %s %s)\n",cc, cc->interface->pinf(1,1,NULL,NULL), cc->interface->pinf(1,2,NULL,NULL), cc->interface->pinf(1,3,NULL,NULL),cc->interface->pinf(1,4,NULL,NULL) );
   int i;
   if(dynamic_cs_mutex) pthread_mutex_lock(dynamic_cs_mutex); 
         
   for(i=1;i<=cc->interface->nvar;i++) if(cc->link[i] &&  ((cc->link[i]-varValues)) < *currentVarPtr)  cc->interface->va[i]=*(cc->link[i]); else 
   {   
     int nin=cc->interface->nin;
     int nout=cc->interface->nout;
     for(int  i=1;i<=nin;i++) printf("%s ", cc->interface->pinf(1,i,NULL,NULL));
     printf("->" );
     for(int  i=nin+1;i<=nin+nout;i++) printf("%s ", cc->interface->pinf(1,i,NULL,NULL));
     printf("\n");
     
     
     printf("Value of variable  '%s' needed for calculation of '%s' is not known yet.\n",
      varNames[*currentVarPtr], cc->interface->varName[i]);  
     cc->interface->va[i]=*(cc->link[i]); FError=1;
     exit(0);
   }
   
   int err=cc->interface->calcFunc();

// if(*currentVarPtr!=53)  printf("  currentVarPtr=%d\n", *currentVarPtr);      

   if(dynamic_cs_mutex) pthread_mutex_unlock(dynamic_cs_mutex);    
   if(err>0) { printf("cannot calculate constraint %s\n",cc->interface->varName[err]); return 1;}
   
   return 0;
}


int slhaDecayPrint(char * name, int dVirt, FILE*f)
{
   txtList all;
   int i,dim; 
   int PDG;
   char N[5][P_NAME_SIZE];
   int id[5];
           
   PDG=pNum(name);
   if(!PDG) return 0;
   fprintf(f,"DECAY %d  %E  # %s\n",PDG,pWidth(name,&all),name);
   for(;all;all=all->next)
   {  
      char pn[20], *chB,*chE;
      double br;
      
      sscanf(all->txt,"%s",pn);
      sscanf(pn,"%lf",&br);
      chB=strstr(all->txt,"->");
      chB+=2;
      for(dim=0,chE=chB ; chE;dim++, chE=strchr(chE+1,',')) 
      {  sscanf(chE+1,"%[^,]",N[dim]); 
         trim(N[dim]); 
         id[dim]=pNum(N[dim]);
      }
      if(dVirt && dim==2 && pMass(name)<= pMass(N[0])+pMass(N[1]))  
      {  
         int v[2],k; 
         for(k=0;k<2;k++) v[k] = (id[k]==23 || abs(id[k])==24);
          
         if(v[0]||v[1])
         { 
            txtList LV;
            if(id[0]!=id[1]) br/=2;   
            for(k=0;k<2;k++) if(v[k])
            { 
               pWidth(N[k],&LV);
               for(;LV;LV=LV->next)
               {  double brV;
                  char* chD=strstr(LV->txt,"->")+2;
                  char name1[20],name2[20];
                  sscanf(chD,"%[^,],%s", name1,name2);
                  trim(name1);
                  sscanf(LV->txt,"%lf",&brV);
                  fprintf(f,"   %e  3  %d  %d  %d # %s,%s->%s\n",br*brV , id[1-k] , pNum(name1),pNum(name2),N[1-k],N[k],chD);
               } 
               if(id[0]==id[1]) break;
            }
            continue;   
         }   
      }     
      
      fprintf(f,"   %s   %d  ",pn,dim);
      for(i=0;i<dim;i++) fprintf(f," %d", id[i] ); 
      chB=strstr(all->txt,"->");
      fprintf(f,"  # %s \n",chB+2);
   } 
   fprintf(f,"\n");
   return PDG;
} 

double pWidthCC(numout*cc,int*err)
{ int ntot,nin,nout;
  double res;

  if(!cc)  {*err=1; return 0;}
  procInfo1(cc,&ntot,&nin,&nout);
  if(ntot<0 || nin>1||nout>4) { *err=2; return 0;}
  *err=0;
  switch(nout)
  { case 2: res=pWidth2(cc,1);     return res;
    case 3: res=width13(cc,1,err); return res;
    case 4: res=width14(cc,err);   return res;          
  }
  return 0; // to avoid warning 
} 
