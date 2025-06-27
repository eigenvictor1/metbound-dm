/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs  
================================*/
      
#define MASSES_INFO      /* Display information about mass spectrum  */

#define OMEGA       /*  Calculate Freeze out relic density and display contribution of  individual channels */
//#define FREEZEIN  /*  Calculate relic density in Freeze-in scenario  */
  
/*===== end of Modules  ======*/

/*===== Options ========*/
/*#define SHOWPLOTS*/
     /* Display  graphical plots on the screen */ 
#define CLEAN
/*===== End of DEFINE  settings ===== */


#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void readScanParams(const char* filename,
    double* YChi_min, double* YChi_max, double* YChi_step,
    double* ThetaY_min, double* ThetaY_max, double* ThetaY_step,
    double* mD_min, double* mD_max, double* mD_step,
    double* mS_min, double* mS_max, double* mS_step)
{
    FILE* f = fopen(filename, "r");
    if (!f) { printf("Could not open parameter file '%s'\n", filename); exit(1); }
    char line[256];
    while (fgets(line, sizeof(line), f)) {
        char key[32];
        double val;
        if (sscanf(line, "%31[^=]=%lf", key, &val) == 2) {
            if (strcmp(key, "YChi_min") == 0) *YChi_min = val;
            else if (strcmp(key, "YChi_max") == 0) *YChi_max = val;
            else if (strcmp(key, "YChi_step") == 0) *YChi_step = val;
            else if (strcmp(key, "ThetaY_min") == 0) *ThetaY_min = val;
            else if (strcmp(key, "ThetaY_max") == 0) *ThetaY_max = val;
            else if (strcmp(key, "ThetaY_step") == 0) *ThetaY_step = val;
            else if (strcmp(key, "mD_min") == 0) *mD_min = val;
            else if (strcmp(key, "mD_max") == 0) *mD_max = val;
            else if (strcmp(key, "mD_step") == 0) *mD_step = val;
            else if (strcmp(key, "mS_min") == 0) *mS_min = val;
            else if (strcmp(key, "mS_max") == 0) *mS_max = val;
            else if (strcmp(key, "mS_step") == 0) *mS_step = val;
        }
    }
    fclose(f);
}


int main(int argc,char** argv)
{  
  // Help message
  if(argc!=3)
  { 
      printf(" Correct usage:  ./main <input file name>  <output file name>\n");
      printf("Example: ./main input.txt output.txt\n");
      exit(1);
  }

  printf("Dark Matter Relic Abundance Calculation for Singlet-Doublet Model\n");

  // CalcHEP options
  // ForceUG: Force Unitary Gauge
  // useSLHAwidth: Use SLHA width
  // VZdecay, VWdecay: Decay modes for Z and W bosons
  // These options are set to default values
  // They can be modified as needed
  int ForceUG=0; // to Force Unitary Gauge assign 1
  VZdecay=0; // Z decay mode
  VWdecay=0; // W decay mode
  //useSLHAwidth=0; // Use SLHA width, set to 1 to use SLHA width


  // Note: No parameter file is read since the parameter are set directly in the code.
  
  // Define error variable
  int err;

  // Parameter scan variables
  double YChi_min, YChi_max, YChi_step;
  double ThetaY_min, ThetaY_max, ThetaY_step;
  double mD_min, mD_max, mD_step;
  double mS_min, mS_max, mS_step;

  // Read scan parameters from file
  readScanParams(argv[1],
    &YChi_min, &YChi_max, &YChi_step,
    &ThetaY_min, &ThetaY_max, &ThetaY_step,
    &mD_min, &mD_max, &mD_step,
    &mS_min, &mS_max, &mS_step);

  // Read constrained parameters
  err = calcMainFunc();
  if(err==-1)
  {
    printf("Error in computing the constrained parameters.\n"); exit(1);
  }

  // Read model parameters
  double vev = findValW("vev");
  
  // Parameter scan

  // Open file for output
  FILE *fp = fopen(argv[2], "w"); // Open file for writing
  if(!fp) {
    printf("Failed to open output file '%s'.\n",argv[2]);
    exit(1);
  }
  // Print header for table
  // Print to file in CSV format
  fprintf(fp, "# Dark Matter Relic Abundance for Singlet-Doublet Model\n");
  fprintf(fp, "#############################################\n");
  fprintf(fp, "# Parameter explanation:\n");
  fprintf(fp, "# mS: mass parameter of the singlet scalar\n");
  fprintf(fp, "# mD: mass parameter of the doublet scalar\n");
  fprintf(fp, "# YChi: Yukawa coupling of the singlet scalar to the doublet scalar\n");
  fprintf(fp, "# ThetaY: mixing angle between the two Yukawa couplings, defined as ArcTan(Y2/Y1)\n");
  fprintf(fp, "# MN1, MN2, MN3: mass eigenvalues of the neutralino-like Majorana fermions\n");
  fprintf(fp, "# MChi: physical mass of the chargino-like Dirac fermion\n");
  fprintf(fp, "# Xf: dimensionless freeze-out temperature parameter\n");
  fprintf(fp, "# Omegah2: dark matter relic density times dimensionless Hubble rate squared\n");
  fprintf(fp, "# Note: The masses are in GeV, ThetaY is in radians\n");
  fprintf(fp, "# Format: mS,mD,YChi,ThetaY,Xf,Omegah2\n");
  fprintf(fp, "#############################################\n");
  fprintf(fp, "mS,mD,YChi,ThetaY,MN1,MN2,MN3,MChi,Xf,Omegah2\n");

  printf("Starting parameter scan...\n");

for (double YChi = YChi_min; YChi <= YChi_max; YChi += YChi_step) {
    for (double ThetaY = ThetaY_min; ThetaY <= ThetaY_max; ThetaY += ThetaY_step) {
      for (double mD = mD_min; mD <= mD_max; mD += mD_step) {
        for (double mS = mS_min; mS <= mS_max; mS += mS_step) {
          // Assign values to model parameters
          assignValW("mS", mS);
          assignValW("mD", mD);
          assignValW("YChi", YChi);
          assignValW("ThetaY", ThetaY);
          // re-compute constrained parameters
          err = calcMainFunc();
          if(err==-1) {printf("Error in computing the constrained parameters.\n"); exit(1);}

          // Variable declaration for mass matrix
          double M11, M12, M13, M22, M23, M33;
          int id;
          // Assign mass matrix elements
          M11 = mS;
          M12 = YChi * vev * cos(ThetaY) / sqrt(2.0);
          M13 = YChi * vev * sin(ThetaY) / sqrt(2.0);
          M22 = 0;
          M23 = -mD;
          M33 = 0.;

          // Diagonalize the mass matrix
          // initialize
          initDiagonal();
          // diagonalize
          id = rDiagonal(3, (REAL)M11, (REAL)M12, (REAL)M13, (REAL)M22, (REAL)M23, (REAL)M33);
          if (id == -1) {
            printf("Error in diagonalization of mass matrix\n");
            printf("Error code: %i\n", id);
            exit(1);
          }

          // Get mass eigenvalues
          double m1, m2, m3;
          m1 = (double)MassArray(id, 1);
          m2 = (double)MassArray(id, 2);
          m3 = (double)MassArray(id, 3);

          // Get mixing matrix
          double R11, R12, R13, R21, R22, R23, R31, R32, R33;
          R11 = (double)MixMatrix(id, 1, 1);
          R12 = (double)MixMatrix(id, 1, 2);
          R13 = (double)MixMatrix(id, 1, 3);
          R21 = (double)MixMatrix(id, 2, 1);
          R22 = (double)MixMatrix(id, 2, 2);
          R23 = (double)MixMatrix(id, 2, 3);
          R31 = (double)MixMatrix(id, 3, 1);
          R32 = (double)MixMatrix(id, 3, 2);
          R33 = (double)MixMatrix(id, 3, 3);

          // Explicit checks
          double element11 = R11*m1*R11 + R21*m2*R21 + R31*m3*R31;
          double element12 = R11*m1*R12 + R21*m2*R22 + R31*m3*R32;
          double element13 = R11*m1*R13 + R21*m2*R23 + R31*m3*R33;
          double element21 = R12*m1*R11 + R22*m2*R21 + R32*m3*R31;
          double element22 = R12*m1*R12 + R22*m2*R22 + R32*m3*R32;
          double element23 = R12*m1*R13 + R22*m2*R23 + R32*m3*R33;
          double element31 = R13*m1*R11 + R23*m2*R21 + R33*m3*R31;
          double element32 = R13*m1*R12 + R23*m2*R22 + R33*m3*R32;
          double element33 = R13*m1*R13 + R23*m2*R23 + R33*m3*R33;
          // Check if diagonalization is correct
          // The elements should match the physical masses
          // Allow a small numerical tolerance
          if (element11-M11> 1e-6 || element12-M12 > 1e-6 || element13-M13 > 1e-6 ||
              element21-M12 > 1e-6 || element22-M22 > 1e-6 || element23-M23 > 1e-6 ||
              element31-M13 > 1e-6 || element32-M23 > 1e-6 || element33-M33 > 1e-6) {
            printf("Error in diagonalization: m1=%.2f, m2=%.2f, m3=%.2f\n", m1, m2, m3);
            printf("Computed elements: %.2f, %.2f, %.2f\n", element11, element22, element33);
            exit(1);
          }

          // Assign physical masses
          // Note: micromegas retuns masses in ascending order of absolute value
          assignValW("MN1", m1);
          assignValW("MN2", m2);
          assignValW("MN3", m3);

          // Assign mixing matrix elements
          // Note: mixing matrix R here and U in FeynRules are inverse (i.e. transpose since orthogonal) of each other
          // See also micrOMEGAs manual and FeynRules model file/notebook
          assignValW("U11", R11);
          assignValW("U12", R21);
          assignValW("U13", R31);
          assignValW("U21", R12);
          assignValW("U22", R22);
          assignValW("U23", R32);
          assignValW("U31", R13);
          assignValW("U32", R23);
          assignValW("U33", R33);

          //printf("m1: %.2f\n", m1);
          //printf("m2: %.2f\n", m2);
          //printf("m3: %.2f\n", m3);

          // Declarations and definitions for dark matter calculation
          int fast=1;
          double Beps=1.E-4;
          double Omegah2, Xf;
          char cdmName[10];
          int spin2, charge3, cdim;
          err = sortOddParticles(cdmName);
          // Catch error
          if(err) {
            printf("Can't calculate %s\n",cdmName); return 1;
          }
          
          // Check if dark matter particle is weakly interacting
          for(int k=1;k<=Ncdm;k++)
          { 
            qNumbers(CDM[k], &spin2, &charge3, &cdim);
            if (charge3 || cdim!= 1){
              printf("Dark matter candidate '%s' is not a weakly interacting particle.\n", CDM[k]);
              return 1;
            }
          }

          // Print parameters
          if(Ncdm==1) 
          {
            // Calculate relic density for 1 dark matter sector
            Omegah2 = darkOmega(&Xf,fast,Beps,&err);
            // Write output in CSV format
            fprintf(fp, "%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2e\n", findValW("mS"), findValW("mD"), findValW("YChi"), findValW("ThetaY"), fabs(m1), fabs(m2), fabs(m3), mD, Xf, Omegah2);
          } else
          {  
            printf("More than one dark matter particle found. Only one is supported in this example.\n");
          } 
        }
      }
    }
  }
  // End of parameter scan

  // Close output file
  fclose(fp);
  printf("Results written to '%s'.\n", argv[2]);
  printf("Finished successfully.\n");

  return 0;
}
