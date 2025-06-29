      SUBROUTINE MCMCSTEPSP(PAR,PROB,NPROB,IFAIL) 

      IMPLICIT NONE 
 
      CHARACTER CHAN*20

      INTEGER I,IFAIL,M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG
      INTEGER AKFLAG,ALFLAG,OMGFLAG,MAFLAG,MOFLAG
      INTEGER NTOT,IDUM,TOTMIN,TOTMAX,NMAX,NPROB
      INTEGER MCFLAG,NSUSY,NGUT,NMES
      PARAMETER (NSUSY=14,NGUT=21,NMES=21)

      DOUBLE PRECISION PAR(*),PROB(*),P,RAN2,XDUM,PP,R
      DOUBLE PRECISION SIG(5,8),S1,S2,ggF13
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5),CB(5),CL(5)
      DOUBLE PRECISION FTSUSY(NSUSY+2),FTGUT(NGUT+2),FTMES(NMES+2)
      DOUBLE PRECISION SMASS(3),PMASS(2),CMASS,SCOMP(3,3),PCOMP(2,2)
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION M0,M12,A0,MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      DOUBLE PRECISION RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
      DOUBLE PRECISION M1INP,M2INP,M3INP,MHDINP,MHUINP,ALINP,AKINP
      DOUBLE PRECISION XIFINP,XISINP,MUPINP,MSPINP,MSINP,M3HINP
      DOUBLE PRECISION M0CEN,M0DEV,M12CEN,M12DEV,TBCEN,TBDEV,A0CEN
      DOUBLE PRECISION A0DEV,M1CEN,M1DEV,M2CEN,M2DEV,M3CEN,M3DEV
      DOUBLE PRECISION MHDCEN,MHDDEV,MHUCEN,MHUDEV,MSCEN,MSDEV,LCEN
      DOUBLE PRECISION LDEV,KCEN,KDEV,ALCEN,ALDEV,AKCEN,AKDEV,MUCEN
      DOUBLE PRECISION MUDEV,XIFCEN,XIFDEV,XISCEN,XISDEV,MUPCEN
      DOUBLE PRECISION MUPDEV,MSPCEN,MSPDEV,M3HCEN,M3HDEV,XCEN,XDEV,X
      DOUBLE PRECISION M0MIN,M12MIN,TBMIN,A0MIN,M1MIN,M2MIN,M3MIN
      DOUBLE PRECISION MHDMIN,MHUMIN,MSMIN,LMIN,KMIN,ALMIN,AKMIN
      DOUBLE PRECISION MUMIN,XIFMIN,XISMIN,MUPMIN,MSPMIN,M3HMIN
      DOUBLE PRECISION M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT
      DOUBLE PRECISION ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT,MQ3GUT,MU3GUT
      DOUBLE PRECISION MD3GUT,MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT
      DOUBLE PRECISION MEGUT,XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      DOUBLE PRECISION DELMB,DELML,DEL1
      DOUBLE PRECISION Q2,Q2MEM,DELMEM,XIFMEM,XISMEM,MUMEM,MDMEM,MSMEM
      DOUBLE PRECISION PARMEM(25)

      COMMON/SOFTGUT/M0,M12,A0
      COMMON/GUTPAR/M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT,
     . ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT,MQ3GUT,MU3GUT,
     . MD3GUT,MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      COMMON/GUTEXT/XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      COMMON/RENSCALE/Q2
      COMMON/DELMB/DELMB,DELML,DEL1
      COMMON/MEM/Q2MEM,DELMEM,XIFMEM,XISMEM,MUMEM,MDMEM,MSMEM,PARMEM
      COMMON/INPPAR/M1INP,M2INP,M3INP,MHDINP,MHUINP,ALINP,AKINP,
     . XIFINP,XISINP,MUPINP,MSPINP,MSINP,M3HINP
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     . MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     . CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      COMMON/RADCOR/RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
      COMMON/MCMCPAR/M0CEN,M0DEV,M12CEN,M12DEV,TBCEN,TBDEV,A0CEN,
     . A0DEV,M1CEN,M1DEV,M2CEN,M2DEV,M3CEN,M3DEV,MHDCEN,MHDDEV,
     . MHUCEN,MHUDEV,MSCEN,MSDEV,LCEN,LDEV,KCEN,KDEV,ALCEN,ALDEV,
     . AKCEN,AKDEV,MUCEN,MUDEV,XIFCEN,XIFDEV,XISCEN,XISDEV,MUPCEN,
     . MUPDEV,MSPCEN,MSPDEV,M3HCEN,M3HDEV,XCEN,XDEV,X,
     . M0MIN,M12MIN,TBMIN,A0MIN,M1MIN,M2MIN,M3MIN,
     . MHDMIN,MHUMIN,MSMIN,LMIN,KMIN,ALMIN,AKMIN,
     . MUMIN,XIFMIN,XISMIN,MUPMIN,MSPMIN,M3HMIN
      COMMON/FINETUN/FTSUSY,FTGUT,FTMES
      COMMON/SCANFLAGS/M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,
     . AKFLAG,ALFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/MCFLAG/MCFLAG
      COMMON/STEPS/NTOT,IDUM,TOTMIN,TOTMAX,NMAX
      COMMON/SMODELS/R,CHAN
      COMMON/LHCSIG/SIG
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     . BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     . BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     . BRSUSY,WIDTH
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/CB/CB,CL

      X=0d0

      IF(IFAIL.EQ.9.OR.IFAIL.GE.11)THEN
       X=1d40
      ENDIF

      IF(IFAIL.EQ.1.OR.IFAIL.EQ.3.OR.IFAIL.EQ.5.OR.IFAIL.EQ.7)THEN
       X=X+(2d0-SMASS(1))*1d20
      ENDIF
      IF(IFAIL.EQ.2.OR.IFAIL.EQ.3.OR.IFAIL.EQ.6.OR.IFAIL.EQ.7)THEN
       X=X+(2d0-PMASS(1))*1d20
      ENDIF
      IF(IFAIL.EQ.4.OR.IFAIL.EQ.5.OR.IFAIL.EQ.6.OR.IFAIL.EQ.7)THEN
       X=X+(2d0-CMASS)*1d20
      ENDIF

      IF(IFAIL.EQ.8)THEN
       X=MIN(X,RMST1)
       X=MIN(X,RMSB1)
       X=MIN(X,MST1)
       X=MIN(X,MST2)
       X=MIN(X,MSB1)
       X=MIN(X,MSB2)
       X=MIN(X,MUL)
       X=MIN(X,MUR)
       X=MIN(X,MDL)
       X=MIN(X,MDR)
       X=MIN(X,MSL1)
       X=MIN(X,MSNT)
       X=MIN(X,MSMU1)
       X=MIN(X,MSNM)
       X=MIN(X,MLR)
       X=MIN(X,MLL)
       X=MIN(X,MNL)
       X=(1d0-X)*1d20
      ENDIF

      IF(IFAIL.EQ.10)THEN
       DO I=1,NPROB
        X=X+DABS(PROB(I))
       ENDDO
       X=(1d0+X)*1d10
      ENDIF

      IF(X.EQ.0d0)X=1d0

      PP=(X-XCEN)/(.28d0*XDEV*MIN(X,XCEN))
      PP=MAX(PP,-8d0*DLOG(10d0))
      PP=MIN(PP,8d0*DLOG(10d0))
      P=1d0/(1d0+DEXP(PP))
      XDUM=RAN2(IDUM)
!      WRITE(0,*)"X",XCEN,X
!      WRITE(0,*)"P",P,XDUM
      IF(P.GE.XDUM)THEN
!       WRITE(0,*)"OK",IFAIL
       XCEN=X
       M0CEN=M0
       M12CEN=M12
       A0CEN=A0
       TBCEN=PAR(3)
       LCEN=PAR(1)
       IF(ALFLAG.EQ.1)THEN
        ALCEN=ALINP
       ENDIF
       IF(AKFLAG.EQ.1)THEN
        AKCEN=AKINP
       ENDIF
       IF(M1FLAG.EQ.1)THEN
        M1CEN=M1INP
       ENDIF
       IF(M2FLAG.EQ.1)THEN
        M2CEN=M2INP
       ENDIF
       IF(M3FLAG.EQ.1)THEN
        M3CEN=M3INP
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-2)THEN
        XIFCEN=XIFINP
       ELSE
        KCEN=PAR(2)
       ENDIF
       IF(MAFLAG.EQ.-1 .OR. MAFLAG.EQ.-3)THEN
        XISCEN=XISINP
       ELSE
        MSCEN=MSINP
       ENDIF
       IF(MAFLAG.EQ.-5)THEN
        MUCEN=PAR(4)
        KCEN=PAR(2)
        XIFCEN=XIFINP
        XISCEN=XISINP
       ELSE
        IF(MHDFLAG.EQ.1)THEN
         MHDCEN=MHDINP
        ENDIF
        IF(MHUFLAG.EQ.1)THEN
         MHUCEN=MHUINP
        ENDIF
       ENDIF
       MUPINP=MUPCEN
       MSPINP=MSPCEN
       M3HINP=M3HCEN
       Q2MEM=Q2
       DELMEM=DELMB
       XIFMEM=XIFGUT
       XISMEM=XISGUT
       MUMEM=MHUGUT
       MDMEM=MDGUT
       MSMEM=MSGUT
       DO I=1,25
        PARMEM(I)=PAR(I)
       ENDDO
!      ELSE
!       WRITE(0,*)"NO",IFAIL
      ENDIF 
!      WRITE(0,*)""
!      WRITE(0,*)""

      END
