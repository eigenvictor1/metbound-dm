      SUBROUTINE MCMCSTEPGM(PAR,PROB,NPROB,IFAIL) 

      IMPLICIT NONE 
 
      INTEGER IFAIL,I,OMGFLAG,MAFLAG,MOFLAG,GMFLAG
      INTEGER NTOT,IDUM,TOTMIN,TOTMAX,NMAX,NPROB,JM,JL
      INTEGER MCFLAG,NSUSY,NGUT,NMES
      PARAMETER (NSUSY=14,NGUT=21,NMES=21)
 
      DOUBLE PRECISION PAR(*),PROB(*),P,RAN2,XDUM,PP
      DOUBLE PRECISION FTSUSY(NSUSY+2),FTGUT(NGUT+2),FTMES(NMES+2)
      DOUBLE PRECISION SMASS(3),PMASS(2),CMASS,SCOMP(3,3),PCOMP(2,2)
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU,Q2
      DOUBLE PRECISION RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
      DOUBLE PRECISION MSUSYEFF,MMESS,N5
      DOUBLE PRECISION ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP
      DOUBLE PRECISION MSUSYEFFCEN,MSUSYEFFDEV,MMESSCEN,MMESSDEV
      DOUBLE PRECISION TBCEN,TBDEV,LCEN,LDEV,KCEN,KDEV,ALCEN,ALDEV
      DOUBLE PRECISION XIFCEN,XIFDEV,XISCEN,XISDEV,MUPCEN,MUPDEV
      DOUBLE PRECISION MSPCEN,MSPDEV,MSCEN,MSDEV,XIUCEN,XIUDEV
      DOUBLE PRECISION LPPCEN,LPPDEV,LTTCEN,LTTDEV,LUCEN,LUDEV
      DOUBLE PRECISION LDCEN,LDDEV, LTCEN,LTDEV,LBCEN,LBDEV
      DOUBLE PRECISION LLCEN,LLDEV,XCEN,XDEV,X
      DOUBLE PRECISION MSUSYEFFMIN,MMESSMIN,TBMIN,LMIN,KMIN,ALMIN
      DOUBLE PRECISION XIFMIN,XISMIN,MUPMIN,MSPMIN,MSMIN,XIUMIN
      DOUBLE PRECISION LPPMIN,LTTMIN,LUMIN,LDMIN,LTMIN,LBMIN,LLMIN
      DOUBLE PRECISION LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,XIU
      DOUBLE PRECISION MSM,MST,LM,LT,KM,KT,HTM,HTT,LPPM,LPPT,LTTM,LTTT
      DOUBLE PRECISION XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      DOUBLE PRECISION MSREF,D,DMIN,DELMB,DELML,DEL1
      DOUBLE PRECISION Q2MEM,DELMEM,XIFMEM,XISMEM,MSMEM
      DOUBLE PRECISION LPPMEM,LTTMEM,PARMEM(25)
      DOUBLE PRECISION neuttot2(5),neuttot(5),neuttot3(5),neuttotrad(5)

      COMMON/INPPAR/ALINP,XIFINP,XISINP,MSINP,MUPINP,MSPINP
      COMMON/MESEXT/XIFMES,XISMES,MSMES,MUPMES,MSPMES,M3HMES
      COMMON/MESGUT/LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,XIU
      COMMON/MESCAL/MSUSYEFF,MMESS,N5
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     . MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     . CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      COMMON/RADCOR/RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
      COMMON/MCMCPAR/MSUSYEFFCEN,MSUSYEFFDEV,MMESSCEN,MMESSDEV,
     . TBCEN,TBDEV,LCEN,LDEV,KCEN,KDEV,ALCEN,ALDEV,XIFCEN,
     . XIFDEV,XISCEN,XISDEV,MUPCEN,MUPDEV,MSPCEN,MSPDEV,
     . MSCEN,MSDEV,XIUCEN,XIUDEV,LPPCEN,LPPDEV,LTTCEN,LTTDEV,
     . LUCEN,LUDEV,LDCEN,LDDEV, LTCEN,LTDEV,LBCEN,LBDEV,
     . LLCEN,LLDEV,XCEN,XDEV,X,
     . MSUSYEFFMIN,MMESSMIN,TBMIN,LMIN,KMIN,ALMIN,
     . XIFMIN,XISMIN,MUPMIN,MSPMIN,MSMIN,XIUMIN,
     . LPPMIN,LTTMIN,LUMIN,LDMIN,LTMIN,LBMIN,LLMIN
      COMMON/FINETUN/FTSUSY,FTGUT,FTMES
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/STEPS/NTOT,IDUM,TOTMIN,TOTMAX,NMAX
      COMMON/GMSCEN/MSREF,D,DMIN,GMFLAG
      COMMON/GMEM/Q2MEM,DELMEM,XIFMEM,XISMEM,MSMEM,
     . LPPMEM,LTTMEM,PARMEM
      COMMON/GMSAVE/MSM,MST,LM,LT,KM,KT,HTM,HTT,
     . LPPM,LPPT,LTTM,LTTT,JM,JL
      COMMON/RENSCALE/Q2
      COMMON/DELMB/DELMB,DELML,DEL1
      COMMON/MCFLAG/MCFLAG
      COMMON/NEUTRALINO_WIDTH/neuttot2,neuttot,neuttot3,neuttotrad

      X=0d0

      IF(IFAIL.EQ.9.OR.(IFAIL.GE.11 .AND. IFAIL.LE.20))THEN
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
       MSUSYEFFCEN=MSUSYEFF
       MMESSCEN=MMESS
       TBCEN=PAR(3)
       LCEN=PAR(1)
       ALCEN=ALINP
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
       MUPCEN=MUPINP
       MSPCEN=MSPINP
       XIUCEN=XIU
       LPPCEN=LPPMES
       LTTCEN=LTTMES
       LUCEN=LUMES
       LDCEN=LDMES
       LTCEN=LTMES
       LBCEN=LBMES
       LLCEN=LLMES
       Q2MEM=Q2
       DELMEM=DELMB
       XIFMEM=XIFMES
       XISMEM=XISMES
       MSMEM=MSM
       LPPMEM=LPPM
       LTTMEM=LTTM
       DO I=1,25
        PARMEM(I)=PAR(I)
       ENDDO
!       WRITE(0,*)"MSUSYEFFCEN=",MSUSYEFFCEN
!       WRITE(0,*)"MMESSCEN=",MMESSCEN
!       WRITE(0,*)"TBCEN=",TBCEN
!       WRITE(0,*)"LCEN=",LCEN
!       WRITE(0,*)"XIUCEN=",XIUCEN
!      ELSE
!       WRITE(0,*)"NO",IFAIL
      ENDIF 
!      WRITE(0,*)""
!      WRITE(0,*)""

      END
