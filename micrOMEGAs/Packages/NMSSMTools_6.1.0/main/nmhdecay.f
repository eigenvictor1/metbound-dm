      PROGRAM MAIN

*  Program to compute the NMSSM Higgs masses, couplings and
*  Branching ratios, with experimental and theoretical constraints
*
*   On input:
*
*      PAR(1) = lambda
*      PAR(2) = kappa
*      PAR(3) = tan(beta)
*      PAR(4) = mu (effective mu term = lambda*s)
*      PAR(5) = Alambda (if MA is not an input)
*      PAR(6) = Akappa
*      PAR(7) = mQ3**2
*      PAR(8) = mU3**2
*      PAR(9) = mD3**2
*      PAR(10) = mL3**2
*      PAR(11) = mE3**2
*      PAR(12) = AU3
*      PAR(13) = AD3
*      PAR(14) = AE3
*      PAR(15) = mQ2**2
*      PAR(16) = mU2**2
*      PAR(17) = mD2**2
*      PAR(18) = mL2**2
*      PAR(19) = mE2**2
*      PAR(20) = M1
*      PAR(21) = M2
*      PAR(22) = M3
*      PAR(23) = MA (diagonal doublet CP-odd mass matrix element)
*      PAR(24) = MP (diagonal singlet CP-odd mass matrix element)
*      PAR(25) = AE2
*
*      All these parameters are assumed to be defined in DRbar
*      at the scale Q2, except for tan(beta) defined at MZ.
*      Q2 is either defined by the user in the input file or
*      computed as Q2 = (2*mQ2+mU2+mD2)/4
*
*   On output:
*
*      SMASS(1-3): CP-even masses (ordered)
*
*      SCOMP(1-3,1-3): Mixing angles: if HB(I) are the bare states,
*        HB(I) = Re(H1), Re(H2), Re(S), and HM(I) are the mass eigenstates,
*        the convention is HB(I) = SUM_(J=1,3) SCOMP(J,I)*HM(J)
*        which is equivalent to HM(I) = SUM_(J=1,3) SCOMP(I,J)*HB(J)
*
*      AMASS(1-2): CP-odd masses (ordered)
*
*      PCOMP(1-2,1-2): Mixing angles: if AB(I) are the bare states,
*        AB(I) = Im(H1), Im(H2), Im(S), and AM(I) are the mass eigenstates,
*        the convention is
*        AM(I) = PCOMP(I,1)*(COSBETA*AB(1)+SINBETA*AB(2))
*              + PCOMP(I,2)*AB(3)
*
*      CMASS: Charged Higgs mass
*
*      CU,CD,CV,CJ,CG(i) Reduced couplings of h1,h2,h3 (i=1,2,3) or
*                        a1,a2 (i=4,5) to up type fermions, down type
*                        fermions, gauge bosons, gluons and photons
*                        Note: CV(4)=CV(5)=0
*      CB(I)             Reduced couplings of h1,h2,h3 (i=1,2,3) or
*                        a1,a2 (i=4,5) to b-quarks including DELMB corrections
*
*      WIDTH(i) Total decay width of h1,h2,h3,a1,a2 (i=1..5)
*               with the following branching ratios:
*      BRJJ(i) h1,h2,h3,a1,a2 -> hadrons
*      BRMM(i)        "       -> mu mu
*      BRLL(i)        "       -> tau tau
*      BRCC(i)        "       -> cc
*      BRBB(i)        "       -> bb
*      BRTT(i)        "       -> tt
*      BRWW(i)        "       -> WW (BRWW(4)=BRWW(5)=0)
*      BRZZ(i)        "       -> ZZ (BRZZ(4)=BRZZ(5)=0)
*      BRGG(i)        "       -> gamma gamma
*      BRZG(i)        "       -> Z gamma
*      BRHIGGS(i)   (i=1..5)  -> other Higgses, including:
*        BRHAA(i,j)   hi      -> a1a1, a1a2, a2a2 (i=1..3, j=1..3)
*        BRHCHC(i)    hi      -> h+h- (i=1..3)
*        BRHAZ(i,j)   hi      -> Zaj  (i=1..3, j=1..2)
*        BRHCW(i)  h1,h2,h3   -> h+W- (i=1..3), a1,a2 -> h+W- (i=4,5)
*        BRHHH(i)     h2      -> h1h1, h3-> h1h1, h1h2, h2h2 (i=1..4)
*        BRAHA(i)     a2      -> a1hi (i=1..3)
*        BRAHZ(i,j)   ai      -> Zhj  (i=1,2, j=1..3)
*      BRSUSY(i)    (i=1..5)  -> susy particles, including:
*        BRNEU(i,j,k)         -> neutralinos j,k (i=1..5, j,k=1..5)
*        BRCHA(i,j)           -> charginos 11, 12, 22 (i=1..5, j=1..3)
*        BRHSQ(i,j)   hi      -> uLuL, uRuR, dLdL, dRdR, t1t1, t2t2,
*                                t1t2, b1b1, b2b2, b1b2 (i=1..3, j=1..10)
*        BRASQ(i,j)   ai      -> t1t2, b1b2 (i=1,2, j=1,2)
*        BRHSL(i,j)   hi      -> lLlL, lRlR, nLnL, l1l1, l2l2, l1l2,
*                                ntnt (i=1..3, j=1..7)
*        BRASL(i)     ai      -> l1l2 (i=1,2)
*
*      HCWIDTH  Total decay width of the charged Higgs
*               with the following branching ratios:
*      HCBRM         h+ -> mu nu_mu
*      HCBRL         "  -> tau nu_tau
*      HCBRSU        "  -> s u
*      HCBRBU        "  -> b u
*      HCBRSC        "  -> s c
*      HCBRBC        "  -> b c
*      HCBRBT        "  -> b t
*      HCBRWHT       "  -> neutral Higgs W+, including:
*        HCBRWH(i)   "  -> H1W+, H2W+, h3W+, a1W+, a2W+ (i=1..5)
*      HCBRSUSY      "  -> susy particles,including
*        HCBRNC(i,j) "  -> neutralino i chargino j (i=1..5, j=1,2)
*        HCBRSQ(i)   "  -> uLdL, t1b1, t1b2, t2b1, t2b2 (i=1..5)
*        HCBRSL(i)   "  -> lLnL, t1nt, t2nt (i=1..3)
*
*      MNEU(i)   Mass of neutralino chi_i (i=1,5, ordered in mass)
*      NEU(i,j)  chi_i components of bino, wino, higgsino u&d, singlino
*                (i,j=1..5)
*
*      MCHA(i)       Chargino masses
*      U(i,j),V(i,j) Chargino mixing matrices
*
*  ERRORS: IFAIL = 0..14
*
*  IFAIL = 0         OK
*          1,3,5,7   m_h1^2 < 0
*          2,3,6,7   m_a1^2 < 0
*          4,5,6,7   m_h+^2 < 0
*          8         m_sfermion^2 < 0
*          10        Violation of phenomenological constraint(s)
*          11,12     Problem in integration of RGEs
*          13,14     Convergence problem
*
*  Phenomenological constraints:
*
*      PROB(I)  = 0, I = 1..89: OK
*
*      PROB(1) =/= 0:  chargino too light
*      PROB(2) =/= 0:  excluded by Z -> neutralinos
*      PROB(3) =/= 0:  charged Higgs too light
*      PROB(4) =/= 0:  excluded by ee -> hZ
*      PROB(5) =/= 0:  excluded by ee -> hZ, h -> bb
*      PROB(6) =/= 0:  excluded by ee -> hZ, h -> tautau
*      PROB(7) =/= 0:  excluded by ee -> hZ, h -> invisible
*      PROB(8) =/= 0:  excluded by ee -> hZ, h -> 2jets
*      PROB(9) =/= 0:  excluded by ee -> hZ, h -> 2photons
*      PROB(10) =/= 0: excluded by ee -> hZ, h -> AA -> 4bs
*      PROB(11) =/= 0: excluded by ee -> hZ, h -> AA -> 4taus
*      PROB(12) =/= 0: excluded by ee -> hZ, h -> AA -> 2bs 2taus
*      PROB(13) =/= 0: excluded by Z -> hA (Z width)
*      PROB(14) =/= 0: excluded by ee -> hA -> 4bs
*      PROB(15) =/= 0: excluded by ee -> hA -> 4taus
*      PROB(16) =/= 0: excluded by ee -> hA -> 2bs 2taus
*      PROB(17) =/= 0: excluded by ee -> hA -> AAA -> 6bs
*      PROB(18) =/= 0: excluded by ee -> hA -> AAA -> 6taus
*      PROB(19) =/= 0: excluded by ee -> Zh -> ZAA -> Z + light pairs
*      PROB(20) =/= 0: excluded by stop -> b l sneutrino
*      PROB(21) =/= 0: excluded by stop -> neutralino c
*      PROB(22) =/= 0: excluded by sbottom -> neutralino b
*      PROB(23) =/= 0: squark/gluino too light
*      PROB(24) =/= 0: selectron/smuon too light
*      PROB(25) =/= 0: stau too light
*      PROB(26) =/= 0: lightest neutralino is not LSP or < 511 keV
*      PROB(27) =/= 0: Landau Pole in l, k, ht, hb below MGUT
*      PROB(28) =/= 0: unphysical global minimum
*      PROB(29) =/= 0: Higgs soft masses >> Msusy
*      PROB(30) =/= 0: excluded by DM relic density (checked only if OMGFLAG=/=0)
*      PROB(31) =/= 0: excluded by DM SI WIMP-nucleon xs (checked if |OMGFLAG|=2 or 4)
*      PROB(32) =/= 0: b->s gamma more than 2 sigma away
*      PROB(33) =/= 0: Delta M_s more than 2 sigma away
*      PROB(34) =/= 0: Delta M_d more than 2 sigma away
*      PROB(35) =/= 0: B_s->mu+mu- more than 2 sigma away
*      PROB(36) =/= 0: B+-> tau+nu_tau more than 2 sigma away
*      PROB(37) =/= 0: (g-2)_muon more than 2 sigma away
*      PROB(38) =/= 0: excluded by Upsilon(1S) -> H/A gamma
*      PROB(39) =/= 0: excluded by eta_b(1S) mass measurement
*      PROB(40) =/= 0: BR(B-->X_s mu+ mu-) more than 2 sigma away
*      PROB(41) =/= 0: excluded by ee -> hZ, h -> AA -> 4taus (ALEPH analysis)
*      PROB(42) =/= 0: excluded by top -> b H+, H+ -> c s (CDF, D0)
*      PROB(43) =/= 0: excluded by top -> b H+, H+ -> tau nu_tau (D0)
*      PROB(44) =/= 0: excluded by top -> b H+, H+ -> W+ A1, A1 -> 2taus (CDF)
*      PROB(45) =/= 0: excluded by top -> bH+ (ATLAS)
*      PROB(46) =/= 0: No Higgs in the MHmin-MHmax GeV range
*      PROB(50) =/= 0: excluded by H->VV (ATLAS+CMS)
*      PROB(51) =/= 0: excluded by H/A->tautau (ATLAS+CMS)
*      PROB(52) =/= 0: excluded by HSM->AA->mus/taus/bs (ATLAS+CMS)
*      PROB(53) =/= 0: excluded by ggF/VBF->H/A->gamgam (ATLAS+CMS)
*      PROB(55) =/= 0: b -> d gamma more than 2 sigma away
*      PROB(56) =/= 0: B_d -> mu+ mu- more than 2 sigma away
*      PROB(57) =/= 0: b -> s nu nubar more than 2 sigma away
*      PROB(58) =/= 0: b -> c tau nu more than 2 sigma away (as SM)
*      PROB(59) =/= 0: K -> pi nu nubar more than 2 sigma away
*      PROB(60) =/= 0: DMK / epsK more than 2 sigma away
*      PROB(61) =/= 0: excluded by DM SD WIMP-neutron xs (checked if |OMGFLAG|=2 or 4)
*      PROB(62) =/= 0: excluded by DM SD WIMP-proton xs (checked if |OMGFLAG|=2 or 4)
*      PROB(63) =/= 0: excluded by HSM->AA->4gam (ATLAS+CMS)
*      PROB(64) =/= 0: excluded by trilepton searches for charg(neutral)inos (CMS)
*      PROB(65) =/= 0: excluded by light mesons or eta_{c,b} decays
*      PROB(66) =/= 0: uncertainty on SM like Higgs mass > 3 GeV
*      PROB(67) =/= 0: k_WZ(HSM) 2 sigma away from LHC measured value
*      PROB(68) =/= 0: k_top(HSM) 2 sigma away from LHC measured value
*      PROB(69) =/= 0: k_bot(HSM) 2 sigma away from LHC measured value
*      PROB(70) =/= 0: k_glu(HSM) 2 sigma away from LHC measured value
*      PROB(71) =/= 0: k_gam(HSM) 2 sigma away from LHC measured value
*      PROB(72) =/= 0: k_tau(HSM) 2 sigma away from LHC measured value
*      PROB(73) =/= 0: B_bsm(HSM) 2 sigma away from LHC measured value
*      PROB(75) =/= 0: excluded by HSM->Z+A (ATLAS)
*      PROB(76) =/= 0: excluded by HSM->Z+A (CMS)
*      PROB(77) =/= 0: excluded by H/A->toptop (CMS)
*      PROB(78) =/= 0: excluded by A->Z+HSM (CMS)
*      PROB(79) =/= 0: excluded by A->Z+HSM (ATLAS)
*      PROB(80) =/= 0: excluded by H/A->Z+A/H (CMS)
*      PROB(81) =/= 0: excluded by H/A->Z+A/H (ATLAS)
*      PROB(82) =/= 0: excluded by H/A->HSM+H/A->2tau2b (CMS)
*      PROB(83) =/= 0: excluded by H/A->HSM+H/A->4b (CMS)
*      PROB(74) =/= 0: excluded by H/A->HSM+H/A->2gam2b (CMS)
*      PROB(84) =/= 0: excluded by VBF->H->HSM+HSM (ATLAS)
*      PROB(85) =/= 0: excluded by ggF->H->HSM+HSM (ATLAS)
*      PROB(86) =/= 0: excluded by ggF->H->HSM+HSM (CMS)
*      PROB(87) =/= 0: excluded by delta_MW
*      PROB(88) =/= 0: excluded by LHC SUSY searches
*      PROB(89) =/= 0: soft SUSY breaking terms beyond 3 TeV
*
************************************************************************

      IMPLICIT NONE

      CHARACTER(200) PRE,SUF,PAT,IFILE,OFILE,EFILE,TFILE,SFILE

      INTEGER NPROB,NPAR
      PARAMETER (NPROB=89,NPAR=25)
      INTEGER IFAIL,DIFAIL,I,IMAX,NMSFLAG,OMGFLAG,MAFLAG
      INTEGER MOFLAG,UNCERTFLAG,PFLAG,GRFLAG,MWFLAG
      INTEGER IND,strlen,CFLAG(6)

      DOUBLE PRECISION PAR(NPAR),PROB(NPROB),SM
      DOUBLE PRECISION M32,CGR,MPL,Q2,DELMB,DELML,DEL1,D0,EPS
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H

      COMMON/SUSYEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/NMSFLAG/NMSFLAG
      COMMON/DELMB/DELMB,DELML,DEL1
      COMMON/M32/M32,CGR,MPL,GRFLAG
      COMMON/RENSCALE/Q2
      COMMON/UNCERTFLAG/UNCERTFLAG
      COMMON/MWFLAG/MWFLAG
      COMMON/CFLAG/CFLAG
      COMMON/PFLAG/PFLAG
      COMMON/FILES/PAT,IFILE,OFILE,EFILE,TFILE,SFILE

      EPS=1d-2
      IMAX=10

* I/O files

      CALL GET_COMMAND_ARGUMENT(1,IFILE)
      IND=index(IFILE,'inp')
      PRE=trim(IFILE(1:IND-1))
      SUF=trim(IFILE(IND+3:strlen(IFILE)))
      OFILE=trim(PRE)//'spectr'//trim(SUF)
      EFILE=trim(PRE)//'omega'//trim(SUF)
      TFILE=trim(PRE)//'slha'//trim(SUF)
      SFILE=trim(PRE)//'smodels'//trim(SUF)
      IND=0
      DO I=1,strlen(IFILE)
       IF(IFILE(I:I).EQ.'/')IND=I
      ENDDO
      IF(IND.EQ.0)THEN
       PAT='./'
      ELSE
       PAT=trim(IFILE(1:IND))
      ENDIF

*   Initialization

      CALL INITIALIZE()

*   Reading of the input parameters

      CALL INPUT(PAR,NPAR)

!      WRITE(0,*)"MAFLAG=",MAFLAG
!      WRITE(0,*)""
!      WRITE(0,*)"TANB =",PAR(3)
!      WRITE(0,*)"M1 =",PAR(20)
!      WRITE(0,*)"M2 =",PAR(21)
!      WRITE(0,*)"M3 =",PAR(22)
!      WRITE(0,*)"LAMBDA =",PAR(1)
!      WRITE(0,*)"KAPPA =",PAR(2)
!      WRITE(0,*)"MUEFF =",PAR(4)
!      WRITE(0,*)"ALAMBDA =",PAR(5)
!      WRITE(0,*)"AKAPPA =",PAR(6)
!      WRITE(0,*)"XIF =",XIF
!      WRITE(0,*)"XIS =",XIS
!      WRITE(0,*)"MUP =",MUP
!      WRITE(0,*)"MSP =",MSP
!      WRITE(0,*)"M3H =",M3H
!      WRITE(0,*)"MA =",PAR(23)
!      WRITE(0,*)"MP =",PAR(24)
!      WRITE(0,*)"AU3 =",PAR(12)
!      WRITE(0,*)"AD3 =",PAR(13)
!      WRITE(0,*)"AE3 =",PAR(14)
!      WRITE(0,*)"AE2 =",PAR(25)
!      WRITE(0,*)"ML3 =",DSQRT(PAR(10))
!      WRITE(0,*)"ML2 =",DSQRT(PAR(18))
!      WRITE(0,*)"ME3 =",DSQRT(PAR(11))
!      WRITE(0,*)"ME2 =",DSQRT(PAR(19))
!      WRITE(0,*)"MQ3 =",DSQRT(PAR(7))
!      WRITE(0,*)"MQ2 =",DSQRT(PAR(15))
!      WRITE(0,*)"MU3 =",DSQRT(PAR(8))
!      WRITE(0,*)"MU2 =",DSQRT(PAR(16))
!      WRITE(0,*)"MD3 =",DSQRT(PAR(9))
!      WRITE(0,*)"MD2 =",DSQRT(PAR(17))
!      WRITE(0,*)""

*   Initialization of PROB and IFAIL

      DO I=1,NPROB
       PROB(I)=0d0
      ENDDO
      IFAIL=0

*   Check for soft terms < 3 TeV

      SM=MAX(DABS(PAR(20)),DABS(PAR(21)),DABS(PAR(22)),DABS(PAR(4)),
     .DABS(PAR(5)),DABS(PAR(6)),DABS(PAR(12)),DABS(PAR(13)),
     .DABS(PAR(14)),DABS(PAR(25)),DSQRT(PAR(10)),DSQRT(PAR(18)),
     .DSQRT(PAR(11)),DSQRT(PAR(19)),DSQRT(PAR(7)),DSQRT(PAR(15)),
     .DSQRT(PAR(8)),DSQRT(PAR(16)),DSQRT(PAR(9)),DSQRT(PAR(17)))
      PROB(89)=DDIM(SM/3d3,1d0)

*   Begin loop to compute DELMB

      UNCERTFLAG=0
      DELMB=.1d0
      DELML=0d0
      DEL1=0d0
!      WRITE(0,*)""
!      WRITE(0,*)"UNCERTFLAG",UNCERTFLAG
!      WRITE(0,*)"DELMB guess",DELMB
      I=0
 1    I=I+1
      IF(I.GT.IMAX)THEN
       IFAIL=14
!       WRITE(0,*)"Exit : IFAIL =",IFAIL
!       WRITE(0,*)""
       GOTO 11
      ENDIF
      D0=DELMB

*   Computation of parameters at QSTSB

      CALL RUNPAR(PAR)

*   Computation of sfermion masses

      CALL MSFERM(PAR,IFAIL,2)
!      WRITE(0,*)"DELMB,IFAIL",DELMB,IFAIL
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"Exit : IFAIL =",IFAIL
!       WRITE(0,*)""
       GOTO 11
      ENDIF

*   End loop to compute DELMB

      IF(2d0*DABS(DELMB-D0)/MAX(1d-3,DABS(DELMB+D0)).GT.EPS)GOTO 1

*   Computation of Higgs masses

      CALL MHIGGS(PAR,PROB,IFAIL)
      IF(IFAIL.EQ.-1)IFAIL=13
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"Exit : IFAIL =",IFAIL
!       WRITE(0,*)""
       GOTO 11
      ENDIF

*   Begin estimate uncertainties

      IF(PFLAG.GT.2)THEN
      DO UNCERTFLAG=1,2

*   Begin loop to compute DELMB

        DIFAIL=0
        DELMB=.1d0
        DELML=0d0
        DEL1=0d0
!        WRITE(0,*)""
!        WRITE(0,*)"UNCERTFLAG",UNCERTFLAG
!        WRITE(0,*)"DELMB guess",DELMB
        I=0
 2      I=I+1
        IF(I.GT.IMAX)THEN
         DIFAIL=14
!         WRITE(0,*)"Exit : DIFAIL =",DIFAIL
         GOTO 3
        ENDIF
        D0=DELMB

*   Computation of parameters at QSTSB

      CALL RUNPAR(PAR)

*   Computation of sfermion masses

        CALL MSFERM(PAR,DIFAIL,2)
!        WRITE(0,*)"DELMB,DIFAIL",DELMB,DIFAIL
        IF(DIFAIL.NE.0)THEN
!         WRITE(0,*)"Exit : DIFAIL =",DIFAIL
         GOTO 3
        ENDIF

*   End loop to compute DELMB

      IF(2d0*DABS(DELMB-D0)/MAX(1d-3,DABS(DELMB+D0)).GT.EPS)GOTO 2

*   Computation of Higgs masses

 3      CALL MHIGGS(PAR,PROB,DIFAIL)

      ENDDO

* Restore original parameters:

      UNCERTFLAG=3

*   Begin loop to compute DELMB

      DELMB=.1d0
      DELML=0d0
      DEL1=0d0
!      WRITE(0,*)""
!      WRITE(0,*)"UNCERTFLAG",UNCERTFLAG
!      WRITE(0,*)"DELMB guess",DELMB
      I=0
 4    I=I+1
      IF(I.GT.IMAX)THEN
       IFAIL=14
!       WRITE(0,*)"Exit : IFAIL =",IFAIL
!       WRITE(0,*)""
       GOTO 11
      ENDIF
      D0=DELMB

*   Computation of parameters at QSTSB

      CALL RUNPAR(PAR)

*   Computation of sfermion masses

      CALL MSFERM(PAR,IFAIL,2)
!      WRITE(0,*)"DELMB,IFAIL",DELMB,IFAIL
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"Exit : IFAIL =",IFAIL
!       WRITE(0,*)""
       GOTO 11
      ENDIF

*   End loop to compute DELMB

      IF(2d0*DABS(DELMB-D0)/MAX(1d-3,DABS(DELMB+D0)).GT.EPS)GOTO 4

*   Computation of Higgs masses

      CALL MHIGGS(PAR,PROB,IFAIL)
      IF(IFAIL.EQ.-1)IFAIL=13
      IF(IFAIL.NE.0)THEN
!       WRITE(0,*)"Exit : IFAIL =",IFAIL
!       WRITE(0,*)""
       GOTO 11
      ENDIF

*   End estimate uncertainties

      ENDIF

*   Computation of gluino mass

      CALL GLUINO(PAR)

*   Computation of chargino masses

      CALL CHARGINO(PAR)

*   Computation of neutralino masses

      CALL NEUTRALINO(PAR)

*   Computation of Higgs + top branching ratios

      CALL DECAY(PAR,PROB)
      CALL TDECAY(PAR)

*   Sparticle decays

      IF(NMSFLAG.NE.0)CALL NMSDECAY(PAR)

*   Exp. constraints on sparticles (LEP, Tevatron)
*   and Higgses (LEP, Tevatron, LHC)

      CALL SUBEXP(PAR,PROB)
      CALL LHCHIG(PAR,PROB)

*   B + K physics

      CALL BOTTOMONIUM(PROB)
      CALL BSG(PAR,PROB)
      CALL KPHYS(PAR,PROB)

*   Anom. magn. moment of the Muon

      CALL MAGNMU(PAR,PROB)

*   Delta_MW

      IF(MWFLAG.NE.0)CALL MWNMSSM(PAR,PROB)

*   Global minimum?

      CALL CHECKMIN(PROB)

*   Landau Pole?

      CALL RGES(PAR,PROB,IFAIL,0)
      IF(IFAIL.NE.0)THEN
       PROB(27)=1d0
       IFAIL=0
      ENDIF

*   RGEs for the soft terms

      CALL RGESOFT(PAR,IFAIL)
      IF(IFAIL.NE.0)THEN
       PROB(27)=1d0
       IFAIL=0
      ENDIF

*   Relic density

      M32=CGR*DSQRT(Q2/3d0)
      CALL RELDEN(PAR,PROB)

*   Exp. constraints on sparticles (LHC)

      IF(CFLAG(5).NE.0)CALL Higgsino_CMS_Trilep(PROB)
      IF(CFLAG(6).NE.0)CALL LHCSUSY(PAR,PROB,1)

*   Check for problems

!      WRITE(0,*)""
      DO I=1,NPROB
       IF(PROB(I).NE.0d0)THEN
!        WRITE(0,*)"PROB",I,PROB(I)
        IFAIL=10
       ENDIF
      ENDDO
!      WRITE(0,*)""

*   Computation of the fine-tuning

      CALL FTPAR(PAR,0)

*   Recording of the results

 11   CALL OUTPUT(PAR,PROB,IFAIL)
!      WRITE(0,*)""

      END


      SUBROUTINE INPUT(PAR,NPAR)

*******************************************************************
*   This subroutine reads SM and NMSSM parameters from input file   .
*******************************************************************

      IMPLICIT NONE

      CHARACTER CHINL*120,CHBLCK*60,CHDUM*120
      CHARACTER(200) PAT,IFILE,OFILE,EFILE,TFILE,SFILE

      INTEGER I,NLINE,INL,ICH,IX,IVAL,Q2FIX
      INTEGER ERR,N0,NLOOP,NBER,NPAR
      INTEGER GMUFLAG,HFLAG,Z3FLAG,OUTFLAG,NMSFLAG,OMGFLAG
      INTEGER MAFLAG,MOFLAG,PFLAG,VFLAG,GRFLAG,MWFLAG,CFLAG(6)

      DOUBLE PRECISION PAR(*),VAL
      DOUBLE PRECISION ACC,XITLA,XLAMBDA,MC0,MB0,MT0
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,Q2,Q2MIN
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
      DOUBLE PRECISION M32,CGR,MPL

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
      COMMON/RENSCALE/Q2
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/SUSYEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/PFLAG/PFLAG
      COMMON/NMSFLAG/NMSFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/VFLAG/VFLAG
      COMMON/OUTFLAG/OUTFLAG
      COMMON/MWFLAG/MWFLAG
      COMMON/CFLAG/CFLAG
      COMMON/M32/M32,CGR,MPL,GRFLAG
      COMMON/FILES/PAT,IFILE,OFILE,EFILE,TFILE,SFILE

* INPUT FILE
      OPEN(15,FILE=IFILE,STATUS='UNKNOWN')

*   INITIALIZATION OF THE SUSY PARAMETERS
      DO I=1,NPAR
       PAR(I)=1d99
      ENDDO
      PAR(2)=0d0
      XIF=1d99
      XIS=1d99
      MUP=0d0
      MSP=0d0
      M3H=0d0
      CGR=1d0
      MPL=2.4d18

*   DEFAULT VALUES FOR FLAGS
      GMUFLAG=1
      PFLAG=2
      OMGFLAG=0
      NMSFLAG=0
      HFLAG=0
      VFLAG=0
      MOFLAG=0
      OUTFLAG=0
      GRFLAG=0
      DO I=1,4
       CFLAG(I)=1
      ENDDO
      CFLAG(5)=0
      CFLAG(6)=0
      MWFLAG=0

*   DEFAULT VALUE FOR THE RENSCALE Q2
      Q2=0d0

*   INITIALIZE READ LOOP
      NLINE=0
      CHBLCK=' '

*   START TO READ NEW LINE INTO CHINL
 21   CHINL=' '

*   LINE NUMBER
      NLINE=NLINE+1

      READ(15,'(A120)',END=29,ERR=999) CHINL

*   CHECK FOR EMPTY OR COMMENT LINES
      IF(CHINL.EQ.' '.OR.CHINL(1:1).EQ.'#'
     .  .OR.CHINL(1:1).EQ.'*') GOTO 21

*   FORCE UPPER CASE LETTERS IN CHINL (AS REQUIRED BELOW)
      INL=0
 22   INL=INL+1
      IF(CHINL(INL:INL).NE.'#')THEN
       DO ICH=97,122
        IF(CHINL(INL:INL).EQ.CHAR(ICH)) CHINL(INL:INL)=CHAR(ICH-32)
       ENDDO
       IF(INL.LT.120) GOTO 22
      ENDIF

*   CHECK FOR BLOCK STATEMENT
      IF(CHINL(1:1).EQ.'B')THEN
       READ(CHINL,'(A6,A)',ERR=999) CHDUM,CHBLCK
       GOTO 21
      ENDIF

*   CHECK FOR NMSSM MODEL IN MODSEL
*   IF THE RELIC DENSITY SHOULD BE COMPUTED
*   THE BLOCK MODSEL MUST CONTAIN THE LINE "  9     1    "
      IF(CHBLCK(1:6).EQ.'MODSEL')THEN
       READ(CHINL,*,ERR=999) IX,IVAL
       IF(IX.EQ.1) Z3FLAG=IVAL
       IF(IX.EQ.8) PFLAG=IVAL
       IF(IX.EQ.9) OMGFLAG=IVAL
       IF(IX.EQ.13) NMSFLAG=IVAL
       IF(IX.EQ.14) VFLAG=IVAL
       IF(IX.EQ.15) MOFLAG=IVAL
       IF(IX.EQ.16) OUTFLAG=IVAL
       IF(IX.EQ.22) CFLAG(5)=IVAL
       IF(IX.EQ.23) MWFLAG=IVAL
       IF(IX.EQ.24) CFLAG(6)=IVAL

*   READ SMINPUTS
      ELSEIF(CHBLCK(1:8).EQ.'SMINPUTS')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.2) GF=VAL
       IF(IX.EQ.3) ALSMZ=VAL
       IF(IX.EQ.4) MZ=VAL
       IF(IX.EQ.5) MB=VAL
       IF(IX.EQ.6) MT=VAL
       IF(IX.EQ.7) MTAU=VAL

*   READ Q2 AND TANBETA
      ELSEIF(CHBLCK(1:6).EQ.'MINPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.0.AND.Q2.EQ.0d0) Q2=VAL**2
       IF(IX.EQ.3) PAR(3)=VAL
       IF(IX.EQ.6) CGR=VAL

*   READ EXTPAR
      ELSEIF(CHBLCK(1:6).EQ.'EXTPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.0) Q2=VAL**2
       IF(IX.EQ.1) PAR(20)=VAL
       IF(IX.EQ.2) PAR(21)=VAL
       IF(IX.EQ.3) PAR(22)=VAL
       IF(IX.EQ.11) PAR(12)=VAL
       IF(IX.EQ.12) PAR(13)=VAL
       IF(IX.EQ.13) PAR(14)=VAL
       IF(IX.EQ.16) PAR(25)=VAL
       IF(IX.EQ.32) PAR(18)=VAL**2
       IF(IX.EQ.33) PAR(10)=VAL**2
       IF(IX.EQ.35) PAR(19)=VAL**2
       IF(IX.EQ.36) PAR(11)=VAL**2
       IF(IX.EQ.42) PAR(15)=VAL**2
       IF(IX.EQ.43) PAR(7)=VAL**2
       IF(IX.EQ.45) PAR(16)=VAL**2
       IF(IX.EQ.46) PAR(8)=VAL**2
       IF(IX.EQ.48) PAR(17)=VAL**2
       IF(IX.EQ.49) PAR(9)=VAL**2
       IF(IX.EQ.61) PAR(1)=VAL
       IF(IX.EQ.62) PAR(2)=VAL
       IF(IX.EQ.63) PAR(5)=VAL
       IF(IX.EQ.64) PAR(6)=VAL
       IF(IX.EQ.65) PAR(4)=VAL
       IF(IX.EQ.66) XIF=VAL
       IF(IX.EQ.67) XIS=VAL
       IF(IX.EQ.68) MUP=VAL
       IF(IX.EQ.69) MSP=VAL
       IF(IX.EQ.72) M3H=VAL
       IF(IX.EQ.124) PAR(23)=VAL
       IF(IX.EQ.125) PAR(24)=VAL

      ENDIF

      GOTO 21

*   END OF READING FROM INPUT FILE

*   Check for errors

 29   CLOSE(15)
      ERR=0
      IF(CFLAG(5).NE.0 .AND. NMSFLAG.EQ.0)THEN
       WRITE(0,2)"LHC CONSTRAINTS CANNOT BE CHECKED IF NMSDECAY IS ",
     .  "NOT CALLED"
       ERR=1
      ENDIF
      IF(CFLAG(6).GT.0 .AND. NMSFLAG.EQ.0)THEN
       WRITE(0,2)"SMODELS CANNOT BE CALLED IF NMSDECAY IS NOT CALLED"
       ERR=1
      ENDIF
      DO I=1,5
       IF(CFLAG(I).LT.0 .OR. CFLAG(I).GT.1)THEN
        WRITE(0,1)"CONSTRAINT FLAGS MUST BE IN [0,1]"
        ERR=1
       ENDIF
      ENDDO
      IF(CFLAG(6).LT.-1 .OR. CFLAG(6).GT.2)THEN
       WRITE(0,1)"SMODELS FLAG MUST BE IN [-1,2]"
       ERR=1
      ENDIF
      IF(VFLAG.LT.0 .OR. VFLAG.GT.1)THEN
       WRITE(0,1)"VFLAG MUST BE IN [0,1]"
       ERR=1
      ENDIF
      IF(MOFLAG.LT.0 .OR. MOFLAG.GT.7)THEN
       WRITE(0,1)"MOFLAG MUST BE IN [0,7]"
       ERR=1
      ENDIF
      IF(OUTFLAG.LT.0 .OR. OUTFLAG.GT.1)THEN
       WRITE(0,1)"OUTFLAG MUST BE IN [0,1]"
       ERR=1
      ENDIF
      IF(MWFLAG.LT.-1 .OR. MWFLAG.GT.1)THEN
       WRITE(0,1)"|MWFLAG| MUST BE IN [0,1]"
       ERR=1
      ENDIF
      IF(PAR(1).EQ.1d99)THEN
       WRITE(0,1)"LAMBDA MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ELSEIF(PAR(1).LE.0d0)THEN
       WRITE(0,1)"LAMBDA MUST BE STRICTLY POSITIVE"
       ERR=1
      ENDIF
      IF(PAR(3).EQ.1d99)THEN
       WRITE(0,1)"TANB MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ELSEIF(PAR(3).LE.0d0)THEN
       WRITE(0,1)"TANB MUST BE STRICTLY POSITIVE"
       ERR=1
      ENDIF
      IF(PAR(4).EQ.1d99)THEN
       WRITE(0,1)"MUEFF MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ELSEIF(PAR(4).EQ.0d0)THEN
       WRITE(0,1)"MUEFF MUST BE NON ZERO"
       ERR=1
      ENDIF
      IF(PAR(21).EQ.1d99)THEN
       WRITE(0,1)"M2 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ELSE
       IF(PAR(20).EQ.1d99)PAR(20)=PAR(21)/2d0
       IF(PAR(22).EQ.1d99)PAR(22)=PAR(21)*3d0
      ENDIF
      IF(PAR(12).EQ.1d99)THEN
       WRITE(0,1)"AU3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(13).EQ.1d99)THEN
       WRITE(0,1)"AD3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(14).EQ.1d99)THEN
       WRITE(0,1)"AE3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(7).EQ.1d99)THEN
       WRITE(0,1)"MQ3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(8).EQ.1d99)THEN
       WRITE(0,1)"MU3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(9).EQ.1d99)THEN
       WRITE(0,1)"MD3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(10).EQ.1d99)THEN
       WRITE(0,1)"ML3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(PAR(11).EQ.1d99)THEN
       WRITE(0,1)"ME3 MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      DO I=15,19
       IF(PAR(I).EQ.1d99)PAR(I)=PAR(I-8)
      ENDDO
      IF(PAR(25).EQ.1d99)PAR(25)=PAR(14)

*   Relations between (ALAMBDA, MA, XIF) and (AKAPPA, MP, XIS)

      IF(PAR(5).NE.1d99 .AND. PAR(23).NE.1d99 .AND. XIF.NE.1d99)THEN
       WRITE(0,1)"AT MOST 2 OF THE 3 PARAMETERS MA, ALAMBDA AND XIF",
     . " CAN BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF

      IF(PAR(6).NE.1d99 .AND. PAR(24).NE.1d99 .AND. XIS.NE.1d99)THEN
       WRITE(0,1)"AT MOST 2 OF THE 3 PARAMETERS MP, AKAPPA AND XIS",
     . " CAN BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF

      IF(PAR(2).EQ.0d0)THEN
       IF(PAR(6).NE.0d0 .AND. PAR(6).NE.1d99)THEN
        WRITE(0,1)"IF KAPPA IS 0, AKAPPA MUST BE 0"
        ERR=1
       ELSE
        PAR(6)=0d0
       ENDIF
       IF(PAR(24).NE.1d99 .AND. XIS.NE.1d99)THEN
        WRITE(0,1)"IF KAPPA IS 0, EITHER MP OR XIS",
     .  " CAN BE GIVEN IN BLOCK EXTPAR"
       ERR=1
       ENDIF
      ENDIF

*   Set default values

      IF(PAR(5).EQ.1d99.AND.PAR(23).EQ.1d99.AND.XIF.EQ.1d99)THEN
       PAR(5)=0d0
       XIF=0d0
      ELSEIF(PAR(5).EQ.1d99.AND.PAR(23).EQ.1d99)THEN
       PAR(5)=0d0
      ELSEIF(PAR(5).EQ.1d99.AND.XIF.EQ.1d99)THEN
       XIF=0d0
      ELSEIF(PAR(23).EQ.1d99.AND.XIF.EQ.1d99)THEN
       XIF=0d0
      ENDIF

      IF(PAR(6).EQ.1d99.AND.PAR(24).EQ.1d99.AND.XIS.EQ.1d99)THEN
       PAR(6)=0d0
       XIS=0d0
      ELSEIF(PAR(6).EQ.1d99.AND.PAR(24).EQ.1d99)THEN
       PAR(6)=0d0
      ELSEIF(PAR(6).EQ.1d99.AND.XIS.EQ.1d99)THEN
       XIS=0d0
      ELSEIF(PAR(24).EQ.1d99.AND.XIS.EQ.1d99)THEN
       XIS=0d0
      ENDIF

*   Set MAFLAG

      IF(PAR(23).EQ.1d99)MAFLAG=0
      IF(PAR(5).EQ.1d99)MAFLAG=1
      IF(XIF.EQ.1d99)MAFLAG=2
      IF(PAR(6).EQ.1d99)MAFLAG=MAFLAG+3
      IF(XIS.EQ.1d99)MAFLAG=MAFLAG+6

*   Check for Z3 breaking terms

      IF(MOD(MAFLAG,3).EQ.2 .OR. MAFLAG/3.EQ.2 .OR.
     . MUP.NE.0d0 .OR. MSP.NE.0d0 .OR. XIF.NE.0d0 .OR.
     . XIS.NE.0d0 .OR. M3H.NE.0d0)THEN
       IF(MOD(PFLAG,3).NE.0)THEN
        WRITE(0,1)
     .  "HIGGS MASS PRECISION = 1, 2, 4, 5, 7 OR 8 FOR Z3-NMSSM"
        ERR=1
       ENDIF
       IF(Z3FLAG.GT.2)THEN
        WRITE(0,1)"PRESENCE OF Z3 BREAKING TERMS"
        ERR=1
       ENDIF
      ENDIF

*   Stop if error

      IF(ERR.EQ.1)THEN
       WRITE(*,1)"ERROR IN INPUT FILE"
       STOP 1
      ENDIF

*   Set Q2MIN, Q2FIX:
      Q2MIN=100d0**2
      Q2FIX=1
      IF(Q2.LE.Q2MIN)THEN
       Q2FIX=0
      ENDIF

*   Initialization for ALPHAS and RUNM (as in hdecay)
*   The bottom quark pole mass MBP is set in INIT and can be changed
*   only there (changing its running mass MB above has no effect
*   on MBP, since one would have to compute alpha_s(MB) first)

      MC0=MC
      MB0=MBP
      MT0=MT
      N0=5
      NLOOP=3
      NBER=18
      ACC=1d-10
      XLAMBDA=XITLA(NLOOP,ALSMZ,ACC)
      CALL ALSINI(ACC)
      CALL BERNINI(NBER)

*    g1,g2  and sin(theta)^2 in the on-shell scheme in terms of
*    GF, MZ(pole) and MW(pole)

      g2=4d0*DSQRT(2d0)*GF*MW**2
      g1=4d0*DSQRT(2d0)*GF*(MZ**2-MW**2)
      S2TW=1d0-(MW/MZ)**2

* LHC constraints on SM Higgs couplings
      CALL COUPLINGFITS(CFLAG(3))

      RETURN

 999  WRITE(0,1)"READ ERROR ON LINE:", NLINE
      WRITE(0,*)CHINL(1:80)
      STOP 1

 1    FORMAT(2A)
 2    FORMAT(A,A)

      END


      SUBROUTINE OUTPUT(PAR,PROB,IFAIL)

*********************************************************************
*   Subroutine writing all the results in the the output files.
*********************************************************************

      IMPLICIT NONE

      CHARACTER(200) PAT,IFILE,OFILE,EFILE,TFILE,SFILE
      CHARACTER CHAN*20

      INTEGER I,NBIN,IFAIL,Q2FIX,NMSFLAG,OMGFLAG,MAFLAG,MOFLAG
      INTEGER PFLAG,VFLAG,OUTFLAG,GRFLAG,MWFLAG,NSUSY,NGUT,NMES
      INTEGER CFLAG(6),PDGH1,PDGH2

      PARAMETER (NSUSY=14,NGUT=21,NMES=21)

      DOUBLE PRECISION PAR(*),PROB(*),SIG(5,8),R,ggF13,CHTBXS
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),AMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC
      DOUBLE PRECISION HCBRBT,HCBRWH(5),HCBRWHT,HCBRNC(5,2)
      DOUBLE PRECISION HCBRSQ(5),HCBRSL(3),HCBRSUSY,HCWIDTH
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5),CB(5),CL(5)
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,TANB,SINB,COSB
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      DOUBLE PRECISION SST,SSB,SSL,Q2,Q2MIN
      DOUBLE PRECISION MGUT,g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION MHUS,MHDS,MSS
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT
      DOUBLE PRECISION HBOTGUT,HTAUGUT,M1GUT,M2GUT,M3GUT,ALGUT,AKGUT
      DOUBLE PRECISION ATGUT,ABGUT,ATAUGUT,AMUGUT,MHUGUT,MHDGUT
      DOUBLE PRECISION MSGUT,MQ3GUT,MU3GUT,MD3GUT,MQGUT,MUGUT
      DOUBLE PRECISION MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      DOUBLE PRECISION XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      DOUBLE PRECISION OMG,OMGMIN,OMGMAX
      DOUBLE PRECISION Xf,sigmaV,vcsll,vcsbb,x(100),dNdx(100),EMIN
      DOUBLE PRECISION sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      DOUBLE PRECISION XSMAX,XENON_SI,XENON_SDn,XENON_SDp,PandaX_SI
      DOUBLE PRECISION LUX_SI,LUX_SDn,LUX_SDp,PICO60_SDp,PandaX_SDn
      DOUBLE PRECISION CRESST_SI,DarkSide50_SI,LZ_SI,LZ_SDn,LZ_SDp
      DOUBLE PRECISION PRINTCHANNELS,omg_
      DOUBLE PRECISION MHUQ,MHDQ,MSX,LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ,QSTSB
      DOUBLE PRECISION BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     . DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     . BRBtaunumax,BRBtaunumin
      DOUBLE PRECISION BRBSll,BRBSllmin,BRBSllMax,
     .             BRBShll,BRBShllmin,BRBShllMax,
     .             BRDG,BRDGmin,BRDGmax,
     .             BRBdMUMU,BRBdMUMUmin,BRBdMUMUmax,
     .             BRBXsnunu,BRBXsnunumin,BRBXsnunumax,
     .             BRBpKpnunu,BRBpKpnunumin,BRBpKpnunuMax,
     .             BRBKsnunu,BRBKsnunumin,BRBKsnunuMax,
     .             RD_taul,RD_taulmin,RD_taulmax,
     .             RDs_taul,RDs_taulmin,RDs_taulmax
      DOUBLE PRECISION BRKp_Pipnunub,BRKp_Pipnunubmin,BRKp_PipnunubMax,
     .       BRKL_Pi0nunub,BRKL_Pi0nunubmin,BRKL_Pi0nunubMax,
     .       DMK,DMKmin,DMKmax,epsK,epsKmin,epsKmax
      DOUBLE PRECISION eps0,epst0,epst1,epst2,epst3,epsts,epstb
      DOUBLE PRECISION epsY32,epsY31,epsY23,epsY13,epscs,epscb
      DOUBLE PRECISION delmagmu,damumin,damumax,amuthmax,amuthmin
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2),toptot
      DOUBLE PRECISION FTSUSY(NSUSY+2),FTGUT(NGUT+2),FTMES(NMES+2)
      DOUBLE PRECISION PX,PA(6),PB(2),PL(7),PK(8)
      DOUBLE PRECISION MHmin,MHmax
      DOUBLE PRECISION M32,CGR,MPL,DELMB,DELML,DEL1
      DOUBLE PRECISION DSMASS(3),DAMASS(2),DCMASS
      DOUBLE PRECISION xsectot,limtrilep
      DOUBLE PRECISION MWNMSSM,dumw,dMW0,DrNMSSM,MWSM,dMWSM,decztt,
     .      deltadecztt,deczee,deltadeczee,BRZTauTau,BRZTTmin,
     .      BRZTTmax,ratio,deltaratio,S2TWeffTau,deltaS2TWeffTau,
     .      S2TWeffSM

      COMMON/EWPO/MWNMSSM,dumw,dMW0,DrNMSSM,MWSM,dMWSM,decztt,
     .      deltadecztt,deczee,deltadeczee,BRZTauTau,BRZTTmin,
     .      BRZTTmax,ratio,deltaratio,S2TWeffTau,deltaS2TWeffTau,
     .      S2TWeffSM
      COMMON/PFLAG/PFLAG
      COMMON/NMSFLAG/NMSFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     . BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     . BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     . BRSUSY,WIDTH
      COMMON/BRC/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,
     . HCBRBT,HCBRWH,HCBRWHT,HCBRNC,HCBRSQ,HCBRSL,
     . HCBRSUSY,HCWIDTH
      COMMON/BRSG/BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     .      DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     .      BRBtaunumax,BRBtaunumin
      COMMON/FLAV2/BRBSll,BRBSllmin,BRBSllMax,
     .      BRBShll,BRBShllmin,BRBShllMax,
     .      BRDG,BRDGmin,BRDGmax,
     .      BRBdMUMU,BRBdMUMUmin,BRBdMUMUmax,
     .      BRBXsnunu,BRBXsnunumin,BRBXsnunumax,
     .      BRBpKpnunu,BRBpKpnunumin,BRBpKpnunuMax,
     .      BRBKsnunu,BRBKsnunumin,BRBKsnunuMax,
     .      RD_taul,RD_taulmin,RD_taulmax,
     .      RDs_taul,RDs_taulmin,RDs_taulmax
      COMMON/FLAV3/BRKp_Pipnunub,BRKp_Pipnunubmin,BRKp_PipnunubMax,
     .             BRKL_Pi0nunub,BRKL_Pi0nunubmin,BRKL_Pi0nunubMax,
     .             DMK,DMKmin,DMKmax,epsK,epsKmin,epsKmax
      COMMON/EPSCOUP/eps0,epst0,epst1,epst2,epst3,epsts,epstb,
     .               epsY32,epsY31,epsY23,epsY13,epscs,epscb
      COMMON/MAGMU/delmagmu,damumin,damumax,amuthmax,amuthmin
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/topwidth/toptot
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/CB/CB,CL
      COMMON/HIGGSPEC/SMASS,SCOMP,AMASS,PCOMP,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     . MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     . CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/MGUT/MGUT
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SUSYMH/MHUS,MHDS,MSS
      COMMON/QMHIGGS/MHUQ,MHDQ,MSX
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT,
     . HBOTGUT,HTAUGUT
      COMMON/GUTPAR/M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT,
     . ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT,MQ3GUT,MU3GUT,MD3GUT,
     . MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      COMMON/GUTEXT/XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      COMMON/MICROMG/OMG,OMGMIN,OMGMAX,Xf,sigmaV,vcsll,vcsbb,
     .      x,dNdx,EMIN,NBIN
      COMMON/MICROMG2/sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      COMMON/FINETUN/FTSUSY,FTGUT,FTMES
      COMMON/EFFCOUP/PX,PA,PB,PL,PK
      COMMON/DELMB/DELMB,DELML,DEL1
      COMMON/LHCSIG/SIG
      COMMON/HIGGSFIT/MHmin,MHmax
      COMMON/VFLAG/VFLAG
      COMMON/OUTFLAG/OUTFLAG
      COMMON/M32/M32,CGR,MPL,GRFLAG
      COMMON/DHIGGSSPEC/DSMASS,DAMASS,DCMASS
      COMMON/XSECTRILEP/xsectot,limtrilep
      COMMON/MWFLAG/MWFLAG
      COMMON/CFLAG/CFLAG
      COMMON/FILES/PAT,IFILE,OFILE,EFILE,TFILE,SFILE
      COMMON/SMODELS/R,CHAN
      COMMON/PDGHIGGS/PDGH1,PDGH2
      COMMON/CHTBXS/CHTBXS

      OPEN(16,FILE=OFILE,STATUS='UNKNOWN')
      IF(OMGFLAG.NE.0)OPEN(17,FILE=EFILE,STATUS='UNKNOWN')

      TANB=PAR(3)
      COSB=1d0/DSQRT(1d0+TANB**2)
      SINB=TANB*COSB

      WRITE(16,899) "# NMSSMTools OUTPUT IN SLHA FORMAT"
      WRITE(16,899) "# Info about spectrum calculator"
      WRITE(16,899) "BLOCK SPINFO   # Program information"
      WRITE(16,930) 1,"NMSSMTools # Spectrum calculator"
      WRITE(16,930) 2,"6.1.0      # Version number"

      IF(PROB(1).NE.0d0)
     . WRITE(16,930) 3,"# Chargino excluded by LEP"
      IF(PROB(2).NE.0d0)
     . WRITE(16,930) 3,"# Neutralinos excluded by LEP"
      IF(PROB(3).NE.0d0)
     . WRITE(16,930) 3,"# Charged Higgs excluded by LEP"
      IF(PROB(4).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ee -> hZ, ind. of h decay"
      IF(PROB(5).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ee -> hZ, h -> bb"
      IF(PROB(6).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ee -> hZ, h -> tautau"
      IF(PROB(7).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ee -> hZ, h -> invisible"
      IF(PROB(8).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ee -> hZ, h -> 2jets"
      IF(PROB(9).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ee -> hZ, h -> 2photons"
      IF(PROB(10).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ee -> hZ, h -> AA -> 4bs"
      IF(PROB(11).NE.0d0 .OR. PROB(41).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ee -> hZ, h -> AA -> 4taus"
      IF(PROB(12).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ee -> hZ, h -> AA -> 2bs 2taus"
      IF(PROB(13).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ee -> Z -> hA (Z width)"
      IF(PROB(14).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ee -> hA -> 4bs"
      IF(PROB(15).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ee -> hA -> 4taus"
      IF(PROB(16).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ee -> hA -> 2bs 2taus"
      IF(PROB(17).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ee -> hA -> AAA -> 6bs"
      IF(PROB(18).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ee -> hA -> AAA -> 6taus"
      IF(PROB(19).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ee -> hZ, h -> AA,A -> light pair"
      IF(PROB(20).NE.0d0)
     . WRITE(16,930) 3,"# excluded by stop -> b l sneutrino"
      IF(PROB(21).NE.0d0)
     . WRITE(16,930) 3,"# excluded by stop -> neutralino c"
      IF(PROB(22).NE.0d0)
     . WRITE(16,930) 3,"# excluded by sbottom -> neutralino b"
      IF(PROB(23).NE.0d0)
     . WRITE(16,930) 3,"# Squark/gluino too light"
      IF(PROB(24).NE.0d0)
     . WRITE(16,930) 3,"# Selectron/smuon too light"
      IF(PROB(25).NE.0d0)
     . WRITE(16,930) 3,"# Stau too light"
      IF(PROB(26).GT.0d0)
     . WRITE(16,930) 3,"# Lightest neutralino is not the LSP"
      IF(PROB(26).LT.0d0)
     . WRITE(16,930) 3,"# Mass of the lightest neutralino < 511 keV"
      IF(PROB(27).NE.0d0)
     . WRITE(16,930) 3,"# Landau Pole below MGUT"
      IF(PROB(28).NE.0d0)
     . WRITE(16,930) 3,"# Unphysical global minimum"
      IF(PROB(29).NE.0d0)
     . WRITE(16,930) 3,"# Higgs soft masses >> Msusy"
      IF(PROB(30).GT.0d0)
     . WRITE(16,930) 3,"# DM relic density too large"
      IF(PROB(30).LT.0d0.AND.PROB(30).GT.-1d0)
     . WRITE(16,930) 3,"# DM relic density too small"
      IF(PROB(30).LE.-1d0)
     . WRITE(16,930) 3,"# Problem in micrOMEGAs"
      IF(PROB(31).NE.0d0)THEN
       WRITE(16,930) 3,"# DM direct detection rate too large (SI)"
       IF(OMGFLAG.LT.0 .AND. OMG.LT.OMGMAX)
     .  WRITE(16,915) "         # Cross sections rescaled with the",
     .  " ratio relic density/Planck observed value"
      ENDIF
      IF(PROB(61).NE.0d0)THEN
       WRITE(16,930) 3,"# DM direct detection rate too large (SD-n)"
       IF(OMGFLAG.LT.0 .AND. OMG.LT.OMGMAX)
     .  WRITE(16,915) "         # Cross sections rescaled with the",
     .  " ratio relic density/Planck observed value"
      ENDIF
      IF(PROB(62).NE.0d0)THEN
       WRITE(16,930) 3,"# DM direct detection rate too large (SD-p)"
       IF(OMGFLAG.LT.0 .AND. OMG.LT.OMGMAX)
     .  WRITE(16,915) "         # Cross sections rescaled with the",
     .  " ratio relic density/Planck observed value"
      ENDIF
      IF(PROB(32).NE.0d0)
     . WRITE(16,930) 3,"# b -> s gamma more than 2 sigma away"
      IF(PROB(33).NE.0d0)
     . WRITE(16,930) 3,"# Delta M_s more than 2 sigma away"
      IF(PROB(34).NE.0d0)
     . WRITE(16,930) 3,"# Delta M_d more than 2 sigma away"
      IF(PROB(35).NE.0d0)
     . WRITE(16,930) 3,"# B_s -> mu+ mu- more than 2 sigma away"
      IF(PROB(36).NE.0d0)
     . WRITE(16,930) 3,"# B+ -> tau nu_tau more than 2 sigma away"
      IF(PROB(37).NE.0d0)
     . WRITE(16,930) 3,"# Muon magn. mom. more than 2 sigma away"
      IF(PROB(38).NE.0d0)
     . WRITE(16,930) 3,"# excluded by Upsilon(1S) -> H/A gamma"
      IF(PROB(38).LT.0d0)
     . WRITE(16,930) 3,"# (but A width> 10 MeV)"
      IF(PROB(39).NE.0d0)
     . WRITE(16,930) 3,
     . "# excluded etab(1S) mass difference (BABAR - theory)"
      IF(PROB(40).NE.0d0)
     . WRITE(16,930) 3,"# excluded by BR(B -> X_s mu +mu-)"
      IF(PROB(42).NE.0d0)
     . WRITE(16,930) 3,"# excluded by top -> b H+, H+ -> c s"
      IF(PROB(43).NE.0d0)
     . WRITE(16,930) 3,"# excluded by top -> b H+, H+ -> tau nu_tau"
      IF(PROB(44).NE.0d0)
     . WRITE(16,930) 3,
     . "# excluded by top -> b H+, H+ -> W+ A1, A1 -> 2taus"
      IF(PROB(45).NE.0d0)
     . WRITE(16,930) 3,"# excluded by top -> bH+ (ATLAS)"
      IF(PROB(46).NE.0d0)
     . WRITE(16,929) 3,"# No Higgs in the",MHMIN,MHMAX," GeV mass range"
      IF(PROB(50).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ggF/VBF->H->VV (ATLAS+CMS)"
      IF(PROB(51).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ggF/bb->H/A->tautau (ATLAS+CMS)"
      IF(PROB(52).NE.0d0)
     . WRITE(16,930) 3,
     ."# excluded by HSM->AA->mus/taus/bs (ATLAS+CMS)"
      IF(PROB(53).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ggF/VBF->H/A->gamgam (ATLAS+CMS)"
      IF(PROB(55).NE.0d0)
     . WRITE(16,930) 3,"# b -> d gamma more than 2 sigma away"
      IF(PROB(56).NE.0d0)
     . WRITE(16,930) 3,"# B_d -> mu+ mu- more than 2 sigma away"
      IF(PROB(57).NE.0d0)
     . WRITE(16,930) 3,"# b -> s nu nubar more than 2 sigma away"
      IF(PROB(58).NE.0d0)
     . WRITE(16,930) 3,"# b -> c tau nu more than 2 sigma away (as SM)"
      IF(PROB(59).NE.0d0)
     . WRITE(16,930) 3,"# K -> pi nu nubar more than 2 sigma away"
      IF(PROB(60).NE.0d0)
     . WRITE(16,930) 3,"# DMK / epsK more than 2 sigma away"
      IF(PROB(63).NE.0d0)
     . WRITE(16,930) 3,"# excluded by HSM->AA->4gam (ATLAS+CMS)"
      IF(PROB(64).NE.0d0)
     . WRITE(16,930) 3,
     . "# excluded by trilepton searches for charg(neutral)inos (CMS)"
      IF(PROB(65).NE.0d0)
     . WRITE(16,930) 3,"# excluded by light mesons or eta_{c,b} decays"
      IF(PROB(66).NE.0d0)
     . WRITE(16,930) 3,"# uncertainty on SM-like Higgs mass > 3 GeV"
      IF(PROB(67).NE.0d0)
     . WRITE(16,930) 3,"# k_W/Z(HSM) more than 2 sigma away"
      IF(PROB(68).NE.0d0)
     . WRITE(16,930) 3,"# k_top(HSM) more than 2 sigma away"
      IF(PROB(69).NE.0d0)
     . WRITE(16,930) 3,"# k_bot(HSM) more than 2 sigma away"
      IF(PROB(70).NE.0d0)
     . WRITE(16,930) 3,"# k_glu(HSM) more than 2 sigma away"
      IF(PROB(71).NE.0d0)
     . WRITE(16,930) 3,"# k_gam(HSM) more than 2 sigma away"
      IF(PROB(72).NE.0d0)
     . WRITE(16,930) 3,"# k_tau(HSM) more than 2 sigma away"
      IF(PROB(73).NE.0d0)
     . WRITE(16,930) 3,"# B_bsm(HSM) more than 2 sigma away"
      IF(PROB(75).NE.0d0)
     . WRITE(16,930) 3,"# excluded by HSM->Z+A (ATLAS)"
      IF(PROB(76).NE.0d0)
     . WRITE(16,930) 3,"# excluded by HSM->Z+A (CMS)"
      IF(PROB(77).NE.0d0)
     . WRITE(16,930) 3,"# excluded by H/A->toptop (CMS)"
      IF(PROB(78).NE.0d0)
     . WRITE(16,930) 3,"# excluded by A->Z+HSM->bb (CMS)"
      IF(PROB(79).NE.0d0)
     . WRITE(16,930) 3,"# excluded by A->Z+HSM->bb (ATLAS)"
      IF(PROB(80).NE.0d0)
     . WRITE(16,930) 3,"# excluded by H/A->Z+A/H (CMS)"
      IF(PROB(81).NE.0d0)
     . WRITE(16,930) 3,"# excluded by H/A->Z+A/H (ATLAS)"
      IF(PROB(82).NE.0d0)
     . WRITE(16,930) 3,"# excluded by H/A->HSM+H/A->2tau2b (CMS)"
      IF(PROB(83).NE.0d0)
     . WRITE(16,930) 3,"# excluded by H/A->HSM+H/A->4b (CMS)"
      IF(PROB(74).NE.0d0)
     . WRITE(16,930) 3,"# excluded by H/A->HSM+H/A->2gam2b (CMS)"
      IF(PROB(84).NE.0d0)
     . WRITE(16,930) 3,"# excluded by VBF->H->HSM+HSM (ATLAS)"
      IF(PROB(85).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ggF->H->HSM+HSM (ATLAS)"
      IF(PROB(86).NE.0d0)
     . WRITE(16,930) 3,"# excluded by ggF->H->HSM+HSM (CMS)"
      IF(PROB(87).LT.0d0)
     . WRITE(16,930) 3,"# MW too small"
      IF(PROB(87).GT.0d0)
     . WRITE(16,930) 3,"# MW too large"
      IF(CFLAG(6).LT.0)THEN
       IF(PROB(88).NE.0d0)
     .  WRITE(16,930) 3,"# sparticles masses < 1 TeV"
      ELSE
       IF(PROB(88).GT.0d0)
     .  WRITE(16,930) 3,"# excluded by LHC SUSY searches"
       IF(PROB(88).LT.0d0)
     .  WRITE(16,930) 3,"# Problem in SmodelS"
      ENDIF
      IF(PROB(89).NE.0d0)
     . WRITE(16,930) 3,"# Soft terms > 3 TeV"

      IF(IFAIL.EQ.1.OR.IFAIL.EQ.3.OR.IFAIL.EQ.5.OR.IFAIL.EQ.7)
     . WRITE(16,930) 4,"# M_H1^2<1"
      IF(IFAIL.EQ.2.OR.IFAIL.EQ.3.OR.IFAIL.EQ.6.OR.IFAIL.EQ.7)
     . WRITE(16,930) 4,"# M_A1^2<1"
      IF(IFAIL.EQ.4.OR.IFAIL.EQ.5.OR.IFAIL.EQ.6.OR.IFAIL.EQ.7)
     . WRITE(16,930) 4,"# M_HC^2<1"
      IF(IFAIL.EQ.8)
     . WRITE(16,930) 4,"# Negative sfermion mass squared"
      IF(IFAIL.EQ.11)
     . WRITE(16,930) 4,"# Integration problem in RGES"
      IF(IFAIL.EQ.12)
     . WRITE(16,930) 4,"# Integration problem in RGESOFT"
      IF(IFAIL.EQ.13)
     . WRITE(16,930) 4,"# Convergence Problem"

      WRITE(16,899) "# Input parameters"
      WRITE(16,899) "BLOCK MODSEL"
      WRITE(16,921) 3,1,"NMSSM particle content"
      WRITE(16,921) 1,0,"IMOD"
      WRITE(16,921) 10,0,"ISCAN"
      WRITE(16,921) 9,OMGFLAG,"Call micrOmegas"
      WRITE(16,921) 8,PFLAG,"Precision for Higgs masses"
      WRITE(16,921) 13,NMSFLAG,"Sparticle decays via NMSDECAY"
      WRITE(16,921) 14,VFLAG,"H-> VV,VV*,(V*V*)"
      WRITE(16,921) 15,MOFLAG,"Precision for micromegas"
      WRITE(16,921) 16,OUTFLAG,"Extra BLOCK's yes/no"
      WRITE(16,921) 23,MWFLAG,"Delta_MW"
      IF(CFLAG(6).NE.0)WRITE(16,921) 24,CFLAG(6),"LHC SUSY constraints"

      WRITE(16,899) "BLOCK SMINPUTS"
      WRITE(16,901) 1,1d0/ALEMMZ,"ALPHA_EM^-1(MZ)"
      WRITE(16,901) 2,GF,"GF"
      WRITE(16,901) 3,ALSMZ,"ALPHA_S(MZ)"
      WRITE(16,901) 4,MZ,"MZ"
      WRITE(16,901) 5,MB,"MB(MB)"
      WRITE(16,901) 6,MT,"MTOP (POLE MASS)"
      WRITE(16,901) 7,MTAU,"MTAU"
      WRITE(16,899) "# SMINPUTS Beyond SLHA:"
      WRITE(16,906) "MW:",MW
      WRITE(16,906) "MS:",MS
      WRITE(16,906) "MC:",MC
      WRITE(16,906) "VUS:",VUS
      WRITE(16,906) "VCB:",VCB
      WRITE(16,906) "VUB:",VUB

      WRITE(16,899) "BLOCK MINPAR"
      IF(Q2FIX.EQ.1)WRITE(16,901) 0,DSQRT(Q2),"REN. SCALE"
      WRITE(16,901) 3,TANB,"TANBETA(MZ)"

      WRITE(16,899) "BLOCK EXTPAR"
      WRITE(16,901) 1,PAR(20), "M1"
      WRITE(16,901) 2,PAR(21), "M2"
      WRITE(16,901) 3,PAR(22), "M3"
      WRITE(16,901) 11,PAR(12), "ATOP"
      WRITE(16,901) 12,PAR(13), "ABOTTOM"
      WRITE(16,901) 13,PAR(14), "ATAU"
      WRITE(16,901) 16,PAR(25), "AMUON"

      WRITE(16,901) 31,DSQRT(PAR(18)),"LEFT SELECTRON"
      WRITE(16,901) 32,DSQRT(PAR(18)),"LEFT SMUON"
      WRITE(16,901) 33,DSQRT(PAR(10)),"LEFT STAU"

      WRITE(16,901) 34,DSQRT(PAR(19)),"RIGHT SELECTRON"
      WRITE(16,901) 35,DSQRT(PAR(19)),"RIGHT SMUON"
      WRITE(16,901) 36,DSQRT(PAR(11)),"RIGHT STAU"

      WRITE(16,901) 41,DSQRT(PAR(15)),"LEFT 1ST GEN. SQUARKS"
      WRITE(16,901) 42,DSQRT(PAR(15)),"LEFT 2ND GEN. SQUARKS"
      WRITE(16,901) 43,DSQRT(PAR(7)),"LEFT 3RD GEN. SQUARKS"

      WRITE(16,901) 44,DSQRT(PAR(16)),"RIGHT U-SQUARKS"
      WRITE(16,901) 45,DSQRT(PAR(16)),"RIGHT C-SQUARKS"
      WRITE(16,901) 46,DSQRT(PAR(8)),"RIGHT T-SQUARKS"

      WRITE(16,901) 47,DSQRT(PAR(17)),"RIGHT D-SQUARKS"
      WRITE(16,901) 48,DSQRT(PAR(17)),"RIGHT S-SQUARKS"
      WRITE(16,901) 49,DSQRT(PAR(9)),"RIGHT B-SQUARKS"

      WRITE(16,901) 61,PAR(1),"LAMBDA"
      WRITE(16,901) 62,PAR(2),"KAPPA"
      IF(MOD(MAFLAG,3).NE.1)THEN
       WRITE(16,901) 63,PAR(5),"ALAMBDA"
      ELSE
       WRITE(16,920) 63,PAR(5),"ALAMBDA"
      ENDIF
      IF(MAFLAG/3.NE.1)THEN
       WRITE(16,901) 64,PAR(6),"AKAPPA"
      ELSE
       WRITE(16,920) 64,PAR(6),"AKAPPA"
      ENDIF
      WRITE(16,901) 65,PAR(4),"MUEFF"
      IF(MOD(MAFLAG,3).NE.2)THEN
       IF(XIFSUSY.NE.0d0)
     .  WRITE(16,901) 66,XIFSUSY,"XIF"
      ELSE
       WRITE(16,920) 66,XIFSUSY,"XIF"
      ENDIF
      IF(MAFLAG/3.NE.2)THEN
       IF(XISSUSY.NE.0d0)
     .  WRITE(16,901) 67,XISSUSY,"XIS"
      ELSE
       WRITE(16,920) 67,XISSUSY,"XIS"
      ENDIF
      IF(MUPSUSY.NE.0d0)
     .  WRITE(16,901) 68,MUPSUSY,"MUP"
      IF(MSPSUSY.NE.0d0)
     .  WRITE(16,901) 69,MSPSUSY,"MSP"
      IF(M3HSUSY.NE.0d0)
     .  WRITE(16,901) 72,M3HSUSY,"M3H"
      IF(MOD(MAFLAG,3).NE.0)THEN
       WRITE(16,901) 124,PAR(23),"MA AT QSTSB"
      ELSE
       WRITE(16,920) 124,PAR(23),"MA AT QSTSB"
      ENDIF
      IF(MAFLAG/3.NE.0)THEN
       WRITE(16,901) 125,PAR(24),"MP AT QSTSB"
      ELSE
       WRITE(16,920) 125,PAR(24),"MP AT QSTSB"
      ENDIF

      IF(IFAIL.GT.0.AND.IFAIL.LT.10) RETURN

      WRITE(16,899) "# "
      WRITE(16,899) "BLOCK MASS   # Mass spectrum "
      WRITE(16,899) "#  PDG Code     mass             particle "
      WRITE(16,902) 5,MB,"MB(MB)"
      WRITE(16,902) 6,MT,"MTOP (POLE MASS)"
      WRITE(16,902) 15,MTAU,"MTAU"
      WRITE(16,902) 23,MZ,"MZ"
      IF(MWFLAG.NE.0)THEN
       WRITE(16,902) 24,MWNMSSM,"MW incl. Delta_MW"
      ELSE
       WRITE(16,902) 24,MW,"MW"
      ENDIF
      WRITE(16,902) PDGH1,SMASS(1),"lightest neutral scalar"
      WRITE(16,902) PDGH2,SMASS(2),"second neutral scalar"
      WRITE(16,902) 45,SMASS(3),"third neutral scalar"
      WRITE(16,902) 36,AMASS(1),"lightest pseudoscalar"
      WRITE(16,902) 46,AMASS(2),"second pseudoscalar"
      WRITE(16,902) 37,CMASS,"charged Higgs"
      WRITE(16,902) 1000001,MDL," ~d_L"
      WRITE(16,902) 2000001,MDR," ~d_R"
      WRITE(16,902) 1000002,MUL," ~u_L"
      WRITE(16,902) 2000002,MUR," ~u_R"
      WRITE(16,902) 1000003,MDL," ~s_L"
      WRITE(16,902) 2000003,MDR," ~s_R"
      WRITE(16,902) 1000004,MUL," ~c_L"
      WRITE(16,902) 2000004,MUR," ~c_R"
      WRITE(16,902) 1000005,MSB1," ~b_1"
      WRITE(16,902) 2000005,MSB2," ~b_2"
      WRITE(16,902) 1000006,MST1," ~t_1"
      WRITE(16,902) 2000006,MST2," ~t_2"
      WRITE(16,902) 1000011,MLL," ~e_L"
      WRITE(16,902) 2000011,MLR," ~e_R"
      WRITE(16,902) 1000012,MNL," ~nue_L"
      WRITE(16,902) 1000013,MSMU1," ~mu_L"
      WRITE(16,902) 2000013,MSMU2," ~mu_R"
      WRITE(16,902) 1000014,MNL," ~numu_L"
      WRITE(16,902) 1000015,MSL1," ~tau_1"
      WRITE(16,902) 2000015,MSL2," ~tau_2"
      WRITE(16,902) 1000016,MSNT," ~nutau_L"
      WRITE(16,902) 1000021,MGL," ~g"
      WRITE(16,902) 1000022,MNEU(1),"neutralino(1)"
      WRITE(16,902) 1000023,MNEU(2),"neutralino(2)"
      WRITE(16,902) 1000025,MNEU(3),"neutralino(3)"
      WRITE(16,902) 1000035,MNEU(4),"neutralino(4)"
      WRITE(16,902) 1000045,MNEU(5),"neutralino(5)"
      WRITE(16,902) 1000024,MCHA(1),"chargino(1)"
      WRITE(16,902) 1000037,MCHA(2),"chargino(2)"
      IF(GRFLAG.NE.0)WRITE(16,902) 1000039,M32,"gravitino"
      WRITE(16,899) "# "

       IF(OUTFLAG.NE.1) THEN

      IF(CFLAG(5).NE.0)THEN
       WRITE(16,899) "BLOCK TRILEP"
       WRITE(16,907) "# ",xsectot,
     .    "[pb]    Higgsino Xsect*BR into trileptons"
       WRITE(16,907) "# ",limtrilep,
     .    "[pb]    Upper limit from CMS, 1801.03957, Fig. 7 & 8a"
       WRITE(16,899) "# "
      ENDIF

      IF(CFLAG(6).GT.0)THEN
       WRITE(16,899) "BLOCK SMODELS"
       WRITE(16,111) "# Analysis: ",CHAN,"   Rmax: ",R
       WRITE(16,899) "# "
      ENDIF

      IF(PFLAG.GT.2)THEN
      WRITE(16,899) "BLOCK DMASS"
      WRITE(16,899)
     .   "# Mass uncertainties from QSTSB -> QSTSB*4 or QSTSB/4: "
      WRITE(16,927) DSMASS(1),"DM(lightest neutral scalar)"
      WRITE(16,927) DSMASS(2),"DM(second neutral scalar)"
      WRITE(16,927) DSMASS(3),"DM(third neutral scalar)"
      WRITE(16,927) DAMASS(1),"DM(lightest pseudoscalar)"
      WRITE(16,927) DAMASS(2),"DM(second pseudoscalar)"
      WRITE(16,927) DCMASS,"DM(charged Higgs)"
      WRITE(16,899) 
     .         "# If 300000.00: Delta(QSTSB) generates Msquark^2 < 0 or 
     .Mhiggs^2 < 0"
      WRITE(16,899) "# For Higgs mass precision = 0:"
      WRITE(16,899) 
     .    "#  2-loop correction assumes degenerate Stops and Sbottoms"
      WRITE(16,899) 
     .    "#  if not true, uncertainties are potentially underestimated"
      WRITE(16,899) "# "
      ENDIF

      WRITE(16,899) "# Low energy observables"
      WRITE(16,899) "BLOCK LOWEN"
      WRITE(16,899)
     .   "# Exp. 2 Sigma: 3.02E-4 < BR(b -> s gamma) < 3.62E-4:"
      WRITE(16,901) 1,BRSG,"BR(b -> s gamma)"
      WRITE(16,901) 11,BRSGMAX,"(BR(b -> s gamma)+Theor.Err.)"
      WRITE(16,901) 12,BRSGMIN,"(BR(b -> s gamma)-Theor.Err.)"
      WRITE(16,899) "# Exp. 2 Sigma: 5.027E-1 < Delta M_d < 5.103E-1:"
      WRITE(16,901) 2,DMD,"Delta M_d in ps^-1"
      WRITE(16,901) 21,DMdmax,"Delta M_d +Theor.Err."
      WRITE(16,901) 22,DMdmin,"Delta M_d -Theor.Err."
      WRITE(16,899) 
     .   "# Exp. 2 Sigma: 1.7715E+1 < Delta Ms < 1.7799E+1:"
      WRITE(16,901) 3,DMS,"Delta M_s in ps^-1"
      WRITE(16,901) 31,DMsmax,"Delta M_s +Theor.Err."
      WRITE(16,901) 32,DMsmin,"Delta M_s -Theor.Err."
      WRITE(16,899) "# Exp. 2 Sigma: 1.7E-9 < BR(Bs->mu+mu-) < 4.5E-9:"
      WRITE(16,901) 4,BRBMUMU,"BR(Bs -> mu+mu-)"
      WRITE(16,901) 41,BRBMUMUmax,"BR(Bs -> mu+mu-)+Theor.Err."
      WRITE(16,901) 42,BRBMUMUmin,"BR(Bs -> mu+mu-)-Theor.Err."
      WRITE(16,899) 
     .   "# Exp. 2 Sigma: 0.78E-4 < BR(B+ > tau+ + nu_tau) < 1.44E-4:"
      WRITE(16,901) 5,BRBtaunu,"BR(B+ -> tau+ + nu_tau)"
      WRITE(16,901) 51,BRBtaunumax,
     .   "BR(B+ -> tau+ + nu_tau) + Theor.Err."
      WRITE(16,901) 52,BRBtaunumin,
     .   "BR(B+ -> tau+ + nu_tau) - Theor.Err."
      WRITE(16,899) 
     .   "# Exp. 2 Sigma: 0.84E-6 < BR(B->Xs mu+mu-)low < 2.32E-6:"
      WRITE(16,907) "# ",BRBSll,"   BR(B->Xs mu+mu-)low"
      WRITE(16,907) "# ",BRBSllmax,"   BR(B->Xs mu+mu-)low+Theor.Err."
      WRITE(16,907) "# ",BRBSllmin,"   BR(B->Xs mu+mu-)low-Theor.Err."
      WRITE(16,899) 
     .   "# Exp. 2 Sigma: 2.8E-7 < BR(B->Xs mu+mu-)high < 6.8E-7:"
      WRITE(16,907) "# ",BRBShll,"   BR(B->Xs mu+mu-)high"
      WRITE(16,907) "# ",BRBShllmax,"   BR(B->Xs mu+mu-)high+Theor.Err."
      WRITE(16,907) "# ",BRBShllmin,"   BR(B->Xs mu+mu-)high-Theor.Err."
      WRITE(16,899)
     .   "# Exp. 2 Sigma: 2.7E-6 < BR(b -> d gamma) < 2.55E-5:"
      WRITE(16,907) "# ",BRDG,"   BR(b -> d gamma)"
      WRITE(16,907) "# ",BRDGMAX,"   BR(b -> d gamma)+Theor.Err.)"
      WRITE(16,907) "# ",BRDGMIN,"   BR(b -> d gamma)-Theor.Err.)"
      WRITE(16,899) 
     .   "# Exp. 2 Sigma: 1.1E-10 < BR(Bd->mu+mu-) < 7.1E-10:"
      WRITE(16,907) "# ",BRBdMUMU,"   BR(Bd -> mu+mu-)"
      WRITE(16,907) "# ",BRBdMUMUmax,"   BR(Bd -> mu+mu-)+Theor.Err."
      WRITE(16,907) "# ",BRBdMUMUmin,"   BR(Bd -> mu+mu-)-Theor.Err."
      WRITE(16,899) 
     .   "# Exp. 2 Sigma: BR(B-> Xs nu nubar) < 6.4E-4:"
      WRITE(16,907) "# ",BRBXsnunu,"   BR(B-> Xs nu nubar)"
      WRITE(16,907)
     .   "# ",BRBXsnunumax,"   BR(B-> Xs nu nubar)+Theor.Err."
      WRITE(16,907)
     .   "# ",BRBXsnunumin,"   BR(B-> Xs nu nubar)-Theor.Err."
      WRITE(16,899) 
     .   "# Exp. 2 Sigma: BR(B+-> K+ nu nubar) < 1.6E-5:"
      WRITE(16,907) "# ",BRBpKpnunu,"   BR(B+-> K+ nu nubar)"
      WRITE(16,907)
     .   "# ",BRBpKpnunumax,"   BR(B+-> K+ nu nubar)+Theor.Err."
      WRITE(16,907)
     .   "# ",BRBpKpnunumin,"   BR(B+-> K+ nu nubar)-Theor.Err."
      WRITE(16,899) 
     .   "# Exp. 2 Sigma: BR(B-> Ks nu nubar) < 5.5E-5:"
      WRITE(16,907) "# ",BRBKsnunu,"   BR(B-> Ks nu nubar)"
      WRITE(16,907)
     .   "# ",BRBKsnunumax,"   BR(B-> Ks nu nubar)+Theor.Err."
      WRITE(16,907)
     .   "# ",BRBKsnunumin,"   BR(B-> Ks nu nubar)-Theor.Err."
      WRITE(16,899) "# RD = BR[B+ -> D tau+ nu_tau]/BR[B+ -> D l+ nu_l]"
      WRITE(16,899) 
     .   "# Exp. 2 Sigma: 2.80E-1 < RD < 4.00E-1:"
      WRITE(16,907) "# ",RD_taul,"   RD"
      WRITE(16,907) "# ",RD_taulmax,"   RD+Theor.Err."
      WRITE(16,907) "# ",RD_taulmin,"   RD-Theor.Err."
      WRITE(16,899)
     .    "# RD* = BR[B+ -> D* tau+ nu_tau]/BR[B+ -> D* l+ nu_l]"
      WRITE(16,899) 
     .   "# Exp. 2 Sigma: 2.67E-1 < RD* < 3.23E-1:"
      WRITE(16,907) "# ",RDs_taul,"   RD*"
      WRITE(16,907) "# ",RDs_taulmax,"   RD*+Theor.Err."
      WRITE(16,907) "# ",RDs_taulmin,"   RD*-Theor.Err."
      WRITE(16,899) 
     .   "# Exp. 2 Sigma: BR(K+-> Pi+ nu nubar) < 4.03E-10:"
      WRITE(16,907) "# ",BRKp_Pipnunub,"   BR(K+-> Pi+ nu nubar)"
      WRITE(16,907)
     .   "# ",BRKp_Pipnunubmax,"   BR(K+-> Pi+ nu nubar)+Theor.Err."
      WRITE(16,907)
     .   "# ",BRKp_Pipnunubmin,"   BR(K+-> Pi+ nu nubar)-Theor.Err."
      WRITE(16,899) 
     .   "# Exp. 2 Sigma: BR(KL-> Pi0 nu nubar) < 2.6E-8:"
      WRITE(16,907) "# ",BRKL_Pi0nunub,"   BR(KL-> Pi0 nu nubar)"
      WRITE(16,907)
     .   "# ",BRKL_Pi0nunubmax,"   BR(KL-> Pi0 nu nubar)+Theor.Err."
      WRITE(16,907)
     .   "# ",BRKL_Pi0nunubmin,"   BR(KL-> Pi0 nu nubar)-Theor.Err."
      WRITE(16,899) 
     .   "# Exp. 2 Sigma: 5.275E-3 < DMK < 5.311E-3:"
      WRITE(16,907) "# ",DMK,"   DMK in ps-1"
      WRITE(16,907)
     .   "# ",DMKmax,"   DMK+Theor.Err."
      WRITE(16,907)
     .   "# ",DMKmin,"   DMK-Theor.Err."
      WRITE(16,899) 
     .   "# Exp. 2 Sigma: 2.206E-3 < epsK < 2.250E-3:"
      WRITE(16,907) "# ",epsK,"   epsK"
      WRITE(16,907)
     .   "# ",epsKmax,"   epsK+Theor.Err."
      WRITE(16,907)
     .   "# ",epsKmin,"   epsK-Theor.Err."
      WRITE(16,899) "# " 
      WRITE(16,899) "# BSM contr. to the muon anomalous magn. moment:"
      WRITE(16,901) 6,delmagmu,"Del_a_mu"
      WRITE(16,901) 61,amuthmax,"Del_a_mu + Theor.Err."
      WRITE(16,901) 62,amuthmin,"Del_a_mu - Theor.Err."
      WRITE(16,907) "# Minimal Exp.-SM (2 sigma):",damumin
      WRITE(16,907) "# Maximal Exp.-SM (2 sigma):",damumax

      IF(OMGFLAG.NE.0)THEN
        WRITE(16,899) "# "
        IF(OMGFLAG.GT.0)WRITE(16,911)
     .   "# Omega h^2 (allowed:",OMGMIN," < Omega h^2 <",OMGMAX,"):"
        IF(OMGFLAG.LT.0)WRITE(16,911)
     .   "# Omega h^2 (allowed: Omega h^2 <",OMGMAX,"):"
        IF(OMG.EQ.0d0)THEN
          WRITE(16,899) "# Cannot compute Omega h^2"
        ELSEIF(OMG.EQ.-1d0)THEN
          WRITE(16,899)
     .      "# Charged LSP"
        ELSEIF(OMG.LE.-2d0)THEN
          WRITE(16,899) "# Problem in micrOMEGAs"
        ELSE
          WRITE(16,901) 10,OMG,"Omega h^2"
          omg_=printChannels(Xf,1d-3,1d-4,1,16)
        ENDIF
      ENDIF

      IF(IABS(OMGFLAG).EQ.2 .OR. IABS(OMGFLAG).EQ.4)THEN

* New June 2019
       XSMAX=MIN(PandaX_SI(DABS(MNEU(1))),LUX_SI(DABS(MNEU(1))),
     .       XENON_SI(DABS(MNEU(1))),CRESST_SI(DABS(MNEU(1))),
     .       DarkSide50_SI(DABS(MNEU(1))),LZ_SI(DABS(MNEU(1))))
       IF(XSMAX.EQ.1d99)THEN
        WRITE(16,899)"# sigma(p)_SI (all values allowed)"
       ELSEIF(OMGFLAG.LT.0 .AND. OMG.LT.OMGMAX)THEN
        WRITE(16,928)
     .   "# allowed sigma(p)_SI after rescaling by OMGMAX/OMG=",
     .   OMGMAX/OMG,": < ",OMGMAX/OMG*XSMAX
       ELSE
        WRITE(16,907)
     .   "# allowed sigma(p)_SI assuming Omega h^2 = 0.1187: <",XSMAX
       ENDIF
       WRITE(16,901) 20,DABS(CSPSI),"sigma_p^SI"

       XSMAX=MIN(LUX_SDn(DABS(MNEU(1))),XENON_SDn(DABS(MNEU(1))),
     .       PandaX_SDn(DABS(MNEU(1))),LZ_SDn(DABS(MNEU(1))))
       IF(XSMAX.EQ.1d99)THEN
        WRITE(16,899)"# sigma(n)_SD (all values allowed)"
       ELSEIF(OMGFLAG.LT.0 .AND. OMG.LT.OMGMAX)THEN
        WRITE(16,928)
     .   "# allowed sigma(n)_SD after rescaling by OMGMAX/OMG=",
     .   OMGMAX/OMG,": < ",OMGMAX/OMG*XSMAX
       ELSE
        WRITE(16,907)
     .   "# allowed sigma(n)_SD assuming Omega h^2 = 0.1187: <",XSMAX
       ENDIF
       WRITE(16,901) 30,DABS(CSNSD),"sigma_n^SD"

       XSMAX=MIN(LUX_SDp(DABS(MNEU(1))),XENON_SDp(DABS(MNEU(1))),
     .        PICO60_SDp(DABS(MNEU(1))),LZ_SDp(DABS(MNEU(1))))
       IF(XSMAX.EQ.1d99)THEN
        WRITE(16,899)"# sigma(p)_SD (all values allowed)"
       ELSEIF(OMGFLAG.LT.0 .AND. OMG.LT.OMGMAX)THEN
        WRITE(16,928)
     .   "# allowed sigma(p)_SD after rescaling by OMGMAX/OMG=",
     .   OMGMAX/OMG,": < ",OMGMAX/OMG*XSMAX
       ELSE
        WRITE(16,907)
     .   "# allowed sigma(p)_SD assuming Omega h^2 = 0.1187: <",XSMAX
       ENDIF
       WRITE(16,901) 40,DABS(CSPSD),"sigma_p^SD"

        WRITE(16,915)"# values used for sigma_piN,sigma_S",
     .  " (strange content of the proton)"
        WRITE(16,901) 50,sigmapiN,"sigma_piN"
        WRITE(16,901) 60,sigmaS,"sigma_S"
      ENDIF
      WRITE(16,899) "# "
*  From IF(OUTFLAG.NE.1)
       ENDIF

      WRITE(16,907) "BLOCK HMIX Q=",DSQRT(QSTSB),
     .    " # (STOP/SBOTTOM MASSES)"
      WRITE(16,901) 1,MUQ,"MUEFF"
      WRITE(16,901) 2,TANBQ,"TAN(BETA)"
      WRITE(16,901) 3,DSQRT(2d0*(H1Q**2+H2Q**2)),"V(Q)"
      WRITE(16,901) 4,PAR(23)**2,"MA^2"
      WRITE(16,901) 5,PAR(24)**2,"MP^2"

      WRITE(16,899) "# "
      WRITE(16,899) "# 3*3 Higgs mixing"
      WRITE(16,899) "BLOCK NMHMIX"
      WRITE(16,903) 1,1,SCOMP(1,2),"S_(1,1)"
      WRITE(16,903) 1,2,SCOMP(1,1),"S_(1,2)"
      WRITE(16,903) 1,3,SCOMP(1,3),"S_(1,3)"
      WRITE(16,903) 2,1,SCOMP(2,2),"S_(2,1)"
      WRITE(16,903) 2,2,SCOMP(2,1),"S_(2,2)"
      WRITE(16,903) 2,3,SCOMP(2,3),"S_(2,3)"
      WRITE(16,903) 3,1,SCOMP(3,2),"S_(3,1)"
      WRITE(16,903) 3,2,SCOMP(3,1),"S_(3,2)"
      WRITE(16,903) 3,3,SCOMP(3,3),"S_(3,3)"

      WRITE(16,899) "# "
      WRITE(16,899) "# 3*3 Pseudoscalar Higgs mixing"
      WRITE(16,899) "BLOCK NMAMIX"
      WRITE(16,903) 1,1,SINB*PCOMP(1,1),"P_(1,1)"
      WRITE(16,903) 1,2,COSB*PCOMP(1,1),"P_(1,2)"
      WRITE(16,903) 1,3,PCOMP(1,2),"P_(1,3)"
      WRITE(16,903) 2,1,SINB*PCOMP(2,1),"P_(2,1)"
      WRITE(16,903) 2,2,COSB*PCOMP(2,1),"P_(2,2)"
      WRITE(16,903) 2,3,PCOMP(2,2),"P_(2,3)"

      SST=DSQRT(1-CST**2)
      SSB=DSQRT(1-CSB**2)
      SSL=DSQRT(1-CSL**2)

      WRITE(16,899) "# "
      WRITE(16,899) "# 3rd generation sfermion mixing"
      WRITE(16,899) "BLOCK STOPMIX  # Stop mixing matrix"
      WRITE(16,903) 1,1,CST,"Rst_(1,1)"
      WRITE(16,903) 1,2,SST,"Rst_(1,2)"
      WRITE(16,903) 2,1,-SST,"Rst_(2,1)"
      WRITE(16,903) 2,2,CST,"Rst_(2,2)"
      WRITE(16,899) "BLOCK SBOTMIX  # Sbottom mixing matrix"
      WRITE(16,903) 1,1,CSB,"Rsb_(1,1)"
      WRITE(16,903) 1,2,SSB,"Rsb_(1,2)"
      WRITE(16,903) 2,1,-SSB,"Rsb_(2,1)"
      WRITE(16,903) 2,2,CSB,"Rsb_(2,2)"
      WRITE(16,899) "BLOCK STAUMIX  # Stau mixing matrix"
      WRITE(16,903) 1,1,CSL,"Rsl_(1,1)"
      WRITE(16,903) 1,2,SSL,"Rsl_(1,2)"
      WRITE(16,903) 2,1,-SSL,"Rsl_(2,1)"
      WRITE(16,903) 2,2,CSL,"Rsl_(2,2)"

      WRITE(16,899) "# "
      WRITE(16,899) "# Gaugino-Higgsino mixing"
      WRITE(16,899) "BLOCK NMNMIX  # 5*5 Neutralino Mixing Matrix"
      WRITE(16,903) 1,1,NEU(1,1),"N_(1,1)"
      WRITE(16,903) 1,2,NEU(1,2),"N_(1,2)"
      WRITE(16,903) 1,3,NEU(1,4),"N_(1,3)"
      WRITE(16,903) 1,4,NEU(1,3),"N_(1,4)"
      WRITE(16,903) 1,5,NEU(1,5),"N_(1,5)"
      WRITE(16,903) 2,1,NEU(2,1),"N_(2,1)"
      WRITE(16,903) 2,2,NEU(2,2),"N_(2,2)"
      WRITE(16,903) 2,3,NEU(2,4),"N_(2,3)"
      WRITE(16,903) 2,4,NEU(2,3),"N_(2,4)"
      WRITE(16,903) 2,5,NEU(2,5),"N_(2,5)"
      WRITE(16,903) 3,1,NEU(3,1),"N_(3,1)"
      WRITE(16,903) 3,2,NEU(3,2),"N_(3,2)"
      WRITE(16,903) 3,3,NEU(3,4),"N_(3,3)"
      WRITE(16,903) 3,4,NEU(3,3),"N_(3,4)"
      WRITE(16,903) 3,5,NEU(3,5),"N_(3,5)"
      WRITE(16,903) 4,1,NEU(4,1),"N_(4,1)"
      WRITE(16,903) 4,2,NEU(4,2),"N_(4,2)"
      WRITE(16,903) 4,3,NEU(4,4),"N_(4,3)"
      WRITE(16,903) 4,4,NEU(4,3),"N_(4,4)"
      WRITE(16,903) 4,5,NEU(4,5),"N_(4,5)"
      WRITE(16,903) 5,1,NEU(5,1),"N_(5,1)"
      WRITE(16,903) 5,2,NEU(5,2),"N_(5,2)"
      WRITE(16,903) 5,3,NEU(5,4),"N_(5,3)"
      WRITE(16,903) 5,4,NEU(5,3),"N_(5,4)"
      WRITE(16,903) 5,5,NEU(5,5),"N_(5,5)"

      WRITE(16,899) "# "
      WRITE(16,899) "BLOCK UMIX  # Chargino U Mixing Matrix"
      WRITE(16,903) 1,1,U(1,1),"U_(1,1)"
      WRITE(16,903) 1,2,U(1,2),"U_(1,2)"
      WRITE(16,903) 2,1,U(2,1),"U_(2,1)"
      WRITE(16,903) 2,2,U(2,2),"U_(2,2)"

      WRITE(16,899) "# "
      WRITE(16,899) "BLOCK VMIX  # Chargino V Mixing Matrix"
      WRITE(16,903) 1,1,V(1,1),"V_(1,1)"
      WRITE(16,903) 1,2,V(1,2),"V_(1,2)"
      WRITE(16,903) 2,1,V(2,1),"V_(2,1)"
      WRITE(16,903) 2,2,V(2,2),"V_(2,2)"
      WRITE(16,899) "# "

       IF(OUTFLAG.NE.1) THEN
      WRITE(16,899) "# SM-Higgs reduced couplings"
      WRITE(16,899) "# (as compared to a SM Higgs with same mass)"
      WRITE(16,899) "# Allowed 2sigma intervals combining"
      WRITE(16,899) "# ATLAS 2207.00092 Fig.6, CMS 2207.00043, Fig.4a"
      WRITE(16,899) "# 0.826 < CU < 1.129 " ! top quark
      WRITE(16,899) "# 0.744 < CB < 1.111 " ! b quark
      WRITE(16,899) "# 0.820 < CL < 1.033 " ! taus
      WRITE(16,899) "# 0.923 < CV < 1.000 " ! electroweak couplings
      WRITE(16,899) "# 0.834 < CJ < 1.044 " ! gluons
      WRITE(16,899) "# 0.946 < CG < 1.143 " ! photons
      WRITE(16,899) "BLOCK REDCOUP"
      WRITE(16,899) "# H1"
      WRITE(16,903) 1,1,CU(1),"U-type fermions"
      WRITE(16,903) 1,2,CD(1),"D-type fermions"
      WRITE(16,903) 1,3,CB(1),"b-quarks"
      WRITE(16,903) 1,4,CL(1),"taus"
      WRITE(16,903) 1,5,CV(1),"W,Z bosons"
      WRITE(16,903) 1,6,CJ(1),"Gluons"
      WRITE(16,903) 1,7,CG(1),"Photons"
      WRITE(16,899) "# H2"
      WRITE(16,903) 2,1,CU(2),"U-type fermions"
      WRITE(16,903) 2,2,CD(2),"D-type fermions"
      WRITE(16,903) 2,3,CB(2),"b-quarks"
      WRITE(16,903) 2,4,CL(2),"taus"
      WRITE(16,903) 2,5,CV(2),"W,Z bosons"
      WRITE(16,903) 2,6,CJ(2),"Gluons"
      WRITE(16,903) 2,7,CG(2),"Photons"
      WRITE(16,899) "# H3"
      WRITE(16,903) 3,1,CU(3),"U-type fermions"
      WRITE(16,903) 3,2,CD(3),"D-type fermions"
      WRITE(16,903) 3,3,CB(3),"b-quarks"
      WRITE(16,903) 3,4,CL(3),"taus"
      WRITE(16,903) 3,5,CV(3),"W,Z bosons"
      WRITE(16,903) 3,6,CJ(3),"Gluons"
      WRITE(16,903) 3,7,CG(3),"Photons"
      WRITE(16,899) "# A1"
      WRITE(16,903) 4,1,CU(4),"U-type fermions"
      WRITE(16,903) 4,2,CD(4),"D-type fermions"
      WRITE(16,903) 4,3,CB(4),"b-quarks"
      WRITE(16,903) 4,4,CL(4),"taus"
      WRITE(16,903) 4,5,0d0,"W,Z bosons"
      WRITE(16,903) 4,6,CJ(4),"Gluons"
      WRITE(16,903) 4,7,CG(4),"Photons"
      WRITE(16,899) "# A2"
      WRITE(16,903) 5,1,CU(5),"U-type fermions"
      WRITE(16,903) 5,2,CD(5),"D-type fermions"
      WRITE(16,903) 5,3,CB(5),"b-quarks"
      WRITE(16,903) 5,4,CL(5),"taus"
      WRITE(16,903) 5,5,0d0,"W,Z bosons"
      WRITE(16,903) 5,6,CJ(5),"Gluons"
      WRITE(16,903) 5,7,CG(5),"Photons"

      WRITE(16,'(a)')'#'
      WRITE(16,'(a)')'Block HiggsCouplingsBosons'
      WRITE(16,'(G16.6,4I6,a)')CV(1),        3,PDGH1,24,24,
     .     ' # Higgs(1)-W-W reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')CV(2),        3,PDGH2,24,24,
     .     ' # Higgs(2)-W-W reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')CV(3),        3,45,24,24,
     .     ' # Higgs(3)-W-W reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')0d0,        3,36,24,24,
     .     ' # CP-odd Higgs(1)-W-W reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')0d0,        3,46,24,24,
     .     ' # CP-odd Higgs(2)-W-W reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')CV(1),        3,PDGH1,23,23,
     .     ' # Higgs(1)-Z-Z reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')CV(2),        3,PDGH2,23,23,
     .     ' # Higgs(2)-Z-Z reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')CV(3),        3,45,23,23,
     .     ' # Higgs(3)-Z-Z reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')0d0,        3,36,23,23,
     .     ' # CP-odd Higgs(1)-Z-Z reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')0d0,        3,46,23,23,
     .     ' # CP-odd Higgs(2)-Z-Z reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')CJ(1),        3,PDGH1,21,21,
     .  ' # Higgs(1)-gluon-gluon reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')CJ(2),        3,PDGH2,21,21,
     .  ' # Higgs(2)-gluon-gluon reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')CJ(3),        3,45,21,21,
     .  ' # Higgs(3)-gluon-gluon reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')CJ(4),        3,36,21,21,
     .  ' # CP-odd Higgs(1)-gluon-gluon reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')CJ(5),        3,46,21,21,
     .  ' # CP-odd Higgs(2)-gluon-gluon reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')CG(1),        3,PDGH1,22,22,
     .  ' # Higgs(1)-gamma-gamma reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')CG(2),        3,PDGH2,22,22,
     .  ' # Higgs(2)-gamma-gamma reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')CG(3),        3,45,22,22,
     .  ' # Higgs(3)-gamma-gamma reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')CG(4),        3,36,22,22,
     .  ' # CP-odd Higgs(1)-gamma-gamma reduced coupling'
      WRITE(16,'(G16.6,4I6,a)')CG(5),        3,46,22,22,
     .  ' # CP-odd Higgs(2)-gamma-gamma reduced coupling'
      WRITE(16,'(a)')'#'
      WRITE(16,'(a)')'Block HiggsCouplingsFermions'
      WRITE(16,'(a)')'#     Scalar       Pseudoscalar'
      WRITE(16,'(G16.6,G16.6,4I6,a)')CB(1),0d0, 3,PDGH1,5,5,
     .   ' # Higgs(1)-b-b red. coupling'
      WRITE(16,'(G16.6,G16.6,4I6,a)')CB(2),0d0, 3,PDGH2,5,5,
     .   ' # Higgs(2)-b-b red. coupling'
      WRITE(16,'(G16.6,G16.6,4I6,a)')CB(3),0d0, 3,45,5,5,
     .   ' # Higgs(3)-b-b red. coupling'
      WRITE(16,'(G16.6,G16.6,4I6,a)')0d0, CB(4), 3,36,5,5,
     .   ' # CP-odd Higgs(1)-b-b red. coupling'
      WRITE(16,'(G16.6,G16.6,4I6,a)')0d0, CB(5),3,46,5,5,
     .   ' # CP-odd Higgs(2)-b-b red. coupling'
      WRITE(16,'(G16.6,G16.6,4I6,a)')CU(1),0d0,3,PDGH1,6,6,
     .   ' # Higgs(1)-top-top red. coupling'
      WRITE(16,'(G16.6,G16.6,4I6,a)')CU(2),0d0, 3,PDGH2,6,6,
     .   ' # Higgs(2)-top-top red. coupling'
      WRITE(16,'(G16.6,G16.6,4I6,a)')CU(3),0d0, 3,45,6,6,
     .   ' # Higgs(3)-top-top red. coupling'
      WRITE(16,'(G16.6,G16.6,4I6,a)')0d0, CU(4), 3,36,6,6,
     .   ' # CP-odd Higgs(1)-top-top red. coupling'
      WRITE(16,'(G16.6,G16.6,4I6,a)')0d0, CU(5), 3,46,6,6,
     .   ' # CP-odd Higgs(2)-top-top red. coupling'
      WRITE(16,'(G16.6,G16.6,4I6,a)')CL(1),0d0, 3,PDGH1,15,15,
     .   ' # Higgs(1)-tau-tau red. coupling'
      WRITE(16,'(G16.6,G16.6,4I6,a)')CL(2),0d0, 3,PDGH2,15,15,
     .   ' # Higgs(2)-tau-tau red. coupling'
      WRITE(16,'(G16.6,G16.6,4I6,a)')CL(3),0d0, 3,45,15,15,
     .   ' # Higgs(3)-tau-tau red. coupling'
      WRITE(16,'(G16.6,G16.6,4I6,a)')0d0, CL(4), 3,36,15,15,
     .   ' # CP-odd Higgs(1)-tau-tau red. coupling'
      WRITE(16,'(G16.6,G16.6,4I6,a)')0d0, CL(5), 3,46,15,15,
     .   ' # CP-odd Higgs(2)-tau-tau red. coupling'
      WRITE(16,899) "# "
*  From IF(OUTFLAG.NE.1)
       ENDIF

      WRITE(16,899) "# GAUGE AND YUKAWA COUPLINGS AT THE SUSY SCALE"
      WRITE(16,907) "BLOCK GAUGE Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(16,901) 1,DSQRT(G1S),"g1(Q,DR_bar)"
      WRITE(16,901) 2,DSQRT(G2S),"g2(Q,DR_bar)"
      WRITE(16,901) 3,DSQRT(G3S),"g3(Q,DR_bar)"

      WRITE(16,907) "BLOCK YU Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(16,903) 3,3,HTOPS,"HTOP(Q,DR_bar)"
      WRITE(16,907) "BLOCK YD Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(16,903) 3,3,HBOTS,"HBOT(Q,DR_bar)"
      WRITE(16,907) "BLOCK YE Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(16,903) 3,3,HTAUS,"HTAU(Q,DR_bar)"

      WRITE(16,899) "# "
      WRITE(16,899) "# SOFT TRILINEAR COUPLINGS AT THE SUSY SCALE"
      WRITE(16,899) "# (BOTH SLHA1 AND SLHA2 FORMAT)"
      WRITE(16,907) "BLOCK AU Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(16,903) 3,3,PAR(12),"ATOP"
      WRITE(16,907) "BLOCK TU Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(16,903) 3,3,PAR(12)*HTOPS,"TU"
      WRITE(16,907) "BLOCK AD Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(16,903) 3,3,PAR(13),"ABOT"
      WRITE(16,907) "BLOCK TD Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(16,903) 3,3,PAR(13)*HBOTS,"TD"
      WRITE(16,907) "BLOCK AE Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(16,903) 2,2,PAR(25),"AMUON"
      WRITE(16,903) 3,3,PAR(14),"ATAU"
      WRITE(16,907) "BLOCK TE Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(16,903) 2,2,PAR(25)*HTAUS*MMUON/MTAU,"TE"
      WRITE(16,903) 3,3,PAR(14)*HTAUS,"TE"

      WRITE(16,899) "# "
      WRITE(16,899) "# SOFT MASSES AT THE SUSY SCALE"
      WRITE(16,907) "BLOCK MSOFT Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(16,901) 1,PAR(20),"M1"
      WRITE(16,901) 2,PAR(21),"M2"
      WRITE(16,901) 3,PAR(22),"M3"
      WRITE(16,901) 21,MHDS,"M_HD^2"
      WRITE(16,901) 22,MHUS,"M_HU^2"
      WRITE(16,901) 31,PAR(18)/DSQRT(DABS(PAR(18))),"M_eL"
      WRITE(16,901) 32,PAR(18)/DSQRT(DABS(PAR(18))),"M_muL"
      WRITE(16,901) 33,PAR(10)/DSQRT(DABS(PAR(10))),"M_tauL"
      WRITE(16,901) 34,PAR(19)/DSQRT(DABS(PAR(19))),"M_eR"
      WRITE(16,901) 35,PAR(19)/DSQRT(DABS(PAR(19))),"M_muR"
      WRITE(16,901) 36,PAR(11)/DSQRT(DABS(PAR(11))),"M_tauR"
      WRITE(16,901) 41,PAR(15)/DSQRT(DABS(PAR(15))),"M_q1L"
      WRITE(16,901) 42,PAR(15)/DSQRT(DABS(PAR(15))),"M_q2L"
      WRITE(16,901) 43,PAR(7)/DSQRT(DABS(PAR(7))),"M_q3L"
      WRITE(16,901) 44,PAR(16)/DSQRT(DABS(PAR(16))),"M_uR"
      WRITE(16,901) 45,PAR(16)/DSQRT(DABS(PAR(16))),"M_cR"
      WRITE(16,901) 46,PAR(8)/DSQRT(DABS(PAR(8))),"M_tR"
      WRITE(16,901) 47,PAR(17)/DSQRT(DABS(PAR(17))),"M_dR"
      WRITE(16,901) 48,PAR(17)/DSQRT(DABS(PAR(17))),"M_sR"
      WRITE(16,901) 49,PAR(9)/DSQRT(DABS(PAR(9))),"M_bR"

      WRITE(16,899) "# "
      WRITE(16,899) "# NMSSM SPECIFIC PARAMETERS THE SUSY SCALE"
      WRITE(16,907) "BLOCK NMSSMRUN Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(16,901) 1,PAR(1),"LAMBDA(Q,DR_bar)"
      WRITE(16,901) 2,PAR(2),"KAPPA(Q,DR_bar)"
      WRITE(16,901) 3,PAR(5),"ALAMBDA"
      WRITE(16,901) 4,PAR(6),"AKAPPA"
      WRITE(16,901) 5,PAR(4),"MUEFF"
      WRITE(16,901) 6,XIFSUSY,"XIF"
      WRITE(16,901) 7,XISSUSY,"XIS"
      WRITE(16,901) 8,MUPSUSY,"MUP"
      WRITE(16,901) 9,MSPSUSY,"MSP"
      WRITE(16,901) 10,MSS,"MS^2"
      WRITE(16,920) 12,M3HSUSY,"M3H"
      WRITE(16,899) "# "
      WRITE(16,899) "BLOCK MSQ2  # Soft l.h. squark masses squared"
      WRITE(16,903) 1,1,PAR(15),"M_q1L"
      WRITE(16,903) 2,2,PAR(15),"M_q2L"
      WRITE(16,903) 3,3,PAR(7),"M_q3L"
      WRITE(16,899) "# "
      WRITE(16,899) "BLOCK MSU2  # Soft r.h. up-squark masses squared"
      WRITE(16,903) 1,1,PAR(16),"M_uR"
      WRITE(16,903) 2,2,PAR(16),"M_cR"
      WRITE(16,903) 3,3,PAR(8),"M_tR"
      WRITE(16,899) "# "
      WRITE(16,899) "BLOCK MSD2  # Soft r.h. down-squark masses squared"
      WRITE(16,903) 1,1,PAR(17),"M_dR"
      WRITE(16,903) 2,2,PAR(17),"M_sR"
      WRITE(16,903) 3,3,PAR(9),"M_bR"
      WRITE(16,899) "# "
      WRITE(16,899) "BLOCK MSL2  # Soft l.h. slepton masses squared"
      WRITE(16,903) 1,1,PAR(18),"M_eL"
      WRITE(16,903) 2,2,PAR(18),"M_muL"
      WRITE(16,903) 3,3,PAR(10),"M_tauL"
      WRITE(16,899) "# "
      WRITE(16,899) "BLOCK MSE2  # Soft r.h. slepton masses squared"
      WRITE(16,903) 1,1,PAR(19),"M_eR"
      WRITE(16,903) 2,2,PAR(19),"M_muR"
      WRITE(16,903) 3,3,PAR(11),"M_tauR"
      WRITE(16,899) "# "
      WRITE(16,899) "BLOCK USQMIX  # Elements of 6x6 up-squark matrix"
      WRITE(16,903) 1,1,1.d0,"R_u_11"
      WRITE(16,903) 2,2,1.d0,"R_u_22"
      WRITE(16,903) 3,3,CST,"R_u_33"
      WRITE(16,903) 3,6,SST,"R_u_36"
      WRITE(16,903) 4,4,1.d0,"R_u_44"
      WRITE(16,903) 5,5,1.d0,"R_u_55"
      WRITE(16,903) 6,6,CST,"R_u_66"
      WRITE(16,903) 6,3,-SST,"R_u_63"
      WRITE(16,899) "# "
      WRITE(16,899) "BLOCK DSQMIX  # Elements of 6x6 down-squark matrix"
      WRITE(16,903) 1,1,1.d0,"R_d_11"
      WRITE(16,903) 2,2,1.d0,"R_d_22"
      WRITE(16,903) 3,3,CSB,"R_d_33"
      WRITE(16,903) 3,6,SSB,"R_u_36"
      WRITE(16,903) 4,4,1.d0,"R_d_44"
      WRITE(16,903) 5,5,1.d0,"R_d_55"
      WRITE(16,903) 6,6,CSB,"R_d_66"
      WRITE(16,903) 6,3,-SSB,"R_u36"
      WRITE(16,899) "# "
      WRITE(16,899) "BLOCK SELMIX  # Elements of 6x6 ch. slepton matrix"
      WRITE(16,903) 1,1,1.d0,"R_e_11"
      WRITE(16,903) 2,2,1.d0,"R_e_22"
      WRITE(16,903) 3,3,CSL,"R_e_33"
      WRITE(16,903) 3,6,SSL,"R_e_36"
      WRITE(16,903) 4,4,1.d0,"R_e_44"
      WRITE(16,903) 5,5,1.d0,"R_e_55"
      WRITE(16,903) 6,6,CSL,"R_e_66"
      WRITE(16,903) 6,3,-SSL,"R_e_63"
      WRITE(16,899) "# "

       IF(OUTFLAG.NE.1) THEN
      WRITE(16,899) "# GAUGE AND YUKAWA COUPLINGS AT THE GUT SCALE"
      WRITE(16,907) "BLOCK GUTGAUGE Q=",MGUT," # (GUT SCALE)"
      WRITE(16,901) 1,DSQRT(5d0/3d0*G1GUT),
     .        "g1(MGUT,DR_bar), GUT normalization"
      WRITE(16,901) 2,DSQRT(G2GUT),"g2(MGUT,DR_bar)"
      WRITE(16,901) 3,DSQRT(G3GUT),"g3(MGUT,DR_bar)"
      WRITE(16,907) "BLOCK GUTYU Q=",MGUT," # (GUT SCALE)"
      WRITE(16,903) 3,3,HTOPGUT/DSQRT(DABS(HTOPGUT)),
     .            "HTOP(MGUT,DR_bar)"
      WRITE(16,907) "BLOCK GUTYD Q=",MGUT," # (GUT SCALE)"
      WRITE(16,903) 3,3,HBOTGUT/DSQRT(DABS(HBOTGUT)),
     .            "HBOT(MGUT,DR_bar)"
      WRITE(16,907) "BLOCK GUTYE Q=",MGUT," # (GUT SCALE)"
      WRITE(16,903) 3,3,HTAUGUT/DSQRT(DABS(HTAUGUT)),
     .        "HTAU(MGUT,DR_bar)"

      WRITE(16,899) "# "
      WRITE(16,899) "# SOFT TRILINEAR COUPLINGS AT THE GUT SCALE"
      WRITE(16,907) "BLOCK GUTAU Q=",MGUT," # (GUT SCALE)"
      WRITE(16,903) 3,3,ATGUT,"ATOP"
      WRITE(16,907) "BLOCK GUTAD Q=",MGUT," # (GUT SCALE)"
      WRITE(16,903) 3,3,ABGUT,"ABOT"
      WRITE(16,907) "BLOCK GUTAE Q=",MGUT," # (GUT SCALE)"
      WRITE(16,903) 2,2,AMUGUT,"AMUON"
      WRITE(16,903) 3,3,ATAUGUT,"ATAU"

      WRITE(16,899) "# "
      WRITE(16,899) "# SOFT MASSES SQUARED AT THE GUT SCALE"
      WRITE(16,907) "BLOCK GUTMSOFT Q=",MGUT," # (GUT SCALE)"
      WRITE(16,901) 1,M1GUT,"M1"
      WRITE(16,901) 2,M2GUT,"M2"
      WRITE(16,901) 3,M3GUT,"M3"
      WRITE(16,901) 21,MHDGUT,"M_HD^2"
      WRITE(16,901) 22,MHUGUT,"M_HU^2"
      WRITE(16,901) 31,MLGUT/DSQRT(DABS(MLGUT)),"M_eL"
      WRITE(16,901) 32,MLGUT/DSQRT(DABS(MLGUT)),"M_muL"
      WRITE(16,901) 33,ML3GUT/DSQRT(DABS(ML3GUT)),"M_tauL"
      WRITE(16,901) 34,MEGUT/DSQRT(DABS(MEGUT)),"M_eR"
      WRITE(16,901) 35,MEGUT/DSQRT(DABS(MEGUT)),"M_muR"
      WRITE(16,901) 36,ME3GUT/DSQRT(DABS(ME3GUT)),"M_tauR"
      WRITE(16,901) 41,MQGUT/DSQRT(DABS(MQGUT)),"M_q1L"
      WRITE(16,901) 42,MQGUT/DSQRT(DABS(MQGUT)),"M_q2L"
      WRITE(16,901) 43,MQ3GUT/DSQRT(DABS(MQ3GUT)),"M_q3L"
      WRITE(16,901) 44,MUGUT/DSQRT(DABS(MUGUT)),"M_uR"
      WRITE(16,901) 45,MUGUT/DSQRT(DABS(MUGUT)),"M_cR"
      WRITE(16,901) 46,MU3GUT/DSQRT(DABS(MU3GUT)),"M_tR"
      WRITE(16,901) 47,MDGUT/DSQRT(DABS(MDGUT)),"M_dR"
      WRITE(16,901) 48,MDGUT/DSQRT(DABS(MDGUT)),"M_sR"
      WRITE(16,901) 49,MD3GUT/DSQRT(DABS(MD3GUT)),"M_bR"

      WRITE(16,899) "# "
      WRITE(16,899) "# NMSSM SPECIFIC PARAMETERS AT THE GUT SCALE"
      WRITE(16,907) "BLOCK GUTNMSSMRUN Q=",MGUT," # (GUT SCALE)"
      WRITE(16,901) 1,LGUT/DSQRT(DABS(LGUT)),"LAMBDA(MGUT,DR_bar)"
      WRITE(16,901) 2,KGUT,"KAPPA(MGUT,DR_bar)"
      WRITE(16,901) 3,ALGUT,"ALAMBDA"
      WRITE(16,901) 4,AKGUT,"AKAPPA"
      WRITE(16,901) 6,XIFGUT,"XIF"
      WRITE(16,901) 7,XISGUT,"XIS"
      WRITE(16,901) 8,MUPGUT,"MUP"
      WRITE(16,901) 9,MSPGUT,"MSP"
      WRITE(16,901) 10,MSGUT,"MS^2"
      WRITE(16,901) 12,M3HGUT,"M3H"

      WRITE(16,899) "# "
      WRITE(16,899) "# FINE-TUNING parameter d(ln Mz^2)/d(ln PG^2)"
      WRITE(16,899) "# BLOCK FINETUNING"
      WRITE(16,901) 1,FTSUSY(1),"PS=MHU"
      WRITE(16,901) 2,FTSUSY(2),"PS=MHD"
      WRITE(16,901) 3,FTSUSY(3),"PS=MS"
      WRITE(16,901) 4,FTSUSY(4),"PS=ALAMBDA"
      WRITE(16,901) 5,FTSUSY(5),"PS=AKAPPA"
      WRITE(16,901) 6,FTSUSY(6),"PS=XIF"
      WRITE(16,901) 7,FTSUSY(7),"PS=XIS"
      WRITE(16,901) 8,FTSUSY(8),"PS=MUP"
      WRITE(16,901) 9,FTSUSY(9),"PS=MSP"
      WRITE(16,901) 10,FTSUSY(10),"PS=M3H"
      WRITE(16,901) 11,FTSUSY(11),"PS=LAMBDA"
      WRITE(16,901) 12,FTSUSY(12),"PS=KAPPA"
      WRITE(16,901) 13,FTSUSY(13),"PS=HTOP"
      WRITE(16,901) 14,FTSUSY(14),"PS=G"
      WRITE(16,901) NSUSY+1,FTSUSY(NSUSY+1),"MAX"
      WRITE(16,914) NSUSY+2,INT(FTSUSY(NSUSY+2)),"IMAX"

      WRITE(16,899) "# "
      WRITE(16,899) "# REDUCED CROSS SECTIONS AT LHC"
      WRITE(16,899) "BLOCK LHCCROSSSECTIONS"
      WRITE(16,901) 11,SIG(1,1),"VBF/VH -> H1 -> tautau"
      WRITE(16,901) 12,SIG(1,2),"ggF -> H1 -> tautau"
      WRITE(16,901) 13,SIG(1,3),"VBF/VH -> H1 -> bb"
      WRITE(16,901) 14,SIG(1,4),"ttH -> H1 -> bb"
      WRITE(16,901) 15,SIG(1,5),"VBF/VH -> H1 -> ZZ/WW"
      WRITE(16,901) 16,SIG(1,6),"ggF -> H1 -> ZZ/WW"
      WRITE(16,901) 17,SIG(1,7),"VBF/VH -> H1 -> gammagamma"
      WRITE(16,901) 18,SIG(1,8),"ggF -> H1 -> gammagamma"
      WRITE(16,901) 21,SIG(2,1),"VBF/VH -> H2 -> tautau"
      WRITE(16,901) 22,SIG(2,2),"ggF -> H2 -> tautau"
      WRITE(16,901) 23,SIG(2,3),"VBF/VH -> H2 -> bb"
      WRITE(16,901) 24,SIG(2,4),"ttH -> H2 -> bb"
      WRITE(16,901) 25,SIG(2,5),"VBF/VH -> H2 -> ZZ/WW"
      WRITE(16,901) 26,SIG(2,6),"ggF -> H2 -> ZZ/WW"
      WRITE(16,901) 27,SIG(2,7),"VBF/VH -> H2 -> gammagamma"
      WRITE(16,901) 28,SIG(2,8),"ggF -> H2 -> gammagamma"
      WRITE(16,901) 31,SIG(3,1),"VBF/VH -> H3 -> tautau"
      WRITE(16,901) 32,SIG(3,2),"ggF -> H3 -> tautau"
      WRITE(16,901) 33,SIG(3,3),"VBF/VH -> H3 -> bb"
      WRITE(16,901) 34,SIG(3,4),"ttH -> H3 -> bb"
      WRITE(16,901) 35,SIG(3,5),"VBF/VH -> H3 -> ZZ/WW"
      WRITE(16,901) 36,SIG(3,6),"ggF -> H3 -> ZZ/WW"
      WRITE(16,901) 37,SIG(3,7),"VBF/VH -> H3 -> gammagamma"
      WRITE(16,901) 38,SIG(3,8),"ggF -> H3 -> gammagamma"
      WRITE(16,901) 42,SIG(4,2),"ggF -> A1 -> tautau"
      WRITE(16,901) 44,SIG(4,4),"ttH -> A1 -> bb"
      WRITE(16,901) 46,SIG(4,6),"ggF -> A1 -> ZZ/WW"
      WRITE(16,901) 48,SIG(4,8),"ggF -> A1 -> gammagamma"
      WRITE(16,901) 52,SIG(5,2),"ggF -> A2 -> tautau"
      WRITE(16,901) 54,SIG(5,4),"ttH -> A2 -> bb"
      WRITE(16,901) 56,SIG(5,6),"ggF -> A2 -> ZZ/WW"
      WRITE(16,901) 58,SIG(5,8),"ggF -> A2 -> gammagamma"

      WRITE(16,899) "# "
      WRITE(16,899)
     . "# PARAMETERS OF THE EFFECTIVE LAGRANGIAN IN THE HIGGS SECTOR"
      WRITE(16,899)
     . "# AS USED IN MICROMEGAS"
      WRITE(16,899) "BLOCK  EFFECTIVE_COUPLINGS"
      WRITE(16,917) " X",PX
      DO I=1,6
      WRITE(16,916) " A",I,PA(I)
      ENDDO
      DO I=1,2
       WRITE(16,916) " B",I,PB(I)
      ENDDO
      DO I=1,7
       WRITE(16,916) " L",I,PL(I)
      ENDDO
      DO I=1,8
       WRITE(16,916) " K",I,PK(I)
      ENDDO
      WRITE(16,917) " XVEV", MUQ/LQ*DSQRT(ZS)
      WRITE(16,917) " DELMB",DELMB
c      WRITE(16,917) " DELML",DELML
c      WRITE(16,917) " X(1)",(CB(1)/CL(1))**2
c      WRITE(16,917) " X(2)",(CB(2)/CL(2))**2
      WRITE(16,899) "# "
*  From IF(OUTFLAG.NE.1)
      ENDIF

      WRITE(16,899) "# "
      WRITE(16,899) "# HIGGS + TOP BRANCHING RATIOS IN SLHA FORMAT"
      WRITE(16,899) "# Allowed BSM branching ratio for the most SM-like"
      WRITE(16,899) "# Higgs state in the 125.1+/-3 GeV window (LHC):"
      WRITE(16,899) "# BR_BSM < 0.194 "
      WRITE(16,899) "# Info about decay package"
      WRITE(16,899) "BLOCK DCINFO   # Program information"
      WRITE(16,900) 1,"NMSSMTools # Decay package"
      WRITE(16,900) 2,"6.1.0      # Version number"

      WRITE(16,899) "#           PDG          Width"
      WRITE(16,904) PDGH1,WIDTH(1),"Lightest neutral Higgs scalar"
      IF(BRJJ(1).GT.0d0)
     .  WRITE(16,905) BRJJ(1),2,21,21,"BR(H_1 -> hadrons)"
      IF(BREE(1).GT.0d0)
     .  WRITE(16,905) BREE(1),2,11,-11,"BR(H_1 -> e- e+)"
      IF(BRMM(1).GT.0d0)
     .  WRITE(16,905) BRMM(1),2,13,-13,"BR(H_1 -> muon muon)"
      IF(BRLL(1).GT.0d0)
     .  WRITE(16,905) BRLL(1),2,15,-15,"BR(H_1 -> tau tau)"
      IF(BRCC(1).GT.0d0)
     .  WRITE(16,905) BRCC(1),2,4,-4,"BR(H_1 -> c cbar)"
      IF(BRBB(1).GT.0d0)
     .  WRITE(16,905) BRBB(1),2,5,-5,"BR(H_1 -> b bbar)"
      IF(BRTT(1).GT.0d0)
     .  WRITE(16,905) BRTT(1),2,6,-6,"BR(H_1 -> t tbar)"
      IF(BRWW(1).GT.0d0)
     .  WRITE(16,905) BRWW(1),2,24,-24,"BR(H_1 -> W+ W-)"
      IF(BRZZ(1).GT.0d0)
     .  WRITE(16,905) BRZZ(1),2,23,23,"BR(H_1 -> Z Z)"
      IF(BRGG(1).GT.0d0)
     .  WRITE(16,905) BRGG(1),2,22,22,"BR(H_1 -> gamma gamma)"
      IF(BRZG(1).GT.0d0)
     .  WRITE(16,905) BRZG(1),2,23,22,"BR(H_1 -> Z gamma)"
      IF(BRHAA(1,1).GT.0d0)
     .  WRITE(16,905) BRHAA(1,1),2,36,36,"BR(H_1 -> A_1 A_1)"
      IF(BRHAA(1,2).GT.0d0)
     .  WRITE(16,905) BRHAA(1,2),2,36,46,"BR(H_1 -> A_1 A_2)"
      IF(BRHAA(1,3).GT.0d0)
     .  WRITE(16,905) BRHAA(1,3),2,46,46,"BR(H_1 -> A_2 A_2)"
      IF(BRHAZ(1,1).GT.0d0)
     .  WRITE(16,905) BRHAZ(1,1),2,23,36,"BR(H_1 -> A_1 Z)"
      IF(BRHAZ(1,2).GT.0d0)
     .  WRITE(16,905) BRHAZ(1,2),2,23,46,"BR(H_1 -> A_2 Z)"
      IF(BRHCHC(1).GT.0d0)
     .  WRITE(16,905) BRHCHC(1),2,37,-37,"BR(H_1 -> H+ H-)"
      IF(BRHCW(1).GT.0d0)
     .  WRITE(16,905) BRHCW(1),2,24,-37,"BR(H_1 -> W+ H-)"
      IF(BRHCW(1).GT.0d0)
     .  WRITE(16,905) BRHCW(1),2,-24,37,"BR(H_1 -> W- H+)"
      IF(BRNEU(1,1,1).GT.0d0)
     .  WRITE(16,905) BRNEU(1,1,1),2,1000022,1000022,
     .    "BR(H_1 -> neu_1 neu_1)"
      IF(BRNEU(1,1,2).GT.0d0)
     .  WRITE(16,905) BRNEU(1,1,2),2,1000022,1000023,
     .    "BR(H_1 -> neu_1 neu_2)"
      IF(BRNEU(1,1,3).GT.0d0)
     .  WRITE(16,905) BRNEU(1,1,3),2,1000022,1000025,
     .    "BR(H_1 -> neu_1 neu_3)"
      IF(BRNEU(1,1,4).GT.0d0)
     .  WRITE(16,905) BRNEU(1,1,4),2,1000022,1000035,
     .    "BR(H_1 -> neu_1 neu_4)"
      IF(BRNEU(1,1,5).GT.0d0)
     .  WRITE(16,905) BRNEU(1,1,5),2,1000022,1000045,
     .    "BR(H_1 -> neu_1 neu_5)"
      IF(BRNEU(1,2,2).GT.0d0)
     .  WRITE(16,905) BRNEU(1,2,2),2,1000023,1000023,
     .    "BR(H_1 -> neu_2 neu_2)"
      IF(BRNEU(1,2,3).GT.0d0)
     .  WRITE(16,905) BRNEU(1,2,3),2,1000023,1000025,
     .    "BR(H_1 -> neu_2 neu_3)"
      IF(BRNEU(1,2,4).GT.0d0)
     .  WRITE(16,905) BRNEU(1,2,4),2,1000023,1000035,
     .    "BR(H_1 -> neu_2 neu_4)"
      IF(BRNEU(1,2,5).GT.0d0)
     .  WRITE(16,905) BRNEU(1,2,5),2,1000023,1000045,
     .    "BR(H_1 -> neu_2 neu_5)"
      IF(BRNEU(1,3,3).GT.0d0)
     .  WRITE(16,905) BRNEU(1,3,3),2,1000025,1000025,
     .    "BR(H_1 -> neu_3 neu_3)"
      IF(BRNEU(1,3,4).GT.0d0)
     .  WRITE(16,905) BRNEU(1,3,4),2,1000025,1000035,
     .    "BR(H_1 -> neu_3 neu_4)"
      IF(BRNEU(1,3,5).GT.0d0)
     .  WRITE(16,905) BRNEU(1,3,5),2,1000025,1000045,
     .    "BR(H_1 -> neu_3 neu_5)"
      IF(BRNEU(1,4,4).GT.0d0)
     .  WRITE(16,905) BRNEU(1,4,4),2,1000035,1000035,
     .    "BR(H_1 -> neu_4 neu_4)"
      IF(BRNEU(1,4,5).GT.0d0)
     .  WRITE(16,905) BRNEU(1,4,5),2,1000035,1000045,
     .    "BR(H_1 -> neu_4 neu_5)"
      IF(BRNEU(1,5,5).GT.0d0)
     .  WRITE(16,905) BRNEU(1,5,5),2,1000045,1000045,
     .    "BR(H_1 -> neu_5 neu_5)"
      IF(BRCHA(1,1).GT.0d0)
     .  WRITE(16,905) BRCHA(1,1),2,1000024,-1000024,
     .    "BR(H_1 -> cha_1 cha_1bar)"
      IF(BRCHA(1,2).GT.0d0)
     .  WRITE(16,905) BRCHA(1,2),2,1000024,-1000037,
     .    "BR(H_1 -> cha_1 cha_2bar)"
      IF(BRCHA(1,2).GT.0d0)
     .  WRITE(16,905) BRCHA(1,2),2,1000037,-1000024,
     .    "BR(H_1 -> cha_2 cha_1bar)"
      IF(BRCHA(1,3).GT.0d0)
     .  WRITE(16,905) BRCHA(1,3),2,1000037,-1000037,
     .    "BR(H_1 -> cha_2 cha_2bar)"
      IF(BRHSQ(1,1).GT.0d0)
     .  WRITE(16,905) BRHSQ(1,1),2,1000002,-1000002,
     .    "BR(H_1 -> ~u_L ~ubar_L)"
      IF(BRHSQ(1,1).GT.0d0)
     .  WRITE(16,905) BRHSQ(1,1),2,1000004,-1000004,
     .    "BR(H_1 -> ~c_L ~cbar_L)"
      IF(BRHSQ(1,2).GT.0d0)
     .  WRITE(16,905) BRHSQ(1,2),2,2000002,-2000002,
     .    "BR(H_1 -> ~u_R ~ubar_R)"
      IF(BRHSQ(1,2).GT.0d0)
     .  WRITE(16,905) BRHSQ(1,2),2,2000004,-2000004,
     .    "BR(H_1 -> ~c_R ~cbar_R)"
      IF(BRHSQ(1,3).GT.0d0)
     .  WRITE(16,905) BRHSQ(1,3),2,1000001,-1000001,
     .    "BR(H_1 -> ~d_L ~dbar_L)"
      IF(BRHSQ(1,3).GT.0d0)
     .  WRITE(16,905) BRHSQ(1,3),2,1000003,-1000003,
     .    "BR(H_1 -> ~s_L ~sbar_L)"
      IF(BRHSQ(1,4).GT.0d0)
     .  WRITE(16,905) BRHSQ(1,4),2,2000001,-2000001,
     .    "BR(H_1 -> ~d_R ~dbar_R)"
      IF(BRHSQ(1,4).GT.0d0)
     .  WRITE(16,905) BRHSQ(1,4),2,2000003,-2000003,
     .    "BR(H_1 -> ~s_R ~sbar_R)"
      IF(BRHSQ(1,5).GT.0d0)
     .  WRITE(16,905) BRHSQ(1,5),2,1000006,-1000006,
     .    "BR(H_1 -> ~t_1 ~tbar_1)"
      IF(BRHSQ(1,6).GT.0d0)
     .  WRITE(16,905) BRHSQ(1,6),2,2000006,-2000006,
     .    "BR(H_1 -> ~t_2 ~tbar_2)"
      IF(BRHSQ(1,7).GT.0d0)
     .  WRITE(16,905) BRHSQ(1,7),2,1000006,-2000006,
     .    "BR(H_1 -> ~t_1 ~tbar_2)"
      IF(BRHSQ(1,7).GT.0d0)
     .  WRITE(16,905) BRHSQ(1,7),2,2000006,-1000006,
     .    "BR(H_1 -> ~t_2 ~tbar_1)"
      IF(BRHSQ(1,8).GT.0d0)
     .  WRITE(16,905) BRHSQ(1,8),2,1000005,-1000005,
     .    "BR(H_1 -> ~b_1 ~bbar_1)"
      IF(BRHSQ(1,9).GT.0d0)
     .  WRITE(16,905) BRHSQ(1,9),2,2000005,-2000005,
     .    "BR(H_1 -> ~b_2 ~bbar_2)"
      IF(BRHSQ(1,10).GT.0d0)
     .  WRITE(16,905) BRHSQ(1,10),2,1000005,-2000005,
     .    "BR(H_1 -> ~b_1 ~bbar_2)"
      IF(BRHSQ(1,10).GT.0d0)
     .  WRITE(16,905) BRHSQ(1,10),2,2000005,-1000005,
     .    "BR(H_1 -> ~b_2 ~bbar_1)"
      IF(BRHSL(1,1).GT.0d0)
     .  WRITE(16,905) BRHSL(1,1),2,1000011,-1000011,
     .    "BR(H_1 -> ~e_L ~ebar_L)"
      IF(BRHSL(1,1).GT.0d0)
     .  WRITE(16,905) BRHSL(1,1),2,1000013,-1000013,
     .    "BR(H_1 -> ~mu_L ~mubar_L)"
      IF(BRHSL(1,2).GT.0d0)
     .  WRITE(16,905) BRHSL(1,2),2,2000011,-2000011,
     .    "BR(H_1 -> ~e_R ~ebar_R)"
      IF(BRHSL(1,2).GT.0d0)
     .  WRITE(16,905) BRHSL(1,2),2,2000013,-2000013,
     .    "BR(H_1 -> ~mu_R ~mubarRL)"
      IF(BRHSL(1,3).GT.0d0)
     .  WRITE(16,905) BRHSL(1,3),2,1000012,-1000012,
     .    "BR(H_1 -> ~nu_e_L ~nu_ebar_L)"
      IF(BRHSL(1,3).GT.0d0)
     .  WRITE(16,905) BRHSL(1,3),2,1000014,-1000014,
     .    "BR(H_1 -> ~nu_mu_L ~nu_mubar_L)"
      IF(BRHSL(1,4).GT.0d0)
     .  WRITE(16,905) BRHSL(1,4),2,1000015,-1000015,
     .    "BR(H_1 -> ~tau_1 ~taubar_1)"
      IF(BRHSL(1,5).GT.0d0)
     .  WRITE(16,905) BRHSL(1,5),2,2000015,-2000015,
     .    "BR(H_1 -> ~tau_2 ~taubar_2)"
      IF(BRHSL(1,6).GT.0d0)
     .  WRITE(16,905) BRHSL(1,6),2,1000015,-2000015,
     .    "BR(H_1 -> ~tau_1 ~taubar_2)"
      IF(BRHSL(1,6).GT.0d0)
     .  WRITE(16,905) BRHSL(1,6),2,2000015,-1000015,
     .    "BR(H_1 -> ~tau_2 ~taubar_1)"
      IF(BRHSL(1,7).GT.0d0)
     .  WRITE(16,905) BRHSL(1,7),2,1000016,-1000016,
     .    "BR(H_1 -> ~nu_tau_L ~nu_taubar_L)"

      WRITE(16,904) PDGH2,WIDTH(2),"2nd neutral Higgs scalar"
      IF(BRJJ(2).GT.0d0)
     .  WRITE(16,905) BRJJ(2),2,21,21,"BR(H_2 -> hadrons)"
      IF(BREE(2).GT.0d0)
     .  WRITE(16,905) BREE(2),2,11,-11,"BR(H_2 -> e- e+)"
      IF(BRMM(2).GT.0d0)
     .  WRITE(16,905) BRMM(2),2,13,-13,"BR(H_2 -> muon muon)"
      IF(BRLL(2).GT.0d0)
     .  WRITE(16,905) BRLL(2),2,15,-15,"BR(H_2 -> tau tau)"
      IF(BRCC(2).GT.0d0)
     .  WRITE(16,905) BRCC(2),2,4,-4,"BR(H_2 -> c cbar)"
      IF(BRBB(2).GT.0d0)
     .  WRITE(16,905) BRBB(2),2,5,-5,"BR(H_2 -> b bbar)"
      IF(BRTT(2).GT.0d0)
     .  WRITE(16,905) BRTT(2),2,6,-6,"BR(H_2 -> t tbar)"
      IF(BRWW(2).GT.0d0)
     .  WRITE(16,905) BRWW(2),2,24,-24,"BR(H_2 -> W+ W-)"
      IF(BRZZ(2).GT.0d0)
     .  WRITE(16,905) BRZZ(2),2,23,23,"BR(H_2 -> Z Z)"
      IF(BRGG(2).GT.0d0)
     .  WRITE(16,905) BRGG(2),2,22,22,"BR(H_2 -> gamma gamma)"
      IF(BRZG(2).GT.0d0)
     .  WRITE(16,905) BRZG(2),2,23,22,"BR(H_2 -> Z gamma)"
      IF(BRHHH(1).GT.0d0)
     .  WRITE(16,905) BRHHH(1),2,PDGH1,PDGH1,"BR(H_2 -> H_1 H_1)"
      IF(BRHAA(2,1).GT.0d0)
     .  WRITE(16,905) BRHAA(2,1),2,36,36,"BR(H_2 -> A_1 A_1)"
      IF(BRHAA(2,2).GT.0d0)
     .  WRITE(16,905) BRHAA(2,2),2,36,46,"BR(H_2 -> A_1 A_2)"
      IF(BRHAA(2,3).GT.0d0)
     .  WRITE(16,905) BRHAA(2,3),2,46,46,"BR(H_2 -> A_2 A_2)"
      IF(BRHAZ(2,1).GT.0d0)
     .  WRITE(16,905) BRHAZ(2,1),2,23,36,"BR(H_2 -> A_1 Z)"
      IF(BRHAZ(2,2).GT.0d0)
     .  WRITE(16,905) BRHAZ(2,2),2,23,46,"BR(H_2 -> A_2 Z)"
      IF(BRHCHC(2).GT.0d0)
     .  WRITE(16,905) BRHCHC(2),2,37,-37,"BR(H_2 -> H+ H-)"
      IF(BRHCW(2).GT.0d0)
     .  WRITE(16,905) BRHCW(2),2,24,-37,"BR(H_2 -> W+ H-)"
      IF(BRHCW(2).GT.0d0)
     .  WRITE(16,905) BRHCW(2),2,-24,37,"BR(H_2 -> W- H+)"
      IF(BRNEU(2,1,1).GT.0d0)
     .  WRITE(16,905) BRNEU(2,1,1),2,1000022,1000022,
     .    "BR(H_2 -> neu_1 neu_1)"
      IF(BRNEU(2,1,2).GT.0d0)
     .  WRITE(16,905) BRNEU(2,1,2),2,1000022,1000023,
     .    "BR(H_2 -> neu_1 neu_2)"
      IF(BRNEU(2,1,3).GT.0d0)
     .  WRITE(16,905) BRNEU(2,1,3),2,1000022,1000025,
     .    "BR(H_2 -> neu_1 neu_3)"
      IF(BRNEU(2,1,4).GT.0d0)
     .  WRITE(16,905) BRNEU(2,1,4),2,1000022,1000035,
     .    "BR(H_2 -> neu_1 neu_4)"
      IF(BRNEU(2,1,5).GT.0d0)
     .  WRITE(16,905) BRNEU(2,1,5),2,1000022,1000045,
     .    "BR(H_2 -> neu_1 neu_5)"
      IF(BRNEU(2,2,2).GT.0d0)
     .  WRITE(16,905) BRNEU(2,2,2),2,1000023,1000023,
     .    "BR(H_2 -> neu_2 neu_2)"
      IF(BRNEU(2,2,3).GT.0d0)
     .  WRITE(16,905) BRNEU(2,2,3),2,1000023,1000025,
     .    "BR(H_2 -> neu_2 neu_3)"
      IF(BRNEU(2,2,4).GT.0d0)
     .  WRITE(16,905) BRNEU(2,2,4),2,1000023,1000035,
     .    "BR(H_2 -> neu_2 neu_4)"
      IF(BRNEU(2,2,5).GT.0d0)
     .  WRITE(16,905) BRNEU(2,2,5),2,1000023,1000045,
     .    "BR(H_2 -> neu_2 neu_5)"
      IF(BRNEU(2,3,3).GT.0d0)
     .  WRITE(16,905) BRNEU(2,3,3),2,1000025,1000025,
     .    "BR(H_2 -> neu_3 neu_3)"
      IF(BRNEU(2,3,4).GT.0d0)
     .  WRITE(16,905) BRNEU(2,3,4),2,1000025,1000035,
     .    "BR(H_2 -> neu_3 neu_4)"
      IF(BRNEU(2,3,5).GT.0d0)
     .  WRITE(16,905) BRNEU(2,3,5),2,1000025,1000045,
     .    "BR(H_2 -> neu_3 neu_5)"
      IF(BRNEU(2,4,4).GT.0d0)
     .  WRITE(16,905) BRNEU(2,4,4),2,1000035,1000035,
     .    "BR(H_2 -> neu_4 neu_4)"
      IF(BRNEU(2,4,5).GT.0d0)
     .  WRITE(16,905) BRNEU(2,4,5),2,1000035,1000045,
     .    "BR(H_2 -> neu_4 neu_5)"
      IF(BRNEU(2,5,5).GT.0d0)
     .  WRITE(16,905) BRNEU(2,5,5),2,1000045,1000045,
     .    "BR(H_2 -> neu_5 neu_5)"
      IF(BRCHA(2,1).GT.0d0)
     .  WRITE(16,905) BRCHA(2,1),2,1000024,-1000024,
     .    "BR(H_2 -> cha_1 cha_1bar)"
      IF(BRCHA(2,2).GT.0d0)
     .  WRITE(16,905) BRCHA(2,2),2,1000024,-1000037,
     .    "BR(H_2 -> cha_1 cha_2bar)"
      IF(BRCHA(2,2).GT.0d0)
     .  WRITE(16,905) BRCHA(2,2),2,1000037,-1000024,
     .    "BR(H_2 -> cha_2 cha_1bar)"
      IF(BRCHA(2,3).GT.0d0)
     .  WRITE(16,905) BRCHA(2,3),2,1000037,-1000037,
     .    "BR(H_2 -> cha_2 cha_2bar)"
      IF(BRHSQ(2,1).GT.0d0)
     .  WRITE(16,905) BRHSQ(2,1),2,1000002,-1000002,
     .    "BR(H_2 -> ~u_L ~ubar_L)"
      IF(BRHSQ(2,1).GT.0d0)
     .  WRITE(16,905) BRHSQ(2,1),2,1000004,-1000004,
     .    "BR(H_2 -> ~c_L ~cbar_L)"
      IF(BRHSQ(2,2).GT.0d0)
     .  WRITE(16,905) BRHSQ(2,2),2,2000002,-2000002,
     .    "BR(H_2 -> ~u_R ~ubar_R)"
      IF(BRHSQ(2,2).GT.0d0)
     .  WRITE(16,905) BRHSQ(2,2),2,2000004,-2000004,
     .    "BR(H_2 -> ~c_R ~cbar_R)"
      IF(BRHSQ(2,3).GT.0d0)
     .  WRITE(16,905) BRHSQ(2,3),2,1000001,-1000001,
     .    "BR(H_2 -> ~d_L ~dbar_L)"
      IF(BRHSQ(2,3).GT.0d0)
     .  WRITE(16,905) BRHSQ(2,3),2,1000003,-1000003,
     .    "BR(H_2 -> ~s_L ~sbar_L)"
      IF(BRHSQ(2,4).GT.0d0)
     .  WRITE(16,905) BRHSQ(2,4),2,2000001,-2000001,
     .    "BR(H_2 -> ~d_R ~dbar_R)"
      IF(BRHSQ(2,4).GT.0d0)
     .  WRITE(16,905) BRHSQ(2,4),2,2000003,-2000003,
     .    "BR(H_2 -> ~s_R ~sbar_R)"
      IF(BRHSQ(2,5).GT.0d0)
     .  WRITE(16,905) BRHSQ(2,5),2,1000006,-1000006,
     .    "BR(H_2 -> ~t_1 ~tbar_1)"
      IF(BRHSQ(2,6).GT.0d0)
     .  WRITE(16,905) BRHSQ(2,6),2,2000006,-2000006,
     .    "BR(H_2 -> ~t_2 ~tbar_2)"
      IF(BRHSQ(2,7).GT.0d0)
     .  WRITE(16,905) BRHSQ(2,7),2,1000006,-2000006,
     .    "BR(H_2 -> ~t_1 ~tbar_2)"
      IF(BRHSQ(2,7).GT.0d0)
     .  WRITE(16,905) BRHSQ(2,7),2,2000006,-1000006,
     .    "BR(H_2 -> ~t_2 ~tbar_1)"
      IF(BRHSQ(2,8).GT.0d0)
     .  WRITE(16,905) BRHSQ(2,8),2,1000005,-1000005,
     .    "BR(H_2 -> ~b_1 ~bbar_1)"
      IF(BRHSQ(2,9).GT.0d0)
     .  WRITE(16,905) BRHSQ(2,9),2,2000005,-2000005,
     .    "BR(H_2 -> ~b_2 ~bbar_2)"
      IF(BRHSQ(2,10).GT.0d0)
     .  WRITE(16,905) BRHSQ(2,10),2,1000005,-2000005,
     .    "BR(H_2 -> ~b_1 ~bbar_2)"
      IF(BRHSQ(2,10).GT.0d0)
     .  WRITE(16,905) BRHSQ(2,10),2,2000005,-1000005,
     .    "BR(H_2 -> ~b_2 ~bbar_1)"
      IF(BRHSL(2,1).GT.0d0)
     .  WRITE(16,905) BRHSL(2,1),2,1000011,-1000011,
     .    "BR(H_2 -> ~e_L ~ebar_L)"
      IF(BRHSL(2,1).GT.0d0)
     .  WRITE(16,905) BRHSL(2,1),2,1000013,-1000013,
     .    "BR(H_2 -> ~mu_L ~mubar_L)"
      IF(BRHSL(2,2).GT.0d0)
     .  WRITE(16,905) BRHSL(2,2),2,2000011,-2000011,
     .    "BR(H_2 -> ~e_R ~ebar_R)"
      IF(BRHSL(2,2).GT.0d0)
     .  WRITE(16,905) BRHSL(2,2),2,2000013,-2000013,
     .    "BR(H_2 -> ~mu_R ~mubarRL)"
      IF(BRHSL(2,3).GT.0d0)
     .  WRITE(16,905) BRHSL(2,3),2,1000012,-1000012,
     .    "BR(H_2 -> ~nu_e_L ~nu_ebar_L)"
      IF(BRHSL(2,3).GT.0d0)
     .  WRITE(16,905) BRHSL(2,3),2,1000014,-1000014,
     .    "BR(H_2 -> ~nu_mu_L ~nu_mubar_L)"
      IF(BRHSL(2,4).GT.0d0)
     .  WRITE(16,905) BRHSL(2,4),2,1000015,-1000015,
     .    "BR(H_2 -> ~tau_1 ~taubar_1)"
      IF(BRHSL(2,5).GT.0d0)
     .  WRITE(16,905) BRHSL(2,5),2,2000015,-2000015,
     .    "BR(H_2 -> ~tau_2 ~taubar_2)"
      IF(BRHSL(2,6).GT.0d0)
     .  WRITE(16,905) BRHSL(2,6),2,1000015,-2000015,
     .    "BR(H_2 -> ~tau_1 ~taubar_2)"
      IF(BRHSL(2,6).GT.0d0)
     .  WRITE(16,905) BRHSL(2,6),2,2000015,-1000015,
     .    "BR(H_2 -> ~tau_2 ~taubar_1)"
      IF(BRHSL(2,7).GT.0d0)
     .  WRITE(16,905) BRHSL(2,7),2,1000016,-1000016,
     .    "BR(H_2 -> ~nu_tau_L ~nu_taubar_L)"

      WRITE(16,904) 45,WIDTH(3),"3rd neutral Higgs scalar"
      IF(BRJJ(3).GT.0d0)
     .  WRITE(16,905) BRJJ(3),2,21,21,"BR(H_3 -> hadrons)"
      IF(BREE(3).GT.0d0)
     .  WRITE(16,905) BREE(3),2,11,-11,"BR(H_3 -> e- e+)"
      IF(BRMM(3).GT.0d0)
     .  WRITE(16,905) BRMM(3),2,13,-13,"BR(H_3 -> muon muon)"
      IF(BRLL(3).GT.0d0)
     .  WRITE(16,905) BRLL(3),2,15,-15,"BR(H_3 -> tau tau)"
      IF(BRCC(3).GT.0d0)
     .  WRITE(16,905) BRCC(3),2,4,-4,"BR(H_3 -> c cbar)"
      IF(BRBB(3).GT.0d0)
     .  WRITE(16,905) BRBB(3),2,5,-5,"BR(H_3 -> b bbar)"
      IF(BRTT(3).GT.0d0)
     .  WRITE(16,905) BRTT(3),2,6,-6,"BR(H_3 -> t tbar)"
      IF(BRWW(3).GT.0d0)
     .  WRITE(16,905) BRWW(3),2,24,-24,"BR(H_3 -> W+ W-)"
      IF(BRZZ(3).GT.0d0)
     .  WRITE(16,905) BRZZ(3),2,23,23,"BR(H_3 -> Z Z)"
      IF(BRGG(3).GT.0d0)
     .  WRITE(16,905) BRGG(3),2,22,22,"BR(H_3 -> gamma gamma)"
      IF(BRZG(3).GT.0d0)
     .  WRITE(16,905) BRZG(3),2,23,22,"BR(H_3 -> Z gamma)"
      IF(BRHHH(2).GT.0d0)
     .  WRITE(16,905) BRHHH(2),2,PDGH1,PDGH1,"BR(H_3 -> H_1 H_1)"
      IF(BRHHH(3).GT.0d0)
     .  WRITE(16,905) BRHHH(3),2,PDGH1,PDGH2,"BR(H_3 -> H_1 H_2)"
      IF(BRHHH(4).GT.0d0)
     .  WRITE(16,905) BRHHH(4),2,PDGH2,PDGH2,"BR(H_3 -> H_2 H_2)"
      IF(BRHAA(3,1).GT.0d0)
     .  WRITE(16,905) BRHAA(3,1),2,36,36,"BR(H_3 -> A_1 A_1)"
      IF(BRHAA(3,2).GT.0d0)
     .  WRITE(16,905) BRHAA(3,2),2,36,46,"BR(H_3 -> A_1 A_2)"
      IF(BRHAA(3,3).GT.0d0)
     .  WRITE(16,905) BRHAA(3,3),2,46,46,"BR(H_3 -> A_2 A_2)"
      IF(BRHAZ(3,1).GT.0d0)
     .  WRITE(16,905) BRHAZ(3,1),2,23,36,"BR(H_3 -> A_1 Z)"
      IF(BRHAZ(3,2).GT.0d0)
     .  WRITE(16,905) BRHAZ(3,2),2,23,46,"BR(H_3 -> A_2 Z)"
      IF(BRHCHC(3).GT.0d0)
     .  WRITE(16,905) BRHCHC(3),2,37,-37,"BR(H_3 -> H+ H-)"
      IF(BRHCW(3).GT.0d0)
     .  WRITE(16,905) BRHCW(3),2,24,-37,"BR(H_3 -> W+ H-)"
      IF(BRHCW(3).GT.0d0)
     .  WRITE(16,905) BRHCW(3),2,-24,37,"BR(H_3 -> W- H+)"
      IF(BRNEU(3,1,1).GT.0d0)
     .  WRITE(16,905) BRNEU(3,1,1),2,1000022,1000022,
     .    "BR(H_3 -> neu_1 neu_1)"
      IF(BRNEU(3,1,2).GT.0d0)
     .  WRITE(16,905) BRNEU(3,1,2),2,1000022,1000023,
     .    "BR(H_3 -> neu_1 neu_2)"
      IF(BRNEU(3,1,3).GT.0d0)
     .  WRITE(16,905) BRNEU(3,1,3),2,1000022,1000025,
     .    "BR(H_3 -> neu_1 neu_3)"
      IF(BRNEU(3,1,4).GT.0d0)
     .  WRITE(16,905) BRNEU(3,1,4),2,1000022,1000035,
     .    "BR(H_3 -> neu_1 neu_4)"
      IF(BRNEU(3,1,5).GT.0d0)
     .  WRITE(16,905) BRNEU(3,1,5),2,1000022,1000045,
     .    "BR(H_3 -> neu_1 neu_5)"
      IF(BRNEU(3,2,2).GT.0d0)
     .  WRITE(16,905) BRNEU(3,2,2),2,1000023,1000023,
     .    "BR(H_3 -> neu_2 neu_2)"
      IF(BRNEU(3,2,3).GT.0d0)
     .  WRITE(16,905) BRNEU(3,2,3),2,1000023,1000025,
     .    "BR(H_3 -> neu_2 neu_3)"
      IF(BRNEU(3,2,4).GT.0d0)
     .  WRITE(16,905) BRNEU(3,2,4),2,1000023,1000035,
     .    "BR(H_3 -> neu_2 neu_4)"
      IF(BRNEU(3,2,5).GT.0d0)
     .  WRITE(16,905) BRNEU(3,2,5),2,1000023,1000045,
     .    "BR(H_3 -> neu_2 neu_5)"
      IF(BRNEU(3,3,3).GT.0d0)
     .  WRITE(16,905) BRNEU(3,3,3),2,1000025,1000025,
     .    "BR(H_3 -> neu_3 neu_3)"
      IF(BRNEU(3,3,4).GT.0d0)
     .  WRITE(16,905) BRNEU(3,3,4),2,1000025,1000035,
     .    "BR(H_3 -> neu_3 neu_4)"
      IF(BRNEU(3,3,5).GT.0d0)
     .  WRITE(16,905) BRNEU(3,3,5),2,1000025,1000045,
     .    "BR(H_3 -> neu_3 neu_5)"
      IF(BRNEU(3,4,4).GT.0d0)
     .  WRITE(16,905) BRNEU(3,4,4),2,1000035,1000035,
     .    "BR(H_3 -> neu_4 neu_4)"
      IF(BRNEU(3,4,5).GT.0d0)
     .  WRITE(16,905) BRNEU(3,4,5),2,1000035,1000045,
     .    "BR(H_3 -> neu_4 neu_5)"
      IF(BRNEU(3,5,5).GT.0d0)
     .  WRITE(16,905) BRNEU(3,5,5),2,1000045,1000045,
     .    "BR(H_3 -> neu_5 neu_5)"
      IF(BRCHA(3,1).GT.0d0)
     .  WRITE(16,905) BRCHA(3,1),2,1000024,-1000024,
     .    "BR(H_3 -> cha_1 cha_1bar)"
      IF(BRCHA(3,2).GT.0d0)
     .  WRITE(16,905) BRCHA(3,2),2,1000024,-1000037,
     .    "BR(H_3 -> cha_1 cha_2bar)"
      IF(BRCHA(3,2).GT.0d0)
     .  WRITE(16,905) BRCHA(3,2),2,1000037,-1000024,
     .    "BR(H_3 -> cha_2 cha_1bar)"
      IF(BRCHA(3,3).GT.0d0)
     .  WRITE(16,905) BRCHA(3,3),2,1000037,-1000037,
     .    "BR(H_3 -> cha_2 cha_2bar)"
      IF(BRHSQ(3,1).GT.0d0)
     .  WRITE(16,905) BRHSQ(3,1),2,1000002,-1000002,
     .    "BR(H_3 -> ~u_L ~ubar_L)"
      IF(BRHSQ(3,1).GT.0d0)
     .  WRITE(16,905) BRHSQ(3,1),2,1000004,-1000004,
     .    "BR(H_3 -> ~c_L ~cbar_L)"
      IF(BRHSQ(3,2).GT.0d0)
     .  WRITE(16,905) BRHSQ(3,2),2,2000002,-2000002,
     .    "BR(H_3 -> ~u_R ~ubar_R)"
      IF(BRHSQ(3,2).GT.0d0)
     .  WRITE(16,905) BRHSQ(3,2),2,2000004,-2000004,
     .    "BR(H_3 -> ~c_R ~cbar_R)"
      IF(BRHSQ(3,3).GT.0d0)
     .  WRITE(16,905) BRHSQ(3,3),2,1000001,-1000001,
     .    "BR(H_3 -> ~d_L ~dbar_L)"
      IF(BRHSQ(3,3).GT.0d0)
     .  WRITE(16,905) BRHSQ(3,3),2,1000003,-1000003,
     .    "BR(H_3 -> ~s_L ~sbar_L)"
      IF(BRHSQ(3,4).GT.0d0)
     .  WRITE(16,905) BRHSQ(3,4),2,2000001,-2000001,
     .    "BR(H_3 -> ~d_R ~dbar_R)"
      IF(BRHSQ(3,4).GT.0d0)
     .  WRITE(16,905) BRHSQ(3,4),2,2000003,-2000003,
     .    "BR(H_3 -> ~s_R ~sbar_R)"
      IF(BRHSQ(3,5).GT.0d0)
     .  WRITE(16,905) BRHSQ(3,5),2,1000006,-1000006,
     .    "BR(H_3 -> ~t_1 ~tbar_1)"
      IF(BRHSQ(3,6).GT.0d0)
     .  WRITE(16,905) BRHSQ(3,6),2,2000006,-2000006,
     .    "BR(H_3 -> ~t_2 ~tbar_2)"
      IF(BRHSQ(3,7).GT.0d0)
     .  WRITE(16,905) BRHSQ(3,7),2,1000006,-2000006,
     .    "BR(H_3 -> ~t_1 ~tbar_2)"
      IF(BRHSQ(3,7).GT.0d0)
     .  WRITE(16,905) BRHSQ(3,7),2,2000006,-1000006,
     .    "BR(H_3 -> ~t_2 ~tbar_1)"
      IF(BRHSQ(3,8).GT.0d0)
     .  WRITE(16,905) BRHSQ(3,8),2,1000005,-1000005,
     .    "BR(H_3 -> ~b_1 ~bbar_1)"
      IF(BRHSQ(3,9).GT.0d0)
     .  WRITE(16,905) BRHSQ(3,9),2,2000005,-2000005,
     .    "BR(H_3 -> ~b_2 ~bbar_2)"
      IF(BRHSQ(3,10).GT.0d0)
     .  WRITE(16,905) BRHSQ(3,10),2,1000005,-2000005,
     .    "BR(H_3 -> ~b_1 ~bbar_2)"
      IF(BRHSQ(3,10).GT.0d0)
     .  WRITE(16,905) BRHSQ(3,10),2,2000005,-1000005,
     .    "BR(H_3 -> ~b_2 ~bbar_1)"
      IF(BRHSL(3,1).GT.0d0)
     .  WRITE(16,905) BRHSL(3,1),2,1000011,-1000011,
     .    "BR(H_3 -> ~e_L ~ebar_L)"
      IF(BRHSL(3,1).GT.0d0)
     .  WRITE(16,905) BRHSL(3,1),2,1000013,-1000013,
     .    "BR(H_3 -> ~mu_L ~mubar_L)"
      IF(BRHSL(3,2).GT.0d0)
     .  WRITE(16,905) BRHSL(3,2),2,2000011,-2000011,
     .    "BR(H_3 -> ~e_R ~ebar_R)"
      IF(BRHSL(3,2).GT.0d0)
     .  WRITE(16,905) BRHSL(3,2),2,2000013,-2000013,
     .    "BR(H_3 -> ~mu_R ~mubarRL)"
      IF(BRHSL(3,3).GT.0d0)
     .  WRITE(16,905) BRHSL(3,3),2,1000012,-1000012,
     .    "BR(H_3 -> ~nu_e_L ~nu_ebar_L)"
      IF(BRHSL(3,3).GT.0d0)
     .  WRITE(16,905) BRHSL(3,3),2,1000014,-1000014,
     .    "BR(H_3 -> ~nu_mu_L ~nu_mubar_L)"
      IF(BRHSL(3,4).GT.0d0)
     .  WRITE(16,905) BRHSL(3,4),2,1000015,-1000015,
     .    "BR(H_3 -> ~tau_1 ~taubar_1)"
      IF(BRHSL(3,5).GT.0d0)
     .  WRITE(16,905) BRHSL(3,5),2,2000015,-2000015,
     .    "BR(H_3 -> ~tau_2 ~taubar_2)"
      IF(BRHSL(3,6).GT.0d0)
     .  WRITE(16,905) BRHSL(3,6),2,1000015,-2000015,
     .    "BR(H_3 -> ~tau_1 ~taubar_2)"
      IF(BRHSL(3,6).GT.0d0)
     .  WRITE(16,905) BRHSL(3,6),2,2000015,-1000015,
     .    "BR(H_3 -> ~tau_2 ~taubar_1)"
      IF(BRHSL(3,7).GT.0d0)
     .  WRITE(16,905) BRHSL(3,7),2,1000016,-1000016,
     .    "BR(H_3 -> ~nu_tau_L ~nu_taubar_L)"

      WRITE(16,904) 36,WIDTH(4),"Lightest pseudoscalar"
      IF(BRJJ(4).GT.0d0)
     .  WRITE(16,905) BRJJ(4),2,21,21,"BR(A_1 -> hadrons)"
      IF(BREE(4).GT.0d0)
     .  WRITE(16,905) BREE(4),2,11,-11,"BR(A_1 -> e- e+)"
      IF(BRMM(4).GT.0d0)
     .  WRITE(16,905) BRMM(4),2,13,-13,"BR(A_1 -> muon muon)"
      IF(BRLL(4).GT.0d0)
     .  WRITE(16,905) BRLL(4),2,15,-15,"BR(A_1 -> tau tau)"
      IF(BRCC(4).GT.0d0)
     .  WRITE(16,905) BRCC(4),2,4,-4,"BR(A_1 -> c cbar)"
      IF(BRBB(4).GT.0d0)
     .  WRITE(16,905) BRBB(4),2,5,-5,"BR(A_1 -> b bbar)"
      IF(BRTT(4).GT.0d0)
     .  WRITE(16,905) BRTT(4),2,6,-6,"BR(A_1 -> t tbar)"
      IF(BRGG(4).GT.0d0)
     .  WRITE(16,905) BRGG(4),2,22,22,"BR(A_1 -> gamma gamma)"
      IF(BRZG(4).GT.0d0)
     .  WRITE(16,905) BRZG(4),2,23,22,"BR(A_1 -> Z gamma)"
      IF(BRAHZ(1,1).GT.0d0)
     .  WRITE(16,905) BRAHZ(1,1),2,23,PDGH1,"BR(A_1 -> Z H_1)"
      IF(BRAHZ(1,2).GT.0d0)
     .  WRITE(16,905) BRAHZ(1,2),2,23,PDGH2,"BR(A_1 -> Z H_2)"
      IF(BRAHZ(1,3).GT.0d0)
     .  WRITE(16,905) BRAHZ(1,3),2,23,45,"BR(A_1 -> Z H_3)"
      IF(BRHCW(4).GT.0d0)
     .  WRITE(16,905) BRHCW(4),2,24,-37,"BR(A_1 -> W+ H-)"
      IF(BRHCW(4).GT.0d0)
     .  WRITE(16,905) BRHCW(4),2,-24,37,"BR(A_1 -> W- H+)"
      IF(BRNEU(4,1,1).GT.0d0)
     .  WRITE(16,905) BRNEU(4,1,1),2,1000022,1000022,
     .    "BR(A_1 -> neu_1 neu_1)"
      IF(BRNEU(4,1,2).GT.0d0)
     .  WRITE(16,905) BRNEU(4,1,2),2,1000022,1000023,
     .    "BR(A_1 -> neu_1 neu_2)"
      IF(BRNEU(4,1,3).GT.0d0)
     .  WRITE(16,905) BRNEU(4,1,3),2,1000022,1000025,
     .    "BR(A_1 -> neu_1 neu_3)"
      IF(BRNEU(4,1,4).GT.0d0)
     .  WRITE(16,905) BRNEU(4,1,4),2,1000022,1000035,
     .    "BR(A_1 -> neu_1 neu_4)"
      IF(BRNEU(4,1,5).GT.0d0)
     .  WRITE(16,905) BRNEU(4,1,5),2,1000022,1000045,
     .    "BR(A_1 -> neu_1 neu_5)"
      IF(BRNEU(4,2,2).GT.0d0)
     .  WRITE(16,905) BRNEU(4,2,2),2,1000023,1000023,
     .    "BR(A_1 -> neu_2 neu_2)"
      IF(BRNEU(4,2,3).GT.0d0)
     .  WRITE(16,905) BRNEU(4,2,3),2,1000023,1000025,
     .    "BR(A_1 -> neu_2 neu_3)"
      IF(BRNEU(4,2,4).GT.0d0)
     .  WRITE(16,905) BRNEU(4,2,4),2,1000023,1000035,
     .    "BR(A_1 -> neu_2 neu_4)"
      IF(BRNEU(4,2,5).GT.0d0)
     .  WRITE(16,905) BRNEU(4,2,5),2,1000023,1000045,
     .    "BR(A_1 -> neu_2 neu_5)"
      IF(BRNEU(4,3,3).GT.0d0)
     .  WRITE(16,905) BRNEU(4,3,3),2,1000025,1000025,
     .    "BR(A_1 -> neu_3 neu_3)"
      IF(BRNEU(4,3,4).GT.0d0)
     .  WRITE(16,905) BRNEU(4,3,4),2,1000025,1000035,
     .    "BR(A_1 -> neu_3 neu_4)"
      IF(BRNEU(4,3,5).GT.0d0)
     .  WRITE(16,905) BRNEU(4,3,5),2,1000025,1000045,
     .    "BR(A_1 -> neu_3 neu_5)"
      IF(BRNEU(4,4,4).GT.0d0)
     .  WRITE(16,905) BRNEU(4,4,4),2,1000035,1000035,
     .    "BR(A_1 -> neu_4 neu_4)"
      IF(BRNEU(4,4,5).GT.0d0)
     .  WRITE(16,905) BRNEU(4,4,5),2,1000035,1000045,
     .    "BR(A_1 -> neu_4 neu_5)"
      IF(BRNEU(4,5,5).GT.0d0)
     .  WRITE(16,905) BRNEU(4,5,5),2,1000045,1000045,
     .    "BR(A_1 -> neu_5 neu_5)"
      IF(BRCHA(4,1).GT.0d0)
     .  WRITE(16,905) BRCHA(4,1),2,1000024,-1000024,
     .    "BR(A_1 -> cha_1 cha_1bar)"
      IF(BRCHA(4,2).GT.0d0)
     .  WRITE(16,905) BRCHA(4,2),2,1000024,-1000037,
     .    "BR(A_1 -> cha_1 cha_2bar)"
      IF(BRCHA(4,2).GT.0d0)
     .  WRITE(16,905) BRCHA(4,2),2,1000037,-1000024,
     .    "BR(A_1 -> cha_2 cha_1bar)"
      IF(BRCHA(4,3).GT.0d0)
     .  WRITE(16,905) BRCHA(4,3),2,1000037,-1000037,
     .    "BR(A_1 -> cha_2 cha_2bar)"
      IF(BRASQ(1,1).GT.0d0)
     .  WRITE(16,905) BRASQ(1,1),2,1000006,-2000006,
     .    "BR(A_1 -> ~t_1 ~tbar_2)"
      IF(BRASQ(1,1).GT.0d0)
     .  WRITE(16,905) BRASQ(1,1),2,2000006,-1000006,
     .    "BR(A_1 -> ~t_2 ~tbar_1)"
      IF(BRASQ(1,2).GT.0d0)
     .  WRITE(16,905) BRASQ(1,2),2,1000005,-2000005,
     .    "BR(A_1 -> ~b_1 ~bbar_2)"
      IF(BRASQ(1,2).GT.0d0)
     .  WRITE(16,905) BRASQ(1,2),2,2000005,-1000005,
     .    "BR(A_1 -> ~b_2 ~bbar_1)"
      IF(BRASL(1).GT.0d0)
     .  WRITE(16,905) BRASL(1),2,1000015,-2000015,
     .    "BR(A_1 -> ~tau_1 ~taubar_2)"
      IF(BRASL(1).GT.0d0)
     .  WRITE(16,905) BRASL(1),2,2000015,-1000015,
     .    "BR(A_1 -> ~tau_2 ~taubar_1)"

      WRITE(16,904) 46,WIDTH(5),"2nd pseudoscalar"
      IF(BRJJ(5).GT.0d0)
     .  WRITE(16,905) BRJJ(5),2,21,21,"BR(A_2 -> hadrons)"
      IF(BREE(5).GT.0d0)
     .  WRITE(16,905) BREE(5),2,11,-11,"BR(A_2 -> e- e+)"
      IF(BRMM(5).GT.0d0)
     .  WRITE(16,905) BRMM(5),2,13,-13,"BR(A_2 -> muon muon)"
      IF(BRLL(5).GT.0d0)
     .  WRITE(16,905) BRLL(5),2,15,-15,"BR(A_2 -> tau tau)"
      IF(BRCC(5).GT.0d0)
     .  WRITE(16,905) BRCC(5),2,4,-4,"BR(A_2 -> c cbar)"
      IF(BRBB(5).GT.0d0)
     .  WRITE(16,905) BRBB(5),2,5,-5,"BR(A_2 -> b bbar)"
      IF(BRTT(5).GT.0d0)
     .  WRITE(16,905) BRTT(5),2,6,-6,"BR(A_2 -> t tbar)"
      IF(BRGG(5).GT.0d0)
     .  WRITE(16,905) BRGG(5),2,22,22,"BR(A_2 -> gamma gamma)"
      IF(BRZG(5).GT.0d0)
     .  WRITE(16,905) BRZG(5),2,23,22,"BR(A_2 -> Z gamma)"
      IF(BRAHA(1).GT.0d0)
     .  WRITE(16,905) BRAHA(1),2,36,PDGH1,"BR(A_2 -> A_1 H_1)"
      IF(BRAHA(2).GT.0d0)
     .  WRITE(16,905) BRAHA(2),2,36,PDGH2,"BR(A_2 -> A_1 H_2)"
      IF(BRAHA(3).GT.0d0)
     .  WRITE(16,905) BRAHA(3),2,36,45,"BR(A_2 -> A_1 H_3)"
      IF(BRAHZ(2,1).GT.0d0)
     .  WRITE(16,905) BRAHZ(2,1),2,23,PDGH1,"BR(A_2 -> Z H_1)"
      IF(BRAHZ(2,2).GT.0d0)
     .  WRITE(16,905) BRAHZ(2,2),2,23,PDGH2,"BR(A_2 -> Z H_2)"
      IF(BRAHZ(2,3).GT.0d0)
     .  WRITE(16,905) BRAHZ(2,3),2,23,45,"BR(A_2 -> Z H_3)"
      IF(BRHCW(5).GT.0d0)
     .  WRITE(16,905) BRHCW(5),2,24,-37,"BR(A_2 -> W+ H-)"
      IF(BRHCW(5).GT.0d0)
     .  WRITE(16,905) BRHCW(5),2,-24,37,"BR(A_2 -> W- H+)"
      IF(BRNEU(5,1,1).GT.0d0)
     .  WRITE(16,905) BRNEU(5,1,1),2,1000022,1000022,
     .    "BR(A_2 -> neu_1 neu_1)"
      IF(BRNEU(5,1,2).GT.0d0)
     .  WRITE(16,905) BRNEU(5,1,2),2,1000022,1000023,
     .    "BR(A_2 -> neu_1 neu_2)"
      IF(BRNEU(5,1,3).GT.0d0)
     .  WRITE(16,905) BRNEU(5,1,3),2,1000022,1000025,
     .    "BR(A_2 -> neu_1 neu_3)"
      IF(BRNEU(5,1,4).GT.0d0)
     .  WRITE(16,905) BRNEU(5,1,4),2,1000022,1000035,
     .    "BR(A_2 -> neu_1 neu_4)"
      IF(BRNEU(5,1,5).GT.0d0)
     .  WRITE(16,905) BRNEU(5,1,5),2,1000022,1000045,
     .    "BR(A_2 -> neu_1 neu_5)"
      IF(BRNEU(5,2,2).GT.0d0)
     .  WRITE(16,905) BRNEU(5,2,2),2,1000023,1000023,
     .    "BR(A_2 -> neu_2 neu_2)"
      IF(BRNEU(5,2,3).GT.0d0)
     .  WRITE(16,905) BRNEU(5,2,3),2,1000023,1000025,
     .    "BR(A_2 -> neu_2 neu_3)"
      IF(BRNEU(5,2,4).GT.0d0)
     .  WRITE(16,905) BRNEU(5,2,4),2,1000023,1000035,
     .    "BR(A_2 -> neu_2 neu_4)"
      IF(BRNEU(5,2,5).GT.0d0)
     .  WRITE(16,905) BRNEU(5,2,5),2,1000023,1000045,
     .    "BR(A_2 -> neu_2 neu_5)"
      IF(BRNEU(5,3,3).GT.0d0)
     .  WRITE(16,905) BRNEU(5,3,3),2,1000025,1000025,
     .    "BR(A_2 -> neu_3 neu_3)"
      IF(BRNEU(5,3,4).GT.0d0)
     .  WRITE(16,905) BRNEU(5,3,4),2,1000025,1000035,
     .    "BR(A_2 -> neu_3 neu_4)"
      IF(BRNEU(5,3,5).GT.0d0)
     .  WRITE(16,905) BRNEU(5,3,5),2,1000025,1000045,
     .    "BR(A_2 -> neu_3 neu_5)"
      IF(BRNEU(5,4,4).GT.0d0)
     .  WRITE(16,905) BRNEU(5,4,4),2,1000035,1000035,
     .    "BR(A_2 -> neu_4 neu_4)"
      IF(BRNEU(5,4,5).GT.0d0)
     .  WRITE(16,905) BRNEU(5,4,5),2,1000035,1000045,
     .    "BR(A_2 -> neu_4 neu_5)"
      IF(BRNEU(5,5,5).GT.0d0)
     .  WRITE(16,905) BRNEU(5,5,5),2,1000045,1000045,
     .    "BR(A_2 -> neu_5 neu_5)"
      IF(BRCHA(5,1).GT.0d0)
     .  WRITE(16,905) BRCHA(5,1),2,1000024,-1000024,
     .    "BR(A_2 -> cha_1 cha_1bar)"
      IF(BRCHA(5,2).GT.0d0)
     .  WRITE(16,905) BRCHA(5,2),2,1000024,-1000037,
     .    "BR(A_2 -> cha_1 cha_2bar)"
      IF(BRCHA(5,2).GT.0d0)
     .  WRITE(16,905) BRCHA(5,2),2,1000037,-1000024,
     .    "BR(A_2 -> cha_2 cha_1bar)"
      IF(BRCHA(5,3).GT.0d0)
     .  WRITE(16,905) BRCHA(5,3),2,1000037,-1000037,
     .    "BR(A_2 -> cha_2 cha_2bar)"
      IF(BRASQ(2,1).GT.0d0)
     .  WRITE(16,905) BRASQ(2,1),2,1000006,-2000006,
     .    "BR(A_2 -> ~t_1 ~tbar_2)"
      IF(BRASQ(2,1).GT.0d0)
     .  WRITE(16,905) BRASQ(2,1),2,2000006,-1000006,
     .    "BR(A_2 -> ~t_2 ~tbar_1)"
      IF(BRASQ(2,2).GT.0d0)
     .  WRITE(16,905) BRASQ(2,2),2,1000005,-2000005,
     .    "BR(A_2 -> ~b_1 ~bbar_2)"
      IF(BRASQ(2,2).GT.0d0)
     .  WRITE(16,905) BRASQ(2,2),2,2000005,-1000005,
     .    "BR(A_2 -> ~b_2 ~bbar_1)"
      IF(BRASL(2).GT.0d0)
     .  WRITE(16,905) BRASL(2),2,1000015,-2000015,
     .    "BR(A_2 -> ~tau_1 ~taubar_2)"
      IF(BRASL(2).GT.0d0)
     .  WRITE(16,905) BRASL(2),2,2000015,-1000015,
     .    "BR(A_2 -> ~tau_2 ~taubar_1)"

      WRITE(16,904) 37,HCWIDTH,"Charged Higgs"
      IF(HCBRM.GT.0d0)
     .  WRITE(16,905) HCBRM,2,-13,14,"BR(H+ -> muon nu_muon)"
      IF(HCBRL.GT.0d0)
     .  WRITE(16,905) HCBRL,2,-15,16,"BR(H+ -> tau nu_tau)"
      IF(HCBRSU.GT.0d0)
     .  WRITE(16,905) HCBRSU,2,2,-3,"BR(H+ -> u sbar)"
      IF(HCBRSC.GT.0d0)
     .  WRITE(16,905) HCBRSC,2,4,-3,"BR(H+ -> c sbar)"
      IF(HCBRBU.GT.0d0)
     .  WRITE(16,905) HCBRBU,2,2,-5,"BR(H+ -> u bbar)"
      IF(HCBRBC.GT.0d0)
     .  WRITE(16,905) HCBRBC,2,4,-5,"BR(H+ -> c bbar)"
      IF(HCBRBT.GT.0d0)
     .  WRITE(16,905) HCBRBT,2,6,-5,"BR(H+ -> t bbar)"
      IF(HCBRWH(1).GT.0d0)
     .  WRITE(16,905) HCBRWH(1),2,24,PDGH1,"BR(H+ -> W+ H_1)"
      IF(HCBRWH(2).GT.0d0)
     .  WRITE(16,905) HCBRWH(2),2,24,PDGH2,"BR(H+ -> W+ H_2)"
      IF(HCBRWH(4).GT.0d0)
     .  WRITE(16,905) HCBRWH(4),2,24,36,"BR(H+ -> W+ A_1)"
      IF(HCBRNC(1,1).GT.0d0)
     .  WRITE(16,905) HCBRNC(1,1),2,1000024,1000022,
     .    "BR(H+ -> cha_1 neu_1)"
      IF(HCBRNC(2,1).GT.0d0)
     .  WRITE(16,905) HCBRNC(2,1),2,1000024,1000023,
     .    "BR(H+ -> cha_1 neu_2)"
      IF(HCBRNC(3,1).GT.0d0)
     .  WRITE(16,905) HCBRNC(3,1),2,1000024,1000025,
     .    "BR(H+ -> cha_1 neu_3)"
      IF(HCBRNC(4,1).GT.0d0)
     .  WRITE(16,905) HCBRNC(4,1),2,1000024,1000035,
     .    "BR(H+ -> cha_1 neu_4)"
      IF(HCBRNC(5,1).GT.0d0)
     .  WRITE(16,905) HCBRNC(5,1),2,1000024,1000045,
     .    "BR(H+ -> cha_1 neu_5)"
      IF(HCBRNC(1,2).GT.0d0)
     .  WRITE(16,905) HCBRNC(1,2),2,1000037,1000022,
     .    "BR(H+ -> cha_2 neu_1)"
      IF(HCBRNC(2,2).GT.0d0)
     .  WRITE(16,905) HCBRNC(2,2),2,1000037,1000023,
     .    "BR(H+ -> cha_2 neu_2)"
      IF(HCBRNC(3,2).GT.0d0)
     .  WRITE(16,905) HCBRNC(3,2),2,1000037,1000025,
     .    "BR(H+ -> cha_2 neu_3)"
      IF(HCBRNC(4,2).GT.0d0)
     .  WRITE(16,905) HCBRNC(4,2),2,1000037,1000035,
     .    "BR(H+ -> cha_2 neu_4)"
      IF(HCBRNC(5,2).GT.0d0)
     .  WRITE(16,905) HCBRNC(5,2),2,1000037,1000045,
     .    "BR(H+ -> cha_2 neu_5)"
      IF(HCBRSQ(1).GT.0d0)
     .  WRITE(16,905) HCBRSQ(1),2,1000002,-1000001,
     .    "BR(H+ -> ~u_L ~dbar_L)"
      IF(HCBRSQ(1).GT.0d0)
     .  WRITE(16,905) HCBRSQ(1),2,1000004,-1000003,
     .    "BR(H+ -> ~c_L ~sbar_L)"
      IF(HCBRSQ(2).GT.0d0)
     .  WRITE(16,905) HCBRSQ(2),2,1000006,-1000005,
     .    "BR(H+ -> ~t_1 ~bbar_1)"
      IF(HCBRSQ(3).GT.0d0)
     .  WRITE(16,905) HCBRSQ(3),2,1000006,-2000005,
     .    "BR(H+ -> ~t_1 ~bbar_2)"
      IF(HCBRSQ(4).GT.0d0)
     .  WRITE(16,905) HCBRSQ(4),2,2000006,-1000005,
     .    "BR(H+ -> ~t_2 ~bbar_1)"
      IF(HCBRSQ(5).GT.0d0)
     .  WRITE(16,905) HCBRSQ(5),2,2000006,-2000005,
     .    "BR(H+ -> ~t_2 ~bbar_2)"
      IF(HCBRSL(1).GT.0d0)
     .  WRITE(16,905) HCBRSL(1),2,1000012,-1000011,
     .    "BR(H+ -> ~nu_e_L ~ebar_L)"
      IF(HCBRSL(1).GT.0d0)
     .  WRITE(16,905) HCBRSL(1),2,1000014,-1000013,
     .    "BR(H+ -> ~nu_mu_L ~mubar_L)"
      IF(HCBRSL(2).GT.0d0)
     .  WRITE(16,905) HCBRSL(2),2,1000016,-1000015,
     .    "BR(H+ -> ~nu_tau_L ~taubar_1)"
      IF(HCBRSL(3).GT.0d0)
     .  WRITE(16,905) HCBRSL(3),2,1000016,-2000015,
     .    "BR(H+ -> ~nu_tau_L ~taubar_2)"

      WRITE(16,904) 6,toptot,'Top Quark'
      IF(brtopbw.ne.0.D0)
     .  WRITE(16,905) brtopbw,2,5,24,'BR(t ->  b    W+)'
      IF(brtopbh.ne.0.D0)
     .  WRITE(16,905) brtopbh,2,5,37,'BR(t ->  b    H+)'
      IF(brtopneutrstop(1,1).ne.0.D0)
     .  WRITE(16,905) brtopneutrstop(1,1),2,1000006,1000022,
     . 'BR(t -> ~t_1 ~chi_10)'
      IF(brtopneutrstop(2,1).ne.0.D0)
     .  WRITE(16,905) brtopneutrstop(2,1),2,1000006,1000023,
     . 'BR(t -> ~t_1 ~chi_20)'
      IF(brtopneutrstop(3,1).ne.0.D0)
     .  WRITE(16,905) brtopneutrstop(3,1),2,1000006,1000025,
     . 'BR(t -> ~t_1 ~chi_30)'
      IF(brtopneutrstop(4,1).ne.0.D0)
     .  WRITE(16,905) brtopneutrstop(4,1),2,1000006,1000025,
     . 'BR(t -> ~t_1 ~chi_40)'
      IF(brtopneutrstop(5,1).ne.0.D0)
     .  WRITE(16,905) brtopneutrstop(5,1),2,1000006,1000025,
     . 'BR(t -> ~t_1 ~chi_50)'
      IF(brtopneutrstop(1,2).ne.0.D0)
     .  WRITE(16,905) brtopneutrstop(1,2),2,2000006,1000022,
     . 'BR(t -> ~t_2 ~chi_10)'
      IF(brtopneutrstop(2,2).ne.0.D0)
     .  WRITE(16,905) brtopneutrstop(2,2),2,2000006,1000023,
     . 'BR(t -> ~t_2 ~chi_20)'
      IF(brtopneutrstop(3,2).ne.0.D0)
     .  WRITE(16,905) brtopneutrstop(3,2),2,2000006,1000025,
     .'BR(t -> ~t_2 ~chi_30)'
      IF(brtopneutrstop(4,2).ne.0.D0)
     .  WRITE(16,905) brtopneutrstop(4,2),2,2000006,1000025,
     . 'BR(t -> ~t_2 ~chi_40)'
      IF(brtopneutrstop(5,2).ne.0.D0)
     .  WRITE(16,905) brtopneutrstop(5,2),2,2000006,1000025,
     . 'BR(t -> ~t_2 ~chi_50)'

      IF(NMSFLAG.NE.0)CALL NS_OUTPUT(16)

      IF(OMGFLAG.EQ.0) RETURN
      WRITE(17,899) "# RELIC DENSITY CALCULATED BY MICROMEGAS"
      WRITE(17,899) "#"
      WRITE(17,899) "BLOCK RDINFO   # Program information"
      WRITE(17,900) 1,"MicrOmegas # Dark matter package"
      WRITE(17,900) 2,"6.1.0      # Version number"
      IF(PROB(30).GT.0d0)
     . WRITE(17,900) 3,"# DM relic density too large"
      IF(PROB(30).LT.0d0.AND.PROB(30).GT.-1d0)
     . WRITE(17,900) 3,"# DM relic density too small"
      IF(PROB(30).LE.-1d0)
     . WRITE(17,900) 3,"# Problem in micrOMEGAs"
      IF(PROB(31).NE.0d0)THEN
       WRITE(17,900) 3,"# DM direct detection rate too large (SI)"
       IF(OMGFLAG.LT.0 .AND. OMG.LT.OMGMAX)
     .  WRITE(17,915) "         # Cross sections rescaled with the",
     .  " ratio relic density/Planck observed value"
      ENDIF
      IF(PROB(61).NE.0d0)THEN
       WRITE(17,900) 3,"# DM direct detection rate too large (SD-n)"
       IF(OMGFLAG.LT.0 .AND. OMG.LT.OMGMAX)
     .  WRITE(17,915) "         # Cross sections rescaled with the",
     .  " ratio relic density/Planck observed value"
      ENDIF
      IF(PROB(62).NE.0d0)THEN
       WRITE(17,900) 3,"# DM direct detection rate too large (SD-p)"
       IF(OMGFLAG.LT.0 .AND. OMG.LT.OMGMAX)
     .  WRITE(17,915) "         # Cross sections rescaled with the",
     .  " ratio relic density/Planck observed value"
      ENDIF
      IF(IFAIL.EQ.0.OR.IFAIL.GE.10)THEN
        CALL printRelDen(17)
      ELSE
        WRITE(17,900) 4,"# Cannot compute Omega h^2 (0<IFAIL<10)"
      ENDIF

      CLOSE(16)
      CLOSE(17)

 111  FORMAT(3A,E10.4)
 899  FORMAT(A)
 900  FORMAT(1X,I5,3X,A)
 901  FORMAT(1X,I5,3X,1P,E16.8,0P,3X,'#',1X,A)
 902  FORMAT(1X,I9,3X,1P,E16.8,0P,3X,'#',1X,A)
 903  FORMAT(1X,I2,1X,I2,3X,1P,E16.8,0P,3X,'#',1X,A)
 904  FORMAT('DECAY',1X,I9,3X,1P,E16.8,0P,3X,'#',1X,A)
 905  FORMAT(3X,1P,E16.8,0P,3X,I2,3X,I9,1X,I9,1X,2X,'#',1X,A)
 906  FORMAT('#',1X,A,3X,E16.8)
 907  FORMAT(A,1P,E16.8,A)
 911  FORMAT(A,F8.5,A,F8.5,A)
 914  FORMAT(1X,I5,3X,1P,I16,0P,3X,'#',1X,A)
 915  FORMAT(A,A)
 916  FORMAT(A,I1,4E16.8)
 917  FORMAT(A,4E16.8)
 918  FORMAT(1X,I5,3X,A,F6.1,'-',F5.1,A)
 920  FORMAT('#',0P,I5,3X,1P,E16.8,0P,3X,'#',1X,A)
 921  FORMAT(1X,I2,1X,I2,3X,'#',1X,A)
 927  FORMAT(13X,F10.2,3X,"#",X,A)
 928  FORMAT(A,1P,E10.4,A,1P,E14.8)
 929  FORMAT('#',1X,I4,3X,A,F6.1,'-',F5.1,A)
 930  FORMAT('#',1X,I4,3X,A)

      END
