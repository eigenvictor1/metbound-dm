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
*      The input file contains lower and upper bounds as well as number
*      of steps for the parameters on which the scan is performed
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
*      PMASS(1-2): CP-odd masses (ordered)
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
*          9         mu = 0 or (kappa=0 and Akappa=/=0)
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

      INTEGER NFL,NPROB,NPAR
      PARAMETER (NFL=14,NPROB=89,NPAR=25)
      INTEGER NFAIL(NFL),IFAIL,I,IMAX,TOT,ITOT,NTOT
      INTEGER M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,MSFLAG,AKFLAG,ALFLAG
      INTEGER NMSFLAG,PFLAG,OMGFLAG,MAFLAG,MOFLAG,GMUFLAG,HFLAG,GRFLAG
      INTEGER NL,NK,NTB,NMU,NMUP,NMSP,NM3H
      INTEGER N1,N2,N3,N4,NM1,NM2,NM3
      INTEGER NAU3,NAD3,NAE3,NAE2,NML3,NML2,NME3
      INTEGER NME2,NMQ3,NMQ2,NMU3,NMU2,NMD3,NMD2
      INTEGER IL,IK,ITB,IMU,IMUP,IMSP,IM3H
      INTEGER I1,I2,I3,I4,IM1,IM2,IM3
      INTEGER IAU3,IAD3,IAE3,IAE2,IML3,IML2,IME3
      INTEGER IME2,IMQ3,IMQ2,IMU3,IMU2,IMD3,IMD2
      INTEGER UNCERTFLAG,DIFAIL,MWFLAG
      INTEGER IND,strlen,CFLAG(6)

      DOUBLE PRECISION PAR(NPAR),PROB(NPROB),SM
      DOUBLE PRECISION LMIN,LMAX,KMIN,KMAX,TBMIN,TBMAX,MUMIN,MUMAX
      DOUBLE PRECISION ALMIN,ALMAX,AKMIN,AKMAX,XIFMIN,XIFMAX
      DOUBLE PRECISION XISMIN,XISMAX,MUPMIN,MUPMAX,MSPMIN,MSPMAX
      DOUBLE PRECISION M3HMIN,M3HMAX,MAMIN,MAMAX,MPMIN,MPMAX
      DOUBLE PRECISION M1MIN,M1MAX,M2MIN,M2MAX,M3MIN,M3MAX
      DOUBLE PRECISION AU3MIN,AU3MAX,AD3MIN,AD3MAX,AE3MIN,AE3MAX
      DOUBLE PRECISION AE2MIN,AE2MAX,ML3MIN,ML3MAX,ML2MIN,ML2MAX
      DOUBLE PRECISION ME3MIN,ME3MAX,ME2MIN,ME2MAX,MQ3MIN,MQ3MAX
      DOUBLE PRECISION MQ2MIN,MQ2MAX,MU3MIN,MU3MAX,MU2MIN,MU2MAX
      DOUBLE PRECISION MD3MIN,MD3MAX,MD2MIN,MD2MAX
      DOUBLE PRECISION LN,LNN,KN,KNN,TBN,TBNN,MUN,MUNN
      DOUBLE PRECISION ALN,ALNN,AKN,AKNN,XIFN,XIFNN
      DOUBLE PRECISION XISN,XISNN,MUPN,MUPNN,MSPN,MSPNN
      DOUBLE PRECISION M3HN,M3HNN,MAN,MANN,MPN,MPNN
      DOUBLE PRECISION M1N,M1NN,M2N,M2NN,M3N,M3NN
      DOUBLE PRECISION AU3N,AU3NN,AD3N,AD3NN,AE3N,AE3NN,AE2N,AE2NN
      DOUBLE PRECISION ML3N,ML3NN,ML2N,ML2NN,ME3N,ME3NN,ME2N,ME2NN
      DOUBLE PRECISION MQ3N,MQ3NN,MQ2N,MQ2NN,MU3N,MU3NN,MU2N,MU2NN
      DOUBLE PRECISION MD3N,MD3NN,MD2N,MD2NN
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
      DOUBLE PRECISION M32,CGR,MPL,Q2,DELMB,DELML,DEL1,D0,EPS

      COMMON/NMSFLAG/NMSFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/SCANFLAGS/M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,
     . MSFLAG,AKFLAG,ALFLAG
      COMMON/STEPS/NTOT,NL,NK,NTB,NMU,NMUP,NMSP,NM3H,
     . N1,N2,N3,N4,NM1,NM2,NM3,NAU3,NAD3,NAE3,NAE2,NML3,
     . NML2,NME3,NME2,NMQ3,NMQ2,NMU3,NMU2,NMD3,NMD2
      COMMON/BOUNDS/LN,LNN,KN,KNN,TBN,TBNN,MUN,MUNN,
     . ALN,ALNN,AKN,AKNN,XIFN,XIFNN,
     . XISN,XISNN,MUPN,MUPNN,MSPN,MSPNN,
     . M3HN,M3HNN,MAN,MANN,MPN,MPNN,
     . M1N,M1NN,M2N,M2NN,M3N,M3NN,
     . AU3N,AU3NN,AD3N,AD3NN,AE3N,AE3NN,AE2N,AE2NN,
     . ML3N,ML3NN,ML2N,ML2NN,ME3N,ME3NN,ME2N,ME2NN,
     . MQ3N,MQ3NN,MQ2N,MQ2NN,MU3N,MU3NN,MU2N,MU2NN,
     . MD3N,MD3NN,MD2N,MD2NN
      COMMON/MINMAX/LMIN,LMAX,KMIN,KMAX,TBMIN,TBMAX,MUMIN,MUMAX,
     . ALMIN,ALMAX,AKMIN,AKMAX,XIFMIN,XIFMAX,XISMIN,XISMAX,
     . MUPMIN,MUPMAX,MSPMIN,MSPMAX,M3HMIN,M3HMAX,MAMIN,MAMAX,
     . MPMIN,MPMAX,M1MIN,M1MAX,M2MIN,M2MAX,M3MIN,M3MAX,
     . AU3MIN,AU3MAX,AD3MIN,AD3MAX,AE3MIN,AE3MAX,
     . AE2MIN,AE2MAX,ML3MIN,ML3MAX,ML2MIN,ML2MAX,
     . ME3MIN,ME3MAX,ME2MIN,ME2MAX,MQ3MIN,MQ3MAX,
     . MQ2MIN,MQ2MAX,MU3MIN,MU3MAX,MU2MIN,MU2MAX,
     . MD3MIN,MD3MAX,MD2MIN,MD2MAX
      COMMON/SUSYEXT/XIF,XIS,MUP,MSP,M3H
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
      OFILE=trim(PRE)//'out'//trim(SUF)
      EFILE=trim(PRE)//'err'//trim(SUF)
      TFILE=trim(PRE)//'slha'//trim(SUF)
      IND=0
      DO I=1,strlen(IFILE)
       IF(IFILE(I:I).EQ.'/')IND=I
      ENDDO
      IF(IND.EQ.0)THEN
       PAT='./'
      ELSE
       PAT=trim(IFILE(1:IND))
      ENDIF
      OPEN(16,FILE=OFILE,STATUS='UNKNOWN')
      OPEN(17,FILE=EFILE,STATUS='UNKNOWN')

*   Initialization

      CALL INITIALIZE()
      DO I=1,NFL
       NFAIL(I)=0
      ENDDO
      TOT=0
      ITOT=0

*   Reading of the input parameters

      CALL INPUT(PAR,NPAR)

*   Initialization of the range of parameters that has passed all tests

      TBN=1d99
      TBNN=-1d99
      M1N=1d99
      M1NN=-1d99
      M2N=1d99
      M2NN=-1d99
      M3N=1d99
      M3NN=-1d99
      LN=1d99
      LNN=-1d99
      KN=1d99
      KNN=-1d99
      MUN=1d99
      MUNN=-1d99
      MAN=1d99
      MANN=-1d99
      ALN=1d99
      ALNN=-1d99
      XIFN=1d99
      XIFNN=-1d99
      MPN=1d99
      MPNN=-1d99
      AKN=1d99
      AKNN=-1d99
      XISN=1d99
      XISNN=-1d99
      MUPN=1d99
      MUPNN=-1d99
      MSPN=1d99
      MSPNN=-1d99
      M3HN=1d99
      M3HNN=-1d99
      AU3N=1d99
      AU3NN=-1d99
      AD3N=1d99
      AD3NN=-1d99
      AE3N=1d99
      AE3NN=-1d99
      AE2N=1d99
      AE2NN=-1d99
      ML3N=1d99
      ML3NN=-1d99
      ML2N=1d99
      ML2NN=-1d99
      ME3N=1d99
      ME3NN=-1d99
      ME2N=1d99
      ME2NN=-1d99
      MQ3N=1d99
      MQ3NN=-1d99
      MQ2N=1d99
      MQ2NN=-1d99
      MU3N=1d99
      MU3NN=-1d99
      MU2N=1d99
      MU2NN=-1d99
      MD3N=1d99
      MD3NN=-1d99
      MD2N=1d99
      MD2NN=-1d99

*   Beginning of the scan

      DO ITB=1,NTB
      IF(NTB.EQ.1)THEN
       PAR(3)=TBMIN
      ELSE
       PAR(3)=TBMIN+(TBMAX-TBMIN)*DFLOAT(ITB-1)/DFLOAT(NTB-1)
      ENDIF

      DO IM2=1,NM2
      IF(NM2.EQ.1)THEN
       PAR(21)=M2MIN
      ELSE
       PAR(21)=M2MIN+(M2MAX-M2MIN)*DFLOAT(IM2-1)/DFLOAT(NM2-1)
      ENDIF

      DO IM1=1,NM1
      IF(M1FLAG.EQ.0)THEN
       PAR(20)=PAR(21)/2d0
      ELSEIF(NM1.EQ.1)THEN
       PAR(20)=M1MIN
      ELSE
       PAR(20)=M1MIN+(M1MAX-M1MIN)*DFLOAT(IM1-1)/DFLOAT(NM1-1)
      ENDIF

      DO IM3=1,NM3
      IF(M3FLAG.EQ.0)THEN
       PAR(22)=PAR(21)*3d0
      ELSEIF(NM3.EQ.1)THEN
       PAR(22)=M3MIN
      ELSE
       PAR(22)=M3MIN+(M3MAX-M3MIN)*DFLOAT(IM3-1)/DFLOAT(NM3-1)
      ENDIF

      DO IL=1,NL
      IF(NL.EQ.1)THEN
       PAR(1)=LMIN
      ELSE
       PAR(1)=LMIN+(LMAX-LMIN)*DFLOAT(IL-1)/DFLOAT(NL-1)
      ENDIF

      DO IK=1,NK
      IF(NK.EQ.1)THEN
       PAR(2)=KMIN
      ELSE
       PAR(2)=KMIN+(KMAX-KMIN)*DFLOAT(IK-1)/DFLOAT(NK-1)
      ENDIF

      DO IMU=1,NMU
      IF(NMU.EQ.1)THEN
       PAR(4)=MUMIN
      ELSE
       PAR(4)=MUMIN+(MUMAX-MUMIN)*DFLOAT(IMU-1)/DFLOAT(NMU-1)
      ENDIF

      DO I1=1,N1
      IF(N1.EQ.1)THEN
       IF(MOD(MAFLAG,3).NE.1)THEN
        PAR(5)=ALMIN
       ELSE
        PAR(23)=MAMIN
       ENDIF
      ELSE
       IF(MOD(MAFLAG,3).NE.1)THEN
        PAR(5)=ALMIN+(ALMAX-ALMIN)*DFLOAT(I1-1)/DFLOAT(N1-1)
       ELSE
        PAR(23)=MAMIN+(MAMAX-MAMIN)*DFLOAT(I1-1)/DFLOAT(N1-1)
       ENDIF
      ENDIF

      DO I2=1,N2
      IF(N2.EQ.1)THEN
       IF(MOD(MAFLAG,3).NE.2)THEN
        XIF=XIFMIN
       ELSE
        PAR(23)=MAMIN
       ENDIF
      ELSE
       IF(MOD(MAFLAG,3).NE.2)THEN
        XIF=XIFMIN+(XIFMAX-XIFMIN)*DFLOAT(I2-1)/DFLOAT(N2-1)
       ELSE
        PAR(23)=MAMIN+(MAMAX-MAMIN)*DFLOAT(I2-1)/DFLOAT(N2-1)
       ENDIF
      ENDIF

      DO I3=1,N3
      IF(N3.EQ.1)THEN
       IF(MAFLAG/3.NE.1)THEN
        PAR(6)=AKMIN
       ELSE
        PAR(24)=MPMIN
       ENDIF
      ELSE
       IF(MAFLAG/3.NE.1)THEN
        PAR(6)=AKMIN+(AKMAX-AKMIN)*DFLOAT(I3-1)/DFLOAT(N3-1)
       ELSE
        PAR(24)=MPMIN+(MPMAX-MPMIN)*DFLOAT(I3-1)/DFLOAT(N3-1)
       ENDIF
      ENDIF

      DO I4=1,N4
      IF(N4.EQ.1)THEN
       IF(MAFLAG/3.NE.2)THEN
        XIS=XISMIN
       ELSE
        PAR(24)=MPMIN
       ENDIF
      ELSE
       IF(MAFLAG/3.NE.2)THEN
        XIS=XISMIN+(XISMAX-XISMIN)*DFLOAT(I4-1)/DFLOAT(N4-1)
       ELSE
        PAR(24)=MPMIN+(MPMAX-MPMIN)*DFLOAT(I4-1)/DFLOAT(N4-1)
       ENDIF
      ENDIF

      DO IMUP=1,NMUP
      IF(NMUP.EQ.1)THEN
       MUP=MUPMIN
      ELSE
       MUP=MUPMIN+(MUPMAX-MUPMIN)*DFLOAT(IMUP-1)/DFLOAT(NMUP-1)
      ENDIF

      DO IMSP=1,NMSP
      IF(NMSP.EQ.1)THEN
       MSP=MSPMIN
      ELSE
       MSP=MSPMIN+(MSPMAX-MSPMIN)*DFLOAT(IMSP-1)/DFLOAT(NMSP-1)
      ENDIF

      DO IM3H=1,NM3H
      IF(NM3H.EQ.1)THEN
       M3H=M3HMIN
      ELSE
       M3H=M3HMIN+(M3HMAX-M3HMIN)*DFLOAT(IM3H-1)/DFLOAT(NM3H-1)
      ENDIF

      DO IAU3=1,NAU3
      IF(NAU3.EQ.1)THEN
       PAR(12)=AU3MIN
      ELSE
       PAR(12)=AU3MIN+(AU3MAX-AU3MIN)*DFLOAT(IAU3-1)/DFLOAT(NAU3-1)
      ENDIF

      DO IAD3=1,NAD3
      IF(NAD3.EQ.1)THEN
       PAR(13)=AD3MIN
      ELSE
       PAR(13)=AD3MIN+(AD3MAX-AD3MIN)*DFLOAT(IAD3-1)/DFLOAT(NAD3-1)
      ENDIF

      DO IAE3=1,NAE3
      IF(NAE3.EQ.1)THEN
       PAR(14)=AE3MIN
      ELSE
       PAR(14)=AE3MIN+(AE3MAX-AE3MIN)*DFLOAT(IAE3-1)/DFLOAT(NAE3-1)
      ENDIF

      DO IAE2=1,NAE2
      IF(NAE2.EQ.1)THEN
       PAR(25)=AE2MIN
      ELSE
       PAR(25)=AE2MIN+(AE2MAX-AE2MIN)*DFLOAT(IAE2-1)/DFLOAT(NAE2-1)
      ENDIF

      DO IML3=1,NML3
      IF(NML3.EQ.1)THEN
       PAR(10)=ML3MIN**2
      ELSE
       PAR(10)=(ML3MIN+(ML3MAX-ML3MIN)*DFLOAT(IML3-1)/DFLOAT(NML3-1))**2
      ENDIF

      DO IML2=1,NML2
      IF(NML2.EQ.1)THEN
       PAR(18)=ML2MIN**2
      ELSE
       PAR(18)=(ML2MIN+(ML2MAX-ML2MIN)*DFLOAT(IML2-1)/DFLOAT(NML2-1))**2
      ENDIF

      DO IME3=1,NME3
      IF(NME3.EQ.1)THEN
       PAR(11)=ME3MIN**2
      ELSE
       PAR(11)=(ME3MIN+(ME3MAX-ME3MIN)*DFLOAT(IME3-1)/DFLOAT(NME3-1))**2
      ENDIF

      DO IME2=1,NME2
      IF(NME2.EQ.1)THEN
       PAR(19)=ME2MIN**2
      ELSE
       PAR(19)=(ME2MIN+(ME2MAX-ME2MIN)*DFLOAT(IME2-1)/DFLOAT(NME2-1))**2
      ENDIF

      DO IMQ3=1,NMQ3
      IF(NMQ3.EQ.1)THEN
       PAR(7)=MQ3MIN**2
      ELSE
       PAR(7)=(MQ3MIN+(MQ3MAX-MQ3MIN)*DFLOAT(IMQ3-1)/DFLOAT(NMQ3-1))**2
      ENDIF

      DO IMQ2=1,NMQ2
      IF(NMQ2.EQ.1)THEN
       PAR(15)=MQ2MIN**2
      ELSE
       PAR(15)=(MQ2MIN+(MQ2MAX-MQ2MIN)*DFLOAT(IMQ2-1)/DFLOAT(NMQ2-1))**2
      ENDIF

      DO IMU3=1,NMU3
      IF(NMU3.EQ.1)THEN
       PAR(8)=MU3MIN**2
      ELSE
       PAR(8)=(MU3MIN+(MU3MAX-MU3MIN)*DFLOAT(IMU3-1)/DFLOAT(NMU3-1))**2
      ENDIF

      DO IMU2=1,NMU2
      IF(NMU2.EQ.1)THEN
       PAR(16)=MU2MIN**2
      ELSE
       PAR(16)=(MU2MIN+(MU2MAX-MU2MIN)*DFLOAT(IMU2-1)/DFLOAT(NMU2-1))**2
      ENDIF

      DO IMD3=1,NMD3
      IF(NMD3.EQ.1)THEN
       PAR(9)=MD3MIN**2
      ELSE
       PAR(9)=(MD3MIN+(MD3MAX-MD3MIN)*DFLOAT(IMD3-1)/DFLOAT(NMD3-1))**2
      ENDIF

      DO IMD2=1,NMD2
      IF(NMD2.EQ.1)THEN
       PAR(17)=MD2MIN**2
      ELSE
       PAR(17)=(MD2MIN+(MD2MAX-MD2MIN)*DFLOAT(IMD2-1)/DFLOAT(NMD2-1))**2
      ENDIF

      ITOT=ITOT+1

!      WRITE(17,*)""
!      WRITE(17,*)"-----------------------------------------------------"
!      WRITE(17,*)""
!      WRITE(17,*)"Point ",ITOT
!      WRITE(17,*)""
!      WRITE(17,*)"MAFLAG=",MAFLAG
!      WRITE(17,*)""
!      WRITE(17,*)"TANB =",PAR(3)
!      WRITE(17,*)"M1 =",PAR(20)
!      WRITE(17,*)"M2 =",PAR(21)
!      WRITE(17,*)"M3 =",PAR(22)
!      WRITE(17,*)"LAMBDA =",PAR(1)
!      WRITE(17,*)"KAPPA =",PAR(2)
!      WRITE(17,*)"MUEFF =",PAR(4)
!      WRITE(17,*)"ALAMBDA =",PAR(5)
!      WRITE(17,*)"AKAPPA =",PAR(6)
!      WRITE(17,*)"XIF =",XIF
!      WRITE(17,*)"XIS =",XIS
!      WRITE(17,*)"MUP =",MUP
!      WRITE(17,*)"MSP =",MSP
!      WRITE(17,*)"M3H =",M3H
!      WRITE(17,*)"MA =",PAR(23)
!      WRITE(17,*)"MP =",PAR(24)
!      WRITE(17,*)"AU3 =",PAR(12)
!      WRITE(17,*)"AD3 =",PAR(13)
!      WRITE(17,*)"AE3 =",PAR(14)
!      WRITE(17,*)"AE2 =",PAR(25)
!      WRITE(17,*)"ML3 =",DSQRT(PAR(10))
!      WRITE(17,*)"ML2 =",DSQRT(PAR(18))
!      WRITE(17,*)"ME3 =",DSQRT(PAR(11))
!      WRITE(17,*)"ME2 =",DSQRT(PAR(19))
!      WRITE(17,*)"MQ3 =",DSQRT(PAR(7))
!      WRITE(17,*)"MQ2 =",DSQRT(PAR(15))
!      WRITE(17,*)"MU3 =",DSQRT(PAR(8))
!      WRITE(17,*)"MU2 =",DSQRT(PAR(16))
!      WRITE(17,*)"MD3 =",DSQRT(PAR(9))
!      WRITE(17,*)"MD2 =",DSQRT(PAR(17))
!      WRITE(17,*)""

*   Initialization of PROB and IFAIL

      DO I=1,NPROB
       PROB(I)=0d0
      ENDDO
      IFAIL=0

*   Check for mu = 0 or (kappa=0 and Akappa=/=0)

      IF(PAR(4).EQ.0d0 .OR. (PAR(2).EQ.0d0 .AND. MAFLAG/3.EQ.1))THEN
       IFAIL=9
       GOTO 11
      ENDIF

*   Check for soft terms < 3 TeV

      SM=MAX(DABS(PAR(20)),DABS(PAR(21)),DABS(PAR(22)),DABS(PAR(4)),
     .DABS(PAR(5)),DABS(PAR(6)),DABS(PAR(12)),DABS(PAR(13)),
     .DABS(PAR(14)),DABS(PAR(25)),DSQRT(PAR(10)),DSQRT(PAR(18)),
     .DSQRT(PAR(11)),DSQRT(PAR(19)),DSQRT(PAR(7)),DSQRT(PAR(15)),
     .DSQRT(PAR(8)),DSQRT(PAR(16)),DSQRT(PAR(9)),DSQRT(PAR(17)))
      PROB(89)=DDIM(SM/3d3,1d0)
      IF(PROB(89).GT.0d0)THEN
!	WRITE(17,*)"PROB",89,PROB(89)
!	WRITE(17,*)""
       IFAIL=10
       GOTO 11
      ENDIF

*   Begin loop to compute DELMB

      UNCERTFLAG=0
      DELMB=.1d0
      DELML=0d0
      DEL1=0d0
!      WRITE(17,*)""
!      WRITE(17,*)"UNCERTFLAG",UNCERTFLAG
!      WRITE(17,*)"DELMB guess",DELMB
      I=0
 1    I=I+1
      IF(I.GT.IMAX)THEN
       IFAIL=14
!       WRITE(17,*)"Exit : IFAIL =",IFAIL
!       WRITE(17,*)""
       GOTO 11
      ENDIF
      D0=DELMB

*   Computation of parameters at QSTSB

      CALL RUNPAR(PAR)

*   Computation of sfermion masses

      CALL MSFERM(PAR,IFAIL,2)
!      WRITE(17,*)"DELMB,IFAIL",DELMB,IFAIL
      IF(IFAIL.NE.0)THEN
!       WRITE(17,*)"Exit : IFAIL =",IFAIL
!       WRITE(17,*)""
       GOTO 11
      ENDIF

*   End loop to compute DELMB

      IF(2d0*DABS(DELMB-D0)/MAX(1d-3,DABS(DELMB+D0)).GT.EPS)GOTO 1

*   Computation of Higgs masses

      CALL MHIGGS(PAR,PROB,IFAIL)
      IF(IFAIL.EQ.-1)IFAIL=13
      IF(IFAIL.NE.0)THEN
!       WRITE(17,*)"Exit : IFAIL =",IFAIL
!       WRITE(17,*)""
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
!        WRITE(17,*)""
!        WRITE(17,*)"UNCERTFLAG",UNCERTFLAG
!        WRITE(17,*)"DELMB guess",DELMB
        I=0
 2      I=I+1
        IF(I.GT.IMAX)THEN
         DIFAIL=14
!         WRITE(17,*)"Exit : DIFAIL =",DIFAIL
         GOTO 3
        ENDIF
        D0=DELMB

*   Computation of parameters at QSTSB

      CALL RUNPAR(PAR)

*   Computation of sfermion masses

        CALL MSFERM(PAR,DIFAIL,2)
!        WRITE(17,*)"DELMB,DIFAIL",DELMB,DIFAIL
        IF(DIFAIL.NE.0)THEN
!         WRITE(17,*)"Exit : DIFAIL =",DIFAIL
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
!      WRITE(17,*)""
!      WRITE(17,*)"UNCERTFLAG",UNCERTFLAG
!      WRITE(17,*)"DELMB guess",DELMB
      I=0
 4    I=I+1
      IF(I.GT.IMAX)THEN
       IFAIL=14
!       WRITE(17,*)"Exit : IFAIL =",IFAIL
!       WRITE(17,*)""
       GOTO 11
      ENDIF
      D0=DELMB

*   Computation of parameters at QSTSB

      CALL RUNPAR(PAR)

*   Computation of sfermion masses

      CALL MSFERM(PAR,IFAIL,2)
!      WRITE(17,*)"DELMB,IFAIL",DELMB,IFAIL
      IF(IFAIL.NE.0)THEN
!       WRITE(17,*)"Exit : IFAIL =",IFAIL
!       WRITE(17,*)""
       GOTO 11
      ENDIF

*   End loop to compute DELMB

      IF(2d0*DABS(DELMB-D0)/MAX(1d-3,DABS(DELMB+D0)).GT.EPS)GOTO 4

*   Computation of Higgs masses

      CALL MHIGGS(PAR,PROB,IFAIL)
      IF(IFAIL.EQ.-1)IFAIL=13
      IF(IFAIL.NE.0)THEN
!       WRITE(17,*)"Exit : IFAIL =",IFAIL
!       WRITE(17,*)""
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
      IF(CFLAG(4).EQ.0)PROB(65)=0d0
      CALL TDECAY(PAR)

*   Sparticle decays

      IF(NMSFLAG.NE.0)CALL NMSDECAY(PAR)

*   Exp. constraints on sparticles (LEP, Tevatron)
*   and Higgses (LEP, Tevatron, LHC)

      IF(CFLAG(2).NE.0)CALL SUBEXP(PAR,PROB)
      IF(CFLAG(3).NE.0)CALL LHCHIG(PAR,PROB)

*   B + K physics

      IF(CFLAG(4).NE.0)THEN
       CALL BOTTOMONIUM(PROB)
       CALL BSG(PAR,PROB)
       CALL KPHYS(PAR,PROB)
       PROB(58)=0d0
      ENDIF

*   Anom. magn. moment of the Muon

      IF(GMUFLAG.NE.0)CALL MAGNMU(PAR,PROB)

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
      IF(CFLAG(6).NE.0)CALL LHCSUSY(PAR,PROB,0)

*   Check for problems

      DO I=1,NPROB
       IF(PROB(I).NE.0d0)THEN
!        WRITE(17,*)I,PROB(I)
        IFAIL=10
       ENDIF
      ENDDO
!      WRITE(17,*)

*   Computation of the fine-tuning

      CALL FTPAR(PAR,0)

*   Recording of the results

 11   CALL OUTPUT(PAR,PROB,IFAIL)
      IF(IFAIL.EQ.0)THEN
       TOT=TOT+1
       LN=MIN(PAR(1),LN)
       LNN=MAX(PAR(1),LNN)
       KN=MIN(PAR(2),KN)
       KNN=MAX(PAR(2),KNN)
       TBN=MIN(PAR(3),TBN)
       TBNN=MAX(PAR(3),TBNN)
       MUN=MIN(PAR(4),MUN)
       MUNN=MAX(PAR(4),MUNN)
       ALN=MIN(PAR(5),ALN)
       ALNN=MAX(PAR(5),ALNN)
       AKN=MIN(PAR(6),AKN)
       AKNN=MAX(PAR(6),AKNN)
       XIFN=MIN(XIF,XIFN)
       XIFNN=MAX(XIF,XIFNN)
       XISN=MIN(XIS,XISN)
       XISNN=MAX(XIS,XISNN)
       MUPN=MIN(MUP,MUPN)
       MUPNN=MAX(MUP,MUPNN)
       MSPN=MIN(MSP,MSPN)
       MSPNN=MAX(MSP,MSPNN)
       M3HN=MIN(M3H,M3HN)
       M3HNN=MAX(M3H,M3HNN)
       M1N=MIN(PAR(20),M1N)
       M1NN=MAX(PAR(20),M1NN)
       M2N=MIN(PAR(21),M2N)
       M2NN=MAX(PAR(21),M2NN)
       M3N=MIN(PAR(22),M3N)
       M3NN=MAX(PAR(22),M3NN)
       MAN=MIN(PAR(23),MAN)
       MANN=MAX(PAR(23),MANN)
       MPN=MIN(PAR(24),MPN)
       MPNN=MAX(PAR(24),MPNN)
       AU3N=MIN(PAR(12),AU3N)
       AU3NN=MAX(PAR(12),AU3NN)
       AD3N=MIN(PAR(13),AD3N)
       AD3NN=MAX(PAR(13),AD3NN)
       AE3N=MIN(PAR(14),AE3N)
       AE3NN=MAX(PAR(14),AE3NN)
       AE2N=MIN(PAR(25),AE2N)
       AE2NN=MAX(PAR(25),AE2NN)
       ML3N=MIN(DSQRT(PAR(10)),ML3N)
       ML3NN=MAX(DSQRT(PAR(10)),ML3NN)
       ML2N=MIN(DSQRT(PAR(18)),ML2N)
       ML2NN=MAX(DSQRT(PAR(18)),ML2NN)
       ME3N=MIN(DSQRT(PAR(11)),ME3N)
       ME3NN=MAX(DSQRT(PAR(11)),ME3NN)
       ME2N=MIN(DSQRT(PAR(19)),ME2N)
       ME2NN=MAX(DSQRT(PAR(19)),ME2NN)
       MQ3N=MIN(DSQRT(PAR(7)),MQ3N)
       MQ3NN=MAX(DSQRT(PAR(7)),MQ3NN)
       MQ2N=MIN(DSQRT(PAR(15)),MQ2N)
       MQ2NN=MAX(DSQRT(PAR(15)),MQ2NN)
       MU3N=MIN(DSQRT(PAR(8)),MU3N)
       MU3NN=MAX(DSQRT(PAR(8)),MU3NN)
       MU2N=MIN(DSQRT(PAR(16)),MU2N)
       MU2NN=MAX(DSQRT(PAR(16)),MU2NN)
       MD3N=MIN(DSQRT(PAR(9)),MD3N)
       MD3NN=MAX(DSQRT(PAR(9)),MD3NN)
       MD2N=MIN(DSQRT(PAR(17)),MD2N)
       MD2NN=MAX(DSQRT(PAR(17)),MD2NN)
      ELSE
       NFAIL(IFAIL)=NFAIL(IFAIL)+1
      ENDIF

      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

*   Summary of the scanning:
*   Number of points that passed/failed the tests
*   and range for scanned parameters

      CALL ERROR(TOT,NTOT,NFAIL)

      CLOSE(16)
      CLOSE(17)

      END


      SUBROUTINE INPUT(PAR,NPAR)

*******************************************************************
*   This subroutine reads SM and NMSSM parameters from input file   .
*******************************************************************

      IMPLICIT NONE

      CHARACTER CHINL*120,CHBLCK*60,CHDUM*120
      CHARACTER(200) PAT,IFILE,OFILE,EFILE,TFILE,SFILE

      INTEGER I,NLINE,INL,ICH,IX,IVAL,Q2FIX,VFLAG,Z3FLAG
      INTEGER N0,NLOOP,NBER,NPAR,ERR,GMUFLAG,HFLAG,GRFLAG
      INTEGER M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,MSFLAG
      INTEGER AKFLAG,ALFLAG,OMGFLAG,MAFLAG,MOFLAG,PFLAG,NMSFLAG
      INTEGER NTOT,NL,NK,NTB,NMU,NMUP,NMSP,NM3H
      INTEGER N1,N2,N3,N4,NM1,NM2,NM3
      INTEGER NAU3,NAD3,NAE3,NAE2,NML3,NML2,NME3
      INTEGER NME2,NMQ3,NMQ2,NMU3,NMU2,NMD3,NMD2
      INTEGER NAL,NMA,NXIF,NAK,NMP,NXIS
      INTEGER MWFLAG,CFLAG(6)

      DOUBLE PRECISION PAR(*),VAL
      DOUBLE PRECISION ACC,XITLA,XLAMBDA,MC0,MB0,MT0
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,Q2,Q2MIN
      DOUBLE PRECISION LMIN,LMAX,KMIN,KMAX,TBMIN,TBMAX,MUMIN,MUMAX
      DOUBLE PRECISION ALMIN,ALMAX,AKMIN,AKMAX,XIFMIN,XIFMAX
      DOUBLE PRECISION XISMIN,XISMAX,MUPMIN,MUPMAX,MSPMIN,MSPMAX
      DOUBLE PRECISION M3HMIN,M3HMAX,MAMIN,MAMAX,MPMIN,MPMAX
      DOUBLE PRECISION M1MIN,M1MAX,M2MIN,M2MAX,M3MIN,M3MAX
      DOUBLE PRECISION AU3MIN,AU3MAX,AD3MIN,AD3MAX,AE3MIN,AE3MAX
      DOUBLE PRECISION AE2MIN,AE2MAX,ML3MIN,ML3MAX,ML2MIN,ML2MAX
      DOUBLE PRECISION ME3MIN,ME3MAX,ME2MIN,ME2MAX,MQ3MIN,MQ3MAX
      DOUBLE PRECISION MQ2MIN,MQ2MAX,MU3MIN,MU3MAX,MU2MIN,MU2MAX
      DOUBLE PRECISION MD3MIN,MD3MAX,MD2MIN,MD2MAX
      DOUBLE PRECISION M32,CGR,MPL

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
      COMMON/RENSCALE/Q2
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/MINMAX/LMIN,LMAX,KMIN,KMAX,TBMIN,TBMAX,MUMIN,MUMAX,
     . ALMIN,ALMAX,AKMIN,AKMAX,XIFMIN,XIFMAX,XISMIN,XISMAX,
     . MUPMIN,MUPMAX,MSPMIN,MSPMAX,M3HMIN,M3HMAX,MAMIN,MAMAX,
     . MPMIN,MPMAX,M1MIN,M1MAX,M2MIN,M2MAX,M3MIN,M3MAX,
     . AU3MIN,AU3MAX,AD3MIN,AD3MAX,AE3MIN,AE3MAX,
     . AE2MIN,AE2MAX,ML3MIN,ML3MAX,ML2MIN,ML2MAX,
     . ME3MIN,ME3MAX,ME2MIN,ME2MAX,MQ3MIN,MQ3MAX,
     . MQ2MIN,MQ2MAX,MU3MIN,MU3MAX,MU2MIN,MU2MAX,
     . MD3MIN,MD3MAX,MD2MIN,MD2MAX
      COMMON/STEPS/NTOT,NL,NK,NTB,NMU,NMUP,NMSP,NM3H,
     . N1,N2,N3,N4,NM1,NM2,NM3,NAU3,NAD3,NAE3,NAE2,NML3,
     . NML2,NME3,NME2,NMQ3,NMQ2,NMU3,NMU2,NMD3,NMD2

      COMMON/SCANFLAGS/M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,
     . MSFLAG,AKFLAG,ALFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/PFLAG/PFLAG
      COMMON/NMSFLAG/NMSFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/VFLAG/VFLAG
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

*   INITIALIZATION OF THE SCANNING PARAMETERS
      TBMIN=1d99
      TBMAX=1d99
      M1MIN=1d99
      M1MAX=1d99
      M2MIN=1d99
      M2MAX=1d99
      M3MIN=1d99
      M3MAX=1d99
      LMIN=1d99
      LMAX=1d99
      KMIN=0d0
      KMAX=1d99
      ALMIN=1d99
      ALMAX=1d99
      AKMIN=1d99
      AKMAX=1d99
      MUMIN=1d99
      MUMAX=1d99
      XIFMIN=1d99
      XIFMAX=1d99
      XISMIN=1d99
      XISMAX=1d99
      MUPMIN=0d0
      MUPMAX=1d99
      MSPMIN=0d0
      MSPMAX=1d99
      M3HMIN=0d0
      M3HMAX=1d99
      MAMIN=1d99
      MAMAX=1d99
      MPMIN=1d99
      MPMAX=1d99
      AU3MIN=1d99
      AU3MAX=1d99
      AD3MIN=1d99
      AD3MAX=1d99
      AE3MIN=1d99
      AE3MAX=1d99
      AE2MIN=1d99
      AE2MAX=1d99
      ML3MIN=1d99
      ML3MAX=1d99
      ML2MIN=1d99
      ML2MAX=1d99
      ME3MIN=1d99
      ME3MAX=1d99
      ME2MIN=1d99
      ME2MAX=1d99
      MQ3MIN=1d99
      MQ3MAX=1d99
      MQ2MIN=1d99
      MQ2MAX=1d99
      MU3MIN=1d99
      MU3MAX=1d99
      MU2MIN=1d99
      MU2MAX=1d99
      MD3MIN=1d99
      MD3MAX=1d99
      MD2MIN=1d99
      MD2MAX=1d99
      NM1=0
      NM2=1
      NM3=0
      NL=1
      NK=1
      NTB=1
      NMU=1
      NAL=0
      NAK=0
      NXIF=0
      NXIS=0
      NMA=0
      NMP=0
      NMUP=1
      NMSP=1
      NM3H=1
      NAU3=1
      NAD3=1
      NAE3=1
      NAE2=1
      NML3=1
      NML2=1
      NME3=1
      NME2=1
      NMQ3=1
      NMQ2=1
      NMU3=1
      NMU2=1
      NMD3=1
      NMD2=1
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
      M1FLAG=0
      M3FLAG=0
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
*   THE BLOCK MODSELMIN MUST CONTAIN THE LINE "  9     1    "
      IF(CHBLCK(1:6).EQ.'MODSEL')THEN
       READ(CHINL,*,ERR=999) IX,IVAL
       IF(IX.EQ.1) Z3FLAG=IVAL
       IF(IX.EQ.8) PFLAG=IVAL
       IF(IX.EQ.9) OMGFLAG=IVAL
       IF(IX.EQ.11) GMUFLAG=IVAL
       IF(IX.EQ.13) NMSFLAG=IVAL
       IF(IX.EQ.14) VFLAG=IVAL
       IF(IX.EQ.15) MOFLAG=IVAL
       IF(IX.EQ.17) CFLAG(1)=IVAL
       IF(IX.EQ.18) CFLAG(2)=IVAL
       IF(IX.EQ.19) CFLAG(3)=IVAL
       IF(IX.EQ.20) CFLAG(4)=IVAL
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
       IF(IX.EQ.37) TBMIN=VAL
       IF(IX.EQ.38) TBMAX=VAL
       IF(IX.EQ.6) CGR=VAL

*   READ EXTPAR
      ELSEIF(CHBLCK(1:6).EQ.'EXTPAR')THEN
       READ(CHINL,*,ERR=999) IX,VAL
       IF(IX.EQ.0) Q2=VAL**2
       IF(IX.EQ.107) M1MIN=VAL
       IF(IX.EQ.108) M1MAX=VAL
       IF(IX.EQ.207) M2MIN=VAL
       IF(IX.EQ.208) M2MAX=VAL
       IF(IX.EQ.307) M3MIN=VAL
       IF(IX.EQ.308) M3MAX=VAL
       IF(IX.EQ.117) AU3MIN=VAL
       IF(IX.EQ.118) AU3MAX=VAL
       IF(IX.EQ.127) AD3MIN=VAL
       IF(IX.EQ.128) AD3MAX=VAL
       IF(IX.EQ.137) AE3MIN=VAL
       IF(IX.EQ.138) AE3MAX=VAL
       IF(IX.EQ.167) AE2MIN=VAL
       IF(IX.EQ.168) AE2MAX=VAL
       IF(IX.EQ.337) ML3MIN=VAL
       IF(IX.EQ.338) ML3MAX=VAL
       IF(IX.EQ.327) ML2MIN=VAL
       IF(IX.EQ.328) ML2MAX=VAL
       IF(IX.EQ.367) ME3MIN=VAL
       IF(IX.EQ.368) ME3MAX=VAL
       IF(IX.EQ.357) ME2MIN=VAL
       IF(IX.EQ.358) ME2MAX=VAL
       IF(IX.EQ.437) MQ3MIN=VAL
       IF(IX.EQ.438) MQ3MAX=VAL
       IF(IX.EQ.427) MQ2MIN=VAL
       IF(IX.EQ.428) MQ2MAX=VAL
       IF(IX.EQ.467) MU3MIN=VAL
       IF(IX.EQ.468) MU3MAX=VAL
       IF(IX.EQ.457) MU2MIN=VAL
       IF(IX.EQ.458) MU2MAX=VAL
       IF(IX.EQ.497) MD3MIN=VAL
       IF(IX.EQ.498) MD3MAX=VAL
       IF(IX.EQ.487) MD2MIN=VAL
       IF(IX.EQ.488) MD2MAX=VAL
       IF(IX.EQ.617) LMIN=VAL
       IF(IX.EQ.618) LMAX=VAL
       IF(IX.EQ.627) KMIN=VAL
       IF(IX.EQ.628) KMAX=VAL
       IF(IX.EQ.637) ALMIN=VAL
       IF(IX.EQ.638) ALMAX=VAL
       IF(IX.EQ.647) AKMIN=VAL
       IF(IX.EQ.648) AKMAX=VAL
       IF(IX.EQ.657) MUMIN=VAL
       IF(IX.EQ.658) MUMAX=VAL
       IF(IX.EQ.667) XIFMIN=VAL
       IF(IX.EQ.668) XIFMAX=VAL
       IF(IX.EQ.677) XISMIN=VAL
       IF(IX.EQ.678) XISMAX=VAL
       IF(IX.EQ.687) MUPMIN=VAL
       IF(IX.EQ.688) MUPMAX=VAL
       IF(IX.EQ.697) MSPMIN=VAL
       IF(IX.EQ.698) MSPMAX=VAL
       IF(IX.EQ.727) M3HMIN=VAL
       IF(IX.EQ.728) M3HMAX=VAL
       IF(IX.EQ.1247) MAMIN=VAL
       IF(IX.EQ.1248) MAMAX=VAL
       IF(IX.EQ.1257) MPMIN=VAL
       IF(IX.EQ.1258) MPMAX=VAL

*   READ STEPS
       ELSEIF(CHBLCK(1:6).EQ.'STEPS')THEN
       READ(CHINL,*,ERR=999) IX,IVAL
       IF(IX.EQ.109) NM1=IVAL
       IF(IX.EQ.209) NM2=IVAL
       IF(IX.EQ.309) NM3=IVAL
       IF(IX.EQ.39) NTB=IVAL
       IF(IX.EQ.619) NL=IVAL
       IF(IX.EQ.629) NK=IVAL
       IF(IX.EQ.639) NAL=IVAL
       IF(IX.EQ.649) NAK=IVAL
       IF(IX.EQ.659) NMU=IVAL
       IF(IX.EQ.669) NXIF=IVAL
       IF(IX.EQ.679) NXIS=IVAL
       IF(IX.EQ.689) NMUP=IVAL
       IF(IX.EQ.699) NMSP=IVAL
       IF(IX.EQ.729) NM3H=IVAL
       IF(IX.EQ.1249) NMA=IVAL
       IF(IX.EQ.1259) NMP=IVAL
       IF(IX.EQ.119) NAU3=IVAL
       IF(IX.EQ.129) NAD3=IVAL
       IF(IX.EQ.139) NAE3=IVAL
       IF(IX.EQ.169) NAE2=IVAL
       IF(IX.EQ.339) NML3=IVAL
       IF(IX.EQ.329) NML2=IVAL
       IF(IX.EQ.369) NME3=IVAL
       IF(IX.EQ.359) NME2=IVAL
       IF(IX.EQ.439) NMQ3=IVAL
       IF(IX.EQ.429) NMQ2=IVAL
       IF(IX.EQ.469) NMU3=IVAL
       IF(IX.EQ.459) NMU2=IVAL
       IF(IX.EQ.499) NMD3=IVAL
       IF(IX.EQ.489) NMD2=IVAL

      ENDIF

      GOTO 21

*   END OF READING FROM INPUT FILE

*   Check for errors

 29   CLOSE(15)
      ERR=0
      IF(CFLAG(5).NE.0 .AND. NMSFLAG.EQ.0)THEN
       WRITE(0,2)"CMS CHARG(NEUTRAL)INO CONSTRAINTS CANNOT BE CHECKED ",
     .  "IF NMSDECAY IS NOT CALLED"
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
       WRITE(0,1)"LHC CONSTRAINT FLAG MUST BE IN [-1,2]"
       ERR=1
      ENDIF
      IF(GMUFLAG.LT.-1 .OR. GMUFLAG.GT.1)THEN
       WRITE(0,1)"|GMUFLAG| MUST BE IN [0,1]"
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
      IF(MWFLAG.LT.-1 .OR. MWFLAG.GT.1)THEN
       WRITE(0,1)"|MWFLAG| MUST BE IN [0,1]"
       ERR=1
      ENDIF
      IF(LMIN.EQ.1d99)THEN
       WRITE(0,1)"LMIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ELSEIF(LMIN.LE.0d0 .OR. LMAX.LE.0d0)THEN
       WRITE(0,1)"LAMBDA MUST BE STRICTLY POSITIVE"
       ERR=1
      ENDIF
      IF(TBMIN.EQ.1d99)THEN
       WRITE(0,1)"TANBMIN MUST BE GIVEN IN BLOCK MINPAR"
       ERR=1
      ELSEIF(TBMIN.LE.0d0 .OR. TBMAX.LE.0d0)THEN
       WRITE(0,1)"TANB MUST BE STRICTLY POSITIVE"
       ERR=1
      ENDIF
      IF(MUMIN.EQ.1d99)THEN
       WRITE(0,1)"MUEFFMIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ELSEIF(MUMIN.EQ.0d0 .OR. MUMAX.EQ.0d0)THEN
       WRITE(0,1)"MUEFF MUST BE NON ZERO"
       ERR=1
      ENDIF
      IF(M1MIN.EQ.1d99 .AND. M1MAX.NE.1d99)THEN
       WRITE(0,1)"M1MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(M2MIN.EQ.1d99)THEN
       WRITE(0,1)"M2MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(M3MIN.EQ.1d99 .AND. M3MAX.NE.1d99)THEN
       WRITE(0,1)"M3MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(AU3MIN.EQ.1d99)THEN
       WRITE(0,1)"AU3MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(AD3MIN.EQ.1d99)THEN
       WRITE(0,1)"AD3MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(AE3MIN.EQ.1d99)THEN
       WRITE(0,1)"AE3MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MQ3MIN.EQ.1d99)THEN
       WRITE(0,1)"MQ3MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MU3MIN.EQ.1d99)THEN
       WRITE(0,1)"MU3MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MD3MIN.EQ.1d99)THEN
       WRITE(0,1)"MD3MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(ML3MIN.EQ.1d99)THEN
       WRITE(0,1)"ML3MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(ME3MIN.EQ.1d99)THEN
       WRITE(0,1)"ME3MIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MAMIN.EQ.1d99 .AND. MAMAX.NE.1d99)THEN
       WRITE(0,1)"MAMIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(MPMIN.EQ.1d99 .AND. MPMAX.NE.1d99)THEN
       WRITE(0,1)"MPMIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(ALMIN.EQ.1d99 .AND. ALMAX.NE.1d99)THEN
       WRITE(0,1)"ALMIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(AKMIN.EQ.1d99 .AND. AKMAX.NE.1d99)THEN
       WRITE(0,1)"AKMIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XIFMIN.EQ.1d99 .AND. XIFMAX.NE.1d99)THEN
       WRITE(0,1)"XIFMIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(XISMIN.EQ.1d99 .AND. XISMAX.NE.1d99)THEN
       WRITE(0,1)"XISMIN MUST BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF
      IF(AE2MIN.EQ.1d99)THEN
       AE2MIN=AE3MIN
      ENDIF
      IF(ML2MIN.EQ.1d99)THEN
       ML2MIN=ML3MIN
      ENDIF
      IF(ME2MIN.EQ.1d99)THEN
       ME2MIN=ME3MIN
      ENDIF
      IF(MQ2MIN.EQ.1d99)THEN
       MQ2MIN=MQ3MIN
      ENDIF
      IF(MU2MIN.EQ.1d99)THEN
       MU2MIN=MU3MIN
      ENDIF
      IF(MD2MIN.EQ.1d99)THEN
       MD2MIN=MD3MIN
      ENDIF

*   Relations between (ALAMBDA, MA, XIF) and (AKAPPA, MP, XIS)

      IF(ALMIN.NE.1d99 .AND. XIFMIN.NE.1d99 .AND. MAMIN.NE.1d99)THEN
       WRITE(0,1)"AT MOST 2 OF THE 3 PARAMETERS ALAMBDA, MA AND XIF",
     . " CAN BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF

      IF(AKMIN.NE.1d99 .AND. XISMIN.NE.1d99 .AND. MPMIN.NE.1d99)THEN
       WRITE(0,1)"AT MOST 2 OF THE 3 PARAMETERS AKAPPA, MP AND XIS",
     . " CAN BE GIVEN IN BLOCK EXTPAR"
       ERR=1
      ENDIF

      IF(KMIN.EQ.0d0 .AND. (KMAX.EQ.0d0 .OR. KMAX.EQ.1d99))THEN
       IF((AKMIN.NE.0d0 .AND. AKMIN.NE.1d99) .OR. (AKMIN.EQ.0d0
     . .AND. AKMAX.NE.0d0 .AND. AKMAX.NE.1d99))THEN
        WRITE(0,1)"IF KAPPA IS 0, AKAPPA MUST BE 0"
        ERR=1
       ELSE
        AKMIN=0d0
        AKMAX=0d0
        IF(XISMIN.NE.1d99 .AND. MPMIN.NE.1d99)THEN
         WRITE(0,1)"IF KAPPA IS 0, EITHER MP OR XIS",
     .   " CAN BE GIVEN IN BLOCK EXTPAR"
         ERR=1
        ENDIF
       ENDIF
      ENDIF

*   Set default values

      IF(ALMIN.EQ.1d99.AND.MAMIN.EQ.1d99.AND.XIFMIN.EQ.1d99)THEN
       ALMIN=0d0
       XIFMIN=0d0
      ELSEIF(ALMIN.EQ.1d99.AND.MAMIN.EQ.1d99)THEN
       ALMIN=0d0
      ELSEIF(ALMIN.EQ.1d99.AND.XIFMIN.EQ.1d99)THEN
       XIFMIN=0d0
      ELSEIF(MAMIN.EQ.1d99.AND.XIFMIN.EQ.1d99)THEN
       XIFMIN=0d0
      ENDIF

      IF(AKMIN.EQ.1d99.AND.MPMIN.EQ.1d99.AND.XISMIN.EQ.1d99)THEN
       AKMIN=0d0
       XISMIN=0d0
      ELSEIF(AKMIN.EQ.1d99.AND.MPMIN.EQ.1d99)THEN
       AKMIN=0d0
      ELSEIF(AKMIN.EQ.1d99.AND.XISMIN.EQ.1d99)THEN
       XISMIN=0d0
      ELSEIF(MPMIN.EQ.1d99.AND.XISMIN.EQ.1d99)THEN
       XISMIN=0d0
      ENDIF

*   Set MAFLAG, SCANFLAGS

      IF(MAMIN.EQ.1d99)MAFLAG=0
      IF(ALMIN.EQ.1d99)MAFLAG=1
      IF(XIFMIN.EQ.1d99)MAFLAG=2
      IF(AKMIN.EQ.1d99)MAFLAG=MAFLAG+3
      IF(XISMIN.EQ.1d99)MAFLAG=MAFLAG+6
      IF(M1MIN.NE.1d99)M1FLAG=1
      IF(M3MIN.NE.1d99)M3FLAG=1

*   Bounds

      IF(TBMAX.EQ.1d99)TBMAX=TBMIN
      IF(M1MAX.EQ.1d99)M1MAX=M1MIN
      IF(M2MAX.EQ.1d99)M2MAX=M2MIN
      IF(M3MAX.EQ.1d99)M3MAX=M3MIN
      IF(LMAX.EQ.1d99)LMAX=LMIN
      IF(KMAX.EQ.1d99)KMAX=KMIN
      IF(ALMAX.EQ.1d99)ALMAX=ALMIN
      IF(AKMAX.EQ.1d99)AKMAX=AKMIN
      IF(MUMAX.EQ.1d99)MUMAX=MUMIN
      IF(XIFMAX.EQ.1d99)XIFMAX=XIFMIN
      IF(XISMAX.EQ.1d99)XISMAX=XISMIN
      IF(MUPMAX.EQ.1d99)MUPMAX=MUPMIN
      IF(MSPMAX.EQ.1d99)MSPMAX=MSPMIN
      IF(M3HMAX.EQ.1d99)M3HMAX=M3HMIN
      IF(MAMAX.EQ.1d99)MAMAX=MAMIN
      IF(MPMAX.EQ.1d99)MPMAX=MPMIN
      IF(AU3MAX.EQ.1d99)AU3MAX=AU3MIN
      IF(AD3MAX.EQ.1d99)AD3MAX=AD3MIN
      IF(AE3MAX.EQ.1d99)AE3MAX=AE3MIN
      IF(AE2MAX.EQ.1d99)AE2MAX=AE2MIN
      IF(ML3MAX.EQ.1d99)ML3MAX=ML3MIN
      IF(ML2MAX.EQ.1d99)ML2MAX=ML2MIN
      IF(ME3MAX.EQ.1d99)ME3MAX=ME3MIN
      IF(ME2MAX.EQ.1d99)ME2MAX=ME2MIN
      IF(MQ3MAX.EQ.1d99)MQ3MAX=MQ3MIN
      IF(MQ2MAX.EQ.1d99)MQ2MAX=MQ2MIN
      IF(MU3MAX.EQ.1d99)MU3MAX=MU3MIN
      IF(MU2MAX.EQ.1d99)MU2MAX=MU2MIN
      IF(MD3MAX.EQ.1d99)MD3MAX=MD3MIN
      IF(MD2MAX.EQ.1d99)MD2MAX=MD2MIN

*   Number of points

      IF(MOD(MAFLAG,3).EQ.0)THEN
       IF(NMA.NE.0)THEN
        WRITE(0,1)"NMA CANNOT BE GIVEN IN BLOCK STEPS"
       ERR=1
       ENDIF
       IF(NAL.EQ.0)NAL=1
       IF(NXIF.EQ.0)NXIF=1
       N1=NAL
       N2=NXIF
      ENDIF

      IF(MOD(MAFLAG,3).EQ.1)THEN
       IF(NAL.NE.0)THEN
        WRITE(0,1)"NAL CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NMA.EQ.0)NMA=1
       IF(NXIF.EQ.0)NXIF=1
       N1=NMA
       N2=NXIF
      ENDIF

      IF(MOD(MAFLAG,3).EQ.2)THEN
       IF(NXIF.NE.0)THEN
        WRITE(0,1)"NXIF CANNOT BE GIVEN IN BLOCK STEPS"
       ERR=1
       ENDIF
       IF(NAL.EQ.0)NAL=1
       IF(NMA.EQ.0)NMA=1
       N1=NAL
       N2=NMA
      ENDIF

      IF(MAFLAG/3.EQ.0)THEN
       IF(NMP.NE.0)THEN
        WRITE(0,1)"NMP CANNOT BE GIVEN IN BLOCK STEPS"
       ERR=1
       ENDIF
       IF(NAK.EQ.0)NAK=1
       IF(NXIS.EQ.0)NXIS=1
       N3=NAK
       N4=NXIS
      ENDIF

      IF(MAFLAG/3.EQ.1)THEN
       IF(NAK.NE.0)THEN
        WRITE(0,1)"NAK CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ENDIF
       IF(NMP.EQ.0)NMP=1
       IF(NXIS.EQ.0)NXIS=1
       N3=NMP
       N4=NXIS
      ENDIF

      IF(MAFLAG/3.EQ.2)THEN
       IF(NXIS.NE.0)THEN
        WRITE(0,1)"NXIS CANNOT BE GIVEN IN BLOCK STEPS"
       ERR=1
       ENDIF
       IF(NAK.EQ.0)NAK=1
       IF(NMP.EQ.0)NMP=1
       N3=NAK
       N4=NMP
      ENDIF

      IF(M1FLAG.EQ.0)THEN
       IF(NM1.NE.0)THEN
        WRITE(0,1)"NM1 CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ELSE
        NM1=1
       ENDIF
      ELSE
       IF(NM1.EQ.0)NM1=1
      ENDIF

      IF(M3FLAG.EQ.0)THEN
       IF(NM3.NE.0)THEN
        WRITE(0,1)"NM3 CANNOT BE GIVEN IN BLOCK STEPS"
        ERR=1
       ELSE
        NM3=1
       ENDIF
      ELSE
       IF(NM3.EQ.0)NM3=1
      ENDIF

      IF(NM1.LE.0 .OR. NM2.LE.0 .OR. NM3.LE.0 .OR. NTB.LE.0 .OR.
     . NL.LE.0 .OR. NK.LE.0 .OR. N1.LE.0 .OR. N2.LE.0 .OR.
     . N3.LE.0 .OR. N4.LE.0 .OR. NMU.LE.0 .OR. NMUP.LE.0 .OR.
     . NMSP.LE.0 .OR. NM3H.LE.0 .OR. NAU3.LE.0 .OR. NAD3.LE.0 .OR.
     . NAE3.LE.0 .OR. NAE2.LE.0 .OR. NML3.LE.0 .OR. NML2.LE.0 .OR. 
     . NME3.LE.0 .OR. NME2.LE.0 .OR. NMQ3.LE.0 .OR. NMQ2.LE.0 .OR. 
     . NMU3.LE.0 .OR. NMU2.LE.0 .OR. NMD3.LE.0 .OR. NMD2.LE.0)THEN
       WRITE(0,1)"WRONG NUMBER OF POINTS IN BLOCK STEPS"
       ERR=1
      ENDIF

*   Check for Z3 breaking terms

      IF(MOD(MAFLAG,3).EQ.2 .OR. MAFLAG/3.EQ.2 .OR.
     . MUPMIN.NE.0d0 .OR. MUPMAX.NE.0d0 .OR. MSPMIN.NE.0d0 .OR. 
     . MSPMAX.NE.0d0 .OR. XIFMIN.NE.0d0 .OR. XIFMAX.NE.0d0 .OR.
     . XISMIN.NE.0d0 .OR. XISMAX.NE.0d0 .OR. M3HMIN.NE.0d0 .OR.
     . M3HMAX.NE.0d0)THEN
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
       WRITE(0,1)"ERROR IN INPUT FILE"
       STOP 1
      ENDIF

*   Warnings

      IF(TBMIN.EQ.TBMAX .AND. NTB.GT.1)THEN
       WRITE(0,1)"WARNING TANBMIN = TANBMAX => NTB=1"
       NTB=1
      ENDIF
      IF(M1MIN.EQ.M1MAX .AND. NM1.GT.1)THEN
       WRITE(0,1)"WARNING M1MIN = M1MAX => NM1=1"
       NM1=1
      ENDIF
      IF(M2MIN.EQ.M2MAX .AND. NM2.GT.1)THEN
       WRITE(0,1)"WARNING M2MIN = M2MAX => NM2=1"
       NM2=1
      ENDIF
      IF(M3MIN.EQ.M3MAX .AND. NM3.GT.1)THEN
       WRITE(0,1)"WARNING M3MIN = M3MAX => NM3=1"
       NM3=1
      ENDIF
      IF(LMIN.EQ.LMAX .AND. NL.GT.1)THEN
       WRITE(0,1)"WARNING LMIN = LMAX => NL=1"
       NL=1
      ENDIF
      IF(KMIN.EQ.KMAX .AND. NK.GT.1)THEN
       WRITE(0,1)"WARNING KMIN = KMAX => NK=1"
       NK=1
      ENDIF
      IF(MUMIN.EQ.MUMAX .AND. NMU.GT.1)THEN
       WRITE(0,1)"WARNING MUEFFMIN = MUEFFMAX => NMU=1"
       NMU=1
      ENDIF
      IF(MUPMIN.EQ.MUPMAX .AND. NMUP.GT.1)THEN
       WRITE(0,1)"WARNING MUPMIN = MUPMAX => NMUP=1"
       NMUP=1
      ENDIF
      IF(MSPMIN.EQ.MSPMAX .AND. NMSP.GT.1)THEN
       WRITE(0,1)"WARNING MSPMIN = MSPMAX => NMSP=1"
       NMSP=1
      ENDIF
      IF(M3HMIN.EQ.M3HMAX .AND. NM3H.GT.1)THEN
       WRITE(0,1)"WARNING M3HMIN = M3HMAX => NM3H=1"
       NM3H=1
      ENDIF
      IF(AU3MIN.EQ.AU3MAX .AND. NAU3.GT.1)THEN
       WRITE(0,1)"WARNING AU3MIN = AU3MAX => NAU3=1"
       NAU3=1
      ENDIF
      IF(AD3MIN.EQ.AD3MAX .AND. NAD3.GT.1)THEN
       WRITE(0,1)"WARNING AD3MIN = AD3MAX => NAD3=1"
       NAD3=1
      ENDIF
      IF(AE3MIN.EQ.AE3MAX .AND. NAE3.GT.1)THEN
       WRITE(0,1)"WARNING AE3MIN = AE3MAX => NAE3=1"
       NAE3=1
      ENDIF
      IF(AE2MIN.EQ.AE2MAX .AND. NAE2.GT.1)THEN
       WRITE(0,1)"WARNING AE2MIN = AE2MAX => NAE2=1"
       NAE2=1
      ENDIF
      IF(ML3MIN.EQ.ML3MAX .AND. NML3.GT.1)THEN
       WRITE(0,1)"WARNING ML3MIN = ML3MAX => NML3=1"
       NML3=1
      ENDIF
      IF(ML2MIN.EQ.ML2MAX .AND. NML2.GT.1)THEN
       WRITE(0,1)"WARNING ML2MIN = ML2MAX => NML2=1"
       NML2=1
      ENDIF
      IF(ME3MIN.EQ.ME3MAX .AND. NME3.GT.1)THEN
       WRITE(0,1)"WARNING ME3MIN = ME3MAX => NME3=1"
       NME3=1
      ENDIF
      IF(ME2MIN.EQ.ME2MAX .AND. NME2.GT.1)THEN
       WRITE(0,1)"WARNING ME2MIN = ME2MAX => NME2=1"
       NME2=1
      ENDIF
      IF(MQ3MIN.EQ.MQ3MAX .AND. NMQ3.GT.1)THEN
       WRITE(0,1)"WARNING MQ3MIN = MQ3MAX => NMQ3=1"
       NMQ3=1
      ENDIF
      IF(MQ2MIN.EQ.MQ2MAX .AND. NMQ2.GT.1)THEN
       WRITE(0,1)"WARNING MQ2MIN = MQ2MAX => NMQ2=1"
       NMQ2=1
      ENDIF
      IF(MU3MIN.EQ.MU3MAX .AND. NMU3.GT.1)THEN
       WRITE(0,1)"WARNING MU3MIN = MU3MAX => NMU3=1"
       NMU3=1
      ENDIF
      IF(MU2MIN.EQ.MU2MAX .AND. NMU2.GT.1)THEN
       WRITE(0,1)"WARNING MU2MIN = MU2MAX => NMU2=1"
       NMU2=1
      ENDIF
      IF(MD3MIN.EQ.MD3MAX .AND. NMD3.GT.1)THEN
       WRITE(0,1)"WARNING MD3MIN = MD3MAX => NMD3=1"
       NMD3=1
      ENDIF
      IF(MD2MIN.EQ.MD2MAX .AND. NMD2.GT.1)THEN
       WRITE(0,1)"WARNING MD2MIN = MD2MAX => NMD2=1"
       NMD2=1
      ENDIF

      IF(TBMIN.NE.TBMAX .AND. NTB.EQ.1)THEN
       WRITE(0,10)"WARNING NTB=1 => TANBMAX=TANBMIN=",TBMIN
       TBMAX=TBMIN
      ENDIF
      IF(M1MIN.NE.M1MAX .AND. NM1.EQ.1)THEN
       WRITE(0,10)"WARNING NM1=1 => M1MAX=M1MIN=",M1MIN
       M1MAX=M1MIN
      ENDIF
      IF(M2MIN.NE.M2MAX .AND. NM2.EQ.1)THEN
       WRITE(0,10)"WARNING NM2=1 => M2MAX=M2MIN=",M2MIN
       M2MAX=M2MIN
      ENDIF
      IF(M3MIN.NE.M3MAX .AND. NM3.EQ.1)THEN
       WRITE(0,10)"WARNING NM3=1 => M3MAX=M3MIN=",M3MIN
       M3MAX=M3MIN
      ENDIF
      IF(LMIN.NE.LMAX .AND. NL.EQ.1)THEN
       WRITE(0,10)"WARNING NL=1 => LMAX=LMIN=",LMIN
       LMAX=LMIN
      ENDIF
      IF(KMIN.NE.KMAX .AND. NK.EQ.1)THEN
       WRITE(0,10)"WARNING NK=1 => KMAX=KMIN=",KMIN
       KMAX=KMIN
      ENDIF
      IF(MUMIN.NE.MUMAX .AND. NMU.EQ.1)THEN
       WRITE(0,10)"WARNING NMU=1 => MUEFFMAX=MUEFFMIN=",MUMIN
       MUMAX=MUMIN
      ENDIF
      IF(MUPMIN.NE.MUPMAX .AND. NMUP.EQ.1)THEN
       WRITE(0,10)"WARNING NMUP=1 => MUPMAX=MUPMIN=",MUPMIN
       MUPMAX=MUPMIN
      ENDIF
      IF(MSPMIN.NE.MSPMAX .AND. NMSP.EQ.1)THEN
       WRITE(0,10)"WARNING NMSP=1 =>  MSPMAX= MSPMIN=",MSPMIN
       MSPMAX=MSPMIN
      ENDIF
      IF(M3HMIN.NE.M3HMAX .AND. NM3H.EQ.1)THEN
       WRITE(0,10)"WARNING NM3H=1 => M3HMAX=M3HMIN=",M3HMIN
       M3HMAX=M3HMIN
      ENDIF
      IF(AU3MIN.NE.AU3MAX .AND. NAU3.EQ.1)THEN
       WRITE(0,10)"WARNING NAU3=1 => AU3MAX=AU3MIN=",AU3MIN
       AU3MAX=AU3MIN
      ENDIF
      IF(AD3MIN.NE.AD3MAX .AND. NAD3.EQ.1)THEN
       WRITE(0,10)"WARNING NAD3=1 => AD3MAX=AD3MIN=",AD3MIN
       AD3MAX=AD3MIN
      ENDIF
      IF(AE3MIN.NE.AE3MAX .AND. NAE3.EQ.1)THEN
       WRITE(0,10)"WARNING NAE3=1 => AE3MAX=AE3MIN=",AE3MIN
       AE3MAX=AE3MIN
      ENDIF
      IF(AE2MIN.NE.AE2MAX .AND. NAE2.EQ.1)THEN
       WRITE(0,10)"WARNING NAE2=1 => AE2MAX=AE2MIN=",AE2MIN
       AE2MAX=AE2MIN
      ENDIF
      IF(ML3MIN.NE.ML3MAX .AND. NML3.EQ.1)THEN
       WRITE(0,10)"WARNING NML3=1 => ML3MAX=ML3MIN=",ML3MIN
       ML3MAX=ML3MIN
      ENDIF
      IF(ML2MIN.NE.ML2MAX .AND. NML2.EQ.1)THEN
       WRITE(0,10)"WARNING NML2=1 => ML2MAX=ML2MIN=",ML2MIN
       ML2MAX=ML2MIN
      ENDIF
      IF(ME3MIN.NE.ME3MAX .AND. NME3.EQ.1)THEN
       WRITE(0,10)"WARNING NME3=1 => ME3MAX=ME3MIN=",ME3MIN
       ME3MAX=ME3MIN
      ENDIF
      IF(ME2MIN.NE.ME2MAX .AND. NME2.EQ.1)THEN
       WRITE(0,10)"WARNING NME2=1 => ME2MAX=ME2MIN=",ME2MIN
       ME2MAX=ME2MIN
      ENDIF
      IF(MQ3MIN.NE.MQ3MAX .AND. NMQ3.EQ.1)THEN
       WRITE(0,10)"WARNING NMQ3=1 => MQ3MAX=MQ3MIN=",MQ3MIN
       MQ3MAX=MQ3MIN
      ENDIF
      IF(MQ2MIN.NE.MQ2MAX .AND. NMQ2.EQ.1)THEN
       WRITE(0,10)"WARNING NMQ2=1 => MQ2MAX=MQ2MIN=",MQ2MIN
       MQ2MAX=MQ2MIN
      ENDIF
      IF(MU3MIN.NE.MU3MAX .AND. NMU3.EQ.1)THEN
       WRITE(0,10)"WARNING NMU3=1 => MU3MAX=MU3MIN=",MU3MIN
       MU3MAX=MU3MIN
      ENDIF
      IF(MU2MIN.NE.MU2MAX .AND. NMU2.EQ.1)THEN
       WRITE(0,10)"WARNING NMU2=1 => MU2MAX=MU2MIN=",MU2MIN
       MU2MAX=MU2MIN
      ENDIF
      IF(MD3MIN.NE.MD3MAX .AND. NMD3.EQ.1)THEN
       WRITE(0,10)"WARNING NMD3=1 => MD3MAX=MD3MIN=",MD3MIN
       MD3MAX=MD3MIN
      ENDIF
      IF(MD2MIN.NE.MD2MAX .AND. NMD2.EQ.1)THEN
       WRITE(0,10)"WARNING NMD2=1 => MD2MAX=MD2MIN=",MD2MIN
       MD2MAX=MD2MIN
      ENDIF

      IF(MOD(MAFLAG,3).EQ.0)THEN
       IF(ALMIN.EQ.ALMAX .AND. NAL.GT.1)THEN
        WRITE(0,1)"WARNING ALMIN=ALMAX => NAL=1"
        NAL=1
        N1=NAL
       ENDIF
       IF(XIFMIN.EQ.XIFMAX .AND. NXIF.GT.1)THEN
        WRITE(0,1)"WARNING XIFMIN=XIFMAX => NXIF=1"
        NXIF=1
        N2=NXIF
       ENDIF
       IF(ALMIN.NE.ALMAX .AND. NAL.EQ.1)THEN
        WRITE(0,10)"WARNING NAL=1 => ALMAX=ALMIN=",ALMIN
        ALMAX=ALMIN
       ENDIF
       IF(XIFMIN.NE.XIFMAX .AND. NXIF.EQ.1)THEN
        WRITE(0,10)"WARNING NXIF=1 => XIFMAX=XIFMIN=",XIFMIN
        XIFMAX=XIFMIN
       ENDIF
      ENDIF

      IF(MOD(MAFLAG,3).EQ.1)THEN
       IF(MAMIN.EQ.MAMAX .AND. NMA.GT.1)THEN
        WRITE(0,1)"WARNING MAMIN=MAMAX => NMA=1"
        NMA=1
        N1=NMA
       ENDIF
       IF(XIFMIN.EQ.XIFMAX .AND. NXIF.GT.1)THEN
        WRITE(0,1)"WARNING XIFMIN=XIFMAX => NXIF=1"
        NXIF=1
        N2=NXIF
       ENDIF
       IF(MAMIN.NE.MAMAX .AND. NMA.EQ.1)THEN
        WRITE(0,10)"WARNING NMA=1 => MAMAX=MAMIN=",MAMIN
        MAMAX=MAMIN
       ENDIF
       IF(XIFMIN.NE.XIFMAX .AND. NXIF.EQ.1)THEN
        WRITE(0,10)"WARNING NXIF=1 => XIFMAX=XIFMIN=",XIFMIN
        XIFMAX=XIFMIN
       ENDIF
      ENDIF

      IF(MOD(MAFLAG,3).EQ.2)THEN
       IF(ALMIN.EQ.ALMAX .AND. NAL.GT.1)THEN
        WRITE(0,1)"WARNING ALMIN=ALMAX => NAL=1"
        NAL=1
        N1=NAL
       ENDIF
       IF(MAMIN.EQ.MAMAX .AND. NMA.GT.1)THEN
        WRITE(0,1)"WARNING MAMIN=MAMAX => NMA=1"
        NMA=1
        N2=NMA
       ENDIF
       IF(ALMIN.NE.ALMAX .AND. NAL.EQ.1)THEN
        WRITE(0,10)"WARNING NAL=1 => ALMAX=ALMIN=",ALMIN
        ALMAX=ALMIN
       ENDIF
       IF(MAMIN.NE.MAMAX .AND. NMA.EQ.1)THEN
        WRITE(0,10)"WARNING NMA=1 => MAMAX=MAMIN=",MAMIN
        MAMAX=MAMIN
       ENDIF
      ENDIF

      IF(MAFLAG/3.EQ.0)THEN
       IF(AKMIN.EQ.AKMAX .AND. NAK.GT.1)THEN
        WRITE(0,1)"WARNING AKMIN=AKMAX => NAK=1"
        NAK=1
        N3=NAK
       ENDIF
       IF(XISMIN.EQ.XISMAX .AND. NXIS.GT.1)THEN
        WRITE(0,1)"WARNING XISMIN=XISMAX => NXIS=1"
        NXIS=1
        N4=NXIS
       ENDIF
       IF(AKMIN.NE.AKMAX .AND. NAK.EQ.1)THEN
        WRITE(0,10)"WARNING NAK=1 => AKMAX=AKMIN=",AKMIN
        AKMAX=AKMIN
       ENDIF
       IF(XISMIN.NE.XISMAX .AND. NXIS.EQ.1)THEN
        WRITE(0,10)"WARNING NXIS=1 => XISMAX=XISMIN=",XISMIN
        XISMAX=XISMIN
       ENDIF
      ENDIF

      IF(MAFLAG/3.EQ.1)THEN
       IF(MPMIN.EQ.MPMAX .AND. NMP.GT.1)THEN
        WRITE(0,1)"WARNING MPMIN=MPMAX => NMP=1"
        NMP=1
        N3=NMP
       ENDIF
       IF(XISMIN.EQ.XISMAX .AND. NXIS.GT.1)THEN
        WRITE(0,1)"WARNING XISMIN=XISMAX => NXIS=1"
        NXIS=1
        N4=NXIS
       ENDIF
       IF(MPMIN.NE.MPMAX .AND. NMP.EQ.1)THEN
        WRITE(0,10)"WARNING NMP=1 => MPMAX=MPMIN=",MPMIN
        MPMAX=MPMIN
       ENDIF
       IF(XISMIN.NE.XISMAX .AND. NXIS.EQ.1)THEN
        WRITE(0,10)"WARNING NXIS=1 => XISMAX=XISMIN=",XISMIN
        XISMAX=XISMIN
       ENDIF
      ENDIF

      IF(MAFLAG/3.EQ.2)THEN
       IF(AKMIN.EQ.AKMAX .AND. NAK.GT.1)THEN
        WRITE(0,1)"WARNING AKMIN=AKMAX => NAK=1"
        NAK=1
        N3=NAK
       ENDIF
       IF(MPMIN.EQ.MPMAX .AND. NMP.GT.1)THEN
        WRITE(0,1)"WARNING MPMIN=MPMAX => NMP=1"
        NMP=1
        N4=NMP
       ENDIF
       IF(AKMIN.NE.AKMAX .AND. NAK.EQ.1)THEN
        WRITE(0,10)"WARNING NAK=1 => AKMAX=AKMIN=",AKMIN
        AKMAX=AKMIN
       ENDIF
       IF(MPMIN.NE.MPMAX .AND. NMP.EQ.1)THEN
        WRITE(0,10)"WARNING NMP=1 => MPMAX=MPMIN=",MPMIN
        MPMAX=MPMIN
       ENDIF
      ENDIF

*   Total number of points

      NTOT=NTB*NM1*NM2*NM3*NL*NMU*NK*NMUP*NMSP*NM3H*N1*N2*N3*N4

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
 10   FORMAT(A,E10.3)

      END


      SUBROUTINE OUTPUT(PAR,PROB,IFAIL)

*********************************************************************
*   Subroutine writing all the results in the the output file.
*********************************************************************

      IMPLICIT NONE

      CHARACTER CHAN*20

      INTEGER IFAIL,NBIN,I,NRES,IRES,GRFLAG,NSUSY,NGUT,NMES,IMAX
      PARAMETER (NSUSY=14,NGUT=21,NMES=21,IMAX=200)

      DOUBLE PRECISION RES(IMAX),PAR(*),PROB(*),SIG(5,8),R,S,ggF13
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
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
      DOUBLE PRECISION VUS,VCB,VUB
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU,Q2
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
*
      DOUBLE PRECISION chartot2(2),chartot(2),chartot3(2)
      DOUBLE PRECISION brcharsel(2),brcharser(2),brcharsmu1(2),
     .         brcharsmu2(2),brcharstau1(2),brcharstau2(2),
     .         brcharsne1(2),brcharsne2(2),brcharsnm1(2),brcharsnm2(2),
     .         brcharsnt1(2),brcharsnt2(2),brcharsupl(2),brcharsupr(2),
     .         brcharsdownl(2),brcharsdownr(2),brcharst1(2),
     .         brcharst2(2),brcharsb1(2),brcharsb2(2),brcharwneut(2,5),
     .         brcharhcneut(2,5),brcharzchic,brcharHchic(3),
     .         brcharAchic(2),brntaunut(2,5),brnelnue(2,5),
     .         brnmunumu(2,5),brnupdb(2,5),brnchsb(2,5),brntopbb(2,5),
     .         brglupdb(2),brglchsb(2),brgltopbb(2),brchee,brchmumu,
     .         brchtautau,brchnene,brchnmunmu,brchntauntau,brchupup,
     .         brchdodo,brchchch,brchstst,brchtoptop,brchbotbot
*
      DOUBLE PRECISION neuttot2(5),neuttot(5),neuttot3(5),neuttotrad(5)
      DOUBLE PRECISION brneutst1(5),brneutst2(5),brneutsb1(5),
     .         brneutsb2(5),
     .         brneutsupl(5),brneutsupr(5),brneutsdownl(5),
     .         brneutsdownr(5),brneutsnel(5),brneutsn1(5),
     .         brneutsn2(5),brneutsell(5),brneutselr(5),
     .         brneutsnmu(5),brneutsmu1(5),brneutsmu2(5),
     .         brneutstau1(5),brneutstau2(5),brneutwchar(5,2),
     .         brneuthcchar(5,2),brneutzneut(5,5),
     .         brneutHneut(5,5,3),brneutAneut(5,5,2),brnraddec(5,5)
      DOUBLE PRECISION brneutup(5,5),brneutdow(5,5),brneutch(5,5),
     .         brneutst(5,5),brneutbot(5,5),brneuttop(5,5),
     .         brneutel(5,5),brneutmu(5,5),brneuttau(5,5),
     .         brneutnue(5,5),brneutnumu(5,5),brneutnutau(5,5),
     .         brchubd(5,2),brchcbs(5,2),brchtbb(5,2),brchelne(5,2),
     .         brchmunmu(5,2),brchtauntau(5,2),brglup(5),brgldo(5),
     .         brglch(5),brglst(5),brgltop(5),brglbot(5)
*
      DOUBLE PRECISION selltot,selltot2,selltot3,selrtot,selrtot2,
     .         selrtot3,smu1tot,smu1tot2,smu1tot3,smu2tot,smu2tot2,
     .         smu2tot3,stau1tot2,stau2tot,stau2tot2,stau2tot3,
     .         snelltot,snelltot2,snelltot3,snmu1tot,snmu1tot2,
     .         snmu1tot3,sntautot,sntautot2,sntautot3
      DOUBLE PRECISION brsellneute(5),brselrneute(5),brsellcharnue(2),
     .         brselrcharnue(2),brsnellneut(5),brsnellchar(5),
     .         brsmu1neutmu(5),brsmu2neutmu(5),brsmu1charnumu(2),
     .         brsmu2charnumu(2),brsnmu1neut(5),brsnmu1char(5),
     .         brstau1neut(5),brstau2neut(5),brstau1char(2),
     .         brstau2char(2),brstau1hcsn(2),brstau2hcsn(2),
     .         brstau1wsn(2),brstau2wsn(2),brstau2ztau,brstau2H(3),
     .         brstau2A(2),brsntauneut(5),brsntauchar(2),
     .         brsntau1hcstau(2),brsntau1wstau(2)
      DOUBLE PRECISION brsellstau1star,brsellstau1,
     .         brsellstau1nutau,brselrstau1star,brselrstau1,
     .         brselrstau1nutau,brsnestau1star,brsnestau1,
     .         brsnestau1nutau,brsmu1stau1star,brsmu1stau1,
     .         brsmu1stau1nutau,brsmu2stau1star,brsmu2stau1,
     .         brsmu2stau1nutau,brsnmustau1star,brsnmustau1,
     .         brsnmustau1nutau,brstau2stau1star,brstau2stau1,
     .         brstau2stau1nn,brsntaustau1star,brsntaustau1,
     .         brsntaustau1nutau
*
      DOUBLE PRECISION supltot2,suprtot2,sdowltot2,sdowrtot2
      DOUBLE PRECISION brsuplnup(5),brsuplcdow(2),brsuplglui,
     .         brsuprnup(5),brsuprcdow(2),brsuprglui,
     .         brsdowlndow(5),brsdowlchup(2),brsdowlglui,
     .         brsdowrndow(5),brsdowrchup(2),brsdowrglui
*
      DOUBLE PRECISION stoptot(2),stoptot2(2),stoptot3(2),stoptotrad(2)
      DOUBLE PRECISION brst1neutt(5),brst2neutt(5),brst1charb(2),
     .         brst2charb(2),brst1hcsb(2),brst2hcsb(2),brst1wsb(2),
     .         brst2wsb(2),brst1glui,brst2glui,brst2H(3),brst2A(2),
     .         brst2ztop,brgamma,brgammaup,brgammagluino
      DOUBLE PRECISION brstopw(2,5),brstoph(2,5),brststau(2,2),
     .         brstsntau(2,2),brstsmu(2,2),brstsnmu(2),brstsel(2,2),
     .         brstsnel(2),brstbsbst(2,2),brstbbsbt(2,2),
     .         brsttausbnu(2,2),brstelsbnu(2,2),brstupsbdow(2,2),
     .         brst2st1tt,brst2st1startt,brst2st1bb,brst2st1uu,
     .         brst2st1dd,brst2st1ee,brst2st1nunu,brst2st1tautau
*
      DOUBLE PRECISION sbottot(2),sbottot2(2),sbottot3(2)
      DOUBLE PRECISION brsb1neutt(5),brsb2neutt(5),brsb1chart(2),
     .         brsb2chart(2),brsb1hcst(2),brsb2hcst(2),
     .         brsb1glui,brsb2glui,brsb1wst(2),
     .         brsb2wst(2),brsb2H(3),brsb2A(2),brsb2zbot
      DOUBLE PRECISION  brsbstau(2,2),brsbsntau(2,2),brsbsel(2,2),
     .         brsbtstsb(2,2),brsbtbstb(2,2),brsbtaustnu(2,2),
     .         brsbelstnu(2,2),brsbupstdow(2,2),brsbsnel(2),
     .         brsb2sb1bb,brsb2sb1starbb,brsb2sb1tt,
     .         brsb2sb1uu,brsb2sb1dd,brsb2sb1ee,brsb2sb1nunu,
     .         brsb2sb1tautau,brsbsmu(2,2),brsbsnmu(2)
*
      DOUBLE PRECISION gluitot,gluitot2,gluitot3,gluitotrad
      DOUBLE PRECISION brgst1,brgst2,brgsb1,brgsb2,brgsupl,brgsupr,
     .         brgsdownl,brgsdownr,brglnjgluon(5)
      DOUBLE PRECISION brgoup(5),brgoch(5),brgodn(5),brgost(5),
     .         brgotp(5),brgobt(5),brgoud(2),brgocs(2),brgotb(2),
     .         brhcst1b,brwst1b
*
      DOUBLE PRECISION MWNMSSM,dumw,dMW0,DrNMSSM,MWSM,dMWSM,decztt,
     .      deltadecztt,deczee,deltadeczee,BRZTauTau,BRZTTmin,
     .      BRZTTmax,ratio,deltaratio,S2TWeffTau,deltaS2TWeffTau,
     .      S2TWeffSM
*
      COMMON/EWPO/MWNMSSM,dumw,dMW0,DrNMSSM,MWSM,dMWSM,decztt,
     .      deltadecztt,deczee,deltadeczee,BRZTauTau,BRZTTmin,
     .      BRZTTmax,ratio,deltaratio,S2TWeffTau,deltaS2TWeffTau,
     .      S2TWeffSM
*
      COMMON/CHARGINO_WIDTH/chartot2,chartot,chartot3
      COMMON/CHARGINO_BR_2BD/brcharsel,brcharser,brcharsmu1,
     .         brcharsmu2,brcharstau1,brcharstau2,
     .         brcharsne1,brcharsne2,brcharsnm1,brcharsnm2,
     .         brcharsnt1,brcharsnt2,brcharsupl,brcharsupr,
     .         brcharsdownl,brcharsdownr,brcharst1,
     .         brcharst2,brcharsb1,brcharsb2,brcharwneut,
     .         brcharhcneut,brcharzchic,brcharHchic,
     .         brcharAchic
      COMMON/CHARGINO_BR_3BD/brntaunut,brnelnue,brnmunumu,
     .         brnupdb,brnchsb,brntopbb,
     .         brglupdb,brglchsb,brgltopbb,
     .         brchee,brchmumu,brchtautau,brchnene,
     .         brchnmunmu,brchntauntau,brchupup,brchdodo,
     .         brchchch,brchstst,brchtoptop,brchbotbot
*
      COMMON/NEUTRALINO_WIDTH/neuttot2,neuttot,neuttot3,neuttotrad
      COMMON/NEUTRALINO_BR_2BD/brneutst1,brneutst2,brneutsb1,brneutsb2,
     .         brneutsupl,brneutsupr,brneutsdownl,brneutsdownr,
     .         brneutsnel,brneutsn1,brneutsn2,brneutsell,brneutselr,
     .         brneutsnmu,brneutsmu1,brneutsmu2,
     .         brneutstau1,brneutstau2,brneutwchar,brneuthcchar,
     .         brneutzneut,brneutHneut,brneutAneut,brnraddec
      COMMON/NEUTRALINO_BR_3BD/brneutup,brneutdow,brneutch,brneutst,
     .         brneutbot,brneuttop,brneutel,brneutmu,brneuttau,
     .         brneutnue,brneutnumu,brneutnutau,brchubd,brchcbs,
     .         brchtbb,brchelne,brchmunmu,brchtauntau,brglup,brgldo,
     .         brglch,brglst,brgltop,brglbot
*
      COMMON/SLEPTON_WIDTH/selltot,selltot2,selltot3,selrtot,
     .         selrtot2,selrtot3,smu1tot,smu1tot2,smu1tot3,smu2tot,
     .         smu2tot2,smu2tot3,stau1tot2,stau2tot,stau2tot2,
     .         stau2tot3,snelltot,snelltot2,snelltot3,snmu1tot,
     .         snmu1tot2,snmu1tot3,sntautot2,sntautot3,sntautot
      COMMON/SLEPTON_BR_2BD/brsellneute,brselrneute,brsellcharnue,
     .         brselrcharnue,brsnellneut,brsnellchar,brsmu1neutmu,
     .         brsmu2neutmu,brsmu1charnumu,brsmu2charnumu,brsnmu1neut,
     .         brsnmu1char,brstau1neut,brstau2neut,brstau1char,
     .         brstau2char,brstau1hcsn,brstau2hcsn,brstau1wsn,
     .         brstau2wsn,brstau2ztau,brstau2H,brstau2A,brsntauneut,
     .         brsntauchar,brsntau1hcstau,brsntau1wstau
      COMMON/SLEPTON_BR_3BD/brsellstau1star,brsellstau1,
     .         brsellstau1nutau,brselrstau1star,brselrstau1,
     .         brselrstau1nutau,brsnestau1star,brsnestau1,
     .         brsnestau1nutau,brsmu1stau1star,brsmu1stau1,
     .         brsmu1stau1nutau,brsmu2stau1star,brsmu2stau1,
     .         brsmu2stau1nutau,brsnmustau1star,brsnmustau1,
     .         brsnmustau1nutau,brstau2stau1star,brstau2stau1,
     .         brstau2stau1nn,brsntaustau1star,brsntaustau1,
     .         brsntaustau1nutau
*
      COMMON/SQUARK_WIDTH/supltot2,suprtot2,sdowltot2,sdowrtot2
      COMMON/SQUARK_BR_2BD/brsuplnup,brsuplcdow,brsuplglui,
     .         brsuprnup,brsuprcdow,brsuprglui,
     .         brsdowlndow,brsdowlchup,brsdowlglui,
     .         brsdowrndow,brsdowrchup,brsdowrglui
*
      COMMON/STOP_WIDTH/stoptot,stoptot2,stoptot3,stoptotrad
      COMMON/STOP_BR_2BD/brst1neutt,brst2neutt,brst1charb,
     .         brst2charb,brst1hcsb,brst2hcsb,brst1wsb,
     .         brst2wsb,brst1glui,brst2glui,brst2H,brst2A,
     .         brst2ztop,brgamma,brgammaup,brgammagluino
      COMMON/STOP_BR_3BD/brstopw,brstoph,brststau,
     .         brstsntau,brstsmu,brstsnmu,brstsel,
     .         brstsnel,brstbsbst,brstbbsbt,
     .         brsttausbnu,brstelsbnu,brstupsbdow,
     .         brst2st1tt,brst2st1startt,brst2st1bb,brst2st1uu,
     .         brst2st1dd,brst2st1ee,brst2st1nunu,brst2st1tautau
*
      COMMON/SBOTTOM_WIDTH/sbottot,sbottot2,sbottot3
      COMMON/SBOTTOM_BR_2BD/brsb1neutt,brsb2neutt,brsb1chart,
     .         brsb2chart,brsb1hcst,brsb2hcst,
     .         brsb1glui,brsb2glui,brsb1wst,
     .         brsb2wst,brsb2H,brsb2A,brsb2zbot
      COMMON/SBOTTOM_BR_3BD/brsbstau,brsbsntau,brsbsel,
     .         brsbtstsb,brsbtbstb,brsbtaustnu,
     .         brsbelstnu,brsbupstdow,brsbsnel,
     .         brsb2sb1bb,brsb2sb1starbb,brsb2sb1tt,
     .         brsb2sb1uu,brsb2sb1dd,brsb2sb1ee,brsb2sb1nunu,
     .         brsb2sb1tautau,brsbsmu,brsbsnmu
*
      COMMON/GLUINO_WIDTH/gluitot,gluitot2,gluitot3,gluitotrad
      COMMON/GLUINO_BR_2BD/brgst1,brgst2,brgsb1,brgsb2,brgsupl,brgsupr,
     .         brgsdownl,brgsdownr,brglnjgluon
      COMMON/GLUINO_BR_3BD/brgoup,brgoch,brgodn,brgost,brgotp,
     .         brgobt,brgoud,brgocs,brgotb,brhcst1b,brwst1b
*
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
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     . MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     . CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
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
      COMMON/M32/M32,CGR,MPL,GRFLAG
      COMMON/SMODELS/R,CHAN

      IF(IFAIL.NE.0)RETURN

      IRES=23
      NRES=5+IRES

      RES(1)=PAR(3)           !TB
      RES(2)=PAR(20)          !M1
      RES(3)=PAR(21)          !M2
      RES(4)=PAR(22)          !M3
      RES(5)=PAR(12)          !AU3
      RES(6)=PAR(13)          !AD3
      RES(7)=PAR(14)          !AE3
      RES(8)=PAR(25)          !AE2
      RES(9)=DSQRT(PAR(10))   !ML3
      RES(10)=DSQRT(PAR(18))  !ML2
      RES(11)=DSQRT(PAR(11))  !ME3
      RES(12)=DSQRT(PAR(19))  !ME2
      RES(13)=DSQRT(PAR(7))   !MQ3
      RES(14)=DSQRT(PAR(15))  !MQ2
      RES(15)=DSQRT(PAR(8))   !MU3
      RES(16)=DSQRT(PAR(16))  !MU2
      RES(17)=DSQRT(PAR(9))   !MD3
      RES(18)=DSQRT(PAR(17))  !MD2
      RES(19)=PAR(1)          !L
      RES(20)=PAR(2)          !K
      RES(21)=PAR(5)          !AL
      RES(22)=PAR(6)          !AK
      RES(23)=PAR(4)          !MU

      DO I=1,3
       RES(IRES+I)=SMASS(I)
      ENDDO
      DO I=1,2
       RES(IRES+3+I)=PMASS(I)
      ENDDO

      WRITE(16,11)(RES(I),I=1,NRES)
 11   FORMAT(200E14.6)

      END


      SUBROUTINE ERROR(TOT,NTOT,NFAIL)

*********************************************************************
*   Subroutine for the error file. It contains a summary of the scan:
*   Number of points that passed/failed the tests
*   and ranges for scanned parameters that passed the tests
*********************************************************************

      IMPLICIT NONE

      INTEGER I,S,TOT,NTOT,NFAIL(*),GMUFLAG,HFLAG,MWFLAG
      INTEGER CFLAG(6),M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG
      INTEGER MSFLAG,AKFLAG,ALFLAG,OMGFLAG,MAFLAG,MOFLAG

      DOUBLE PRECISION LN,LNN,KN,KNN,TBN,TBNN,MUN,MUNN
      DOUBLE PRECISION ALN,ALNN,AKN,AKNN,XIFN,XIFNN
      DOUBLE PRECISION XISN,XISNN,MUPN,MUPNN,MSPN,MSPNN
      DOUBLE PRECISION M3HN,M3HNN,MAN,MANN,MPN,MPNN
      DOUBLE PRECISION M1N,M1NN,M2N,M2NN,M3N,M3NN
      DOUBLE PRECISION AU3N,AU3NN,AD3N,AD3NN,AE3N,AE3NN,AE2N,AE2NN
      DOUBLE PRECISION ML3N,ML3NN,ML2N,ML2NN,ME3N,ME3NN,ME2N,ME2NN
      DOUBLE PRECISION MQ3N,MQ3NN,MQ2N,MQ2NN,MU3N,MU3NN,MU2N,MU2NN
      DOUBLE PRECISION MD3N,MD3NN,MD2N,MD2NN,DEV

      COMMON/BOUNDS/LN,LNN,KN,KNN,TBN,TBNN,MUN,MUNN,
     . ALN,ALNN,AKN,AKNN,XIFN,XIFNN,
     . XISN,XISNN,MUPN,MUPNN,MSPN,MSPNN,
     . M3HN,M3HNN,MAN,MANN,MPN,MPNN,
     . M1N,M1NN,M2N,M2NN,M3N,M3NN,
     . AU3N,AU3NN,AD3N,AD3NN,AE3N,AE3NN,AE2N,AE2NN,
     . ML3N,ML3NN,ML2N,ML2NN,ME3N,ME3NN,ME2N,ME2NN,
     . MQ3N,MQ3NN,MQ2N,MQ2NN,MU3N,MU3NN,MU2N,MU2NN,
     . MD3N,MD3NN,MD2N,MD2NN
      COMMON/SCANFLAGS/M1FLAG,M2FLAG,M3FLAG,MHDFLAG,MHUFLAG,
     . MSFLAG,AKFLAG,ALFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/MWFLAG/MWFLAG
      COMMON/CFLAG/CFLAG

      WRITE(17,10)"NMSSMTools scan info               "
      WRITE(17,10)"Version number: 6.1.0              "
      WRITE(17,*)
      WRITE(17,*)
      WRITE(17,10)"Number of points:                  "
      WRITE(17,*)
      WRITE(17,10)"  scanned                          ",NTOT
      WRITE(17,10)"  mu=0 or (kappa=0 and Akappa=/=0) ",NFAIL(9)
      S=0
      DO I=1,7
       S=S+NFAIL(I)
      ENDDO
      WRITE(17,10)"  with mh1^2 or ma1^2 or mhc^2 < 0 ",S
      WRITE(17,10)"  with m_sfermion^2 < 0            ",NFAIL(8)
      WRITE(17,10)"  violating constraints            ",NFAIL(10)
      S=NFAIL(11)+NFAIL(12)
      WRITE(17,10)"  RGE integration problem          ",S
      S=NFAIL(13)+NFAIL(14)
      WRITE(17,10)"  convergence problem              ",S
      WRITE(17,*)
      WRITE(17,10)"Remaining good points              ",TOT
      WRITE(17,*)
      WRITE(17,*)

      IF(OMGFLAG.EQ.0 .AND. CFLAG(1).EQ.0 .AND. CFLAG(2).EQ.0 .AND.
     .   CFLAG(3).EQ.0 .AND. CFLAG(4).EQ.0 .AND. CFLAG(5).EQ.0.AND. 
     .   CFLAG(6).EQ.0 .AND. GMUFLAG.EQ.0)THEN
       WRITE(17,20)"Contraints taken into account: none               "
      ELSE
       WRITE(17,20)"Contraints taken into account:                    "
       IF(OMGFLAG.GT.0)
     .  WRITE(17,20)" - Relic density from Planck +/- 10% [0.107,0.131]"
       IF(OMGFLAG.LT.0)
     .  WRITE(17,20)" - Relic density from Planck upper bound < 0.131  "
       IF(IABS(OMGFLAG).EQ.2 .OR. IABS(OMGFLAG).EQ.4)
     .  WRITE(17,20)" - DM direct detection                            "
       IF(IABS(OMGFLAG).EQ.3 .OR. IABS(OMGFLAG).EQ.4)
     .  WRITE(17,20)" - DM indirect detection                          "
       IF(CFLAG(1).NE.0)
     .  WRITE(17,20)" - Landau poles and false minima                  "
       IF(CFLAG(2).NE.0)
     .  WRITE(17,20)" - LEP/Tevatron Higgs+sparticle                   "
       IF(CFLAG(3).NE.0)
     .  WRITE(17,20)" - LHC Higgs                                      "
       IF(CFLAG(4).NE.0)
     .  WRITE(17,20)" - Upsilon, B and K decays                        "
       IF(CFLAG(5).NE.0)
     .  WRITE(17,20)" - CMS charg(neutal)ino                           "
       IF(GMUFLAG.EQ.1)
     .  WRITE(17,20)" - (g-2)_muon                                     "
       IF(MWFLAG.EQ.1)
     .  WRITE(17,20)" - Delta_MW                                       "
       IF(CFLAG(6).LT.0)
     .  WRITE(17,20)" - Squarks and gluinos > 1 TeV                    "
       IF(CFLAG(6).GT.0)
     .  WRITE(17,20)" - LHC SUSY constraints via SmodelS               "
      ENDIF

      IF(TOT.GT.0)THEN

       WRITE(17,*)
       WRITE(17,*)
       WRITE(17,20)"Parameter ranges for good points:"
       WRITE(17,*)
       WRITE(17,30)" TANB: ",TBN,TBNN,DEV(TBN,TBNN)
       IF(M1FLAG.NE.0)THEN
        WRITE(17,30)" M1: ",M1N,M1NN,DEV(M1N,M1NN)
       ENDIF
       WRITE(17,30)" M2: ",M2N,M2NN,DEV(M2N,M2NN)
       IF(M3FLAG.NE.0)THEN
        WRITE(17,30)" M3: ",M3N,M3NN,DEV(M3N,M3NN)
       ENDIF
       WRITE(17,30)" LAMBDA: ",LN,LNN,DEV(LN,LNN)
       WRITE(17,30)" KAPPA: ",KN,KNN,DEV(KN,KNN)
       WRITE(17,30)" MUEFF: ",MUN,MUNN,DEV(MUN,MUNN)
       IF(MOD(MAFLAG,3).EQ.0)THEN
        WRITE(17,30)" ALAMBDA: ",ALN,ALNN,DEV(ALN,ALNN)
        WRITE(17,30)" MA: ",MAN,MANN
        WRITE(17,40)"(MA is not an input parameter)"
        IF(XIFN.NE.0d0 .OR. XIFNN.NE.0d0)
     .   WRITE(17,30)" XIF: ",XIFN,XIFNN,DEV(XIFN,XIFNN)
       ELSEIF(MOD(MAFLAG,3).EQ.1)THEN
        WRITE(17,30)" ALAMBDA: ",ALN,ALNN
        WRITE(17,40)"(ALAMBDA is not an input parameter)"
        WRITE(17,30)" MA: ",MAN,MANN,DEV(MAN,MANN)
        IF(XIFN.NE.0d0 .OR. XIFNN.NE.0d0)
     .   WRITE(17,30)" XIF: ",XIFN,XIFNN,DEV(MAN,MANN)
       ELSE
        WRITE(17,30)" ALAMBDA: ",ALN,ALNN,DEV(ALN,ALNN)
        WRITE(17,30)" MA: ",MAN,MANN,DEV(MAN,MANN)
        WRITE(17,30)" XIF: ",XIFN,XIFNN
        WRITE(17,40)"(XIF is not an input parameter)"
       ENDIF
       IF(MAFLAG/3.EQ.0)THEN
        WRITE(17,30)" AKAPPA: ",AKN,AKNN,DEV(AKN,AKNN)
        WRITE(17,30)" MP: ",MPN,MPNN
        WRITE(17,40)"(MP is not an input parameter)"
        IF(XISN.NE.0d0 .OR. XISNN.NE.0d0)
     .   WRITE(17,30)" XIS: ",XISN,XISNN,DEV(XISN,XISNN)
       ELSEIF(MAFLAG/3.EQ.1)THEN
        WRITE(17,30)" AKAPPA: ",AKN,AKNN
        WRITE(17,40)"(AKAPPA is not an input parameter)"
        WRITE(17,30)" MP: ",MPN,MPNN,DEV(MPN,MPNN)
        IF(XISN.NE.0d0 .OR. XISNN.NE.0d0)
     .   WRITE(17,30)" XIS: ",XISN,XISNN,DEV(XISN,XISNN)
       ELSE
        WRITE(17,30)" AKAPPA: ",AKN,AKNN,DEV(AKN,AKNN)
        WRITE(17,30)" MP: ",MPN,MPNN,DEV(MPN,MPNN)
        WRITE(17,30)" XIS: ",XISN,XISNN
        WRITE(17,40)"(XIS is not an input parameter)"
       ENDIF
       IF(MUPN.NE.0d0 .OR. MUPNN.NE.0d0)
     .  WRITE(17,30)" MUP: ",MUPN,MUPNN,DEV(MUPN,MUPNN)
       IF(MSPN.NE.0d0 .OR. MSPNN.NE.0d0)
     .  WRITE(17,30)" MSP: ",MSPN,MSPNN,DEV(MSPN,MSPNN)
       IF(M3HN.NE.0d0 .OR. M3HNN.NE.0d0)
     .  WRITE(17,30)" M3H: ",M3HN,M3HNN,DEV(M3HN,M3HNN)
       WRITE(17,30)" AU3: ",AU3N,AU3NN,DEV(AU3N,AU3NN)
       WRITE(17,30)" AD3: ",AD3N,AD3NN,DEV(AD3N,AD3NN)
       WRITE(17,30)" AE3: ",AE3N,AE3NN,DEV(AE3N,AE3NN)
       WRITE(17,30)" AE2: ",AE2N,AE2NN,DEV(AE2N,AE2NN)
       WRITE(17,30)" ML3: ",ML3N,ML3NN,DEV(ML3N,ML3NN)
       WRITE(17,30)" ML2: ",ML2N,ML2NN,DEV(ML2N,ML2NN)
       WRITE(17,30)" ME3: ",ME3N,ME3NN,DEV(ME3N,ME3NN)
       WRITE(17,30)" ME2: ",ME2N,ME2NN,DEV(ME2N,ME2NN)
       WRITE(17,30)" MQ3: ",MQ3N,MQ3NN,DEV(MQ3N,MQ3NN)
       WRITE(17,30)" MQ2: ",MQ2N,MQ2NN,DEV(MQ2N,MQ2NN)
       WRITE(17,30)" MU3: ",MU3N,MU3NN,DEV(MU3N,MU3NN)
       WRITE(17,30)" MU2: ",MU2N,MU2NN,DEV(MU2N,MU2NN)
       WRITE(17,30)" MD3: ",MD3N,MD3NN,DEV(MD3N,MD3NN)
       WRITE(17,30)" MD2: ",MD2N,MD2NN,DEV(MD2N,MD2NN)

      ENDIF

 10   FORMAT(A35,I10)
 20   FORMAT(A50)
 30   FORMAT(A11,3E15.4)
 40   FORMAT(A36)

      END
