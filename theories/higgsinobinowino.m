(* ::Package:: *)

paramsHiggsinoBino = { L[x], YT[x], YB[x], YTA[x], g[x], gp[x], gs[x], yL[x], yR[x], yW[x] };
paramsHiggsinoBinoSubs =  { L -> L[x], YT -> YT[x], YB -> YB[x], YTA -> YTA[x], g -> g[x], gp -> gp[x], gs -> gs[x], yL -> yL[x], yR -> yR[x], yW -> yW[x] };


modelname = "higgsinobinowino";


(*Generate beta functions for specific model - requires RGBeta package*)
ResetModel[]
(*STANDARD MODEL*)

AddGaugeGroup[gp,U1Y,U1]
AddGaugeGroup[g,SU2L,SU[2]]
AddGaugeGroup[gs,SU3c,SU[3]]

Dim[gen]=ng(*3*);

AddFermion[ q, GaugeRep->{U1Y[1/6],SU2L[fund],SU3c[fund]},FlavorIndices->{gen}]
AddFermion[ u, GaugeRep->{U1Y[-2/3],Bar@SU3c[fund]},FlavorIndices->{gen}]
AddFermion[ d, GaugeRep->{U1Y[1/3],Bar@SU3c[fund]},FlavorIndices->{gen}]
AddFermion[ l, GaugeRep->{U1Y[-1/2],SU2L[fund]},FlavorIndices->{gen}]
AddFermion[ e, GaugeRep->{U1Y[1]},FlavorIndices->{gen}]

AddScalar[H,GaugeRep->{U1Y[1/2],SU2L[fund]}]

AddYukawa[yu,{H,q,u},
GroupInvariant->(del[SU3c@fund,#2,#3]eps[SU2L@fund,#2,#1]&),
CouplingIndices->({gen[#2],gen[#3]}&),
Chirality->Right]
AddYukawa[yd,{Bar@H,q,d},
GroupInvariant->(del[SU3c@fund,#2,#3]del[SU2L@fund,#1,#2]&),
CouplingIndices->({gen[#2],gen[#3]}&),
Chirality->Right]
AddYukawa[ye,{Bar@H,l,e},
GroupInvariant->(del[SU2L@fund,#1,#2]&),
CouplingIndices->({gen[#2],gen[#3]}&),
Chirality->Right]

matrixSubs={yu->DiagonalMatrix@{0,0,YT}, yd->DiagonalMatrix@{0,0,YB}, ye->DiagonalMatrix@{0,0,YTA},ng->3};

AddQuartic[L,{Bar@H,H,Bar@H,H},
GroupInvariant->((del[SU2L@fund,#1,#2]del[SU2L@fund,#3,#4])&)]

AddScalarMass[MH2,{Bar@H,H},
GroupInvariant->(del[SU2L@fund,#1,#2]&)]

(*NEW PHYSICS*)
AddFermion[psiL]
AddFermion[chiL,GaugeRep->{U1Y[1/2],SU2L[fund]}]
AddFermion[chiRBar,GaugeRep->{U1Y[-1/2],SU2L[fund]}]
AddFermion[chiW,GaugeRep->{U1Y[-1/2],SU2L[fund]}]

AddYukawa[yL,{Bar@H,chiL,psiL},
GroupInvariant->(del[SU2L@fund,#1,#2]&),
Chirality->Right]
AddYukawa[yR,{H,chiRBar,psiL},
GroupInvariant->(del[SU2L@fund,#1,#2]&),
Chirality->Right]
AddYukawa[yW,{H,chiW,psiL},
GroupInvariant->(del[SU2L@fund,#1,#2]&),
Chirality->Right]
AddFermionMass[mS,{psiL,psiL}]
AddFermionMass[mD,{chiL,chiRBar},
GroupInvariant->(del[SU2L@fund,#1,#2]&),
Chirality->Right]
AddFermionMass[mW,{chiL,chiRBar},
GroupInvariant->(del[SU2L@fund,#1,#2]&),
Chirality->Right]

SetReal[YT,YB,YTA,yL,yR, yW]


Keys@$couplings
CheckProjection/@%


betasHiggsinoBino =
{ 
	Finalize[BetaFunction[ L,  2 ], Parameterizations->matrixSubs], 
	Finalize[BetaFunction[ yu, 2 ], Parameterizations->matrixSubs][[3,3]], 
	Finalize[BetaFunction[ yd, 2 ], Parameterizations->matrixSubs][[3,3]], 
	Finalize[BetaFunction[ ye, 2 ], Parameterizations->matrixSubs][[3,3]], 
	Finalize[BetaFunction[ g,  3 ], Parameterizations->matrixSubs], 
	Finalize[BetaFunction[ gp, 3 ], Parameterizations->matrixSubs], 
	Finalize[BetaFunction[ gs, 3 ], Parameterizations->matrixSubs], 
	Finalize[BetaFunction[ yL, 2 ], Parameterizations->matrixSubs], 
	Finalize[BetaFunction[ yR, 2 ], Parameterizations->matrixSubs],
	Finalize[BetaFunction[ yW, 2 ], Parameterizations->matrixSubs]
}/.paramsHiggsinoBinoSubs;


ThresholdHiggsinoBino[{L_, YT_, YB_, YTA_, g_, gp_, gs_, yL_, yR_, yW_, mD_, mS_}]:=
	Module[{x},
	x = Min[{mD, mS}]/Max[{mD, mS}];
	{
		L+1/(16 \[Pi]^2 (-1+x^2)^3) (x^2 (-3+2 x^2+x^4-x^2 (-5+x^2) Log[1/x^2]) yL^4+8 x (-1+x^4+2 x^2 Log[1/x^2]) yL^3 yR+L(-1+7 x^2-7 x^4+x^6+2 x^4 (-3+x^2) Log[1/x^2]) yR^2+x^2 (-3+2 x^2+x^4-x^2 (-5+x^2) Log[1/x^2]) yR^4-4 x (-1+x^4+2 x^2 Log[1/x^2]) yL yR (L-2 yR^2)+yL^2 (L (-1+7 x^2-7 x^4+x^6+2 x^4 (-3+x^2) Log[1/x^2])-2 (2+7 x^2-8 x^4-x^6+x^2 (-6-7 x^2+x^4) Log[1/x^2]) yR^2)),
		YT -YT/(64 \[Pi]^2 (-1+x^2)^3) ((-1+7 x^2-7 x^4+x^6+2 x^4 (-3+x^2) Log[1/x^2]) yL^2-4 x (-1+x^4+2 x^2 Log[1/x^2]) yL yR+(-1+7 x^2-7 x^4+x^6+2 x^4 (-3+x^2) Log[1/x^2]) yR^2) ,
		YB, 
		YTA, 
		g, 
		gp, 
		gs, 
		yL,
		yR
	}//N
	];


massvals[md_,ms_,yL_,yR_,x_]:=(Eigenvalues[{{0, 0, - bet Exp[x]},{0, -0, - alph Exp[x]}, {-bet Exp[x], -  alph Exp[x], -0}}/.{alph->(yL+yR)/2,bet->(yL-yR)/2}]//N)^2; 
LOneLoopHiggsinoBino[{L_, YT_, YB_, YTA_, g_, gp_, gs_, yL_, yR_, yW_, mD_, mS_, x_}]:=(4/(2048 Pi^2) (-5 gp^4+6 (gp^2+g^2)^2 Log[(gp^2+g^2)/4]-10 gp^2 g^2-15 gp^4+12 g^4 Log[g^2/4]+144 YT^4-96 YT^4 Log[YT^2/2])-1/(8\[Pi]^2) Sum[(massvals[mD,mS,yL,yR,x][[i]])^2/Exp[4x](Log[(massvals[mD,mS,yL,yR,x][[i]]+0.0000000001)/Exp[2x]]-3/2),{i,1,3}])


OrderedEigenvalsDoubletSinglet[{y_, theta_, mD_, mS_}]:=
Module[
{alph, bet, M, m},
alph=N[y/2(Sin[theta]+Cos[theta])];
bet=N[y/2(Sin[theta]-Cos[theta])];
M = {{mD, 0, - bet v},{0, -mD, -alph v}, {-bet v, - alph v, -mS}}/.v->vH;
m=Eigenvalues[M];
Abs[m[[Ordering[Abs[m]]]]]
]


paramsHiggsinoBino::usage = "running parameters: SM + Higgsino + Bino";
betasHiggsinoBino::usage = "beta functions: SM + Higgsino + Bino";
ThresholdHiggsinoBino::usage = "ThresholdHiggsinoBino[{L_, YT_, YB_, YTA_, g_, gp_, gs_, yL_, yR_, mS_, mD_}]: adds threshold corrections on top of integrated SM variables";
LOneLoopHiggsinoBino::usage = "LOneLoopHiggsinoBino[{L_, YT_, YB_, YTA_, g_, gp_, gs_, yL_, yR_, mD_, mS_, x_}]: adds one loop corrections IMPORTANT AROUND INST SCALE";
