(* ::Package:: *)

(*Standard Model*)
paramsSM = {L[x], YT[x], YB[x], YTA[x], g[x], gp[x], gs[x]};
inisSM = {Lt, YTt, YBt, YTAt, gt, gpt, gst};
paramsSMSubs =  {L -> L[x], YT -> YT[x], YB -> YB[x], YTA -> YTA[x], g -> g[x], gp -> gp[x], gs -> gs[x] };


(*Generate beta functions for specific model - requires RGBeta package*)
ResetModel[]
(*STANDARD MODEL*)

AddGaugeGroup[gp,U1Y,U1]
AddGaugeGroup[g,SU2L,SU[2]]
AddGaugeGroup[gs,SU3c,SU[3]]

Dim[gen]=ng(*3*);

AddFermion[q,GaugeRep->{U1Y[1/6],SU2L[fund],SU3c[fund]},FlavorIndices->{gen}]
AddFermion[u,GaugeRep->{U1Y[-2/3],Bar@SU3c[fund]},FlavorIndices->{gen}]
AddFermion[d,GaugeRep->{U1Y[1/3],Bar@SU3c[fund]},FlavorIndices->{gen}]
AddFermion[l,GaugeRep->{U1Y[-1/2],SU2L[fund]},FlavorIndices->{gen}]
AddFermion[e,GaugeRep->{U1Y[1]},FlavorIndices->{gen}]

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

SetReal[YT,YB,YTA]


Keys@$couplings
CheckProjection/@%


betasSM =
{ 
	Finalize[BetaFunction[L,2], Parameterizations->matrixSubs], 
	Finalize[BetaFunction[yu,2], Parameterizations->matrixSubs][[3,3]], 
	Finalize[BetaFunction[yd,2], Parameterizations->matrixSubs][[3,3]], 
	Finalize[BetaFunction[ye,2], Parameterizations->matrixSubs][[3,3]], 
	Finalize[BetaFunction[g,3], Parameterizations->matrixSubs], 
	Finalize[BetaFunction[gp,3], Parameterizations->matrixSubs], 
	Finalize[BetaFunction[gs,3], Parameterizations->matrixSubs]
}/.paramsHiggsinoBinoSubs;


paramsSM::usage = "running parameters in the SM";
inisSM::usage = "values of the SM running parameters at x = mtop";
betasSM::usage = "beta functions in the SM";
