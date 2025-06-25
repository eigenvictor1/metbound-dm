(* ::Package:: *)

(* ::Text:: *)
(*RGE related functions*)


RGEconstructor[params_,betas_]:= Map[D[#[[1]], x] == #[[2]]&, Thread@{params, betas}];
RGEconstructor::usage = "RGEconstructor[params_,betas_]: builds RGEs from list of parameters and beta functions";


IniCondsConstructor[params_, inis_, xini_]:= Map[#[[1]] == #[[2]]/.{x->xini}&, Thread@{params, inis}];
IniCondsConstructor::usage = "IniCondsConstructor[params_, inis_, xini_]: builds list enforcing initial conditions at x = xini";


RGEintegrate[params_, betas_,inis_, xini_, xend_, conds_]:=
Module[{eqs},
eqs = Join[
			RGEconstructor[params, betas],               (*RGEs*)
			IniCondsConstructor[params, inis, xini],    (*ini conds*)
			conds                                        (*extra conds*)
		];
(*Print[IniCondsConstructor[params, inis, xini]];
Print[params];
Print[eqs];*)
	Flatten@params/.NDSolve[eqs, params, {x, xini, xend}]
];
RGEintegrate::usage = "RGEintegrate[params_, betas_,inis_, xini_, xend_, conds_]: outputs RGE evolution from xini to xend for given inputs";


RGEmerge[pieces_, listthres_]:= 
Module[{ntotalparams, paddedpieces},
ntotalparams = TakeLargest[Map[Length, pieces],1][[1]];
paddedpieces = Map[Join[#,ConstantArray[-10^3, ntotalparams-Length[#]]]&, pieces];
Sum[paddedpieces[[i]]*HeavisideTheta[x-listthres[[i]]]HeavisideTheta[listthres[[i+1]]-x],{i, 1, Length[paddedpieces]}]
]
RGEmerge::usage = "RGEmerge[pieces_, listthres_]: gives total evolution from n pieces and n-1 thresholds. First in listthres is xini, last is xend";


CalculateInstabilityScale[Lwith1Loop_, xsearch_]:= x/.FindRoot[Lwith1Loop == 0, {x, xsearch}];
CalculateInstabilityScale::usage = "CalculateInstabilityScale[int_, xsearch_] does what the name says using result form NDSolve including 1 loop corrections";


CalculateMetastabilityBound[Lwith1Loop_,xsearch_]:=Module[
{muIlog},
muIlog = CalculateInstabilityScale[Lwith1Loop, xsearch];
Sqrt@(Abs[D[Lwith1Loop,x]/.{x->muIlog}] E^{-3/2} (E^muIlog)^2)//N
][[1]];
CalculateMetastabilityBound::usage = "CalculateMetastabilityBound[Lwith1Loop_, xend_] does what the name says";


(* ::Code:: *)
(*(*DEBUG*)*)


(* ::Code:: *)
(*paramsbelow = {L[x], Y[x]}*)
(*newparams = {t2[x]}*)
(*Join[paramsbelow, newparams]*)
(*funcs[L_] := {L[[1]], L[[2]](1+(4Pi)^-2 L[[3]]), L[[3]]};*)
(*AddThresholdCorrections[paramsbelow, newparams, funcs]*)


(* ::Code:: *)
(*params = {L[x], Y[x]}*)
(*betas = {(4 Pi)^(-2)(24 L[x]^2), (4 Pi)^(-2)Y[x]}*)
(*inis = {1, 0};*)
(*conds = {};*)
(*xini = 1; *)
(*xend = 10;*)
(*RGEconstructor[params, betas]*)
(*IniCondsConstructor[params, inis, xini] *)
(*aux = RGEintegrate[params, betas, inis,  xini, xend, conds]*)
(*y[xM]==YSM[xM](1+(4Pi)^-2 t2in),b[xM]==BSM[xM],t[xM]==TSM[xM],L[xM]==lSM[xM](1+(4Pi)^-2 t2in)+(4Pi)^-2 t2in^2,g[xM]==GSM[xM],h[xM]==HSM[xM],s[xM]==SSM[xM],t2[xM]==t2in,*)


(* ::Code:: *)
(*Plot[aux[[1]], {x, xini, xend}]*)
