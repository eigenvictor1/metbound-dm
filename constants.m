(* ::Package:: *)

(* ::Text:: *)
(*Known  physical  quantities  to  be  used  globally*)


(*scales*)
mp = 1.22*10^19;             mp::usage = "plank mass";
xp =  Log@1.22*10^19;        xp::usage  = "log of the plank mass";
mt = 173.1;                  mt::usage  = "mass of the top";
xt = Log@mt;                 xt::usage  = "log of the mass of the top" ; 
vH = 246;                    vH::usage  = "mass of the higgs";
(*gauge couplings*)
gpt = 0.349946(*0.35830 0.358545*);              gpt::usage  = "hypercharge gauge coupling at x = mt" ;
gt  = 0.652825(*0.64779 0.64765*);               gt::usage = "electroweak gauge coupling at x = mt";
gst = 1.1650(*1.1666 1.1618*);                gst::usage = "strong gauge coupling at x = mt";
(*Higgs quartic coupling*)
Lt  = 0.12987(*0.12604 0.12607*);               Lt::usage  = "higgs quartic coupling at x = mt";
(*yukawas*)
YTt  = 0.9373(* 0.93690 161.98/246 Sqrt[2]*);   YTt::usage = "yukawa of the top at x = mt";
YBt  = 0.10206(*4.2*Sqrt[2]/246.2 2.702/246 Sqrt[2]*);   YBt::usage = "yukawa of the bottom at x = mt";
YTAt = 0.024026 (*1.777*Sqrt[2]/246.2 1.78412/246 Sqrt[2]*);  YTAt::usage = "yukawa of the tau at x = mt";
