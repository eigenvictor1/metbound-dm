(* ::Package:: *)

(* ::Text:: *)
(*Fancy pre-formatted versions of the function plot*)


Clear[FancyPlot]
FancyPlot[func_, rangex_, rangey_, labels_, ticks_:Automatic, aspectratio_:1, imagesize_:{400, Automatic}]:=
Plot[func, Flatten@{x, rangex},
PlotRange->{rangex, rangey}, 
Frame->True,
FrameLabel-> labels,
FrameTicks->ticks,
FrameTicksStyle->Directive[Black, FontSize->12],
AspectRatio->aspectratio,
ImageSize->imagesize
]
FancyPlot::usage = "FancyPlot[func_, rangex_, rangey_, labels_, ticks_:Automatic, aspectratio_:1, imagesize_:350]: nicer plot";



FancyListContourPlot[list_, contours_,colorfunction_, xrange_, yrange_, labels_, ticks_:Automatic, aspectratio_:1, imagesize_:{400, Automatic}] := 
ListContourPlot[list, 
Contours-> contours,
ColorFunction-> colorfunction,
PlotRange-> {xrange, yrange},
PlotRangeClipping-> True,
ImagePadding-> Automatic,
PlotRangePadding-> 0,
PerformanceGoal-> "Quality",
Frame-> True, 
FrameLabel-> labels,
FrameTicksStyle-> Directive[Black, FontSize->13],
FrameTicks-> ticks,
AspectRatio-> aspectratio,
ImageSize-> imagesize
]
FancyListContourPlot::usage = "FancyListContourPlot[list_, contours_, colorfunction_, xrange_, yrange_, labels_, ticks_:Automatic, aspectratio_:1, imagesize_:{400, Automatic}]: pre-formated contour plot";


(* ::Code:: *)
(*(*DEBUG*)*)


(* ::Code:: *)
(*func = x^2;*)
(*rangex = {0, 1};*)
(*rangey = {0, 1};*)
(*labels = {"x", "y"};*)
(*FancyPlot[func, rangex, rangey, labels]*)
