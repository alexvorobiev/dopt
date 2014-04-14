(* ::Package:: *)

Off[General::compat];
Needs["Combinatorica`"];

Tsp[file_String] := Module[
  {data = Import[file, "Table"][[2;;-1]]},
  FindShortestTour[data]];

c = Tsp[$CommandLine[[-1]]];
WriteString[$Output, First@c];
WriteString[$Output, " 0\n"];
WriteString[$Output, #, " "]& /@ Rest@c;
WriteString["\n"];
