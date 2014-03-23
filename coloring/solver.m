(* ::Package:: *)

Off[General::compat];
Needs["Combinatorica`"];

LoadAndColor[file_String] := Module[
  {data = Import[file, "Table"][[2;;-1]], g, am},
  am = Normal@SparseArray[data + 1 -> Table[1, {i, 1, Length@data}],{Max@data + 1,Max@data + 1}];
  g = FromAdjacencyMatrix[am, Type -> Undirected];
  VertexColoring[g, Algorithm -> Optimum]];

c = LoadAndColor[$CommandLine[[-1]]];
WriteString[$Output, Length@DeleteDuplicates@c];
WriteString[$Output, " 0\n"];
WriteString[$Output, #, " "]& /@ c;
WriteString["\n"];
