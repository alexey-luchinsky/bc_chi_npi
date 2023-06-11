<< MyPackages`HistV2`

Clear[extractFF];
extractFF::usage = "extractFF[json, num]";
extractFF[json_, n_, fact_ : 1] := Module[{name, data},
  name = json[[3, 2, n, 1, 2]];
  data = json[[3, 2, n, 5, 2, All, 3, 2]];
  data[[All, 2]] = fact*data[[All, 2]];
  Return[<|"name" -> name, "data" -> data|>]
  ]

Clear[extractNames];
extractNames[json_] := 
  Module[{names = json[[3, 2, All, 1, 2]]}, 
   Table[{i, names[[i]]}, {i, Length[names]}]];

fitFunc = a0 + a1 q2 + a2 q2^2 + a3 q2^3;

fitParams = DeleteCases[Variables[fitFunc] , q2];


Clear[getFileNameName];
Options[getFileName] = {Var -> "m2_0m1"};
getFileName[out_, mode_, ff_, OptionsPattern[]] := 
  "../c++/results/B_c+_" <> out <> "_" <> mode <> "_ff" <> ff <> "/" <>
    OptionValue[Var] <> ".txt";
getFileName[out_, mode_, opt : OptionsPattern[]] := 
  getFileName[out, mode, "1", opt];


Clear[readROOT];
Options[readROOT] = {Norm -> False, Color -> Red};
readROOT[fileName_, OptionsPattern[]] := Module[{root},
  root = ReadList[fileName, {Number, Number, Number}];
  If[Head[root] =!= List,
   Print["File not found",  fileName];
   Abort[];
   ];
  root = HST2D[{#[[1]], #[[2]] \[PlusMinus] #[[3]]} & /@ root, 
    PlotStyle -> OptionValue[Color]];
  If[OptionValue[Norm] =!= False,
   root *= OptionValue[Norm]/IntegralHST[root];
   ];
  root
  ]

Clear[NormalizeHST];
Options[NormalizeHST] = {Norm -> 1};
NormalizeHST[hist_, OptionsPattern[]] := 
  OptionValue[Norm]*hist/IntegralHST[hist];

Clear[SaveVar];
Options[SaveVar] = {Overwrite -> True};
SaveVar[fileName_, var_, OptionsPattern[]] := Module[{file},
  If[OptionValue[Overwrite] === True,
   file = OpenWrite[fileName];,
   file = OpenAppend[fileName];
   ];
  WriteString[file, ToString[Definition[var], InputForm] <> "\n"];
  Close[file];
  ]

