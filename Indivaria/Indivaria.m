(* Mathematica Package *)

(* Author: Sompob Saralamba *)

(* Mathematical Modelling Team, Mahidol-Oxford Tropical Medicine Research Unit *)


BeginPackage["Indivaria`",{"Indivaria`BasicFunc`"}]
(* Exported symbols added here with SymbolName::usage *) 

IDVL::usage= "IDVL[input.in] Indivaria main function"
RI::usage="RI[ec50,emax,gamma,1/alpha] calculate the resistance index."
ListPar::usage="ListPar[modelID,outform] shows the list of the parameters needed by the model."
ParaFitD::usage = "ParaFitD[parafile_, concfile_, initN_, PMR_, mu_, sigma_,
	LC_, KZ_, everyH_, Ndrug_, gamma_, ec50_, emin_, emax_, T_, dorfrac_, dortime_, output_]"

Sexy8::usage="[initN_, PMR_, mu_, sigma_, hours_Integer, KZ_List, concdatafile_, 
   everyH_, Ndrug_, gamma_List, ec50_List, emin_List, emax_List, T_, runMax_, 
   dorfrac_, dortime_, outfile_]"


Begin["`Private`"]
(* Implementation of the package *)

(*
DistributeN[initN_, LC_Integer, mu_, sigma_] :=
(*DistributeN[initN,LifeCycle,mu,sigma] distributes the initial number over the life-cycle (LC)
with the mean mu and SD sigma hours.*)
  Module[{distr, x},
   distr = Table[PDF[NormalDistribution[mu, sigma], x]//N, {x, 1, LC}];
   Developer`ToPackedArray[(initN/Total[distr])*distr]
];
   
   
Shiftonehour[ls_List, PMR_Integer] := Module[{tmp},
   (*shift the list of number to the right with a multiplication pmr*)
   tmp = ls;
   tmp = RotateRight[tmp];
   tmp = ReplacePart[tmp, 1 -> tmp[[1]]*PMR];
   Developer`ToPackedArray[tmp]
   ];


Eff[c_, gamma_, ec50_, emin_, emax_] := 
	(*calculate the percentage of killing*)
 emin + (emax - emin) *c^gamma/(c^gamma + ec50^gamma)


PRingFunc[i_Integer, a1_Integer, a2_Integer] := 
	(* inferred from Fig2 of "Febrile temperatures induce cytoadherence of 
	ring-stage Plasmodium falciparum-infected erythrocytes" by
	Rachanee Udomsangpetch et al. (http://www.pnas.org/content/99/18/11825)*)
  Piecewise[{{1, i < a1}, {Exp[Log[0.5]*(i - a1)/(a2 - a1)], i >= a1}}];


LsDot[ls1_List, ls2_List] := Module[{tmp},
   (* dot vector *)
   If[Length@ls1 == Length@ls2, 
   		tmp = ls1*ls2;
   		, 
    	(*Print["Please check your input."];*)
    	tmp = {};
    ];
   Developer`ToPackedArray[tmp]
   ];
*)

(*
CountRing[ls_List, Age_List, fn_Integer] := 
	(*count parasites at Rings*)
  Module[{AGEbegin, AGEend, tot, pls, i},
   (** fn = 0 not using the partition function, 
   		  = 1 use the partition function (PRingFunc) 
     Age = {AGEbegin, AGEend} age range of ring parasites
   		  ***)
    AGEbegin = Age[[1]];
    AGEend = Age[[2]];
    
    pls = Table[PRingFunc[i, AGEbegin, AGEend], {i, 1, Length@ls}];
   
   If[fn == 1,
    tot = ls.pls;
    ,
    tot = ls[[AGEbegin ;; AGEend]] // Total;
    ];
   tot
   ];
*)

ConcMod[datname_String, everyH_Integer, Ndrug_Integer, gamma_List, ec50_List, emin_List, emax_List] := 
  (* calculate the drug concentration (from fitting) and its efficacy for Ring, Trophozoite and Schizont stages *)
  Module[{i, j, nd, dat, maxpoint, part1, part2, fitke, modka, modke, 
    ke, a, fn, endtime, t, conc, concls, TIMEMAX},
   
   TIMEMAX = 1000; (**maximum time for running the models*)
   dat = Import[datname];
   
   (*position of Cmax*)
   maxpoint = Position[dat, Max[dat]];
   endtime = Last[dat];
   
   (*absorbtion*)
   part1 = dat[[1 ;; maxpoint[[1, 1]]]];
   
   (*elimination*)
   part2 = dat[[maxpoint[[1, 1]] ;; Length[dat]]];
   
   (**modka=Fit[part1,{1,t,t^2},t];**)
   (*modka=Fit[part1,{1,t},t];*)
   (*using straight line*)
   modka = Fit[part1, {t}, t];
   
   fitke = FindFit[part2, a*Exp[-ke*t], {a, ke}, t];
   modke = Function[{t}, Evaluate[a*Exp[-ke*t] /. fitke]];
   
   (*output equations*)
   fn = Piecewise[{{modka, 
       0 <= t < dat[[maxpoint[[1, 1]], 1]]}, {modke[t], 
       t >= dat[[maxpoint[[1, 1]], 1]]}}];
   
   (*Table[{i,fn/.t->i},{i,0,Ceiling[endtime[[1]]]}]*)
   
   j = 0;
   nd = 1;
   concls = {};
   
   conc = 0; (*init concentration*)
   
   For[i = 0, i < Ndrug*everyH, i = i + 1,
    If[j != everyH,
     conc = fn /. t -> j;
     
     AppendTo[concls, {i, conc,
       Eff[conc, gamma[[1]], ec50[[1]], emin[[1]], emax[[1]]], 
       Eff[conc, gamma[[2]], ec50[[2]], emin[[2]], emax[[2]]], 
       Eff[conc, gamma[[3]], ec50[[3]], emin[[3]], emax[[3]]]}];
     ];
    
    If[j == everyH && nd < Ndrug,
     nd = nd + 1;
     conc = conc + fn /. t -> j;
     AppendTo[concls, {i, conc,
       Eff[conc, gamma[[1]], ec50[[1]], emin[[1]], emax[[1]]], 
       Eff[conc, gamma[[2]], ec50[[2]], emin[[2]], emax[[2]]], 
       Eff[conc, gamma[[3]], ec50[[3]], emin[[3]], emax[[3]]]}];
     j = 0;
     ];
    (*update ndrug*)
    j = j + 1; 
    ];
   For[i = Ndrug*everyH; j = everyH, i <= TIMEMAX, i = i + 1,
    conc = fn /. t -> j;
    AppendTo[concls, {i, conc,
      Eff[conc, gamma[[1]], ec50[[1]], emin[[1]], emax[[1]]], 
      Eff[conc, gamma[[2]], ec50[[2]], emin[[2]], emax[[2]]], 
      Eff[conc, gamma[[3]], ec50[[3]], emin[[3]], emax[[3]]]}
     ];
    j = j + 1;
    ];

   Developer`ToPackedArray[concls]
   
   ];


ConcMod2[datname_String, everyH_Integer, Ndrug_Integer, ec50_Real] := 
  (* calculate the drug concentration (from fitting) and its efficacy (for Hoshen's model) *)
  Module[{i, j, nd, dat, maxpoint, part1, part2, fitke, modka, modke, 
    ke, a, fn, endtime, t, conc, concls, TIMEMAX},
   
   TIMEMAX = 1000; (**maximum time for running the models*)
   dat = Import[datname];
   
   (*position of Cmax*)
   maxpoint = Position[dat, Max[dat]];
   endtime = Last[dat];
   
   (*absorbtion*)
   part1 = dat[[1 ;; maxpoint[[1, 1]]]];
   
   (*elimination*)
   part2 = dat[[maxpoint[[1, 1]] ;; Length[dat]]];
   
   (**modka=Fit[part1,{1,t,t^2},t];**)
   (*modka=Fit[part1,{1,t},t];*)
   (*using straight line*)
   modka = Fit[part1, {t}, t];
   
   fitke = FindFit[part2, a*Exp[-ke*t], {a, ke}, t];
   modke = Function[{t}, Evaluate[a*Exp[-ke*t] /. fitke]];
   
   (*output equations*)
   fn = Piecewise[{{modka, 
       0 <= t < dat[[maxpoint[[1, 1]], 1]]}, {modke[t], 
       t >= dat[[maxpoint[[1, 1]], 1]]}}];
   
   (*Table[{i,fn/.t->i},{i,0,Ceiling[endtime[[1]]]}]*)
   
   j = 0;
   nd = 1;
   concls = {};
   
   conc = 0; (*init concentration*)
   
   For[i = 0, i < Ndrug*everyH, i = i + 1,
    If[j != everyH,
     conc = fn /. t -> j;
     
     AppendTo[concls, {i, conc, HoshenPDEff[conc, ec50]}];
     ];
    
    If[j == everyH && nd < Ndrug,
     nd = nd + 1;
     conc = conc + fn /. t -> j;
     AppendTo[concls, {i, conc, HoshenPDEff[conc, ec50]}];
     j = 0;
     ];
    (*update ndrug*)
    j = j + 1; 
    ];

   For[i = Ndrug*everyH; j = everyH, i <= TIMEMAX, i = i + 1,
    conc = fn /. t -> j;
    AppendTo[concls, {i, conc, HoshenPDEff[conc, ec50]}
     ];
    j = j + 1;
    ];

   Developer`ToPackedArray[concls]
   
];


Ki[T_Real, concls_List] := Module[{i, lsk},
(*k_i (t) generate the list of decay constant (beta in the paper)   
	concls = the list of concentration and effects from	the function ConcMod
	T = 1/alpha 
*)
   lsk = Table[{(1/T)*Log[100./(100.-concls[[i, 3]])], 
      (1/T)*Log[100./(100.-concls[[i, 4]])], 
      (1/T)*Log[100./(100.- concls[[i, 5]])]}, {i, 1, Length[concls]}];
   lsk
   ];


(*decay function*)
(*ages = the ages of parasites*)
Fdecay[ages_, lst_, attime_, lsk_, stages_] := Module[{tmp, i},
   If[ages < 1 || ages > Length[lst], 
    Print["Please check your input!"]];
   i = ages;
   tmp = lst[[i]];
   If[stages[[i]] == 1, tmp = lst[[i]]*Exp[-lsk[[attime, 1]]]];
   If[stages[[i]] == 2, tmp = lst[[i]]*Exp[-lsk[[attime, 2]]]];
   If[stages[[i]] == 3, tmp = lst[[i]]*Exp[-lsk[[attime, 3]]]];
   tmp
   ];

(* decay function with dormancy *)
(*decay function with dormancy*)
FdecayDor[lst_, attime_, lsk_, stages_, dorfrac_] := 
  Module[{tmp, tmpdor, i},
   tmp = lst;
   tmp = Table[
     Which[
      stages[[i]] == 0, tmp[[i]],
      stages[[i]] == 1, tmp[[i]]*Exp[-lsk[[attime, 1]]],
      stages[[i]] == 2, tmp[[i]]*Exp[-lsk[[attime, 2]]], 
      stages[[i]] == 3, tmp[[i]]*Exp[-lsk[[attime, 3]]]
      ], {i, 1, Length@tmp}];
   
   tmpdor = tmp[[6 ;; 26]]*dorfrac;
   tmp[[6 ;; 26]] = tmp[[6 ;; 26]] - tmpdor;
   {tmp, tmpdor}
   ];


WhichRTS[lst_List, KZ_List] := Module[{i, tmp},
	(*WhichRTS defines the stages of the parasites from the kill zone and their ages.*)
	(*0=not in the kill zone 1=Rings 2=Trophozoites 3=Schizonts*)
	(*KZ={{BeginRing,EndRing},{BeginTrop,EndTrop},{BeginSchi,EndSchi}}*)
   tmp = Table[0, {Length[lst]}];(*initial values*)
   If[Length[KZ] != 3, Print["Please check your kill zone!\n"];];
   For[i = 1, i <= Length[lst], i = i + 1,
    If[IntervalMemberQ[Interval[KZ[[1]]], i] == True, 
     tmp = ReplacePart[tmp, i -> 1];];
    If[IntervalMemberQ[Interval[KZ[[2]]], i] == True, 
     tmp = ReplacePart[tmp, i -> 2];];
    If[IntervalMemberQ[Interval[KZ[[3]]], i] == True, 
     tmp = ReplacePart[tmp, i -> 3];];
    ];
   tmp
   ];


(*This is the main function for running the model.  Sexy5 reads the \
concentration data and generates concentration profile followings the \
dose regime. it uses CountRing with partition function*)
Sexy5[initN_, PMR_Integer, mu_, sigma_, hours_Integer, KZ_List, datafile_String, 
   everyH_Integer, Ndrug_Integer, gamma_List, ec50_List, emin_List, emax_List, 
   T_Real, runMax_Integer, outfile_Integer] := 
  Module[{runs, stages, concls, i, j, lst, lsk, time, output},
   
   (*input forms for KZ,gamma,ec50,emin and emax
   KZ={{rb,re},{tb,te},{sb,
   se}}  gamma={gammar,gammat,gammas}  ec50={r50,t50,s50}  
   emin={rmin,tmin,smin}  emax={rmax,tmax,smax}*)
   (**outfile = 1 write the output to modgen.csv ; 
   			  =	0 don't write output file.*)
   
   (*running status*)
   runs = True;
   
   (*time for calculating the concentrations (hr)..it is the 
   maximum time that the program can run!*)
   time = runMax;  
   
   (*list of the parasite distribution at each time step*)
   output = {}; 

   (*template of the drug and its effects*)
   concls = ConcMod[datafile, everyH, Ndrug, gamma, ec50, emin, emax];
   
   (*k_i (t)*)
   lsk = Ki[T, concls];
   
   (*initial parasite load*)
   lst = DistributeN[initN, hours, mu, sigma];
   
   (*tempplate for using k_i*)
   stages = WhichRTS[lst, KZ];
   
   (*(* sum all 1-11 + 11-48 * decay fn  *)
   tmp1 = Log10[CountRing[lst, {11, 14}, 1]] // Abs; 
   pmin = tmp1;*)
   
   (*adding the values*)
   AppendTo[output, lst];
   
   i = 0;
   While[runs == True && i < time,
    (*evolving the system*)
    i = i + 1;
    
    (*Parasites are growing. Feed them!*)
    lst = Shiftonehour[lst, PMR]; 
   
    (*a time to kill.*)
    For[j = 1, j <= Length[lst], j = j + 1,
     lst = ReplacePart[lst, j -> Fdecay[j, lst, i, lsk, stages]];
     ];
    
    (*adding a point for ploting*)
    Developer`ToPackedArray[AppendTo[output, lst]];
    
    ]; (* end While loop *)
    
   (* write the output to file*)
   If[outfile == 1,
    Export["modgen.csv", 
     Table[{i - 1, Log10@Total[output[[i]]]}, {i, 1, 
       Length[output]}], "CSV"];
    Print["The output has been written to modgen.csv.\n"];
    ,
    If[outfile == 2, 
     Developer`ToPackedArray[
      Table[{i - 1, Log10@Total[output[[i]]]}, {i, 1, 
        Length[output]}]]
     , If[outfile == 0, output]
     ]
    ]
       
   ];

Sexy52[initN_, PMR_Integer, mu_, sigma_, hours_Integer, KZ_List, datafile_String, 
   everyH_Integer, Ndrug_Integer, gamma_List, ec50_List, emin_List, emax_List, 
   T_Real, runMax_Integer, outfile_Integer] := 
  Module[{runs, stages, concls, i, j, lst, lsk, time, output},
   
   (*input forms for KZ,gamma,ec50,emin and emax
   KZ={{rb,re},{tb,te},{sb,
   se}}  gamma={gammar,gammat,gammas}  ec50={r50,t50,s50}  
   emin={rmin,tmin,smin}  emax={rmax,tmax,smax}*)
   (**outfile = 1 write the output to modgen.csv ; 
   			  =	0 don't write output file.*)
   
   (*running status*)
   runs = True;
   
   (*time for calculating the concentrations (hr)..it is the 
   maximum time that the program can run!*)
   time = runMax;  
   
   (*list of the parasite distribution at each time step*)
   output = {}; 

   (*template of the drug and its effects*)
   concls = ConcMod[datafile, everyH, Ndrug, gamma, ec50, emin, emax];
   
   (*k_i (t)*)
   lsk = Ki[T, concls];
   
   (*initial parasite load*)
   lst = DistributeN2[initN, hours, mu, sigma, PMR];
   
   (*tempplate for using k_i*)
   stages = WhichRTS[lst, KZ];
   
   (*(* sum all 1-11 + 11-48 * decay fn  *)
   tmp1 = Log10[CountRing[lst, {11, 14}, 1]] // Abs; 
   pmin = tmp1;*)
   
   (*adding the values*)
   AppendTo[output, lst];
   
   i = 0;
   While[runs == True && i < time,
    (*evolving the system*)
    i = i + 1;
    
    (*Parasites are growing. Feed them!*)
    lst = Shiftonehour[lst, PMR]; 
   
    (*a time to kill.*)
    For[j = 1, j <= Length[lst], j = j + 1,
     lst = ReplacePart[lst, j -> Fdecay[j, lst, i, lsk, stages]];
     ];
    
    (*adding a point for ploting*)
    Developer`ToPackedArray[AppendTo[output, lst]];
    
    ]; (* end While loop *)
    
   (* write the output to file*)
   If[outfile == 1,
    Export["modgen.csv", 
     Table[{i - 1, Log10@Total[output[[i]]]}, {i, 1, 
       Length[output]}], "CSV"];
    Print["The output has been written to modgen.csv.\n"];
    ,
    If[outfile == 2, 
     Developer`ToPackedArray[
      Table[{i - 1, Log10@Total[output[[i]]]}, {i, 1, 
        Length[output]}]]
     , If[outfile == 0, output]
     ]
    ]
       
   ];


(*This is the main function for running the model. Sexy5 reads the concentration 
data and generates concentration profile followings the dose regime. 
It uses CountRing with partition function*)
RunSexy5[initN_, PMR_Integer, mu_, sigma_, hours_Integer, KZ_List, datafile_String, 
  everyH_Integer, Ndrug_Integer, gamma_List, ec50_List, emin_List, emax_List, T_Real, runMax_Integer, 
  lim_Real] := Module[{runs, p1, p24, p48, tmp1, tmp2, stages, 
   pmin, prr24, prr48, concls, ct, rc, i, j, lst, lsk, cs, rcs, time, dec, inc, period, clcheck, output},
  
  (*
  datafile = name of the DHA concentration file
  input forms for KZ,gamma,ec50,emin and emax
  KZ={{rb,re},{tb,te},{sb, se}}  
  gamma={gammar,gammat,gammas}  
  ec50={r50,t50,s50}  
  emin={rmin,tmin,smin}  
  emax={rmax,tmax,smax}
  *)
  
  (* outfile = 1 write the output to modgen.csv ; 
  			 = 0 don't write output file. *) 
  runs = True;(*running status*)
  
  (*time for calculating the concentrations (hr)..it is 
  the maximum time that the program can run!*)
  time = runMax;
  
  (*period of blood sampling*)
  period = 6; 
  
  (*parasite clearance status: 0=not clear 1=clear*)
  cs = 0;
  
  (*0 = no recrudescence  1= recrudescence*)
  rcs = 0; 
  dec = {}; 
  inc = {}; 
  clcheck = {};
  
  (*initial clearance time*)
  ct = 6;
  
  (*list of the parasite distribution at each time step*)
  output = {}; 

  (*template of the drug and its effects*)
  concls = ConcMod[datafile, everyH, Ndrug, gamma, ec50, emin, emax];
  
  (*k_i (t)*)
  lsk = Ki[T, concls];
  
  (*initial parasite load*)
  lst = DistributeN[initN, hours, mu, sigma];
  
  (*tempplate for using k_i*)
  stages = WhichRTS[lst, KZ];
    
  (* sum all 1-11 + 11-48 * decay fn  *)
  tmp1 = Log10[CountRing[lst, {11, 14}, 1]] // Abs; 
  pmin = tmp1;
  
  (*adding the values*)
  AppendTo[output, lst];
  
  i = 0;
  While[runs == True && i < time,
   (*evolving the system*)
   i = i + 1;
   
   (*Parasites are growing. Feed them!*)
   lst = Shiftonehour[lst, PMR]; 
  
   (*a time to kill.*)
   For[j = 1, j <= Length[lst], j = j + 1,
    lst = ReplacePart[lst, j -> Fdecay[j, lst, i, lsk, stages]];
    ];
   
   AppendTo[output, lst];(*adding a point for ploting*)
  
   (* sum all 1-11 + 11-48 * decay fn *) 
   tmp2 = Log10[CountRing[lst, {11, 14}, 1]] // Abs; 
   
   (*finding minimum parasites*)
   If[tmp2 < pmin, pmin = tmp2;];
   
   (*parasite reduction ratio
   calculate from the ring parasites (using CountRing with decay function)
   *)
   If[i == 1, p1 = CountRing[lst, {11, 14}, 1];];
   If[i == 24, p24 = CountRing[lst, {11, 14}, 1];];
   If[i == 48, p48 = CountRing[lst, {11, 14}, 1];];
   prr24 = p1/p24; (*reduction ratio at 24*)
   prr48 = p1/p48; (*reduction ratio at 48*)
   
   (*finding the cross section point at the limit line*)
   If[tmp2 <=  lim && tmp1 > lim, AppendTo[dec, i - 1];];
   If[tmp2 > lim && tmp1 <=  lim, AppendTo[inc, i - 1];];
   
   (*update initial value for total parasites*)
   tmp1 = tmp2; 
   
   (*===clinical check===*)
   If[Mod[i, period] == 0,
    If[tmp2 > lim, 
    	AppendTo[clcheck, True];
    	, 
      	AppendTo[clcheck, False];
      ];
   ];
   
   (*clearance time*)
   If[Length[clcheck] >= 3,
    If[clcheck[[Length[clcheck]]] == False && 
       clcheck[[Length[clcheck] - 1]] == False && 
       clcheck[[Length[clcheck] - 2]] == True
       , 
       ct = (Length[clcheck] - 1)*period;
      ];
   ];
   
   (*recrudescence time*)
   If[Length[clcheck] >= 3,
    If[clcheck[[Length[clcheck]]] == True && 
       clcheck[[Length[clcheck] - 1]] == False && 
       clcheck[[Length[clcheck] - 2]] == False
       , 
       rc = Length[clcheck]*period; 
       rcs = 1; 
       runs = False;
      ];
   ];
   
   (*stop when the parasites have gone*)
   If[pmin < 0, 
   	runs = False; 
   	cs = 1; 
   	rcs = 0;
   	, 
   	cs = 0; 
   	rcs = 1;
   	];
   
   (*in the case it has no ct&rc....not below the detection limit*)
   If[Length[clcheck] >=  40 && clcheck[[Length[clcheck]]] == True && 
     clcheck[[Length[clcheck] - 1]] == True && 
     clcheck[[Length[clcheck] - 2]] == True, 
     cs = 0; 
     rc = "NO RC!"; 
     ct = "NO CT!"; 
     rcs = 0; 
     runs = False;
     ];
   
   ]; (* end While loop *)
  
  (**)
  If[ct == 6 && cs == 1,
   If[Length[Position[clcheck, False]] == 0, 
     ct = (Length[clcheck] + 1)*6;
     ];
   ];
  
  (*===output===*)
  Print["Total time: ", i];
  Print["intersection: ", {dec, inc}];
  If[i >= 24, Print["prr24: ", prr24];, prr24 = "NO PRR24";];
  If[i >= 48, Print["prr48: ", prr48];, prr48 = "NO PRR48";];
  Print["clearance time: ", ct];
  If[cs == 1 && rcs == 0, 
  		Print["Parasites have gone!"]; 
   		rc = "NO RC";
   	];
  If[rcs == 1 && cs == 0, 
  		Print["pmin(log10): ", pmin]; 
   		Print["recrudescence time: ", rc];
   	];
  
  {i, {dec, inc}, prr24, prr48, ct, pmin, rc}

];



(* comparing the parasite count data and model *)
ParaFit[parafile_, concfile_, initN_, PMR_, mu_, sigma_,
	LC_, KZ_, everyH_, Ndrug_, gamma_, ec50_, emin_, emax_, T_, output_] := 
  Module[{junk, tmp, modlist, parals, paralsln, rmsd, runMax, factor, i},
   (*output	=	0  return the RMSD;
   			= 	1 show the rmsd and graph
   *)
   parals = Import[parafile];
   paralsln = Table[{parals[[i, 1]], Log10@parals[[i, 3]] // N}, {i, 2, Length[parals]}];
   
   (**maximum time from the data**)
   runMax = paralsln[[All, 1]] // Max;
   Print["runmax: ",runMax];
   
   (** detection limit in log10 scale **)
   factor = Log10@parals[[1, 1]] // N;
   
   (* Sexy5[initN_, PMR_Integer, mu_, sigma_, hours_Integer, KZ_List, datafile_String, 
   everyH_Integer, Ndrug_Integer, gamma_List, ec50_List, emin_List, emax_List, 
   T_Real, runMax_Integer, outfile_Integer] *)
    junk = Sexy5[initN,PMR,mu,sigma,LC,KZ,concfile,everyH,Ndrug,gamma,ec50,emin,emax,T,runMax,0];
   
   (** use only the rings for comparison with the data **)
   junk = LsDot[#,Table[PRingFunc[i, 11, 14], {i, 1, LC}]]&/@junk;
   
   tmp = Table[{i,Log10[junk[[i + 1]]//Total]}, {i, 0, Length[junk] - 1}];
      
   junk =.;
   
      (**extract points*)
      modlist = MapDat[tmp, paralsln];
   (*modlist = MapDat[junk, datalist];*)
   
   (*
   If[Length[modlist] != Length[paralsln],
   (*Print["model: ",modlist];
   Print["data: ",paralsln];*)
    rmsd = "Bad input!!";];
   If[Length[modlist] == Length[paralsln],
       rmsd = RootMeanSquare[modlist[[All, 2]]-paralsln[[All, 2]]]];
   *)
   
   rmsd = RMSD[paralsln, modlist, factor];
       
   If[output == 1,
    (*Print["RMSD: ",rmsd]*)
    (*
    ListPlot[{modlist,paralsln},Epilog->Line[{{0,8},{80,8}}],
    PlotRange->{{0,Last[paralsln][[1]]},{0,13}},PlotStyle->{Red,Blue},
    Frame->True,Joined->{False,False},PlotMarkers->Automatic,
    Filling->{1->{2}},FrameLabel->{{"Subscript[log, 10]
     Number of parasites",""},{"time(hrs)",ToString[parafile]<>
    " rmsd:"<>ToString@rmsd}}]
    *)
    ListPlot[{modlist, paralsln}, Epilog -> Line[{{0, 8}, {runMax, 8}}], 
     PlotRange -> {{0, Last[paralsln][[1]]}, {0, 13}}, 
     PlotStyle -> {Red, Blue}, Frame -> True, Joined -> {True, False},
      FrameLabel -> {{"log10 Number of parasites", ""}, {"time(hrs)", 
        ToString[parafile] <> " rmsd:" <> ToString@rmsd}}]
     , 
    If[output == 0, rmsd]
    ]
   
   ];


(* comparing the parasite count data and model *)
ParaFit2[parafile_, concfile_, initN_, PMR_, mu_, sigma_,
	LC_, KZ_, everyH_, Ndrug_, gamma_, ec50_, emin_, emax_, T_, output_] := 
  Module[{junk, tmp, modlist, parals, paralsln, rmsd, runMax, factor, i},
   (*output	=	0  return the RMSD;
   			= 	1 show the rmsd and graph
   *)
   parals = Import[parafile];
   paralsln = Table[{parals[[i, 1]], Log10@parals[[i, 3]] // N}, {i, 2, Length[parals]}];
   
   (**maximum time from the data**)
   runMax = paralsln[[All, 1]] // Max;
   Print["runmax: ",runMax];
   
   (** detection limit in log10 scale **)
   factor = Log10@parals[[1, 1]] // N;
   
   (* Sexy5[initN_, PMR_Integer, mu_, sigma_, hours_Integer, KZ_List, datafile_String, 
   everyH_Integer, Ndrug_Integer, gamma_List, ec50_List, emin_List, emax_List, 
   T_Real, runMax_Integer, outfile_Integer] *)
    junk = Sexy52[initN,PMR,mu,sigma,LC,KZ,concfile,everyH,Ndrug,gamma,ec50,emin,emax,T,runMax,0];
   
   (** use only the rings for comparison with the data **)
   junk = LsDot[#,Table[PRingFunc[i, 11, 14], {i, 1, LC}]]&/@junk;
   
   tmp = Table[{i,Log10[junk[[i + 1]]//Total]}, {i, 0, Length[junk] - 1}];
      
   junk =.;
   
      (**extract points*)
      modlist = MapDat[tmp, paralsln];
   (*modlist = MapDat[junk, datalist];*)
   
   (*
   If[Length[modlist] != Length[paralsln],
   (*Print["model: ",modlist];
   Print["data: ",paralsln];*)
    rmsd = "Bad input!!";];
   If[Length[modlist] == Length[paralsln],
       rmsd = RootMeanSquare[modlist[[All, 2]]-paralsln[[All, 2]]]];
   *)
   
   rmsd = RMSD[paralsln, modlist, factor];
       
   If[output == 1,
    (*Print["RMSD: ",rmsd]*)
    (*
    ListPlot[{modlist,paralsln},Epilog->Line[{{0,8},{80,8}}],
    PlotRange->{{0,Last[paralsln][[1]]},{0,13}},PlotStyle->{Red,Blue},
    Frame->True,Joined->{False,False},PlotMarkers->Automatic,
    Filling->{1->{2}},FrameLabel->{{"Subscript[log, 10]
     Number of parasites",""},{"time(hrs)",ToString[parafile]<>
    " rmsd:"<>ToString@rmsd}}]
    *)
    ListPlot[{modlist, paralsln}, Epilog -> Line[{{0, 8}, {runMax, 8}}], 
     PlotRange -> {{0, Last[paralsln][[1]]}, {0, 13}}, 
     PlotStyle -> {Red, Blue}, Frame -> True, Joined -> {True, False},
      FrameLabel -> {{"log10 Number of parasites", ""}, {"time(hrs)", 
        ToString[parafile] <> " rmsd:" <> ToString@rmsd}}]
     , 
    If[output == 0, rmsd]
    ]
   
   ];



(* comparing the parasite count data and model *)
ParaFitHoshen[parafile_, concfile_, initN_, PMR_, mu_, sigma_,
	LC_, everyH_, Ndrug_, ec50_, dorfrac_, dortime_, output_] := 
  Module[{junk, modlist, parals, paralsln, rmsd, runMax, factor, i},
   (*output	=	0  return the RMSD;
   			= 	1 show the rmsd and graph
   *)
   parals = Import[parafile];
   paralsln = Table[{parals[[i, 1]], Log10@parals[[i, 3]] // N}, {i, 2, Length[parals]}];
   
   (**maximum time from the data**)
   runMax = paralsln[[All, 1]] // Max;
   Print["runmax: ",runMax];
   
   (** detection limit in log10 scale **)
   factor = Log10@parals[[1, 1]] // N;
   
   (* Sexy5[initN_, PMR_Integer, mu_, sigma_, hours_Integer, KZ_List, datafile_String, 
   everyH_Integer, Ndrug_Integer, gamma_List, ec50_List, emin_List, emax_List, 
   T_Real, runMax_Integer, outfile_Integer] *)
    junk = Hoshen[initN,PMR,mu,sigma,LC,concfile,everyH,Ndrug,ec50,runMax,dorfrac,dortime,2];
   
   (** use only the rings for comparison with the data **)
   (*junk = LsDot[#,Table[PRingFunc[i, 11, 14], {i, 1, LC}]]&/@junk;
   
   tmp = Table[{i,Log10[junk[[i + 1]]//Total]}, {i, 0, Length[junk] - 1}];
      
   junk =.;
   *)
      (**extract points*)
      modlist = MapDat[junk, paralsln];
   (*modlist = MapDat[junk, datalist];*)
   
   (*
   If[Length[modlist] != Length[paralsln],
   (*Print["model: ",modlist];
   Print["data: ",paralsln];*)
    rmsd = "Bad input!!";];
   If[Length[modlist] == Length[paralsln],
       rmsd = RootMeanSquare[modlist[[All, 2]]-paralsln[[All, 2]]]];
   *)
   
   rmsd = RMSD[paralsln, modlist, factor];
       
   If[output == 1,
    (*Print["RMSD: ",rmsd]*)
    (*
    ListPlot[{modlist,paralsln},Epilog->Line[{{0,8},{80,8}}],
    PlotRange->{{0,Last[paralsln][[1]]},{0,13}},PlotStyle->{Red,Blue},
    Frame->True,Joined->{False,False},PlotMarkers->Automatic,
    Filling->{1->{2}},FrameLabel->{{"Subscript[log, 10]
     Number of parasites",""},{"time(hrs)",ToString[parafile]<>
    " rmsd:"<>ToString@rmsd}}]
    *)
    ListPlot[{modlist, paralsln}, Epilog -> Line[{{0, 8}, {runMax, 8}}], 
     PlotRange -> {{0, Last[paralsln][[1]]}, {0, 13}}, 
     PlotStyle -> {Red, Blue}, Frame -> True, Joined -> {True, False},
      FrameLabel -> {{"log10 Number of parasites", ""}, {"time(hrs)", 
        ToString[parafile] <> " rmsd:" <> ToString@rmsd}}]
     , 
    If[output == 0, rmsd]
    ]
   
   ];

ParaFitHoshenNew[parafile_, concfile_, initN_, PMR_, mu_, sigma_,
	LC_, everyH_, Ndrug_, ec50_, dorfrac_, dortime_, output_] := 
  Module[{junk, modlist, parals, paralsln, rmsd, runMax, factor, i},
   (*output	=	0  return the RMSD;
   			= 	1 show the rmsd and graph
   *)
   parals = Import[parafile];
   paralsln = Table[{parals[[i, 1]], Log10@parals[[i, 3]] // N}, {i, 2, Length[parals]}];
   
   (**maximum time from the data**)
   runMax = paralsln[[All, 1]] // Max;
   Print["runmax: ",runMax];
   
   (** detection limit in log10 scale **)
   factor = Log10@parals[[1, 1]] // N;
   
   (* Sexy5[initN_, PMR_Integer, mu_, sigma_, hours_Integer, KZ_List, datafile_String, 
   everyH_Integer, Ndrug_Integer, gamma_List, ec50_List, emin_List, emax_List, 
   T_Real, runMax_Integer, outfile_Integer] *)
    junk = HoshenNew[initN,PMR,mu,sigma,LC,concfile,everyH,Ndrug,ec50,runMax,dorfrac,dortime,2];
   
   (** use only the rings for comparison with the data **)
   (*junk = LsDot[#,Table[PRingFunc[i, 11, 14], {i, 1, LC}]]&/@junk;
   
   tmp = Table[{i,Log10[junk[[i + 1]]//Total]}, {i, 0, Length[junk] - 1}];
      
   junk =.;
   *)
      (**extract points*)
      modlist = MapDat[junk, paralsln];
   (*modlist = MapDat[junk, datalist];*)
   
   (*
   If[Length[modlist] != Length[paralsln],
   (*Print["model: ",modlist];
   Print["data: ",paralsln];*)
    rmsd = "Bad input!!";];
   If[Length[modlist] == Length[paralsln],
       rmsd = RootMeanSquare[modlist[[All, 2]]-paralsln[[All, 2]]]];
   *)
   
   rmsd = RMSD[paralsln, modlist, factor];
       
   If[output == 1,
    (*Print["RMSD: ",rmsd]*)
    (*
    ListPlot[{modlist,paralsln},Epilog->Line[{{0,8},{80,8}}],
    PlotRange->{{0,Last[paralsln][[1]]},{0,13}},PlotStyle->{Red,Blue},
    Frame->True,Joined->{False,False},PlotMarkers->Automatic,
    Filling->{1->{2}},FrameLabel->{{"Subscript[log, 10]
     Number of parasites",""},{"time(hrs)",ToString[parafile]<>
    " rmsd:"<>ToString@rmsd}}]
    *)
    ListPlot[{modlist, paralsln}, Epilog -> Line[{{0, 8}, {runMax, 8}}], 
     PlotRange -> {{0, Last[paralsln][[1]]}, {0, 13}}, 
     PlotStyle -> {Red, Blue}, Frame -> True, Joined -> {True, False},
      FrameLabel -> {{"log10 Number of parasites", ""}, {"time(hrs)", 
        ToString[parafile] <> " rmsd:" <> ToString@rmsd}}]
     , 
    If[output == 0, rmsd]
    ]
   
   ];



(* compare dormancy model (sexy6->8) with the data *)
ParaFitD[parafile_, concfile_, initN_, PMR_, mu_, sigma_,
	LC_, KZ_, everyH_, Ndrug_, gamma_, ec50_, emin_, emax_, T_, dorfrac_, dortime_, output_] := 
  Module[{junk, tmp, modlist, parals, paralsln, rmsd, runMax, factor, i},
   (*output	=	0  return the RMSD;
   			= 	1 show the rmsd and graph
   *)
   parals = Import[parafile];
   paralsln = Table[{parals[[i, 1]], Log10@parals[[i, 3]] // N}, {i, 2, Length[parals]}];
   
   (**maximum time from the data**)
   runMax = paralsln[[All, 1]] // Max;
   (*Print["runmax: ",runMax];*)
   
   (** detection limit in log10 scale **)
   factor = Log10@parals[[1, 1]] // N;

(*Sexy6[initN_, PMR_, mu_, sigma_, hours_Integer, KZ_List, datafile_, 
   everyH_, Ndrug_, gamma_List, ec50_List, emin_List, emax_List, T_, runMax_, 
   dorfrac_, dortime_, outfile_]*)   
   
   (*junk = Sexy6[initN,PMR,mu,sigma,LC,KZ,concfile,everyH,Ndrug,gamma,ec50,emin,emax,T,runMax,dorfrac,dortime,0];*)
   
   junk = Sexy8[initN,PMR,mu,sigma,LC,KZ,concfile,everyH,Ndrug,gamma,ec50,emin,emax,T,runMax,dorfrac,dortime,0];
   
   
   (** use only the rings for comparison with the data **)
   junk = LsDot[#,Table[PRingFunc[i, 11, 14], {i, 1, LC}]]&/@junk;
   
   tmp = Table[{i,Log10[junk[[i + 1]]//Total]}, {i, 0, Length[junk] - 1}];
      
   junk =.;
   
      (**extract points*)
      modlist = MapDat[tmp, paralsln];
   (*modlist = MapDat[junk, datalist];*)
   
   (*
   If[Length[modlist] != Length[paralsln],
   (*Print["model: ",modlist];
   Print["data: ",paralsln];*)
    rmsd = "Bad input!!";];
   If[Length[modlist] == Length[paralsln],
       rmsd = RootMeanSquare[modlist[[All, 2]]-paralsln[[All, 2]]]];
   *)
   
   rmsd = RMSD[paralsln, modlist, factor];
       
   If[output == 1,
    (*Print["RMSD: ",rmsd]*)
    (*
    ListPlot[{modlist,paralsln},Epilog->Line[{{0,8},{80,8}}],
    PlotRange->{{0,Last[paralsln][[1]]},{0,13}},PlotStyle->{Red,Blue},
    Frame->True,Joined->{False,False},PlotMarkers->Automatic,
    Filling->{1->{2}},FrameLabel->{{"Subscript[log, 10]
     Number of parasites",""},{"time(hrs)",ToString[parafile]<>
    " rmsd:"<>ToString@rmsd}}]
    *)
    ListPlot[{modlist, paralsln}, Epilog -> Line[{{0, 8}, {runMax, 8}}], 
     PlotRange -> {{0, Last[paralsln][[1]]}, {0, 13}}, 
     PlotStyle -> {Red, Blue}, Frame -> True, Joined -> {True, False},
      FrameLabel -> {{"log10 Number of parasites", ""}, {"time(hrs)", 
        StringDrop[StringDrop[parafile,4],-4] <> " rmsd:" <> ToString@rmsd}}]
     , 
    If[output == 0, rmsd]
    ]
   
   ];


(**calculate the root mean square deviation of the data and model*)
RMSD[datls_List, modls_List, factor_Real] := 
  Module[{datlength, modlength, zeropos, rmsdout, dattemp, i},
   
   (*datls=list of {time , log10 (parasite)}*)
   
   rmsdout = 0;
   dattemp = datls;
   
   datlength = Length[dattemp];
   modlength = Length[modls];
   
   (*position of zero in data *)
   (** BUG !!  the data must be 0 not 0.0 **)
   (**
   zeropos=Position[dattemp[[All,2]],-Infinity]//Flatten;
   **)
   zeropos = Position[dattemp[[All, 2]], Indeterminate]~Join~ Position[dattemp[[All, 2]], -Infinity] // Flatten // Sort;
   
   (*If[Length[zeropos]!=0, Print["We found Zero!! ",zeropos];];*)
   
   (** if no zero **)
   If[Length[zeropos] == 0,
    If[datlength != modlength,
      (*not the same length*)
      (**rmsdout="NSL"; or 10^10 (large value)**) 
      rmsdout = 10^10;       
      ,
      rmsdout = RootMeanSquare[modls[[1;;Length[modls],2]]-dattemp[[1;;Length[dattemp],2]]];
      ];
    ];
   
   (** if the data has zero points **)
   If[Length[zeropos] != 0,
    If[datlength == modlength,
      For[i = 1, i <= Length[zeropos], i++,
       If[modls[[zeropos[[i]], 2]] >=  factor,       	
         	(*it will be (model (@zero)-factor)^2 when calculate rmsd. *)
         	dattemp[[zeropos[[i]], 2]] = factor;
         ,
         	(* it will be (model (@zero)-factor)^2 = 0 *)
         	dattemp[[zeropos[[i]], 2]] = modls[[zeropos[[i]], 2]];
         ];
       ];
      rmsdout = RootMeanSquare[modls[[1;;Length[modls],2]]-dattemp[[1;;Length[dattemp],2]]];
      ,
      (*not the same length and then return the large value (10^10).*)
      (*rmsdout="NSL";*)
      rmsdout = 10^10;
      
      ];
    ];
   
   rmsdout
];


MapDat[modlist_List, datalist_List] := Module[{i, j, outls},
(*make the model data has the same first column with the input data*)
   (**datalist = log10 real data**)
     outls = {};
     For[i = 1, i <= Length[datalist], i++,
    	For[j = 1, j <= Length[modlist], j++,
      		If[datalist[[i, 1]] == modlist[[j, 1]], 
       			AppendTo[outls, modlist[[j]]];
       		];
      	];
    ];
    outls
 ];

ExtractVal[txtin_List, keyword_String] := Module[{pos, val},
   (* extract the value of the keyword from the input file *)
   (* txtin = Import[inputfile,"Table"]; *)
   	pos = Select[Position[ToLowerCase /@ txtin, keyword], #[[2]] == 1 &];
   	If[Length@pos == 1, val = txtin[[pos[[1, 1]], 2]];];
	If[pos=={},
		Print["Please enter the value of ",keyword];
		Abort[];
		,
		Print[keyword,": ",val];
	];
   	
   ToExpression@val
];


(* const = T =1/alpha*)
RI[ec50_Real,emax_Real,gamma_Real,const_Real]:=Module[{ri},
	(*Print["Resistance Index = (Const*Ec50)/(Emax*Gamma)"];*)
	ri=const*ec50/(emax*gamma);
	ri
];


BatchRun[parafile_String, concfile_String, runsteps_Integer, outfilename_String] := 
  Module[{outfile, outstream, runid, initN, PMR, 
    Mu, Sigma, LC, rKZb, rKZe, tKZb, tKZe, sKZb, sKZe, everyH, Ndrug, 
    GammaR, GammaT, GammaS, ec50R, ec50T, ec50S, eminR, eminT, eminS, 
    emaxR, emaxT, emaxS, T, ec, rn, em, rmsd,riR,riT,riS},
   
   outfile = outfilename <> ".csv";
   (*parafile = "para" <> id <> ".csv";
   concfile = id <> ".csv";
   *)
   (***outstep=0 don't print run steps,1 print run steps***)
   
   (*Print["Data Id:",id ,"\tInput files: ",parafile,"\t",concfile];*)

   If[MemberQ[FileNames[], parafile] == False,
    	Print[parafile, " is not in the working folder."];
    	Abort[];
    ];
   If[MemberQ[FileNames[], concfile] == False,
    	Print[concfile, " is not in the working folder."];
    	Abort[];
    ];
   
   If[MemberQ[FileNames[], outfile] == False,
    	Print["Creating an output file.\n"];
    
    	(****write header of the output file****)
    	outstream = OpenWrite[outfile, PageWidth -> 350];
    
    	WriteString[outstream, "no,initN,PMR,Mu,Sigma,LC,rKZb,rKZe,tKZb,tKZe,sKZb,sKZe,everyH,Ndrug,GammaR,GammaT,GammaS,ec50R,ec50T,ec50S,eminR,eminT,eminS,emaxR,emaxT,emaxS,1/alpha,rmsd,riR,riT,riS\n"];
    
    	Close[outstream];
    	,
    	Print[outfile, " exists already."];
    	(*Print[outfile," has been copied to ","BAK-"<>outfile];
    	CopyFile[outfile,"BAK-"<>outfile];*)
    	Abort[];
    ];
   
   	For[runid = 1, runid <= runsteps, runid++, 
		Print["\n===Run no. " <> ToString[runid] <> "==="];
    
    	(**parasite load*)
    	initN = 10^RandomReal[NormalDistribution[Log10@(2.8*10^11), 0.5]] // Abs;
    
    	(**PMR**)
    	PMR = RandomInteger[{6, 12}];
    
    	(***Mu*)
    	Mu = RandomInteger[{4, 16}];
    
    	(**Sigma**)
    	Sigma = RandomInteger[{2, 8}];
    
    	(**life cycle**)
    	LC = 48;
    
    	(**kill zone*)
    	rKZb = 6;
    	rKZe = RandomInteger[{20, 30}];
    	tKZb = rKZe + 1;
    	tKZe = 38;
    	sKZb = tKZe + 1;
    	sKZe = 44;
    
    	ec = RandomChoice[{0, 1, 2, 3}];
    
    Which[
    	ec == 0, 
    	ec50R = RandomReal[{0, 100}];
		ec50T = ec50R;
     	ec50S = ec50T;
     	rn = RandomReal[{-5, 0}];
     	em = (1 - 10^rn)*100;
     	emaxR = em;
     	emaxT = emaxR;
     	emaxS = emaxT;
     	GammaR = RandomReal[{1.5, 9.5}];
     	GammaT = GammaR;
     	GammaS = GammaT;
     ,
     	ec == 1, ec50R = RandomReal[{0, 100}];
     	ec50T = RandomReal[{0, 100}];
     	ec50S = RandomReal[{0, 100}];
     	rn = RandomReal[{-5, 0}];
     	em = (1 - 10^rn)*100;
     	emaxR = em;
     	rn = RandomReal[{-5, 0}];
     	em = (1 - 10^rn)*100;
     	emaxT = em;
     	rn = RandomReal[{-5, 0}];
     	em = (1 - 10^rn)*100;
     	emaxS = em;
     	GammaR = RandomReal[{1.5, 9.5}];
     	GammaT = RandomReal[{1.5, 9.5}];
     	GammaS = RandomReal[{1.5, 9.5}];
     ,	
     	ec == 2, ec50R = RandomReal[{10, 100}];
     	ec50T = RandomReal[{0, 40}];
     	ec50S = ec50T;
     	rn = RandomReal[{-5, 0}];
     	em = (1 - 10^rn)*100;
     	emaxR = em;
     	rn = RandomReal[{-5, 0}];
     	em = (1 - 10^rn)*100;
     	emaxT = em;
     	emaxS = emaxT;
     	GammaR = RandomReal[{1.5, 9.5}];
     	GammaT = RandomReal[{1.5, 9.5}];
     	GammaS = GammaT;
    ,	
    	ec == 3, 
    	ec50R = RandomReal[{0, 100}];
     	ec50T = RandomReal[{0, 10}];
     	ec50S = ec50T;
     	rn = RandomReal[{-5, 0}];
     	em = (1 - 10^rn)*100;
     	emaxR = em;
     	emaxT = 99.99;
     	emaxS = emaxT;
     	GammaR = RandomReal[{1.5, 9.5}];
     	GammaT = RandomReal[{1.5, 9.5}];
     	GammaS = GammaT;
     ];
    	eminR = 0.0;
    	eminT = eminR;
    	eminS = eminT;
    	T = RandomReal[NormalDistribution[5, 2.5]] // Abs;
    	everyH = 24; Ndrug = 7;
      
    	rmsd = ParaFit[parafile, concfile, initN, PMR, Mu, Sigma, 
      		LC, {{rKZb, rKZe}, {tKZb, tKZe}, {sKZb, sKZe}}, everyH, 
      		Ndrug, {GammaR, GammaT, GammaS}, {ec50R, ec50T, ec50S}, {eminR, 
       		eminT, eminS}, {emaxR, emaxT, emaxS}, T, 0];
    	
    	riR=RI[ec50R,emaxR,GammaR,T]; (*  T*ec50R/(GammaR*emaxR)*)
    	riT=RI[ec50T,emaxT,GammaT,T]; (*T*ec50T/(GammaT*emaxT)*)
    	riS=RI[ec50S,emaxS,GammaS,T]; (*T*ec50S/(GammaS*emaxS)*)
    	
    Print["RMSD : ", rmsd];
    
    If[(rmsd != "NSL" || rmsd != 10^10) && rmsd <= 1.0,
    	outstream = OpenAppend[outfile, PageWidth -> 350];
    	WriteString[outstream, runid, ",", initN // FortranForm, ",", PMR, 
      		",", Mu, ",", Sigma, ",", LC, ",", rKZb, ",", rKZe, ",", tKZb, 
      		",", tKZe, ",", sKZb, ",", sKZe, ",", everyH, ",", Ndrug, ",", 
      		GammaR, ",", GammaT, ",", GammaS, ",", ec50R, ",", ec50T, ",", 
      		ec50S, ",", eminR, ",", eminT, ",", eminS, ",", emaxR, ",", 
     	 	emaxT, ",", emaxS, ",", T, ",", rmsd,",",riR,",",riT,",",riS,",", "\n"];
     
     Close[outstream];
     Print["The set of parameters is collected."];
     ,
     Print["The set of parameters is rejected."];
     ];
    ];
	Print["The program is finised."]
];



BatchRun2[parafile_String, concfile_String, runsteps_Integer, outfilename_String] := 
  Module[{outfile, outstream, runid, initN, PMR, 
    Mu, Sigma, LC, rKZb, rKZe, tKZb, tKZe, sKZb, sKZe, everyH, Ndrug, 
    GammaR, GammaT, GammaS, ec50R, ec50T, ec50S, eminR, eminT, eminS, 
    emaxR, emaxT, emaxS, T, ec, rn, em, rmsd,riR,riT,riS},
   
   outfile = outfilename <> ".csv";
   (*parafile = "para" <> id <> ".csv";
   concfile = id <> ".csv";
   *)
   (***outstep=0 don't print run steps,1 print run steps***)
   
   (*Print["Data Id:",id ,"\tInput files: ",parafile,"\t",concfile];*)

   If[MemberQ[FileNames[], parafile] == False,
    	Print[parafile, " is not in the working folder."];
    	Abort[];
    ];
   If[MemberQ[FileNames[], concfile] == False,
    	Print[concfile, " is not in the working folder."];
    	Abort[];
    ];
   
   If[MemberQ[FileNames[], outfile] == False,
    	Print["Creating an output file.\n"];
    
    	(****write header of the output file****)
    	outstream = OpenWrite[outfile, PageWidth -> 350];
    
    	WriteString[outstream, "no,initN,PMR,Mu,Sigma,LC,rKZb,rKZe,tKZb,tKZe,sKZb,sKZe,everyH,Ndrug,GammaR,GammaT,GammaS,ec50R,ec50T,ec50S,eminR,eminT,eminS,emaxR,emaxT,emaxS,1/alpha,rmsd,riR,riT,riS\n"];
    
    	Close[outstream];
    	,
    	Print[outfile, " exists already."];
    	(*Print[outfile," has been copied to ","BAK-"<>outfile];
    	CopyFile[outfile,"BAK-"<>outfile];*)
    	Abort[];
    ];
   
   	For[runid = 1, runid <= runsteps, runid++, 
		Print["\n===Run no. " <> ToString[runid] <> "==="];
    
    	(**parasite load*)
    	initN = 10^RandomReal[NormalDistribution[Log10@(2.8*10^11), 0.5]] // Abs;
    
    	(**PMR**)
    	PMR = RandomInteger[{6, 12}];
    
    	(***Mu*)
    	Mu = RandomInteger[{4, 16}];
    
    	(**Sigma**)
    	Sigma = RandomInteger[{2, 8}];
    
    	(**life cycle**)
    	LC = 48;
    
    	(**kill zone*)
    	rKZb = 6;
    	rKZe = RandomInteger[{20, 30}];
    	tKZb = rKZe + 1;
    	tKZe = 38;
    	sKZb = tKZe + 1;
    	sKZe = 44;
    
    	ec = RandomChoice[{0, 1, 2, 3}];
    
    Which[
    	ec == 0, 
    	ec50R = RandomReal[{0, 100}];
		ec50T = ec50R;
     	ec50S = ec50T;
     	rn = RandomReal[{-5, 0}];
     	em = (1 - 10^rn)*100;
     	emaxR = em;
     	emaxT = emaxR;
     	emaxS = emaxT;
     	GammaR = RandomReal[{1.5, 9.5}];
     	GammaT = GammaR;
     	GammaS = GammaT;
     ,
     	ec == 1, ec50R = RandomReal[{0, 100}];
     	ec50T = RandomReal[{0, 100}];
     	ec50S = RandomReal[{0, 100}];
     	rn = RandomReal[{-5, 0}];
     	em = (1 - 10^rn)*100;
     	emaxR = em;
     	rn = RandomReal[{-5, 0}];
     	em = (1 - 10^rn)*100;
     	emaxT = em;
     	rn = RandomReal[{-5, 0}];
     	em = (1 - 10^rn)*100;
     	emaxS = em;
     	GammaR = RandomReal[{1.5, 9.5}];
     	GammaT = RandomReal[{1.5, 9.5}];
     	GammaS = RandomReal[{1.5, 9.5}];
     ,	
     	ec == 2, ec50R = RandomReal[{10, 100}];
     	ec50T = RandomReal[{0, 40}];
     	ec50S = ec50T;
     	rn = RandomReal[{-5, 0}];
     	em = (1 - 10^rn)*100;
     	emaxR = em;
     	rn = RandomReal[{-5, 0}];
     	em = (1 - 10^rn)*100;
     	emaxT = em;
     	emaxS = emaxT;
     	GammaR = RandomReal[{1.5, 9.5}];
     	GammaT = RandomReal[{1.5, 9.5}];
     	GammaS = GammaT;
    ,	
    	ec == 3, 
    	ec50R = RandomReal[{0, 100}];
     	ec50T = RandomReal[{0, 10}];
     	ec50S = ec50T;
     	rn = RandomReal[{-5, 0}];
     	em = (1 - 10^rn)*100;
     	emaxR = em;
     	emaxT = 99.99;
     	emaxS = emaxT;
     	GammaR = RandomReal[{1.5, 9.5}];
     	GammaT = RandomReal[{1.5, 9.5}];
     	GammaS = GammaT;
     ];
    	eminR = 0.0;
    	eminT = eminR;
    	eminS = eminT;
    	T = RandomReal[NormalDistribution[5, 2.5]] // Abs;
    	everyH = 24; Ndrug = 7;
      
    	rmsd = ParaFit2[parafile, concfile, initN, PMR, Mu, Sigma, 
      		LC, {{rKZb, rKZe}, {tKZb, tKZe}, {sKZb, sKZe}}, everyH, 
      		Ndrug, {GammaR, GammaT, GammaS}, {ec50R, ec50T, ec50S}, {eminR, 
       		eminT, eminS}, {emaxR, emaxT, emaxS}, T, 0];
    	
    	riR=RI[ec50R,emaxR,GammaR,T]; (*  T*ec50R/(GammaR*emaxR)*)
    	riT=RI[ec50T,emaxT,GammaT,T]; (*T*ec50T/(GammaT*emaxT)*)
    	riS=RI[ec50S,emaxS,GammaS,T]; (*T*ec50S/(GammaS*emaxS)*)
    	
    Print["RMSD : ", rmsd];
    
    If[(rmsd != "NSL" || rmsd != 10^10) && rmsd <= 1.0,
    	outstream = OpenAppend[outfile, PageWidth -> 350];
    	WriteString[outstream, runid, ",", initN // FortranForm, ",", PMR, 
      		",", Mu, ",", Sigma, ",", LC, ",", rKZb, ",", rKZe, ",", tKZb, 
      		",", tKZe, ",", sKZb, ",", sKZe, ",", everyH, ",", Ndrug, ",", 
      		GammaR, ",", GammaT, ",", GammaS, ",", ec50R, ",", ec50T, ",", 
      		ec50S, ",", eminR, ",", eminT, ",", eminS, ",", emaxR, ",", 
     	 	emaxT, ",", emaxS, ",", T, ",", rmsd,",",riR,",",riT,",",riS,",", "\n"];
     
     Close[outstream];
     Print["The set of parameters is collected."];
     ,
     Print["The set of parameters is rejected."];
     ];
    ];
	Print["The program is finised."]
];



BatchRunHoshen[parafile_String, concfile_String, runsteps_Integer, outfilename_String] := 
  Module[{outfile, outstream, runid, initN, PMR, 
    Mu, Sigma, LC, everyH, Ndrug, ec50, dorfrac, dortime, rmsd},
   
   (*ParaFitHoshen[parafile_, concfile_, initN_, PMR_, mu_, sigma_,
	LC_, everyH_, Ndrug_, ec50_, dorfrac_, dortime_, output_]*)
   
   outfile = outfilename <> ".csv";
   (*parafile = "para" <> id <> ".csv";
   concfile = id <> ".csv";
   *)
   (***outstep=0 don't print run steps,1 print run steps***)
   
   (*Print["Data Id:",id ,"\tInput files: ",parafile,"\t",concfile];*)

   If[MemberQ[FileNames[], parafile] == False,
    	Print[parafile, " is not in the working folder."];
    	Abort[];
    ];
   If[MemberQ[FileNames[], concfile] == False,
    	Print[concfile, " is not in the working folder."];
    	Abort[];
    ];
   
   If[MemberQ[FileNames[], outfile] == False,
    	Print["Creating an output file.\n"];
    
    	(****write header of the output file****)
    	outstream = OpenWrite[outfile, PageWidth -> 350];
    
    	WriteString[outstream, "no,initN,PMR,Mu,Sigma,LC,everyH,Ndrug,ec50,dorfrac,dortime,rmsd\n"];
    
    	Close[outstream];
    	,
    	Print[outfile, " exists already."];
    	(*Print[outfile," has been copied to ","BAK-"<>outfile];
    	CopyFile[outfile,"BAK-"<>outfile];*)
    	Abort[];
    ];
   
   	For[runid = 1, runid <= runsteps, runid++, 
		Print["\n===Run no. " <> ToString[runid] <> "==="];
    
    	(**parasite load*)
    	initN = 10^RandomReal[NormalDistribution[Log10@(2.8*10^11), 0.5]] // Abs;
    
    	(**PMR**)
    	PMR = RandomInteger[{6, 12}];
    
    	(***Mu*)
    	Mu = RandomInteger[{4, 16}];
    
    	(**Sigma**)
    	Sigma = RandomInteger[{2, 8}];
    
    	(**life cycle**)
    	LC = 48;
        
    	ec50 = RandomReal[{1., 100.}];
 
    	everyH = 24; Ndrug = 7;
      	
      	dorfrac=RandomReal[0.001];
      	dortime=RandomInteger[{1,480}];
      	(*ParaFitHoshen[parafile_, concfile_, initN_, PMR_, mu_, sigma_,
	LC_, everyH_, Ndrug_, ec50_, dorfrac_, dortime_, output_] *)
	
    	rmsd = ParaFitHoshen[parafile, concfile, initN, PMR, Mu, Sigma, 
      		LC, everyH, Ndrug,ec50,dorfrac,dortime, 0];
    	    	
    Print["RMSD : ", rmsd];
    
    If[(rmsd != "NSL" || rmsd != 10^10) && rmsd <= 1.0,
    	outstream = OpenAppend[outfile, PageWidth -> 350];
    	WriteString[outstream, runid, ",", initN // FortranForm, ",", PMR, 
      		",", Mu, ",", Sigma, ",", LC, ",", everyH, ",", Ndrug, ",", 
      		ec50, ",", dorfrac//FortranForm,",",dortime,",",rmsd,"\n"];
     
     Close[outstream];
     Print["The set of parameters is collected."];
     ,
     Print["The set of parameters is rejected."];
     ];
    ];
	Print["The program is finised."]
];


BatchRunHoshenNew[parafile_String, concfile_String, runsteps_Integer, outfilename_String] := 
  Module[{outfile, outstream, runid, initN, PMR, 
    Mu, Sigma, LC, everyH, Ndrug, ec50, dorfrac, dortime, rmsd},
   
   (*ParaFitHoshen[parafile_, concfile_, initN_, PMR_, mu_, sigma_,
	LC_, everyH_, Ndrug_, ec50_, dorfrac_, dortime_, output_]*)
   
   outfile = outfilename <> ".csv";
   (*parafile = "para" <> id <> ".csv";
   concfile = id <> ".csv";
   *)
   (***outstep=0 don't print run steps,1 print run steps***)
   
   (*Print["Data Id:",id ,"\tInput files: ",parafile,"\t",concfile];*)

   If[MemberQ[FileNames[], parafile] == False,
    	Print[parafile, " is not in the working folder."];
    	Abort[];
    ];
   If[MemberQ[FileNames[], concfile] == False,
    	Print[concfile, " is not in the working folder."];
    	Abort[];
    ];
   
   If[MemberQ[FileNames[], outfile] == False,
    	Print["Creating an output file.\n"];
    
    	(****write header of the output file****)
    	outstream = OpenWrite[outfile, PageWidth -> 350];
    
    	WriteString[outstream, "no,initN,PMR,Mu,Sigma,LC,everyH,Ndrug,ec50,dorfrac,dortime,rmsd\n"];
    
    	Close[outstream];
    	,
    	Print[outfile, " exists already."];
    	(*Print[outfile," has been copied to ","BAK-"<>outfile];
    	CopyFile[outfile,"BAK-"<>outfile];*)
    	Abort[];
    ];
   
   	For[runid = 1, runid <= runsteps, runid++, 
		Print["\n===Run no. " <> ToString[runid] <> "==="];
    
    	(**parasite load*)
    	initN = 10^RandomReal[NormalDistribution[Log10@(2.8*10^11), 0.5]] // Abs;
    
    	(**PMR**)
    	PMR = RandomInteger[{6, 12}];
    
    	(***Mu*)
    	Mu = RandomInteger[{4, 16}];
    
    	(**Sigma**)
    	Sigma = RandomInteger[{2, 8}];
    
    	(**life cycle**)
    	LC = 48;
        
    	ec50 = RandomReal[{1., 100.}];
 
    	everyH = 24; Ndrug = 7;
      	
      	dorfrac=RandomReal[0.001];
      	dortime=RandomInteger[{1,120}];
      	(*ParaFitHoshen[parafile_, concfile_, initN_, PMR_, mu_, sigma_,
	LC_, everyH_, Ndrug_, ec50_, dorfrac_, dortime_, output_] *)
	
    	rmsd = ParaFitHoshenNew[parafile, concfile, initN, PMR, Mu, Sigma, 
      		LC, everyH, Ndrug,ec50,dorfrac,dortime, 0];
    	    	
    Print["RMSD : ", rmsd];
    
    If[(rmsd != "NSL" || rmsd != 10^10) && rmsd <= 1.0,
    	outstream = OpenAppend[outfile, PageWidth -> 350];
    	WriteString[outstream, runid, ",", initN // FortranForm, ",", PMR, 
      		",", Mu, ",", Sigma, ",", LC, ",", everyH, ",", Ndrug, ",", 
      		ec50, ",", dorfrac//FortranForm,",",dortime,",",rmsd,"\n"];
     
     Close[outstream];
     Print["The set of parameters is collected."];
     ,
     Print["The set of parameters is rejected."];
     ];
    ];
	Print["The program is finised."]
];


(* model 2.7 ==>  do the batch run of dormancy model (model 1.9 & 2.6) *)
BatchRunD[parafile_String, concfile_String, runsteps_Integer, outfilename_String] := 
  Module[{outfile, outstream, runid, initN, PMR, 
    Mu, Sigma, LC, rKZb, rKZe, tKZb, tKZe, sKZb, sKZe, everyH, Ndrug, 
    GammaR, GammaT, GammaS, ec50R, ec50T, ec50S, eminR, eminT, eminS, 
    emaxR, emaxT, emaxS, T, (*ec, em*) rn, rmsd,riR,riT,riS,dorfrac, dortime},
   
   outfile = outfilename <> ".csv";
   (*parafile = "para" <> id <> ".csv";
   concfile = id <> ".csv";
   *)
   (***outstep=0 don't print run steps,1 print run steps***)
   
   (*Print["Data Id:",id ,"\tInput files: ",parafile,"\t",concfile];*)

   If[MemberQ[FileNames[], parafile] == False,
    	Print[parafile, " is not in the working folder."];
    	Abort[];
    ];
   If[MemberQ[FileNames[], concfile] == False,
    	Print[concfile, " is not in the working folder."];
    	Abort[];
    ];
   
   If[MemberQ[FileNames[], outfile] == False,
    	Print["Creating an output file.\n"];
    
    	(****write header of the output file****)
    	outstream = OpenWrite[outfile, PageWidth -> 350];
    
    	WriteString[outstream, "no,initN,PMR,Mu,Sigma,LC,rKZb,rKZe,tKZb,tKZe,sKZb,sKZe,everyH,Ndrug,GammaR,GammaT,GammaS,ec50R,ec50T,ec50S,eminR,eminT,eminS,emaxR,emaxT,emaxS,1/alpha,rmsd,riR,riT,riS,dorfrac,dortime\n"];
    
    	Close[outstream];
    	,
    	Print[outfile, " exists already."];
    	(*Print[outfile," has been copied to ","BAK-"<>outfile];
    	CopyFile[outfile,"BAK-"<>outfile];*)
    	Abort[];
    ];
   
   	For[runid = 1, runid <= runsteps, runid++, 
		Print["\n===Run no. " <> ToString[runid] <> "==="];
    
    	(**parasite load*)
    	initN = 10^RandomReal[NormalDistribution[Log10@(2.8*10^11), 0.5]] // Abs;
    
    	(**PMR**)
    	PMR = RandomInteger[{6, 12}];
    
    	(***Mu*)
    	Mu = RandomInteger[{4, 16}];
    
    	(**Sigma**)
    	Sigma = RandomInteger[{2, 8}];
    
    	(**life cycle**)
    	LC = 48;
    
    	(**kill zone*)
    	rKZb = 6;
    	(*rKZe = RandomInteger[{20, 30}];*)
    	rKZe= 26;
    	(*tKZb = rKZe + 1;*)
    	tKZb = 27;
    	tKZe = 38;
    	(*sKZb = tKZe + 1;*)
    	sKZb = 39;
    	sKZe = 44;


		GammaR = RandomReal[{1.5, 9.5}];
     	GammaT = RandomReal[{1.5, 9.5}];
     	GammaS = RandomReal[{1.5, 9.5}];
     	
     	ec50R = RandomReal[{0, 100}];
     	ec50T = RandomReal[{0, 100}];
		ec50S = RandomReal[{0, 100}];
		
		eminR=0.0; eminT=0.0; eminS=0.0;
		rn = RandomReal[{-5, 0},3];
     	{emaxR,emaxT,emaxS} = (1 - 10^rn)*100;
		
		T = RandomReal[NormalDistribution[5, 2.5]] // Abs;
		

    	everyH = 24; Ndrug = 7;
        
        (* about 1 %  *) 
      	dorfrac=RandomReal[0.01];
      	(* about a week 168 hours *)
      	dortime=RandomInteger[{1,168}];
      	
    	rmsd = ParaFitD[parafile, concfile, initN, PMR, Mu, Sigma, 
      		LC, {{rKZb, rKZe}, {tKZb, tKZe}, {sKZb, sKZe}}, everyH, 
      		Ndrug, {GammaR, GammaT, GammaS}, {ec50R, ec50T, ec50S}, {eminR, 
       		eminT, eminS}, {emaxR, emaxT, emaxS}, T, dorfrac, dortime, 0];
    	
    	riR=RI[ec50R,emaxR,GammaR,T]; (*  T*ec50R/(GammaR*emaxR)*)
    	riT=RI[ec50T,emaxT,GammaT,T]; (*T*ec50T/(GammaT*emaxT)*)
    	riS=RI[ec50S,emaxS,GammaS,T]; (*T*ec50S/(GammaS*emaxS)*)
    	
    Print["RMSD : ", rmsd];
    
    If[(rmsd != "NSL" || rmsd != 10^10) && rmsd <= 1.0,
    	outstream = OpenAppend[outfile, PageWidth -> 350];
    	WriteString[outstream, runid, ",", initN // FortranForm, ",", PMR, 
      		",", Mu, ",", Sigma, ",", LC, ",", rKZb, ",", rKZe, ",", tKZb, 
      		",", tKZe, ",", sKZb, ",", sKZe, ",", everyH, ",", Ndrug, ",", 
      		GammaR, ",", GammaT, ",", GammaS, ",", ec50R, ",", ec50T, ",", 
      		ec50S, ",", eminR, ",", eminT, ",", eminS, ",", emaxR, ",", 
     	 	emaxT, ",", emaxS, ",", T, ",", rmsd,",",riR,",",riT,",",riS,",",dorfrac//FortranForm,",",dortime, "\n"];
     
     Close[outstream];
     Print["The set of parameters is collected."];
     ,
     Print["The set of parameters is rejected."];
     ];
    ];
	Print["The program is finised."]
];

(* randall-27052013 *)
WPDAT = {{8.17159, 5.36593, 2.40103, 13.2581, 36.7942, 80.5593, 40.8794, 
  99.6815, 99.96, 1.44036}, {4.95536, 8.34241, 4.36846, 71.3956, 
  61.2569, 84.4491, 92.7302, 99.998, 99.9552, 1.59431}, {3.5729, 
  2.56308, 3.02016, 96.1771, 9.35341, 27.1718, 99.9972, 99.9988, 
  99.8338, 3.42768}, {6.25071, 8.17665, 4.29768, 75.4635, 84.138, 
  74.5975, 99.9958, 99.9875, 73.1494, 3.69538}, {3.82955, 8.92304, 
  7.99278, 23.5114, 1.83057, 9.94115, 7.12565, 99.9419, 99.9521, 
  5.37209}, {8.74277, 6.43023, 6.07184, 86.0671, 30.5329, 50.3848, 
  99.7973, 99.9687, 99.3983, 2.07667}, {6.13732, 3.43439, 8.02635, 
  57.7876, 84.4255, 96.7605, 99.9586, 99.9954, 99.9337, 
  2.49827}, {2.83287, 3.52711, 7.0418, 87.0532, 54.0949, 59.9887, 
  99.9941, 99.8069, 95.8823, 2.73589}, {2.78463, 9.17345, 2.68464, 
  65.3424, 38.4376, 14.1801, 99.9895, 99.9987, 99.1744, 
  4.2841}, {8.06267, 9.31438, 5.16169, 53.7065, 1.06551, 66.0359, 
  99.9216, 99.3586, 99.0276, 2.66782}, {5.33608, 9.42678, 5.9529, 
  59.3989, 57.0834, 56.9265, 99.7972, 99.9599, 17.2647, 
  3.35787}, {8.44489, 5.796, 4.74902, 48.212, 77.8883, 66.446, 
  97.0709, 99.9915, 99.9434, 4.03325}, {4.17846, 4.06008, 2.63439, 
  77.971, 98.5652, 60.0894, 99.9623, 99.9975, 74.0285, 
  1.81848}, {4.13513, 3.61146, 4.15068, 53.048, 5.36826, 60.261, 
  99.0128, 99.9974, 96.715, 2.6736}, {4.63587, 4.50722, 6.86977, 
  96.674, 57.9288, 28.1495, 99.9501, 99.1498, 99.7076, 
  0.640381}, {9.17426, 2.8915, 8.53325, 50.1065, 99.1786, 83.7381, 
  99.9208, 99.997, 99.9057, 6.24037}, {5.18681, 3.10254, 3.31382, 
  62.4808, 39.5133, 40.7548, 99.4355, 99.9868, 33.1732, 
  2.29856}, {6.39937, 5.14185, 3.64495, 33.0319, 65.4909, 51.1381, 
  99.9527, 26.2934, 99.8867, 6.26167}, {3.65218, 5.80987, 7.8867, 
  56.8687, 37.7537, 49.288, 99.9089, 99.984, 99.7716, 
  1.45185}, {7.27038, 7.92816, 8.00021, 31.8276, 1.20662, 66.1422, 
  98.8886, 99.9886, 99.9984, 3.21963}, {5.1327, 3.05238, 5.39086, 
  74.1372, 4.83492, 91.5842, 71.5181, 98.7166, 97.6114, 
  1.61217}, {7.67913, 4.83743, 8.41432, 37.7506, 49.9714, 47.0959, 
  58.4835, 99.9848, 98.2689, 3.41445}, {8.32902, 7.20122, 3.79352, 
  0.468523, 1.82454, 80.5256, 63.5577, 99.9944, 99.9989, 
  6.09665}, {5.11015, 2.26752, 6.35988, 22.4655, 3.16507, 37.622, 
  99.6686, 99.9988, 99.9055, 5.25434}, {5.0826, 4.66639, 3.75638, 
  48.7149, 70.9417, 40.1057, 90.3881, 99.9975, 99.9969, 
  2.79889}, {4.49386, 6.46556, 4.57484, 64.2482, 70.7738, 89.3048, 
  99.7916, 99.9895, 99.985, 3.98762}, {4.84466, 4.94479, 3.17098, 
  18.2377, 25.172, 1.34331, 93.1352, 99.3342, 99.9899, 
  5.31699}, {2.0915, 6.06529, 2.3645, 61.6902, 11.2562, 17.8373, 
  99.9566, 99.8703, 99.9787, 4.89333}, {2.02792, 6.35801, 2.14702, 
  26.7839, 5.75638, 30.4114, 99.3265, 99.9988, 99.997, 
  3.0153}, {3.95097, 1.5547, 5.13588, 66.4936, 1.84877, 16.624, 
  99.8632, 97.2839, 95.6912, 3.6343}, {3.40614, 5.49948, 2.09142, 
  29.71, 34.5661, 18.587, 92.312, 99.9231, 99.996, 2.72057}, {7.93797,
   2.13583, 4.22515, 78.7755, 78.0107, 88.3944, 22.3061, 93.2726, 
  93.3089, 0.439619}, {5.70515, 8.82963, 5.24869, 86.1309, 61.5892, 
  3.54161, 99.9758, 96.9036, 99.7888, 2.00155}, {8.55301, 6.06864, 
  3.9556, 40.9389, 48.884, 36.8049, 99.9621, 99.9166, 97.0352, 
  3.95529}, {3.08177, 6.45961, 6.09323, 22.2209, 40.3739, 92.768, 
  99.9441, 99.994, 98.9784, 5.37458}, {2.95495, 7.55884, 5.06722, 
  4.90849, 47.9136, 16.0792, 95.7009, 99.9872, 99.872, 
  3.27054}, {7.56739, 1.89577, 5.51855, 94.332, 3.04507, 40.8891, 
  99.9509, 99.9652, 99.7864, 5.41104}, {8.73375, 4.16594, 7.69404, 
  80.6439, 66.8618, 12.9582, 99.9498, 99.0696, 99.951, 
  2.2517}, {7.15995, 6.07758, 8.88232, 24.3732, 66.557, 3.48642, 
  99.9959, 99.9985, 98.5624, 4.27949}, {1.59385, 6.4461, 2.5331, 
  16.7328, 80.3903, 18.4052, 99.2645, 99.9534, 99.8963, 
  3.12186}, {5.73706, 9.13062, 7.56655, 98.8064, 15.5981, 70.7789, 
  99.9111, 99.9103, 99.3152, 3.40204}, {5.83923, 8.90475, 4.99659, 
  80.3864, 41.7393, 64.0852, 99.9989, 99.9515, 99.9839, 
  3.24704}, {3.40614, 5.49948, 2.09142, 29.71, 34.5661, 18.587, 
  92.312, 99.9231, 99.996, 2.72057}, {7.56739, 1.89577, 5.51855, 
  94.332, 3.04507, 40.8891, 99.9509, 99.9652, 99.7864, 
  5.41104}, {8.56064, 4.0798, 6.23229, 68.4971, 3.07562, 81.7657, 
  39.0724, 99.9104, 99.933, 1.5255}, {8.31968, 4.87243, 5.45809, 
  90.1482, 45.6475, 0.00888459, 99.9975, 91.7872, 99.9721, 
  3.79542}, {8.62941, 7.48201, 5.64701, 85.9279, 31.1717, 33.9174, 
  99.5665, 99.9806, 99.9843, 5.3}, {2.53688, 2.11991, 4.06095, 
  45.3182, 8.45816, 77.063, 99.9919, 99.9976, 99.9892, 
  3.07049}, {2.70571, 7.5683, 8.36938, 53.7243, 11.4417, 41.5471, 
  99.9032, 97.54, 99.9929, 2.90316}, {7.93797, 2.13583, 4.22515, 
  78.7755, 78.0107, 88.3944, 22.3061, 93.2726, 93.3089, 
  0.439619}, {2.55174, 2.7422, 7.52042, 15.258, 43.5598, 54.6368, 
  89.779, 99.9706, 99.9983, 3.61505}, {8.61137, 6.21647, 6.8243, 
  45.2802, 4.76797, 11.2627, 80.3174, 85.8542, 98.6798, 
  2.88269}, {7.5004, 3.29065, 4.9553, 57.5345, 25.8202, 40.4459, 
  99.821, 99.5656, 99.9965, 4.6851}, {5.04987, 4.95548, 6.65274, 
  97.1297, 99.6683, 54.5529, 92.0771, 99.998, 99.9979, 
  1.71289}, {8.02895, 6.54929, 6.11998, 71.8937, 26.6532, 39.7071, 
  99.9931, 99.998, 99.9886, 6.00061}, {4.59092, 2.64302, 7.4616, 
  39.0572, 36.6232, 38.5731, 99.9569, 99.9441, 99.9976, 
  6.05951}, {8.79266, 4.36402, 1.52192, 79.3538, 23.3831, 15.8973, 
  99.2874, 99.3115, 99.9561, 3.54502}, {5.31975, 9.49133, 8.9164, 
  11.4241, 52.4424, 68.735, 97.3428, 99.9989, 99.9703, 
  5.82214}, {1.55774, 7.29011, 6.98743, 69.2623, 73.4909, 87.3938, 
  99.9521, 99.987, 99.7647, 1.37633}, {8.75326, 7.23149, 5.94295, 
  45.5134, 1.21615, 34.6691, 96.0117, 99.9533, 99.9816, 
  4.54101}, {2.10286, 7.4609, 7.87806, 24.3854, 3.75073, 5.56342, 
  15.2494, 90.8364, 99.807, 4.38845}, {9.16186, 5.86724, 1.90294, 
  19.67, 1.60062, 43.1174, 53.4679, 99.9982, 99.8227, 
  5.20774}, {3.52058, 1.583, 4.89188, 86.8855, 3.3733, 43.7047, 
  97.4038, 99.8581, 99.9969, 5.75182}, {5.27486, 7.50424, 6.70286, 
  25.9414, 13.214, 9.77231, 38.8266, 99.5509, 99.9974, 
  3.5927}, {2.10815, 7.79961, 1.67748, 95.9313, 75.6882, 69.0133, 
  73.9512, 99.9341, 99.955, 3.489}, {4.25224, 9.29299, 8.03188, 
  40.3278, 0.316924, 43.2878, 33.6809, 99.9924, 98.4207, 
  5.13252}, {7.32832, 6.76179, 5.24326, 55.2096, 20.3594, 67.1554, 
  89.6494, 99.9925, 99.9848, 8.72502}, {6.23698, 4.52677, 9.1236, 
  69.4198, 0.514033, 26.6689, 95.4592, 99.9961, 99.9965, 
  11.4131}, {3.22243, 1.81218, 7.61006, 88.1414, 2.4508, 73.4654, 
  10.7209, 99.9971, 99.9849, 4.3851}, {8.17508, 2.34137, 8.07461, 
  41.4143, 23.45, 9.19868, 43.1932, 97.8863, 99.9989, 
  3.55399}, {4.18296, 4.01782, 8.92793, 96.085, 31.6022, 84.0935, 
  73.4458, 99.9977, 99.8578, 0.501713}, {6.29123, 9.41658, 5.44269, 
  57.8927, 9.8449, 6.32406, 99.5925, 99.9958, 99.8966, 
  5.9371}, {6.33912, 7.19431, 5.96699, 90.517, 0.265973, 85.4556, 
  99.9295, 99.7114, 99.6267, 3.12036}, {1.54183, 4.05336, 9.30767, 
  43.0725, 23.8976, 43.5793, 96.6807, 99.3555, 99.9812, 
  1.2868}, {4.44636, 7.53874, 4.64896, 80.517, 47.5257, 46.9692, 
  99.6971, 99.7545, 98.9655, 1.49608}, {3.85371, 6.10772, 4.61011, 
  80.4676, 1.10796, 31.2201, 98.3393, 99.3778, 98.1158, 
  1.29833}, {2.40203, 5.23382, 5.94657, 86.684, 5.49382, 83.7561, 
  8.26938, 99.9659, 99.3921, 0.480626}, {8.10511, 7.41246, 1.84512, 
  64.5034, 14.0499, 12.0411, 99.9912, 99.9205, 99.9838, 
  3.50452}, {4.21313, 6.37235, 4.02966, 74.7741, 40.9294, 54.9772, 
  96.4724, 99.8585, 92.8554, 0.239513}, {6.3682, 8.48475, 9.09448, 
  14.9998, 11.2494, 64.7907, 91.8403, 99.9958, 99.8473, 
  5.67984}, {3.57622, 2.10501, 2.31041, 31.713, 5.02368, 65.0744, 
  99.9956, 99.9598, 71.4918, 1.44109}, {7.19988, 3.35608, 1.87687, 
  78.207, 26.8646, 58.6329, 99.9526, 99.6672, 99.9946, 
  0.861573}, {8.04852, 6.87097, 6.41331, 60.8064, 34.9159, 2.51863, 
  99.9941, 99.3997, 90.2107, 0.638113}, {4.18296, 4.01782, 8.92793, 
  96.085, 31.6022, 84.0935, 73.4458, 99.9977, 99.8578, 
  0.501713}, {4.5915, 6.61179, 5.71736, 22.7116, 2.04282, 66.2641, 
  99.9905, 99.9918, 99.9799, 1.10231}, {4.08263, 1.81886, 4.68386, 
  72.7076, 14.2189, 89.0024, 93.9262, 99.9085, 99.9694, 
  1.63112}, {4.21313, 6.37235, 4.02966, 74.7741, 40.9294, 54.9772, 
  96.4724, 99.8585, 92.8554, 0.239513}, {7.89187, 6.41637, 5.26262, 
  95.1501, 94.5765, 41.2555, 99.9988, 99.7735, 79.7334, 
  0.520329}, {4.30137, 3.91935, 4.54562, 64.0219, 50.314, 24.4864, 
  29.4455, 99.9915, 99.9982, 0.673102}, {3.9708, 6.45849, 9.07293, 
  63.7537, 75.6585, 44.1156, 99.9967, 99.9417, 98.0813, 
  1.57786}, {4.04281, 2.89641, 7.28808, 26.959, 1.38694, 2.57197, 
  41.1733, 99.9937, 99.6898, 1.75819}, {8.1186, 7.11612, 8.46622, 
  64.8289, 64.4101, 17.9903, 12.4893, 99.9981, 97.5647, 
  2.56602}, {5.93666, 5.82491, 2.00109, 53.428, 89.7787, 55.4647, 
  64.0043, 99.9979, 99.8208, 1.65672}, {5.75574, 6.94877, 9.18209, 
  82.5037, 4.60125, 23.1514, 96.2799, 99.9979, 99.9989, 
  6.20389}, {9.42128, 8.11236, 9.0945, 90.043, 24.4238, 73.0078, 
  99.4015, 99.7692, 99.9926, 2.95099}, {4.51027, 8.66742, 3.92109, 
  71.5847, 96.7984, 69.5328, 99.958, 99.762, 99.9594, 
  2.50903}, {3.86584, 6.37478, 6.16191, 81.5083, 0.758943, 23.3043, 
  97.3344, 99.9901, 99.6073, 3.73415}, {4.84126, 5.55246, 2.8502, 
  38.9013, 20.8126, 47.6756, 99.6433, 99.9961, 99.861, 
  4.06967}, {3.1403, 9.05126, 4.84558, 28.0323, 14.8991, 74.9112, 
  23.9802, 99.9983, 98.6783, 1.67307}, {8.09042, 4.86911, 5.02845, 
  87.9317, 19.9094, 40.4908, 96.0024, 99.9957, 84.6074, 
  0.476204}, {4.95586, 3.62765, 8.51698, 22.7759, 67.156, 45.7975, 
  66.0947, 99.9947, 99.9983, 1.90454}, {9.41473, 2.66972, 4.41192, 
  18.3297, 7.01418, 32.925, 97.7745, 99.9983, 99.9593, 
  2.6264}, {2.15615, 6.85465, 7.43165, 78.6018, 9.09551, 99.1833, 
  99.0762, 99.2573, 99.9971, 1.85394}, {3.56619, 3.62185, 7.34852, 
  46.3311, 98.7194, 85.6877, 96.6561, 97.3831, 93.0697, 
  0.254493}, {3.00693, 8.65827, 7.14155, 53.5078, 33.974, 30.9188, 
  99.9866, 99.9988, 99.9914, 2.91214}, {3.06525, 5.16252, 8.29072, 
  20.7764, 57.5311, 10.6992, 99.9727, 98.5274, 99.5453, 
  0.891177}, {7.63223, 4.07721, 1.96054, 39.0547, 40.4339, 95.4157, 
  99.9859, 99.9652, 99.2232, 1.77152}, {7.29558, 5.64949, 2.8811, 
  5.44578, 27.9093, 16.0614, 98.8066, 99.998, 99.9933, 
  4.58653}, {4.28815, 7.41513, 1.94484, 7.74768, 2.67541, 0.351201, 
  99.9864, 99.9916, 99.9911, 5.79955}, {1.80885, 3.64285, 6.68724, 
  86.1387, 46.5206, 48.9943, 99.9986, 99.9635, 99.9314, 
  0.886434}, {4.59179, 3.23396, 6.21096, 20.4324, 5.56389, 0.827754, 
  49.087, 55.1668, 99.9984, 5.52731}, {2.34635, 6.08839, 5.47863, 
  93.4553, 84.8885, 1.44247, 63.3604, 20.7211, 99.9963, 
  3.67058}, {8.67904, 9.0719, 8.94064, 22.3645, 0.788655, 61.1947, 
  8.8939, 99.8468, 99.8533, 2.25064}, {4.37183, 4.13533, 7.17843, 
  97.5977, 82.8214, 45.5433, 2.65429, 99.9979, 71.3277, 
  0.0907113}, {6.03657, 3.38797, 5.11441, 87.9552, 13.5933, 2.60125, 
  99.1886, 99.8365, 99.989, 5.48895}, {3.40663, 7.40894, 8.59269, 
  41.7428, 24.0868, 21.4381, 41.6435, 99.9931, 99.9941, 
  3.76347}, {8.44154, 2.1391, 2.8496, 8.90238, 6.85588, 4.35837, 
  21.2213, 99.9637, 99.6076, 3.42671}, {7.60202, 2.68719, 4.30726, 
  97.6691, 20.5688, 89.1342, 91.5111, 99.9652, 99.844, 
  5.15461}, {6.12289, 5.97563, 9.48301, 78.9959, 18.4541, 22.7614, 
  60.2269, 99.9862, 99.9499, 4.3972}, {7.44264, 7.71285, 6.08001, 
  64.2758, 51.7434, 2.93313, 88.88, 93.4001, 99.2909, 
  3.55898}, {7.37178, 8.87779, 6.6448, 81.2205, 15.8779, 4.18643, 
  41.2351, 99.9601, 92.2717, 1.32915}, {7.59748, 9.32757, 7.12076, 
  14.2005, 0.245115, 57.5699, 91.9871, 99.7947, 99.9979, 
  4.31516}, {3.01153, 6.83957, 7.7081, 32.5739, 73.9266, 46.4212, 
  99.9064, 99.9819, 38.8691, 2.09092}, {6.40414, 7.00064, 6.88734, 
  22.5594, 2.59197, 81.8061, 88.9638, 99.9866, 99.9358, 
  2.18183}, {7.0047, 2.39931, 5.67804, 66.788, 48.8115, 33.007, 
  97.0608, 99.9418, 90.5477, 1.83165}, {1.70618, 8.15528, 6.2876, 
  58.9402, 82.8962, 28.1009, 99.6348, 99.9299, 99.9634, 
  1.6998}, {2.62647, 6.60293, 5.364, 58.5981, 17.3428, 40.7356, 
  93.5446, 99.9284, 99.7406, 4.36042}, {2.50892, 4.50689, 6.13296, 
  48.4049, 0.0676907, 14.9763, 44.071, 99.9956, 99.9834, 
  1.37487}, {2.29687, 2.17293, 7.53671, 33.2556, 69.2853, 13.503, 
  99.3481, 99.9715, 99.1699, 1.8859}, {3.40663, 7.40894, 8.59269, 
  41.7428, 24.0868, 21.4381, 41.6435, 99.9931, 99.9941, 
  3.76347}, {7.60149, 7.92392, 9.43736, 33.6019, 48.586, 44.7875, 
  99.9986, 99.999, 98.791, 5.92025}, {9.35331, 8.78697, 2.85835, 
  83.0237, 59.4489, 17.4363, 96.855, 74.3414, 99.9904, 
  0.474819}, {1.99748, 5.57573, 7.28023, 85.944, 30.2311, 88.1763, 
  92.6333, 88.7521, 99.9908, 1.14854}, {8.06659, 7.5091, 5.12076, 
  32.0864, 46.0467, 97.5903, 99.9471, 99.9805, 99.6783, 
  3.98772}, {8.62156, 3.99597, 5.03718, 71.9549, 96.5722, 4.33995, 
  99.8785, 99.5954, 99.896, 1.82725}, {8.21073, 4.29281, 2.48654, 
  92.1154, 42.2122, 33.7003, 99.9731, 99.8616, 70.2982, 
  3.02936}, {2.12293, 8.64597, 1.50045, 52.7673, 7.03622, 85.2054, 
  99.0512, 97.6231, 99.9989, 2.52192}, {4.70209, 9.49757, 8.0795, 
  26.1104, 79.0581, 25.9698, 99.9555, 99.9204, 90.3236, 
  3.50739}, {8.40554, 4.13941, 6.26039, 75.4673, 14.7746, 30.3462, 
  99.918, 99.9542, 99.992, 3.89145}, {8.02054, 8.55722, 4.62459, 
  19.093, 63.1166, 80.7548, 99.9874, 99.9988, 99.6507, 
  5.08589}, {5.00931, 5.27433, 3.54218, 3.2995, 78.7094, 75.6799, 
  99.9816, 99.4067, 99.9629, 6.20381}, {4.54218, 9.38192, 2.30052, 
  7.31615, 23.4484, 95.269, 99.8918, 99.8577, 3.73309, 
  3.4914}, {8.09059, 1.63799, 2.17583, 7.20271, 66.9868, 60.829, 
  90.8138, 99.8832, 99.9311, 1.37417}, {8.10848, 7.55569, 5.49387, 
  8.51411, 18.3697, 59.1282, 93.9798, 99.9634, 97.2866, 
  2.57122}, {3.35748, 6.73568, 3.00201, 5.92528, 49.1173, 36.5895, 
  99.9674, 99.9981, 99.9689, 5.57522}, {4.24574, 3.96488, 8.27671, 
  3.88158, 41.0361, 71.6264, 99.99, 91.8686, 87.3365, 
  2.53283}, {9.14054, 9.46186, 4.59029, 12.2438, 76.6261, 67.2662, 
  99.9583, 89.2674, 97.5297, 0.596405}, {3.98249, 6.52713, 7.38191, 
  22.6486, 27.9321, 0.173861, 99.9973, 97.3058, 97.59, 
  1.36853}, {4.92255, 3.97742, 7.30116, 14.7848, 27.5801, 9.78051, 
  99.9957, 99.9826, 94.0391, 1.37089}, {2.44611, 8.8879, 1.64598, 
  28.7825, 76.3194, 36.7545, 99.767, 99.9068, 88.4087, 
  0.304613}, {9.00868, 6.08745, 9.46507, 42.7424, 99.418, 78.5111, 
  98.6082, 99.9942, 99.9832, 5.66273}, {2.61167, 3.3706, 8.0893, 
  76.3216, 4.96106, 55.8986, 95.1788, 99.9943, 29.7518, 
  1.8122}, {4.33232, 1.97586, 4.0075, 59.7703, 9.87178, 83.4126, 
  99.6145, 99.9454, 97.033, 4.01024}, {6.34695, 9.09703, 4.2607, 
  69.1521, 0.474443, 61.8155, 98.9945, 99.9981, 99.9837, 
  5.44662}, {7.74035, 5.43369, 4.7712, 13.4951, 21.6419, 65.6951, 
  99.9309, 99.8925, 99.9961, 7.62992}, {3.5673, 7.1409, 3.01119, 
  64.7896, 71.7061, 91.3339, 99.9694, 96.9165, 99.9884, 
  3.23233}, {2.24763, 6.66019, 5.22456, 41.9074, 97.9161, 15.0735, 
  94.9915, 99.0848, 99.7698, 2.33747}, {2.36144, 2.34403, 2.70347, 
  46.9886, 75.7711, 19.6808, 99.7677, 99.9966, 99.8866, 
  4.2401}, {5.74194, 2.74304, 6.58683, 53.1441, 97.6277, 80.8848, 
  92.1167, 93.8893, 73.7756, 1.97752}, {4.77949, 8.47201, 4.9314, 
  67.1359, 67.516, 53.5443, 99.7951, 99.998, 99.995, 
  5.06378}, {7.96353, 6.92866, 5.86909, 78.1488, 5.01099, 39.2779, 
  99.3583, 99.9844, 97.9908, 2.57833}, {2.59691, 7.75894, 4.68982, 
  23.3573, 55.6155, 44.6937, 39.0912, 99.9951, 47.5557, 
  1.77592}, {2.97662, 7.93488, 5.50166, 13.0154, 52.4564, 57.557, 
  94.0594, 99.9974, 99.7756, 2.91149}, {3.10634, 8.45944, 8.85597, 
  61.8798, 2.99506, 5.47342, 96.6643, 99.9849, 99.8912, 
  3.33391}, {4.10993, 5.55887, 2.39008, 26.0869, 13.5394, 80.319, 
  35.3354, 99.9446, 99.8603, 2.3829}, {8.29929, 8.43892, 2.75217, 
  92.384, 23.177, 20.4054, 13.9011, 76.9172, 47.232, 
  0.317526}, {1.74279, 7.06368, 2.12305, 93.6107, 59.5389, 21.3121, 
  99.9987, 99.9957, 67.0618, 2.7333}, {8.07755, 4.79852, 9.3047, 
  88.6784, 9.82053, 10.6528, 98.0251, 99.9826, 99.9984, 
  1.98553}, {7.35769, 6.13095, 3.16914, 42.4178, 39.6633, 43.0004, 
  91.012, 99.9986, 99.9969, 3.36847}, {2.22681, 4.04375, 5.15608, 
  76.7037, 16.0346, 65.8873, 99.4531, 99.8588, 99.7911, 
  2.80963}, {5.19194, 4.40342, 4.63089, 14.1441, 56.0962, 37.3248, 
  61.4923, 94.5646, 99.9988, 3.5357}, {3.36683, 7.96797, 2.14212, 
  18.9348, 43.8551, 7.4812, 79.8588, 99.8319, 99.9961, 
  4.59124}, {5.23371, 2.33771, 3.56199, 2.69287, 63.9188, 48.9992, 
  98.617, 85.1346, 94.3789, 7.9415}, {9.27604, 2.18206, 9.11544, 
  68.8873, 2.27447, 68.9287, 86.6299, 98.8338, 99.9315, 
  3.93646}, {8.41856, 1.63889, 7.40906, 72.1709, 39.1747, 18.3598, 
  46.3754, 49.09, 99.9369, 1.87566}, {8.18875, 4.77217, 2.57405, 
  5.69575, 33.2058, 19.1148, 92.4681, 99.6173, 99.0641, 
  6.78586}, {5.49926, 3.33896, 7.01447, 25.4559, 50.4373, 1.29872, 
  98.3124, 99.5204, 99.1888, 11.3243}, {1.69227, 4.09975, 7.39314, 
  62.564, 45.7367, 34.8869, 99.8993, 99.2479, 99.9965, 
  3.43916}, {5.7313, 3.2032, 6.68056, 99.4135, 84.5362, 54.7399, 
  97.3012, 99.9954, 95.9113, 1.66782}, {8.21561, 8.95252, 1.50047, 
  13.3891, 76.7185, 22.8239, 90.0687, 99.6195, 99.9751, 
  4.92659}, {3.55233, 4.11058, 8.19377, 86.4736, 38.507, 56.9374, 
  14.9708, 99.1806, 99.8982, 0.719405}, {8.2179, 8.33442, 6.15574, 
  58.1628, 37.5233, 7.95139, 99.9985, 99.998, 99.9967, 
  4.27566}, {3.79495, 8.45828, 3.85569, 93.7994, 22.8551, 84.2106, 
  99.9977, 99.9971, 99.3894, 3.59364}, {6.3049, 4.39136, 6.28408, 
  80.3853, 3.059, 51.3386, 99.9982, 99.9961, 99.994, 
  7.77428}, {5.20796, 3.12107, 3.37735, 36.8342, 11.7245, 15.4064, 
  99.9934, 99.8573, 70.1734, 4.46193}, {2.5561, 6.93066, 1.98418, 
  52.341, 27.2714, 59.9739, 99.2437, 99.9863, 91.4549, 
  1.40394}, {9.20581, 7.82841, 1.91565, 73.4332, 6.09666, 47.8882, 
  99.99, 99.6677, 68.8957, 2.84018}, {8.49142, 4.8934, 2.54901, 
  27.736, 26.5347, 0.851412, 99.1393, 99.9986, 94.8133, 
  1.84519}, {1.50443, 1.51369, 5.15044, 80.0924, 82.0089, 30.9841, 
  66.5867, 92.2235, 99.9744, 0.601006}, {2.28371, 8.38251, 4.82825, 
  81.5278, 45.5205, 78.6621, 99.9384, 96.7052, 58.987, 
  1.56813}, {9.47872, 7.31259, 9.00488, 50.6832, 1.92908, 40.5155, 
  51.4801, 99.8201, 99.995, 4.50164}, {5.23493, 9.10444, 8.5701, 
  52.2922, 17.7846, 42.1136, 59.9657, 99.996, 99.941, 
  4.63387}, {2.89093, 5.92061, 4.01907, 77.6319, 34.011, 67.7192, 
  99.7205, 99.9746, 79.5243, 3.86807}, {8.2179, 8.33442, 6.15574, 
  58.1628, 37.5233, 7.95139, 99.9985, 99.998, 99.9967, 
  4.27566}, {7.88067, 8.34754, 3.94803, 87.797, 31.6479, 94.5544, 
  67.3695, 99.8122, 99.2328, 3.90403}, {6.3049, 4.39136, 6.28408, 
  80.3853, 3.059, 51.3386, 99.9982, 99.9961, 99.994, 
  7.77428}, {3.19807, 4.55317, 5.5089, 38.5646, 6.26737, 28.9509, 
  99.2664, 99.9667, 99.9901, 5.0656}, {1.59107, 5.85093, 3.42501, 
  44.9876, 16.5302, 25.9782, 77.5608, 99.9988, 86.5532, 
  5.14966}, {9.02554, 2.64676, 4.97729, 93.5845, 5.44662, 76.0264, 
  5.20263, 99.9839, 99.5946, 1.01089}, {7.09657, 3.93928, 5.12548, 
  55.2294, 1.19625, 97.3552, 85.6011, 99.9974, 99.8243, 5.34807}};

(*randall-29042013*)
(* gamma 1,2,3 ec50 4,5,6 emax 7,8,9 1/alpha 10 *)
WPDATold = {{5.42343, 7.82418, 9.05988, 44.9941, 0.329823, 54.3449, 95.158, 
  99.6025, 99.2221, 3.61581}, {8.40389, 9.33133, 3.36183, 76.5061, 
  2.9751, 62.2868, 99.9972, 99.9827, 99.9812, 5.71557}, {4.68182, 
  6.23802, 2.89062, 56.2508, 0.72896, 82.9546, 53.6156, 99.5572, 
  99.0838, 3.24602}, {9.36408, 8.86381, 4.07153, 24.9403, 2.09547, 
  90.724, 99.4914, 99.3872, 99.993, 1.05936}, {2.15284, 7.43136, 
  3.09007, 43.5992, 33.8943, 26.0492, 99.9454, 99.995, 99.9564, 
  5.28075}, {2.63991, 8.50424, 5.64276, 11.1708, 22.4289, 61.7505, 
  99.3079, 74.7955, 99.9926, 0.0542775}, {8.23716, 5.05332, 2.83279, 
  90.9142, 4.29998, 1.63821, 99.9731, 99.927, 98.9149, 
  4.86457}, {6.28761, 2.22787, 2.09095, 69.4104, 2.3817, 25.4591, 
  99.9829, 99.9285, 99.9989, 5.30785}, {6.37916, 7.56168, 8.59668, 
  3.05894, 1.65549, 19.6139, 66.6134, 91.1744, 99.966, 
  1.70654}, {6.43796, 3.14232, 7.63836, 52.8518, 35.555, 70.0107, 
  99.9976, 99.9894, 99.9975, 3.41022}, {2.40293, 2.48239, 1.94417, 
  59.0533, 2.77191, 73.2739, 99.7906, 89.4005, 99.8075, 
  1.83032}, {9.17904, 8.94836, 3.85801, 80.7798, 10.9455, 44.5071, 
  99.9961, 99.9895, 99.3898, 4.83266}, {7.93373, 7.83238, 6.80784, 
  53.6554, 32.6302, 11.7288, 99.9936, 98.1161, 99.6155, 
  3.4422}, {7.20397, 3.98519, 6.82427, 72.253, 61.3176, 56.7029, 
  99.9984, 99.7429, 96.4066, 4.80005}, {4.05175, 9.02749, 5.61823, 
  34.1391, 80.7539, 99.2389, 99.3971, 99.7841, 98.0312, 
  4.61868}, {1.89905, 7.48888, 2.47309, 5.98643, 2.03184, 73.4825, 
  99.9805, 99.6706, 99.998, 4.97719}, {2.46403, 3.78763, 5.41791, 
  62.3111, 84.9618, 62.0889, 99.8738, 99.9963, 99.5831, 
  2.96787}, {9.00293, 7.24281, 5.82232, 45.514, 49.1486, 29.4127, 
  99.9655, 99.97, 99.6605, 5.7785}, {6.1322, 9.38718, 7.1769, 50.6748,
   62.2098, 74.7998, 99.9822, 99.864, 99.9973, 6.69555}, {7.8471, 
  4.66555, 2.15832, 44.9533, 5.84541, 20.1962, 99.9935, 93.1985, 
  98.9265, 5.99352}, {8.60401, 8.28301, 4.94222, 90.5093, 12.9309, 
  28.3358, 99.7351, 99.998, 99.9986, 2.78887}, {8.93782, 4.87389, 
  8.02019, 86.8627, 3.46723, 58.8603, 90.88, 99.961, 97.6201, 
  2.32972}, {2.50417, 3.26907, 3.15748, 36.7908, 4.6422, 27.2744, 
  84.9317, 99.881, 99.9965, 2.46448}, {5.01878, 5.09261, 7.51171, 
  82.7745, 75.8291, 51.4349, 99.9931, 99.9842, 99.9988, 
  1.95533}, {7.29261, 8.4522, 3.41311, 80.8308, 3.04444, 56.2353, 
  99.746, 99.9574, 99.6133, 4.08017}, {6.1926, 7.0979, 6.72984, 
  64.7488, 18.2829, 97.9806, 99.984, 99.9931, 99.9894, 
  0.0534263}, {1.57285, 9.3591, 4.24472, 28.2199, 6.33901, 39.882, 
  99.8544, 97.3813, 99.9986, 3.73951}, {5.96438, 7.94073, 6.05144, 
  73.1119, 14.901, 12.7774, 99.7186, 99.979, 93.8408, 
  2.0753}, {1.98106, 7.89009, 6.97999, 48.5951, 41.5299, 85.8338, 
  86.0575, 99.9976, 99.9972, 3.80189}, {1.56712, 7.93547, 7.212, 
  94.752, 88.6017, 36.1364, 94.769, 97.2303, 99.9334, 
  0.583036}, {3.23321, 5.22801, 5.60171, 85.4328, 34.9539, 58.7829, 
  99.8787, 99.7268, 99.9456, 0.717242}, {1.69176, 4.01458, 2.88695, 
  8.61557, 34.2252, 42.496, 99.969, 99.9928, 99.9984, 
  4.34743}, {7.23777, 2.7477, 5.39018, 85.4008, 46.4168, 82.9589, 
  99.313, 99.0554, 99.9909, 1.46495}, {9.34768, 9.42247, 3.68075, 
  52.6607, 6.90526, 1.14363, 99.9986, 99.9828, 99.9851, 
  3.88758}, {5.91552, 8.47278, 5.36452, 6.10466, 23.7263, 46.3227, 
  95.6691, 99.9957, 98.4888, 0.388615}, {4.34451, 4.11515, 5.18221, 
  95.7664, 23.4506, 1.29784, 59.4963, 99.8489, 99.9953, 
  0.232086}, {5.01878, 5.09261, 7.51171, 82.7745, 75.8291, 51.4349, 
  99.9931, 99.9842, 99.9988, 1.95533}, {3.95388, 5.89469, 5.84469, 
  59.9984, 54.0551, 24.7657, 99.8237, 99.9968, 99.9978, 
  3.94891}, {8.3698, 9.258, 7.16862, 80.794, 25.081, 1.93789, 99.9881,
   90.3457, 99.9591, 0.603529}, {6.53897, 5.34137, 1.88803, 73.6159, 
  95.1085, 90.5397, 98.4606, 99.9719, 99.9817, 2.10078}, {9.23804, 
  8.53276, 4.41656, 52.128, 0.946719, 22.4389, 97.4109, 99.9464, 
  99.9971, 3.81918}, {6.83858, 6.7606, 8.16094, 67.1166, 34.3897, 
  55.3815, 99.9508, 99.9478, 99.2022, 1.70534}, {1.52106, 7.8609, 
  7.60857, 90.5692, 56.5728, 46.2518, 99.9274, 99.9894, 99.0432, 
  1.83884}, {2.04415, 4.27399, 3.97745, 48.1142, 2.42805, 27.7945, 
  99.9156, 99.9189, 99.9871, 1.34293}, {6.7205, 2.87401, 8.70776, 
  82.0308, 4.3483, 94.4633, 99.9445, 99.5408, 98.1459, 
  1.23331}, {6.64341, 6.85558, 5.33014, 90.8618, 52.562, 24.624, 
  99.5644, 98.5547, 99.9957, 0.463091}, {2.49251, 5.87076, 6.78298, 
  30.4349, 47.3536, 95.7116, 99.9447, 99.9866, 99.6159, 
  3.82425}, {7.9626, 4.28262, 5.48584, 35.7327, 3.74632, 46.8289, 
  98.1774, 93.7177, 99.9928, 0.915217}, {7.99134, 4.52003, 5.73742, 
  76.4878, 10.9316, 51.3915, 99.9867, 98.4705, 86.7383, 
  1.53983}, {2.81144, 9.47506, 7.3062, 18.0721, 40.5614, 61.1026, 
  99.9227, 99.9956, 99.9989, 3.84399}, {7.43219, 5.26849, 8.26024, 
  78.0182, 65.0353, 98.3549, 56.6865, 99.2992, 99.9428, 
  0.767176}, {6.34963, 2.45528, 2.80127, 42.7313, 30.7003, 64.5111, 
  99.9983, 99.8444, 98.9138, 6.38402}, {2.14677, 9.47672, 6.04446, 
  28.7658, 65.7724, 34.8623, 99.528, 99.9976, 98.0719, 
  4.05928}, {4.39556, 8.31346, 6.30642, 74.8422, 80.1779, 60.4736, 
  91.7028, 94.7157, 85.3805, 2.20132}, {3.66443, 9.04537, 9.36893, 
  15.1057, 47.3205, 98.407, 78.8529, 99.9965, 99.9988, 
  1.92904}, {4.902, 5.23429, 8.4436, 21.1572, 51.1348, 69.1025, 
  97.318, 99.9695, 99.4061, 3.76183}, {6.77994, 6.40822, 5.53998, 
  14.978, 59.2321, 40.4844, 96.2849, 99.7207, 99.6855, 
  4.80017}, {2.45859, 8.97124, 4.78341, 50.1152, 85.026, 11.6042, 
  99.8322, 99.9836, 99.198, 1.43748}, {4.50554, 7.28306, 4.86755, 
  74.5551, 37.4872, 86.8104, 99.9982, 96.2161, 99.9768, 
  1.43761}, {4.45464, 4.84792, 6.68465, 42.4105, 59.3338, 58.0091, 
  99.9979, 99.9835, 99.9858, 4.79899}, {3.00102, 6.45503, 3.52136, 
  91.132, 8.3456, 48.4053, 21.679, 99.9782, 99.9889, 
  3.14912}, {2.35963, 6.41743, 6.90313, 72.0153, 2.22392, 17.5749, 
  94.2446, 99.972, 99.6617, 5.56934}, {7.64067, 5.29415, 5.78688, 
  93.3128, 23.0279, 38.733, 91.0267, 98.4828, 99.535, 
  1.83562}, {4.83775, 7.85894, 2.63096, 92.0844, 8.0906, 51.8109, 
  67.4413, 99.9965, 99.8087, 2.73796}, {9.13208, 5.52067, 4.71086, 
  89.9639, 2.66402, 51.502, 96.5982, 99.966, 99.4009, 
  3.22931}, {4.62731, 7.07077, 2.32192, 81.4595, 1.41984, 21.0625, 
  94.6722, 99.997, 99.998, 4.14672}, {8.21748, 3.53658, 3.97596, 
  88.0734, 13.3493, 5.61567, 99.518, 99.9384, 99.9966, 
  2.50853}, {6.38629, 6.74361, 4.14486, 66.5871, 12.6955, 20.8768, 
  99.0902, 99.9976, 99.9923, 3.97614}, {2.73145, 6.72499, 4.39451, 
  57.4823, 26.9154, 36.2258, 99.9985, 99.9821, 99.9988, 
  4.42135}, {5.16475, 3.5118, 9.49883, 89.3371, 17.3497, 46.8377, 
  99.999, 99.9848, 99.9989, 6.0013}, {6.35701, 8.66867, 3.7594, 
  84.7732, 41.6905, 58.5386, 99.5269, 99.9981, 99.9928, 
  2.14449}, {7.93829, 5.19263, 3.79075, 61.6727, 1.74552, 45.9986, 
  99.9902, 97.8178, 99.2622, 1.71552}, {4.6779, 6.48259, 5.93278, 
  62.3418, 0.688125, 24.1431, 99.7348, 98.7988, 99.928, 
  1.34453}, {6.89058, 6.65165, 4.8541, 94.4137, 82.6484, 0.695372, 
  98.9401, 99.981, 99.9844, 2.14615}, {7.51863, 5.58594, 6.44074, 
  95.0198, 23.4933, 60.8957, 99.9896, 99.9679, 99.9954, 
  3.67415}, {6.26284, 9.37439, 9.33853, 79.7853, 23.4731, 33.8574, 
  99.7785, 99.6771, 99.9977, 2.52965}, {4.86378, 3.24177, 2.51183, 
  76.8885, 5.77778, 13.782, 99.8586, 99.9983, 98.9277, 
  0.180044}, {3.90305, 4.05545, 4.12654, 89.833, 46.5645, 46.3876, 
  99.9963, 98.6823, 99.2598, 1.58412}, {9.11891, 6.72612, 8.75484, 
  70.2636, 73.9755, 48.7303, 38.9251, 99.8468, 99.8794, 
  1.85572}, {2.25362, 6.55268, 5.19854, 50.9922, 1.113, 22.059, 
  99.9972, 99.9864, 99.9628, 3.58654}, {6.35701, 8.66867, 3.7594, 
  84.7732, 41.6905, 58.5386, 99.5269, 99.9981, 99.9928, 
  2.14449}, {8.7567, 4.48745, 2.48314, 75.3618, 0.0442428, 20.1048, 
  99.3038, 99.9839, 98.0646, 2.79286}, {3.79522, 8.00724, 3.15755, 
  44.3838, 35.439, 59.0231, 99.9977, 99.9356, 99.9978, 
  2.7691}, {2.82075, 4.84785, 2.89161, 85.7013, 16.8734, 71.1759, 
  99.9542, 99.9657, 99.9629, 1.46431}, {6.47079, 8.12274, 5.93411, 
  51.1797, 81.5165, 93.2315, 99.3372, 99.7311, 99.9843, 
  1.76112}, {6.89058, 6.65165, 4.8541, 94.4137, 82.6484, 0.695372, 
  98.9401, 99.981, 99.9844, 2.14615}, {8.62872, 9.11419, 7.80186, 
  39.9582, 99.3153, 65.6166, 58.2105, 98.4461, 84.0971, 
  0.19217}, {1.75094, 3.8592, 2.86701, 22.1581, 84.5036, 94.3184, 
  99.967, 99.9989, 93.4128, 1.56146}, {8.25598, 6.15145, 7.78887, 
  24.8734, 63.0095, 12.7558, 83.6026, 99.9882, 99.8508, 
  0.0850392}, {6.54842, 6.63239, 3.24717, 64.3909, 77.796, 73.8905, 
  99.9933, 99.9979, 96.9658, 1.66135}, {3.73107, 6.45405, 7.02576, 
  82.7048, 17.5249, 93.7181, 99.8992, 99.7539, 99.9989, 
  3.97013}, {3.96261, 9.10479, 9.11612, 21.0002, 66.05, 33.5966, 
  43.3521, 99.9038, 99.9967, 1.18668}, {2.94227, 1.97971, 8.73854, 
  72.7833, 40.5635, 8.27333, 99.5548, 99.9967, 97.571, 
  2.57591}, {2.26605, 4.59335, 9.09374, 56.9462, 35.6634, 62.2665, 
  99.9918, 99.9984, 99.7006, 2.34548}, {9.40753, 8.82972, 9.29934, 
  90.1814, 50.4726, 15.1575, 99.8528, 99.9753, 97.4834, 
  2.69906}, {9.45803, 8.03571, 7.02649, 83.3349, 29.3323, 76.511, 
  99.4903, 99.862, 92.7691, 1.92114}, {6.47079, 8.12274, 5.93411, 
  51.1797, 81.5165, 93.2315, 99.3372, 99.7311, 99.9843, 
  1.76112}, {3.32349, 1.51273, 7.84947, 45.425, 21.6592, 42.8984, 
  95.8074, 99.8865, 99.9981, 1.83987}, {4.32035, 3.56004, 4.49681, 
  64.4378, 93.0023, 11.7024, 99.9392, 98.6126, 97.7887, 
  1.17359}, {3.62244, 9.24814, 2.99819, 93.8775, 21.6367, 83.5531, 
  92.9575, 99.9983, 99.9985, 4.45703}, {7.60953, 6.50124, 6.12339, 
  16.386, 13.1496, 22.9608, 43.8544, 99.9916, 99.9972, 
  1.83647}, {6.19244, 4.49547, 4.36785, 31.7549, 72.3001, 41.2945, 
  99.9987, 99.9898, 99.9933, 1.70228}, {2.90121, 6.99485, 9.33199, 
  27.6062, 41.7231, 26.2713, 99.9095, 99.9972, 99.9924, 
  2.18624}, {6.32523, 4.25935, 6.166, 42.3668, 0.255913, 30.8892, 
  99.9505, 76.9382, 93.2737, 0.851386}, {7.95035, 9.13364, 6.07099, 
  56.4888, 40.3266, 39.0374, 98.468, 99.9827, 99.9554, 
  3.35964}, {6.77924, 9.00788, 9.13704, 20.0236, 22.6341, 17.2244, 
  99.9235, 99.9807, 99.6731, 3.40328}, {2.18236, 2.54817, 4.4725, 
  55.5286, 39.3764, 65.5646, 97.6736, 99.953, 98.338, 
  0.462636}, {9.22636, 5.32741, 5.52227, 18.0654, 1.5229, 0.531061, 
  98.6428, 99.9158, 99.9723, 0.563297}, {3.48495, 1.63719, 9.33577, 
  36.2037, 28.419, 19.312, 99.9914, 86.2158, 99.9439, 
  0.700042}, {2.86886, 3.24682, 9.36568, 4.33985, 0.407499, 92.1026, 
  77.4528, 99.8814, 99.9955, 2.52442}, {5.25515, 4.95936, 1.5152, 
  50.2664, 44.2234, 81.2728, 73.4899, 99.9582, 99.6168, 
  4.28028}, {1.64663, 1.77012, 7.00267, 89.994, 18.3872, 76.3029, 
  39.1403, 99.9452, 93.9624, 0.892627}, {4.47179, 9.28405, 4.7659, 
  75.0451, 7.10452, 20.2265, 94.6835, 99.6092, 99.9625, 
  3.07874}, {4.22218, 6.37458, 7.30086, 78.2668, 29.9816, 14.9599, 
  99.8475, 99.027, 99.9848, 4.55211}, {8.37158, 3.26551, 9.18756, 
  39.6585, 48.1139, 2.32836, 37.5451, 96.8814, 98.7667, 
  2.70629}, {4.3508, 7.56421, 7.2492, 97.4252, 70.6212, 48.5639, 
  60.0043, 99.9832, 96.0945, 4.03795}, {7.56608, 8.83598, 9.27765, 
  42.557, 46.9555, 41.4823, 92.189, 99.9669, 99.2967, 
  4.17528}, {2.0152, 2.07668, 6.13796, 81.8178, 9.62361, 97.8863, 
  94.0946, 92.932, 99.3247, 2.72121}, {1.9114, 9.39761, 7.60389, 
  85.1817, 0.470713, 90.7364, 90.4545, 99.9795, 99.9956, 
  3.99393}, {4.5972, 6.36282, 9.08005, 66.9569, 36.5554, 73.2513, 
  85.3943, 99.9971, 98.7804, 3.73896}, {5.6771, 8.39232, 8.80232, 
  24.7797, 8.64293, 80.2261, 98.1567, 99.9987, 98.8606, 
  3.19322}, {6.35158, 9.18929, 2.83867, 45.4077, 91.7641, 31.6159, 
  99.9489, 99.9974, 95.1902, 2.32026}, {5.48057, 6.07737, 5.72096, 
  50.5491, 95.1263, 25.1588, 97.7547, 99.9978, 93.8402, 
  2.34373}, {7.56321, 5.01104, 2.72818, 19.3553, 35.0302, 54.0338, 
  99.9958, 99.9987, 87.9803, 3.73357}, {7.85592, 4.25077, 3.10268, 
  67.3977, 58.1893, 45.433, 99.9719, 99.9961, 99.7203, 
  1.34138}, {6.08366, 7.92359, 5.06594, 34.6899, 38.1664, 20.8079, 
  99.9289, 99.8814, 99.9491, 4.1804}, {4.89048, 4.46818, 6.47577, 
  16.8482, 22.184, 9.34782, 94.7376, 99.9881, 58.1808, 
  1.00194}, {7.87581, 9.44385, 2.49903, 80.1039, 21.1369, 28.422, 
  74.0076, 99.9924, 82.3836, 0.667942}, {3.86502, 8.45364, 8.69792, 
  74.0578, 71.1443, 99.5087, 93.9054, 99.9354, 99.9987, 
  0.903332}, {3.01075, 4.29474, 2.73233, 43.0842, 58.3853, 88.5179, 
  99.8805, 99.9341, 99.997, 2.68768}, {9.36564, 7.7445, 3.27081, 
  41.1263, 26.7761, 62.1124, 76.8347, 97.1327, 99.9286, 
  2.11509}, {5.68369, 2.75343, 2.33696, 93.4856, 16.3405, 45.1025, 
  95.0213, 94.8543, 99.9962, 1.38817}, {3.10244, 8.83732, 8.90439, 
  47.4893, 42.5321, 10.3193, 99.923, 99.9826, 96.5674, 
  4.0955}, {8.58146, 4.5631, 9.12266, 80.1707, 72.5799, 66.8751, 
  99.9984, 99.8956, 99.9862, 2.58321}, {5.36628, 4.86763, 3.34799, 
  78.4057, 13.2018, 48.413, 99.9965, 99.985, 99.9899, 
  4.20505}, {9.04233, 8.08101, 6.34039, 67.9415, 83.456, 57.0291, 
  99.9959, 99.9989, 12.4836, 3.46473}, {7.93007, 8.8259, 8.16275, 
  88.0465, 97.0801, 31.1458, 99.9575, 99.9461, 99.9735, 
  2.42972}, {3.21243, 8.34212, 9.15711, 47.4382, 53.8211, 38.6439, 
  97.1474, 99.9979, 95.8236, 2.35657}, {7.13696, 8.14008, 8.47027, 
  66.5283, 18.1746, 84.3462, 99.999, 99.9924, 99.9736, 
  4.08685}, {4.68916, 6.40536, 1.90461, 33.0123, 36.0667, 72.4098, 
  99.9974, 99.9973, 99.3211, 4.26156}, {8.05472, 3.15277, 6.47182, 
  7.10111, 8.89058, 17.6682, 99.9725, 99.8438, 99.9841, 
  5.35136}, {3.69629, 4.6906, 3.14472, 77.2106, 6.22775, 73.8975, 
  60.2139, 88.815, 99.9979, 0.0340163}, {4.85523, 9.40565, 7.65178, 
  5.50199, 25.3638, 97.9044, 99.9662, 99.9803, 99.9872, 
  6.47291}, {8.27038, 5.35301, 7.14323, 21.2227, 64.0171, 17.0228, 
  99.9979, 99.999, 99.9905, 2.90807}, {8.20636, 2.03795, 4.84133, 
  40.1357, 90.1485, 45.5459, 99.9497, 84.6377, 60.1836, 
  0.118246}, {5.42016, 8.75501, 9.34528, 14.1758, 63.2995, 63.9034, 
  99.991, 99.9839, 95.3669, 4.48359}, {6.51582, 8.31866, 8.89431, 
  1.74372, 19.0975, 35.2394, 98.7138, 99.9989, 99.5777, 
  5.60351}, {5.01579, 9.21205, 7.25092, 28.592, 57.4504, 12.7989, 
  92.4383, 99.6987, 96.0684, 0.870299}, {7.42205, 4.8246, 7.15327, 
  4.40471, 22.9466, 88.8273, 99.9666, 99.9913, 99.9983, 
  5.83477}, {8.30752, 7.9834, 4.88726, 10.9442, 6.95718, 72.494, 
  97.034, 99.8374, 99.9213, 2.35092}, {6.29518, 7.85602, 4.97815, 
  64.824, 20.4998, 32.3557, 99.7595, 99.6576, 99.997, 
  5.82594}, {9.31072, 4.52377, 6.08765, 24.581, 26.9969, 1.4734, 
  97.7352, 99.9984, 99.9963, 4.60628}, {3.78646, 7.54138, 8.76269, 
  42.1995, 83.901, 32.9972, 98.8175, 99.6753, 99.969, 
  2.51534}, {6.79339, 7.9225, 2.75267, 29.7641, 39.1797, 29.1703, 
  99.979, 99.9057, 99.9739, 5.06601}, {2.41502, 7.32489, 3.66117, 
  24.3674, 50.0984, 53.1745, 99.665, 99.9954, 99.9565, 
  2.66135}, {8.11059, 2.02733, 4.83048, 77.7158, 29.8623, 21.6056, 
  99.9806, 99.9899, 99.7809, 4.83107}, {9.22216, 2.39271, 4.15387, 
  78.5062, 23.4382, 89.2905, 99.8238, 99.0475, 99.9958, 
  2.93372}, {3.89283, 5.73695, 6.83283, 34.5102, 89.6937, 80.0263, 
  90.1741, 99.7804, 99.7444, 3.97131}, {1.51349, 8.31518, 8.59835, 
  26.6777, 97.3064, 94.5339, 99.9298, 99.9947, 98.4046, 
  3.15783}, {5.91388, 3.14472, 7.00511, 58.4374, 25.5233, 17.8801, 
  36.7254, 99.9969, 99.9989, 0.529528}, {7.37883, 7.80582, 3.56486, 
  27.0608, 32.0299, 25.4145, 22.2008, 99.9881, 99.9814, 
  2.05338}, {7.15425, 4.31503, 4.49657, 87.8414, 26.4461, 15.4222, 
  99.9945, 99.9942, 59.9875, 1.6374}, {6.78933, 7.71021, 5.42738, 
  51.0644, 23.4346, 10.4709, 99.9597, 99.9942, 96.573, 
  3.56561}, {5.42726, 8.06661, 4.94872, 47.2003, 3.80055, 6.07947, 
  99.8839, 99.9843, 99.9746, 3.43985}, {7.29454, 9.20782, 6.61482, 
  23.0903, 10.3791, 57.0153, 99.7848, 99.9809, 99.7742, 
  2.58648}, {6.85647, 9.30392, 8.29602, 86.3079, 40.0127, 97.7631, 
  99.9914, 99.9146, 99.9532, 1.11187}, {3.38192, 9.06491, 9.07175, 
  73.6998, 0.129058, 71.9719, 99.993, 99.999, 99.5598, 
  0.0887269}, {3.89079, 7.68934, 4.96513, 34.1374, 4.84671, 41.3519, 
  97.4671, 98.582, 91.0073, 0.486984}, {6.0131, 6.35753, 9.01698, 
  93.9025, 25.7201, 76.7946, 99.9885, 99.991, 99.9875, 
  1.74044}, {3.90313, 2.85983, 6.8574, 28.871, 48.4695, 71.095, 
  98.9007, 99.9908, 99.998, 1.58404}, {8.19047, 3.59348, 7.31048, 
  91.1416, 8.01978, 47.8396, 56.2002, 98.9869, 99.9846, 
  3.18693}, {1.88723, 6.95135, 6.41864, 98.1718, 5.37623, 22.921, 
  95.7612, 98.1193, 99.9987, 3.57909}, {2.1885, 4.54648, 9.05606, 
  67.9703, 10.2835, 4.80778, 99.9688, 99.3363, 99.9885, 
  4.4185}, {5.18293, 4.14445, 8.64303, 17.786, 0.816339, 87.0541, 
  85.5608, 99.4432, 99.995, 5.15254}, {3.20009, 8.31918, 3.41134, 
  87.0026, 1.8533, 24.9227, 66.4883, 99.5117, 99.8458, 
  5.52355}, {2.44529, 7.52169, 4.61081, 96.2277, 66.4129, 0.0710521, 
  99.9101, 99.9668, 99.963, 5.98205}, {2.46964, 7.12046, 6.3912, 
  82.3847, 1.76468, 18.2442, 99.9987, 99.9882, 99.9895, 
  8.30305}, {8.6039, 7.14484, 7.09388, 93.5039, 7.16031, 0.604178, 
  29.7391, 99.9635, 99.9671, 1.7561}, {5.13981, 5.6487, 1.68684, 
  1.13975, 33.5399, 60.8159, 98.3015, 84.0276, 99.8863, 
  5.33229}, {4.50732, 5.12167, 7.61528, 82.1747, 6.5602, 22.0956, 
  98.858, 99.9956, 99.9918, 5.10145}, {2.37703, 8.22347, 7.94758, 
  31.5176, 1.2071, 82.782, 99.6394, 99.9533, 99.953, 
  3.32898}, {8.19883, 7.78783, 3.62287, 84.7512, 93.4215, 82.9446, 
  99.2796, 89.8296, 99.8157, 1.52342}, {8.8476, 3.21167, 7.56194, 
  92.9835, 85.1061, 92.9659, 98.9746, 99.9926, 99.9918, 
  2.55486}, {4.57418, 2.7745, 3.76513, 91.8839, 24.575, 40.2422, 
  99.9471, 99.9387, 99.8233, 1.90642}, {6.81387, 7.4903, 7.21084, 
  47.9711, 0.0938269, 55.1364, 96.6685, 99.9936, 99.9902, 
  5.78368}, {8.50724, 9.36438, 7.99663, 67.809, 46.8395, 0.60067, 
  99.993, 99.7875, 99.9919, 0.831683}, {6.38024, 8.05635, 3.06172, 
  26.9246, 4.45453, 66.2476, 98.5592, 99.9931, 99.7649, 
  3.14441}, {3.3114, 8.57459, 9.18162, 56.0049, 90.8517, 99.9895, 
  99.9974, 99.9814, 79.5952, 4.17873}, {7.14249, 7.02208, 4.39958, 
  5.50075, 42.0972, 94.0981, 99.9729, 89.3775, 99.9974, 
  3.06181}, {6.13684, 2.6087, 2.83233, 62.3323, 19.9714, 77.5092, 
  99.4711, 99.093, 99.8356, 3.22332}, {3.78802, 1.8289, 8.27538, 
  19.5744, 14.1619, 58.4756, 2.12028, 99.9762, 77.2782, 
  0.191023}, {2.20214, 9.24315, 3.35792, 45.7033, 14.7732, 20.5733, 
  81.2384, 99.9821, 85.5638, 0.114677}, {8.07996, 8.14634, 8.65129, 
  76.9516, 33.1731, 3.8808, 97.9998, 99.4556, 89.1094, 
  2.98124}, {5.28265, 9.47566, 5.49567, 76.3101, 4.69308, 35.3974, 
  99.9937, 99.9354, 99.9858, 2.47232}, {4.6494, 4.8952, 4.37379, 
  98.3046, 2.55692, 11.653, 95.8449, 99.8779, 99.9983, 
  8.6116}, {8.82059, 4.98947, 7.24892, 85.5897, 0.851099, 29.478, 
  94.6528, 99.8826, 99.7717, 2.43585}, {1.95306, 5.7252, 7.4486, 
  45.6626, 14.5315, 87.1724, 99.9933, 99.9978, 99.9788, 
  7.05918}, {7.14211, 5.15312, 3.2553, 43.7696, 4.92801, 82.8677, 
  76.1672, 99.9766, 99.9936, 3.58517}, {7.54449, 3.05092, 4.69485, 
  88.804, 0.88617, 56.6909, 99.9904, 99.9974, 99.617, 
  1.47596}, {6.12855, 7.51343, 4.55836, 25.0907, 5.05831, 22.385, 
  99.3742, 99.9294, 97.9959, 0.140223}};



(* model 2.71 to make RI of Pailin the same as in WangPha 
runsteps = run steps per RI {RI_R,RI_T,RI_S} 
 *)
BatchRunDfixedRI[parafile_String, concfile_String, runsteps_Integer, outfilename_String] := 
  Module[{outfile, outstream, runid, initN, PMR, 
    Mu, Sigma, LC, rKZb, rKZe, tKZb, tKZe, sKZb, sKZe, everyH, Ndrug, 
    GammaR, GammaT, GammaS, ec50R, ec50T, ec50S, eminR, eminT, eminS, 
    emaxR, emaxT, emaxS, T,rmsd,riR,riT,riS,dorfrac, dortime,countrun,
    j,dx},
   
   outfile = outfilename <> ".csv";
   (*parafile = "para" <> id <> ".csv";
   concfile = id <> ".csv";
   *)
   (***outstep=0 don't print run steps,1 print run steps***)
   
   (*Print["Data Id:",id ,"\tInput files: ",parafile,"\t",concfile];*)

   If[MemberQ[FileNames[], parafile] == False,
    	Print[parafile, " is not in the working folder."];
    	Abort[];
    ];
   If[MemberQ[FileNames[], concfile] == False,
    	Print[concfile, " is not in the working folder."];
    	Abort[];
    ];
   
   If[MemberQ[FileNames[], outfile] == False,
    	Print["Creating an output file.\n"];
    
    	(****write header of the output file****)
    	outstream = OpenWrite[outfile, PageWidth -> 350];
    
    	WriteString[outstream, "no,initN,PMR,Mu,Sigma,LC,rKZb,rKZe,tKZb,tKZe,sKZb,sKZe,everyH,Ndrug,GammaR,GammaT,GammaS,ec50R,ec50T,ec50S,eminR,eminT,eminS,emaxR,emaxT,emaxS,1/alpha,rmsd,riR,riT,riS,dorfrac,dortime\n"];
    
    	Close[outstream];
    	,
    	Print[outfile, " exists already."];
    	(*Print[outfile," has been copied to ","BAK-"<>outfile];
    	CopyFile[outfile,"BAK-"<>outfile];*)
    	Abort[];
    ];
    
    dx=0.001;
    countrun=1;
For[j=1,j<=Length@WPDAT,j++,    
   	For[runid = 1, runid <= runsteps, runid++, 
		Print["\n===Run no. " <> ToString[countrun] <> "==="];
        countrun=countrun+1;  
    	(**parasite load*)
    	initN = 10^RandomReal[NormalDistribution[Log10@(2.8*10^11), 0.5]] // Abs;
    
    	(**PMR**)
    	PMR = RandomInteger[{6, 12}];
    
    	(***Mu*)
    	Mu = RandomInteger[{4, 16}];
    
    	(**Sigma**)
    	Sigma = RandomInteger[{2, 8}];
    
    	(**life cycle**)
    	LC = 48;
    
    	(**kill zone*)
    	rKZb = 6;
    	(*rKZe = RandomInteger[{20, 30}];*)
    	rKZe= 26;
    	(*tKZb = rKZe + 1;*)
    	tKZb = 27;
    	tKZe = 38;
    	(*sKZb = tKZe + 1;*)
    	sKZb = 39;
    	sKZe = 44;

   (* gamma 1,2,3 ec50 4,5,6 emax 7,8,9 1/alpha 10 *)
   
		GammaR = RandomReal[{WPDAT[[j,1]]-dx,WPDAT[[j,1]]+dx}];
     	GammaT = RandomReal[{WPDAT[[j,2]]-dx,WPDAT[[j,2]]+dx}];
     	GammaS = RandomReal[{WPDAT[[j,3]]-dx,WPDAT[[j,3]]+dx}];
     	
     	ec50R = RandomReal[{WPDAT[[j,4]]-dx,WPDAT[[j,4]]+dx}];
     	ec50T = RandomReal[{WPDAT[[j,5]]-dx,WPDAT[[j,5]]+dx}];
		ec50S = RandomReal[{WPDAT[[j,6]]-dx,WPDAT[[j,6]]+dx}];
		
		{eminR,eminT,eminS}={0.0,0.0,0.0};
		emaxR = RandomReal[{WPDAT[[j,7]]-dx,WPDAT[[j,7]]+dx}];
		emaxT = RandomReal[{WPDAT[[j,8]]-dx,WPDAT[[j,8]]+dx}];
		emaxS = RandomReal[{WPDAT[[j,9]]-dx,WPDAT[[j,9]]+dx}];
		
		T = RandomReal[{WPDAT[[j,10]]-dx,WPDAT[[j,10]]+dx}]//Abs;
		
    	everyH = 24; Ndrug = 7;
        
        (* about 1 %  *) 
      	dorfrac=RandomReal[0.01];
      	(* about a week 168 hours *)
      	dortime=RandomInteger[{1,168}];
      	
    	rmsd = ParaFitD[parafile, concfile, initN, PMR, Mu, Sigma, 
      		LC, {{rKZb, rKZe}, {tKZb, tKZe}, {sKZb, sKZe}}, everyH, 
      		Ndrug, {GammaR, GammaT, GammaS}, {ec50R, ec50T, ec50S}, {eminR, 
       		eminT, eminS}, {emaxR, emaxT, emaxS}, T, dorfrac, dortime, 0];
    	
    	riR=RI[ec50R,emaxR,GammaR,T]; (*  T*ec50R/(GammaR*emaxR)*)
    	riT=RI[ec50T,emaxT,GammaT,T]; (*T*ec50T/(GammaT*emaxT)*)
    	riS=RI[ec50S,emaxS,GammaS,T]; (*T*ec50S/(GammaS*emaxS)*)
    	
    Print["RMSD : ", rmsd];
    
    If[(rmsd != "NSL" || rmsd != 10^10) && rmsd <= 1.0,
    	outstream = OpenAppend[outfile, PageWidth -> 350];
    	WriteString[outstream, runid, ",", initN // FortranForm, ",", PMR, 
      		",", Mu, ",", Sigma, ",", LC, ",", rKZb, ",", rKZe, ",", tKZb, 
      		",", tKZe, ",", sKZb, ",", sKZe, ",", everyH, ",", Ndrug, ",", 
      		GammaR, ",", GammaT, ",", GammaS, ",", ec50R, ",", ec50T, ",", 
      		ec50S, ",", eminR, ",", eminT, ",", eminS, ",", emaxR, ",", 
     	 	emaxT, ",", emaxS, ",", T, ",", rmsd,",",riR,",",riT,",",riS,",",dorfrac//FortranForm,",",dortime, "\n"];
     
     Close[outstream];
     Print["The set of parameters is collected."];
     ,
     Print["The set of parameters is rejected."];
     ];
    ];
];
	Print["The program is finised."]
];


(* DHA efficacy for each stage from Hoshen et al. Parasitology (2000). *)
HoshenDHAEff[age_Integer]:=Module[{fn},
	fn=Piecewise[{
		{0., 0 <= age < 3}, 
		{(age - 3.)/8., 3 <= age <= 11}, 
		{1., 11 < age <= 32}, 
		{(-age/8.) + 5., 32 < age <= 40}, 
		{0., 40 < age <= 48}}
	];
	fn	
];

(* calculate the stage-specific time-dependent kill-rate *)
HoshenPDEff[conc_,ec50_]:= Module[{pd},
	pd=conc/(conc+ec50);
	pd
];

(* model 1.3 Hoshen et al. Parasitology(2000) *)
Hoshen[initN_, PMR_Integer, mu_, sigma_, hours_Integer, datafile_, everyH_Integer, 
	Ndrug_Integer, ec50_, runMax_Integer, dorfrac_, dortime_, outfile_] := 
  	Module[{runs, concls, i, j, lst, time, output, dorls, dorperd},
   
   (*  
   dorfrac  = fraction of parasites to be dormant
   dortime = dormancy period (hours) 
   	*)
   
   (**outfile = 1 write the output to modgen.csv ; 0 = 
   don't write output file.*)
   
   	runs = True;(*running status*)
   	time = runMax;(*time for calculating the concentrations(hr)..it is the maximum time that the program can run!*)
   	dorperd = 0; (*dormancy period counter*)
   	dorls = {};(* dormant parasites at 6-26 hours*)
   
   	output = {}; (*list of the parasite distribution at each time step*)

   	(*template of the drug and its effects*)
  	concls = ConcMod2[datafile, everyH, Ndrug, ec50];
   
   	(*initial parasite load*)
   	lst = DistributeN[initN, hours, mu, sigma];
    
   	(*tmp1=Log[10,Total[lst]];*) (*count all*)
   	(*tmp1=Log[10,Total[lst,{1,38}]]//Abs; (*count up to Schizonts!!!*)*)

   	(* sum all 1-24 Hosen et al. (2000) *)
   	(*tmp1 = Log10[CountRing[lst, {1, 24}, 0]] // Abs;*) 
   	   
   	(*adding the values*)
   	AppendTo[output, lst];
   
   	i = 0;
   
   	While[runs == True && i < time,
    (*evolving the system*)
    	i = i + 1;
    
    	lst = Shiftonehour[lst, PMR]; (*Parasites are growing. Feed them!*)    
    
    	(*become dormant*)
		dorls = DorCollect[dorls, dorfrac*lst~Join~{0}];
    	lst = (1. - dorfrac)*lst;
    
        (* dormancy parasites can't be killed and can not be observed in blood sample Hoshen et al Parasitology (2000)*)
    	If[dorls[[1, Length@dorls[[1]]]] == Round@RandomReal[NormalDistribution[dortime,2]],
     		lst = lst + dorls[[1, 1 ;; Length@dorls[[1]] - 1]];
     		dorls = Drop[dorls, 1];
     	];
        
    	(*a time to kill.*)
    	For[j = 1, j <= Length[lst], j = j + 1,
     		lst = ReplacePart[lst, j -> (1.-HoshenDHAEff[j]*concls[[i,3]])*lst[[j]]];
     	];
    
    	(* from Hoshen et al. Parasitology(2000) dormancy parasites can not be observed from blood sample *)
    	(*lst[[6 ;; 26]] = lst[[6 ;; 26]] + dorfrac*lst[[6 ;; 26]];*)
       
    	(*adding a point for ploting*)
    	Developer`ToPackedArray[AppendTo[output, lst]];
    
	];
   
   	(* write the output to file*)
   	If[outfile == 1,
    	Export["modgen.csv", 
     	Table[{i - 1, Log10@CountRing[output[[i]],{1,24},0]}, {i, 1, 
       	Length[output]}], "CSV"];
    	Print["The output has been written to modgen.csv.\n"];
    ,
    	If[outfile == 2, 
     	Developer`ToPackedArray[
      	Table[{i - 1, Log10@CountRing[output[[i]],{1,24},0]}, {i, 1, 
        Length[output]}]]
     , 	
     	If[outfile == 0, output]
     	]
    ]   
    
];

(* modified Hoshen et al. Parasitology(2000) #dormant parasites can be observed in the blood sample *)
HoshenNew[initN_, PMR_Integer, mu_, sigma_, hours_Integer, datafile_, everyH_Integer, 
	Ndrug_Integer, ec50_, runMax_Integer, dorfrac_, dortime_, outfile_] := 
  	Module[{runs, concls, i, j, lst, time, output, dorls, dorperd},
   
   (*  
   dorfrac  = fraction of parasites to be dormant
   dortime = dormancy period (hours) 
   	*)
   
   (**outfile = 1 write the output to modgen.csv ; 0 = 
   don't write output file.*)
   
   	runs = True;(*running status*)
   	time = runMax;(*time for calculating the concentrations(hr)..it is the maximum time that the program can run!*)
   	dorperd = 0; (*dormancy period counter*)
   	dorls = {};
   
   	output = {}; (*list of the parasite distribution at each time step*)

   	(*template of the drug and its effects*)
  	concls = ConcMod2[datafile, everyH, Ndrug, ec50];
   
   	(*initial parasite load*)
   	lst = DistributeN[initN, hours, mu, sigma];
    
   	(*tmp1=Log[10,Total[lst]];*) (*count all*)
   	(*tmp1=Log[10,Total[lst,{1,38}]]//Abs; (*count up to Schizonts!!!*)*)

   	(* sum all 1-24 Hosen et al. (2000) *)
   	(*tmp1 = Log10[CountRing[lst, {1, 24}, 0]] // Abs;*) 
   	   
   	(*adding the values*)
   	AppendTo[output, lst];
   
   	i = 0;
   
   	While[runs == True && i < time,
    (*evolving the system*)
    	i = i + 1;
    
    	lst = Shiftonehour[lst, PMR]; (*Parasites are growing. Feed them!*)    
    
    	(*become dormant*)
		dorls = DorCollect[dorls, (dorfrac*lst)~Join~{0}];
    	lst = lst - lst*dorfrac;
    
        (* the wake up time of dormancy parasites is from a normal distribution (original)*)
        (* the modified version won't use random number *)
    	If[Last@dorls[[1]] == dortime ,
     		lst = lst + dorls[[1,1;;Length@dorls[[1]]-1]];
     		dorls = Drop[dorls, {1}];
     	];
        
                
    	(*a time to kill.*)
    	For[j = 1, j <= Length[lst], j = j + 1,
     		lst = ReplacePart[lst,j->(1.-HoshenDHAEff[j]*concls[[i,3]])*lst[[j]]];
     	];
    
    	(* modified Hoshen et al. Parasitology(2000) that dormancy parasites can be observed from blood sample *)
    	lst = lst + Total@dorls[[All,1;;Length@dorls[[1]]-1]];
        
    	(*adding a point for ploting*)
    	Developer`ToPackedArray[AppendTo[output, lst]];
    
	];
   
   	(* write the output to file*)
   	Which[outfile == 1,
    	Export["modgen.csv", 
     	Table[{i - 1, Log10@CountRing[output[[i]],{1,24},0]}, {i, 1, 
       	Length[output]}], "CSV"];
    	Print["The output has been written to modgen.csv.\n"];
    ,
    	outfile == 2, 
     	Developer`ToPackedArray[
      	Table[{i - 1, Log10@CountRing[output[[i]],{1,24},0]}, {i, 1, 
        Length[output]}]]
     , 	
     	outfile == 0, output
     	
     	]   
    
];


DorCollect[dorls_, doraddls_] := 
  Module[{newdorls, dorlslength},
   	newdorls = dorls;
   
   	AppendTo[newdorls, doraddls];
   
   	dorlslength = Length@newdorls;
   	If[Length@newdorls > 0,
   		newdorls[[1 ;; Length@newdorls - 1, Length@newdorls[[1]]]] = 
   		newdorls[[1 ;; Length@newdorls - 1, Length@newdorls[[1]]]] + 1;
    ];
   newdorls
];

(* model 1.4 *)
NJWhite[initN_,LC_,mu_,sigma_,pmr_,endAge_,runT_,outform_]:=Module[{lst,temp,outls},
	(* generate the parasite in a patient without treatment following NJW's model.
		outform = 0 list of the parasite numbers at the ages from 1 - endAge end at time 0 -> runT hours.
				= 1 list of the total parasite from age 1-endAge
				= 2 list of log10 of total parasite at every time step. 
	 *)
	outls={};
	temp = DistributeN[initN,LC,mu,sigma];
	lst = NestList[Shiftonehour[#,pmr]&,temp,runT];
	
	Which[
		outform==0,
			outls=lst[[All,1;;endAge]];,
		outform==1,
			outls=Total[#[[1;;endAge]]]&/@lst;,
		outform==2,
			outls=Log10@(Total[#[[1;;endAge]]]&/@lst);
		
	];
	temp=.;lst=.;
	outls
];

(* model 1.5 AMG's model*)
(* original model Anderson, May, and Gupta Parasitology (1989),99,S59-S79 *)
AMG[lambda_,mu_,beta_,alpha_,r_,d_,initx_,inity_,inits_,runT_,dt_,outform_Integer]:=Module[{sol,x,y,s,t,outls},
	(* outform 	= 0 list of the value of each parameter at t=0 -> t=runT with step dt  
	*)
	sol= NDSolve[{	x'[t] == lambda-mu*x[t]-beta*x[t]*s[t], 
  					y'[t] == beta*x[t]*s[t]-alpha*y[t], 
  					s'[t] == alpha*r*y[t]-d*s[t]-beta*s[t]*x[t], 
  					x[0] == initx, y[0] == inity, s[0] == inits}, {x, y, s}, {t, 0, runT}, 
 					MaxSteps -> Infinity];
	Which[
		outform==0,
			outls=Table[{t, x[t], y[t], s[t]} /. First[sol], {t, 0.0, runT, dt}];
	];
 	outls
 					
];


(* model 1.6 Kwiatkowski and Nowak 2 stages*)
DominicNowak1[initxy_List, lsd_List, r_, lsPS_List, runT_,outform_] := 
  Module[{d1, d2, i, outtemp, tot,Fx,p,s,x0,y0,xold,yold,xnew,ynew},
   outtemp = {};
   
   p=lsPS[[1]]; s=lsPS[[2]];
   
   x0 = initxy[[1]]; y0 = initxy[[2]]; 
   d1 = lsd[[1]]; d2 = lsd[[2]]; 
   tot = x0+y0;
   Fx=36.5+4(1-Exp[-0.0003*x0]);
   i = 0;
   xold=x0;yold=y0;
   AppendTo[outtemp, {i,x0,y0,tot,Fx}];
   While[i <= runT,
    i = i + 1;
    ynew=d1*xold*Exp[-p*xold];
    xnew=r*d2*yold*Exp[-s*xold];
    
    Fx = 36.5 + 4*(1 - Exp[-0.0003*xnew]);
    tot = xnew+ynew;
    AppendTo[outtemp, {i,xnew,ynew,tot,Fx}];
    xold=xnew;yold=ynew;
   ];
   If[outform==0,
   		outtemp
    ]
];

(* model 1.7 Kwiatkowski and Nowak 4 stages*)
DominicNowak2[initx_List, r_, h_, lsd_List, lss_List, runT_,outform_] := 
  Module[{x1, x2, x3, x4, d1, d2, d3, d4, s1, s2, s3, s4, f1, f2, f3, 
    f4, i, outtemp, x1h, x2h, x3h, x4h, tot, Fx1},
   outtemp = {};
   
   x1 = initx[[1]]; x2 = initx[[2]]; x3 = initx[[3]]; x4 = initx[[4]];
   d1 = lsd[[1]]; d2 = lsd[[2]]; d3 = lsd[[3]]; d4 = lsd[[4]];
   s1 = lss[[1]]; s2 = lss[[2]]; s3 = lss[[3]]; s4 = lss[[4]];
   Fx1 = 36.5 + 4*(1 - Exp[-0.0003*x1]);
   tot = x1 + x2 + x3 + x4;
   i = 0;
   AppendTo[outtemp, {i, x1, tot, Fx1}];
   While[i <= runT,
    i = i + 0.5;
    f1 = Exp[-s1*x1]; f2 = Exp[-s2*x1]; f3 = Exp[-s3*x1]; 
    f4 = Exp[-s4*x1];
    
    x1h = (1 - h)*r*d4*x4*f4 + h*d1*x1*f1;
    x2h = (1 - h)*d1*x1*f1 + h*d2*x2*f2;
    x3h = (1 - h)*d2*x2*f2 + h*d3*x3*f3;
    x4h = (1 - h)*d3*x3*f3 + h*d4*x4*f4;

    Fx1 = 36.5 + 4*(1 - Exp[-0.0003*x1h]);
    tot = x1h + x2h + x3h + x4h;
    AppendTo[outtemp, {i, x1, tot, Fx1}];
    x1 = x1h; x2 = x2h; x3 = x3h; x4 = x4h;
    ];
    If[outform==0,
   		outtemp
    ]
];

(* model 1.8
Estimating sequestered parasite population dynamics in cerebral 
malaria 1. Michael B.Gravenor, 2. Michael Boele van Hensbroek,and 3. Dominic 
Kwiatkowski PNAS June 23,1998 vol.95 no.13 7620-7624 *)
Seques[initn_List,lambda_List,mu_List,r_,runT_,dt_,outform_]:=
Module[{sol,outls,n1,n2,t,lambda1,lambda2,mu1,mu2,initn1,initn2},
	(* outform 	= 0 list of the value of each parameter at t=0 -> t=runT with step dt  
				= 1	
	*)

	initn1=initn[[1]];initn2=initn[[2]];
	lambda1=lambda[[1]];lambda2=lambda[[2]];
	mu1=mu[[1]];mu2=mu[[2]];
		
	sol = NDSolve[{
		n1'[t]==r*lambda2*n2[t]-(lambda1 + mu1)*n1[t], 
    	n2'[t]==lambda1*n1[t] - (lambda2 + mu2)*n2[t], 
    	n1[0]==initn1,n2[0]==initn2}, {n1, n2}, {t, 0.0, runT},MaxSteps -> Infinity
    	];

	Which[
		outform==0,
			outls=Table[{t, n1[t], n2[t]} /. First[sol], {t, 0, runT, dt}];
	];
 	outls    	
    
];


(*Sexy5 plus dormancy assume 6-26 hours parasites to be dormant

BUT the counting is wrong.
*)
Sexy6[initN_, PMR_, mu_, sigma_, hours_Integer, KZ_List, datafile_, 
   everyH_, Ndrug_, gamma_List, ec50_List, emin_List, emax_List, T_, runMax_, 
   dorfrac_, dortime_, outfile_] := 
  Module[{runs, tmp1, stages, concls, i, j, lst, lsk, 
    time, period, output, dorls, dorperd},
   
   (*input forms for KZ,gamma,ec50,emin and emax
   KZ={{rb,re},{tb,te},{sb,se}}  
   gamma={gammar,gammat,gammas}  
   ec50={r50,t50,s50}  
   emin={rmin,tmin,smin}  
   emax={rmax,tmax,smax}*)
   	(*  dorfrac  = fraction of parasites to be dormant
           dortime = dormancy period (hours) 
   	*)
   
   (**outfile = 1 write the output to modgen.csv ; 0 = 
   don't write output file.*)
   
   runs = True;(*running status*)
   
   (*time for calculating the concentrations(hr)..it is the maximum time that the program can run!*)
    time = runMax;
	
   period = 6; (*period of blood sampling *)
  
   dorperd = 0; (*dormancy period counter*)
   dorls = {};(* dormant parasites at 6-26 hours*)
   
   output = {}; (*list of the parasite distribution at each time step*)

   (*template of the drug and its effects*)
   concls = ConcMod[datafile, everyH, Ndrug,gamma, ec50, emin, emax];
 	
   (*k_i(t)*)
   lsk = Ki[T, concls];
   
   (*initial parasite load*)
   lst = DistributeN[initN, hours, mu, sigma];
   
   (*template for using k_i*)
   stages = WhichRTS[lst, KZ];
      
   tmp1 = Log10[CountRing[lst, {11, 14}, 1]] // Abs; (* 
   sum all 1-11 + 11-48 * decay fn  *)
   
   
   (*adding the values*)
   AppendTo[output, lst];
   
   i = 0;
   
   While[runs == True && i < time,
    (*evolving the system*)
    i = i + 1;
    
    lst = Shiftonehour[lst, PMR]; (*Parasites are growing. 
    Feed them!*)
     
    (*become dormant*)
    dorls = DorCollect[dorls, dorfrac*lst[[6 ;; 26]]~Join~{0}];
    lst[[6 ;; 26]] = (1. - dorfrac)*lst[[6 ;; 26]];
      
    (* dormancy parasites can't be killed but can be observed in blood sample*)
    If[dorls[[1, Length@dorls[[1]]]] == dortime,
     lst[[6 ;; 26]] = 
      lst[[6 ;; 26]] + dorls[[1, 1 ;; Length@dorls[[1]] - 1]];
     dorls = Drop[dorls, 1];
     ];
       
    (*a time to kill.*)
    For[j = 1, j <= Length[lst], j = j + 1,
     lst = ReplacePart[lst, j -> Fdecay[j, lst, i, lsk, stages]];
     ];
    
    (* dormancy parasites can be observed from blood sample *)
    (* old version was WRONG!!*)
    (*lst[[6 ;; 26]] = lst[[6 ;; 26]] + dorfrac*lst[[6 ;; 26]];  *)
    (* the total parasites = present + all dormant parasties *)
    lst[[6;;26]] = lst[[6;;26]] + Total@dorls[[All,1;;Length@dorls[[1]]-1]]; 
   
    (*adding a point for ploting*)
    Developer`ToPackedArray[AppendTo[output, lst]];
    
    ];

    
   (* write the output to file*)
   Which[
   	outfile==0,
   		output,
   	outfile==2,
   		Developer`ToPackedArray[
      Table[{i - 1, Log10@Total[output[[i]]]}, {i, 1, 
        Length[output]}]],
    outfile==1,
    	Export["modgen.csv", 
     Table[{i - 1, Log10@Total[output[[i]]]}, {i, 1, 
       Length[output]}], "CSV"];
    Print["The output has been written to modgen.csv.\n"];
   ]
 
];


(*Sexy6 plus dormancy, dormancy age distribution. *)
Sexy7[initN_, PMR_, mu_, sigma_, hours_, KZ_, datafile_, 
   everyH_, Ndrug_, gamma_, ec50_, emin_, emax_, T_, runMax_,dorfrac_, 
   dortime_, dormu_, dorsig_, outfile_] := 
  Module[{runs, tmp1, stages, concls, i, j, lst, lsk, time,period, output, dorls, dorperd, dordist, 
    nlst, totaldor},
   
   (*input forms for KZ,\[Gamma],ec50,emin and emax
   KZ={{rb,re},{tb,te},{sb,
   se}}  \[Gamma]={\[Gamma]r,\[Gamma]t,\[Gamma]s}  ec50={r50,t50,s50}  
   emin={rmin,tmin,smin}  emax={rmax,tmax,smax}*)
   	(*  dorfrac  = 
   fraction of parasites to be dormant
           dortime = dormancy period (hours) 
   dormu,dorsig = mean age and sd of dormant parasites
   	*)
   
   runs = True;(*running status*)
   
   time = runMax;
   
   period = 6; (*period of blood sampling *)
   
   (*lim=8.0;(*detection limit is set at 8 (log10 scale)*)*)
     
   dorperd = 0; (*dormancy period counter*)
   dorls = {};
   
   dordist = (If[# <= 10^-6, 0.0, #]&/@DistributeN[dorfrac*100., hours, dormu, dorsig])*0.01;
	   
   output = {}; (*list of the parasite distribution at each time step*)

   (*template of the drug and its effects*)
   
   concls =ConcMod[datafile, everyH, Ndrug, gamma, ec50, emin, emax];
   
   (*k_i(t)*)
   lsk = Ki[T, concls];
   
   (*initial parasite load*)
   
   lst = DistributeN[initN, hours, mu, sigma];
   nlst = lst;
   
   (*tempplate for using k_i*)
   stages = WhichRTS[lst, KZ];
   
   (*tmp1=Log[10,Total[lst]];*) (*count all*)
   (*tmp1=Log[10,Total[
   lst,{1,38}]]//Abs; (*count up to Schizonts!!!*)*)
   
   tmp1 = Log10[CountRing[lst, {11, 14}, 1]] // Abs; (* 
   sum all 1-11 + 11-48 * decay fn  *)
   
   (*adding the values*)
   AppendTo[output, lst];
   
   i = 0;
      
   While[runs == True && i < time,
    (*evolving the system*)
    i = i + 1;
    
    nlst = Shiftonehour[nlst, PMR];
    
    (*become dormant*)
    dorls = DorCollect2[dorls, (dordist*nlst)~Join~{0}, 0];
    
    nlst = nlst-(dordist*nlst);
    
    
    (*a time to kill.*)
    For[j = 1, j <= Length[lst], j = j + 1,
     nlst = ReplacePart[nlst, j -> Fdecay[j, nlst, i, lsk, stages]];
     ];
      
    (* dormancy parasites can't be killed but can be observed in blood sample*)
    If[Last[dorls[[1]]] == dortime,
     nlst = nlst + dorls[[1, 1 ;; hours]];
     dorls = Drop[dorls, 1];
     ];
    
    totaldor = Total@dorls[[All, 1;;Length@dorls[[1]]-1]];
    
    (* dormancy parasites can be observed from blood sample *)  
    lst = nlst + totaldor;
       
    Developer`ToPackedArray[AppendTo[output, lst]];
    
    ];
   (* write the output to file*)
   If[outfile == 1,
    Export["modgen.csv", 
     Table[{i - 1, Log10[Total[output[[i]]]]}, {i, 1, 
       Length[output]}], "CSV"];
    Print["The output has been written to modgen.csv.\n"];
    ,
    If[outfile == 2, 
     Developer`ToPackedArray[
      Table[{i - 1, Log10[Total[output[[i]]]]}, {i, 1, 
        Length[output]}]]
     , If[outfile == 0, output]
     ]
    ]
];

(* the dormant parasites are from killing percentage * no * dormancy fraction
   see FdecayDor
 *)
Sexy8[initN_, PMR_, mu_, sigma_, hours_Integer, KZ_List, concdatafile_, 
   everyH_, Ndrug_, gamma_List, ec50_List, emin_List, emax_List, T_, runMax_, 
   dorfrac_, dortime_, outfile_] := 
  Module[{runs, obs, stages, concls, i, lst, lsk, 
    time, period, output, dorls, dorperd, tmpdor, totdor, outfile1,outfile2,junk},
   
   (*input forms for KZ,gamma,ec50,emin and emax
   KZ={{rb,re},{tb,te},{sb,se}}  
   gamma={gammar,gammat,gammas}  
   ec50={r50,t50,s50}  
   emin={rmin,tmin,smin}  
   emax={rmax,tmax,smax}*)
   	(*  dorfrac  = fraction of parasites to be dormant
           dortime = dormancy period (hours) 
   	*)
   
   (**outfile = 1 write the output to modgen.csv ; 0 = 
   don't write output file.*)
   
   runs = True;(*running status*)
   
   (*time for calculating the concentrations(hr)..it is the maximum time that the program can run!*)
    time = runMax;
	
   period = 6; (*period of blood sampling *)
  
   dorperd = 0; (*dormancy period counter*)
   dorls = {};(* dormant parasites at 6-26 hours*)
   
   totdor = {}; (* total of dormant parasites *)
   
   output = {}; (*list of the parasite distribution at each time step*)

   (*template of the drug and its effects*)
   concls = ConcMod[concdatafile, everyH, Ndrug,gamma, ec50, emin, emax];
 	
   (*k_i(t)*)
   lsk = Ki[T, concls];
   
   (* initial parasite distribution *)
   lst = DistributeN[initN, hours, mu, sigma];
   
   (*template for using k_i*)
   stages = WhichRTS[lst, KZ];
   
   (*adding the values*)
   AppendTo[output, lst];
   
   i = 0;
   
	While[runs == True && i < time,
	    (*evolving the system*)
	    i = i + 1;
	    
	    lst = Shiftonehour[lst, PMR]; (*Parasites are growing. Feed them!*)
	    
	    (* time to kill; lst is the list of sensitive parasites *)
	    {lst, tmpdor} = FdecayDor[lst, i, lsk, stages, dorfrac];   
	
	     
	    (*become dormant*)
	    dorls = DorCollect[dorls, tmpdor~Join~{0}];
	    
	      
	    (* the dormant parasites wake up*)
	    If[dorls[[1, Length@dorls[[1]]]] == dortime,
	     lst[[6 ;; 26]] = 
	      lst[[6 ;; 26]] + dorls[[1, 1 ;; Length@dorls[[1]] - 1]];
	     dorls = Drop[dorls, 1];
	     ];
	        
	    totdor = {totdor, Total@dorls[[All,1;;Length@(dorls[[1]])-1]]//Total};
	
	    (* obs = parasites that could be counted ; sensitive + dormancy *)
	    obs = lst;
	       
	    (* dormancy parasites can be observed from blood sample *)
	    (* the total parasites = present + all dormant parasties *)
	    obs[[6;;26]] = obs[[6;;26]] + Total@dorls[[All,1;;Length@(dorls[[1]])-1]]; 
	   
	    (*adding a point for ploting*)
	    AppendTo[output, obs];    
	];

    
   (* write the output to file*)
   junk = LsDot[#, Table[PRingFunc[i, 11, 14], {i, 1, hours}]] & /@ output;

	outfile1 = Table[{i - 1, Log10@Total[output[[i]]]}, {i, 1, Length[output]}];
	outfile2 = Table[{i, Log10[junk[[i + 1]] // Total]}, {i, 0, Length[junk] - 1}];

	(*write the output to file*)
	Which[
	 outfile == 0, output,
	 outfile == 1, outfile1,
	 outfile == 2, outfile2,
	 outfile == 3, {output, outfile1, outfile2, totdor // Flatten}
	 ]

];




(* for using with Sexy7 *)
DorCollect2[dorls_, doraddls_, outfile_] := 
  Module[{newdorls,dorlslength, outname},
   (*outfile = 0 no export , = 1 export dorls.csv*)
   
   outname = "dorls.csv";
   If[outfile == 0,
    newdorls = dorls;
    
    AppendTo[newdorls, doraddls];
    
    dorlslength = Length@newdorls;
    If[Length@newdorls > 0,
     newdorls[[1 ;; Length@newdorls - 1, Length@newdorls[[1]]]] = 
       newdorls[[1 ;; Length@newdorls - 1, Length@newdorls[[1]]]] + 
        1;
     ];
    newdorls
    ,
    newdorls = dorls;
    
    AppendTo[newdorls, doraddls];
    
    dorlslength = Length@newdorls;
    If[Length@newdorls > 0,
     newdorls[[1 ;; Length@newdorls - 1, Length@newdorls[[1]]]] = 
       newdorls[[1 ;; Length@newdorls - 1, Length@newdorls[[1]]]] + 
        1;
     ];
    Export[outname, newdorls, "CSV"];
    newdorls
    
    ]
   ];

(** model 2.1 outform 0**)
SimpsonMef1[p0_,a_,k_,k1_,c0_,gamma_,c50_]:=Module[{ls,i},
	ls=Table[{i,Log10[10^p0 Exp[a i] ((c50^gamma + (c0 Exp[-k i])^gamma)/(
    c50^gamma + c0^gamma))^(k1/(k gamma))]},{i,0,80,0.1}];
    ls
];

(** model 2.2 outform 0**)
SimpsonMef2[p0_,a_,k_,k1_,c0_,gamma_,c50_,mic_]:=Module[{ls,i},
	ls=Table[{i,Log10[10^p0 Exp[a i] ((((k1 - a)/a) (c0/mic)^(-gamma) + 
     Exp[-k i]^gamma)/(((k1 - a)/a) (c0/mic)^(-gamma) + 
     1))^(k1/(k gamma))]},{i,0,80,0.01}];
	ls
];


(* show the list of parameters for the insterested model *)
ListPar[modelID_,outform_]:=Module[{parms},

   (* lists of the parameters for the interested functions *)
   Which[
   	modelID==1.0&&outform==0, (*Sexy5*)
   	parms = {"model","initn", "pmr", "mu", "sigma", "lifecycle", "killzone", "concfile", "everyh", 
    "ndrug", "gamma", "ec50", "emin", "emax", "1/alpha", "runmax", "outform"};
   	,
   	modelID==1.0&&outform==1, (*Sexy5*)
   	parms = {"model","initn", "pmr", "mu", "sigma", "lifecycle", "killzone", "concfile", "everyh", 
    "ndrug", "gamma", "ec50", "emin", "emax", "1/alpha", "runmax", "outform"};
   	,
   	modelID==1.0&&outform==2, (*Sexy5*)
   	parms = {"model","initn", "pmr", "mu", "sigma", "lifecycle", "killzone", "concfile", "everyh", 
    "ndrug", "gamma", "ec50", "emin", "emax", "1/alpha", "runmax", "outform"};
   	,
   	modelID==1.0&&outform==3, (*RunSexy5*)
   	parms = {"model","initn", "pmr", "mu", "sigma", "lifecycle", "killzone", "concfile", "everyh", 
    "ndrug", "gamma", "ec50", "emin", "emax", "1/alpha", "lim", "outform"};
   	,
   	modelID==1.1&&(outform==0||outform==1), (*ParaFit*)
	parms = {"model","parafile","concfile","initn","pmr","mu","sigma","lifecycle","killzone","everyh","ndrug","gamma"
		,"ec50","emin","emax","1/alpha","outform"};
	,
	modelID==1.2&&outform==0, (*BatchRun*)   
    (*BatchRun[parafile_String, concfile_String, runsteps_Integer, outfilename_String]*)
    parms = {"model","parafile","concfile","runsteps","outname","outform"};
   ,
 	modelID==1.21&&outform==0, (*BatchRun*)   
    (*BatchRun[parafile_String, concfile_String, runsteps_Integer, outfilename_String]*)
    parms = {"model","parafile","concfile","runsteps","outname","outform"};
   ,
   	modelID==1.3&&(outform==0||outform==1||outform==2), (* Hoshen et al. Prasitology (2000)*) 
   	parms = {"model","initn","pmr","mu","sigma","lifecycle","concfile","everyh","ndrug","ec50","runmax",
   		"dorfrac","dortime","outform"};
   ,
   	modelID==1.4&&(outform==0||outform==1||outform==2),
   	parms = {"model","initn","lifecycle","mu","sigma","pmr","endage","runmax","outform"};	   
   ,
   	modelID==1.5&&(outform==0),
   	(* lambda_,mu_,beta_,alpha_,r_,d_,initx_,inity_,inits_,runT_Integer,dt_,outform_Integer *)
   	parms = {"model","lambda","mu","beta","alpha","r","d","initx","inity","inits","runsteps","stepsize","outform"};	   
   ,
    modelID==1.6&&(outform==0),
    (* initxy_List, lsd_List, r_, lsPS_List, runT_,outform_ *)
   	parms = {"model","initxy","initd","r","ps","runsteps","outform"};
   ,
   	modelID==1.7&&(outform==0),
   	(* initx_List, r_, h_, lsd_List, lss_List, runT_,outform_*)
   	parms = {"model","initx","r","h","initd","inits","runsteps","outform"};	   
   ,
   	modelID==1.8&&(outform==0),
   	(* Sques[initn_List,lambda_List,mu_List,r_,runT_,dt_,outform_] *)
   	parms = {"model","initn","lambdals","muls","r","runsteps","stepsize","outform"};	  
   ,
   	modelID==1.9&&(outform==0||outform==1||outform==2),
   	(* initN_, PMR_, mu_, sigma_, hours_Integer, KZ_List, datafile_, 
   everyH_, Ndrug_, gamma_List, ec50_List, emin_List, emax_List, T_, runMax_, 
   dorfrac_, dortime_, outfile_ *)
   	parms = {"model","initn","pmr","mu","sigma","lifecycle","killzone","concfile","everyh","ndrug","gamma","ec50"
   	,"emin","emax","1/alpha","runmax","dorfrac","dortime","outform"};
   ,
   	modelID==2.0&&(outform==0||outform==1||outform==2),
   	parms = {"model","initn","pmr","mu","sigma","lifecycle","killzone","concfile","everyh","ndrug","gamma","ec50",
   		"emin","emax","1/alpha","runmax", "dorfrac","dortime", "dormu","dorsig","outform"};	  
   ,
   	modelID==2.1&&(outform==0),
   	(* p0_,a_,k_,k1_,c0_,gamma_,c50_*)
   	parms = {"model","p0","a","k","k1","c0","gamma","c50","outform"};	   
   ,
   	modelID==2.2&&(outform==0),
   	(* p0_,a_,k_,k1_,c0_,gamma_,c50_, mic_*)
   	parms = {"model","p0","a","k","k1","c0","gamma","c50","mic","outform"};	   
   ,
    modelID==2.3&&(outform==0||outform==1),
    (* ParaFitHoshen[parafile_, concfile_, initN_, PMR_, mu_, sigma_,
	LC_, everyH_, Ndrug_, ec50_, dorfrac_, dortime_, output_] *)
	parms = {"model","parafile","concfile","initn","pmr","mu","sigma","lifecycle","everyh","ndrug","ec50","dorfrac","dortime","outform"};  	 
   ,
   	modelID==2.4&&(outform==0),
   	(*BatchRunHoshen[parafile_String, concfile_String, runsteps_Integer, outfilename_String]*)
	parms = {"model","parafile","concfile","runsteps","outname","outform"};
   ,
  	modelID==2.5&&(outform==0||outform==1||outform==2), (* modified Hoshen et al. Prasitology (2000)// dormant can be observed*) 
   	parms = {"model","initn","pmr","mu","sigma","lifecycle","concfile","everyh","ndrug","ec50","runmax",
   		"dorfrac","dortime","outform"};
   ,
   	modelID==2.6&&(outform==0||outform==1),
   	(*ParaFitD[parafile_, concfile_, initN_, PMR_, mu_, sigma_,
	LC_, KZ_, everyH_, Ndrug_, gamma_, ec50_, emin_, emax_, T_, dorfrac_, dortime_, output_]*)
	parms = {"model","parafile","concfile","initn","pmr","mu","sigma","lifecycle","killzone","everyh","ndrug","gamma","ec50","emin","emax","1/alpha","dorfrac","dortime","outform"};
   ,
    modelID==2.7&&(outform==0),
	parms = {"model","parafile","concfile","runsteps","outname","outform"};
   ,
    modelID==2.71&&(outform==0),
	parms = {"model","parafile","concfile","runsteps","outname","outform"};
   ,
    modelID==2.8&&(outform==0),
   	parms = {"model","parafile","concfile","runsteps","outname","outform"};  	 
   ,
    modelID==2.9&&(outform==0||outform==1),
	parms = {"model","parafile","concfile","initn","pmr","mu","sigma","lifecycle","everyh","ndrug","ec50","dorfrac","dortime","outform"};  	 
     
	
   ];

	parms	

];


ReadPar[inputfile_String] := Module[{txtin, outls, parms,modelID,outform},
   txtin = Import[inputfile, "Table"];
	Print["******INDIVARIA*****"];
   	modelID=ExtractVal[txtin,"model"];
   	outform=ExtractVal[txtin,"outform"];
   	Print["********************"];
   (*Print["Model ID: ",modelID];*)
   
   outls = {};
   
   (* lists of the parameters for the interested functions *)
   Which[
   	modelID==1.0&&outform==0, (*Sexy5*)
   	parms = {"model","initn", "pmr", "mu", "sigma", "lifecycle", "killzone", "concfile", "everyh", 
    "ndrug", "gamma", "ec50", "emin", "emax", "1/alpha", "runmax", "outform"};
   	,
   	modelID==1.0&&outform==1, (*Sexy5*)
   	parms = {"model","initn", "pmr", "mu", "sigma", "lifecycle", "killzone", "concfile", "everyh", 
    "ndrug", "gamma", "ec50", "emin", "emax", "1/alpha", "runmax", "outform"};
   	,
   	modelID==1.0&&outform==2, (*Sexy5*)
   	parms = {"model","initn", "pmr", "mu", "sigma", "lifecycle", "killzone", "concfile", "everyh", 
    "ndrug", "gamma", "ec50", "emin", "emax", "1/alpha", "runmax", "outform"};
   	,
   	modelID==1.0&&outform==3, (*RunSexy5*)
   	parms = {"model","initn", "pmr", "mu", "sigma", "lifecycle", "killzone", "concfile", "everyh", 
    "ndrug", "gamma", "ec50", "emin", "emax", "1/alpha", "lim", "outform"};
   	,
   	modelID==1.1&&(outform==0||outform==1), (*ParaFit*)
	parms = {"model","parafile","concfile","initn","pmr","mu","sigma","lifecycle","killzone","everyh","ndrug","gamma"
		,"ec50","emin","emax","1/alpha","outform"};
	,
	modelID==1.2&&outform==0, (*BatchRun*)   
    (*BatchRun[parafile_String, concfile_String, runsteps_Integer, outfilename_String]*)
    parms = {"model","parafile","concfile","runsteps","outname","outform"};
   ,
	modelID==1.21&&outform==0, (*BatchRun2*)   
    parms = {"model","parafile","concfile","runsteps","outname","outform"};
   ,
   	modelID==1.3&&(outform==0||outform==1||outform==2), (* Hoshen et al. Prasitology (2000)*) 
   	parms = {"model","initn","pmr","mu","sigma","lifecycle","concfile","everyh","ndrug","ec50","runmax",
   		"dorfrac","dortime","outform"};
   ,
   	modelID==1.4&&(outform==0||outform==1||outform==2),  (*NJW model*)
   	parms = {"model","initn","lifecycle","mu","sigma","pmr","endage","runmax","outform"};	   
   ,
   	modelID==1.5&&(outform==0),  (*AMG the original model*)
   	(* lambda_,mu_,beta_,alpha_,r_,d_,initx_,inity_,inits_,runT_Integer,dt_,outform_Integer *)
   	parms = {"model","lambda","mu","beta","alpha","r","d","initx","inity","inits","runsteps","stepsize","outform"};	   
   ,
    modelID==1.6&&(outform==0),  (* Kwaitkowski et al. 2 stages *)
    (* initxy_List, lsd_List, r_, lsPS_List, runT_,outform_ *)
   	parms = {"model","initxy","initd","r","ps","runsteps","outform"};
   ,
   	modelID==1.7&&(outform==0),  (* Kwaitkowski et al. 4 stages *)
   	(* initx_List, r_, h_, lsd_List, lss_List, runT_,outform_*)
   	parms = {"model","initx","r","h","initd","inits","runsteps","outform"};	   
   ,
   	modelID==1.8&&(outform==0),  (* Gravenor's model (Estimating sequestered parasites) *)
   	(* Sques[initn_List,lambda_List,mu_List,r_,runT_,dt_,outform_] *)
   	parms = {"model","initn","lambdals","muls","r","runsteps","stepsize","outform"};	  
   ,
   	modelID==1.9&&(outform==0||outform==1||outform==2),  (* Saralamba dormancy model 1 *)
   	(* initN_, PMR_, mu_, sigma_, hours_Integer, KZ_List, datafile_, 
   everyH_, Ndrug_, gamma_List, ec50_List, emin_List, emax_List, T_, runMax_, 
   dorfrac_, dortime_, outfile_ *)
   	parms = {"model","initn","pmr","mu","sigma","lifecycle","killzone","concfile","everyh","ndrug","gamma","ec50"
   	,"emin","emax","1/alpha","runmax","dorfrac","dortime","outform"};
   ,
   	modelID==2.0&&(outform==0||outform==1||outform==2),   (* Saralamba dormancy model 2 *)
   	parms = {"model","initn","pmr","mu","sigma","lifecycle","killzone","concfile","everyh","ndrug","gamma","ec50",
   		"emin","emax","1/alpha","runmax", "dorfrac","dortime", "dormu","dorsig","outform"};	  	  	 
   ,
   	modelID==2.1&&(outform==0),  (* SimpsonMef1 model *)
   	parms = {"model","p0","a","k","k1","c0","gamma","c50","outform"};	  
   ,
   	modelID==2.2&&(outform==0),  (* (* SimpsonMef2 model *) *)
    parms = {"model","p0","a","k","k1","c0","gamma","c50","mic","outform"};	   
   ,
    modelID==2.3&&(outform==0||outform==1),
    (* ParaFitHoshen[parafile_, concfile_, initN_, PMR_, mu_, sigma_,
	LC_, everyH_, Ndrug_, ec50_, dorfrac_, dortime_, output_] *)
	parms = {"model","parafile","concfile","initn","pmr","mu","sigma","lifecycle","everyh","ndrug","ec50","dorfrac","dortime","outform"};  	 
   ,
   	modelID==2.4&&(outform==0),
   	parms = {"model","parafile","concfile","runsteps","outname","outform"};  	 
   ,
  	modelID==2.5&&(outform==0||outform==1||outform==2), (* modified Hoshen et al. Prasitology (2000)// dormant can be observed*) 
   	parms = {"model","initn","pmr","mu","sigma","lifecycle","concfile","everyh","ndrug","ec50","runmax",
   		"dorfrac","dortime","outform"};
   ,
    modelID==2.6&&(outform==0||outform==1),
   	(*ParaFitD[parafile_, concfile_, initN_, PMR_, mu_, sigma_,
	LC_, KZ_, everyH_, Ndrug_, gamma_, ec50_, emin_, emax_, T_, dorfrac_, dortime_, output_]*)
	parms = {"model","parafile","concfile","initn","pmr","mu","sigma","lifecycle","killzone","everyh","ndrug","gamma","ec50","emin","emax","1/alpha","dorfrac","dortime","outform"};
   ,
    modelID==2.7&&(outform==0),
   	parms = {"model","parafile","concfile","runsteps","outname","outform"};
   ,
    modelID==2.71&&(outform==0),
   	parms = {"model","parafile","concfile","runsteps","outname","outform"};
   ,
   	modelID==2.8&&(outform==0),
   	parms = {"model","parafile","concfile","runsteps","outname","outform"};  	 
   ,
    modelID==2.9&&(outform==0||outform==1),
	parms = {"model","parafile","concfile","initn","pmr","mu","sigma","lifecycle","everyh","ndrug","ec50","dorfrac","dortime","outform"};  	 
    
   
   ];
       
   outls = ExtractVal[txtin, #] & /@ parms;
   
   (*If[Length@outls!=Length@parms,
   	Print["Please check your input file. Some values are missing."];
   	];
   *)
   
   DeleteCases[outls,Null]
   ];



(** the main engine  **)
IDVL[inputfile_String]:=Module[{par,modelID,output,outform,runmax},

	output=Null;
	par=ReadPar[inputfile];
	
	(* modelID *)
	modelID=First[par];
	outform=Last[par];
	
	Which[		
		modelID==1.0&&outform==0,
		(*parms={"model","initn", "pmr", "mu", "sigma", "lifecycle", "killzone", "concfile", "everyh", 
    	"ndrug", "gamma", "ec50", "emin", "emax", "1/alpha", "runmax", "outform"};*)
		output= Sexy5[par[[2]],par[[3]],par[[4]],par[[5]],par[[6]],par[[7]],StringJoin@DeleteCases[Characters[ToString@par[[8]]], " "],
		par[[9]],par[[10]],par[[11]],par[[12]],par[[13]],par[[14]],par[[15]],par[[16]],par[[17]]];
	,
		modelID==1.0&&outform==1,
		(*parms={"model","initn", "pmr", "mu", "sigma", "lifecycle", "killzone", "concfile", "everyh", 
    	"ndrug", "gamma", "ec50", "emin", "emax", "1/alpha", "runmax", "outform"};*)
		output= Sexy5[par[[2]],par[[3]],par[[4]],par[[5]],par[[6]],par[[7]],StringJoin@DeleteCases[Characters[ToString@par[[8]]], " "],
		par[[9]],par[[10]],par[[11]],par[[12]],par[[13]],par[[14]],par[[15]],par[[16]],1];
	,	
		modelID==1.0&&outform==2,
		(*parms={"model","initn", "pmr", "mu", "sigma", "lifecycle", "killzone", "concfile", "everyh", 
    	"ndrug", "gamma", "ec50", "emin", "emax", "1/alpha", "runmax", "outform"};*)
		output= Sexy5[par[[2]],par[[3]],par[[4]],par[[5]],par[[6]],par[[7]],StringJoin@DeleteCases[Characters[ToString@par[[8]]], " "],
		par[[9]],par[[10]],par[[11]],par[[12]],par[[13]],par[[14]],par[[15]],par[[16]],2];
	,
		modelID==1.0&&outform==3, 
		(* parms = {"model","initn", "pmr", "mu", "sigma", "lifecycle", "killzone", "concfile", "everyh", 
    	"ndrug", "gamma", "ec50", "emin", "emax", "1/alpha", "lim","output"}; *)
		runmax=1000;
		output=RunSexy5[par[[2]],par[[3]],par[[4]],par[[5]],par[[6]],par[[7]],StringJoin@DeleteCases[Characters[ToString@par[[8]]], " "],par[[9]],
		par[[10]],par[[11]],par[[12]],par[[13]],par[[14]],par[[15]],runmax,par[[16]]];
	,
		modelID==1.1&&(outform==0||outform==1),
		(* parms = {"model","parafile","concfile","initn","pmr","mu","sigma","lifecycle","killzone","everyh","ndrug","gamma"
		,"ec50","emin","emax","1/alpha","outform"} *)
		output=ParaFit[StringJoin@DeleteCases[Characters[ToString@par[[2]]], " "],StringJoin@DeleteCases[Characters[ToString@par[[3]]], " "],par[[4]],par[[5]],par[[6]],par[[7]],par[[8]],
			par[[9]],par[[10]],par[[11]],par[[12]],par[[13]],par[[14]],par[[15]],par[[16]],par[[17]]];
	,
		modelID==1.2&&outform==0,
		(* parms = {"model","parafile","concfile","runsteps","outname","outform"}; *)
		output=BatchRun[StringJoin@DeleteCases[Characters[ToString@par[[2]]], " "],StringJoin@DeleteCases[Characters[ToString@par[[3]]], " "],par[[4]],ToString@par[[5]]];
	,
		modelID==1.21&&outform==0,
		(* parms = {"model","parafile","concfile","runsteps","outname","outform"}; *)
		output=BatchRun2[StringJoin@DeleteCases[Characters[ToString@par[[2]]], " "],StringJoin@DeleteCases[Characters[ToString@par[[3]]], " "],par[[4]],ToString@par[[5]]];
	,
		modelID==1.3&&(outform==0||outform==1||outform==2),
		(* parms = {"model","initn","pmr","mu","sigma","lifecycle","concfile","everyh","ndrug","ec50","runmax","dorfrac","dortime","outform"} *)	
		(*Hoshen[initN_, PMR_Integer, mu_, sigma_, hours_Integer, datafile_, everyH_Integer, 
		Ndrug_Integer, ec50_, runMax_Integer, dorfrac_, dortime_, outfile_]*)
		output=Hoshen[par[[2]],par[[3]],par[[4]],par[[5]],par[[6]],StringJoin@DeleteCases[Characters[ToString@par[[7]]], " "],par[[8]],par[[9]],par[[10]],par[[11]],par[[12]],par[[13]],par[[14]]];
	,
		modelID==1.4&&(outform==0||outform==1||outform==2),
		(*parms = {"model","initn","lifecycle","mu","sigma","pmr","endage","runmax","outform"}*)
		output=NJWhite[par[[2]],par[[3]],par[[4]],par[[5]],par[[6]],par[[7]],par[[8]],par[[9]]];
	,
		modelID==1.5&&(outform==0),
		(* parms = {"model","lambda","mu","beta","alpha","r","d","initx","inity","inits","runsteps","stepsize","outform"}*)
		output=AMG[par[[2]],par[[3]],par[[4]],par[[5]],par[[6]],par[[7]],par[[8]],par[[9]],par[[10]],par[[11]],par[[12]],par[[13]]];
	,
		modelID==1.6&&(outform==0),
    	(* initxy_List, lsd_List, r_, lsPS_List, runT_,outform_ *)
   		(*parms = {"model","initxy","initd","r","ps","runsteps","outform"};*)
   		output = DominicNowak1[par[[2]],par[[3]],par[[4]],par[[5]],par[[6]],par[[7]]];		
	,
		modelID==1.7&&(outform==0),
		(* parms = {"model","initx","r","h","initd","inits","runsteps","outform"}*)
		output=DominicNowak2[par[[2]], par[[3]],par[[4]],par[[5]],par[[6]],par[[7]],par[[8]]];
	,
		modelID==1.8&&(outform==0),
		(* parms = {"model","initn","lambdals","muls","r","runsteps","stepsize","outform"};*)
		output=Seques[par[[2]],par[[3]],par[[4]],par[[5]],par[[6]],par[[7]],par[[8]]];
	,	
		modelID==1.9&&(outform==0||outform==1||outform==2),
   	   	(* parms = {"model","initn","pmr","mu","sigma","lifecycle","killzone","concfile","everyh","ndrug","gamma","ec50"
   		,"emin","emax","1/alpha","runmax","dorfrac","dortime","outform"}; *)
   		output=Sexy6[par[[2]],par[[3]],par[[4]], par[[5]], par[[6]],par[[7]],StringJoin@DeleteCases[Characters[ToString@par[[8]]], " "],par[[9]],par[[10]], 
   		par[[11]],par[[12]],par[[13]],par[[14]],par[[15]],par[[16]],par[[17]],par[[18]],outform];
	,	
		modelID==2.0&&(outform==0||outform==1||outform==2),
   	   	(*parms = {"model","initn","pmr","mu","sigma","lifecycle","killzone","concfile","everyh","ndrug","gamma","ec50",
   		"emin","emax","1/alpha","runmax", "dorfrac","dortime", "dormu","dorsig","outform"};*)
   		output=Sexy7[par[[2]],par[[3]],par[[4]], par[[5]], par[[6]],par[[7]],StringJoin@DeleteCases[Characters[ToString@par[[8]]], " "],par[[9]],par[[10]], 
   		par[[11]],par[[12]],par[[13]],par[[14]],par[[15]],par[[16]],par[[17]],par[[18]],par[[19]],par[[20]],outform];

	,
		modelID==2.1&&outform==0,
		(*parms = {"model","p0","a","k","k1","c0","gamma","c50","outform"};	 *)
		output=SimpsonMef1[par[[2]],par[[3]],par[[4]],par[[5]],par[[6]],par[[7]],par[[8]]];		
	,
		modelID==2.2&&outform==0,
		(*parms = {"model","p0","a","k","k1","c0","gamma","c50","mic","outform"};	 *)
		output=SimpsonMef2[par[[2]],par[[3]],par[[4]],par[[5]],par[[6]],par[[7]],par[[8]],par[[9]]];		
	,
		modelID==2.3&&(outform==0||outform==1),
		(* 	parms = {"model","parafile","concfile","initn","pmr","mu","sigma","lifecycle","everyh","ndrug","ec50","dorfrac","dortime","outform"}; *)
		output=ParaFitHoshen[StringJoin@DeleteCases[Characters[ToString@par[[2]]], " "],StringJoin@DeleteCases[Characters[ToString@par[[3]]], " "],par[[4]],
		par[[5]],par[[6]],par[[7]],par[[8]],par[[9]],par[[10]],par[[11]],par[[12]],par[[13]],par[[14]]];
	,
		modelID==2.4&&(outform==0),
		(*BatchRunHoshen[parafile_String, concfile_String, runsteps_Integer, outfilename_String]*)
		output=BatchRunHoshen[StringJoin@DeleteCases[Characters[ToString@par[[2]]], " "],StringJoin@DeleteCases[Characters[ToString@par[[3]]], " "],par[[4]],
		ToString@par[[5]]];
	,
		modelID==2.5&&(outform==0||outform==1||outform==2),
		(* parms = {"model","initn","pmr","mu","sigma","lifecycle","concfile","everyh","ndrug","ec50","runmax","dorfrac","dortime","outform"} *)	
		(*Hoshen[initN_, PMR_Integer, mu_, sigma_, hours_Integer, datafile_, everyH_Integer, 
		Ndrug_Integer, ec50_, runMax_Integer, dorfrac_, dortime_, outfile_]*)
		output=HoshenNew[par[[2]],par[[3]],par[[4]],par[[5]],par[[6]],StringJoin@DeleteCases[Characters[ToString@par[[7]]], " "],par[[8]],par[[9]],par[[10]],par[[11]],par[[12]],par[[13]],par[[14]]];
	,
		modelID==2.6&&(outform==0||outform==1),
		(*	parms = {"model","parafile","concfile","initn","pmr","mu","sigma","lifecycle","killzone","everyh","ndrug","gamma","ec50","emin","emax","1/alpha","dorfrac","dortime","outform"};*)		
		output=ParaFitD[StringJoin@DeleteCases[Characters[ToString@par[[2]]], " "],StringJoin@DeleteCases[Characters[ToString@par[[3]]], " "],par[[4]],par[[5]],par[[6]],par[[7]],par[[8]],par[[9]],par[[10]],par[[11]],par[[12]],par[[13]],par[[14]],par[[15]],par[[16]],par[[17]],par[[18]],par[[19]]];
	,
		modelID==2.7&&(outform==0),
		output=output=BatchRunD[StringJoin@DeleteCases[Characters[ToString@par[[2]]], " "],StringJoin@DeleteCases[Characters[ToString@par[[3]]], " "],par[[4]],
		ToString@par[[5]]];
	,
		modelID==2.71&&(outform==0),
		output=output=BatchRunDfixedRI[StringJoin@DeleteCases[Characters[ToString@par[[2]]], " "],StringJoin@DeleteCases[Characters[ToString@par[[3]]], " "],par[[4]],
		ToString@par[[5]]];
	,
		modelID==2.8&&(outform==0),
		(*BatchRunHoshen[parafile_String, concfile_String, runsteps_Integer, outfilename_String]*)
		output=BatchRunHoshenNew[StringJoin@DeleteCases[Characters[ToString@par[[2]]], " "],StringJoin@DeleteCases[Characters[ToString@par[[3]]], " "],par[[4]],
		ToString@par[[5]]];
	,
		modelID==2.9&&(outform==0||outform==1),
		(* 	parms = {"model","parafile","concfile","initn","pmr","mu","sigma","lifecycle","everyh","ndrug","ec50","dorfrac","dortime","outform"}; *)
		output=ParaFitHoshenNew[StringJoin@DeleteCases[Characters[ToString@par[[2]]], " "],StringJoin@DeleteCases[Characters[ToString@par[[3]]], " "],par[[4]],
		par[[5]],par[[6]],par[[7]],par[[8]],par[[9]],par[[10]],par[[11]],par[[12]],par[[13]],par[[14]]];


	];

	output
];


End[]

(*SetAttributes[{IDVL,RI},{Locked,ReadProtected,Protected}];*)

EndPackage[]

