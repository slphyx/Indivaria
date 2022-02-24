(* Mathematica package *)

(* Author: Sompob Saralamba *)

(* Mathematical Modelling Team, Mahidol-Oxford Tropical Medicine Research Unit *)


BeginPackage["Indivaria`PfHRP2`",{"Indivaria`BasicFunc`"}]

(* Exported symbols added here with SymbolName::usage *) 
PfHRP2Conc::usage="PfHRP2Conc[hrp2perRBC,initN,PMR,mu,sigma,LC,KZ,concdat,everyH,Ndrug,gamma,ec50,emin,emax,T,runMax]"


Begin["`Private`"]


Ki[T_Real, concls_List] := 
  Module[{i, lsk}, 
   lsk = Table[{1/T Log[100./(100. - concls[[i, 3]])], 
      1/T Log[100./(100. - concls[[i, 4]])], 
      1/T Log[100./(100. - concls[[i, 5]])]}, {i, 1, Length[concls]}];
   lsk];

(*decay function*)
Fdecay[ages_, lst_, attime_, lsk_, stages_] := Module[{tmp, i},
   i = ages;
   tmp = lst[[i]];
   If[stages[[i]] == 1, tmp = lst[[i]]*Exp[-lsk[[attime, 1]]]];
   If[stages[[i]] == 2, tmp = lst[[i]]*Exp[-lsk[[attime, 2]]]];
   If[stages[[i]] == 3, tmp = lst[[i]]*Exp[-lsk[[attime, 3]]]];
   tmp];

WhichRTS[lst_List, KZ_List] := Module[{i, tmp},
   tmp = Table[0, {Length[lst]}];
   For[i = 1, i <= Length[lst], i = i + 1, 
    If[IntervalMemberQ[Interval[KZ[[1]]], i] == True, 
     tmp = ReplacePart[tmp, i -> 1];];
    If[IntervalMemberQ[Interval[KZ[[2]]], i] == True, 
     tmp = ReplacePart[tmp, i -> 2];];
    If[IntervalMemberQ[Interval[KZ[[3]]], i] == True, 
     tmp = ReplacePart[tmp, i -> 3];];];
   tmp];

(* gives the list of the number of ruptured schizonts over time *)
RupturedSchizonts[initN_, PMR_Integer, mu_, sigma_, hours_Integer, 
   KZ_List, concdat_List, everyH_Integer, Ndrug_Integer, gamma_List, 
   ec50_List, emin_List, emax_List, T_Real, runMax_Integer] := 
  Module[{runs = True, stages, concls, i, j, lst, lsk, output, 
    efffn,efffnparms},
  
   (* 
   (*template of the drug and its effects*)
   xm = concdat[[1]];
   ym = concdat[[2]];
   ke = concdat[[3]];
   *)
   (* EffConc[Ndrug,everyH,efffn,efffnparms,runtime]
   
   (* generates the list of the drug concentration and its efficacy over time following the dose regimen  
   ndrug = number of the drugs; everh = time for taking next dose; fn = the efficacy function; fnparms = the list of parameters of fn 

	efffn = 0 , efffnparms = {"dataname",gamma_List,ec50_List,emin_List,emax_List}
	efffn = 1 , efffnparms = {gamma_List,ec50_List,emin_List,emax_List, concparms_List} , concparms = {xm,ym,ke}
    efffn = 2 , efffnparms = {gamma,ec50,emin,emax, concparm_List}, concparms = {xm,ym,ke}
*)
  *) 
   (*concls = ConcMod[xm, ym, ke, everyH, Ndrug, gamma, ec50, emin, emax];*)
   efffn = 1;
   efffnparms = {gamma,ec50,emin,emax, concdat};
   concls = EffConc[Ndrug,everyH,efffn,efffnparms,runMax];
   
   (*k_i (t)*)
   lsk = Ki[T, concls];
   (*initial parasite load*)
   
   output=Reap[
   	
   lst = DistributeN[initN, hours, mu, sigma];
   Sow[lst,dist];
   Sow[Last@lst,schi];
   
   stages = WhichRTS[lst, KZ];
   
   (*AppendTo[output, lst];*)
   
   i = 0;
   While[runs == True && i < runMax,
    i = i + 1;
    
    (*Parasites are growing.Feed them!*)
    lst = Shiftonehour[lst, PMR];
    Sow[lst,dist];
    Sow[Last@lst,schi];
    
    (*a time to kill.*)
    For[j = 1, j <= Length[lst], j = j + 1, 
     lst = ReplacePart[lst, j -> Fdecay[j, lst, i, lsk, stages]];
     ];
    
    (*adding a point for ploting*)
    (*AppendTo[output, lst];*)
    
    ];
    ,{dist,schi}][[2]];
    
   (*
   (**non sequestered parasites**)
   junk = 
    LsDot[#, Table[PRingFunc[i, 11, 14], {i, 1, hours}]] & /@ output;
   onlyring = 
    Table[{i - 1, Log10[junk[[i]] // Total]}, {i, 1, Length[junk]}];
   
   (*total of parasites over time*)
   tot = Table[{i - 1, Log10@Total[output[[i]]]}, {i, 1, 
      Length[output]}];
   
   {onlyring, tot}
   
   *)
   output
 ];

PfHRP2Conc[hrp2perRBC_,initN_, PMR_Integer, mu_, sigma_, hours_Integer, 
   KZ_List, concdat_List, everyH_Integer, Ndrug_Integer, gamma_List, 
   ec50_List, emin_List, emax_List, T_Real, runMax_Integer]:=Module[{schils,temp,paradist},
   
   temp = RupturedSchizonts[initN, PMR, mu, sigma, hours, KZ, concdat, everyH, Ndrug, gamma, 
    ec50, emin, emax, T, runMax];
   
   schils = temp[[2,1]];
   paradist = temp[[1,1]];
   
   hrp2perRBC*schils
 ];


End[] (* End Private Context *)

EndPackage[]