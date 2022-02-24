(* Mathematica Package *)

BeginPackage["Indivaria`Game`", { "Indivaria`BasicFunc`"}]
(* Exported symbols added here with SymbolName::usage *)  
Game::usage = "Game[initN, mu, sigma, LC, PMR, gcons, grate, gtriggerT, runtime, outform]"
Game2::usage="Game2[initN, mu, sigma, LC, PMR, gcons, gtriggerT, gages, runtime, outform]" 
Game3::usage="Game3[initN, mu, sigma, LC, PMR, gcons, gages, runtime, outform]"
Game4::usage="Game4[initN, mu, sigma, LC, PMR, gcons, gages, runtime, outform]"
Game5::usage="Game5[initN, mu, sigma, LC, PMR, gcons, gages, pkpdparms, runtime, outform]"
Game6::usage="Game6[initN, mu, sigma, LC, PMR, gcons, gages, seqparms, pkpdparms, runtime, outform]"
Game7::usage="Game7[initN, mu, sigma, LC, PMR, gcons, gages, seqparms, pkpdparms, runtime, outform]"
Game8::usage="Game8[initN, mu, sigma, LC, PMR, gconstants, gages, immparms, trigger, immuneFunc, runtime, outform]" 
Game9::usage="Game9[initN, mu, sigma, LC, PMR, gconstants, gages, immparms,trigger, immuneFunc, drugAct, pkpdparms, runtime, outform]"
Game10::usage="Game10[initN, mu, sigma, LC, PMR, gconstants, gages, immparms,trigger, immuneFunc, drugAct, pkpdparms, runtime, outform]"
Game11::usage="Game11[initN, mu, sigma, LC, PMR, gconstants, gages, immparms,trigger, immuneFunc, drugAct, pkpdparms, runtime, outform]"


Begin["`Private`"] (* Begin Private Context *) 


 (* the main model *)
Game[initN_, mu_, sigma_, LC_Integer, PMR_, gcons_List, gtriggerT_Integer, runtime_, outform_] := 
	Module[{ls, i,grate, tmp, gls, ngam, nring, nring2, k, ygls, ogls, obsex, obasex},
   
   
   
   (*initial parasites after coming out from the liver*)
   tmp = DistributeN[initN, LC, mu, sigma]; 
   ls = {tmp}; 
   gls = {0};
   ygls = {0};
   ogls = {0};
   ngam = 0;
   nring = {};

   nring2 = {0}; (* gametocytes look like rings *)
   
  (*gcons {ratio become gametocytes , ration look like rings} *)
   grate = gcons[[1]];
  (* observation ratio (look like ring ratio) *)
   k = gcons[[2]];
   
  {{ls},{nring2},{ygls},{ogls}}= Reap[
   Do[
     tmp = Shiftonehour[tmp, PMR];        
	 (* become gametocytes after gtrigger+LC *)
     If[i >= gtriggerT+LC,
     	ngam = grate*First@tmp;
     	tmp = ReplacePart[tmp,1->(1-grate)*First@tmp];
     ];
          
     (*AppendTo[ls, tmp];*)
     Sow[tmp,0];
     PrependTo[gls,ngam];
     
     (*
     (* look like rings *)
     AppendTo[nring2,CountObservedGame[gls,{1,48},k]];
	 (* young observed gametocytes *)
	 AppendTo[ygls,CountObservedGame[gls,{1,48},(1-k)]];
	 (* old observed gametocytes *)
	 AppendTo[ogls,CountObservedGame[gls,{144,240},1]]
     *)
	 Sow[CountObservedGame[gls,{1,48},k],1];
	 Sow[CountObservedGame[gls,{1,48},(1-k)],2];
	 Sow[CountObservedGame[gls,{144,240},1],3];
	 
	 
     , {i, 1, runtime}]
     ,{0,1,2,3}
     ][[2]];

   (* observed rings (no gametocytes) *)
   nring = CountRing[#,{11,14},1]&/@ls;
   
   obasex = nring + nring2;
   obsex = ygls + ogls;

   Which[
   	(* show age distribution over time*)
   	outform==0,ls  
    (* show the last age distribution*)
   ,outform==1,tmp
    (* show the number of gametocyte over time*)
   ,outform==2,gls
    (* big results! *)
   ,outform==3,{ls,gls,nring,obasex,obsex}
    
   ]

];

(* k percentage that gametocytes look like ring *)
CountObservedGame[gls_List,AgeRange_List,k_]:=Module[{glsLength,obgam},
	glsLength = gls//Length;
	
	If[(glsLength<=AgeRange[[2]]) && (AgeRange[[1]] < glsLength),
		obgam = k*Total[gls[[AgeRange[[1]];;glsLength]]];
	,	
		If[(glsLength > AgeRange[[2]]) && (AgeRange[[1]] < glsLength),
		obgam = k*Total[gls[[AgeRange[[1]];;AgeRange[[2]]]]];
		,
		obgam = 0;
		];	
	];
	
	
	obgam
];


(* Game2 has the input for the age ranges of the gametocytes *)
Game2[initN_, mu_, sigma_, LC_Integer, PMR_, gcons_List, gtriggerT_Integer, gages_List, runtime_, outform_] := 
	Module[{ls, i,grate, tmp, gls, ngam, nring, nring2, k, ygls, ogls,infgls,obsex, obasex},
   
   (*grate {ratio become gametocytes, ration look like rings} *)
   (*gages {ag1,ag2,ag3,ag4,ag5} gametocyte age ranges*)
   (*initial parasites after coming out from the liver*)
   tmp = DistributeN[initN, LC, mu, sigma]; 
   ls = {tmp}; 
   gls = {0};
   ygls = {0};
   ogls = {0};
   infgls={0};
   ngam = 0;
   nring = {};

   nring2 = {0}; (* gametocytes look like rings *)
   
   grate = gcons[[1]];
  (* observation ratio (look like ring ratio) *)
   k = gcons[[2]];
   
   Do[
     tmp = Shiftonehour[tmp, PMR];        
	 (* become gametocytes after gtrigger+LC *)
     If[i >= gtriggerT+LC,
     ngam = grate*First@tmp;
     tmp = ReplacePart[tmp,1->(1-grate)*First@tmp];
     ];
          
     AppendTo[ls, tmp];
     
     PrependTo[gls,ngam];
     (* reset ngam *)
     ngam = 0;
     
     (* look like rings *)
     AppendTo[nring2,CountObservedGame[gls,{gages[[1]],gages[[2]]},k]];
	 (* young observed gametocytes *)
	 AppendTo[ygls,CountObservedGame[gls,{gages[[1]],gages[[2]]},(1-k)]];
	 (* old observed gametocytes *)
	 AppendTo[ogls,CountObservedGame[gls,{gages[[3]],gages[[4]]},1]];
	 (* infectiousness gametocytes *)
	 AppendTo[infgls,CountObservedGame[gls,{gages[[4]],gtriggerT+LC+240},1]];
	 
     , {i, 1, runtime}];

   (* observed rings (no gametocytes) *)
   nring = CountRing[#,{11,14},1]&/@ls;
   
   obasex = nring + nring2;
   obsex = ygls + ogls;
      
   Which[
   	(* show age distribution over time*)
   	outform==0,ls  
    (* show the last age distribution*)
   ,outform==1,tmp
    (* show the number of gametocyte over time*)
   ,outform==2,gls
    (* big results! *)
   ,outform==3,{ls,gls,nring,obasex,obsex,infgls}
    
   ]

];

(* Game3 has the input for the age ranges of the gametocytes, no trigger time  *)
Game3[initN_, mu_, sigma_, LC_Integer, PMR_, gcons_List, gages_List, runtime_, outform_] := 
	Module[{ls, i,grate, tmp, gls, ngam, nring, nring2, k, ygls, 
		ogls,infgls,obsex, obasex,threshold,thresholdT, ag1, ag2, ag3, ag4,
		gotit
		
		},
   
   (*grate {ratio become gametocytes, ration look like rings} *)
   (*gages {ag1,ag2,ag3,ag4} gametocyte age ranges
   I:0-ag1,II:ag1-ag2,III:ag2-ag3,IV:ag3-ag4,V:ag4-ag5 
   ###life cycle of gametocytes = ag5  ###
   
   obg1 = look like rings (0-obg1)
   obg1 - obg2 = sequestered
   obg3 - ag5 = observed gams
    *)
   ag1=gages[[1]]; 
   ag2=gages[[2]];  
   ag3=gages[[3]];
   ag4=gages[[4]];
   
   (*initial parasites after coming out from the liver*)
   tmp = DistributeN[initN, LC, mu, sigma]; 
   ls = {tmp}; 

   gls = {0};
   ygls = {0};
   ogls = {0};
   
   thresholdT = 0;

   infgls={0};
   ngam = 0;
   nring = {};

   nring2 = {0}; (* gametocytes look like rings *)
   
   (* gcons = List of the constants i.e ,threshold (number of asexual) for switching to gams 
      gcons = {% of being gams/100 (1/h), threshold, look like rings ratio}
   *)
   (* % of becoming gams/100 (1/h) *)   
   grate = gcons[[1]];
  (* observation ratio (look like ring ratio) *)
   k = gcons[[2]];
   (* threshold for switching to gams *)
   threshold = gcons[[3]];
   
   gotit = 0;
   
   Do[
     tmp = Shiftonehour[tmp, PMR];
             
	 (* become gametocytes if npara > threshold *)
     If[Total@tmp > threshold && gotit==0, 
     	thresholdT = i;
     	gotit = 1;
     ];
	 (* gametocytes are generated LC hours after  Npara reaches the threshold *) 
     If[i >= thresholdT+LC,
     	ngam = grate*First@tmp;
     	tmp = ReplacePart[tmp,1->(1-grate)*First@tmp];
     	,
     	ngam = 0;
     ];
          
     AppendTo[ls, tmp];     
     
     PrependTo[gls,ngam];
     (* reset ngam *)
     ngam = 0;
     
     (* look like rings *)
     AppendTo[nring2,CountObservedGame[gls,{1,ag1},k]];
	 (* young observed gametocytes *)
	 AppendTo[ygls,CountObservedGame[gls,{1,ag1},(1-k)]];
	 (* old observed gametocytes *)
	 AppendTo[ogls,CountObservedGame[gls,{ag2+1,ag3},1]];
	 (* infectiousness gametocytes *)
	 AppendTo[infgls,CountObservedGame[gls,{ag3+1,ag4},1]];
	 
     , {i, 1, runtime}];

   (* observed rings (no gametocytes) *)
   nring = CountRing[#,{11,14},1]&/@ls;
   
   obasex = nring + nring2;
   obsex = ygls + ogls;
      
 (*** output forms ***)
   Which[
   	(* show age distribution over time*)
   	outform==0,ls  
    (* show the last age distribution*)
   ,outform==1,tmp
    (* show the number of gametocyte over time*)
   ,outform==2,gls
    (* big results! *)
   ,outform==3,{ls,gls,nring,obasex,obsex,infgls,thresholdT}
    (* threshold time *)
   ,outform==4,thresholdT 
    
   ]

];


(* Game4 has the input for the age ranges of the gametocytes, no trigger time, modified game3, no killing effect on gametocyte  *)
Game4[initN_, mu_, sigma_, LC_Integer, PMR_, gcons_List, gages_List, runtime_, outform_] := 
	Module[{ls, i,grate, tmp, gls, ngam, nring, nring2, k, ygls, 
		ogls,infgls,obsex, obasex,threshold,thresholdT, ag1, ag2, ag3, ag4	
		},
   
   (*grate {ratio become gametocytes, ration look like rings} *)
   (*gages {ag1,ag2,ag3,ag4} gametocyte age ranges
   I:0-ag1,II:ag1-ag2,III:ag2-ag3,IV:ag3-ag4,V:ag4-ag5 
   ###life cycle of gametocytes = ag5  ###
   
   obg1 = look like rings (0-obg1)
   obg1 - obg2 = sequestered
   obg3 - ag5 = observed gams
    *)
   ag1=gages[[1]]; 
   ag2=gages[[2]];  
   ag3=gages[[3]];
   ag4=gages[[4]];
   
   (*initial parasites after coming out from the liver*)
   tmp = DistributeN[initN, LC, mu, sigma]; 
   ls = {tmp}; 

   gls = {0};
   ygls = {0};
   ogls = {0};
   
   thresholdT = {0};

   infgls={0};
   ngam = 0;
   nring = {};

   nring2 = {0}; (* gametocytes look like rings *)
   
   (* gcons = List of the constants i.e ,threshold (number of asexual) for switching to gams 
      gcons = {% of being gams/100 (1/h), threshold, look like rings ratio}
   *)
   (* % of becoming gams/100 (1/h) *)   
   grate = gcons[[1]];
  (* observation ratio (look like ring ratio) *)
   k = gcons[[2]];
   (* threshold for switching to gams *)
   threshold = gcons[[3]];
   
   
   Do[
     tmp = Shiftonehour[tmp, PMR];
             
	 (* become gametocytes if npara > threshold *)
     If[Total@tmp > threshold ,      	
     	ngam = grate*(First@tmp)/PMR;
     	tmp = ReplacePart[tmp,1->(1-grate)*First@tmp];
     	AppendTo[thresholdT,i];
     	,
     	ngam = 0;
     ];
          
     AppendTo[ls, tmp];     
     
     PrependTo[gls,ngam];
     (* reset ngam *)
     ngam = 0;
     
     (* look like rings *)
     AppendTo[nring2,CountObservedGame[gls,{1,ag1},k]];
	 (* young observed gametocytes *)
	 AppendTo[ygls,CountObservedGame[gls,{1,ag1},(1-k)]];
	 (* old observed gametocytes *)
	 AppendTo[ogls,CountObservedGame[gls,{ag2+1,ag3},1]];
	 (* infectiousness gametocytes *)
	 AppendTo[infgls,CountObservedGame[gls,{ag3+1,ag4},1]];
	 
     , {i, 1, runtime}];

   (* observed rings (no gametocytes) *)
   nring = CountRing[#,{11,14},1]&/@ls;
   
   obasex = nring + nring2;
   obsex = ygls + ogls;
      
 (*** output forms ***)
   Which[
   	(* show age distribution over time*)
   	outform==0,ls  
    (* show the last age distribution*)
   ,outform==1,tmp
    (* show the number of gametocyte over time*)
   ,outform==2,gls
    (* big results! *)
   ,outform==3,{ls,gls,nring,obasex,obsex,infgls,thresholdT}
    (* threshold time *)
   ,outform==4,thresholdT 
    
   ]

];



(* Game5 has the input for the age ranges of the gametocytes, no trigger time, modified game4, add one EC-curve for gametocytes  *)
Game5[initN_, mu_, sigma_, LC_Integer, PMR_, gcons_List, gages_List, pkpdparms_List, runtime_Integer, outform_Integer] := 
	Module[{ls, i,grate, tmp, gls, ngam, nring, nring2, k, ygls, 
		ogls,infgls,obsex, obasex,threshold,thresholdT, ag1, ag2, ag3, ag4,
		gkillzones,ndrug,everyh,takeT,gamma,ec50,emin,emax,xm,ym,ke,drugeff	
		},
   
   (*grate {ratio become gametocytes, ration look like rings} *)
   (*gages {ag1,ag2,ag3,ag4} gametocyte age ranges
   I:0-ag1,II:ag1-ag2,III:ag2-ag3,IV:ag3-ag4,V:ag4-ag5 
   ###life cycle of gametocytes = ag5  ###
   
   obg1 = look like rings (0-obg1)
   obg1 - obg2 = sequestered
   obg3 - ag5 = observed gams
    *)
   ag1=gages[[1]]; 
   ag2=gages[[2]];  
   ag3=gages[[3]];
   ag4=gages[[4]];
   
   (*initial parasites after coming out from the liver*)
   tmp = DistributeN[initN, LC, mu, sigma]; 
   ls = {tmp}; 

   gls = {0};
   ygls = {0};
   ogls = {0};
   
   thresholdT = {0};

   infgls={0};
   ngam = 0;
   nring = {};

   nring2 = {0}; (* gametocytes look like rings *)
   
   (* gcons = List of the constants i.e ,threshold (number of asexual) for switching to gams 
      gcons = {% of being gams/100 (1/h), threshold, look like rings ratio}
   *)
   (* % of becoming gams/100 (1/h) *)   
   grate = gcons[[1]];
  (* observation ratio (look like ring ratio) *)
   k = gcons[[2]];
   (* threshold for switching to gams *)
   threshold = gcons[[3]];
  
  (* there is only one killzone.
     pkpdparms = {killzones_List,ndrug,everyh,takeT,{gamma,ec50,emin,emax},concparms_List}
     concparms = {xm,ym,ke}
     takeT = time that the patient taking the first dose (0 -> runtime)
   *) 
   
   {gkillzones,ndrug,everyh,takeT,{gamma,ec50,emin,emax},{xm,ym,ke}} = pkpdparms;
   
   (* calculate the concentration and its efficacy from t0 = takeT --> runtime + (takeT-t0) *)
   drugeff = EffConc[ndrug,everyh, 2, {gamma,ec50,emin,emax, {xm,ym,ke}} ,runtime];
   (* change time to be the same as the system *)
   drugeff[[All,1]] = drugeff[[All,1]]+takeT;
   
(*   
   gkillzone = pkpdparms[[1]];
   ndrug = pkpdparms[[2]];
   everyh = pkpdparms[[3]];
   gamma = pkpdparms[[4,1]];
   ec50 = plpdparms[[4,2]];
   emin = pkpdparms[[4,3]];
   emax = pkpdparms[[4,4]];
   {xm,ym,ke} = concparms;
*)
    
   
   Do[
     tmp = Shiftonehour[tmp, PMR];
             
	 (* become gametocytes if npara > threshold *)
     If[Total@tmp > threshold ,      	
     	ngam = grate*(First@tmp);
     	tmp = ReplacePart[tmp,1->(1.0-grate)*First@tmp];
     	AppendTo[thresholdT,i];
     	,
     	ngam = 0;
     ];
     
     (* list ot total parasites; (time) past --> future  *)     
     AppendTo[ls, tmp];     
     
     (* list of total gametocytes; (time) young --> old *)
     PrependTo[gls,ngam];
         
     (*drug action for gametocyte*)
     gls = drugActionGam[gls,gkillzones,drugeff,i];
     
     (* reset ngam *)
     ngam = 0;
     
     (* look like rings *)
     AppendTo[nring2,CountObservedGame[gls,{1,ag1},k]];
	 
	 (* young observed gametocytes *)
	 AppendTo[ygls,CountObservedGame[gls,{1,ag1},(1-k)]];
	
	 (* old observed gametocytes *)
	 AppendTo[ogls,CountObservedGame[gls,{ag2+1,ag3},1]];
	
	 (* infectiousness gametocytes *)
	 AppendTo[infgls,CountObservedGame[gls,{ag3+1,ag4},1]];
	 
     , {i, 1, runtime}];

   (* observed rings (no gametocytes) *)
   nring = CountRing[#,{11,14},1]&/@ls;
   
   obasex = nring + nring2;
   obsex = ygls + ogls;
      
 (*** output forms ***)
   Which[
   	(* show age distribution over time*)
   	outform==0,ls  
    (* show the last age distribution*)
   ,outform==1,tmp
    (* show the number of gametocyte over time*)
   ,outform==2,gls
    (* big results! *)
   ,outform==3,{ls,gls,nring,obasex,obsex,infgls,thresholdT,drugeff}
    (* threshold time *)
   ,outform==4,thresholdT 
    (* drug effect *)
   ,outform==5,drugeff 
   ]

];

(*****)

(* Game6 has the input for the age ranges of the gametocytes, no trigger time, modified game5, add one EC-curve for gametocytes 
   
   sequestered sexual + asexual parasites are killed by the immune system
 *)
Game6[initN_, mu_, sigma_, LC_Integer, PMR_, gcons_List, gages_List, seqparms_List,pkpdparms_List, runtime_Integer, outform_Integer] := 
	Module[{ls, i,grate, tmp, gls, ngam, nring, nring2, k, ygls, 
		ogls,infgls,obsex, obasex,threshold,thresholdT, ag1, ag2, ag3, ag4,
		gkillzones,ndrug,everyh,takeT,gamma,ec50,emin,emax,xm,ym,ke,drugeff,
		seqK1,seqK2,seqK3,seqK4,seqAge1,seqAge2,seqAge3,seqAge4	
		},
   
   (*grate {ratio become gametocytes, ration look like rings} *)
   (*gages {ag1,ag2,ag3,ag4} gametocyte age ranges
   I:0-ag1,II:ag1-ag2,III:ag2-ag3,IV:ag3-ag4,V:ag4-ag5 
   ###life cycle of gametocytes = ag5  ###
   
   obg1 = look like rings (0-obg1)
   obg1 - obg2 = sequestered
   obg3 - ag5 = observed gams
    *)
   ag1=gages[[1]]; 
   ag2=gages[[2]];  
   ag3=gages[[3]];
   ag4=gages[[4]];
   
   (*initial parasites after coming out from the liver*)
   tmp = DistributeN[initN, LC, mu, sigma]; 
   ls = {tmp}; 

   gls = {0};
   ygls = {0};
   ogls = {0};
   
   thresholdT = {0};

   infgls={0};
   ngam = 0;
   nring = {};

   nring2 = {0}; (* gametocytes look like rings *)
   
   (* gcons = List of the constants i.e ,threshold (number of asexual) for switching to gams 
      gcons = {% of being gams/100 (1/h), threshold, look like rings ratio}
   *)
   (* % of becoming gams/100 (1/h) *)   
   grate = gcons[[1]];
  (* observation ratio (look like ring ratio) *)
   k = gcons[[2]];
   (* threshold for switching to gams *)
   threshold = gcons[[3]];
  
  (* there is only one killzone.
     pkpdparms = {killzones_List,ndrug,everyh,takeT,{gamma,ec50,emin,emax},concparms_List}
     concparms = {xm,ym,ke}
     takeT = time that the patient taking the first dose (0 -> runtime)
   *) 
   
   {gkillzones,ndrug,everyh,takeT,{gamma,ec50,emin,emax},{xm,ym,ke}} = pkpdparms;
   
   (* calculate the concentration and its efficacy from t0 = takeT --> runtime + (takeT-t0) *)
   drugeff = EffConc[ndrug,everyh, 2, {gamma,ec50,emin,emax, {xm,ym,ke}} ,runtime];
  
   (* change time to be the same as the system *)
   drugeff[[All,1]] = drugeff[[All,1]]+takeT;
   
(*   
   gkillzone = pkpdparms[[1]];
   ndrug = pkpdparms[[2]];
   everyh = pkpdparms[[3]];
   gamma = pkpdparms[[4,1]];
   ec50 = plpdparms[[4,2]];
   emin = pkpdparms[[4,3]];
   emax = pkpdparms[[4,4]];
   {xm,ym,ke} = concparms;
*)
  (* sequestered asexual + sequestered gametocytes are killed by immune system ???  
    seqparms = {asexual k1,asexual k2,asexual agerange1,asexual agerange2,
    			gam k1,gam k2,gam agerange1,gam agerange2}; 
    
     *)
   {seqK1,seqK2,seqAge1,seqAge2,seqK3,seqK4,seqAge3,seqAge4} = seqparms;
   
   Do[
     tmp = Shiftonehour[tmp,PMR];
     (* kill sequester asexual [parls_,k1_,k2_,agerange_List,t_]*)
   	 tmp = KillSequestered[tmp,seqK1,seqK2,{seqAge1,seqAge2},i];
   	         
	 (* become gametocytes if npara > threshold *)
     If[Total@tmp > threshold ,      	
     	ngam = grate*(First@tmp);
     	tmp = ReplacePart[tmp,1->(1.0-grate)*First@tmp];
     	AppendTo[thresholdT,i];
     	,
     	ngam = 0;
     ];
     
     (* list ot total parasites; (time) past --> future  *)     
     AppendTo[ls, tmp];     
     
     (* list of total gametocytes; (time) young --> old *)
     PrependTo[gls,ngam];
         
     (*drug action for gametocyte*)
     gls = drugActionGam[gls,gkillzones,drugeff,i];
 
     (* immune system kills sequestered gametocytes *)
     gls = KillSequestered[gls,seqK3,seqK4,{seqAge3,seqAge4},i];
  
     (* reset ngam *)
     ngam = 0;
     
     (* look like rings *)
     AppendTo[nring2,CountObservedGame[gls,{1,ag1},k]];
	 
	 (* young observed gametocytes *)
	 AppendTo[ygls,CountObservedGame[gls,{1,ag1},(1-k)]];
	
	 (* old observed gametocytes *)
	 AppendTo[ogls,CountObservedGame[gls,{ag2+1,ag3},1]];
	
	 (* infectiousness gametocytes *)
	 AppendTo[infgls,CountObservedGame[gls,{ag3+1,ag4},1]];
	 
     , {i, 1, runtime}];

   (* observed rings (no gametocytes) *)
   nring = CountRing[#,{11,14},1]&/@ls;
   
   obasex = nring + nring2;
   obsex = ygls + ogls;
      
 (*** output forms ***)
   Which[
   	(* show age distribution over time*)
   	outform==0,ls  
    (* show the last age distribution*)
   ,outform==1,tmp
    (* show the number of gametocyte over time*)
   ,outform==2,gls
    (* big results! *)
   ,outform==3,{ls,gls,nring,obasex,obsex,infgls,thresholdT,drugeff}
    (* threshold time *)
   ,outform==4,thresholdT 
    (* drug effect *)
   ,outform==5,drugeff 
   ]

];

(* immune kill sequestered, drug kill specific age of gam  *)
Game7[initN_, mu_, sigma_, LC_Integer, PMR_, gcons_List, gages_List, seqparms_List,pkpdparms_List, runtime_Integer, outform_Integer] := 
	Module[{ls, i,grate, tmp, gls, ngam, nring, nring2, k, ygls, 
		ogls,infgls,obsex, obasex,threshold,thresholdT, ag1, ag2, ag3, ag4,
		gkillzones,ndrug,everyh,takeT,gamma,ec50,emin,emax,xm,ym,ke,drugeff,
		rho,delta,kappa	
		},
   
   (*grate {ratio become gametocytes, ration look like rings} *)
   (*gages {ag1,ag2,ag3,ag4} gametocyte age ranges; sequestered ag1->ag2
   I:0-ag1,II:ag1-ag2,III:ag2-ag3,IV:ag3-ag4,V:ag4-ag5 
   ###life cycle of gametocytes = ag5  ###
   
   obg1 = look like rings (0-obg1)
   obg1 - obg2 = sequestered
   obg3 - ag5 = observed gams
    *)
   ag1=gages[[1]]; 
   ag2=gages[[2]];  
   ag3=gages[[3]];
   ag4=gages[[4]];
   
   (*initial parasites after coming out from the liver*)
   tmp = DistributeN[initN, LC, mu, sigma]; 
   ls = {tmp}; 

   gls = {0};
   ygls = {0};
   ogls = {0};
   
   thresholdT = {0};

   infgls={0};
   ngam = 0;
   nring = {};

   nring2 = {0}; (* gametocytes look like rings *)
   
   (* gcons = List of the constants i.e ,threshold (number of asexual) for switching to gams 
      gcons = {% of being gams/100 (1/h), threshold, look like rings ratio}
   *)
   (* % of becoming gams/100 (1/h) *)   
   grate = gcons[[1]];
  (* observation ratio (look like ring ratio) *)
   k = gcons[[2]];
   (* threshold for switching to gams *)
   threshold = gcons[[3]];
  
  (* there is only one killzone.
     pkpdparms = {killzones_List,ndrug,everyh,takeT,{gamma,ec50,emin,emax},concparms_List}
     concparms = {xm,ym,ke}
     takeT = time that the patient taking the first dose (0 -> runtime)
   *) 
   
   {gkillzones,ndrug,everyh,takeT,{gamma,ec50,emin,emax},{xm,ym,ke}} = pkpdparms;
   
   (* calculate the concentration and its efficacy from t0 = takeT --> runtime + (takeT-t0) *)
   drugeff = EffConc[ndrug,everyh, 2, {gamma,ec50,emin,emax, {xm,ym,ke}} ,runtime];
 
   (* change time to be the same as the system *)
   drugeff[[All,1]] = drugeff[[All,1]]+takeT;


  (* sequestered asexual + sequestered gametocytes are killed by immune system ???  
    seqparms = {rho,delta,kappa}; 
    		y'(t)=rho S - delta y(t); S=sequestered parasites
    	% of kill = kappa*y[t]/S 	
    		
     *)
   {rho,delta,kappa} = seqparms;
 
   
   
   
   Do[
   	
     tmp = Shiftonehour[tmp,PMR];
     
     (* kill sequester asexual [parls_,k1_,k2_,agerange_List,t_]*)
   	 (*tmp = KillSequestered[tmp,seqK1,seqK2,{seqAge1,seqAge2},i];*)
   	 
     (* here if the number of parasites less that 10^-10, it'll be replaced by 0 *)   	 
   	 tmp = Chop@ImmuneKillSequesteredAX[tmp,rho,delta,kappa,i];
   	         
	 (* become gametocytes if npara > threshold *)
     If[Total@tmp > threshold ,      	
     	ngam = grate*(First@tmp);
     	tmp = ReplacePart[tmp,1->(1.0-grate)*First@tmp];
     	AppendTo[thresholdT,i];
     	,
     	ngam = 0;
     ];
     
     (* list ot total parasites; (time) past --> future  *)     
     AppendTo[ls, tmp];     
     
     (* list of total gametocytes; (time) young --> old *)
     PrependTo[gls,ngam];
         
     (*drug action for gametocyte*)
     gls = drugActionGam[gls,gkillzones,drugeff,i];
     
     (* immune system kills sequestered gametocytes *)
  	 gls = ImmuneKillSequesteredX[gls,rho,delta,kappa,{ag1,ag2},i];
      		
     (* reset ngam *)
     ngam = 0;
     
     (* look like rings *)
     AppendTo[nring2,CountObservedGame[gls,{1,ag1},k]];
	 
	 (* young observed gametocytes *)
	 AppendTo[ygls,CountObservedGame[gls,{1,ag1},(1-k)]];
	
	 (* old observed gametocytes *)
	 AppendTo[ogls,CountObservedGame[gls,{ag2+1,ag3},1]];
	
	 (* infectiousness gametocytes *)
	 AppendTo[infgls,CountObservedGame[gls,{ag3+1,ag4},1]];
	 
     , {i, 1, runtime}];

   (* observed rings (no gametocytes) *)
   nring = CountRing[#,{11,14},1]&/@ls;
   
   obasex = nring + nring2;
   obsex = ygls + ogls;
      
 (*** output forms ***)
   Which[
   	(* show age distribution over time*)
   	outform==0,ls  
    (* show the last age distribution*)
   ,outform==1,tmp
    (* show the number of gametocyte over time*)
   ,outform==2,gls
    (* big results! *)
   ,outform==3,{ls,gls,nring,obasex,obsex,infgls,thresholdT,drugeff}
    (* threshold time *)
   ,outform==4,thresholdT 
    (* drug effect *)
   ,outform==5,drugeff 
   ]

];

(** modified from Game7 , no drug and only immnue  **)
Game8[initN_, mu_, sigma_, LC_Integer, PMR_, gconstants_List, gages_List, immparms_List,trigger_Integer, immuneFunc_Integer, runtime_Integer, outform_Integer] := 
	Module[{ls, i,grate, tmp, gls, ngam, nring, nring2, k, ygls, 
		ogls,infgls,obsex, obasex,threshold,thresholdT, ag1, ag2, ag3, ag4,
		rho,delta,kappa,seqK1,seqK2,seqAge1,seqAge2,seqK3,seqK4,seqAge3,seqAge4,gtriggerT	
		},
   (*trigger = 0 "time" from the beginning of the model
   			 = 1 "threshold"  or parasite density   
     *)
   
   (*
   		immuneFunc = 0 "no immune"
   				   = 1 " k1(1-e^(-k2 t)) "
   				   = 2 " y'(t)=rho S - delta y(t) "	
   *)
   
   (*grate {ratio become gametocytes, ration look like rings} *)
   (*gages {ag1,ag2,ag3,ag4} gametocyte age ranges; sequestered ag1->ag2
   I:0-ag1,II:ag1-ag2,III:ag2-ag3,IV:ag3-ag4,V:ag4-ag5 
   ###life cycle of gametocytes = ag5  ###
   
   obg1 = look like rings (0-obg1)
   obg1 - obg2 = sequestered
   obg3 - ag5 = observed gams ; so in counting gams , ag4 can replaced ag5. 
    *)
   ag1=gages[[1]]; 
   ag2=gages[[2]];  
   ag3=gages[[3]];
   ag4=gages[[4]];
   
   (*initial parasites after coming out from the liver*)
   tmp = DistributeN[initN, LC, mu, sigma]; 
   ls = {tmp}; 

   gls = {0};
   ygls = {0};
   ogls = {0};
   
   (* time that the number of parasites > threshold *)
   thresholdT = {0};

   infgls={0};
   ngam = 0;
   nring = {};

   nring2 = {0}; (* gametocytes look like rings *)
   
   (* gconstants = List of the constants i.e ,threshold (number of asexual) for switching to gams 
      gconstants = {% of being gams/100 (1/h), threshold, look like rings ratio}
      !! threshold will be used if trigger = 1 !!
   *)
   (* % of becoming gams/100 (1/h) *)   
   grate = gconstants[[1]];
  (* observation ratio (look like ring ratio) *)
   k = gconstants[[2]];
   Which[trigger==1,
   (* threshold for switching to gams *)
   threshold = gconstants[[3]];
   ,trigger==0,
   (* trigger time *)
   gtriggerT=gconstants[[3]];
   ];
  (* there is only one killzone.
     pkpdparms = {killzones_List,ndrug,everyh,takeT,{gamma,ec50,emin,emax},concparms_List}
     concparms = {xm,ym,ke}
     takeT = time that the patient taking the first dose (0 -> runtime)
   *) 
   
  (* {gkillzones,ndrug,everyh,takeT,{gamma,ec50,emin,emax},{xm,ym,ke}} = pkpdparms;*)
   
   (* calculate the concentration and its efficacy from t0 = takeT --> runtime + (takeT-t0) *)
   (*drugeff = EffConc[ndrug,everyh, 2, {gamma,ec50,emin,emax, {xm,ym,ke}} ,runtime];*)
 
   (* change time to be the same as the system *)
   (*drugeff[[All,1]] = drugeff[[All,1]]+takeT;*)


  (* sequestered asexual parasites + sequestered gametocytes are killed by immune system (ref.???)  
    seqparms = {rho,delta,kappa}; 
    		y'(t)=rho S - delta y(t); S=sequestered parasites
    	% of kill = kappa*y[t]/S 	
   
   
   killing function for sequestered asexual and sexual stages
           fn=k1*(1.-Exp[-k2*t]);
    seqparms = {k1,k2,k3,k4}   
    k1,k2 for asexual  k3,k4 for games 
 
    		
     *)
   Which[immuneFunc==1&&Length@immparms==3,  
   {rho,delta,kappa} = immparms;,
   immuneFunc==2&&Length@immparms==8,
   {seqK1,seqK2,seqAge1,seqAge2,seqK3,seqK4,seqAge3,seqAge4} = immparms;
   ]; 
 
   
   Do[
   	
     tmp = Shiftonehour[tmp,PMR];
     
     (* kill sequester asexual [parls_,k1_,k2_,agerange_List,t_]*)
   	 (*tmp = KillSequestered[tmp,seqK1,seqK2,{seqAge1,seqAge2},i];*)	 
     (* if the number of parasites less that 10^-10, it'll be replaced by 0 *)   	 
   	 (*tmp = Chop@ImmuneKillSequesteredAX[tmp,rho,delta,kappa,i];*)
   	 
    Which[immuneFunc==1,tmp = Chop@ImmuneKillSequesteredAX[tmp,rho,delta,kappa,i],immuneFunc==2, tmp = KillSequestered[tmp,seqK1,seqK2,{seqAge1,seqAge2},i] ];
   	
	Which[trigger==0,
	(* become gametocytes after gtrigger + LC *)
     If[i >= gtriggerT+LC,
     ngam = grate*First@tmp;
     tmp = ReplacePart[tmp,1->(1-grate)*First@tmp];
     ,
     ngam = 0;
     
     ];
	,trigger==1,
	 (* become gametocytes if npara > threshold *)
     If[Total@tmp > threshold ,      	
     	ngam = grate*(First@tmp);
     	tmp = ReplacePart[tmp,1->(1.0-grate)*First@tmp];
     	AppendTo[thresholdT,i];
     	,
     	ngam = 0;
     ];
	];
	
     (* list ot total parasites; (time) past --> future  *)     
     AppendTo[ls, tmp];     
     
     (* list of total gametocytes; (time) young --> old *)
     PrependTo[gls,ngam];
         
     (*drug action for gametocyte*)
     (*gls = drugActionGam[gls,gkillzones,drugeff,i];*)
     
     (* immune system kills sequestered gametocytes *)
  	 (*gls = ImmuneKillSequesteredX[gls,rho,delta,kappa,{ag1,ag2},i];*)
     (* immune system kills sequestered gametocytes *)
     (*gls = KillSequestered[gls,seqK3,seqK4,{seqAge3,seqAge4},i];*)
     Which[immuneFunc==1,gls = ImmuneKillSequesteredX[gls,rho,delta,kappa,{ag1,ag2},i],immuneFunc==2, gls = KillSequestered[gls,seqK3,seqK4,{seqAge3,seqAge4},i]];
     		
     		 		
     (* reset ngam *)
     ngam = 0;
     
     (* look like rings *)
     AppendTo[nring2,CountObservedGame[gls,{1,ag1},k]];
	 
	 (* young observed gametocytes *)
	 AppendTo[ygls,CountObservedGame[gls,{1,ag1},(1-k)]];
	
	 (* old observed gametocytes *)
	 AppendTo[ogls,CountObservedGame[gls,{ag2+1,ag3},1]];
	
	 (* infectiousness gametocytes *)
	 AppendTo[infgls,CountObservedGame[gls,{ag3+1,ag4},1]];
	 
     , {i, 1, runtime}];

   (* observed rings (no gametocytes) *)
   nring = CountRing[#,{11,14},1]&/@ls;
   
   obasex = nring + nring2;
   obsex = ygls + ogls;
      
 (*** output forms ***)
   Which[
   	(* show age distribution over time*)
   	outform==0,ls  
    (* show the last age distribution*)
   ,outform==1,tmp
    (* show the number of gametocyte over time*)
   ,outform==2,gls
    (* big results! *)
   ,outform==3,{ls,gls,nring,obasex,obsex,infgls,thresholdT}
    (* threshold time *)
   ,outform==4,thresholdT 
   ]

];


(***)


(* kill gametocytes by drug *)
drugActionSexual[gls_List, killzones_List, drugeff_List, time_Integer]:=Module[{glength,tempgls,effattime,bKZ,eKZ,CEattime},
	{bKZ,eKZ} = killzones;
	glength = Length@gls;
	tempgls = gls;
    CEattime = Select[drugeff,#[[1]]==time&]; (* "%" of killing *)
	
	If[(glength < eKZ) && (glength > bKZ) && CEattime != {},
		effattime = 0.01*CEattime[[1,3]];
		tempgls[[bKZ;;glength]] = (1.0-effattime)*gls[[bKZ;;glength]];
		,
		If[(eKZ <= glength) && CEattime != {},
			effattime = 0.01*CEattime[[1,3]];
			tempgls[[bKZ;;eKZ]] = (1.0-effattime)*gls[[bKZ;;eKZ]];
		];
	];
	
	tempgls
];

(**kill parasites (1 kill zone) **)
drugActionASexual[parals_List,parakillzone_List,drugeff_List,time_Integer]:=Module[
	{paralength,bKZ,eKZ,CEattime,tmpparals},
	
	(* parals = parasite distribution *)
	tmpparals = parals;
	{bKZ,eKZ}=parakillzone;
	
	(* force the input of killzone not lower or greater than the length of parals *)
	Which[bKZ<1,bKZ=1;
	,
	eKZ>Length@parals,eKZ=Length@parals;
	];
	
	paralength=Length@parals;
	CEattime=Select[drugeff,#[[1]]==time&];
	
	If[CEattime!={},	
	tmpparals[[bKZ;;eKZ]] = (1-CEattime[[1,3]]*0.01)*parals[[bKZ;;eKZ]];
	];
	tmpparals 
];


(** modified from Game8 , add drug which kill both gams and asexual (depends on drugAct), one set of PD parameters for both gams and asexual  **)
Game9[initN_, mu_, sigma_, LC_Integer, PMR_, gconstants_List, gages_List, immparms_List,trigger_Integer, immuneFunc_Integer, drugAct_List, pkpdparms_List, runtime_Integer, outform_Integer] := 
	Module[{ls, i,grate, tmp, gls, ngam, nring, nring2, k, ygls, 
		ogls,infgls,obsex, obasex,threshold,thresholdT, ag1, ag2, ag3, ag4,
		rho,delta,kappa,seqK1,seqK2,seqAge1,seqAge2,seqK3,seqK4,seqAge3,seqAge4,gtriggerT,
		parakillzone,gkillzone,ndrug,everyh,takeT,gamma,ec50,emin,emax,xm,ym,ke,drugeff	
		},
   (*trigger = 0 "time" from the beginning of the model
   			 = 1 "threshold"  or parasite density   
     *)
   
   (*
   		immuneFunc = 0 "no immune"
   				   = 1 " k1(1-e^(-k2 t)) "
   				   = 2 " y'(t)=rho S - delta y(t) "	
   *)
   
   (*grate {ratio become gametocytes, ration look like rings} *)
   (*gages {ag1,ag2,ag3,ag4} gametocyte age ranges; sequestered ag1->ag2
   I:0-ag1,II:ag1-ag2,III:ag2-ag3,IV:ag3-ag4,V:ag4-ag5 
   ###life cycle of gametocytes = ag5  ###
   
   obg1 = look like rings (0-obg1)
   obg1 - obg2 = sequestered
   obg3 - ag5 = observed gams ; so in counting gams , ag4 can replaced ag5. 
    *)
   ag1=gages[[1]]; 
   ag2=gages[[2]];  
   ag3=gages[[3]];
   ag4=gages[[4]];
   
   (*initial parasites after coming out from the liver*)
   tmp = DistributeN[initN, LC, mu, sigma]; 
   ls = {tmp}; 

   gls = {0};
   ygls = {0};
   ogls = {0};
   
   (* time that the number of parasites > threshold *)
   thresholdT = {0};

   infgls={0};
   ngam = 0;
   nring = {};

   nring2 = {0}; (* gametocytes look like rings *)
   
   (* gconstants = List of the constants i.e ,threshold (number of asexual) for switching to gams 
      gconstants = {% of being gams/100 (1/h), threshold, look like rings ratio}
      !! threshold will be used if trigger = 1 !!
   *)
   (* % of becoming gams/100 (1/h) *)   
   grate = gconstants[[1]];
  (* observation ratio (look like ring ratio) *)
   k = gconstants[[2]];
   Which[trigger==1,
   (* threshold for switching to gams *)
   threshold = gconstants[[3]];
   ,trigger==0,
   (* trigger time *)
   gtriggerT=gconstants[[3]];
   ];
  
  (* drugAct = {asexual, sexual} 
  		asexual/sexual = 0 not kill, 1 kill by the drug
  *) 
  (* there is only one killzone.
     pkpdparms = {parakillzone_List,gkillzones_List,ndrug,everyh,takeT,{gamma,ec50,emin,emax},concparms_List}
     concparms = {xm,ym,ke}
     takeT = time that the patient taking the first dose (0 -> runtime)
     
     !!pkpdparms will be read when Total@drugAct > 0.     
   *) 
   
   If[Total@drugAct>0,
   	
   		{parakillzone,gkillzone,ndrug,everyh,takeT,{gamma,ec50,emin,emax},{xm,ym,ke}} = pkpdparms;
   
   		(* calculate the concentration and its efficacy from t0 = takeT --> runtime + (takeT-t0) *)
   		drugeff = EffConc1[ndrug,everyh,{gamma,ec50,emin,emax, {xm,ym,ke}} ,runtime];
 
   		(* change time to be the same as the system *)
   		drugeff[[All,1]] = drugeff[[All,1]]+takeT;
   	];
	

  (* sequestered asexual parasites + sequestered gametocytes are killed by immune system (ref.???)  
    seqparms = {rho,delta,kappa}; 
    		y'(t)=rho S - delta y(t); S=sequestered parasites
    	% of kill = kappa*y[t]/S 	
   killing function for sequestered asexual and sexual stages
           fn=k1*(1.-Exp[-k2*t]);
    seqparms = {k1,k2,k3,k4}   
    k1,k2 for asexual  k3,k4 for games    		
     *)
     
   Which[immuneFunc==1&&Length@immparms==3,  
   		{rho,delta,kappa} = immparms;
   		,
   		immuneFunc==2&&Length@immparms==8,
   		{seqK1,seqK2,seqAge1,seqAge2,seqK3,seqK4,seqAge3,seqAge4} = immparms;
   ]; 
 
   
   Do[
   	
     tmp = Shiftonehour[tmp,PMR];
     
    If[drugAct[[1]]==1,
     (* drug action for asexual parasites *)
     tmp = drugActionASexual[tmp,parakillzone,drugeff,i];
    ];     
     (* kill sequester asexual [parls_,k1_,k2_,agerange_List,t_]*)
   	 (*tmp = KillSequestered[tmp,seqK1,seqK2,{seqAge1,seqAge2},i];*)	 
     (* if the number of parasites less that 10^-10, it'll be replaced by 0 *)   	 
   	 (*tmp = Chop@ImmuneKillSequesteredAX[tmp,rho,delta,kappa,i];*)
   	 
    Which[immuneFunc==1,tmp = Chop@ImmuneKillSequesteredAX[tmp,rho,delta,kappa,i],immuneFunc==2, tmp = KillSequestered[tmp,seqK1,seqK2,{seqAge1,seqAge2},i] ];
   	
	Which[
		trigger==0,
	(* become gametocytes after gtrigger + LC *)
     If[i >= gtriggerT+LC,
     	ngam = grate*First@tmp;
     	tmp = ReplacePart[tmp,1->(1-grate)*First@tmp];
     	,
     	ngam = 0;
     	];
	,
	trigger==1,
	 (* become gametocytes if npara > threshold *)
     If[Total@tmp > threshold ,      	
     	ngam = grate*(First@tmp);
     	tmp = ReplacePart[tmp,1->(1.0-grate)*First@tmp];
     	AppendTo[thresholdT,i];
     	,
     	ngam = 0;
     ];
	];
	
     (* list ot total parasites; (time) past --> future  *)     
     AppendTo[ls, tmp];     
     
     (* list of total gametocytes; (time) young --> old *)
     PrependTo[gls,ngam];
         
     If[drugAct[[2]]==1,
     	(*drug action for gametocyte*)
     	gls = drugActionSexual[gls,gkillzone,drugeff,i];
     ];
     (* immune system kills sequestered gametocytes *)
  	 (*gls = ImmuneKillSequesteredX[gls,rho,delta,kappa,{ag1,ag2},i];*)
     (* immune system kills sequestered gametocytes *)
     (*gls = KillSequestered[gls,seqK3,seqK4,{seqAge3,seqAge4},i];*)
     Which[immuneFunc==1,gls = ImmuneKillSequesteredX[gls,rho,delta,kappa,{ag1,ag2},i],immuneFunc==2, gls = KillSequestered[gls,seqK3,seqK4,{seqAge3,seqAge4},i]];
     		
     		 		
     (* reset ngam *)
     ngam = 0;
     
     (* look like rings *)
     AppendTo[nring2,CountObservedGame[gls,{1,ag1},k]];
	 
	 (* young observed gametocytes *)
	 AppendTo[ygls,CountObservedGame[gls,{1,ag1},(1-k)]];
	
	 (* old observed gametocytes ; from ag2 to ag4 *)
	 AppendTo[ogls,CountObservedGame[gls,{ag2+1,ag4},1]];
	
	 (* infectiousness gametocytes *)
	 AppendTo[infgls,CountObservedGame[gls,{ag3+1,ag4},1]];
	 
     , {i, 1, runtime}];

   (* observed rings (no gametocytes) *)
   nring = CountRing[#,{11,14},1]&/@ls;
   
   (* replace the approximate numbers that are close to zero (10^-10) by zero *)
   obasex = Chop@(nring + nring2);
   obsex = Chop@(ygls + ogls);
      
 (*** output forms ***)
   Which[
   	(* show age distribution over time*)
   	outform==0,ls  
    (* show the last age distribution*)
   ,outform==1,tmp
    (* show the number of gametocyte over time*)
   ,outform==2,gls
    (* big results! *)
   ,outform==3,{ls,gls,nring,obasex,obsex,infgls,thresholdT}
    (* threshold time *)
   ,outform==4,thresholdT 
   ]

];

(*********)
(** modified from Game9 , gams and asexual have their own EC-curves for the drug effect  **)
Game10[initN_, mu_, sigma_, LC_Integer, PMR_, gconstants_List, gages_List, immparms_List,trigger_Integer, immuneFunc_Integer, drugAct_List, pkpdparms_List, runtime_Integer, outform_Integer] := 
	Module[{ls, i,grate, tmp, gls, ngam, nring, nring2, k, ygls, 
		ogls,infgls,obsex, obasex,threshold,thresholdT, ag1, ag2, ag3, ag4,
		rho,delta,kappa,seqK1,seqK2,seqAge1,seqAge2,seqK3,seqK4,seqAge3,seqAge4,gtriggerT,
		parakillzone,gkillzone,ndrug,everyh,takeT,gamma,ec50,emin,emax,xm,ym,ke,drugeffAS,drugeffS	
		},
   (*trigger = 0 "time" from the beginning of the model
   			 = 1 "threshold"  or parasite density   
     *)
   
   (*
   		immuneFunc = 0 "no immune"
   				   = 1 " k1(1-e^(-k2 t)) "
   				   = 2 " y'(t)=rho S - delta y(t) "	
   *)
   
   (*grate {ratio become gametocytes, ration look like rings} *)
   (*gages {ag1,ag2,ag3,ag4} gametocyte age ranges; sequestered ag1->ag2
   I:0-ag1,II:ag1-ag2,III:ag2-ag3,IV:ag3-ag4,V:ag4-ag5 
   ###life cycle of gametocytes = ag5  ###
   
   obg1 = look like rings (0-obg1)
   obg1 - obg2 = sequestered
   obg3 - ag5 = observed gams ; so in counting gams , ag4 can replaced ag5. 
    *)
   ag1=gages[[1]]; 
   ag2=gages[[2]];  
   ag3=gages[[3]];
   ag4=gages[[4]];
   
   (*initial parasites after coming out from the liver*)
   tmp = DistributeN[initN, LC, mu, sigma]; 
   ls = {tmp}; 

   gls = {0};
   ygls = {0};
   ogls = {0};
   
   (* time that the number of parasites > threshold *)
   thresholdT = {0};

   infgls={0};
   ngam = 0;
   nring = {};

   nring2 = {0}; (* gametocytes look like rings *)
   
   (* gconstants = List of the constants i.e ,threshold (number of asexual) for switching to gams 
      gconstants = {% of being gams/100 (1/h), threshold, look like rings ratio}
      !! threshold will be used if trigger = 1 !!
   *)
   (* % of becoming gams/100 (1/h) *)   
   grate = gconstants[[1]];
  (* observation ratio (look like ring ratio) *)
   k = gconstants[[2]];
   Which[trigger==1,
   (* threshold for switching to gams *)
   threshold = gconstants[[3]];
   ,trigger==0,
   (* trigger time *)
   gtriggerT=gconstants[[3]];
   ];
  
  (* drugAct = {asexual, sexual} 
  		asexual/sexual = 0 not kill, 1 kill by the drug
  *) 
  (* there is only one killzone.
     pkpdparms = {parakillzone_List,gkillzones_List,ndrug,everyh,takeT,{gamma,ec50,emin,emax},concparms_List}
     concparms = {xm,ym,ke}
     takeT = time that the patient taking the first dose (0 -> runtime)
     
     !!pkpdparms will be read when Total@drugAct > 0.    
     
     PD parameters are in the list format.
     {gamma_List,ec50_List,emin_List,emax_List} = {{gammaAS,gammaS},{ec50AS,ec50S},{eminAS,eminS},{emaxAS,emaxS}}  
      
   *) 
   
   If[Total@drugAct>0,
   	
   		{parakillzone,gkillzone,ndrug,everyh,takeT,{gamma,ec50,emin,emax},{xm,ym,ke}} = pkpdparms;
   
   		(* calculate the concentration and its efficacy from t0 = takeT --> runtime + (takeT-t0) *)
   		drugeffAS = EffConc1[ndrug,everyh,{gamma[[1]],ec50[[1]],emin[[1]],emax[[1]], {xm,ym,ke}} ,runtime];
   		drugeffS = EffConc1[ndrug,everyh,{gamma[[2]],ec50[[2]],emin[[2]],emax[[2]], {xm,ym,ke}} ,runtime];
 
   		(* change time to be the same as the system *)
   		drugeffAS[[All,1]] = drugeffAS[[All,1]]+takeT;
   		drugeffS[[All,1]] = drugeffS[[All,1]] + takeT;
   	];
	

  (* sequestered asexual parasites + sequestered gametocytes are killed by immune system (ref.???)  
    seqparms = {rho,delta,kappa}; 
    		y'(t)=rho S - delta y(t); S=sequestered parasites
    	% of kill = kappa*y[t]/S 	
   killing function for sequestered asexual and sexual stages
           fn=k1*(1.-Exp[-k2*t]);
    seqparms = {k1,k2,k3,k4}   
    k1,k2 for asexual  k3,k4 for games    		
     *)
     
   Which[immuneFunc==1&&Length@immparms==3,  
   		{rho,delta,kappa} = immparms;
   		,
   		immuneFunc==2&&Length@immparms==8,
   		{seqK1,seqK2,seqAge1,seqAge2,seqK3,seqK4,seqAge3,seqAge4} = immparms;
   ]; 
 
   
   Do[
   	
     tmp = Shiftonehour[tmp,PMR];
     
    If[drugAct[[1]]==1,
     (* drug action for asexual parasites *)
     tmp = drugActionASexual[tmp,parakillzone,drugeffAS,i];
    ];     
     (* kill sequester asexual [parls_,k1_,k2_,agerange_List,t_]*)
   	 (*tmp = KillSequestered[tmp,seqK1,seqK2,{seqAge1,seqAge2},i];*)	 
     (* if the number of parasites less that 10^-10, it'll be replaced by 0 *)   	 
   	 (*tmp = Chop@ImmuneKillSequesteredAX[tmp,rho,delta,kappa,i];*)
   	 
    Which[immuneFunc==1,tmp = Chop@ImmuneKillSequesteredAX[tmp,rho,delta,kappa,i],immuneFunc==2, tmp = KillSequestered[tmp,seqK1,seqK2,{seqAge1,seqAge2},i] ];
   	
	Which[
		trigger==0,
	(* become gametocytes after gtrigger + LC *)
     If[i >= gtriggerT+LC,
     	ngam = grate*First@tmp;
     	tmp = ReplacePart[tmp,1->(1-grate)*First@tmp];
     	,
     	ngam = 0;
     	];
	,
	trigger==1,
	 (* become gametocytes if npara > threshold *)
     If[Total@tmp > threshold ,      	
     	ngam = grate*(First@tmp);
     	tmp = ReplacePart[tmp,1->(1.0-grate)*First@tmp];
     	AppendTo[thresholdT,i];
     	,
     	ngam = 0;
     ];
	];
	
     (* list ot total parasites; (time) past --> future  *)     
     AppendTo[ls, tmp];     
     
     (* list of total gametocytes; (time) young --> old *)
     PrependTo[gls,ngam];
         
     If[drugAct[[2]]==1,
     	(*drug action for gametocyte*)
     	gls = drugActionSexual[gls,gkillzone,drugeffS,i];
     ];
     (* immune system kills sequestered gametocytes *)
  	 (*gls = ImmuneKillSequesteredX[gls,rho,delta,kappa,{ag1,ag2},i];*)
     (* immune system kills sequestered gametocytes *)
     (*gls = KillSequestered[gls,seqK3,seqK4,{seqAge3,seqAge4},i];*)
     Which[immuneFunc==1,gls = ImmuneKillSequesteredX[gls,rho,delta,kappa,{ag1,ag2},i],immuneFunc==2, gls = KillSequestered[gls,seqK3,seqK4,{seqAge3,seqAge4},i]];
     		
     		 		
     (* reset ngam *)
     ngam = 0;
     
     (* look like rings *)
     AppendTo[nring2,CountObservedGame[gls,{1,ag1},k]];
	 
	 (* young observed gametocytes *)
	 AppendTo[ygls,CountObservedGame[gls,{1,ag1},(1-k)]];
	
	 (* old observed gametocytes ; from ag2 to ag4 *)
	 AppendTo[ogls,CountObservedGame[gls,{ag2+1,ag4},1]];
	
	 (* infectiousness gametocytes *)
	 AppendTo[infgls,CountObservedGame[gls,{ag3+1,ag4},1]];
	 
     , {i, 1, runtime}];

   (* observed rings (no gametocytes) *)
   nring = CountRing[#,{11,14},1]&/@ls;
   
   (* replace the approximate numbers that are close to zero (10^-10) by zero *)
   obasex = Chop@(nring + nring2);
   obsex = Chop@(ygls + ogls);
      
 (*** output forms ***)
   Which[
   	(* show age distribution over time*)
   	outform==0,ls  
    (* show the last age distribution*)
   ,outform==1,tmp
    (* show the number of gametocyte over time*)
   ,outform==2,gls
    (* big results! *)
   ,outform==3,{ls,gls,nring,obasex,obsex,infgls,thresholdT}
    (* threshold time *)
   ,outform==4,thresholdT 
   ]

];

(**  **)

(** modified from Game10 , the threshold and time triggers are combined.  **)
Game11[initN_, mu_, sigma_, LC_Integer, PMR_, gconstants_List, gages_List, immparms_List,trigger_Integer, immuneFunc_Integer, drugAct_List, pkpdparms_List, runtime_Integer, outform_Integer] := 
	Module[{ls, i,grate, tmp, gls, ngam, nring, nring2, k, ygls, 
		ogls,infgls,obsex, obasex,threshold,thresholdT, ag1, ag2, ag3, ag4,
		rho,delta,kappa,seqK1,seqK2,seqAge1,seqAge2,seqK3,seqK4,seqAge3,seqAge4,gtriggerT,
		parakillzone,gkillzone,ndrug,everyh,takeT,gamma,ec50,emin,emax,xm,ym,ke,drugeffAS,drugeffS
		},
   (*trigger = 0 "time" from the beginning of the model
   			 = 1 "threshold"  or parasite density   
     *)
   
   (*
   		immuneFunc = 0 "no immune"
   				   = 1 " k1(1-e^(-k2 t)) "
   				   = 2 " y'(t)=rho S - delta y(t) "	
   *)
   
   (*grate {ratio become gametocytes, ratio look like rings} *)
   (*gages {ag1,ag2,ag3,ag4} gametocyte age ranges; sequestered ag1->ag2
   I:0-ag1,II:ag1-ag2,III:ag2-ag3,IV:ag3-ag4,V:ag4-ag5 
   ###life cycle of gametocytes = ag5  ###
   
   obg1 = look like rings (0-obg1)
   obg1 - obg2 = sequestered
   obg3 - ag5 = observed gams ; so in counting gams , ag4 can replaced ag5. 
    *)
   ag1=gages[[1]]; 
   ag2=gages[[2]];  
   ag3=gages[[3]];
   ag4=gages[[4]];
   
   (*initial parasites after coming out from the liver*)
   tmp = DistributeN[initN, LC, mu, sigma]; 
   ls = {tmp}; 

   gls = {0};
   ygls = {0};
   ogls = {0};
   
   (* time that the number of parasites > threshold *)
   thresholdT = {0};

   infgls={0};
   ngam = 0;
   nring = {};

   nring2 = {0}; (* gametocytes look like rings *)
   
   (* gconstants = List of the constants i.e ,threshold (number of asexual) for switching to gams 
      if trigger = 0 
      gconstants = {% of being gams/100 (1/h), look like rings ratio, gtriggerT}
      if trigger = 1
      gconstants = {% of being gams/100 (1/h), look like rings ratio, threshold}
	  if trigger = 2
      gconstants = {% of being gams/100 (1/h), look like rings ratio, gtriggerT, threshold}      
   *)
   (* % of becoming gams/100 (1/h) *)   
   grate = gconstants[[1]];
  (* observation ratio (look like ring ratio) *)
   k = gconstants[[2]];
   Which[trigger==1,
   (* threshold for switching to gams *)
   threshold = gconstants[[3]];
   ,trigger==0,
   (* trigger time *)
   gtriggerT=gconstants[[3]];
   ,trigger==2,
   (* time and threshold *)
   gtriggerT = gconstants[[3]];
   threshold = gconstants[[4]];
   ];
  
  (* drugAct = {asexual, sexual} 
  		asexual/sexual = 0 not kill, 1 kill by the drug
  *) 
  (* there is only one killzone.
     pkpdparms = {parakillzone_List,gkillzones_List,ndrug,everyh,takeT,{gamma,ec50,emin,emax},concparms_List}
     concparms = {xm,ym,ke}
     takeT = time that the patient taking the first dose (0 -> runtime)
     
     !!pkpdparms will be read when Total@drugAct > 0.    
     
     PD parameters are in the list format.
     {gamma_List,ec50_List,emin_List,emax_List} = {{gammaAS,gammaS},{ec50AS,ec50S},{eminAS,eminS},{emaxAS,emaxS}}  
      
   *) 
   
   If[Total@drugAct>0,
   	
   		{parakillzone,gkillzone,ndrug,everyh,takeT,{gamma,ec50,emin,emax},{xm,ym,ke}} = pkpdparms;
   
   		(* calculate the concentration and its efficacy from t0 = takeT --> runtime + (takeT-t0) *)
   		drugeffAS = EffConc1[ndrug,everyh,{gamma[[1]],ec50[[1]],emin[[1]],emax[[1]], {xm,ym,ke}} ,runtime];
   		drugeffS = EffConc1[ndrug,everyh,{gamma[[2]],ec50[[2]],emin[[2]],emax[[2]], {xm,ym,ke}} ,runtime];
 
   		(* change time to be the same as the system *)
   		drugeffAS[[All,1]] = drugeffAS[[All,1]]+takeT;
   		drugeffS[[All,1]] = drugeffS[[All,1]] + takeT;
   	];
	

  (* sequestered asexual parasites + sequestered gametocytes are killed by immune system (ref.???)  
    seqparms = {rho,delta,kappa}; 
    		y'(t)=rho S - delta y(t); S=sequestered parasites
    	% of kill = kappa*y[t]/S 	
   killing function for sequestered asexual and sexual stages
           fn=k1*(1.-Exp[-k2*t]);
    seqparms = {k1,k2,k3,k4}   
    k1,k2 for asexual  k3,k4 for games    		
     *)
     
   Which[immuneFunc==1&&Length@immparms==3,  
   		{rho,delta,kappa} = immparms;
   		,
   		immuneFunc==2&&Length@immparms==8,
   		{seqK1,seqK2,seqAge1,seqAge2,seqK3,seqK4,seqAge3,seqAge4} = immparms;
   ]; 
 
   
   Do[
   	
     tmp = Shiftonehour[tmp,PMR];
     
    If[drugAct[[1]]==1,
     (* drug action for asexual parasites *)
     tmp = drugActionASexual[tmp,parakillzone,drugeffAS,i];
    ];     
     (* kill sequester asexual [parls_,k1_,k2_,agerange_List,t_]*)
   	 (*tmp = KillSequestered[tmp,seqK1,seqK2,{seqAge1,seqAge2},i];*)	 
     (* if the number of parasites less that 10^-10, it'll be replaced by 0 *)   	 
   	 (*tmp = Chop@ImmuneKillSequesteredAX[tmp,rho,delta,kappa,i];*)
   	 
    Which[
    	immuneFunc==1,tmp = Chop@ImmuneKillSequesteredAX[tmp,rho,delta,kappa,i],
    	immuneFunc==2, tmp = KillSequestered[tmp,seqK1,seqK2,{seqAge1,seqAge2},i] 
    ];
   	
   	Print[tmp];
   	
	Which[
		trigger==0,
	(* become gametocytes after gtrigger + LC *)
     If[i >= gtriggerT+LC,
     	ngam = grate*First@tmp;
     	tmp = ReplacePart[tmp,1->(1.0-grate)*First@tmp];
     	,
     	ngam = 0;
     ];
	,
	trigger==1,
	 (* become gametocytes if npara > threshold *)
     If[Total@tmp > threshold ,      	
     	ngam = grate*(First@tmp);
     	tmp = ReplacePart[tmp,1->(1.0-grate)*First@tmp];
     	AppendTo[thresholdT,i];
     	,
     	ngam = 0;
     ];
     ,
     trigger==2,
     (* become gametocytes if npara > threshold && i >=gtriggerT+LC *)
      If[i>=gtriggerT+LC && Total@tmp > threshold,
      	ngam = grate*First@tmp;
      	tmp = ReplacePart[tmp,1->(1.0-grate)*First@tmp];
      	AppendTo[thresholdT,i];
      	,
      	ngam = 0;
      ];	
     	
	];
	
     (* list ot total parasites; (time) past --> future  *)     
     AppendTo[ls, tmp];     
     
     (* list of total gametocytes; (time) young --> old *)
     PrependTo[gls,ngam];
         
     If[drugAct[[2]]==1,
     	(*drug action for gametocyte*)
     	gls = drugActionSexual[gls,gkillzone,drugeffS,i];
     ];
     (* immune system kills sequestered gametocytes *)
  	 (*gls = ImmuneKillSequesteredX[gls,rho,delta,kappa,{ag1,ag2},i];*)
     (* immune system kills sequestered gametocytes *)
     (*gls = KillSequestered[gls,seqK3,seqK4,{seqAge3,seqAge4},i];*)
     Which[immuneFunc==1,gls = ImmuneKillSequesteredX[gls,rho,delta,kappa,{ag1,ag2},i],immuneFunc==2, gls = KillSequestered[gls,seqK3,seqK4,{seqAge3,seqAge4},i]];
     		
     		 		
     (* reset ngam *)
     ngam = 0;
     
     (* look like rings *)
     AppendTo[nring2,CountObservedGame[gls,{1,ag1},k]];
	 
	 (* young observed gametocytes *)
	 AppendTo[ygls,CountObservedGame[gls,{1,ag1},(1-k)]];
	
	 (* old observed gametocytes ; from ag2 to ag4 *)
	 AppendTo[ogls,CountObservedGame[gls,{ag2+1,ag4},1]];
	
	 (* infectiousness gametocytes *)
	 AppendTo[infgls,CountObservedGame[gls,{ag3+1,ag4},1]];
	 
     , {i, 1, runtime}];

   (* observed rings (no gametocytes) *)
   nring = CountRing[#,{11,14},1]&/@ls;
   
   (* replace the approximate numbers that are close to zero (10^-10) by zero *)
   obasex = Chop@(nring + nring2);
   obsex = Chop@(ygls + ogls);
      
 (*** output forms ***)
   Which[
   	(* show age distribution over time*)
   	outform==0,ls  
    (* show the last age distribution*)
   ,outform==1,tmp
    (* show the number of gametocyte over time*)
   ,outform==2,gls
    (* big results! *)
   ,outform==3,{ls,gls,nring,obasex,obsex,infgls,thresholdT}
    (* threshold time *)
   ,outform==4,thresholdT 
   ]

];


Game12[initN_, mu_, sigma_, LC_Integer, PMR_, gconstants_List, gages_List, immparms_List,trigger_Integer, immuneFunc_Integer, drugAct_List, pkpdparms_List, runtime_Integer, outform_Integer] := 
	Module[{ls, i,grate, tmp, gls, ngam, nring, nring2, k, ygls, 
		ogls,infgls,obsex, obasex,threshold,thresholdT, ag1, ag2, ag3, ag4,
		rho,delta,kappa,seqK1,seqK2,seqAge1,seqAge2,seqK3,seqK4,seqAge3,seqAge4,gtriggerT,
		parakillzone,gkillzone,ndrug,everyh,takeT,gamma,ec50,emin,emax,xm,ym,ke,drugeffAS,drugeffS,
		Ima1,Ima2
		},
   (*trigger = 0 "time" from the beginning of the model
   			 = 1 "threshold"  or parasite density   
     *)
   
   (*
   		immuneFunc = 0 "no immune"
   				   = 1 " k1(1-e^(-k2 t)) "
   				   = 2 " y'(t)=rho S - delta y(t) "	
   *)
   
   (*grate {ratio become gametocytes, ratio look like rings} *)
   (*gages {ag1,ag2,ag3,ag4} gametocyte age ranges; sequestered ag1->ag2
   I:0-ag1,II:ag1-ag2,III:ag2-ag3,IV:ag3-ag4,V:ag4-ag5 
   ###life cycle of gametocytes = ag5  ###
   
   obg1 = look like rings (0-obg1)
   obg1 - obg2 = sequestered
   obg3 - ag5 = observed gams ; so in counting gams , ag4 can replaced ag5. 
    *)
   ag1=gages[[1]]; 
   ag2=gages[[2]];  
   ag3=gages[[3]];
   ag4=gages[[4]];
   
   (*initial parasites after coming out from the liver*)
   tmp = DistributeN2[initN, LC, mu, sigma, PMR]; 
   ls = {tmp}; 

   gls = {0};
   ygls = {0};
   ogls = {0};
   
   (* time that the number of parasites > threshold *)
   thresholdT = {0};

   infgls={0};
   ngam = 0;
   nring = {};

   nring2 = {0}; (* gametocytes look like rings *)
   
   (* gconstants = List of the constants i.e ,threshold (number of asexual) for switching to gams 
      if trigger = 0 
      gconstants = {% of being gams/100 (1/h), look like rings ratio, gtriggerT}
      if trigger = 1
      gconstants = {% of being gams/100 (1/h), look like rings ratio, threshold}
	  if trigger = 2
      gconstants = {% of being gams/100 (1/h), look like rings ratio, gtriggerT, threshold}      
   *)
   (* % of becoming gams/100 (1/h) *)   
   grate = gconstants[[1]];
  (* observation ratio (look like ring ratio) *)
   k = gconstants[[2]];
   Which[trigger==1,
   (* threshold for switching to gams *)
   threshold = gconstants[[3]];
   ,trigger==0,
   (* trigger time *)
   gtriggerT=gconstants[[3]];
   ,trigger==2,
   (* time and threshold *)
   gtriggerT = gconstants[[3]];
   threshold = gconstants[[4]];
   ];
  
  (* drugAct = {asexual, sexual} 
  		asexual/sexual = 0 not kill, 1 kill by the drug
  *) 
  (* there is only one killzone.
     pkpdparms = {parakillzone_List,gkillzones_List,ndrug,everyh,takeT,{gamma,ec50,emin,emax},concparms_List}
     concparms = {xm,ym,ke}
     takeT = time that the patient taking the first dose (0 -> runtime)
     
     !!pkpdparms will be read when Total@drugAct > 0.    
     
     PD parameters are in the list format.
     {gamma_List,ec50_List,emin_List,emax_List} = {{gammaAS,gammaS},{ec50AS,ec50S},{eminAS,eminS},{emaxAS,emaxS}}  
      
   *) 
   
   If[Total@drugAct>0,
   	
   		{parakillzone,gkillzone,ndrug,everyh,takeT,{gamma,ec50,emin,emax},{xm,ym,ke}} = pkpdparms;
   
   		(* calculate the concentration and its efficacy from t0 = takeT --> runtime + (takeT-t0) *)
   		drugeffAS = EffConc1[ndrug,everyh,{gamma[[1]],ec50[[1]],emin[[1]],emax[[1]], {xm,ym,ke}} ,runtime];
   		drugeffS = EffConc1[ndrug,everyh,{gamma[[2]],ec50[[2]],emin[[2]],emax[[2]], {xm,ym,ke}} ,runtime];
 
   		(* change time to be the same as the system *)
   		drugeffAS[[All,1]] = drugeffAS[[All,1]]+takeT;
   		drugeffS[[All,1]] = drugeffS[[All,1]] + takeT;
   	];
	

  (* sequestered asexual parasites + sequestered gametocytes are killed by immune system (ref.???)  
    seqparms = {rho,delta,kappa}; 
    		y'(t)=rho S - delta y(t); S=sequestered parasites
    	% of kill = kappa*y[t]/S 	
   killing function for sequestered asexual and sexual stages
           fn=k1*(1.-Exp[-k2*t]);
    seqparms = {k1,k2,k3,k4}   
    k1,k2 for asexual  k3,k4 for games    		
     *)
     
   Which[immuneFunc==1&&Length@immparms==3,  
   		{rho,delta,kappa} = immparms;
   		,
   		immuneFunc==2&&Length@immparms==8,
   		{seqK1,seqK2,seqAge1,seqAge2,seqK3,seqK4,seqAge3,seqAge4} = immparms;
   		,
   		immuneFunc==3&&Length@immparms==2,
   		{Ima1,Ima2} = immparms;
   ]; 
 
   
   Do[
   	
     tmp = Shiftonehour[tmp,PMR];
     
    If[drugAct[[1]]==1,
     (* drug action for asexual parasites *)
     tmp = drugActionASexual[tmp,parakillzone,drugeffAS,i];
    ];     
     (* kill sequester asexual [parls_,k1_,k2_,agerange_List,t_]*)
   	 (*tmp = KillSequestered[tmp,seqK1,seqK2,{seqAge1,seqAge2},i];*)	 
     (* if the number of parasites less that 10^-10, it'll be replaced by 0 *)   	 
   	 (*tmp = Chop@ImmuneKillSequesteredAX[tmp,rho,delta,kappa,i];*)
   	 
    Which[
    	immuneFunc==1,tmp = Chop@ImmuneKillSequesteredAX[tmp,rho,delta,kappa,i],
    	immuneFunc==2, tmp = KillSequestered[tmp,seqK1,seqK2,{seqAge1,seqAge2},i],
    	immuneFunc==3, tmp = tmp*Exp[-((Ima1*i)+(i^2)*(Ima2-Ima1)/runtime)] 
    ];
   	
   	Print[tmp];
   	
	Which[
		trigger==0,
	(* become gametocytes after gtrigger + LC *)
     If[i >= gtriggerT+LC,
     	ngam = grate*First[tmp];  (* S->G *)
     	tmp = ReplacePart[tmp,1->(1.0-grate)*First@tmp];
     	,
     	ngam = 0;
     ];
	,
	trigger==1,
	 (* become gametocytes if npara > threshold *)
     If[Total@tmp > threshold ,      	
     	ngam = grate*First[tmp];
     	tmp = ReplacePart[tmp,1->(1.0-grate)*First@tmp];
     	AppendTo[thresholdT,i];
     	,
     	ngam = 0;
     ];
     ,
     trigger==2,
     (* become gametocytes if npara > threshold && i >=gtriggerT+LC *)
      If[i>=gtriggerT+LC && Total@tmp > threshold,
      	ngam = grate*First[tmp];
      	tmp = ReplacePart[tmp,1->(1.0-grate)*First@tmp];
      	AppendTo[thresholdT,i];
      	,
      	ngam = 0;
      ];	
     	
	];
	
     (* list ot total parasites; (time) past --> future  *)     
     AppendTo[ls, tmp];     
     
     (* list of total gametocytes; (time) young --> old *)
     PrependTo[gls,ngam];
         
     If[drugAct[[2]]==1,
     	(*drug action for gametocyte*)
     	gls = drugActionSexual[gls,gkillzone,drugeffS,i];
     ];
     (* immune system kills sequestered gametocytes *)
  	 (*gls = ImmuneKillSequesteredX[gls,rho,delta,kappa,{ag1,ag2},i];*)
     (* immune system kills sequestered gametocytes *)
     (*gls = KillSequestered[gls,seqK3,seqK4,{seqAge3,seqAge4},i];*)
     Which[
     	immuneFunc==1,gls = ImmuneKillSequesteredX[gls,rho,delta,kappa,{ag1,ag2},i],
     	immuneFunc==2, gls = KillSequestered[gls,seqK3,seqK4,{seqAge3,seqAge4},i]
     ];
     		
     		 		
     (* reset ngam; why!! *)
     ngam = 0;
     
     (* look like rings *)
     AppendTo[nring2,CountObservedGame[gls,{1,ag1},k]];
	 
	 (* young observed gametocytes *)
	 AppendTo[ygls,CountObservedGame[gls,{1,ag1},(1-k)]];
	
	 (* old observed gametocytes ; from ag2 to ag4 *)
	 AppendTo[ogls,CountObservedGame[gls,{ag2+1,ag4},1]];
	
	 (* infectiousness gametocytes *)
	 AppendTo[infgls,CountObservedGame[gls,{ag3+1,ag4},1]];
	 
     , {i, 1, runtime}]; (* Do loop *)

   (* observed rings (no gametocytes) *)
   nring = CountRing[#,{11,14},1]&/@ls;
   
   (* replace the approximate numbers that are close to zero (10^-10) by zero *)
   obasex = Chop@(nring + nring2);
   obsex = Chop@(ygls + ogls);
      
 (*** output forms ***)
   Which[
   	(* show age distribution over time*)
   	outform==0,ls  
    (* show the last age distribution*)
   ,outform==1,tmp
    (* show the number of gametocyte over time*)
   ,outform==2,gls
    (* big results! *)
   ,outform==3,{ls,gls,nring,obasex,obsex,infgls,thresholdT}
    (* threshold time *)
   ,outform==4,thresholdT 
   ]

];



(** this function has never be used! **)
(* decay function for PMR   *)
(*pmrFn[pmr0_,pmrmin_,k_,t_]:=(pmr0-pmrmin)*Exp[-k*t]+pmrmin*)

(* the gametocytes are generated if the total number of asexual parasites is below a threshold 
this function will return 0 if the difference of the total number at the last and the previous is less than the threshold
and 1 if the difference greater than the threshold.  
*)
stressDifN[totls_List,stressthreshold_,position4compare_Integer]:= Module[{last, comp, output, difN, inputlength},
  inputlength = Length@totls;
  last = Last@totls;
  output=0; (* default *)
  If[position4compare < inputlength,
   comp = totls[[inputlength - position4compare]];
   difN = comp - last;,
   difN = 0;

   If[difN > stressthreshold,
    output = 1, output = 0
   ];

   ];
  
  output
];


End[] (* End Private Context *)

EndPackage[]