(* Mathematica Package *)
(* the collection of basics functions required by subpackages of Indivaria*)

BeginPackage["Indivaria`BasicFunc`"]
(* Exported symbols added here with SymbolName::usage *)  

DistributeN::usage="DistributeN[initN,lifecycle, mu, sigma] generates the age distribution of the parasites."
DistributeN2::usage="DistributeN2[initN,lifecycle, mu, sigma, pmr] generates the age distribution of the parasites."

Shiftonehour::usage="Shiftonehour[ls,pmr] rotate the list of numbers to the right 1 step (an hour) with a multiplication rate pmr."
Eff::usage="Eff[c, gamma, ec50, emin, emax] calculates the efficacy."
PRingFunc::usage="PRingFunc[i, a1, a2]"
PSequestered::usage="PSequestered[i, a1, a2]"
CountRing::usage="CountRing[ls, Age, fn]"
LsDot::usage="LsDot[ls1, ls2]"
EffConc::usage="EffConc[Ndrug,everyH,efffn,efffnparms,runtime]" 
EffConc1::usage="EffConc1[Ndrug, everyH, efffnparms, runtime]"
KillSequestered::usage="KillSequestered[parls,k1,k2,agerange_List,t]"
ImmuneKillSequestered::usage="ImmuneKillSequestered[parls,rho,delta,kappa,killrange,t]"
ImmuneKillSequesteredAX::usage="ImmuneKillSequesteredAX[parls,rho,delta,kappa,t]"
ImmuneKillSequesteredX::usage="ImmuneKillSequesteredX[parls,rho,delta,kappa,killrange,t]"

Begin["`Private`"] (* Begin Private Context *) 

(* original version that was used in PNAS *)
DistributeN[initN_, LC_Integer, mu_, sigma_] :=
(*DistributeN[initN,LifeCycle,mu,sigma] distributes the initial number over the life-cycle (LC)
with the mean mu and SD sigma hours.*)
  Module[{distr, x},
   distr = Table[PDF[NormalDistribution[mu, sigma], x]//N, {x, 1, LC}];
   Developer`ToPackedArray[Chop[(initN/Total[distr])*distr]]
];

(*based on NJW*)
DistributeN2[initN_, LC_Integer, mu_, sigma_,pmr_] :=
(*DistributeN2[initN,LifeCycle,mu,sigma,pmr] distributes the initial number over the life-cycle (LC)
with the mean mu and SD sigma hours.*)
  Module[{distr, x},
   distr = Table[PDF[NormalDistribution[mu, sigma],x-LC]/pmr+PDF[NormalDistribution[mu, sigma], x]+(PDF[NormalDistribution[mu, sigma],x+LC]*pmr), {x, 1, LC}]//N;
   Developer`ToPackedArray[Chop[(initN/Total[distr])*distr]]
];
   

Shiftonehour[ls_List, PMR_] := Module[{tmp},
   (*shift the list of number to the right with a multiplication pmr*)
   tmp = ls;
   tmp = RotateRight[tmp];
   tmp = ReplacePart[tmp, 1 -> tmp[[1]]*PMR];
   tmp
 ];


(*calculate the percentage for killing*)
Eff[c_, gamma_, ec50_, emin_, emax_] :=  emin + (emax - emin) *c^gamma/(c^gamma + ec50^gamma)

(* inferred from Fig2 of "Febrile temperatures induce cytoadherence of 
	ring-stage Plasmodium falciparum-infected erythrocytes" by
	Rachanee Udomsangpetch et al. (http://www.pnas.org/content/99/18/11825)
	
	here it's estimated that 50% circulating
	*)
PRingFunc[i_Integer, a1_Integer, a2_Integer] := 
  Piecewise[{{1, i < a1}, {Exp[Log[0.5]*(i - a1)/(a2 - a1)], i >= a1}}];

(* probability of becoming sequestered *)
PSequestered[i_Integer,a1_Integer,a2_Integer]:=1.0-PRingFunc[i, a1, a2]

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


LsDot[ls1_List, ls2_List] := Module[{tmp},
   (* dot vector *)
   If[Length@ls1 == Length@ls2, 
   		tmp = ls1*ls2;
   		, 
    	(*Print["Please check your input."];*)
    	tmp = {};
    ];
   tmp
   ];

(* generates the list of the drug concentration and its efficacy over time following the dose regimen  
   ndrug = number of the drugs; everh = time for taking next dose; fn = the efficacy function; fnparms = the list of parameters of fn 

	efffn = 0 , efffnparms = {"dataname",gamma_List,ec50_List,emin_List,emax_List}
	efffn = 1 , efffnparms = {gamma_List,ec50_List,emin_List,emax_List, concparm_List} , concparms = {xm,ym,ke}
    efffn = 2 , efffnparms = {gamma,ec50,emin,emax, concparm_List}, concparms = {xm,ym,ke}
*)
EffConc[Ndrug_Integer,everyH_Integer,efffn_Integer,efffnparms_List,runtime_Integer] := 
Module[{concls,nd,i,j,conc,t,fn,dat,datname,maxpoint,endtime,part1,part2,modka,fitke,modke,a,ke,emax,emin,gamma,ec50,
	TIMEMAX,xm,ym,x},

   TIMEMAX = runtime;
   j = 0;
   nd = 1;
   concls = {};
   
   conc = 0; (*init concentration*)
 
 Which[efffn==0, (* PNAS three stages EC-Curves *)
   (* efffnparms = {"dataname",gamma_List,ec50_List,emin_List,emax_List} *)
   datname = efffnparms[[1]];
   gamma = efffnparms[[2]];
   ec50 = efffnparms[[3]];
   emin = efffnparms[[4]];
   emax = efffnparms[[5]];
   
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
   
   
   For[i = 0, i < Ndrug*everyH, i = i + 1,
    If[j != everyH,
     conc = Chop@(fn /. t -> j);
     
     AppendTo[concls, {i, conc,
       Eff[conc, gamma[[1]], ec50[[1]], emin[[1]], emax[[1]]], 
       Eff[conc, gamma[[2]], ec50[[2]], emin[[2]], emax[[2]]], 
       Eff[conc, gamma[[3]], ec50[[3]], emin[[3]], emax[[3]]]}];
     ];
    
    If[j == everyH && nd < Ndrug,
     nd = nd + 1;
     conc = Chop@(conc + fn /. t -> j);
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
    conc = Chop@(fn /. t -> j);
    AppendTo[concls, {i, conc,
      Eff[conc, gamma[[1]], ec50[[1]], emin[[1]], emax[[1]]], 
      Eff[conc, gamma[[2]], ec50[[2]], emin[[2]], emax[[2]]], 
      Eff[conc, gamma[[3]], ec50[[3]], emin[[3]], emax[[3]]]}
     ];
    j = j + 1;
    ];
  
  ,efffn==1, (* three stages without conc data; use straight line and exponential decay for absorbtion and eliminatio phases.*)
  (* efffnparms = {gamma_List,ec50_List,emin_List,emax_List, concparm_List} 
     concparms = {xm,ym,ke}
  *)
   gamma = efffnparms[[1]];
   ec50 = efffnparms[[2]];
   emin = efffnparms[[3]];
   emax = efffnparms[[4]];
   
   xm = efffnparms[[5,1]];
   ym = efffnparms[[5,2]];
   ke = efffnparms[[5,3]];
    
  For[i = 0, i < Ndrug*everyH, 
  	i = i + 1, 
 	If[j != everyH, conc = Chop@(genpoints[xm, ym, ke, x] /. x -> j);
 	AppendTo[concls, {i, conc, 
    Eff[conc, gamma[[1]], ec50[[1]], emin[[1]], emax[[1]]], 
    Eff[conc, gamma[[2]], ec50[[2]], emin[[2]], emax[[2]]], 
    Eff[conc, gamma[[3]], ec50[[3]], emin[[3]], emax[[3]]]}];
     ];
  
 	If[j == everyH && nd < Ndrug, nd = nd + 1;
  	 conc = Chop@(conc + genpoints[xm, ym, ke, x] /. x -> j);
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
 		conc = Chop@(genpoints[xm, ym, ke, x] /. x -> j);
 		AppendTo[concls, {i, conc, 
   		Eff[conc, gamma[[1]], ec50[[1]], emin[[1]], emax[[1]]], 
   		Eff[conc, gamma[[2]], ec50[[2]], emin[[2]], emax[[2]]], 
   		Eff[conc, gamma[[3]], ec50[[3]], emin[[3]], emax[[3]]]}];
 		j = j + 1;
 	];
  
  ,efffn==2,  (* one killzone *)
    (* efffnparms = {gamma_,ec50_,emin_,emax_, concparm_List} 
     concparms = {xm,ym,ke}
  *)
   gamma = efffnparms[[1]];
   ec50 = efffnparms[[2]];
   emin = efffnparms[[3]];
   emax = efffnparms[[4]];
   
   xm = efffnparms[[5,1]];
   ym = efffnparms[[5,2]];
   ke = efffnparms[[5,3]];
    
  For[i = 0, i < Ndrug*everyH, 
  	i = i + 1, 
 	If[j != everyH, conc = Chop@(genpoints[xm, ym, ke, x] /. x -> j);
 	AppendTo[concls, {i, conc, Eff[conc, gamma, ec50, emin, emax]}];
     ];
  
 	If[j == everyH && nd < Ndrug, nd = nd + 1;
  	 conc = Chop@(conc + genpoints[xm, ym, ke, x] /. x -> j);
  	 AppendTo[concls, {i, conc, Eff[conc, gamma, ec50, emin, emax]}];
      j = 0;
     ];
 
 (*update ndrug*)
 	 j = j + 1;
 	 ];

    For[i = Ndrug*everyH; j = everyH, i <= TIMEMAX, i = i + 1, 
 		conc = Chop@(genpoints[xm, ym, ke, x] /. x -> j);
 		AppendTo[concls, {i, conc, Eff[conc, gamma, ec50, emin, emax]}];
 		j = j + 1;
 	];   
    
     
  ];
  
   concls
   


];

(***)

(* drug concentration and efficacy for one kill zone *)
EffConc1[Ndrug_Integer, everyH_Integer, efffnparms_List, runtime_Integer] := 
  Module[{nd, i, j, conc, ke, emax, emin, gamma, ec50, TIMEMAX, xm, 
    ym, x, concls},(*efffnparms={gamma_,ec50_,emin_,emax_,
   concparm_List} concparms={xm,ym,ke}*)TIMEMAX = runtime;
   j = 0;
   nd = 1;
   
   conc = 0;
   gamma = efffnparms[[1]];
   ec50 = efffnparms[[2]];
   emin = efffnparms[[3]];
   emax = efffnparms[[4]];
   xm = efffnparms[[5, 1]];
   ym = efffnparms[[5, 2]];
   ke = efffnparms[[5, 3]];
   concls = Reap[
     For[i = 0, i < Ndrug*everyH, i = i + 1, 
      If[j != everyH, conc = Chop@(genpoints[xm, ym, ke, x] /. x -> j);
       Sow[{i, conc, Eff[conc, gamma, ec50, emin, emax]}];
       ];
      If[j == everyH && nd < Ndrug, nd = nd + 1;
       conc = Chop@(conc + genpoints[xm, ym, ke, x] /. x -> j);
       Sow[{i, conc, Chop@Eff[conc, gamma, ec50, emin, emax]}];
       
       j = 0;];
      (*update ndrug*)j = j + 1;];
     For[i = Ndrug*everyH; j = everyH, i <= TIMEMAX, i = i + 1, 
      conc = Chop@(genpoints[xm, ym, ke, x] /. x -> j);
      Sow[{i, conc, Chop@Eff[conc, gamma, ec50, emin, emax]}];
      j = j + 1;];
     ];
   concls[[2, 1]]
   ];



(* concentration of DHA *)
genpoints[xm_, ym_, ke_, x_] := 
  Piecewise[{{ym - (ym/xm)*(xm - x), x <= xm}, {ym*Exp[-ke*(x - xm)], 
     x > xm}}];



(* killing function for sequestered asexual and sexual stages
  fn=k1*(1.-Exp[-k2*t]);
 *)
KillSequestered[parls_,k1_,k2_,agerange_List,t_]:=Module[{fn,tmp,beginage,endage,lengthls},
	
	lengthls=Length@parls;
	
	(* agerange={beginage,endage} *)
	{beginage,endage}=agerange;
	
	fn=k1*(1.-Exp[-k2*t]);
	tmp=parls;
	
	Which[
		beginage<=lengthls&&endage<=lengthls,
		tmp[[beginage;;endage]]=(1.-fn)*parls[[beginage;;endage]];
	,
		beginage <= lengthls &&lengthls<=endage,
		tmp[[beginage;;lengthls]] = (1.-fn)*parls[[beginage;;lengthls]];
	,
		lengthls<=beginage,
		tmp = 1.*tmp;	
	];
	
	If[#<0,0,#]&/@tmp
];


(*
immune function for killing sequestered asexual and sexual parasites
y'(t)=rho S - delta y(t); S=sequestered parasites
*)
ImmuneKillSequestered[parls_,rho_,delta_,kappa_,killrange_List,runtime_Integer,t_]:=Module[{S,tmp,lengthls,beginkill,endkill,y,i,fn},
	lengthls=Length@parls;	
	tmp=parls;

	{beginkill,endkill}=killrange;
	
	Which[
		beginkill<=lengthls&&endkill<=lengthls,
		S=Total[tmp[[beginkill;;endkill]]];
		fn=NDSolve[{y'[i] == rho*S - delta*y[i], y[0] == 0},y, {i, 0., runtime}];
		tmp[[beginkill;;endkill]]=(1.-(kappa*y[t]/S)/.fn)[[1]]*tmp[[beginkill;;endkill]];
	,
		beginkill < lengthls &&lengthls<endkill,
		S=Total[tmp[[beginkill;;lengthls]]];
		fn=NDSolve[{y'[i] == rho*S - delta*y[i], y[0] == 0},y, {i, 0., runtime}];
		tmp[[beginkill;;lengthls]] = (1.-(kappa*y[t]/S)/.fn)[[1]]*tmp[[beginkill;;lengthls]];
	,
		lengthls<=beginkill,
		tmp = 1*tmp;	
	];
	
	If[#<0,0,#]&/@tmp
	
];


(* killrange = sequestered from ag1 to ag2 *)
ImmuneKillSequesteredX[parls_List,rho_,delta_,kappa_,killrange_List,t_]:=Module[{S,tmp,lengthls,beginkill,endkill,y,n,imm,output},
	lengthls=Length@parls;	
	tmp=parls;

	{beginkill,endkill}=killrange;
	
	Which[
		beginkill<=lengthls&&endkill<=lengthls,
		S=Total[tmp[[beginkill;;endkill]]];
		imm=Last@RecurrenceTable[{y[n + 1] == y[n]+(rho*S-delta*y[n]),y[0]==.0},y,{n,1,t}];
		tmp[[beginkill;;endkill]]=tmp[[beginkill;;endkill]]-kappa*imm*tmp[[beginkill;;endkill]];
	,
		beginkill < lengthls &&lengthls<endkill,
		S=Total[tmp[[beginkill;;lengthls]]];
		imm=Last@RecurrenceTable[{y[n + 1] == y[n]+(rho*S-delta*y[n]),y[0]==0.},y,{n,1,t}];
		tmp[[beginkill;;lengthls]] = tmp[[beginkill;;lengthls]]-kappa*imm*tmp[[beginkill;;lengthls]];
	,
		lengthls<=beginkill,
		tmp = 1*tmp;	
	];
	
	(* it the number is negative then replace it with 0 *)
	output=If[#<0,0,#]&/@tmp;
	
	output
];


(* y[n+1]=y[n]+(rho S[n] - delta y[n]) 
S = total@sequestered
sequestered = PSequestered * pardist
*)
ImmuneKillSequesteredAX[parls_,rho_,delta_,kappa_,t_]:=Module[{S,tmp,lengthls,y,i,imm,n,pseq,output},
	lengthls=Length@parls;	
	tmp=parls;
	pseq=Table[PSequestered[i,11,14],{i,1,Length@parls}];
	
	S=pseq.tmp;
	imm=Last@RecurrenceTable[{y[n + 1] == y[n]+(rho*S-delta*y[n]),y[0]==0},y,{n,1,t}];
	
	tmp=tmp-kappa*imm*tmp*pseq;
	output=If[#<0,0,#]&/@tmp;
	
	output
	
];



End[] (* End Private Context *)

EndPackage[]