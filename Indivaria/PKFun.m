(* Mathematica Package *)

BeginPackage["Indivaria`PKFun`"]
(* Exported symbols added here with SymbolName::usage *)  
DrugConc1::usage = "DrugConc1[C0_, k_, t_] := C0*E^(-k*t)"


Begin["`Private`"] (* Begin Private Context *) 

(**)
DrugConc1[C0_, k_, t_] := C0*E^(-k*t)

(**)
DrugConc2[Bio_, dose_, Vd_, Kab_, Kel_, t_] := ((Bio*dose)/Vd)*(Kab/(Kab - Kel))*(E^(-Kel*t) - E^(-Kab*t))


(**time to maximum concentration*)
Tmax[Kab_, Kel_] := (Log[Kel] - Log[Kab])/(Kel - Kab)



End[] (* End Private Context *)

EndPackage[]