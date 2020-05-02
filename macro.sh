#!/bin/bash
# este es un script que genera un macro para geant4 avanzado
#clear
#T="$(date +%s)"
$n = 0
for ((x1 = 100 ; x1 <= 300 ; x1+=1)); 
do
	for ((x2 = 100 ; x2 <= 300 ; x2+=1)); 
	do

	r1=$(echo "scale=2; $x1/10.0" | bc)
	r2=$(echo "scale=2; $x2/10.0" | bc)

	##  imprimimos el macro
	echo -e "/control/verbose 0 \n/run/verbose 0 \n/vis/disable \n/det/field 1 \n/det/target/material G4_lH2 \n/det/target/Z 1 \n/det/target/A 2 \n/det/target/width 0.001 \n/det/target/pos 0. 0. 301. \n/det/primary/energy 20 MeV \n/det/primary/Z 5 \n/det/primary/A 8 \n/det/primary/pos 0. 0. -115. \n/det/recoil/A 3 \n/det/recoil/Z 2 \n/det/ejectile/A 7 \n/det/ejectile/Z 4" "\n/det/currentValue $r1" "\n/det/currentValue2 $r2" "\n/det/update" "\n/run/beamOn 100" > tmp.mac

	## imprimimos un temporal de los datos usados en la simulacion
	echo $r1 $r2 > tmp.dat 

	##  rodamos la simulacion
	./arquivobinario tmp.mac

	## limpiamos todo
	mv tree_run_0.root tree_run_$r2.root
	mv tree_run_$r2.root ROOT/
	
	done
done
#T="$(($(date +%s)-T))"
#echo "Time in seconds: ${T}"
