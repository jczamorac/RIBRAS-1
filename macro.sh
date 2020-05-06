T="$(date +%s)"

for((x3 = 900 ; x3 <= 1150 ; x3+=50));	
do
	for ((x1 = 300 ; x1 <= 500 ; x1+=5)); 
	do
		for ((x2 = 300 ; x2 <= 500 ; x2+=5)); 
		do
			r1=$(echo "scale=2; $x1/10.0" | bc)
			r2=$(echo "scale=2; $x2/10.0" | bc)
			r3=$(echo "scale=2; $x3/10.0" | bc)

			## Create vis.mac file
			echo -e "/control/verbose 0 \n/run/verbose 0 \n/vis/disable \n/det/field 1 \n/det/target/material G4_lH2 \n/det/target/Z 1 \n/det/target/A 2 \n/det/target/width 0.001 \n/det/target/pos 0. 0. 301. \n/det/primary/energy 30 MeV \n/det/primary/Z 5 \n/det/primary/A 8 \n/det/primary/pos 0. 0. -$r3 \n/det/recoil/A 3 \n/det/recoil/Z 2 \n/det/ejectile/A 7 \n/det/ejectile/Z 4" "\n/det/currentValue $r1" "\n/det/currentValue2 $r2" "\n/det/update" "\n/run/beamOn 100" > tmp.mac

			## Printing current
			echo $r1 $r2 $r3 > tmp.dat 

			## Run simulation
			./arquivobinario tmp.mac

			## Moving and renaming everything
			mv tree_run_0.root tree_run_{$r1}_{$r2}_{$r3}.root
			mv tree_run_{$r1}_{$r2}_{$r3}.root ROOT/
		done
	done
done

T="$(($(date +%s)-T))"
echo "Time in seconds: ${T}"
