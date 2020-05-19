T="$(date +%s)"

# Delete every .root file first
cd ROOT
rm *.root
cd ..

# Initial value for current 1 and 2 (int value)
iCurrent1=240
iCurrent2=330

# Final value for current 1 and 2 (int value)
fCurrent1=250
fCurrent2=340

# Primary beam position
Pos=-80

# Changing analyse macro
sed -i "s/^    for(Current1.*/    for(Current1 = $((iCurrent1/10)); Current1<=$((fCurrent1/10)); Current1 += 0.2)/" macros/analise3.mac
sed -i "s/^        for(Current2.*/        for(Current2 = $((iCurrent2/10)); Current2<=$((fCurrent2/10)); Current2 += 0.2)/" macros/analise3.mac

# Initializing simulation
echo '---------------------------------------'
echo 'Initializing simulation'
echo "Current 1 interval: $((iCurrent1/10))A to $((fCurrent1/10))A"
echo "Current 2 interval: $((iCurrent2/10))A to $((fCurrent2/10))A"
echo "Primary beam position: $Pos m"
echo '---------------------------------------'
sleep 5

for ((x1 = $iCurrent1 ; x1 <= $fCurrent1 ; x1+=2)); 
do
    # Getting value of current 1
    r1=$(echo "scale=2; $x1/10.0" | bc)
    r2i=$(echo "scale=2; $iCurrent2/10.0" | bc)
    r2f=$(echo "scale=2; $fCurrent2/10.0" | bc)

    # Create macro ready for loop
    echo -e "/control/verbose 0 \n/run/verbose 0 \n/vis/disable \n/det/field 1 \n/det/target/material G4_lH2 \n/det/target/Z 1 \n/det/target/A 2 \n/det/target/width 0.001 \n/det/target/pos 0. 0. 301. \n/det/primary/energy 30 MeV \n/det/primary/Z 5 \n/det/primary/A 8 \n/det/primary/pos 0. 0. $Pos \n/det/recoil/A 3 \n/det/recoil/Z 2 \n/det/ejectile/A 7 \n/det/ejectile/Z 4" "\n/det/currentValue $r1" "\n/det/currentValue2 {current}" "\n/det/update" "\n/run/beamOn 100" > tmp.mac

    # Create looping macro
    echo -e "/control/loop tmp.mac current $r2i $r2f 0.2" > looping.mac

    # Run looping macro
    ./arquivobinario looping.mac
done

# Delete everything when done
rm looping.mac
rm tmp.mac
echo 'Simulation finished'

# Then, analyse everything
echo ' '
echo 'Now, analysing...'

# Running macro
cd macros
root -l <<-EOF
.x analise3.mac
.q
EOF

cd ..

# Print simulation time
T="$(($(date +%s)-T))"
echo "Duration: ${T}s"