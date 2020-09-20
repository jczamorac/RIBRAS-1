T="$(date +%s)"

# Delete every .root file first
cd ROOT
rm *.root
cd ..

# Interval for Current1
iCurrent1=150
fCurrent1=180

# Interval for Current2
iCurrent2=500
fCurrent2=700

# Primary beam position (cm)
Primary_pos=-115

# Primary beam energy (MeV)
Primary_e=30

# Target definition
Target_Z=3
Target_A=6

# Primary beam particles
Primary_Z=9
Primary_A=17

# Recoil particles
Recoil_Z=1
Recoil_A=2

# Ejectile particles
Ejectile_Z=11
Ejectile_A=21

# Target z position (cm)
Target_z=0

# Changing analyse macro
sed -i "s/^    for(Current1.*/    for(Current1 = $((iCurrent1/10)); Current1<=$((fCurrent1/10)); Current1 ++)/" macros/analise4.mac
sed -i "s/^        for(Current2.*/        for(Current2 = $((iCurrent2/10)); Current2<=$((fCurrent2/10)); Current2 ++)/" macros/analise4.mac

# Initializing simulation
echo 'Initializing simulation:'
echo ' '
echo "Current 1 interval: $((iCurrent1/10))A to $((fCurrent1/10))A"
echo "Current 2 interval: $((iCurrent2/10))A to $((fCurrent2/10))A"
echo "Primary beam position: $Primary_pos m"
echo "Primary beam energy: $Primary_e"
echo "Target position: $Target_z m"                            
sleep 5

# Running simulation with a interval of 1.0
for ((x1 = $iCurrent1 ; x1 <= $fCurrent1 ; x1+=10)); 
do
    # Getting value of current 1
    r1=$(echo "scale=2; $x1/10.0" | bc)
    r2i=$(echo "scale=2; $iCurrent2/10.0" | bc)
    r2f=$(echo "scale=2; $fCurrent2/10.0" | bc)

    # Create macro ready for loop
    echo -e "/control/verbose 0 \n/run/verbose 0 \n/vis/disable \n/det/field 1 \n/det/target/material G4_POLYETHYLENE \n/det/target/Z $Target_Z \n/det/target/A $Target_A \n/det/target/width 1.89e-3 \n/det/target/pos 0. 0. $Target_z \n/det/primary/energy $Primary_e MeV \n/det/primary/Z $Primary_Z \n/det/primary/A $Primary_A \n/det/primary/pos 0. 0. $Primary_pos \n/det/recoil/A $Recoil_A \n/det/recoil/Z $Recoil_Z \n/det/ejectile/A $Ejectile_A \n/det/ejectile/Z $Ejectile_Z" "\n/det/currentValue $r1" "\n/det/currentValue2 {current}" "\n/det/update" "\n/run/beamOn 100" > tmp.mac

    # Create looping macro
    echo -e "/control/loop tmp.mac current $r2i $r2f 1.0" > looping.mac

    # Run looping macro
    ./arquivobinario looping.mac
done

# Analyse everything
cd macros
root -l <<-EOF
.x analise4.mac
.q
EOF

# Now, getting currents
lines=(`cat "Current.txt"`)
fCurrent1=${lines[0]}
fCurrent2=${lines[1]}
iCurrent1=${lines[2]}
iCurrent2=${lines[3]}

if [ $fCurrent1 -ne -1000 ]
then
    # Delete everything for the next run
    cd ..
    cd ROOT
    rm *.root
    cd ..

    # Changing analyse macro
    sed -i "s/^    for(Current1.*/    for(Current1 = $((iCurrent1/10)); Current1<=$((fCurrent1/10)); Current1 += 0.2)/" macros/analise3.mac
    sed -i "s/^        for(Current2.*/        for(Current2 = $((iCurrent2/10)); Current2<=$((fCurrent2/10)); Current2 += 0.2)/" macros/analise3.mac

    # Running simulation with a interval of 0.2
    for ((x1 = $iCurrent1 ; x1 <= $fCurrent1 ; x1+=2)); 
    do
        # Getting value of current 1
        r1=$(echo "scale=2; $x1/10.0" | bc)
        r2i=$(echo "scale=2; $iCurrent2/10.0" | bc)
        r2f=$(echo "scale=2; $fCurrent2/10.0" | bc)

        # Create macro ready for loop
        echo -e "/control/verbose 0 \n/run/verbose 0 \n/vis/disable \n/det/field 1 \n/det/target/material G4_POLYETHYLENE \n/det/target/Z $Target_Z \n/det/target/A $Target_A \n/det/target/width 1.89e-3 \n/det/target/pos 0. 0. $Target_z \n/det/primary/energy $Primary_e MeV \n/det/primary/Z $Primary_Z \n/det/primary/A $Primary_A \n/det/primary/pos 0. 0. $Primary_pos \n/det/recoil/A $Recoil_A \n/det/recoil/Z $Recoil_Z \n/det/ejectile/A $Ejectile_A \n/det/ejectile/Z $Ejectile_Z" "\n/det/currentValue $r1" "\n/det/currentValue2 {current}" "\n/det/update" "\n/run/beamOn 100" > tmp.mac

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
    cd ROOT
    rm *.root
    cd ..
else
echo " "
echo "-----------------------------"
echo "-- Couldn't find any focus --"
echo "-----------------------------"
echo " "
cd ..
cd ROOT
rm *.root
cd ..
fi

# Print simulation time
T="$(($(date +%s)-T))"
echo "Duration: $T s"