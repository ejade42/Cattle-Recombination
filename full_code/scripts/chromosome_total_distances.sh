# Makes a file in inputs with the total genetic distance of each chromosome.
# Requires the .markers files in each chromosome folder (i.e. must be run after LinkphaseInput)
# Runs from the scripts directory (i.e. can be executed directly from here)

cd ../qc
echo > ../files/chromosome_total_distances.txt
for i in {1..29}
do
    cd chr_$i
    tail -1 chr_$i".markers" >> ../../files/chromosome_total_distances.txt
    cd ..
done