# a master script to run all of the windows in the somd simulation
# this will start with vacuum then solvated legs

echo "STARTING SOMD RUN"
echo "STARTING VACUUM PHASE"
for i in vacuum free; do
cd $i
# first run dis
cd discharge
bash fep_discharge.sh
cd ../vanish
bash fep_vanish.sh
cd ../../
done