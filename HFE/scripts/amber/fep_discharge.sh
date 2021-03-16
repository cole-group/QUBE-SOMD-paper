# a util script to run all of the discharge stages.
# the lambda values are hard set to match the cfg

for i in 0.0000 0.0625 0.1250 0.1875 0.2500 0.3125 0.3750 0.4375 0.5000 0.5625 0.6250 0.6875 0.7500 0.8125 0.8750 0.9375 1; do
mkdir lambda-$i
cd lambda-$i
# move the files from the input
cp ../../input/MORPH.pert.discharge MORPH.pert
cp ../../input/sim.cfg .
cp ../../input/SYSTEM.* .

somd-freenrg -C sim.cfg -t SYSTEM.top -c SYSTEM.crd -m MORPH.pert -l $i -p CUDA
cd ..
done
