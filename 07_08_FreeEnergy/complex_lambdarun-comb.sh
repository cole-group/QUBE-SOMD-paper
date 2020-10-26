Home="$PWD"


alias somd-freenrg=/home/sire-combRuleUpdate.app/bin/somd-freenrg
export SIRE_DONT_PHONEHOME=1
source /home/Documents/smirnoff/miniconda3/etc/profile.d/conda.sh
conda activate fep

for i in $(cat $Home/pertlist_3); do 
a=${i:0:`expr index "$i" -`-1};
b=${i:`expr index "$i" -`};
c=L${a}_to_L${b};

cd $Home/Perturbations/$c/bound; 
echo $PWD;
	for s in $(seq 0.0 0.1 1.0);
	do 
		mkdir lambda-${s} ;
		cd lambda-${s} ;
		somd-freenrg -t ../$c.prm7 -c ../$c.rst7 -m ../$c.pert -C $Home/Parameters/lambda.cfg -l $s ;
		echo "		=======================
				done $c  lambda $s 	
			======================="
		cd ../ ;
	done
   
done

cd $Home
conda deactivate
