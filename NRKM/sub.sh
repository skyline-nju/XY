#! /bin/bash

T=0.1

sigma=0.2
h=0.1
n_step=2000000
snap_dt=10000
IC=rand

for seed in 1002 1003 1004 1005 1006 1007 1008 1009
do
    job_name="KM512_s${sigma}_T${T}_h${h}_S${seed}"
    cat <<EOF > ${job_name}.sh
#!/bin/bash
#SBAtCH --job-name=${job_name}
#SBATCH --partition=localhost
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --hint=nomultithread

./a.out $T $sigma $h $seed $snap_dt $n_step $IC
exit 0
EOF
sleep 0.25
sbatch ${job_name}.sh
done
