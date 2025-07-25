#! /bin/bash

T=0.1
h=0.1
n_step=10000000
snap_dt=100000
seed=3640
IC=ordered

for sigma in 0.2
do
    job_name="KM_L4096_T${T}_s${sigma}_h${h}"
    cat <<EOF > ${job_name}.sh
#!/bin/bash
#SBAtCH --job-name=${job_name}
#SBATCH --partition=localhost
#SBATCH -N 1
#SBATCH --ntasks=64
#SBATCH --hint=nomultithread

mpirun --bind-to none ./a8_8.out $T $sigma $h $seed $snap_dt $n_step $IC
exit 0
EOF
sleep 0.25
sbatch ${job_name}.sh
done
