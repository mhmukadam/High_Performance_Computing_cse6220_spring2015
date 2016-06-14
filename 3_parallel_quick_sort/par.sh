#PBS -q class
#PBS -l nodes=jinx12+jinx13+jinx21+jinx24
#PBS -l walltime=02:00:00
#PBS -N test

export PBS_O_WORKDIR=$HOME/prog3
EXE=$PBS_O_WORKDIR/sort

for n in 50000 100000 250000 500000 750000 1000000
do
echo "----------Parallel--n-fixed---------" >&2
echo "n=$n" >&2
for p in 1 3 6 9 12 15 18 21 24
do
echo "----------Parallel--p-fixed---------" >&2
echo "p=$p" >&2
OMPI_MCA_mpi_yield_when_idle=0 mpirun --hostfile $PBS_NODEFILE -np $p $EXE -r -n $n
done
done
