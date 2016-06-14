#PBS -q class
#PBS -l nodes=6:sixcore
#PBS -l walltime=02:00:00
#PBS -N test

export PBS_O_WORKDIR=$HOME/prog2
EXE=$PBS_O_WORKDIR/jacobi

for n in 2000 4000 6000 8000 10000 12000
do
echo "----------Parallel--n-fixed---------" >&2
echo "n=$n" >&2
for p in 1 4 9 16 25 36 49
do
echo "----------Parallel--p-fixed---------" >&2
echo "p=$p" >&2
for d in 0.2 0.4 0.5 0.6 0.8 1.0
do
echo "----------Parallel--d-varying-------" >&2
echo "d=$d" >&2
OMPI_MCA_mpi_yield_when_idle=0 mpirun --hostfile $PBS_NODEFILE -np $p $EXE  -n $n -d $d
done
done
done
