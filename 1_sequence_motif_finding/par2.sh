#PBS -q class
#PBS -l nodes=4:sixcore
#PBS -l walltime=02:00:00
#PBS -N 10_34_8

export PBS_O_WORKDIR=$HOME/prog1

INPUT_FILE=$PBS_O_WORKDIR/input/10_34_8.txt
OUTPUT_FILE=$PBS_O_WORKDIR/result/10_34_8.txt
EXE=$PBS_O_WORKDIR/findmotifs
MASTER_DEPTH=1
p=1

echo "---------------Serial---------------" >&2
OMPI_MCA_mpi_yield_when_idle=0 mpirun --hostfile $PBS_NODEFILE -np $p $EXE $MASTER_DEPTH $INPUT_FILE $OUTPUT_FILE

echo "----------Parallel--p-fixed---------" >&2
p=15
for MASTER_DEPTH in 3 4 5 6 7 8
do
echo "p=$p, k=$MASTER_DEPTH" >&2
OMPI_MCA_mpi_yield_when_idle=0 mpirun --hostfile $PBS_NODEFILE -np $p $EXE $MASTER_DEPTH $INPUT_FILE $OUTPUT_FILE
done

echo "----------Parallel--k-fixed---------" >&2
MASTER_DEPTH=6
for p in 3 6 9 12 15 18 21 24
do
echo "p=$p, k=$MASTER_DEPTH" >&2
OMPI_MCA_mpi_yield_when_idle=0 mpirun --hostfile $PBS_NODEFILE -np $p $EXE $MASTER_DEPTH $INPUT_FILE $OUTPUT_FILE
done