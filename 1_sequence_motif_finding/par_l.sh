#PBS -q class
#PBS -l nodes=4:sixcore
#PBS -l walltime=02:00:00
#PBS -N 10_32_l

export PBS_O_WORKDIR=$HOME/prog1
EXE=$PBS_O_WORKDIR/findmotifs
MASTER_DEPTH=4
p=6

echo "----------Parallel-----l-varied-----" >&2
echo "10_32_5" >&2
INPUT_FILE=$PBS_O_WORKDIR/input/10_32_5.txt
OUTPUT_FILE=$PBS_O_WORKDIR/result/10_32_5.txt
OMPI_MCA_mpi_yield_when_idle=0 mpirun --hostfile $PBS_NODEFILE -np $p $EXE $MASTER_DEPTH $INPUT_FILE $OUTPUT_FILE

echo "10_32_6" >&2
INPUT_FILE=$PBS_O_WORKDIR/input/10_32_6.txt
OUTPUT_FILE=$PBS_O_WORKDIR/result/10_32_6.txt
OMPI_MCA_mpi_yield_when_idle=0 mpirun --hostfile $PBS_NODEFILE -np $p $EXE $MASTER_DEPTH $INPUT_FILE $OUTPUT_FILE

echo "10_32_7" >&2
INPUT_FILE=$PBS_O_WORKDIR/input/10_32_7.txt
OUTPUT_FILE=$PBS_O_WORKDIR/result/10_32_7.txt
OMPI_MCA_mpi_yield_when_idle=0 mpirun --hostfile $PBS_NODEFILE -np $p $EXE $MASTER_DEPTH $INPUT_FILE $OUTPUT_FILE

echo "10_32_8" >&2
INPUT_FILE=$PBS_O_WORKDIR/input/10_32_8.txt
OUTPUT_FILE=$PBS_O_WORKDIR/result/10_32_8.txt
OMPI_MCA_mpi_yield_when_idle=0 mpirun --hostfile $PBS_NODEFILE -np $p $EXE $MASTER_DEPTH $INPUT_FILE $OUTPUT_FILE

echo "10_32_9" >&2
INPUT_FILE=$PBS_O_WORKDIR/input/10_32_9.txt
OUTPUT_FILE=$PBS_O_WORKDIR/result/10_32_9.txt
OMPI_MCA_mpi_yield_when_idle=0 mpirun --hostfile $PBS_NODEFILE -np $p $EXE $MASTER_DEPTH $INPUT_FILE $OUTPUT_FILE

echo "10_32_10" >&2
INPUT_FILE=$PBS_O_WORKDIR/input/10_32_10.txt
OUTPUT_FILE=$PBS_O_WORKDIR/result/10_32_10.txt
OMPI_MCA_mpi_yield_when_idle=0 mpirun --hostfile $PBS_NODEFILE -np $p $EXE $MASTER_DEPTH $INPUT_FILE $OUTPUT_FILE
