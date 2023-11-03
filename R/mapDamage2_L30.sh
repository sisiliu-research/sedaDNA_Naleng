#!/bin/bash
#
# contact: sisi.liu@awi.de
#
# slurm options and variables under >set required variables< 
#================== HPC setting: user dependent ==============

#SBATCH --account=
#SBATCH --job-name=
#SBATCH --partition=
#SBATCH -t 00:00:00
#SBATCH --qos=
#SBATCH --array=
#SBATCH --cpus-per-task=
#SBATCH --mem=
#SBATCH --mail-type=
#SBATCH --mail-user=


# set required variables (adapt according to your own requirements)
#===================================================================
CPU=${SLURM_CPUS_PER_TASK}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

BWA="bwa/0.7.17"
VIEW="samtools/1.16.1"
MAPD="mapdamage2/2.2.1"

#===================================================================
PARENT="/path"
WORK="${PARENT}/naleng"
INDIR="${PARENT}/naleng/ngs/out.fastp.extract"
DBP="${PARENT}/naleng/mapdamage/fasta"

cd ${DBP}
ACC=$(basename $(find . -mindepth 1 -maxdepth 1 -type d | sed -n ${SLURM_ARRAY_TASK_ID}p))

tax2acc_file="${PARENT}/naleng/mapdamage/tax2acc2.txt"
DB=$(awk -v org="$ACC" '$1 == org {print $1}' "$tax2acc_file")
ORG=$(awk -v org="$ACC" '$1 == org {print $2}' "$tax2acc_file")

DBi="${DBP}/${DB}"
DIR_PATH="${INDIR}/${ORG}"
FILE_END="_${ORG}.fq.gz"

cd ${DIR_PATH}
INFILE=$(ls *${FILE_END})

# Load required modules
module load ${BWA}
module load ${VIEW}
module load ${MAPD}

echo ${DB}
echo ${ORG}
echo ${DBi}
echo ${DIR_PATH}

# Loop through FASTQ files in the ORG directory
for i in ${DIR_PATH}/${INFILE}; do
	echo $i
	bname=$(basename $i)
	FILEBASE=${bname%$FILE_END}
	
	echo "Aligning $i $FILEBASE"
	srun bwa mem ${DBi}/${DB}.fasta ${i} -t ${SLURM_CPUS_PER_TASK} > ${DBi}/${FILEBASE}.sam
	srun samtools view -Sb ${DBi}/${FILEBASE}.sam > ${DBi}/${FILEBASE}.bam
	
	echo "Sorting BAM for $FILEBASE"
	srun samtools sort ${DBi}/${FILEBASE}.bam -o ${DBi}/${FILEBASE}.sort.bam
	
	echo "Index sorted BAM file"
	srun samtools index ${DBi}/${FILEBASE}.sort.bam
	
	echo "Run mapDamage"
	srun mapDamage -i ${DBi}/${FILEBASE}.sort.bam -r ${DBi}/${DB}.fasta --rescale --single-stranded --rescale-out=${DBi}/${FILEBASE}.rescale.bam -d ${DBi}/${FILEBASE}
done

# Unload modules after processing
module unload ${BWA}
module unload ${VIEW}
module unload ${MAPD}
#== END ==
