#!/bin/bash

# by Sisi Liu: sisi.liu@awi.de; sisi.liu.research@gmail.com
#== HPC setting: dependent ==

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


#== START ==
INDIR="/PARENT/naleng/ngs/out.fastp"
OUTDIR="/PARENT/naleng/ngs/out.fastp.extract"
DBP="/PARENT/naleng/mapdamage"

# input fastq.gz
END_MERGED="_fastp_merged_R2.fq.gz"

cd ${INDIR}
FILE_R1=$(ls *${END_MERGED} | sed -n ${SLURM_ARRAY_TASK_ID}p)
FILEBASE=${FILE_R1%${END_MERGED}}
OUT_MERGED="${FILEBASE}_fastp_merged_R2.fq.gz"

for ti in ${DBP}/*.txt; do\
	bname=$(basename "$ti")
	IFS="_" read -ra parts <<< "$bname"
	TAXID="${parts[0]}"
	ORG="${parts[1]}"
	RANK="${parts[3]}"
	SEQID="${TAXID}_${ORG}_${RANK}"
	mkdir -p ${OUTDIR}/${ORG}
	SEQOUT="${OUTDIR}/${ORG}/${FILEBASE}_${ORG}.fq.gz"
	srun /PARENT/naleng/programmes/seqtk/seqtk subseq ${INDIR}/${OUT_MERGED} ${DBP}/${SEQID} > ${SEQOUT}
done
#== END ==
