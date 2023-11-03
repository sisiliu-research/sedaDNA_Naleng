# sedaDNA analysis of Lake Naleng on the Tibetan plateau
This repository contains source codes used to analysis shotgun sequencing/metagenomic data of Liu et al., 2023. [Tibetan terrestrial and aquatic ecosystems collapsed with cryosphere loss inferred from sedimentary ancient metagenomics] 

The codes for bioinformatics are done by Sisi Liu (sisi.liu@awi.de/sisi.liu.research@gmail.com) and Lars Harms (lars.harms@awi.de). Other scripts are done by Sisi Liu.

## Bioinformatic analysis

0. Slrum option setting: [Slrum manual](https://slurm.schedmd.com/sbatch.html)

1. Quality check on raw sequencing data (*_R1.fastq.gz and *_R2.fastq.gz): [fastqc manual](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
```
srun fastqc -q -o ${OUTDIR}/${OUT_FQC}/${PRE} -t 2 ${INDIR}/${FILE_R1} ${INDIR}/${FILE_R2}

```
2. Remove duplicate raw reads: [clumpify.sh](https://github.com/BioInfoTools/BBMap/blob/master/sh/clumpify.sh)
```
srun clumpify.sh in=${INDIR}/${FILE_R1} in2=${INDIR}/${FILE_R2} out=${OUTDIR}/${OUT_DEDUP}/${OUT_R1_CL} out2=${OUTDIR}/${OUT_DEDUP}/${OUT_R2_CL} dedupe=t

```
3. Adapter trimming and merging of paired-end reads in parallel: [fastp manual](https://github.com/OpenGene/fastp#merge-pe-reads)
```
srun fastp --in1 ${OUTDIR}/${OUT_DEDUP}/${OUT_R1_CL} --in2 ${OUTDIR}/${OUT_DEDUP}/${OUT_R2_CL} --out1 ${OUTDIR}/${OUT_FASTP}/${OUT_R1} --out2 ${OUTDIR}/${OUT_FASTP}/${OUT_R2}\
-m --merged_out ${OUTDIR}/${OUT_FASTP}/${OUT_MERGED}\
-a auto --poly_g_min_len 10 --poly_x_min_len 10 -q 15 -u 40 -n 5 -l 30 --low_complexity_filter 30\
-c --overlap_len_require 30 --overlap_diff_limit 5 --overlap_diff_percent_limit 20\
-w ${CPU} --verbose --json=${OUTDIR}/${OUT_FASTP}/${FILEBASE}.json --html=${OUTDIR}/${OUT_FASTP}/${FILEBASE}.html
```
4. Quality check on merged reads
```
srun fastqc -q -o ${OUTDIR}/${OUT_FQC}/${POST} -t 3 ${OUTDIR}/${OUT_FASTP}/${OUT_R1} ${OUTDIR}/${OUT_FASTP}/${OUT_R2} ${OUTDIR}/${OUT_FASTP}/${OUT_MERGED}

```
5. Taxonomic classification
    
5.1. Refseq database establishment: [script source](https://github.com/miwipe/KapCopenhagen)
  
5.2. End-to-end alignment using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml): 
```
#== alignment against archaea, fungi, viral refseq (same settings for other indexed refseq DBs)
DB="/PARENT/naleng/bowtie2db/afv/afv"
srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${OUTDIR}/${OUT_ALIGN}/${FILEBASE}/${FILEBASE}.$(basename $DB).sam
#== sam to bam
cd ${OUTDIR}/${OUT_ALIGN}/${FILEBASE}
srun samtools view -bS ${FILEBASE}.$(basename $DB).sam > ${FILEBASE}.$(basename $DB).bam
rm ${FILEBASE}.$(basename $DB).sam
```
5.3 Merge all bam files, convert to *.sam.gz file, and sorting: [script source](https://github.com/miwipe/KapCopenhagen)

5.4. Taxonomic classification

Tool is [ngsLCA](https://github.com/miwipe/ngsLCA).

Citation: [Wang et al., 2022](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.14006)
```
#==taxonomy
NAMES="/PARENT/naleng/T2023_08_db/taxonomy/names.dmp"
NODES="/PARENT/naleng/T2023_08_db/taxonomy/nodes.dmp"
ACC="/PARENT/naleng/T2023_08_db/taxonomy/nucl.accession2taxid.gz"

#== classification
srun ngsLCA -simscorelow 0.95 -simscorehigh 1.0 -names ${NAMES} -nodes ${NODES} -acc2tax ${ACC} -bam ${WORKFOLDER}/${FILEBASE}_L30.sorted.sam.gz -outnames ${WORKFOLDER}/${FILEBASE}.L30.ss095to1
```
Other parameters are not modified.

## Ancient DNA (aDNA) damage pattern analysis (Fig. S3 and Fig. S4)
Tool is [mapDamage2](https://ginolhac.github.io/mapDamage/).

1. Download the Refseq sequnces for 26 major taxa from NCBI
2. Indexing Refseq
```
INDIR="/PARENT/mapdamage/fasta"
cd ${INDIR}
FILEBASE=$(basename $(find . -mindepth 1 -maxdepth 1 -type d | sed -n ${SLURM_ARRAY_TASK_ID}p))
# for *.fasta < 2GB
srun bwa index -a is ${FILEBASE}.fasta
srun samtools faidx ${FILEBASE}.fasta
#== END

INDIR="/PARENT/mapdamage/fna0"
cd ${INDIR}
FILEBASE=$(basename $(find . -mindepth 1 -maxdepth 1 -type d | sed -n ${SLURM_ARRAY_TASK_ID}p))
# for *.fna < 2GB
srun bwa index -a is ${FILEBASE}.fna
srun samtools faidx ${FILEBASE}.fna
#== END

INDIR="/PARENT/mapdamage/fna1"
cd ${INDIR}
FILEBASE=$(basename $(find . -mindepth 1 -maxdepth 1 -type d | sed -n ${SLURM_ARRAY_TASK_ID}p))
# for *.fna > 2GB
srun bwa index -a bwtsw ${FILEBASE}.fna
srun samtools faidx ${FILEBASE}.fna
#== END
```
3. Extract classified read IDs from output of ngsLCA/*.lca using [extract_seqid_from_lca.py](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/scripts/extract_seqid_from_lca.py) and merge them into single txt file using [seqid_merge.R](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/scripts/seqid_merge.R)

4. Subset the QC fastq files: [subset_qc-reads_fastp.sh](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/scripts/subset_qc-reads_fastp.sh)

5. mapDamage: [mapDamage2_L30.sh](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/scripts/mapDamage2_L30.sh)
```
srun mapDamage -i ${DBi}/${FILEBASE}.sort.bam -r ${DBi}/${DB}.fasta --rescale --single-stranded --rescale-out=${DBi}/${FILEBASE}.rescale.bam -d ${DBi}/${FILEBASE}
```
## Normolization prior to Redundancy Analysis (RDA)

Script: [normalization.R](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/scripts/normalization.R). \n

The key functions are rlog and varianceStabilizingTransformation from [DESeq2 package](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). 

Citation: [Love et al., 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)


## Redundancy Analysis (RDA) and variance partitioning analysis (Fig. 2 K and Fig. 3 E-F)
0. Terrestrial vegetation: [rda_venn_terrestrialPlants.R](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/scripts/rda_venn_terrestrialPlants.R)
1. Terrestrial mammalian: [rda_venn_terrestrialMamm.R](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/scripts/rda_venn_terrestrialMamm.R)
2. Aquatic communities: [rda_venn_aquatic.R](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/scripts/rda_venn_aquatic.R)

## Network analysis (Fig. 2 L and Fig. 3 G)
Tool is ecoCopula (https://github.com/gordy2x/ecoCopula). It is detailed explained by [Popovic et al., 2018](https://www.sciencedirect.com/science/article/pii/S0047259X17307522?via%3Dihub) and [2019](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13247)

0. network analysis for terrestrial ecosystem: [network_terrestrial.R](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/scripts/network_terrestrial.R)
1. network analysis for aquatic ecosystem: [network_aquatic.R](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/scripts/network_aquatic.R)

## Past permafrost simulation (Fig. 1 B and Fig. S1)
0. Download the source data: [palaeo proxies-based temperature](https://github.com/StefanKruse/R_PastElevationChange), [BIO1](https://www.worldclim.org/data/worldclim21.html), [present-day permafrost distribution in the Tibetan Plateau](https://tc.copernicus.org/articles/11/2527/2017/), and [SRTM 30 m digital elevation data with SRTM-Downloader plugin in the QGIS software](https://qgis.org/de/site/). 
1. Merging SRTM raster files in the [QGIS software](https://qgis.org/de/site/): Raster > Miscellaneous > Merge
2. Prepare the raster files for permafrost simulation using r-script: [prepare_raster_files.R](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/scripts/prepare_raster_files.R).
3. Permafrost simulation using r-script: [permafrost_predict.R](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/scripts/permafrost_predict.R)





