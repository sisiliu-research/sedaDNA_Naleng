# sedaDNA analysis of Lake Naleng
This repository contains source codes used to analysis shotgun sequencing/metagenomic data of Liu et al., 2023. [Tibetan terrestrial and aquatic ecosystems collapsed with cryosphere loss inferred from sedimentary ancient metagenomics] 

The script for bioinformatics is done by Lars Harms (lars.harms@awi.de) and modified by Sisi Liu (sisi.liu@awi.de/sisi.liu.research@gmail.com). Other scripts are done by Sisi Liu.

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
-a auto -g --poly_g_min_len 10 -q 15 -u 40 -n 5 -l 15 --low_complexity_filter\
-c --overlap_len_require 30 --overlap_diff_limit 5 --overlap_diff_percent_limit 20\
-w ${CPU} --verbose --json=${OUTDIR}/${OUT_FASTP}/${FILEBASE}.json --html=${OUTDIR}/${OUT_FASTP}/${FILEBASE}.html
```
4. Quality check on merged reads
```
srun fastqc -q -o ${OUTDIR}/${OUT_FQC}/${POST} -t 3 ${OUTDIR}/${OUT_FASTP}/${OUT_R1} ${OUTDIR}/${OUT_FASTP}/${OUT_R2} ${OUTDIR}/${OUT_FASTP}/${OUT_MERGED}

```
5. Taxonomic classification 
Tools: [HOPS-malt](https://github.com/rhuebler/HOPS) and [MEGAN software](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/megan6/) 

5.1. Classification using Megan Alignment Tool (malt) with default setting: 
- top alignments = 100 (unique accession hits)
- min score = 50
- max expected (e-value) = 1
- minimal percent identity = 90 
- m = BlastN (model for nucleotide)
- lowest common ancestor (LCA) = naïve
- top percent = 1 (taking 1% of the top hits into consideration for LCA assignment) 
- minimal support read = 1 (taking a taxon having one read at least for LCA assignment)
```
srun hops -Xmx800G -input ${INPUT}/${INRMA} -output ${EXOUT}/${FILE_BASE} -m malt -c ${CONFIG};
```
5.2. Re-classification in MEGAN software with strict parameters (open the *.rma6 generated by HOPS-malt in the MEGAN software): 
- max expected (e-value) = 1E-5
- minimal percent identity = 95
- top percent = 20
- minimal support read = 3
- lowest common ancestor (LCA) = weighted
- Percent to cover = 80

Other parameters are not modified.

## Ancient DNA (aDNA) damage pattern analysis (Fig. S3 and Fig. S4)
Tool is [HOPS-MaltExtract](https://github.com/rhuebler/MaltExtract). It is detailed explained by [Hübler et al., 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1903-0)

```
srun hops -Xmx60G -input ${INPUT}/${INMEGANRMA} -output ${EXOUT}/${FILE_BASE} -m me_po -c ${CONFIG};
```
Eight main taxa with sufficient read counts are authenticated. They represent terrestrial mammals (Bos mutus), terrestrial plants (Saxifraga sinomontana, Asteroideae, Asteraceae, Salix, and Saliceae), aquatic plants (Potamogeton perfoliatus), and aquatic microbes (N. limnetica). The rates of C>T substitutions are assessed to determine whether the age-dependent signal and climate-dependent are present by integrating the deamination results of these taxa showed the damage patterns across most time slides spanning 17 to 3 ka. The related r-script is [01_MEGAN_id95_aDNA_damage.R](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/R/01_MEGAN_id95_aDNA_damage.R).

## Compositional analysis (Fig. 2 A-E and Fig. 3 A-D)
0. After data cleaning and filtering, two terrestrial datasets (plants and mammals) and one aquatic dataset are rarefied by considering the taxa coverage. Then, the common taxa (have the maximum relative rarefied abundance >= 1% and occur in 5 samples at least) are selected for ordination and network analysis. The r-script is: [02_MEGAN_id95_rarefied.R](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/R/02_MEGAN_id95_rarefied.R).


## Ordination analysis (Fig. 2 K and Fig. 3 E-F)
0. Redundancy Analysis (RDA) and variance partitioning analysis for terrestrial vegetation: [03_MEGAN_id95_rda_venn_terrestrialPlants.R](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/R/03_MEGAN_id95_rda_venn_terrestrialPlants.R)
1. Redundancy Analysis (RDA) and variance partitioning analysis for terrestrial mammalian: [04_MEGAN_id95_rda_venn_terrestrialMamm.R](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/R/04_MEGAN_id95_rda_venn_terrestrialMamm.R)
2. Redundancy Analysis (RDA) and variance partitioning analysis for aquatic communities: [05_MEGAN_id95_rda_venn_aquatic.R](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/R/05_MEGAN_id95_rda_venn_aquatic.R)

## Network analysis (Fig. 2 L and Fig. 3 G)
Tool is ecoCopula (https://github.com/gordy2x/ecoCopula). It is detailed explained by [Popovic et al., 2018](https://www.sciencedirect.com/science/article/pii/S0047259X17307522?via%3Dihub) and [2019](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13247)

0. network analysis for terrestrial ecosystem: [06_MEGAN_id95_network_terrestrial.R](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/R/06_MEGAN_id95_Network_terrestrial.R)
1. network analysis for aquatic ecosystem: [07_MEGAN_id95_network_aquatic.R](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/R/07_MEGAN_id95_Network_aquatic.R)

## Past permafrost simulation (Fig. 1 B and Fig. S1)
0. Download the source data: [palaeo proxies-based temperature](https://github.com/StefanKruse/R_PastElevationChange), [BIO1](https://www.worldclim.org/data/worldclim21.html), [present-day permafrost distribution in the Tibetan Plateau](https://tc.copernicus.org/articles/11/2527/2017/), and [SRTM 30 m digital elevation data with SRTM-Downloader plugin in the QGIS software](https://qgis.org/de/site/). 
1. Merging SRTM raster files in the QGIS software: Raster > Miscellaneous > Merge
2. Prepare the raster files for permafrost simulation using r-script: [08_prepare_raster_files.R](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/R/08_prepare_raster_files.R).
3. Permafrost simulation using r-script: [09_permafrost_predict.R](https://github.com/sisiliu-research/sedaDNA_Naleng/blob/master/R/09_permafrost_predict.R)





