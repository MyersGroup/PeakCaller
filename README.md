[![DOI](https://img.shields.io/badge/DOI-10.7554%2FeLife.28383-brightgreen)](https://doi.org/10.7554/eLife.28383)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3783600.svg)](https://doi.org/10.5281/zenodo.3783600)

# ChIPseq Peak Caller

By Simon Myers, Nicolas Altemose & Daniel Wells

This pipeline implements the algorithms for calling peaks from ChIPseq data described in [Altemose et al. eLife 2017](https://elifesciences.org/articles/28383).

## Example Usage
### Convert BAM to fragment BED
First we need to summarise the PE BAM files as BED files of physical DNA fragments.
The BAM files should already have been QC'd (e.g. low quality & duplicates removed).
```{bash}
for sample in Chip Input
do
(
bedtools bamtobed -i ${sample}.bam -bedpe | \
  awk '{print $1, $2, $6}' OFS="\t" | \
  LC_ALL=C sort -k1,1V -k2,2n -S5G --parallel=5 \
  > Fragment_Position_${sample}.sorted.bed
) &
done
```

### (Create Pseudoreplicates)
The method requires 2 replicates per ChIP sample, if you don't have this you can split your sample into two pseudoreplicates.
```{bash}
awk '{print > (rand()<0.5 ? (FILENAME".PR1") : (FILENAME".PR2"))}' Fragment_Position_Chip.sorted.bed
```

### Call Peaks
Then we can run the script to fit the parameters of the model to our data and then calculate coverage at each base pair and find the peaks.
```{bash}
sh MAPeakCaller.sh \
	--intermediates intermediates_folder/ \
	--outconstfile Constants_SampleName1.tsv \
	--outpeakfile SingleBasePairPeaks_SampleName1.bed \
	--chrsizes hg38.sizes \
	-a Fragment_Position_Chip.sorted.bed.PR1 \
	-b Fragment_Position_Chip.sorted.bed.PR2 \
	-i Fragment_Position_Input.sorted.bed \
	--autosomes 22 \
	--pthresh 0.000001 \
	--peakminsep 250
```

This will produce two files:
 - Constants_SampleName1.tsv
 - SingleBasePairPeaks_SampleName1.bed

Total runtime ~ 20 mins on 16 core server.


### Force Calling
We can also 'force call' p-values and enrichments at abitrary locations,
In the example below we're force calling Sample2 at the peak 
postions from Sample1 (you can use any BED file for --forcepositions):

```
sh MAPeakCaller.sh \
	--intermediates intermediates_folder/ \
	--chrsizes hg38.sizes \
	--outconstfile Constants_SampleName2.tsv \
	--outpeakfile ForceCalledPeaks_SampleName2.bed \
	-a Fragment_Position_Chip.sorted.bed.PR1 \
	-b Fragment_Position_Chip.sorted.bed.PR2 \
	-i Fragment_Position_Input.sorted.bed \
	--autosomes 22 \
	--forcepositions SingleBasePairPeaks_SampleName1.bed

```

## Dependencies
- R (3.5.1)
- R packages "data.table" (1.11.4), "IRanges", and "parallel"
- Bedtools (2.27.1)
- Samtools (1.7)

This code is ported from the original at [https://github.com/altemose/PRDM9-map](https://github.com/altemose/PRDM9-map)

If you use this program, please cite [Altemose et al. eLife 2017](https://elifesciences.org/articles/28383).

This is free software shared in the hope it may be of use; no warranty is given or implied.
