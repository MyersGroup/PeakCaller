# ChIPseq Peak Caller

By Simon Myers, Nicolas Altemose & Daniel Wells

This pipeline implements the algorithms for calling peaks from ChIPseq data described in [Altemose et al. eLife 2017](https://elifesciences.org/articles/28383).

## Dependencies
- R (3.5.1)
- R packages "data.table" (1.11.4) and "parallel"
- Bedtools (2.27.1)
- Samtools (1.7)
- Perl (5.26.2)

## Example Usage
First we need to summarise the BAM files as fragment position BED files
The BAM files should already have been QC'd (e.g. duplicates removed)
```{bash}
for sample in Chip, Input
do
(
samtools view -F12 -q1 path_to_bams/${sample}.bam | \
	perl GetFragmentPositions.pl path_to_files/Fragment_Position_${sample} 75
) &
done
```

The method requires 2 replicates per ChIP sample, if you don't have this you can split your sample into two pseudoreplicates
```{bash}
perl MakePseudoreplicates.pl path_to_files/Fragment_Position_Chip.sorted.bed
```

Then we need to fit the parameters of the model to our data
```{bash}
Rscript EstimateConstants.R \
	path_to_files/ \
	hg38.sizes \
	sample_name \
	Fragment_Position_Chip.sorted.bed.PR1 \
	Fragment_Position_Chip.sorted.bed.PR2 \
	Fragment_Position_Input.sorted.bed \
	22
```

Finally we can calculate coverage at each base pair and find the peaks
```{bash}
Rscript DeNovoPeakCalling-SingleBase.R \
	path_to_files/ \
	sample_name \
	Fragment_Position_Chip.sorted.bed.PR1 \
	Fragment_Position_Chip.sorted.bed.PR2 \
	Fragment_Position_Input.sorted.bed \
	Constants.sample_name.tsv \
	0.000001 250 "1,1,1" 1 1 \
	22 \
	hg38.sizes
```

Total runtime ~ 20 mins on 16 core server

This code is ported from the original at [https://github.com/altemose/PRDM9-map](https://github.com/altemose/PRDM9-map)

If you use this program, please cite [Altemose et al. eLife 2017](https://elifesciences.org/articles/28383).

This is free software shared in the hope it may be of use; no warranty is given or implied.