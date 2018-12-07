# ChIPseq Peak Caller

By Simon Myers, Nicolas Altemose & Daniel Wells

This pipeline implements the algorithms for calling peaks from ChIPseq data described in [Altemose et al. eLife 2017](https://elifesciences.org/articles/28383).

## Requirements
- R (3.5.1)
- R packages "data.table" (1.11.4) and "parallel"
- Bedtools (>2.27.1)
- Perl (5.26.2)

## Example Usage
First we need to convert the BAM files to fragment position BED files
The BAM files should already have been QC'd (e.g. duplicates removed)
```{bash}
for sample in Chip, Input
do
(
samtools view -F12 -q1 path_to_bams/${sample}.bam | \
	perl ../GetFragmentDepth.pl \
		Fragment_Depth_${sample} \
		Fragment_Position_${sample} \
		75 \
		10000 \
		path/to/bedtools2/bin/ \
		../hg38.sizes
) &
done
```

The method requires 2 replicates per ChIP sample, if you don't have this you can split your sample into two pseudoreplicates
```{bash}
perl ../MakePseudoreplicates.pl \
	Fragment_Position_Chip.sorted.bed \
	0.5 \
	../hg38.sizes \
	path/to/bedtools2/bin/
```

Then we need to fit the parameters of the model to our data
```{bash}
Rscript EstimateConstants.R \
	path_to_files/ \
	hg38.sizes \
	sample_name \
	Fragment_Position_Chip.sorted.bed.PR1.sorted.bed \
	Fragment_Position_Chip.sorted.bed.PR2.sorted.bed \
	Fragment_Position_Input.sorted.bed \
	22
```

Finally we can calculate coverage at each base pair and find the peaks
```{bash}
Rscript DeNovoPeakCalling-SingleBase.R \
	path_to_files/ \
	sample_name \
	Fragment_Position_Chip.sorted.bed.PR1.sorted.bed \
	Fragment_Position_Chip.sorted.bed.PR2.sorted.bed \
	Fragment_Position_Input.sorted.bed \
	Constants.sample_name.tsv \
	0.000001 250 "1,1,1" 1 1 \
	22 \
	hg38.sizes
```

This code is ported from the original at [https://github.com/altemose/PRDM9-map](https://github.com/altemose/PRDM9-map)

If you use this program, please cite [Altemose et al. eLife 2017](https://elifesciences.org/articles/28383).

This is free software shared in the hope it may be of use; no warranty is given or implied.