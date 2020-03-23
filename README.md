# SMAC-seq-scripts

**1. Tombo extraction** 

Run the `TomboSingleReadsExtract-tombo_de_novo.py` (for Tombo versions prior to 1.5) or the `TomboSingleReadsExtract-tombo_de_novo-1.5.py` (for Tombo version 1.5) scripts in order to convert Tombo `per_read_stats` files into text files. The script has multiple options for different sequence contexts, excluding certain sequence contexts, etc.:

```
python TomboSingleReadsExtract-tombo_de_novo-1.5.py tombo.per_read_stats genome.fa outfile_prefix 
[-m5C-only] [-m6A-only] [-CG-only] [-CG-CG-only] [-GC-only] [-m6A-CG-only] [-m6A-GC-only] 
[-m6A-GC-CG-only] [-doT] [-T-only] [-generic bases(comma-separated)]
[-excludeContext string(,string2,string3,...,stringN) radius]
[-excludeChr chr1[,chr2,...,chrN]] [-chrPrefix string]
```

Example for `A` positions:

```
python TomboSingleReadsExtract-tombo_de_novo.py 0.tombo.per_read_stats genome.fa 0.tombo.m6A-only -m6A-only
```

Run the script for each individual `tombo.per_read_stats` file.

**2. Merging and indexing** 

Merge  the converted files into a single file, sorted by coordinates:

```
cat *.m6A-only.reads.tsv | sort -k1,1 -k2,2n -k3,3n | bgzip > merged.m6A-only.reads.tsv.bgz &
```

`tabix`-index the file:

```
tabix -s 1 -b 2 -e 3 merged.m6A-only.reads.tsv.bgz
```

**3. Calculate mapping statistics**

```
python NanoporeTSVMappingStats.py merged.m6A-only.reads.tsv.bgz NanoporeTSVMappingStats-merged.m6A-only
```

**4. Create coverage file**

Use the `methylation-reads-tsv-to_coverage.py` script to create a coverage file. 

```
python methylation_reads_all.tsv threshold outfile 
[-stranded +|-] [-minAbsLogLike float] [-minAbsPValue float]
[-BayesianIntegration window(bp) step alpha beta pseudosamplesize] 
[-N6mAweight pseudosamplesize genome.fa] [-saveNewSingleMoleculeFile filename]
```

Example, using 0.5 as the threshold:

```
python methylation-reads-tsv-to_coverage.py merged.m6A-only.reads.tsv.bgz 0.5 merged.m6A-only.cutoff_0.5.coverage
```

Convert to a `.bgz` file:

```
cat merged.m6A-only.cutoff_0.5.coverage | bgzip > merged.m6A-only.cutoff_0.5.coverage.bgz
```

Then `tabix`-index:

```
tabix -s 1 -b 2 -e 3 merged.m6A-only.cutoff_0.5.coverage.bgz
```

**5. Bayseian integration**

The Bayseian integration calculation is also carried out using the `methylation-reads-tsv-to_coverage.py` script. For efficiency of calculation, compute it on the individual converted tombo files, as follows (for a 10-bp context and (10,10) prior):

```
python methylation-reads-tsv-to_coverage.py 0.tombo.m6A-only 0.5 
0.tombo.m6A-only.all0.min_p_val0.4.cutoff_0.5.coverage.BI_w10_a10_b10 
-minAbsPValue 0.4 -BayesianIntegration 10 1 10 10 50 
-saveNewSingleMoleculeFile 0.tombo.m6A-only.BI_w10_a10_b10.reads.tsv
```

Merge the files:

```
cat *tombo.m6A-only.BI_w10_a10_b10.reads.tsv | sort -k1,1 -k2,2n -k3,3n | 
bgzip > merged.BI_w10_a10_b10.reads.tsv.bgz
```

Then `tabix`-index:

```
tabix -s 1 -b 2 -e 3 merged.BI_w10_a10_b10.reads.tsv.bgz
```

**6. Filtering fully methylated reads**

This operation can be done using the `filterFullyMethylatedReads.py` script:

```
python filterFullyMethylatedReads.py methylation_reads_all.tsv WindowSize minFraction 
[-keepShort] [-missingBasesFilter genome.fa basecontexts(comma-separated) minFraction [-doMBFSet]]
```

**7. Create genome browser tracks**

```
python coverage_to_wig.py coverage.bgz window step chrField MfieldID UfieldID chrom.sizes outprefix [-minCov N_reads]
```

Where the `M` and the `U` fields indicate the column IDs of the numbers of methyllated and unmethylated reads, respectively, and the `window` and `step` parameters specify the width and the stride of the averaging. 

This script will output two `bedGraph` files -- a `coverage.wig` one (which contains the number of reads covering a position) and a `meth.wig` one (which contains the fraction of methylated reads). These can then be converted into `bigWig` files that can in turn be displayed on a genome browser.

**8. Making metaplots around a position**

The `coverage.bgz` can be used to make metaplots around a set of positions, as follows, with a variety of parameters (window size, minimal coverage, etc.):

```
python signalAroundPeaks-nano.py inputfilename chrFieldID posField strandField radius window coverage.bgz outputfilename 
[-bismark.cov] [-bed] [-minCov N] [-unstranded] [-ERANGE_hts] [-narrowPeak] [-first number]
        Input format: <fields .. tabs> chr <tab> position <tab> strandField
        This script outputs the average signal over all regions within the given radius
        if the the -bed option is used, the middle point of a bed region will be used; specifiy the posField as the left coordinate of the region
        if the the -narrowPeak option is used, the posField will be ignored and strand will be assumed to be +
        Note: the script will normalize only against the number of windows that have a CpG ot other signal in the methylation file
        use the [-bismark.cov] option if you want to use the script for a Bismark bedcov output file
```

**9. Making single molecule plots**

Single-molecule plots can be generated over a list of regions (one plot per region will be generated) using the `SMAC-footprints-from-methylation-reads-tsv-tabix.py` and `SMAC-footprints-from-methylation-reads-tsv-tabix-kmeans.py` scripts. The first script will apply hierarchical clustering while the second one will use k-means. The commands are otherwise the same. There is a wide variety of options regarding display, subsampling, etc.:

```
python methylation_reads_all.tsv peak_list chrFieldID leftFieldID rightFieldID strandFieldID tabix_path outfile_prefix 
[-resize factor] [-subset N] [-label fieldID] [-minCov fraction]
[-minPassingBases fraction] [-minReads N] [-unstranded] [-minAbsLogLike float]
[-scatterPlot colorscheme minScore maxScore color|none] [-window bp] [-readStrand +|-] 
[-printMatrix] [-deleteMatrix] [-binarize threshold]' 
        Use the [-subset] option if you want only N of the fragments; the script will pick the N fragments best covering each region, and will discard regions with fewer than N covering fragments
        Use the [-label] option if you want regions to be labeled with something other than their coordinates'
        The [-heatmap] option will generate png heatmaps instead of text file matrices'
        The [-minCov] option will remove all fragments that cover the region at less than the specified fraction'
```

Example:

```
python SMAC-footprints-from-methylation-reads-tsv-tabix-kmeans.py 
2019_01_16_60min_Diamide-rep2.all.BI_w10_a10_b10.reads.filtered_1kb_0.75.tsv.bgz 
AAD6.TSS-600bp.bed 0 1 2 3 tabix AAD6.TSS-600bp.binary-0.5-gist_heat.2019_01_16_60min_Diamide-rep2.BI.filt.10bp.resize0.5
-window 10 -minCov 1 -deleteMatrix -binarize 0.5 -scatterPlot gist_heat 0 1.1 w -resize 0.5 -unstranded
```

There are  also analogous scripts, `SMAC-footprints-from-methylation-reads-tsv-tabix-all-sites.py`  and `SMAC-footprints-from-methylation-reads-tsv-tabix-kmeans-all-sites.py` that will create single-molecule plots combining reads covering multiple regions. 

**10. Calculating single-molecule correlations**

To estimate coaccessibility between all pairs of regions (with sufficient coverage), use the `SingleMoleculeCorrelation-empirical-quantiles.py` script:

```
python SingleMoleculeCorrelation-empirical-quantiles.py methylation_reads_all.tsv peaks chrFieldID leftFiled RightFieldID
minCoverage maxDist N_samplings tabix_location outfile [-subsample N] [-quantiles N]
```

Example:

```
python SingleMoleculeCorrelation-empirical-quantiles.py 
20180515_Yeast_Run-tombo_denovo_1.3.reads.filtered_1kb_0.75.tsv.bgz
Saccharomyces_cerevisiae.SacCer_Apr2011.20.TSS-100bp.bed 0 1 2 100 20000 1000 tabix 
SMCorrEQ.20180515_Yeast_Run-tombo_denovo_1.3.reads.filtered_1kb_0.75.TSS-100bp.q5.ss50 
-quantiles 5 -subsample 50
```

**11. Calculating NMI matrices**

To calculate NMI matrices, the `SingleMoleculeCorrelation-NMI-matrix.py` script can be used. 

```
python SingleMoleculeCorrelation-NMI-matrix.py methylation_reads_all.tsv region.bed
chrFieldID leftField rightFieldID minCoverage windowsize stepsize tabix_location outfileprefix 
[-subsample N] [-expectedMaxDist bp] [-label fieldID]
```

Example:

```
python SingleMoleculeCorrelation-NMI-matrix.py 
2018_07_05_Diamide_0min.all.BI_w10_a10_b10.reads.filtered_1kb_0.75.tsv.bgz 
CTT1.TSS-600bp.bed 0 1 2 50 1 1200 tabix 
NMI.min50cov.1bp.TIF-seq-updated.CTT1.TSS-600bp.2018-07-05_Diamide_0min -expectedMaxDist 1500
```

If running genome-wide, split the genome into overlapping bins for paralellization efficiency, e.g. 50-kbp in size with a 10-kbp stride, and calculate a separate matrix for each, then take the average NMI values for each pair of coordinates for downstream analyses. 

