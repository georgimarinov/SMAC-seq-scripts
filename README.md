# SMAC-seq-scripts

**1. Tombo extraction** 

Run the `TomboSingleReadsExtract-tombo_de_novo.py` (for Tombo versions prior to 1.5) or the `TomboSingleReadsExtract-tombo_de_novo-1.5.py` (for Tombo version 1.5) scripts in order to convert Tombo per_read_stats files into text files. The script has multiple options for different sequence contexts, excluding certain sequence contexts, etc.:

```
python TomboSingleReadsExtract-tombo_de_novo-1.5.py tombo.per_read_stats genome.fa outfile_prefix 
[-m5C-only] [-m6A-only] [-CG-only] [-CG-CG-only] [-GC-only] [-m6A-CG-only] [-m6A-GC-only] [-m6A-GC-CG-only] 
[-doT] [-T-only] [-generic bases(comma-separated)] [-excludeContext string(,string2,string3,...,stringN) radius]
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
python NanoporeTSVMappingStats.py merged.m6A-only.reads.tsv.bgz NanoporeTSVMappingStats-merged.m6A-only &
```

**4. Create coverage file**

Use the `methylation-reads-tsv-to_coverage.py` script to create a coverage file. 

```
python methylation_reads_all.tsv threshold outfile [-stranded +|-] [-minAbsLogLike float] [-minAbsPValue float] [-BayesianIntegration window(bp) step alpha beta pseudosamplesize] [-N6mAweight pseudosamplesize genome.fa] [-saveNewSingleMoleculeFile filename]
```

Example, using 0.5 as the threshold:

```
python methylation-reads-tsv-to_coverage.py merged.m6A-only.reads.tsv.bgz 0.5 merged.m6A-only.cutoff_0.5.coverage &
```
