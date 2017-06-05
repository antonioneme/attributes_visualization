# attributes_visualization
This program creates a circle plot of the provided genomic attributes, such as (differential) gene expression (RNA-seq), transcription factor strength (ChIP-seq), chromatin  accessibility (FAIRE-seq), and any other transcriptomic or epigenetic attribute. It also allows the display of TADs (topological association domains). It also computes a bunch of statistics, such as average scores of each attribute per window, it applys a basic hotspot detection (anomaly detection). To use it:

python interactions_visualization_CTCF.py -i chrall_TAD.csv -chr all -tss -res 4000 -img visual/CTCF_chrall_4.png -indep tracks_indep_4.csv -motifs tracks_motifs_1.csv -out visual/CTCF_chrall_out_1.csv -outliers visual/outliers_chrall_1 -stats visual/CTCF_chrall_stats_1 -CTCFloc visual/CTCF_chrall_loc_1.csv

There are two levels of resolution: the whole genome (-chr all), or a specific chromosome (-chr chr1 f.e.). The number of pixels (windows) for the visualization is specified by -res. The output file is indicated by -img.

The "independent" tracks are specified in the file following the flag "-indep". It is a tab separated file containing the name and path to the files. The format of the independet track files is:

#TAD	chrall_TAD_gwv.csv	0.06	100,100,200	20	TAD	no	0.1	105,120,190	20.0
TG_prim	TG_primary_gwv.csv	0.2	250,0,0	20.0	TG_prim	no	0.25	30,230,25	20.0
TG_sec	TG_secondary_gwv.csv	0.4	0,0,250	20.0	TG_sec	no	0.45	70,140,55	20.0
VDR_P	VDR_persistent_gwv_2.csv	0.55	255,140,0	30.0	VDR_P	no	0.6	240,23,30	30.0

The first line is ignored since a "#" appears there. Then, a name for each track is provided, followed by the source of that track, and followed by formatting aspects. The format of the tracks is, for example, for TG_primary_gwv.csv, like this:

chr1	3371147	3371375	0.549554364812071
chr1	13910252	13910595	0.574686903057089
chr1	24194774	24194776	0.583020612934639
chr1	27153201	27153630	0.380095938435784
chr1	28199055	28199176	0.799576133078772
chr1	28286504	28286873	0.408865748111838
chr1	28503103	28503455	0.663798269210502
chr1	31882412	31882699	0.299360215495836
chr1	39456916	39457272	0.497658688635961
chr1	41327797	41327799	0.237445078237307

The first three columns specify the loci, the last one is the value of the attribute at that location
All files are tab separated files.

The motifs track is again a file containig the source and formatting information for motifs. A motif file is of the form:

chr1	910163	910493	1.0
chr1	970776	971553	1.0
chr1	973750	974636	0.0
chr1	1080524	1081204	1.0
chr1	1373437	1374156	0.0
chr1	1476856	1477144	1.0
chr1	1481261	1481715	1.0
chr1	1706807	1707447	0.0
chr1	1711189	1712313	1.0
chr1	2165249	2166412	1.0

A 1 in dicates that a motif was found, 0 otherwise. There is a reason for maintaiing motifs in a separate file from the independent attributes, to be discussed later.

You can 
