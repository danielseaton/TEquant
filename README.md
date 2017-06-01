# TEquant
Tools for manipulating quantification of transposable element reads from RNA-seq data.

## Workflow

* Map reads to reference

* Create a GTF
  * Curated TE gtf files are available from TEToolkit, here: http://labshare.cshl.edu/shares/mhammelllab/www-data/TEToolkit/ (github page: https://github.com/mhammell-laboratory/tetoolkit)
  * Can combine transcript and repeatmasker GTF files, filtering out repeats that lie in annotated transcripts

* Create a table relating different levels of the TE heirarchy
  * Sample code for this (that works for the hg19_rmsk_TE.gtf file) is in [parse_TE_heirarchies.py](./bin/parse_TE_heirarchies.py).

* Run [featurecounts](http://bioinf.wehi.edu.au/featureCounts/) on each sample
  * Settings depend on library type (e.g. stranded, unstranded)
  * Can ignore or keep multimapping reads
  * Note that featurecounts will by default throw away reads that overlap multiple features by default

* Merge featurecounts outputs across samples
  * Using [merge_featurecounts_lowmem.py](./bin/merge_featurecounts_lowmem.py)
  * This generates a file with counts across all samples, and some metadata (e.g. total number of mapped reads in each sample), indicated by a '_' prefix (e.g. '_assigned' for reads assigned to a feature - see featurecounts documentation for details).

* Subselect and sum data across transcript types:
  * Using [split_and_sum_TE_counts](./bin/split_and_sum_TE_counts.py)
  * Creates a count matrix for each type of transcript:
    * gene_id
    * TEclass_id (most general)
    * TEfamily_id
    * TEgene_id
    * TEtranscript_id (most specific - individual loci)
