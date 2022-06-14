# DIVE

## The algorithm

DIVE is a purely statistical and completely annotation-free algorithm that proposes a new conceptual approach to discovering k-mer sequences associated with high rates of sequence diversification. DIVE is an efficient algorithm designed to identify sequences that may mechanistically cause sequence diversification (e.g., CRISPR repeat or transposon end) and the variable sequences near them, such as an insertion site. The identified sequences are assigned statistical scores for biologists to prioritize them. For more details, see [1].

## Installation

### pip

To install DIVE simply run the following pip command on the terminal:

```bash
pip install biodive
```

### github

To install DIVE directly from the repository simply run the following commands:

```python
git clone https://github.com/jordiabante/biodive.git
cd biodive
conda create -n biodive python=3.6.8
conda activate biodive
pip install -e .
```

## Usage

To run a single-sample analysis

```python
# import bio module
from biodive import bio

# define input file and output dir
outdir = "/path/to/outdir/"
fqfile = "/path/to/fastq.gz"

# configure run
config = bio.Config(
    outdir=outdir,              # directory where output files will be stored
    kmer_size=27,               # k-mer size used in the analysis
    min_smp_sz=5,               # minimum sample size to compute p-value
    max_smp_sz=50,              # maximum number of sequences sampled per
    lmer_size=7,                # l-mer size used to compute jaccard similarity between k-mers
    jsthrsh=0.25                # jaccard similarity threshold used to collapse the observed sequences
    max_fastq_reads=5000000,    # maximum number of FASTQ reads to process
    annot_fasta=[]              # array containing fasta files to use with blast
)

# run analysis
bio.biodive_single_sample_analysis(fqfile,config)
```

If `len(annot_fasta)>0`, then `blast` must be available on the path.

## Output files

### Anchor sequences table

A table with suffix `_anchors.txt.gz` is produced containing information about the interesting anchors detected (keys in old convention). The file contains the following columns:

```bash
    sequence id | assembly of {anchor1,anchor2,...} | max_c_up | max_n_up | max_efct_sz_up | max_efct_sz_qval_up | max_kmer_up | max_c_dn | max_n_dn | max_efct_sz_dn | max_efct_sz_qval_dn | max_kmer_dn | A% | C% | G% | T% | {anchor1,anchor2,...}
```

where `up/dn` indicate the position of the HVR with respect to the anchor and:

* `max_c_*`: number of clusters formed for the maximizing anchor in the set in `*` direction.
* `max_n_*`: corresponding number of target sequences observed.
* `max_efct_sz_*`: corresponding effect size.
* `max_efct_sz_qval_*`: corresponding adjusted p-value.
* `max_kmer_*`: corresponding k-mer sequence.

If the `len(annot_fasta)>0`, then two extra columns will be added to the previous table, for each direction (upstream, downstream) and for each FASTA in `annot_fasta`, containing the lowest e-value and the corresponding hit in the FASTA (sequence in FASTA resulting in lowest e-value), and the output will be stored in a new table with suffix `_anchors_annot.txt.gz` (NA will be assigned when e-value>1). For example, if we pass `annot_fasta=[fasta1]` we will see four extra columns:

```bash
    sequence id | ... | {anchor1,anchor2,...} | best_eval_up_fasta1 | best_hit_up_fasta1 | best_eval_dn_fasta1 | best_hit_dn_fasta1 
```

The intermediate XML files produced by blast are also stored for further analysis.

### Re-running annotation

In some cases we might want to update the set of FASTA files we want to blast the results against. Say, for example, that we want to re-run the annotation with FASTA files `f1.fasta`, `f2.fasta`, and `f3.fasta`, with our output `SRRXYZ_anchors.txt.gz` (note the suffix is `_anchors.txt.gz`). In that case, we can use the following python code:

```python
from biodive import bio

anchorfile = "/path/to/SRRXYZ_anchors.txt.gz"
annot_fasta = ["/path/to/annotations/f1.fa", "/path/to/annotations/f2.fa", "/path/to/annotations/f3.fa"]
config = bio.Config(annot_fasta=annot_fasta)

bio.biodive_single_sample_analysis_annotation(anchorfile,config)
```

### Anchor sequences FASTA

Three FASTA files are produced:

1. FASTA file with suffix `_assemb_anchors.fasta`: assembled anchor sequences.
2. FASTA file with suffix `_max_anchor_up.fasta`: maximizing anchor sequence upstream.
3. FASTA file with suffix `_max_anchor_dn.fasta`: maximizing anchor sequence downstream.

Note that not all anchor sequences in 2 and 3 are necessarily significant.

### Target sequences table

For each anchor in the set `{anchor1,anchor2,...}`, the target sequences are stored in a file with suffix `_targets.txt.gz` containing the following columns:

```bash
    anchor | upstream/downstream | distance | target | number of instances observed
```

## References

[1] J. Abante, P.L. Wang, J. Salzman. *DIVE: a reference-free statistical approach to diversity-generating & mobile genetic element discovery*, bioarxiv (2022).
