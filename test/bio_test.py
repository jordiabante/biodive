from biodive import bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

## Process a record
ok_record = SeqRecord(
    Seq("ATAGATAGTAGATAGATAGATAGTCATGATTAAAA"),
    id="test_rec_ok",
    name="Test record OK",
    description="Test record for testing functions in bio with no N",
)
not_ok_record = SeqRecord(
    Seq("ATAGATAGTAGATAGATAGATAGTCATGNTTAAAA"),
    id="test_rec_not_ok",
    name="Test record not OK",
    description="Test record for testing functions in bio with N",
)

# Init config
config = bio.Config(kmer_size=5,rand_lkp=True)

# count reads in fastq
from biodive import bio
fastqfile = "./test_data/test_data.fastq.gz"
bio.count_nreads_fastq(fastqfile)

# Pre-process FASTQ file
from biodive import bio
fqfile = "./test_data/test_data.fastq.gz"
config = bio.Config(outdir="./test_data/out/", kmer_size=35, min_smp_sz=10, max_smp_sz=50, lmer_size=7, jsthrsh=0.25)
bio.biodive_single_sample_analysis(fqfile,config)

# Pre-process FASTQ file (with annotation)
from biodive import bio
fqfile = "./test_data/test_data.fastq.gz"
annot_fasta = ["./test_data/test_annot.fa"]
config = bio.Config(outdir="./test_data/out/", kmer_size=35, min_smp_sz=10, max_smp_sz=50, lmer_size=7, jsthrsh=0.25, annot_fasta=annot_fasta)
bio.biodive_single_sample_analysis(fqfile,config)

# add annotation to processed file
from biodive import bio
anchorfile = "./test_data/out/test_data_anchors.txt.gz"
annot_fasta = ["./test_data/test_annot.fa"]
config = bio.Config(annot_fasta=annot_fasta)
bio.biodive_single_sample_analysis_annotation(anchorfile,config)
