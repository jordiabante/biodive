from biodive import io

# class for configuration
class Config:

    """
    Class for configuration object.
    """

    def __init__(self, kmer_size=25, rand_lkp=False, dist=0, outdir="./", min_smp_sz=20, max_smp_sz=100, no_new_min_p=0.90,
                 lmer_size=7, jsthrsh=0.25, q_thresh=0.05, batch_sz_poibin=10000, max_fastq_reads=0, annot_fasta=[]):
        
        self.dist = dist                                    # fixed lookup distance (only matters if rand_lkp is False)
        self.outdir = outdir                                # output directory
        self.jsthrsh = jsthrsh                              # jaccard similarity threshold for value collapsing
        self.q_thresh = q_thresh                            # q-value threshold
        self.rand_lkp = rand_lkp                            # true if we randomize lookup distance
        self.lmer_size = lmer_size                          # l-mer size (k-mer is divided into shorter l-mers of size l)
        self.kmer_size = kmer_size                          # k-mer analysis size
        self.max_smp_sz = max_smp_sz                        # maximum sample size required for testing
        self.min_smp_sz = min_smp_sz                        # minimum sample size required for testing
        self.annot_fasta = annot_fasta                      # array containing fasta files to use for annotation
        self.no_new_min_p = no_new_min_p                    # minimum prob that we will have observed at least 1 instance by now
        self.batch_sz_poibin = batch_sz_poibin              # batch size for Poisson binomial model
        self.max_fastq_reads = max_fastq_reads              # maximum number of fastq records to be processed (all if 0)

    def report(self):

        io.print_mess("******************* CONFIGURATION *******************")
        io.print_mess(f"Output directory: {self.outdir}")
        io.print_mess(f"Maximum reads to be processed (0 means no limit): {self.max_fastq_reads}")
        io.print_mess(f"Min probability of observing at least 1 count by now: {self.no_new_min_p}")
        io.print_mess(f"Lookahead distance: {self.dist}")
        io.print_mess(f"k-mer size: {self.kmer_size}")
        io.print_mess(f"l-mer size: {self.lmer_size}")
        io.print_mess(f"Minimum sample size: {self.min_smp_sz}")
        io.print_mess(f"Maximum sample size: {self.max_smp_sz}")
        io.print_mess(f"Jaccard similarity threshold: {self.jsthrsh}")
        io.print_mess(f"Batch size Poisson binomial null: {self.batch_sz_poibin}")
        io.print_mess(f"Q-value threshold: {self.q_thresh}")
        io.print_mess("********************* ANALYSIS **********************")
