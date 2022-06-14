## DBSCAN
import numpy as np
from biodive import mystats

# Similar k-mers 
kmer_1 = 'ATGGACCAGATATAGGGAGAGCCAGGTAGG'
kmer_2 = 'ATGGACCAGAATAGGGAGAGCCACGTAGGG'
# diff    ----------*-------------*--------

# Clearly different k-mer
kmer_3 = 'GCGTATTGGTAGCAGACTAGATTAGAGGAG'
kmer_4 = 'GCGTATTGGATAGCAGACTAGATTAAGGAT'
kmer_5 = 'GCGTATTGGATAGCAGACTAGATTAAGGAT'

kmers = [kmer_1, kmer_2, kmer_3, kmer_4, kmer_5]

mystats.dbscan_stats(kmers)

## Key Oracle

# Normal approximation
import timeit
import numpy as np
from biodive import mystats
from scipy.stats import binom,norm

x = 1
min_smp_sz = 20
n_tot = 1000000
n_so_far = 10000
phat = x/n_so_far

1.0-binom.cdf(min_smp_sz, n_tot, phat)
1.0-norm.cdf((min_smp_sz-n_tot*phat)/np.sqrt(n_tot*phat*(1.0-phat)),0,1)

mystats.kp_kmer_oracle_says_slow(x,n_so_far,n_tot,min_smp_sz)
mystats.kp_kmer_oracle_says(x,n_so_far,n_tot,min_smp_sz)
timeit.timeit('mystats.kp_kmer_oracle_says_slow(x,n_so_far,n_tot,min_smp_sz)', 'from __main__ import mystats,x,n_so_far,n_tot,min_smp_sz', number=10000)
timeit.timeit('mystats.kp_kmer_oracle_says(x,n_so_far,n_tot,min_smp_sz)', 'from __main__ import mystats,x,n_so_far,n_tot,min_smp_sz', number=10000)

## Poibin model
from biodive import mystats, simulations, bio

# generate simulated data
hvr_len = 30
rd_len = 150
lmer_size = 7
num_rds = 500
num_gnm = 1000
gnm_len = 1000
kmer_size = 28
gnms, hvr_rng = simulations.sim_gnms_sngl_species(num_gnm, gnm_len, hvr_len)
gnm_ids, rds, pos = simulations.samp_reads(gnms, rd_len, num_rds, sigma=0.02)
seq_dct = simulations.preproc_rds(rds, pos, hvr_rng, kmer_size, lmer_size, dist=1, jsthrsh=0.25, quiet=True)

# check with oracle
config = bio.Config(min_smp_sz=10)
mystats.check_w_oracle(seq_dct, 500, 500, config)

# train data
null_ind = mystats.poibin_train_data(len(list()), 1000)

# fit Poisson model
n_max = 50
pnvec = mystats.fit_poibin_model(seq_dct, n_max, null_ind)

# compute p-values
n_min = 10
mystats.comp_poibin_pval(seq_dct, pnvec, n_min)

# correct for multiple hypothesis
mystats.poibin_mult_hyp_corr(seq_dct)

# get Poisson binomial p-value
seq_dct = simulations.preproc_rds(rds, pos, hvr_rng, kmer_size, lmer_size, dist=1, jsthrsh=0.25, quiet=True)
mystats.poibin_test(seq_dct, m_min=1000, n_min=10, n_max=50)

# retain only significant ones
mystats.drop_nonsig_anchors(seq_dct, 0.1)

# collapse k-mers (should return two sequences: one upstream one downstream)
kmer_stats = mystats.cllps_anchors(seq_dct)

# get range of positive sequences
rng_fs1 = range(gnms[0].find(kmer_stats[0][0]), gnms[0].find(kmer_stats[0][0])+len(kmer_stats[0][0]))
rng_fs2 = range(gnms[0].find(kmer_stats[1][0]), gnms[0].find(kmer_stats[1][0])+len(kmer_stats[0][0]))
print(f"Range of first FS detected: {rng_fs1}")
print(f"Range of second FS detected: {rng_fs2}")
print(f"Range of HVR: {hvr_rng}")
