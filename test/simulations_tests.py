from biodive import simulations, io

plt_bool = False

# generate random genomes
num_gnm = 5
hvr_len = 20
gnm_len = 200
gnms = simulations.sim_gnm(num_gnm, gnm_len, hvr_len)

# check sizes are right
[len(x) for x in gnms]

# sample reads from it
rd_len = 112
num_rds = 100
gnm_ids, rds, pos = simulations.samp_reads(gnms, rd_len, num_rds, sigma=0.0)

# get sequence dictionary
seq_dct = simulations.preproc_rds(rds, 28, 7, False)

# store reads in fastq format
from biodive import simulations, io
hvr_len = 30
rd_len = 150
num_rds = 500
num_gnm = 1000
gnm_len = 1000
outfile = "./test_data/test_data.fastq.gz"
gnms, hvr_rng = simulations.sim_gnms_sngl_species(num_gnm, gnm_len, hvr_len)
gnm_ids, rds, pos = simulations.samp_reads(gnms, rd_len, num_rds, sigma=0.02)
io.write_fastq(rds, outfile)
