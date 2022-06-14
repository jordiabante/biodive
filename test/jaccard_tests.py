## Jaccard similarity
from biodive import jaccard

# Similar k-mers 
kmer_1 = 'ATGGACCAGATATAGGGAGAGCCAGGTAGG'
kmer_2 = 'ATGGACCAGAATATAGGGAGAGCCACGTAG'
# diff    ----------i--------------*-----

# Clearly different k-mer
kmer_3 = 'GCGTATTGGATAGCAGACTAGATTAGAGGA'

assert len(kmer_1)==len(kmer_2)==len(kmer_3)
k = len(kmer_1)

# Create l-mers
lmer_size = 4
lmers_1 = jaccard.get_lmers(kmer_1, lmer_size)
lmers_2 = jaccard.get_lmers(kmer_2, lmer_size)
lmers_3 = jaccard.get_lmers(kmer_3, lmer_size)

# Compute JS
jaccard.comp_js_exact(set(lmers_1), set(lmers_2))
jaccard.comp_js_exact(set(lmers_1), set(lmers_3))
jaccard.comp_js_exact(set(lmers_2), set(lmers_3))

# Create hashes
hash_1 = jaccard.hash_kmers(lmers_1, 0)
hash_2 = jaccard.hash_kmers(lmers_2, 0)
hash_3 = jaccard.hash_kmers(lmers_3, 0)

# Compute exact JS with hashes (it slightly changes b/c of canonical k-mer calculation!)
jaccard.comp_js_exact(hash_1, hash_2)
jaccard.comp_js_exact(hash_1, hash_3)
jaccard.comp_js_exact(hash_2, hash_3)

# Create k-mer sketches
s = 5
mysetup = f"from biodive import jaccard; s = {s}; lmers_1={lmers_1}"
s1 = jaccard.kmer_sketch(s, lmers_1, 0)
s2 = jaccard.kmer_sketch(s, lmers_2, 0)
s3 = jaccard.kmer_sketch(s, lmers_3, 0)

# Compute approximate JS with sketches
jaccard.comp_js(s1, s2)
jaccard.comp_js(s1, s3)
jaccard.comp_js(s2, s3)

# Test addition of new value to dictionary
l = 3
vals = {'ACGTG':1,'AGTTA':2}
new_vals = ['ACGTA', 'CCCCC']

for new_val in new_vals:
    vals = jaccard.add_new_val(vals, new_val, l)

vals