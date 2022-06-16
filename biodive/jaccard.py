from Bio.Seq import Seq
from random import random
from scipy.stats import binom, hypergeom

# returns list of l-mers in k-mer
def get_lmers(kmer, l):
    
    """
    Creates list of l-mers from k-mer.
    """

    # Init
    lmers = []
    n_lmers = len(kmer) - l + 1

    # Loop over k-mer
    for i in range(n_lmers):
        lmer = kmer[i:i + l]
        lmers.append(lmer)

    return lmers

# Computes Jaccard similarity exactly
def comp_js_exact(a, b):

    """
    Computes exact Jaccard similarity between two lists of l-mers.
    
    Arguments:
        - a: set of l-mers A.
        - b: set of l-mers B.

    Returns:
        - exact jaccard similarity

    """

    # Return |A int B| / |A union B|
    return len(a.intersection(b)) / len(a.union(b))

# Computes exact p-value for a given number of matches 
def comp_pval_js_exact(x, y, sx, sy):

    """
    Computes exact p-value under null hypothesis that observed approximated JS can be explained by the 
    variability in JS derived from a random distribution of k-mers. The hypergeometric distribution is 
    used in this case. See Ondov et al 2016 for details.

    Arguments:
        - x: set of unique l-mers in k-mer X
        - y: set of unique l-mers in k-mer Y
        - sx: sketch of k-mer X
        - sy: sketch of k-mer Y

    Returns:
        - p-value

    """

    # Get quantities
    s = len(sx)
    m = len(x.union(y))
    w = len(x.intersection(y))
    z = len(sx.intersection(sy))
        # print(f's={s}')
        # print(f'm={m}')
        # print(f'w={w}')
        # print(f'z={z}')

    # Return pval
    return 1.0-hypergeom.cdf(z, m, w, s) if w>0 else 1.0

# collapse keys based on JS
def add_new_val_js(vals, new_val, l, jsthrsh):

    """
    Adds new value to dictionary. If a key with large enough JS is found, then the count 
    is updated. Otherwise, a new key is created and initialized at 1. Returns true if a 
    new key had to be created.
    """

    # init
    x = False
    merged = False
    key_match = ""
    lmers_val = set(get_lmers(new_val, l))

    # if key is present update and return
    if new_val in vals:
        vals[new_val] += 1
        return x

    # loop over keys
    for key in vals:
        
        # if similarity is good enough add to key
        if comp_js_exact(set(get_lmers(key, l)), lmers_val)>=jsthrsh:
            
            # update entry
            merged = True
            vals[key] += 1
            key_match = key
            
            # no need to check further
            break

    # if merged...
    if merged:

        # coin flip to decide if we change key. note: needs to be different, otherwise we remove the entry
        if new_val!=key_match and random()>=0.5:
            vals[new_val] = vals.pop(key_match)

    else:

        # create key if didn't merge
        x = True
        vals[new_val] = 1

    # return bool
    return x
