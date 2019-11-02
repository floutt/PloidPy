import numpy as np
import pysam
import os
import scipy.stats as stats


# in the case that a value produces a probability of 0 (or a probability lower
# than what can be represented with floating-point numbers, we replace it with
# the lowest representable number for a 64bit operating system. We do this in
# order to prevent any complications arising from NaN values (i.e. log(0) is
# replaced with log(EPS)
EPS = np.finfo(np.float64).tiny


def get_biallelic_coverage(bamfile, outfile, bed = False, map_quality = 15):
    bam = pysam.AlignmentFile(bamfile, "rb")
    out = open(outfile, "w+")
    qual_num = 0
    qual_dnm = 0
    nuc_map = {"a":0, "t":1, "g":2, "c":3}

    def pcol_iter(pcol):
        ATGC = np.zeros(4)
        bases = pcol.get_query_sequences()
        for b in bases:
            if b == '': continue
            ATGC[nuc_map[b.lower()]] += 1
        ATGC = ATGC[ATGC > 0]
        if not len(ATGC) == 2:
            return
        out.write("%d %d\n" % (np.min(ATGC), np.sum(ATGC)))

    if bed:
        f = open(bed)
        for line in f:
            br = line.split()
            for pcol in bam.pileup(br[0], int(br[1]), int(br[2]), min_base_quality = map_quality):
                pcol_iter(pcol)
                qual_num += np.sum(10 ** -(np.array(pcol.get_query_qualities())/10))
                qual_dnm += pcol.get_num_aligned()
    else:
        for pcol in bam.pileup(min_base_quality = map_quality):
            pcol_iter(pcol)
            qual_num += np.sum(10 ** -(np.array(pcol.get_query_qualities())/10))
            qual_dnm += pcol.get_num_aligned()
    print("Base quality error probability:")
    print(qual_num/qual_dnm)
    out.close()


# denoises read count file generated from get_biallelic_coverage by removing
# removing the putative false positive biallelic sites. This is done by
# comparing the data to a given binomial error model. A normal distribution is
# used to represent the "true" data - not because it is necessarily
# representative
def denoise_reads(readfile, p_err):
    # calculates the log likelihood value of an array of likelihood values
    def log_lh(mat):
        return np.sum(np.log(mat))

    reads = np.loadtxt(readfile)

    x = reads[:,0]
    error_model = stats.binom(np.mean(reads[:,1]), p_err)
    em_lh = error_model.pmf(x)
    em_lh[em_lh < EPS] = EPS  # replace 0s with EPS
    # set prior values
    nm_mean = np.mean(x[np.random.randint(len(x), size = 100)])
    nm_std = np.std(x[np.random.randint(len(x), size = 100)])
    nm_lh_old = np.ones_like(x) * EPS  # initial old value assumes 0 probability
    nm_lh = stats.norm.pdf(x, nm_mean, nm_std)
    nm_lh[nm_lh <  EPS] = EPS  # replace 0s with EPS
    posterior = None

    iters = 0
    # run till it converges upon a maximum likelihood value
    while log_lh(nm_lh_old) < log_lh(nm_lh):
        print(log_lh(nm_lh))
        iters += 1
        posterior = nm_lh / (nm_lh + em_lh)
        nm_mean = np.dot(posterior, x) / np.sum(posterior)
        nm_std = np.sqrt(np.dot(posterior, (x - nm_mean) ** 2) / np.sum(posterior))
        nm_lh_old = nm_lh
        nm_lh = stats.norm.pdf(x, nm_mean, nm_std)
        nm_lh[nm_lh <  EPS] = EPS  # replace 0s with EPS
    print(log_lh(nm_lh))

    print("Performed %s iterations" % iters)
    print(nm_mean, nm_std)
    return posterior, reads[posterior == 1]
