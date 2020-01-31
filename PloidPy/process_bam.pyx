from cpython cimport array
from pysam.libcalignmentfile cimport AlignmentFile
import numpy as np
import pysam
import scipy.stats as stats
import time
import gzip
import os.path as path

# in the case that a value produces a probability of 0 (or a probability lower
# than what can be represented with floating-point numbers, we replace it with
# the lowest representable number for a 64bit operating system. We do this in
# order to prevent any complications arising from NaN values (i.e. log(0) is
# replaced with log(EPS)
EPS = np.finfo(np.float64).tiny


def get_biallelic_coverage(bamfile, outfile, bed=False, map_quality=15):
    start = time.time()
    cdef AlignmentFile bam
    bam = pysam.AlignmentFile(bamfile, "rb")
    outf = outfile if outfile[-3:] == ".gz" else outfile + ".gz"  # compressed
    info_file = (outfile[:-3] + ".info" if outfile[-3:] == ".gz"
                 else outfile + ".info")
    # check if file exists
    if path.exists(outf) or path.exists(outfile) or path.exists(info_file):
        IOError("Output file %s exists! Will not overwrite." % outfile)
    out = gzip.open(outf, "wt")
    qual_num = 0.0
    qual_dnm = 0
    nuc_map = {"a": 0, "t": 1, "g": 2, "c": 3, "A": 0, "T": 1, "G": 2, "C": 3}
    cdef array.array allele_num = array.array('i', [0, 0, 0, 0])

    def pcol_iter(pcol):
        cdef array.array ATGC = array.array('i', [0, 0, 0, 0])
        bases = pcol.get_query_sequences()
        cdef i
        for i in range(len(bases)):
            b = bases[i]
            if b == '' or b == 'N' or b == 'n':
                continue
            ATGC[nuc_map[b]] += 1
        cdef a_num, min_cnt, max_cnt
        a_num = 0
        min_cnt = 0
        max_cnt = 0
        first = True
        for i in range(len(ATGC)):
            cnt = ATGC[i]
            if cnt == 0:
                continue
            elif first:
                first = False
                min_cnt = cnt
            if cnt > 0:
                a_num += 1
            if cnt > max_cnt:
                max_cnt = ATGC[i]
            elif cnt < min_cnt:
                min_cnt = ATGC[i]

        allele_num[a_num - 1] += 1
        if not a_num == 2:
            return False
        else:
            out.write("%i %i\n" % (min_cnt, min_cnt + max_cnt))
            return True

    if bed:
        f = open(bed)
        for line in f:
            br = line.split()
            for pcol in bam.pileup(br[0], int(br[1]), int(br[2]),
                                   min_mapping_quality=map_quality):
                if pcol_iter(pcol):
                    qual_num += sum(map(lambda x: 10 ** -(x/10),
                                        pcol.get_query_qualities()))
                    qual_dnm += pcol.get_num_aligned()
    else:
        for pcol in bam.pileup(min_mapping_quality=map_quality):
            if pcol_iter(pcol):
                qual_num += sum(map(lambda x: 10 ** -(x/10),
                                    pcol.get_query_qualities()))
                qual_dnm += pcol.get_num_aligned()
    info = open(info_file, "w+")
    info.write("p_err\t%s\n" % (qual_num/qual_dnm))
    info.write("1-count\t%s\n" % allele_num[0])
    info.write("2-count\t%s\n" % allele_num[1])
    info.write("3-count\t%s\n" % allele_num[2])
    info.write("4-count\t%s\n" % allele_num[3])
    info.close()
    out.close()
    print("Count data stored in %s" % outf)
    print("Secondary information stored in %s" % info_file)
    end = time.time()
    print("Process finished in %s minutes" % ((end - start)/60))


# denoises read count file generated from get_biallelic_coverage by removing
# removing the putative false positive biallelic sites. This is done by
# comparing the data to a given binomial error model. A normal distribution is
# used to represent the "true" data - not because it is necessarily
# representative
def denoise_reads(readfile):
    # calculates the log likelihood value of an array of likelihood values
    def log_lh(mat):
        return np.sum(np.log(mat))

    info_file = (readfile[:-3] + ".info" if readfile[-3:] == ".gz"
                 else readfile + ".info")
    info = open(info_file)
    p_err = float(info.readline().strip().split("\t")[1])
    info.close()
    reads = np.loadtxt(readfile)
    original_len = len(reads)
    x = reads[:, 0]
    error_model = stats.binom(np.mean(reads[:, 1]), p_err)
    em_lh = error_model.pmf(x)
    em_lh[em_lh < EPS] = EPS  # replace 0s with EPS
    # set prior values
    nm_mean = np.mean(x[np.random.randint(len(x), size=100)])
    nm_std = np.std(x[np.random.randint(len(x), size=100)])
    # initial old value assumes 0 probability
    nm_lh_old = np.ones_like(x) * EPS
    nm_lh = stats.norm.pdf(x, nm_mean, nm_std)
    nm_lh[nm_lh < EPS] = EPS  # replace 0s with EPS
    posterior = None

    iters = 0
    # run till it converges upon a maximum likelihood value
    while log_lh(nm_lh_old) < log_lh(nm_lh):
        print(log_lh(nm_lh))
        iters += 1
        posterior = nm_lh / (nm_lh + em_lh)
        nm_mean = np.dot(posterior, x) / np.sum(posterior)
        nm_std = np.sqrt(np.dot(posterior, (x - nm_mean) ** 2) /
                         np.sum(posterior))
        nm_lh_old = nm_lh
        nm_lh = stats.norm.pdf(x, nm_mean, nm_std)
        nm_lh[nm_lh < EPS] = EPS  # replace 0s with EPS
    print(log_lh(nm_lh))

    print("Performed %s iterations" % iters)
    print(nm_mean, nm_std)
    print("")
    print("%s percent of data removed" % (
        (1 - (sum(posterior == 1)/original_len)) * 100))
    return posterior, reads[posterior == 1]
