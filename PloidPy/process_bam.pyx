from cpython cimport array
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
import numpy as np
cimport numpy as np
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

cdef float g_div = 0
cdef int g_denom = 0


cdef np.ndarray get_count(AlignmentFile bam, str contig,
                                            int start=-1, end=-1,
                                            int mapq_thresh=15,
                                            int base_qual=13):
    if not ((start == end == -1) or (isinstance(start, int) and
                                       isinstance(end, int))):
        raise ValueError("""Invalid configuration of contig, start, and end
                         parameters""")
    cdef int length
    length = (bam.get_reference_length(contig) if (start == end == -1) else
              end - start)
    # count arrays
    cdef np.ndarray ACGT
    ACGT = np.zeros((4, length))
    cdef float numer
    cdef long denom
    numer = 0
    denom = 0
    cdef AlignedSegment read
    for read in bam.fetch(contig=contig, start=(None if start == -1 else start),
                          end=(None if end == -1 else end)):
        if (read.flag & (0x4 | 0x100 | 0x200 | 0x400)):
            continue
        if read.is_duplicate or read.is_qcfail or read.is_secondary or read.is_unmapped:
            continue
        if read.mapping_quality < mapq_thresh:
            continue

        seq = read.seq

        for r_pos, ref_pos in read.get_aligned_pairs(True):
            if read.query_qualities[r_pos] < base_qual:
                continue
            else:
                numer += 10 ** -(read.query_qualities[r_pos]/10)
                denom += 1
            start = (0 if start == -1 else start)
            end = (length - 1 if end == -1 else end)
            if r_pos is None or ref_pos is None:
                continue
            elif ref_pos > end or ref_pos < start:
                continue

            if seq[r_pos] == 'A':
                ACGT[0][ref_pos - start] += 1
            elif seq[r_pos] == 'C':
                ACGT[1][ref_pos - start] += 1
            elif seq[r_pos] == 'G':
                ACGT[2][ref_pos - start] += 1
            elif seq[r_pos] == 'T':
                ACGT[3][ref_pos - start] += 1
    global g_div
    global g_denom
    g_div = numer/denom
    g_denom = denom
    return ACGT


def get_MAC_TRC(ACGT, f):
    length = len(ACGT[0])
    # number of (i+1) allele sites
    cnt = array.array('l', [0, 0, 0, 0])
    for i in range(length):
        num_allele = np.sum(ACGT[:,i] > 0)
        cnt[num_allele - 1] += 1
        if num_allele != 2:
            continue
        else:
            min_val = float('inf')
            summation = 0
            if ACGT[0][i] != 0:
                summation += ACGT[0][i]
                if min_val > ACGT[0][i]:
                    min_val = ACGT[0][i]
            if ACGT[1][i] != 0:
                summation += ACGT[1][i]
                if min_val > ACGT[1][i]:
                    min_val = ACGT[1][i]
            if ACGT[2][i] != 0:
                summation += ACGT[2][i]
                if min_val > ACGT[2][i]:
                    min_val = ACGT[2][i]
            if ACGT[3][i] != 0:
                summation += ACGT[3][i]
                if min_val > ACGT[3][i]:
                    min_val = ACGT[3][i]
            f.write("%i %i\n" % (min_val, summation))
    return cnt


def get_biallelic_coverage(bamfile, outfile, bed=False, map_quality=15,
                           base_quality=13):
    start = time.time()
    cdef AlignmentFile bam
    bam = pysam.AlignmentFile(bamfile, "rb")
    outf = outfile if outfile[-3:] == ".gz" else outfile + ".gz"  # compressed
    info_file = (outfile[:-3] + ".info" if outfile[-3:] == ".gz"
                 else outfile + ".info")
    # check if file exists
    if path.exists(outf) or path.exists(outfile) or path.exists(info_file):
        raise IOError("Output file %s exists! Will not overwrite." % outfile)
    out = gzip.open(outf, "wt")
    cnt = [0, 0, 0, 0]
    error_list = []
    if bed:
        f = open(bed)
        for line in f:
            br = line.split()
            contig = br[0]
            start = -1
            end = -1
            if len(br) >= 3:
                start = br[1]
                end = br[2]
            ACGT = get_count(bam, contig, start, end, map_quality, base_quality)
            cnt0 = get_MAC_TRC(ACGT, out)
            for i in range(4): cnt[i] += cnt0[i]
            global g_div
            global g_denom
            error_list.append((g_div, g_denom))
        f.close()
    else:
        for contig in bam.header.references:
            ACGT = get_count(bam, contig, -1, -1, map_quality, base_quality)
            cnt0 = get_MAC_TRC(ACGT, out)
            for i in range(4): cnt[i] += cnt0[i]
            global g_div
            global g_denom
            error_list.append((g_div, g_denom))
    numer = 0
    denom = 0
    for etpl in error_list:
        numer += etpl[0] * etpl[1]
        denom += etpl[1]
    info = open(info_file, "w+")
    info.write("p_err\t%8.10f\n" % (numer/denom))
    info.write("1-count\t%i\n" % cnt[0])
    info.write("2-count\t%i\n" % cnt[1])
    info.write("3-count\t%i\n" % cnt[2])
    info.write("4-count\t%i\n" % cnt[3])
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
