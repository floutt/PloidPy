import numpy as np
cimport numpy as np
import scipy.stats as stats
import time
import gzip
from math import ceil
import os.path as path
from process_bam cimport *

# in the case that a value produces a probability of 0 (or a probability lower
# than what can be represented with floating-point numbers, we replace it with
# the lowest representable number for a 64bit operating system. We do this in
# order to prevent any complications arising from NaN values (i.e. log(0) is
# replaced with log(EPS)
EPS = np.finfo(np.float64).tiny


cdef long double numer = 0
cdef unsigned long long denom = 0

# calculates the base-specific coverage of a specific region and stores data in
# a NumPy array
# b - string with location of BAM file
# bamfile - pointer to C htsFile structure
# c - contig/chromosome name
# bamHdr - C header file struct for the corrosponding BAM file
# start - start position in contig (0-index)
# end - end position in contig
# mapq_thresh - minimum mapping quality (Phred scale)
# base_qual - minimum base quality (Phred scale)
cdef np.ndarray get_count(char *b, htsFile *bamfile, char *c,
                          sam_hdr_t *bamHdr, int start, int end,
                          int mapq_thresh, int base_qual):
    cdef int ref_len = get_ref_length(bamHdr, c)
    np_arr = np.zeros((ref_len, 4), dtype=np.dtype("i"))
    cdef int[:,:] cm = np_arr

    cmap(b, bamfile, c, start, end, mapq_thresh, base_qual, &numer, &denom,
         &cm[0,0], bamHdr)
    return np_arr


# writes biallelic MAC and TRC values onto the file object f. Input ACGT coming
# from get_count
# ACGT - base coverage matrix
# f - file object
def get_MAC_TRC(ACGT, f):
    allele_cnt = np.sum(ACGT > 0, axis=1)
    cnt = np.bincount(allele_cnt)[1:]
    if len(cnt) < 4:
        cnt = np.append(cnt, [0] * (4 - len(cnt)))
    two_sites = ACGT[np.sum(ACGT > 0, axis=1) == 2, :]
    bi_mac = np.quantile(two_sites, 0.75, axis=1, interpolation="lower")
    bi_trc = np.sum(two_sites, axis=1)
    [f.write("%i %i\n" % (bi_mac[i], bi_trc[i])) for i in range(cnt[1])]
    return cnt


# writes biallelic MAC and TRC values to a file, option for bed file is given.
# bam - location of bam file
# outfile - output file location
# bed - False if not bed file given, else the variable should be set to the
# location of the bed file
# map_quality - minimum mapping quality
# base_quality - minimum base quality
def get_biallelic_coverage(bam, outfile, bed=False, map_quality=15,
                           base_quality=13, MAX_NP=200000000):
    start = time.time()

    b = bam.encode('UTF-8')
    if not path.exists("%s.bai" % bam):
        print("BAM index not detected. Building index...")
        status = sam_index_build(b, 0)
        if status == 0:
            print("\tSuccess!")
        else:
            raise Exception("An error occurred while creating the index.")

    cdef htsFile *bamfile = hts_open(b, b'r')
    cdef sam_hdr_t *hdr = sam_hdr_read(bamfile)
    cdef int num_ctg = hdr.n_targets

    outf = outfile if outfile[-3:] == ".gz" else outfile + ".gz"  # compressed
    info_file = (outfile[:-3] + ".info" if outfile[-3:] == ".gz"
                 else outfile + ".info")
    # check if file exists
    if path.exists(outf) or path.exists(outfile) or path.exists(info_file):
        raise IOError("Output file %s exists! Will not overwrite." % outfile)
    out = gzip.open(outf, "wt")
    cnt = [0, 0, 0, 0]
    l_nmr = 0
    l_dnm = 0
    if bed:
        f = open(bed)
        for line in f:
            br = line.split()
            contig = br[0].encode('UTF-8')
            ref_len = get_ref_length(hdr, contig)
            start = 0
            end = ref_len
            if len(br) >= 3:
                start = int(br[1])
                end = int(br[2])
            if (end - start) > MAX_NP:
                for idx in range(ceil((end - start) / MAX_NP)):
                    s0 = idx * MAX_NP
                    e0 = min(end - start, (idx + 1) * MAX_NP)
                    ACGT = get_count(b, bamfile, contig, hdr, s0, e0,
                                     map_quality, base_quality)
                    cnt0 = get_MAC_TRC(ACGT, out)
                    for i in range(4): cnt[i] += cnt0[i]
                    l_nmr += numer
                    l_dnm += denom
            else:
                ACGT = get_count(b, bamfile, contig, hdr, start, end,
                                 map_quality, base_quality)
                cnt0 = get_MAC_TRC(ACGT, out)
                for i in range(4): cnt[i] += cnt0[i]
                l_nmr += numer
                l_dnm += denom
        f.close()
    else:
        for i in range(num_ctg):
            contig = sam_hdr_tid2name(hdr, i)
            ref_len = get_ref_length(hdr, contig)
            if ref_len > MAX_NP:
                for idx in range(ceil(ref_len / MAX_NP)):
                    s0 = idx * MAX_NP
                    e0 = min(ref_len, (idx + 1) * MAX_NP)
                    ACGT = get_count(b, bamfile, contig, hdr, s0, e0,
                                     map_quality, base_quality)
                    cnt0 = get_MAC_TRC(ACGT, out)
                    for i in range(4): cnt[i] += cnt0[i]
                    l_nmr += numer
                    l_dnm += denom
            else:
                ACGT = get_count(b, bamfile, contig, hdr, 0, ref_len,
                                 map_quality, base_quality)
                cnt0 = get_MAC_TRC(ACGT, out)
                for i in range(4): cnt[i] += cnt0[i]
                l_nmr += numer
                l_dnm += denom
    info = open(info_file, "w+")
    info.write("p_err\t%8.10f\n" % (l_nmr/l_dnm))
    info.write("1-count\t%i\n" % cnt[0])
    info.write("2-count\t%i\n" % cnt[1])
    info.write("3-count\t%i\n" % cnt[2])
    info.write("4-count\t%i\n" % cnt[3])
    info.close()
    out.close()
    hts_close(bamfile)
    print("Count data stored in %s" % outf)
    print("Secondary information stored in %s" % info_file)
    end = time.time()
    print("Process finished in %s minutes" % ((end - start)/60))


# denoises read count file generated from get_biallelic_coverage by removing
# removing the putative false positive biallelic sites. This is done by
# comparing the data to a given binomial error model. A normal distribution is
# used to represent the "true" data - not because it is necessarily
# representative
def denoise_reads(readfile, post_thresh=1.0):
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
        (1 - (sum(posterior == post_thresh)/original_len)) * 100))
    return posterior, reads[posterior == post_thresh]
