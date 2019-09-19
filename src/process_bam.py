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


def get_biallelic_coverage(bamfile, outfile, bed = False, quality = 15):
    bam = pysam.AlignmentFile(bamfile, "rb")
    f = None
    if bed:
        f = open(bed)
    else:
        # temporary storage file for contig names
        temp = "temp." + str(np.random.randint(100))
        # command to store data in temporary file
        bashcmd = ("samtools view -H %s | grep SQ | cut -f2 | cut -c 4- > " +
                   temp)
        os.system(bashcmd % bamfile)
        f = open(temp)

    # counts for the total number of sites with a certain amount of unique
    # nucleotide variants
    counts = np.array([0, 0, 0, 0])

    # variables used to calculate the mean coverage of the data AS A WHOLE
    for line in f:
        br = line.split()
        nuc_cov = None
        if bed:
            nuc_cov = np.array(bam.count_coverage(br[0], int(br[1]),
                                                  int(br[2]),
                                                  quality_threshold = quality))
        else:
            nuc_cov = np.array(bam.count_coverage(br[0],
                                                  quality_threshold = quality))

        # categorize sites by unique nucleotide count
        pop_sites = np.sum((nuc_cov == 0) == False, axis = 0)
        site_counts = np.unique(pop_sites, return_counts = True)
        counts += site_counts[1][np.where(site_counts[0] > 0)[0]]

        # only get biallelic site
        nuc_cov = nuc_cov[:, pop_sites == 2].T
        nuc_cov = nuc_cov[(nuc_cov == 0) == False]
        nuc_cov = np.reshape(nuc_cov, (int(len(nuc_cov) / 2), 2)).T

        # get and save minimum and total allele count numbers
        svf = open(outfile, 'ab')
        np.savetxt(svf, np.array([np.min(nuc_cov, axis = 0),
                                  np.sum(nuc_cov, axis = 0)]).T, fmt='%d')
        svf.close()
    for i in range(4):
        print("%s %s-allele sites detected" % (counts[i], i + 1))
    f.close()
    # remove temporary file
    os.system("rm %s" % temp)


# denoises read count file generated from get_biallelic_coverage by removing
# removing the putative false positive biallelic sites. This is done by
# comparing the data to a given binomial error model. A normal distribution is
# used to represent the "true" data - not because it is necessarily
# representative
def denoise_reads(readfile, total_mean, p_err):
    # calculates the log likelihood value of an array of likelihood values
    def log_lh(mat):
        return np.sum(np.log(mat))

    x = np.loadtxt(readfile)[:,0]
    error_model = stats.binom(total_mean, p_err)
    em_lh = error_model.pmf(x)
    em_lh[em_lh < EPS] = EPS  # replace 0s with EPS
    # set prior values
    nm_mean = np.mean(np.max(x))
    nm_std = 1
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
    return posterior
