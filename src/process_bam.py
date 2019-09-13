import numpy as np
import pysam
import os


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
        nuc_cov = np.reshape(nuc_cov, (2, int(len(nuc_cov) / 2))).T

        # get and save minimum and total allele count numbers
        svf = open(outfile, 'ab')
        np.savetxt(svf, np.array([np.min(nuc_cov, axis = 1),
                                  np.sum(nuc_cov, axis = 1)]).T, fmt='%d')
        svf.close(i)
    for i in range(4):
        print("%s %s-allele sites detected" % (count[i], i + 1))
    f.close()
    # remove temporary file
    os.system("rm %s" % temp)


# denoises read count file generated from get_biallelic_coverage by removing
# clustering the data into two parts. One which should capture the low coverage
# minor allele reads which often plague the data
def denoise_reads(readfile, iter = 3):
    reads = np.loadtxt(readfile)
    for i in range(iter):
        km = kmeans(reads, 2)
        centroids = km[0]
        degen = np.argmin(centroids[:,0])
        reads = reads[(vq(reads, centroids)[0] == degen) == False]
    return reads
