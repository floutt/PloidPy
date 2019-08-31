import numpy as np
import pysam
import os


def get_biallelic_coverage(bamfile, bed = False, outfile):
    bam = pysam.AlignmentFile(bamfile, "rb")
    f = None
    if bed:
        f = open(bed)
    else:
        # temporary storage file for contig names
        temp = "temp.file"
        # command to store data in temporary file
        bashcmd = ("samtools view -H %s | grep SQ | cut -f2 | cut -c 4- > " +
                   temp)
        os.system(bashcmd % bamfile)
        f = open(temp)

    for line in f:
        bedrow = line.split()
        nuc_cov = None
        if bed:
            nuc_cov = np.arrray(bam.count_coverage(line[0], int(line[1]),
                                                   int(line[2])))
        else:
            nuc_cov = np.arrray(bam.count_coverage(line[0]))

        # only get biallelic site
        nuc_cov = nuc_cov[:, np.sum((nuc_cov == 0) == False, axis = 0) == 2].T
        nuc_cov = nuc_cov[(nuc_cov == 0) == False]
        nuc_cov = np.reshape(nuc_cov, (2, int(len(nuc_cov) / 2))).T

        # get and save minimum and total allele count numbers
        svf = pen(outfile, 'ab')
        np.savetxt(svf, np.array([np.min(nuc_cov, axis = 0),
                                  np.sum(nuc_cov, axis=0)]).T, fmt='%d')
        svf.close()
    f.close()
