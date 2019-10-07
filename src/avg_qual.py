import argparse
import sys
import pysam
import numpy as np


def avg_quality(bamfile, out, quality):
    nmrtr = 0  # numerator for mean calculation
    denom = 0  # denominator for mean calculation
    bam = pysam.AlignmentFile(bamfile, "rb")
    for pcol in bam.pileup(min_base_quality = quality):
        qual = 10 ** -(np.array(pcol.get_query_qualities())/10)
        cov = len(qual)
        if cov < 2:
            continue
        nmrtr += sum(qual)
        denom += len(qual)

    out.write("%s\n" % nmrtr/denom)
    out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam")
    parser.add_argument("--out")
    parser.add_argument("--quality", type = int, default = 15)
    args = parser.parse_args()

    outbuf = None
    if args.out:
        outbuf = open(args.out, "w+")
    else:
        outbuf = sys.stdout

    avg_quality(args.bam, outbuf, args.quality)
