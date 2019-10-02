import argparse
import sys
import pysam


def avg_quality(bamfile, out, quality):
    nmrtr = {}  # numerator for mean calculation
    denom = {}  # denominator for mean calculation
    bam = pysam.AlignmentFile(bamfile, "rb")
    for pcol in bam.pileup(min_base_quality = quality):
        qual = 10 ** -(np.array(pcol.get_mapping_qualities())/10)
        cov = len(qual)
        if cov < 2:
            continue
        if cov not in nmrtr:
            nmrtr[cov] = sum(qual)
            denom[cov] = len(qual)
        else:
            nmrtr[cov] += sum(qual)
            denom[cov] += len(qual)

    for cov in sorted(nmrtr.keys()):
        out.write("%s\t%s\n" % (cov, nmrtr[cov]/denom[cov]))

    fbuf.close()
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

    avg_quality(pufbuf, outbuf, args.quality)
