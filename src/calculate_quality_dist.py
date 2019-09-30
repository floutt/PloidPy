import argparse
import sys
import os


def avg_quality(fbuf, out):
    nmrtr = {}  # numerator for mean calculation
    denom = {}  # denominator for mean calculation
    for line in fbuf:
        lsplt = line.split("\t")
        cov = int(lsplt[3])
        if cov < 2:
            continue
        qstr = list(filter(lambda x: not (x in [">", "<", "^", "~", "*"]), lsplt[5]))
        qual = list(map(lambda x: 10 ** (-(ord(x))/10), qstr))
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
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--bam")
    group.add_argument("--pileup")
    parser.add_argument("--out")
    parser.add_argument("--quality", type = int, default = 15)
    args = parser.parse_args()

    tmp = "tmp/mpu"
    outbuf = None
    pufbuf = None

    if args.out:
        outbuf = open(args.out, "w+")
    else:
        outbuf = sys.stdout

    if args.pileup:
        pufbuf = open(args.pileup)
    elif args.bam:
        os.system("mkdir -p tmp")
        os.system("samtools mpileup -Q %s -q 0 -A %s > %s" % (args.quality, args.bam, tmp))
        pufbuf = open(tmp)
    else:
        pufbuf = sys.stdin

    avg_quality(pufbuf, outbuf)

    if args.bam:
        os.system("rm %s" % tmp)
