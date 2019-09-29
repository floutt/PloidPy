import argparse


def avg_quality(fbuf, outfile):
    out = open(outfile, "w+")
    nmrtr = {}  # numerator for mean calculation
    denom = {}  # denominator for mean calculation
    for line in fbuf:
        lsplt = line.split("\t")
        cov = int(lsplt[3])
        if cov == 0:
            continue
        qstr = list(filter(lambda x: not (x in [">", "<", "^", "~", "*"]), lsplt[5]))
        qual = list(map(lambda x: 10 ** (-(ord(x))/10), qstr))  # phred 33 offset
        if len(qual[qual > 1]) > 0:
            print("ERROR")
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
