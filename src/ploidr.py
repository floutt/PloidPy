import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='subparser')

    process_bam = subparsers.add_parser("process_bam")
    process_bam.add_argument("--bam", required = True)
    process_bam.add_argument("--out", required = True)


    denoise = subparsers.add_parser("denoise")
    denoise.add_argument("--count_file", required = True)
    denoise.add_argument("--out")
    denoise.add_argument("iter", type = int)

    histo = subparsers.add_parser("histo")
    histo.add_argument("--count_file", required = True)
    histo.add_argument("--out", required = True)

    assess = subparsers.add_parser("assess")
    assess.add_argument("--count_file", required = True)
    assess.add_argument("--ploidies", nargs='+', type = int, required = True)

    args = parser.parse_args()
