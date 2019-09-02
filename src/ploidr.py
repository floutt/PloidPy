import argparse
import numpy as np
import process_bam as pb
import ploidy_model as pm


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='subparser')

    process_bam = subparsers.add_parser("process_bam")
    process_bam.add_argument("--bam", required = True)
    process_bam.add_argument("--out", required = True)
    process_bam.add_argument("--bed", default = False)

    denoise = subparsers.add_parser("denoise")
    denoise.add_argument("--count_file", required = True)
    denoise.add_argument("--out", required = True)
    denoise.add_argument("iter", type = int, default = 3, required = True)

    histo = subparsers.add_parser("histo")
    histo.add_argument("--count_file", required = True)
    histo.add_argument("--out", required = True)

    assess = subparsers.add_parser("assess")
    assess.add_argument("--count_file", required = True)
    assess.add_argument("--ploidies", nargs='+', type = int, required = True)

    args = parser.parse_args()

    if args.subparser == 'process_bam':
        subargs = process_bam.parse_args()
        print("Now processing BAM file %s..." % subargs.bam)
        if subargs.bed:
            print("\tprocessing subset found in %s" % subargs.bed)
        pb.get_biallelic_coverage(subargs.bam, subargs.out, subargs.bed)
        print("Success!")
    elif args.subparser == 'denoise':
        subargs = denoise.parse_args()
        print("Denoising count file %s..." % subargs.count_file)
        np.savetxt(subargs.out, pb.denoise_reads(subargs.count_file,
                                                 subargs.iter),
                   fmt='%d')
        print("Filtered data sucessfully saved in %s" % subargs.out)
    elif args.subparser == 'histo':
        # FILL IN TODO
        None
    elif args.subparser == 'assess':
       # FILL IN TODO
       None
