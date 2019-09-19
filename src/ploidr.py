import argparse
import numpy as np
import process_bam as pb
import ploidy_model as pm
import plot_read_data as plot


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='subparser')

    process_bam = subparsers.add_parser("process_bam")
    process_bam.add_argument("--bam", required = True)
    process_bam.add_argument("--out", required = True)
    process_bam.add_argument("--bed", default = False)
    process_bam.add_argument("--quality", default = 15, type = int)

    denoise = subparsers.add_parser("denoise")
    denoise.add_argument("--count_file", required = True)
    denoise.add_argument("--out", required = True)
    denoise.add_argument("--iter", type = int, default = 3, required = True)

    histo = subparsers.add_parser("histo")
    histo.add_argument("--count_file", required = True)
    histo.add_argument("--out", required = True)

    assess = subparsers.add_parser("assess")
    assess.add_argument("--count_file", required = True)
    assess.add_argument("--ploidies", nargs='+', type = int, required = True)

    args = parser.parse_args()

    if args.subparser == 'process_bam':
        print("Now processing BAM file %s..." % args.bam)
        if args.bed:
            print("\tprocessing subset found in %s" % args.bed)
        pb.get_biallelic_coverage(args.bam, args.out, args.bed, args.quality)
        print("Success!")
    elif args.subparser == 'denoise':
        print("Denoising count file %s..." % args.count_file)
        np.savetxt(args.out, pb.denoise_reads(args.count_file,
                                                 args.iter),
                   fmt='%d')
        print("Filtered data sucessfully saved in %s" % args.out)
    elif args.subparser == 'histo':
        cnts = np.loadtxt(args.count_file)
        plot.plot_read_hexbin(cnts, args.out + ".hexbin.pdf")
        plot.plot_read_histogram(cnts, args.out + ".histo.pdf")
        print("Files saved in %s.hexbin.pdf and %s.histo.pdf!" % (args.count_file, args.count_file))
    elif args.subparser == 'assess':
        cnts = np.loadtxt(args.count_file)
        uniq = np.unique(cnts[:,1], return_counts = True)
        dens = uniq[1] / np.sum(uniq[1])
        n_val = uniq[0]
        pld = np.array(args.ploidies)
        aicllh = pm.get_Log_Likelihood_AIC(cnts[:,0], pld, dens, n_val)
        print(aicllh[0])
        print(aicllh[1])
        print(pld)
        print("The most likely model is %s-ploid." % pm.select_model(aicllh[1], pld))
