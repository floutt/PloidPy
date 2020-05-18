#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <htslib/sam.h>


void cmap(char *bam, samFile *fp_in, char *chrom, int start, int end,
				int min_mapq, int min_bseq, double *numer,
				long *denom, int *count_mat, sam_hdr_t *bamHdr);

int get_ref_length(sam_hdr_t *bamHdr, char *chrom);
