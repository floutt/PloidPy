#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <htslib/sam.h>
#include <htslib/hts.h>

#define ROWSTRIDE 4

void cmap(char *bam, samFile *bamfile, char *chrom, int start, int end,
		  int min_mapq, int min_bseq, long double *numer, unsigned long long *denom,
		  int *count_mat, sam_hdr_t *hdr, hts_idx_t *idx) {
	int tid = bam_name2id(hdr, chrom);
	bam1_t *read = bam_init1();
	hts_itr_t *itr = sam_itr_queryi(idx, tid, start, end);

	double n = 0;
	unsigned long long d = 0;	
	while(sam_itr_next(bamfile, itr, read) > 0) {
		
		int32_t pos = read->core.pos;
		int read_pos = 0;
		int ref_pos = 0;
		uint32_t cigar_len = read->core.n_cigar;
		
		uint8_t *q = bam_get_seq(read);
		uint32_t q2 = read->core.qual;
		if ( min_mapq > q2 ) { continue; }
		if (read->core.flag & (BAM_FUNMAP | BAM_FDUP | BAM_FSECONDARY | BAM_FQCFAIL)) { continue; }
		uint32_t *cigar = bam_get_cigar(read);
		int i;
		for(i=0; i < cigar_len ; i++) {
			int op = cigar[i] & BAM_CIGAR_MASK;
			int l = cigar[i] >> BAM_CIGAR_SHIFT;
			if (op  == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
				int j;
				for(j = 0; j < l; j++) {
					int read_p = read_pos;
					int ref_p = ref_pos;
					read_pos += 1;
					ref_pos += 1;
					if ( (pos + ref_p) > end || (pos + ref_p) < start) { continue; }
					if ( min_bseq > q[read_p] ) { continue; }
					char base = (char) seq_nt16_str[bam_seqi(q, read_p)];
					n = n + pow(10, - q[read_p] / 10);
					d = d + 1;
					if (base == 'A') { count_mat[(((pos + ref_p) - start) * ROWSTRIDE)] += 1; }
					if (base == 'C') { count_mat[(((pos + ref_p) - start) * ROWSTRIDE) + 1] += 1; }
					if (base == 'G') { count_mat[(((pos + ref_p) - start) * ROWSTRIDE) + 2] += 1; }
					if (base == 'T') { count_mat[(((pos + ref_p) - start) * ROWSTRIDE) + 3] += 1; }
				}
			}
			else if (op == BAM_CREF_SKIP || op == BAM_CDEL) { ref_pos += l; }
			else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) { read_pos += l; }
		}
	}
	
	*numer = n;
	*denom = d;
	bam_destroy1(read);
	sam_itr_destroy(itr);
}

int get_ref_length(sam_hdr_t *hdr, char *chrom) {
		int tid = bam_name2id(hdr, chrom);
		int len = sam_hdr_tid2len(hdr, tid);
		return len;
}
