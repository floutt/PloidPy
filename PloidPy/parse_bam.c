#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <htslib/sam.h>

#define ROWSTRIDE 4

void cmap(char *bam, samFile *bamfile, char *chrom, int start, int end,
		  int min_mapq, int min_bseq, double *numer, long *denom,
		  int *count_mat, sam_hdr_t *hdr) {
	int tid = bam_name2id(hdr, chrom);
	bam1_t *read = bam_init1();
	hts_idx_t *idx = sam_index_load(bamfile, bam);
	hts_itr_t *itr = sam_itr_queryi(idx, tid, start, end);

	double n = 0;
	int d = 0;	
	while(sam_itr_next(bamfile, itr, read) > 0) {
		
		int32_t pos = read->core.pos;
		uint32_t len = read->core.l_qseq;
		
		uint8_t *q = bam_get_seq(read);
		uint32_t q2 = read->core.qual;
		if ( min_mapq > q2 ) { continue; }
		if (read->core.flag & (BAM_FUNMAP | BAM_FDUP | BAM_FSECONDARY | BAM_FQCFAIL)) { continue; }
			
		char *qseq = (char *)malloc(len);
		for(int i=0; i< len ; i++){
			qseq[i] = seq_nt16_str[bam_seqi(q,i)];
			if ( (pos + i) > end || (pos + i) < start) { continue; }
			if ( min_bseq > q[i] ) { continue; }
			n = n + pow(10, - q[i] / 10);
			d = d + 1;
			if (qseq[i] == 'A') { count_mat[(((pos+i) - start)* ROWSTRIDE)] += 1; }
			if (qseq[i] == 'C') { count_mat[(((pos+i) - start)* ROWSTRIDE) + 1] += 1; }
			if (qseq[i] == 'G') { count_mat[(((pos+i) - start)* ROWSTRIDE) + 2] += 1; }
			if (qseq[i] == 'T') { count_mat[(((pos+i) - start)* ROWSTRIDE) + 3] += 1; }
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
