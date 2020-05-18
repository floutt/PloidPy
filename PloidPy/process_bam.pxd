cdef extern from "parse_bam.h":
    void cmap(char *bam, htsFile *bamfile, char *chrom, int start, int end,
              int min_mapq, int min_bseq, double *numer,
              long *denom, int *count_mat, sam_hdr_t *hdr)
    int get_ref_length(sam_hdr_t *hdr, char *chrom)

cdef extern from "htslib/hts.h":
    ctypedef struct htsFile:
        pass
    htsFile *hts_open(const char *fn, const char *mode)
    int hts_close(htsFile *fp)

cdef extern from "htslib/sam.h":
    ctypedef struct sam_hdr_t:
        int n_targets
    sam_hdr_t *sam_hdr_read(htsFile *fp)
    char *sam_hdr_tid2name(const sam_hdr_t *h, int tid)
    int sam_index_build(const char *fn, int min_shift)
