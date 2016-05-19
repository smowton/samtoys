
#include <htslib/hts.h>
#include <htslib/sam.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv) {

  if(argc < 2) {
    fprintf(stderr, "Usage: filter_match_ratio match_proportion (e.g. 0.5)\n");
    exit(1);
  }

  double required_prop = atof(argv[1]);
  if(required_prop <= 0 || required_prop > 1) {
    fprintf(stderr, "Match proportion must be a real number > 0 and <= 1\n");
    exit(1);
  }

  char filter_message_buf[1024];

  htsFile* hf = hts_open("-", "r");
  if(!hf) {
    fprintf(stderr, "Failed to open stdin\n");
    exit(1);
  }

  bam_hdr_t* header = sam_hdr_read(hf);

  htsFile* hfo = hts_open("-", "wb");
  if(!hfo) {
    fprintf(stderr, "Failed to open stdout\n");
    exit(1);
  }

  sam_hdr_write(hfo, header);

  bam1_t *rec = bam_init1();

  while(sam_read1(hf, header, rec) >= 0) {

    int32_t total_bases = rec->core.l_qseq;
    int32_t matched_bases = 0;
    
    uint32_t* cigar_ops = bam_get_cigar(rec);
    
    for(uint32_t i = 0; i != rec->core.n_cigar; ++i) {

      uint32_t cigar_op = bam_cigar_op(cigar_ops[i]);
      uint32_t cigar_len = bam_cigar_oplen(cigar_ops[i]);

      if(cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL)
	matched_bases += cigar_len;
      else if(cigar_op == BAM_CHARD_CLIP)
	total_bases += cigar_len;

    }

    double match_prop = ((double)matched_bases) / total_bases;
    if(match_prop < required_prop) {

      // Force unmapped:
      rec->core.flag |= BAM_FUNMAP;
      // Note how it got that way:
      sprintf(filter_message_buf, "Filtered by filter_match_ratio (threshold match %g; actual %g)", required_prop, match_prop);
      bam_aux_append(rec, "rf", 'Z', strlen(filter_message_buf) + 1, (uint8_t*)filter_message_buf);

    }

    sam_write1(hfo, header, rec);

  }

  hts_close(hf);
  hts_close(hfo);

  return 0;

}
