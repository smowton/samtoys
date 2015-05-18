
#include <htslib/hts.h>
#include <htslib/sam.h>

#include <string.h>

#include <map>
#include <string>
#include <iostream>
#include <algorithm>

int main(int argc, char** argv) {

  htsFile* hf = hts_open(argv[1], "r");
  if(!hf) {
    fprintf(stderr, "Failed to open %s\n", argv[1]);
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

    // Note l_qname includes the trailing null
    if(rec->core.l_qname > 3) {

      char* qn = bam_get_qname(rec);
      char* qnsuffix = qn + (rec->core.l_qname - 3);
      if(qnsuffix[0] == '/' && isdigit(qnsuffix[1])) {
	
	// Remove trailing chars.
	// Move to (2-before-qname-terminator), from (qname-terminator), length (everything except the qname, plus its null terminator)
	memmove(qn + (rec->core.l_qname - 3), qn + (rec->core.l_qname - 1), rec->l_data - (rec->core.l_qname - 1));
	rec->core.l_qname -= 2;
	rec->l_data -= 2;

      }

    }

    sam_write1(hfo, header, rec);

  }

  hts_close(hf);
  hts_close(hfo);

}
