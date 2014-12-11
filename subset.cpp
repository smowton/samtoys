#include <iostream>
#include <fstream>
#include <unordered_set>

#include <stdlib.h>
#include <string.h>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>

int main(int argc, char** argv) {

  // Filter SAM/BAM file on stdin, producing an uncompressed BAM on stdout featuring only those QNAMEs in argv[1]

  if(argc < 2) {
    std::cerr << "Usage: subset filterfile [thread_count] <samorbam >bam\n";
    exit(1);
  }

  std::unordered_set<std::string> keep_qnames;

  std::cerr << "Reading Qnames to keep...\n";

  std::ifstream qname_file(argv[1]);
  std::string qname;

  while(std::getline(qname_file, qname)) {

    size_t end = qname.find_last_not_of(" \n\r\t");
    if(end == std::string::npos)
      continue;
    qname.erase(end + 1);
    qname.erase(0, qname.find_first_not_of(" \n\r\t"));

    keep_qnames.insert(qname);

  }

  std::cerr << "Read " << keep_qnames.size() << " Qnames\n";

  htsFile* hfi = hts_open("-", "r");
  if(argc >= 3) {

    int nthreads = strtol(argv[2], 0, 0);
    if(!hfi->is_bin)
      std::cerr << "Thread count ignored (non-BAM input)\n";
    else {
      hts_set_threads(hfi, nthreads);
      bgzf_set_cache_size(hfi->fp.bgzf, BGZF_MAX_BLOCK_SIZE * nthreads * 256);
    }
    
  }

  htsFile* hfo = hts_open("-", "wb0");

  bam_hdr_t* header = sam_hdr_read(hfi);
  if(!header) {
    std::cerr << "Failed to read SAM/BAM header from stdin\n";
    exit(1);
  }

  if(sam_hdr_write(hfo, header)) {
    std::cerr << "Failed to write SAM/BAM header to stdout\n";
    exit(1);
  }

  unsigned long total = 0, kept = 0;

  bam1_t rec;
  memset(&rec, 0, sizeof(rec));

  while(sam_read1(hfi, header, &rec) >= 0) {

    ++total;

    qname = bam_get_qname(&rec);
    if(!keep_qnames.count(qname))
      continue;

    ++kept;

    if(sam_write1(hfo, header, &rec) < 0) {
      std::cerr << "Failed to write BAM record\n";
      exit(1);
    }

  }

  hts_close(hfi);
  hts_close(hfo);

  std::cerr << "Kept " << kept << " of " << total << " records\n";
  return 0;

}
