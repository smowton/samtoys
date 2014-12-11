
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>

#include <vector>
#include <iostream>
#include <string>
#include <cstdint>
#include <random>
#include <functional>

#include <stdlib.h>
#include <unistd.h>
#include <string.h>

int main(int argc, char** argv) {

  if(argc < 2) {
    std::cerr << "Usage: seektest samorbamfile [readahead]\n";
    exit(1);
  }

  htsFile* hfi = hts_open(argv[1], "r");
  if(argc >= 3 && std::string(argv[2]) == "readahead") {

    if(!hfi->is_bin) {
      std::cerr << "Readahead currently only usable on BAM files\n";
      exit(1);
    }

    hts_set_threads(hfi, 2);
    bgzf_set_cache_size(hfi->fp.bgzf, BGZF_MAX_BLOCK_SIZE * 2 * 256);

  }

  bam_hdr_t* header = sam_hdr_read(hfi);
  if(!header) {
    std::cerr << "Failed to read input header\nn";
    exit(1);
  }

  bam1_t rec;
  memset(&rec, 0, sizeof(rec));

  std::vector<std::string> qnames;
  std::vector<int64_t> offsets;

  offsets.push_back(bgzf_tell(hfi->fp.bgzf));

  while(sam_read1(hfi, header, &rec) >= 0) {
   
    qnames.push_back(bam_get_qname(&rec));
    offsets.push_back(bgzf_tell(hfi->fp.bgzf));    
    
  }

  std::cerr << "Read " << qnames.size() << " records\n";

  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution(0, qnames.size() - 1);
  auto pickrec = std::bind(distribution, generator);

  for(int i = 0; i < 10000; ++i) {

    int getrec = pickrec();
    bgzf_seek(hfi->fp.bgzf, offsets[getrec], SEEK_SET);
    sam_read1(hfi, header, &rec);
    std::string qname(bam_get_qname(&rec));
    if(qname != qnames[getrec]) {

      std::cerr << "Test failed at try " << i << ": expected " << qnames[getrec] << " but got " << qname << "\n";
      exit(1);

    }

  }

  std::cerr << "All tests passed\n";

  hts_close(hfi);

}
